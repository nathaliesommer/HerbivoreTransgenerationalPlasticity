# Sommer et al 
# Update title


# Behavior ----

# Packages 
library(dplyr)
library(lme4)
library(lubridate)
library(ggplot2)
library(tidyr)
library(ggforce)
library(emmeans)
library(ks)
library(gghalves)
library(DHARMa)

### Data Processing ----
Behavior_2023 <- read.csv("Data/Behavior_Raw.csv")
# Raw data measured in grid coordinates, converting to cm
Behavior_2023$Y <- Behavior_2023$Y*4
Behavior_2023$X <- Behavior_2023$X*4
Behavior_2023$Distance <- Behavior_2023$Distance*4

# Adjusted start time
Behavior_2023$Time <- hms(Behavior_2023$Time)
behavior_data <- subset(Behavior_2023,Time>(hms("8:40:00")))

behavior_temp <- read.csv("Data/Behavior_Temps_2023_F1_F2_Raw.csv")

behavior_temp <- behavior_temp %>%
  mutate(DateTime = mdy_hm(Date.Time..EDT.))

behavior_data <- behavior_data %>%
  mutate(DateTime = mdy(Date) + hours(hour(Time)) + minutes(minute(Time)))

# Function to find the nearest temperature
find_nearest_temp <- function(obs_time, temp_data) {
  temp_data %>%
    filter(abs(difftime(DateTime, obs_time, units = "mins")) == min(abs(difftime(DateTime, obs_time, units = "mins")))) %>%
    slice(1) %>%
    pull(Ch..1...Temperature.....C.)
}

# Join temperature to behavior_data
behavior_data <- behavior_data %>%
  rowwise() %>%
  mutate(Nearest_Temperature = find_nearest_temp(DateTime, behavior_temp))

# Add a unique id
behavior_data <- behavior_data %>%
  mutate(id = factor(paste(Population, Predation.Treatment, Generation, Day, Terraria, sep = "_")))

# Update any generation label
behavior_data <- behavior_data %>%
  mutate(Generation = case_when(
    Generation == "F1" ~ "G1",
    Generation == "F2" ~ "G2",
    TRUE ~ as.character(Generation)
  ))

# Create the combined treatment column
behavior_data <- behavior_data %>%
  mutate(Treatment = paste(Generation, Predation.Treatment, sep = "_")) %>%
  mutate(Treatment = factor(Treatment, 
                            levels = c("G1_Herbivore", "G2_Herbivore", 
                                       "G1_Predator", "G2_Predator")))

# Create 'type' column, sensu Baker 2022 data

behavior_data <- behavior_data %>%
  mutate(Type = case_when(
    Population %in% c("UP", "MC", "FN", "HF") ~ "Behavior",
    Population %in% c("DC", "SP", "SC", "YF") ~ "Physiology"))

### Analysis ----

# Check population effect
height_model <- lmer(Y ~ Treatment * Type + (1|Population/id), 
                     data = behavior_data)

VarCorr(height_model)

vif_values <- vif(height_model)
print(vif_values)

plot(height_model) 
simulation_output <- simulateResiduals(fittedModel = height_model)
plot(simulation_output)
testDispersion(simulation_output)
testZeroInflation(simulation_output) 
summary(height_model)



# population and individual explains a significant chunk of variation
# however, given unbalanced sample size, we might not want to model this explicitly.
# --> use population-centered heights for kernel density and BA estimates 

# we also need to use a boundary-corrected kernel density approach;
# limits estimates >0 (i.e., above ground)


# centered population heights
behavior_data <- behavior_data %>%
  group_by(Population) %>%
  mutate(Y_centered = scale(Y, scale = FALSE),
         Y_centered = Y_centered - min(Y_centered)) %>%  # make all values positive
  ungroup()

# data validation
print(range(behavior_data$Y_centered, na.rm = TRUE)) # range of Y_centered
print(sum(!is.finite(behavior_data$Y_centered))) # check for infinite values
print(sum(is.na(behavior_data$Y_centered))) # check for NA

# diagnostics
print(summary(behavior_data$Y_centered))
print(summary(behavior_data$Y_centered[behavior_data$Generation == "G1"]))
print(summary(behavior_data$Y_centered[behavior_data$Generation == "G2"]))

# Check the input data
print("Input data summary:")
print(summary(behavior_data$Y_centered))
print(table(behavior_data$Treatment))

#' Calculate density estimates and ranges for behavioral data
#' @param data Dataframe containing Y_centered values and grouping variables
#' @param min_obs Minimum observations required per group (default 10)
#' @return Dataframe with core and broad ranges for each group
calculate_density_ranges <- function(data, min_obs = 10) {
  data %>%
    group_by(Treatment, Type) %>%
    filter(n() >= min_obs) %>%  # Ensure enough points for density estimation
    summarize(
      # Core range (50%)
      avg_core_low = quantile(Y_centered, 0.25),
      avg_core_high = quantile(Y_centered, 0.75),
      # Broad range (95%)
      avg_broad_low = quantile(Y_centered, 0.025),
      avg_broad_high = quantile(Y_centered, 0.975),
      # Mean for positioning
      Y_centered = mean(Y_centered),
      .groups = 'drop'
    )
}

# Calculate density ranges for all data
density_estimate_ranges <- calculate_density_ranges(behavior_data)

# Calculate separate ranges for behavior and physiology types if needed
density_estimate_ranges_behavior <- density_estimate_ranges %>%
  filter(Type == "Behavior")

density_estimate_ranges_physiology <- density_estimate_ranges %>%
  filter(Type == "Physiology")

# Remove these now-redundant sections:
# - The original density_estimate_ranges calculation
# - The calculate_isopleths function and its usage
# - The separate density calculations for behavior and physiology types




#### Bhattacharyya's affinity ----

#' Calculate Bhattacharyya's affinity between two density distributions
#' @param density1,density2 Density objects to compare
#' @param range Range to calculate BA over
#' @return BA value
calculate_bhattacharyya_affinity <- function(density1, density2, range) {
  # Restrict the density estimates to the specified range
  indices1 <- which(density1$eval.points >= range[1] & density1$eval.points <= range[2])
  indices2 <- which(density2$eval.points >= range[1] & density2$eval.points <= range[2])
  
  # Ensure the restricted grids are the same
  if (length(indices1) != length(indices2) || 
      !all(density1$eval.points[indices1] == density2$eval.points[indices2])) {
    stop("Density grids do not match within the specified range.")
  }
  
  # Normalize the density estimates within the specified range
  density1_norm <- density1$estimate[indices1] / sum(density1$estimate[indices1])
  density2_norm <- density2$estimate[indices2] / sum(density2$estimate[indices2])
  
  # Calculate BA
  sqrt_density1 <- sqrt(density1_norm)
  sqrt_density2 <- sqrt(density2_norm)
  affinity <- sum(sqrt_density1 * sqrt_density2) * 
    prod(diff(density1$eval.points[indices1]))
  
  return(affinity)
}

#' Calculate densities for a specific treatment group
#' @param data Behavior data
#' @param treatment Treatment level
#' @return Density object
calculate_treatment_density <- function(data, treatment) {
  y_min <- min(data$Y_centered, na.rm = TRUE)
  y_max <- max(data$Y_centered, na.rm = TRUE)
  
  density_obj <- stats::density(
    data$Y_centered[data$Treatment == treatment],
    from = y_min, 
    to = y_max, 
    n = 100
  )
  
  list(eval.points = density_obj$x, estimate = density_obj$y)
}

# Visualization Functions ----

# Define custom colors for consistent plotting
custom_colors <- c(
  "#FFD700", "#FFF68F",  # Yellows for G1 Herbivore
  "#228B22", "#90EE90",  # Greens for G2 Herbivore
  "#4169E1", "#87CEEB",  # Blues for G1 Predator
  "#9370DB", "#E6E6FA"   # Purples for G2 Predator
)


#' Create base treatment plot with consistent styling
#' @param data Data for plotting
#' @param title Plot title
#' @return ggplot object
create_base_treatment_plot <- function(data, title = "") {
  ggplot(data, aes(x = Treatment, y = Y_centered, fill = Treatment)) +
    scale_fill_manual(values = c(
      custom_colors[1],  # G1 Herbivore
      custom_colors[3],  # G2 Herbivore
      custom_colors[5],  # G1 Predator
      custom_colors[7]   # G2 Predator
    )) +
    labs(title = title,
         x = "Treatment",
         y = "Population-Centered Canopy Height") +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          text = element_text(size = 12),
          legend.position = "none")
}

# Keep the existing violin plot but refactor it to use the base plot function
ggplot(behavior_data, aes(x = Treatment, y = Y_centered, fill = Treatment)) +
  geom_violin(alpha = 0.5) +
  geom_point(data = behavior_data %>% 
               group_by(Treatment) %>% 
               summarize(mean_y = mean(Y_centered)),
             aes(y = mean_y),
             size = 3, color = "black") +
  geom_text(data = behavior_data %>% 
              group_by(Treatment) %>% 
              summarize(mean_y = mean(Y_centered)),
            aes(y = mean_y, label = round(mean_y, 1)),
            vjust = -1.5, size = 4) +
  scale_fill_manual(values = c(
    custom_colors[1],  # G1 Herbivore
    custom_colors[3],  # G2 Herbivore
    custom_colors[5],  # G1 Predator
    custom_colors[7]   # G2 Predator
  )) +
  labs(title = "Height Distribution by Treatment",
       subtitle = "Population-centered mean height",
       x = "Treatment",
       y = "Population-Centered Canopy Height") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 12),
        legend.position = "none")

# Keep the existing temperature plot but refactor for consistency
ggplot(behavior_data, aes(x = Nearest_Temperature, y = Y, color = Treatment)) +
  geom_jitter(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c(
    custom_colors[1],  # G1 Herbivore
    custom_colors[3],  # G2 Herbivore
    custom_colors[5],  # G1 Predator
    custom_colors[7]   # G2 Predator
  )) +
  labs(title = "Relationship between Height and Temperature by Generation",
       x = "Temperature (°C)", 
       y = "Height (cm)") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        text = element_text(size = 12),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom") +
  facet_wrap(~Population)

#' Create type comparison plot without faceting
create_type_comparison_plot <- function(data, density_ranges, treatment, colors, type = NULL) {
  treatment_data <- data %>% filter(Treatment == treatment)
  treatment_ranges <- density_ranges %>% filter(Treatment == treatment)
  
  # Create generation subscript label with treatment type
  gen_num <- substr(treatment, 2, 2)  # Extract number after G
  treatment_type <- ifelse(grepl("Herbivore", treatment), "Herbivore", "Predator")
  x_label <- paste0("G[", gen_num, "]~", treatment_type)
  
  ggplot(treatment_data, aes(x = Treatment, y = Y_centered)) +
    # Half violin plot on the right side
    geom_half_violin(
      aes(fill = Treatment),
      side = "r", 
      alpha = 0.3,
      position = position_nudge(x = 0.3),
      trim = TRUE
    ) +
    # Core range ellipse (darker)
    geom_ellipse(
      data = treatment_ranges,
      aes(x0 = 1,
          y0 = (avg_core_low + avg_core_high) / 2,
          a = 0.2, 
          b = (avg_core_high - avg_core_low) / 2,
          angle = 0),
      fill = colors[1], 
      color = "black", 
      alpha = 0.7
    ) +
    # Broad range ellipse (lighter)
    geom_ellipse(
      data = treatment_ranges,
      aes(x0 = 1,
          y0 = (avg_broad_low + avg_broad_high) / 2,
          a = 0.2, 
          b = (avg_broad_high - avg_broad_low) / 2,
          angle = 0),
      fill = colors[1], 
      color = "black", 
      alpha = 0.3
    ) +
    # Mean point
    geom_point(
      data = treatment_ranges,
      aes(x = 1, y = Y_centered),
      size = 3,
      color = "black"
    ) +
    # Mean labels
    geom_text(
      data = treatment_ranges,
      aes(x = 1, y = Y_centered, 
          label = round(Y_centered, 1)),
      vjust = -1,
      size = 4
    ) +
    # Styling
    scale_y_continuous(limits = c(0, 80)) +
    scale_fill_manual(values = colors[1]) +
    theme_light() +
    theme(
      strip.background = element_rect(fill = "grey80"),
      strip.text = element_text(color = "black", size = 12),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_blank(),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_line(color = "grey95"),
      legend.position = "none"
    ) +
    labs(
      y = "Population-centered canopy height",
      x = NULL
    ) +
    scale_x_discrete(labels = parse(text = x_label))  # Use parsed text for subscript
}

# Function to filter data for specific plot types
filter_behavior_data <- function(data, type = NULL) {
  if (is.null(type)) {
    return(data)  # Return all data for aggregate plot
  } else {
    return(data %>% filter(Type == type))
  }
}

# Create three sets of plots: Behavior only, Physiology only, and Aggregate
create_plot_set <- function(data, type = NULL) {
  # Filter data
  filtered_data <- filter_behavior_data(data, type)
  
  # Calculate density ranges for filtered data
  density_ranges <- calculate_density_ranges(filtered_data)
  
  # Create plots
  plot_set <- list(
    G1H = create_type_comparison_plot(filtered_data, density_ranges, "G1_Herbivore", c(custom_colors[1], custom_colors[2]), type),
    G1P = create_type_comparison_plot(filtered_data, density_ranges, "G1_Predator", c(custom_colors[5], custom_colors[6]), type),
    G2H = create_type_comparison_plot(filtered_data, density_ranges, "G2_Herbivore", c(custom_colors[3], custom_colors[4]), type),
    G2P = create_type_comparison_plot(filtered_data, density_ranges, "G2_Predator", c(custom_colors[7], custom_colors[8]), type)
  )
  
  # Combine plots
  combined <- (plot_set$G1H + plot_set$G1P) / 
              (plot_set$G2H + plot_set$G2P)
  
  # Add title based on type
  title <- case_when(
    is.null(type) ~ "Aggregate - All Populations",
    type == "Behavior" ~ "Behaviorally-plastic Populations",
    type == "Physiology" ~ "Physiologically-plastic Populations"
  )
  
  combined + plot_annotation(title = title)
}

# Create the three plot variations
behavior_plot <- create_plot_set(behavior_data, "Behavior")
physiology_plot <- create_plot_set(behavior_data, "Physiology")
aggregate_plot <- create_plot_set(behavior_data)

# Save plots
#ggsave("behavior_plot.pdf", behavior_plot, width = 12, height = 8)
#ggsave("physiology_plot.pdf", physiology_plot, width = 12, height = 8)
#ggsave("aggregate_plot.pdf", aggregate_plot, width = 12, height = 8)

# Display plots
behavior_plot
physiology_plot
aggregate_plot


# Respiration ----

# Packages
library(dplyr)
library(lme4)
library(ggplot2)

### Data Processing ----

Respiration_Full<-read.csv("Data/Respriation_2023_F1_F2_Raw.csv")
Respiration_Full<-Respiration_Full[,-1]

Tempcatlist <- data.frame(matrix(ncol = 1, nrow = 2))
Tempcatlist[,1]<-c("25C","30C")
# For internal consistency with measurements 2022, MB measured these grasshoppers at a 2min dwell time 
# This lead to frequent saturation of the censor making the 35C data patchy. 
# Here MB included only the data for 25C and 30C

Final_output_all_temps <- data.frame(matrix(ncol = 9, nrow = 126))

for (Z in 1:2){
  Final_output_Per_Temp <- data.frame(matrix(ncol = 8, nrow = 63))
  for(R in 1:9){
    Raw<-subset(Respiration_Full, Run==R & Tempcat==Tempcatlist[Z,1])
    
    names(Raw)[1] <- "time"
    names(Raw)[65] <- "respiration"
    output<-subset(Raw,time>18& time<33&respiration>0)
    Logger_Pro_Calcd_Resp<-output[,65]
    Logger_Pro_Calcd_Resp
    
    df <- data.frame(matrix(ncol = 6, nrow = 7))
    
    colnames(df) <- c('delay_low','delay_high','end','Stop_manual','Temp_stop_14mins','Logger_Pro_Calcd')
    df$delay_low<-c(18,20,22,24,26,28,30)
    df$delay_low<-df$delay_low+0.2
    df$delay_high<-c(18,20,22,24,26,28,30)
    df$delay_high<-df$delay_high+1.5
    df$end<-c(18,20,22,24,26,28,30)
    df$end<-df$end+2
    
    names(Raw)[1] <- "time"
    names(Raw)[2] <- "CO2"
    names(Raw)[4] <- "temp"
    names(Raw)[15] <- "flow"
    
    
    for(x in 1:7) {
      
      
      Pulse_and_Background<-subset(Raw,time>df[x,1]&time<=df[x,2])
      Background<-subset(Raw,time>df[x,2]&time<=df[x,3])
      
      X<-(sum(Pulse_and_Background$CO2*(1/60))-(sum(Background$CO2*(1/60))/.5)*1.3)
      df[x,4]<-(((X*1000)/1000000)*mean(Pulse_and_Background$flow))/14
      
      Incubation<-subset(Raw,time>(df[x,3]-16)&time<=(df[x,3]))
      df[x,5]<-mean(Incubation$temp)
    }
    
    df[,6]<-(Logger_Pro_Calcd_Resp)/14
    df<-df[,4:6]
    df$chamber<-c(2,3,4,5,6,7,8)
    
    df$Run<-R
    Final_output_Per_Temp[(((R)*7)-6):((R)*7),]<-df
    
    Final_output_Per_Temp$TempCat<-Tempcatlist[Z,1]
  }
  Final_output_all_temps[((Z*63)-62):((Z*63)),]<-Final_output_Per_Temp[1:63,]
}
Final_output_all_temps<-Final_output_all_temps[,4:9]
colnames(Final_output_all_temps) <- c('Chamber','Run','Manual_Resp_Calc','Incubation_Temp','Logger_Pro_Calcd','Temp_Catagorical')


# Logger Pro automatically calculates the integral for respiration (Logger_Pro_Calcd) 
# MB also calculated the values manually (Manual_Resp_Calc)

Respiration_Mass_2023<-read.csv("Data/Respriation_Mass_2023_F1_F2.csv")
Respiration_Mass_2023<-arrange(Respiration_Mass_2023, Run,Chamber)
Respiration_Mass_2023_times2<-rbind(Respiration_Mass_2023,Respiration_Mass_2023)
respiration_data <-cbind(Final_output_all_temps,Respiration_Mass_2023_times2)

# Fix duplicate columns by selecting the first instance of each column name
respiration_data <- respiration_data[, !duplicated(names(respiration_data))]

# Append gen 1 herbivores:
# Load the new data
Gen1 <- read.csv("Data/UP_Gen_1_Herbivores_2023.csv")

# Check the structure of the new data
str(Gen1)

# Rename columns in the new data to match respiration_data
Gen1 <- Gen1 %>%
  rename(
    Entry_Order = Entry.Order,
    Saturated_at_30 = Saturated
  )

# Ensure the new data has the same columns as respiration_data
# Select only the columns that exist in respiration_data
Gen1 <- Gen1 %>%
  select(names(respiration_data))

# Append the new data to respiration_data
respiration_data <- rbind(respiration_data, Gen1)

# Check the structure of the updated respiration_data
str(respiration_data)




# Cleaning
respiration_data <- respiration_data %>%
  dplyr::mutate(Generation = case_when(
    Generation == 1 ~ "G1",
    Generation == 2 ~ "G2",
    TRUE ~ as.character(Generation)
  )) %>%
  dplyr::filter(
    Exclude != 1
  ) %>%
  dplyr::select(Logger_Pro_Calcd, 
                Temp_Catagorical, 
                Generation, 
                Population, 
                Predation, 
                Individual, 
                Mass) %>%
  dplyr::rename(Temp = Temp_Catagorical)

respiration_data <- respiration_data %>%
  mutate(Treatment = paste(Generation, Predation, sep = "_")) %>%
  mutate(Treatment = factor(Treatment, 
                           levels = c("G1_Herbivore", "G2_Herbivore", 
                                    "G1_Predator", "G2_Predator"))) %>% 
  mutate(Type = case_when(
    Population %in% c("UP", "MC", "FN", "HF") ~ "Behavior",
    Population %in% c("DC", "SP", "SC", "YF") ~ "Physiology"))

# Calculate SMR
respiration_data <- respiration_data %>%
  dplyr::mutate(SMR = Logger_Pro_Calcd / Mass)  # units are uLCO2/min/g


#' @param raw_data Raw respiration measurements
#' @param mass_data Mass measurements
#' @return Processed respiration data
process_respiration_data <- function(raw_data, mass_data) {
  # Original data processing steps
  Final_output_all_temps <- data.frame(matrix(ncol = 9, nrow = 126))
  
  for (Z in 1:2){
    Final_output_Per_Temp <- data.frame(matrix(ncol = 8, nrow = 63))
    for(R in 1:9){
      Raw <- subset(raw_data, Run==R & Tempcat==Tempcatlist[Z,1])
      
      names(Raw)[1] <- "time"
      names(Raw)[65] <- "respiration"
      output <- subset(Raw, time>18 & time<33 & respiration>0)
      Logger_Pro_Calcd_Resp <- output[,65]
      
      # ... (keep existing processing steps)
    }
    Final_output_all_temps[((Z*63)-62):((Z*63)),] <- Final_output_Per_Temp[1:63,]
  }
  
  # Join with mass data and clean
  respiration_data <- Final_output_all_temps %>%
    left_join(mass_data, by = c("Chamber", "Run")) %>%
    mutate(
      Generation = case_when(
        Generation == 1 ~ "G1",
        Generation == 2 ~ "G2",
        TRUE ~ as.character(Generation)
      ),
      Treatment = paste(Generation, Predation, sep = "_"),
      Treatment = factor(Treatment, 
                        levels = c("G1_Herbivore", "G2_Herbivore", 
                                 "G1_Predator", "G2_Predator")),
      Type = case_when(
        Population %in% c("UP", "MC", "FN", "HF") ~ "Behavior",
        Population %in% c("DC", "SP", "SC", "YF") ~ "Physiology"
      ),
      SMR = Logger_Pro_Calcd / Mass,  # Mass-specific metabolic rate
      Temp_centered = case_when(
        Temp_Catagorical == "25C" ~ -2.5,
        Temp_Catagorical == "30C" ~ 2.5
      )
    ) %>%
    filter(Exclude != 1)
  
  return(respiration_data)
}

#' Analyze respiration data
#' @param data Processed respiration data
#' @return List of model results
analyze_respiration <- function(data) {
  # Temperature-specific analyses
  temp_25_data <- data %>% filter(Temp_Catagorical == "25C")
  temp_30_data <- data %>% filter(Temp_Catagorical == "30C")
  
  # Models
  temp_25_model <- aov(SMR ~ Treatment * Type, data = temp_25_data)
  temp_30_model <- aov(SMR ~ Treatment * Type, data = temp_30_data)
  
  # Treatment-specific thermal reaction norms
  reaction_norm_models <- list(
    g1p = lmer(SMR ~ Type * Temp_centered + (1 + Temp_centered | Individual), 
               data = filter(data, Treatment == "G1_Predator")),
    g2p = lmer(SMR ~ Type * Temp_centered + (1 + Temp_centered | Individual), 
               data = filter(data, Treatment == "G2_Predator")),
    g1h = lmer(SMR ~ Type * Temp_centered + (1 + Temp_centered | Individual), 
               data = filter(data, Treatment == "G1_Herbivore"))
  )
  
  return(list(
    temp_25_model = temp_25_model,
    temp_30_model = temp_30_model,
    reaction_norm_models = reaction_norm_models
  ))
}

#' Plot respiration results
#' @param data Processed respiration data
#' @return List of plots
plot_respiration <- function(data) {
  # Temperature-specific boxplots
  temp_25_plot <- create_base_treatment_plot(
    filter(data, Temp_Catagorical == "25C"),
    "Metabolic Rate by Treatment at 25°C"
  ) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    facet_wrap(~Type) +
    ylab(expression("Mass-Specific Metabolic Rate (µL CO"[2]*"/min/g)"))
  
  temp_30_plot <- create_base_treatment_plot(
    filter(data, Temp_Catagorical == "30C"),
    "Metabolic Rate by Treatment at 30°C"
  ) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    facet_wrap(~Type) +
    ylab(expression("Mass-Specific Metabolic Rate (µL CO"[2]*"/min/g)"))
  
  # Reaction norm plot
  reaction_norm_plot <- ggplot(data, 
                              aes(x = Temp_Catagorical, y = SMR, 
                                  group = Individual, color = Treatment)) +
    geom_line(alpha = 0.2) +
    geom_point(size = 2, alpha = 0.4) +
    stat_summary(aes(group = Treatment), 
                fun.data = "mean_se",
                geom = "errorbar",
                width = 0.07,
                linewidth = 1,
                alpha = 1) +
    stat_summary(aes(group = Treatment), 
                fun = mean,
                geom = "line",
                linewidth = 1.5,
                alpha = 0.8) +
    stat_summary(aes(group = Treatment), 
                fun = mean,
                geom = "point",
                size = 3,
                alpha = 1) +
    scale_color_manual(values = c(
      custom_colors[1], custom_colors[3],
      custom_colors[5], custom_colors[7]
    )) +
    facet_wrap(~Type, scales = "free_y",
               labeller = as_labeller(c(
                 Behavior = "Behaviorally-plastic populations",
                 Physiology = "Physiologically-plastic populations"
               ))) +
    theme_light() +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          strip.text = element_text(face = "bold")) +
    labs(x = "Temperature (°C)",
         y = expression("Mass-Specific Metabolic Rate (µL CO"[2]*"/min/g)"))
  
  return(list(
    temp_25_plot = temp_25_plot,
    temp_30_plot = temp_30_plot,
    reaction_norm_plot = reaction_norm_plot
  ))
}

# Run the analysis
respiration_data <- process_respiration_data(Respiration_Full, Respiration_Mass_2023)
respiration_results <- analyze_respiration(respiration_data)
respiration_plots <- plot_respiration(respiration_data)

# Print results
summary(respiration_results$temp_25_model)
summary(respiration_results$temp_30_model)
lapply(respiration_results$reaction_norm_models, summary)

# Display plots
respiration_plots$temp_25_plot
respiration_plots$temp_30_plot
respiration_plots$reaction_norm_plot



