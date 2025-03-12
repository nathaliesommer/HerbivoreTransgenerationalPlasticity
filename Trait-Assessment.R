# Sommer et al 
# Update title


# Behavior ----

# Packages 
library(dplyr)
library(car)
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

behavior_data <- behavior_data %>%
  filter(Population != "UP")

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
    Population %in% c("MC", "FN", "HF") ~ "Cool-site origin",
    Population %in% c("DC", "SP", "SC") ~ "Warm-site origin"))


# Center population heights
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
    filter(n() >= min_obs) %>%
    summarize(
      # Calculate quantiles directly
      avg_core_low = quantile(Y_centered, probs = 0.25),
      avg_core_high = quantile(Y_centered, probs = 0.75),
      avg_broad_low = quantile(Y_centered, probs = 0.025),
      avg_broad_high = quantile(Y_centered, probs = 0.975),
      
      # Mean for positioning
      Y_centered = mean(Y_centered),
      .groups = 'drop'
    )
}

#' Print diagnostic information for density ranges
#' @param data The density range data
print_diagnostic <- function(data) {
  data %>%
    group_by(Treatment, Type) %>%
    summarize(
      core_range = paste(round(avg_core_low, 1), "to", round(avg_core_high, 1)),
      core_width = round(avg_core_high - avg_core_low, 1),
      broad_range = paste(round(avg_broad_low, 1), "to", round(avg_broad_high, 1)),
      broad_width = round(avg_broad_high - avg_broad_low, 1),
      mean_height = round(Y_centered, 1)
    ) %>%
    print(n = Inf)
}

# Calculate ranges and print diagnostics
density_estimate_ranges <- calculate_density_ranges(behavior_data)
print("\nDensity Range Diagnostics:")
print_diagnostic(density_estimate_ranges)

# Calculate separate ranges for cool-site and warm-site origins
density_estimate_ranges_behavior <- density_estimate_ranges %>%
  filter(Type == "Cool-site origin")

density_estimate_ranges_physiology <- density_estimate_ranges %>%
  filter(Type == "Warm-site origin")

### Analysis ----


#### ANOVA on cool-site populations data ----
behavior_aov <- aov(Y_centered ~ Treatment, 
                    data = subset(behavior_data, Type == "Cool-site origin"))

# Tukey's HSD
behavior_tukey <- TukeyHSD(behavior_aov)

# Extract comparisons of interest
selected_comparisons <- behavior_tukey$Treatment[c(
  "G2_Herbivore-G1_Herbivore",  # H1 vs H2
  "G1_Predator-G1_Herbivore",   # H1 vs P1
  "G2_Predator-G1_Predator",     # P1 vs P2
  "G2_Predator-G2_Herbivore"  #H2 vs P2
), ]

# Print results
print("Selected Treatment Comparisons for Cool-site origin Type:")
print(selected_comparisons)


#### ANOVA on warm-site populations data ----
physiology_aov <- aov(Y_centered ~ Treatment, 
                    data = subset(behavior_data, Type == "Warm-site origin"))

# Run Tukey's HSD
physiology_tukey <- TukeyHSD(physiology_aov)

# Extract comparisons of interest
selected_comparisons <- physiology_tukey$Treatment[c(
  "G2_Herbivore-G1_Herbivore",  # H1 vs H2
  "G1_Predator-G1_Herbivore",   # H1 vs P1
  "G2_Predator-G1_Predator", # P1 vs P2
  "G2_Predator-G2_Herbivore"  #H2 vs P2
), ]

# Print results
print("Selected Treatment Comparisons for Warm-site origin Type:")
print(selected_comparisons)

#### ANOVA for all data ----
all_data_aov <- aov(Y_centered ~ Treatment, data = behavior_data)
summary(all_data_aov)

# Tukey's HSD for all data
all_data_tukey <- TukeyHSD(all_data_aov)

# Extract comparisons of interest
selected_comparisons_all <- all_data_tukey$Treatment[c(
  "G2_Herbivore-G1_Herbivore",  # H1 vs H2
  "G1_Predator-G1_Herbivore",   # H1 vs P1
  "G2_Predator-G1_Predator",
  "G2_Predator-G2_Herbivore" 
), ]

# Print results
print("Selected Treatment Comparisons for All Data:")
print(selected_comparisons_all)

## Distribution test: Wilcoxon rank-sum test ----

# Update Type labels in the data
behavior_data <- behavior_data %>%
  mutate(Type = case_when(
    Type == "Behavior" ~ "Cool-site origin",
    Type == "Physiology" ~ "Warm-site origin",
    TRUE ~ Type
  ))

#' Run distribution comparisons using Wilcoxon rank-sum test
#' @param data The behavior dataset
#' @return Data frame of comparison results
compare_distributions <- function(data) {
  results <- data.frame(
    Comparison_Group = character(),
    Comparison = character(),
    W_statistic = numeric(),
    P_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Function to run Wilcoxon test and add results
  add_wilcox_result <- function(data1, data2, group_name, comparison_name) {
    w_result <- wilcox.test(data1, data2)
    new_row <- data.frame(
      Comparison_Group = group_name,
      Comparison = comparison_name,
      W_statistic = w_result$statistic,
      P_value = w_result$p.value
    )
    return(new_row)
  }
  
  # Cool-site origin comparisons
  cool_site_data <- subset(data, Type == "Cool-site origin")
  results <- rbind(results,
    add_wilcox_result(
      cool_site_data$Y_centered[cool_site_data$Treatment == "G1_Herbivore"],
      cool_site_data$Y_centered[cool_site_data$Treatment == "G2_Herbivore"],
      "Cool-site origin", "G1 Herb vs G2 Herb"
    ),
    add_wilcox_result(
      cool_site_data$Y_centered[cool_site_data$Treatment == "G1_Predator"],
      cool_site_data$Y_centered[cool_site_data$Treatment == "G1_Herbivore"],
      "Cool-site origin", "G1 Pred vs G1 Herb"
    ),
    add_wilcox_result(
      cool_site_data$Y_centered[cool_site_data$Treatment == "G1_Predator"],
      cool_site_data$Y_centered[cool_site_data$Treatment == "G2_Predator"],
      "Cool-site origin", "G1 Pred vs G2 Pred"
    ),
    add_wilcox_result(
      cool_site_data$Y_centered[cool_site_data$Treatment == "G2_Herbivore"],
      cool_site_data$Y_centered[cool_site_data$Treatment == "G2_Predator"],
      "Cool-site origin", "G2 Herb vs G2 Pred"
    )
  )
  
  # Warm-site origin comparisons
  warm_site_data <- subset(data, Type == "Warm-site origin")
  results <- rbind(results,
    add_wilcox_result(
      warm_site_data$Y_centered[warm_site_data$Treatment == "G1_Herbivore"],
      warm_site_data$Y_centered[warm_site_data$Treatment == "G2_Herbivore"],
      "Warm-site origin", "G1 Herb vs G2 Herb"
    ),
    add_wilcox_result(
      warm_site_data$Y_centered[warm_site_data$Treatment == "G1_Predator"],
      warm_site_data$Y_centered[warm_site_data$Treatment == "G1_Herbivore"],
      "Warm-site origin", "G1 Pred vs G1 Herb"
    ),
    add_wilcox_result(
      warm_site_data$Y_centered[warm_site_data$Treatment == "G1_Predator"],
      warm_site_data$Y_centered[warm_site_data$Treatment == "G2_Predator"],
      "Warm-site origin", "G1 Pred vs G2 Pred"
    ),
    add_wilcox_result(
      warm_site_data$Y_centered[warm_site_data$Treatment == "G2_Herbivore"],
      warm_site_data$Y_centered[warm_site_data$Treatment == "G2_Predator"],
      "Warm-site origin", "G2 Herb vs G2 Pred"
    )
  )
  
  # Type comparisons for each treatment
  results <- rbind(results,
    add_wilcox_result(
      data$Y_centered[data$Type == "Warm-site origin" & data$Treatment == "G1_Herbivore"],
      data$Y_centered[data$Type == "Cool-site origin" & data$Treatment == "G1_Herbivore"],
      "G1 Herbivores", "Warm-site vs Cool-site origin"
    ),
    add_wilcox_result(
      data$Y_centered[data$Type == "Warm-site origin" & data$Treatment == "G2_Herbivore"],
      data$Y_centered[data$Type == "Cool-site origin" & data$Treatment == "G2_Herbivore"],
      "G2 Herbivores", "Warm-site vs Cool-site origin"
    ),
    add_wilcox_result(
      data$Y_centered[data$Type == "Warm-site origin" & data$Treatment == "G1_Predator"],
      data$Y_centered[data$Type == "Cool-site origin" & data$Treatment == "G1_Predator"],
      "G1 Predators", "Warm-site vs Cool-site origin"
    ),
    add_wilcox_result(
      data$Y_centered[data$Type == "Warm-site origin" & data$Treatment == "G2_Predator"],
      data$Y_centered[data$Type == "Cool-site origin" & data$Treatment == "G2_Predator"],
      "G2 Predators", "Warm-site vs Cool-site origin"
    )
  )
  
  # Add significance indicators
  results$Significance <- ifelse(results$P_value < 0.001, "***",
                               ifelse(results$P_value < 0.01, "**",
                                     ifelse(results$P_value < 0.05, "*", "ns")))
  
  return(results)
}

# Run the comparisons
wilcox_results <- compare_distributions(behavior_data)

# Display results in a nice format
print("Distribution Comparison Results:")
print(knitr::kable(wilcox_results, digits = 3))



## Figures ----

# Define custom colors
custom_colors <- c(
  "#484A5A", "#9fa1ac",  # G1 Herbivore (dark, light) - slate blue
  "#C69648", "#ddc095",  # G1 Predator (dark, light) - golden
  "#9BA48C", "#c5cabb",  # G2 Herbivore (dark, light) - sage
  "#6B4E71", "#a895ad"   # G2 Predator (dark, light) - purple
)

# Define treatment labels for plots
treatment_labels <- c(
  "G1_Herbivore" = "G[1] Herbivore",
  "G1_Predator" = "G[1] Predator",
  "G2_Herbivore" = "G[2] Herbivore",
  "G2_Predator" = "G[2] Predator"
)

#' Create single panel plot for cool-site/warm-site specific plots
create_single_panel_plot <- function(data, type) {
  # Calculate means for the plot
  means <- data %>%
    filter(Type == type) %>%
    group_by(Treatment) %>%
    summarise(
      Y_centered = mean(Y_centered),
      .groups = 'drop'
    )
  
  ggplot(data %>% filter(Type == type), aes(x = Treatment, y = Y_centered)) +
    # Violin plots
    geom_violin(
      aes(fill = Treatment),
      alpha = 0.5,
      trim = TRUE
    ) +
    # Add jittered points
    geom_jitter(
      aes(color = Treatment),
      alpha = 0.2,
      size = 0.8,
      width = 0.15
    ) +
    # Mean points
    geom_point(
      data = means,
      size = 3,
      color = "black"
    ) +
    # Mean labels with background
    geom_label(
      data = means,
      aes(label = round(Y_centered, 1)),
      vjust = -0.5,
      size = 3.5,
      label.padding = unit(0.2, "lines"),
      label.size = 0,
      fill = "white",
      alpha = 0.8
    ) +
    # Styling
    scale_y_continuous(limits = c(0, 80)) +
    scale_x_discrete(
      labels = parse(text = c(
        "G[1]~Herbivore",
        "G[2]~Herbivore",
        "G[1]~Predator",
        "G[2]~Predator"
      ))
    ) +
    scale_fill_manual(values = c(
      "G1_Herbivore" = custom_colors[1],
      "G2_Herbivore" = custom_colors[3],
      "G1_Predator" = custom_colors[5],
      "G2_Predator" = custom_colors[7]
    )) +
    scale_color_manual(values = c(
      "G1_Herbivore" = custom_colors[1],
      "G2_Herbivore" = custom_colors[3],
      "G1_Predator" = custom_colors[5],
      "G2_Predator" = custom_colors[7]
    )) +
    theme_light() +
    theme(
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_line(color = "grey95"),
      legend.position = "none",
      axis.text.x = element_text(size = 10)
    ) +
    labs(
      title = paste0(type, " populations"),
      y = "Population-centered canopy height",
      x = NULL
    )
}


# Create the aggregated plot
aggregated_plot <- ggplot(behavior_data, aes(x = Treatment, y = Y_centered)) +
    # Violin plots
    geom_violin(
        aes(fill = Treatment),
        alpha = 0.5,
        trim = TRUE
    ) +
    # Add jittered points
    geom_jitter(
        aes(color = Treatment),
        alpha = 0.2,
        size = 0.8,
        width = 0.15
    ) +
    # Mean points and labels
    stat_summary(
        fun = mean,
        geom = "point",
        size = 3,
        color = "black"
    ) +
    # Add single mean value label per group with background
    stat_summary(
        fun = mean,
        geom = "label",
        aes(label = round(..y.., 1)),
        vjust = -0.5,
        size = 3.5,
        label.padding = unit(0.2, "lines"),
        fill = "white",
        alpha = 0.8
    ) +
    scale_x_discrete(
        labels = parse(text = c(
            "G[1]~Herbivore",
            "G[2]~Herbivore",
            "G[1]~Predator",
            "G[2]~Predator"
        ))
    ) +
    scale_fill_manual(values = c(
        "G1_Herbivore" = custom_colors[1],
        "G2_Herbivore" = custom_colors[3],
        "G1_Predator" = custom_colors[5],
        "G2_Predator" = custom_colors[7]
    )) +
    scale_color_manual(values = c(
        "G1_Herbivore" = custom_colors[1],
        "G2_Herbivore" = custom_colors[3],
        "G1_Predator" = custom_colors[5],
        "G2_Predator" = custom_colors[7]
    )) +
    theme_light() +
    theme(panel.grid.major.x = element_blank(),
        legend.position = "none"
    ) +
    labs(
        x = "Treatment",
        y = "Population-centered canopy height",
        title = "Aggregated Population Heights"
    )

# Temperature plot
temperature_plot <- ggplot(behavior_data, aes(x = Nearest_Temperature, y = Y, color = Treatment)) +
  geom_jitter(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c(
    custom_colors[1],  # G1 Herbivore
    custom_colors[3],  # G2 Herbivore
    custom_colors[5],  # G1 Predator
    custom_colors[7]   # G2 Predator
  ),
  labels = expression(
    G[1]~Herbivore,
    G[2]~Herbivore, 
    G[1]~Predator,
    G[2]~Predator
  )) +
  labs(title = "No relationship between ambient temperature and canopy height",
       x = "Temperature (°C)", 
       y = "Height (cm)") +
  theme_light() +
  theme(plot.title = element_text(face = "bold"),
        text = element_text(size = 10),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom") +
  facet_wrap(~Treatment)

# Display all plots
aggregated_plot

# Create the individual plots for cool-site and warm-site origins
behavior_plot <- create_single_panel_plot(behavior_data, "Cool-site origin")
physiology_plot <- create_single_panel_plot(behavior_data, "Warm-site origin")

# Display the plots
behavior_plot
physiology_plot
temperature_plot


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

# Append gen 1 herbivores data set:
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

# Drop UP
respiration_data <- respiration_data %>%
  filter(Population != "UP")


# Cleaning
respiration_data <- respiration_data %>%
  dplyr::mutate(
    # Convert Generation to G1/G2 format
    Generation = case_when(
      Generation == "1" ~ "G1",
      Generation == "2" ~ "G2",
      Generation == 1 ~ "G1",  # Handle numeric values
      Generation == 2 ~ "G2",  # Handle numeric values
      TRUE ~ as.character(Generation)
    )
  ) %>%
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

# Ensure proper factor levels and data types
respiration_data <- respiration_data %>%
  mutate(
    Generation = factor(Generation, levels = c("G1", "G2")),
    Treatment = paste(Generation, Predation, sep = "_"),
    Treatment = factor(Treatment, 
                      levels = c("G1_Herbivore", "G2_Herbivore", 
                                "G1_Predator", "G2_Predator")),
    Type = case_when(
      Population %in% c("MC", "FN", "HF") ~ "Cool-site origin",
      Population %in% c("DC", "SP", "SC") ~ "Warm-site origin"
    ),
    Type = factor(Type),
    Predation = factor(Predation)
  )

# Calculate SMR
respiration_data <- respiration_data %>%
  dplyr::mutate(SMR = Logger_Pro_Calcd / Mass)  # units are uLCO2/min/g

# Print data structure for verification
print("Respiration data structure after preparation:")
print(str(respiration_data))
print("\nUnique values in Generation column:")
print(table(respiration_data$Generation))
print("\nUnique values in Treatment column:")
print(table(respiration_data$Treatment))
print("\nUnique values in Type column:")
print(table(respiration_data$Type))

### Analysis ----



# Center temperature to make intercepts more interpretable
respiration_data <- respiration_data %>%
  mutate(Temp_centered = case_when(
    Temp == "25C" ~ -2.5,
    Temp == "30C" ~ 2.5
  ))


# Test for differences between treatments at 25C
respiration_25C <- respiration_data %>%
  filter(Temp == "25C")

aov_25C <- aov(SMR ~ Treatment*Type, data = respiration_25C)
summary(aov_25C)


# Plot for 25C
ggplot(respiration_25C, aes(x = Treatment, y = SMR, fill = Treatment)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_fill_manual(values = c(custom_colors[1], custom_colors[3], custom_colors[5], custom_colors[7])) +
  labs(title = "Metabolic Rate by Treatment at 25°C", 
       x = "Treatment", 
       y = expression("Mass-Specific Metabolic Rate (µL CO"[2]*"/min/g)")) +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        text = element_text(size = 12),
        legend.position = "none") +
  facet_wrap(~Type)

# Test for differences between treatments at 30C
respiration_30C <- respiration_data %>%
  filter(Temp == "30C")

# ANOVA instead of t-test for three groups
aov_30C <- aov(SMR ~ Treatment*Type, data = respiration_30C)
summary(aov_30C)

# Plot for 30C
ggplot(respiration_30C, aes(x = Treatment, y = SMR, fill = Treatment)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_fill_manual(values = c(custom_colors[1], custom_colors[3], custom_colors[5], custom_colors[7])) +
  labs(title = "Metabolic Rate by Treatment at 30°C", 
       x = "Treatment", 
       y = expression("Mass-Specific Metabolic Rate (µL CO"[2]*"/min/g)")) +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        text = element_text(size = 12),
        legend.position = "none") +
  facet_wrap(~Type)



# Test for differences in thermal reaction norms between Types within each Treatment

# Separate models for each Treatment
# G1 Predator
g1p_data <- respiration_data %>% 
  filter(Treatment == "G1_Predator")
g1p_model <- lmer(SMR ~ Type * Temp_centered + (1 | Individual), data = g1p_data)

# G2 Predator
g2p_data <- respiration_data %>% 
  filter(Treatment == "G2_Predator")
g2p_model <- lmer(SMR ~ Type * Temp_centered + (1 | Individual), data = g2p_data)

# G1 Herbivore
g1h_data <- respiration_data %>% 
  filter(Treatment == "G1_Herbivore")
g1h_model <- lmer(SMR ~ Type * Temp_centered + (1 | Individual), data = g1h_data)

# G2 Herbivore: Not applicable; no G2 herbivores of behavior

# Print interpretable results for each comparison
cat("\nG1 Predator - Behavior vs Physiology:\n")
cat("Main effect of Type:", fixef(g1p_model)["TypePhysiology"], "\n")
cat("Interaction (difference in slopes):", fixef(g1p_model)["TypePhysiology:Temp_centered"], "\n")
cat("P-values from ANOVA:\n")
print(anova(g1p_model))

cat("\nG2 Predator - Behavior vs Physiology:\n")
cat("Main effect of Type:", fixef(g2p_model)["TypePhysiology"], "\n")
cat("Interaction (difference in slopes):", fixef(g2p_model)["TypePhysiology:Temp_centered"], "\n")
cat("P-values from ANOVA:\n")
print(anova(g2p_model))

cat("\nG1 Herbivore - Behavior vs Physiology:\n")
cat("Main effect of Type:", fixef(g1h_model)["TypePhysiology"], "\n")
cat("Interaction (difference in slopes):", fixef(g1h_model)["TypePhysiology:Temp_centered"], "\n")
cat("P-values from ANOVA:\n")
print(anova(g1h_model))




## Reaction Norm Plot by Type ----
# Reaction Norm Plot
ggplot(respiration_data, aes(x = Temp, y = SMR, group = Treatment, color = Treatment)) +
  # Add individual reaction norms
  geom_line(aes(group = Individual), alpha = 0.2) +
  # Add mean reaction norms
  stat_summary(fun = mean, geom = "line", linewidth = 1.5) +
  # Add mean points
  stat_summary(fun = mean, geom = "point", size = 3) +
  # Styling
  scale_color_manual(values = c(
    "G1_Herbivore" = custom_colors[1],
    "G2_Herbivore" = custom_colors[3],
    "G1_Predator" = custom_colors[5],
    "G2_Predator" = custom_colors[7]
  ),
  labels = expression(
    G[1]~Herbivore,
    G[2]~Herbivore, 
    G[1]~Predator,
    G[2]~Predator
  )) +
  scale_x_discrete(labels = c("25°C", "30°C")) +
  facet_wrap(~factor(Type, levels = c("Cool-site origin", "Warm-site origin")), 
             ncol = 1) +
  theme_light() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 12)
  ) +
  labs(
    x = "Temperature (°C)",
    y = expression(paste("Mass-specific metabolic rate (", mu, "L ", CO[2], "/min/g)")))




### Aggregated reaction norm plot ----
# Reaction Norm Plot
ggplot(respiration_data, aes(x = Temp, y = SMR, group = Treatment, color = Treatment)) +
  # Add individual reaction norms
  geom_line(aes(group = Individual), alpha = 0.2) +
  # Add mean reaction norms
  stat_summary(fun = mean, geom = "line", linewidth = 1.5) +
  # Add mean points
  stat_summary(fun = mean, geom = "point", size = 3) +
  # Styling
  scale_color_manual(values = c(
    "G1_Herbivore" = custom_colors[1],
    "G2_Herbivore" = custom_colors[3],
    "G1_Predator" = custom_colors[5],
    "G2_Predator" = custom_colors[7]
  ),
  labels = expression(
    G[1]~Herbivore,
    G[2]~Herbivore, 
    G[1]~Predator,
    G[2]~Predator
  )) +
  scale_x_discrete(labels = c("25°C", "30°C")) +
  theme_light() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 12)
  ) +
  labs(
    x = "Temperature (°C)",
    y = expression(paste("Mass-specific metabolic rate (", mu, "L ", CO[2], "/min/g)")))



# Effect Size Comparisons ----

# Run Temperatures.R first to generate env_changes and mi_results_14day objects

#' Calculate generational trait changes
#' @param data Dataset containing trait measurements
#' @param trait_col Name of the trait column
#' @return Dataframe of trait changes between G1 and G2
calc_trait_changes <- function(data, trait_col) {
  result <- data %>%
    mutate(
      Generation = factor(Generation, levels = c("G1", "G2"))
    ) %>%
    group_by(Population, Type, Predation) %>%
    summarise(
      g1_trait = mean(get(trait_col)[Generation == "G1"], na.rm = TRUE),
      g2_trait = mean(get(trait_col)[Generation == "G2"], na.rm = TRUE),
      n_g1 = sum(Generation == "G1"),
      n_g2 = sum(Generation == "G2"),
      var_g1 = if(sum(Generation == "G1") > 1) 
                 var(get(trait_col)[Generation == "G1"], na.rm = TRUE) 
               else NA_real_,
      var_g2 = if(sum(Generation == "G2") > 1) 
                 var(get(trait_col)[Generation == "G2"], na.rm = TRUE) 
               else NA_real_,
      trait_change = g2_trait - g1_trait,
      .groups = 'drop'
    ) %>%
    filter(!is.na(trait_change))
  
  return(result)
}

# Prepare behavior data
behavior_changes <- behavior_data %>%
  rename(Predation = Predation.Treatment) %>%
  mutate(
    Generation = factor(Generation, levels = c("G1", "G2")),
    Type = factor(Type),
    Predation = factor(Predation)
  ) %>%
  calc_trait_changes("Y_centered")

#' Calculate predation effect sizes
#' @param changes Dataframe of trait changes
#' @param effect_name Name for the effect type
#' @return Dataframe of standardized effect sizes
calc_predation_effects <- function(changes, effect_name) {
  changes %>%
    group_by(Type) %>%
    summarise(
      pred_effect = case_when(
        Type == "Cool-site origin" ~ mean(trait_change[Predation == "Predator"], na.rm = TRUE),
        Type == "Warm-site origin" ~ mean(trait_change[Predation == "Predator"], na.rm = TRUE) - 
                                    mean(trait_change[Predation == "Herbivore"], na.rm = TRUE)
      ),
      n_pred = sum(Predation == "Predator"),
      n_herb = sum(Predation == "Herbivore"),
      var_pred = if(sum(Predation == "Predator") > 1) 
                   var(trait_change[Predation == "Predator"], na.rm = TRUE) 
                 else NA_real_,
      var_herb = case_when(
        Type == "Cool-site origin" ~ 0,
        Type == "Warm-site origin" & sum(Predation == "Herbivore") > 1 ~ 
          var(trait_change[Predation == "Herbivore"], na.rm = TRUE),
        TRUE ~ NA_real_
      ),
      pooled_var = case_when(
        Type == "Cool-site origin" & n_pred >= 2 ~ var_pred/n_pred,
        Type == "Warm-site origin" & n_pred >= 2 & n_herb >= 1 ~ 
          (var_pred/n_pred + var_herb/n_herb),
        TRUE ~ NA_real_
      ),
      SE = sqrt(pmax(pooled_var, 0)),
      std_effect = if_else(!is.na(pooled_var) & pooled_var > 0,
                          pred_effect/sqrt(pooled_var),
                          NA_real_),
      effect_type = effect_name,
      .groups = 'drop'
    )
}

#' Calculate environmental effects (CV, autocorrelation, and MI)
#' @param changes Dataframe of trait changes
#' @param env_var Environmental variable to correlate with
#' @param effect_name Name for the effect type
#' @return Dataframe of standardized effect sizes
calc_env_effects <- function(changes, env_var, effect_name) {
  changes %>%
    left_join(env_changes, by = c("Population" = "Field")) %>%
    group_by(Type) %>%
    summarise(
      pred_effect = {
        x <- trait_change
        y <- get(env_var)
        complete_cases <- !is.na(x) & !is.na(y)
        x <- x[complete_cases]
        y <- y[complete_cases]
        
        if(length(x) >= 3 && sd(x) > 0 && sd(y) > 0) {
          if(env_var == "mean_mi") {
            mean(y, na.rm = TRUE) * 2 - 1
          } else {
            cor(x, y)
          }
        } else {
          NA_real_
        }
      },
      n_pred = sum(!is.na(get(env_var)) & !is.na(trait_change)),
      n_herb = n_pred,
      pooled_var = case_when(
        env_var == "mean_mi" ~ if(n_pred >= 3) 1/n_pred else NA_real_,
        TRUE ~ if(n_pred >= 4) 1/(n_pred - 3) else NA_real_
      ),
      SE = sqrt(pmax(pooled_var, 0)),
      std_effect = case_when(
        env_var == "mean_mi" ~ pred_effect,
        !is.na(pred_effect) && !is.na(pooled_var) ~ atanh(pred_effect),
        TRUE ~ NA_real_
      ),
      effect_type = effect_name,
      .groups = 'drop'
    )
}

# Calculate effects
effects_list <- list(
  # Predation effects
  height_pred = calc_predation_effects(behavior_changes, "Predation"),
  
  # Temperature CV effects
  height_cv = calc_env_effects(behavior_changes, "cv_diff", "CV (°C)"),
  
  # Temperature autocorrelation effects
  height_acf = calc_env_effects(behavior_changes, "auto_diff", "Autocorrelation (°C)"),
  
  # Mutual Information effects
  height_mi = calc_env_effects(behavior_changes, "mean_mi", "Mutual Information (°C)")
)

# Combine effects and set factor levels for proper ordering
combined_effects <- bind_rows(effects_list) %>%
  filter(!is.na(std_effect), !is.na(SE)) %>%
  mutate(effect_type = factor(effect_type, 
                             levels = c("Predation",
                                      "CV (°C)",
                                      "Autocorrelation (°C)",
                                      "Mutual Information (°C)")))

# Create forest plot
plot_effects <- function(data) {
  ggplot(data, aes(y = effect_type, x = std_effect, color = Type)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_errorbarh(aes(xmin = std_effect - 1.96*SE, 
                      xmax = std_effect + 1.96*SE),
                  position = position_dodge(width = 0.5),
                  height = 0.2) +
    scale_color_manual(values = c("Cool-site origin" = "#5D74A5FF", 
                                "Warm-site origin" = "#A8554EFF")) +
    theme_light() +
    labs(x = "Standardized Effect Size",
         y = NULL) +
    theme(legend.position = "right")
}

# Create and display plot
height_plot <- plot_effects(combined_effects)
print(height_plot)

# Correlation Assessment ----

# Calculate correlations between trait changes and environmental variables
trait_env_correlations <- behavior_changes %>%
  left_join(env_changes, by = c("Population" = "Field")) %>%
  group_by(Type) %>%
  summarise(
    # CV correlations
    cv_cor = cor(trait_change, cv_diff, use = "complete.obs"),
    cv_n = sum(complete.cases(trait_change, cv_diff)),
    cv_p = cor.test(trait_change, cv_diff, use = "complete.obs")$p.value,
    
    # Autocorrelation correlations
    acf_cor = cor(trait_change, auto_diff, use = "complete.obs"),
    acf_n = sum(complete.cases(trait_change, auto_diff)),
    acf_p = cor.test(trait_change, auto_diff, use = "complete.obs")$p.value,
    
    # Mutual Information correlations
    mi_cor = cor(trait_change, mean_mi, use = "complete.obs"),
    mi_n = sum(complete.cases(trait_change, mean_mi)),
    mi_p = cor.test(trait_change, mean_mi, use = "complete.obs")$p.value,
    .groups = 'drop'
  )

# Print correlation results
print("\nCorrelation Results:")
print(trait_env_correlations)

# Create correlation plots
cv_plot <- ggplot(behavior_changes %>% 
                    left_join(env_changes, by = c("Population" = "Field")),
                  aes(x = cv_diff, y = trait_change, color = Type)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Change in CV", y = "Trait Change") +
  theme_light()

acf_plot <- ggplot(behavior_changes %>% 
                     left_join(env_changes, by = c("Population" = "Field")),
                   aes(x = auto_diff, y = trait_change, color = Type)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Change in Autocorrelation", y = "Trait Change") +
  theme_light()

mi_plot <- ggplot(behavior_changes %>% 
                    left_join(env_changes, by = c("Population" = "Field")),
                  aes(x = mean_mi, y = trait_change, color = Type)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Mutual Information", y = "Trait Change") +
  theme_light()

# Display correlation plots
print(cv_plot)
print(acf_plot)
print(mi_plot)