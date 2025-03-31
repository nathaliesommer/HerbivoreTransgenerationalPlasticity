# Temperature data
# Sommer et al. (in review)

# Load libraries
library(dplyr)
library(lubridate)
library(readr)
library(ggplot2)
library(zoo)
library(tidyr)
library(entropy)

# Import temperature data
read_temperature_files <- function(dir_path = "Data/Temperatures") {
    files <- list.files(path = dir_path, pattern = "*.csv", full.names = TRUE)
    expected_cols <- c("Event", "Date_Time", "TemperatureC", "Field", 
                      "Pop", "Cage", "Height", "Treatment")
    all_data <- do.call(rbind, lapply(files, function(file) {
        data <- read_csv(file, show_col_types = FALSE)
        
        missing_cols <- setdiff(expected_cols, names(data))
        if (length(missing_cols) > 0) {
            warning(sprintf("Missing columns in %s: %s", 
                          basename(file), 
                          paste(missing_cols, collapse = ", ")))
            for (col in missing_cols) {
                data[[col]] <- NA
            }
        }
        
        data <- data[, expected_cols]
        data$filename <- basename(file)
        return(data)
    }))
    return(all_data)
}

temp_data <- read_temperature_files()

# Date conversion - warning expected
temp_data <- temp_data %>%
    filter(!is.na(Date_Time)) %>%
    mutate(
        Date_Time = parse_date_time(Date_Time, 
                                  orders = c("mdy HM", "mdy HMS"),
                                  truncated = 3)
    )

# Add temperature diagnostics
temp_diagnostics <- temp_data %>%
    mutate(Date = as.Date(Date_Time)) %>%
    # Filter for June through October
    filter(month(Date) >= 6 & month(Date) <= 10) %>%
    # Filter for Height = 40 or 60
    filter(Height %in% c(40, 60)) %>%
    # Filter for specified field/year combinations
    filter(
        (Field == "DC" & year(Date) == 2022) |
        (Field == "FN" & year(Date) == 2022) |
        (Field == "MC" & year(Date) == 2022) |
        (Field == "HF" & year(Date) == 2022) |
        (Field == "UP" & year(Date) %in% c(2022, 2023)) |
        (Field == "SC" & year(Date) == 2022) |
        (Field == "SP" & year(Date) == 2022)
    ) %>%
    # Remove true temperature outliers
    filter(TemperatureC <= 47, TemperatureC >= -5) %>%
    group_by(Field, year(Date), Height) %>%
    summarise(
        min_temp = min(TemperatureC, na.rm = TRUE),
        max_temp = max(TemperatureC, na.rm = TRUE),
        mean_temp = mean(TemperatureC, na.rm = TRUE),
        q25 = quantile(TemperatureC, 0.25, na.rm = TRUE),
        q75 = quantile(TemperatureC, 0.75, na.rm = TRUE),
        n_outliers = sum(TemperatureC > 50 | TemperatureC < -10, na.rm = TRUE)
    ) %>%
    rename(Year = `year(Date)`)

# Check diagnostics
cat("\nTemperature Diagnostics by Field, Year, and Height:\n")
print(temp_diagnostics)

# Calculate overall means, combined height
overall_means <- temp_data %>%
    # Create date column without time
    mutate(Date = as.Date(Date_Time)) %>%
    # Filter for June through October
    filter(month(Date) >= 6 & month(Date) <= 10) %>%
    # Include all heights
    filter(Height %in% c(20, 40, 60)) %>%
    # Exclude YF - no trait data in 2023
    filter(Field != "YF") %>% 
    # Filter for specified field/year combinations
    filter(
        (Field == "DC" & year(Date) == 2022) |
        (Field == "FN" & year(Date) == 2022) |
        (Field == "MC" & year(Date) == 2022) |
        (Field == "HF" & year(Date) == 2022) |
        (Field == "UP" & year(Date) %in% c(2022, 2023)) |
        (Field == "SC" & year(Date) == 2022) |
        (Field == "SP" & year(Date) == 2022)
    ) %>%
    # Remove temperature outliers
    filter(TemperatureC <= 47, TemperatureC >= -5) %>%
    # Drop loggers with known issues
    filter(!(Field == "MC" & Height == 60) & 
           !(Field == "UP" & Height == 40 & year(Date) == 2022) &
           !(Field == "SC" & Height == 20)) %>%
    # Group by Field, Year, and Date (not Height)
    group_by(Field, year(Date), Date) %>%
    # Calculate daily maximum
    summarise(daily_max = max(TemperatureC, na.rm = TRUE)) %>%
    # Calculate mean of daily maximums and count of days
    group_by(Field, `year(Date)`) %>%
    summarise(
        mean_daily_max = mean(daily_max, na.rm = TRUE),
        n_days = n()
    ) %>%
    # Rename year column
    rename(Year = `year(Date)`)

# Print overall means
cat("\nOverall Mean Daily Maximum Temperatures (all heights):\n")
print(overall_means)

plot_data <- temp_data %>%
  mutate(Date = as.Date(Date_Time)) %>%
  filter(month(Date) >= 6 & month(Date) <= 10) %>%
  filter(Height %in% c(20, 40, 60)) %>%
  filter(TemperatureC <= 47, TemperatureC >= -5) %>%
  filter(!(Field == "MC" & Height == 60) &
           !(Field == "UP" & Height == 40 & year(Date) == 2022) &
           !(Field == "SC" & Height == 20)) %>%
  filter(
    (Field == "DC" & year(Date) == 2022) |
      (Field == "FN" & year(Date) == 2022) |
      (Field == "MC" & year(Date) == 2022) |
      (Field == "HF" & year(Date) == 2022) |
      (Field == "UP" & year(Date) %in% c(2022, 2023)) |
      (Field == "SC" & year(Date) == 2022) |
      (Field == "SP" & year(Date) == 2022)
  ) %>%
  group_by(Field, year(Date), Date) %>%
  summarise(daily_max = max(TemperatureC, na.rm = TRUE), .groups = "drop") %>%
  group_by(Field, `year(Date)`) %>%
  summarise(
    mean_daily_max = mean(daily_max, na.rm = TRUE),
    se = sd(daily_max, na.rm = TRUE) / sqrt(n()),
    n = n()
  ) %>%
  rename(Year = `year(Date)`) %>%
  mutate(
    site_group = factor(case_when(
      Field %in% c("FN", "MC", "HF") ~ "Cool-site origin",
      Field %in% c("SC", "DC", "SP") ~ "Warm-site origin",
      Field == "UP" ~ "Common garden"
    ), levels = c("Cool-site origin", "Common garden", "Warm-site origin")),
    field_group = case_when(
      Field %in% c("FN", "MC", "HF") ~ "northern",
      Field %in% c("SC", "DC", "SP") ~ "southern",
      Field == "UP" ~ "garden"
    )
  )

# Create temperature plot
ggplot(plot_data, aes(x = site_group, y = mean_daily_max)) +
  geom_point(aes(color = interaction(field_group, Year),
                 shape = field_group,
                 group = interaction(Field, Year)), 
             size = 3, 
             position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(color = interaction(field_group, Year),
                    ymin = mean_daily_max - se, 
                    ymax = mean_daily_max + se,
                    group = interaction(Field, Year)),
                width = 0.2,
                position = position_dodge(width = 0.5)) +
  geom_text(data = subset(plot_data, Field == "UP"),
            aes(label = paste("Year", ifelse(Year == 2022, "1", "2")),
                color = interaction(field_group, Year),
                group = interaction(Field, Year)),
            hjust = -0.5,
            size = 3,
            position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c(
    "northern.2022" = "#5D74A5FF",
    "northern.2023" = "#5D74A5FF",
    "southern.2022" = "#A8554EFF",
    "southern.2023" = "#A8554EFF",
    "garden.2022" = "darkgray",
    "garden.2023" = "black"
  )) +
  scale_shape_manual(values = c(
    "northern" = 16,
    "southern" = 17,
    "garden" = 15
  )) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  ) +
  labs(
    x = "",
    y = "Mean Daily Maximum Temperature (°C)",
    title = ""
  )

# Temperature Variability-Persistence Analysis ----

# Function to calculate metrics for 14-day periods
calc_14day_metrics <- function(data) {
    # Create 14-day periods
    data <- data %>%
        mutate(
            Date = as.Date(Date_Time),
            period = as.numeric(format(Date, "%j")) %/% 14  # Create 14-day periods
        )
    
    # Calculate metrics by period
    results <- data %>%
        group_by(Field, year(Date), period) %>%
        summarise(
            CV = sd(TemperatureC, na.rm = TRUE) / mean(TemperatureC, na.rm = TRUE),
            Autocorrelation = {
                if(n() >= 15) {  # Ensure enough data points
                    acf_result <- acf(TemperatureC, lag.max = 14, plot = FALSE)
                    acf_result$acf[15]
                } else {
                    NA_real_
                }
            },
            n_measurements = n(),
            mean_temp = mean(TemperatureC, na.rm = TRUE),
            start_date = min(Date),
            end_date = max(Date),
            .groups = "drop"
        ) %>%
        rename(Year = `year(Date)`)
    
    return(results)
}

# Prepare data for 14-day period analysis
filtered_metrics <- temp_data %>%
    mutate(Date = as.Date(Date_Time)) %>%
    filter(month(Date) >= 7 & month(Date) <= 10) %>%
    filter(Height %in% c(20, 40, 60)) %>%
    filter(
        (Field == "DC" & year(Date) == 2022) |
        (Field == "FN" & year(Date) == 2022) |
        (Field == "MC" & year(Date) == 2022) |
        (Field == "HF" & year(Date) == 2022) |
        (Field == "UP" & year(Date) %in% c(2022, 2023)) |
        (Field == "SC" & year(Date) == 2022) |
        (Field == "SP" & year(Date) == 2022)
    ) %>%
    filter(TemperatureC <= 47, TemperatureC >= -5) %>%
    filter(!(Field == "MC" & Height == 60) & 
           !(Field == "UP" & Height == 40 & year(Date) == 2022) &
           !(Field == "SC" & Height == 20))

# Apply to the filtered data
results_14day <- filtered_metrics %>%
    calc_14day_metrics()

# Add site type classification
results_14day_with_type <- results_14day %>%
    mutate(
        site_type = case_when(
            Field %in% c("FN", "MC", "HF") ~ "Cool-site origin",
            Field %in% c("SC", "DC", "SP") ~ "Warm-site origin",
            Field == "UP" ~ "Common garden"
        ),
        site_type = factor(site_type, 
                          levels = c("Cool-site origin", "Warm-site origin", "Common garden"))
    )

# Calculate mutual information for 14-day periods
calc_mutual_information <- function(data1, data2, n_bins = 8) {
    # Only calculate if we have enough data points
    if(length(data1) < 5 || length(data2) < 5) return(NA_real_)
    
    # Remove any NA pairs
    complete_cases <- complete.cases(data1, data2)
    data1 <- data1[complete_cases]
    data2 <- data2[complete_cases]
    
    # Standardize the data
    data1_std <- scale(data1)
    data2_std <- scale(data2)
    
    # Create bins using equal-width binning on standardized data
    breaks1 <- seq(min(data1_std), max(data1_std), length.out = n_bins + 1)
    breaks2 <- seq(min(data2_std), max(data2_std), length.out = n_bins + 1)
    
    # Ensure unique break points
    breaks1 <- unique(breaks1)
    breaks2 <- unique(breaks2)
    
    # Only proceed if we have enough unique break points
    if(length(breaks1) < 3 || length(breaks2) < 3) return(NA_real_)
    
    # Bin the data
    bins1 <- cut(data1_std, breaks = breaks1, include.lowest = TRUE)
    bins2 <- cut(data2_std, breaks = breaks2, include.lowest = TRUE)
    
    # Create joint frequency table
    joint_freq <- table(bins1, bins2)
    
    # Add small constant and normalize
    joint_freq <- joint_freq + 1e-10  # Very small constant to avoid log(0)
    joint_freq <- joint_freq / sum(joint_freq)
    
    # Calculate marginal probabilities
    p1 <- rowSums(joint_freq)
    p2 <- colSums(joint_freq)
    
    # Calculate mutual information directly
    mi <- 0
    for(i in 1:nrow(joint_freq)) {
        for(j in 1:ncol(joint_freq)) {
            if(joint_freq[i,j] > 0) {
                mi <- mi + joint_freq[i,j] * log2(joint_freq[i,j]/(p1[i]*p2[j]))
            }
        }
    }
    
    # Calculate entropy of each variable for normalization
    h1 <- -sum(p1 * log2(p1))
    h2 <- -sum(p2 * log2(p2))
    
    # Normalize by minimum entropy
    mi_normalized <- mi / min(h1, h2)
    
    return(mi_normalized)
}

# Calculate mutual information between each site and UP 2023
mi_results_14day <- results_14day_with_type %>%
    filter(Field != "UP" | (Field == "UP" & Year == 2023)) %>%
    select(Field, Year, period, mean_temp) %>%
    pivot_wider(
        names_from = c(Field, Year),
        values_from = mean_temp,
        names_sep = "_"
    ) %>%
    pivot_longer(
        cols = -c(period, UP_2023),
        names_to = c("site_2022", "year"),
        names_sep = "_",
        values_to = "site_temp"
    ) %>%
    group_by(site_2022) %>%
    summarise(
        mutual_information = calc_mutual_information(UP_2023, site_temp),
        comparison_type = "G1",
        .groups = 'drop'
    )

# Calculate environmental changes for effect size analysis
# Get UP 2023 values as reference
up_2023_values <- results_14day_with_type %>%
  filter(Field == "UP", Year == 2023) %>%
  summarise(
    up_cv = mean(CV, na.rm = TRUE),
    up_acf = mean(Autocorrelation, na.rm = TRUE)
  )

# Create environmental changes summary
env_changes <- results_14day_with_type %>%
  filter(Year == 2022) %>%  # Use 2022 values for all sites
  group_by(Field) %>%
  summarise(
    source_cv = mean(CV, na.rm = TRUE),
    source_acf = mean(Autocorrelation, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  # Calculate differences (UP 2023 - source 2022)
  mutate(
    cv_diff = up_2023_values$up_cv - source_cv,
    auto_diff = up_2023_values$up_acf - source_acf
  )

# Add mean MI by site
site_mi <- mi_results_14day %>%
  filter(comparison_type == "G1") %>%
  group_by(site_2022) %>%
  summarise(
    mean_mi = mean(mutual_information, na.rm = TRUE)
  ) %>%
  rename(Field = site_2022)

# Combine all environmental metrics
env_changes <- env_changes %>%
  left_join(site_mi, by = "Field")




