# Temperatures

# Load libraries
library(dplyr)
library(lubridate)
library(readr)
library(ggplot2)

# Read all CSV files in the temperatures folder
read_temperature_files <- function(dir_path = "Data/Temperatures") {
    files <- list.files(path = dir_path, pattern = "*.csv", full.names = TRUE)
    expected_cols <- c("Event", "Date_Time", "TemperatureC", "Field", 
                      "Pop", "Cage", "Height", "Treatment")
    all_data <- do.call(rbind, lapply(files, function(file) {

        data <- read_csv(file, show_col_types = FALSE)
        
        # Check columns
        missing_cols <- setdiff(expected_cols, names(data))
        if (length(missing_cols) > 0) {
            warning(sprintf("Missing columns in %s: %s", 
                          basename(file), 
                          paste(missing_cols, collapse = ", ")))
            # Add missing columns with NA values
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

# Read all temperature data
temp_data <- read_temperature_files()

# Date conversion - warning expected
temp_data <- temp_data %>%
    # Remove rows with NA dates
    filter(!is.na(Date_Time)) %>%
    mutate(
        Date_Time = parse_date_time(Date_Time, 
                                  orders = c("mdy HM", "mdy HMS"),
                                  truncated = 3)
    )

# Add temperature diagnostics
temp_diagnostics <- temp_data %>%
    # Create date column without time
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
    # Remove temperature outliers
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


# Calculate overall means (combining heights)
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
cat("\nOverall Mean Daily Maximum Temperatures (combining heights):\n")
print(overall_means)






# Create plot data
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

# Create the plot
ggplot(plot_data, aes(x = site_group, y = mean_daily_max)) +
  geom_point(aes(color = interaction(field_group, Year),
                 group = interaction(Field, Year)),  # This ensures jittering for all sites
             size = 3, 
             position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(color = interaction(field_group, Year),
                    ymin = mean_daily_max - se, 
                    ymax = mean_daily_max + se),
                width = 0.1,
                position = position_dodge(width = 0.3)) +
  geom_text(data = subset(plot_data, Field == "UP"),
            aes(label = paste("Year", ifelse(Year == 2022, "1", "2")),
                color = interaction(field_group, Year),
                group = interaction(Field, Year)),
            hjust = -0.5,
            size = 3,
            position = position_dodge(width = 0.3)) +
  scale_color_manual(values = c(
    "northern.2022" = "#5D74A5FF",
    "northern.2023" = "#5D74A5FF",
    "southern.2022" = "#A8554EFF",
    "southern.2023" = "#A8554EFF",
    "garden.2022" = "darkgray",
    "garden.2023" = "black"
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

