# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)

# --- 1. Load and Prepare the Data ---

# Create a directory to store the new plots
if (!dir.exists("plots/regional_plots")) {
  dir.create("plots/regional_plots", recursive = TRUE)
}

# Load the wide-format data
disturbance_data_wide <- read_csv("exports/n2k_disturbance_r_results_wide_format.csv")

# Reshape the data from wide to long format to make it plottable
# We are only interested in the columns containing absolute hectares ('haYYYY')
disturbance_data_long <- disturbance_data_wide %>%
  # Select only the necessary columns
  select(country = `country.x`, region_name, starts_with("ha2")) %>%
  # Pivot the 'haYYYY' columns into two new columns: 'year' and 'hectares'
  pivot_longer(
    cols = starts_with("ha2"),
    names_to = "year",
    values_to = "hectares_disturbed"
  ) %>%
  # The 'year' column currently looks like 'ha2001', 'ha2002', etc.
  # We need to extract just the numbers.
  mutate(
    year = as.numeric(str_extract(year, "\\d{4}"))
  ) %>%
  # For this plot, we are interested in the total disturbance per region per year.
  # We group by country, region, and year, and sum the hectares.
  group_by(country, region_name, year) %>%
  summarise(total_hectares_disturbed = sum(hectares_disturbed, na.rm = TRUE), .groups = 'drop') %>%
  # Remove any rows where region name is missing
  filter(!is.na(region_name))


# --- 2. Loop Through Each Country and Generate Plots ---

# Get a unique list of countries from the dataset
countries <- unique(disturbance_data_long$country)

# Loop through each country in the list
for (current_country in countries) {

  # Filter the data for only the current country
  country_data <- disturbance_data_long %>%
    filter(country == current_country)

  # Create the plot for the current country
  regional_plot <- ggplot(country_data, aes(x = year, y = total_hectares_disturbed, color = region_name)) +
    geom_line(linewidth = 1) +
    labs(
      title = paste("Annual Forest Disturbance in", current_country),
      subtitle = "Trends by region (2001-2023)",
      x = "Year",
      y = "Area of Forest Disturbance (hectares)",
      color = "Region" # Legend title
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

  # Define the filename for the plot
  file_name <- paste0("plots/regional_plots/disturbance_trend_", current_country, ".jpg")

  # Save the plot
  ggsave(file_name, plot = regional_plot, width = 11, height = 7, dpi = 300)

  # Print a message to track progress
  message(paste("Saved plot for", current_country, "to", file_name))
}

message("\nAll regional plots have been generated and saved.")