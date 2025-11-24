# member_state_stats.R
# 
# Purpose: Calculate detailed disturbance statistics by Member State and Year,
#          and generate country-level trend graphs with detailed summaries.

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(scales) # For nice number formatting

# --- 1. Load Data ---

# Load the site-level summary (contains Country/Member State info and total forest area)
site_summary_path <- "analysis_outputs/full_disturbance_and_forest_area_by_site.csv"
if (!file.exists(site_summary_path)) {
  stop("Site summary file not found. Please check the path.")
}
site_summary <- read_csv(site_summary_path, show_col_types = FALSE)

# Load the yearly disturbance data (contains Site, Year, and Disturbed Hectares)
yearly_data_path <- "exports/n2k_disturbance_by_year_no_fire.csv"
if (!file.exists(yearly_data_path)) {
  stop("Yearly disturbance file not found. Please check the path.")
}
yearly_data <- read_csv(yearly_data_path, show_col_types = FALSE)

# Create output directories
dir.create("analysis_outputs/stats", recursive = TRUE, showWarnings = FALSE)
dir.create("plots/country_trends", recursive = TRUE, showWarnings = FALSE)

# --- 2. Data Preparation ---

# Ensure we have the member_state for every row in the yearly data
# We join the yearly data with the site_summary to get the 'member_state' column
yearly_data_joined <- yearly_data %>%
  left_join(site_summary %>% select(SITECODE, member_state), by = "SITECODE") %>%
  filter(!is.na(member_state)) # Remove rows if sitecode didn't match

# --- 3. Calculate General Member State Statistics (Totals) ---

# Calculate totals for the entire period/country
country_totals <- site_summary %>%
  group_by(member_state) %>%
  summarise(
    total_n2k_sites = n(),
    sites_with_disturbance = sum(disturbed_ha_inside > 0, na.rm = TRUE),
    total_forest_ha = sum(forest_ha_inside, na.rm = TRUE),
    total_disturbed_ha = sum(disturbed_ha_inside, na.rm = TRUE)
  ) %>%
  mutate(
    pct_sites_disturbed = (sites_with_disturbance / total_n2k_sites) * 100,
    pct_forest_disturbed = (total_disturbed_ha / total_forest_ha) * 100
  ) %>%
  arrange(member_state)

# Save the total stats
write_csv(country_totals, "analysis_outputs/stats/member_state_totals.csv")
message("Saved 'member_state_totals.csv'")

# --- 4. Calculate Yearly Statistics by Member State ---

country_yearly_stats <- yearly_data_joined %>%
  group_by(member_state, year) %>%
  summarise(
    # Count distinct sites that had disturbance recorded in this specific year
    sites_experiencing_clearcut = n_distinct(SITECODE),
    # Sum the disturbed hectares for this year
    disturbed_ha_in_year = sum(disturbed_ha, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(member_state, year)

# Save the yearly stats
write_csv(country_yearly_stats, "analysis_outputs/stats/member_state_yearly_stats.csv")
message("Saved 'member_state_yearly_stats.csv'")

# --- 5. Generate Plots per Country ---

unique_countries <- unique(country_yearly_stats$member_state)

for (country in unique_countries) {
  
  # Filter data for the specific country
  plot_data <- country_yearly_stats %>% 
    filter(member_state == country)
  
  # Retrieve summary stats for this country
  stats <- country_totals %>% filter(member_state == country)
  
  # Format numbers for the subtitle (adding commas, rounding)
  site_text <- paste0(label_comma(stats$sites_with_disturbance), " sites disturbed(", 
                      label_comma(stats$total_n2k_sites), " total sites)")
  
  area_text <- paste0(comma(round(stats$total_disturbed_ha)), " ha forest disturbed (", 
                      comma(round(stats$total_forest_ha)), " ha total)")
  
  subtitle_text <- paste0("N2K Stats: ", site_text, "\nArea Stats: ", area_text)
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = year, y = disturbed_ha_in_year)) +
    geom_col(fill = "#E69F00", alpha = 0.8) + # Bar chart for annual values
    geom_line(color = "#D55E00", linewidth = 1) + # Line to show trend
    #geom_smooth(method = "loess", se = FALSE, color = "#56B4E9", linetype = "dashed", linewidth = 0.8) + # Smoothed trend
    labs(
      title = paste("Annual Forest Disturbance:", country),
      subtitle = subtitle_text,
      x = "Year",
      y = "Disturbed Area (ha)",
      caption = "PFPI Analysis based on GFC Data"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray30"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_y_continuous(labels = comma) + # Add commas to y-axis numbers
    scale_x_continuous(breaks = seq(min(plot_data$year), max(plot_data$year), by = 2))
  
  # Save the plot
  ggsave(
    filename = paste0("plots/country_trends/disturbance_trend_", country, ".png"),
    plot = p,
    width = 9, # Slightly wider to accommodate text
    height = 6,
    bg = "white"
  )
}

message("Plots generated in 'plots/country_trends/'")