# member_state_stats.R
# Purpose: Generate country-level trend graphs matching the style of Figure 2
# (Comparing Inside N2K vs Buffer with lines, points, and trends)

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(scales)

# --- 1. Load Data ---

# Load Site Info (to map sites to Member States)
site_summary_path <- "imports/full_disturbance_and_forest_area_by_site.csv"
if (!file.exists(site_summary_path)) stop("Site summary file not found.")
site_summary <- read_csv(site_summary_path, show_col_types = FALSE)

# Load N2K Disturbance
n2k_path <- "imports/n2k_disturbance_by_year_no_fire.csv"
if (!file.exists(n2k_path)) stop("N2K disturbance file not found.")
n2k_data <- read_csv(n2k_path, show_col_types = FALSE)

# Load Buffer Disturbance (Required for the Orange/Blue comparison)
buffer_path <- "imports/n2k_buffer_disturbance_by_year_no_fire.csv"
if (!file.exists(buffer_path)) stop("Buffer disturbance file not found.")
buffer_data <- read_csv(buffer_path, show_col_types = FALSE)

# Create output directories
dir.create("plots/country_trends", recursive = TRUE, showWarnings = FALSE)

# --- 2. Data Preparation & Aggregation ---

# Helper function to attach Member State and sum by Year
process_disturbance <- function(dist_data, sites, col_name) {
  dist_data %>%
    left_join(sites %>% select(SITECODE, member_state), by = "SITECODE") %>%
    filter(!is.na(member_state)) %>%
    group_by(member_state, year) %>%
    summarise(!!col_name := sum(disturbed_ha, na.rm = TRUE), .groups = "drop")
}

# Process both datasets
n2k_aggregated <- process_disturbance(n2k_data, site_summary, "Inside Natura 2000")
buffer_aggregated <- process_disturbance(buffer_data, site_summary, "1-km Buffer Zone")

# Join them into one master dataset
full_country_stats <- full_join(n2k_aggregated, buffer_aggregated, by = c("member_state", "year")) %>%
  arrange(member_state, year)

# --- 3. Generate Plots per Country ---

unique_countries <- unique(full_country_stats$member_state)

for (country in unique_countries) {
  
  # Filter for the specific country
  country_data_wide <- full_country_stats %>% 
    filter(member_state == country)
  
  # Pivot longer for ggplot (essential for the multi-line style)
  country_data_long <- country_data_wide %>%
    pivot_longer(
      cols = c("Inside Natura 2000", "1-km Buffer Zone"),
      names_to = "Location",
      values_to = "Area_ha"
    )
  
  # Create the plot matching Figure 2 styling
  p <- ggplot(country_data_long, aes(x = year, y = Area_ha, color = Location)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) + 
    # Dashed trend line (Linear Model)
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
    
    # Custom Colors: Blue for Inside, Orange for Buffer
    scale_color_manual(values = c("Inside Natura 2000" = "#0072B2", "1-km Buffer Zone" = "#D55E00")) +
    
    # Format Y axis with "k" for thousands
    scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "k")) +
    
    labs(
      # title = paste(country, "- Annual Forest Disturbance"),
      # subtitle = "Comparison of Disturbance Inside N2K vs 1-km Buffer",
      x = "Year",
      y = "Area of Forest Disturbance (ha)",
      color = "Location"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # Save the plot
  ggsave(
    filename = paste0("plots/country_trends/trend_", country, ".jpg"),
    plot = p,
    width = 10,
    height = 6,
    dpi = 300
  )
}

message("Country trend plots (styled like Figure 2) generated in 'plots/country_trends/'")