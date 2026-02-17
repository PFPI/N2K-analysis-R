# area_ratio_analysis.R
# 
# Purpose: Calculate the average ratio of N2K site forest area to buffer forest area,
#          and generate a distribution graph of this ratio against N2K forest area.

library(dplyr)
library(readr)
library(ggplot2)
library(scales)

# --- 1. Load Data ---

data_path <- "imports/full_disturbance_and_forest_area_by_site.csv"
if (!file.exists(data_path)) {
  stop("Data file not found. Please check 'imports/full_disturbance_and_forest_area_by_site.csv'")
}

site_data <- read_csv(data_path, show_col_types = FALSE)

# --- 2. Calculate Ratios ---

# Calculate ratio of Inside Area / Buffer Area
# Note: Using 'forest_ha' columns as these are the available area metrics.
ratio_data <- site_data %>%
  mutate(
    ratio_inside_to_buffer = forest_ha_inside / forest_ha_buffer
  ) %>%
  # Filter out infinite or NA ratios (e.g., if buffer has 0 forest)
  filter(is.finite(ratio_inside_to_buffer), !is.na(ratio_inside_to_buffer))

# --- 3. Statistical Summary ---

average_ratio <- mean(ratio_data$ratio_inside_to_buffer, na.rm = TRUE)
median_ratio <- median(ratio_data$ratio_inside_to_buffer, na.rm = TRUE)

message(paste0("Analysis complete for ", nrow(ratio_data), " sites."))
message(paste0("Average Ratio (Inside/Buffer): ", round(average_ratio, 4)))
message(paste0("Median Ratio (Inside/Buffer):  ", round(median_ratio, 4)))

# Save summary to text file
writeLines(
  c(
    paste0("Analysis Date: ", Sys.Date()),
    paste0("Total Sites: ", nrow(ratio_data)),
    paste0("Average Ratio (Inside/Buffer): ", round(average_ratio, 4)),
    paste0("Median Ratio (Inside/Buffer): ", round(median_ratio, 4))
  ),
  "summaries/ratio_summary.txt"
)

# --- 4. Generate Distribution Graph ---

# Create directory for plot if it doesn't exist
dir.create("plots", recursive = TRUE, showWarnings = FALSE)

p <- ggplot(ratio_data, aes(x = forest_ha_inside, y = ratio_inside_to_buffer)) +
  geom_point(alpha = 0.5, color = "#2c7fb8") +
  geom_hline(yintercept = average_ratio, color = "#d95f02", linetype = "dashed", size = 1) +
  scale_x_log10(labels = comma) + # Log scale for area as it varies significantly
  scale_y_continuous(labels = comma) +
  labs(
    title = "Distribution of N2K/Buffer Area Ratio by Site Size",
    subtitle = paste0("Mean Ratio: ", round(average_ratio, 2), " (Orange Dashed Line)"),
    x = "N2K Forest Area (ha) [Log Scale]",
    y = "Ratio (Inside Area / Buffer Area)",
    caption = "PFPI Analysis | Source: Global Forest Change & EEA"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text = element_text(color = "gray20")
  )

# Save the plot
output_plot_path <- "plots/n2k_buffer_ratio_distribution.png"
ggsave(output_plot_path, plot = p, width = 8, height = 6, bg = "white")

message(paste0("Graph saved to: ", output_plot_path))