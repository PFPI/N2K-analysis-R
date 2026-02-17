# ratio_regression.R
# 
# Purpose: Fit trend lines to the N2K/Buffer ratio distribution and calculate
#          regression statistics to quantify the "hockey stick" effect.

library(dplyr)
library(readr)
library(ggplot2)
library(scales)
library(broom) # For tidy model outputs

# --- 1. Load and Prepare Data ---

data_path <- "analysis_outputs/full_disturbance_and_forest_area_by_site.csv"
if (!file.exists(data_path)) stop("Data file not found.")

site_data <- read_csv(data_path, show_col_types = FALSE)

# Calculate ratio and filter valid data
regression_data <- site_data %>%
  mutate(
    ratio_inside_to_buffer = forest_ha_inside / forest_ha_buffer,
    log_forest_area = log10(forest_ha_inside)
  ) %>%
  filter(
    is.finite(ratio_inside_to_buffer),
    !is.na(ratio_inside_to_buffer),
    forest_ha_inside > 0 # Ensure log calculation is valid
  )

# --- 2. Statistical Modeling ---

# Model 1: Linear Regression on Log(Area)
# Tests if the ratio increases linearly with orders of magnitude of area
linear_model <- lm(ratio_inside_to_buffer ~ log_forest_area, data = regression_data)

# Extract stats
model_summary <- tidy(linear_model)
r_squared <- summary(linear_model)$r.squared

message("--- Linear Model Summary (Ratio ~ Log10(Area)) ---")
print(model_summary)
message(paste0("R-squared: ", round(r_squared, 4)))

# Save stats to file
sink("analysis_outputs/stats/ratio_regression_stats.txt")
cat("Regression Analysis: N2K/Buffer Ratio vs Site Area\n")
cat("==================================================\n\n")
cat("Model: Ratio ~ Log10(Forest Area)\n")
cat(paste0("R-squared: ", round(r_squared, 4), "\n"))
cat(paste0("P-value: ", format.pval(summary(linear_model)$coefficients[2,4]), "\n\n"))
print(summary(linear_model))
sink()

# --- 3. Generate Plot with Fits ---

p <- ggplot(regression_data, aes(x = forest_ha_inside, y = ratio_inside_to_buffer)) +
  geom_point(alpha = 0.3, color = "#2c7fb8", size = 1) +
  
  # Add the 'Hockey Stick' smooth fit (GAM)
  # This allows the line to curve naturally with the data
  geom_smooth(
    method = "gam", 
    formula = y ~ s(x, bs = "cs"), # Cubic regression spline
    color = "#d95f02", 
    linewidth = 1,
    fill = "#d95f02",
    alpha = 0.2
  ) +
  
  scale_x_log10(labels = comma) + 
  scale_y_continuous(labels = comma) +
  labs(
    title = "Trend Analysis: N2K/Buffer Ratio vs Site Size",
    subtitle = "Orange line: GAM Smooth Fit (Visualizing the 'Hockey Stick')",
    x = "N2K Forest Area (ha) [Log Scale]",
    y = "Ratio (Inside Area / Buffer Area)",
    caption = paste0("Linear Fit R2: ", round(r_squared, 3), " | Source: GFC & EEA")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text = element_text(color = "gray20")
  )

# Save the plot
ggsave("plots/supporting/ratio_trend_fit.png", plot = p, width = 8, height = 6, bg = "white")
message("Plot saved to 'plots/supporting/ratio_trend_fit.png'")