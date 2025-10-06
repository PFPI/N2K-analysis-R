# figures.R (Corrected Again)
# Load necessary libraries
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(scales)

# --- Figure 2: Annual Disturbance Trend (Time Series Plot) ---
# This plot's code was correct as it used different source files. No changes needed here.

# --- checking figure 2 --- #
# Load the two datasets used for the plot
n2k_disturbance <- read_csv("exports/n2k_disturbance_by_year_no_fire.csv")
buffer_disturbance <- read_csv("exports/n2k_buffer_disturbance_by_year_no_fire.csv")

n2k_disturbance_summary <- n2k_disturbance %>%
  group_by(year) %>%
  summarise(Disturbance_Inside_N2K_ha = sum(disturbed_ha, na.rm = TRUE))

buffer_disturbance_summary <- buffer_disturbance %>%
  group_by(year) %>%
  summarise(Disturbance_in_Buffer_ha = sum(disturbed_ha, na.rm = TRUE))

# Now, the join will work correctly because 'year' is unique in both tables
disturbance_summary_table <- full_join(
    n2k_disturbance_summary,
    buffer_disturbance_summary,
    by = "year"
  ) %>%
  # Arrange by year to ensure correct calculations
  arrange(year) %>%
  # Calculate year-over-year change
  mutate(
    Pct_Change_Inside = ((Disturbance_Inside_N2K_ha - lag(Disturbance_Inside_N2K_ha)) / lag(Disturbance_Inside_N2K_ha)) * 100,
    Pct_Change_in_Buffer = ((Disturbance_in_Buffer_ha - lag(Disturbance_in_Buffer_ha)) / lag(Disturbance_in_Buffer_ha)) * 100
  ) %>%
  # Select and reorder columns for the final table
  select(
    Year = year,
    `Disturbance Inside N2K (ha)` = Disturbance_Inside_N2K_ha,
    `YoY Change Inside (%)` = Pct_Change_Inside,
    `Disturbance in Buffer (ha)` = Disturbance_in_Buffer_ha,
    `YoY Change in Buffer (%)` = Pct_Change_in_Buffer
  )
# Use the 'gt' package to create a publication-quality table
summary_gt <- disturbance_summary_table %>%
  gt() %>%
  tab_header(
    title = "Annual Forest Disturbance Summary",
    subtitle = "Year-over-year analysis for Natura 2000 sites and 1-km buffer zones."
  ) %>%
  # Format numbers for readability
  fmt_number(
    columns = c(2, 4),
    decimals = 0
  ) %>%
  fmt_number(
    columns = c(3, 5),
    decimals = 1
  ) %>%
  # Highlight the year 2018 for easy reference
  tab_style(
    style = list(cell_fill(color = "yellow")),
    locations = cells_body(rows = Year == 2018)
  ) %>%
  # Replace NA values from the first year of change calculation
  fmt_missing(
    columns = everything(),
    missing_text = "---"
  )

# Print the final table
print(summary_gt)

disturbance_summary_wide  <- full_join(
    n2k_disturbance_summary,
    buffer_disturbance_summary,
    by = "year"
  )

# To plot both lines with ggplot, we need to pivot the data to a 'long' format
disturbance_summary_long <- disturbance_summary_wide %>%
  pivot_longer(
    cols = c("Disturbance_Inside_N2K_ha", "Disturbance_in_Buffer_ha"),
    names_to = "Location",
    values_to = "Area_ha"
  ) %>%
  # Make the location names more readable for the plot legend
  mutate(Location = case_when(
    Location == "Disturbance_Inside_N2K_ha" ~ "Inside Natura 2000",
    Location == "Disturbance_in_Buffer_ha"  ~ "1-km Buffer Zone"
  ))

# Create the plot using our new long-format summary table
figure2_plot <- ggplot(disturbance_summary_long, aes(x = year, y = Area_ha, color = Location)) +
  geom_line(linewidth = 1.2) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
  scale_color_manual(values = c("Inside Natura 2000" = "#0072B2", "1-km Buffer Zone" = "#D55E00")) +
  
  scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "k")) +
  labs(
    title = "Figure 2: Annual Forest Disturbance Inside and Outside Natura 2000 Sites",
    subtitle = "Disturbance has shown a significant positive trend in both protected and unprotected areas.",
    x = "Year",
    y = "Area of Forest Disturbance (hectares)",
    color = "Location"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top", 
        axis.line.x = element_line(color = "black", linewidth = 0.5))

# Save the plot to a file
ggsave("plots/Figure_2_Annual_Disturbance_Trend.jpg", plot = figure2_plot, width = 10, height = 6, dpi = 300)

# Add a confirmation message
message("Figure 2 has been successfully generated from the summary table and saved to the 'figures' directory.")

# --- Figure 3: Country Comparison Bar Chart ---
# This chart's code is updated to reflect the correct column names.

# Load the data again for this figure
disturbance_by_site <- read_csv("analysis_outputs/full_disturbance_and_forest_area_by_site.csv")

# Perform the same summary as in the tables script
country_summary_for_plot <- disturbance_by_site %>%
  group_by(member_state) %>%
  summarise(
    total_n2k_forest_area = sum(forest_ha_inside, na.rm = TRUE),
    total_n2k_disturbance_area = sum(disturbed_ha_inside, na.rm = TRUE)
  ) %>%
  mutate(
    disturbance_as_percent_of_forest = (total_n2k_disturbance_area / total_n2k_forest_area) * 100
  )

figure3_plot <- country_summary_for_plot %>%
  # Use reorder() to sort the countries by the disturbance percentage for a clean look
  mutate(member_state = reorder(member_state, disturbance_as_percent_of_forest)) %>%
  ggplot(aes(x = member_state, y = disturbance_as_percent_of_forest)) +
  geom_bar(stat = "identity", fill = "#009E73") +
  coord_flip() + # Flip coordinates to make country names readable
  labs(
    title = "Figure 3: Forest Disturbance in N2K Sites as a Percentage of Total N2K Forest Area",
    subtitle = "Ranking by proportional impact highlights the most-affected member states.",
    x = "Member State",
    y = "Disturbance as % of Natura 2000 Forest Area"
  ) +
  theme_minimal(base_size = 12)

print(figure3_plot)
ggsave("plots/figure3.jpg", figure3_plot)


