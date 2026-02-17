# tables.R (Corrected Again)
# Load necessary libraries
library(dplyr)
library(gt)
library(readr)

# --- Table 1: Summary of Forest Disturbance by Member State ---

# Load the necessary files with correct column names
# This file contains the detailed site-by-site data.
disturbance_by_site <- read_csv("imports/full_disturbance_and_forest_area_by_site.csv")
# This file contains the pre-calculated statistical trend data (Sen's Slope).
member_state_stats <- read_csv("imports/member_state_statistical_summary.csv")

# Summarize the site data by member state
country_summary <- disturbance_by_site %>%
  group_by(member_state) %>%
  summarise(
    # Sum up the forest and disturbance areas for each state
    total_n2k_forest_area = sum(forest_ha_inside, na.rm = TRUE),
    total_n2k_disturbance_area = sum(disturbed_ha_inside, na.rm = TRUE)
  ) %>%
  # Calculate disturbance as a percentage of the total forest area
  mutate(
    disturbance_as_percent_of_forest = (total_n2k_disturbance_area / total_n2k_forest_area) * 100
  )

# Join the summary with the statistical trend data
# Note: The column in the stats file is also 'member_state', so the join is direct.
final_summary <- country_summary %>%
  left_join(
    member_state_stats %>% select(member_state, sens_slope_inside),
    by = "member_state"
  )

table1 <- final_summary %>%
  # Select and rename columns for the final presentation table
  select(
    `Member State` = member_state,
    `Total N2K Forest Area (ha)` = total_n2k_forest_area,
    `Gross Disturbance in N2K Sites (ha, 2001-2023)` = total_n2k_disturbance_area,
    `Annual Rate of Change (ha/yr)` = sens_slope_inside,
    `Disturbance as % of Total N2K Forest` = disturbance_as_percent_of_forest
  ) %>%
  # Arrange by the highest percentage of disturbance for impact
  arrange(desc(`Disturbance as % of Total N2K Forest`)) %>%
  gt() %>%
  tab_header(
    title = "Table 1: Forest Disturbance in Natura 2000 Sites by Member State (2001-2023)",
    subtitle = "Disturbance totals and rates of change, ranked by proportional impact."
  ) %>%
  fmt_number(columns = c(2:4), decimals = 0) %>%
  fmt_percent(columns = 5, decimals = 2, scale_values = FALSE) %>%
  tab_source_note(source_note = "Data derived from Hansen et al. Global Forest Change and EEA Natura 2000 boundaries.")

print(table1)
gtsave(table1, "tables/t1-forest_disturbance_by_member_state.docx")

# --- Table 2: Visual Verification Results ---
# This table remains unchanged as it is derived from the manuscript text.

verification_data <- tibble(
  `Classification Category` = c("Probable Clearcut", "Probable Thinning", "False Positive", "Undeterminable"),
  `Percentage of Sample` = c(NA, NA, 5, 25)
)

table2 <- verification_data %>%
  gt() %>%
  tab_header(
    title = "Table 2: Summary of Visual Verification of Potential Disturbance Events",
    subtitle = "A total of 976 points were examined using Google Earth Pro."
  ) %>%
  fmt_missing(columns = everything(), missing_text = "Not specified") %>%
  tab_source_note(source_note = "Based on manual classification of high-resolution satellite imagery.")

print(table2)
gtsave(table2, "tables/t2-visual_verification_disturbance_events.docx")

# --- Table 3: Top 10 Most Disturbed Natura 2000 Sites ---
# This table is corrected to use the available SITECODE instead of SITENAME.

table3 <- disturbance_by_site %>%
  # Arrange by the total disturbance area in descending order
  arrange(desc(disturbed_ha_inside)) %>%
  head(10) %>%
  # CORRECTED: Select SITECODE and member_state
  select(
    `Site ID` = SITECODE,
    `Member State` = member_state,
    `Total Disturbance (ha, 2001-2023)` = disturbed_ha_inside
  ) %>%
  gt() %>%
  tab_header(
    title = "Table 3: Top 10 Natura 2000 Sites by Gross Forest Disturbance Area",
    subtitle = "These sites represent the largest cumulative disturbance from 2001-2023."
  ) %>%
  fmt_number(columns = 3, decimals = 0)

print(table3)
gtsave(table3, "tables/t3-top_10_n2k_by_gross_disturbed_area.docx")