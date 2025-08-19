# new_results_transform.R
#
# Purpose: Transform the R-based analysis results (long format) into the
#          wide format of the original ArcGIS Pro export for consistency.

# --- Load libraries --- #
library(dplyr)
library(readr)
library(tidyr)
library(sf)

# --- Step 1: Load required datasets --- #

cat("--- Loading required data ---\n")

# Load the new R-based results (long format)
new_results_path <- "exports/n2k_disturbance_by_year_no_fire.csv"
new_results <- read_csv(new_results_path, show_col_types = FALSE)

# Load the old ArcGIS results to get metadata and target column structure
old_results_path <- "imports/export_n2kloggedareas_20241030.csv"
old_results <- read_csv(old_results_path, show_col_types = FALSE)

# Load the authoritative site areas from the Natura 2000 geopackage
n2k_gpkg_path <- "D:/GIS/Data/Raw Files/2024-05 Natura 2000 Europe/eea_v_3035_100_k_natura2000_p_2022_v01_r00/Natura2000_end2022.gpkg"
n2k_layername <- "NATURA2000SITES"

site_areas <- st_read(n2k_gpkg_path, layer = n2k_layername, quiet = TRUE) %>%
  st_drop_geometry() %>%
  select(SITECODE, totalha = AREAHA) %>%
  # Ensure there is only one area per site
  distinct(SITECODE, .keep_all = TRUE)

# --- Step 2: Prepare and Reshape R-based results --- #

cat("--- Reshaping new R results to wide format ---\n")

# Clean the new results: filter for relevant years and select columns
new_results_clean <- new_results %>%
  filter(year > 2000) %>%
  select(SITECODE = SITECODE, year, disturbed_ha_r = disturbed_ha)

# Ensure every site has a row for every year (2001-2023) for a complete pivot
new_results_complete <- new_results_clean %>%
  tidyr::complete(SITECODE, year = 2001:2023, fill = list(disturbed_ha_r = 0))

# Pivot the disturbance data from long to wide format (ha2001, ha2002, etc.)
new_results_wide_ha <- new_results_complete %>%
  mutate(year_col = paste0("ha", year)) %>%
  pivot_wider(
    id_cols = SITECODE,
    names_from = year_col,
    values_from = disturbed_ha_r
  )

# --- Step 3: Join data and calculate derived columns --- #

cat("--- Joining data and calculating percentage columns ---\n")

# Join the disturbance data with the site area data
# Use a right_join to ensure we only keep sites present in the area file
recreated_data <- right_join(new_results_wide_ha, site_areas, by = "SITECODE") %>%
  # Replace NAs in hectare columns with 0 (for sites with no disturbance)
  mutate(across(starts_with("ha20"), ~replace_na(., 0)))

# To calculate percentages efficiently, pivot longer, calculate, then pivot wider
recreated_long_perc <- recreated_data %>%
  select(SITECODE, totalha, starts_with("ha20")) %>%
  pivot_longer(
    cols = starts_with("ha20"),
    names_to = "year_col",
    values_to = "disturbed_ha"
  ) %>%
  # Calculate percentage, handling division by zero
  mutate(
    perc = if_else(totalha > 0, (disturbed_ha / totalha), 0),
    year = readr::parse_number(year_col)
  )

# Pivot the new percentage columns back to a wide format
recreated_wide_perc <- recreated_long_perc %>%
  select(SITECODE, year, perc) %>%
  mutate(year_col = paste0("perc", year)) %>%
  pivot_wider(
    id_cols = SITECODE,
    names_from = year_col,
    values_from = perc
  )

# Calculate the final 'perclogged' summary column
perclogged_summary <- recreated_long_perc %>%
  group_by(SITECODE) %>%
  summarise(perclogged = sum(perc, na.rm = TRUE), .groups = "drop")

# --- Step 4: Assemble and export the final table --- #

cat("--- Assembling, ordering, and exporting final table ---\n")

# Extract metadata from the old results file
old_results_meta <- old_results %>%
  select(sitecode, country.x, country.y, region_gid1, region_name) %>%
  distinct(sitecode, .keep_all = TRUE)

# Join all the pieces together and select columns in the target order
final_output <- recreated_data %>%
  left_join(recreated_wide_perc, by = "SITECODE") %>%
  left_join(perclogged_summary, by = "SITECODE") %>%
  rename(sitecode = SITECODE) %>%
  left_join(old_results_meta, by = "sitecode") %>%
  mutate(`...1` = row_number()) %>%
  # rename_with(~ if_else(. == "...1", "", .)) %>%
  select(any_of(names(old_results)))

# Export to a new CSV file
output_path <- "exports/n2k_disturbance_r_results_wide_format.csv"
write_csv(final_output, output_path, na = "")
cat(paste("\n--- Transformation complete. File saved to:", output_path, "---\n"))
