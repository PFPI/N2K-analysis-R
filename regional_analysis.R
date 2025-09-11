# regional_analysis.R
#
# Purpose: This script breaks down forest disturbance by GADM Level 1 
#          (regional) boundaries within each EU Member State. It identifies
#          the most disturbed regions and the most disturbed N2K sites
#          within those regions.
#
# The process is as follows:
# 1. Load necessary data:
#    - Pre-calculated disturbance stats per N2K site.
#    - N2K site geometries.
#    - GADM Level 1 administrative boundaries for all EU countries.
# 2. Spatially join N2K sites to their corresponding administrative region.
# 3. Aggregate disturbance data by region to find regional totals.
# 4. Identify and rank the top 3 most disturbed regions per country.
# 5. Identify and rank the top 3 most disturbed sites within each region.
# 6. Save the results to both CSV and a formatted text file.

# --- Load Libraries --- #
library(sf)
library(dplyr)
library(readr)
library(purrr) # For map() function to iterate over files

# --- Configuration --- #

# Input data paths
site_disturbance_data_path <- "analysis_outputs/full_disturbance_and_forest_area_by_site.csv"
n2k_gpkg_path <- "D:/GIS/Data/Raw Files/2024-05 Natura 2000 Europe/eea_v_3035_100_k_natura2000_p_2022_v01_r00/Natura2000_end2022.gpkg"
n2k_layername <- "NaturaSite_polygon"
gadm_base_path <- "D:/GIS/Data/Political Boundaries/EU Subdivisions" # The main folder

# Target CRS
target_crs <- "EPSG:3035"

# Output paths
output_dir <- "analysis_outputs"
region_text_output_path <- "exports/region_output.txt"
regional_summary_csv_path <- file.path(output_dir, "regional_disturbance_summary.csv")
top_sites_csv_path <- file.path(output_dir, "top_disturbed_sites_by_region.csv")

# --- Part 1: Load All Necessary Data --- #
message("--- Part 1: Loading Data ---")

# 1a: Load the pre-calculated site-level disturbance data
message("Loading site-level disturbance data...")
site_data <- read_csv(site_disturbance_data_path, show_col_types = FALSE)

# 1b: Load N2K site geometries. We'll use their centroids for the join.
message("Loading Natura 2000 site geometries...")
n2k_sites_sf <- st_read(n2k_gpkg_path, layer = n2k_layername) %>%
  select(SITECODE) %>%
  st_transform(crs = target_crs) %>%
  st_make_valid()

# Create centroids for a clean spatial join
n2k_centroids <- st_centroid(n2k_sites_sf)

# 1c: Find all GADM Level 1 shapefiles
message("Finding all GADM Level 1 regional boundaries...")
country_folders <- list.dirs(gadm_base_path, recursive = FALSE, full.names = TRUE)

shapefile_paths <- map(country_folders, ~ list.files(
    path = .x,
    pattern = "gadm41_.*_1\\.shp$",
    full.names = TRUE,
    recursive = TRUE
  )) %>%
  list_flatten() %>%
  unlist()

if (length(shapefile_paths) == 0) {
  stop("No GADM Level 1 shapefiles found! Check 'gadm_base_path' and file pattern.")
}

message(paste("Found", length(shapefile_paths), "regional shapefiles. Processing one country at a time..."))

# --- Part 2: Spatially Join N2K Sites to Regions (Iteratively) --- #
message("\n--- Part 2: Joining N2K sites to regions ---")

# Create a function to process a single country's shapefile
process_country <- function(shp_path, all_sites_centroids) {
  # Load the regions for just one country
  country_regions <- st_read(shp_path, quiet = TRUE) %>%
    st_transform(crs = target_crs) %>%
    st_make_valid() %>%
    select(COUNTRY, NAME_1) # Keep only necessary columns

  # Get the country code from the data to filter N2K sites
  country_name_from_gadm <- country_regions$COUNTRY[1]
  
  # A small lookup list to match GADM country names to N2K country codes (SITECODE prefix)
  # This may need to be expanded!
  country_lookup <- c("France" = "FR", "Spain" = "ES", "Germany" = "DE", 
                      "Italy" = "IT", "Poland" = "PL", "Romania" = "RO",
                      "Sweden" = "SE", "Greece" = "GR", "Portugal" = "PT",
                      "Austria" = "AT", "Belgium" = "BE", "Bulgaria" = "BG",
                      "Croatia" = "HR", "Cyprus" = "CY", "Czechia" = "CZ",
                      "Denmark" = "DK", "Estonia" = "EE", "Finland" = "FI",
                      "Hungary" = "HU", "Ireland" = "IE", "Latvia" = "LV",
                      "Lithuania" = "LT", "Luxembourg" = "LU", "Malta" = "MT",
                      "Netherlands" = "NL", "Slovakia" = "SK", "Slovenia" = "SI")

  country_code <- names(country_lookup[which(country_lookup == country_name_from_gadm)])[1]
  if (is.na(country_code)) {
       country_code <- country_lookup[country_name_from_gadm]
  }


  if (is.na(country_code)) {
    message(paste("Warning: No country code found for", country_name_from_gadm, "- skipping."))
    return(NULL)
  }
  
  message(paste("Processing:", country_name_from_gadm))

  # Filter N2K sites for the current country
  sites_in_country <- all_sites_centroids %>%
    filter(substr(SITECODE, 1, 2) == country_code)
    
  if (nrow(sites_in_country) == 0) {
    message("No N2K sites found for this country. Skipping.")
    return(NULL)
  }

  # Perform the spatial join
  sites_with_region <- st_join(sites_in_country, country_regions, join = st_intersects) %>%
    st_drop_geometry() # We only need the table of SITECODE and NAME_1

  return(sites_with_region)
}

# Apply this function to every shapefile path
# map_df will automatically combine the results into a single data frame
all_sites_with_regions <- map_df(shapefile_paths, ~ process_country(.x, n2k_centroids))

message("Spatial join complete for all countries.")
glimpse(all_sites_with_regions)

# --- Part 3: Regional Analysis and Ranking --- #
message("\n--- Part 3: Aggregating disturbance data by region ---")

# First, ensure there are no duplicate SITECODEs in our region lookup table
# This can happen if a site centroid falls exactly on a boundary
sites_with_region_unique <- all_sites_with_regions %>%
  distinct(SITECODE, .keep_all = TRUE)

# Join the regional information back to the main site disturbance data
full_regional_data <- site_data %>%
  left_join(sites_with_region_unique, by = "SITECODE") %>%
  # Filter out any sites we couldn't match to a region
  filter(!is.na(NAME_1))

# (a) Which regions are the most disturbed?
# Aggregate the data by country and region to get total disturbance
regional_summary <- full_regional_data %>%
  group_by(COUNTRY, NAME_1) %>%
  summarise(
    total_forest_ha = sum(forest_ha_inside, na.rm = TRUE),
    total_disturbed_ha = sum(disturbed_ha_inside, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Calculate the percentage of forest that was disturbed
  mutate(
    percent_disturbed = (total_disturbed_ha / total_forest_ha) * 100
  ) %>%
  # Handle cases with zero forest area to avoid NaN results
  mutate(percent_disturbed = ifelse(total_forest_ha == 0, 0, percent_disturbed))

message("Regional aggregation complete.")
glimpse(regional_summary)

# Now, let's get the rankings you asked for.
message("\n--- Identifying top disturbed regions and sites ---")

# Top 3 regions per country by GROSS disturbance (total_disturbed_ha)
top_regions_by_gross <- regional_summary %>%
  group_by(COUNTRY) %>%
  slice_max(order_by = total_disturbed_ha, n = 3) %>%
  ungroup()

# Top 3 regions per country by PERCENT disturbance (percent_disturbed)
top_regions_by_percent <- regional_summary %>%
  group_by(COUNTRY) %>%
  slice_max(order_by = percent_disturbed, n = 3) %>%
  ungroup()

# (b) Within each region, which Natura 2000 sites are the most disturbed?
# Calculate percent disturbance at the SITE level first
full_regional_data <- full_regional_data %>%
    mutate(percent_disturbed_site = (disturbed_ha_inside / forest_ha_inside) * 100) %>%
    mutate(percent_disturbed_site = ifelse(forest_ha_inside == 0, 0, percent_disturbed_site))


# Top 3 sites per region by gross disturbance
top_sites_by_gross <- full_regional_data %>%
  group_by(COUNTRY, NAME_1) %>%
  slice_max(order_by = disturbed_ha_inside, n = 3) %>%
  ungroup() %>%
  select(COUNTRY, NAME_1, SITECODE, disturbed_ha_inside, percent_disturbed_site)

message("Ranking complete.")

# --- Part 4: Save All Outputs --- #
message("\n--- Part 4: Saving analysis outputs ---")

# First, save the main data frames to CSV files for easy access later.
message("Saving regional summary to CSV...")
write_csv(regional_summary, regional_summary_csv_path)

message("Saving top disturbed sites summary to CSV...")
# We'll save the more detailed top_sites_by_gross table
write_csv(top_sites_by_gross, top_sites_csv_path)


# Now, create the formatted text file output.
# This is a helper function to make writing to the file cleaner.
save_regional_output <- function(text_to_save) {
  cat(text_to_save, file = region_text_output_path, append = TRUE)
}

# Clear the file if it already exists to start fresh
if (file.exists(region_text_output_path)) file.remove(region_text_output_path)

# Loop through each country to write its summary to the text file
countries <- unique(regional_summary$COUNTRY)

for (country in countries) {
  
  header <- paste0("\n\n======================================================\n",
                   "      ANALYSIS FOR: ", toupper(country), "\n",
                   "======================================================\n")
  save_regional_output(header)
  
  # --- Top Regions by Gross Disturbance ---
  save_regional_output("\n--- Top 3 Regions by Gross Disturbance (Hectares) ---\n")
  
  regions_gross <- top_regions_by_gross %>%
    filter(COUNTRY == country) %>%
    arrange(desc(total_disturbed_ha))
    
  for (i in 1:nrow(regions_gross)) {
    row <- regions_gross[i,]
    output_line <- paste0(
      i, ". ", row$NAME_1, ": ", round(row$total_disturbed_ha, 2), " ha disturbed\n"
    )
    save_regional_output(output_line)
  }
  
  # --- Top Regions by Percent Disturbance ---
  save_regional_output("\n--- Top 3 Regions by Percent of Forest Disturbed ---\n")

  regions_percent <- top_regions_by_percent %>%
    filter(COUNTRY == country) %>%
    arrange(desc(percent_disturbed))
    
  for (i in 1:nrow(regions_percent)) {
    row <- regions_percent[i,]
    output_line <- paste0(
      i, ". ", row$NAME_1, ": ", round(row$percent_disturbed, 2), "% of forest area disturbed\n"
    )
    save_regional_output(output_line)
  }

  # --- Top Sites within those top regions ---
  save_regional_output("\n--- Top Disturbed Sites Within These Regions ---\n")
  
  # Get the list of top regions to report on
  top_region_names <- unique(c(regions_gross$NAME_1, regions_percent$NAME_1))
  
  for (region_name in top_region_names) {
    sites_in_region <- top_sites_by_gross %>%
      filter(NAME_1 == region_name) %>%
      arrange(desc(disturbed_ha_inside))
      
    if(nrow(sites_in_region) > 0) {
      save_regional_output(paste0("  In ", region_name, ":\n"))
      
      for(j in 1:nrow(sites_in_region)) {
        site_row <- sites_in_region[j,]
        site_line <- paste0(
          "    ", j, ". Site ", site_row$SITECODE, ": ",
          round(site_row$disturbed_ha_inside, 2), " ha disturbed (",
          round(site_row$percent_disturbed_site, 1), "%)\n"
        )
        save_regional_output(site_line)
      }
    }
  }
}

message(paste("\n--- Regional analysis complete. Results saved to:", output_dir, "and", region_text_output_path, "---"))

# --- PAN-EUROPEAN ANALYSIS ---
pan_eu_header <- "\n\n======================================================\n"
pan_eu_header <- paste0(pan_eu_header, "      OVERALL TOP 10 REGIONS (PAN-EUROPEAN)\n")
pan_eu_header <- paste0(pan_eu_header, "======================================================\n")
save_regional_output(pan_eu_header)

# Top 10 Regions across all of Europe by Gross Disturbance
save_regional_output("\n--- Top 10 Regions by Gross Disturbance (Hectares) ---\n")
top_10_gross_eu <- regional_summary %>%
  slice_max(order_by = total_disturbed_ha, n = 10)

for (i in 1:nrow(top_10_gross_eu)) {
  row <- top_10_gross_eu[i,]
  output_line <- paste0(
    i, ". ", row$NAME_1, " (", row$COUNTRY, "): ", 
    round(row$total_disturbed_ha, 2), " ha disturbed\n"
  )
  save_regional_output(output_line)
}

# Top 10 Regions across all of Europe by Percent Disturbance
save_regional_output("\n--- Top 10 Regions by Percent of Forest Disturbed ---\n")
top_10_percent_eu <- regional_summary %>%
  # We should add a filter here to ensure we are only looking at regions
  # with a meaningful amount of forest to avoid misleading high percentages
  # from very small areas. Let's set a threshold of 1000 ha of forest.
  filter(total_forest_ha > 1000) %>%
  slice_max(order_by = percent_disturbed, n = 10)

for (i in 1:nrow(top_10_percent_eu)) {
  row <- top_10_percent_eu[i,]
  output_line <- paste0(
    i, ". ", row$NAME_1, " (", row$COUNTRY, "): ", 
    round(row$percent_disturbed, 2), "% of forest area disturbed\n"
  )
  save_regional_output(output_line)
}

