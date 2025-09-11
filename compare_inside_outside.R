# compare_inside_outside.R
#
# Purpose: This script compares forest disturbance rates inside Natura 2000
#          sites versus in the 1km buffer zones surrounding them. It directly
#          implements the statistical tests described in the project's
#          methods document.
#
# The process is as follows:
# 1. Calculate the total forested area (c. 2000) for each N2K site and its
#    1km buffer ring. This is a necessary, one-time calculation.
# 2. Load the annual disturbance data for both zones (from start.R and
#    buffer_analysis.R).
# 3. Merge all data sources into a single analysis-ready data frame.
# 4. Perform pan-European analysis:
#    - Chi-square test on disturbed vs. undisturbed forest area.
#    - Mann-Kendall trend tests on the annual disturbance time series.
# 5. Perform Member State level analysis, running the same tests for each
#    country.
# 6. Save and print summary results.

# --- Load Libraries --- #
library(sf)
library(terra)
library(dplyr)
library(readr)
library(tidyr)
library(exactextractr)
library(trend) # For mk.test and sens.slope

# --- Configuration --- #

# Input data paths
n2k_disturbance_path <- "exports/n2k_disturbance_by_year_no_fire.csv"
buffer_disturbance_path <- "exports/n2k_buffer_disturbance_by_year_no_fire.csv"
n2k_gpkg_path <- "D:/GIS/Data/Raw Files/2024-05 Natura 2000 Europe/eea_v_3035_100_k_natura2000_p_2022_v01_r00/Natura2000_end2022.gpkg" # nolint
n2k_layername <- "NaturaSite_polygon"

# Path to the GFC lossyear raster. This will be used to define forest area.
# Any non-NA pixel (values 0-23) is considered forest c. 2000.
gfc_lossyear_path <- "D:/GIS/Projects/N2k-EU/Europe_LossYear_2023_v1_1.tif"

target_crs <- "EPSG:3035"

# Output directories
output_dir <- "analysis_outputs"
if (!dir.exists(output_dir)) dir.create(output_dir)
rdata_dir <- "rdata"
if (!dir.exists(rdata_dir)) dir.create(rdata_dir)

# Path for the intermediate forest area file
forest_area_output_path <- file.path(output_dir, "n2k_and_buffer_forest_area.csv")


# --- Part 1: Calculate Total Forest Area (The Missing Piece) --- #

# This section is refactored to be resumable and to save progress incrementally.
message("--- Part 1: Calculating/Verifying Total Forest Area ---")

# Load N2K sites first, as we'll need them regardless.
message("Loading N2K sites for area calculation...")
n2k_sites_raw <- st_read(n2k_gpkg_path, layer = n2k_layername)
## Transform to the same CRS and validate
n2k_sites_all <- st_transform(n2k_sites_raw, crs = target_crs) %>%
  st_make_valid() %>%
  select(SITECODE) # Only keep the SITECODE and geometry

# Determine which sites need processing.
sites_to_process <- n2k_sites_all
if (file.exists(forest_area_output_path)) {
  message("Found existing forest area file. Checking for unprocessed sites...")
  processed_areas <- read_csv(forest_area_output_path, show_col_types = FALSE)
  processed_sitecodes <- unique(processed_areas$SITECODE)

  sites_to_process <- n2k_sites_all %>%
    filter(!SITECODE %in% processed_sitecodes)

  message(paste(length(processed_sitecodes), "sites already processed."))
  message(paste(length(sites_to_process), "sites need to be processed."))
} else {
  message("Forest area file not found. Will process all sites from scratch.")
  # Write the header for the new CSV file
  write_csv(
    tibble(SITECODE = character(), forest_ha_inside = double(), forest_ha_buffer = double()),
    forest_area_output_path
  )
}

total_to_process <- nrow(sites_to_process)

if (total_to_process > 0) {
  message(paste("Processing", total_to_process, "remaining sites. This may take a long time."))
  message("Using lossyear raster to define forest: any non-NA pixel is considered forest.")

  # Load the large raster data only if we need to process something.
  gfc_lossyear <- rast(gfc_lossyear_path)
  gfc_crs <- crs(gfc_lossyear)

  # Define a function to calculate forest area for a given polygon
  # This function uses the lossyear raster, where any non-NA pixel is forest.
  calculate_forest_area <- function(polygon, raster, target_crs) {
    tryCatch({
      # Crop raster to polygon extent for efficiency
      poly_bbox_native_crs <- st_bbox(st_transform(st_as_sfc(st_bbox(polygon)), crs = crs(raster))) # nolint
      raster_cropped <- crop(raster, ext(poly_bbox_native_crs))
      if (ncell(raster_cropped) == 0) return(0)

      # Reproject the small cropped piece using nearest neighbor for categorical data
      raster_proj <- project(raster_cropped, target_crs, method = "near")

      # Create a binary forest mask (1 = forest, 0/NA = not forest)
      # Any pixel that is not NA in the lossyear raster was forest in 2000.
      forest_mask <- !is.na(raster_proj)

      # Use exact_extract to count forested pixels covered by the polygon
      # 'sum' will sum the values of the mask (1s) weighted by coverage
      pixel_sum <- exact_extract(forest_mask, polygon, fun = "sum")

      # Convert pixel count to hectares (assuming 30m GFC resolution)
      forest_ha <- pixel_sum * (30 * 30) / 10000
      return(forest_ha)
    }, error = function(e) {
      message("Error calculating forest area for site ", polygon$SITECODE, ": ", e$message) # nolint
      return(0)
    })
  }

  # Loop through only the sites that need processing
  for (i in 1:total_to_process) {
    current_site <- sites_to_process[i, ]
    sitecode <- current_site$SITECODE

    # More frequent progress reporting
    if (i %% 25 == 0 || i == 1 || i == total_to_process) {
      message(paste0("Processing site ", i, " of ", total_to_process, " (", sitecode, ")"))
    }

    # Calculate forest area INSIDE the site
    area_inside <- calculate_forest_area(current_site, gfc_lossyear, target_crs)

    # Create buffer ring and calculate forest area in the BUFFER
    area_buffer <- 0 # Default value
    tryCatch({
        current_site_buffered <- st_buffer(current_site, dist = 1000)
        current_buffer_ring <- st_difference(current_site_buffered, st_geometry(current_site))

        if (nrow(current_buffer_ring) > 0 && !st_is_empty(current_buffer_ring)) {
          area_buffer <- calculate_forest_area(current_buffer_ring, gfc_lossyear, target_crs)
        }
    }, error = function(e) {
        message(paste("Warning: Could not create or process buffer for site", sitecode, "-", e$message))
        # area_buffer remains 0
    })

    # Create a one-row tibble with the results
    result_tibble <- tibble(
      SITECODE = sitecode,
      forest_ha_inside = area_inside,
      forest_ha_buffer = area_buffer
    )

    # Append the result to the CSV file immediately
    write_csv(result_tibble, forest_area_output_path, append = TRUE)
  }
  message("--- Forest area calculation complete. Results saved. ---")
} else {
  message("All sites have been processed. Moving to analysis.")
}


### IN HERE I HAVE TO MODIFY # Aggregate total disturbance per site


# --- Part 2: Load and Aggregate All Data --- #
message("\n--- Part 2: Loading and preparing disturbance data. ---")

# Load the complete data for the rest of the script.
forest_areas <- read_csv(forest_area_output_path, show_col_types = FALSE)

# Load disturbance data
disturbance_inside <- read_csv(n2k_disturbance_path, show_col_types = FALSE)
disturbance_buffer <- read_csv(buffer_disturbance_path, show_col_types = FALSE)

# Aggregate total disturbance per site
total_disturbance_inside <- disturbance_inside %>%
  filter(loss_year != 0) %>%
  group_by(SITECODE) %>%
  summarise(disturbed_ha_inside = sum(disturbed_ha, na.rm = TRUE),
            .groups = "drop")

# Aggregate total disturbance in 1km buffer zones
total_disturbance_buffer <- disturbance_buffer %>%
  filter(year != 2000) %>% ## new, untested. remove if problematic. 
  group_by(SITECODE) %>%
  summarise(disturbed_ha_buffer = sum(disturbed_ha, na.rm = TRUE),
            .groups = "drop")

# Join everything together
analysis_df <- forest_areas %>%
  left_join(total_disturbance_inside, by = "SITECODE") %>%
  left_join(total_disturbance_buffer, by = "SITECODE") %>%
  # Replace NA disturbance with 0
  mutate(
    disturbed_ha_inside = replace_na(disturbed_ha_inside, 0),
    disturbed_ha_buffer = replace_na(disturbed_ha_buffer, 0)
  ) %>%
  # Calculate undisturbed area
  mutate(
    undisturbed_ha_inside = forest_ha_inside - disturbed_ha_inside,
    undisturbed_ha_buffer = forest_ha_buffer - disturbed_ha_buffer
  ) %>%
  # Add member state for aggregation
  mutate(member_state = substr(SITECODE, 1, 2)) %>%
  # Ensure undisturbed is not negative due to minor geometric mismatches
  mutate(
    undisturbed_ha_inside = pmax(0, undisturbed_ha_inside),
    undisturbed_ha_buffer = pmax(0, undisturbed_ha_buffer)
  )

glimpse(analysis_df)
write.csv(analysis_df, "exports/intermediate_cio2-3_analysis.csv")

# --- Part 3: Pan-European Statistical Analysis --- #
message("\n--- Part 3: Performing Pan-European Analysis ---")

analysis_df <- read.csv("exports/intermediate_cio2-3_analysis.csv")
# Aggregate all data to the pan-European level
pan_europe_summary <- analysis_df %>%
  summarise(
    across(where(is.numeric), sum, na.rm = TRUE)
  )

# 1. Chi-square Test
message("\n--- Chi-square Test (Pan-European) ---")
contingency_table_eu <- matrix(
  c(
    pan_europe_summary$disturbed_ha_inside,
    pan_europe_summary$disturbed_ha_buffer,
    pan_europe_summary$undisturbed_ha_inside,
    pan_europe_summary$undisturbed_ha_buffer
  ),
  nrow = 2,
  dimnames = list(Zone = c("Inside", "Buffer"),
                  Status = c("Disturbed", "Undisturbed"))
)

print(contingency_table_eu)
chi_test_eu <- chisq.test(contingency_table_eu)
print(chi_test_eu)

# Calculate and print Odds Ratio
a <- contingency_table_eu[1, 1] # Disturbed Inside
b <- contingency_table_eu[1, 2] # Undisturbed Inside
c <- contingency_table_eu[2, 1] # Disturbed Buffer
d <- contingency_table_eu[2, 2] # Undisturbed Buffer
odds_ratio_eu <- (a / b) / (c / d)
odds_ratio_disturbed <- (c / d) / (a / b)


save_output <- function(msg_string, file_path = "exports/output.txt") {
  cat(msg_string, file = file_path, append = TRUE)
}

save_output("Pan European Analysis: Odds Ratios\n")

save_output(paste0("\tOdds Ratio (Inside vs. Buffer): ",
                   round(odds_ratio_eu, 3), "\n"))

save_output(paste0("\tOdds Ratio (Buffer vs. Inside): ",
                   round(odds_ratio_disturbed, 3), "\n"))

save_output(paste0("Interpretation: The odds of a protected forest ",
                   "being disturbed are ", round(odds_ratio_eu, 3),
                   ". Values less than one indicate a lowered probability.\n"))

if (odds_ratio_eu < 1) {
  save_output(paste0("\tA protected forest is ",
                     round(odds_ratio_eu / (odds_ratio_eu + 1), 3) * 100,
                     " percent less likely to be disturbed",
                     " compared to an unprotected forest.\n"))
} else {
  save_output(paste0("\tA protected forest is ",
                     round(odds_ratio_eu / (odds_ratio_eu + 1), 3) * 100,
                     " percent more likely to be disturbed",
                     " compared to an unprotected forest.\n"))
}


save_output(paste0("Interpretation: The odds of an uprotected forest ",
                   "being disturbed are ", round(odds_ratio_disturbed, 3),
                   ". Values less than one indicate a lowered probability.\n"))

if (odds_ratio_disturbed < 1) {
  save_output(paste0("\tAn unprotected forest is ",
                     round(odds_ratio_disturbed /
                             (odds_ratio_disturbed + 1), 3) * 100,
                     " percent less likely to be disturbed",
                     " compared to a protected forest.\n"))
} else {
  save_output(paste0("\tAn unprotected forest is ",
                     round(odds_ratio_disturbed /
                             (odds_ratio_disturbed + 1), 3) * 100,
                     " percent more likely to be disturbed",
                     " compared to a protected forest.\n"))
}
save_output("\n\n")

message("Odds Ratio and explanation can be found in exports/output.txt")

# 2. Trend Analysis
message("\n--- Trend Analysis (Pan-European) ---")
save_output("Pan European Analysis: Slope Tests\n")


# disturbance_inside <- read_csv(n2k_disturbance_path, show_col_types = FALSE)
# disturbance_buffer <- read_csv(buffer_disturbance_path, show_col_types = FALSE)

# Prepare time series data
ts_inside <- disturbance_inside %>%
  group_by(year) %>%
  summarise(total_ha = sum(disturbed_ha))
ts_buffer <- disturbance_buffer %>%
  group_by(year) %>%
  summarise(total_ha = sum(disturbed_ha))

# Trend test for INSIDE
mk_inside <- mk.test(ts_inside$total_ha)
sen_inside <- sens.slope(ts_inside$total_ha)
save_output("Trend for Disturbance INSIDE N2K sites:\n")
save_output(paste("\tMann-Kendall p-value:", round(mk_inside$p.value, 4)))
save_output(paste("\n\tSen's Slope (ha/year):", round(sen_inside$estimates, 2)))

save_output("\n")
# Trend test for BUFFER
mk_buffer <- mk.test(ts_buffer$total_ha)
sen_buffer <- sens.slope(ts_buffer$total_ha)
save_output("Trend for Disturbance in BUFFER zones:\n")
save_output(paste("\tMann-Kendall p-value:", round(mk_buffer$p.value, 4)))
save_output(paste("\n\tSen's Slope (ha/year):", round(sen_buffer$estimates, 2)))



# --- Part 4: Member State Statistical Analysis --- #
message("\n\n--- Part 4: Performing Member State Analysis ---")
save_output("\n\nMember State Analysis\n")

# Aggregate data by member state
ms_summary <- analysis_df %>%
  group_by(member_state) %>%
  summarise(
    across(where(is.numeric), sum, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Only test states with area in both zones
  filter(forest_ha_inside > 0, forest_ha_buffer > 0) 

 write.csv(ms_summary, "exports/intermediate_step4.csv")
## Originally showing no disturbed HA in buffers in AT or BE.
## Made a few edits to buffer_analysis.R and re-ran to get corrections.

# For efficiency, pre-aggregate the time series data for trend analysis
ms_ts_inside <- disturbance_inside %>%
  mutate(member_state = substr(SITECODE, 1, 2)) %>%
  filter(year != 2000) %>%
  group_by(member_state, year) %>%
  summarise(total_ha = sum(disturbed_ha, na.rm = TRUE), .groups = "drop")

ms_ts_buffer <- disturbance_buffer %>%
  mutate(member_state = substr(SITECODE, 1, 2)) %>%
  group_by(member_state, year) %>%
  summarise(total_ha = sum(disturbed_ha, na.rm = TRUE), .groups = "drop")


# Function to perform tests for one member state
run_ms_analysis <- function(df_row, ts_data_in, ts_data_out) {
  # Chi-square
  tbl <- matrix(c(df_row$disturbed_ha_inside, df_row$disturbed_ha_buffer,
                  df_row$undisturbed_ha_inside, df_row$undisturbed_ha_buffer), nrow = 2)
  # Use suppressWarnings for cases with low counts
  chi_test <- suppressWarnings(chisq.test(tbl))
  
  # Odds Ratio
  a <- tbl[1, 1]; b <- tbl[1, 2]; c <- tbl[2, 1]; d <- tbl[2, 2]
  print(c(df_row$member_state, a, b, c, d))
  odds_ratio <- if (b != 0 && c != 0) (a * d) / (b * c) else NA
  
  # Trend Analysis
  ms_code <- df_row$member_state
  ts_in <- ts_data_in %>% filter(member_state == ms_code)
  ts_out <- ts_data_out %>% filter(member_state == ms_code)
    
  mk_in_p <- if(nrow(ts_in) > 1) mk.test(ts_in$total_ha)$p.value else NA
  sen_in_slope <- if(nrow(ts_in) > 1) sens.slope(ts_in$total_ha)$estimates else NA
  mk_out_p <- if(nrow(ts_out) > 1) mk.test(ts_out$total_ha)$p.value else NA
  sen_out_slope <- if(nrow(ts_out) > 1) sens.slope(ts_out$total_ha)$estimates else NA
  
  return(tibble(
    member_state = ms_code,
    chi2_p_value = chi_test$p.value,
    odds_ratio = odds_ratio,
    mk_p_inside = mk_in_p,
    sens_slope_inside = sen_in_slope,
    mk_p_buffer = mk_out_p,
    sens_slope_buffer = sen_out_slope
  ))
}

# Apply the function to each member state
ms_results <- bind_rows(lapply(
  1:nrow(ms_summary),
  function(i) run_ms_analysis(ms_summary[i, ], ms_ts_inside, ms_ts_buffer)
))

message("\n--- Member State Results Summary ---")
print(ms_results, n = 30)
save_output(paste(capture.output(print(ms_results, n = 30)), collapse = "\n"))


# Add these lines
for (i in 1:nrow(ms_results)) {
  row <- ms_results[i, ]

  # Format the output string
  output_string <- paste0(
    "\n\n--- Analysis for ", row$member_state, " ---\n",
    "Chi-square p-value: ", round(row$chi2_p_value, 4), "\n",
    "Odds Ratio (Inside/Buffer): ", round(row$odds_ratio, 3), "\n",
    "  - Interpretation: For every one instance of disturbance in the buffer, there are ", round(row$odds_ratio, 3), " instances in the site.\n",
    "Trend Inside - p-value: ", round(row$mk_p_inside, 4), ", slope: ", round(row$sens_slope_inside, 2), " ha/yr\n",
    "Trend Buffer - p-value: ", round(row$mk_p_buffer, 4), ", slope: ", round(row$sens_slope_buffer, 2), " ha/yr\n"
  )

  # Save the formatted string to the output file
  save_output(output_string)
}

# --- Part 5: Save Outputs --- #
message("\n\n--- Part 5: Saving analysis outputs ---")

# Save the main analysis dataframe
write_csv(analysis_df, file.path(output_dir, "full_disturbance_and_forest_area_by_site.csv"))
message("Full site-level data saved.")

# Save the member state statistical results
write_csv(ms_results, file.path(output_dir, "member_state_statistical_summary.csv"))
message("Member state statistical summary saved.")

# Save workspace
save.image(file = file.path(rdata_dir, "comparison_analysis_checkpoint.RData"))
message("R workspace saved.")

message("\n--- Comparison analysis complete. ---")