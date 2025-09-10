# buffer_analysis.R
#
# Purpose: This script performs a disturbance analysis on 1km buffers created
#          AROUND Natura 2000 sites. It calculates the area of forest loss
#          (from the Global Forest Change dataset) within these buffer "rings",
#          excluding areas affected by fire (from EFFIS).
#
#          The process is as follows:
#          1. Load Natura 2000 site polygons.
#          2. Create a 1km buffer around each site.
#          3. Subtract the original site area from the buffer to create a "ring".
#          4. For each buffer ring, extract non-fire forest loss pixels.
#          5. Summarize the disturbed area (in hectares) by year for each site's buffer.
#          6. Export the results to a CSV file.

# --- Load libraries --- #
library(sf)            # For vector data manipulation
library(terra)         # For raster data processing
library(dplyr)         # For data wrangling
library(exactextractr) # For fast zonal statistics

# --- Configuration: File Paths --- #
# Ensure these paths are correct for your system.
n2k_filepath <- "D:/GIS/Data/Raw Files/2024-05 Natura 2000 Europe/eea_v_3035_100_k_natura2000_p_2022_v01_r00/Natura2000_end2022.gpkg" # nolint: line_length_linter.
n2k_layername <- "NaturaSite_polygon"

gfc_filepath <- "D:/GIS/Projects/N2k-EU/Europe_LossYear_2023_v1_1.tif"

effis_filepath <- "D:/GIS/Data/Raw Files/2024-05 Natura 2000 Europe/effis-2024-05/modis.ba.poly.shp" # nolint: line_length_linter.

output_dir <- "exports"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
output_csv_path <- file.path(output_dir, "n2k_buffer_disturbance_by_year_no_fire.csv") #nolint

rdata_dir <- "rdata"
if (!dir.exists(rdata_dir)) {
  dir.create(rdata_dir)
}
output_rdata_path <- file.path(rdata_dir, "buffer_analysis_checkpoint.RData")

# --- Load vector data and set up CRS --- #
message("--- Loading vector data ---")
n2k_sites_raw <- st_read(n2k_filepath, layer = n2k_layername)
effis_fires <- st_read(effis_filepath)

# Set the target CRS for the analysis (an equal-area projection for Europe)
target_crs <- "EPSG:3035"
n2k_sites <- st_transform(n2k_sites_raw, crs = target_crs)
effis_fires <- st_transform(effis_fires, crs = target_crs)

# --- Prepare for Analysis --- #
# Ensure geometries are valid before performing buffer and difference operations
# This is an important step to prevent errors during geometric operations.
n2k_sites <- st_make_valid(n2k_sites)


# --- Open large raster without loading into memory --- #
message("--- Loading GFC raster data ---")
gfc_lossyear_merged <- rast(gfc_filepath)
gfc_crs <- crs(gfc_lossyear_merged) # Get the raster's original CRS

# --- Define a function to process a single N2K site buffer --- #
# This function is reused from the original analysis (start.R). It performs
# all heavy processing on a small chunk of data (a single polygon).
process_site_buffer <- function(site_buffer_polygon, gfc_raster, fires_data, target_crs) {
  tryCatch({
    # 1. Get the buffer's bounding box and transform it to the raster's CRS.
    site_bbox_target_crs <- st_bbox(site_buffer_polygon)
    site_bbox_gfc_crs <- st_bbox(
                                 st_transform(
                                              st_as_sfc(site_bbox_target_crs),
                                              crs = gfc_crs))

    # 2. Crop the large GFC raster to the buffer's extent.
    gfc_cropped <- crop(gfc_raster, ext(site_bbox_gfc_crs))

    if (ncell(gfc_cropped) == 0) return(NULL)

    # 3. Reproject only the small, cropped raster.
    gfc_cropped_proj <- project(gfc_cropped, target_crs, method = "near")

    # 4. Find fires that intersect with the current buffer.
    intersecting_fires <- st_filter(
                                    fires_data,
                                    site_buffer_polygon,
                                    .predicate = st_intersects)

    # 5. Mask out fire disturbances from the cropped & projected raster.
    gfc_no_fire <- gfc_cropped_proj
    if (nrow(intersecting_fires) > 0) {
      gfc_no_fire <- mask(gfc_no_fire, intersecting_fires,
                          inverse = TRUE, updatevalue = NA)
    }

    # 6. Extract raw pixel values for the single buffer ring.
    disturbance_pixels <- exact_extract(gfc_no_fire, site_buffer_polygon,
                                        include_cols = "SITECODE")

    if (length(disturbance_pixels) == 0 || nrow(disturbance_pixels[[1]]) == 0) {
        return(NULL)
    }

    # 7. Summarize pixel counts for each disturbance year.
    disturbance_stats <- disturbance_pixels[[1]] %>%
        group_by(SITECODE, value) %>%
        summarize(pixel_count = sum(coverage_fraction), .groups = 'drop')
    return(disturbance_stats)
  }, error = function(e) {
    message("Error processing buffer for site ", site_buffer_polygon$SITECODE, ": ", e$message)
    return(NULL)
  })
}

# --- Loop through all sites, create buffer, and process one by one --- #
# This approach avoids creating all buffers in memory at once by writing
# results to a CSV file incrementally.

# 1. Create the header for our output file first.
# This initializes/overwrites the file.
header <- "SITECODE,year,disturbed_ha"
writeLines(header, output_csv_path)

total_sites <- nrow(n2k_sites)
message(paste("\n--- Starting analysis for", total_sites, "site buffers ---"))
# The results_list is no longer needed.

for (i in 1:total_sites) {
  # Print progress less frequently to avoid cluttering the console
  if (i %% 100 == 0 || i == 1 || i == total_sites) {
      message(paste0("Processing site ", i,
                     " of ", total_sites,
                     " (", n2k_sites$SITECODE[i], ")"))
  }

  # Get the single site polygon
  current_site <- n2k_sites[i, ]

  # --- Create Buffer Ring for the CURRENT site --- #
  # This is the memory-efficient approach.
  # 1. Create a 1km buffer (distance is in meters, as per the CRS)
  current_site_buffered <- st_buffer(current_site, dist = 1000)

  # 2. Create the buffer "ring" by subtracting the original site geometry
  current_buffer_ring <- st_difference(current_site_buffered, st_geometry(current_site))

  # 3. Skip if the resulting ring is empty or invalid
  if (nrow(current_buffer_ring) == 0 || st_is_empty(current_buffer_ring)) {
    next
  }

  site_result <- process_site_buffer(current_buffer_ring,
                                     gfc_lossyear_merged,
                                     effis_fires,
                                     target_crs)

  # Process this single result and append it to the CSV
  if (!is.null(site_result) && nrow(site_result) > 0) {
    # Transform the data for the current site
    processed_result <- site_result %>%
      rename(loss_year = value) %>%
      filter(loss_year > 0) %>%
      mutate(year = loss_year + 2000) %>%
      mutate(disturbed_ha = pixel_count * (30 * 30) / 10000) %>%
      select(SITECODE, year, disturbed_ha)

    # Append the processed result to the CSV file without writing the header again
    if (nrow(processed_result) > 0) {
      write.table(processed_result,
                  output_csv_path,
                  sep = ",",
                  append = TRUE,
                  row.names = FALSE,
                  col.names = FALSE)
    }
  }
}

# --- Save workspace ---
message(paste("--- Saving R workspace to:", output_rdata_path, "---"))
save.image(file = output_rdata_path)

message("\n--- Buffer analysis complete. ---")