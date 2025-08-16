# Install packages if you haven't already
# Load libraries
library(sf)            # For vector data manipulation
library(terra)         # For raster data processing
library(dplyr)         # For data wrangling
library(exactextractr) # For fast zonal statistics
library(trend)         # For Mann-Kendall and Sen's Slope tests
library(tidyr)         # For tidyr::unnest

# --- Load datasets --- #
n2k_filepath <- "D:/GIS/Data/Raw Files/2024-05 Natura 2000 Europe/eea_v_3035_100_k_natura2000_p_2022_v01_r00/Natura2000_end2022.gpkg" # nolint: line_length_linter.
st_layers(n2k_filepath)
n2k_layername <- "NaturaSite_polygon"

gfc_filepath <- "D:/GIS/Projects/N2k-EU/Europe_LossYear_2023_v1_1.tif"

effis_filepath <- "D:/GIS/Data/Raw Files/2024-05 Natura 2000 Europe/effis-2024-05/modis.ba.poly.shp" # nolint: line_length_linter.

# --- Load vector data and set up CRS --- #
n2k_sites <- st_read(n2k_filepath, layer = n2k_layername)
effis_fires <- st_read(effis_filepath)

# Set the target CRS for the analysis (an equal-area projection for Europe)
target_crs <- "EPSG:3035"
n2k_sites <- st_transform(n2k_sites, crs = target_crs)
effis_fires <- st_transform(effis_fires, crs = target_crs)

# --- Open large raster without loading into memory --- #
# This creates a pointer to the file on disk, keeping memory usage low.
gfc_lossyear_merged <- rast(gfc_filepath)
gfc_crs <- crs(gfc_lossyear_merged) # Get the raster's original CRS

# --- Define a function to process a single N2K site --- #
# This function will perform all heavy processing on a small chunk of data.
process_site <- function(site_polygon, gfc_raster, fires_data, target_crs) {
  tryCatch({
    # 1. Get the site's bounding box and transform it to the raster's CRS.
    #    This defines the small area we need to read from the large raster file.
    site_bbox_target_crs <- st_bbox(site_polygon)
    site_bbox_gfc_crs <- st_bbox(
                                 st_transform(
                                              st_as_sfc(site_bbox_target_crs),
                                              crs = gfc_crs))

    # 2. Crop the large GFC raster to the site's extent. This reads a small
    #    chunk from disk into memory.
    gfc_cropped <- crop(gfc_raster, ext(site_bbox_gfc_crs))

    # If the site is outside the raster's extent, the crop will be empty.
    if (ncell(gfc_cropped) == 0) {
      return(NULL)
    }

    # 3. Reproject only the small, cropped raster. This is now a fast,
    #    memory-friendly operation.
    gfc_cropped_proj <- project(gfc_cropped, target_crs, method = "near")

    # 4. Find only the fires that intersect with the current site.
    intersecting_fires <- st_filter(
                                    fires_data,
                                    site_polygon,
                                    .predicate = st_intersects)

    # 5. Mask out fire disturbances from the cropped & projected raster.
    gfc_no_fire <- gfc_cropped_proj
    if (nrow(intersecting_fires) > 0) {
      # Mask using all relevant fires at once, more efficient than looping.
      gfc_no_fire <- mask(gfc_no_fire, intersecting_fires, # nolint: object_name_linter.
                          inverse = TRUE, updatevalue = NA)
    }

    # 6. Extract raw pixel values for the single site.
    # This returns a list containing one data frame of pixel values.
    disturbance_pixels <- exact_extract(gfc_no_fire, site_polygon,
                                        include_cols = "SITECODE")

    # If no pixels intersect, return NULL
    if (length(disturbance_pixels) == 0 || nrow(disturbance_pixels[[1]]) == 0) {
        return(NULL)
    }

    # Extract the data frame and summarize to get pixel counts for each
    # event year (value=0 represents undisturbed pixels).
    disturbance_stats <- disturbance_pixels[[1]] %>%
        group_by(SITECODE, value) %>%
        summarize(pixel_count = sum(coverage_fraction), .groups = 'drop')
    return(disturbance_stats)
  }, error = function(e) {
    # If any error occurs for a site, print it and return NULL
    message("Error processing site ", site_polygon$SITECODE, ": ", e$message)
    return(NULL)
  })
}

# --- Loop through all sites and process them one by one --- #
# Pre-allocating the list is more efficient. Set to total_sites for the full run.
results_list <- vector("list", length = 10)
total_sites <- nrow(n2k_sites)

for (i in 1:total_sites) {
  # Print progress
  message(paste0("Processing site ", i, 
                 " of ", total_sites, 
                 " (", n2k_sites$SITECODE[i], ")"))

  # Get the single site polygon
  current_site <- n2k_sites[i, ]

  # Process the site
  site_result <- process_site(current_site,
                              gfc_lossyear_merged,
                              effis_fires,
                              target_crs)

  # Add the result to our list
  if (!is.null(site_result)) {
    results_list[[i]] <- site_result
  }
}


# --- Process the results into a clean table ---
# Combine the list of results into a single data frame
disturbance_df_n2k <- bind_rows(results_list)

disturbance_df_n2k_final <- disturbance_df_n2k %>%
  rename(loss_year = value) %>%
  # The GFC dataset codes years as 1-23 for 2001-2023
  mutate(year = loss_year + 2000) %>%
  # Convert pixel count to hectares (assuming 30m x 30m Hansen pixels)
  mutate(disturbed_ha = pixel_count * (30 * 30) / 10000)

# View the final tidy data frame
head(disturbance_df_n2k_final)
str(disturbance_df_n2k_final)
write.csv(disturbance_df_n2k_final, 
          "n2k_disturbance_by_year_no_fire.csv", 
          row.names = FALSE)