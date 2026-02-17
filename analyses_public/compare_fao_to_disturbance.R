# analyze_forestry_trends.R
#
# Purpose: Analyze broad forestry trends (Production, Trade) for EU Member States
#          and the EU-27 aggregate using FAOSTAT data.
#          Calculates Sen's Slope to quantify intensification/mobilization rates.

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(trend)   # For Sen's slope and Mann-Kendall test
library(scales)  # For number formatting

# --- 1. Load and Clean Data ---

input_path <- "imports/FAOSTAT_data_en_12-5-2025.csv"

if (!file.exists(input_path)) {
  stop(paste("File not found:", input_path))
}

# Load data
fao_data <- read_csv(input_path, show_col_types = FALSE)

# Basic cleaning
clean_data <- fao_data %>%
  select(Area, Element, Item, Year, Unit, Value) %>%
  # Filter out non-EU aggregates if they exist (keep specific countries + EU27)
  filter(Area != "World") %>% 
  # Ensure Value is numeric
  mutate(Value = as.numeric(Value)) %>%
  filter(!is.na(Value)) %>%
  # IMPORTANT: Arrange by Year to ensure time series calculations are correct
  arrange(Area, Item, Element, Year)

# --- Dynamic Area Detection ---
# Identify exactly how "European Union" is spelled in this specific file
# This handles variations like "European Union (27)" vs "European Union (28)" or hidden spaces.
eu_candidates <- unique(clean_data$Area[grepl("European Union", clean_data$Area)])

if (length(eu_candidates) > 0) {
  # Use the first match found (usually "European Union (27)")
  eu_area_name <- eu_candidates[1] 
  message(paste("Detected EU Aggregate Name as:", eu_area_name))
} else {
  eu_area_name <- "European Union (27)" # Fallback
  message("WARNING: No Area matching 'European Union' found. Using default.")
}

# Debug: Print found areas to confirm what is in the file
found_areas <- unique(clean_data$Area)
message("--- Data Loaded ---")
message(paste("Found", length(found_areas), "unique Areas."))
message("Areas: ", paste(head(found_areas, 10), collapse = ", "), 
        ifelse(length(found_areas) > 10, "...", ""))

# Create Output Directory
dir.create("plots/fao_trends", recursive = TRUE, showWarnings = FALSE)

# --- 2. Define Analysis Functions ---

# Function to calculate Stats (Mann-Kendall P-value AND Sen's Slope)
calc_trend_stats <- function(years, values) {
  # Need at least 3 points for a meaningful trend test
  if (length(values) < 3) {
    return(tibble(MK_p_val = NA_real_, SS_p_val = NA_real_, Sens_Slope = NA_real_))
  }
  
  # Create a time series object
  ts_data <- ts(values, start = min(years), frequency = 1)
  
  # 1. Run Mann-Kendall Trend Test
  mk_res <- tryCatch(mk.test(ts_data), error = function(e) NULL)
  
  if (is.null(mk_res)) {
     return(tibble(MK_p_val = NA_real_, SS_p_val = NA_real_, Sens_Slope = NA_real_))
  }
  
  mk_p <- mk_res$p.value
  
  # 2. Run Sen's Slope ONLY if Mann-Kendall is positive/significant (p < 0.05)
  ss_p <- NA_real_
  slope <- NA_real_
  
  # Using 0.05 as the standard alpha for significance
  if (!is.na(mk_p) && mk_p < 0.05) {
    ss_res <- tryCatch(sens.slope(ts_data), error = function(e) NULL)
    if (!is.null(ss_res)) {
      ss_p <- ss_res$p.value
      slope <- ss_res$estimates
    }
  }
  
  # Return a single-row tibble with all stats
  return(tibble(MK_p_val = mk_p, SS_p_val = ss_p, Sens_Slope = slope))
}

# --- 3. Run Analysis by Item and Element ---

# Get unique combinations of Item (Product) and Element (Metric: Production/Import/Export)
unique_combos <- clean_data %>%
  distinct(Item, Element, Unit)

# Initialize lists to store results
stats_results <- list()

message("\nStarting analysis of forestry trends...")

for (i in 1:nrow(unique_combos)) {
  
  curr_item <- unique_combos$Item[i]
  curr_element <- unique_combos$Element[i]
  curr_unit <- unique_combos$Unit[i]
  
  message(paste("Processing:", curr_item, "-", curr_element))
  
  # Filter data for this specific combination
  plot_data <- clean_data %>%
    filter(Item == curr_item, Element == curr_element)
  
  if (nrow(plot_data) == 0) next
  
  # --- A. EU-27 Analysis (Individual Plots) ---
  eu_data <- plot_data %>% filter(Area == eu_area_name)
  
  if (nrow(eu_data) > 0) {
    p_eu <- ggplot(eu_data, aes(x = Year, y = Value)) +
      geom_line(color = "#0072B2", linewidth = 1.2) +
      geom_point(size = 2, color = "#0072B2") +
      geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "gray50") +
      scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
      labs(
        title = paste("EU-27 Trend:", curr_item),
        subtitle = paste(curr_element, "(", curr_unit, ")"),
        x = "Year",
        y = paste("Value", curr_unit)
      ) +
      theme_minimal()
    
    # Save EU plot
    safe_filename <- gsub("[^A-Za-z0-9]", "_", paste("EU", curr_item, curr_element, sep = "_"))
    ggsave(paste0("plots/fao_trends/", safe_filename, ".jpg"), p_eu, width = 8, height = 5)
  }
  
  # --- B. Member State Analysis (Facet Plot) ---
  # Exclude EU aggregate for the country plots
  ms_data <- plot_data %>% 
    filter(Area != eu_area_name) %>%
    filter(!is.na(Value))
  
  if (nrow(ms_data) > 0) {
    p_ms <- ggplot(ms_data, aes(x = Year, y = Value)) +
      geom_line(color = "#D55E00") +
      facet_wrap(~Area, scales = "free_y") + # free_y allows each country to have its own scale
      scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
      labs(
        title = paste("Member State Trends:", curr_item),
        subtitle = paste(curr_element, "- Note: Y-axis scales differ by country"),
        x = "Year",
        y = curr_unit
      ) +
      theme_minimal() +
      theme(
        strip.text = element_text(size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6)
      )
    
    # Save MS plot
    safe_filename_ms <- gsub("[^A-Za-z0-9]", "_", paste("MS", curr_item, curr_element, sep = "_"))
    ggsave(paste0("plots/fao_trends/", safe_filename_ms, ".jpg"), p_ms, width = 12, height = 10)
  }
  
  # --- C. Calculate Statistics (Mann-Kendall + Sen's Slope) ---
  # We calculate this for ALL areas (including EU-27)
  
  stats <- plot_data %>%
    group_by(Area) %>%
    summarise(
      Item = curr_item,
      Element = curr_element,
      Unit = curr_unit,
      Start_Year = min(Year),
      End_Year = max(Year),
      Min_Value = min(Value), 
      Max_Value = max(Value), 
      Total_Change_Percent = (last(Value) - first(Value)) / first(Value) * 100,
      # Call the custom function which returns a tibble of stats
      trend_stats = list(calc_trend_stats(Year, Value)),
      .groups = "drop"
    ) %>%
    # Unnest the stats columns (MK_p_val, SS_p_val, Sens_Slope)
    unnest(trend_stats)
  
  stats_results[[length(stats_results) + 1]] <- stats
}

# --- 4. Generate Master EU Combination Plot ---
# We calculate this directly from clean_data to ensure it works even if the loop filtering is strict

message("\nGenerating Master EU Combination Plot...")

# Filter specifically for the EU area detected earlier
eu_master_data <- clean_data %>% 
  filter(Area == eu_area_name)

if (nrow(eu_master_data) > 0) {
  
  combined_eu_df <- eu_master_data %>%
    # Create a nice label for the facet strips combining Item and Element
    mutate(Facet_Label = paste0(Item, "\n(", Element, " - ", Unit, ")"))
  
  p_combined <- ggplot(combined_eu_df, aes(x = Year, y = Value)) +
    geom_line(color = "#0072B2", linewidth = 1) +
    # Facet by the Item+Element label. "scales = free_y" is critical here.
    facet_wrap(~Facet_Label, scales = "free_y", ncol = 3) + 
    scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
    labs(
      title = "All EU-27 Forestry Trends",
      subtitle = paste("Overview of Production and Trade Metrics (2000-2023) -", eu_area_name),
      x = "Year",
      y = "Value (See panel title for units)"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 8, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(size = 7)
    )
  
  ggsave("plots/fao_trends/MASTER_EU_TRENDS.jpg", p_combined, width = 14, height = 12)
  message("Master plot saved successfully.")
  
} else {
  message("WARNING: No data found for Master EU Plot. Check area name detection.")
}

# --- 5. Save Statistical Summary ---

if (length(stats_results) > 0) {
  full_stats <- bind_rows(stats_results) %>%
    # Sort by Slope descending (highest growth), putting NAs at the bottom
    arrange(Item, Element, desc(Sens_Slope)) %>%
    # Round all numeric columns to 2 decimal places
    mutate(across(where(is.numeric), ~ round(., 2)))
  
  write_csv(full_stats, "csvs/forestry_trends_summary.csv")
  message("Statistical summary saved to 'csvs/forestry_trends_summary.csv'")
  
  # Print top 10 fastest growing trends (where slope is significant)
  message("\n--- Top 10 Fastest Growing Trends (Significant Only) ---")
  print(head(full_stats %>% 
               filter(!is.na(Sens_Slope)) %>% 
               select(Area, Item, Element, MK_p_val, Sens_Slope), 10))
} else {
  message("No data found to analyze.")
}

message("\nAnalysis Complete. Check plots or csvs for graphs.")