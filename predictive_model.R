# Title: Predictive Model for Forest Disturbance in Natura 2000 Sites
# Author: Gemini
# Date: 2025-10-06
# Description: This script builds a Random Forest model to predict high-disturbance
#              Natura 2000 sites based on geographic and site-specific factors.

# --- 1. SETUP: Load necessary libraries ---
# If you don't have these packages installed, run: install.packages(c("dplyr", "randomForest", "caret"))
library(dplyr)
library(randomForest)
library(caret)
library(knitr)

# --- 2. DATA LOADING AND PREPARATION ---

# Load the core datasets from your analysis outputs.
# This file contains the disturbance data for each N2K site.
disturbance_by_site <- read.csv("analysis_outputs/full_disturbance_and_forest_area_by_site.csv")

# This file contains data on the 1-km buffer zones around the sites.
buffer_data <- read.csv("analysis_outputs/n2k_and_buffer_forest_area.csv")

# This file contains the region information
region_data <- read.csv("imports/NATURA2000_v2022_WITH_REGIONS.csv")

# Merge the two datasets to combine site and buffer information.
# We'll join them by the SITECODE, which is the unique identifier for each Natura 2000 site.
model_data <- disturbance_by_site %>%
  left_join(buffer_data, by = "SITECODE") %>%
  left_join(region_data, by = "SITECODE") %>%
  rename(n2k_forest_area = forest_ha_inside.x) %>%
  rename(buffer_forest_area = forest_ha_buffer.x) %>%
  rename(total_disturbance_in_site = disturbed_ha_inside) %>%
  rename(buffer_disturbance_area = disturbed_ha_buffer) %>%
  rename(COUNTRY.x = member_state) %>%
  rename(REGION.x = GID_1) %>%
  select(-MS)

str(model_data)
# --- 3. FEATURE ENGINEERING ---

# Create the target variable: 'high_disturbance'
# First, calculate the disturbance rate to normalize by the size of the forest area in the site.
# This prevents larger sites from always being flagged as high disturbance.
model_data <- model_data %>%
  mutate(
    disturbance_rate = ifelse(n2k_forest_area > 0, total_disturbance_in_site / n2k_forest_area, 0)
  )

# Define "high disturbance" sites. We'll use the 75th percentile as the threshold.
# This means sites in the top 25% of disturbance rates will be classified as "high disturbance".
disturbance_threshold <- quantile(model_data$disturbance_rate, 0.75, na.rm = TRUE)

model_data <- model_data %>%
  mutate(
    high_disturbance = as.factor(ifelse(disturbance_rate >= disturbance_threshold, "Yes", "No"))
  )

# Prepare the final dataset for the model
# Select the predictor variables and the target variable.
# We're also filtering out sites with no forest area and removing any rows with missing data.
final_model_data <- model_data %>%
  select(
    high_disturbance,
    SITECODE,
    COUNTRY.x,
    REGION.x,
    n2k_forest_area,
    buffer_disturbance_area
  ) %>%
  rename(
    country = COUNTRY.x,
    region = REGION.x,
    site_forest_area = n2k_forest_area,
    buffer_disturbance = buffer_disturbance_area
  ) %>%
  na.omit() # Remove rows with any missing values

# --- 4. MODEL TRAINING AND EVALUATION ---

# Split the data into training (80%) and testing (20%) sets.
set.seed(123) # for reproducibility
trainIndex <- createDataPartition(final_model_data$high_disturbance, p = .8,
                                  list = FALSE,
                                  times = 1)
training_set <- final_model_data[ trainIndex,]
testing_set  <- final_model_data[-trainIndex,]

# --- MODIFICATION TO ADDRESS CLASS IMBALANCE ---
# Calculate class weights to penalize the model for misclassifying the minority class ('Yes')
# We'll give the 'Yes' class a higher weight.
# A simple way is to use the inverse of the class frequencies.
n_no <- sum(training_set$high_disturbance == "No")
n_yes <- sum(training_set$high_disturbance == "Yes")

# Train the Random Forest model.
# We are predicting 'high_disturbance' based on the other variables.
# ntree is the number of trees to grow in the forest. 500 is a good default.
rf_model <- randomForest(high_disturbance ~ .,
                         data = training_set,
                         ntree = 500,
                         importance = TRUE, 
                         sampsize = c(Yes = n_yes, No = n_yes),
                         strata = training_set$high_disturbance)

# Print the model summary
print(rf_model)

# Make predictions on the testing set
predictions <- predict(rf_model, testing_set)

# Evaluate model performance using a confusion matrix
confusionMatrix <- confusionMatrix(predictions, testing_set$high_disturbance)
print(confusionMatrix)

# --- 5. INTERPRET RESULTS: FEATURE IMPORTANCE ---

# Plot the importance of each variable in the model.
# This helps us understand which factors are the most powerful predictors of high disturbance.
# The 'MeanDecreaseGini' indicates how much a variable contributes to the purity of the nodes in the random forest.
# Higher values mean the variable is more important.
png("plots/feature_importance_plot.png", width = 800, height = 600)
varImpPlot(rf_model, main = "Feature Importance for Predicting High Disturbance")
dev.off()

cat("\nAnalysis complete. A plot named 'feature_importance_plot.png' has been saved to your directory.\n")
cat("This plot shows the most important variables for predicting forest disturbance.\n")


### REDO PART 4 WITH FULL DATA
# --- 4. MODEL TRAINING (using the entire dataset for the final model) ---

# For generating a final prediction list, it's best practice to train the model
# on ALL available data so it learns from every possible example.
set.seed(123) # for reproducibility

# Use balanced sampling
n_yes <- sum(final_model_data$high_disturbance == "Yes")

# --- MODIFICATION TO PRESERVE SITECODE ---
# We will use the x, y interface for randomForest instead of the formula.
# This prevents non-predictor columns (like SITECODE) from being dropped.

# 1. Separate predictors (x) from the response variable (y)
predictors <- final_model_data %>%
  select(country, region, site_forest_area, buffer_disturbance)

response <- final_model_data$high_disturbance

# 2. Train the model using the x, y interface
final_rf_model <- randomForest(x = predictors,
                               y = response,
                               ntree = 500,
                               importance = TRUE,
                               sampsize = c(Yes = n_yes, No = n_yes),
                               strata = response) # Use the response vector for strata

cat("Final model has been trained on the full dataset.\n")


# Print the model summary
print(final_rf_model)
importance(final_rf_model)
# Make predictions on the testing set
predictions <- predict(final_rf_model, testing_set)

# Evaluate model performance using a confusion matrix
confusionMatrix <- confusionMatrix(predictions, testing_set$high_disturbance)
print(confusionMatrix)


# --- 5. GENERATE PREDICTIONS FOR ALL SITES ---

# Use the trained model to predict on the 'predictors' data frame ONLY.
# This ensures we only use the columns the model was trained on.
risk_scores <- predict(final_rf_model, predictors, type = "prob")
predicted_class <- predict(final_rf_model, predictors)

# Combine the original data (which has SITECODE) with the new prediction columns.
predictions_output <- final_model_data %>%
  mutate(
    predicted_risk_score = risk_scores[, "Yes"],
    predicted_class = predicted_class
  )

# Arrange the final list from highest risk to lowest risk.
prioritized_risk_list <- predictions_output %>%
  arrange(desc(predicted_risk_score))


head(prioritized_risk_list)

  # --- 6. SAVE THE RESULTS ---

# Save this prioritized list to a new CSV file.
write.csv(prioritized_risk_list, "analysis_outputs/prioritized_disturbance_risk_list.csv", row.names = FALSE)

cat("\nPrediction complete.\n")
cat("A new file named 'prioritized_disturbance_risk_list.csv' has been created.\n")
cat("This file contains all Natura 2000 sites, ranked by their predicted risk of high disturbance.\n")


# --- 5. SUMMARIZE BY REGION (NEW SECTION) ---
# We can also group by region to get a more granular view.
region_summary <- prioritized_risk_list %>%
  group_by(country, region) %>%
  summarise(
    total_sites = n(),
    high_risk_sites = sum(predicted_class == "Yes"),
    .groups = 'drop'
  ) %>%
  mutate(
    percent_high_risk = (high_risk_sites / total_sites) * 100
  ) %>%
  filter(total_sites > 10) %>% # Optional: filter for regions with a meaningful number of sites
  arrange(desc(percent_high_risk))

  print(region_summary)

  # --- 6. DISPLAY AND SAVE REGIONAL RESULTS (NEW SECTION) ---

# Print the regional summary table to the console.
cat("\n\nSummary of Predicted High-Risk Sites by Region (top 20):\n\n")
print(head(region_summary, 20))

# Save the regional summary table to a new CSV file.
write.csv(region_summary, "analysis_outputs/region_risk_summary.csv", row.names = FALSE)

cat("\n\nRegional summary complete.\n")
cat("A new file named 'region_risk_summary.csv' has been created.\n")



### Trying to understand RF outputs
library(randomForestExplainer)
explain_forest(final_rf_model)