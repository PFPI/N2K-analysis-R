# Title: LASSO Predictive Model for Forest Disturbance in Natura 2000 Sites
# Author: Gemini (adapted from user's Random Forest script)
# Date: 2025-10-09
# Description: This script builds a LASSO regularized logistic regression model to predict
#              high-disturbance Natura 2000 sites and creates a prioritized risk list.

# --- 1. SETUP: Load necessary libraries ---
# If you don't have these packages installed, run: install.packages(c("dplyr", "glmnet", "caret"))
library(dplyr)
library(glmnet) # For LASSO modeling
library(caret)  # For the confusion matrix and data splitting

# --- 2. DATA LOADING AND PREPARATION (This section is identical to your script) ---

# Load the core datasets from your analysis outputs.
disturbance_by_site <- read.csv("analysis_outputs/full_disturbance_and_forest_area_by_site.csv")
buffer_data <- read.csv("analysis_outputs/n2k_and_buffer_forest_area.csv")
region_data <- read.csv("imports/NATURA2000_v2022_WITH_REGIONS.csv")

# Merge the datasets
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

# --- 3. FEATURE ENGINEERING (This section is identical to your script) ---

# Create the target variable: 'high_disturbance'
model_data <- model_data %>%
  mutate(
    disturbance_rate = ifelse(n2k_forest_area > 0, total_disturbance_in_site / n2k_forest_area, 0)
  )

disturbance_threshold <- quantile(model_data$disturbance_rate, 0.75, na.rm = TRUE)

model_data <- model_data %>%
  mutate(
    high_disturbance = as.factor(ifelse(disturbance_rate >= disturbance_threshold, "Yes", "No"))
  )

# Prepare the final dataset
final_model_data <- model_data %>%
  select(
    high_disturbance,
    total_disturbance_in_site,
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
  na.omit()

# --- 4. LASSO MODEL TRAINING AND EVALUATION ---

### 4a. Data Preparation for glmnet ###
# LASSO requires a numeric matrix for predictors (x) and a vector for the response (y).
# model.matrix() will automatically create dummy variables for factors like 'country' and 'region'.
predictors_x <- model.matrix(high_disturbance ~ country + region + site_forest_area + buffer_disturbance, data = final_model_data)[, -1]
response_y <- final_model_data$high_disturbance

### 4b. Train/Test Split ###
set.seed(123)
trainIndex <- createDataPartition(response_y, p = .8, list = FALSE, times = 1)

# Split both the predictor matrix and the response vector
training_x <- predictors_x[trainIndex, ]
testing_x  <- predictors_x[-trainIndex, ]

training_y <- response_y[trainIndex]
testing_y  <- response_y[-trainIndex]

### 4c. Cross-Validation to Find Best Lambda ###
# We use cv.glmnet to find the optimal penalty strength (lambda) for the model.
# alpha = 1 specifies LASSO. family = "binomial" is for binary classification.
set.seed(123)
cv_lasso_model <- cv.glmnet(training_x, training_y, alpha = 1, family = "binomial")

# Plot the cross-validation curve
png("plots/lasso_cv_plot.png", width = 800, height = 600)
plot(cv_lasso_model)
dev.off()
cat("\nAnalysis running... A plot named 'lasso_cv_plot.png' has been saved.\n")


### 4d. Model Evaluation on Testing Set ###
# Get the best lambda value (lambda.1se is often preferred for a simpler model)
best_lambda <- cv_lasso_model$lambda.1se

# Make predictions on the testing set
# type = "class" gives the "Yes"/"No" prediction
predictions <- predict(cv_lasso_model, s = best_lambda, newx = testing_x, type = "class")

# Evaluate performance with a confusion matrix
confusionMatrix <- confusionMatrix(as.factor(predictions), testing_y)
print(confusionMatrix)


# --- 5. INTERPRET RESULTS: VIEWING COEFFICIENTS ---

# The "feature importance" in LASSO comes from which variables were NOT shrunk to zero.
# This tells you which predictors the model found most useful.
lasso_coeffs <- coef(cv_lasso_model, s = best_lambda)
print("--- LASSO Model Coefficients (Non-Zero Predictors) ---")
print(lasso_coeffs[which(lasso_coeffs != 0),])


# --- 6. RETRAIN ON FULL DATA & GENERATE FINAL RISK LIST ---

# Train the final model on ALL available data to get the best possible coefficients.
set.seed(123)
final_lasso_model <- cv.glmnet(predictors_x, response_y, alpha = 1, family = "binomial")
final_best_lambda <- final_lasso_model$lambda.1se

cat("Final LASSO model has been trained on the full dataset.\n")

# Generate predictions for ALL sites using the final model
# type = "response" gives the probability of the "Yes" class (our risk score)
risk_scores <- predict(final_lasso_model, s = final_best_lambda, newx = predictors_x, type = "response")
# type = "class" gives the final predicted class
predicted_class <- predict(final_lasso_model, s = final_best_lambda, newx = predictors_x, type = "class")

# Combine the original data (with SITECODE) with the new predictions
predictions_output <- final_model_data %>%
  mutate(
    predicted_risk_score = as.vector(risk_scores),
    predicted_class = as.vector(predicted_class)
  )

# Arrange the final list from highest risk to lowest risk
prioritized_risk_list <- predictions_output %>%
  arrange(desc(predicted_risk_score))

# --- 7. SUMMARIZE AND SAVE RESULTS (This section is identical to your script) ---

# Save the prioritized list to a new CSV file
write.csv(prioritized_risk_list, "analysis_outputs/prioritized_disturbance_risk_list_LASSO.csv", row.names = FALSE)
cat("\nPrediction complete.\n")
cat("A new file named 'prioritized_disturbance_risk_list_LASSO.csv' has been created.\n")

# Summarize by region
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
  filter(total_sites > 10) %>%
  arrange(desc(percent_high_risk))

# Print and save the regional summary
cat("\n\nSummary of Predicted High-Risk Sites by Region (top 20) using LASSO:\n\n")
print(head(region_summary, 20))

write.csv(region_summary, "analysis_outputs/region_risk_summary_LASSO.csv", row.names = FALSE)
cat("\n\nRegional summary complete.\n")
cat("A new file named 'region_risk_summary_LASSO.csv' has been created.\n")