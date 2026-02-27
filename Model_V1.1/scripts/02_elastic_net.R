#!/usr/bin/env Rscript

# Elastic Net script for NeuroPathPredict Pipeline V1.1

# Comment out potential initial standalone execution logic if present (e.g., arg parsing, direct calls)
# This script seems to primarily define functions called by run_pipeline.R.
# However, to be safe and consistent, we ensure no top-level execution is triggered on source.

# Example of what might be commented out if it existed at the top or bottom:
# args <- commandArgs(trailingOnly = TRUE)
# if (length(args) != 1) {
#   stop("Usage: Rscript 02_elastic_net.R <subject_id>")
# }
# subject_id_arg <- args[1]
# 
# # Source common functions (if any, assuming handled by run_pipeline or internalized)
# # source(file.path(dirname(sys.frame(1)$ofile), "common_functions.R"))
# 
# # Config loading and logging (if any, assuming handled by run_pipeline or within functions)
# # config_obj <- load_config()
# # init_logging(...)
# 
# # Top-level call (if any)
# # tryCatch({
# #   result <- run_elastic_net(subject_id_arg, config_obj) # This function is defined below
# # }, error = function(e) {
# #   ...
# # })

# Function to load transformed data
load_transformed_data <- function(subject_id, config) {
  # Remove any extra spaces from subject_id
  subject_id <- trimws(subject_id)
  
  # Construct path components
  base_dir <- trimws(config$output$base_dir)
  preprocessing_dir <- trimws(config$output$preprocessing)
  filename <- sprintf("%s_step3_transformed_standardized.csv", subject_id)
  
  # Construct full path
  input_file <- file.path(base_dir, subject_id, preprocessing_dir, filename)
  
  if (!file.exists(input_file)) {
    stop(sprintf("Transformed data file not found: %s", input_file))
  }
  
  data <- readr::read_csv(input_file, show_col_types = FALSE)
  logging::loginfo("Loaded transformed data: %d rows, %d columns", 
                  nrow(data), ncol(data))
  return(data)
}

# Function to prepare data for Elastic Net
prepare_data <- function(data) {
  # Get all BZ and EDT columns
  bz_cols <- grep("^bz_", names(data), value = TRUE)
  edt_cols <- grep("^v1_edt_", names(data), value = TRUE)
  
  # Separate features and target
  X <- as.matrix(data[, c(bz_cols, edt_cols)])
  y <- data$AT8_intensity
  
  # Create train-test split (80-20)
  set.seed(42)  # For reproducibility
  train_idx <- sample(1:nrow(X), size = floor(0.8 * nrow(X)))
  
  X_train <- X[train_idx, ]
  y_train <- y[train_idx]
  X_test <- X[-train_idx, ]
  y_test <- y[-train_idx]
  
  # Center and scale features
  X_train_scaled <- scale(X_train)
  X_test_scaled <- scale(X_test, 
                        center = attr(X_train_scaled, "scaled:center"),
                        scale = attr(X_train_scaled, "scaled:scale"))
  
  return(list(
    X_train = X_train_scaled,
    y_train = y_train,
    X_test = X_test_scaled,
    y_test = y_test,
    feature_names = c(bz_cols, edt_cols)
  ))
}

# Function to calculate model metrics
calculate_metrics <- function(y_true, y_pred) {
  metrics <- list(
    mse = mean((y_true - y_pred)^2),
    rmse = sqrt(mean((y_true - y_pred)^2)),
    mae = mean(abs(y_true - y_pred)),
    r2 = 1 - sum((y_true - y_pred)^2) / sum((y_true - mean(y_true))^2),
    adj_r2 = 1 - (1 - (1 - sum((y_true - y_pred)^2) / sum((y_true - mean(y_true))^2))) * 
             (length(y_true) - 1) / (length(y_true) - 2)  # Using 2 for intercept and slope
  )
  return(metrics)
}

# Function to train Elastic Net model
train_elastic_net <- function(data, config) {
  logging::loginfo("Training Elastic Net model...")
  
  # Prepare data
  prepared_data <- prepare_data(data)
  X_train <- prepared_data$X_train
  y_train <- prepared_data$y_train
  X_test <- prepared_data$X_test
  y_test <- prepared_data$y_test
  feature_names <- prepared_data$feature_names
  
  # Set up cross-validation with progress bar
  logging::loginfo("Performing cross-validation...")
  pb <- txtProgressBar(min = 0, max = config$elastic_net$cv_folds, style = 3)
  cv_fit <- glmnet::cv.glmnet(
    x = X_train,
    y = y_train,
    alpha = config$elastic_net$alpha_grid[1],  # Use first alpha value
    nfolds = config$elastic_net$cv_folds,
    trace.it = 1,
    set.seed = 42
  )
  close(pb)
  
  # Get best lambda
  best_lambda <- cv_fit$lambda.min
  
  # Fit final model
  logging::loginfo("Fitting final model...")
  model <- glmnet::glmnet(
    x = X_train,
    y = y_train,
    alpha = config$elastic_net$alpha_grid[1],  # Use first alpha value
    lambda = best_lambda
  )
  
  # Calculate predictions
  train_predictions <- predict(model, X_train)
  test_predictions <- predict(model, X_test)
  
  # Calculate metrics for both train and test sets
  train_metrics <- calculate_metrics(y_train, train_predictions)
  test_metrics <- calculate_metrics(y_test, test_predictions)
  
  return(list(
    model = model,
    cv_fit = cv_fit,
    train_predictions = train_predictions,
    test_predictions = test_predictions,
    train_metrics = train_metrics,
    test_metrics = test_metrics,
    feature_names = feature_names
  ))
}

# Function to analyze feature importance
analyze_features <- function(model, data) {
  # Get coefficients
  coefs <- as.matrix(coef(model))
  coefs <- coefs[-1, ]  # Remove intercept
  
  # Create feature importance data frame
  importance <- data.frame(
    feature = names(coefs),
    coefficient = as.vector(coefs),
    abs_coefficient = abs(as.vector(coefs))
  )
  
  # Sort by absolute coefficient
  importance <- importance[order(-importance$abs_coefficient), ]
  
  return(importance)
}

# Function to save significant predictors
save_significant_predictors <- function(importance, subject_id, config, threshold = 30) {
  output_dir <- file.path(config$output$base_dir, subject_id, 
                         config$output$elastic_net)
  
  # Get top predictors
  top_predictors <- head(importance[importance$coefficient != 0, ], threshold)
  
  # Save to text file
  predictors_file <- file.path(output_dir, sprintf("%s_predictors_to_review.txt", subject_id))
  sink(predictors_file)
  cat("Significant Predictors for Elastic Net Model\n")
  cat("=========================================\n\n")
  for (i in 1:nrow(top_predictors)) {
    cat(sprintf("%d. %s (Coefficient: %.4f)\n", 
                i, top_predictors$feature[i], top_predictors$coefficient[i]))
  }
  sink()
  
  # Save predictors list as CSV for use in factor analysis
  predictors_csv <- file.path(output_dir, sprintf("%s_significant_predictors.csv", subject_id))
  write.csv(top_predictors, predictors_csv, row.names = FALSE)
  
  return(top_predictors)
}

# Function to create SigPred visualizations
create_sigpred_visualizations <- function(data, top_predictors, subject_id, config) {
  # Create SigPred directory
  sigpred_dir <- file.path(config$output$base_dir, subject_id, 
                          config$output$elastic_net, "SigPred")
  dir.create(sigpred_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Prepare data for visualizations
  pred_data <- data[, c("AT8_intensity", top_predictors$feature)]
  
  # Correlation plot
  cor_matrix <- cor(pred_data, use = "complete.obs")
  png(file.path(sigpred_dir, "correlation_plot.png"), 
      width = 1200, height = 1200, res = 100)
  corrplot::corrplot(cor_matrix, method = "color", type = "upper", 
                    tl.col = "black", tl.srt = 45)
  dev.off()
  
  # Dendrogram
  hc <- hclust(dist(cor_matrix))
  png(file.path(sigpred_dir, "dendrogram.png"), 
      width = 1200, height = 800, res = 100)
  plot(hc, main = "Hierarchical Clustering of Predictors",
       xlab = "", sub = "", cex = 0.8)
  dev.off()
  
  # Factor Analysis
  fa_result <- psych::fa(cor_matrix, nfactors = 5, rotate = "varimax")
  png(file.path(sigpred_dir, "factor_analysis.png"), 
      width = 1200, height = 800, res = 100)
  psych::fa.diagram(fa_result)
  dev.off()
  
  # Scree Plot
  png(file.path(sigpred_dir, "scree_plot.png"), 
      width = 800, height = 600, res = 100)
  psych::scree(cor_matrix, main = "Scree Plot")
  dev.off()
  
  # Save factor analysis results
  sink(file.path(sigpred_dir, "factor_analysis_results.txt"))
  print(fa_result)
  sink()
}

# Function to save model results
save_model_results <- function(results, data, subject_id, config) {
  output_dir <- file.path(config$output$base_dir, subject_id, 
                         config$output$elastic_net)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save model summary
  summary_file <- file.path(output_dir, sprintf("%s_model_summary.txt", subject_id))
  sink(summary_file)
  cat("Elastic Net Model Summary\n")
  cat("=======================\n\n")
  
  cat("Model Parameters:\n")
  cat(sprintf("Alpha: %.2f\n", config$elastic_net$alpha_grid[1]))
  cat(sprintf("Lambda: %.4f\n", results$cv_fit$lambda.min))
  cat("\n")
  
  cat("Training Set Metrics:\n")
  cat(sprintf("MSE: %.4f\n", results$train_metrics$mse))
  cat(sprintf("RMSE: %.4f\n", results$train_metrics$rmse))
  cat(sprintf("MAE: %.4f\n", results$train_metrics$mae))
  cat(sprintf("R-squared: %.4f\n", results$train_metrics$r2))
  cat(sprintf("Adjusted R-squared: %.4f\n", results$train_metrics$adj_r2))
  cat("\n")
  
  cat("Test Set Metrics:\n")
  cat(sprintf("MSE: %.4f\n", results$test_metrics$mse))
  cat(sprintf("RMSE: %.4f\n", results$test_metrics$rmse))
  cat(sprintf("MAE: %.4f\n", results$test_metrics$mae))
  cat(sprintf("R-squared: %.4f\n", results$test_metrics$r2))
  cat(sprintf("Adjusted R-squared: %.4f\n", results$test_metrics$adj_r2))
  cat("\n")
  
  cat("Cross-Validation Results:\n")
  print(results$cv_fit)
  cat("\n")
  
  cat("Feature Importance:\n")
  importance <- analyze_features(results$model, data)
  print(importance)
  sink()
  
  # Save significant predictors
  save_significant_predictors(importance, subject_id, config)
  
  # Save predictions
  predictions_file <- file.path(output_dir, sprintf("%s_predictions.csv", subject_id))
  readr::write_csv(
    data.frame(
      set = c(rep("train", length(results$train_predictions)),
              rep("test", length(results$test_predictions))),
      actual = c(data$AT8_intensity[1:length(results$train_predictions)],
                data$AT8_intensity[(length(results$train_predictions) + 1):nrow(data)]),
      predicted = c(results$train_predictions, results$test_predictions)
    ),
    predictions_file
  )
  
  # Save model
  model_file <- file.path(output_dir, sprintf("%s_model.rds", subject_id))
  saveRDS(results$model, model_file)
  
  logging::loginfo("Saved model results to %s", output_dir)
}

# Main Elastic Net function
run_elastic_net <- function(subject_id, config) {
  logging::loginfo("Starting Elastic Net analysis for subject %s", subject_id)
  
  # Check if model already exists
  model_file <- file.path(config$output$base_dir, subject_id, 
                         config$output$elastic_net,
                         sprintf("%s_model.rds", subject_id))
  
  if (file.exists(model_file)) {
    logging::loginfo("Model already exists for subject %s. Loading from file.", subject_id)
    results <- readRDS(model_file)
    
    # Load predictions if they exist
    predictions_file <- file.path(config$output$base_dir, subject_id,
                                config$output$elastic_net,
                                sprintf("%s_predictions.csv", subject_id))
    if (file.exists(predictions_file)) {
      predictions <- readr::read_csv(predictions_file, show_col_types = FALSE)
      results$train_predictions <- predictions$predicted[predictions$set == "train"]
      results$test_predictions <- predictions$predicted[predictions$set == "test"]
    }
    
    return(results)
  }
  
  # Load transformed data
  data <- load_transformed_data(subject_id, config)
  
  # Train model
  results <- train_elastic_net(data, config)
  
  # Save results
  save_model_results(results, data, subject_id, config)
  
  logging::loginfo("Completed Elastic Net analysis for subject %s", subject_id)
  return(results)
} 