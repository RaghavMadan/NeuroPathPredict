#!/usr/bin/env Rscript

# Factor Analysis script for NeuroPathPredict Pipeline V1.1

# Comment out command line argument processing
# args <- commandArgs(trailingOnly = TRUE)
# if (length(args) != 1) {
#   stop("Usage: Rscript 03_factor_analysis.R <subject_id>")
# }
# subject_id_arg <- args[1] # Renamed

# Comment out sourcing of common_functions and other setup
# source(file.path(dirname(sys.frame(1)$ofile), "common_functions.R"))
# env <- setup_environment()
# source(file.path("scripts", "setup.R"))
# source_functions()
# config_obj <- load_config() # Renamed
# log_file_path <- file.path("output", "pipeline.log") # Renamed
# init_logging(log_file_path)

# Comment out main execution block
# tryCatch({
#     logging::loginfo("Starting Factor Analysis for subject %s", subject_id_arg)
#     result <- run_factor_analysis(subject_id_arg, config_obj) # This function is defined below
#     
#     # Check if predictors need review (this logic might need to be in run_factor_analysis or run_pipeline)
#     output_dir_val <- file.path(config_obj$output$base_dir, subject_id_arg, config_obj$output$factor_analysis) # Renamed
#     predictors_reviewed_file <- file.path(output_dir_val, sprintf("%s_predictors_reviewed.txt", subject_id_arg))
#     predictors_to_review_file <- file.path(output_dir_val, sprintf("%s_predictors_to_review.txt", subject_id_arg))
#     
#     if (!file.exists(predictors_reviewed_file)) {
#         message("\n")
#         message(paste(rep("=", 80), collapse = ""))
#         message(sprintf("REVIEW REQUIRED: Review %s_predictors_to_review.txt and save it as %s_predictors_reviewed.txt before universal kriging can be implemented.", subject_id_arg, subject_id_arg))
#         message(paste(rep("=", 80), collapse = ""))
#         message("\n")
#         message(sprintf("Location of file to review: %s", predictors_to_review_file))
#         message(sprintf("Save reviewed file as: %s", predictors_reviewed_file))
#         message("\n")
#     }
#     
#     logging::loginfo("Completed Factor Analysis for subject %s", subject_id_arg)
# }, error = function(e) {
#     logging::logerror("Error in Factor Analysis for subject %s: %s", subject_id_arg, e$message)
#     quit(status = 1)
# })

# Function to load significant predictors
load_significant_predictors <- function(subject_id, config) {
  predictors_file <- file.path(config$output$base_dir, subject_id,
                              config$output$elastic_net,
                              sprintf("%s_significant_predictors.csv", subject_id))
  
  if (!file.exists(predictors_file)) {
    stop(sprintf("Significant predictors file not found: %s", predictors_file))
  }
  
  predictors <- read.csv(predictors_file)
  logging::loginfo("Loaded %d significant predictors", nrow(predictors))
  return(predictors)
}

# Function to load transformed data
load_transformed_data <- function(subject_id, config) {
  input_file <- file.path(config$output$base_dir, subject_id,
                         config$output$preprocessing,
                         sprintf("%s_step3_transformed_standardized.csv", subject_id))
  
  if (!file.exists(input_file)) {
    stop(sprintf("Transformed data file not found: %s", input_file))
  }
  
  data <- readr::read_csv(input_file, show_col_types = FALSE)
  logging::loginfo("Loaded transformed data: %d rows, %d columns", 
                  nrow(data), ncol(data))
  return(data)
}

# Function to perform factor analysis
perform_factor_analysis <- function(data, predictors, subject_id, config) {
  # Create output directory
  output_dir <- file.path(config$output$base_dir, subject_id, 
                         config$output$factor_analysis)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save predictors to review in factor analysis folder - simple list format
  predictors_file <- file.path(output_dir, sprintf("%s_predictors_to_review.txt", subject_id))
  writeLines(predictors$feature, predictors_file)
  
  # Prepare data for analysis
  pred_data <- data[, c("AT8_intensity", predictors$feature)]
  
  # Calculate correlation matrix
  cor_matrix <- cor(pred_data, use = "complete.obs")
  
  # Create visualizations
  # 1. Correlation plot with correlation values
  png(file.path(output_dir, sprintf("%s_correlation_plot.png", subject_id)), 
      width = 1200, height = 1200, res = 100)
  corrplot::corrplot(cor_matrix, method = "color", type = "upper", 
                    tl.col = "black", tl.srt = 45,
                    addCoef.col = "black", number.cex = 0.7,
                    tl.cex = 0.65)  # Reduced text size by ~35%
  dev.off()
  
  # 2. Dendrogram
  hc <- hclust(dist(cor_matrix))
  png(file.path(output_dir, sprintf("%s_dendrogram.png", subject_id)), 
      width = 1200, height = 800, res = 100)
  plot(hc, main = "Hierarchical Clustering of Predictors",
       xlab = "", sub = "", cex = 0.8)
  dev.off()
  
  # 3. Factor Analysis with improved settings
  # Try different factor analysis methods and rotations
  methods <- c("ml", "minres", "pa")
  rotations <- c("varimax", "promax", "oblimin")
  fa_result <- NULL
  
  # Get eigenvalues to determine number of factors
  ev <- eigen(cor_matrix)$values
  n_factors <- min(sum(ev > 1), min(6, ncol(pred_data) - 1))  # Use Kaiser criterion but cap at 6 factors
  
  if (n_factors < 2) n_factors <- 2  # Ensure at least 2 factors
  
  for (method in methods) {
    for (rotation in rotations) {
      tryCatch({
        fa_result <- psych::fa(cor_matrix, 
                             nfactors = n_factors, 
                             rotate = rotation, 
                             fm = method,
                             scores = "regression",  # Use regression scores
                             n.obs = nrow(pred_data))
        # If we get here, the factor analysis succeeded
        break
      }, error = function(e) {
        # Continue to next iteration if this combination fails
        NULL
      })
    }
    if (!is.null(fa_result)) break
  }
  
  if (is.null(fa_result)) {
    # If all combinations failed, try PCA as a fallback
    tryCatch({
      pca_result <- prcomp(pred_data, scale. = TRUE)
      n_comp <- min(which(cumsum(pca_result$sdev^2/sum(pca_result$sdev^2)) > 0.8))
      fa_result <- list(
        loadings = pca_result$rotation[, 1:min(n_comp, ncol(pca_result$rotation))],
        scores = pca_result$x[, 1:min(n_comp, ncol(pca_result$x))],
        method = "PCA fallback"
      )
    }, error = function(e) {
      stop("Both factor analysis and PCA fallback failed")
    })
  }
  
  # Create factor analysis diagram
  png(file.path(output_dir, sprintf("%s_factor_analysis.png", subject_id)), 
      width = 1200, height = 800, res = 100)
  psych::fa.diagram(fa_result)
  dev.off()
  
  # 4. Scree Plot
  png(file.path(output_dir, sprintf("%s_scree_plot.png", subject_id)), 
      width = 800, height = 600, res = 100)
  psych::scree(cor_matrix, main = "Scree Plot")
  dev.off()
  
  # Save factor analysis results
  sink(file.path(output_dir, sprintf("%s_factor_analysis_results.txt", subject_id)))
  print(fa_result)
  sink()
  
  # Save factor loadings
  loadings_file <- file.path(output_dir, sprintf("%s_factor_loadings.csv", subject_id))
  write.csv(fa_result$loadings, loadings_file)
  
  # Save factor scores if available
  if (!is.null(fa_result$scores)) {
    scores_file <- file.path(output_dir, sprintf("%s_factor_scores.csv", subject_id))
    write.csv(fa_result$scores, scores_file)
  }
  
  logging::loginfo("Saved factor analysis results to %s", output_dir)
  return(fa_result)
}

# Function to find highly correlated pairs
find_high_correlations <- function(cor_matrix, threshold = 0.8) {
  # Get upper triangle of correlation matrix
  cor_matrix[upper.tri(cor_matrix, diag = TRUE)] <- NA
  
  # Find pairs with absolute correlation > threshold
  high_cor_pairs <- which(abs(cor_matrix) > threshold, arr.ind = TRUE)
  
  # Create data frame with variable pairs and their correlations
  if (nrow(high_cor_pairs) > 0) {
    pairs_df <- data.frame(
      Var1 = rownames(cor_matrix)[high_cor_pairs[, 1]],
      Var2 = colnames(cor_matrix)[high_cor_pairs[, 2]],
      Correlation = cor_matrix[high_cor_pairs]
    )
    # Sort by absolute correlation
    pairs_df <- pairs_df[order(abs(pairs_df$Correlation), decreasing = TRUE), ]
    return(pairs_df)
  } else {
    return(NULL)
  }
}

# Function to create SigPred visualizations
create_sigpred_visualizations <- function(data, predictors, subject_id, config) {
  # Create SigPred directory
  sigpred_dir <- file.path(config$output$base_dir, subject_id, 
                          config$output$factor_analysis, "SigPred")
  dir.create(sigpred_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Prepare data for visualizations
  pred_data <- data[, c("AT8_intensity", predictors)]
  
  # Calculate correlation matrix
  cor_matrix <- cor(pred_data, use = "complete.obs")
  
  # Correlation plot with correlation values
  png(file.path(sigpred_dir, sprintf("%s_correlation_plot.png", subject_id)), 
      width = 1200, height = 1200, res = 100)
  corrplot::corrplot(cor_matrix, method = "color", type = "upper", 
                    tl.col = "black", tl.srt = 45,
                    addCoef.col = "black", number.cex = 0.7,
                    tl.cex = 0.65)  # Reduced text size by ~35%
  dev.off()
  
  # Dendrogram
  hc <- hclust(dist(cor_matrix))
  png(file.path(sigpred_dir, sprintf("%s_dendrogram.png", subject_id)), 
      width = 1200, height = 800, res = 100)
  plot(hc, main = "Hierarchical Clustering of Predictors",
       xlab = "", sub = "", cex = 0.8)
  dev.off()
  
  # Find and save highly correlated pairs
  high_cor_pairs <- find_high_correlations(cor_matrix, threshold = 0.8)
  if (!is.null(high_cor_pairs)) {
    # Save to CSV
    write.csv(high_cor_pairs, 
              file.path(sigpred_dir, sprintf("%s_high_correlations.csv", subject_id)),
              row.names = FALSE)
    
    # Save to text file with formatted output
    sink(file.path(sigpred_dir, sprintf("%s_high_correlations.txt", subject_id)))
    cat("Highly Correlated Variable Pairs (|correlation| > 0.8)\n")
    cat("===================================================\n\n")
    for (i in 1:nrow(high_cor_pairs)) {
      cat(sprintf("%d. %s - %s (correlation = %.3f)\n",
                 i, high_cor_pairs$Var1[i], high_cor_pairs$Var2[i],
                 high_cor_pairs$Correlation[i]))
    }
    sink()
  } else {
    # Create empty files to indicate no high correlations found
    write.csv(data.frame(Var1 = character(), Var2 = character(), Correlation = numeric()),
              file.path(sigpred_dir, sprintf("%s_high_correlations.csv", subject_id)),
              row.names = FALSE)
    
    sink(file.path(sigpred_dir, sprintf("%s_high_correlations.txt", subject_id)))
    cat("No variable pairs found with |correlation| > 0.8\n")
    sink()
  }
}

# Main factor analysis function
run_factor_analysis <- function(subject_id, config) {
  logging::loginfo("Starting factor analysis for subject %s", subject_id)

  output_dir <- file.path(config$output$base_dir, subject_id, config$output$factor_analysis)
  predictors_reviewed_file <- file.path(output_dir, sprintf("%s_predictors_reviewed.txt", subject_id))

  # If predictors_reviewed.txt exists, only run SigPred analysis
  if (file.exists(predictors_reviewed_file)) {
    predictors <- data.frame(feature = readLines(predictors_reviewed_file))
    data <- load_transformed_data(subject_id, config)
    create_sigpred_visualizations(data, predictors$feature, subject_id, config)
    logging::loginfo("Created SigPred visualizations for subject %s", subject_id)
    return(invisible(NULL))
  }

  # Otherwise, run the main (pre-review) factor analysis as before
  predictors <- load_significant_predictors(subject_id, config)
  data <- load_transformed_data(subject_id, config)
  fa_result <- perform_factor_analysis(data, predictors, subject_id, config)

  # Add subject ID to factor loadings file
  loadings_file <- file.path(output_dir, sprintf("%s_factor_loadings.csv", subject_id))
  if (file.exists(loadings_file)) {
    loadings <- read.csv(loadings_file)
    loadings$subject_id <- subject_id
    write.csv(loadings, loadings_file, row.names = FALSE)
  }

  logging::loginfo("Completed factor analysis for subject %s", subject_id)
  return(fa_result)
} 