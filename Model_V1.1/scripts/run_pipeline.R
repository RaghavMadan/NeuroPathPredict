#!/usr/bin/env Rscript

# NeuroPathPredict Pipeline V1.1
# Main pipeline script

# Get the directory where the script is located
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script_dir <- dirname(script.name)

# Set working directory to the project root
setwd(dirname(script_dir))

# Source setup script
source(file.path("scripts", "setup.R"))

# Source pipeline components
source(file.path("scripts", "01_preprocessing.R"))
source(file.path("scripts", "02_elastic_net.R"))
source(file.path("scripts", "03_factor_analysis.R"))
source(file.path("scripts", "04_universal_kriging.R"))

# Read configuration
config <- config::get(file = "config/input_variables.txt")
subjects <- readLines("config/subject_list.txt")
subjects <- subjects[!grepl("^#", subjects)]  # Remove comments
subjects <- subjects[nzchar(trimws(subjects))]  # Remove blank lines

# Initialize logging
log_file <- file.path("output", "pipeline.log")
# Set logging level to DEBUG to capture detailed variogram info
logging::basicConfig(level = "DEBUG")
logging::addHandler(logging::writeToFile, file = log_file)

# Function to create subject-specific UK config file
create_uk_config <- function(subject_id, config) {
  # Define path for subject-specific config
  uk_config_dir <- file.path(config$output$base_dir, subject_id, 
                            config$universal_kriging$output_dir)
  uk_config_file <- file.path(uk_config_dir, sprintf("%s_uk_config.txt", subject_id))
  
  # Skip if config already exists
  if (file.exists(uk_config_file)) {
    logging::loginfo("Subject-specific UK config already exists for %s", subject_id)
    return(uk_config_file)
  }
  
  # Create directory if it doesn't exist
  dir.create(uk_config_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Create config file with default values
  sink(uk_config_file)
  cat("# Universal Kriging Configuration for Subject", subject_id, "\n")
  cat("# Generated on:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  
  # Write each parameter in YAML format
  for (param in names(config$universal_kriging$default_config)) {
    value <- config$universal_kriging$default_config[[param]]
    if (is.character(value)) {
      cat(sprintf("%s: '%s'\n", param, value))
    } else if (param == "n_splits") {
      cat(sprintf("%s: [%s]\n", param, paste(value, collapse = ", ")))
    } else {
      cat(sprintf("%s: %s\n", param, value))
    }
  }
  sink()
  
  logging::loginfo("Created subject-specific UK config for %s", subject_id)
  return(uk_config_file)
}

# Process each subject
for (subject_id in subjects) {
  logging::loginfo("Processing subject %s", subject_id)
  
  # Step 1: Preprocessing (including transformations and standardization)
  logging::loginfo("Starting preprocessing for subject %s", subject_id)
  preprocessed_data <- preprocess_data(subject_id, config)
  
  # Step 2: Elastic Net Analysis
  logging::loginfo("Starting Elastic Net analysis for subject %s", subject_id)
  elastic_net_results <- run_elastic_net(subject_id, config)
  
  # Step 3: Factor Analysis
  logging::loginfo("Starting factor analysis for subject %s", subject_id)
  factor_analysis_results <- run_factor_analysis(subject_id, config)
  
  # Step 4: Universal Kriging (check predictors_reviewed.txt, then enabled parameter)
  predictors_reviewed_file <- file.path(config$output$base_dir, subject_id, config$output$factor_analysis, sprintf("%s_predictors_reviewed.txt", subject_id))
  if (!file.exists(predictors_reviewed_file)) {
    logging::logwarn("Universal Kriging skipped: %s_predictors_reviewed.txt does not exist.", subject_id)
    logging::logwarn("Please review the predictors and create the reviewed file before proceeding with Universal Kriging.")
    next
  }
  uk_config_file <- file.path(config$output$base_dir, subject_id, config$universal_kriging$output_dir, sprintf("%s_uk_config.txt", subject_id))
  
  # Create UK config file if it doesn't exist
  if (!file.exists(uk_config_file)) {
    logging::loginfo("Creating UK config file for subject %s", subject_id)
    uk_config_file <- create_uk_config(subject_id, config)
  }
  
  uk_enabled <- TRUE
  if (file.exists(uk_config_file)) {
    uk_lines <- readLines(uk_config_file)
    enabled_line <- grep("^enabled:", uk_lines, value = TRUE)
    if (length(enabled_line) > 0) {
      uk_enabled <- tolower(trimws(sub("enabled:", "", enabled_line))) == "true"
    }
  } else if (!is.null(config$universal_kriging$default_config$enabled)) {
    uk_enabled <- isTRUE(config$universal_kriging$default_config$enabled)
  }
  if (uk_enabled) {
    logging::loginfo("Universal kriging enabled for subject %s", subject_id)
    tryCatch({
      kriging_results <- run_universal_kriging(subject_id, config)
    }, error = function(e) {
      logging::logerror("Error processing subject %s: %s", subject_id, e$message)
    })
  } else {
    logging::loginfo("Universal kriging disabled for subject %s (enabled parameter is FALSE)", subject_id)
  }
}

logging::loginfo("Pipeline completed") 