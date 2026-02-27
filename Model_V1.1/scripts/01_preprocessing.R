#!/usr/bin/env Rscript

# NeuroPathPredict Pipeline V1.1
# Preprocessing script for data cleaning and validation

# Helper Functions
standardize <- function(x) {
  tm <- mean(x, na.rm = TRUE)
  ts <- sd(x, na.rm = TRUE)
  round((x - tm) / ts, 3)
}

cbrt_transform <- function(x) {
  sign(x) * abs(x)^(1/3)
}

clean_colnames <- function(df) {
  strings_to_remove <- c("_proc_0.5_MFG", "_tal_nlin_sym", 
                        "2011_17Networks_MNI152_FreeSurferConformed1mm_LooseMask",
                        "2011_17Networks_MNI152_FreeSurferConformed1mm_LiberalMask")
  
  # Clean column names
  clean_names <- colnames(df)
  for (string in strings_to_remove) {
    clean_names <- gsub(string, "", clean_names)
  }
  
  # Remove any repeated underscores
  clean_names <- gsub("_+", "_", clean_names)
  # Remove leading/trailing underscores
  clean_names <- gsub("^_|_$", "", clean_names)
  
  colnames(df) <- clean_names
  return(df)
}

check_missing_values <- function(df, df_name) {
  missing_count <- sum(is.na(df))
  if (missing_count > 0) {
    logging::logwarn("Found %d missing values in %s", missing_count, df_name)
    missing_cols <- colnames(df)[colSums(is.na(df)) > 0]
    logging::logwarn("Columns with missing values: %s", paste(missing_cols, collapse = ", "))
  } else {
    logging::loginfo("No missing values found in %s", df_name)
  }
}

# Additional transformation functions
sqrt_transform <- function(x) {
  sqrt(x)
}

inverse_transform <- function(x) {
  1/(x + 1e-10)  # Add small constant to avoid division by zero
}

squared_transform <- function(x) {
  x^2
}

# Function to calculate distribution metrics
calculate_distribution_metrics <- function(x) {
  # Remove any NA values
  x <- x[!is.na(x)]
  
  # Handle empty data
  if (length(x) == 0) {
    return(list(
      mean = NA,
      median = NA,
      sd = NA,
      skewness = NA,
      kurtosis = NA
    ))
  }
  
  # Calculate metrics and round to 3 decimal places
  metrics <- list(
    mean = round(mean(x), 3),
    median = round(median(x), 3),
    sd = round(sd(x), 3),
    skewness = round(moments::skewness(x), 3),
    kurtosis = round(moments::kurtosis(x), 3)
  )
  
  return(metrics)
}

# STEP 1: Merge and Subset Data
merge_and_subset_data <- function(config, subject_id) {
  # Convert subject_id to integer
  subject_id <- as.integer(subject_id)
  logging::loginfo("STEP 1: Starting data merge and subset for subject %d", subject_id)
  
  # Load input data
  logging::loginfo("Loading input data...")
  qnp_data <- data.table::as.data.table(data.table::fread(config$input$raw$qnp_data))
  cov_data <- data.table::as.data.table(data.table::fread(config$input$raw$covariates))
  
  # Ensure sub.id is integer
  qnp_data[, sub.id := as.integer(sub.id)]
  
  # Log initial data about subject
  subject_rows <- nrow(qnp_data[sub.id == subject_id])
  logging::loginfo("Subject %d data in QNP file: %d rows", subject_id, subject_rows)
  
  if (subject_rows == 0) {
    logging::logwarn("Available subjects in QNP data: %s", 
                    paste(unique(qnp_data$sub.id), collapse = ", "))
    stop(sprintf("Subject %d not found in QNP data", subject_id))
  }
  
  # Log data dimensions and column names
  logging::loginfo("QNP data: %d rows, %d columns", nrow(qnp_data), ncol(qnp_data))
  logging::loginfo("QNP columns: %s", paste(names(qnp_data), collapse = ", "))
  
  logging::loginfo("Covariate data: %d rows, %d columns", nrow(cov_data), ncol(cov_data))
  logging::loginfo("First few covariate columns: %s", paste(names(cov_data)[1:10], collapse = ", "))
  
  # Set keys for faster merging
  data.table::setkey(qnp_data, X, Y, Z)
  data.table::setkey(cov_data, X, Y, Z)
  
  # Validate input data
  check_missing_values(qnp_data, "QNP data")
  check_missing_values(cov_data, "Covariate data")
  
  # Merging QNP data with covariates
  logging::loginfo("Merging QNP data with covariates...")
  merged_data <- qnp_data[cov_data, nomatch=0]
  logging::loginfo("After merge: %d rows", nrow(merged_data))
  
  subject_rows <- nrow(merged_data[sub.id == subject_id])
  logging::loginfo("Subject %d rows after merge: %d", subject_id, subject_rows)
  
  if (subject_rows == 0) {
    logging::logwarn("Lost subject data during merge. Checking coordinate overlap...")
    qnp_coords <- qnp_data[sub.id == subject_id, .(X, Y, Z)]
    matching_coords <- cov_data[.(qnp_coords), nomatch=0]
    logging::logwarn("Found %d matching coordinates in covariates", nrow(matching_coords))
    stop("Lost subject data during merge")
  }
  
  # Rename columns if needed
  if ("AT8 intensity" %in% names(merged_data)) {
    data.table::setnames(merged_data, "AT8 intensity", "AT8_intensity")
  }
  
  # Subset for subject
  logging::loginfo("Subsetting data for subject %d", subject_id)
  subject_data <- merged_data[sub.id == subject_id]
  
  if (nrow(subject_data) == 0) {
    logging::logwarn("Available subjects in merged data: %s", 
                    paste(unique(merged_data$sub.id), collapse = ", "))
    stop(sprintf("No data found for subject %d", subject_id))
  }
  logging::loginfo("Final subject data: %d rows", nrow(subject_data))
  
  # Save intermediate output
  output_dir <- file.path(config$output$base_dir, subject_id, config$output$preprocessing)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_file <- file.path(output_dir, sprintf("%d_step1_merged.csv", subject_id))
  data.table::fwrite(subject_data, output_file)
  logging::loginfo("Saved Step 1 output to: %s", output_file)
  
  return(subject_data)
}

# STEP 2: Clean Zero Values
clean_zero_values <- function(data, subject_id, config) {
  # Convert subject_id to integer for consistent handling
  subject_id <- as.integer(subject_id)
  logging::loginfo("STEP 2: Cleaning zero values for subject %d", subject_id)
  
  # Remove specified columns
  cols_to_remove <- c("i", "j", "k", "i.i", "i.j", "i.k")
  existing_cols <- cols_to_remove[cols_to_remove %in% names(data)]
  if (length(existing_cols) > 0) {
    logging::loginfo("Removing columns: %s", paste(existing_cols, collapse = ", "))
    data[, (existing_cols) := NULL]
  }
  
  # Calculate zero proportions for each column
  zero_props <- colMeans(data == 0, na.rm = TRUE)
  high_zero_cols <- names(zero_props)[zero_props > 0.9]
  
  if (length(high_zero_cols) > 0) {
    logging::loginfo("Removing %d columns with >90%% zero values: %s", 
                    length(high_zero_cols), 
                    paste(high_zero_cols, collapse = ", "))
    data[, (high_zero_cols) := NULL]
  } else {
    logging::loginfo("No columns found with >90%% zero values")
  }
  
  # Clean column names
  data <- clean_colnames(data)
  
  # Save intermediate output
  output_dir <- file.path(config$output$base_dir, as.character(subject_id), config$output$preprocessing)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_file <- file.path(output_dir, sprintf("%d_step2_cleaned.csv", subject_id))
  data.table::fwrite(data, output_file)
  logging::loginfo("Saved Step 2 output to: %s", output_file)
  
  return(data)
}

# STEP 3: Transform and Standardize Data
transform_and_standardize <- function(data, subject_id, config) {
  # Convert subject_id to integer for consistent handling
  subject_id <- as.integer(subject_id)
  logging::loginfo("STEP 3: Transforming and standardizing data for subject %d", subject_id)
  
  # Create a copy of the data
  transformed_data <- data
  
  # Initialize metrics list
  metrics_list <- list()
  
  # Get numeric columns (excluding X, Y, Z, sub.id)
  exclude_cols <- c("X", "Y", "Z", "sub.id")
  numeric_cols <- names(data)[sapply(data, is.numeric)]
  numeric_cols <- setdiff(numeric_cols, exclude_cols)
  
  # Calculate pre-transformation metrics
  pre_metrics <- list()
  for (col in numeric_cols) {
    pre_metrics[[col]] <- calculate_distribution_metrics(data[[col]])
  }
  
  # Apply transformations based on variable type
  for (col in numeric_cols) {
    if (grepl("^BZ", col)) {
      # Apply cube root transformation to BZ variables
      transformed_data[[col]] <- round(cbrt_transform(data[[col]]), 3)
    } else if (col == "AT8_intensity") {
      # Apply log transformation to AT8_intensity
      transformed_data[[col]] <- round(log1p(data[[col]]), 3)
    }
    # EDT and Var variables remain unchanged but still rounded
    else {
      transformed_data[[col]] <- round(data[[col]], 3)
    }
  }
  
  # Calculate post-transformation metrics
  post_transform_metrics <- list()
  for (col in numeric_cols) {
    post_transform_metrics[[col]] <- calculate_distribution_metrics(transformed_data[[col]])
  }
  
  # Apply standardization to all numeric columns except X, Y, Z, sub.id
  for (col in numeric_cols) {
    transformed_data[[col]] <- round(standardize(transformed_data[[col]]), 3)
  }
  
  # Calculate post-standardization metrics
  post_standardize_metrics <- list()
  for (col in numeric_cols) {
    post_standardize_metrics[[col]] <- calculate_distribution_metrics(transformed_data[[col]])
  }
  
  # Create metrics data frame
  metrics_df <- data.frame(
    variable = character(),
    stage = character(),
    mean = numeric(),
    median = numeric(),
    sd = numeric(),
    skewness = numeric(),
    kurtosis = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Combine all metrics
  for (col in numeric_cols) {
    # Pre-transformation metrics
    metrics_df <- rbind(metrics_df, data.frame(
      variable = col,
      stage = "pre_transformation",
      mean = round(pre_metrics[[col]]$mean, 3),
      median = round(pre_metrics[[col]]$median, 3),
      sd = round(pre_metrics[[col]]$sd, 3),
      skewness = round(pre_metrics[[col]]$skewness, 3),
      kurtosis = round(pre_metrics[[col]]$kurtosis, 3),
      stringsAsFactors = FALSE
    ))
    
    # Post-transformation metrics
    metrics_df <- rbind(metrics_df, data.frame(
      variable = col,
      stage = "post_transformation",
      mean = round(post_transform_metrics[[col]]$mean, 3),
      median = round(post_transform_metrics[[col]]$median, 3),
      sd = round(post_transform_metrics[[col]]$sd, 3),
      skewness = round(post_transform_metrics[[col]]$skewness, 3),
      kurtosis = round(post_transform_metrics[[col]]$kurtosis, 3),
      stringsAsFactors = FALSE
    ))
    
    # Post-standardization metrics
    metrics_df <- rbind(metrics_df, data.frame(
      variable = col,
      stage = "post_standardization",
      mean = round(post_standardize_metrics[[col]]$mean, 3),
      median = round(post_standardize_metrics[[col]]$median, 3),
      sd = round(post_standardize_metrics[[col]]$sd, 3),
      skewness = round(post_standardize_metrics[[col]]$skewness, 3),
      kurtosis = round(post_standardize_metrics[[col]]$kurtosis, 3),
      stringsAsFactors = FALSE
    ))
  }
  
  # Save metrics
  output_dir <- file.path(config$output$base_dir, as.character(subject_id), config$output$preprocessing)
  metrics_file <- file.path(output_dir, sprintf("%d_covariate_metrics.csv", subject_id))
  
  # Round all numeric columns in metrics_df to 3 decimal places
  numeric_cols <- sapply(metrics_df, is.numeric)
  metrics_df[numeric_cols] <- round(metrics_df[numeric_cols], 3)
  
  write.csv(metrics_df, metrics_file, row.names = FALSE)
  logging::loginfo("Saved covariate metrics to: %s", metrics_file)
  
  # Save transformed and standardized data
  output_file <- file.path(output_dir, sprintf("%d_step3_transformed_standardized.csv", subject_id))
  data.table::fwrite(transformed_data, output_file)
  logging::loginfo("Saved transformed and standardized data to: %s", output_file)
  
  return(transformed_data)
}

# Main preprocessing function that runs all steps
preprocess_data <- function(subject_id, config) {
  # Convert subject_id to integer for consistent handling
  subject_id <- as.integer(subject_id)
  logging::loginfo("Starting preprocessing pipeline for subject %d", subject_id)
  
  # Define output paths
  output_dir <- file.path(config$output$base_dir, as.character(subject_id), config$output$preprocessing)
  step1_file <- file.path(output_dir, sprintf("%d_step1_merged.csv", subject_id))
  step2_file <- file.path(output_dir, sprintf("%d_step2_cleaned.csv", subject_id))
  step3_file <- file.path(output_dir, sprintf("%d_step3_transformed_standardized.csv", subject_id))
  
  # Check if step 3 output exists
  if (file.exists(step3_file)) {
    logging::loginfo("Transformed and standardized data already exists for subject %d. Skipping preprocessing.", subject_id)
    return(data.table::fread(step3_file))
  }
  
  # Run each step sequentially, checking for intermediate files
  if (file.exists(step1_file)) {
    logging::loginfo("Step 1 output exists. Loading from: %s", step1_file)
    step1_data <- data.table::fread(step1_file)
  } else {
    logging::loginfo("Running Step 1: Merge and subset data")
    step1_data <- merge_and_subset_data(config, subject_id)
  }
  
  if (file.exists(step2_file)) {
    logging::loginfo("Step 2 output exists. Loading from: %s", step2_file)
    step2_data <- data.table::fread(step2_file)
  } else {
    logging::loginfo("Running Step 2: Clean zero values")
    step2_data <- clean_zero_values(step1_data, subject_id, config)
  }
  
  if (file.exists(step3_file)) {
    logging::loginfo("Step 3 output exists. Loading from: %s", step3_file)
    final_data <- data.table::fread(step3_file)
  } else {
    logging::loginfo("Running Step 3: Transform and standardize data")
    final_data <- transform_and_standardize(step2_data, subject_id, config)
  }
  
  logging::loginfo("Completed preprocessing steps for subject %d", subject_id)
  return(final_data)
} 