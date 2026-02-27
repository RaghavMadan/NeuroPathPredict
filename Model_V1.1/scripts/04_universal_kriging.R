#!/usr/bin/env Rscript

# Universal Kriging script for NeuroPathPredict Pipeline V1.1

# Add ggplot2 for plotting
library(ggplot2)

# Function to prepare spatial data
prepare_spatial_data <- function(data, coords) {
  # No spatial object creation, just ensure numeric
  data <- as.data.frame(data)
  data[] <- lapply(data, function(x) as.numeric(as.character(x)))
  return(data)
}

# Function to load subject-specific UK config
load_uk_config <- function(subject_id, config) {
  uk_config_file <- file.path(config$output$base_dir, subject_id,
                             config$universal_kriging$output_dir,
                             sprintf("%s_uk_config.txt", subject_id))
  
  if (!file.exists(uk_config_file)) {
    stop(sprintf("Subject-specific UK config not found: %s", uk_config_file))
  }
  
  uk_config <- yaml::read_yaml(uk_config_file)
  return(uk_config)
}

# Function to plot initial variogram
plot_initial_variogram <- function(v, subject_id, uk_config, config) {
  # Create output directory for plots
  plot_dir <- file.path(config$output$base_dir, subject_id, config$universal_kriging$output_dir, "plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save empirical variogram data for review
  write.csv(v, file.path(plot_dir, "empirical_variogram.csv"), row.names = FALSE)
  
  # Create and save variogram plot
  v_df <- as.data.frame(v)
  initial_plot_path <- file.path(plot_dir, "initial_variogram.png")
  p_init <- ggplot2::ggplot(v_df, ggplot2::aes(x = dist, y = gamma)) +
    ggplot2::geom_point(color = 'blue', size = 2) +
    ggplot2::labs(title = paste("Initial Variogram - Subject", subject_id), x = "Distance", y = "Semivariance") +
    ggplot2::theme_minimal()
  ggplot2::ggsave(initial_plot_path, plot = p_init, width = 8, height = 6, dpi = 300)
  if (!file.exists(initial_plot_path)) {
    logging::logwarn(sprintf("initial_variogram.png was not created at %s", initial_plot_path))
  }
  logging::loginfo(paste("Initial variogram plot saved to", plot_dir))
}

# Function to fit variogram
fit_variogram <- function(df, uk_config, subject_id, config, variogram_formula) {
  logging::loginfo("Computing experimental variogram with formula: %s", deparse(variogram_formula))
  v <- gstat::variogram(variogram_formula, data = df, locations = ~X+Y+Z, cutoff = uk_config$cutoff, width = uk_config$width)
  #logging::loginfo("Empirical variogram (v) structure:")
  #logging::logdebug(capture.output(str(v)))
  #logging::loginfo("Empirical variogram (v) head:")
  #logging::logdebug(capture.output(head(v)))
  logging::loginfo("Fitting variogram model: %s", uk_config$variogram_model)
  v_fit <- gstat::fit.variogram(
    v,
    model = gstat::vgm(
      psill = uk_config$psill,
      model = uk_config$variogram_model,
      range = uk_config$range,
      nugget = uk_config$nugget
    )
  )
  logging::loginfo("Fitted variogram model (v_fit):")
  logging::logdebug(capture.output(print(v_fit)))
  # Save only the fitted variogram plot with annotation and white background
  plot_dir <- file.path(config$output$base_dir, subject_id, config$universal_kriging$output_dir, "plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  fit_plot_path <- file.path(plot_dir, "fitted_variogram.png")
  v_df <- as.data.frame(v)
  fit_line <- gstat::variogramLine(v_fit, maxdist = max(v_df$dist), n = 200)
  fit_df <- data.frame(dist = fit_line$dist, gamma = fit_line$gamma)
  # Extract fitted parameters
  nugget <- v_fit$psill[1]
  psill <- v_fit$psill[2]
  range <- v_fit$range[2]
  model <- as.character(v_fit$model[2])
  param_text <- sprintf("Nugget: %.3f\nPartial Sill: %.3f\nRange: %.3f\nModel: %s", nugget, psill, range, model)
  p_fit <- ggplot2::ggplot() +
    ggplot2::geom_point(data = v_df, ggplot2::aes(x = dist, y = gamma), color = 'blue', size = 2) +
    ggplot2::geom_line(data = fit_df, ggplot2::aes(x = dist, y = gamma), color = 'red', linewidth = 1) +
    ggplot2::labs(title = "Fitted Variogram Model", x = "Distance", y = "Semivariance") +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    ) +
    ggplot2::annotate(
      "text",
      x = max(v_df$dist) * 0.98,
      y = min(v_df$gamma) + (max(v_df$gamma) - min(v_df$gamma)) * 0.1,
      label = param_text,
      hjust = 1, vjust = 0, size = 5, family = "mono", fontface = "bold"
    )
  ggplot2::ggsave(fit_plot_path, plot = p_fit, width = 8, height = 6, dpi = 300)
  if (!file.exists(fit_plot_path)) {
    logging::logwarn(sprintf("fitted_variogram.png was not created at %s", fit_plot_path))
  }
  logging::loginfo("Saved fitted variogram to: %s", fit_plot_path)
  # --- Update config file with fitted parameters ---
  uk_config_file <- file.path(config$output$base_dir, subject_id, config$universal_kriging$output_dir, sprintf("%s_uk_config.txt", subject_id))
  if (file.exists(uk_config_file)) {
    uk_lines <- readLines(uk_config_file)
    uk_lines <- gsub("^psill:.*", sprintf("psill: %.6f", psill), uk_lines)
    uk_lines <- gsub("^nugget:.*", sprintf("nugget: %.6f", nugget), uk_lines)
    uk_lines <- gsub("^range:.*", sprintf("range: %.6f", range), uk_lines)
    writeLines(uk_lines, uk_config_file)
    logging::loginfo("Updated %s with fitted variogram parameters.", uk_config_file)
  } else {
    logging::logwarn("Could not find uk_config.txt to update fitted parameters.")
  }
  return(list(empirical = v, fitted = v_fit))
}

# Updated function to create/load spatial block assignments
create_spatial_blocks <- function(sp_data, current_n_split_value, subject_id, config) {
  logging::loginfo("Managing spatial blocks for %d splits for subject %s", current_n_split_value, subject_id)

  # Define paths
  uk_data_path <- file.path(config$output$base_dir, subject_id, 
                            config$universal_kriging$output_dir, 
                            sprintf("%s_uk_data.csv", subject_id))
  python_script_path <- file.path("scripts", "05_block_visualization.py")  # Fixed path
  output_vis_dir <- file.path(config$output$base_dir, subject_id, 
                              config$universal_kriging$output_dir, 
                              "block_visualizations")
  
  # Load subject-specific UK config to get all n_splits values intended for the pipeline run
  uk_config <- load_uk_config(subject_id, config)
  all_pipeline_n_splits <- uk_config$n_splits # This should be a vector like c(2, 4, 8, 16)
  if (is.null(all_pipeline_n_splits) || length(all_pipeline_n_splits) == 0) {
    stop(sprintf("n_splits not defined or empty in uk_config for subject %s", subject_id))
  }

  # Ensure output directory for visualizations exists for Python script
  dir.create(output_vis_dir, recursive = TRUE, showWarnings = FALSE)

  # Check if the primary _uk_data.csv file exists
  if (!file.exists(uk_data_path)) {
    stop(sprintf("_uk_data.csv file not found at %s. It should be created by prepare_uk_data.", uk_data_path))
  }
  
  # Load the _uk_data.csv file
  logging::loginfo("Loading data from %s to check/add block assignments.", uk_data_path)
  subject_uk_data <- readr::read_csv(uk_data_path, show_col_types = FALSE)
  
  # Determine which n_splits columns are missing
  missing_n_splits_cols <- c()
  for (n_val in all_pipeline_n_splits) {
    col_name <- sprintf("%d_block_splits", n_val)
    if (!col_name %in% names(subject_uk_data)) {
      missing_n_splits_cols <- c(missing_n_splits_cols, n_val)
    }
  }
  
  # If any block assignment columns are missing, run the Python script
  if (length(missing_n_splits_cols) > 0) {
    logging::loginfo("Missing block assignment columns for n_splits: %s. Running Python script to generate them.", paste(missing_n_splits_cols, collapse=", "))
    n_splits_str_for_python <- paste(all_pipeline_n_splits, collapse = ",")
    cmd <- sprintf("python3 \"%s\" \"%s\" \"%s\" \"%s\" \"%s\"",
                   python_script_path,
                   uk_data_path,
                   n_splits_str_for_python,
                   output_vis_dir,
                   subject_id)
    logging::loginfo("Executing Python script: %s", cmd)
    status <- system(cmd)
    if (status != 0) {
      stop(sprintf("Python script for block creation/assignment failed with status %d.", status))
    }
    logging::loginfo("Reloading data from %s after Python script execution.", uk_data_path)
    subject_uk_data <- readr::read_csv(uk_data_path, show_col_types = FALSE)
  } else {
    logging::loginfo("All required block assignment columns already exist in %s. Python script not called.", uk_data_path)
  }
  
  # Extract the specific block assignment vector for the current_n_split_value
  target_col_name <- sprintf("%d_block_splits", current_n_split_value)
  if (!target_col_name %in% names(subject_uk_data)) {
    stop(sprintf("Target column %s not found in %s even after Python script execution. Check Python script logic and n_splits configuration.", 
                 target_col_name, uk_data_path))
  }
  
  block_assignments_vector <- subject_uk_data[[target_col_name]]
  
  if (is.null(block_assignments_vector) || length(block_assignments_vector) != nrow(subject_uk_data)){
      stop(sprintf("Block assignments for %d splits are null or have incorrect length.", current_n_split_value))
  }

  logging::loginfo("Successfully retrieved block assignments for %d splits.", current_n_split_value)
  return(block_assignments_vector)
}

# Function to assign a point to a block
assign_block <- function(x, y, z, blocks) {
  # Input validation
  if (!is.numeric(x) || !is.numeric(y) || !is.numeric(z)) {
    logging::logerror("Non-numeric coordinates provided: x=%s, y=%s, z=%s", x, y, z)
    return(NA)
  }
  
  if (length(blocks) == 0) {
    logging::logerror("Empty blocks list provided")
    return(NA)
  }
  
  # Track which blocks the point could potentially belong to
  potential_blocks <- c()
  min_dist <- Inf
  closest_block <- NA
  
  # Calculate epsilon based on coordinate scale
  coord_range <- max(abs(c(x, y, z)))
  eps <- coord_range * 1e-6  # Dynamic epsilon based on coordinate scale
  
  for (i in seq_along(blocks)) {
    b <- blocks[[i]]
    
    # Validate block structure
    if (!all(c("x", "y", "z") %in% names(b))) {
      logging::logerror("Invalid block structure for block %d", i)
      next
    }
    
    # Check if point is within block boundaries (with epsilon tolerance)
    if (x >= (b$x[1] - eps) && x <= (b$x[2] + eps) &&
        y >= (b$y[1] - eps) && y <= (b$y[2] + eps) &&
        z >= (b$z[1] - eps) && z <= (b$z[2] + eps)) {
      potential_blocks <- c(potential_blocks, i)
    }
    
    # Calculate distance to block center for tie-breaking
    center_x <- mean(b$x)
    center_y <- mean(b$y)
    center_z <- mean(b$z)
    dist <- sqrt((x - center_x)^2 + (y - center_y)^2 + (z - center_z)^2)
    
    if (dist < min_dist) {
      min_dist <- dist
      closest_block <- i
    }
  }
  
  if (length(potential_blocks) == 0) {
    logging::logwarn("Point (%.6f, %.6f, %.6f) not assigned to any block (closest block: %d, distance: %.6f)", 
                    x, y, z, closest_block, min_dist)
    # If point is very close to closest block, assign it anyway
    if (min_dist < eps * 10) {
      logging::loginfo("Assigning point to closest block %d due to small distance", closest_block)
      return(closest_block)
    }
    return(NA)
  } else if (length(potential_blocks) > 1) {
    # If multiple blocks are possible, choose the closest one
    block_distances <- sapply(potential_blocks, function(i) {
      b <- blocks[[i]]
      center_x <- mean(b$x)
      center_y <- mean(b$y)
      center_z <- mean(b$z)
      sqrt((x - center_x)^2 + (y - center_y)^2 + (z - center_z)^2)
    })
    chosen_block <- potential_blocks[which.min(block_distances)]
    logging::loginfo("Point (%.6f, %.6f, %.6f) assigned to closest block %d among %d candidates", 
                    x, y, z, chosen_block, length(potential_blocks))
    return(chosen_block)
  } else {
    logging::loginfo("Point (%.6f, %.6f, %.6f) assigned to block %d", x, y, z, potential_blocks[1])
    return(potential_blocks[1])
  }
}

# Function to plot observed vs predicted and predicted vs residuals
plot_obs_pred_residuals <- function(df, subject_id, n_split, config) {
  output_dir <- file.path(config$output$base_dir, subject_id, config$universal_kriging$output_dir, "plots")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  df <- as.data.frame(df)
  df <- df[!is.na(df$predicted) & !is.na(df$observed), ]
  if (nrow(df) == 0) return()
  
  # Filter out predictions outside the -10 to +10 range
  valid_pred_mask <- df$predicted >= -10 & df$predicted <= 10
  df_filtered <- df[valid_pred_mask, ]
  
  if (nrow(df_filtered) == 0) {
    logging::logwarn("No predictions within -10 to +10 range for plotting")
    return()
  }
  
  df_filtered$residual <- df_filtered$observed - df_filtered$predicted
  
  # Log exclusion statistics
  n_excluded <- sum(!valid_pred_mask)
  n_total <- nrow(df)
  if (n_total > 0) {
    logging::loginfo("Plot filtering: %d/%d predictions excluded (%.2f%%) for being outside -10 to +10 range", 
                    n_excluded, n_total, 100 * n_excluded / n_total)
  }
  
  # Observed vs Predicted
  p1 <- ggplot2::ggplot(df_filtered, ggplot2::aes(x = observed, y = predicted)) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = 'red', linetype = 'dashed') +
    ggplot2::labs(title = sprintf("Observed vs Predicted (n_splits=%d)", n_split), x = "Observed", y = "Predicted") +
    ggplot2::theme_classic() +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", color = NA),
                   plot.background = ggplot2::element_rect(fill = "white", color = NA))
  ggplot2::ggsave(filename = file.path(output_dir, sprintf("obs_vs_pred_%d_splits.png", n_split)), plot = p1, width = 6, height = 5, dpi = 150)
  
  # Predicted vs Residuals
  p2 <- ggplot2::ggplot(df_filtered, ggplot2::aes(x = predicted, y = residual)) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::geom_hline(yintercept = 0, color = 'red', linetype = 'dashed') +
    ggplot2::labs(title = sprintf("Predicted vs Residuals (n_splits=%d)", n_split), x = "Predicted", y = "Residual") +
    ggplot2::theme_classic() +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", color = NA),
                   plot.background = ggplot2::element_rect(fill = "white", color = NA))
  ggplot2::ggsave(filename = file.path(output_dir, sprintf("pred_vs_residual_%d_splits.png", n_split)), plot = p2, width = 6, height = 5, dpi = 150)
}

# Function to perform cross-validation using blocks with optimizations
perform_cross_validation <- function(df, block_assignments, formula, uk_config, current_n_split_value,
 subject_id = NULL, config = NULL, n_split = NULL) {
  if (!requireNamespace("future", quietly = TRUE)) install.packages("future")
  if (!requireNamespace("future.apply", quietly = TRUE)) install.packages("future.apply")
  library(future)
  library(future.apply)
  block_assignments <- as.numeric(block_assignments)
  n_blocks <- max(block_assignments)
  n_cores <- min(parallel::detectCores() - 1, n_blocks)
  plan(multicore, workers = n_cores)
  results <- future_lapply(1:n_blocks, function(block) {
    train_idx <- block_assignments != block
    test_idx <- block_assignments == block
    train_data <- df[train_idx, ]
    test_data <- df[test_idx, ]
    
    # Store original row indices for coordinate mapping
    original_indices <- which(test_idx)
    
    # Fit variogram (external drift, block-based)
    v <- gstat::variogram(formula, data = train_data, locations = ~X+Y+Z)
    v_fit <- gstat::fit.variogram(v, gstat::vgm(uk_config$psill, uk_config$variogram_model, uk_config$range, uk_config$nugget))
    # Log fitted variogram parameters for this block
    logging::loginfo(
      "Block %d variogram fit: psill=%.4f, range=%.4f, nugget=%.4f",
      block,
      v_fit$psill[2],
      v_fit$range[2],
      v_fit$psill[1]
    )
    # Create gstat object (external drift)
    g <- gstat::gstat(
      formula = formula, # external drift, no X/Y/Z in trend
      data = train_data,
      model = v_fit,
      locations = ~X+Y+Z,
      maxdist = as.numeric(uk_config$maxdist),
      nmin = as.numeric(uk_config$nmin),
      nmax = as.numeric(uk_config$nmax),
      #force = TRUE
    )

    # Predict in chunks
    chunk_size <- 500
    n_chunks <- ceiling(nrow(test_data) / chunk_size)
    predictions <- vector("list", n_chunks)
    for(i in 1:n_chunks) {
      start_idx <- (i-1) * chunk_size + 1
      end_idx <- min(i * chunk_size, nrow(test_data))
      chunk_data <- test_data[start_idx:end_idx, ]
      chunk_original_indices <- original_indices[start_idx:end_idx]
      tryCatch({
        pred <- predict(g, newdata = chunk_data)
        predictions[[i]] <- data.frame(
          observed = as.numeric(chunk_data$response),
          predicted = as.numeric(pred$var1.pred),
          variance = as.numeric(pred$var1.var),
          X = chunk_data$X,
          Y = chunk_data$Y,
          Z = chunk_data$Z,
          original_row_index = chunk_original_indices
        )
      }, error = function(e) {
        logging::logwarn("Kriging failed for chunk %d in block %d: %s", i, block, e$message)
        return(data.frame(
          observed = as.numeric(chunk_data$response),
          predicted = NA_real_,
          variance = NA_real_,
          X = chunk_data$X,
          Y = chunk_data$Y,
          Z = chunk_data$Z,
          original_row_index = chunk_original_indices
        ))
      })
    }
    block_pred <- do.call(rbind, predictions)
    block_pred$block <- block
    # Exclude predictions outside the range -10 to +10
    valid_pred <- !is.na(block_pred$predicted)
    outlier_mask <- block_pred$predicted < -10 | block_pred$predicted > 10
    exclude_mask <- valid_pred & outlier_mask
    block_pred$excluded <- FALSE
    block_pred$excluded[exclude_mask] <- TRUE
    # For metrics, use only non-excluded predictions
    metric_mask <- valid_pred & !block_pred$excluded
    n_excluded <- sum(block_pred$excluded, na.rm = TRUE)
    n_total <- nrow(block_pred)
    pct_excluded <- ifelse(n_total > 0, 100 * n_excluded / n_total, 0)
    # Skip metrics if all predictions failed or all excluded
    if (sum(metric_mask) == 0) {
      n_split_val <- if(is.null(n_split)) NA_integer_ else as.integer(n_split)
      block_val <- as.integer(block)
      block_metrics <- data.frame(
        n_split = n_split_val,
        block = block_val,
        mse = NA_real_,
        rmse = NA_real_,
        mae = NA_real_,
        r2 = NA_real_,
        obs_var = NA_real_,
        n_total = n_total,
        n_valid = 0L,
        n_excluded = n_excluded,
        pct_excluded = pct_excluded
      )
      logging::loginfo("Block %d: No valid predictions after exclusion", block)
      empty_pred <- block_pred[FALSE, ]
      return(list(pred = empty_pred, metrics = block_metrics))
    }
    # Calculate metrics for blocks with valid predictions
    mse <- mean((block_pred$observed[metric_mask] - block_pred$predicted[metric_mask])^2)
    rmse <- sqrt(mse)
    mae <- mean(abs(block_pred$observed[metric_mask] - block_pred$predicted[metric_mask]))
    obs_var <- var(block_pred$observed[metric_mask])
    ss_res <- sum((block_pred$observed[metric_mask] - block_pred$predicted[metric_mask])^2)
    ss_tot <- sum((block_pred$observed[metric_mask] - mean(block_pred$observed[metric_mask]))^2)
    r2 <- ifelse(ss_tot == 0, NA, 1 - ss_res / ss_tot)
    n_valid <- sum(metric_mask)
    block_metrics <- data.frame(
      n_split = n_split,
      block = block,
      mse = mse,
      rmse = rmse,
      mae = mae,
      r2 = r2,
      obs_var = obs_var,
      n_total = n_total,
      n_valid = n_valid,
      n_excluded = n_excluded,
      pct_excluded = pct_excluded
    )
    logging::loginfo("Block %d: RMSE=%.4f, R2=%.4f, Excluded=%d (%.2f%%)", block, rmse, r2, n_excluded, pct_excluded)
    return(list(pred = block_pred, metrics = block_metrics))
  }, future.seed = TRUE)
  valid_results <- !sapply(results, is.null)
  if (sum(valid_results) == 0) {
    logging::logwarn("No valid predictions in any block")
    return(data.frame(
      observed = numeric(0),
      predicted = numeric(0),
      variance = numeric(0),
      block = integer(0),
      excluded = logical(0)
    ))
  }
  valid_preds <- lapply(results[valid_results], function(x) x$pred)
  valid_metrics <- lapply(results[valid_results], function(x) x$metrics)
  all_pred <- do.call(rbind, valid_preds)
  all_block_metrics <- do.call(rbind, valid_metrics)
  if (!is.null(subject_id) && !is.null(config) && !is.null(n_split)) {
    combined_metrics_file <- file.path(config$output$base_dir, subject_id, config$universal_kriging$output_dir, sprintf("%s_combined_block_metrics.csv", subject_id))
    if (file.exists(combined_metrics_file)) {
      write.table(all_block_metrics, combined_metrics_file, append = TRUE, sep = ",", row.names = FALSE, col.names = !file.exists(combined_metrics_file))
    } else {
      write.csv(all_block_metrics, combined_metrics_file, row.names = FALSE)
    }
    if (nrow(all_block_metrics) > 0) {
      mean_rmse <- mean(all_block_metrics$rmse, na.rm = TRUE)
      mean_r2 <- mean(all_block_metrics$r2, na.rm = TRUE)
      logging::loginfo("CV split %d summary: mean RMSE=%.4f, mean R2=%.4f", n_split, mean_rmse, mean_r2)
    }
  }
  # Return combined predictions and block metrics
  return(list(pred = all_pred, metrics = all_block_metrics))
}

# Function to create prediction grid
create_prediction_grid <- function(df, uk_config) {
  bbox <- apply(df[, c("X", "Y", "Z")], 2, range)
  grid_size <- as.numeric(uk_config$grid_size)
  grid <- expand.grid(
    X = seq(bbox[1, 1], bbox[2, 1], by = grid_size),
    Y = seq(bbox[1, 2], bbox[2, 2], by = grid_size),
    Z = seq(bbox[1, 3], bbox[2, 3], by = grid_size)
  )
  coordinates(grid) <- ~X+Y+Z
  return(grid)
}

# Function to perform kriging
perform_kriging <- function(df, grid_data, variogram, uk_config, formula_str) {
  # grid_data is a SpatialPoints object, df is a data.frame
  g <- gstat::gstat(
    formula = as.formula(formula_str),
    data = df,
    model = variogram,
    locations = ~X+Y+Z,
    nmax = as.numeric(uk_config$nmax),
    nmin = as.numeric(uk_config$nmin),
    maxdist = as.numeric(uk_config$maxdist),
    #force = TRUE
  )
  krige_result <- predict(g, newdata = grid_data)
  return(krige_result)
}

# Add new function for creating diagnostic plots
create_diagnostic_plots <- function(predictions, subject_id, config, n_splits, plot_dir) {
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Remove NA values and handle infinite values
  valid_data <- predictions %>%
    dplyr::filter(!is.na(observed) & !is.na(predicted) & 
                  !is.infinite(observed) & !is.infinite(predicted))
  
  if(nrow(valid_data) == 0) {
    logging::logwarn("No valid data points for diagnostic plots")
    return()
  }
  
  # Filter out predictions outside the -10 to +10 range
  valid_pred_mask <- valid_data$predicted >= -10 & valid_data$predicted <= 10
  valid_data_filtered <- valid_data[valid_pred_mask, ]
  
  if(nrow(valid_data_filtered) == 0) {
    logging::logwarn("No predictions within -10 to +10 range for diagnostic plots")
    return()
  }
  
  # Log exclusion statistics
  n_excluded <- sum(!valid_pred_mask)
  n_total <- nrow(valid_data)
  if (n_total > 0) {
    logging::loginfo("Diagnostic plot filtering: %d/%d predictions excluded (%.2f%%) for being outside -10 to +10 range", 
                    n_excluded, n_total, 100 * n_excluded / n_total)
  }
  
  # Ensure residual column exists
  if (!"residual" %in% names(valid_data_filtered) && nrow(valid_data_filtered) > 0) {
    res_vec <- valid_data_filtered$observed - valid_data_filtered$predicted
    if (length(res_vec) == nrow(valid_data_filtered)) {
      valid_data_filtered$residual <- res_vec
    }
  }
  
  # Calculate plot limits
  x_range <- range(valid_data_filtered$observed, na.rm = TRUE)
  y_range <- range(valid_data_filtered$predicted, na.rm = TRUE)
  plot_limits <- range(c(x_range, y_range))
  
  # Create observation vs predicted plot
  obs_pred_plot <- ggplot2::ggplot(valid_data_filtered, ggplot2::aes(x = observed, y = predicted)) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    ggplot2::labs(title = sprintf("Observed vs Predicted Values (n_splits=%d)", n_splits),
                 x = "Observed",
                 y = "Predicted") +
    ggplot2::coord_equal(xlim = plot_limits, ylim = plot_limits) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", color = NA),
                   plot.background = ggplot2::element_rect(fill = "white", color = NA))
  
  # Create residual plot
  residual_plot <- ggplot2::ggplot(valid_data_filtered, ggplot2::aes(x = predicted, y = residual)) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    ggplot2::labs(title = sprintf("Residuals vs Predicted Values (n_splits=%d)", n_splits),
                 x = "Predicted",
                 y = "Residuals") +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", color = NA),
                   plot.background = ggplot2::element_rect(fill = "white", color = NA))
  
  # Save plots
  ggplot2::ggsave(
    filename = file.path(plot_dir, sprintf("obs_vs_pred_%d_splits.png", n_splits)),
    plot = obs_pred_plot,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  ggplot2::ggsave(
    filename = file.path(plot_dir, sprintf("residuals_%d_splits.png", n_splits)),
    plot = residual_plot,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  logging::loginfo("Saved diagnostic plots to %s", plot_dir)
}

# Update save_kriging_results to report total/percentage excluded
save_kriging_results <- function(results, subject_id, config, uk_config, formula_str, sp_data) {
  output_dir <- file.path(config$output$base_dir, subject_id, config$universal_kriging$output_dir)
  uk_data_files_dir <- file.path(output_dir, "UK_data_files")
  dir.create(uk_data_files_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Create a file to store all kriging covariate coefficients
  coef_file <- file.path(uk_data_files_dir, sprintf("%s_all_kriging_coefficients.txt", subject_id))
  sink(coef_file)
  cat("Kriging Covariate Coefficients for All Split Levels\n")
  cat("================================================\n\n")
  
  # Save summary
  summary_file <- file.path(output_dir, sprintf("%s_kriging_summary.txt", subject_id))
  sink(summary_file)
  cat("Universal Kriging Summary\n")
  cat("=======================\n\n")
  
  cat("Kriging Formula (trend model):\n")
  cat(sprintf("%s\n\n", formula_str))

  cat("Variogram Formula (detrending):\n")
  cat("response ~ X + Y + Z\n\n")
  
  cat("Variogram Parameters:\n")
  cat(sprintf("Model: %s\n", uk_config$variogram_model))
  cat(sprintf("Nugget: %.4f\n", uk_config$nugget))
  cat(sprintf("Partial Sill: %.4f\n", uk_config$psill))
  cat(sprintf("Range: %.4f\n", uk_config$range))
  cat(sprintf("Cutoff: %.4f\n", uk_config$cutoff))
  cat(sprintf("Width: %.4f\n", uk_config$width))
  cat("\n")
  
  cat("Kriging Parameters:\n")
  cat(sprintf("Maximum Neighbors (nmax): %d\n", uk_config$nmax))
  cat(sprintf("Minimum Neighbors (nmin): %d\n", uk_config$nmin))
  cat(sprintf("Maximum Distance: %.4f\n", uk_config$maxdist))
  cat(sprintf("Grid Size: %.4f\n", uk_config$grid_size))
  cat("\n")
  
  cat("Cross-validation Metrics:\n")
  total_excluded <- 0
  total_points <- 0
  
  # Store all coefficients for comparison
  all_coefficients <- list()
  
  for (n_splits in names(results$cv$metrics)) {
    block_metrics <- results$cv$metrics[[n_splits]]
    cat(sprintf("\nResults for %s splits:\n", n_splits))
    cat(sprintf("MSE: %.4f\n", mean(block_metrics$mse, na.rm=TRUE)))
    cat(sprintf("RMSE: %.4f\n", mean(block_metrics$rmse, na.rm=TRUE)))
    cat(sprintf("MAE: %.4f\n", mean(block_metrics$mae, na.rm=TRUE)))
    cat(sprintf("R-squared: %.4f\n", mean(block_metrics$r2, na.rm=TRUE)))
    cat(sprintf("Total predictions: %d\n", sum(block_metrics$n_total, na.rm=TRUE)))
    cat(sprintf("Predictions excluded (>2SD): %d (%.2f%%)\n", sum(block_metrics$n_excluded, na.rm=TRUE),
                100 * sum(block_metrics$n_excluded, na.rm=TRUE) / max(1, sum(block_metrics$n_total, na.rm=TRUE))))
    total_excluded <- total_excluded + sum(block_metrics$n_excluded, na.rm=TRUE)
    total_points <- total_points + sum(block_metrics$n_total, na.rm=TRUE)
  }
  cat("\nOverall Exclusion Summary:\n")
  cat(sprintf("Total predictions excluded (>2SD): %d\n", total_excluded))
  cat(sprintf("Percentage excluded: %.2f%%\n", 100 * total_excluded / max(1, total_points)))
  sink()
  
  # Save cross-validation results and create diagnostic plots
  for (n_splits in names(results$cv$cv)) {
    cv_results <- results$cv$cv[[n_splits]]
    
    # Coordinates are now included directly in cv_results, no need to add them
    # Remove the problematic coordinate assignment that was causing mismatches
    
    # Save CV results with X, Y, Z columns
    cv_file <- file.path(uk_data_files_dir, sprintf("%s_cv_%s_splits.csv", subject_id, n_splits))
    write.csv(cv_results, cv_file, row.names = FALSE)
    
    # Create diagnostic plots for this split
    if (nrow(cv_results) > 0 && sum(!is.na(cv_results$predicted)) > 0) {
      create_diagnostic_plots(cv_results, subject_id, config, as.integer(n_splits), plot_dir = file.path(uk_data_files_dir, "plots"))
    }
    
    # Save metrics
    metrics_file <- file.path(uk_data_files_dir, sprintf("%s_cv_%s_splits_metrics.csv", subject_id, n_splits))
    write.csv(results$cv$metrics[[n_splits]], metrics_file, row.names = FALSE)

    # --- Outlier analysis for predictions ---
    pred_df <- as.data.frame(cv_results)
    # Calculate residual if not present
    if (!"residual" %in% names(pred_df) && nrow(pred_df) > 0) {
      res_vec <- pred_df$observed - pred_df$predicted
      if (length(res_vec) == nrow(pred_df)) {
        pred_df$residual <- res_vec
      }
    }
    # Identify outliers: predicted outside -10 to +10 (for outlier CSVs)
    valid_pred <- !is.na(pred_df$predicted)
    outlier_mask <- pred_df$predicted < -10 | pred_df$predicted > 10
    exclude_mask <- valid_pred & outlier_mask
    pred_df$excluded <- FALSE
    pred_df$excluded[exclude_mask] <- TRUE
    outlier_df <- pred_df[outlier_mask & valid_pred, , drop = FALSE]
    
    if (nrow(outlier_df) > 0) {
      # Add all covariates from sp_data if not present
      for (col in setdiff(names(sp_data), names(outlier_df))) {
        outlier_df[[col]] <- sp_data[[col]][as.numeric(rownames(outlier_df))]
      }
      # Ensure observed == AT8_intensity for all outliers
      if ("AT8_intensity" %in% names(sp_data)) {
        outlier_df$observed <- sp_data$AT8_intensity[as.numeric(rownames(outlier_df))]
      }
      # Calculate residual if not present
      if (!"residual" %in% names(outlier_df) && nrow(outlier_df) > 0) {
        res_vec <- outlier_df$observed - outlier_df$predicted
        if (length(res_vec) == nrow(outlier_df)) {
          outlier_df$residual <- res_vec
        }
      }
      # Save outlier CSV
      outlier_csv <- file.path(uk_data_files_dir, sprintf("%s_outlier_predictions_%s_splits.csv", subject_id, n_splits))
      write.csv(outlier_df, outlier_csv, row.names = FALSE)
    }

    # --- Fit final model for this split and extract trend coefficients ---
    # Use the full data and block assignments for this split
    # Fit variogram on full data
    v_full <- gstat::variogram(as.formula(formula_str), data = sp_data, locations = ~X+Y+Z)
    v_fit_full <- gstat::fit.variogram(v_full, gstat::vgm(uk_config$psill, uk_config$variogram_model, uk_config$range, uk_config$nugget))
    g_full <- gstat::gstat(
      formula = as.formula(formula_str),
      data = sp_data,
      model = v_fit_full,
      locations = ~X+Y+Z,
      maxdist = as.numeric(uk_config$maxdist),
      nmax = as.numeric(uk_config$nmax),
      nmin = as.numeric(uk_config$nmin)
    )
    
    # Predict on full data to get trend coefficients
    krige_result_full <- predict(g_full, newdata = sp_data, BLUE = TRUE)
    trend_coefs <- attr(krige_result_full, "beta")
    # Debugging output
    print("DEBUG: formula_str")
    print(formula_str)
    print("DEBUG: sp_data columns")
    print(names(sp_data))
    print("DEBUG: trend_coefs content:")
    print(trend_coefs)
    # Log coefficients for debugging
    logging::loginfo(sprintf("Coefficients for split %s: %s", n_splits, paste(trend_coefs, collapse=", ")))
    # Define coef_file in this scope
    coef_file <- file.path(uk_data_files_dir, sprintf("%s_all_kriging_coefficients.txt", subject_id))
    # Write coefficients to the combined file only (no individual CSVs)
    if (!is.null(trend_coefs) && length(trend_coefs) > 0) {
      cat(sprintf("\nCoefficients for %s splits:\n", n_splits), file=coef_file, append=TRUE)
      cat("--------------------------------\n", file=coef_file, append=TRUE)
      for (i in seq_along(trend_coefs)) {
        cat(sprintf("%s: %.6f\n", names(trend_coefs)[i], trend_coefs[i]), file=coef_file, append=TRUE)
      }
      cat("\n", file=coef_file, append=TRUE)
    } else {
      cat(sprintf("\nNo coefficients extracted for %s splits\n", n_splits), file=coef_file, append=TRUE)
    }

    # --- Kriging Trend Model Metrics (BLUE) ---
    krige_pred <- as.numeric(krige_result_full$var1.pred)
    krige_obs <- sp_data$response
    # Remove NA predictions for metrics
    valid_krige <- !is.na(krige_pred)
    krige_pred <- krige_pred[valid_krige]
    krige_obs <- krige_obs[valid_krige]
    krige_rmse <- sqrt(mean((krige_obs - krige_pred)^2, na.rm = TRUE))
    krige_mae <- mean(abs(krige_obs - krige_pred), na.rm = TRUE)
    krige_r2 <- 1 - sum((krige_obs - krige_pred)^2, na.rm = TRUE) / sum((krige_obs - mean(krige_obs, na.rm = TRUE))^2, na.rm = TRUE)
    krige_n <- length(krige_obs)
    krige_p <- length(trend_coefs)
    krige_sigma2 <- mean((krige_obs - krige_pred)^2, na.rm = TRUE)
    krige_rss <- sum((krige_obs - krige_pred)^2, na.rm = TRUE)
    krige_aic <- krige_n * log(krige_rss / krige_n) + 2 * krige_p
    krige_bic <- krige_n * log(krige_rss / krige_n) + log(krige_n) * krige_p
    krige_cp <- krige_rss / krige_sigma2 - krige_n + 2 * krige_p

    # --- Linear Regression Model Metrics ---
    lm_formula <- as.formula(formula_str)
    lm_fit <- lm(lm_formula, data = sp_data)
    lm_pred <- predict(lm_fit, newdata = sp_data)
    lm_obs <- sp_data$response
    valid_lm <- !is.na(lm_pred)
    lm_pred <- lm_pred[valid_lm]
    lm_obs <- lm_obs[valid_lm]
    lm_rmse <- sqrt(mean((lm_obs - lm_pred)^2, na.rm = TRUE))
    lm_mae <- mean(abs(lm_obs - lm_pred), na.rm = TRUE)
    lm_r2 <- 1 - sum((lm_obs - lm_pred)^2, na.rm = TRUE) / sum((lm_obs - mean(lm_obs, na.rm = TRUE))^2, na.rm = TRUE)
    lm_n <- length(lm_obs)
    lm_p <- length(coef(lm_fit))
    lm_sigma2 <- mean((lm_obs - lm_pred)^2, na.rm = TRUE)
    lm_rss <- sum((lm_obs - lm_pred)^2, na.rm = TRUE)
    lm_aic <- lm_n * log(lm_rss / lm_n) + 2 * lm_p
    lm_bic <- lm_n * log(lm_rss / lm_n) + log(lm_n) * lm_p
    lm_cp <- lm_rss / lm_sigma2 - lm_n + 2 * lm_p

    # Save kriging and linear regression metrics to summary
    sink(summary_file, append = TRUE)
    cat("\nKriging Trend Model (BLUE) Metrics:\n")
    cat(sprintf("RMSE: %.4f\n", krige_rmse))
    cat(sprintf("MAE: %.4f\n", krige_mae))
    cat(sprintf("R-squared: %.4f\n", krige_r2))
    cat(sprintf("AIC: %.4f\n", krige_aic))
    cat(sprintf("BIC: %.4f\n", krige_bic))
    cat(sprintf("Cp: %.4f\n", krige_cp))
    cat("\nLinear Regression Metrics:\n")
    cat(sprintf("RMSE: %.4f\n", lm_rmse))
    cat(sprintf("MAE: %.4f\n", lm_mae))
    cat(sprintf("R-squared: %.4f\n", lm_r2))
    cat(sprintf("AIC: %.4f\n", lm_aic))
    cat(sprintf("BIC: %.4f\n", lm_bic))
    cat(sprintf("Cp: %.4f\n", lm_cp))
    sink()
  }
  
  # Add coefficient comparison to the coefficients file
  sink(coef_file, append = TRUE)
  cat("\nCoefficient Comparison Across Split Levels\n")
  cat("=========================================\n\n")
  
  # Get all unique coefficient names
  all_coef_names <- unique(unlist(lapply(all_coefficients, names)))
  
  # Create comparison table
  cat("Split Level |")
  for (coef_name in all_coef_names) {
    cat(sprintf(" %-15s |", coef_name))
  }
  cat("\n")
  
  cat(paste(rep("-", 16 + 17 * length(all_coef_names)), collapse = ""))
  cat("\n")
  
  for (n_splits in names(all_coefficients)) {
    cat(sprintf("%-11s |", n_splits))
    for (coef_name in all_coef_names) {
      if (coef_name %in% names(all_coefficients[[n_splits]])) {
        cat(sprintf(" %-15.6f |", all_coefficients[[n_splits]][coef_name]))
      } else {
        cat(" NA              |")
      }
    }
    cat("\n")
  }
  sink()
  
  # Plot and save variogram if enabled
  if (uk_config$plot_variogram) {
    plot_dir <- file.path(output_dir, "plots")
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    
    tryCatch({
      v_emp <- results$variogram$empirical
      v_fit <- results$variogram$fitted
      # Check for valid data before plotting
      valid_emp <- !is.null(v_emp) && nrow(v_emp) > 0 && all(is.finite(v_emp$gamma)) && all(is.finite(v_emp$dist))
      valid_fit <- !is.null(v_fit) && all(is.finite(unlist(v_fit)))
      if (valid_emp && valid_fit) {
        # Prepare fitted line using gstat::variogramLine
        fit_line <- gstat::variogramLine(v_fit, maxdist = max(v_emp$dist), n = 200)
        fit_df <- data.frame(dist = fit_line$dist, gamma = fit_line$gamma)
        emp_df <- data.frame(dist = v_emp$dist, gamma = v_emp$gamma)
        # Extract variogram parameters
        nugget <- v_fit$psill[1]
        psill <- v_fit$psill[2]
        range <- v_fit$range[2]
        # Create annotation text
        param_text <- sprintf("Nugget: %.3f\nPartial Sill: %.3f\nRange: %.3f", 
                            nugget, psill, range)
        p <- ggplot2::ggplot() +
          ggplot2::geom_point(data = emp_df, ggplot2::aes(x = dist, y = gamma), 
                             color = 'blue', size = 2) +
          ggplot2::geom_line(data = fit_df, ggplot2::aes(x = dist, y = gamma), 
                            color = 'red', linewidth = 1) +
          ggplot2::labs(title = 'Empirical and Fitted Variogram', 
                       x = 'Distance', y = 'Semivariance') +
          ggplot2::theme_minimal() +
          ggplot2::theme(
            panel.background = ggplot2::element_rect(fill = "white", color = NA),
            plot.background = ggplot2::element_rect(fill = "white", color = NA)
          ) +
          ggplot2::annotate("text", 
                           x = max(emp_df$dist) * 0.95, 
                           y = min(emp_df$gamma) * 1.2,
                           label = param_text, 
                           hjust = 1, 
                           vjust = 0,
                           size = 4.8,
                           family = "mono",
                           fontface = "bold")
        # Save as PNG with high resolution
        ggplot2::ggsave(filename = file.path(plot_dir, sprintf("%s_variogram.png", subject_id)), 
                       plot = p, width = 8, height = 6, dpi = 300)
      } else {
        logging::logwarn("Skipping variogram plot: invalid or non-finite data.")
      }
    }, error = function(e) {
      logging::logwarn("Failed to plot variogram: %s", e$message)
    })
  }
  
  logging::loginfo("Saved kriging results to %s", output_dir)
}

# Function to prepare data for Universal Kriging
prepare_uk_data <- function(subject_id, config) {
  # Paths for input files
  transformed_data_path <- file.path(config$output$base_dir, subject_id, "preprocessing", sprintf("%s_step3_transformed_standardized.csv", subject_id))
  predictors_path <- file.path(config$output$base_dir, subject_id, "factor_analysis", sprintf("%s_predictors_reviewed.txt", subject_id))
  uk_data_path <- file.path(config$output$base_dir, subject_id, "universal_kriging", sprintf("%s_uk_data.csv", subject_id))

  # If the uk_data.csv file already exists, load and return it
  if (file.exists(uk_data_path)) {
    logging::loginfo("%s already exists, loading instead of recreating.", uk_data_path)
    uk_data <- readr::read_csv(uk_data_path, show_col_types = FALSE)
    return(uk_data)
  }

  # Check if input files exist
  if (!file.exists(transformed_data_path)) {
    stop(sprintf("Transformed data file not found at %s", transformed_data_path))
  }
  if (!file.exists(predictors_path)) {
    stop(sprintf("Predictors file not found at %s", predictors_path))
  }

  # Load transformed data
  transformed_data <- readr::read_csv(transformed_data_path, show_col_types = FALSE)

  # Read predictor names as a character vector (one per line, no header)
  predictor_names <- readr::read_lines(predictors_path)
  predictor_names <- predictor_names[predictor_names != ""] # Remove empty lines

  # Columns to keep: default + predictors
  default_cols <- c("X", "Y", "Z", "AT8_intensity", "sub.id")
  all_cols <- unique(c(default_cols, predictor_names))

  # Subset the columns from transformed_data
  uk_data <- transformed_data[, all_cols, drop = FALSE]

  # Add a 'response' column as a copy of AT8_intensity
  uk_data$response <- uk_data$AT8_intensity

  # Save combined data
  readr::write_csv(uk_data, uk_data_path)
  logging::loginfo("Created uk_data.csv at %s", uk_data_path)
  return(uk_data)
}

# Main function to run Universal Kriging
run_universal_kriging <- function(subject_id, config) {
  logging::loginfo("Starting Universal Kriging for subject %s", subject_id)
  
  # Load subject-specific UK config
  uk_config <- load_uk_config(subject_id, config)
  
  # Prepare data
  data <- prepare_uk_data(subject_id, config)
  
  # Define output directories
  output_dir <- file.path(config$output$base_dir, subject_id, config$universal_kriging$output_dir)
  uk_data_files_dir <- file.path(output_dir, "UK_data_files")
  dir.create(uk_data_files_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Extract coordinates and create spatial data
  coords <- data[, c("X", "Y", "Z")]
  sp_data <- prepare_spatial_data(data, coords)
  
  # Ensure all relevant columns are numeric in sp_data
  sp_data[] <- lapply(sp_data, function(x) as.numeric(as.character(x)))

  # Update predictor columns to exclude X, Y, Z
  all_data_cols <- names(sp_data)
  cols_to_exclude_from_predictors <- c("response", "AT8_intensity", "sub.id", "X", "Y", "Z")
  potential_predictor_cols <- setdiff(all_data_cols, cols_to_exclude_from_predictors)
  actual_predictor_cols <- potential_predictor_cols[!grepl("block_", potential_predictor_cols)]
  # Remove X_std, Y_std, Z_std if present
  actual_predictor_cols <- setdiff(actual_predictor_cols, c("X_std", "Y_std", "Z_std"))
  # Insert predictor checks here
  # Remove predictors with zero or near-zero variance
  nzv <- sapply(sp_data[, actual_predictor_cols, drop=FALSE], function(x) sd(x, na.rm=TRUE) < 1e-8)
  if (any(nzv)) {
    logging::logwarn("Removing predictors with near-zero variance: %s", paste(actual_predictor_cols[nzv], collapse=","))
    actual_predictor_cols <- actual_predictor_cols[!nzv]
  }
  # Remove predictors with NA values
  na_cols <- sapply(sp_data[, actual_predictor_cols, drop=FALSE], function(x) any(is.na(x)))
  if (any(na_cols)) {
    logging::logwarn("Removing predictors with NA values: %s", paste(actual_predictor_cols[na_cols], collapse=","))
    actual_predictor_cols <- actual_predictor_cols[!na_cols]
  }
  # Remove linearly dependent predictors using caret
  if (!requireNamespace("caret", quietly = TRUE)) install.packages("caret")
  library(caret)
  mm <- model.matrix(as.formula(paste("response ~", paste(actual_predictor_cols, collapse = " + "))), data = sp_data)
  combos <- findLinearCombos(mm)
  if (!is.null(combos$remove) && length(combos$remove) > 0) {
    logging::logwarn("Removing linearly dependent predictors: %s", paste(colnames(mm)[combos$remove], collapse=","))
    actual_predictor_cols <- actual_predictor_cols[-(combos$remove - 1)] # -1 for intercept
  }
  # Rebuild formula
  formula_str <- paste("response ~", paste(actual_predictor_cols, collapse = " + "))
  formula <- as.formula(formula_str)
  logging::loginfo("Final predictors used for kriging: %s", paste(actual_predictor_cols, collapse=","))
  variogram_formula <- formula
  
  # Fit initial variogram using the full dataset (no covariates)
  variogram <- fit_variogram(sp_data, uk_config, subject_id, config, variogram_formula)
  
  # Plot initial variogram if enabled
  if (uk_config$plot_variogram) {
    plot_initial_variogram(variogram$empirical, subject_id, uk_config, config)
  }
  
  # Initialize lists to store results for each n_split value
  all_cv_results_list <- list()
  all_cv_metrics_list <- list()
  all_cv_block_metrics_list <- list()
  
  # Process each n_split value
  if (length(uk_config$n_splits) == 0) {
    logging::logwarn("No n_splits values defined in uk_config. Skipping cross-validation.")
  } else {
    ordered_n_splits <- sort(uk_config$n_splits, decreasing = FALSE)
    logging::loginfo("Starting cross-validation loop for n_splits: %s (ordered)", 
                    paste(ordered_n_splits, collapse=", "))
    for (current_n_split_value in ordered_n_splits) {
      tryCatch({
        logging::loginfo("Processing %d splits", current_n_split_value)
        # Get block assignments for current split value
        block_assignments <- create_spatial_blocks(sp_data, current_n_split_value, subject_id, config)
        # Ensure block_assignments is numeric
        block_assignments <- as.numeric(block_assignments)
        # Perform cross-validation
        cv_results <- perform_cross_validation(sp_data, block_assignments, formula, uk_config, current_n_split_value, subject_id, config, current_n_split_value)
        # Store results
        split_key <- as.character(current_n_split_value)
        all_cv_results_list[[split_key]] <- cv_results$pred
        all_cv_metrics_list[[split_key]] <- cv_results$metrics
        all_cv_block_metrics_list[[split_key]] <- cv_results$metrics
        # Log progress
        logging::loginfo("Completed CV for %d splits", current_n_split_value)
        # Only log metrics if they exist
        if (!is.null(cv_results$metrics) && 
            any(!is.na(cv_results$metrics$rmse)) && 
            any(!is.na(cv_results$metrics$r2))) {
          logging::loginfo(
            "RMSE (mean across blocks): %.4f, R2 (mean across blocks): %.4f",
            mean(cv_results$metrics$rmse, na.rm = TRUE),
            mean(cv_results$metrics$r2, na.rm = TRUE)
          )
        } else {
          logging::loginfo("No valid metrics available (all predictions failed)")
        }
        # Call plot_obs_pred_residuals, but exclude rows with any predictor >3SD from mean
        # Compute mask for rows with any predictor >3SD from mean
        predictor_cols <- actual_predictor_cols
        pred_mask <- rep(FALSE, nrow(sp_data))
        for (col in predictor_cols) {
          vals <- sp_data[[col]]
          mu <- mean(vals, na.rm = TRUE)
          sigma <- sd(vals, na.rm = TRUE)
          pred_mask <- pred_mask | (abs(sp_data[[col]] - mu) > 3 * sigma)
        }
        # Exclude these rows from plotting
        plot_df <- as.data.frame(cv_results$pred)
        if ("row.names" %in% colnames(plot_df)) rownames(plot_df) <- plot_df$row.names
        # Try to align indices if possible
        if (nrow(plot_df) == nrow(sp_data)) {
          plot_df <- plot_df[!pred_mask, , drop = FALSE]
        } else if ("row.names" %in% colnames(plot_df)) {
          plot_df <- plot_df[!pred_mask[as.numeric(rownames(plot_df))], , drop = FALSE]
        }
        plot_obs_pred_residuals(plot_df, subject_id, current_n_split_value, config)

        # --- Outlier analysis for predictions ---
        pred_df <- as.data.frame(cv_results$pred)
        # Calculate residual if not present
        if (!"residual" %in% names(pred_df) && nrow(pred_df) > 0) {
          res_vec <- pred_df$observed - pred_df$predicted
          if (length(res_vec) == nrow(pred_df)) {
            pred_df$residual <- res_vec
          }
        }
        # Identify outliers: predicted outside -10 to +10 (for outlier CSVs)
        valid_pred <- !is.na(pred_df$predicted)
        outlier_mask <- pred_df$predicted < -10 | pred_df$predicted > 10
        exclude_mask <- valid_pred & outlier_mask
        pred_df$excluded <- FALSE
        pred_df$excluded[exclude_mask] <- TRUE
        outlier_df <- pred_df[outlier_mask & valid_pred, , drop = FALSE]
        # Add X, Y, Z, and all covariates if available
        # Try to merge with sp_data by row index if needed
        if (nrow(outlier_df) > 0) {
          # If X, Y, Z not present, try to add from sp_data
          for (coord in c("X", "Y", "Z")) {
            if (!coord %in% names(outlier_df) && coord %in% names(sp_data)) {
              outlier_df[[coord]] <- sp_data[[coord]][as.numeric(rownames(outlier_df))]
            }
          }
          # Add all covariates from sp_data if not present
          for (col in setdiff(names(sp_data), names(outlier_df))) {
            outlier_df[[col]] <- sp_data[[col]][as.numeric(rownames(outlier_df))]
          }
          # Ensure observed == AT8_intensity for all outliers
          if ("AT8_intensity" %in% names(sp_data)) {
            outlier_df$observed <- sp_data$AT8_intensity[as.numeric(rownames(outlier_df))]
          }
          # Calculate residual if not present
          if (!"residual" %in% names(outlier_df) && nrow(outlier_df) > 0) {
            res_vec <- outlier_df$observed - outlier_df$predicted
            if (length(res_vec) == nrow(outlier_df)) {
              outlier_df$residual <- res_vec
            }
          }
          # Save outlier CSV
          outlier_csv <- file.path(uk_data_files_dir, sprintf("%s_outlier_predictions_%s_splits.csv", subject_id, current_n_split_value))
          write.csv(outlier_df, outlier_csv, row.names = FALSE)
        }

        # --- Fit final model for this split and extract trend coefficients ---
        # Use the full data and block assignments for this split
        # Fit variogram on full data
        v_full <- gstat::variogram(as.formula(formula_str), data = sp_data, locations = ~X+Y+Z)
        v_fit_full <- gstat::fit.variogram(v_full, gstat::vgm(uk_config$psill, uk_config$variogram_model, uk_config$range, uk_config$nugget))
        g_full <- gstat::gstat(
          formula = as.formula(formula_str),
          data = sp_data,
          model = v_fit_full,
          locations = ~X+Y+Z,
          maxdist = as.numeric(uk_config$maxdist),
          nmax = as.numeric(uk_config$nmax),
          nmin = as.numeric(uk_config$nmin)
        )
        # Predict on full data to get trend coefficients
        krige_result_full <- predict(g_full, newdata = sp_data, BLUE = TRUE)
        trend_coefs <- attr(krige_result_full, "beta")
        # Debugging output
        print("DEBUG: formula_str")
        print(formula_str)
        print("DEBUG: sp_data columns")
        print(names(sp_data))
        print("DEBUG: trend_coefs content:")
        print(trend_coefs)
        # Log coefficients for debugging
        logging::loginfo(sprintf("Coefficients for split %s: %s", current_n_split_value, paste(trend_coefs, collapse=", ")))
        # Define coef_file in this scope
        coef_file <- file.path(uk_data_files_dir, sprintf("%s_all_kriging_coefficients.txt", subject_id))
        # Write coefficients to the combined file only (no individual CSVs)
        if (!is.null(trend_coefs) && length(trend_coefs) > 0) {
          cat(sprintf("\nCoefficients for %s splits:\n", current_n_split_value), file=coef_file, append=TRUE)
          cat("--------------------------------\n", file=coef_file, append=TRUE)
          for (i in seq_along(trend_coefs)) {
            cat(sprintf("%s: %.6f\n", names(trend_coefs)[i], trend_coefs[i]), file=coef_file, append=TRUE)
          }
          cat("\n", file=coef_file, append=TRUE)
        } else {
          cat(sprintf("\nNo coefficients extracted for %s splits\n", current_n_split_value), file=coef_file, append=TRUE)
        }
      }, error = function(e) {
        logging::logerror("Error in split %d for subject %s: %s", current_n_split_value, subject_id, e$message)
        split_key <- as.character(current_n_split_value)
        all_cv_results_list[[split_key]] <- data.frame()
        all_cv_metrics_list[[split_key]] <- data.frame()
        all_cv_block_metrics_list[[split_key]] <- data.frame()
      })
    }
  }
  
  # Save only cross-validation results (no final kriging prediction)
  consolidated_cv_results <- list(
    cv = all_cv_results_list,
    metrics = all_cv_metrics_list,
    block_metrics = all_cv_block_metrics_list
  )
  
  # Save results (no grid or krige)
  results <- list(
    variogram = variogram,
    cv = consolidated_cv_results
  )
  # Save cross-validation metrics and summary
  save_kriging_results(results, subject_id, config, uk_config, formula_str, sp_data)
  
  logging::loginfo("Completed Universal Kriging for subject %s", subject_id)
  return(results)
}

#' Visualize Split Blocks as NIfTI
#' 
#' This function takes split blocks and converts them into a NIfTI file where each block
#' is represented by a unique label (1 to n_blocks).
#' 
#' @param blocks List of split blocks, each containing spatial coordinates
#' @param n_blocks Number of blocks (2, 4, 8, or 16)
#' @param output_dir Directory to save the NIfTI file
#' @param subject_id Subject ID for the output filename
#' @return Path to the created NIfTI file
visualize_split_blocks <- function(blocks, n_blocks, output_dir, subject_id) {
  if (!require("RNifti")) {
    install.packages("RNifti")
    library(RNifti)
  }
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Try to load reference NIfTI, use default if not found
  # Path relative to project root (set by run_pipeline.R)
  ref_path <- file.path("input", "processed", "mni_icbm152_t1_nlin_sym_09b_hires_stripped.nii.gz")
  
  if (file.exists(ref_path)) {
    # Load reference NIfTI
    ref_nifti <- tryCatch({
      RNifti::readNifti(ref_path)
    }, error = function(e) {
      logging::logwarn("Failed to read reference NIfTI: %s. Using default parameters.", e$message)
      NULL
    })
    
    if (!is.null(ref_nifti)) {
      # Get affine matrix and dimensions
      affine <- RNifti::xform(ref_nifti)
      ref_dim <- dim(ref_nifti)
      logging::loginfo("Using reference NIfTI dimensions: %s", paste(ref_dim, collapse="x"))
    } else {
      ref_nifti <- NULL
    }
  } else {
    logging::logwarn("Reference NIfTI not found at %s. Using default parameters.", ref_path)
    ref_nifti <- NULL
  }
  
  # If reference NIfTI not available, use default parameters
  if (is.null(ref_nifti)) {
    # Create default affine matrix (1mm isotropic, centered)
    affine <- diag(4)
    affine[1:3, 4] <- c(-90, -126, -72)  # Standard MNI space origin
    ref_dim <- c(182, 218, 182)  # Standard MNI dimensions
    logging::loginfo("Using default dimensions: %s", paste(ref_dim, collapse="x"))
  }
  
  # Validate affine matrix
  if (!is.matrix(affine) || nrow(affine) != 4 || ncol(affine) != 4) {
    stop("Invalid affine matrix")
  }
  
  # Calculate inverse affine for world to voxel transformation
  inv_affine <- try(solve(affine))
  if (inherits(inv_affine, "try-error")) {
    stop("Failed to compute inverse affine matrix")
  }

  # Create array with reference dimensions
  block_array <- array(0L, dim = ref_dim)  # Integer array
  
  # Function to convert world coordinates to voxel indices
  world_to_voxel <- function(coords) {
    # Add homogeneous coordinate
    coords_hom <- cbind(coords, 1)
    # Transform to voxel space
    ijk <- t(inv_affine %*% t(coords_hom))
    # Round to nearest voxel
    ijk <- round(ijk[, 1:3])
    return(ijk)
  }

  # Process each block
  for (i in seq_along(blocks)) {
    block <- blocks[[i]]
    
    # Generate grid of points within block bounds
    x_seq <- seq(block$x[1], block$x[2], length.out = 50)
    y_seq <- seq(block$y[1], block$y[2], length.out = 50)
    z_seq <- seq(block$z[1], block$z[2], length.out = 50)
    
    # Create meshgrid
    grid <- expand.grid(x = x_seq, y = y_seq, z = z_seq)
    
    # Convert to voxel coordinates
    voxel_coords <- world_to_voxel(as.matrix(grid))
    
    # Filter valid voxel coordinates
    valid_idx <- voxel_coords[,1] >= 1 & voxel_coords[,1] <= ref_dim[1] &
                 voxel_coords[,2] >= 1 & voxel_coords[,2] <= ref_dim[2] &
                 voxel_coords[,3] >= 1 & voxel_coords[,3] <= ref_dim[3]
    
    if (any(valid_idx)) {
      valid_coords <- voxel_coords[valid_idx, , drop = FALSE]
      
      # Assign block number to valid voxels
      for (j in seq_len(nrow(valid_coords))) {
        x <- valid_coords[j, 1]
        y <- valid_coords[j, 2]
        z <- valid_coords[j, 3]
        block_array[x, y, z] <- i
      }
      
      logging::loginfo("Block %d: Filled %d voxels", i, sum(valid_idx))
    } else {
      logging::logwarn("Block %d: No valid voxels found", i)
    }
  }

  # Create NIfTI object
  nifti_img <- RNifti::asNifti(block_array)
  
  # Set the transformation matrix
  RNifti::sform(nifti_img) <- affine
  
  # Save NIfTI file
  output_file <- file.path(output_dir, sprintf("%s_split_blocks_%d.nii.gz", subject_id, n_blocks))
  tryCatch({
    RNifti::writeNifti(nifti_img, output_file)
    logging::loginfo("Successfully saved block visualization to %s", output_file)
  }, error = function(e) {
    stop(sprintf("Failed to save NIfTI file: %s", e$message))
  })

  return(output_file)
} 