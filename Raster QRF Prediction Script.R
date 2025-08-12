# Clear workspace
rm(list = ls())

library(raster)
library(caret)
library(quantregForest)
library(stringr)

setwd("G:/")

regions <- c("China")
sampling_data_template <- "XX/{region}_BD{depth}_combined_data.csv"
raster_stack_template <- "XXX/Block/Forests_Raster_stack_item/" #The folder where the grid stack is located
best_features_template <- "XXXX/{region}_BD{depth}_best_features.rds"
checkpoint_template <- "XXXXX/{region}_BD{depth}_checkpoint.rds"
depths <- c("05") #Depth of layer 05 515 1530 3060 60100

region_info <- lapply(regions, function(region) {
  raster_stack_path <- gsub("\\{region\\}", region, raster_stack_template)
  sampling_data <- sapply(depths, function(depth) {
    gsub("\\{region\\}", region, gsub("\\{depth\\}", depth, sampling_data_template))
  })
  best_features <- sapply(depths, function(depth) {
    gsub("\\{region\\}", region, gsub("\\{depth\\}", depth, best_features_template))
  })
  checkpoints <- sapply(depths, function(depth) {
    gsub("\\{region\\}", region, gsub("\\{depth\\}", depth, checkpoint_template))
  })
  list(
    region = region,
    sampling_data = sampling_data,
    raster_stack_path = raster_stack_path,
    best_features = best_features,
    checkpoints = checkpoints
  )
})

calculate_model_performance <- function(actual, predicted) {
  rsq <- cor(actual, predicted)^2
  rmse <- sqrt(mean((actual - predicted)^2))
  mae <- mean(abs(actual - predicted))
  return(list(R2 = rsq, RMSE = rmse, MAE = mae))
}

merge_rasters <- function(raster_list) {
  if (length(raster_list) == 0) stop("Raster list is empty")
  if (length(raster_list) == 1) return(raster_list[[1]])
  merged_raster <- raster_list[[1]]
  for (i in 2:length(raster_list)) {
    merged_raster <- tryCatch({
      mosaic(merged_raster, raster_list[[i]], fun = mean, tolerance = 0.1, na.rm = TRUE)
    }, error = function(e) {
      warning(paste("Error merging raster", i, ":", e$message))
      return(merged_raster)
    })
  }
  return(merged_raster)
}

save_block_results <- function(pred_raster, lower_raster, upper_raster, cache_dir, output_prefix, file_suffix) {
  writeRaster(pred_raster, file.path(cache_dir, paste0(output_prefix, "_", file_suffix, "_pred.tif")), format = "GTiff", overwrite = TRUE)
  writeRaster(lower_raster, file.path(cache_dir, paste0(output_prefix, "_", file_suffix, "_lower.tif")), format = "GTiff", overwrite = TRUE)
  writeRaster(upper_raster, file.path(cache_dir, paste0(output_prefix, "_", file_suffix, "_upper.tif")), format = "GTiff", overwrite = TRUE)
  cat("Saved block results:", paste0(output_prefix, "_", file_suffix, "_pred.tif"), paste0(output_prefix, "_", file_suffix, "_lower.tif"), paste0(output_prefix, "_", file_suffix, "_upper.tif"), "\n")
}

initialize_checkpoint <- function(checkpoint_file) {
  if (file.exists(checkpoint_file)) {
    cat("Loading checkpoint file:", checkpoint_file, "\n")
    checkpoint <- tryCatch({
      readRDS(checkpoint_file)
    }, error = function(e) {
      cat("Error reading checkpoint file:", e$message, "\n")
      cat("Deleting corrupted checkpoint:", checkpoint_file, "\n")
      file.remove(checkpoint_file)
      checkpoint <- list()
    })
  } else {
    checkpoint <- list()
  }
  if (!is.null(checkpoint$processed_files) && !all(sapply(checkpoint$processed_files, is.character))) {
    checkpoint$processed_files <- lapply(checkpoint$processed_files, as.character)
  }
  return(checkpoint)
}

optimized_model_results <- read.csv("G:/pH_BD_maaping/3.QRF/2.mosaic/total_optimized_model_results.csv", stringsAsFactors = FALSE)

for (region_data in region_info) {
  current_region <- region_data$region
  cat("Processing region:", current_region, "\n")
  pred_rasters <- list()
  lower_rasters <- list()
  upper_rasters <- list()
  for (i in 1:length(region_data$sampling_data)) {
    current_sampling_file <- region_data$sampling_data[i]
    output_prefix <- sub("_combined_data.csv", "", basename(current_sampling_file))
    samples <- read.csv(current_sampling_file)
    best_features <- readRDS(region_data$best_features[i])
    optimized_model_results$Depth <- sprintf("%02d", as.numeric(optimized_model_results$Depth))
    current_depth <- depths[i]
    response_variable <- "BD"
    optimized_params <- optimized_model_results[
      optimized_model_results$Region == current_region &
        optimized_model_results$Depth == current_depth &
        optimized_model_results$Response == response_variable, ]
    if (nrow(optimized_params) == 0) stop(paste("No optimized parameter found for", current_region, current_depth, response_variable))
    best_mtry <- optimized_params$Best_mtry
    best_num_trees <- optimized_params$Best_num.trees
    best_min_node_size <- optimized_params$Best_min.node.size
    response_variable <- optimized_params$Response
    final_qrf_model <- quantregForest(
      x = samples[, best_features],
      y = samples[["y"]],
      ntree = best_num_trees,
      mtry = best_mtry,
      nodesize = best_min_node_size
    )
    processed_files <- list.files(
      path = "Z:/ccc_Modeling_Regional_90m/Block/Forests_Raster_stack_item/",
      pattern = "^raster_stack_Forests_block_\\d+\\.tif$",
      full.names = TRUE
    )
    if (length(processed_files) == 0) stop(paste("No raster stack file found at:", region_data$raster_stack_path))
    cat("Found raster stack files:", processed_files, "\n")
    checkpoint_file <- region_data$checkpoints[i]
    checkpoint <- initialize_checkpoint(checkpoint_file)
    for (processed_file in processed_files) {
      if (!is.null(checkpoint$processed_files) && basename(processed_file) %in% checkpoint$processed_files) next
      file_suffix <- gsub("raster_stack_|\\.tif$", "", basename(processed_file))
      pred_file <- file.path(cache_dir, paste0(output_prefix, "_", file_suffix, "_pred.tif"))
      if (pred_file %in% processed_files) {
        cat("File exists, skipping:", pred_file, "\n")
        pred_rasters <- append(pred_rasters, list(raster(pred_file)))
        lower_rasters <- append(lower_rasters, list(raster(file.path(cache_dir, paste0(output_prefix, "_", file_suffix, "_lower.tif")))))
        upper_rasters <- append(upper_rasters, list(raster(file.path(cache_dir, paste0(output_prefix, "_", file_suffix, "_upper.tif")))))
        checkpoint$processed_files <- append(checkpoint$processed_files, basename(processed_file))
        next
      }
      checkpoint$processed_file <- basename(processed_file)
      cat("Processing file:", processed_file, "\n")
      rst <- stack(processed_file)
      rds_file <- gsub("raster_stack_Forests_block_(\\d+)\\.tif$", "raster_stack_names_Forests_block_\\1.rds", processed_file)
      if (file.exists(rds_file)) {
        cat("Loading layer names:", rds_file, "\n")
        layer_names <- readRDS(rds_file)
        names(rst) <- layer_names
      } else {
        stop(paste("No corresponding .rds file found:", rds_file))
      }
      missing_features <- setdiff(best_features, names(rst))
      if (length(missing_features) > 0) stop(paste("Missing features:", paste(missing_features, collapse = ", ")))
      rst_selected <- rst[[best_features]]
      set.seed(666)
      predicted_raster <- predict(rst_selected, final_qrf_model, what = 0.5)
      predicted_raster_10 <- predict(rst_selected, final_qrf_model, what = 0.05)
      predicted_raster_90 <- predict(rst_selected, final_qrf_model, what = 0.95)
      save_block_results(predicted_raster, predicted_raster_10, predicted_raster_90, cache_dir, output_prefix, file_suffix)
      pred_rasters <- append(pred_rasters, list(predicted_raster))
      lower_rasters <- append(lower_rasters, list(predicted_raster_10))
      upper_rasters <- append(upper_rasters, list(predicted_raster_90))
      checkpoint$pred_rasters <- pred_rasters
      checkpoint$lower_rasters <- lower_rasters
      checkpoint$upper_rasters <- upper_rasters
      if (is.null(checkpoint$processed_files)) checkpoint$processed_files <- list()
      checkpoint$processed_files <- append(checkpoint$processed_files, basename(processed_file))
      saveRDS(checkpoint, checkpoint_file)
      checkpoint$processed_file <- NULL
    }
    if (!dir.exists(new_output_dir)) dir.create(new_output_dir, recursive = TRUE)
    if (length(pred_rasters) > 0) {
      pred_mosaic <- merge_rasters(pred_rasters)
      lower_mosaic <- merge_rasters(lower_rasters)
      upper_mosaic <- merge_rasters(upper_rasters)
      writeRaster(pred_mosaic, file.path(new_output_dir, paste0(output_prefix, "_pred.tif")), format = "GTiff", overwrite = TRUE)
      writeRaster(lower_mosaic, file.path(new_output_dir, paste0(output_prefix, "_lower.tif")), format = "GTiff", overwrite = TRUE)
      writeRaster(upper_mosaic, file.path(new_output_dir, paste0(output_prefix, "_upper.tif")), format = "GTiff", overwrite = TRUE)
      cat("Merged files saved:", file.path(new_output_dir, paste0(output_prefix, "_pred.tif")),
          file.path(new_output_dir, paste0(output_prefix, "_lower.tif")),
          file.path(new_output_dir, paste0(output_prefix, "_upper.tif")), "\n")
    }
  }
}
