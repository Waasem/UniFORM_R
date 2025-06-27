#' UniFORM Normalization Functions
#' 
#' R implementation of the UniFORM normalization pipeline

library(signal)
library(ggplot2)

correlation_based_normalization <- function(ref_hist, hist_list) {
  shifts_direct <- c()
  shifts_fft <- c()
  
  for (i in seq_along(hist_list)) {
    hist <- hist_list[[i]]
    
    # Cross-correlation using signal package
    correlation_direct <- ccf(as.vector(hist), as.vector(ref_hist), 
                             plot = FALSE, lag.max = length(ref_hist) - 1)
    correlation_fft <- correlation_direct  # R's ccf uses FFT by default
    
    shift_direct <- which.max(correlation_direct$acf) - length(ref_hist)
    shift_fft <- shift_direct
    
    shifts_direct <- c(shifts_direct, shift_direct)
    shifts_fft <- c(shifts_fft, shift_fft)
  }
  
  return(list(
    shifts_direct = shifts_direct,
    shifts_fft = shifts_fft
  ))
}

landmark_shift <- function(fd_list, landmarks, location = NULL) {
  landmarks <- as.numeric(landmarks)
  
  if (length(landmarks) != length(fd_list)) {
    stop("Landmark list must have the same length as the number of samples")
  }
  
  if (is.null(location)) {
    loc_array <- mean(landmarks)
  } else {
    loc_array <- as.numeric(location)
  }
  
  return(landmarks - loc_array)
}

normalize_and_plot_distributions <- function(adata_list, histograms, markers, 
                                           reference_sample, landmarks = NULL,
                                           bin_counts = 1024, xlim = NULL, 
                                           colormap = "Set1", figsize = c(13, 5)) {
  
  sample_names <- unique(adata_list$obs$sample_id)
  t <- seq(0, bin_counts - 1, length.out = bin_counts)
  shifts_fft_dict <- list()
  
  for (i in seq_along(markers)) {
    marker <- markers[i]
    reference_index <- which(sample_names == reference_sample[[marker]])
    
    cat("Processing marker", marker, "\n")
    cat("Reference sample for", marker, "is", reference_sample[[marker]], "\n")
    
    hist_data <- histograms[[marker]]$hist_list
    ref_hist <- hist_data[[reference_index]]
    
    if (!is.null(landmarks) && marker %in% names(landmarks)) {
      marker_landmark <- landmarks[[marker]]
      if (!is.null(marker_landmark)) {
        cat("Performing landmark finetuning...\n")
        shift <- landmark_shift(hist_data, marker_landmark, marker_landmark[reference_index])
      } else {
        cat("Performing automatic normalization...\n")
        corr_result <- correlation_based_normalization(ref_hist, hist_data)
        shift <- corr_result$shifts_fft
      }
    } else {
      cat("Performing automatic normalization...\n")
      corr_result <- correlation_based_normalization(ref_hist, hist_data)
      shift <- corr_result$shifts_fft
    }
    
    shifts_fft_dict[[marker]] <- shift
    cat("Shifts for", marker, ":", paste(shift, collapse = ", "), "\n\n")
  }
  
  return(shifts_fft_dict)
}

calculate_shift_in_log_pixels <- function(data_range, shifts_fft_dict, bin_counts = 1024) {
  shift_in_log_pixels_dict <- list()
  
  for (key in names(shifts_fft_dict)) {
    min_val <- data_range[[key]]$global_min
    max_val <- data_range[[key]]$global_max
    increment <- (max_val - min_val) / (bin_counts - 1)
    shifts <- shifts_fft_dict[[key]]
    shift_in_log_pixels <- shifts * increment
    shift_in_log_pixels_dict[[key]] <- shift_in_log_pixels
    cat("shift_in_log_pixels for", key, "is", paste(shift_in_log_pixels, collapse = ", "), "\n")
  }
  
  return(shift_in_log_pixels_dict)
}

generate_normalized_feature <- function(adata_list, shift_in_log_pixels_dict, 
                                      reference_sample, output_directory, 
                                      num_bins = 1024, plot_dist = FALSE, 
                                      save_normalized_features = FALSE) {
  
  # Generate timestamp for layer key
  timestamp <- format(Sys.time(), "%Y/%m/%d/%H/%M")
  layer_key <- paste0(timestamp, "_normalized")
  
  # Negate the shifts for normalization
  negated_dict <- lapply(shift_in_log_pixels_dict, function(x) -x)
  
  cat("Performing Feature Normalization\n")
  
  # Extract data components
  feature_data <- adata_list$X
  marker_list <- adata_list$var$marker_name
  markers_to_normalize <- names(negated_dict)
  sample_names <- unique(adata_list$obs$sample_id)
  
  # Initialize normalized data
  normalized_data <- feature_data
  
  for (marker in markers_to_normalize) {
    cat("Processing marker =", marker, ", reference =", reference_sample[[marker]], "\n")
    
    reference_index <- which(sample_names == reference_sample[[marker]])
    marker_index <- which(marker_list == marker)
    
    # Get shift values for current marker
    negated_list <- negated_dict[[marker]]
    
    for (sample_index in seq_along(sample_names)) {
      sample_name <- sample_names[sample_index]
      
      # Extract raw intensities for current marker and sample
      sample_mask <- adata_list$obs$sample_id == sample_name
      marker_raw <- feature_data[sample_mask, marker_index]
      shift_scaling_factor <- 10^(negated_list[sample_index])
      
      # Apply normalization
      marker_shifted <- marker_raw * shift_scaling_factor
      
      # Store normalized data
      normalized_data[sample_mask, marker_index] <- marker_shifted
      
      # Optional plotting
      if (plot_dist) {
        reference_mask <- adata_list$obs$sample_id == reference_sample[[marker]]
        reference_marker_raw <- feature_data[reference_mask, marker_index]
        
        # Create comparison plot
        plot_data <- data.frame(
          value = c(log10(marker_raw + 1), log10(marker_shifted + 1), log10(reference_marker_raw + 1)),
          type = rep(c("Original", "Normalized", "Reference"), 
                    c(length(marker_raw), length(marker_shifted), length(reference_marker_raw)))
        )
        
        p <- ggplot(plot_data, aes(x = value, fill = type)) +
          geom_histogram(alpha = 0.7, bins = num_bins, position = "identity") +
          labs(title = paste(sample_name, marker, "Cell Mean Intensity Distribution"),
               x = "Log10 Cell Mean Intensity",
               y = "Frequency") +
          theme_minimal()
        
        print(p)
      }
    }
  }
  
  # Save normalized data
  if (!exists("layers", where = adata_list)) {
    adata_list$layers <- list()
  }
  adata_list$layers[[layer_key]] <- normalized_data
  
  if (save_normalized_features) {
    if (!dir.exists(output_directory)) {
      dir.create(output_directory, recursive = TRUE)
    }
    
    output_path <- file.path(output_directory, "normalized_adata.rds")
    saveRDS(adata_list, output_path)
    cat("Normalized data saved in adata$layers under key '", layer_key, "'\n")
    cat("Updated AnnData object saved at", output_path, "\n")
  }
  
  cat("Feature Normalization Done\n\n")
  
  return(adata_list)
}