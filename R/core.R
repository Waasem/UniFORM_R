#' Calculate Histograms and GMM Analysis for UniFORM
#' 
#' @param adata_list List containing X (data matrix), obs (observations), var (variables)
#' @param marker_to_plot Character vector of markers to analyze
#' @param rep Character, either "X" or "normalized"
#' @param bin_counts Integer, number of histogram bins
#' @param xlims Named list of xlim ranges for each marker
#' @param ylims Named list of ylim ranges for each marker
#' @param colormap Character, color palette name
#' @param plot_local_threshold Logical, whether to plot local thresholds
#' @param plot_global_threshold Logical, whether to plot global threshold
#' @return List containing results_range, results_hist, and gmm_curves

library(mixtools)
library(ggplot2)
library(dplyr)
library(tidyr)

preprocess_raw <- function(image) {
  log10(image + 1)
}

plot_global_gmm <- function(X) {
  X <- as.vector(X)
  X <- X[!is.na(X) & is.finite(X)]
  
  if (length(X) < 10) {
    return(mean(X, na.rm = TRUE))
  }
  
  tryCatch({
    gmm_result <- normalmixEM(X, k = 2, maxit = 1000)
    
    # Determine which component is negative/positive based on means
    means <- gmm_result$mu
    if (means[1] > means[2]) {
      neg_component <- 2
      pos_component <- 1
    } else {
      neg_component <- 1
      pos_component <- 2
    }
    
    # Find the cutoff as the maximum of the negative component
    cut <- means[neg_component] + 2 * sqrt(gmm_result$sigma[neg_component])
    return(cut)
  }, error = function(e) {
    return(quantile(X, 0.75, na.rm = TRUE))
  })
}

plot_local_gmm <- function(X, bin_counts, global_min, global_max) {
  X <- as.vector(X)
  X <- X[!is.na(X) & is.finite(X)]
  
  if (length(X) < 10) {
    x_axis <- seq(global_min, global_max, length.out = bin_counts)
    return(list(
      cut = mean(X, na.rm = TRUE),
      x_axis = x_axis,
      y_axis0 = rep(0, bin_counts),
      y_axis1 = rep(0, bin_counts)
    ))
  }
  
  tryCatch({
    gmm_result <- normalmixEM(X, k = 2, maxit = 1000)
    
    x_axis <- seq(global_min, global_max, length.out = bin_counts)
    
    # Calculate probability density for each component
    y_axis0 <- gmm_result$lambda[1] * dnorm(x_axis, gmm_result$mu[1], sqrt(gmm_result$sigma[1]))
    y_axis1 <- gmm_result$lambda[2] * dnorm(x_axis, gmm_result$mu[2], sqrt(gmm_result$sigma[2]))
    
    # Determine negative component
    means <- gmm_result$mu
    if (means[1] > means[2]) {
      neg_component <- 2
    } else {
      neg_component <- 1
    }
    
    cut <- gmm_result$mu[neg_component] + 2 * sqrt(gmm_result$sigma[neg_component])
    
    return(list(
      cut = cut,
      x_axis = x_axis,
      y_axis0 = y_axis0,
      y_axis1 = y_axis1
    ))
  }, error = function(e) {
    x_axis <- seq(global_min, global_max, length.out = bin_counts)
    return(list(
      cut = quantile(X, 0.75, na.rm = TRUE),
      x_axis = x_axis,
      y_axis0 = rep(0, bin_counts),
      y_axis1 = rep(0, bin_counts)
    ))
  })
}

uniform_calculate_histogram <- function(adata_list, marker_to_plot, rep = "X", 
                                      bin_counts = 1024, xlims = NULL, ylims = NULL,
                                      colormap = "Set1", plot_local_threshold = TRUE,
                                      plot_global_threshold = TRUE, save_filename = NULL) {
  
  # Extract data components
  if (rep == "X") {
    data_matrix <- adata_list$X
  } else if (rep == "normalized" && "normalized" %in% names(adata_list$layers)) {
    data_matrix <- adata_list$layers$normalized
  } else {
    stop("Invalid rep value. Use 'X' or 'normalized'.")
  }
  
  sample_names <- unique(adata_list$obs$sample_id)
  marker_list <- adata_list$var$marker_name
  
  results_range <- list()
  results_hist <- list()
  gmm_curves <- list()
  
  # Create plots
  plot_list <- list()
  
  for (marker_name in marker_to_plot) {
    if (!marker_name %in% marker_list) next
    
    marker_index <- which(marker_list == marker_name)
    cat("Processing marker:", marker_name, "\n")
    
    min_list <- c()
    max_list <- c()
    global_min <- Inf
    global_max <- -Inf
    combined_data <- c()
    
    # First pass: get global min/max
    for (sample_name in sample_names) {
      sample_mask <- adata_list$obs$sample_id == sample_name
      marker_data <- data_matrix[sample_mask, marker_index]
      marker_mean_intensity_scaled <- preprocess_raw(marker_data)
      combined_data <- c(combined_data, marker_mean_intensity_scaled)
      
      min_val <- min(marker_mean_intensity_scaled, na.rm = TRUE)
      max_val <- max(marker_mean_intensity_scaled, na.rm = TRUE)
      min_list <- c(min_list, min_val)
      max_list <- c(max_list, max_val)
      
      global_min <- min(global_min, min_val)
      global_max <- max(global_max, max_val)
    }
    
    results_range[[marker_name]] <- list(
      min_list = min_list,
      max_list = max_list,
      global_min = global_min,
      global_max = global_max
    )
    
    # Second pass: create histograms and GMM analysis
    hist_list <- list()
    bin_edge_list <- list()
    gmm_curves[[marker_name]] <- list()
    
    plot_data <- data.frame()
    
    for (i in seq_along(sample_names)) {
      sample_name <- sample_names[i]
      sample_mask <- adata_list$obs$sample_id == sample_name
      marker_data <- data_matrix[sample_mask, marker_index]
      marker_mean_intensity_scaled <- preprocess_raw(marker_data)
      
      # Create histogram
      hist_result <- hist(marker_mean_intensity_scaled, 
                         breaks = bin_counts,
                         range = c(global_min, global_max),
                         plot = FALSE)
      
      hist_list[[i]] <- hist_result$counts
      bin_edge_list[[i]] <- hist_result$breaks
      
      # GMM analysis
      local_gmm <- plot_local_gmm(marker_mean_intensity_scaled, bin_counts, global_min, global_max)
      gmm_curves[[marker_name]][[sample_name]] <- list(
        cut_means = local_gmm$cut,
        x_axis = local_gmm$x_axis,
        y_axis0 = local_gmm$y_axis0,
        y_axis1 = local_gmm$y_axis1
      )
      
      # Prepare plot data
      temp_data <- data.frame(
        bin_centers = hist_result$mids,
        counts = hist_result$counts,
        sample = sample_name
      )
      plot_data <- rbind(plot_data, temp_data)
    }
    
    results_hist[[marker_name]] <- list(
      hist_list = hist_list,
      bin_edge_list = bin_edge_list
    )
    
    # Create ggplot
    p <- ggplot(plot_data, aes(x = bin_centers, y = counts, color = sample)) +
      geom_line(linewidth = 1, alpha = 0.7) +
      labs(title = marker_name,
           x = "Log10 Cell Mean Intensity",
           y = "Frequency") +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 16),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 2)
      )
    
    if (!is.null(xlims) && marker_name %in% names(xlims)) {
      p <- p + xlim(xlims[[marker_name]])
    }
    if (!is.null(ylims) && marker_name %in% names(ylims)) {
      p <- p + ylim(ylims[[marker_name]])
    }
    
    # Add global threshold
    if (plot_global_threshold) {
      global_threshold <- plot_global_gmm(combined_data)
      p <- p + geom_vline(xintercept = global_threshold, linetype = "dashed", color = "black")
    }
    
    plot_list[[marker_name]] <- p
  }
  
  return(list(
    results_range = results_range,
    results_hist = results_hist,
    gmm_curves = gmm_curves,
    plots = plot_list
  ))
}