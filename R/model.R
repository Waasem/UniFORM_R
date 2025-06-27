#' UniFORM Model Class
#' 
#' R6 class implementation of the UniFORM model

library(R6)

UniFORM <- R6Class("UniFORM",
  public = list(
    adata = NULL,
    
    initialize = function(adata_list) {
      self$adata <- adata_list
    },
    
    uniform_calculate_histogram = function(marker_to_plot, ...) {
      uniform_calculate_histogram(self$adata, marker_to_plot, ...)
    },
    
    normalize_and_plot_distributions = function(histograms, markers, reference_sample, ...) {
      normalize_and_plot_distributions(self$adata, histograms, markers, reference_sample, ...)
    },
    
    calculate_shift_in_log_pixels = function(data_range, shifts_fft_dict, ...) {
      calculate_shift_in_log_pixels(data_range, shifts_fft_dict, ...)
    },
    
    generate_normalized_feature = function(shift_in_log_pixels_dict, reference_sample, output_directory, ...) {
      generate_normalized_feature(self$adata, shift_in_log_pixels_dict, reference_sample, output_directory, ...)
    }
  )
)