#' Basic Usage Example for UniFORM-R
#' 
#' This example demonstrates how to use the UniFORM normalization pipeline in R

library(UniFORM)

# Load your data
data_dir <- "path/to/your/data"  # Update this path
adata <- load_example_data(data_dir)

# Initialize UniFORM model
uniform_model <- UniFORM$new(adata)

# Step 1: Calculate histograms and perform exploratory data analysis
markers_to_analyze <- c("CD3 (Opal 650) Cytoplasm Intensity", 
                       "CD4 (Opal 570) Cytoplasm Intensity",
                       "CD8 (Opal 520) Cytoplasm Intensity")

results <- uniform_model$uniform_calculate_histogram(
  marker_to_plot = markers_to_analyze,
  bin_counts = 1024,
  plot_local_threshold = TRUE,
  plot_global_threshold = TRUE
)

# Extract results
results_range <- results$results_range
results_hist <- results$results_hist
gmm_curves <- results$gmm_curves

# Step 2: Define reference samples for each marker
# Update these based on your analysis from Step 1
reference_samples <- list(
  "CD3 (Opal 650) Cytoplasm Intensity" = "sample1",
  "CD4 (Opal 570) Cytoplasm Intensity" = "sample1", 
  "CD8 (Opal 520) Cytoplasm Intensity" = "sample2"
)

# Step 3: Perform normalization and calculate shifts
shifts_fft_dict <- uniform_model$normalize_and_plot_distributions(
  histograms = results_hist,
  markers = names(reference_samples),
  reference_sample = reference_samples
)

# Step 4: Calculate shifts in log pixel scale
shift_in_log_pixels_dict <- uniform_model$calculate_shift_in_log_pixels(
  data_range = results_range,
  shifts_fft_dict = shifts_fft_dict
)

# Step 5: Generate normalized features
output_directory <- "path/to/output"  # Update this path
normalized_adata <- uniform_model$generate_normalized_feature(
  shift_in_log_pixels_dict = shift_in_log_pixels_dict,
  reference_sample = reference_samples,
  output_directory = output_directory,
  plot_dist = TRUE,
  save_normalized_features = TRUE
)

cat("Normalization complete! Check", output_directory, "for results.\n")