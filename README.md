# UniFORM_R: Universal ImmunoFluorescence nORMalization

[![R](https://img.shields.io/badge/R-%3E%3D4.0.0-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

UniFORM_R is an R implementation of the UniFORM (Universal ImmunoFluorescence nORMalization) pipeline for normalizing multiplex immunofluorescence data. This package provides robust normalization methods to harmonize intensity measurements across different samples, batches, and experimental conditions.

## Overview

UniFORM addresses the critical need for standardization in multiplex immunofluorescence imaging by:
- **Histogram-based normalization**: Uses cross-correlation to align intensity distributions
- **Gaussian Mixture Model (GMM) analysis**: Identifies positive and negative cell populations
- **Feature-level normalization**: Works directly with extracted cell features
- **Batch effect correction**: Harmonizes data across different experimental batches

## Installation

### Prerequisites

Ensure you have R version 4.0.0 or higher installed on your system.

### Install Required Packages

```r
# Install required CRAN packages
install.packages(c(
  "R6",           # Object-oriented programming
  "ggplot2",      # Plotting and visualization
  "dplyr",        # Data manipulation
  "tidyr",        # Data tidying
  "mixtools",     # Gaussian mixture models
  "signal",       # Signal processing for cross-correlation
  "testthat"      # Testing framework
))
```

### Install UniFORM_R

Since this is a local package, you can install it using:

```r
# Install from local directory
devtools::install("path/to/UniFORM_R")

# Or if you're in the UniFORM_R directory
devtools::install(".")
```

## Quick Start

### 1. Load the Package

```r
library(UniFORM)
```

### 2. Prepare Your Data

UniFORM_R expects CSV files with cell-level intensity measurements. Your data should contain:
- Cell intensity measurements for each marker
- Sample identifiers
- Optional metadata (cell type, area, etc.)

Example data format:
```
Image Location, Object Id, CD3 (Opal 650) Cytoplasm Intensity, CD4 (Opal 570) Cytoplasm Intensity, ...
sample1.tif, 1, 1.23, 0.45, ...
sample1.tif, 2, 2.10, 0.89, ...
```

### 3. Load Your Data

```r
# For the provided example data
data_dir <- "/UniFORM/data"
adata <- load_example_data(data_dir)

# For custom data
csv_files <- c("sample1_data.csv", "sample2_data.csv", "sample3_data.csv")
sample_names <- c("Sample1", "Sample2", "Sample3")
marker_columns <- c(
  "CD3 (Opal 650) Cytoplasm Intensity",
  "CD4 (Opal 570) Cytoplasm Intensity",
  "CD8 (Opal 520) Cytoplasm Intensity"
)

adata <- load_csv_to_adata(csv_files, sample_names, marker_columns)
```

### 4. Initialize UniFORM Model

```r
uniform_model <- UniFORM$new(adata)
```

### 5. Run the Normalization Pipeline

```r
# Step 1: Calculate histograms and perform exploratory analysis
markers_to_analyze <- c(
  "CD3 (Opal 650) Cytoplasm Intensity",
  "CD4 (Opal 570) Cytoplasm Intensity", 
  "CD8 (Opal 520) Cytoplasm Intensity"
)

results <- uniform_model$uniform_calculate_histogram(
  marker_to_plot = markers_to_analyze,
  bin_counts = 1024,
  plot_local_threshold = TRUE,
  plot_global_threshold = TRUE
)

# Step 2: Define reference samples for each marker
# Choose reference samples based on the histogram analysis
reference_samples <- list(
  "CD3 (Opal 650) Cytoplasm Intensity" = "Sample1",
  "CD4 (Opal 570) Cytoplasm Intensity" = "Sample1",
  "CD8 (Opal 520) Cytoplasm Intensity" = "Sample2"
)

# Step 3: Calculate normalization shifts
shifts_fft_dict <- uniform_model$normalize_and_plot_distributions(
  histograms = results$results_hist,
  markers = names(reference_samples),
  reference_sample = reference_samples
)

# Step 4: Convert shifts to log pixel scale
shift_in_log_pixels_dict <- uniform_model$calculate_shift_in_log_pixels(
  data_range = results$results_range,
  shifts_fft_dict = shifts_fft_dict
)

# Step 5: Apply normalization
output_directory <- "./normalized_output"
normalized_adata <- uniform_model$generate_normalized_feature(
  shift_in_log_pixels_dict = shift_in_log_pixels_dict,
  reference_sample = reference_samples,
  output_directory = output_directory,
  plot_dist = TRUE,
  save_normalized_features = TRUE
)
```

## Data Structure

UniFORM_R uses an AnnData-like structure with the following components:

```r
adata <- list(
  X = matrix,           # Cell x Marker intensity matrix
  obs = data.frame,     # Cell metadata (sample_id, cell_type, etc.)
  var = data.frame,     # Marker metadata (marker_name, etc.)
  layers = list()       # Additional data layers (normalized, etc.)
)
```

## Key Functions

### Core Functions

- **`uniform_calculate_histogram()`**: Calculates intensity histograms and performs GMM analysis
- **`normalize_and_plot_distributions()`**: Computes normalization shifts using cross-correlation
- **`calculate_shift_in_log_pixels()`**: Converts integer shifts to log pixel scale
- **`generate_normalized_feature()`**: Applies normalization to feature data

### Data Loading Functions

- **`load_csv_to_adata()`**: Load CSV files into UniFORM-compatible format
- **`load_example_data()`**: Load the provided example dataset

### Model Class

- **`UniFORM`**: R6 class that encapsulates the entire normalization pipeline

## Example Workflow

See the complete example in [`examples/basic_usage.R`](examples/basic_usage.R):

```r
# Load package
library(UniFORM)

# Load data
data_dir <- "path/to/your/data"
adata <- load_example_data(data_dir)

# Initialize model
uniform_model <- UniFORM$new(adata)

# Run pipeline
results <- uniform_model$uniform_calculate_histogram(markers_to_analyze)
shifts <- uniform_model$normalize_and_plot_distributions(...)
normalized_data <- uniform_model$generate_normalized_feature(...)
```

## Output

The normalization pipeline produces:

1. **Histogram plots**: Showing intensity distributions before and after normalization
2. **GMM analysis plots**: Displaying positive/negative population identification
3. **Normalized data**: Stored in `adata$layers` with timestamp
4. **Shift parameters**: Quantifying the normalization adjustments
5. **Quality control plots**: Comparing original vs. normalized distributions

## Advanced Usage

### Custom Reference Selection

```r
# You can specify different reference samples for each marker
reference_samples <- list(
  "CD3 (Opal 650) Cytoplasm Intensity" = "HighQualitySample",
  "CD4 (Opal 570) Cytoplasm Intensity" = "BestCD4Sample",
  "CD8 (Opal 520) Cytoplasm Intensity" = "OptimalCD8Sample"
)
```

### Landmark-based Fine-tuning

```r
# For advanced users: provide custom landmarks for each marker
landmarks <- list(
  "CD3 (Opal 650) Cytoplasm Intensity" = c(1.2, 1.5, 1.8),
  "CD4 (Opal 570) Cytoplasm Intensity" = c(0.8, 1.0, 1.3)
)

shifts <- uniform_model$normalize_and_plot_distributions(
  histograms = results$results_hist,
  markers = names(reference_samples),
  reference_sample = reference_samples,
  landmarks = landmarks
)
```

## Testing

Run the test suite to verify installation:

```r
# Run all tests
testthat::test_dir("tests/")

# Run specific test file
testthat::test_file("tests/basic_func_test.R")
```

## Troubleshooting

### Common Issues

1. **Package installation errors**:
   ```r
   # Try installing dependencies individually
   install.packages("mixtools")
   install.packages("signal")
   ```

2. **Data loading issues**:
   ```r
   # Check your CSV file structure
   data <- read.csv("your_file.csv")
   head(data)
   colnames(data)
   ```

3. **Memory issues with large datasets**:
   ```r
   # Process markers individually or reduce bin_counts
   results <- uniform_model$uniform_calculate_histogram(
     marker_to_plot = markers_to_analyze[1:3],  # Process subset
     bin_counts = 512  # Reduce from default 1024
   )
   ```

### Getting Help

1. Check function documentation: `?uniform_calculate_histogram`
2. Review example usage: `examples/basic_usage.R`
3. Run tests to verify installation: `testthat::test_dir("tests/")`

## Comparison with Python Version

This R implementation maintains full compatibility with the original Python UniFORM pipeline:

| Feature | Python | R |
|---------|--------|---|
| Data Structure | AnnData | List (AnnData-like) |
| GMM Analysis | scikit-learn | mixtools |
| Cross-correlation | scipy.signal | signal |
| Plotting | matplotlib | ggplot2 |
| Object-Oriented | Python classes | R6 classes |

## Contributing

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Ensure all tests pass
5. Submit a pull request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use UniFORM_R in your research, please cite:

```
Chang, M.K. et al. (2024). UniFORM: Universal ImmunoFluorescence Normalization 
for Multiplex Tissue Imaging. [Journal/Preprint].
```

## Authors

- **Mark Kunlun Wang** - Original Python implementation - wangmar@ohsu.edu
- **R Implementation Team** - R conversion and maintenance - asim.waqas@moffitt.org

## Acknowledgments

- Original UniFORM Python implementation by Chang Lab
- UniFORM_R implementation by Schabath Lab
- R package structure follows Bioconductor guidelines
- GMM analysis adapted from mixtools package
- Cross-correlation implementation uses signal package
