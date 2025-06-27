#' Data Loading Functions
#' 
#' Functions to load and format data for UniFORM

#' Load CSV data and convert to AnnData-like format
#' @param csv_files Character vector of CSV file paths
#' @param sample_names Character vector of sample names
#' @param marker_columns Character vector of marker column names
#' @return List with X, obs, var components similar to Python AnnData

load_csv_to_adata <- function(csv_files, sample_names, marker_columns) {
  
  all_data <- list()
  obs_data <- data.frame()
  
  for (i in seq_along(csv_files)) {
    cat("Loading", csv_files[i], "\n")
    
    # Read CSV file
    data <- read.csv(csv_files[i], stringsAsFactors = FALSE)
    
    # Extract marker intensity columns
    if (missing(marker_columns)) {
      # Auto-detect intensity columns (assuming they contain "Intensity" or similar patterns)
      intensity_cols <- grep("Intensity|Mean|Positive", names(data), value = TRUE)
      # Filter to get the main intensity measurements
      marker_cols <- intensity_cols[!grepl("Nucleus|Cytoplasm", intensity_cols)]
    } else {
      marker_cols <- marker_columns
    }
    
    # Extract intensity data
    intensity_data <- data[, marker_cols, drop = FALSE]
    
    # Create observation metadata
    n_cells <- nrow(intensity_data)
    sample_obs <- data.frame(
      sample_id = rep(sample_names[i], n_cells),
      cell_id = paste0(sample_names[i], "_", seq_len(n_cells)),
      stringsAsFactors = FALSE
    )
    
    # Add additional metadata if available
    if ("Classifier Label" %in% names(data)) {
      sample_obs$cell_type <- data[["Classifier Label"]]
    }
    if ("Cell Area (µm²)" %in% names(data)) {
      sample_obs$cell_area <- data[["Cell Area (µm²)"]]
    }
    
    all_data[[i]] <- intensity_data
    obs_data <- rbind(obs_data, sample_obs)
  }
  
  # Combine all intensity data
  combined_data <- do.call(rbind, all_data)
  
  # Create variable metadata
  var_data <- data.frame(
    marker_name = names(combined_data),
    stringsAsFactors = FALSE
  )
  
  # Create AnnData-like structure
  adata_list <- list(
    X = as.matrix(combined_data),
    obs = obs_data,
    var = var_data,
    layers = list()
  )
  
  return(adata_list)
}

#' Load example data from the provided CSV file
load_example_data <- function(data_dir) {
  # Get all CSV files in the data directory
  csv_files <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)
  
  # Extract sample names from filenames
  sample_names <- gsub(".*?([^/]+)\\.csv$", "\\1", csv_files)
  sample_names <- gsub("_object_Data", "", sample_names)
  sample_names <- gsub(".*_", "", sample_names)
  
  # Define marker columns based on the provided data structure
  marker_columns <- c(
    "DAPI (DAPI) Nucleus Intensity",
    "CD3 (Opal 650) Cytoplasm Intensity", 
    "CD4 (Opal 570) Cytoplasm Intensity",
    "CD8 (Opal 520) Cytoplasm Intensity", 
    "CD19 (Opal 540) Cytoplasm Intensity",
    "CD138 (Opal 620) Cytoplasm Intensity"
  )
  
  return(load_csv_to_adata(csv_files, sample_names, marker_columns))
}