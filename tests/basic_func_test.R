#' Test basic UniFORM functionality

library(testthat)
library(UniFORM)

test_that("Data loading works correctly", {
  # Create sample test data
  test_data <- data.frame(
    "DAPI (DAPI) Nucleus Intensity" = rnorm(100, 1000, 200),
    "CD3 (Opal 650) Cytoplasm Intensity" = rnorm(100, 500, 100),
    "Classifier Label" = sample(c("Tumor", "Stroma"), 100, replace = TRUE),
    stringsAsFactors = FALSE
  )
  
  # Write test CSV
  temp_file <- tempfile(fileext = ".csv")
  write.csv(test_data, temp_file, row.names = FALSE)
  
  # Load data
  adata <- load_csv_to_adata(
    csv_files = temp_file,
    sample_names = "test_sample", 
    marker_columns = c("DAPI (DAPI) Nucleus Intensity", "CD3 (Opal 650) Cytoplasm Intensity")
  )
  
  expect_true(is.list(adata))
  expect_true("X" %in% names(adata))
  expect_true("obs" %in% names(adata))
  expect_true("var" %in% names(adata))
  expect_equal(nrow(adata$X), 100)
  expect_equal(ncol(adata$X), 2)
  
  # Clean up
  unlink(temp_file)
})

test_that("UniFORM model initialization works", {
  # Create minimal test data
  adata <- list(
    X = matrix(rnorm(200), nrow = 100, ncol = 2),
    obs = data.frame(sample_id = rep(c("sample1", "sample2"), each = 50)),
    var = data.frame(marker_name = c("marker1", "marker2"))
  )
  
  model <- UniFORM$new(adata)
  expect_true(is(model, "UniFORM"))
  expect_equal(model$adata, adata)
})