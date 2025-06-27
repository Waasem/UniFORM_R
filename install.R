#!/usr/bin/env Rscript

# UniFORM-R Installation Script
# Run this script to install UniFORM-R and all dependencies

cat("=== UniFORM-R Installation Script ===\n\n")

# Check R version
r_version <- R.Version()
cat("R version:", r_version$version.string, "\n")

if (as.numeric(paste(r_version$major, r_version$minor, sep = ".")) < 4.0) {
  stop("UniFORM-R requires R version 4.0.0 or higher")
}

# Install required packages
required_packages <- c(
  "R6",
  "ggplot2", 
  "dplyr",
  "tidyr",
  "mixtools",
  "signal",
  "testthat",
  "devtools"
)

cat("\nChecking and installing required packages...\n")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, dependencies = TRUE)
  } else {
    cat(pkg, "is already installed\n")
  }
}

# Install UniFORM-R package
cat("\nInstalling UniFORM-R package...\n")

if (require("devtools", quietly = TRUE)) {
  devtools::install(".", dependencies = TRUE)
} else {
  stop("devtools package is required for installation")
}

# Run tests
cat("\nRunning package tests...\n")
if (require("testthat", quietly = TRUE)) {
  testthat::test_dir("tests/")
} else {
  cat("Warning: testthat not available, skipping tests\n")
}

cat("\n=== Installation Complete! ===\n")
cat("Load the package with: library(UniFORM)\n")
cat("See examples/basic_usage.R for usage examples\n")