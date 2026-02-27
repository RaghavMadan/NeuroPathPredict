#!/usr/bin/env Rscript

# Setup script for NeuroPathPredict Pipeline V1.1

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Required R packages
required_packages <- c(
  "dplyr", "tidyr", "readr", "stringr", "ggplot2",
  "caret", "glmnet", "gstat", "sp", "raster",
  "logging", "config", "purrr", "furrr", "tictoc",
  "moments", "MASS", "yaml", "RNifti"  # Added RNifti package
)

# Required Python packages
required_python_packages <- c(
  "numpy", "pandas", "scikit-learn", "nibabel"
)

# Function to install and load R packages
install_and_load_packages <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    }
  }
}

# Function to check Python installation and pip
check_python <- function() {
  # Check Python
  python_version <- system("python3 --version", intern = TRUE)
  if (length(python_version) == 0) {
    stop("Python 3 not found. Please install Python 3.")
  }
  cat(sprintf("Found %s\n", python_version))
  
  # Check pip
  pip_version <- system("pip3 --version", intern = TRUE)
  if (length(pip_version) == 0) {
    stop("pip3 not found. Please install pip for Python 3.")
  }
  cat(sprintf("Found pip: %s\n", pip_version))
}

# Function to install Python packages
install_python_packages <- function(packages) {
  cat("\nChecking Python packages...\n")
  for (pkg in packages) {
    # Check if package is installed
    check_cmd <- sprintf("pip3 show %s > /dev/null 2>&1", pkg)
    if (system(check_cmd) != 0) {
      cat(sprintf("Installing %s...\n", pkg))
      install_cmd <- sprintf("pip3 install %s", pkg)
      if (system(install_cmd) != 0) {
        warning(sprintf("Failed to install %s", pkg))
      }
    } else {
      cat(sprintf("✓ %s\n", pkg))
    }
  }
}

# Install and load all required R packages
cat("Installing and loading required R packages...\n")
install_and_load_packages(required_packages)

# Verify R package installation
cat("\nVerifying R package installation:\n")
for (pkg in required_packages) {
  if (require(pkg, character.only = TRUE)) {
    version <- packageVersion(pkg)
    cat(sprintf("✓ %s (v%s)\n", pkg, version))
  } else {
    cat(sprintf("✗ %s\n", pkg))
  }
}

# Check Python and install required Python packages
cat("\nChecking Python installation...\n")
tryCatch({
  check_python()
  install_python_packages(required_python_packages)
  
  # Print Python package versions
  cat("\nPython package versions:\n")
  for (pkg in required_python_packages) {
    version_cmd <- sprintf("python3 -c 'import %s; print(\"%s: \" + %s.__version__)'", 
                         pkg, pkg, pkg)
    system(version_cmd)
  }
}, error = function(e) {
  cat(sprintf("\nError: %s\n", e$message))
  cat("Please ensure Python 3 and pip are installed correctly.\n")
})

cat("\nSetup completed!\n") 