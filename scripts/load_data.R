# Script Title: Data Download and Extraction
# Author: Alexander J. Keeley   
# Date: 20/02/2025 (Updated 15/05/2025)
# 
# Description:
# This script automates the download and extraction of data to reproduce analyses 
# that will be publicly available under a Creative Commons (CC) license 
# upon publication. The dataset is hosted on Zenodo.
#
# Package Version
# devtools     devtools   2.4.5
# inborutils inborutils   0.4.0
# librarian   librarian   1.8.1


# The data will be distributed under the Creative Commons Attribution 
# 4.0 International (CC BY-NC-ND 4.0) license.
# ---------------------------------------------------------

# Install librarian if not already installed
if (!requireNamespace("librarian", quietly = TRUE)) install.packages("librarian")

library(librarian)


# Install and load required packages
if (!requireNamespace("inborutils", quietly = TRUE)) {
    # load devtools 
    shelf(devtools)
    # install.packages("devtools")
    devtools::install_github("inbo/inborutils")
}

shelf(inborutils)

# Define the Zenodo DOI and target directory
doi <- "10.5281/zenodo.14887949"
data_dir <- "data"

# Create the data directory if it doesn't exist
if (!dir.exists(data_dir)) {
    dir.create(data_dir)
}

# Download all files from the Zenodo record to the data directory
message("Downloading dataset using inborutils::download_zenodo()...")
download_zenodo(doi = doi, path = data_dir, quiet = TRUE)

# List files for verification
message("Downloaded files:")
print(list.files(data_dir, recursive = TRUE))



# Detach all non-base packages
detach_all_packages <- function() {
    base_pkgs <- c("stats", "graphics", "grDevices", "utils", "datasets", "methods", "base")
    loaded_pkgs <- setdiff(loadedNamespaces(), base_pkgs)
    for (pkg in loaded_pkgs) {
        try(detach(paste0("package:", pkg), character.only = TRUE, unload = TRUE), silent = TRUE)
    }
}

detach_all_packages()
