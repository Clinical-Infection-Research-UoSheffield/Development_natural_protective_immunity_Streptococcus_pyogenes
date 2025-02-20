# ---------------------------------------------------------
# Script Title: Data Download and Extraction
# Author: Alex Keeley   
# Date: 20/02/2025
# 
# Description:
# This script automates the download and extraction of data to reproduce analyses 
# that will be publicly available under a Creative Commons (CC) license 
# upon publication. The dataset is hosted on Zenodo.
#
# The data will be distributed under the Creative Commons Attribution 
# 4.0 International (CC BY-NC-ND 4.0) license.
# ---------------------------------------------------------


# Load required packages

library(librarian)
shelf(httr)

# Define URL and file paths
zip_url <- "<insert once public>"
zip_file <- "data.zip"
data_dir <- "data"

# Ensure the 'data/' directory exists
if (!dir.exists(data_dir)) {
    dir.create(data_dir)
}

# Download ZIP file if it doesn't exist
if (!file.exists(zip_file)) {
    message("Downloading dataset from Zenodo...")
    download.file(zip_url, zip_file, mode = "wb")
} else {
    message("ZIP file already exists, skipping download.")
}

# Extract the ZIP file into 'data/' directory
message("Extracting dataset...")
unzip(zip_file, exdir = data_dir)

message("Dataset successfully extracted to 'data/' directory.")