# Title: Setup Environments for Analysis Scripts
# Version: 1.0
# Date: 2024-07-23
# Author: Dr. Alexander J. Keeley
# Inputs:  sets global options and loads necessary libraries
# Outputs: age.df, read.titres function. 

# Description:
# This script sets up the environment for analysis scripts


# Install required packages, including librarian if not already installed
if (!requireNamespace("librarian", quietly = TRUE)) install.packages("librarian")

# Create output directory if it doesn't exist 
if (!dir.exists("R_output")) {
    dir.create("R_output")
}

if (!dir.exists("R_output/final_submission")) {
    dir.create("R_output/final_submission")
}
# Load librarian
library(librarian)

shelf(tidyverse)
select <- dplyr::select

# Set global options for R

options(
    digits = 3,  # Limit the number of significant digits printed in output to 3
    scipen = 5   # Bias towards fixed notation rather than scientific notation
)

# Define a function to convert log10 values to their exponential form
log10_to_exp <- function(x) {
    10^x
}

# set levels for age groups 
level_order <- c("< 2 years", "2-4 years", "5-11 years","12-18 years", "Over 18 years")


final_dems <- readRDS("data/final_dems.RDS") 
age <- final_dems



#final_dems <- readRDS("R_objects/final_dems.RDS") %>%
#    mutate(
#        age_grp = case_when(
#            age < 2 ~ "< 2 years",
#            age < 5 & age >= 2 ~ "2-4 years",
#            age < 12 & age >= 5 ~ "5-11 years",
#            age < 19 & age >= 12 ~ "12-18 years",
#            TRUE ~ age_grp
#        ),
#        age_grp = factor(age_grp, levels = c("< 2 years", "2-4 years", "5-11 years","12-18 years", "Over 18 years"))
#    )   %>%
#    mutate(age_cat = 
#               case_when(
#                   age >= 6  & age < 10 ~ 6,
#                   age >= 10 & age < 15 ~ 7,
#                   age >= 15 & age < 20 ~ 8,
#                   age >= 20 & age <= 30 ~ 9,
#                   age >= 30 & age <= 40 ~ 10,
#                   age >= 40 ~ 11,
#                   T ~ age
#               ))
#
#### Would have to remove data of birth from the main dataframe and the script  
#
#age <- final_dems %>% 
#    select(pid, hid,
#           age,
#           age_grp,
#           sex,
#           hhsize = median_hhsize,
#           # dob,
#           age_cat)
#

read.titres <- function(titre.df, var = "titre")    {
    
    fun_titre <- {
        if (var == "log_RLU_titres") {
            readRDS(titre.df) %>%
                rename(titre = log_RLU_titres,
                       visit_date = Date) %>%
                select(-log_AUC, -RLU, -AUC)
            
        } 
        else if (var == "AUC") {
            readRDS(titre.df) %>%
                rename(titre = log_AUC,
                       visit_date = Date) %>%
                select(-log_RLU_titres, -RLU, -AUC)
        } 
        else {
            readRDS(titre.df) 
        }
    }
    
    return(fun_titre)
}

############
############
############

library(ggplot2)

# Universal theme function
theme_universal <- function(base_size = 8, base_family = "") {
    theme_minimal(base_size = base_size, base_family = base_family) +
        theme(
            # Title settings
            plot.title = element_text(
                size = base_size,
                hjust = 0.5,
                face = "bold",
                margin = margin(b = base_size / 2)
            ),
            # Subtitle settings
            plot.subtitle = element_text(
                size = base_size * 1.2,
                hjust = 0.5,
                margin = margin(b = base_size / 2)
            ),
            # Facet title settings
            strip.text = element_text(
                size = base_size,
                face = "bold",
                hjust = 0.5
            ),
            # Axis title settings
            axis.title.x = element_text(
                size = base_size,
                face = "bold",
                margin = margin(t = base_size / 2)  # Add space above the x-axis title
            ),
            axis.title.y = element_text(
                size = base_size,
                face = "bold" 
            ),
            # Axis text settings
            axis.text = element_text(
                size = base_size * 0.8
            ),
            # Legend title and text
            legend.title = element_text(
                size = base_size,
                face = "bold"
            ),
            legend.text = element_text(
                size = base_size
            ),
            # Margins
            plot.margin = margin(t = base_size, r = base_size, b = base_size, l = base_size)
        )
}


plot_basesize = 14


############
############

StrepA_colscheme <- c(
    "SpyCEP" = "#FDC086",  # Original colors
    "SpyAD" = "#d19c2f",
    "SLO" = "#386CB0",
    "GAC" = "#7FC97F",
    "DNAseB" = "#BEAED4",
    # E3 m peptides distinct shades of blue for M103, M113, M25, M44, M82, M87
    "M103" = "#1F78B4",  # Bright blue
    "M113" = "#6BAED6",  # Sky blue
    "M25"  = "#08519C",  # Dark blue
    "M44"  = "#4292C6",  # Medium blue
    "M82"  = "#08306B",  # Navy blue
    "M87"  = "#3182BD",  # Steel blue
    # Conserves M for P17, J8, K4S2
    "P17"  = "#D73027",  # Deep red
    "J8"   = "#FC8D59",  # Lighter red-orange
    "K4S2" = "#B2182B",  # Dark crimson red
    # Distinct colors for other Antigens
    "M1"   = "#66C2A5",
    "M18"  = "#FC8D62",
    "M3"   = "#8DA0CB",
    "M4"   = "#E78AC3",
    "M53"  = "#A6D854",
    "M55"  = "#FFD92F",
    "M6"   = "#E5C494",
    "M71"  = "#B3B3B3",
    "M74"  = "#1F78B4",
    "M75"  = "#33A02C",
    "M76"  = "#FB9A99",
    "M89"  = "#E31A1C",
    "M97"  = "#FDBF6F"
)

event_colours <-
    c("skin carriage" = "#3B9AB2",
    "skin disease" = "#78B7C5",
    "throat carraige" = "#EBCC2A",
    "throat disease" = "#E1AF00",
    "other" = "#F21A00") 

