# Import baseline dataframes

# Title: Analysis of M peptide responses
# Version: 1.1
# Date: 2024-09-25
# Author: Dr. Alexander J. Keeley
# Inputs: 
# Outputs: 

# Description:

# This script performs a baseline analysis of IgG titres in participants without Strep A disease events. The analysis focuses on:
# 1. Extracting baseline titres from participants with no Strep A positive disease events.
# 2. Fitting fractional polynomial models to predict IgG titres based on age and generating residual and fitted values.
# 3. Calculating centile values and adding prediction bands based on these centiles.
# 4. Generating plots for each antigen showing the distribution of titres and prediction bands.
# 5. Import raw MFI values from the same samples, inorder to visualise relative 


# Requirements:

#            Package Version
#CorrMixed CorrMixed     1.1
#gridExtra gridExtra     2.3
#mfp             mfp 1.5.4.1
#tidyverse tidyverse   2.0.0

# Setup environment
source("scripts/setup_environment.R")
source("scripts/load_functions.R")

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
                size = base_size * 1.2
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



plot_basesize = 16



# Load required packages using librarian
shelf(tidyverse, mfp, CorrMixed, gridExtra, lme4, gtsummary) 
shelf(FSA)

#1. make a wide dataframe with age

df_titre <- readRDS("data/M_CRC_baseline_no_disease.RDS")

df.text1 <- df_titre %>%
    select(pid, Antigen, titre) 


df <- readRDS("data/baseline_blood_no_disease_titres.RDS")
df.text2 <- df %>%
    select(pid, Antigen, titre)

df.text <- rbind(df.text1, df.text2) %>%
    spread(key = Antigen, value = titre) %>%
    ungroup() %>%
    select(-pid)

age <- readRDS("data/final_dems.RDS") %>%
    mutate(
        age_grp = case_when(
            age < 2 ~ "< 2 years",
            age < 5 & age >= 2 ~ "2-4 years",
            age < 12 & age >= 5 ~ "5-11 years",
            age < 19 & age >= 12 ~ "12-18 years",
            TRUE ~ age_grp
        ),
        age_grp = factor(age_grp, levels = c("< 2 years", "2-4 years", "5-11 years","12-18 years", "Over 18 years"))
    ) %>%
    filter(!is.na(age)) %>%
    select(pid, age_grp)


wide_age_df <- bind_rows(df.text1,df.text2) %>%
    left_join(age) %>%
    spread(Antigen,titre)

simple_df <- wide_age_df %>%
    ungroup() %>%
    select(-pid)



# Desired fixed order for non-M antigens
fixed_order <- c("age_grp", "GAC", "SLO", "SpyAD", "SpyCEP", "DNAseB")

# Extract M-type columns dynamically
m_columns <- setdiff(colnames(df.text), fixed_order)
m_columns <- sort(m_columns)  # Optional: sort the M-types alphabetically

# Reorder dataframe columns
simple_df <- simple_df[, c(fixed_order, m_columns)]




# Initialize an empty list to store the correlation matrices with p-values
correlation_results <- list()

# Iterate over each age group
for (ag in unique(simple_df$age_grp)) {
    
    print(ag)
    
    # Filter the dataframe for the current age group
    for_df <- simple_df %>% 
        filter(age_grp == ag) %>%
        ungroup() %>%
        select(-age_grp)
    
    # Compute the correlation matrix and p-values for the current age group
    corr_result <- Hmisc::rcorr(as.matrix(for_df), type = "pearson")
    corr_matrix <- corr_result$r  # Correlation coefficients
    p_matrix <- corr_result$P     # P-values
    
    # Convert the correlation matrix into a long format
    corr_long <- as.data.frame(as.table(corr_matrix))
    names(corr_long) <- c("xName", "yName", "corr")  # Rename columns
    
    # Convert the p-value matrix into a long format
    p_long <- as.data.frame(as.table(p_matrix))
    names(p_long) <- c("xName", "yName", "p_value")  # Rename columns
    
    # Merge the correlation and p-value dataframes
    corr_long <- merge(corr_long, p_long, by = c("xName", "yName"))
    
    # Add the current age group as a new column
    corr_long$age_grp <- ag
    
    # Append the result to the list
    correlation_results[[ag]] <- corr_long
}

# Combine all the correlation matrices into one data frame
final_correlation_df <- do.call(rbind, correlation_results)


final_correlation_df$age_grp %>% as.factor() %>% levels()
# set levels for age groups 
level_order <- c("< 2 years", "2-4 years", "5-11 years","12-18 years", "Over 18 years")

final_correlation_df$age_grp <- factor(final_correlation_df$age_grp,levels = level_order)


final_correlation_df$age_grp %>% as.factor() %>% levels()


pal <- wesanderson::wes_palette("Zissou1", 100, type = "continuous")

corr_plot <- final_correlation_df %>%
    filter(!is.na(age_grp), xName != yName) %>%
    ggplot(
        aes(
            x = age_grp, y = 1, fill = corr
        )
    ) +
    geom_tile() +
    facet_grid(xName ~ yName) +
    scale_fill_gradientn(colors = pal) +
    theme_minimal() +
    labs(
        x = "Age group",
        y = NULL,
        fill = "Spearman's rho"
    ) +
    theme(
        axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        strip.text.y = element_text(angle = 0))


corr_plot

ggsave("R_output/supp_r2r_M_conserved_age_corr.png",plot = corr_plot, width = 1400 / 96, height = 700/ 96, dpi = 300, bg = "white")


final_correlation_df %>%
    filter(!is.na(age_grp), xName != yName,
           age_grp == "Over 18 years") %>%
    ggplot(
        aes(
            x = 1, y = 1, fill = corr
        )
    ) +
    geom_tile() +
    facet_grid(xName ~ yName) +
    scale_fill_gradientn(colors = pal) +
    theme_minimal() +
    labs(
        x = NULL,
        y = NULL,
        fill = "Spearman's rho"
    ) +
    theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        strip.text.y = element_text(angle = 0))

