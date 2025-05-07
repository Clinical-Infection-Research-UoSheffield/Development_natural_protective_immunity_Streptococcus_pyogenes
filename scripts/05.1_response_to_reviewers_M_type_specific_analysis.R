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


#########



df_titre <- readRDS("R_objects/M_CRC_baseline_no_disease.RDS")

df.text1 <- df_titre %>%
    select(pid, Antigen, titre) 


df <- readRDS("R_objects/baseline_blood_no_disease_titres.RDS")
df.text2 <- df %>%
    select(pid, Antigen, titre)

df.text <- rbind(df.text1, df.text2) %>%
    spread(key = Antigen, value = titre) %>%
    ungroup() %>%
    select(-pid)


##### ##### ##### ##### ##### 
##### correlation matrix #### 
##### ##### ##### ##### ##### 

##### Tidy data to allow a correlation matrix to be draw for each individual at baseline

colnames(df.text)

# Desired fixed order for non-M antigens
fixed_order <- c("GAC", "SLO", "SpyAD", "SpyCEP", "DNAseB")

# Extract M-type columns dynamically
m_columns <- setdiff(colnames(df.text), fixed_order)
m_columns <- sort(m_columns)  # Optional: sort the M-types alphabetically

# Reorder dataframe columns
df.text <- df.text[, c(fixed_order, m_columns)]


# Calculate the correlation matrix
cor_matrix <- cor(df.text, use = "pairwise.complete.obs")

# Turn the correlation matrix into a long dataframe
cor_data_long <- as.data.frame(as.table(cor_matrix))

# Remove vales where both antigens are the same
cor_data_filtered <- cor_data_long %>%
    filter(Var1 != Var2)

# Calculate the range of correlation coefficients for each group
range_df <- cor_data_filtered %>%
    summarise(min_coef = min(Freq, na.rm = TRUE), max_coef = max(Freq, na.rm = TRUE))

# Construct the sentence
sentence <- sprintf(
    "When comparing baseline anti-M antibody titres within individuals, a highly variable correlation coefficient was observed, ranging from %.2f to %.2f",
    range_df$min_coef, range_df$max_coef
)

# Print the sentence
print(sentence)





# Function to calculate correlation coefficients and p-values
get_cor_and_pval <- function(df) {
    # Initialize an empty dataframe to store results
    results <- tibble(Var1 = character(), Var2 = character(),
                      Correlation = numeric(), PValue = numeric())
    
    # Loop through all combinations of variables
    for (i in 1:(ncol(df) - 1)) {
        for (j in (i + 1):ncol(df)) {
            # Perform correlation test
            test <- cor.test(df[[i]], df[[j]], use = "pairwise.complete.obs")
            
            # Store results in the dataframe
            results <- results %>%
                add_row(Var1 = colnames(df)[i],
                        Var2 = colnames(df)[j],
                        Correlation = test$estimate,
                        PValue = test$p.value)
        }
    }
    return(results)
}

# Apply the function to your dataframe
cor_pval_results <- get_cor_and_pval(df.text)

# View the results
print(cor_pval_results)

# Calculate the range of correlation coefficients for each group
maxp <- cor_pval_results %>%
    summarise(max_p = max(PValue, na.rm = TRUE))

# Construct the sentence with the dynamic condition
if (maxp$max_p < 0.0001) {
    sentence <- sprintf(
        "When comparing baseline anti-M antibody titres within individuals, a highly variable correlation coefficient was observed, ranging from %.2f to %.2f (p < 0.0001 for all comparisons).",
        min(cor_pval_results$Correlation, na.rm = TRUE),
        max(cor_pval_results$Correlation, na.rm = TRUE)
    )
} else {
    sentence <- sprintf(
        "When comparing baseline anti-M antibody titres within individuals, a highly variable correlation coefficient was observed, ranging from %.2f to %.2f (p < %.4f for all comparisons).",
        min(cor_pval_results$Correlation, na.rm = TRUE),
        max(cor_pval_results$Correlation, na.rm = TRUE),
        maxp$max_p
    )
}

# Print the sentence
print(sentence)

# Plot the correlation matrix for interactions
corrplot::corrplot.mixed(cor_matrix, upper = "color", lower = "number",
                         tl.col = "black", tl.srt = 45, lower.col = "black", number.cex = 1.3,
                         tl.pos = "lt",tl.cex = 1.5)  # "lt" for left and top labels




############ Now add age to this: 

#1. make a wide dataframe with age

df_titre <- readRDS("R_objects/M_CRC_baseline_no_disease.RDS")

df.text1 <- df_titre %>%
    select(pid, Antigen, titre) 


df <- readRDS("R_objects/baseline_blood_no_disease_titres.RDS")
df.text2 <- df %>%
    select(pid, Antigen, titre)

age <- readRDS("R_objects/final_dems.RDS") %>%
    mutate(
    age_grp = case_when(
        age < 2 ~ "< 2 years",
        age < 5 & age >= 2 ~ "2-4 years",
        age < 12 & age >= 5 ~ "5-11 years",
        age < 19 & age >= 12 ~ "12-18 years",
        TRUE ~ age_grp
    ),
    age_grp = factor(age_grp, levels = c("< 2 years", "2-4 years", "5-11 years","12-18 years", "Over 18 years"))
)   %>%
    mutate(age_cat = 
               case_when(
                   age >= 6  & age < 10 ~ 6,
                   age >= 10 & age < 15 ~ 7,
                   age >= 15 & age < 20 ~ 8,
                   age >= 20 & age <= 30 ~ 9,
                   age >= 30 & age <= 40 ~ 10,
                   age >= 40 ~ 11,
                   T ~ age
               )) %>%
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


# Initialize an empty list to store the correlation matrices
correlation_results <- list()

# Iterate over each age group
for (ag in unique(simple_df$age_grp)) {
    
    print(ag)
    
    # Filter the dataframe for the current age group
    for_df <- simple_df %>% 
        filter(age_grp == ag) %>%
        ungroup() %>%
        select(-age_grp)
    
    # Compute the correlation matrix for the current age group
    corr_matrix_interaction_for <- for_df %>%
        cor(use = "pairwise.complete.obs", method = "spearman")
    
    # Convert the correlation matrix into a long format
    corr_long <- as.data.frame(as.table(corr_matrix_interaction_for))
    names(corr_long) <- c("xName", "yName", "corr")  # Rename columns
    
    # Add the current age group as a new column
    corr_long$age_grp <- ag
    
    # Append the result to the list
    correlation_results[[ag]] <- corr_long
}

# Combine all the correlation matrices into one data frame
final_correlation_df <- do.call(rbind, correlation_results)

levels(as.factor(final_correlation_df$age_grp))
# set levels for age groups 
level_order <- c("< 2 years", "2-4 years", "5-11 years","12-17 years", "Over 18 years")

final_correlation_df$age_grp <- factor(final_correlation_df$age_grp,levels = level_order)



pal <- wesanderson::wes_palette("Zissou1", 100, type = "continuous")

final_correlation_df %>%
    filter(! is.na(age_grp),
           xName != yName) %>%
    ggplot(
        aes(
            x = age_grp, y = corr, col = corr
        )
    ) +
    geom_point() + 
    facet_grid(xName ~ yName) +
    scale_color_gradientn(colors = pal) +
    theme_minimal() + 
    labs(
        title = "Correlation coefficients within individuals at baseline (n = 431)",
        x = "Age group",
        y = "pearsons correlation coefficient"
    ) +
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1))

final_correlation_df %>%
    filter(! is.na(age_grp),
           xName != yName) %>%
    ggplot(
        aes(
            x = age_grp, y ="1", fill = corr
        )
    ) +
    geom_tile() + 
    facet_grid(xName ~ yName) +
    scale_fill_gradientn(colors = pal) +
    theme_minimal() + 
    labs(
        title = "Correlation coefficients within individuals at baseline (n = 431)",
        x = "Age group",
    ) +
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1))

final_correlation_df %>%
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
        title = "Correlation coefficients within individuals at baseline (n = 431)",
        x = "Age group",
        y = NULL
    ) +
    theme(
        axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()
    )


final_correlation_df %>% colnames()

#### now get p values: 

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
        select(GAC:DNAseB)
    
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

# set levels for age groups 
level_order <- c("< 2 years", "2-4 years", "5-11 years","12-18 years", "Over 18 years")

final_correlation_df$age_grp <- factor(final_correlation_df$age_grp,levels = level_order)


final_correlation_df %>%
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
        title = "Correlation coefficients within individuals at baseline (n = 431)",
        x = "Age group",
        y = NULL
    ) +
    theme(
        axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()
    )




#1. make a wide dataframe with age

df_titre <- readRDS("R_objects/M_CRC_baseline_no_disease.RDS")

df.text1 <- df_titre %>%
    select(pid, Antigen, titre) 


df <- readRDS("R_objects/baseline_blood_no_disease_titres.RDS")
df.text2 <- df %>%
    select(pid, Antigen, titre)

age <- readRDS("R_objects/final_dems.RDS") %>%
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

