# Title: Baseline Analysis of IgG Titres in Participants Without Strep A Disease Events from SpyCATS household cohort
# Version: 1.1
# Date: 2024-12-13
# Author: Dr Alexander J Keeley
# Inputs: data/baseline_blood_no_disease_titres.RDS, data/bloodIgG_protective_threshold_df.RDS
# Outputs: Plots and summary statistics describing baseline antibody titres across age groups and antigens

# Description:

# This script performs a baseline analysis of IgG titres in participants without Strep A disease events. The analysis focuses on:
# 1. Extracting baseline IgG levels from participants with no Strep A positive disease events.
# 2. Fitting fractional polynomial models to predict IgG titres based on age and generating residual and fitted values.
# 3. Calculating centile values and adding prediction bands based on these centiles.
# 4. Generating plots for each antigen showing the distribution of titres and prediction bands.
# 5. Performing correlation matrix and analysis of baseline IgG titres 
# 6. Calculate and visualise the proportion of participants who have baseline titres above the putative 50% protective thresholds 

# Requirements:

# Package Version

# CorrMixed CorrMixed     1.1
# corrplot   corrplot    0.92
# gridExtra gridExtra     2.3
# gtsummary gtsummary   2.0.0
# mfp             mfp 1.5.4.1
# tidyverse tidyverse   2.0.0

# Setup environment
source("scripts/setup_environment.R")

# Load required packages using librarian
shelf(tidyverse, mfp, CorrMixed, gridExtra) 


###### breakdown of SpyCATS participants by age groups: 

final_dems %>%
    group_by(age_grp) %>%
    summarise(n = n())

# import baseline titres in those participants with no Strep A positive disease events at baseline

### IMPORTANT: Data Anonymization for Public Data Sharing

# To protect the confidentiality of individual participants, 
# this dataframe has been modified before public release.
# Specifically, we have replaced the actual ages with a "pseudo-age".
# For each participant, a random age is generated within their corresponding age group.
# This approach maintains the overall age group distribution and supports reproducibility of analyses,
# while significantly reducing the risk of re-identifying individual participants.

df <- readRDS("data/baseline_blood_no_disease_titres.RDS")  

# create an empty dataframe to append the output of your for loop
centile_df <- data.frame()

# Create an empty list to store plots
plots <- list()
plots2 <- list()
# loop through each antigen

for (a in unique(df$Antigen)) {        
    
    
    # create a datframe for that antigen
    ctrl <- df %>%
        filter(Antigen == a) %>%
        arrange(age)                                
    
    common_y_limits <- c(min(df$titre), max(df$titre))
    
    
    # fit a fractional polynomial model 
    mfp_model <- mfp(titre ~ fp(age), data = ctrl, family =  gaussian())  
    # add the residuals to filtered dataframe 
    ctrl$.resid <- residuals(mfp_model)             
    # add the fitted values to the dataframe
    ctrl$.fitted <- fitted(mfp_model)              
    
    # Calculate RMSE
    rmse <- sqrt(mean(ctrl$.resid^2))       
    
    # Calculate centiles
    centile_values <- quantile(ctrl$.resid / rmse, c(0.025, 0.50, 0.80, 0.975))
    
    # Add the upper prediction bands based on centiles and update centile dataframe
    ctrl <- ctrl %>%
        mutate(
            upp25 = .fitted + centile_values[1] * rmse,
            upp50 = .fitted + centile_values[2] * rmse,
            upp80 = .fitted + centile_values[3] * rmse,
            upp975 = .fitted + centile_values[4] * rmse
        )
    
    # Append the centile data for this antigen to the centile dataframe
    centile_df <- rbind(centile_df, ctrl)
    
    # Plot
    p <- ggplot(ctrl, aes(x = age, y = titre, col = Antigen)) +
        scale_colour_manual(values = c("SpyCEP" = "#FDC086", "SpyAD" = "#d19c2f", "SLO" = "#386CB0", "GAC" = "#7FC97F", "DNAseB" = "#BEAED4")) +
        guides(color = "none") +
        geom_point(alpha = 0.5) +
        geom_line(aes(y = .fitted), linetype = "dashed", color = "black") +
        geom_line(aes(y = upp25), linetype = "dashed", color = "red", alpha = 0.5) +
        geom_line(aes(y = upp80), linetype = "dashed", color = "blue") +
        geom_line(aes(y = upp975), linetype = "dashed", color = "red", alpha = 0.5) +
        labs(title = paste(a), x = "", y = "IgG level (RLU/mL)") +
        theme_minimal() +
        scale_y_continuous(limits =common_y_limits, breaks = c(1:6), labels = log10_to_exp) +
        theme_universal(base_size = plot_basesize)
    
    # Store the plot in the list with the antigen name as the key
    
    plots[[a]] <- p 
    
    
    
}

# Arrange the plots into a grid layout

plot_02_main01_V1.0 <- do.call(grid.arrange, c(plots, ncol = 5))


#repeat with an age cut off of under 16s. 

for (a in unique(df$Antigen)) {        
    
    
    # create a datframe for that antigen
    ctrl <- df %>%
        filter(Antigen == a) %>%
        arrange(age)                                
    
    
    common_y_limits <- c(min(df$titre), max(df$titre))
    
    
    # fit a fractional polynomial model 
    mfp_model <- mfp(titre ~ fp(age), data = ctrl, family =  gaussian())  
    # add the residuals to filtered dataframe 
    ctrl$.resid <- residuals(mfp_model)             
    # add the fitted values to the dataframe
    ctrl$.fitted <- fitted(mfp_model)              
    
    # Calculate RMSE
    rmse <- sqrt(mean(ctrl$.resid^2))       
    
    # Calculate centiles
    centile_values <- quantile(ctrl$.resid / rmse, c(0.025, 0.50, 0.80, 0.975))
    
    # Add the upper prediction bands based on centiles and update centile dataframe
    ctrl <- ctrl %>%
        mutate(
            upp25 = .fitted + centile_values[1] * rmse,
            upp50 = .fitted + centile_values[2] * rmse,
            upp80 = .fitted + centile_values[3] * rmse,
            upp975 = .fitted + centile_values[4] * rmse
        )
    
    # Append the centile data for this antigen to the centile dataframe
    centile_df <- rbind(centile_df, ctrl)
    
    
    # Plot
    p2 <- 
        
        ctrl %>%
        filter(age <=15) %>%
        ggplot(aes(x = age, y = titre, col = Antigen)) +
        scale_colour_manual(values = c("SpyCEP" = "#FDC086", "SpyAD" = "#d19c2f", "SLO" = "#386CB0", "GAC" = "#7FC97F", "DNAseB" = "#BEAED4")) +
        guides(color = "none") +
        geom_point(alpha = 0.5) +
        geom_line(aes(y = .fitted), linetype = "dashed", color = "black") +
        geom_line(aes(y = upp25), linetype = "dashed", color = "red", alpha = 0.5) +
        geom_line(aes(y = upp80), linetype = "dashed", color = "blue") +
        geom_line(aes(y = upp975), linetype = "dashed", color = "red", alpha = 0.5) +
        labs(x = "Age", y = "IgG level (RLU/mL)") +
        theme_minimal() +
        scale_y_continuous(limits =common_y_limits,breaks = c(1:6), labels = log10_to_exp) +
        theme_universal(base_size = plot_basesize)
    
    # Store the plot in the list with the antigen name as the key
    
    plots2[[a]] <- p2 
    
    
    
}

plot_02_main02_V1.0 <- do.call(grid.arrange, c(plots2, ncol = 5))

plot_02_main02_V1.0

##### ##### ##### ##### ##### 
##### correlation matrix #### 
##### ##### ##### ##### ##### 

##### Tidy data to allow a correlation matrix to be draw for each individual at baseline

df.text <- df %>%
    select(pid, Antigen, titre) %>%
    spread(key = Antigen, value = titre) %>%
    ungroup() %>%
    select(-pid)

head(df.text)

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
    "When comparing baseline antibody titres within individuals, a strong correlation was observed between antibodies to all antigens, ranging from %.2f to %.2f",
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
        "When comparing baseline antibody titres within individuals, a strong correlation was observed between antibodies to all antigens, ranging from %.2f to %.2f (p < 0.0001 for all comparisons).",
        min(cor_pval_results$Correlation, na.rm = TRUE),
        max(cor_pval_results$Correlation, na.rm = TRUE)
    )
} else {
    sentence <- sprintf(
        "When comparing baseline antibody titres within individuals, a strong correlation was observed between antibodies to all antigens, ranging from %.2f to %.2f (p < %.4f for all comparisons).",
        min(cor_pval_results$Correlation, na.rm = TRUE),
        max(cor_pval_results$Correlation, na.rm = TRUE),
        maxp$max_p
    )
}

# Print the sentence
print(sentence)

# Save the corrplot as an image
png("R_output/Supp_fig_baselineIgG_corplot_211124.png", width = 1500, height = 1500, res = 300)

# Plot the correlation matrix for interactions
corrplot::corrplot.mixed(cor_matrix, upper = "color", lower = "number",
                         tl.col = "black", tl.srt = 45, lower.col = "black", number.cex = 1.3,
                         tl.pos = "lt",tl.cex = 1.5)  # "lt" for left and top labels

dev.off()

##############################################################################
#### Describe baseline titres in relation to putative 50% protective thresholds 
##############################################################################

# Import 50% putative thresholds: generated in 07_protection.R scripts 

protective_threshold <- readRDS("R_objects/bloodIgG_protective_threshold_df.RDS")



# Calculate percentages of participants at study baseline with IgG above
# putative 50% threshold for each Antigen

df %>%
    filter(Antigen %in% c("SLO", "SpyAD", "SpyCEP")) %>%
    left_join(protective_threshold, by = join_by(Antigen)) %>%  # Correct join
    mutate(above_threshold = ifelse(titre > Threshold, 1, 0)) %>%
    ungroup() %>%
    select(Antigen, above_threshold) %>%
    gtsummary::tbl_summary(by = Antigen)

percentage_labels <- df %>%
    filter(Antigen %in% c("SLO", "SpyAD", "SpyCEP")) %>%
    left_join(protective_threshold, by = join_by(Antigen)) %>%  # Correct join
    mutate(above_threshold = ifelse(titre > Threshold, 1, 0)) %>%
    group_by(Antigen) %>%
    summarise(
        n = n(),
        n_above = sum(above_threshold),  # Sum the `above_threshold` column to count occurrences
        percentage_above = mean(above_threshold) * 100  # Calculate the mean percentage
    ) %>%
    mutate(label = paste0(round(percentage_above, 0), "%"))  # Create a label for annotation


# Dynamically create summary sentence
sentence <- sprintf(
    "At baseline, antibody levels from the cohort showed that %d individuals (%s) for SLO, %d (%s) for SpyAD, and %d (%s) for SpyCEP had IgG levels above the protective thresholds (Figure 5A).",
    percentage_labels$n_above[percentage_labels$Antigen == "SLO"], 
    percentage_labels$label[percentage_labels$Antigen == "SLO"], 
    percentage_labels$n_above[percentage_labels$Antigen == "SpyAD"], 
    percentage_labels$label[percentage_labels$Antigen == "SpyAD"], 
    percentage_labels$n_above[percentage_labels$Antigen == "SpyCEP"], 
    percentage_labels$label[percentage_labels$Antigen == "SpyCEP"]
)

# Print the sentence
print(sentence)

# Demonstrate the percentage above 50% thresholds at each age group for each antigen 

df %>%
    filter(Antigen %in% c("SLO", "SpyAD", "SpyCEP")) %>%
    left_join(protective_threshold) %>%
    mutate(above_threshold = ifelse(titre > Threshold, 1,0)) %>%
    group_by(age_grp, Antigen) %>%
    summarise(
        proportion_protected = round(mean(above_threshold)*100, 0)) %>%
    spread(Antigen, proportion_protected)


# Plot age vs titre and the proportion above the putatuve 50% protective threshold

df %>%
    filter(Antigen %in% c("SLO", "SpyAD", "SpyCEP")) %>%
    left_join(protective_threshold, by = join_by(Antigen)) %>%  # Correct join
    mutate(above_threshold = ifelse(titre > Threshold, 1, 0)) %>%
    ggplot(aes(x = age, y = titre, col = factor(above_threshold))) +  # Use factor to color by above/below threshold
    guides(color = "none") +
    facet_wrap(~Antigen) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = c("#d73027","#7570b3")) +
    labs(
        y = "IgG level (log10 RLU/mL)",
        x = 'Age'
        # title = "Blood IgG titres above 50% protective threshold by age group"
    ) +
    geom_hline(aes(yintercept = Threshold), linetype = "dashed", color = "red") +  # Reference Threshold correctly
    theme_minimal() +
    theme_universal(base_size = plot_basesize) +
    # Annotate percentages
    geom_text(
        data = percentage_labels,
        aes(x = Inf, y = Inf, label = label),  # Place in top-right corner
        inherit.aes = FALSE, hjust = 1.2, vjust = 1.2, 
        size = plot_basesize
    )


# Plot age group vs titre and the proportion above the putatuve 50% protective threshold



df %>%
    filter(Antigen %in% c("SLO", "SpyAD", "SpyCEP")) %>%
    left_join(protective_threshold) %>%
    mutate(above_threshold = ifelse(titre > Threshold, 1,0)) %>%
    group_by(age_grp, Antigen) %>%
    summarise(
        proportion_protected = mean(above_threshold)) %>%
    ggplot(aes(x = age_grp, y = proportion_protected, fill = Antigen)) +
    scale_fill_manual(values = c("SpyCEP" = "#FDC086", "SpyAD" = "#d19c2f", "SLO" = "#386CB0", "GAC" = "#7FC97F", "DNAseB" = "#BEAED4")) +
    guides(fill = "none") +
    facet_wrap(~Antigen) +
    geom_col(alpha = 0.7) +
    labs(
        title = "Proportion of participants with baseline titres above 50% protected threshold",
        x = "Age Group",
        y = "Proportion"
    ) +
    theme_minimal() 





# Define a function to keep only plot-related objects

# List all objects in the environment
all_objects <- ls()

# Identify objects that contain "plot" in their name
plot_objects <- grep("plot_", all_objects, value = TRUE, ignore.case = TRUE)

# Remove all objects that do not contain "plot" in their name
rm(list = setdiff(all_objects, plot_objects), envir = .GlobalEnv)

rm(all_objects,plot_objects,keep_plot_objects)




