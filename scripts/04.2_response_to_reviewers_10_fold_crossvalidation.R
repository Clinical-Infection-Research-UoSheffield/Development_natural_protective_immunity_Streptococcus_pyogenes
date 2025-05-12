# Title: Cross-Validated Estimation of Protective IgG Thresholds
# Version: 1.0
# Date: 2025-03-26
# Author: Alexander J Keeley

# Reviewer Comment:
#   For Figure 4 thresholds, a 10-fold cross-validation would help build enhanced 
#   statistical confidence in the established threshold.
#
# Purpose:
#   This script performs 10-fold cross-validation to estimate antigen-specific protective 
#   IgG titre thresholds with greater statistical robustness. It addresses reviewer requests 
#   by quantifying uncertainty around the 50% protection threshold and defining 95% confidence 
#   intervals for each antigen.
#
# Key Objectives:
#   - Fit logistic regression models to estimate protective thresholds for each antigen.
#   - Implement 10-fold cross-validation to produce bootstrapped thresholds.
#   - Exponentiate log10 thresholds into RLU/mL for interpretability.
#   - Quantify the proportion of the baseline cohort with titres above the protective threshold.
#   - Stratify protective proportions by age group and visualize titres against age.
#
# Inputs:
#   - data/blood_IgG_titres.RDS
#   - data/baseline_blood_no_disease_titres.RDS
#   - data/SpyCATS_incidence_df.RDS
#   - data/incidence_start_dates.RDS
#
# Outputs:
#   - Formatted table of mean thresholds and 95% CI for each antigen
#   - Age-stratified summaries of protection proportions
#   - ggplot figure (Figure 4E): titres vs age, with protection threshold overlayed
#   - Sentence for manuscript text summarizing thresholds and protection estimates


# Setup the environement and load functions
source("scripts/setup_environment.R")
source("scripts/load_functions.R")

# Set output directory

output_dir <- "R_output/"

### LOAD DATA

# Load reference dates for incidence window alignment
start_dates <- readRDS("data/incidence_start_dates.RDS")
titre_breakpoint_df <- tibble(Antigen = c("SLO","SpyAD", "SpyCEP", "GAC","DNAseB"), transition_point = c(4.3,4.1,4.3,3,3))


# load packages:

shelf(officer,flextable, lme4,dplyr,tibble,lubridate,lme4,readr)


# Define function to estimate 50% protection threshold from training data

estimate_protective_threshold <- function(data_subset, antigen, next_event_window = 45) {
    
    
    # Read titre data and merge with start dates and age, adjust visit dates, and categorize titres
    fun_titres <- data_subset %>%
        left_join(start_dates) %>%
        left_join(age) %>%
        mutate(visit_date = as.numeric(ymd(as.character(visit_date))) - entry_1) %>%
        group_by(Antigen) %>%
        ungroup()
    
    # Filter for specific antigen
    antigen_df <- fun_titres %>% filter(Antigen == antigen & !is.na(pid)) 
    
    pos_incidence_zero <- readRDS("data/SpyCATS_incidence_df.RDS")
    
    # Prepare event data create a binary indicator of whether GAS event occurs within 45 days
    fun_df <- pos_incidence_zero %>%
        select(pid, date, gas_event) %>%
        arrange(pid, date) %>%
        group_by(pid) %>%
        mutate(
            next_date = lead(date),
            check = next_date - date,
            event_next_n = case_when(
                next_date - date <= next_event_window & lead(gas_event) == 1 ~ 1,
                next_date  - date <= next_event_window & lead(gas_event) == 0 ~ 0,
                next_date -date > next_event_window ~ NA
            )) %>%
        select(-next_date) %>%
        ungroup()
    
    # Extract the predefined titre breakpoint for piecewise regression
    titre_breakpoint <- titre_breakpoint_df %>%
        filter(Antigen == antigen) %>%
        pull(transition_point)
    
    # Merge antibody data with event data, fill forward missing covariates within each individual - assuming titres remain constant between meaasurements
    
    final_df_3  <-fun_df %>%
        left_join(antigen_df %>% 
                      select(pid, date = visit_date, titre, age)) %>%
        group_by(pid) %>%
        fill(titre) %>%
        ungroup() %>%
        mutate(hid = substring(pid, 1, 3))
    
    final_df_5 <- final_df_3 %>%
        mutate(titre_below_threshold = ifelse(titre <= titre_breakpoint, titre, titre_breakpoint),
               titre_above_threshold = ifelse(titre > titre_breakpoint, titre - titre_breakpoint, 0))
    
    # Fit logistic model to above-threshold data
    model <- glmer(event_next_n ~ titre_above_threshold + (1 | pid) + (1 | hid),
                   data = final_df_5, family = binomial)
    
    # Compute maximum predicted probability (i.e., reference point)
    max_prob <- predict(model, newdata = data.frame(
        titre_below_threshold = titre_breakpoint, 
        titre_above_threshold = 0), type = "response", re.form = NA)
    
    # Define 50% threshold of max predicted protection
    threshold_prob <- 0.5 * max_prob
    
    # Search for titre at which predicted probability reaches this threshold
    titre_seq <- seq(0, max(final_df_5$titre_above_threshold, na.rm = TRUE), length.out = 1000)
    probs_for_threshold <- predict(model, newdata = data.frame(
        titre_below_threshold = titre_breakpoint,
        titre_above_threshold = titre_seq), type = "response", re.form = NA)
    
    titre_50 <- titre_seq[which.min(abs(probs_for_threshold - threshold_prob))]
    
    return(titre_50 + titre_breakpoint)
}

# Perform 10-fold cross-validation across all antigens
cross_validate_threshold <- function(full_data, antigen, k = 10) {
    set.seed(342)  # Fixed for reproducibility
    folds <- sample(rep(1:k, length.out = nrow(full_data)))
    
    threshold_estimates <- numeric(k)
    
    for (i in 1:k) {
        train_data <- full_data[folds != i, ]
        
        threshold_estimates[i] <- estimate_protective_threshold(
            data_subset = train_data,
            antigen = antigen
        )
    }
    
    summary_df <- tibble(
        Antigen = antigen,
        Mean_Threshold = mean(threshold_estimates, na.rm = TRUE),
        SD = sd(threshold_estimates, na.rm = TRUE),
        CI_Lower = quantile(threshold_estimates, 0.025, na.rm = TRUE),
        CI_Upper = quantile(threshold_estimates, 0.975, na.rm = TRUE),
        Thresholds = list(threshold_estimates)
    )
    
    return(summary_df)
}

antigen_list <- c("SLO", "SpyAD", "SpyCEP")
cv_results_list <- list()

# Loop through antigens and estimate cross validation-based thresholds
for (ag in antigen_list) {
    res <- cross_validate_threshold(
        full_data = readRDS("data/blood_IgG_titres.RDS"),
        antigen = ag
    )
    
    cv_results_list[[ag]] <- res
}

# Combine results from all antigens
cv_thresholds_df <- bind_rows(cv_results_list)
print(cv_thresholds_df)

# Exponentiate and format the thresholds
formatted_thresholds <- cv_thresholds_df %>%
    mutate(
        Mean_exp = 10^Mean_Threshold,
        CI_Lower_exp = 10^CI_Lower,
        CI_Upper_exp = 10^CI_Upper,
        sentence = sprintf("%.0f (CI %.0f–%.0f) RLU/mL", Mean_exp, CI_Lower_exp, CI_Upper_exp)
    )

formatted_thresholds

# construct sentence for reporting in manuscript text
final_sentence <- sprintf(
    "These thresholds were %s for SLO, %s for SpyAD, and %s for SpyCEP.",
    formatted_thresholds$sentence[formatted_thresholds$Antigen == "SLO"],
    formatted_thresholds$sentence[formatted_thresholds$Antigen == "SpyAD"],
    formatted_thresholds$sentence[formatted_thresholds$Antigen == "SpyCEP"]
)

cat(final_sentence)

# Print sentence for reporting in manuscript text
protective_threshold_df <- 
    formatted_thresholds %>%
    select(Antigen, Threshold = Mean_Threshold)

df <- readRDS("data/baseline_blood_no_disease_titres.RDS")

protective_threshold <-protective_threshold_df

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

# reubttal # 1 = update this to figure 4 - panel E
fig_4E_age_thresholds <- df %>%
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
    geom_hline(aes(yintercept = Threshold), linetype = "dashed", color = "blue") +  # Reference Threshold correctly
    theme_minimal() +
    theme_universal(base_size = plot_basesize) +
    # Annotate percentages
    geom_text(
        data = percentage_labels,
        aes(x = Inf, y = Inf, label = label),  # Place in top-right corner
        inherit.aes = FALSE, hjust = 1.2, vjust = 1.2, 
        size = plot_basesize
    )


fig_4E_age_thresholds










