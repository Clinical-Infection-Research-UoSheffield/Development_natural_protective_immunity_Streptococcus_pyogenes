# Title: 10 fold threshold validation 
# Version: 1.0
# Date: 2025-03_26
# Author: Alexander J Keeley

# Reivewer comment: 
# For figure 4 thresholds, a 10-fold cross-validation would help build enhanced statistical confidence in the established threshold. 


# Setup the environement and load functions
source("scripts/setup_environment.R")
source("scripts/load_functions.R")

# Set output directory

output_dir <- "R_output/"


# load packages:

shelf(officer,flextable, lme4)


shelf("dplyr", "tibble", "lubridate", "lme4", "readr")

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
    
    # Prepare event data
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
    
    titre_breakpoint <- titre_breakpoint_df %>%
        filter(Antigen == antigen) %>%
        pull(transition_point)
    
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
    
    # Predict max and threshold probability
    max_prob <- predict(model, newdata = data.frame(
        titre_below_threshold = titre_breakpoint, 
        titre_above_threshold = 0), type = "response", re.form = NA)
    
    threshold_prob <- 0.5 * max_prob
    
    titre_seq <- seq(0, max(final_df_5$titre_above_threshold, na.rm = TRUE), length.out = 1000)
    probs_for_threshold <- predict(model, newdata = data.frame(
        titre_below_threshold = titre_breakpoint,
        titre_above_threshold = titre_seq), type = "response", re.form = NA)
    
    titre_50 <- titre_seq[which.min(abs(probs_for_threshold - threshold_prob))]
    
    return(titre_50 + titre_breakpoint)
}

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

for (ag in antigen_list) {
    res <- cross_validate_threshold(
        full_data = readRDS("data/blood_IgG_titres.RDS"),
        antigen = ag
    )
    
    cv_results_list[[ag]] <- res
}

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

# Construct the sentence
final_sentence <- sprintf(
    "These thresholds were %s for SLO, %s for SpyAD, and %s for SpyCEP.",
    formatted_thresholds$sentence[formatted_thresholds$Antigen == "SLO"],
    formatted_thresholds$sentence[formatted_thresholds$Antigen == "SpyAD"],
    formatted_thresholds$sentence[formatted_thresholds$Antigen == "SpyCEP"]
)

cat(final_sentence)

protective_threshold_df <- 
    formatted_thresholds %>%
    select(Antigen, Threshold = Mean_Threshold)

df <- readRDS("R_objects/baseline_blood_no_disease_titres.RDS")


# protective_threshold <- readRDS("R_objects/bloodIgG_protective_threshold_df.RDS")

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




############# updating plot piece wise regression 



plot_probabilites_piecewise <- function(path_to_titre.df, sample, class, next_event_window = 45, var_name = "titre", antigen, df = 1, slice_window = 0.1) {
    
    # Load and process data
    fun_titres <- read.titres(path_to_titre.df, var_name) %>%
        left_join(start_dates) %>%
        left_join(age) %>%
        mutate(visit_date = as.numeric(ymd(as.character(visit_date))) - entry_1) %>%
        group_by(Antigen) %>%
        ungroup()
    
    antigen_df <- fun_titres %>% filter(Antigen == antigen & !is.na(pid))
    pos_incidence_zero <- readRDS("R_objects/SpyCATS_incidence_df.RDS")
    
    # Define helper: estimate_protective_threshold (same as your function)
    estimate_protective_threshold <- function(data_subset, antigen, next_event_window = 45) {
        titre_breakpoint <- titre_breakpoint_df %>%
            filter(Antigen == antigen) %>%
            pull(transition_point)
        
        fun_titres <- data_subset %>%
            left_join(start_dates) %>%
            left_join(age) %>%
            mutate(visit_date = as.numeric(ymd(as.character(visit_date))) - entry_1) %>%
            group_by(Antigen) %>%
            ungroup()
        
        antigen_df <- fun_titres %>% filter(Antigen == antigen & !is.na(pid)) 
        
        fun_df <- pos_incidence_zero %>%
            select(pid, date, gas_event) %>%
            arrange(pid, date) %>%
            group_by(pid) %>%
            mutate(
                next_date = lead(date),
                event_next_n = case_when(
                    next_date - date <= next_event_window & lead(gas_event) == 1 ~ 1,
                    next_date  - date <= next_event_window & lead(gas_event) == 0 ~ 0,
                    TRUE ~ NA
                )) %>%
            select(-next_date) %>%
            ungroup()
        
        final_df_3 <- fun_df %>%
            left_join(antigen_df %>% select(pid, date = visit_date, titre, age)) %>%
            group_by(pid) %>%
            fill(titre) %>%
            ungroup() %>%
            mutate(hid = substring(pid, 1, 3))
        
        final_df_5 <- final_df_3 %>%
            mutate(titre_below_threshold = ifelse(titre <= titre_breakpoint, titre, titre_breakpoint),
                   titre_above_threshold = ifelse(titre > titre_breakpoint, titre - titre_breakpoint, 0))
        
        model <- glmer(event_next_n ~ titre_above_threshold + (1 | pid) + (1 | hid), data = final_df_5, family = binomial)
        
        max_prob <- predict(model, newdata = data.frame(
            titre_below_threshold = titre_breakpoint, titre_above_threshold = 0), type = "response", re.form = NA)
        
        threshold_prob <- 0.5 * max_prob
        
        titre_seq <- seq(0, max(final_df_5$titre_above_threshold, na.rm = TRUE), length.out = 1000)
        probs_for_threshold <- predict(model, newdata = data.frame(
            titre_below_threshold = titre_breakpoint, titre_above_threshold = titre_seq), type = "response", re.form = NA)
        
        titre_50 <- titre_seq[which.min(abs(probs_for_threshold - threshold_prob))]
        return(titre_50 + titre_breakpoint)
    }
    
    # Define helper: cross_validate_threshold
    cross_validate_threshold <- function(full_data, antigen, k = 10) {
        set.seed(342)
        folds <- sample(rep(1:k, length.out = nrow(full_data)))
        threshold_estimates <- numeric(k)
        
        for (i in 1:k) {
            train_data <- full_data[folds != i, ]
            threshold_estimates[i] <- estimate_protective_threshold(train_data, antigen)
        }
        
        tibble(
            Antigen = antigen,
            Mean_Threshold = mean(threshold_estimates, na.rm = TRUE),
            SD = sd(threshold_estimates, na.rm = TRUE),
            CI_Lower = quantile(threshold_estimates, 0.025, na.rm = TRUE),
            CI_Upper = quantile(threshold_estimates, 0.975, na.rm = TRUE),
            Thresholds = list(threshold_estimates)
        )
    }
    
    # Run cross-validation for this antigen
    cv_result <- cross_validate_threshold(antigen_df, antigen)
    titre_mean <- cv_result$Mean_Threshold
    titre_CI_low <- cv_result$CI_Lower
    titre_CI_high <- cv_result$CI_Upper
    
    
    # Extract Odds Ratios, CIs, and p-value
    model_summary <- summary(model)
    ORs <- exp(model_summary$coefficients[, "Estimate"])
    CIs <- exp(confint(model, parm = "titre", method = "Wald"))
    p_value <- coef(summary(model))["titre", "Pr(>|z|)"]  # Extract p-value for titre
    
    # Check if p-value is significant (less than 0.05) and format accordingly
    p_label <- if (p_value < 0.05) {
        paste0("p = ", sprintf("%.3f", p_value), "*")
    } else {
        paste0("p = ", sprintf("%.3f", p_value))
    }
    
    # Create a sequence of titre values for predictions
    new_data <- data.frame(titre = seq(min(final_df_3$titre, na.rm = TRUE), 
                                       max(final_df_3$titre, na.rm = TRUE), length.out = 100))
    new_data$prob <- predict(model, newdata = new_data, type = "response", re.form = NA)
    
    conf_int <- predict(model, newdata = new_data, type = "link", re.form = NA, se.fit = TRUE)
    new_data$lower <- plogis(conf_int$fit - 1.96 * conf_int$se.fit)
    new_data$upper <- plogis(conf_int$fit + 1.96 * conf_int$se.fit)
    
  
    final_df_5 <-
        final_df_3 %>%
        mutate(titre_below_threshold = ifelse(titre <= titre_breakpoint, titre, titre_breakpoint),
               titre_above_threshold = ifelse(titre > titre_breakpoint, titre - titre_breakpoint, 0))
    
    
    
    # Fit the logistic regression model
    model <- lme4::glmer(event_next_n ~ titre_below_threshold+titre_above_threshold + (1 | pid) + (1 | hid), data = final_df_5, family = binomial)
    
    tb1 <- model %>%
        tbl_regression(exponentiate = TRUE) %>%
        add_nevent() %>%
        bold_p(t = 0.05) %>%
        modify_header(list(label ~ paste(sample, class, "titre"))) 
    
    tb1
    
    
    
    
    # Extract Odds Ratios, CIs, and p-value for titre_above_threshold
    model_summary <- summary(model)
    OR_above <- exp(model_summary$coefficients["titre_above_threshold", "Estimate"])
    CI_above <- exp(confint(model, parm = "titre_above_threshold", method = "Wald"))
    p_value_above <- model_summary$coefficients["titre_above_threshold", "Pr(>|z|)"]
    
    p_label_above <- if (p_value_above < 0.05) {
        paste0("p = ", sprintf("%.3f", p_value_above), "*")
    } else {
        paste0("p = ", sprintf("%.3f", p_value_above))
    }
    
    
    
    
    OR_below <- exp(model_summary$coefficients["titre_below_threshold", "Estimate"])
    CI_below <- exp(confint(model, parm = "titre_below_threshold", method = "Wald"))
    p_value_below <- model_summary$coefficients["titre_below_threshold", "Pr(>|z|)"]
    
    p_label_below <- if (p_value_below < 0.05) {
        paste0("p = ", sprintf("%.3f", p_value_below), "*")
    } else {
        paste0("p = ", sprintf("%.3f", p_value_below))
    }
    
    
    # Create a sequence of titre values for predictions above threshold
    new_data <- data.frame(
        titre = seq(titre_breakpoint, max(final_df_5$titre, na.rm = TRUE), length.out = 100))
    
    new_data$titre_below_threshold <- ifelse(new_data$titre < titre_breakpoint, new_data$titre, titre_breakpoint)
    new_data$titre_above_threshold <- new_data$titre - titre_breakpoint
    
    new_data$prob <- predict(model, newdata = new_data, type = "response", re.form = NA)
    
    conf_int <- predict(model, newdata = new_data, type = "link", re.form = NA, se.fit = TRUE)
    new_data$lower <- plogis(conf_int$fit - 1.96 * conf_int$se.fit)
    new_data$upper <- plogis(conf_int$fit + 1.96 * conf_int$se.fit)
    
    
    ########## Plot it #############
    
    uxlm <- max(final_df_3$titre, na.rm = T)
    lxlm <- 3 # min(final_df_3$titre, na.rm = T)
    
    
    plot <- ggplot(new_data, aes(x = titre, y = prob)) +
        geom_col(data = probabilities, aes(x = Threshold, y = Probability),alpha = 0.5, fill = "#138aa8") +
        # geom_line(color = "blue") +
        geom_vline(xintercept = titre_breakpoint, color = "red", linetype = "dashed") +
        #    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
        labs(title = paste(antigen),
             x = paste0("IgG threshold (log10 RLU/mL)"),
             y = paste0("Proportion with event within ", next_event_window, " days")) +
        theme_minimal() +
        xlim(lxlm,uxlm) +
        guides(fill = "none", alpha = "none")  +
        theme_universal(base_size = plot_basesize)
    #   annotate("text", x = 4.1 + 0.2, y = max(new_data$prob) * 1.6, 
    #           label = paste0("OR (above breakpoint): ", round(OR_above, 2), 
    #                         "\n95% CI: [", round(CI_above[1], 2), ", ", round(CI_above[2], 2), "]", 
    #                        "\n", p_label_above), 
    #        hjust = 0)
    
    
    
    plot
    
    
    
    # Use titre_mean and CIs for annotation and plotting
    ### now modify your plot2 block like this:
    plot2 <- ggplot(new_data, aes(x = titre, y = prob)) +
        geom_line(color = "blue") +
        geom_vline(xintercept = titre_breakpoint, color = "red", linetype = "dashed") +
        geom_vline(xintercept = titre_mean, color = "blue", linetype = "dashed") +
        geom_errorbarh(aes(xmin = titre_CI_low, xmax = titre_CI_high, y = 0.025), color = "blue", height = 0.01) +
        theme_minimal() +
        theme_universal(base_size = plot_basesize)
    
    # Add back the rest of your plotting logic (plot1, plot3, etc.), annotations, and return objects
    # And make sure the output includes cv_result
    
    
    plot3 <- final_df_3 %>%
        ggplot(
            aes(
                x = titre
            )
        ) +
        geom_density(fill = StrepA_colscheme[antigen], alpha = 0.8) +
        geom_vline(xintercept = titre_breakpoint, color = "red", linetype = "dashed") +
        xlim(lxlm,uxlm) +
        labs(
            x = "IgG level(log10 RLU/mL)") +
        
        theme_minimal() +
        theme_universal(base_size = plot_basesize)
    
    plot3
    
    #combi_plot <- patchwork::wrap_plots(plot,plot2,plot3, ncol = 1, heights = c(1, 1, 0.3))
    #
    #
    ## combi_plot <- cowplot::plot_grid(plot,plot2,plot3, ncol = 1, rel_heights = c(1, 1, 0.3), labels = c("A","B","C"))
    #
    #results_list <- list("plot" = combi_plot)
    #
    #return(results_list)
    #
    #}
    
    
    return(list(
        plot1 = plot,
        plot2 = plot2,
        plot3 = plot3,
        tb1 = tb1,
        ptdf = cv_result
    ))
}



# Blood IgG

Antigen_list <- c("SLO" , "SpyAD",  "SpyCEP")

plot_list1 <- list()
plot_list2 <- list()
plot_list3 <- list()
tbl_list <- list()
protective_threshold_df <- tibble()


# Loop through each Antigen and apply the function: mixed_effects_protection_glmer

for (ag in Antigen_list) {
    
    results <- plot_probabilites_piecewise(path_to_titre.df = "data/blood_IgG_titres.RDS",
                                           sample = "Blood",
                                           class = "IgG",
                                           var_name = "titre",
                                           antigen = ag,
                                           slice_window = 0.1)
    
    
    
    # Extract and store the table and plot
    plot_list1[[ag]] <- results$plot1
    plot_list2[[ag]] <- results$plot2
    plot_list3[[ag]] <- results$plot3
    tbl_list[[ag]] <- results$tb1
    
    protective_threshold_df <- bind_rows(protective_threshold_df, results$ptdf)
    
}
# Combine the gtsummary tables
combined_table <- tbl_stack(tbl_list)
combined_table
library(cowplot)

# Combine the plots for each row
row1 <- cowplot::plot_grid(plotlist = plot_list1, labels = "A", nrow = 1) # Unpacks plot_list1
row2 <- cowplot::plot_grid(plotlist = plot_list2, labels = "B", nrow = 1) # Unpacks plot_list2
row3 <- cowplot::plot_grid(plotlist = plot_list3, labels = "C", nrow = 1) # Unpacks plot_list3

row2
# Combine the rows into one grid
combined_plot <- cowplot::plot_grid(row1, row2, row3, ncol = 1,rel_heights = c(1, 1, 0.5), align = "v")



print(combined_plot)

main_07_fig4A_V3 <- combined_plot

main_07_fig4A_V3


###### get 50% putatuve protective thresholds

pt<- protective_threshold_df  %>%
    mutate(Threshold_RLU = 10^Threshold)

protective_threshold_df  %>%
    saveRDS("data/bloodIgG_protective_threshold_df.RDS")

sprintf("The 50%% protective thresholds were %d RLU/mL for SLO, %d RLU/mL for SpyAD, and %d RLU/mL for SpyCEP",
        as.integer(pt$Threshold_RLU[pt$Antigen == "SLO"]),
        as.integer(pt$Threshold_RLU[pt$Antigen == "SpyAD"]),
        as.integer(pt$Threshold_RLU[pt$Antigen == "SpyCEP"]))

