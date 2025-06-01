# Title: The association between IgG level and culture-confirmed S pyogenes events 
# Version: 2.0
# Date: 2025-01-12
# Author: Alexander J Keeley

# Inputs: 
#   - data/blood_IgG_titres.RDS
#   - data/incidence_start_dates.RDS
#   - data/SpyCATS_incidence_df.RDS
#   - data/survivial_df_blood_IgG.RDS

# Outputs: 
#   - Plots visualizing IgG titre distributions and event probabilities
#   - Statistical models assessing relationships between IgG levels and event risk
#   - Regression tables summarizing Cox PH survival analyses and logistic regression models

# Description:
#
# This script examines the relationship between blood IgG titres and the risk of S pyogenes-related events.
# Key components of the analysis include:
#
# 1. Loading appropriate functions and dataframes
# 2. Calculation of event probabilities at different IgG level thresholds.
# 3. Mixed-effects regression models to assess impact of IgG levels on event occurrence within 45 days.
# 4. Determination of optimal titre thresholds for a piecewise regression using AIC-based model comparisons.
# 5. Piecewise logistic regression analysis including estimatation of putative 50% protective thresholds.
# 6. Orthogonal Anderson Gill extension of Cox proportional hazards modeling to analyze time-to-event relationships.

### Date Anonymisation for Public Data Sharing
#
# To protect the confidentiality of individual participants, all date values have been anonymised.
# This is achieved by adding a constant offset (offset_days) to each date field (e.g., visit_date, incidence date, enrollment date).
# The offset preserves the relative intervals between dates (so the timing and spacing of events remains accurate),
# but it prevents the disclosure of the actual calendar dates.
# This approach ensures that while the temporal relationships within the data are maintained, individual identification is minimized.


# Requirements:

# Package  Version
# cowplot     cowplot    1.1.3
# flextable flextable    0.9.6
# forcats     forcats    1.0.0
# lme4           lme4 1.1-35.5
# officer     officer    0.6.6
# patchwork patchwork    1.2.0


# Setup the environement and load functions
source("scripts/setup_environment.R")
source("scripts/load_functions.R")

# Set output directory

output_dir <- "R_output/"


# load packages:

shelf(officer,flextable,lme4)



# The number of individuals with titre measured 

readRDS("data/blood_IgG_titres.RDS") %>%
    pull(pid) %>%
    unique() %>%
    length()

### Import a dataframe with individual study start dates: 
start_dates <- readRDS("data/incidence_start_dates.RDS")


################################################
########## plot probabilities  #############
################################################

# Input:
# - Antibody titre data (e.g., "data/blood_IgG_titres.RDS") including participant IDs, dates, and titres.
# - Event incidence data ("data/SpyCATS_incidence_df.RDS") with dates of events and participant IDs.
# - Demographic data

# Description:
# - Data cleaning and preparation: Merges antibody titre data with demographic and event data, filters by antigen, and calculates variable (next at next visit within a given time period).
# - Analysis: None

# Output:
# - Visualization plot showing the relationship between titre levels and event probability at each threshold


plot_probabilites <- function(path_to_titre.df, sample, class, next_event_window = 45, var_name = "titre", antigen, df = 1, slice_window = 0.1, protective_threshold = 0.67) {
    
    # Read titre data and merge with start dates and age, adjust visit dates
    fun_titres <- read.titres(path_to_titre.df, var_name) %>%
        left_join(start_dates) %>%
        left_join(age) %>%
        mutate(visit_date = as.numeric(ymd(as.character(visit_date))) - entry_1)
    
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
    
    # Merge antibody data with event data, fill forward missing covariates within each individual - assuming titres remain constant between meaasurements
    
    final_df_3  <-fun_df %>%
        left_join(antigen_df %>% 
                      select(pid, date = visit_date, titre, age)) %>%
        group_by(pid) %>%
        fill(titre) %>%
        ungroup() %>%
        mutate(hid = substring(pid, 1, 3))
    
    # Remove rows with NA in titre or event_next_day
    final_df_4 <- na.omit(final_df_3, cols = c("titre", "event_next_day"))
    
    # Calculate and adjust thresholds for analysis
    min_titre <- min(final_df_4$titre, na.rm = TRUE)
    max_titre <- max(final_df_4$titre, na.rm = TRUE)
    
    # Round down/up to the nearest multiple of slice_window
    rounded_min <- floor(min_titre / slice_window) * slice_window
    rounded_max <- ceiling(max_titre / slice_window) * slice_window
    
    thresholds <- seq(from = rounded_min, to = rounded_max, by = slice_window)
    
    # Initialize a data frame to store results
    probabilities <- data.frame(Threshold = numeric(), Probability = numeric(), percent = numeric())
    
    for (threshold in thresholds) {
        # Subset the data for rows where titre >= threshold
        subset_df <- final_df_4 %>% filter(titre >= threshold)
        
        if (nrow(subset_df) > 0) {
            prob <- sum(subset_df$event_next_n == 1) / nrow(subset_df)
        } else {
            prob <- NA
        }
        
        n_obs <- nrow(subset_df) / nrow(final_df_4) * 100
        probabilities <- rbind(probabilities, data.frame(Threshold = threshold, Probability = prob, percent= n_obs))
    }
    
    max_threshold <- max(probabilities$Threshold[probabilities$percent >= 97.5])
    
    
    ########## Plot it #############
    plot <- ggplot() +
        geom_col(data = probabilities, aes(x = Threshold, y = Probability, alpha = percent), fill = "#138aa8") +
        #    ylim(0, 0.15) +
        xlim(max_threshold, max(thresholds)) +
        labs(title = paste(antigen),
             x = "IgG level (log10 RLU/mL)",
             y = paste0("Proportion experiencing event within ", next_event_window, " days")) +
        guides(fill = "none", alpha = "none") +
        theme_minimal() +
        theme_universal(base_size = plot_basesize)
    
    results_list <- list("plot" = plot)
    
    return(results_list)
}



# Blood IgG

Antigen_list <- c("GAC", "SLO" , "SpyAD",  "SpyCEP","DNAseB")

plot_list <- list()

# Loop through each Anitgen and apply the function: mixed_effects_protection_glmer
for (ag in Antigen_list) {
    
    results <- plot_probabilites(path_to_titre.df = "data/blood_IgG_titres.RDS",
                                 sample = "Blood",
                                 class = "IgG",
                                 var_name = "titre",
                                 antigen = ag,
                                 slice_window = 0.1)
    
    
    
    # Extract and store the table and plot
    plot_list[[ag]] <- results$plot
}

# Combine the plots 
combined_plot <- patchwork::wrap_plots(plot_list, ncol = 5) # adjust ncol as needed

# Print the combined plot

print(combined_plot)


################################################
########## mixed effect regression #############
################################################


# This function performs a mixed effects regression to explore the relationship between
# IgG level and event within the next n days (set to 45)


# Input:
# - Antibody titre data (e.g., "data/blood_IgG_titres.RDS") including participant IDs, dates, and titres.
# - Event incidence data ("data/SpyCATS_incidence_df.RDS") with dates of events and participant IDs.
# - Demographic data

# Description:
# - Data cleaning and preparation: Merges antibody titre data with demographic and event data, filters by antigen, and calculates variable (next at next visit within a given time period).
# - Analysis: Calculates titre thresholds, fits a generalized linear mixed-effects model to evaluate the effect of titres on event_next_n. 

# Output:
# - Regression table summarizing the effects of titre levels on event probability.
# - Visualization plot showing the relationship between titre levels and event probability, with mixed effects regression + confidence intervals plotted.

mixed_effects_protection_glmer <- function(path_to_titre.df, sample, class, next_event_window = 45, var_name = "titre", antigen, df = 1, slice_window = 0.1, protective_threshold = 0.67) {
    
    # Read titre data and merge with start dates and age, adjust visit dates, and categorize titres
    fun_titres <- read.titres(path_to_titre.df, var_name) %>%
        left_join(start_dates) %>%
        left_join(age) %>%
        mutate(visit_date = as.numeric(ymd(as.character(visit_date))) - entry_1)
    
    # Filter for specific antigen
    antigen_df <- fun_titres %>% filter(Antigen == antigen & !is.na(pid)) 
    
    #load events incidence datdframe
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
    
    # Merge antibody data with event data, fill forward missing covariates within each individual - assuming titres remain constant between meaasurements
    
    final_df_3  <-fun_df %>%
        left_join(antigen_df %>% 
                      select(pid, date = visit_date, titre, age)) %>%
        group_by(pid) %>%
        fill(titre) %>%
        ungroup() %>%
        mutate(hid = substring(pid, 1, 3))
    
    # Remove rows with NA in titre or event_next_day
    final_df_4 <- na.omit(final_df_3, cols = c("titre", "event_next_day"))
    
    # Calculate and adjust thresholds for analysis
    min_titre <- min(final_df_4$titre, na.rm = TRUE)
    max_titre <- max(final_df_4$titre, na.rm = TRUE)
    
    # Round down/up to the nearest multiple of slice_window
    rounded_min <- floor(min_titre / slice_window) * slice_window
    rounded_max <- ceiling(max_titre / slice_window) * slice_window
    
    thresholds <- seq(from = rounded_min, to = rounded_max, by = slice_window)
    
    # Initialize a data frame to store results
    probabilities <- data.frame(Threshold = numeric(), Probability = numeric(), percent = numeric())
    
    for (threshold in thresholds) {
        # Subset the data for rows where titre >= threshold
        subset_df <- final_df_4 %>% filter(titre >= threshold)
        
        if (nrow(subset_df) > 0) {
            prob <- sum(subset_df$event_next_n == 1) / nrow(subset_df)
        } else {
            prob <- NA
        }
        
        n_obs <- nrow(subset_df) / nrow(final_df_4) * 100
        probabilities <- rbind(probabilities, data.frame(Threshold = threshold, Probability = prob, percent= n_obs))
    }
    
    max_threshold <- max(probabilities$Threshold[probabilities$percent >= 97.5])
    
    # Fit the logistic regression model
    model <- lme4::glmer(event_next_n ~ titre + (1 | pid) + (1 | hid), data = final_df_3, family = binomial)
    
    # Extract AIC
    model_aic <- AIC(model)
    
    print(paste("sample size:",dim(final_df_3)[1]))
    print(paste(antigen,":", AIC(model)))
    
    tb1 <- model %>%
        tbl_regression(exponentiate = TRUE) %>%
        add_nevent() %>%
        bold_p(t = 0.05) %>%
        modify_header(list(label ~ paste(sample, class, "titre"))) 
    
    tb1$table_body$label <- paste0(antigen)
    
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
    
    ########## Plot it #############
    plot1 <- ggplot(data = probabilities, aes(x = Threshold, y = Probability)) +
        geom_col(fill = "#138aa8") +
        xlim(max_threshold, max(thresholds)) +
        labs(title = paste(antigen),
             x = paste0("IgG threshold (log10 RLU/mL)"),
             y = paste0("Proportion with event within ", next_event_window, " days")) +
        guides(fill = "none", alpha = "none") +
        theme_minimal() +
        theme_universal(base_size = plot_basesize)
    
    
    
    
    plot2 <- ggplot(new_data, aes(x = titre, y = prob)) +
        geom_line(color = "blue") +
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
        ylim(0, 0.15) +
        xlim(max_threshold, max(thresholds)) +
        labs(title = paste(antigen),
             x = paste0("IgG level (log10 RLU/mL)"),
             y = paste0("Probability of event in ", next_event_window, " days")) +
        guides(fill = "none", alpha = "none") +
        theme_minimal() +
        annotate("text", x = max_threshold + 0.5, y = 0.10, 
                 label = paste0("OR: ", round(ORs[2], 2), 
                                "\n95% CI: [", round(CIs[1, 1], 2), ", ", round(CIs[1, 2], 2), "]", 
                                "\n", p_label), 
                 hjust = 0,
                 size = plot_basesize / 2.8) +
        theme_universal(base_size = plot_basesize)
    
    
    results_list <- list("table" = tb1, "plot1" = plot1, "plot2" = plot2, "AIC_glmer" = model_aic)
    
    return(results_list)
}



# Blood IgG 

Antigen_list <- c("GAC", "SLO" , "SpyAD",  "SpyCEP","DNAseB")

table_list <- list()
plot_list1 <- list()
plot_list2 <- list()


# Initialize a dataframe to store AIC values
AIC_glmer <- data.frame(Antigen = character(), AIC_glmer = numeric(), stringsAsFactors = FALSE)



# Loop through each Anitgen and apply the function: mixed_effects_protection_glmer
for (ag in Antigen_list) {
    
    results <- mixed_effects_protection_glmer(path_to_titre.df = "data/blood_IgG_titres.RDS",
                                              sample = "Blood",
                                              class = "IgG",
                                              var_name = "titre",
                                              antigen = ag,
                                              slice_window = 0.1)
    
    
    
    # Extract and store the table and plot
    table_list[[ag]] <- results$table
    plot_list1[[ag]] <- results$plot1
    plot_list2[[ag]] <- results$plot2
    # Append the AIC to the dataframe
    AIC_glmer <- rbind(AIC_glmer, data.frame(Antigen = ag, AIC_glmer = results$AIC_glmer))
    
}

# Combine the gtsummary tables
combined_table <- tbl_stack(table_list)


row1 <- cowplot::plot_grid(plotlist = plot_list1, labels = "A", nrow = 1, label_size = 8) # Unpacks plot_list1
row2 <- cowplot::plot_grid(plotlist = plot_list2, labels = "B", nrow = 1, label_size = 8) # Unpacks plot_list2

# Combine the plots 
combined_plot <- cowplot::plot_grid(row1, row2, ncol = 1,rel_heights = c(1, 1), align = "v")


# Print the combined table and plot
print(combined_table)
print(combined_plot)
print(AIC_glmer)

supp_fig3_protection_V3.0_261124 <- combined_plot
supp_fig3_protection_V3.0_261124 




##### ##### ##### ##### ##### ##### ##### ##### ### 
##### determine optimal breakpoint threshold: ##### 
##### ##### ##### ##### ##### ##### ##### ##### ###

# identify_piecewise_transition: Determines the optimal titre threshold for subsequent piecewise regression.

# Inputs: 
#   - path_to_titre.df: File path to titre data; sample, class, antigen, next_event_window, var_name,
#     slice_window, and protective_threshold as parameters.

# Process:
#   1. Merge titre data with start dates and age; adjust visit dates relative to enrolment.
#   2. Filter by antigen and prepare event data (computing event_next_n for events within next_event_window).
#   3. Iterate over candidate thresolds, split titre into below/above components,
#      fit mixed-effects logistic models (glmer), and record AIC values for each threshold.
# Output: Returns a data frame of candidate breakpoints with corresponding AIC values for optimal threshold selection.


identify_piecewise_transition<- function(path_to_titre.df, sample, class, next_event_window = 45, var_name = "titre", antigen, df = 1, slice_window = 0.1, protective_threshold = 0.67) {
    
    
    # Read titre data and merge with start dates and age, adjust visit dates, and categorize titres
    fun_titres <- read.titres(path_to_titre.df, var_name) %>%
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
    
    
    titre_breakpoint_df <- tibble(Antigen = c("SLO","SpyAD", "SpyCEP"), titre_breakpoint = c(4.3,4.1,4.3))
    
    titre_breakpoint <- titre_breakpoint_df %>%
        filter(Antigen == antigen) %>%
        pull(titre_breakpoint)
    
    # Merge antibody data with event data, fill forward missing covariates within each individual - assuming titres remain constant between meaasurements
    
    final_df_3  <-fun_df %>%
        left_join(antigen_df %>% 
                      select(pid, date = visit_date, titre, age)) %>%
        group_by(pid) %>%
        fill(titre) %>%
        ungroup() %>%
        mutate(hid = substring(pid, 1, 3))
    
    
    AICs <- data.frame(titre_breakpoint = numeric(), AIC = numeric())
    
    for (bp in seq(3, 5, by = 0.1)) {
        final_df_3 <- final_df_3 %>%
            mutate(
                titre_below_threshold = ifelse(titre <= bp, titre, bp),
                titre_above_threshold = ifelse(titre > bp, titre - bp, 0)
            )
        
        model <- tryCatch(
            {
                lme4::glmer(event_next_n ~ titre_below_threshold + titre_above_threshold + (1 | pid) + (1 | hid),
                            data = final_df_3, family = binomial)
            },
            error = function(e) NA
        )
        
        aic <- if (!is.na(model)) AIC(model) else NA
        
        AICs <- rbind(AICs, data.frame(titre_breakpoint = bp, AIC = aic))
    }
    
    
    
    results_list <- list("AICs" = AICs)
    
    return(results_list)
    
}

#all_aics <- data.frame(Antigen = character(), titre_breakpoint = numeric(), AIC = numeric())


# Loop through each Anitgen and apply the function:identify_piecewise_transition
#for (ag in Antigen_list) {
#    
#    results <- identify_piecewise_transition(path_to_titre.df = "data/blood_IgG_titres.RDS",
#                                           sample = "Blood",
#                                           class = "IgG",
#                                           var_name = "titre",
#                                           antigen = ag,
#                                           slice_window = 0.1)
#    
#
#    # Append AIC results to the master dataframe
#    aic_df <- results$AICs  # Extract AIC dataframe from the function output
#    aic_df$Antigen <- ag    # Add the antigen as a column
#    all_aics <- rbind(all_aics, aic_df)  # Append to the master dataframe
#    
#    
#}
#
## create a breakpoint dataframe with visualised breakpoint

titre_breakpoint_df <- tibble(Antigen = c("SLO","SpyAD", "SpyCEP", "GAC","DNAseB"), transition_point = c(4.3,4.1,4.3,3,3))

#all_aics <-all_aics %>%
#    left_join(AIC_glmer) %>%
#    left_join(titre_breakpoint_df)


# Plot AIC values for the model at each iteratnion of titre threhold. Values within 2 are considered comparable

# all_aics %>%
#     ggplot(
#         aes(x = titre_breakpoint, y = AIC)
#     ) +
#     geom_point() +
#     geom_hline(aes(yintercept = AIC_glmer), linetype = "dashed", color = "blue") +
#     geom_hline(aes(yintercept = AIC_glmer -2), linetype = "dashed", color = "green") +
#     geom_vline(aes(xintercept = transition_point), linetype = "dashed", color = "red") +
#     facet_wrap(~ Antigen, scales = "free") +
#     labs(
#         title = "AIC across Titre Breakpoints by Antigen",
#         x = "Titre Breakpoint",
#         y = "AIC"
#     ) +
#    theme_minimal() 


#############################################################################
###### Estimate 50% protective thresholds via k-fold validation model #######
#############################################################################

# estimate_protective_threshold:
# Computes the IgG titre at which the predicted risk of GAS event is 50% of the maximum probability, 
# using a piecewise logistic mixed-effects model fitted above a previsouly defined transition point

# cross_validate_threshold:
# Performs 10-fold cross-validation of the threshold estimation process.
# For each fold, recomputes the 50% protective threshold on the training data.
# Returns the cross-validated mean threshold, standard deviation, and 95% CI.

# Output:
# cv_thresholds_df – a data frame with cross-validated mean 50% protective thresholds and uncertainty bounds for each antigen.

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
    
    # identify the antigen specific transition point 
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
    set.seed(342)  # This is a random number - set for reproducibility
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

saveRDS(protective_threshold_df, "data/bloodIgG_protective_threshold_df.RDS")

###############################################################
########### perform piecewise logistic regression #############
###############################################################

# plot_probabilites_piecewise:
# Evaluates the relationship between IgG titre and subsequent GAS event risk using piecewise logistic regression.

# Inputs: titre data file (path_to_titre.df), incidence dataframe 

# Process:
#   1. Merges IgG titre data with enrolment and demographic information, aligns visit dates relative to baseline.
#   2. Constructs a time-to-event outcome variable (event_next_n) within `next_event_window` using GAS event data.
#   3. Fits piecewise logistic regression with titres segmented below and above a fixed transition point.
#   4. Applies a previously computed 10-fold cross-validated threshold to highlight estimated 50% protective titre.
#   5. Generates:
#      - Plot 1: Empirical probability of subsequent infection vs titre
#      - Plot 2: Model-predicted probability curve with confidence bands, odds ratios, and the cross-validated threshold with CI
#      - Plot 3: Distribution of titre values with the transition point overlayed
#   6. Outputs a results list containing model summaries, plots, and metadata.
#
# Output:
#   A named list with elements:
#   - "plot1": Empirical titre-event relationship
#   - "plot2": Model-predicted protection curve with threshold annotation
#   - "plot3": Titre distribution density
#   - "tb1": gtsummary output of ORs from piecewise model


plot_probabilites_piecewise <- function(path_to_titre.df, sample, class, next_event_window = 45, var_name = "titre", antigen, df = 1, slice_window = 0.1) {
    
    
    # Read titre data and merge with start dates and age, adjust visit dates, and categorize titres
    fun_titres <- read.titres(path_to_titre.df, var_name) %>%
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
    
    
    titre_breakpoint_df <- tibble(Antigen = c("SLO","SpyAD", "SpyCEP"), titre_breakpoint = c(4.3,4.1,4.3))
    
    titre_breakpoint <- titre_breakpoint_df %>%
        filter(Antigen == antigen) %>%
        pull(titre_breakpoint)
    
    # Merge antibody data with event data, fill forward missing covariates within each individual - assuming titres remain constant between meaasurements
    
    final_df_3  <-fun_df %>%
        left_join(antigen_df %>% 
                      select(pid, date = visit_date, titre, age)) %>%
        group_by(pid) %>%
        fill(titre) %>%
        ungroup() %>%
        mutate(hid = substring(pid, 1, 3))
    
    
    
    # Remove rows with NA in titre or event_next_day
    final_df_4 <- na.omit(final_df_3, cols = c("titre", "event_next_day"))
    
    # Calculate and adjust thresholds for analysis
    min_titre <- min(final_df_4$titre, na.rm = TRUE)
    max_titre <- max(final_df_4$titre, na.rm = TRUE)
    
    # Round down/up to the nearest multiple of slice_window
    rounded_min <- floor(min_titre / slice_window) * slice_window
    rounded_max <- ceiling(max_titre / slice_window) * slice_window
    
    thresholds <- seq(from = rounded_min, to = rounded_max, by = slice_window)
    
    # Initialize a data frame to store results
    probabilities <- data.frame(Threshold = numeric(), Probability = numeric(), percent = numeric())
    
    for (threshold in thresholds) {
        # Subset the data for rows where titre >= threshold
        subset_df <- final_df_4 %>% filter(titre >= threshold)
        
        if (nrow(subset_df) > 0) {
            prob <- sum(subset_df$event_next_n == 1) / nrow(subset_df)
        } else {
            prob <- NA
        }
        
        n_obs <- nrow(subset_df) / nrow(final_df_4) * 100
        probabilities <- rbind(probabilities, data.frame(Threshold = threshold, Probability = prob, percent= n_obs))
    }
    
    max_threshold <- max(probabilities$Threshold[probabilities$percent >= 97.5])
    
    
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
    
    
    ########## Plot 2 #############
    
    titre_mean <- cv_thresholds_df$Mean_Threshold[cv_thresholds_df$Antigen == antigen]
    titre_CI_low <- cv_thresholds_df$CI_Lower[cv_thresholds_df$Antigen == antigen]
    titre_CI_high <- cv_thresholds_df$CI_Upper[cv_thresholds_df$Antigen == antigen]
    
    # Predict the event probability at the cross-validated mean threshold
    predicted_prob_mean <- predict(model, newdata = data.frame(
        titre_below_threshold = titre_breakpoint, 
        titre_above_threshold = titre_mean - titre_breakpoint
    ), type = "response", re.form = NA)
    
    
    plot2 <- ggplot(new_data, aes(x = titre, y = prob)) +
        geom_line(color = "blue") +
        geom_vline(xintercept = titre_breakpoint, color = "red", linetype = "dashed") +
        geom_errorbarh(aes(xmin = titre_CI_low, xmax = titre_CI_high, y = 0), color = "blue", height = 0.005) +
        
        geom_segment(aes(x = titre_mean, 
                         xend = titre_mean, 
                         y = 0, 
                         yend = predicted_prob_mean), 
                     color = "blue", linetype = "dashed") +
        
        geom_segment(aes(x = 3, 
                         xend = titre_mean, 
                         y = predicted_prob_mean, 
                         yend = predicted_prob_mean), 
                     color = "blue", linetype = "dashed") +
        
        
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
        labs(
            x = paste0("IgG level (log10 RLU/mL)"),
            y = paste0("Probability of event in ", next_event_window, " days")) +
        theme_minimal() +
        guides(fill = "none", alpha = "none") +
        xlim(lxlm,uxlm) +
        annotate("text", x = titre_breakpoint + 0.15, y = max(new_data$prob) +0.005, 
                 size = plot_basesize /2.8,
                 label = paste0("Above transition:",
                                "\nOR : ", round(OR_above, 2), 
                                "\nCI:", round(CI_above[1], 2), ",", round(CI_above[2], 2), 
                                "\n", p_label_above), 
                 hjust = 0)  +
        annotate("text", x = lxlm + 0.1  , y = max(new_data$prob)+0.005, 
                 size = plot_basesize /2.8,
                 label = paste0("Below transition:",
                                "\nOR : ", round(OR_below, 2), 
                                "\nCI:", round(CI_below[1], 2), ",", round(CI_below[2], 2), 
                                "\n", p_label_below), 
                 hjust = 0) +
        theme_universal(base_size = plot_basesize)
    
    
    
    plot2
    
    ########## Plot 3 #############
    
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
    
    results_list <- list("plot1" = plot, "plot2" = plot2, "plot3" = plot3, "tb1" = tb1, "ptdf" = protective_threshold_df)
    
    return(results_list)
    
}

# Blood IgG

Antigen_list <- c("SLO" , "SpyAD",  "SpyCEP")

plot_list1 <- list()
plot_list2 <- list()
plot_list3 <- list()
tbl_list <- list()
#protective_threshold_df <- tibble()


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

}
# Combine the gtsummary tables
combined_table <- tbl_merge(tbl_list)
combined_table

shelf(cowplot)

# Combine the plots for each row
row1 <- cowplot::plot_grid(plotlist = plot_list1, labels = "A", nrow = 1, label_size = 8) # Unpacks plot_list1
row2 <- cowplot::plot_grid(plotlist = plot_list2, labels = "B", nrow = 1, label_size = 8) # Unpacks plot_list2
row3 <- cowplot::plot_grid(plotlist = plot_list3, labels = "C", nrow = 1, label_size = 8) # Unpacks plot_list3

# Combine the rows into one grid
combined_plot <- cowplot::plot_grid(row1, row2, row3, ncol = 1,rel_heights = c(1, 1, 0.5), align = "v")



print(combined_plot)

main_07_fig4A_V3 <- combined_plot

main_07_fig4A_V3


################################
################################

# explore the relationship with the IgG level above transition point in regression analysis, adjusting for age group, sex and household size and event within next (default =45) days 

glmer_adjusted_piecewise <- function(path_to_titre.df, sample, class, var_name = "titre", antigen, next_event_window = 45) {
    
    # Read titre data and merge with start dates and age, adjust visit dates, and categorize titres
    fun_titres <- read.titres(path_to_titre.df, var_name) %>%
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
    
    
    titre_breakpoint_df <- tibble(Antigen = c("SLO", "SpyAD", "SpyCEP"), titre_breakpoint = c(4.3,4.1,4.3))
    
    titre_breakpoint <- titre_breakpoint_df %>%
        filter(Antigen == antigen) %>%
        pull(titre_breakpoint)
    
    # Merge antibody data with event data, fill forward missing covariates within each individual - assuming titres remain constant between meaasurements
    
    final_df_3  <-fun_df %>%
        left_join(antigen_df %>% 
                      select(pid, date = visit_date, titre, age_grp, sex, hhsize)) %>%
        group_by(pid) %>%
        fill(titre,age_grp,sex,hhsize) %>%
        ungroup() %>%
        mutate(hid = substring(pid, 1, 3))
    
    final_df_5 <-
        final_df_3 %>%
        mutate(titre_below_threshold = ifelse(titre <= titre_breakpoint, titre, titre_breakpoint),
               titre_above_threshold = ifelse(titre > titre_breakpoint, titre - titre_breakpoint, 0))
    
    final_df_5$age_grp <- relevel(final_df_5$age_grp, ref = "Over 18 years")
    
    # Fit the logistic regression model
    model_above_only <- lme4::glmer(event_next_n ~ titre_above_threshold 
                                    +age_grp 
                                    +sex
                                    +hhsize
                                    + (1 | pid) + (1 | hid), data = final_df_5, family = binomial)
    
    tb1 <- model_above_only %>%
        tbl_regression(exponentiate = TRUE,
                       pvalue_fun = format_p_value) %>%
        bold_p(t = 0.05) %>%
        modify_header(list(label ~ paste(sample, class, "titre"))) 
    
    
    # Extract the summary of the model
    model_summary <- summary(model_above_only)
    
    # Calculate confidence intervals
    conf_intervals <- confint(model_above_only, parm = "beta_", method = "Wald")  # Wald CIs 
    
    # Create a dataframe with hazard ratios (exponentiated coefficients), confidence intervals, and p-values
    or_df <- data.frame(
        term = rownames(model_summary$coefficients),  # Variable names
        estimate = exp(model_summary$coefficients[, "Estimate"]),  # OR (exp(coef))
        conf.low = exp(conf_intervals[, 1]),  # Lower bound of 95% CI (exponentiated)
        conf.high = exp(conf_intervals[, 2]),  # Upper bound of 95% CI (exponentiated)
        p.value = model_summary$coefficients[, "Pr(>|z|)"]  # P-values
    ) %>%
        mutate(
            # Add descriptive labels for the variables
            variable_label = case_when(
                grepl("titre_above_threshold", term) ~ "Titre",
                grepl("age_grp", term) ~ "Age group",
                grepl("hhsize", term) ~ "Household size",
                grepl("sexFemale", term) ~ "Sex",
                TRUE ~ NA_character_  # Default to NA if no match
            ),
            # Simplify term names for clarity
            term = case_when(
                grepl("titre_above_threshold", term) ~ "Titre",
                grepl("age_grp< 2 years", term) ~ "< 2 years",
                grepl("age_grp12-18 years", term) ~ "12-18 years",
                grepl("age_grp2-4 years", term) ~ "2-4 years",
                grepl("age_grp5-11 years", term) ~ "5-11 years",
                grepl("sexFemale", term) ~ "Female",
                grepl("hhsize", term) ~ "Household size",
                TRUE ~ term  # Retain term if it doesn't match any pattern
            )
        )
    
    # Add an antigen column to the dataframe
    or_df$antigen <- antigen
    
    
    results_list <- list("tbl" = tb1, "or_df" = or_df)
    
    return(results_list)
}


tbl_list <- list()
df_list <- list()
adjusted_model_output <- tibble()

# Loop through each Anitgen and apply the function: glmer_adjusted_piecewise

for (ag in Antigen_list) {
    
    results <- glmer_adjusted_piecewise(path_to_titre.df = "data/blood_IgG_titres.RDS",
                                        sample = "Blood",
                                        class = "IgG",
                                        var_name = "titre",
                                        antigen = ag)
    
    
    
    # Extract and store the table and plot
    tbl_list[[ag]] <- results$tbl
    df_list[[ag]] <- results$or_df
    
    
    adjusted_model_output <- bind_rows(adjusted_model_output,  results$or_df)
}


adjusted_model_output


# Combine the gtsummary tables
combined_table <- tbl_merge(tbl_list,
                            tab_spanner = names(tbl_list))

combined_table



#######################################################
###### Plot forest plots ############

# create reference points in dataframe for subsequent forest plots
new_rows <- unique(adjusted_model_output$antigen) %>%
    expand.grid(antigen = ., term = c("Over 18 years (ref)", "Male (ref)")) %>%
    mutate(
        estimate = NA, 
        conf.low = NA, 
        conf.high = NA,
        p.value = NA,
        variable_label = case_when(
            term == "Below threshold (ref)" ~ "Relation to titre threshold",
            term == "Over 18 years (ref)" ~ "Age group",
            term == "Male (ref)" ~ "Sex"
        )
    )


# 
final_fp_df<- bind_rows(new_rows,adjusted_model_output)

levels(final_fp_df$term) <- gsub("hhsize", "Household size", levels(final_fp_df$term))


final_fp_df <-final_fp_df %>%
    mutate(
        term = factor(term, levels = c(
            "Titre", 
            "Over 18 years (ref)",
            "< 2 years",
            "2-4 years", 
            "5-11 years", 
            "12-18 years", 
            "Male (ref)",
            "Female",
            "Household size" 
        ))
    )

final_fp_df  %>%
    filter(term != "(Intercept)") %>% # Exclude the intercept
    ggplot(aes(x = estimate, y = factor(term, levels = c(
        "Titre", 
        "Over 18 years (ref)",
        "< 2 years",
        "2-4 years", 
        "5-11 years", 
        "12-18 years", 
        "Male (ref)",
        "Female",
        "Household size"
    )))) +
    geom_point() +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
    #  facet_wrap(~ antigen, scales = "free") + # Facet by antigen
    labs(
        title = "",
        x = "Odds ratios with 95% CIs",
        y = ""
    ) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    # Add a point for the reference level at x = 1, y = position of reference variable
    
    theme_minimal() +
    geom_point(data = new_rows %>% filter(!antigen %in% c("DNAseB", "GAC")),
               aes(x = 1, y = term), color = "blue", size = 3) +
    # Separate age_grp and event_type under subheadings on y-axis, including reference level
    facet_grid(cols = vars(antigen), rows = vars(variable_label), scales = "free", space = "free_y") +
    scale_x_log10() + 
    theme(
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        #panel.background = element_blank(),
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(size = 18, hjust = 0.5),
        strip.text.x = element_text(size = 18),
        strip.text.y = element_text(size = 18, angle = 0), 
        legend.title=element_text(size=16),
        legend.text=element_text(size=16),
        axis.text.x =element_text(size=12, angle = 90),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18))


# Ensure 'term' and 'variable_label' are ordered as intended
final_fp_df <- final_fp_df %>%
    mutate(
        term = factor(term, levels = c(
            "Titre", "Over 18 years (ref)", 
            "12-18 years", "5-11 years", "2-4 years", "< 2 years", 
            "Male (ref)", "Female", "Household size"
        )),
        variable_label = forcats::fct_relevel(variable_label, 
                                              "Titre", "Age group", "Sex", "Household size")
    )

# Plot

main_07_fig4B_V3 <- final_fp_df %>%
    filter(term != "(Intercept)") %>% # Exclude the intercept
    ggplot(aes(x = estimate, y = term)) +  # term already defined with levels
    geom_point() +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
    labs(
        title = "",
        x = "Odds ratios with 95% CIs",
        y = ""
    ) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    theme_minimal() +
    geom_point(
        data = new_rows %>% filter(!antigen %in% c("DNAseB", "GAC")),
        aes(x = 1, y = term), color = "blue", size = 1.5
    ) +
    facet_grid(
        cols = vars(antigen),
        rows = vars(forcats::fct_relevel(variable_label, 
                                         "Titre", "Age group", "Sex", "Household size")), 
        scales = "free", 
        space = "free_y"
    ) +
    scale_x_log10() +  # Log scale for odds ratios
    theme(
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size = 18, hjust = 0.5),
        strip.text.x = element_text(size = 18),
        strip.text.y = element_text(size = 18, angle = 0),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.text.x = element_text(size = 12, angle = 90),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18)
    ) +
    theme_universal(base_size = plot_basesize)

main_07_fig4B_V3


##############################################################################
#### Describe baseline titres in relation to putative 50% protective thresholds 
##############################################################################


# run this line to source annoymised publically available data 
#df <- readRDS("data/baseline_blood_no_disease_titres.RDS")  

# run this line to source original data -> not annonymised only for manuscript revisions 
df <- readRDS("R_objects/baseline_blood_no_disease_titres.RDS")
protective_threshold <- readRDS("data/bloodIgG_protective_threshold_df.RDS")



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
    "At baseline, antibody levels from the cohort showed that %d individuals (%s) for SLO, %d (%s) for SpyAD, and %d (%s) for SpyCEP had IgG levels above the protective thresholds (Figure 4E).",
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
    geom_point(alpha = 0.5, size = dot_size) +
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






################################################
########## Cox AG survival analysis ############
################################################

# coxag_cont

# Dynamically fits a Cox proportional hazards model and returns a formatted regression table.
# Inputs: outcome and predictor (as expressions), optional covariates, a label for the predictor, and the dataset (survival dataframe).
# Process:
#   1. Constructs a survival model formula using tstart, tstop, and the provided outcome, predictor, and covariates.
#   2. Fits a Cox model (with clustering by 'hid' and participant id 'pid').
#   3. Formats the model output into a regression table with exponentiated coefficients (hazard ratios), p-values, and event counts.


# Output: A regression table summarizing the effect of the predictor on the specified outcome.

coxag_cont <- function(outcome, predictor, covariates = NULL, label = as.character(), data) {
    
    outcome <- substitute(outcome)
    predictor <- substitute(predictor)
    
    # Create the formula for the model dynamically
    covariate_formula <- if (!is.null(covariates) && length(covariates) > 0) {
        paste(covariates, collapse = " + ")
    } else {
        ""
    }
    
    full_formula <- as.formula(paste("Surv(tstart, tstop, eval(outcome, data)) ~ eval(predictor, data)", 
                                     covariate_formula, sep = ifelse(covariate_formula != "", " + ", "")))
    
    model <- coxph(full_formula, cluster = hid, id = pid, data = data)
    
    tb1 <- model %>%
        tbl_regression(exponentiate = TRUE,
                       label = list("eval(predictor, data)" ~ label),
                       include = "eval(predictor, data)",
                       pvalue_fun = format_p_value) %>%
        add_nevent() %>%
        bold_labels() %>%
        bold_p(t = 0.05) %>%
        modify_header(list(label ~ "**Antibody titre**"))
    
    return(tb1)
}

# draw_coxag_cont_tables_cont_piecewise: Generates merged Cox regression tables for different event outcomes.
# Inputs: Data frame, sample, class, event_breakdown flag, and optional covariates.

# Process:
#   1. Defines a helper to create a caption based on sample, class, and covariate adjustments.
#   2. Runs piecewise Cox models (via coxag_cont) for multiple outcomes (e.g., all events, disease, carriage, etc.).
#   3. Stacks and merges the resulting tables into one composite table.
#   4. Modifies the caption and returns the final merged table.

# Output: A table of hazard ratios, confidence intervals, and p-values for the specified outcomes.



draw_coxag_cont_tables_cont_piecewise <- function(df, sample, class, event_breakdown = F, covariates = NULL) {
    
    
    
    # Helper function to create the caption string
    create_caption <- function(sample, class, covariates) {
        covariate_str <- if (!is.null(covariates) && length(covariates) > 0) {
            paste(", adjusted for", paste(covariates, collapse = ", "),".")
        } else {
            ", unadjusted."
        }
        paste("Table demonstrating hazards ratios of incidence of Strep A events by titre of", sample, class, covariate_str)
    }
    
    
    # Run the models and store the output tables in a list, ordered by DNAseB, GAC, SLO, SpyAD, SpyCEP
    table_list1 <- list(
        
        coxag_cont(gas_event, slo_above_threshold, "SLO", df, covariates = covariates),
        coxag_cont(gas_event, spyad_above_threshold, "SpyAD", df, covariates = covariates),
        coxag_cont(gas_event, spycep_above_threshold, "SpyCEP", df, covariates = covariates)
    )
    
    # Combine the tables into one
    merged_table1 <- tbl_stack(table_list1)
    
    # Repeat for other tables, keeping the same order
    table_list2 <- list(
        
        coxag_cont(gas_infection, slo_above_threshold, "SLO", df, covariates = covariates),
        coxag_cont(gas_infection, spyad_above_threshold, "SpyAD", df, covariates = covariates),
        coxag_cont(gas_infection, spycep_above_threshold, "SpyCEP", df, covariates = covariates)
    )
    
    merged_table2 <- tbl_stack(table_list2)
    
    table_list3 <- list(
        coxag_cont(gas_carriage, slo_above_threshold, "SLO", df, covariates = covariates),
        coxag_cont(gas_carriage, spyad_above_threshold, "SpyAD", df, covariates = covariates),
        coxag_cont(gas_carriage, spycep_above_threshold, "SpyCEP", df, covariates = covariates)
    )
    
    merged_table3 <- tbl_stack(table_list3)
    
    table_list4 <- list(
        coxag_cont(gas_pharyngitis, slo_above_threshold, "SLO", df, covariates = covariates),
        coxag_cont(gas_pharyngitis, spyad_above_threshold, "SpyAD", df, covariates = covariates),
        coxag_cont(gas_pharyngitis, spycep_above_threshold, "SpyCEP", df, covariates = covariates)
    )
    
    merged_table4 <- tbl_stack(table_list4)
    
    table_list5 <- list(
        coxag_cont(gas_pyoderma, slo_above_threshold, "SLO", df, covariates = covariates),
        coxag_cont(gas_pyoderma, spyad_above_threshold, "SpyAD", df, covariates = covariates),
        coxag_cont(gas_pyoderma, spycep_above_threshold, "SpyCEP", df, covariates = covariates)
    )
    
    merged_table5 <- tbl_stack(table_list5)
    
    table_list6 <- list(
        coxag_cont(gas_throat_carriage, slo_above_threshold, "SLO", df, covariates = covariates),
        coxag_cont(gas_throat_carriage, spyad_above_threshold, "SpyAD", df, covariates = covariates),
        coxag_cont(gas_throat_carriage, spycep_above_threshold, "SpyCEP", df, covariates = covariates)
    )
    
    merged_table6 <- tbl_stack(table_list6)
    
    table_list7 <- list(
        coxag_cont(gas_skin_carriage, slo_above_threshold, "SLO", df, covariates = covariates),
        coxag_cont(gas_skin_carriage, spyad_above_threshold, "SpyAD", df, covariates = covariates),
        coxag_cont(gas_skin_carriage, spycep_above_threshold, "SpyCEP", df, covariates = covariates)
    )
    
    merged_table7 <- tbl_stack(table_list7)
    
    
    
    # Event breakdown based on flags
    
    if (event_breakdown) {
        
        table_2 <- tbl_merge(list(merged_table4, merged_table5, merged_table6, merged_table7),
                             tab_spanner = c("Pharyngitis", "Pyoderma", "Throat carriage", "Skin carriage")) %>%
            modify_caption(create_caption(sample, class, covariates))
        
        return(table_2)
        
    }  else {
        
        table_1 <- tbl_merge(list(merged_table1, merged_table2, merged_table3),
                             tab_spanner = c("All events", "Disease events", "Carriage events")) %>%
            modify_caption(create_caption(sample, class, covariates))
        
        return(table_1)
    }
}


#### Import survival dataframe. 

# A detailed description of the generation of this dataframe can be found at https://edwinarmitage.github.io/SpyCATS_primary_analysis.html

survivial_df_blood_IgG <-  readRDS("data/survivial_df_blood_IgG.RDS")

cox_PH1 <- draw_coxag_cont_tables_cont_piecewise(survivial_df_blood_IgG, sample = "blood", class = "IgG", covariates = c("sex", "age_grp", "hhsize"))
cox_PH2 <- draw_coxag_cont_tables_cont_piecewise(survivial_df_blood_IgG, sample = "blood", class = "IgG", event_breakdown = T, covariates = c("sex", "age_grp", "hhsize"))

cox_PH1
cox_PH2

doc <- read_docx()
# Add regression table to word document
doc <- add_table_to_doc(doc, cox_PH1)
doc <- add_table_to_doc(doc, cox_PH2)

print(doc, target = "R_output/CoxPH_supplemen_tables_V3.docx")

