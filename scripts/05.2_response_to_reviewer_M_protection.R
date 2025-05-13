# Title: Respondin to reviewer: unrelated M protein and protection 
# Version: 1.0
# Date: 2025-03-24
# Author: Alexander J Keeley

# Inputs: 
#   - data/M_protection_dataframe.RDS
#   - data/incidence_start_dates.RDS
#   - data/SpyCATS_incidence_df.RDS
#   - data/all_events_long_incidence_wgs.RData

# Outputs: 
#   - Logistic regression models assessing Z-scored IgG titres (and unrelated comparisons) for cluster relatied anti-M and protection
#   - Combined figure comparing event proportions and logistic predictions
#   - AIC values for model comparison

# Description:
#
# This script models the association between M-protein-specific IgG titres (expressed as Z-scores) and the short-term risk 
# of a culture confirmed event within 45 days. It compares homologous responses to unrelated antigens using mixed-effects logistic regression.
#


# Setup environment
source("scripts/setup_environment.R")
source("scripts/load_functions.R")

# load datframes
M_protection_df <- readRDS("data/M_protection_dataframe.RDS")
start_dates <- readRDS("data/incidence_start_dates.RDS")

# A function to add the unrelated M titre and compare AICs of the various models when exploring unrelated titre vs protection

create_regression_dataframe_M_unrelated<- function(path_to_titre.df,sample,class, next_event_window = 45, var_name = "titre", n_tile = 4, antigen, df = 1,slice_window = 0.1, protective_threshold = 0.67) {
    
    
    
    # Read titre data and merge with start dates and age, adjust visit dates, and categorize titres
    fun_titres <- path_to_titre.df %>%
        left_join(start_dates) %>%
        left_join(age) %>%
        mutate(visit_date = as.numeric(ymd(as.character(visit_date))) - entry_1)
    
    
    #import incidence dataframe
    pos_incidence_zero <- readRDS("data/SpyCATS_incidence_df.RDS")
    
    # Filter for specific antigen
    antigen_df <- fun_titres 
    
    # Prepare event data
    fun_df <- pos_incidence_zero %>%
        select(pid,date, gas_event)
    
    # Calculate time until next event and categorize events
    fun_df <- fun_df %>%
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
    
    
    # Create buckets for titre
    bucket_size <- 0.2   
    
    final_df_3  <-fun_df %>%
        left_join(antigen_df %>% 
                      select(pid,date = visit_date,titre, unrelated_comparison, age)) %>%
        mutate(titre_bucket = floor(titre / bucket_size) * bucket_size)  %>%
        group_by(pid) %>%
        #   fill(titre) %>%
        ungroup() %>%
        mutate(hid = substring(pid,1,3)) %>%
        ## This is key ##- to deduplicate identiacal rows of data which articifically raise the power of subsequent models 
        unique()
    
    
    
    # Remove rows with NA in titre or event_next_day
    final_df_4 <- na.omit(final_df_3, cols = c("titre", "event_next_day"))
    
    # Calculate and adjust thresholds for analysis
    min_titre <- min(final_df_4$titre, na.rm = TRUE)
    max_titre <- max(final_df_4$titre, na.rm = TRUE)
    
    
    # Round down/up to the nearest multiple of slice_window
    rounded_min <- floor(min_titre / slice_window) * slice_window
    rounded_max <- ceiling(max_titre / slice_window) * slice_window
    
    # Adjust if necessary to align with multiples of slice window
    rounded_min_multiplier <- 1 / slice_window
    rounded_max_multiplier <- 1 / slice_window
    
    if (rounded_min * rounded_min_multiplier %% 1 != 0) {
        rounded_min <- rounded_min - (rounded_min %% slice_window)
    }
    if (rounded_max * rounded_max_multiplier %% 1 != 0) {
        rounded_max <- rounded_max + (slice_window - (rounded_max %% slice_window))
    }
    
    thresholds <- seq(from = rounded_min, to = rounded_max, by = slice_window)
    
    
    # Initialize a data frame to store results
    probabilities <- data.frame(Threshold = numeric(), Probability = numeric(), percent = numeric())
    
    for(threshold in thresholds) {
        
        # Subset the data for rows where titre >= threshold
        subset_df <- final_df_4 %>%
            filter(titre >= threshold)
        
        # Check if the subset is not empty
        if(nrow(subset_df) > 0) {
            # Calculate the probability
            prob <- sum(subset_df$event_next_n == 1) / nrow(subset_df)
        } else {
            # Set probability to NA if no rows meet the threshold
            prob <- NA
        }
        
        n_obs <- nrow(subset_df) / nrow(final_df_4) * 100
        
        # Store the result
        probabilities <- rbind(probabilities, data.frame(Threshold = threshold, Probability = prob, percent= n_obs))
    }
    
    
    # Find the maximum threshold where Number_obs is at least 97.5%
    max_threshold <- max(probabilities$Threshold[probabilities$percent >= 97.5])
    
    
    # then do logistic regression: 
    
    library(splines)
    
    stopifnot("pid" %in% names(final_df_3))
    
    model <- lme4::glmer(event_next_n ~ titre + (1 | pid) + (1 | hid), data=final_df_3, family=binomial)
    
    print(paste("AIC titre:",AIC(model)))
    
    tb1 <- model %>%
        tbl_regression(exponentiate = TRUE) %>%
        add_nevent() %>%
        bold_p(t = 0.05) %>%
        modify_header(list(label ~ paste(sample, class, "titre"))) # Modify header here 
    
    # Update the row labels with the antigen name
    tb1$table_body$label <- paste0(antigen)
    
    # Extract AIC
    model_aic <- AIC(model)
    print(paste(antigen,":", AIC(model)))
    
    
    
    model2 <- lme4::glmer(event_next_n ~ titre + unrelated_comparison + (1 | pid) + (1 | hid), data=final_df_3, family=binomial)
    
    print(paste("AIC titre + unrelated titre:",AIC(model2)))
    
    tb2 <- model2 %>%
        tbl_regression(exponentiate = TRUE) %>%
        add_nevent() %>%
        bold_p(t = 0.05) %>%
        modify_header(list(label ~ paste(sample, class, "titre"))) # Modify header here 
    
    model3 <- lme4::glmer(event_next_n ~ unrelated_comparison + (1 | pid) + (1 | hid), data=final_df_3, family=binomial)
    
    print(paste("AIC  unrelated titre:",AIC(model3)))
    
    tb3 <- model3 %>%
        tbl_regression(exponentiate = TRUE) %>%
        add_nevent() %>%
        bold_p(t = 0.05) %>%
        modify_header(list(label ~ paste(sample, class, "titre"))) # Modify header here 
    
    
    
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
    
    
    # Create a sequence of titre values
    new_data <- data.frame(titre = seq(min(final_df_3$titre, na.rm=TRUE), 
                                       max(final_df_3$titre, na.rm=TRUE), length.out=100))
    
    # Predict probabilities
    new_data$prob <- predict(model, newdata=new_data, type="response", re.form=NA)
    
    # Calculate confidence intervals
    conf_int <- predict(model, newdata=new_data, type="link",re.form=NA, se.fit=TRUE)
    new_data$lower <- plogis(conf_int$fit - 1.96 * conf_int$se.fit)
    new_data$upper <- plogis(conf_int$fit + 1.96 * conf_int$se.fit)
    
    
    
    
    
    # Step 3: Modify Return Value
    
    
    
    ########## plot it 
    
    options(digits = 3, scipen = 9)
    
    # Define a function to convert log10 values to their exponential form
    log10_to_exp <- function(x) {
        10^x
    }
    
    
    plot2 <- ggplot(new_data, aes(x=titre, y=prob)) +
        geom_col(data =  probabilities, aes(x = Threshold, y = Probability, alpha = percent), fill = "#138aa8") +
        #   geom_smooth(data =  probabilities, aes(x = Threshold, y = Probability, fill = percent), col = "#f5bff0") +
        geom_line(color="blue") +
        geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
        
        #scale_x_continuous(limits = c(max_threshold, max(probabilities$Threshold)),labels = trans_format("identity", log10_to_exp)) +
        #    ylim(0,0.1) +
        # xlim(max_threshold, max(thresholds)) +
        labs(title= paste(sample, class, antigen),
             x=paste0("Titre (Z score) to cluster homologous event"),
             y=paste0("Proportion experiencing event within ", next_event_window, " days")) +
        guides(fill = "none",
               alpha = "none")+
        theme_minimal()      +
        annotate("text", x = max_threshold + 0.5, y = 0.15, 
                 label = paste0("OR: ", round(ORs[2], 2), 
                                "\n95% CI: [", round(CIs[1, 1], 2), ", ", round(CIs[1, 2], 2), "]", 
                                "\n", p_label), 
                 hjust = 0) +
        theme_universal(base_size = plot_basesize)
    
    plot1 <- ggplot(data = probabilities, aes(x = Threshold, y = Probability)) +
        geom_col(fill = "#138aa8") +
        labs(
            x=paste0("IgG (Z score) to cluster homologous event"),
            y = paste0("Proportion with event within ", next_event_window, " days")) +
        guides(fill = "none", alpha = "none") +
        theme_minimal() +
        theme_universal(base_size = plot_basesize)
    
    
    
    
    plot2 <- ggplot(new_data, aes(x = titre, y = prob)) +
        geom_line(color = "blue") +
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
        labs(
            x=paste0("IgG (Z score) to cluster homologous event"),
            y = paste0("Probability of event in ", next_event_window, " days")) +
        guides(fill = "none", alpha = "none") +
        theme_minimal() +
        annotate("text", x = max_threshold + 1, y = 0.25, 
                 label = paste0("OR: ", round(ORs[2], 2), 
                                "\n95% CI: [", round(CIs[1, 1], 2), ", ", round(CIs[1, 2], 2), "]", 
                                "\n", p_label), 
                 hjust = 0) +
        theme_universal(base_size = plot_basesize)
    
    plot <- cowplot::plot_grid(plot1,plot2,ncol = 1)
    
    
    # Instead of printing, add the table and plot to a list
    results_list <- list("table" = tb1,"table2" = tb2,"table3" = tb3, "plot" = plot, "AIC_glmer" = model_aic)
    
    
    results_list$probabilities2 <- new_data
    
    
    
    
    # Return the list
    
    
    return(results_list)
    
    
}

load("data/all_events_long_incidence_wgs.RData")

AIC_glmer <- data.frame(Antigen = character(), AIC_glmer = numeric(), stringsAsFactors = FALSE)

result <- create_regression_dataframe_M_unrelated(path_to_titre.df =  M_protection_df,
                                                  sample = "Blood",
                                                  class = "IgG",
                                                  var_name = "titre",
                                                  antigen = "to cluster homologous M peptide by Z score",
                                                  slice_window = 0.1)



plot_09_fig06_panelE <- result$plot
plot_09_fig06_panelE 


AIC_glmer <- rbind(AIC_glmer, data.frame(Antigen = "M", AIC_glmer = result$AIC_glmer))

AIC_glmer

result$table
result$table2
result$table3

