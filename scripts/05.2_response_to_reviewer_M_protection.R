




create_regression_dataframe_M_unrelated<- function(path_to_titre.df,sample,class, next_event_window = 45, var_name = "titre", n_tile = 4, antigen, df = 1,slice_window = 0.1, protective_threshold = 0.67) {
    
    
    
    # Read titre data and merge with start dates and age, adjust visit dates, and categorize titres
    fun_titres <- path_to_titre.df %>%
        left_join(start_dates) %>%
        left_join(age) %>%
        mutate(visit_date = as.numeric(ymd(as.character(visit_date))) - entry_1)
    
    
    #start and stop dates from f/u periods. Set all start dates at zero.
    # Adjust follow-up periods, setting start dates to zero and calculating periods and total participant years
    
    incidence_zero <- follow_up_dates_incidence %>%
        select(!entry_6:exit_7) %>%
        mutate(across(!pid, ~ .x - entry_1)) %>%   # for each row - take the entry 1 date away from the date variable (except pid )
        filter(!exit_1 == 0) %>% # removes the participant who had no follow up time.
        rowwise() %>%
        mutate(entry_1 = entry_1,
               period_1 = (exit_1 - entry_1),
               period_2 = (exit_2 - entry_2),
               period_3 = (exit_3 - entry_3),
               period_4 = (exit_4 - entry_4),
               period_5 = (exit_5 - entry_5),
               total_pyears = sum(c_across(period_1:period_5), na.rm = T) / 365.25,
               start = min(c_across(entry_1:entry_5), na.rm = T),
               stop = max(c_across(exit_1:exit_5), na.rm = T)) %>%
        select(-contains("period"))
    
    # Pivot entry and exit dates to long and calculate "gap" status - whether there was a gap in follow up.
    incidence_entry_exit_long_zero <- incidence_zero %>%
        pivot_longer(entry_1:exit_5, names_to = "date", values_to = c("timepoint"), values_drop_na = T) %>% 
        mutate(
            gap_status = case_when(
                grepl("entry_",date) ~ 0,
                grepl("exit_",date) ~ 1)) %>%
        select(-total_pyears)
    
    # Merge demographic data
    
    incidence_zero <- right_join(final_dems,incidence_zero)
    
    
    
    # Create df for positive events incidence and zero dates (to be dependent variables)
    pos_incidence_zero <- left_join(all_events_long_incidence_wgs,final_dems)
    pos_incidence_zero <- left_join(pos_incidence_zero,start_dates)
    pos_incidence_zero <- pos_incidence_zero %>%
        mutate(date = as.numeric(date - entry_1)) %>%
        filter(date >= 0)
    pos_incidence_zero <- pos_incidence_zero %>%
        mutate(across(contains("gas") | contains("sdse"), ~ ifelse(date == 0 & . == 1, NA, .)))  # This removes the event if it occureed in the first visit 
    
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

load("R_objects/all_events_long_incidence_wgs.RData")
#load("../../../SpyCATS_github/R outputs/all_events_long_incidence_pcr.RData")

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
