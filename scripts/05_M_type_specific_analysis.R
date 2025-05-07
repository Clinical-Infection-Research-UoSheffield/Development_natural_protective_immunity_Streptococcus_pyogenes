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

# Import baseline dataframes

df_titre <- readRDS("R_objects/M_CRC_baseline_no_disease.RDS")

n1 <- df_titre %>% pull(pid) %>% unique() %>% length()

sprintf(
    "Basline titres to cluster representative M peptide hyper variable regions by age in SpyCATS study baseline (n = %d)", n1
)



# create an empty dataframe to append the output of your for loop
centile_df <- data.frame()


# Create an empty list to store plots
plots <- list()
plots2 <- list()
plots3 <- list()

# Add emm cluster as a variable 

cluster_rep.df <- readRDS("R_objects/cluster_rep.df.RDS")

df_titre2 <- df_titre %>%
    left_join(cluster_rep.df) %>%
    mutate(
        Antigen_cluster = 
            case_when(Antigen %in% c("M6","M74", "M55", "M18") ~ "Singleton",
                      T ~ Antigen_cluster)
    )

# loop through each antigen, fit polynomial model to the data, and plot age vs antibody titre

for (a in unique(df_titre$Antigen)) {        
    
    
    # create a datframe for that antigen
    ctrl <- df_titre2 %>%
        filter(Antigen == a) %>%
        arrange(age)                                
    
    cluster_label <- as.character(ctrl$Antigen_cluster[1])
    
    common_y_limits <- c(min(df_titre$titre), max(df_titre$titre))
    
    
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
        scale_color_manual(values = StrepA_colscheme) +
        guides(color = "none") +
        geom_point(alpha = 0.5) +
        geom_line(aes(y = .fitted), linetype = "dashed", color = "black") +
        geom_line(aes(y = upp25), linetype = "dashed", color = "red", alpha = 0.5) +
        geom_line(aes(y = upp80), linetype = "dashed", color = "blue") +
        geom_line(aes(y = upp975), linetype = "dashed", color = "red", alpha = 0.5) +
        labs(title = paste0(a," \n(",cluster_label,")"), x = "Age", y = "IgG level (RLU/mL)") +
        theme_minimal() +
        scale_y_continuous(limits =common_y_limits, breaks = c(1:6), labels = log10_to_exp) +
        theme_universal(base_size = plot_basesize)
    
    # Store the plot in the list with the antigen name as the key
    
    plots[[a]] <- p 
    
    
    
}


# Arrange the plots into a grid layout

plot_09_fig06_panelA <- do.call(grid.arrange, c(plots, ncol = 7))
print(plot_09_fig06_panelA)


########## plot 2 ####################

#### Baseline raw MFI ################


# Use this when annonymised data loaded: 
# MFI <- readRDS("data/M_CRC_baseline_MFI_no_disease.RDS")

MFI <- readRDS("R_objects/M_CRC_baseline_MFI_no_disease.RDS")

plot_09_fig06_panelB <- MFI %>%
    # Reorder Antigen based on the mean MFI for each group, with highest first
    group_by(Antigen) %>%
    mutate(mean_MFI = mean(MFI)) %>%
    ungroup() %>%
    mutate(Antigen = fct_reorder(Antigen, mean_MFI, .desc = TRUE)) %>%
    filter(!Antigen %in% c("J8","P17", "K4S2")) %>%
    ggplot(aes(x = Antigen, y = MFI, col = Antigen)) +
    scale_color_manual(values = StrepA_colscheme) +
    geom_violin() +
    geom_jitter(alpha = 0.1) +
    guides(col = "none") +
    labs(
        y = "log10 MFI * sample dilution",
        #  title = "IgG Antibodies as raw median flourescene from from baseline samples demonstrating relative abundance of specific IgG in pariticpant sera"
    ) +
    #  stat_compare_means() +
    theme_minimal() +
    theme_universal(base_size = plot_basesize)


plot_09_fig06_panelB


##############################################
#### Changes around events #####
##############################################

# Import a dataframe of events with their emm type, cluster and titres

##### each row is a set of measurements for each antigen arounda n event.

# Import the antiM, absolute changes in Z score around events data


cleaned_data <- readRDS("R_objects/M_CRC_E3_around_events.RDS") %>%
    group_by(event_id) %>%
    filter(n() >1) # This removes two readings where only a single antigen has titre data - not true values 

#### Number of events 

cleaned_data %>%
    pull(event_id) %>%
    unique() %>%
    length()

# number of paired measurements
cleaned_data %>%
    dim()

# number of paired measurements by category
cleaned_data %>%
    group_by(relation_to_event) %>%
    summarise(n = n())

# Convert relation_to_event to a factor with the specified order
cleaned_data <- cleaned_data %>%
    mutate(relation_to_event = factor(relation_to_event, 
                                      levels = c("homologous", "cluster_homologous", "unrelated")))

# Perform Dunn's test
dunn_test_results <- dunnTest(pre_post ~ relation_to_event, data = cleaned_data, method = "bonferroni")

# Extract the p-values and comparisons
dunn_results_df <- as.data.frame(dunn_test_results$res)
dunn_results_df <- dunn_results_df %>%
    mutate(comparison = gsub(" ", " vs ", Comparison),
           p.adj = P.adj) %>%
    separate(Comparison, into = c("group1", "group2"), sep = " - ") %>%
    mutate(p.adj = P.adj) %>%
    mutate(significance = sapply(p.adj, function(p) {
        if (p < 0.0001) {
            "****"
        } else if (p < 0.001) {
            "***"
        } else if (p < 0.01) {
            "**"
        } else if (p < 0.05) {
            "*"
        } else {
            "ns"  # Non-significant
        }
    }))

dunn_results_df

# Create the plot
plot <- cleaned_data %>%
    ggplot(aes(x = relation_to_event, y = pre_post)) +
    geom_boxplot(width = 0.05, outlier.size = 0.3) +
    ggdist::geom_dots(
        width = 0.3,
        position = position_nudge(x = 0.05), 
        aes(col = relation_to_event, fill = relation_to_event)
    ) +
    ggdist::stat_slab(
        side = 'bottom', 
        position = position_nudge(x = -0.05), 
        width = 0.3, 
        aes(col = relation_to_event, fill = relation_to_event)
    ) +
    labs(
        x = "Relationship of antibody to emm type of event",
        y = "Absolute change in titre Z score around event"
    ) +
    scale_x_discrete(labels = function(x) str_replace_all(x, "_", " ")) + # Remove underscores from labels
    scale_color_brewer() +
    scale_fill_brewer() +
    guides(col = "none", fill = "none") +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5)  
    ) +
    theme_universal(base_size = plot_basesize)

# Add Dunn's test results to the plot
plot_09_fig06_panelC <- plot + ggpubr::stat_pvalue_manual(dunn_results_df, 
                                                          label = "significance", 
                                                          y.position = max(cleaned_data$pre_post, na.rm = TRUE) + 0.5,
                                                          tip.length = 0.01, 
                                                          step.increase = 0.1,
                                                          size = 8)

# Display the final plot
plot_09_fig06_panelC


############################################################
############### sensitivity analyses #######################
############################################################

# option 2: keeping the best data in the dataframe: 


##### each row is a set of measurements for each antigen arounda n event.

cleaned_data3 <- readRDS("R_objects/M_CRC_E3_around_events.RDS") %>%
    filter(Antigen %in% c("M1", "M4", "M44", "M75", "M89", "M25", "M87")) %>%
    group_by(event_id) %>%
    filter(n() >1)

#### Number of events 

cleaned_data3 %>%
    pull(event_id) %>%
    unique() %>%
    length()

# number of paired measurements
cleaned_data3 %>%
    dim()

# number of paired measurements by category
cleaned_data3 %>%
    group_by(relation_to_event) %>%
    summarise(n = n())

# Convert relation_to_event to a factor with the specified order
cleaned_data3 <- cleaned_data3 %>%
    mutate(relation_to_event = factor(relation_to_event, 
                                      levels = c("homologous", "cluster_homologous", "unrelated")))

# Perform Dunn's test
dunn_test_results <- dunnTest(pre_post ~ relation_to_event, data = cleaned_data3, method = "bonferroni")

# Extract the p-values and comparisons
dunn_results_df <- as.data.frame(dunn_test_results$res)
dunn_results_df <- dunn_results_df %>%
    mutate(comparison = gsub(" ", " vs ", Comparison),
           p.adj = P.adj) %>%
    separate(Comparison, into = c("group1", "group2"), sep = " - ") %>%
    mutate(p.adj = P.adj) %>%
    mutate(significance = sapply(p.adj, function(p) {
        if (p < 0.0001) {
            "****"
        } else if (p < 0.001) {
            "***"
        } else if (p < 0.01) {
            "**"
        } else if (p < 0.05) {
            "*"
        } else {
            "ns"  # Non-significant
        }
    }))

dunn_results_df

# Create the plot
plot <- cleaned_data3 %>%
    ggplot(aes(x = relation_to_event, y = pre_post)) +
    geom_boxplot(width = 0.05, outlier.size = 0.3) +
    ggdist::geom_dots(
        width = 0.3,
        position = position_nudge(x = 0.05), 
        aes(col = relation_to_event, fill = relation_to_event)
    ) +
    ggdist::stat_slab(
        side = 'bottom', 
        position = position_nudge(x = -0.05), 
        width = 0.3, 
        aes(col = relation_to_event, fill = relation_to_event)
    ) +
    labs(
        x = "Relationship of antibody to emm type of event",
        y = "Absolute change in titre Z score around event"
    ) +
    scale_x_discrete(labels = function(x) str_replace_all(x, "_", " ")) + # Remove underscores from labels
    scale_color_brewer() +
    scale_fill_brewer() +
    guides(col = "none", fill = "none") +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5)  
    ) +
    theme_universal(base_size = plot_basesize)

# Add Dunn's test results to the plot
plot_09_fig06_panelC_sensitivity_2 <- plot + ggpubr::stat_pvalue_manual(dunn_results_df, 
                                                                        label = "significance", 
                                                                        y.position = max(cleaned_data3$pre_post, na.rm = TRUE) + 0.5,
                                                                        tip.length = 0.01, 
                                                                        step.increase = 0.1,
                                                                        size = 8)

# Display the final plot
plot_09_fig06_panelC_sensitivity_2



####

############################################################
############### Protection analyses ########################
############################################################


M_protection_df <- readRDS("R_objects/M_protection_dataframe.RDS")

#### The number of emm-typed events
n1 <- M_protection_df %>%
    select(event_id) %>%
    unique() %>%
    dim()

# number of measurements in this dataframe
n2 <- M_protection_df %>%    
    select(pid, visit_date, Antigen, titre) %>%
    unique() %>%
    dim()

# number of timepoints 
n3 <- M_protection_df %>%
    select(pid, visit_date) %>%
    unique() %>%
    dim()

# number of individuals 
n4 <- M_protection_df %>%
    select(pid) %>%
    unique() %>%
    dim()

sprintf("We next explored whether anti-M IgG was associated with protection. For each emm-typed event (n=%d), the homologous or, if unavailable, the cluster-homologous anti-M IgG Z-score (herein called cluster-related) was identified in cases before, during, and after the event, and in household contacts before and after the event. %d measurements from %d timepoints and %d individuals were included",
        n1[1],n2[1],n3[1],n4[1])

### Join conserved antigen titres to M protection dataframe

join1 <- M_protection_df %>%
    select(pid, visit_date, M = titre) %>%
    unique()

conserved <- readRDS("data/blood_IgG_titres.RDS") %>%
    spread(Antigen, titre)

joined.df <- join1 %>%
    left_join(conserved) %>%
    unique()

# Return the number of measurements in this dataframe 
joined.df %>%
    select(pid, visit_date,M) %>%
    unique() %>% dim()


#### Perform a Spearman correlation coefficient between the two dataframes 
corr.df <-
    joined.df %>%
    select(-pid, - visit_date)


corr_matrix_interaction <- corr.df %>%
    cor(use = "pairwise.complete.obs", method = "spearman") 

# Save the corrplot as an image
png("R_output/Fig6_panel_E.png", width = 1200, height = 1200, res = 300)



corrplot::corrplot.mixed(corr_matrix_interaction, upper = "color", lower = "number",
                         tl.col = "black", tl.srt = 45, lower.col = "black", number.cex = .8,
                         tl.pos = "lt") 

dev.off()

# return the number of timepoints in this dataframe
joined.df %>%
    select(pid, visit_date) %>%
    unique()%>%
    dim()

##### produce the sensitivity analysis keeping only the 7 M protins with best specificty data

M_protection_sensitivity <- readRDS("R_objects/M_protection_sensitivity_dataframe.RDS")

# Number of measurements in sensitivity dataframe 
M_protection_sensitivity %>%    
    select(pid, visit_date, Antigen, titre) %>%
    unique() %>%
    dim()

# Number of visits in sensitivity datafrmae 
M_protection_sensitivity %>%
    select(pid, visit_date) %>%
    unique() %>%
    dim()

# Number of individuals in sensitivity protection dataframe 
M_protection_sensitivity %>%
    select(pid) %>%
    unique() %>%
    dim()

# Number of events 
M_protection_sensitivity %>%
    select(event_id) %>%
    unique() %>%
    dim()

join1 <- M_protection_sensitivity %>%
    select(pid, visit_date, M = titre) %>%
    unique()

conserved <- readRDS("R_objects/blood_IgG_titres.RDS") %>%
    spread(Antigen, titre)


joined.df <- join1 %>%
    left_join(conserved) %>%
    unique()

joined.df %>%
    select(pid, visit_date,M) %>%
    unique() %>% dim()

corr.df <-
    joined.df %>%
    select(-pid, - visit_date)


corr_matrix_interaction <- corr.df %>%
    cor(use = "pairwise.complete.obs", method = "spearman") 


# Save the corrplot as an image
png("R_output/Fig6_panel_E_sensitivity.png", width = 1200, height = 1200, res = 300)


corrplot::corrplot.mixed(corr_matrix_interaction, upper = "color", lower = "number",
                         tl.col = "black", tl.srt = 45, lower.col = "black", number.cex = .8,
                         tl.pos = "lt") 

dev.off()



################################################
########## mixed effect regression #############
################################################

# This function performs a mixed effects regression to explore the relationship between
# IgG level and event within the next n days (set to 45)


# Input:
# - Antibody titre data including participant IDs, dates, and titres.
# - Event incidence data with dates of events and participant IDs.
# - Demographic data

# Description:
# - Data cleaning and preparation: Merges antibody titre data with demographic and event data, filters by antigen, and calculates variable (next at next visit within a given time period).
# - Analysis: Calculates titre thresholds, fits a generalized linear mixed-effects model to evaluate the effect of titres on event_next_n. 

# Output:
# - regression summary, visualization, AIC, and predicted probability data.


load("data/all_events_long_incidence_wgs.Rdata")

create_regression_dataframe_glmer_test <- function(path_to_titre.df,sample,class, next_event_window = 45, var_name = "titre", antigen, df = 1,slice_window = 0.1) {
    
    
    
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
                      select(pid,date = visit_date,titre, age)) %>%
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
    
    print(paste("AIC:",AIC(model)))
    
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
            x=paste0("IgG (Z score) to cluster related event"),
            y = paste0("Proportion with event \nwithin ", next_event_window, " days")) +
        guides(fill = "none", alpha = "none") +
        theme_minimal() +
        theme_universal(base_size = plot_basesize)+ 
        theme(axis.title.y = element_text(size = plot_basesize))
    
    
    
    
    plot2 <- ggplot(new_data, aes(x = titre, y = prob)) +
        geom_line(color = "blue") +
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
        labs(
            x=paste0("IgG (Z score) to cluster related event"),
            y = paste0("Probability of event \nwithin ", next_event_window, " days")) +
        guides(fill = "none", alpha = "none") +
        theme_minimal() +
        annotate("text", x = max_threshold + 1, y = 0.35, 
                 label = paste0("OR: ", round(ORs[2], 2), 
                                "\nCI:", round(CIs[1, 1], 2), ", ", round(CIs[1, 2], 2), 
                                "\n", p_label), 
                 size = 5,
                 hjust = 0) +
        theme_universal(base_size = plot_basesize) + 
        theme(axis.title.y = element_text(size = plot_basesize))
    
    plot <- cowplot::plot_grid(plot1,plot2,ncol = 1)
    
    
    # Instead of printing, add the table and plot to a list
    results_list <- list("table" = tb1, "plot" = plot, "AIC_glmer" = model_aic)
    
    
    results_list$probabilities2 <- new_data
    
    
    
    
    # Return the list
    
    
    return(results_list)
    
    
}

AIC_glmer <- data.frame(Antigen = character(), AIC_glmer = numeric(), stringsAsFactors = FALSE)

result <- create_regression_dataframe_glmer_test(path_to_titre.df =  M_protection_df,
                                                 sample = "Blood",
                                                 class = "IgG",
                                                 var_name = "titre",
                                                 antigen = "to cluster homologous M peptide by Z score",
                                                 slice_window = 0.1)



plot_09_fig06_panelE <- result$plot
plot_09_fig06_panelE 


AIC_glmer <- rbind(AIC_glmer, data.frame(Antigen = "M", AIC_glmer = result$AIC_glmer))


##### repeat this performing a sensitivity analysis with only top 7 performing anti-M IgG measurements 

result3 <- create_regression_dataframe_glmer_test(path_to_titre.df =  M_protection_sensitivity,
                                                  sample = "Blood",
                                                  class = "IgG",
                                                  var_name = "titre",
                                                  antigen = "to cluster homologous M peptide by Z score",
                                                  slice_window = 0.1)



plot_09_fig06_panelE_sensitivity_v2 <- result3$plot
plot_09_fig06_panelE_sensitivity_v2


###########################################
###########################################


# Repeat this analysis adjusting for age group, household size and sex

create_regression_dataframe_glmer_adjusted <- function(path_to_titre.df,sample,class, next_event_window = 45, var_name = "titre", antigen, df = 1,slice_window = 0.1) {
    
    
    
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
                      select(pid,date = visit_date,titre, age_grp,sex,hhsize)) %>%
        mutate(titre_bucket = floor(titre / bucket_size) * bucket_size)  %>%
        group_by(pid) %>%
        #   fill(titre) %>%
        ungroup() %>%
        mutate(hid = substring(pid,1,3)) %>%
        unique()
    
    
    stopifnot("pid" %in% names(final_df_3))
    
    final_df_3$age_grp <- relevel(final_df_3$age_grp, ref = "Over 18 years")
    
    model <- lme4::glmer(event_next_n ~ titre + age_grp +sex + hhsize + (1 | pid) + (1 | hid), data=final_df_3, family=binomial)
    
    tb1 <- model %>%
        tbl_regression(exponentiate = TRUE) %>%
        #  add_nevent() %>%
        bold_p(t = 0.05) %>%
        modify_header(list(label ~ paste(sample, class, "titre"))) # Modify header here 
    
    print(paste("AIC:",AIC(model)))
    
    # Instead of printing, add the table and plot to a list
    results_list <- list("table" = tb1)
    
    
    
    
    # Return the list
    # return(results_list)
    
    return(results_list )
    
    
}

output <- create_regression_dataframe_glmer_adjusted(path_to_titre.df =  M_protection_df,
                                                     sample = "Blood",
                                                     class = "IgG",
                                                     var_name = "titre",
                                                     antigen = "to cluster homologous M peptide by Z score",
                                                     slice_window = 0.1)


output$table

shelf(officer,flextable)

combined_table <- output$table

### Draw tables for  manscript ### 

# Create a new Word document
doc <- read_docx()
# Add regression table to word document
doc <- add_table_to_doc(doc, combined_table, "Blood IgG M titres and probability of event in next 45 days")


# Save the Word document in landscape orientation
print(doc, target ="R_output/Supplementary_M_adjusted_GLMER.docx")

rm(doc)

#### sensitivity 2 : only keep the best performing M types in the dataframes. 

output3 <- create_regression_dataframe_glmer_adjusted(path_to_titre.df = M_protection_sensitivity,
                                                      sample = "Blood",
                                                      class = "IgG",
                                                      var_name = "titre",
                                                      antigen = "to cluster homologous M peptide by Z score",
                                                      slice_window = 0.1)


output3$table

###########################################################################
####### combining GLMER with conserved antigen above threshold ############
##########################################################################


create_regression_glmer_conserved_M <- function(M_titres, conserved_titres,sample,class, next_event_window = 45, var_name = "titre", n_tile = 4, antigen) {
    
    
    titre_breakpoint_df <- tibble(Antigen = c("SLO","SpyAD", "SpyCEP"), titre_breakpoint = c(4.3,4.1,4.3))
    
    titre_breakpoint <- titre_breakpoint_df %>%
        filter(Antigen == antigen) %>%
        pull(titre_breakpoint)
    
    # Filter for specific antigen
    antigen_df <- conserved_titres %>% 
        left_join(age) %>%
        filter(Antigen == antigen & !is.na(pid)) %>%
        mutate(titre_below_threshold = ifelse(titre <= titre_breakpoint, titre, titre_breakpoint),
               titre_above_threshold = ifelse(titre > titre_breakpoint, titre - titre_breakpoint, 0)) %>%
        select(pid, visit_date,conserved_titre = titre, titre_above_threshold)
    
    # Read titre data and merge with start dates and age, adjust visit dates, and categorize titres
    fun_titres <- M_titres %>%
        left_join(antigen_df) %>%
        left_join(start_dates) %>%
        left_join(age) %>%
        mutate(visit_date = as.numeric(ymd(as.character(visit_date))) - entry_1)
    
    
    
    
    
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
    
    
    final_df_3  <-fun_df %>%
        left_join(antigen_df %>% 
                      select(pid,date = visit_date,titre,conserved_titre,titre_above_threshold, age_grp, hhsize, sex)) %>%
        group_by(pid) %>%
        #   fill(titre) %>%
        ungroup() %>%
        mutate(hid = substring(pid,1,3)) %>%
        unique()
    
    # then do logistic regression: 
    
    
    stopifnot("pid" %in% names(final_df_3))
    
    
    
    model2 <- lme4::glmer(event_next_n ~ titre + titre_above_threshold + age_grp + sex + hhsize+  (1 | pid) + (1 | hid), data=final_df_3, family=binomial)
    
    print(paste0("model2 (M + ", antigen, " above transitionpoint -fully adjusted):", AIC(model2)))
    
    tb2 <- model2 %>%
        tbl_regression(exponentiate = TRUE) %>%
        add_nevent() %>%
        bold_p(t = 0.05) 
    
    
    
    # Use broom::tidy to extract model2 results
    model2_tidy <- broom::tidy(model2, conf.int = TRUE, exponentiate = TRUE) %>%
        rename(estimate = estimate, conf.low = conf.low, conf.high = conf.high, p.value = p.value) %>%
        mutate(antigen = antigen)  # Add Antigen column for identification
    
    
    # Create a range of values for M titre and the continuous antigen variable
    M_range <- seq(min(final_df_3$titre, na.rm = TRUE), max(final_df_3$titre, na.rm = TRUE), length.out = 100)
    cont_var_range <- seq(min(final_df_3$titre_above_threshold, na.rm = TRUE), max(final_df_3$titre_above_threshold, na.rm = TRUE), length.out = 100)
    
    # Create a grid of data with all combinations of M and the continuous antigen variable
    grid_data <- expand.grid(titre = M_range,titre_above_threshold = cont_var_range)
    
    # Predict probabilities using the model for the grid data
    grid_data$prob <- predict(model, newdata = grid_data, type = "response", re.form = NA)
    
    # Append the current grid_data to the combined_df using bind_rows
    grid_data$antigen <- antigen  # Add antigen name
    
    
    
    plot <- ggplot(grid_data, aes(x = titre_above_threshold + titre_breakpoint, y = titre, fill = prob)) +
        geom_tile()  + # Creates the heatmap tiles based on probability
        scale_fill_gradient(name = "Probability", low = "white", high = "navy") +
        labs(
            title = paste(antigen),  # Add antigen to title
            x = paste0("IgG above ", titre_breakpoint, " (log10 RLU/mL)"),
            y = "Anti-M IgG (Z score)"
        ) +
        theme_minimal()+
        theme_universal(base_size = plot_basesize)
    
    
    results_list <- list("tb1" = tb1, "tb2" = tb2,"tb3" = tb3, "plot" = plot, "model2_results" = model2_tidy )
    
    # Return the list
    
    
    return(results_list)
    
    
}

################################################################
###### Deploy this function on the protection dataframe ########
################################################################

tb2_list <- list()
plot_list <- list()
plot_row_list <- list()

results_df <- data.frame(antigen = character(), term = character(), estimate = numeric(),
                         conf.low = numeric(), conf.high = numeric(), p.value = numeric(),
                         stringsAsFactors = FALSE)


Antigen_list = c("SLO","SpyAD","SpyCEP")
for (ag in Antigen_list) {
    results <- create_regression_glmer_conserved_M(
        M_titres = M_protection_df,
        conserved_titres = readRDS("R_objects/blood_IgG_titres.RDS"),
        sample = "Blood",
        class = "IgG",
        var_name = "titre",
        antigen = ag
    )
    
    tb2_list[[ag]] <- results$tb2
    plot_list[[ag]] <- results$plot
    
    # Combine three plots into one row (shared fill scale is automatic if using ggplot)
    plot_row_list[[ag]] <- cowplot::plot_grid(
        results$plot,
        align = "hv", nrow = 1  # Arrange plots in one row
    )
    
    # Append model2_results to the master dataframe
    results_df <- rbind(results_df, results$model2_results)
}

# Combine all antigen rows into a single plot
combined_plot <- cowplot::plot_grid(plotlist = plot_row_list, ncol = 3, align = "v")

plot_09_fig06_panelF <- combined_plot
plot_09_fig06_panelF


tbl_merge( tb2_list,
           tab_spanner = Antigen_list)

# Display all plots
print(combined_plot)

# Example: Display merged tables for all antigens

# Display all plots
print(combined_plot)
as.vector(levels(as.factor(results_df$term)))


final_fp_df <-results_df %>%
    mutate(
        # Add descriptive labels for the variables
        variable_label = case_when(
            grepl("titre_above_threshold", term) ~ "IgG level",
            grepl("titre", term) ~ "IgG level",
            grepl("age_grp", term) ~ "Age group",
            grepl("hhsize", term) ~ "Household size",
            grepl("sexFemale", term) ~ "Sex",
            TRUE ~ NA_character_  # Default to NA if no match
        ),
        # Simplify term names for clarity
        term = case_when(
            grepl("titre_above_threshold", term) ~ "Conserved IgG",
            grepl("titre", term) ~ "M IgG",
            grepl("age_grp< 2 years", term) ~ "< 2 years",
            grepl("age_grp12-18 years", term) ~ "12-18 years",
            grepl("age_grp2-4 years", term) ~ "2-4 years",
            grepl("age_grp5-11 years", term) ~ "5-11 years",
            grepl("sexFemale", term) ~ "Female",
            grepl("hhsize", term) ~ "Household size",
            TRUE ~ term  # Retain term if it doesn't match any pattern
        )) %>%
    mutate(
        term = factor(term, levels = c(
            "Conserved IgG",
            "M IgG",
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



# Ensure 'term' and 'variable_label' are ordered as intended
final_fp_df <- final_fp_df %>%
    mutate(
        term = factor(term, levels = c(
            "Conserved IgG",
            "M IgG",
            "Over 18 years (ref)",
            "< 2 years",
            "2-4 years", 
            "5-11 years", 
            "12-18 years", 
            "Male (ref)",
            "Female",
            "Household size" 
        )),
        variable_label = forcats::fct_relevel(variable_label, 
                                              "IgG level", "Age group", "Sex", "Household size")
    )

# Plot

new_rows <- unique(results_df$antigen) %>%
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


M_forrest <- final_fp_df %>%
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
        aes(x = 1, y = term), color = "blue", size = 3
    ) +
    facet_grid(
        cols = vars(antigen),
        rows = vars(forcats::fct_relevel(variable_label, 
                                         "IgG level", "Age group", "Sex", "Household size")), 
        scales = "free", 
        space = "free_y"
    ) +
    scale_x_log10() +  # Log scale for odds ratios
    theme_universal(base_size = plot_basesize) +
    theme(strip.text.y = element_text(angle = 0))

plot_09_fig06_panelG <- M_forrest
plot_09_fig06_panelG




#  only using "M1", "M4", "M44", "M75", "M89", "M25", "M87" titres 

tb1_list <- list()
tb2_list <- list()
tb3_list <- list()
plot_list <- list()
plot_row_list <- list()

results_df <- data.frame(antigen = character(), term = character(), estimate = numeric(),
                         conf.low = numeric(), conf.high = numeric(), p.value = numeric(),
                         stringsAsFactors = FALSE)



for (ag in Antigen_list) {
    results <- create_regression_glmer_conserved_M(
        M_titres = M_protection_sensitivity,
        conserved_titres = readRDS("R_objects/blood_IgG_titres.RDS"),
        sample = "Blood",
        class = "IgG",
        var_name = "titre",
        antigen = ag
    )
    
    tb1_list[[ag]] <- results$tb1
    tb2_list[[ag]] <- results$tb2
    tb3_list[[ag]] <- results$tb3
    plot_list[[ag]] <- results$plot
    
    # Combine three plots into one row (shared fill scale is automatic if using ggplot)
    plot_row_list[[ag]] <- cowplot::plot_grid(
        results$plot,
        align = "hv", nrow = 1  # Arrange plots in one row
    )
    
    # Append model2_results to the master dataframe
    results_df <- rbind(results_df, results$model2_results)
}

# Combine all antigen rows into a single plot
combined_plot <- cowplot::plot_grid(plotlist = plot_row_list, ncol = 3, align = "v")

plot_09_fig06_panelF <- combined_plot
plot_09_fig06_panelF
tbl_merge( tb1_list,
           tab_spanner = Antigen_list)

tbl_merge( tb2_list,
           tab_spanner = Antigen_list)

tbl_merge( tb3_list,
           tab_spanner = Antigen_list)

# Display all plots
print(combined_plot)

# Example: Display merged tables for all antigens

# Display all plots
print(combined_plot)
as.vector(levels(as.factor(results_df$term)))


final_fp_df <-results_df %>%
    mutate(
        # Add descriptive labels for the variables
        variable_label = case_when(
            grepl("titre_above_threshold", term) ~ "IgG level",
            grepl("titre", term) ~ "IgG level",
            grepl("age_grp", term) ~ "Age group",
            grepl("hhsize", term) ~ "Household size",
            grepl("sexFemale", term) ~ "Sex",
            TRUE ~ NA_character_  # Default to NA if no match
        ),
        # Simplify term names for clarity
        term = case_when(
            grepl("titre_above_threshold", term) ~ "Conserved IgG",
            grepl("titre", term) ~ "M IgG",
            grepl("age_grp< 2 years", term) ~ "< 2 years",
            grepl("age_grp12-18 years", term) ~ "12-18 years",
            grepl("age_grp2-4 years", term) ~ "2-4 years",
            grepl("age_grp5-11 years", term) ~ "5-11 years",
            grepl("sexFemale", term) ~ "Female",
            grepl("hhsize", term) ~ "Household size",
            TRUE ~ term  # Retain term if it doesn't match any pattern
        )) %>%
    mutate(
        term = factor(term, levels = c(
            "Conserved IgG",
            "M IgG",
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



# Ensure 'term' and 'variable_label' are ordered as intended
final_fp_df <- final_fp_df %>%
    mutate(
        term = factor(term, levels = c(
            "Conserved IgG",
            "M IgG",
            "Over 18 years (ref)",
            "< 2 years",
            "2-4 years", 
            "5-11 years", 
            "12-18 years", 
            "Male (ref)",
            "Female",
            "Household size" 
        )),
        variable_label = forcats::fct_relevel(variable_label, 
                                              "IgG level", "Age group", "Sex", "Household size")
    )

# Plot

new_rows <- unique(results_df$antigen) %>%
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


M_forrest <- final_fp_df %>%
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
        aes(x = 1, y = term), color = "blue", size = 3
    ) +
    facet_grid(
        cols = vars(antigen),
        rows = vars(forcats::fct_relevel(variable_label, 
                                         "IgG level", "Age group", "Sex", "Household size")), 
        scales = "free", 
        space = "free_y"
    ) +
    scale_x_log10() +  # Log scale for odds ratios
    theme_universal(base_size = plot_basesize) +
    theme(strip.text.y = element_text(angle = 0))

plot_09_fig06_panelG_sensitivity_v2 <- M_forrest
plot_09_fig06_panelG_sensitivity_v2



