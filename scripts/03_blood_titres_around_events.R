# Title: Analysis of IgG Responses in Blood Samples around Strep A events 
# Version: 1.0
# Date: 2024-07-31
# Author : Dr. Alexander J. Keeley
# Inputs: data/blood_IgG_titres.RDS, data/all_events_long_immunology.Rdata, data/all_events_long_incidence_wgs.RData
# Outputs: Summaries of events numbers and timing to events. Plot files. 


# Description:
# This script analyzes the IgG  in blood aamples from a study of humoral immunity to Streptococcus pyogenes (Strep A). The analysis includes:
# 1. Calculating the total number of microbiologically confirmed events and the availability of paired blood IgG
# 2. Summarizing the timing of pre- and post-event antibody titre measurements.
# 3. Generating descriptive statistics and plots for titres relative to microbiologically confirmed events.

# Requirements:

# Package Versions
# librarian     librarian   1.8.1
# ggdist           ggdist   3.3.2
# ggpubr           ggpubr   0.6.0
# tidyverse     tidyverse   2.0.0
# wesanderson wesanderson   0.3.7


##################
# Setup environment
source("scripts/setup_environment.R")

# load functions for datacleaning and plotting:
source("scripts/load_functions.R")

# load packages using librarian 
shelf(tidyverse,ggdist,ggpubr,wesanderson, gtsummary,lme4, flextable,officer, tidyr, broom)


##### create "per event dataframes" 

IgG_blood<- create_events_dataframe("data/blood_IgG_titres.RDS", var = "titre", class = "IgG", sample = "Blood")

#######  Describe the numbers of events used in analyses: 

#### Number of events: 

IgG_blood$event_titres_df_IgG_Blood %>%
    select(event_id,event_type) %>%
    unique() %>%
    ungroup() %>%
    select(event_type) %>%
    gtsummary::tbl_summary()


IgG_blood$fold_changes_IgG_Blood %>% 
    filter(!is.na(pre_post)) %>%
    select(event_id, event_type) %>%
    unique() %>%
    ungroup() %>%
    select(event_type) %>%
    gtsummary::tbl_summary()


# Number of people with 
IgG_blood$event_titres_df_IgG_Blood %>%
    select(event_id,pid) %>%
    unique() %>%
    ungroup() %>%
    pull(pid) %>%
    unique() %>%
    length()

# total number of microbiological confirmed events in analysis 
n1 <- IgG_blood$event_titres_df_IgG_Blood %>%
    pull(event_id) %>%
    unique() %>% length

#### Number of events with paired DBS titres: 

# total number of microbiological confirmed events with blood IgG measure before and after 
n2 <- IgG_blood$relative_changes_IgG_Blood %>%
    filter(!is.na(pre_post)) %>%
    pull(event_id) %>%
    unique() %>% length ()


IgG_event_sentence <- 
    sprintf(
    "From %d microbiologically confirmed events, paired pre- and post-event blood IgG titres were available from %d for analysis.",
    n1, n2)

print(IgG_event_sentence)
rm(n2, IgG_event_sentence)


#######  Describe the timing around events of titre measurements: 

IgG_blood$fold_changes_IgG_Blood %>%
    as.data.frame() %>%
    mutate(time_between_measurements = 
               abs(diff_Post) + abs(diff_Pre)) %>%
    select(event_id, time_between_measurements) %>%
    unique() %>%
    summarise(n = n(),
              median = median(abs(time_between_measurements), na.rm=T),
              Q1 = quantile(abs(time_between_measurements), probs = 0.25, na.rm = T),
              Q3 = quantile(abs(time_between_measurements), probs = 0.75, na.rm = T),
              IQR = IQR(abs(time_between_measurements), na.rm=T),
              mean =  mean(abs(time_between_measurements), na.rm=T),
              min =  min(abs(time_between_measurements), na.rm=T),
              max =  max(abs(time_between_measurements), na.rm=T))






#####################################################################
####### paired antibody titres:pre -post events #####################
#####################################################################

# IgG blood 

output <- generate_paired_analysis(IgG_blood$event_titres_df_IgG_Blood,
                                   include_event_as_paired = T, 
                                   sample = "Blood",
                                   class = "IgG"
)
print(output$results)

blood_IgG_for_results <- output$results %>%
    ungroup() %>%
    summarise(min = min(mean_difference),
              max = max(mean_difference))

sprintf(" In paired DBS samples there was a significant rise in geometric mean IgG titre between the pre- and post-event sample to all antigens, ranging from %.2f to %.2f log10 RLU/mL (p<0.0001 for all comparisons, supplementary figure xx).", 
        blood_IgG_for_results$min,
        blood_IgG_for_results$max)


plot_05_supp01_v1.0 <- output$plot + labs(y = "IgG level (RLU/mL)")
plot_05_supp01_v1.0


rm(output, blood_IgG_for_results)


#####################################################################
####### produce plots: absolute titres around events ################
#####################################################################


plot_05_main01_v1.0 <- Longitudinal_event_plot_by_age(IgG_blood$titres_relative_to_baseline, sample = "blood", class = "IgG") + labs(y = "Absolute change from pre event IgG level (log10 RLU/mL)")
plot_05_main01_v1.0 


plot_05_main02_v1.0 <- longitudinal_event_site_plot(IgG_blood$titres_relative_to_baseline_IgG_Blood, sample = "blood", class = "IgG",age_cut = 100) + labs(y = "Absolute change from pre event IgG level (log10 RLU/mL)")
plot_05_main02_v1.0


#####################################################################
####### produce plots:        absolute changes by type  ###########
#####################################################################


df_to_adjust_event_type <- 
    IgG_blood$relative_changes_IgG_Blood %>%
    as.data.frame() %>%
    mutate(Antigen = factor(Antigen, levels = c("GAC", "SLO", "SpyAD", "SpyCEP", "DNAseB"))) %>%
    left_join(age) %>%
    group_by(Antigen) %>%
    ##
    # Use this chuck of code to include the pre to event changes in dataframe if no pre_post titre is available. 
    
    mutate(pre_post = 
               case_when(
                   is.na(pre_post) & ! is.na(pre_event) ~ pre_event,
                   T ~ pre_post
               )) %>%
    filter(event_type != "other",
        #    age < 18,
           !is.na(pre_post)) %>%
    select(event_type, Antigen, pre_post)



# Define the antigens you're interested in

events <- unique(df_to_adjust_event_type$event_type)
antigens <- c("GAC", "SLO", "SpyAD", "SpyCEP","DNAseB")

# Initialize a list to store Dunn test results

dunn_results <- list()

# Loop over Antigen 
for (antigen in antigens) {
    # Filter the data for the current Antigen type
    data_subset <- df_to_adjust_event_type %>%
        filter(Antigen == antigen)
    
    # Perform Dunn's test
    test_result <- dunn.test::dunn.test(x = data_subset$pre_post,
                                        g = data_subset$event_type,
                                        method = "bonferroni")
    
    # Store the results with a unique name for each Antigen type
    dunn_results[[antigen]] <- test_result
}

# Create an empty summary dataframe
summary_df <- data.frame(Antigen = character(),
                         adjusted_p_value = numeric(),
                         group1 = character(),
                         group2 = character(),
                         label = character(),
                         y.position = numeric(),
                         stringsAsFactors = FALSE)

# Extract data for each result stored in the list and create the summary dataframe
for (name in names(dunn_results)) {
    # Extract the adjusted p-values and comparisons
    adjusted_p_values <- dunn_results[[name]]$P.adjusted
    comparisons <- dunn_results[[name]]$comparisons
    
    # Split the comparisons to extract the groups
    split_comparisons <- strsplit(comparisons, " - ")
    group1 <- sapply(split_comparisons, `[`, 1)
    group2 <- sapply(split_comparisons, `[`, 2)
    
    # Create labels based on p-value significance
    labels <- sapply(adjusted_p_values, function(p) {
        if (p < 0.0001) {
            "****"
        } else if (p < 0.001) {
            "***"
        } else if (p < 0.01) {
            "**"
        } else if (p < 0.05) {
            "*"
        } else {
            "NA"  # Indicate non-significant p-values
        }
    })
    
    # Calculate the maximum log10 transformed titre cahnge for plotting
    y_position <- max(df_to_adjust_event_type$pre_post, na.rm = TRUE)
    
    # Combine the individual results into the summary dataframe
    temp_df <- data.frame(Antigen = rep(name, length(adjusted_p_values)),
                          adjusted_p_value = adjusted_p_values,
                          group1 = group1,
                          group2 = group2,
                          label = labels,
                          y.position = y_position,
                          stringsAsFactors = FALSE)
    
    # Bind the temp dataframe to the summary dataframe
    summary_df <- rbind(summary_df, temp_df)
}

# Print the summary dataframe
print(summary_df)


# Print the plot

plot_05_supp02_v1.0 <- df_to_adjust_event_type %>%
    ggplot(
        aes(x = event_type, y = pre_post)) +
    
    geom_boxplot(width = 0.05, outlier.size = 0.3, fill = NA, aes(fill = event_type , col = event_type)) +
    
    ggdist::geom_dots(
        aes(fill = event_type , col = event_type),
        width = 0.3,
        position = position_nudge(x = 0.05)
    ) +
    
    ggdist::stat_slab(
        aes(fill = event_type , col = event_type),
        side = 'bottom', 
        position = position_nudge(x = -0.05), 
        width = 0.3
    ) +
    
    facet_wrap(~Antigen, ncol = 5) +
    theme_minimal() +
    scale_fill_manual(name = "Event Type", values = wesanderson::wes_palette("Zissou1")) + 
    scale_color_manual(name = "Event Type", values = wesanderson::wes_palette("Zissou1")) + 
    labs(y = "Absolute change from pre event IgG level (log10 RLU/mL)",
         x = "Event type") +
  #  ggtitle("Comparison of relative antibody changes between microbiologically confirmed event types in participants under 12") +
    # ggpubr::stat_compare_means()  +
    theme_universal(base_size = plot_basesize) +
   # scale_x_discrete() + 
    summary_df %>%
    filter(!is.na(label) & label != "NA") %>%
    ggpubr::stat_pvalue_manual(
        label = "label",
        tip.length = 0.01,
        size = 12 / (14/5),
        step.increase = 0.035,
        step.group.by = "Antigen",
        position = position_dodge(width = 0.01)
    ) 

plot_05_supp02_v1.0 <- plot_05_supp02_v1.0 + theme(axis.text.x = element_blank())
plot_05_supp02_v1.0

df_to_adjust_event_type <- 
    IgG_blood$relative_changes_IgG_Blood %>%
    as.data.frame() %>%
    mutate(Antigen = factor(Antigen, levels = c("GAC", "SLO", "SpyAD", "SpyCEP", "DNAseB"))) %>%
    left_join(age) %>%
    group_by(Antigen) %>%
    ##
    # Use this chuck of code to include the pre to event changes in dataframe if no pre_post titre is available. 
    
    mutate(pre_post = 
               case_when(
                   is.na(pre_post) & ! is.na(pre_event) ~ pre_event,
                   T ~ pre_post
               )) %>%
    filter(event_type != "other",
           #    age < 18,
           !is.na(pre_post)) %>%
    select(pid, age_grp, event_type, Antigen, pre_post)


supplementary_response_changes_age_facet <- df_to_adjust_event_type %>%
    
    ggplot(
        aes(x = event_type, y = pre_post)) +
    geom_boxplot(aes(fill = event_type , col = event_type), alpha = 0.8) +
    geom_point(aes(fill = event_type), shape = 21, color = "black", stroke = 0.2) +
    facet_grid(Antigen ~ age_grp) +
    theme_minimal() +
    scale_fill_manual(name = "Event Type", values = wesanderson::wes_palette("Zissou1")) + 
    scale_color_manual(name = "Event Type", values = wesanderson::wes_palette("Zissou1")) + 
    labs(y = "Absolute change from pre event IgG level (log10 RLU/mL)",
         x = "Event type") +
    #  ggtitle("Comparison of relative antibody changes between microbiologically confirmed event types in participants under 12") +
    # ggpubr::stat_compare_means()  +
    theme_universal(base_size = plot_basesize) +
    # scale_x_discrete() + 
    summary_df %>%
    filter(!is.na(label) & label != "NA") %>%
    ggpubr::stat_pvalue_manual(
        label = "label",
        tip.length = 0.01,
        size = 12 / (14/5),
        step.increase = 0.035,
        step.group.by = "Antigen",
        position = position_dodge(width = 0.01)
    ) + 
    theme(axis.text.x = element_blank())


plot_05_supp02_v1.0 <- supplementary_response_changes_age_facet

rm(dunn_results,summary_df, df_to_adjust_event_type)

#####################################################################
####### Generalised linear mixed effects regression  ################
#####################################################################


# Assosiation between age and event type on serological responses to Strep A events 


## Repeat the analysis removing interaction term given no significnat interactions 
#### This code returns the AIC of the models #####


# List to store regression tables
tables_list <- list()

# Initialize an empty data frame to store the results
results_df <- data.frame(
    Antigen = character(),
    Variables = character(),
    AIC = numeric(),
    stringsAsFactors = FALSE
)


for (Ag in unique(IgG_blood$fold_changes_IgG_Blood$Antigen)) {
    
    model.df <- IgG_blood$fold_changes_IgG_Blood %>%
        as.data.frame() %>%
        left_join(age) %>%
        filter(Antigen == Ag) %>%
        filter(!is.na(pre_post),
               event_type != "other") %>%
        select(
            pid, Antigen, pre_post, age_grp, event_type
        )   %>%
        mutate(
            # Set reference levels for factors
            age_grp = relevel(as.factor(age_grp), ref = "Over 18 years"),
            event_type = relevel(as.factor(event_type), ref = "throat disease")
        )
    
    
    #model <- lme4::glmer(pre_post ~ age_grp + event_type + (1|pid), data = model.df, family = gaussian)
    
    model <- lmerTest::lmer(pre_post ~ age_grp + event_type + (1|pid), data = model.df)
    
    
    tb1 <- model %>%
        tbl_regression(exponentiate = F,
                       pvalue_fun = format_p_value) %>%
        modify_header(estimate = "**Log10 Antibody Titer Change**")
    
    
    # Tidy the model
    model_summary <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE)
    
    model_summary <- model_summary %>%
        mutate(
            term = recode(term, 
                          "age_grpOver 18 years" = "Age Group: Over 18",
                          "event_typeSkin disease" = "Event Type: Skin Disease"),
            significant = ifelse(p.value < 0.05, "Yes", "No"),
           p.value = format_p_value(p.value) # Apply the custom formatting function
        )
    
    
    
    # Enhanced forest plot
   plot <-  ggplot(model_summary, aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high, color = significant)) +
        geom_pointrange() +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
        coord_flip() +
        scale_color_manual(values = c("Yes" = "blue", "No" = "black")) +
        labs(
            title = paste(Ag),
            x = "Predictor",
            y = "Log10 Titer Change (Estimate)",
            color = "Significant"
        ) +
        theme_minimal()
    
    
    
    # Create a regression table and add it to the list
    
    print(tb1)
    print(plot)
    print(model_summary)
    
    tables_list[[Ag]] <- tb1
    
    
    
    
}

tbl_merge(tables_list)

# Merge tables side by side with custom headers for each antigen
combined_table <- tbl_merge(
    tables_list,
    tab_spanner = names(tables_list)
)


# Print the combined table
print(combined_table)


# Convert the combined gtsummary table to a flextable object for Word export
ft <- as_flex_table(combined_table)


# Add title and legend 
ft <- ft %>%
    add_header_lines(values = "Assosiation between age and event type on serological responses to Strep A events ") %>%
    add_footer_lines(values = "The absolute change in log10 blood IgG titres (RLU/mL) around microbiologically confirmed Strep A events.A generalized linear mixed model (GLMM) was used to assess the impact of age group and event type on absolute titre changes (log10 RLU/mL) around events. The model included the predictors age group and event type, as well as their interaction, with participant identifier included as a random effect to account for repeated measure from individuals.")


# Create a Word document and add the table
doc <- read_docx() %>%
    body_add_flextable(ft) %>%
    body_add_par(" ")

# Save the table as a Word document
print(doc, target = "R_output/table_supp01.docx")


#Plot forests plots for this data 

# Initialize an empty dataframe to store the results for all antigens
combined_results_df <- data.frame()

for (Ag in unique(IgG_blood$fold_changes_IgG_Blood$Antigen)) {
    
    model.df <- IgG_blood$fold_changes_IgG_Blood %>%
        as.data.frame() %>%
        left_join(age) %>%
        filter(Antigen == Ag) %>%
        filter(!is.na(pre_post),
               event_type != "other") %>%
        select(
            pid, Antigen, pre_post, age_grp, event_type
        )   %>%
        mutate(
            # Set reference levels for factors
            age_grp = relevel(as.factor(age_grp), ref = "Over 18 years"),
            event_type = relevel(as.factor(event_type), ref = "throat disease")
        )
    
    # Fit the generalized linear mixed model
    model <- lme4::glmer(pre_post ~ age_grp + event_type + (1|pid), data = model.df, family = gaussian)
   # model <- lmerTest::lmer(pre_post ~ age_grp + event_type + (1|pid), data = model.df)
    
    
    # Extract the model summary, including coefficients and confidence intervals
    model_summary <- broom::tidy(model, conf.int = TRUE, exponentiate = F)
    
    # Exclude sd_Observation and sd_Intercept from the results
    model_summary <- model_summary %>%
        filter(!grepl("^sd_", term)) %>%
        mutate(
            # Label the terms with their reference level
            variable_label = case_when(
                grepl("^age_grp", term) ~ paste0("Age Group"),
                grepl("^event_type", term) ~ paste0("Event Type")
            ),
            antigen = Ag # Add the antigen as a column
        )
    
    # Append the model summary to the combined dataframe
    combined_results_df <- rbind(combined_results_df, model_summary)
}

combined_results_df %>%
    mutate(term = gsub("^age_grp|^event_type", "", term))


# Now create a forest plot with facets by antigen and display reference levels on the x-axis
test <- combined_results_df %>%
    mutate(term = gsub("^age_grp|^event_type", "", term))
   
levels(as.factor(test$term))

new_rows <- unique(test$antigen) %>%
    expand.grid(antigen = ., term = c("Over 18 years (ref)", "throat disease (ref)")) %>%
    mutate(
        effect = NA, 
        group = NA, 
        estimate = NA, 
        std.error = NA, 
        statistic = NA, 
        conf.low = NA, 
        conf.high = NA,
          variable_label = ifelse(term == "Over 18 years (ref)","Age Group", "Event Type")
    )

# Combine the new rows with the original dataframe
combined_results_df <- bind_rows(test, new_rows)
 
as.factor(combined_results_df$term) %>% levels()

combined_results_df <- 
    combined_results_df %>%  mutate(
        term = factor(term, levels = c(
            "(Intercept)", 
            "< 2 years", 
            "2-4 years", 
            "5-11 years", 
            "12-18 years", 
            "Over 18 years (ref)", 
            "throat carriage",
            "skin carriage", 
            "skin disease", 
            "throat disease (ref)"
        ))
    )

levels(as.factor(combined_results_df$term))

levels(combined_results_df$term)

forest_plot <- combined_results_df %>%
    mutate(antigen = factor(antigen, levels = c("GAC", "SLO", "SpyAD", "SpyCEP", "DNAseB"))) %>%
    filter(term != "(Intercept)") %>% # Exclude the intercept
    ggplot(aes(x = estimate, y = term)) +
    geom_point() +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  #  facet_wrap(~ antigen, scales = "free") + # Facet by antigen
    labs(
     #   title = "Association between age and event type and absolute log10 RLU/mL titre changes around event",
        x = "IgG level change (log10 RLU/mL) + 95% CIs",
        y = ""
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    # Add a point for the reference level at x = 1, y = position of reference variable

    theme_minimal() +
    geom_point(data = new_rows ,
               aes(x = 0, y = term), color = "blue", size = 3) +
    # Separate age_grp and event_type under subheadings on y-axis, including reference level
    facet_grid(cols = vars(antigen), rows = vars(variable_label), scales = "free", space = "free_y") +
    theme_universal(base_size = plot_basesize)

# Print the forest plot
print(forest_plot)


plot_05_main02.2 <- forest_plot




#####################################################################
####### Assosiation between baseline titre and pre-post response  ###
#####################################################################


temp1 <- IgG_blood$relative_changes_IgG_Blood %>%
    as.data.frame() %>%
    left_join(age) %>% 
    group_by(Antigen) %>%
    ##
    # Use this chuck of code to include the pre to event changes in dataframe if no pre_post titre is available. 
    
    mutate(pre_post = 
               case_when(
                   is.na(pre_post) & ! is.na(pre_event) ~ pre_event,
                   T ~ pre_post
               )) %>%
    ##
    
    select(event_id,Antigen, pre_post)

temp2 <- IgG_blood$event_titres_df_IgG_Blood %>%
    as.data.frame() %>%
    filter(status == "Pre") %>%
    left_join(age) %>% 
    select(event_id, age, Antigen, titre)




plot_05_supp03_v1.0 <-
    temp1 %>%
    left_join(temp2) %>%
    ggplot(
        aes(x = titre, 
            y = pre_post)
    ) +
    
    geom_smooth(method = "lm", se = FALSE, col = "black") +  # Linear regression line
    geom_point(alpha = 0.8 ,aes(
        col = Antigen)) + 
    facet_wrap(~Antigen, ncol = 5) +
    scale_colour_manual(values = c("SpyCEP" = "#FDC086", "SpyAD" = "#d19c2f", "SLO" = "#386CB0", "GAC" = "#7FC97F", "DNAseB" = "#BEAED4")) +
    theme_minimal() +
    theme(
        text = element_text(color = "black"))+
    guides(col = "none") +
    labs(
        x = "Baseline IgG level (log10 RLU/mL)",
        y = "Absolute change from pre event IgG level (log10 RLU/mL)",
    #    title = "Change in blood IgG titre from baseline around events relative to baseline titres",
     #   caption = "Coefficient and p-value determined using Pearson’s method with linear regression line plotted"
    ) +
    ggpubr::stat_cor(
        aes(label = paste0(
            "italic(r) == ", ..r..,  # Correlation coefficient
            ifelse(..p.. < 0.0001, 
                   "~~italic(p) < 0.0001",  # p < 0.0001 case
                   paste0("~~italic(p) == ", ..p..))  # Exact p-value case
        )),
        label.x = 0.5,  # Adjust the x-position of the label
        label.y = 2.2,  # Set y-position of the label
        digits = 3, 
        size = 5.5,# Number of digits for correlation coefficients
        parse = TRUE    # Enable parsing for math notation
    ) +
    theme_universal(base_size = plot_basesize)


plot_05_supp03_v1.0

rm(temp1,temp2)



#####################################################################
####### Longitudinal titres in those with no events ##############
#####################################################################


## return list of individuals with an event 
event_pid <-
    IgG_blood$event_titres_df_IgG_Blood %>%
    pull(pid) %>%
    unique()

titres <-readRDS("data/blood_IgG_titres.RDS")
## titre profiles in those with no event
no_event_titres <- titres %>%
    filter(
        ! pid %in% event_pid
    )

### create "relative titres" compared to baseline titre variable by for each Antigen, 
no_event_titres_rel <- 
    no_event_titres %>%
    arrange(visit_date) %>%
    group_by(pid, Antigen) %>%
    mutate(rel_titre = titre - first(titre))


geo_mean_monthly_df <- no_event_titres_rel %>%
    filter(visit_date <= min(visit_date) + 410) %>%
    left_join(age) %>%
    mutate(month = floor_date(visit_date, unit = "month")) %>%
    group_by(month, Antigen, age_grp) %>%
    summarise(
        mean = mean(rel_titre, na.rm = T),
        geo_mean = log10(psych::geometric.mean(10^rel_titre, na.rm = TRUE)),
        .groups = "drop"
    )


# Plot longitduinal titres relative to baseline 
plot_05_main03_v1.0 <- no_event_titres_rel %>%
    filter(visit_date <= min(visit_date) + 410) %>%
    left_join(age) %>%
    ggplot(aes(
        x = visit_date,
        y = rel_titre,
        col = Antigen
    )) +
    geom_smooth(data = geo_mean_monthly_df, 
                aes(x = month, y = geo_mean, group = 1), 
                se = F,
                #inherit.aes = FALSE,
                color = "black", size = 0.8, linetype = "solid") +
    geom_point(alpha = 0.8) +
    scale_colour_manual(values = c("SpyCEP" = "#FDC086", "SpyAD" = "#d19c2f", "SLO" = "#386CB0", "GAC" = "#7FC97F", "DNAseB" = "#BEAED4")) +  
    geom_line(aes(y = rel_titre, x = visit_date, group = pid),color = "grey", size = 0.3, alpha = 0.5) +
    geom_smooth(data = geo_mean_monthly_df, 
                aes(x = month, y = geo_mean, group = 1), 
                se = F,
                #inherit.aes = FALSE,
                color = "black", size = 0.8, linetype = "solid", alpha = 0.5) +
    facet_grid(Antigen ~factor(age_grp)) +
    theme_minimal() +
    labs(
        y = "Absolute change from baseline IgG level (log10 RLU/mL)",
        x = "Date",
   #     title = "Longitudinal blood IgG profiles in participants with no microbiologically confirmed events, relative to baseline titre"
    ) +
    guides(col = "none") +
    theme_universal(base_size = plot_basesize)

plot_05_main03_v1.0

n5 <- no_event_titres_rel %>% pull(pid) %>% unique() %>% length()

sprintf("Longitudinal blood IgG profiles in participants (n=%i) without microbiologically confirmed Strep A events during the study period.", n5)


 
######################
######################


# Upset plots to visualse the events distributions and visualise differences between 
# protection focused events and response focused events 


# Load datframe defining response focused events 
load("data/all_events_long_immunology.Rdata")

# Load dataframe defining protection focused events: 
load("data/all_events_long_incidence_wgs.RData")

shelf(gridExtra,UpSetR, wesanderson, grid)

# Summarize the number of events for each category in the immunology dataset
immunology_summary <- all_events_long_immunology %>%
    summarise(
        total_events = sum(gas_event, na.rm = TRUE),
        total_carriage = sum(gas_carriage, na.rm = TRUE),
        total_disease = sum(gas_infection, na.rm = TRUE),
        skin_carriage = sum(gas_skin_carriage, na.rm = TRUE),
        throat_carriage = sum(gas_throat_carriage, na.rm = TRUE),
        pyoderma = sum(gas_pyoderma, na.rm = TRUE),
        pharyngitis = sum(gas_pharyngitis, na.rm = TRUE)
    )

# Summarize the number of events for each category in the incidence dataset
incidence_summary <- all_events_long_incidence_wgs %>%
    summarise(
        total_events = sum(gas_event, na.rm = TRUE),
        total_carriage = sum(gas_carriage, na.rm = TRUE),
        total_disease = sum(gas_infection, na.rm = TRUE),
        skin_carriage = sum(gas_skin_carriage, na.rm = TRUE),
        throat_carriage = sum(gas_throat_carriage, na.rm = TRUE),
        pyoderma = sum(gas_pyoderma, na.rm = TRUE),
        pharyngitis = sum(gas_pharyngitis, na.rm = TRUE)
    )

# View the summaries
print(immunology_summary)
print(incidence_summary)


# Create unique identifiers for each visit
all_events_long_immunology <- all_events_long_immunology %>%
    mutate(unique_id = paste(pid, date, sep = "_"))


all_events_long_incidence_wgs <- all_events_long_incidence_wgs %>%
    mutate(unique_id = paste(pid, date, sep = "_"))

# Ensure the selected data frame contains only binary columns
upset_data_immunology <- all_events_long_immunology %>%
    select(gas_carriage, gas_infection, gas_skin_carriage, gas_throat_carriage, gas_pyoderma, gas_pharyngitis)

# Rename columns for clarity
colnames(upset_data_immunology) <- c("Carriage", "Disease", "Skin Carriage", "Throat Carriage", "Skin Disease", "Throat Disease")

# Convert numeric columns to integer (0/1)
upset_data_immunology <- upset_data_immunology %>%
    mutate(across(everything(), ~ as.integer(.)))

#

# Convert tibble to a regular data frame
upset_data_immunology <- as.data.frame(upset_data_immunology)

# Set the output to a PNG file with 300 DPI
png("R_output/supp_05_fig04.png", width = 1020 /96, height = 530 / 96, units = "in", res = 300)


upset(
    upset_data_immunology, 
    nsets = 6, 
    order.by = "freq",
    sets = c("Carriage", "Disease", "Skin Carriage", "Throat Carriage", "Skin Disease", "Throat Disease"),
    sets.x.label = "Event type",
    keep.order = TRUE,
    text.scale = 2,
  #  set_size.show = T,
    main.bar.color = "lightgreen"
) 

grid.text("Response focused events",x = 0.65, y=0.95, gp=gpar(fontsize=20))

# Close the PNG device to save the file
dev.off()

# Ensure the selected data frame contains only binary columns
upset_data_incidence <- all_events_long_incidence_wgs %>%
    select(gas_carriage, gas_infection, gas_skin_carriage, gas_throat_carriage, gas_pyoderma, gas_pharyngitis)

# Rename columns for clarity
colnames(upset_data_incidence) <- c("Carriage", "Disease", "Skin Carriage", "Throat Carriage", "Skin Disease", "Throat Disease")

# Convert numeric columns to integer (0/1)
upset_data_incidence <- upset_data_incidence %>%
    mutate(across(everything(), ~ as.integer(.)))

#

# Convert tibble to a regular data frame
upset_data_incidence <- as.data.frame(upset_data_incidence)
# Set the output to a PNG file with 300 DPI
png("R_output/supp_05_fig05.png", width = 1020 /96, height = 530 / 96, units = "in", res = 300)

upset(upset_data_incidence, 
      nsets = 6, 
      order.by = "freq",
      sets = c("Carriage", "Disease", "Skin Carriage", "Throat Carriage", "Skin Disease", "Throat Disease"),
      sets.x.label = "Event type",
      keep.order = TRUE,
      text.scale = 2,
    #  set_size.show = T,
      main.bar.color = "blue")

grid.text("Protection focused events",x = 0.65, y=0.95, gp=gpar(fontsize=20))

# Close the PNG device to save the file
dev.off()


table(all_events_long_immunology$visit, all_events_long_immunology$gas_carriage)

####################




event_type_age <- IgG_blood$fold_changes_IgG_Blood %>% 
    filter(!is.na(pre_post)) %>%
    left_join(age) %>%
    select(age_grp, event_id, event_type) %>%
    unique() %>%
    ungroup() %>%
    select(age_grp, event_type) %>%
    gtsummary::tbl_summary(by = age_grp,
                           percent = "row") %>%
    add_overall()


doc <- read_docx()
# Add regression table to word document
doc <- add_table_to_doc(doc, event_type_age)

print(doc, target = "R_output/supp_reviewer_response_event_type_age.docx")





######################
####### fin  #########
######################

# keep only plot objects

# List all objects in the environment
all_objects <- ls()
# Identify objects that contain "plot" in their name
plot_objects <- grep("plot_", all_objects, value = TRUE, ignore.case = TRUE)
# Remove all objects that do not contain "plot" in their name
rm(list = setdiff(all_objects, plot_objects), envir = .GlobalEnv)
rm(all_objects,plot_objects,keep_plot_objects)
