


################################################################
################################################################
#################### create events dataframes ##################

## For each compartment (IgG blood, IgA oral fluid, IgG oral fluid) creates a labelled list object with the following dataframes: 

# titres dataframe - the input dataframe 
# event_titres dataframe: titres with their relationship in days to an event, with identification of the pre-event, event and post event measurement for the individual around each event. 
# fold_change: absolute titre changes and fold changes in antibody titre around events 
# titres_relative_to_baseline: titres in relation to events normalized to the pre event titre
# relative_changes: absolute titre changes from baseline around events 



library(librarian)

shelf(tidyverse,survival,gtsummary)

# Setup environment
source("scripts/setup_environment.R")
#### load the events dataframe 

load(file = "data/all_events_long_immunology.Rdata")


### read titres dynamically

# Function to read and standardize the titres dataframe
read_titres <- function(titre_input, var = "titre", class = NULL, sample = NULL) {
    if (is.character(titre_input)) {
        titre_df <- readRDS(titre_input)
    } else if (is.data.frame(titre_input)) {
        titre_df <- titre_input
    } else {
        stop("Input must be either a file path or a dataframe.")
    }
    
    fun_titre <- {
        if (var == "log_RLU_titres") {
            titre_df %>%
                rename(titre = log_RLU_titres, visit_date = Date) %>%
                select(-log_AUC, -RLU, -AUC)
        } else if (var == "AUC") {
            titre_df %>%
                rename(titre = log_AUC, visit_date = Date) %>%
                select(-log_RLU_titres, -RLU, -AUC)
        } else if (var != "titre") {
            titre_df %>%
                rename(titre = !!sym(var))
        } else {
            titre_df
        }
    }
    
    if (!"class" %in% colnames(titre_df) & !is.null(class)) {
        fun_titre <- fun_titre %>% mutate(class = class)
    }
    if (!"sample" %in% colnames(fun_titre) & !is.null(sample)) {
        fun_titre <- fun_titre %>% mutate(sample = sample)
    }
    
    return(fun_titre)
}

# Function to process events dataframe


process_events <- function(events_df) {
    events_df %>%
        filter(gas_event == 1) %>%
        arrange(pid, date) %>%
        group_by(pid) %>%
        mutate(event_no_per_pid = row_number()) %>%
        mutate(event_id = paste(pid, event_no_per_pid, sep = "_")) %>%
        mutate(skin = case_when(
            gas_pyoderma == 1 | gas_skin_carriage == 1 ~ 1,
            TRUE ~ 0),
            throat = case_when(
                gas_pharyngitis == 1 | gas_throat_carriage == 1 ~ 1,
                TRUE ~ 0),
            site = case_when(
                skin == 1 & throat == 0 ~ "skin",
                skin == 0 & throat == 1 ~ "throat",
                skin == 1 & throat == 1 ~ "mixed"
            ),
               event_type = 
                   case_when(
                       gas_pyoderma ==1 & gas_pharyngitis == 0 ~ "skin disease",
                       gas_infection ==0 & site == "skin" ~ "skin carriage",
                       gas_pyoderma == 0 & gas_pharyngitis == 1 ~ "throat disease",
                       gas_infection ==0 & site == "throat" ~ "throat carriage",
                       T ~ "other"
                   )
            #    ,
           # event_type = 
           #     case_when(
           #         gas_infection ==1 & site == "skin" ~ "skin disease",
           #         gas_infection ==0 & site == "skin" ~ "skin carriage",
           #         gas_infection ==1 & site == "throat" ~ "throat disease",
           #         gas_infection ==0 & site == "throat" ~ "throat carriage",
           #         T ~ "other"
           #     )
            ) %>%
        select(pid, event_visit = visit, event_date = date, event_no_per_pid, event_id, site, gas_infection, event_type)
}



create_visits <- function(titres_df) {
    visits <- titres_df %>%
        select(pid, visit_date, Antigen, titre) %>%
        mutate(Antigen = factor(Antigen, levels = c("GAC", "SLO", "SpyAD", "SpyCEP", "DNAseB"))) %>%
        spread(Antigen, titre) 
    
    unique_visits <- visits %>%
        mutate(id = paste0(pid, visit_date)) %>%
        select(id) %>%
        unique()
    
    return(list(visits = visits, unique_visits_dim = dim(unique_visits)))
}

##### This function:
# 1 takes titres in idividuals are allocates them a relationship to an event in the study
# 2 creates a variable to identify if the sample is the pre sample, the event sample or the post event sample 



process_event_visits <- function(events_df, visits_df,sample) {
    events_visits <- left_join(events_df, visits_df, by = "pid", relationship = "many-to-many") %>%
        arrange(pid, event_date, visit_date) %>%
        group_by(event_id) %>%
        
        
        mutate(diff = as.numeric(difftime(visit_date, event_date, units = "days")),
               post_event = case_when(
                   diff >= 1 ~ case_when(
                       if (sample == "OF") {
                           diff == min(diff[diff >= 1], na.rm = TRUE) ~ 1
                       } else {
                           diff == min(diff[diff >= 14], na.rm = TRUE) ~ 1
                       },
                       TRUE ~ 0
                   ),
                   TRUE ~ 0
               ),
               pre_event = case_when(
                   
                   if (sample == "OF") {
                       
                       diff <= -7 ~ case_when(
                           diff == max(diff[diff <= -7], na.rm = TRUE) ~ 1,
                           TRUE ~ 0
                       )
                       
                       
                   } else {
                       diff <= -14 ~ case_when(
                           diff == max(diff[diff <= -14], na.rm = TRUE) ~ 1,
                           TRUE ~ 0
                       )
                       
                   },
                   TRUE ~ 0
                   
                   
                   
               ),
               event = case_when(event_date == visit_date ~ 1, TRUE ~ 0))
    
    
    
    #
    ### take the titre dataframe and make it wide 
    
    event_titres <- events_visits %>%
        pivot_longer(cols = GAC:DNAseB, names_to = "Antigen", values_to = "titre") %>%
        mutate(Antigen = str_remove(Antigen, "_log")) %>%
        mutate(status = case_when(
            pre_event == 1 ~ "Pre",
            post_event == 1 ~ "Post",
            event == 1 ~ "Event"
        )) %>%
        select(pid, event_id, event_no_per_pid, visit_date, diff, status, Antigen, titre, site, gas_infection, event_type) %>%
        mutate(Antigen = factor(Antigen, levels = c("GAC", "SLO", "SpyAD", "SpyCEP", "DNAseB")))
    
    return(event_titres)
}

# Function to create event-specific dataframes

create_event_specific_dataframes <- function(event_titres) {
    df_filtered <- event_titres %>%
        filter(!is.na(status))
    
    df.rate <- df_filtered %>%
        arrange(event_id, Antigen, diff) %>%
        select(-visit_date) %>%
        group_by(event_id, Antigen) %>%
        pivot_wider(names_from = status, values_from = c(titre, diff)) %>%
        mutate(#pre_event_rate = (titre_Event - titre_Pre) / abs(diff_Pre),
            #event_post_rate = (titre_Post - titre_Event) / abs(diff_Post),
            #pre_post_rate = (titre_Post - titre_Pre) / (abs(diff_Pre) + abs(diff_Post)),
            pre_post = titre_Post - titre_Pre,
            event_post = titre_Post - titre_Event,
            pre_event = titre_Event - titre_Pre,
            pre_post_fold = (10^titre_Post / 10^titre_Pre),
            event_post_fold = (10^titre_Post / 10^titre_Event),
            pre_event_fold = (10^titre_Event / 10^titre_Pre),
            two_fold_increase = ifelse(pre_post_fold > 2, "Yes", "No")) %>%
        ungroup() %>%
        mutate(Antigen = factor(Antigen, levels = c("GAC", "SLO", "SpyAD", "SpyCEP", "DNAseB")))
    
    rel_titres2 <- event_titres %>%
        filter(!is.na(status)) %>%
        group_by(pid, Antigen, event_id) %>%
        filter(sum(status == "Pre") == 1) %>%
        mutate(pre_value = first(titre[status == "Pre"])) %>%
        mutate(rel_titre = titre - pre_value) %>%
        select(-titre, -pre_value) %>%
        ungroup() %>%
        mutate(Antigen = factor(Antigen, levels = c("GAC", "SLO", "SpyAD", "SpyCEP", "DNAseB")))
    
    rel2_df_filtered <- rel_titres2 %>%
        filter(!is.na(status)) %>%
        mutate(Antigen = factor(Antigen, levels = c("GAC", "SLO", "SpyAD", "SpyCEP", "DNAseB")))
    
    rel2_df.rate <- rel2_df_filtered %>%
        arrange(event_id, Antigen, diff) %>%
        select(-visit_date) %>%
        group_by(event_id, Antigen) %>%
        pivot_wider(names_from = status, values_from = c(rel_titre, diff)) %>%
        mutate(pre_post = rel_titre_Post - rel_titre_Pre,
               event_post = rel_titre_Post - rel_titre_Event,
               pre_event = rel_titre_Event - rel_titre_Pre) %>%
        ungroup() %>%
        mutate(Antigen = factor(Antigen, levels = c("GAC", "SLO", "SpyAD", "SpyCEP", "DNAseB")))
    
    return(list(fold_changes = df.rate, titres_relative_to_baseline = rel_titres2, relative_changes = rel2_df.rate))
}

# Main function to create events dataframe and tag output variables
create_events_dataframe <- function(titre_input, var = "titre", class = NULL, sample = NULL) {
    
    titres_df <- read_titres(titre_input, var, class, sample)
    
    events_df <- process_events(all_events_long_immunology)
    
    visits_result <- create_visits(titres_df)
    
    visits_df <- visits_result$visits
    
    event_titres_df <- process_event_visits(events_df, visits_df,sample)
    
    event_specific_dataframes <- create_event_specific_dataframes(event_titres_df)
    
    # Add event_titres_df to the list before tagging
    event_specific_dataframes$event_titres_df <- event_titres_df
    event_specific_dataframes$titres_df <-  titres_df
    
    class_sample_tag <- paste0(class, "_", sample)
    
    names(event_specific_dataframes) <- paste0(names(event_specific_dataframes), "_", class_sample_tag)
    
    return(event_specific_dataframes)
}





#######################################################################
#### plot drawing functions to analyse paired titres around event ####

create_paired_dataframe <- function(input_df, include_event_as_paired = T){
    
    
    paired_df <- 
        
        if (include_event_as_paired) {
            
            input_df %>%
                as.data.frame() %>%
                select(pid, event_id, Antigen, status, titre) %>%
                group_by(event_id, Antigen) %>%
                filter(
                    (any(status == "Pre") & any(status == "Post")) |
                        (any(status == "Pre") & any(status == "Event"))
                ) %>%
                filter(
                    !is.na(status), 
                    status %in% c("Pre", "Post", "Event")
                ) %>%
                mutate(status = if_else(status == "Event" & !any(status == "Post"), "Post", status)) %>%
                filter(status == "Pre" | status == "Post") %>%
                ungroup()
            
        } else {
            
            input_df %>%
                as.data.frame() %>%
                select(pid,event_id, Antigen, status, titre) %>%
                group_by(event_id, Antigen) %>%
                filter(any(status == "Pre") & any(status == "Post")) %>%
                filter(!is.na(status)) %>%
                filter(status == "Pre" | status == "Post") %>%
                ungroup()
            
        }
    
    return(paired_df)
    
}


generate_paired_analysis <- function(input_df, include_event_as_paired = T, sample = "", class = "") {
    # Create a paired dataframe
    
    paired_df <- create_paired_dataframe(input_df, include_event_as_paired = include_event_as_paired)
    # Perform paired Wilcoxon test within each Antigen group
    result <- paired_df %>%
        group_by(Antigen) %>%
        summarise(
            mean_pre = log10(psych::geometric.mean(10^titre[status == "Pre"])),
            mean_post = log10(psych::geometric.mean(10^titre[status == "Post"])),
            p = {
                pre_titre <- titre[status == "Pre"]
                post_titre <- titre[status == "Post"]
                wilcox.test(pre_titre, post_titre, paired = TRUE)$p.value
            }
        ) %>%
        mutate(
            mean_difference = mean_post - mean_pre
        )
    
    # Adjust p-values using FDR correction
    result <- result %>%
        mutate(p_adjusted = p.adjust(p, method = "fdr"))
    
    # Create annotations for p-values
    create_annotations <- function(pvals) {
        annotations <- pvals %>%
            mutate(label = case_when(p_adjusted < 0.0001 ~ "p.adj < 0.0001",
                                     p_adjusted > 0.25 ~ "p.adj > 0.25",
                                     T ~ sprintf("p.adj = %.4f", p_adjusted))) %>%
            select(Antigen, label)
        
        return(annotations)
    }
    
    annotations <- create_annotations(result)
    
    # Generate the plot
    plot <- paired_df %>%
        ggplot(aes(x = factor(status, level = c("Pre", "Post")), y = titre)) +
        geom_violin(aes(fill = Antigen), alpha = 0.1) +
        geom_line(aes(group = event_id), col = "grey", linewidth = 0.075) +
        geom_point(size = dot_size, aes(col = Antigen), width = 0.2, alpha = 0.8) +
        stat_summary(fun = median, geom = "crossbar", width = 0.5, col = "black") +
        scale_colour_manual(values = c("SpyCEP" = "#FDC086", "SpyAD" = "#d19c2f", "SLO" = "#386CB0", "GAC" = "#7FC97F", "DNAseB" = "#BEAED4")) +
        scale_fill_manual(values = c("SpyCEP" = "#FDC086", "SpyAD" = "#d19c2f", "SLO" = "#386CB0", "GAC" = "#7FC97F", "DNAseB" = "#BEAED4")) +
        facet_wrap(~Antigen, nrow = 1) +
        theme_minimal() +
        theme(legend.position = "none") +
        labs(
            #   title = paste(sample, class, "antibody titres in paired samples before and after any Strep A event", sep = " "),
            x = "Relationship to event",
            y = "Titre (RLU/mL)"
        ) +
        scale_y_continuous(limits = c(0, 6), breaks = c(1, 2, 3, 4, 5, 6), labels = log10_to_exp) +
        geom_text(data = annotations, aes(x = 1.5, y = 6, label = label), color = "black", size = plot_basesize/2.8, inherit.aes = FALSE) +
        theme_universal(base_size = plot_basesize)
    
    return(list(results = result, plot = plot))
}




#######################################################################
#### plot drawing functions to analyse titre changes around events ####

# plot by age group 


Longitudinal_event_plot_by_age <- function(data, sample, class, var_label = "titre", Ag = c("GAC", "SLO" , "SpyAD",  "SpyCEP","DNAseB"), age_cut = 100, simplify = F , diff_n = 43, centile = F, normal_dist = T, zscore = F) {
    
    
    titre_label <- {
        
        if (centile == T) {
            "age-specific centile score"
        }
        
        else if (zscore == T) {
            "Z score"
        }
        
        else if (var_label == "AUC") {
            "log10 AUC"
        }
        
        else {
            "log10 RLU"}
    }
    
    data_modified <- {
        if(simplify) {
            data %>%
                left_join(age) %>%
                mutate(
                    age_grp =
                        case_when(
                            age >= 0 & age < 12 ~ "0-11 years",
                            age >= 12 ~ "12+ years"
                        ),
                    age_grp = factor(age_grp, levels = c("0-11 years", "12+ years"))
                )
        } else {
            data %>%
                left_join(age)
        }
    }
    
    
    bootstrap_ci <- function(data, n_bootstrap = 1000, alpha = 0.05) {
        # Generate bootstrap samples and compute the median for each sample
        boot_medians <- map_dbl(1:n_bootstrap, function(i) {
            sample_data <- sample(data$rel_titre, replace = TRUE)
            median(sample_data, na.rm = TRUE)
        })
        
        # Calculate the lower and upper bounds of the confidence interval
        ci_lower <- quantile(boot_medians, probs = alpha / 2)
        ci_upper <- quantile(boot_medians, probs = 1 - alpha / 2)
        
        list(median = median(boot_medians), LL = ci_lower, UL = ci_upper)
    }
    
    joining_df <- data_modified %>%
        left_join(age) %>%
        filter(diff > - diff_n & diff < diff_n,
               Antigen %in% Ag,
               age <= age_cut) %>%
        group_by(status, Antigen, age_grp) %>%
        do({
            boot_results <- bootstrap_ci(.)
            data.frame(
                median_rel_titre = boot_results$median,
                LM = boot_results$LL,
                UM = boot_results$UL
            )
        }) %>%
        ungroup()
    
    summary_df <- data_modified %>%
        left_join(age) %>%
        filter(diff > - diff_n & diff < diff_n,
               Antigen %in% Ag,
               age <= age_cut) %>%
        group_by(status, Antigen, age_grp) %>%
        mutate(
            mean_log_RLU_titres = mean(rel_titre, na.rm = TRUE),
            sd = sd(rel_titre, na.rm = TRUE),
            n = n(),
            SE = sd / sqrt(n),
            LL = mean_log_RLU_titres - (1.96 * SE),
            UL = mean_log_RLU_titres + (1.96 * SE),
            .groups = 'drop'
            
        ) %>%
        left_join(joining_df)
    
    
    
    plot1 <- summary_df %>%
        ggplot(aes(x = diff, y = median_rel_titre)) +
        facet_grid(Antigen ~ age_grp) +
        geom_ribbon(aes(ymin = LM, ymax = UM, fill = Antigen), alpha = 0.5) +
        geom_line(aes(y = rel_titre, x = diff, group = event_id),color = "grey", size = 0.3, alpha = 0.5) +
        geom_point(size = dot_size, width = 0.1, alpha = 0.8, aes(y = rel_titre, x = diff, col = Antigen)) +
        geom_line() +
              labs(x = "Days from event (n)",
             y = paste0("Change from baseline titre around event /n (",titre_label,")")
        ) + 
        guides(col = "none", fill = "none") +
        scale_colour_manual(values = c("SpyCEP" = "#FDC086", "SpyAD" = "#8B8000", "SLO" = "#386CB0", "GAC" = "#7FC97F", "DNAseB" = "#BEAED4")) +  
        scale_fill_manual(values = c("SpyCEP" = "#FDC086", "SpyAD" = "#8B8000", "SLO" = "#386CB0", "GAC" = "#7FC97F", "DNAseB" = "purple")) +  
        # ggtitle(paste0("Relative ", sample, " ", class, " antibody changes around Strep A events in the study")) +
        theme_minimal() +
        theme_universal(base_size = plot_basesize)
    
    
    
    plot2 <- summary_df %>%
        ggplot(aes(x = diff, y = mean_log_RLU_titres)) +
        facet_grid(Antigen ~ age_grp) +
        geom_ribbon(aes(ymin = LL, ymax = UL, fill = Antigen), alpha = 0.5) +
        geom_line(data = summary_df, aes(y = rel_titre, x = diff, group = event_id),color = "grey", size = 0.3, alpha = 0.5) +
        geom_point(size = dot_size,width = 0.1, alpha = 0.8, aes(y = rel_titre, x = diff, col = Antigen)) +
        geom_line() +
        labs(x = "Days from event (n)",
             y = paste0("Change from baseline titre around event /n (",titre_label,")")
        ) + 
        guides(col = "none", fill = "none") +
        scale_colour_manual(values = c("SpyCEP" = "#FDC086", "SpyAD" = "#8B8000", "SLO" = "#386CB0", "GAC" = "#7FC97F", "DNAseB" = "#BEAED4")) +  
        scale_fill_manual(values = c("SpyCEP" = "#FDC086", "SpyAD" = "#8B8000", "SLO" = "#386CB0", "GAC" = "#7FC97F", "DNAseB" = "purple")) +  
        # ggtitle(paste0("Relative ", sample, " ", class, " antibody changes around Strep A events in the study")) +
        theme_minimal() +
        theme_universal(base_size = plot_basesize)
    
    
    
    if (normal_dist) {
        return(plot2)
    }
    else {
        return(plot1)
    }
    
    
}


# plot by event type 


longitudinal_event_site_plot <- function(data, sample, class,var_label = "titre",Ag = c( "GAC", "SLO" , "SpyAD",  "SpyCEP","DNAseB"), age_cut = 100, diff_n = 43, centile = F, zscore = F, normal_dist =T) {
    
    
    age_label <- 
        {
            
            if (age_cut == 100) {
                "all participants"
            }
            
            else {
                paste("participants under", age_cut, sep = " ")}
        }
    
    
    titre_label <- 
        {
            
            if (centile == T) {
                "age-specific centile score"
            }
            
            
            else if (zscore == T) {
                "Z score"
            }
            
            else if (var_label == "AUC") {
                "log10 AUC"
            }
            
            else {
                "log10 RLU"}
        }
    
    
    
    
    ### This is the dataframe of titres relative to the pre event titre: 
    
    
    bootstrap_ci <- function(data, n_bootstrap = 1000, alpha = 0.05) {
        # Generate bootstrap samples and compute the median for each sample
        boot_medians <- map_dbl(1:n_bootstrap, function(i) {
            sample_data <- sample(data$rel_titre, replace = TRUE)
            median(sample_data, na.rm = TRUE)
        })
        
        # Calculate the lower and upper bounds of the confidence interval
        ci_lower <- quantile(boot_medians, probs = alpha / 2)
        ci_upper <- quantile(boot_medians, probs = 1 - alpha / 2)
        
        list(median = median(boot_medians), LL = ci_lower, UL = ci_upper)
    }
    
    joining_df <- data %>%
        filter(event_type != "other",
               !is.na(status)) %>% 
        left_join(age) %>%
        filter(diff > - diff_n & diff < diff_n,
               Antigen %in% Ag,
               age <= age_cut) %>%
        group_by(status, Antigen, event_type) %>%
        do({
            boot_results <- bootstrap_ci(.)
            data.frame(
                median_rel_titre = boot_results$median,
                LM = boot_results$LL,
                UM = boot_results$UL
            )
        }) %>%
        ungroup()
    
    summary_df <- data %>%
        filter(event_type != "other",
               !is.na(status)) %>% 
        # filter(event_no_per_pid == 1) %>%
        filter(diff > - diff_n & diff < diff_n) %>%
        left_join(age) %>%
        filter(Antigen %in% Ag,
               age <= age_cut) %>%
        group_by(status, Antigen, event_type) %>%
        mutate(
            mean_titre = mean(rel_titre, na.rm = TRUE),
            sd = sd(rel_titre, na.rm = TRUE),
            n = n(),
            SE = sd / sqrt(n),
            LL = mean_titre - (1.96 * SE),
            UL = mean_titre + (1.96 * SE)
        ) %>%
        left_join(joining_df)
    
    
    
    plot1 <- summary_df %>%
        ggplot(aes(x = diff, y = median_rel_titre)) +
        facet_grid(Antigen~event_type, scales = "free") +
        geom_ribbon(aes(ymin = LM, ymax = UM, fill = factor(event_type)), alpha = 0.8) +        geom_line(aes(y = rel_titre, x = diff, group = event_id),color = "grey", size = 0.3, alpha = 0.5) +
        geom_point(size = dot_size,width = 0.1, alpha = 0.5, 
                   aes(y = rel_titre, x = diff, fill = factor(event_type), col = factor(event_type), group = event_type), 
                   position = position_dodge(width = 5)) +  # Dodging the points by site
        geom_line() +
        labs(x = "Days from event (n)",
             y = paste0("Change from baseline titre around event /n(",titre_label,")")
        ) + 
        theme_minimal() +
        #  ggtitle(paste0("Relative ", sample, " ", class, " antibody changes around Strep A events in ", age_label, " by event type")) +
        scale_fill_manual(values = wes_palette("Zissou1")) +
        guides(fill = "none",
               col = "none") +
        scale_color_manual(values = wes_palette("Zissou1")) +
        theme_universal(base_size = plot_basesize)
    
    plot2 <- summary_df %>%
        ggplot(aes(x = diff, y = mean_titre)) +
        facet_grid(Antigen~event_type, scales = "free") +
        geom_ribbon(aes(ymin = LL, ymax = UL, fill = factor(event_type)), alpha = 0.8) +
  
        geom_line(aes(y = rel_titre, x = diff, group = event_id),color = "grey", size = 0.3, alpha = 0.5) +
        geom_point(size = dot_size,width = 0.1, alpha = 0.5, 
                   aes(y = rel_titre, x = diff, fill = factor(event_type), col = factor(event_type), group = event_type), 
                   position = position_dodge(width = 5)) +
        geom_line() +
        labs(x = "Days from event (n)",
             y = paste0("Change from baseline titre around event /n (",titre_label,")")
        ) + 
        theme_minimal() +
        #   ggtitle(paste0("Relative ", sample, " ", class, " antibody changes around Strep A events in ", age_label, " by event type")) +
        scale_fill_manual(values = wes_palette("Zissou1")) +
        guides(fill = "none",
               col = "none") +
        scale_color_manual(values = wes_palette("Zissou1")) +
        theme_universal(base_size = plot_basesize)
    
    
    
    if (normal_dist) {
        return(plot2)
    }
    else {
        return(plot1)
    }
    
}







# function to add analysis tables to word docs: 


add_table_to_doc <- function(doc, table, title = NULL) {
    
    
    
    ft <- as_flex_table(table)
    
    # Scale the table to fit the page
    ft <- set_table_properties(ft, layout = "autofit")
    
    if (!is.null(title)) {
        doc <- doc %>%
            body_add_par(title, style = "Table Caption") %>%
            body_add_flextable(ft) %>%
            body_add_par("", style = "Normal") %>%
            body_end_section_landscape()
        
    } else {
        doc <- doc %>%
            body_add_flextable(ft) %>%
            body_add_par("", style = "Normal") %>%
            body_end_section_landscape()
    }
    
    return(doc)
}


##### Format P values to 2 significant figures unless <0.0001 

format_p_value <- function(p) {
    ifelse(p < 0.0001, "<0.0001", formatC(p, format = "g", digits = 2))
}

