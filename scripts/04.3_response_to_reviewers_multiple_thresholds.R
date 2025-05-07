

source("scripts/setup_environment.R")
source("scripts/load_functions.R")
base_size = 14
# Set output directory

output_dir <- "R_output/"


# load packages:

shelf(officer,flextable)



#### Step 1 load constant data frame: 

load("R_objects/follow_up_dates_incidence.RData")

# Find enrolment dates for each person
start_dates <- follow_up_dates_incidence %>%
    select(pid, entry_1)



df <- readRDS("R_objects/baseline_blood_no_disease_titres.RDS")
protective_threshold_df <- readRDS("R_objects/bloodIgG_protective_threshold_df.RDS")
protective_threshold <-protective_threshold_df




   path_to_titre.df <- "R_objects/blood_IgG_titres.RDS"
   sample = "Blood"
   class = "IgG"
   var_name = "titre"
   next_event_window = 45
   titre_breakpoint_df <- tibble(Antigen = c("SLO","SpyAD", "SpyCEP", "GAC","DNAseB"), transition_point = c(4.3,4.1,4.3,3,3))
   
  # protective_threshold_df <- readRDS("R_objects/bloodIgG_protective_threshold_df.RDS")
   
   # Read titre data and merge with start dates and age, adjust visit dates, and categorize titres
   fun_titres <- read.titres(path_to_titre.df, var_name) %>%
       left_join(start_dates) %>%
       left_join(age) %>%
       mutate(visit_date = as.numeric(ymd(as.character(visit_date))) - entry_1) %>%
       group_by(Antigen) %>%
       ungroup() %>%
       left_join(titre_breakpoint_df) %>%
       mutate(titre_below_threshold = ifelse(titre <= transition_point, titre, transition_point),
              titre_above_threshold = ifelse(titre > transition_point, titre - transition_point, 0)) %>%
       filter(Antigen %in% c("SLO", "SpyAD", "SpyCEP")) %>%
       left_join(protective_threshold_df) %>%
       mutate(above_protective_threshold = ifelse(titre > Threshold, 1, 0 )
       )
   
   antigen_df <- fun_titres %>%
       group_by(pid,visit_date, age_grp,sex,hhsize) %>%
       summarise(n_above = sum(above_protective_threshold))
   
   # Filter for specific antigen
   
  
   pos_incidence_zero <- readRDS("R_objects/SpyCATS_incidence_df.RDS")
   
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
   
   final_df_3  <-fun_df %>%
       left_join(antigen_df %>% 
                     rename(date = visit_date)) %>%
       group_by(pid) %>%
       fill(age_grp:n_above) %>%
       ungroup() %>%
       mutate(hid = substring(pid, 1, 3),
              n_above = factor(n_above))



final_df_3$n_above <- case_when(final_df_3$n_above == "3" ~ as.factor("2+"),
                                final_df_3$n_above == "2" ~ as.factor("2+"),
                                T ~ as.factor(final_df_3$n_above))


final_df_3$n_above <- relevel(final_df_3$n_above, ref = "1")


final_df_3$pid <- factor(final_df_3$pid)
final_df_3$hid <- factor(final_df_3$hid)

# Fit the logistic regression model
model_8.2 <- lme4::glmer(event_next_n ~ n_above +  age_grp + sex + hhsize + (1 |pid) + (1 |hid), data = final_df_3, family = binomial)

# Extract AIC
model_aic <- AIC(model_8.2)
print(paste("Model:event_next_n ~ n_above +  age_grp + sex + hhsize + (1 | pid) + (1 | hid) 
                AIC:"
            , AIC(model_8.2)))

tb8.2 <- model_8.2 %>%
    tbl_regression(exponentiate = TRUE) %>%
    bold_p(t = 0.05)

tb8.2


############### reduce to 33% protective threhsolds

estimate_protective_threshold_33 <- function(data_subset, antigen, next_event_window = 45) {
    
    
    # Read titre data and merge with start dates and age, adjust visit dates, and categorize titres
    fun_titres <- data_subset %>%
        left_join(start_dates) %>%
        left_join(age) %>%
        mutate(visit_date = as.numeric(ymd(as.character(visit_date))) - entry_1) %>%
        group_by(Antigen) %>%
        ungroup()
    
    # Filter for specific antigen
    antigen_df <- fun_titres %>% filter(Antigen == antigen & !is.na(pid)) 
    
    pos_incidence_zero <- readRDS("R_objects/SpyCATS_incidence_df.RDS")
    
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
    
    threshold_prob <- 0.67 * max_prob
    
    titre_seq <- seq(0, max(final_df_5$titre_above_threshold, na.rm = TRUE), length.out = 1000)
    probs_for_threshold <- predict(model, newdata = data.frame(
        titre_below_threshold = titre_breakpoint,
        titre_above_threshold = titre_seq), type = "response", re.form = NA)
    
    titre_50 <- titre_seq[which.min(abs(probs_for_threshold - threshold_prob))]
    
    return(titre_50 + titre_breakpoint)
}




cross_validate_threshold_33 <- function(full_data, antigen, k = 10) {
    set.seed(342)  # Reproducibility
    folds <- sample(rep(1:k, length.out = nrow(full_data)))
    
    threshold_estimates <- numeric(k)
    
    for (i in 1:k) {
        train_data <- full_data[folds != i, ]
        
        threshold_estimates[i] <- estimate_protective_threshold_33(
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
    res <- cross_validate_threshold_33(
        full_data = readRDS("R_objects/blood_IgG_titres.RDS"),
        antigen = ag
    )
    
    cv_results_list[[ag]] <- res
}

cv_thresholds_df_33 <- bind_rows(cv_results_list)
print(cv_thresholds_df_33)
# Exponentiate and format the thresholds
formatted_thresholds_33 <- cv_thresholds_df_33 %>%
    mutate(
        Mean_exp = 10^Mean_Threshold,
        CI_Lower_exp = 10^CI_Lower,
        CI_Upper_exp = 10^CI_Upper,
        sentence = sprintf("%.0f (CI %.0f–%.0f) RLU/mL", Mean_exp, CI_Lower_exp, CI_Upper_exp)
    )

formatted_thresholds_33

# Construct the sentence
final_sentence <- sprintf(
    "Thresholds for 33 percent protection were %s for SLO, %s for SpyAD, and %s for SpyCEP.",
    formatted_thresholds_33$sentence[formatted_thresholds_33$Antigen == "SLO"],
    formatted_thresholds_33$sentence[formatted_thresholds_33$Antigen == "SpyAD"],
    formatted_thresholds_33$sentence[formatted_thresholds_33$Antigen == "SpyCEP"]
)

final_sentence

protective_threshold_df_33 <- 
    formatted_thresholds_33 %>%
    select(Antigen, Threshold = Mean_Threshold)




path_to_titre.df <- "R_objects/blood_IgG_titres.RDS"
sample = "Blood"
class = "IgG"
var_name = "titre"

titre_breakpoint_df <- tibble(Antigen = c("SLO","SpyAD", "SpyCEP", "GAC","DNAseB"), transition_point = c(4.3,4.1,4.3,3,3))

protective_threshold_df_33

# Read titre data and merge with start dates and age, adjust visit dates, and categorize titres
fun_titres <- read.titres(path_to_titre.df, var_name) %>%
    left_join(start_dates) %>%
    left_join(age) %>%
    mutate(visit_date = as.numeric(ymd(as.character(visit_date))) - entry_1) %>%
    group_by(Antigen) %>%
    ungroup() %>%
    left_join(titre_breakpoint_df) %>%
    mutate(titre_below_threshold = ifelse(titre <= transition_point, titre, transition_point),
           titre_above_threshold = ifelse(titre > transition_point, titre - transition_point, 0)) %>%
    filter(Antigen %in% c("SLO", "SpyAD", "SpyCEP")) %>%
    left_join(protective_threshold_df_33) %>%
    mutate(above_protective_threshold = ifelse(titre > Threshold, 1, 0 )
    )

antigen_df <- fun_titres %>%
    group_by(pid,hid,visit_date, age_grp,sex,hhsize) %>%
    summarise(n_above = sum(above_protective_threshold))

# Filter for specific antigen


pos_incidence_zero <- readRDS("R_objects/SpyCATS_incidence_df.RDS")

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

final_df_3  <-fun_df %>%
    left_join(antigen_df %>% 
                  rename(date = visit_date)) %>%
    group_by(pid) %>%
    fill(age_grp:n_above) %>%
    ungroup() %>%
    mutate(hid = substring(pid, 1, 3),
           n_above = factor(n_above))





final_df_3$n_above <- case_when(final_df_3$n_above == "3" ~ as.factor("2+"),
                                final_df_3$n_above == "2" ~ as.factor("2+"),
                                T ~ as.factor(final_df_3$n_above))
levels(final_df_3$n_above)

final_df_3$n_above <- relevel(final_df_3$n_above, ref = "1")


# Fit the logistic regression model
model_8.3 <- lme4::glmer(event_next_n ~ n_above +  age_grp + sex + hhsize + (1 | pid) + (1 | hid), data = final_df_3, family = binomial)

# Extract AIC
model_aic <- AIC(model_8.3)
print(paste("Model:event_next_n ~ n_above +  age_grp + sex + hhsize + (1 | pid) + (1 | hid) 
                AIC:"
            , AIC(model_8.3)))

tb8.3 <- model_8.3 %>%
    tbl_regression(exponentiate = TRUE) %>%
    bold_p(t = 0.05)

tb8.3

tbl_merge(tbls = list(tb8.2, tb8.3))

tb8.2_nabove <- model_8.2 %>%
    tbl_regression(exponentiate = TRUE, include = n_above, label = list(n_above = "Number of IgG levels above protective threshold")) %>%
    bold_p(t = 0.05) %>%
    modify_table_styling(columns = "label", 
                         label = "")

tb8.2_nabove

tb8.3_nabove <- model_8.3 %>%
    tbl_regression(exponentiate = TRUE, include = n_above,  label = list(n_above = "Number of IgG levels above protective threshold")) %>%
    bold_p(t = 0.05) %>%
    modify_table_styling(columns = "label", 
                         label = "")

tb8.3_nabove

tbl8.4 <- tbl_merge(tbls = list(tb8.2_nabove, tb8.3_nabove),
          tab_spanner = c("**50% threshold**", "**33% threshold**"))

tbl8.4



######################################################################
######################################################################
######################################################################


# Generate n_above counts from titres
n_above_summary <- fun_titres %>%
    filter(Antigen %in% c("SLO", "SpyAD", "SpyCEP")) %>%
    left_join(protective_threshold_df) %>%
    mutate(above_protective_threshold = ifelse(titre > Threshold, 1, 0)) %>%
    group_by(pid, visit_date) %>%
    summarise(n_above_50 = sum(above_protective_threshold)) %>%
    ungroup()

n_above_summary_33 <- fun_titres %>%
    filter(Antigen %in% c("SLO", "SpyAD", "SpyCEP")) %>%
    left_join(protective_threshold_df_33) %>%
    mutate(above_protective_threshold = ifelse(titre > Threshold, 1, 0)) %>%
    group_by(pid, visit_date) %>%
    summarise(n_above_33 = sum(above_protective_threshold)) %>%
    ungroup()




create_survival_dataframe <- function(titre.df, var = "titre", split = 4) {
    
    
    titre_breakpoint_df <- tibble(Antigen = c("SLO","SpyAD", "SpyCEP"), titre_breakpoint = c(4.3,4.1,4.3))
    
    fun_titre <- {
        if (var == "log_RLU_titres") {
            readRDS(titre.df) %>%
                rename(titre = log_RLU_titres,
                       visit_date = Date) %>%
                select(-log_AUC, -RLU, -AUC)
            
        } 
        else if (var == "AUC") {
            readRDS(titre.df) %>%
                rename(titre = log_AUC,
                       visit_date = Date) %>%
                select(-log_RLU_titres, -RLU, -AUC)
        } 
        else {
            readRDS(titre.df) 
        }
    }
    
    fun_titre <-
        fun_titre %>%
        left_join(start_dates) %>%
        left_join(age) %>%
        mutate(visit_date = as.numeric(ymd(as.character(visit_date))) - entry_1) %>%
        group_by(Antigen) %>%
        mutate(titre_cat = ntile(titre, split)) %>%
        ungroup() %>%
        left_join(titre_breakpoint_df) %>%
        mutate(hid = substring(pid, 1, 3)) %>%
        mutate(titre_below_threshold = ifelse(titre <= titre_breakpoint, titre, titre_breakpoint),
               titre_above_threshold = ifelse(titre > titre_breakpoint, titre - titre_breakpoint, 0))
    
    
    n_above_summary <- fun_titre %>%
        filter(Antigen %in% c("SLO", "SpyAD", "SpyCEP")) %>%
        left_join(protective_threshold_df) %>%
        mutate(above_protective_threshold = ifelse(titre > Threshold, 1, 0)) %>%
        group_by(pid, visit_date) %>%
        summarise(n_above_50 = sum(above_protective_threshold)) %>%
        ungroup()
    
    n_above_summary_33 <- fun_titre %>%
        filter(Antigen %in% c("SLO", "SpyAD", "SpyCEP")) %>%
        left_join(protective_threshold_df_33) %>%
        mutate(above_protective_threshold = ifelse(titre > Threshold, 1, 0)) %>%
        group_by(pid, visit_date) %>%
        summarise(n_above_33 = sum(above_protective_threshold)) %>%
        ungroup()
    
    # incidence zero 
    
    # Get start and stop dates from f/u periods. Set all start dates at zero.
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
    
    # Join with demographics
    incidence_zero <- right_join(age,incidence_zero)
    
    
    #  incidence_zero[c(1,4,6,7,18)] %>% colnames()
    
    # Remove old household size variables as about to add new ones
    incidence_zero <- incidence_zero %>%
        select(!contains("hhsize_"))
    
    
    # load zero datad dynamic household size dataframe 
    hh_size_long_zero <- readRDS("data/edited/hh_size_long_zero.RDS")
    
    
    # Create df for positive events incidence and zero dates (to be dependent variables)
    
    pos_incidence_zero <- readRDS("R_objects/SpyCATS_incidence_df.RDS")
    
    
    
    
    
    
    
    
    
    
    slo <- fun_titre %>% filter(Antigen == "SLO" & !is.na(pid)) 
    
    
    spyad <- fun_titre %>% filter(Antigen == "SpyAD" & !is.na(pid)) 
    
    spycep <- fun_titre %>% filter(Antigen == "SpyCEP" & !is.na(pid)) 
    
    incidence_data_long <- tmerge(incidence_zero, incidence_zero, id=pid, 
                                  tstart = start, tstop = stop)
    incidence_data_long <- tmerge(incidence_data_long, incidence_entry_exit_long_zero, id=pid,
                                  gap = tdc(timepoint, gap_status))
    incidence_data_long <- tmerge(incidence_data_long, hh_size_long_zero, id=pid,
                                  hhsize = tdc(date, hhsize),
                                  rain = tdc(date, rain))
    
    
    incidence_data_long <- tmerge(incidence_data_long, slo, id=pid,
                                  slo = tdc(visit_date, titre),
                                  slo_cat = tdc(visit_date, titre_cat),
                                  slo_above_threshold = tdc(visit_date, titre_above_threshold),
                                  slo_below_threshold = tdc(visit_date, titre_below_threshold))
    incidence_data_long <- tmerge(incidence_data_long, spyad, id=pid,
                                  spyad = tdc(visit_date, titre),
                                  spyad_cat = tdc(visit_date, titre_cat),
                                  spyad_above_threshold = tdc(visit_date, titre_above_threshold),
                                  spyad_below_threshold = tdc(visit_date, titre_below_threshold))
    incidence_data_long <- tmerge(incidence_data_long, spycep, id=pid,
                                  spycep = tdc(visit_date, titre),
                                  spycep_cat = tdc(visit_date, titre_cat),
                                  spycep_above_threshold = tdc(visit_date, titre_above_threshold),
                                  sypcep_below_threshold = tdc(visit_date, titre_below_threshold))
    incidence_data_long <- tmerge(incidence_data_long, pos_incidence_zero, id=pid,
                                  gas_pharyngitis = event(date, gas_pharyngitis),
                                  gas_pyoderma = event(date, gas_pyoderma),
                                  gas_throat_carriage = event(date, gas_throat_carriage),
                                  gas_skin_carriage = event(date, gas_skin_carriage),
                                  pharyngitis = event(date, pharyngitis),
                                  pyoderma = event(date, pyoderma),
                                  pcr_pharyngitis = event(date, pcr_pharyngitis),
                                  pcr_pyoderma = event(date, pcr_pyoderma),
                                  pcr_infection = event(date,pcr_infection),
                                  gas_event = event(date, gas_event), 
                                  gas_infection = event(date,gas_infection),
                                  gas_carriage = event(date,gas_carriage))
    
    
    incidence_data_long <- tmerge(incidence_data_long, n_above_summary, id = pid,
                                  n_above_50 = tdc(visit_date, n_above_50))
    
    incidence_data_long <- tmerge(incidence_data_long, n_above_summary_33, id = pid,
                                  n_above_33 = tdc(visit_date, n_above_33))
                                  
    
    incidence_data_long2 <- incidence_data_long %>%
        filter(gap == 0) %>%
        left_join(age) %>%
        mutate(hid = substring(pid, 1, 3))
    
    return(incidence_data_long2)
    
}

load("R_objects/follow_up_dates_incidence.RData")

# Find enrolment dates for each person
start_dates <- follow_up_dates_incidence %>%
    select(pid, entry_1)

survivial_df_blood_IgG <- create_survival_dataframe("R_objects/blood_IgG_titres.RDS",split = 4)

# Merge into survival dataframe


colnames(survivial_df_blood_IgG
         )


library(dplyr)

# Summarise person-time for each n_above_50 category
person_time_50 <- survivial_df_blood_IgG %>%
    filter(!is.na(n_above_50)) %>%
    group_by(n_above_50) %>%
    summarise(
        person_time_days = sum(tstop - tstart, na.rm = TRUE),
        person_time_years = person_time_days / 365.25,
        .groups = "drop"
    )


person_time_50
# Summarise person-time for each n_above_33 category
person_time_33 <- survivial_df_blood_IgG %>%
    filter(!is.na(n_above_33)) %>%
    group_by(n_above_33) %>%
    summarise(
        person_time_days = sum(tstop - tstart, na.rm = TRUE),
        person_time_years = person_time_days / 365.25,
        .groups = "drop"
    )



# For 50% threshold
event_summary_50 <- survivial_df_blood_IgG %>%
    filter(!is.na(n_above_50)) %>%
    group_by(n_above_50) %>%
    summarise(
        person_time_days = sum(tstop - tstart, na.rm = TRUE),
        person_time_years = person_time_days / 365.25,
        events = sum(gas_event == 1, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    rename(n_above = n_above_50) %>%
    mutate(threshold = "50%")


event_summary_50
# For 33% threshold
event_summary_33 <- survivial_df_blood_IgG %>%
    filter(!is.na(n_above_33)) %>%
    group_by(n_above_33) %>%
    summarise(
        person_time_days = sum(tstop - tstart, na.rm = TRUE),
        person_time_years = person_time_days / 365.25,
        events = sum(gas_event == 1, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    rename(n_above = n_above_33) %>%
    mutate(threshold = "33%")

event_summary_33

# Combine
combined_summary <- bind_rows(event_summary_50, event_summary_33) %>%
    select(threshold, n_above, person_time_years, events)


combined_summary


# Collapse function
collapse_n_above <- function(df) {
    df %>%
        mutate(n_above_collapsed = case_when(
            n_above %in% c(0, 1) ~ as.character(n_above),
            n_above >= 2 ~ "2+"
        )) %>%
        group_by(threshold, n_above_collapsed) %>%
        summarise(
            person_time_years = sum(person_time_years, na.rm = TRUE),
            events = sum(events, na.rm = TRUE),
            .groups = "drop"
        )
}

# Apply to previously created combined_summary
collapsed_summary <- collapse_n_above(combined_summary)


collapsed_summary 




# 1. Collapse summaries
collapsed_50 <- collapse_n_above(event_summary_50) %>% select(-threshold)
collapsed_33 <- collapse_n_above(event_summary_33) %>% select(-threshold)


collapsed_33





######## combine the table: 

# Convert gtsummary tables to data frames
df_50 <- as_tibble(tb8.2_nabove$table_body)
df_33 <- as_tibble(tb8.3_nabove$table_body)

# Add n_above labels to align with your summary tables
df_50$n_above_collapsed <- gsub(".*n_above", "", df_50$label)
df_33$n_above_collapsed <- gsub(".*n_above", "", df_33$label)

# Join with person-time and event data
df_50 <- left_join(df_50, collapsed_50, by = "n_above_collapsed")
df_33 <- left_join(df_33, collapsed_33, by = "n_above_collapsed")

# Format a combined table
shelf(gt, officer)

combined_table <- bind_rows(
    df_50 %>% mutate(Threshold = "50%"),
    df_33 %>% mutate(Threshold = "33%")
) %>%
    select(Threshold, n_above_collapsed,person_time_years, events, estimate, conf.low, conf.high, p.value
           ) %>%
    mutate(across(estimate:p.value, ~ round(., 2))) %>%
    gt() %>%
    cols_label(
        n_above_collapsed = "Threshold group",
        person_time_years = "Person-years",
        events = "Events",
        estimate = "OR",
        conf.low = "CI low",
        conf.high = "CI high",
        p.value = "p-value",

    ) %>%
    tab_spanner(label = "Follow-up Summary", columns = person_time_years:events) %>%
    tab_spanner(label = "Model Results", columns = estimate:p.value)


combined_table

class(combined_table)
 
# Save gt table to a Word document
gtsave(data = combined_table, filename = "R_output/combined_thresholds_combined_table_output.docx")

