# Title: response to reviewer - impact of combination of thresholds 
# Version: 1.0
# Date: 2025-03-24
# Author: Alexander J Keeley

# Inputs: 
#   - data/blood_IgG_titres.RDS
#   - data/incidence_start_dates.RDS
#   - data/SpyCATS_incidence_df.RDS
#   - data/baseline_blood_no_disease_titres.RDS
#   - data/bloodIgG_protective_threshold_df.RDS
#   - data/survivial_df_blood_IgG_n_above.RDS

# Outputs: 
#   - Cross-validated protective titre thresholds for 33% protection
#   - Logistic regression models for 50% and 33% titre thresholds
#   - Merged results table with ORs, confidence intervals, and event counts
#   - Word document with summary tables of logistic model outputs and person-time by antibody profile

# Description:
#
# This script estimates and compares protective blood IgG thresholds against S pyogenes-related events. 
# Key analysis components include:
#
# 1. Loading and preprocessing titre and event data, including merging with demographics and anonymised dates.
# 2. Estimating antigen-specific 33% protection thresholds using logistic regression and 10-fold cross-validation.
# 3. Fitting logistic regression models with random effects (GLMM) to assess associations between protection level (number of titres above threshold) and event risk, using both 50% and 33% thresholds.
# 4. Summarising model performance using AIC, odds ratios, and p-values for number of antigens above threshold.
# 5. Merging model results with person-time and event counts across antibody threshold groups (0, 1, 2+).
# 6. Exporting final comparative results as a publication-ready table in Word format.

### Date Anonymisation for Public Data Sharing
#
# All dates are offset by a fixed constant to preserve relative timing while removing identifiable calendar dates. This protects participant confidentiality while maintaining the internal validity of the time-based analyses.



# Load environment and helper functions
source("scripts/setup_environment.R")
source("scripts/load_functions.R")

# Set output directory
output_dir <- "R_output/"

# load packages:
shelf(officer,flextable, lme4)


### LOAD DATA

# Load reference dates for incidence window alignment
start_dates <- readRDS("data/incidence_start_dates.RDS")

# Load previously defined thresholds for 50% protection
protective_threshold_df <- readRDS("data/bloodIgG_protective_threshold_df.RDS")
protective_threshold <-protective_threshold_df

# File path to raw titres
path_to_titre.df <- "data/blood_IgG_titres.RDS"
sample = "Blood"
class = "IgG"
var_name = "titre"
next_event_window = 45

# Set breakpoint titres to split slope in piecewise modeling
titre_breakpoint_df <- tibble(Antigen = c("SLO","SpyAD", "SpyCEP", "GAC","DNAseB"), transition_point = c(4.3,4.1,4.3,3,3))

### PROCESS TITRES DATA

# Preprocess titres: merge with metadata, calculate transformed variables
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

# Summarise number of antigens above protective threshold
antigen_df <- fun_titres %>%
    group_by(pid,visit_date, age_grp,sex,hhsize) %>%
    summarise(n_above = sum(above_protective_threshold))

### LOAD EVENT DATA

# Load incidence data and calculate future event status within time window
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

# Merge titre summary with event data
final_df_3  <-fun_df %>%
    left_join(antigen_df %>% 
                  rename(date = visit_date)) %>%
    group_by(pid) %>%
    fill(age_grp:n_above) %>%
    ungroup() %>%
    mutate(hid = substring(pid, 1, 3),
           n_above = factor(n_above))


# Collapse 2 and 3 into "2+" category
final_df_3$n_above <- case_when(final_df_3$n_above == "3" ~ as.factor("2+"),
                                final_df_3$n_above == "2" ~ as.factor("2+"),
                                T ~ as.factor(final_df_3$n_above))


final_df_3$n_above <- relevel(final_df_3$n_above, ref = "1")


final_df_3$pid <- factor(final_df_3$pid)
final_df_3$hid <- factor(final_df_3$hid)


### FIT LOGISTIC MODEL FOR 50% THRESHOLD

# Fit the logistic regression model - fully adjusted
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


### ESTIMATE 33% PROTECTIVE THRESHOLDS

# Function to estimate antigen-specific titre thresholds associated with 33% protection
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
    
    # Identify breakpoint titre to split slope
    titre_breakpoint <- titre_breakpoint_df %>%
        filter(Antigen == antigen) %>%
        pull(transition_point)
    
    # Merge titre with event data
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



# Function to apply k-fold cross-validation to estimate 33% protective threhsolds 

cross_validate_threshold_33 <- function(full_data, antigen, k = 10) {
    set.seed(342)  # Fix this random number for reproducibility
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
        full_data = readRDS("data/blood_IgG_titres.RDS"),
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


########## recreate event and titre dataframes with 33% protective thresholds 

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
    mutate(above_protective_threshold = ifelse(titre > Threshold, 1, 0 ),
           hid = substring(pid, 1, 3)
    )

antigen_df <- fun_titres %>%
    group_by(pid,hid,visit_date, age_grp,sex,hhsize) %>%
    summarise(n_above = sum(above_protective_threshold))

# Filter for specific antigen


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

### FIT LOGISTIC MODEL FOR 50% THRESHOLD

# Fit the logistic regression model - fully adjusted
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


#### 

tbl_merge(tbls = list(tb8.2, tb8.3))


#### Now just take the Number of IgG levels above protective threshold data for outputs. 

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


### describe time above n thresholds 


# Generate n_above counts from titres dataframe at each time of measurement
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

survivial_df_blood_IgG_n <- readRDS("data/survivial_df_blood_IgG_n_above.RDS")

# Summarise person-time for each n_above_50 category
person_time_50 <- survivial_df_blood_IgG_n %>%
    filter(!is.na(n_above_50)) %>%
    group_by(n_above_50) %>%
    summarise(
        person_time_days = sum(tstop - tstart, na.rm = TRUE),
        person_time_years = person_time_days / 365.25,
        .groups = "drop"
    )


person_time_50
# Summarise person-time for each n_above_33 category
person_time_33 <- survivial_df_blood_IgG_n %>%
    filter(!is.na(n_above_33)) %>%
    group_by(n_above_33) %>%
    summarise(
        person_time_days = sum(tstop - tstart, na.rm = TRUE),
        person_time_years = person_time_days / 365.25,
        .groups = "drop"
    )



# For 50% threshold
event_summary_50 <- survivial_df_blood_IgG_n %>%
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
event_summary_33 <- survivial_df_blood_IgG_n %>%
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


######## combine the tables: 

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

