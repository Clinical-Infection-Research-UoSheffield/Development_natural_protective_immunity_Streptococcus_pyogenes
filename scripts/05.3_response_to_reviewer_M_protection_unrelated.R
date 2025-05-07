
#### Script to respond to reviewers 

# M immunity 


source("scripts/setup_environment.R")
source("scripts/load_functions.R")
base_size = 14
# Set output directory

output_dir <- "R_output/"


#
# Universal theme function
theme_universal <- function(base_size = 14, base_family = "") {
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
                size = base_size * 0.8
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


plot_basesize = 14



# load packages:

shelf(officer,flextable)



#### Step 1 load constant data frame: 

load("R_objects/follow_up_dates_incidence.RData")

# Find enrolment dates for each person
start_dates <- follow_up_dates_incidence %>%
    select(pid, entry_1)


#####  centered titres above threshold with covariats

# To centre data (a method of attempting to handle multicollinearity in we could consider the approach to centre the titre data. 

path_to_titre.df <- "R_objects/blood_IgG_titres.RDS"
sample = "Blood"
class = "IgG"
var_name = "titre"
next_event_window = 45

titre_breakpoint_df <- tibble(Antigen = c("SLO","SpyAD", "SpyCEP", "GAC","DNAseB"), transition_point = c(4.3,4.1,4.3,3,3))


# Read titre data and merge with start dates and age, adjust visit dates, and categorize titres
fun_titres <- read.titres(path_to_titre.df, var_name) %>%
    left_join(start_dates) %>%
    left_join(age) %>%
    mutate(visit_date = as.numeric(ymd(as.character(visit_date))) - entry_1) %>%
    group_by(Antigen) %>%
    ungroup() %>%
    left_join(titre_breakpoint_df) %>%
    mutate(titre_below_threshold = ifelse(titre <= transition_point, titre, transition_point),
           IgG = ifelse(titre > transition_point, titre - transition_point, 0))

fun_titres <- 
    fun_titres %>%
    group_by(Antigen) %>%
    mutate(c_titre = titre - mean(titre,na.rm = T),
           c_IgG = IgG - mean(IgG, na.rm = T))


antigen_df <- fun_titres %>%
    select(pid, hid, visit_date, Antigen, age_grp, sex, hhsize, titre = c_titre, IgG = c_IgG) %>%
    pivot_wider(
        names_from = Antigen,
        values_from = c(titre, IgG),
        names_glue = "{.value}_{Antigen}"
    )

M_titres <- readRDS("../../M_immunity/M_peptide_testing/data/edited/blood_IgG_M_complete_Z_scores.RDS") %>%
    filter(!Antigen %in% c("J8", "K4S2", "P17","M82","M113","M25" ,"M58","M103","M49")) 

M_av <- M_titres %>%
    left_join(start_dates) %>%
    left_join(age) %>%
    mutate(visit_date = as.numeric(ymd(as.character(visit_date))) - entry_1) %>%
    group_by(pid, visit_date) %>%
    summarise(M_average = mean(titre, na.rm = T))

M_wide <- 
    M_titres %>%
    left_join(start_dates) %>%
    left_join(age) %>%
    mutate(visit_date = as.numeric(ymd(as.character(visit_date))) - entry_1) %>%
    spread(key = Antigen, value = titre)
    
### option 1 - the mean 
antigen_df <- antigen_df %>%
    left_join(M_av) %>%
    fill(M_average)

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
    fill(age_grp:IgG_DNAseB) %>%
    ungroup() %>%
    mutate(hid = substring(pid, 1, 3))


final_df_3$age_grp <- relevel(final_df_3$age_grp, ref = "Over 18 years")


model_M1 <- lme4::glmer(event_next_n ~ M_average + age_grp + sex + hhsize + (1 | pid) + (1 | hid), data = final_df_3, family = binomial)

print(paste("Model: event_next_n ~ M_average + age_grp + sex + hhsize + (1 | pid) + (1 | hid). 
                AIC:"
            , AIC(model_M1)))

tb_M1 <- model_M1 %>%
    tbl_regression(exponentiate = TRUE) %>%
    bold_p(t = 0.05) %>%
    modify_header(list(label ~ paste(sample, class, "titre"))) 

tb_M1



#### and antigens: 


model_M1SLO <- lme4::glmer(event_next_n ~ M_average + IgG_SLO + age_grp + sex + hhsize + (1 | pid) + (1 | hid), data = final_df_3, family = binomial)

print(paste("Model: event_next_n ~ M_average + age_grp + sex + hhsize + (1 | pid) + (1 | hid). 
                AIC:"
            , AIC(model_M1)))

tb_M1SLO <- model_M1SLO %>%
    tbl_regression(exponentiate = TRUE) %>%
    bold_p(t = 0.05) %>%
    modify_header(list(label ~ paste(sample, class, "titre"))) 

tb_M1SLO


#

model_M1SpyCEP <- lme4::glmer(event_next_n ~ M_average + IgG_SpyCEP + age_grp + sex + hhsize + (1 | pid) + (1 | hid), data = final_df_3, family = binomial)

print(paste("Model: event_next_n ~ M_average + age_grp + sex + hhsize + (1 | pid) + (1 | hid). 
                AIC:"
            , AIC(model_M1)))

tb_M1SpyCEP <- model_M1SpyCEP %>%
    tbl_regression(exponentiate = TRUE) %>%
    bold_p(t = 0.05) %>%
    modify_header(list(label ~ paste(sample, class, "titre"))) 

tb_M1SpyCEP

# 

model_M1SpyAD <- lme4::glmer(event_next_n ~ M_average + IgG_SpyAD + age_grp + sex + hhsize + (1 | pid) + (1 | hid), data = final_df_3, family = binomial)

print(paste("Model: event_next_n ~ M_average + age_grp + sex + hhsize + (1 | pid) + (1 | hid). 
                AIC:"
            , AIC(model_M1)))

tb_M1SpyAD <- model_M1SpyAD %>%
    tbl_regression(exponentiate = TRUE) %>%
    bold_p(t = 0.05) %>%
    modify_header(list(label ~ paste(sample, class, "titre"))) 

tb_M1SpyAD


### 



model_M2 <- lme4::glmer(event_next_n ~ M_average + (1 | pid) + (1 | hid), data = final_df_3, family = binomial)

print(paste("Model: event_next_n ~ M_average + (1 | pid) + (1 | hid). 
                AIC:"
            , AIC(model_M2)))

tb_M2 <- model_M2 %>%
    tbl_regression(exponentiate = TRUE) %>%
    bold_p(t = 0.05) %>%
    modify_header(list(label ~ paste(sample, class, "titre"))) 

tb_M2



model_M2SLO <- lme4::glmer(event_next_n ~ M_average + IgG_SLO + (1 | pid) + (1 | hid), data = final_df_3, family = binomial)

print(paste("Model: event_next_n ~ M_average + age_grp + sex + hhsize + (1 | pid) + (1 | hid). 
                AIC:"
            , AIC(model_M2SLO)))

tb_M2SLO <- model_M2SLO %>%
    tbl_regression(exponentiate = TRUE) %>%
    bold_p(t = 0.05) %>%
    modify_header(list(label ~ paste(sample, class, "titre"))) 

tb_M2SLO


#

model_M2SpyCEP <- lme4::glmer(event_next_n ~ M_average + IgG_SpyCEP + (1 | pid) + (1 | hid), data = final_df_3, family = binomial)

print(paste("Model: event_next_n ~ M_average + age_grp + sex + hhsize + (1 | pid) + (1 | hid). 
                AIC:"
            , AIC(model_M2SpyCEP)))

tb_M2SpyCEP <- model_M2SpyCEP %>%
    tbl_regression(exponentiate = TRUE) %>%
    bold_p(t = 0.05) %>%
    modify_header(list(label ~ paste(sample, class, "titre"))) 

tb_M2SpyCEP

# 

model_M2SpyAD <- lme4::glmer(event_next_n ~ M_average + IgG_SpyAD + (1 | pid) + (1 | hid), data = final_df_3, family = binomial)

print(paste("Model: event_next_n ~ M_average + IgG_SpyAD +  (1 | pid) + (1 | hid). 
                AIC:"
            , AIC(model_M2SpyAD)))

tb_M2SpyAD <- model_M2SpyAD %>%
    tbl_regression(exponentiate = TRUE) %>%
    bold_p(t = 0.05) %>%
    modify_header(list(label ~ paste(sample, class, "titre"))) 

tb_M2SpyAD





### Adding each titre in 

# Read titre data and merge with start dates and age, adjust visit dates, and categorize titres
fun_titres <- read.titres(path_to_titre.df, var_name) %>%
    left_join(start_dates) %>%
    left_join(age) %>%
    mutate(visit_date = as.numeric(ymd(as.character(visit_date))) - entry_1) %>%
    group_by(Antigen) %>%
    ungroup() %>%
    left_join(titre_breakpoint_df) %>%
    mutate(titre_below_threshold = ifelse(titre <= transition_point, titre, transition_point),
           IgG = ifelse(titre > transition_point, titre - transition_point, 0))

fun_titres <- 
    fun_titres %>%
    group_by(Antigen) %>%
    mutate(c_titre = titre - mean(titre,na.rm = T),
           c_IgG = IgG - mean(IgG, na.rm = T))


antigen_df <- fun_titres %>%
    select(pid, hid, visit_date, Antigen, age_grp, sex, hhsize, titre = c_titre, IgG = c_IgG) %>%
    pivot_wider(
        names_from = Antigen,
        values_from = c(titre, IgG),
        names_glue = "{.value}_{Antigen}"
    )

M_titres <- readRDS("../../M_immunity/M_peptide_testing/data/edited/blood_IgG_M_complete_Z_scores.RDS") %>%
    filter(!Antigen %in% c("J8", "K4S2", "P17","M82","M87","M113","M25" ,"M58","M103","M49")) 

M_wide <- 
    M_titres %>%
    left_join(start_dates) %>%
    left_join(age) %>%
    mutate(visit_date = as.numeric(ymd(as.character(visit_date))) - entry_1) %>%
    spread(key = Antigen, value = titre)

### option 1 - the mean 
antigen_df <- antigen_df %>%
    left_join(M_wide) %>%
    fill(M1:M97)

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
    fill(age_grp:IgG_DNAseB) %>%
    ungroup() %>%
    mutate(hid = substring(pid, 1, 3))


final_df_3$age_grp <- relevel(final_df_3$age_grp, ref = "Over 18 years")




colnames(final_df_3)


# Fit the logistic regression model
model_M3 <- lme4::glmer(event_next_n ~ M1 + M18 + M3 + M4 + M44 + M53 + M55 + M6 + M71 + M74 + M75 + M76 + M89 + M97 + age_grp + sex + hhsize + (1 | pid) + (1 | hid), data = final_df_3, family = binomial)
# Extract AIC
model_aic <- AIC(model_M3)
print(paste("Model: event_next_n ~ IgG_SpyCEP+IgG_SpyAD+IgG_SLO+IgG_SpyCEP:IgG_SpyAD+IgG_SpyCEP:IgG_SLO+IgG_SpyAD:IgG_SLO + age_grp + sex + hhsize + (1 | pid) + (1 | hid). 
                AIC:"
            , AIC(model_M3)))

tbM3 <- model_M3 %>%
    tbl_regression(exponentiate = TRUE) %>%
    bold_p(t = 0.05) %>%
    modify_header(list(label ~ paste(sample, class, "titre"))) 

tbM3

car::vif(model_M3)




###### option 2: 
#create a dataframe that includes: M protection df at the time of events, 


# Read titre data and merge with start dates and age, adjust visit dates, and categorize titres
fun_titres <- read.titres(path_to_titre.df, var_name) %>%
    left_join(start_dates) %>%
    left_join(age) %>%
    mutate(visit_date = as.numeric(ymd(as.character(visit_date))) - entry_1) %>%
    group_by(Antigen) %>%
    ungroup() %>%
    left_join(titre_breakpoint_df) %>%
    mutate(titre_below_threshold = ifelse(titre <= transition_point, titre, transition_point),
           IgG = ifelse(titre > transition_point, titre - transition_point, 0))

fun_titres <- 
    fun_titres %>%
    group_by(Antigen) %>%
    mutate(c_titre = titre - mean(titre,na.rm = T),
           c_IgG = IgG - mean(IgG, na.rm = T))


antigen_df <- fun_titres %>%
    select(pid, hid, visit_date, Antigen, age_grp, sex, hhsize, titre = c_titre, IgG = c_IgG) %>%
    pivot_wider(
        names_from = Antigen,
        values_from = c(titre, IgG),
        names_glue = "{.value}_{Antigen}"
    )

M_av <- M_titres %>%
    left_join(start_dates) %>%
    left_join(age) %>%
    mutate(visit_date = as.numeric(ymd(as.character(visit_date))) - entry_1) %>%
    group_by(pid, visit_date) %>%
    summarise(M_average = mean(titre, na.rm = T))

protection <- readRDS("R_objects/M_protection_dataframe.RDS") %>%
    left_join(start_dates) %>%
    left_join(age) %>%
    mutate(visit_date = as.numeric(ymd(as.character(visit_date))) - entry_1) %>%
    select(pid, visit_date, titre) %>%
    unique

# First, filter M_av for rows not in protection
M_av_filtered <- M_av %>%
    anti_join(protection, by = c("pid", "visit_date")) %>%
    rename(titre = M_average)

# Then bind them together
protection_extended <- bind_rows(protection, M_av_filtered) %>%
    group_by(pid, visit_date) %>%
    summarise(titre = mean(titre))


### option 1 - the mean 
antigen_df <- antigen_df %>%
    left_join(protection_extended)

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
  #  fill(age_grp:IgG_DNAseB) %>%
    fill(age_grp:titre) %>%
    ungroup() %>%
    mutate(hid = substring(pid, 1, 3))




final_df_3$age_grp <- relevel(final_df_3$age_grp, ref = "Over 18 years")




########################



model_M1 <- lme4::glmer(event_next_n ~ titre + age_grp + sex + hhsize + (1 | pid) + (1 | hid), data = final_df_3, family = binomial)

print(AIC(model_M1))

tb_M1 <- model_M1 %>%
    tbl_regression(exponentiate = TRUE) %>%
    bold_p(t = 0.05) %>%
    modify_header(list(label ~ paste(sample, class, "titre"))) 

tb_M1



#### and antigens: 


model_M1SLO <- lme4::glmer(event_next_n ~ titre + IgG_SLO + age_grp + sex + hhsize + (1 | pid) + (1 | hid), data = final_df_3, family = binomial)

print(paste("Model: event_next_n ~ titre + IgG_SLO + + age_grp + sex + hhsize + (1 | pid) + (1 | hid). 
                AIC:"
            , AIC(model_M1)))

tb_M1SLO <- model_M1SLO %>%
    tbl_regression(exponentiate = TRUE) %>%
    bold_p(t = 0.05) %>%
    modify_header(list(label ~ paste(sample, class, "titre"))) 

tb_M1SLO


## 

model_M1SpyAD <- lme4::glmer(event_next_n ~ titre + IgG_SpyAD + age_grp + sex + hhsize + (1 | pid) + (1 | hid), data = final_df_3, family = binomial)

print(paste("Model: event_next_n ~ M_average + IgG_SpyAD + age_grp + sex + hhsize + (1 | pid) + (1 | hid). 
                AIC:"
            , AIC(model_M1)))

tb_M1SpyAD <- model_M1SpyAD %>%
    tbl_regression(exponentiate = TRUE) %>%
    bold_p(t = 0.05) %>%
    modify_header(list(label ~ paste(sample, class, "titre"))) 

tb_M1SpyAD

# 

model_M1SpyCEP <- lme4::glmer(event_next_n ~ titre + IgG_SpyCEP + age_grp + sex + hhsize + (1 | pid) + (1 | hid), data = final_df_3, family = binomial)

print(paste("Model: event_next_n ~ titre + IgG_SpyCEP + age_grp + sex + hhsize + (1 | pid) + (1 | hid). 
                AIC:"
            , AIC(model_M1)))

tb_M1SpyCEP <- model_M1SpyCEP %>%
    tbl_regression(exponentiate = TRUE) %>%
    bold_p(t = 0.05) %>%
    modify_header(list(label ~ paste(sample, class, "titre"))) 

tb_M1SpyCEP




### without age: 


model_M2 <- lme4::glmer(event_next_n ~ titre + (1 | pid) + (1 | hid), data = final_df_3, family = binomial)

print(paste("Model: event_next_n ~ titre + (1 | pid) + (1 | hid). 
                AIC:"
            , AIC(model_M2)))

tb_M2 <- model_M2 %>%
    tbl_regression(exponentiate = TRUE) %>%
    bold_p(t = 0.05) %>%
    modify_header(list(label ~ paste(sample, class, "titre"))) 

tb_M2



model_M2SLO <- lme4::glmer(event_next_n ~ titre + IgG_SLO + (1 | pid) + (1 | hid), data = final_df_3, family = binomial)

print(paste("Model: event_next_n ~ titre + IgG_SLO + age_grp + sex + hhsize + (1 | pid) + (1 | hid). 
                AIC:"
            , AIC(model_M2SLO)))

tb_M2SLO <- model_M2SLO %>%
    tbl_regression(exponentiate = TRUE) %>%
    bold_p(t = 0.05) %>%
    modify_header(list(label ~ paste(sample, class, "titre"))) 

tb_M2SLO


#

model_M2SpyCEP <- lme4::glmer(event_next_n ~ titre + IgG_SpyCEP + (1 | pid) + (1 | hid), data = final_df_3, family = binomial)

print(paste("Model: event_next_n ~ titre + age_grp + IgG_SpyCEP sex + hhsize + (1 | pid) + (1 | hid). 
                AIC:"
            , AIC(model_M2SpyCEP)))

tb_M2SpyCEP <- model_M2SpyCEP %>%
    tbl_regression(exponentiate = TRUE) %>%
    bold_p(t = 0.05) %>%
    modify_header(list(label ~ paste(sample, class, "titre"))) 

tb_M2SpyCEP

# 

model_M2SpyAD <- lme4::glmer(event_next_n ~ titre + IgG_SpyAD + (1 | pid) + (1 | hid), data = final_df_3, family = binomial)

print(paste("Model: event_next_n ~ titre + IgG_SpyAD +  (1 | pid) + (1 | hid). 
                AIC:"
            , AIC(model_M2SpyAD)))

tb_M2SpyAD <- model_M2SpyAD %>%
    tbl_regression(exponentiate = TRUE) %>%
    bold_p(t = 0.05) %>%
    modify_header(list(label ~ paste(sample, class, "titre"))) 


tb_M2SpyAD
tbl_merge(tbls = list(tb_M1SLO, tb_M1SpyAD, tb_M1SpyCEP))



#######################

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
       # fill(titre) %>%
        ungroup() %>%
        mutate(hid = substring(pid,1,3)) %>%
        unique()
    
    # then do logistic regression: 
    
    
    stopifnot("pid" %in% names(final_df_3))
    
    
    model <- lme4::glmer(event_next_n ~ titre + titre_above_threshold + (1 | pid) + (1 | hid), data=final_df_3, family=binomial)
    
    print(paste0("model1 (M + ", antigen, " above transitionpoint):", AIC(model)))
    
    model2 <- lme4::glmer(event_next_n ~ titre + titre_above_threshold + age_grp + sex + hhsize+  (1 | pid) + (1 | hid), data=final_df_3, family=binomial)
    
    print(paste0("model2 (M + ", antigen, " above transitionpoint -fully adjusted):", AIC(model2)))
    
    model3 <- lme4::glmer(event_next_n ~ titre + titre_above_threshold + age_grp + (1 | pid) + (1 | hid), data=final_df_3, family=binomial)
    
    print(paste0("model3 (M + ", antigen, " above transitionpoint + age_grp):", AIC(model3)))
    
    model4 <- lme4::glmer(event_next_n ~ titre_above_threshold  + age_grp + sex + hhsize+  (1 | pid) + (1 | hid), data=final_df_3, family=binomial)
    
    print(paste0("model4 (", antigen, " above transitionpoint  -fully adjusted only - No M):", AIC(model4)))
    
    model5 <- lme4::glmer(event_next_n ~ titre_above_threshold  + (1 | pid) + (1 | hid), data=final_df_3, family=binomial)
    
    print(paste0("model5 (", antigen, " above transitionpoint  -fully adjusted only - No M):", AIC(model5)))
    
    model6 <- lme4::glmer(event_next_n ~ titre + (1 | pid) + (1 | hid), data=final_df_3, family=binomial)
    
    print(paste0("model6 (M only):", AIC(model6)))
    
    model7 <- lme4::glmer(event_next_n ~ titre + age_grp + sex + hhsize+  (1 | pid) + (1 | hid), data=final_df_3, family=binomial)
    
    print(paste0("model7 (M only -fully adjusted):", AIC(model7)))
    
    model8 <- lme4::glmer(event_next_n ~ titre * titre_above_threshold + age_grp + sex + hhsize+  (1 | pid) + (1 | hid), data=final_df_3, family=binomial)
    
    print(paste0("model8 (M + ", antigen, " above transitionpoint incl. interaction -fully adjusted):", AIC(model8)))
    
    model9 <- lme4::glmer(event_next_n ~ titre + titre_above_threshold + age_grp + sex +  (1 | pid) + (1 | hid), data=final_df_3, family=binomial)
    
    print(paste0("model9 (M + ", antigen, " above transitionpoint + age + sex):", AIC(model9)))
    
    
    
    print(model8%>%
              tbl_regression(exponentiate = TRUE) %>%
              add_nevent() %>%
              bold_p(t = 0.05)) 
    
    tb1 <- model %>%
        tbl_regression(exponentiate = TRUE) %>%
        add_nevent() %>%
        bold_p(t = 0.05)  # Modify header here 
    
    tb2 <- model2 %>%
        tbl_regression(exponentiate = TRUE) %>%
        add_nevent() %>%
        bold_p(t = 0.05) 
    
    tb3 <- model3 %>%
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


tb1_list <- list()
tb2_list <- list()
tb3_list <- list()
plot_list <- list()
plot_row_list <- list()

results_df <- data.frame(antigen = character(), term = character(), estimate = numeric(),
                         conf.low = numeric(), conf.high = numeric(), p.value = numeric(),
                         stringsAsFactors = FALSE)


M_av <- M_titres %>%
    group_by(pid, visit_date) %>%
    summarise(M_average = mean(titre, na.rm = T))

protection <- readRDS("R_objects/M_protection_dataframe.RDS") %>%
    select(pid, visit_date, titre) %>%
    unique

# First, filter M_av for rows not in protection
M_av_filtered <- M_av %>%
    anti_join(protection, by = c("pid", "visit_date")) %>%
    rename(titre = M_average)

# Then bind them together
protection_extended <- bind_rows(protection, M_av_filtered) %>%
    group_by(pid, visit_date) %>%
    summarise(titre = mean(titre))



Antigen_list = c("SLO","SpyAD","SpyCEP")
for (ag in Antigen_list) {
    results <- create_regression_glmer_conserved_M(
        M_titres = protection_extended,
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

plot_09_fig06_panelG <- M_forrest
plot_09_fig06_panelG




#### 

final_dems %>%
    ggplot(aes(x = age)) +
    geom_histogram(fill = "#4D759A") +
    theme_minimal()
