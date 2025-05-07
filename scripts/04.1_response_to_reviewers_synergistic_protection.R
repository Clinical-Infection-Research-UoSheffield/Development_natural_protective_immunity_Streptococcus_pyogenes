
# Title: Responding to reviewer comments – multivariate analysis of IgG titres
# Author: Alexander J Keeley
# Date: 2025-01-22
# Description:
#   This script was written in response to reviewer comments requesting:
#   - Assessment of the independent effects of SLO, SpyAD, and SpyCEP
#   - Evaluation of interaction effects between responses
#   - Visualization of odds ratios and model confidence
#   - Commentary on potential collinearity
#
# Inputs:
#   - R_objects/blood_IgG_titres.RDS
#   - R_objects/follow_up_dates_incidence.RData
#   - R_objects/SpyCATS_incidence_df.RDS
#
# Outputs:
#   - Forest plots (linear and log scale)
#   - Model summary tables with ORs and 95% CIs
#   - VIF diagnostics
#   - AIC value for model comparison




source("scripts/setup_environment.R")
source("scripts/load_functions.R")

# Set output directory

output_dir <- "R_output/"


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

# Fit the logistic regression model
model_3 <- lme4::glmer(event_next_n ~ IgG_SpyCEP+IgG_SpyAD+IgG_SLO+IgG_SpyCEP:IgG_SpyAD+IgG_SpyCEP:IgG_SLO+IgG_SpyAD:IgG_SLO + age_grp + sex + hhsize + (1 | pid) + (1 | hid), data = final_df_3, family = binomial)

# Extract AIC
model_aic <- AIC(model_3)
print(paste("Model: event_next_n ~ IgG_SpyCEP+IgG_SpyAD+IgG_SLO+IgG_SpyCEP:IgG_SpyAD+IgG_SpyCEP:IgG_SLO+IgG_SpyAD:IgG_SLO + age_grp + sex + hhsize + (1 | pid) + (1 | hid). 
                AIC:"
            , AIC(model_3)))

tb3 <- model_3 %>%
    tbl_regression(exponentiate = TRUE) %>%
    bold_p(t = 0.05) %>%
    modify_header(list(label ~ paste(sample, class, "titre"))) 

tb3

car::vif(model_3)


shelf(broom.mixed)


model_df <- tidy(model_3, effects = "fixed", conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term != "(Intercept)") %>%
    mutate(
        conf.high.trunc = pmin(conf.high, 5),
        above_limit = conf.high > 5
    )

ggplot(model_df, aes(x = estimate, y = term)) +
    geom_pointrange(aes(xmin = conf.low, xmax = conf.high.trunc)) +
    geom_text(
        data = filter(model_df, above_limit),
        aes(x = 5, label = "→"),  # You can change to "*" or any symbol
        vjust = -0.5, size = 4
    ) +
    coord_cartesian(xlim = c(0, 5)) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    theme_minimal() +
    labs(
        title = "Forest Plot for Logistic Regression Model",
        x = "Odds Ratio (95% CI)",
        y = ""
    )



# Extract the summary of the model
model_summary <- summary(model_3)

# Calculate confidence intervals
conf_intervals <- confint(model_3, parm = "beta_", method = "Wald")  # Wald CIs are common for GLMMs


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
            grepl("IgG", term) ~ "Titre",
            grepl("age_grp", term) ~ "Age group",
            grepl("hhsize", term) ~ "Household size",
            grepl("sexFemale", term) ~ "Sex",
            TRUE ~ NA_character_  # Default to NA if no match
        ),
        # Simplify term names for clarity
        term = case_when(
            term == "IgG_SpyCEP" ~ "SpyCEP",
            term == "IgG_SpyAD" ~ "SpyAD",
            term == "IgG_SLO" ~ "SLO",
            term == "IgG_SpyCEP:IgG_SpyAD" ~ "SpyCEP:SpyAD",
            term == "IgG_SpyCEP:IgG_SLO" ~ "SpyCEP:SLO",
            term == "IgG_SpyAD:IgG_SLO" ~ "SpyAD:SLO",
            grepl("age_grp< 2 years", term) ~ "< 2 years",
            grepl("age_grp12-18 years", term) ~ "12-18 years",
            grepl("age_grp2-4 years", term) ~ "2-4 years",
            grepl("age_grp5-11 years", term) ~ "5-11 years",
            grepl("sexFemale", term) ~ "Female",
            grepl("hhsize", term) ~ "Household size",
            TRUE ~ term  # Retain term if it doesn't match any pattern
        )
    ) %>%            mutate(
        conf.high.trunc = pmin(conf.high, 5),
        above_limit = conf.high > 5
    )




new_rows <- data.frame() %>%
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

new_rows <- data.frame(
    term = c("Below threshold (ref)", "Over 18 years (ref)", "Male (ref)")
) %>%
    mutate(
        estimate = NA_real_,
        conf.low = NA_real_,
        conf.high = NA_real_,
        p.value = NA_real_,
        variable_label = case_when(
            term == "Below threshold (ref)" ~ "Relation to titre threshold",
            term == "Over 18 years (ref)" ~ "Age group",
            term == "Male (ref)" ~ "Sex",
            TRUE ~ NA_character_
        )
    )





# Plot

or_df %>%
    mutate(
        conf.high.trunc = pmin(conf.high, 5),
        above_limit = conf.high > 5
    ) %>%
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
                                         "Titre", "Age group", "Sex", "Household size")), 
        scales = "free", 
        space = "free_y"
    ) +
    
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


unique(or_df$term)


# Ensure 'term' and 'variable_label' are ordered as intended
or_df <- or_df %>%
    mutate(
        term = factor(term, levels = c(
            "SpyAD:SLO", "SpyCEP:SLO", "SpyCEP:SpyAD","SpyCEP", "SpyAD", "SLO",
            "Titre", "Over 18 years (ref)", 
            "12-18 years", "5-11 years", "2-4 years", "< 2 years", 
            "Male (ref)", "Female", "Household size"
        )),
        variable_label = forcats::fct_relevel(variable_label, 
                                              "Titre", "Age group", "Sex", "Household size")
    )

or_df %>%
    mutate(
        conf.high.trunc = pmin(conf.high, 5),
        above_limit = conf.high > 5
    ) %>%
    filter(term != "(Intercept)") %>% # Exclude the intercept
    ggplot(aes(x = estimate, y = term)) +  # term already defined with levels
    geom_pointrange(aes(xmin = conf.low, xmax = conf.high.trunc)) +
    geom_point(
        data = new_rows %>% filter(!antigen %in% c("DNAseB", "GAC")),
        aes(x = 1, y = term), color = "blue", size = 3
    ) +
    facet_grid(
        rows = vars(forcats::fct_relevel(variable_label, 
                                         "Titre", "Age group", "Sex", "Household size")), 
        scales = "free", 
        space = "free_y"
    ) +
    
    geom_text(
        data = filter(or_df, above_limit),
        aes(x = 5, label = "→"),  # You can change to "*" or any symbol
        vjust = -0.5, size = 4
    ) +
    coord_cartesian(xlim = c(0, 5)) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    theme_minimal() +
    labs(
        title = "Mixed effect logistic regression including all tires and their interactions",
        x = "Odds Ratio (95% CI)",
        y = ""
    )+ theme(strip.text.y = element_text(angle = 0))



####### now plot on log 10: 

# Prepare or_df with transformed values for log scale plotting
or_df_plot <- or_df %>%
    mutate(
        conf.high.trunc = pmin(conf.high, 100),
        above_limit = conf.high > 100,
        estimate_log = log10(estimate),
        conf.low_log = log10(conf.low),
        conf.high_log = log10(conf.high.trunc)
    )

# Prepare new_rows for log10 x-axis (1 becomes log10(1) = 0)
new_rows_plot <- new_rows %>%
    mutate(
        estimate_log = 0  # log10(1)
    )

# Final plot with log10-transformed x-axis
or_df_plot %>%
    filter(term != "(Intercept)") %>%
    ggplot(aes(x = estimate_log, y = term)) +
    geom_pointrange(aes(xmin = conf.low_log, xmax = conf.high_log)) +
    geom_point(
        data = new_rows_plot,
        aes(x = estimate_log, y = term),
        color = "blue", size = 3
    ) +
    facet_grid(
        rows = vars(forcats::fct_relevel(variable_label, 
                                         "Titre", "Age group", "Sex", "Household size")), 
        scales = "free", 
        space = "free_y"
    ) +
    geom_text(
        data = filter(or_df_plot, above_limit),
        aes(x = log10(100), label = "→"),  # adjust to match new axis
        vjust = -0.5, size = 4
    ) +
    scale_x_continuous(
        trans = "identity",  # Already log-transformed manually
        limits = c(log10(0.01), log10(100)),
        breaks = log10(c(0.01, 0.1, 1, 10, 100)),
        labels = c("0.01", "0.1", "1", "10", "100")
    ) +
    geom_vline(xintercept = 0, linetype = "dashed") +  # log10(1) = 0
    theme_minimal() +
    labs(
        title = "Mixed effect logistic regression including all titres and their interactions",
        x = "log10 Odds Ratio (95% CI)",
        y = ""
    ) +
    theme(strip.text.y = element_text(angle = 0))
