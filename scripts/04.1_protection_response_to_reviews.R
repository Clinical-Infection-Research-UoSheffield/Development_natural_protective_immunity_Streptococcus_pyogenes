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
#   - data/blood_IgG_titres.RDS
#   - data/incidence_start_dates.RDS"
#   - data/SpyCATS_incidence_df.RDS
#
# Outputs:
#   - Model summary tables with ORs and 95% CIs
#   - VIF diagnostics
#   - AIC value for model comparison
#   - Forest plot of selected model to explore synergy.

# Comment: 
# The statistical model in Figure 4 elegantly explores non-linear relationships in antibody mediated protection. 
# It would be helpful for the authors to clarify the independent contribution of SLO, SpyAD and SpyCEP-responses to the protective signature. 
# This is important as some of these responses are correlated in the overall response. 
# Were interaction effects explored? Do particular combinations provide more predictive accuracy? 
# Do functions predict better than titers? 
# All of these data could be very important for vaccine development efforts.





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



### LOAD DATA

# Load reference dates for incidence window alignment

start_dates <- readRDS("data/incidence_start_dates.RDS")


#####  centered titres above threshold with covariats

# Define path to IgG titre dataset and some key parameters
path_to_titre.df <- "data/blood_IgG_titres.RDS"
sample = "Blood"
class = "IgG"
var_name = "titre"
next_event_window = 45

# Set breakpoints for piecewise regression
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
           IgG = ifelse(titre > transition_point, titre - transition_point, 0),
           hid = substring(pid, 1, 3))


# Center each antigen's titre and IgG values to reduce collinearity
fun_titres <- 
    fun_titres %>%
    group_by(Antigen) %>%
    mutate(c_titre = titre - mean(titre,na.rm = T),
           c_IgG = IgG - mean(IgG, na.rm = T))

# Reshape data: one row per visit per participant, with wide-format IgG and titre columns
antigen_df <- fun_titres %>%
    select(pid, hid, visit_date, Antigen, age_grp, sex, hhsize, titre = c_titre, IgG = c_IgG) %>%
    pivot_wider(
        names_from = Antigen,
        values_from = c(titre, IgG),
        names_glue = "{.value}_{Antigen}"
    )



# Load culture confirmed event incidence data 
pos_incidence_zero <- readRDS("data/SpyCATS_incidence_df.RDS")

# Create a binary indicator of whether another GAS event occurs within 45 days
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
                  rename(date = visit_date)) %>%
    group_by(pid) %>%
    fill(age_grp:IgG_DNAseB) %>%
    ungroup() %>%
    mutate(hid = substring(pid, 1, 3))

# Set reference category for age group to "Over 18 years"
final_df_3$age_grp <- relevel(final_df_3$age_grp, ref = "Over 18 years")

# Fit the logistic regression model including interaction terms 
model_3 <- lme4::glmer(event_next_n ~ IgG_SpyCEP+IgG_SpyAD+IgG_SLO+IgG_SpyCEP:IgG_SpyAD+IgG_SpyCEP:IgG_SLO+IgG_SpyAD:IgG_SLO + age_grp + sex + hhsize + (1 | pid) + (1 | hid), data = final_df_3, family = binomial)

# Extract AIC
model_aic <- AIC(model_3)
print(paste("Model: event_next_n ~ IgG_SpyCEP+IgG_SpyAD+IgG_SLO+IgG_SpyCEP:IgG_SpyAD+IgG_SpyCEP:IgG_SLO+IgG_SpyAD:IgG_SLO + age_grp + sex + hhsize + (1 | pid) + (1 | hid). 
                AIC:"
            , AIC(model_3)))
# Format regression table for presentation
tb3 <- model_3 %>%
    tbl_regression(exponentiate = TRUE) %>%  # Show odds ratios instead of log-odds
    bold_p(t = 0.05) %>%  # Bold significant results
    modify_header(list(label ~ paste(sample, class, "titre")))  # Custom table header

tb3 <- model_3 %>%
    tbl_regression(exponentiate = TRUE) %>%
    bold_p(t = 0.05) %>%
    modify_header(list(label ~ paste(sample, class, "titre"))) 

tb3

# Check for multicollinearity among predictors using Variance Inflation Factor (VIF)
car::vif(model_3)


shelf(broom.mixed)

# Create a tidy data frame from model_3, including exponentiated coefficients (i.e. odds ratios)

model_df <- tidy(model_3, effects = "fixed", conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term != "(Intercept)") %>%
    mutate(
        conf.high.trunc = pmin(conf.high, 5),
        above_limit = conf.high > 5
    )

# Basic forest plot for model_3 showing truncated confidence intervals

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



# Extract coefficients and confidence intervals manually (Wald method)

model_summary <- summary(model_3)
conf_intervals <- confint(model_3, parm = "beta_", method = "Wald")  


# Build a data frame of results including relabeling terms for interpretability and presenation 

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
        conf.high.trunc = pmin(conf.high, 100),
        above_limit = conf.high > 100,
        conf.low.trunc = pmax(conf.low, 0.01),
        below_limit =  conf.low < 0.01
    )



# Reference rows for categories not explicitly modeled

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



# Prepare factor levels for ordered display
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



# Final annotated and faceted forest plot
or_df %>%
    filter(term != "(Intercept)") %>%
    ggplot(aes(x = estimate, y = term)) +  # term already defined with levels
    geom_pointrange(aes(xmin = conf.low.trunc, xmax = conf.high.trunc)) +
    geom_point(
        data = new_rows,
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
        aes(x = 100, label = "→"),  # You can change to "*" or any symbol
        vjust = -0.5, size = 4
    ) +
    geom_text(
        data = filter(or_df, below_limit),
        aes(x = 0.01, label = "←"),  # You can change to "*" or any symbol
        vjust = -0.5, size = 4
    ) +
    scale_x_log10(
        limits = c(0.01, 100),
        breaks = c(0.01, 0.1, 1, 10, 100),
        labels = c("0.01", "0.1", "1", "10", "100")
    ) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    theme_minimal() +
    labs(
        title = paste0("Mixed effect logistic regression including all tires and their interactions. AIC:",round(AIC(model_3),2)),
        x = "Odds Ratio (95% CI)",
        y = ""
    )+ theme(strip.text.y = element_text(angle = 0)) 





##########################################################################
##########################################################################


shelf(lme4)


# Define the shared covariates used across all models
base_covariates <- "+ age_grp + sex + hhsize + (1 | pid) + (1 | hid)"

# Antigen list
antigens <- c("IgG_SpyCEP", "IgG_SpyAD", "IgG_SLO")

# Make a helper function to standardize model info extraction

get_model_info <- function(name, predictors, model) {
    tibble(
        Model = name,
        Predictors = predictors,
        AIC = formatC(AIC(model), format = "f", digits = 2)
    )
}
# -----------------------------
#  Fit models with various combinations of antigens
# -----------------------------

model_3 <- glmer(event_next_n ~ IgG_SpyCEP + IgG_SpyAD + IgG_SLO +
                     IgG_SpyCEP:IgG_SpyAD + IgG_SpyCEP:IgG_SLO + IgG_SpyAD:IgG_SLO +
                     age_grp + sex + hhsize + (1 | pid) + (1 | hid),
                 data = final_df_3, family = binomial)

model_4 <- glmer(event_next_n ~ IgG_SpyCEP + IgG_SpyAD + IgG_SLO +
                     age_grp + sex + hhsize + (1 | pid) + (1 | hid),
                 data = final_df_3, family = binomial)

model_5 <- glmer(event_next_n ~ IgG_SpyCEP +
                     age_grp + sex + hhsize + (1 | pid) + (1 | hid),
                 data = final_df_3, family = binomial)

model_6 <- glmer(event_next_n ~ IgG_SLO +
                     age_grp + sex + hhsize + (1 | pid) + (1 | hid),
                 data = final_df_3, family = binomial)

model_7 <- glmer(event_next_n ~ IgG_SpyAD +
                     age_grp + sex + hhsize + (1 | pid) + (1 | hid),
                 data = final_df_3, family = binomial)

# summarise models 
original_models <- list(
    get_model_info("Model 3", "SpyCEP + SpyAD + SLO + Interactions + covariates", model_3),
    get_model_info("Model 4", "SpyCEP + SpyAD + SLO + covariates", model_4),
    get_model_info("Model 5", "SpyCEP + covariates", model_5),
    get_model_info("Model 6", "SLO + covariates", model_6),
    get_model_info("Model 7", "SpyAD + covariates", model_7)
)

original_models

# -----------------------------
# Pairwise models (no interactions)
# -----------------------------

pairwise_main_models <- combn(antigens, 2, simplify = FALSE) %>%
    map(~ {
        formula_str <- paste0("event_next_n ~ ", paste(.x, collapse = " + "), base_covariates)
        model <- glmer(as.formula(formula_str), data = final_df_3, family = binomial)
        predictors_clean <- gsub("IgG_", "", paste(.x, collapse = " + "))
        get_model_info(paste0("Pairwise: ", predictors_clean), predictors_clean, model)
    })

# -----------------------------
# Models with interaction terms
# -----------------------------

interaction_models <- combn(antigens, 2, simplify = FALSE) %>%
    map(~ {
        interaction_term <- paste(.x, collapse = ":")
        formula_str <- paste0("event_next_n ~ ", paste(.x, collapse = " + "), " + ", interaction_term, base_covariates)
        model <- glmer(as.formula(formula_str), data = final_df_3, family = binomial)
        predictors_clean <- gsub("IgG_", "", paste(.x, collapse = " + "))
        get_model_info(paste0("Interaction: ", gsub("IgG_", "", gsub(":", " * ", interaction_term))),
                       paste0(predictors_clean, " + interaction"), model)
    })

# -----------------------------
# Combine all models into one summary table
# -----------------------------

all_models_table <- bind_rows(
    pairwise_main_models,
    interaction_models
) %>%
    mutate(
        Predictors = paste0(Predictors, " + covariates")
    ) %>%
    bind_rows(original_models) %>%
    arrange(AIC)


# View result
print(all_models_table)

all_models_table %>% write.csv("R_output/aic_synergistic_mdoels_reponse_reviewer.csv")

# -----------------------------
# Visualize AIC comparison
# -----------------------------

base_size = 14

shelf(ggtext)

plot_1 <- all_models_table %>%
    mutate(
        AIC = as.numeric(AIC),
        Predictors = if_else(
            Predictors == "SpyCEP + SpyAD + SLO + covariates",
            "**SpyCEP + SpyAD + SLO + covariates**",
            Predictors
        ),
        Predictors = forcats::fct_reorder(Predictors, AIC)
    ) %>%
    ggplot(aes(x = AIC, y = Predictors)) +
    geom_segment(aes(x = AIC - 2, xend = AIC + 2, yend = Predictors), 
                 linetype = "dotted", color = "grey50") +
    geom_point(size = 2.5, color = "steelblue") +
    theme_minimal() +
    labs(
        x = "AIC",
        y = "Model Predictors"
    ) +

    theme_universal() +
    theme(
        axis.text.y = element_markdown()
    ) 

plot_1

######################## plotting this: 


model_df <- tidy(model_4, effects = "fixed", conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term != "(Intercept)") %>%
    mutate(
        conf.high.trunc = pmin(conf.high, 100),
        above_limit = conf.high > 100
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
model_summary <- summary(model_4)

# -----------------------------
# Format and relabel Model 4 (event_next_n ~ IgG_SpyCEP + IgG_SpyAD + IgG_SLO) output for faceted forest plot
# -----------------------------

# Calculate confidence intervals
conf_intervals <- confint(model_4, parm = "beta_", method = "Wald")  # Wald CIs are common for GLMMs


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
            grepl("age_grp< 2 years", term) ~ "< 2 years",
            grepl("age_grp12-18 years", term) ~ "12-18 years",
            grepl("age_grp2-4 years", term) ~ "2-4 years",
            grepl("age_grp5-11 years", term) ~ "5-11 years",
            grepl("sexFemale", term) ~ "Female",
            grepl("hhsize", term) ~ "Household size",
            TRUE ~ term  # Retain term if it doesn't match any pattern
        )
    ) %>%            mutate(
        conf.high.trunc = pmin(conf.high, 100),
        above_limit = conf.high > 100,
        conf.low.trunc = pmax(conf.low, 0.01),
        below_limit =  conf.low < 0.01
    ) %>%
    filter(term != "(Intercept)")




# Reference labels for plot
    
new_rows <- data.frame(
    term = c(
        "Over 18 years (ref)", "Male (ref)")
) %>%
    mutate(
        estimate = NA_real_,
        conf.low = NA_real_,
        conf.high = NA_real_,
        p.value = NA_real_,
        variable_label = case_when(
            term == "Over 18 years (ref)" ~ "Age group",
            term == "Male (ref)" ~ "Sex",
            TRUE ~ NA_character_
        )
    )



# Ensure 'term' and 'variable_label' are ordered as intended
or_df <- or_df %>%
    mutate(
        term = factor(term, levels = c(
            "SpyCEP", "SpyAD", "SLO",
            "Titre", "Over 18 years (ref)", 
            "12-18 years", "5-11 years", "2-4 years", "< 2 years", 
            "Male (ref)", "Female", "Household size"
        )),
        variable_label = forcats::fct_relevel(variable_label, 
                                              "Titre", "Age group", "Sex", "Household size")
    )


# Build final forest plot
plot_2 <- or_df %>%
    
    ggplot(aes(x = estimate, y = term)) +  # term already defined with levels
    geom_pointrange(aes(xmin = conf.low.trunc, xmax = conf.high.trunc)) +
    geom_point(
        data = new_rows,
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
        aes(x = 100, label = "→"),  # You can change to "*" or any symbol
        vjust = -0.5, size = 4
    ) +
    geom_text(
        data = filter(or_df, below_limit),
        aes(x = 0.01, label = "←"),  # You can change to "*" or any symbol
        vjust = -0.5, size = 4
    ) +
    scale_x_log10(
        limits = c(0.01, 100),
        breaks = c(0.01, 0.1, 1, 10, 100),
        labels = c("0.01", "0.1", "1", "10", "100")
    ) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    theme_minimal() +
    labs(
    #    title = paste0("Mixed effect logistic regression including all protective tires. AIC:",round(AIC(model_4),2)),
        x = "Odds Ratio (95% CI)",
        y = ""
    ) +  
    theme_universal() + 
    theme(strip.text.y = element_text(angle = 0)) 

plot_2

# -----------------------------
# Combine both plots into one figure
# -----------------------------

shelf(cowplot)
synergistic <- plot_grid(
   plot_1, plot_2,
    ncol = 1,
   nrow = 2,
   rel_heights = c(0.4,0.6),
   labels = c("A", "B")
   # Arrange rows vertically
                 # Equal heights for rows
)

synergistic

ggsave("R_output/response_reviewers_combined_effects_protection.png",plot = synergistic, width = 1200 / 96, height = 800/ 96, dpi = 300, bg = "white")


# Save as pdf for final submission 
ggsave("R_output/final_submission/figE_synergistic_protection.pdf", plot = synergistic, width = 1200 / 96, height = 800/ 96, dpi = 300, bg = "white", device = "pdf")









