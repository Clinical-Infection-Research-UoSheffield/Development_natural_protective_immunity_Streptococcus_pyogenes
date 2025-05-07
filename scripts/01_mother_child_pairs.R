# Title: Analysis of mother child pairs from study of humoral immunity in The Gambia
# Version: 1.0
# Date: 2024-07-18
# Author: Dr Alexander J. Keeley
# Inputs: data/mother_child_pairs_IgG_titres.RDS
# Outputs: Plots to visualise anibody titres in mother child pairs over the first years of life 

# Description:

# This script explores the influence of maternal IgG on the serological profile of Streptococcus pyogenes (Strep A) antigens in the early years of life.
# IgG titres from mother/child pairs were measured at delivery and during the first year of the child's life. The analysis includes:

# 1. Determining the boost in titres between the 6-month visit and subsequent visits for each individual and each antigen.
# 2. Determining whether each individual had serological evidence of exposure 
# 3. Analyzing the titres between paired neonatal cord and maternal serum at delivery.
# 4. Creating summary plots to visualize the serological profiles.

# Requirements:

# Package Versions
# ggdist           ggdist   3.3.2
# ggpubr           ggpubr   0.6.0
# pheatmap       pheatmap  1.0.12
# tidyverse     tidyverse   2.0.0
# wesanderson wesanderson   0.3.7


##############################
# 1. Setup and Data Loading  #
##############################

# Setup environment
source("scripts/setup_environment.R")

# Load required packages using librarian
shelf(tidyverse, ggdist, ggpubr, wesanderson, pheatmap)

# Import titre dataframe from mother child pairs 
titres <- readRDS("data/mother_child_pairs_IgG_titres.RDS")

### Identifying data has been annonymised and individual ids recoded for public availability


# Calculate the number of unique mother/child pairs in the dataset
num_pairs <- titres %>%
    pull(id) %>%
    unique() %>%
    length()

# Print the number of unique pairs
print(paste("Number of unique mother/child pairs in the dataset:", num_pairs))

###############################
# 2. Data Preprocessing for Plotting
###############################


# Determine for each group of id and antigen whether a rise (boost ==1) in titre occured between 6m and next visit

boost <- titres %>%
    filter(participant == "C") %>%     # Filter for newborns only (participant == "C")
    mutate(post_seq = case_when(       # Create sequence order for visits
        participant == "C" & time_point == 0 ~ 0, # Baseline visit for newborns (time_point == 0) is sequence 0
        time_point == 6 ~ 1,                      # 6-month visit is sequence 1
        time_point %in% c(9, 10, 11) ~ 2,         # Visits at 9, 10, or 11 months are sequence 2
        TRUE ~ NA_real_                           # All other cases are set to NA
    )) %>%
    group_by(id, participant, Antigen) %>%
    summarise(
        titre_post_seq_1 = mean(log_RLU_titres[post_seq == 1], na.rm = TRUE),   # Return titre for post_seq 1 (6m) - there is only one value per group
        titre_post_seq_2 = mean(log_RLU_titres[post_seq == 2], na.rm = TRUE)   # Mean titre for post_seq 2 (9, 10, 11m) - there is only one value per group
    ) %>%
    mutate(boost = case_when(
        is.nan(titre_post_seq_1) | is.nan(titre_post_seq_2) ~ 0,   # If any titre is NaN, no boost
        titre_post_seq_2 > titre_post_seq_1 ~ 1,                   # If titre at sequence 2 is greater than sequence 1, boost occurred
        TRUE ~ 0                                                   # Otherwise, no boost
    )) %>%
    ungroup() %>% 
    mutate(change =  titre_post_seq_2 - titre_post_seq_1) %>%    # Calculate the change in titres
    select(id, participant, Antigen, boost, change)

    
### Determine a vaule relative to the birth titre for each group of participant, id and Antigen. 

titres_relative <- titres %>%
    group_by(id, participant, Antigen) %>%
    mutate(base_titre = log_RLU_titres[time_point == 0]) %>%    # return the titre value at time_point == 0 for each group
    mutate(relative_titre =  log_RLU_titres -base_titre ) %>%     # Calculate the relative titre for each row
    filter(participant == "C") %>%  # Filter for newborns only (participant == "C")
    left_join(boost) 

#### Repeat previous two steps but using mutate to create each a value for each row  - which is needed produce the figure 
# Note: boost2 includes per-row values for plotting; boost is a summary table

boost2 <- titres %>%
    filter(participant == "C") %>%
    mutate(post_seq = case_when(
        participant == "C" & time_point == 0 ~ 0,
        time_point == 6 ~ 1,
        time_point %in% c(9, 10, 11) ~ 2,
        TRUE ~ NA_real_
    )) %>%
    group_by(id, participant, Antigen) %>%
    mutate(
        titre_post_seq_1 = mean(log_RLU_titres[post_seq == 1], na.rm = TRUE),
        titre_post_seq_2 = mean(log_RLU_titres[post_seq == 2], na.rm = TRUE)
    ) %>%
    mutate(boost2 = case_when(
        is.nan(titre_post_seq_1) | is.nan(titre_post_seq_2) ~ 0,
        post_seq == 2 & (titre_post_seq_2 > titre_post_seq_1) ~ 1,
        TRUE ~ 0
    )) %>%
    ungroup() %>% 
    mutate(change =  titre_post_seq_2 - titre_post_seq_1)

titres_relative2 <- titres %>%
    group_by(id, participant, Antigen) %>%
    # Compute the titre value at time_point == 0 for each group
    mutate(base_titre = log_RLU_titres[time_point == 0]) %>%
    # Calculate the relative titre for each row
    mutate(relative_titre =  log_RLU_titres -base_titre ) %>%
    filter(participant == "C") %>%
    left_join(boost2)

# Create a variable "serological evidence of exposure" exposure for each ID 
# Summary variable 'exposure' indicates serological evidence of exposure between 6m and follow-up

seroconversion_summary <- boost %>%
    group_by(id) %>%
    mutate(seroconverted_count = sum(boost)) %>%
    mutate(exposure = case_when(
        seroconverted_count >= 2 ~ 1, # returns 1 if there is a rise in titre to two or more Anitgens
        seroconverted_count == 1 & any(boost == 1 & change > 0.5) ~ 1, # returns 1 if there is a rise in titre of more than 0.5 log10 RLU/mL to a single antigen
        TRUE ~ 0 # otherwise returns 0 
    )) %>%
    ungroup() %>%
    select(id, exposure) %>% 
    unique() # only keep the variables id and exposure. 


# Calculate the number and percentage of individuals with serological exposure
total_individuals <- nrow(seroconversion_summary)
exposed_individuals <- sum(seroconversion_summary$exposure)
exposed_percentage <- (exposed_individuals / total_individuals) * 100

# Generate a text summary of those with seroloigcal evidence of exposure
sprintf("Between 6m and the subsequent tested sample (at either 9, 10, or 11 months), %d (%.0f%%) demonstrated serological evidence of exposure defined as a boost in response of IgG titres to two or more antigens, or to one antigen of at least 0.5 log10 RLU/mL.", exposed_individuals, exposed_percentage)

### Ammed dataframe for visualisation of data with different colour schemes and a jittered start point to differentiate mothers from children

plot_df <- titres %>%
    left_join(seroconversion_summary) %>% 
    filter(time_point >= 0) %>%
    mutate(time_point_jittered = case_when(time_point == 0 & participant == "M" ~ time_point - 0.15, 
                                           time_point == 0 & participant == "C" ~ time_point + 0.15,
                                           T ~ time_point)) %>%
    
    mutate(tboost = case_when
           (
               exposure == 1 & visit %in% c("9M","10M","11M") ~ 1,
               T ~ 0
           ))

# create dataframe of just new borns required for plotting individual tranjectories
plot_df2 <- plot_df %>% filter(participant == "C")


#############################
# 3. Plot Longitudinal Data #
#############################

## longitudinal titre figures:
## Points that are colored by green are from individuals with serological evidence of exposure.

plot_01_main01_v1.0 <- 
    plot_df %>%
    ggplot(
        aes(
            x = time_point_jittered,y = log_RLU_titres)) +
    geom_line(alpha = 0.4, col = "grey", aes( x = time_point,y = log_RLU_titres, group = id), data = plot_df2) +
    geom_point(aes(shape = participant,
                   col = interaction(participant, factor(tboost, levels = c(1,0), labels = c("Yes", "No"))))) + 
    geom_smooth(data = titres_relative2 %>% left_join(seroconversion_summary) %>% filter(exposure == 0),
                aes(
                    x = time_point,y = log_RLU_titres, 
                    #group = factor(boost, levels = c(0,1))
                ), method = "loess",se = T, colour = "2b2d42") +
    facet_wrap(~Antigen, ncol = 5) +
    labs(
        x = "Time from delivery (months)",
        y = "IgG level (RLU/mL)",
        col = "Boosted"
    ) +
    scale_colour_manual(values =  wesanderson::wes_palette("FrenchDispatch", n = 4)) +
    theme_minimal() +
    scale_y_continuous(breaks = c(1:6), labels = log10_to_exp) +
    guides(col = "none", shape = "none") +
    theme_universal(base_size = plot_basesize)

plot_01_main01_v1.0


#############################
# 4. Comparison of IgG levels at birth #
#############################

###############################################
#      Extra comments from reviewers #
# calculalate Fetal:materal transfer ratio.  #
###############################################

# calculate F:M ratios on non-transofmred data 

mftr <- titres %>%
    filter(visit == "DEL") %>%
    group_by(Antigen, id) %>%
    select(id, participant, RLU_antibody_titre) %>%
    spread(key = participant, value = RLU_antibody_titre) %>%
    left_join(seroconversion_summary) %>%
    mutate(mftr = C/M,
           exposure = as.character(exposure)) %>%
    ungroup()

# summarise the M:F ratios 
mftr_summary <- mftr %>%
    group_by(Antigen) %>%
    summarise(
        FMR_mean = mean(mftr, na.rm = TRUE),
        FMR_sd = sd(mftr, na.rm = TRUE),
        n = sum(!is.na(mftr)),
        FMR_se = FMR_sd / sqrt(n),
        CI_lower = FMR_mean - 1.96 * FMR_se,
        CI_upper = FMR_mean + 1.96 * FMR_se
    ) %>%
    select(Antigen, FMR_mean, CI_lower, CI_upper)

mftr_summary


#Pairwise comparison between maternal an new born with Paired Wilcoxon test and FDR correction ######

# Filter the dataframe to allow p-value calculation
df2 <- titres %>%
    filter(visit == "DEL") %>%  # To delivery
    mutate(titre = log_RLU_titres) %>%
    select(Antigen, titre, participant, id)  # Include 'id' for pairing

##### calculate the geometric mean for each 

df2 %>%
    group_by(Antigen, participant) %>%
    summarise(
        GMT = log10(psych::geometric.mean(10^titre)))

# Create a function to calculate paired p-values and correct with FDR method
calculate_pvalues_paired <- function(df) {
    pvals <- df %>%
        group_by(Antigen) %>%
        summarise(
            p.value = wilcox.test(
                titre[participant == "M"], 
                titre[participant == "C"], 
                paired = TRUE
            )$p.value
        )
    
    pvals <- pvals %>%
        mutate(p.adjusted = p.adjust(p.value, method = "fdr"))
    
    return(pvals)
}

# Extract adjusted p-values
pvals <- calculate_pvalues_paired(df2)

# Merge p-values back into the original dataframe
df2 <- df2 %>%
    left_join(pvals, by = "Antigen")

# create a function to make annotations of p values for plotting

create_annotations <- function(pvals) {
    annotations <- pvals %>%
        mutate(
            label = ifelse(
                p.adjusted < 0.0001,
                "p.adj < 0.0001",  # Display "< 0.0001" for very small values
                sprintf("p.adj = %.2f", p.adjusted)  # Otherwise, format to 3 decimal places
            )
        ) %>%
        select(Antigen, label)
    
    return(annotations)
}

annotations <- create_annotations(pvals)

fmr_annotations <- mftr_summary %>%
    mutate(
        fmr_label = sprintf("F:M ratio = %.2f\n(%.2f–%.2f)", FMR_mean, CI_lower, CI_upper)
    ) %>%
    select(Antigen, fmr_label)

# Set the colour scheme for the plot

palette <- wes_palette("FrenchDispatch")
custom_palette <- palette[c(2, 3)]

plot_01_main02_v1.0 <- titres %>%
    filter(visit == "DEL") %>%
    ggplot(
        aes(
            x = factor(participant, levels = c("M","C"), labels = c("Mother", "Child")), y = log_RLU_titres, fill = participant, col = participant)) +

    facet_wrap(~Antigen, ncol = 5) +
    geom_line(aes(group = id), col = "grey", alpha = 0.2) +
    geom_boxplot(alpha = 0.1) +
    geom_point(alpha = 0.3) +
    scale_fill_manual(values =  custom_palette) +
    scale_colour_manual(values =  custom_palette) +
    labs(
        x = "Participant",
        y = "IgG level (RLU/mL)"
    ) +
    guides(col = "none", fill = "none") +
    theme_minimal()+ 
    geom_text(data = annotations, aes(x = 1.5, y = 6, label = label), 
              color = "black", size = 3, inherit.aes = FALSE) +
    geom_text(data = fmr_annotations, aes(x = 1.5, y = 3.2, label = fmr_label), 
              color = "black", size = 3, inherit.aes = FALSE) +
    scale_y_continuous(limits = c(3,6), breaks = c(3,4,5,6),labels = log10_to_exp) +
    theme_universal(base_size = plot_basesize)

plot_01_main02_v1.0

# Confirm the numbers in this analysis 
titres %>%
    filter(visit == "DEL") %>%
    pull(id) %>%
    unique %>%
    length()


############################
# Are there any differences between maternal fetal transfer ratios between those who did and did not have serological evidence of exposure 
# Pairwise comparison between M:FR between those with and without serological evidence of exposure. 

# Create a function to calculate paired p-values
calculate_pvalues <- function(df) {
    pvals <- df %>%
        group_by(Antigen) %>%
        summarise(
            p.value = wilcox.test(
                mftr[exposure == "1"], 
                mftr[exposure == "0"], 
                paired = F
            )$p.value
        )
    
    pvals <- pvals %>%
        mutate(p.adjusted = p.adjust(p.value, method = "fdr"))
    
    return(pvals)
}

# Extract adjusted p-values
pvals <- calculate_pvalues(mftr)

# Merge p-values back into the original dataframe
mftr<- mftr %>%
    left_join(pvals, by = "Antigen")





create_annotations <- function(pvals) {
    annotations <- pvals %>%
        mutate(
            label = ifelse(
                p.adjusted < 0.0001,
                "p.adj < 0.0001",  # Display "< 0.0001" for very small values
                sprintf("p.adj = %.2f", p.adjusted)  # Otherwise, format to 3 decimal places
            )
        ) %>%
        select(Antigen, label)
    
    return(annotations)
}

annotations <- create_annotations(pvals)


# Make colour scheme consistent 
palette <- wes_palette("FrenchDispatch")
custom_palette <- palette[c(2, 1)]

mfttr_plot <- mftr %>% 
    ggplot(
        aes(x = as.factor(exposure), y = mftr, colour = as.factor(exposure))) +
    geom_jitter() +
    facet_wrap(~Antigen) +
    geom_text(data = annotations, aes(x = "0", y = 3.5, label = label), 
              color = "black", size = 3, inherit.aes = FALSE) +
    theme_minimal() +
    #scale_colour_manual(values =  wesanderson::wes_palette("FrenchDispatch", n = 2)) +
    scale_colour_manual(values =  custom_palette) +
    theme_universal() +
    guides(colour = "none") +
    labs(y = "Fetal:Maternal IgG transfer ratio",
         x = "Serological evidence of exposure between 6m and subsequent visit")


mfttr_plot 
ggsave("R_output/supp_fetal_maternal_ratios_vs_exposure_V1.0.png", mfttr_plot , dpi = 600, width = 600 / 96, height = 400/96, bg = "white")




 #################################
##### 5. heatmap figures ###########
#################################


# Prepare dataframe for heatmaps figures 

# Extract only the columns 'id', 'Antigen', and 'change' from the boost dataset,
# and remove any rows where the 'change' value is NaN.

boost3 <- boost %>% select(id, Antigen, change) %>%
    filter(!is.nan(change))

# Clustering on 'id':

# Convert the data to wide format where each row represents an observation and
# each column corresponds to an antigen. The cell values are the 'change' values.
wide_boost3_id <- boost3 %>%
    pivot_wider(
        names_from = Antigen,           # Each antigen becomes a column.
        values_from = change,           # Fill these columns with the 'change' values.
        values_fill = list(change = 0)  # Fill missing values with 0.
    ) %>%
    select(-id)  # Remove the 'id' column as it is not needed for clustering.

# Compute Euclidean distances between rows (different ids) and perform hierarchical clustering.
d_id <- dist(as.matrix(wide_boost3_id))
fit_id <- hclust(d_id)
order_id <- fit_id$order  # Extract the clustering order for ids.


# Clustering on 'Antigen':

# Pivot the data so that each row represents an antigen and each column corresponds to an id.
wide_boost3_antigen <- boost3 %>%
    pivot_wider(
        names_from = id,              # Each id becomes a column.
        values_from = change,         # Fill these columns with the 'change' values.
        values_fill = list(change = 0)  # Fill missing values with 0.
    ) %>%
    select(-Antigen)  # Remove the 'Antigen' column; row names will serve as antigen identifiers.

# Compute Euclidean distances between antigens and perform hierarchical clustering.
d_antigen <- dist(as.matrix(wide_boost3_antigen))
fit_antigen <- hclust(d_antigen)
order_antigen <- fit_antigen$order  # Extract the clustering order for antigens.


# Create Ordered Factors for 'id' and 'Antigen':

# Map 'id' to an ordered factor based on the clustering order.
id_order_df <- data.frame(
    id = unique(boost3$id),
    order = order_id
)

# Create the ordered factor for 'id'
boost3$ordered_id <- factor(boost3$id, levels = id_order_df$id[order_id])

# Create a data frame to map 'Antigen' to 'order_antigen'
antigen_order_df <- data.frame(Antigen = unique(boost3$Antigen), order = order_antigen)

# Create the ordered factor for 'Antigen'
boost3$ordered_antigen <- factor(boost3$Antigen, levels = antigen_order_df$Antigen[order_antigen])

# Make it look nice: Define a Continuous Colour Palette for the Heatmap:
# Use a continuous palette from the Wes Anderson package.
pal <- wesanderson::wes_palette("Zissou1", 100, type = "continuous")

# Prepare Matrix for Heatmap (Continuous 'change' Values):
# Convert boost3 to wide format: rows are ordered antigens and columns are ordered ids.
# rows = antigens, columns = individuals, values = change in titre

mat <- boost3 %>%
    select(-id, -Antigen) %>%  # Remove original id and antigen columns.
    pivot_wider(
        names_from = ordered_id,   # Use the ordered id factor as column names.
        values_from = change,
        values_fill = list(change = NA)  # Use NA for missing values.
    ) %>%
    column_to_rownames("ordered_antigen") %>%  # Set the row names to be the ordered antigens.
    as.matrix()

# Set font size for plots
base_size = 10

# Create Heatmap with Dendrograms (Continuous 'change' Values):
plot_01_supp02_V1.0 <- pheatmap(
    mat,
    cluster_rows = TRUE,                 # Cluster the rows (antigens).
    cluster_cols = TRUE,                 # Cluster the columns (ids).
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    display_numbers = FALSE,             # Do not display numbers in the cells.
    show_colnames = FALSE,               # Do not show column names.
    color = pal,                         # Use the continuous colour palette.
    fontsize = base_size,                # Universal font size.
    fontsize_row = base_size * 0.8,        # Slightly smaller font for row labels.
    fontsize_col = base_size * 0.8,        # Slightly smaller font for column labels.
    legend = TRUE,                       # Display a legend.
    cellheight = 20                      # Set the height of each cell.
)
# Display the heatmap with dendrograms.
plot_01_supp02_V1.0


# Create Binary Heatmap of 'change' (Rise vs. No Rise) but keeping the same order as the heatmap above

# Extract the order of rows and columns from the dendrogram
row_order <- plot_01_supp02_V1.0$tree_row$order
col_order <- plot_01_supp02_V1.0$tree_col$order

# Convert matrix to binary for visualising presence/absence of response

binary_mat <- ifelse(mat > 0, 1, 0)

# Reorder binary_mat based on the dendrogram order
binary_mat <- binary_mat[row_order, col_order]

# Remove the column names
##colnames(binary_mat) <- NULL

# Create a custom color palette: white for 0, blue for 1
binary_colors <- c("white", "blue")

# make next panel for dendrogram plot 

plot_01_supp03_V1.0 <- pheatmap(
    binary_mat, 
    cluster_rows = FALSE, 
    cluster_cols = FALSE, 
    display_numbers = FALSE, 
    show_colnames = FALSE,
    color = binary_colors,
    legend_breaks = c(0, 1),
    legend_labels = c("No rise", "Rise"),
    fontsize = base_size,
    fontsize_row = base_size * 0.8,
    fontsize_col = base_size * 0.8,
    cellheight = 20
)


plot_01_supp03_V1.0

## Create Heatmap for serological evidence of exposure ## 

# Ensure 'seroconversion_summary' has the same 'id' order as in 'boost3'
seroconversion_summary_ordered <- seroconversion_summary %>%
    mutate(ordered_id = factor(id, levels = levels(boost3$ordered_id))) %>%
    arrange(ordered_id) %>%
    filter(!is.na(ordered_id))

# Convert to matrix for heatmap, with 'exposure' as the values
seroconversion_mat <-
    
    seroconversion_summary_ordered %>%
    select(ordered_id, exposure) %>%
    pivot_wider(names_from = ordered_id, values_from = exposure, values_fill = list(exposure = 0)) %>%
    as.matrix()

# Create a custom color palette for exposure: white for 0, red for 1
exposure_colors <- c("white", "red")

# Create the heatmap for exposure, maintaining the same order for 'id' without legend and id names

plot_01_supp03B_V1.0 <- pheatmap(seroconversion_mat,
                    cluster_rows = FALSE, 
                    cluster_cols = FALSE, 
                    display_numbers = FALSE, 
                    show_colnames = FALSE,
                    color = exposure_colors,
                    legend_breaks = c(0, 1),
                    legend_labels = c("No", "Yes"),
                    fontsize = base_size,
         cellheight = 10,
                    fontsize_row = base_size * 0.8,
                    fontsize_col = base_size * 0.8)


# Return the plot object
plot_01_supp03B_V1.0

##### Combined output used to make figure 2 panel C. 


### Fin ### 

# keep only plot objects for later compliation

    # List all objects in the environment
all_objects <- ls()
    
    # Identify objects that contain "plot" in their name
plot_objects <- grep("plot_", all_objects, value = TRUE, ignore.case = TRUE)
    
    # Remove all objects that do not contain "plot" in their name
rm(list = setdiff(all_objects, plot_objects), envir = .GlobalEnv)
    
rm(all_objects,plot_objects)





