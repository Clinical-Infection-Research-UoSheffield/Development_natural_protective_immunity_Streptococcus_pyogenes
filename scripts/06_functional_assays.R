
# Title: Functional Immunoassay Analysis and Opsonophagocytic Activity in selected sera from the SpyCATS Cohort

# Version: 1.0
# Date: 2025-02-05
# Author: Alexander J Keeley 

# Inputs: 
#  - "scripts/Compile data.R" - to compile functional data frame 
#   - data/baseline_blood_no_disease_titres.RDS
#   - data/bloodIgG_protective_threshold_df.RDS

# Outputs: 
#   - Plots and summary statistics describing functional assay results in relation to binding IgG levels and correlation matrices. 

# Description:
#
# This script performs an analysis of IgG titres and their functional relationships within the SpyCATS cohort. 
# The analysis includes:
#
# 1. Extracting IgG levels from participants and comparing titres to putative protective thresholds.
# 2. Generating statistical summaries of functional responses and their relationship to protective thresholds.
# 3. Performing correlation matrix analysis for functional assay data, including opsonophagocytic activity (OPA) and binding IgG levels.
# 4. Visualizing the relationship between IgG titres and bacterial opsonophagocytosis using statistical comparisons.
# 5. Comparing the proportions of participants exhibiting functional opsonophagocytosis across between those above and below IgG protective thresholds.
# 6. Creating visual summaries of IgG titre distributions, correlation matrices, and functional outcomes through multiple panels.

# Data Availability:
    
# The data used in this analysis is owned by GSK Vaccines and is not publicly available.
# However, requests for access to the data can be made to the corresponding authors and will be considered.

# Requirements:

# Package Version

# Package  Version
# corrplot   corrplot     0.92
# cowplot     cowplot    1.1.3
# ggdist       ggdist    3.3.2
# ggplot2     ggplot2    3.5.1
# ggpubr       ggpubr    0.6.0
# librarian librarian    1.8.1
# magick       magick    2.8.5
# psych         psych 2.4.6.26
# tidyr         tidyr    1.3.1
# tidyverse tidyverse    2.0.0

# Load dependancies: 

library(librarian)

shelf(tidyverse, ggplot2, ggpubr, magick, cowplot)

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


plot_basesize = 16


### source the dataframes - please note this is not publically available 

source("scripts/Compile data.R")



#######################################################
# Panel A - titres in relation to protective thresholds 
# (replicated from "scripts/02_blood_IgG_baseline.R")
#######################################################

df2 <- readRDS("R_objects/baseline_blood_no_disease_titres.RDS")  
protective_threshold <- readRDS("R_objects/bloodIgG_protective_threshold_df.RDS")

protective_threshold %>%
    mutate(RLU = 10^Threshold)

# Calculate percentages for each Antigen

percentage_labels <- df2 %>%
    filter(Antigen %in% c("SLO", "SpyAD", "SpyCEP")) %>%
    left_join(protective_threshold, by = join_by(Antigen)) %>%  # Correct join
    mutate(above_threshold = ifelse(titre > Threshold, 1, 0)) %>%
    group_by(Antigen) %>%
    summarise(
        n = n(),
        n_above = sum(above_threshold),  # Sum the `above_threshold` column to count occurrences
        percentage_above = mean(above_threshold) * 100  # Calculate the mean percentage
    ) %>%
    mutate(label = paste0(round(percentage_above, 0), "%"))  # Create a label for annotation

percentage_labels


# Dynamically create the improved sentence using sprintf
sentence <- sprintf(
    "At baseline, antibody levels from the cohort showed that %d individuals (%s) for SLO, %d (%s) for SpyAD, and %d (%s) for SpyCEP had IgG levels above the protective thresholds (Figure 5A).",
    percentage_labels$n_above[percentage_labels$Antigen == "SLO"], 
    percentage_labels$label[percentage_labels$Antigen == "SLO"], 
    percentage_labels$n_above[percentage_labels$Antigen == "SpyAD"], 
    percentage_labels$label[percentage_labels$Antigen == "SpyAD"], 
    percentage_labels$n_above[percentage_labels$Antigen == "SpyCEP"], 
    percentage_labels$label[percentage_labels$Antigen == "SpyCEP"]
)

# Print the sentence
print(sentence)



# Create the plot with annotations
panelA <- df2 %>%
    filter(Antigen %in% c("SLO", "SpyAD", "SpyCEP")) %>%
    left_join(protective_threshold, by = join_by(Antigen)) %>%  # Correct join
    mutate(above_threshold = ifelse(titre > Threshold, 1, 0)) %>%
    ggplot(aes(x = age, y = titre, col = factor(above_threshold))) +  # Use factor to color by above/below threshold
    guides(color = "none") +
    facet_wrap(~Antigen) +
    geom_point(alpha = 0.8, size = 3) +
    scale_color_manual(values = c("#d73027","#7570b3")) +
    labs(
        y = "IgG level (log10 RLU/mL)",
        x = 'Age'
        # title = "Blood IgG titres above 50% protective threshold by age group"
    ) +
    geom_hline(aes(yintercept = Threshold), linetype = "dashed", color = "red") +  # Reference Threshold correctly
    theme_minimal() +
    theme_universal(base_size = plot_basesize) +
    # Annotate percentages
    geom_text(
        data = percentage_labels,
        aes(x = Inf, y = Inf, label = label),  # Place in top-right corner
        inherit.aes = FALSE, hjust = 1.2, vjust = 1.2, 
        size = plot_basesize
    )

panelA


###### 

############################################################
# Panel 2: correlation matric for all functional IC50 data #
############################################################


# Take the value to make a correlation matrix 
lim_corr_df <-
    corr_df %>%
    select(DNAseB_IgG:M1_THP1_IC50)

# Update column names for plotting 
updated_colnames <- colnames(lim_corr_df) %>%
    # Remove underscores
    gsub("_", " ", .) %>%
    # Replace "THP1 IC50" with "phagocytosis"
    gsub("THP1 IC50", "phagocytosis", .) %>%
    # Replace "GAC phagocytosis" and "SpyAD phagocytosis" with "<GAC/SpyAD> beads: phagocytosis"
    gsub("GAC phagocytosis", "GAC beads: phagocytosis", .) %>%
    gsub("SpyAD phagocytosis", "SpyAD beads: phagocytosis", .) %>%
    # Replace "M1 phagocytosis" with "M1 bacteria phagocytosis"
    gsub("M1 phagocytosis", "M1 bacteria: phagocytosis", .)


# Create a Spearman's correlation matrix 

corr_matrix_interaction <- lim_corr_df %>%
    cor(use = "pairwise.complete.obs", method = "spearman") 

# Visualise this dataframe
lim_corr_df %>% head()

#####f unction to calculate correlation coefficients and p-values
get_cor_and_pval <- function(df) {
    # Initialize an empty dataframe to store results
    results <- tibble(Var1 = character(), Var2 = character(),
                      Correlation = numeric(), PValue = numeric())
    
    # Loop through all combinations of variables
    for (i in 1:(ncol(df) - 1)) {
        for (j in (i + 1):ncol(df)) {
            # Perform correlation test
            test <- cor.test(df[[i]], df[[j]], use = "pairwise.complete.obs",method = "spearman")
            
            # Store results in the dataframe
            results <- results %>%
                add_row(Var1 = colnames(df)[i],
                        Var2 = colnames(df)[j],
                        Correlation = test$estimate,
                        PValue = test$p.value)
        }
    }
    return(results)
}

# Apply the function to your dataframe
cor_pval_results <- get_cor_and_pval(lim_corr_df)

# Get the p values for each element of the matrix
p_values_matrix <- psych::corr.test(lim_corr_df, use = "pairwise.complete.obs", method = "spearman")$p

# Set updated column names to the dataframe and correlation matrix
colnames(lim_corr_df) <- updated_colnames
colnames(corr_matrix_interaction) <- updated_colnames
rownames(corr_matrix_interaction) <- updated_colnames

# Save the corrplot as an image

png("R_output/Fig_5B_110225.png", width = 4000, height = 4000, res = 300)

corrplot::corrplot.mixed(corr_matrix_interaction, 
                         upper = "color", 
                         lower = "number",
                         tl.col = "black", 
                         tl.srt = 45,            # Rotate top labels to 90 degrees
                         lower.col = "black", 
                         number.cex = 2,
                         tl.pos = "lt", 
                         tl.cex = 2, 
                         cl.cex = 1.5,           # Increase legend number size
                         p.mat = p_values_matrix,
                         insig = "p-value") # "lt" for left and top labels
dev.off()


### Combine Panel A and Panel B for the figure 

panelB <-image_read("R_output/Fig_5B_110225.png")

panelB <- ggdraw() + draw_image(panelB)
panelB
row1 <- plot_grid(
    panelA, panelB,
    nrow = 1,
    labels = c("A","B"),
    rel_widths = c(0.55,0.45)
) 

row1

###############################################################
# Panel 3: IC50 plots in relation to the protective threshold #
###############################################################

# For each functional assay the IC50 are plotted as rain plots in relation to putative protective thresholds 
# and compared with Mann Whitney U test 

A1<- 
    SLO.df %>%
    filter(!is.na(above_threshold)) %>%
    ggplot(
        aes(
            x = factor(above_threshold, labels = c("No", "Yes")), y = log10(haemolysisIC50), col = as.factor(above_threshold), fill = as.factor(above_threshold)
        )
    ) +
    labs(
        x = "IgG above 50% protective threshold",
        y = "Inhibition of haemolysis (log10 IC50)", 
        title = "SLO"
    ) +
    theme_minimal() +
    geom_boxplot(width = 0.05, outliers = F) +
    #  geom_violin(trim = FALSE) +  # Add violin plot
    ggdist::geom_dots(

        width = 0.3,
        size = 4,
        position = position_nudge(x = 0.075)
    ) +
    ggdist::stat_slab(
        alpha = 0.5,
        side = 'bottom', 
        position = position_nudge(x = -0.075), 
        width = 0.3
    ) +
 #   ggpubr::stat_compare_means() +
    ggpubr::stat_compare_means(
        aes(label = ifelse(..p.. < 0.0001, "p < 0.0001", sprintf("p = %.4f", ..p..))),
        method = "wilcox.test", # Replace with your method if needed
        size = 5
        ) +

    scale_color_manual(values = c("#d73027","#7570b3")) +
    scale_fill_manual(values = c("#d73027","#7570b3")) +
    #  scale_color_manual(values = c(StrepA_colscheme_light[["SLO"]], StrepA_colscheme_dark[["SLO"]])) +
    # scale_fill_manual(values = c(StrepA_colscheme_light[["SLO"]], StrepA_colscheme_dark[["SLO"]])) +
    guides(col = "none", fill = "none")  +
    theme_universal(base_size = plot_basesize)

A1


B1 <-
    SpyAD.df %>%
    filter(!is.na(above_threshold)) %>%
    #   filter(haemolysis_activity == 1) %>%
    ggplot(
        aes(
            x = factor(above_threshold, labels = c("No", "Yes")), y = log10(SpyAD_THP1_IC50), col = as.factor(above_threshold), fill = as.factor(above_threshold)
        )
    ) +
    labs(
        x = "IgG above 50% protective threshold",
        y = "Phagocytosis of SpyAD beads (log10 IC50)", 
        title = "SpyAD"
    ) +
    theme_minimal() +
    geom_boxplot(width = 0.05, outliers = F) +
    #  geom_violin(trim = FALSE) +  # Add violin plot
    ggdist::geom_dots(
        width = 0.3,
        size = 4,
        position = position_nudge(x = 0.075)
    ) +
    ggdist::stat_slab(
        alpha = 0.5,
        side = 'bottom', 
        position = position_nudge(x = -0.075), 
        width = 0.3
    ) +
    scale_color_manual(values = c("#d73027","#7570b3")) +
    scale_fill_manual(values = c("#d73027","#7570b3")) +
    ggpubr::stat_compare_means(
        aes(label = ifelse(..p.. < 0.0001, "p < 0.0001", sprintf("p = %.4f", ..p..))),
        method = "wilcox.test", # Replace with your method if needed
        size = 5
        ) +
    #  scale_color_manual(values = c(StrepA_colscheme_light[["SpyAD"]], StrepA_colscheme_dark[["SpyAD"]])) +
    #  scale_fill_manual(values = c(StrepA_colscheme_light[["SpyAD"]], StrepA_colscheme_dark[["SpyAD"]])) +
    guides(col = "none", fill = "none")  +
    theme_universal(base_size = plot_basesize)
    

B1



C1 <- SpyCEP.df %>%
    filter(!is.na(above_threshold)) %>%
    #   filter(haemolysis_activity == 1) %>%
    ggplot(
        aes(
            x = factor(above_threshold, labels = c("No", "Yes")), y = log10(IL8_IC50), col = as.factor(above_threshold), fill = as.factor(above_threshold)
        )
    ) +
    labs(
        x = "IgG above 50% protective threshold",
        y = "Inhibition of IL8 cleavage (log10 IC50)", 
        title = "SpyCEP"
    ) +
    theme_minimal() +
    geom_boxplot(width = 0.05,outliers = F) +
    #  geom_violin(trim = FALSE) +  # Add violin plot
    ggdist::geom_dots(
        
        width = 0.3,
        size = 4,
        position = position_nudge(x = 0.075)
    ) +
    ggdist::stat_slab(
        alpha = 0.5,
        side = 'bottom', 
        position = position_nudge(x = -0.075), 
        width = 0.3
    ) +
    scale_color_manual(values = c("#d73027","#7570b3")) +
    scale_fill_manual(values = c("#d73027","#7570b3")) +
    ggpubr::stat_compare_means(
        aes(label = ifelse(..p.. < 0.0001, "p < 0.0001", sprintf("p = %.4f", ..p..))),
        method = "wilcox.test", # Replace with your method if needed
        size = 5
        ) +
    #  scale_color_manual(values = c(StrepA_colscheme_light[["SpyCEP"]], StrepA_colscheme_dark[["SpyCEP"]])) +
    #  scale_fill_manual(values = c(StrepA_colscheme_light[["SpyCEP"]], StrepA_colscheme_dark[["SpyCEP"]])) +
    guides(col = "none", fill = "none")  +
    theme_universal(base_size = plot_basesize)

C1

adj <-theme(axis.title = element_text(
    size = plot_basesize))

row2 <- 
    cowplot::plot_grid(
        A1 + adj,
        B1+ adj,
        C1+ adj,
        nrow = 1, 
        labels =  c("C","D","E")
    )

row2


###############################################################
# Panel 4: M1 opsonophagocytosis vs binding antibody titre  #
###############################################################


####################################################
####################################################

# Analysis 1: Compare the proportion of participants with observable OPA of M1 between 
# those above and below protective thresholds with Fisher's exact test 

# Prepare data for plotting by selecting relevant columns
data_plot <- df %>%
    # Select relevant columns for the analysis
    select(M1_THP1_activity, SLO_above_threshold, SpyAD_above_threshold, SpyCEP_above_threshold) %>%
    # Pivot longer for easier aggregation
    tidyr::pivot_longer(cols = c(SLO_above_threshold, SpyAD_above_threshold, SpyCEP_above_threshold),
                        names_to = "Antigen",
                        values_to = "Above_Threshold") %>%
    filter(!is.na(Above_Threshold)) %>%
    # Group by antigen and threshold level
    group_by(Antigen, Above_Threshold) %>%
    # Calculate the percentage of M1_THP1_activity == 1
    summarise(Percentage = mean(M1_THP1_activity == 1, na.rm = TRUE) * 100, .groups = "drop") %>%
    mutate(Antigen = str_remove(Antigen, "_above_threshold"))


# The above transformation ensures that we have a summary of the percentage 
# of participants exhibiting opsonophagocytosis within each antigen group.

# Calculate Fisher's exact test p-values
p_values <- df %>%
    # Pivot data to long format
    tidyr::pivot_longer(cols = c(SLO_above_threshold, SpyAD_above_threshold, SpyCEP_above_threshold),
                        names_to = "Antigen",
                        values_to = "Above_Threshold") %>%
    # Remove rows where Above_Threshold variable is missing
    filter(!is.na(Above_Threshold)) %>%
    # Perform Fisher's exact test for each antigen group
    group_by(Antigen) %>%
    summarise(
        P_value = fisher.test(table(M1_THP1_activity, Above_Threshold))$p.value,
        .groups = "drop"
    ) %>%
    # Clean antigen names for consistency
    mutate(Antigen = str_remove(Antigen, "_above_threshold"))

# Merge percentages and p-values
data_plot <- data_plot %>%
    left_join(p_values, by = "Antigen")

# Display the final summary dataset
data_plot

# Extract the values for the sentence
antigens <- c("SLO", "SpyAD", "SpyCEP")

# Extract percentage values for participants above and below protective thresholds
above_threshold <- data_plot %>%
    filter(Above_Threshold == 1) %>%
    arrange(match(Antigen, antigens)) %>%
    pull(Percentage) %>%
    round()

below_threshold <- data_plot %>%
    filter(Above_Threshold == 0) %>%
    arrange(match(Antigen, antigens)) %>%
    pull(Percentage) %>%
    round()

p_values <- data_plot %>%
    filter(Above_Threshold == 1) %>%
    arrange(match(Antigen, antigens)) %>%
    pull(P_value)

# Format p-values
formatted_p_values <- sapply(p_values, function(p) ifelse(p < 0.0001, "< 0.0001", sprintf("%.5f", p)))

# Construct summary sentence dynamically
sentence <- sprintf(
    "M1 bacterial opsonophagocytosis was observed in %d%%, %d%%, %d%% serum samples with IgG levels above 50%% protective thresholds, compared with %d%%, %d%%, %d%% below protective thresholds for SLO, SpyAD, and SpyCEP respectively (Fisher's exact test, p = %s, %s, %s).",
    above_threshold[1], above_threshold[2], above_threshold[3], 
    below_threshold[1], below_threshold[2], below_threshold[3], 
    formatted_p_values[1], formatted_p_values[2], formatted_p_values[3]
)

# Print the sentence
print(sentence)

##

# Extract values for the sentence
antigens <- c("SLO", "SpyAD", "SpyCEP")
above_threshold <- data_plot %>%
    filter(Above_Threshold == 1) %>%
    arrange(match(Antigen, antigens)) %>%
    pull(Percentage) %>%
    round()

below_threshold <- data_plot %>%
    filter(Above_Threshold == 0) %>%
    arrange(match(Antigen, antigens)) %>%
    pull(Percentage) %>%
    round()

p_values <- data_plot %>%
    filter(Above_Threshold == 1) %>%
    arrange(match(Antigen, antigens)) %>%
    pull(P_value)

# Format p-values
formatted_p_values <- sapply(p_values, function(p) ifelse(p < 0.0001, "< 0.0001", sprintf("%.5f", p)))

# Additional dynamic summary with antigen-specific comparisons
antigen_comparisons <- sapply(1:length(antigens), function(i) {
    sprintf(
        "%s (%d%% vs %d%%, p = %s)",
        antigens[i], above_threshold[i], below_threshold[i], formatted_p_values[i]
    )
})

# Construct a second summary sentence with formatted antigen-wise comparisons
sentence <- sprintf(
    "Opsonophagocytosis of M1 bacteria was observed in a greater proportion of samples above the 50%% protective threshold compared to below it for %s.",
    paste(antigen_comparisons, collapse = "; ")
)

# Print the sentence
print(sentence)


################

## Analysis 2: Visualise the relationship between binding IgG titre to SLO SpyAD and SpyCEP and IC50 from M1 OPA assay. 

plot <- ggplot(data_plot, aes(x = Antigen, y = Percentage, fill = factor(Above_Threshold))) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    labs(
      #  title = "Percentage of M1_THP1 Activity by Antigen Threshold",
        x = "Binding IgG",
        y = "Percentage with M1 opsonophagocytosis (%)",
        fill = "Above 50% Protective Threshold"
    ) +
    geom_text(
        aes(label = ifelse(
            Above_Threshold == 1,
            ifelse(P_value < 0.0001, "p < 0.0001", paste0("p = ", round(P_value, 3))),
            ""
        )),
        size = 5,
        position = position_dodge(width = 0.9),
        vjust = -0.5
    ) +
    scale_fill_manual(values = c("#d73027", "#7570b3")) +
    theme_minimal()+
    theme_universal(base_size = plot_basesize)

plot
E1 <- df %>%
    ggplot(
        aes(
            x = SLO_IgG, y = log10(M1_THP1_IC50))
    )  +
    geom_point(aes(col = as.factor(SLO_above_threshold)), alpha = 0.8,size = 3) +
   # ggpubr::stat_cor(method = "spearman")    +
    theme_minimal() +
    scale_color_manual(values = c("#d73027","#7570b3")) +
    labs(x = "IgG level (log10 RLU/mL)",
         y= "Phagocytosis of M1 bacteria (log10 IC50)",
         title = "SLO"
    ) +
    guides(col = "none")+
    theme_universal(base_size = plot_basesize)
E1

E2 <- df %>%
    ggplot(
        aes(
            x = SpyAD_IgG, y = log10(M1_THP1_IC50))
    )  +
    geom_point(aes(col = as.factor(SpyAD_above_threshold)), alpha = 0.8,size = 3) +
    # ggpubr::stat_cor(method = "spearman")    +
    scale_color_manual(values = c("#d73027","#7570b3")) +
    theme_minimal() +
    labs(
        x = "IgG level (log10 RLU/mL)",
        y="",
    #    y= "Phagocytosis of M1 bacteria (log10 IC50)",
        title = "SpyAD"
    ) +
    guides(col = "none")+
    theme_universal(base_size = plot_basesize)


E3 <- df %>%
    ggplot(
        aes(
            x = SpyCEP_IgG, y = log10(M1_THP1_IC50))
    )  +
    geom_point(aes(col = as.factor(SpyCEP_above_threshold)), alpha = 0.8,size = 3) +
  #  ggpubr::stat_cor(method = "spearman")    +
    scale_color_manual(values = c("#d73027","#7570b3")) +
    theme_minimal() +
    labs(x = "IgG level (log10 RLU/mL)",
         y = "",
      #   y= "Phagocytosis of M1 bacteria (log10 IC50)",
         title = "SpyCEP"
    ) +
    guides(col = "none")+
    theme_universal(base_size = plot_basesize)


panelF <- plot_grid(E1, E2, E3, ncol= 3)

row3 <- plot_grid(plot,panelF, ncol= 2,
                  rel_widths = c(0.4,0.6),
                  labels = c("F","G"))

row3



##### Combine all panels 

main_plot5_V3.0_261124 <- plot_grid(
    row1,
    row2,
    row3,
    ncol = 1,
    rel_heights = c(1,0.8,1)
)

ggsave(filename = "R_output/Main_Fig5_V3.3_120225.png", 
       plot = main_plot5_V3.0_261124,
       width = 2800 / 96, 
       height = 1600/ 96, 
       dpi = 300, 
       bg = "white")
    
main_plot5_V3.0_261124


