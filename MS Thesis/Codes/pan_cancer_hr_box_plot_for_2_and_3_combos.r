# Author: Anjaney J Pandey
# This script is for genrating box plots showing the distribution of hazard ratios for different signatures in a pan-cancer manner

library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(tidyr)
library(forcats)
library(RColorBrewer)
library(scales)
library(ggpubr)

pan_cancer_results_dir <- "./Transcriptomics/TCGA_DATA/Combinatorial_Survival_AllNegRef/pan_cancer_combinations"
output_dir <- pan_cancer_results_dir

# Input files
two_combo_file <- file.path(pan_cancer_results_dir, "PanCancer_Two_Signature_Combinations_AllNegRef_HR.tsv")
three_combo_file <- file.path(pan_cancer_results_dir, "PanCancer_Three_Signature_Combinations_AllNegRef_HR.tsv")

# Output filenames
output_plot_file_2combo <- file.path(output_dir, "PanCancer_TwoCombo_SigHR_BoxPlot_Comparisons.png")
output_plot_file_3combo <- file.path(output_dir, "PanCancer_ThreeCombo_SigHR_BoxPlot_Comparisons.png")

# Filtering thresholds
p_value_threshold <- 0.05
hr_upper_threshold <- 10

# Signature
sig_short_forms <- c("EMT", "Gly", "OXP", "FAO", "Hyp", "Fer")

# Define Colors
primary_colors <- c("red", "blue", "black", "#CCCC00", "darkgreen", "#FF69B4")
color_palette <- setNames(primary_colors, sig_short_forms)

# Function to parse POSITIVE signatures from the comparison term label such as EMT+/Gly+/Hyp-
parse_positive_signatures <- function(term_label, all_short_forms) {
  if (is.na(term_label)) return(list(NA_character_)[-1])
  parts <- str_split(term_label, "/", simplify = TRUE)[1,]
  positive_sigs <- parts[str_detect(parts, "\\+$")]
  positive_sigs_clean <- gsub("\\+$", "", positive_sigs)
  intersect(positive_sigs_clean, all_short_forms)
}

# Main function for plotting
create_combo_boxplot_with_comparisons <- function(input_data,
                                                  plot_title,
                                                  comparisons_list,
                                                  output_filename,
                                                  p_filter = p_value_threshold,
                                                  hr_limit = hr_upper_threshold)
 {

  # Filter data
  plot_data_long <- input_data %>%
    filter( p.value < p_filter, is.finite(HR), HR <= hr_limit ) %>%
    mutate( plot_labels = gsub("^Combination_Level", "", term),
            plot_labels = ifelse(plot_labels == "", term, plot_labels) ) %>%
    filter(!is.na(plot_labels) & plot_labels != "") %>%
    rowwise() %>%
    mutate( positive_signatures_list = list(parse_positive_signatures(plot_labels, sig_short_forms)) ) %>%
    ungroup() %>%
    tidyr::unnest(positive_signatures_list, keep_empty = FALSE) %>%
    rename(positive_signature = positive_signatures_list) %>%
    mutate( positive_signature = factor(positive_signature, levels = sig_short_forms),
            cancer_type = factor(cancer_type) ) %>%
    filter(!is.na(positive_signature))

  # Calculate max HR
  max_hr_plot <- max(plot_data_long$HR, na.rm = TRUE)
  num_comparisons <- length(comparisons_list)
  y_positions <- seq(from = max_hr_plot * 1.05,
                     by = max_hr_plot * 0.10,
                     length.out = num_comparisons)


  # Create Box Plot
  p_boxplot <- ggplot(plot_data_long, aes(x = positive_signature, y = HR)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.6) +
    geom_boxplot( aes(fill = positive_signature), outlier.shape = NA, alpha = 0.5, width = 0.6 ) +
    geom_jitter( aes(color = positive_signature), size = 1.5, width = 0.25, alpha = 0.7 ) +
    scale_fill_manual( values = color_palette, name = "Signature", drop = FALSE ) +
    scale_color_manual( values = color_palette, name = "Signature", drop = FALSE ) +
    guides(fill = guide_legend(override.aes = list(alpha=1)), color = "none") +

    # Add Pairwise Comparisons
    stat_compare_means(
        comparisons = comparisons_list,
        method = "wilcox.test",
        label = "p.format",
        label.y = y_positions,
        bracket.size = 0.4,
        step.increase = 0,
        tip.length = 0.015
        ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face="bold"),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 11),
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.title = element_text(size=10, face="bold"),
      legend.text = element_text(size=9),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = plot_title,
      x = "Positive Signature in Comparison Term",
      y = "Hazard Ratio (HR)"
    ) +
    coord_cartesian(ylim = c(NA, max(y_positions) * 1.1))


  # Save Plot
  plot_height <- 6 + (num_comparisons * 0.5)
  plot_width <- 9

  ggsave( filename = output_filename, plot = p_boxplot, width = plot_width,
              height = plot_height, dpi = 300, bg = "white", limitsize = FALSE)
}

# Generate Plot for 2-Signature Combinations
    combo2_data <- read_tsv(two_combo_file, show_col_types = FALSE)

    create_combo_boxplot_with_comparisons(
        input_data = combo2_data,
        plot_title = "Pan-Cancer HR Distrib. for Sig. 2-Signature Combos",
        comparisons_list = list(c("Gly", "OXP"), c("Gly", "FAO")),
        output_filename = output_plot_file_2combo,
        p_filter = p_value_threshold,
        hr_limit = hr_upper_threshold
    )

# Generate Plot for 3-Signature Combinations
    combo3_data <- read_tsv(three_combo_file, show_col_types = FALSE)

    create_combo_boxplot_with_comparisons(
        input_data = combo3_data,
        plot_title = "Pan-Cancer HR Distrib. for Sig. 3-Signature Combos",
        comparisons_list = list(c("Gly", "OXP"), c("Gly", "FAO")),
        output_filename = output_plot_file_3combo,
        p_filter = p_value_threshold,
        hr_limit = hr_upper_threshold
    )
