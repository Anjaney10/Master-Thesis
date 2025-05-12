# Author: Anjaney J Pandey
# This code is used for generating box plots for the scenarios when only one signature was positive (meaning the activity of that particular pathway was high) but the other signatures had low activity (negative sign).
# This is done to look at the impact of the individual signatures on the hazard ratios when they are in combinations, which couldn't be captured otherwise

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

two_combo_file <- file.path(pan_cancer_results_dir, "PanCancer_Two_Signature_Combinations_AllNegRef_HR.tsv")
three_combo_file <- file.path(pan_cancer_results_dir, "PanCancer_Three_Signature_Combinations_AllNegRef_HR.tsv")

output_plot_file_2combo <- file.path(output_dir, "PanCancer_TwoCombo_SinglePosSigHR_BoxPlot_Comparisons.png")
output_plot_file_3combo <- file.path(output_dir, "PanCancer_ThreeCombo_SinglePosSigHR_BoxPlot_Comparisons.png")

p_value_threshold <- 0.05
hr_upper_threshold <- 5
sig_short_forms <- c("EMT", "Gly", "OXP", "FAO", "Hyp", "Fer")
primary_colors <- c("red", "blue", "black", "#CCCC00", "darkgreen", "#FF69B4")
color_palette <- setNames(primary_colors, sig_short_forms)

parse_positive_signatures <- function(term_label, all_short_forms) {
  if (is.na(term_label)) return(list(NA_character_)[-1])
  parts <- str_split(term_label, "/", simplify = TRUE)[1,]; positive_sigs <- parts[str_detect(parts, "\\+$")]
  positive_sigs_clean <- gsub("\\+$", "", positive_sigs); intersect(positive_sigs_clean, all_short_forms) }

# Main Plotting Function
create_combo_boxplot_with_comparisons <- function(input_data,
                                                  plot_title,
                                                  comparisons_list,
                                                  output_filename,
                                                  p_filter = p_value_threshold,
                                                  hr_limit = hr_upper_threshold)
{

    # Calculate max HR
    max_hr_plot <- max(input_data$HR, na.rm = TRUE)
    num_comparisons <- length(comparisons_list)
    y_positions <- seq(from = max_hr_plot * 1.05, by = max_hr_plot * 0.10, length.out = num_comparisons)

    # Create Box Plot
    p_boxplot <- ggplot(input_data, aes(x = positive_signature, y = HR)) +
        geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.6) +
        geom_boxplot( aes(fill = positive_signature), outlier.shape = NA, alpha = 0.5, width = 0.6 ) +
        geom_jitter( aes(color = positive_signature), size = 1.5, width = 0.25, alpha = 0.7 ) +
        scale_fill_manual( values = color_palette, name = "Signature", drop = FALSE ) +
        scale_color_manual( values = color_palette, name = "Signature", drop = FALSE ) +
        guides(fill = guide_legend(override.aes = list(alpha=1)), color = "none") +
        stat_compare_means( comparisons = comparisons_list, method = "wilcox.test", label = "p.signif", size = 3.5,
                            label.y = y_positions, bracket.size = 0.4, step.increase = 0, tip.length = 0.015 ) +
        theme_bw(base_size = 11) +
        theme( axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face="bold"), axis.text.y = element_text(size = 10),
               axis.title = element_text(size = 11), plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
               legend.position = "right", legend.title = element_text(size=10, face="bold"), legend.text = element_text(size=9),
               panel.grid.major.x = element_blank(), panel.grid.minor = element_blank() ) +
        labs( title = plot_title, x = "Signature", y = "Hazard Ratio (HR)" ) +
        coord_cartesian(ylim = c(NA, max(y_positions) * 1.1))

    # Save Plot
    plot_height <- 6 + (num_comparisons * 0.5)
    plot_width <- 9
    ggsave( filename = output_filename, plot = p_boxplot, width = plot_width, height = plot_height,
                dpi = 300, bg = "white", limitsize = FALSE )
}


# Define comparisons
comparisons_to_make <- list(c("Gly", "OXP"), c("Gly", "FAO"), c(("Gly", "Fer")))

# FOr 2-Combo File
    combo2_data_raw <- read_tsv(two_combo_file, show_col_types = FALSE)

    plot_data_2combo_single_pos <- combo2_data_raw %>%
            filter( p.value < p_value_threshold, is.finite(HR), HR <= hr_upper_threshold ) %>%
            mutate( plot_labels = gsub("^Combination_Level", "", term),
                    plot_labels = ifelse(plot_labels == "", term, plot_labels) ) %>%
            filter(!is.na(plot_labels) & plot_labels != "") %>%
            rowwise() %>%
            mutate( positive_signatures = list(parse_positive_signatures(plot_labels, sig_short_forms)) ) %>%
            ungroup() %>%

            # Important: Keep ONLY rows where EXACTLY ONE signature was positive

            filter(lengths(positive_signatures) == 1) %>%
            mutate( positive_signature = map_chr(positive_signatures, ~ .x[1]),
                    positive_signature = factor(positive_signature, levels = sig_short_forms)) %>%
            filter(!is.na(positive_signature))

    create_combo_boxplot_with_comparisons(
                input_data = plot_data_2combo_single_pos,
                plot_title = "Pan-Cancer HR (2-signature combinations)",
                comparisons_list = comparisons_to_make,
                output_filename = output_plot_file_2combo,
                p_filter = p_value_threshold,
                hr_limit = hr_upper_threshold
            )

# For 3-Combo File
    combo3_data_raw <- read_tsv(three_combo_file, show_col_types = FALSE)

    plot_data_3combo_single_pos <- combo3_data_raw %>%
            filter( p.value < p_value_threshold, is.finite(HR), HR <= hr_upper_threshold ) %>%
            mutate( plot_labels = gsub("^Combination_Level", "", term),
                    plot_labels = ifelse(plot_labels == "", term, plot_labels) ) %>%
            filter(!is.na(plot_labels) & plot_labels != "") %>%
            rowwise() %>%
            mutate( positive_signatures = list(parse_positive_signatures(plot_labels, sig_short_forms)) ) %>%
            ungroup() %>%

            # Important: Keep ONLY rows where EXACTLY ONE signature was positive

            filter(lengths(positive_signatures) == 1) %>%
            mutate( positive_signature = map_chr(positive_signatures, ~ .x[1]),
                    positive_signature = factor(positive_signature, levels = sig_short_forms)) %>%
            filter(!is.na(positive_signature))

    create_combo_boxplot_with_comparisons(
                input_data = plot_data_3combo_single_pos,
                plot_title = "Pan-Cancer HR (3-signature combinations)",
                comparisons_list = comparisons_to_make,
                output_filename = output_plot_file_3combo,
                p_filter = p_value_threshold,
                hr_limit = hr_upper_threshold
            )
