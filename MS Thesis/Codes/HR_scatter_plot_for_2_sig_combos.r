# Author: Anjaney J Pandey
# Scatter Plot of Hazard Ratio for PAN CANCER SURVIVAL METRICS

library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(tidyr)
library(forcats)
library(RColorBrewer)

pan_cancer_results_dir <- "./Transcriptomics/TCGA_DATA/MR_Combination_Survival_Analysis/Combinatorial_Survival_Analysis/pan_cancer_combinations"
output_dir <- pan_cancer_results_dir
two_combo_file <- file.path(pan_cancer_results_dir, "PanCancer_Two_Signature_Combinations_AllNegRef_HR.tsv")
output_plot_file <- file.path(output_dir, "PanCancer_TwoCombo_SigHR_ScatterPlot.png")
p_value_threshold <- 0.05
hr_upper_threshold <- 10
sig_short_forms <- c("EMT", "Gly", "OXP", "FAO", "Hyp", "Fer")

primary_colors <- c("red", "blue", "black", "pink", "orange", "yellow")
color_palette <- setNames(primary_colors, sig_short_forms)

if (!file.exists(two_combo_file)) { stop("Required input file not found: ", two_combo_file) }
print(paste("Loading data from:", two_combo_file))
combined_results <- tryCatch(
    readr::read_tsv(two_combo_file, show_col_types = FALSE) %>% mutate(combo_size = 2),
    error = function(e) { stop("Error reading file: ", two_combo_file, "\nSpecific Error: ", conditionMessage(e)) } )
if (nrow(combined_results) == 0) { stop("Loaded results file contains no data.") }

get_first_positive_sig <- function(term_label, all_short_forms) {
  parts <- str_split(term_label, "/", simplify = TRUE)[1,]; positive_sigs <- parts[str_detect(parts, "\\+$")]
  if (length(positive_sigs) > 0) { first_pos <- gsub("\\+$", "", positive_sigs[1])
    if (first_pos %in% all_short_forms) { return(first_pos) } }; return(NA_character_) }

plot_data <- combined_results %>%
  filter( p.value < p_value_threshold, is.finite(HR), 1 < HR,  HR < hr_upper_threshold ) %>%
  mutate( plot_labels = gsub("^Combination_Level", "", term), plot_labels = ifelse(plot_labels == "", term, plot_labels) ) %>%
  filter(!is.na(plot_labels) & plot_labels != "") %>%
  rowwise() %>%
  mutate( color_sig = get_first_positive_sig(plot_labels, sig_short_forms) ) %>%
  ungroup() %>%
  filter(!is.na(color_sig)) %>%
  mutate( color_sig = factor(color_sig, levels = sig_short_forms), cancer_type = factor(cancer_type) )

pan_cancer_scatter_2combo <- ggplot(plot_data, aes(x = cancer_type, y = HR)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey70", linewidth = 0.5) +
  geom_jitter(
    aes(color = color_sig),
    size = 2,
    width = 0.3,
    shape = 16
  ) +
  scale_color_manual(
    values = color_palette,
    name = "Positive Signature",
    na.translate = FALSE,
    drop = FALSE
  ) +
  guides(color = guide_legend(override.aes = list(size=3))) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(size=9),
    legend.text = element_text(size=8),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = paste("Pan-Cancer Significant 2-Signature HRs (p <", p_value_threshold, ", HR <", hr_upper_threshold, ")"),
    subtitle = "Color indicates the first positive signature found in the comparison",
    x = "Cancer Type",
    y = "Hazard Ratio (HR)"
  )

# Save Plot
plot_height <- max(6, 0.2 * n_distinct(plot_data$cancer_type))
plot_width <- 10
tryCatch({
    ggsave( filename = output_plot_file, plot = pan_cancer_scatter_2combo, width = plot_width, height = plot_height, dpi = 300, bg = "white", limitsize = FALSE )
    print(paste("Scatter plot saved to:", output_plot_file))
    }, error = function(e) { warning(paste("Failed to save scatter plot:", conditionMessage(e))) })

print("--- Pan-Cancer 2-Combo Scatter Plot Script Finished ---")
