# Author: Anjaney J Pandey
# Code for single signature box plot with statistical comparisons between signatures

##############################################################

library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(tidyr)
library(forcats)
library(RColorBrewer)
library(scales)
library(ggpubr)

pan_cancer_results_dir <- "./Transcriptomics/TCGA_DATA/Single_Signature_Survival_Analysis"
output_dir <- pan_cancer_results_dir
single_sig_results_file <- file.path(input_dir, "PanCancer_Single_Signature_OS_HR.tsv")
output_plot_file <- file.path(output_dir, "PanCancer_SingleSig_HR_BoxPlot_withComparisons.png")
p_value_threshold <- 0.1
hr_upper_threshold <- 5
sig_short_forms <- c("EMT", "Gly", "OXP", "FAO", "Hyp", "Fer")
primary_colors <- c("red", "blue", "black", "#CCCC00", "darkgreen", "#FF69B4")
color_palette <- setNames(primary_colors, sig_short_forms)

if (!file.exists(single_sig_results_file)) { stop("Required input file not found: ", single_sig_results_file) }
print(paste("Loading data from:", single_sig_results_file))
single_sig_results <- tryCatch(
    readr::read_tsv(single_sig_results_file, show_col_types = FALSE),
    error = function(e) { stop("Error reading file: ", single_sig_results_file, "\nError: ", conditionMessage(e)) } )
if (nrow(single_sig_results) == 0) { stop("Loaded results file contains no data.") }

plot_data_box <- single_sig_results %>%
  filter( p.value < p_value_threshold, is.finite(HR), HR <= hr_upper_threshold ) %>%
  mutate( signature = factor(signature, levels = sig_short_forms) ) %>%
  filter(!is.na(signature))

max_hr_plot <- max(plot_data_box$HR, na.rm = TRUE)
y_pos_1 <- max_hr_plot * 1.05
y_pos_2 <- max_hr_plot * 1.15

# Box Plot
my_comparisons <- list( c("Gly", "OXP"), c("Gly", "FAO") )

pan_cancer_boxplot_compare <- ggplot(plot_data_box, aes(x = signature, y = HR)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.6) +
  geom_boxplot( aes(fill = signature), outlier.shape = NA, alpha = 0.5, width = 0.6 ) +
  geom_jitter( aes(color = signature), size = 1.5, width = 0.25, alpha = 0.7 ) +
  scale_fill_manual( values = color_palette, name = "Signature", drop = FALSE ) +
  scale_color_manual( values = color_palette, name = "Signature", drop = FALSE ) +
  guides(fill = guide_legend(override.aes = list(alpha=1)), color = "none") +
  stat_compare_means(
      comparisons = my_comparisons,
      method = "wilcox.test",
      label = "p.signif",
      label.y = c(y_pos_1, y_pos_2),
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
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(size=10, face="bold"),
    legend.text = element_text(size=9),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = paste("Pan-Cancer HR Distribution for Individual Signatures"),
    x = "Signature",
    y = "Hazard Ratio (HR)"
  ) +
  coord_cartesian(ylim = c(NA, max_hr_plot * 1.25))


tryCatch({
    ggsave(
        filename = output_plot_file,
        plot = pan_cancer_boxplot_compare,
        width = 9,
        height = 7,
        dpi = 300,
        bg = "white",
        limitsize = FALSE
    )
    print(paste("Box plot with comparisons saved to:", output_plot_file))
    }, error = function(e) {
         warning(paste("Failed to save box plot:", conditionMessage(e)))
    })

print("--- Pan-Cancer Single Signature Box Plot w/ Comparisons Script Finished ---")
