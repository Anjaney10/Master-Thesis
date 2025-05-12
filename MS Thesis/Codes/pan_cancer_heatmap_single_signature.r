# Author: Anjaney J Pandey
# This script is for generating the pan-cancer heatmap when each signature is considered individually for analysis

library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(forcats)
library(RColorBrewer)
library(scales)

input_dir <- "./Transcriptomics/TCGA_DATA/Combinatorial_Survival_AllNegRef/Single_Signature_Survival_Analysis"
output_dir <- input_dir

# Input file
single_sig_results_file <- file.path(input_dir, "PanCancer_Single_Signature_OS_HR.tsv")

# Output filename
output_heatmap_file <- file.path(output_dir, "PanCancer_SingleSig_HR_Heatmap_Filtered.png")

# Filtering thresholds
hr_upper_threshold <- 5

# Signature
sig_short_forms <- c("EMT", "Gly", "OXP", "FAO", "Hyp", "Fer")

# Define Significance Stars

get_significance_stars <- function(p_values) {
  case_when( is.na(p_values) ~ "", p_values < 0.001 ~ "***", p_values < 0.01 ~ "**", p_values < 0.05 ~ "*", TRUE ~ "" )
}

# Load Data
if (!file.exists(single_sig_results_file)) {
  stop("Required input file not found: ", single_sig_results_file)
}
print(paste("Loading data from:", single_sig_results_file))
single_sig_results <- tryCatch(
    readr::read_tsv(single_sig_results_file, show_col_types = FALSE),
    error = function(e) { stop("Error reading file: ", single_sig_results_file, "\nError: ", conditionMessage(e)) } )

if (nrow(single_sig_results) == 0) { stop("Loaded results file contains no data.") }

# Filter data
plot_data_heatmap <- single_sig_results %>%
  filter(
    is.finite(HR),
    HR <= hr_upper_threshold
  ) %>%
  mutate(
    signature = factor(signature, levels = rev(sig_short_forms)),
    cancer_type = factor(cancer_type),
    significance = get_significance_stars(p.value)
  ) %>%
  filter(!is.na(signature) & !is.na(cancer_type))

# Define Color Limits
max_abs_dev <- max(abs(plot_data_heatmap$HR - 1), 0.1, na.rm = TRUE)
color_limits <- c(max(0, 1 - max_abs_dev), 1 + max_abs_dev)


# Create Heatmap
pan_cancer_sig_heatmap <- ggplot(plot_data_heatmap, aes(x = cancer_type, y = signature, fill = HR)) +
  geom_tile(color = "grey85", size=0.2) +
  geom_text(aes(label = significance), color = "black", size = 3.5) +
  scale_fill_gradient2(
        name = "Hazard Ratio (HR)",
        low = muted("blue"),
        mid = "white",
        high = muted("red"),
        midpoint = 1,
        limit = color_limits,
        oob = scales::squish
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9, face = "bold"),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(size=9, face="bold"),
    legend.text = element_text(size=8),
    legend.key.height = unit(1.2, "cm")
  ) +
  labs(
    title = paste("Pan-Cancer Single Signature HR")
  )

# Save Heatmap

plot_width <- max(8, 2 + n_distinct(plot_data_heatmap$cancer_type) * 0.4)
plot_height <- max(5, 1.5 + n_distinct(plot_data_heatmap$signature) * 0.5)

ggsave(
        filename = output_heatmap_file,
        plot = pan_cancer_sig_heatmap,
        width = plot_width,
        height = plot_height,
        dpi = 300,
        bg = "white",
        limitsize = FALSE
    )
