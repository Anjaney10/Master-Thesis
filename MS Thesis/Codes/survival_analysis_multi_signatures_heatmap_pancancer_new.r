# Author: Anjaney J Pandey
# FINAL CODE FOR PAN CANCER SURVIVAL METRICS using multiple signatures - Heatmap from already calculated Hazard Ratios for 2 and 3 signature combinations
########################################################################

library(dplyr)
library(ggplot2)
library(scales)
library(readr)
library(stringr)

get_significance_stars <- function(p_values) {
  case_when(
    is.na(p_values) ~ "",
    p_values < 0.001 ~ "***",
    p_values < 0.01 ~ "**",
    p_values < 0.05 ~ "*",
    TRUE ~ ""
  )
}

# Main Function to Generate Heatmaps
generate_pan_cancer_heatmaps_from_file <- function(combined_results_file,
                                                   output_dir,
                                                   conf_high_limit = 7,
                                                   survival_metric_name = "OS")
{
  if (!file.exists(combined_results_file)) {
    stop("Combined results file not found: ", combined_results_file)
  }
  required_cols <- c("combination", "cancer_type", "term", "HR", "conf.high", "p.value")
  results_data_full <- tryCatch(
       readr::read_tsv(combined_results_file, show_col_types = FALSE),
       error = function(e) { stop("Failed to read results file: ", conditionMessage(e)) }
   )
  if (!all(required_cols %in% colnames(results_data_full))) {
      missing_cols <- setdiff(required_cols, colnames(results_data_full))
      stop("Results file missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  results_data_filtered <- results_data_full %>%
    filter(is.finite(HR) & is.finite(conf.high) & conf.high <= conf_high_limit)

  if (nrow(results_data_filtered) == 0) {
    warning("No data remaining after filtering conf.high <= ", conf_high_limit, " in file: ", basename(combined_results_file))
    return(invisible(NULL))
  }

  results_data_filtered <- results_data_filtered %>%
      mutate(significance = get_significance_stars(p.value))

  unique_combinations <- unique(results_data_filtered$combination)

  for (current_combo_name in unique_combinations) {
    print(paste("Generating heatmap for Combination:", current_combo_name))

    combo_data <- results_data_filtered %>%
      filter(combination == current_combo_name)

    if (nrow(combo_data) == 0) {
        warning("No data for combination '", current_combo_name, "' after filtering. Skipping heatmap.")
        next
    }

    combo_data <- combo_data %>%
        mutate(
            plot_labels = gsub("^Combination_Level", "", term),
            plot_labels = ifelse(plot_labels == "", term, plot_labels),
            plot_labels = factor(plot_labels, levels = rev(sort(unique(plot_labels)))),
            cancer_type = factor(cancer_type)
        ) %>%
        filter(!is.na(plot_labels) & plot_labels != "")
    if(nrow(combo_data) == 0) {
        warning("No data remaining for '", current_combo_name, "' after cleaning labels. Skipping heatmap.")
        next
    }

    max_abs_dev <- max(abs(combo_data$HR - 1), 0.1, na.rm = TRUE)
    color_limits <- c(max(0, 1 - max_abs_dev), 1 + max_abs_dev)


    # Create Heatmap
    heatmap_plot <- ggplot(combo_data, aes(x = cancer_type, y = plot_labels, fill = HR)) +
      geom_tile(color = "grey90", size=0.2) +
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
        axis.text.y = element_text(size = 9),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(size=9),
        legend.text = element_text(size=8),
        legend.key.height = unit(1.2, "cm")
      ) +
      labs(
        title = paste("Pan-Cancer HR:", current_combo_name)
      )

    # Save Heatmap
    heatmap_filename <- file.path(output_dir, paste0("PanCancer_", current_combo_name, "_", survival_metric_name, "_Heatmap.png"))

    tryCatch({
        ggsave(
            plot = heatmap_plot,
            filename = heatmap_filename,
            width = 12,
            height = 5,
            dpi = 300,
            bg = "white",
            limitsize = FALSE
        )
        print(paste(".... Heatmap saved:", basename(heatmap_filename)))
    }, error = function(e) {
        warning(paste("Failed to save heatmap for", current_combo_name, ":", conditionMessage(e)))
    })

  }

  print(paste("Heatmap generation complete for file:", basename(combined_results_file)))
  return(invisible(NULL))
}


pan_cancer_results_dir <- "./Transcriptomics/TCGA_DATA/MR_Combination_Survival_Analysis/Combinatorial_Survival_Analysis_New/pan_cancer_combinations"
pan_cancer_heatmap_dir <- pan_cancer_results_dir

two_combo_file <- file.path(pan_cancer_results_dir, "PanCancer_Two_Signature_Combinations_HR.tsv")
three_combo_file <- file.path(pan_cancer_results_dir, "PanCancer_Three_Signature_Combinations_HR.tsv")

# Generate heatmaps for 2-signature combinations
if (file.exists(two_combo_file)) {
  generate_pan_cancer_heatmaps_from_file(
    combined_results_file = two_combo_file,
    output_dir = pan_cancer_heatmap_dir,
    conf_high_limit = 10,
    survival_metric_name = "OS"
  )
} else {
  warning("File not found, cannot generate 2-signature heatmaps: ", two_combo_file)
}

# Generate heatmaps for 3-signature combinations
if (file.exists(three_combo_file)) {
  generate_pan_cancer_heatmaps_from_file(
    combined_results_file = three_combo_file,
    output_dir = pan_cancer_heatmap_dir,
    conf_high_limit = 20,
    survival_metric_name = "OS"
  )
} else {
  warning("File not found, cannot generate 3-signature heatmaps: ", three_combo_file)
}

print("--- Pan-Cancer Heatmap Script Finished")
