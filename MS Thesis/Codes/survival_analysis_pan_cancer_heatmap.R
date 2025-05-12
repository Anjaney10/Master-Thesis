# Author: Anjaney J Pandey
# FINAL CODE FOR PAN CANCER SURVIVAL METRICS - Heatmap

######################################################### Heatmap plot ########################################################
####################################################################

library(dplyr)
library(ggplot2)
library(scales)

generate_survival_heatmaps <- function(results_file_path,
                                       base_output_dir,
                                       surv_metrics)
{

  if (!file.exists(results_file_path)) {
    stop("Results file not found: ", results_file_path)
  }
  required_cols <- c("metric", "cancer_type", "labels", "HR", "p.value", "HR.confint.lower", "HR.confint.upper")
  results_data <- tryCatch(
       read.delim(results_file_path, sep = "\t", stringsAsFactors = FALSE),
       error = function(e) { stop("Failed to read results file: ", conditionMessage(e)) }
   )
  if (!all(required_cols %in% colnames(results_data))) {
      missing_cols <- setdiff(required_cols, colnames(results_data))
      stop("Results file missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  #  Data Pre-processing
  results_data$plot_labels <- gsub("^Type", "", results_data$labels)
  comparison_levels <- c("RPMS+/BPMS+", "RPMS+/BPMS-", "RPMS-/BPMS+")
  results_data$plot_labels <- factor(results_data$plot_labels, levels = rev(comparison_levels))
  results_data$cancer_type <- factor(results_data$cancer_type)
  results_data$significance <- sapply(results_data$p.value, get_significance_stars)


  for (metric_id in names(surv_metrics)) {
    metric_info <- surv_metrics[[metric_id]]
    metric_pretty_name <- metric_info$pretty_name
    metric_dir_name <- metric_info$dir_name

    print(paste("Generating heatmap for Metric:", metric_pretty_name))

    metric_data_initial <- results_data %>%
      filter(metric == metric_id, !is.na(plot_labels))

    if (nrow(metric_data_initial) == 0) {
      warning("No data found for metric: ", metric_pretty_name, ". Skipping heatmap.")
      next
    }
    print(paste("Initial rows for", metric_pretty_name, ":", nrow(metric_data_initial)))

    cancers_with_invalid <- metric_data_initial %>%
        filter(!is.finite(HR) | !is.finite(HR.confint.lower) | !is.finite(HR.confint.upper)) %>%
        distinct(cancer_type) %>%
        pull(cancer_type)

    if (length(cancers_with_invalid) > 0) {
        print(paste("Removing cancer types with non-finite HR/CI for", metric_pretty_name, ":", paste(cancers_with_invalid, collapse=", ")))
        metric_data_filtered <- metric_data_initial %>%
            filter(!(cancer_type %in% cancers_with_invalid))
    } else {
        print("No cancer types found with non-finite HR/CI values.")
        metric_data_filtered <- metric_data_initial
    }
    metric_data_filtered$cancer_type <- droplevels(metric_data_filtered$cancer_type)


    #  Create Heatmap using metric_data_filtered
    max_hr_val <- max(metric_data_filtered$HR, na.rm = TRUE)
    min_hr_val <- min(metric_data_filtered$HR, na.rm = TRUE)
    color_limit <- max(abs(1 - min_hr_val), abs(max_hr_val - 1), 0.1, na.rm = TRUE)

      heatmap_plot <- ggplot(metric_data_filtered, aes(x = cancer_type, y = plot_labels, fill = HR)) +
      geom_tile(color = "grey85", size=0.3) +
      geom_text(aes(label = significance), color = "black", size = 5, vjust = 0.6) +
      scale_fill_gradient2(
            name = "Hazard Ratio (HR)",
            low = muted("blue"),
            mid = "white",
            high = muted("red"),
            midpoint = 1,
            limit = c(max(0, 1 - color_limit), 1 + color_limit),
            oob = scales::squish
      ) +
      guides(fill = guide_colorbar(
          title.position = "top",
          title.hjust = 0.5,
          barheight = unit(15, "lines"),
          barwidth = unit(1.5, "lines")
      )) +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text.align = 0.5,
        legend.text = element_text(size = 10)
      ) +
      labs(
        title = metric_pretty_name
      )


    #  Save Heatmap
    output_subdir <- file.path(base_output_dir, metric_dir_name)
    dir.create(output_subdir, showWarnings = FALSE, recursive = TRUE)
    heatmap_filename <- file.path(output_subdir, paste0(metric_dir_name, "_HR_heatmap_new.png"))
    plot_height <- max(7, 1 + nlevels(metric_data_filtered$cancer_type) * 0.8)
    plot_width <- max(5, 1.5 + nlevels(metric_data_filtered$plot_labels) * 0.6)

    print(paste("Saving heatmap to:", heatmap_filename))
    tryCatch(
        ggsave(
            plot = heatmap_plot,
            filename = heatmap_filename,
            width = 12,
            height = 4,
            dpi = 300,
            bg = "white"
        ),
        error = function(e) {
             warning("Failed to save heatmap for ", metric_pretty_name, ": ", conditionMessage(e))
        }
    )

  }

  return(invisible(NULL))
}

#  Function Call Arguments
surv_metrics_config <- list(
   OS = list(status = "OS", time = "OS.time", pretty_name = "Overall Survival", dir_name = "overall_survival"),
   PFI = list(status = "PFI", time = "PFI.time", pretty_name = "Progression Free Interval", dir_name = "progression_free_interval"),
   DFI = list(status = "DFI", time = "DFI.time", pretty_name = "Disease Free Interval", dir_name = "disease_free_interval"),
   DSS = list(status = "DSS", time = "DSS.time", pretty_name = "Disease Specific Survival", dir_name = "disease_specific_survival")
)

results_file <- "./Transcriptomics/TCGA_DATA/PAN_CANCER_Survival_Analysis_RPMS_BPMS/tpm_HR_all_metrics_results.tsv"
heatmap_output_dir <- "./Transcriptomics/TCGA_DATA/PAN_CANCER_Survival_Analysis_RPMS_BPMS"

# Run the function
generate_survival_heatmaps(
   results_file_path = results_file,
   base_output_dir = heatmap_output_dir,
   surv_metrics = surv_metrics_config
)
