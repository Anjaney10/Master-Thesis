# Author: Anjaney J Pandey
# FINAL CODE FOR PAN CANCER SURVIVAL METRICS - Forest Plot
# We are plotting the forest plot in batch size of 3 as there are 33 cancer types in TCGA

##############################################################################################

library(dplyr)
library(ggplot2)
library(patchwork)
library(survival)
library(scales)

# Define the significance stars function
get_significance_stars <- function(p_value) {
  if (is.na(p_value)) return("")
  if (p_value < 0.001) return("***")
  if (p_value < 0.01) return("**")
  if (p_value < 0.05) return("*")
  return("")
}

# Function to generate plots for multiple survival metrics
plot_survival_metrics <- function(cancer_types,
                                  batch_size = 3,
                                  base_output_dir = "./Transcriptomics/TCGA_DATA/tpm_all_survival_plots",
                                  results_file = "./Transcriptomics/TCGA_DATA/tpm_HR_all_metrics_results.tsv",
                                  survival_data_dir = "./Transcriptomics/tcga/survival_data/",
                                  ssgsea_data_dir = "./Transcriptomics/TCGA_DATA/tpm_ssgsea_results/") {

  # Survival Metrics
  surv_metrics <- list(
    OS = list(status = "OS", time = "OS.time", pretty_name = "Overall Survival", dir_name = "overall_survival"),
    PFI = list(status = "PFI", time = "PFI.time", pretty_name = "Progression Free Interval", dir_name = "progression_free_interval"),
    DFI = list(status = "DFI", time = "DFI.time", pretty_name = "Disease Free Interval", dir_name = "disease_free_interval"),
    DSS = list(status = "DSS", time = "DSS.time", pretty_name = "Disease Specific Survival", dir_name = "disease_specific_survival")
  )

  # Store ALL results
  results_list_all <- list()

  for (i in seq(1, length(cancer_types), by = batch_size)) {
    end_idx <- min(i + batch_size - 1, length(cancer_types))
    print(paste("Processing cancers", i, "to", end_idx, ":", paste(cancer_types[i:end_idx], collapse=", ")))
    cancer_types_subset <- cancer_types[i:end_idx]

    plot_lists_by_metric_batch <- list()

    for (cancer in cancer_types_subset) {
      print(paste("  Processing Cancer:", cancer))
      survival_file <- file.path(survival_data_dir, paste0(cancer, "_survival.txt"))
      ssgsea_file <- file.path(ssgsea_data_dir, paste0("TCGA-", cancer, "_ssgsea_results.tsv"))

      if (!file.exists(survival_file)) {
        warning(paste("Survival file missing for", cancer, ":", survival_file))
        next
      }
      if (!file.exists(ssgsea_file)) {
        warning(paste("ssGSEA file missing for", cancer, ":", ssgsea_file))
        next
      }

      merged_data <- tryCatch({
          survival_data <- read.table(survival_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

          ssgsea_data_temp <- read.delim(ssgsea_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL)
          sample_col_name <- colnames(ssgsea_data_temp)[1]


               ssgsea_data_temp <- read.delim(ssgsea_file, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL)
               header_line <- readLines(ssgsea_file, n=1)
               actual_headers <- strsplit(header_line, "\t")[[1]]
               colnames(ssgsea_data_temp) <- c("SampleID", actual_headers[-1])
               sample_col_name <- "SampleID"

          ssgsea_data_temp[[sample_col_name]] <- gsub("\\.", "-", ssgsea_data_temp[[sample_col_name]])
          survival_data$sample <- gsub("\\.", "-", survival_data$sample)

          ssgsea_data_temp <- ssgsea_data_temp[!duplicated(ssgsea_data_temp[[sample_col_name]]), ]
          rownames(ssgsea_data_temp) <- ssgsea_data_temp[[sample_col_name]]
          ssgsea_data <- ssgsea_data_temp[, !(colnames(ssgsea_data_temp) %in% sample_col_name), drop = FALSE]

          survival_data <- survival_data[!duplicated(survival_data$sample), ]

          # Merge
          merge(survival_data, ssgsea_data, by.x = "sample", by.y = "row.names", all = FALSE)

      }, error = function(e) {
          warning(paste("Error loading/merging data for", cancer, ":", conditionMessage(e)))
          return(NULL)
      })

      merged_data$Type <- as.character(merged_data$Type)
      expected_levels <- c("RPMS-/BPMS-", "RPMS+/BPMS+", "RPMS+/BPMS-", "RPMS-/BPMS+")
      merged_data <- merged_data[merged_data$Type %in% expected_levels, ]
      if (nrow(merged_data) < 2) {
          warning(paste("Not enough data after filtering by Type for", cancer))
          next
      }
      merged_data$Type <- factor(merged_data$Type, levels = expected_levels)
      if (length(levels(droplevels(merged_data$Type))) < 2) {
          warning(paste("Not enough levels in 'Type' after filtering for Cox model in", cancer))
          next
      }

      for (metric_id in names(surv_metrics)) {
          metric_info <- surv_metrics[[metric_id]]
          metric_status_col <- metric_info$status
          metric_time_col <- metric_info$time
          metric_pretty_name <- metric_info$pretty_name
          metric_dir_name <- metric_info$dir_name

          print(paste("    Metric:", metric_pretty_name))

          results_metric <- NULL
          tryCatch({
              if (!all(c(metric_status_col, metric_time_col) %in% colnames(merged_data))) {
                  warning(paste("Survival columns", metric_status_col, "or", metric_time_col, "missing for", cancer))
                  next
              }
              relevant_cols <- c(metric_time_col, metric_status_col, "Type")
              complete_cases_idx <- complete.cases(merged_data[, relevant_cols])
              if(sum(complete_cases_idx) < 2) {
                  warning(paste("Not enough complete cases for", metric_pretty_name, "in", cancer))
                  next
              }
              merged_data_cox <- merged_data[complete_cases_idx, ]

              if (length(unique(merged_data_cox[[metric_status_col]])) < 2 || length(unique(merged_data_cox[[metric_time_col]])) < 2) {
                  warning(paste("Insufficient variability in data for", metric_pretty_name, "in", cancer))
                  next
              }
              merged_data_cox$Type <- droplevels(merged_data_cox$Type)
              if (length(levels(merged_data_cox$Type)) < 2) {
                  warning(paste("Not enough levels in 'Type' after NA removal for", metric_pretty_name, "in", cancer))
                  next
              }

              # Fit Cox model
              formula <- as.formula(paste0("Surv(", metric_time_col, ", ", metric_status_col, ") ~ Type"))
              fit.coxph <- coxph(formula, data = merged_data_cox)
              s <- summary(fit.coxph)

              if(is.null(s$coefficients) || nrow(s$coefficients) == 0) {
                  warning(paste("Cox model summary empty for", metric_pretty_name, "in", cancer))
                  next
              }

              # Extract results
              results_metric <- data.frame(
                  metric = metric_id,
                  HR = signif(s$coefficients[, "exp(coef)"], 3),
                  HR.confint.lower = signif(s$conf.int[, "lower .95"], 3),
                  HR.confint.upper = signif(s$conf.int[, "upper .95"], 3),
                  p.value = signif(s$coefficients[, "Pr(>|z|)"], 3),
                  labels = rownames(s$coefficients),
                  Log_Rank_pvalue = signif(s$logtest[['pvalue']], 3),
                  n = s$n,
                  cancer_type = cancer
              results_metric$significance <- sapply(results_metric$p.value, get_significance_stars)

              results_list_all[[length(results_list_all) + 1]] <- results_metric

          }, error = function(e) {
              warning(paste("Error during Cox model for", metric_pretty_name, "in", cancer, ":", conditionMessage(e)))
          })

          if (!is.null(results_metric) && nrow(results_metric) > 0) {

              plot_title_string <- cancer
              log_rank_p <- results_metric$Log_Rank_pvalue[1]
              if (!is.na(log_rank_p)) {
                    plot_title_string <- paste0(cancer," ", metric_pretty_name, "\nLogRank p = ", signif(log_rank_p, 3))
              }
              results_metric$plot_labels <- gsub("^Type", "", results_metric$labels)
              label_order <- unique(results_metric$plot_labels)
              results_metric$plot_labels <- factor(results_metric$plot_labels, levels = rev(label_order))

              p <- ggplot(results_metric, aes(x = plot_labels, y = HR)) +
                  geom_pointrange(aes(ymin = HR.confint.lower, ymax = HR.confint.upper),
                                  color = "black",
                                  size = 0.8) +
                  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
                  geom_text(aes(label = significance, y = HR.confint.upper),
                            color = "black", vjust = -0.5, hjust = 0.5, size = 5) +
                  coord_flip() +
                  theme_bw(base_size = 14) +
                  theme(
                      legend.position = "none",
                      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, lineheight = 1.1),
                      axis.text.x = element_text(size = 14),
                      axis.text.y = element_text(size = 14),
                      axis.title.x = element_text(size = 16),
                      axis.title.y = element_blank(),
                      plot.margin = margin(t = 15, r = 10, b = 5, l = 10, unit = "pt")
                  ) +
                  labs(
                      title = plot_title_string,
                      x = NULL,
                      y = "Hazard Ratio (HR)"
                  )
              if (!metric_id %in% names(plot_lists_by_metric_batch)) {
                  plot_lists_by_metric_batch[[metric_id]] <- list()
              }
              plot_lists_by_metric_batch[[metric_id]][[cancer]] <- p

          }

      }
    }


    #  Combine and Save Plots
    if (length(plot_lists_by_metric_batch) > 0) {
        for(metric_id in names(plot_lists_by_metric_batch)) {
            metric_info <- surv_metrics[[metric_id]]
            current_metric_plot_list <- plot_lists_by_metric_batch[[metric_id]]

            if (length(current_metric_plot_list) > 0) {
                print(paste("  Combining and saving plots for Metric:", metric_info$pretty_name))

                final_plot <- wrap_plots(current_metric_plot_list, ncol = min(3, length(current_metric_plot_list)))

                output_subdir <- file.path(base_output_dir, metric_info$dir_name)
                dir.create(output_subdir, showWarnings = FALSE, recursive = TRUE)

                plot_filename <- file.path(output_subdir,
                                          paste0(paste(cancer_types_subset, collapse = "_"), "_",
                                          metric_info$dir_name, "_plot.png"))

                # Save the plot
                ggsave(plot = final_plot,
                       filename = plot_filename,
                       width = 7 * min(3, length(current_metric_plot_list)),
                       height = 8.5,
                       dpi = 300,
                       device = 'png',
                       bg = "white")
            }
        }
    }

  }

  # Combine all results

  if (length(results_list_all) > 0) {
      results_all <- bind_rows(results_list_all)
      dir.create(dirname(results_file), showWarnings = FALSE, recursive = TRUE)
      write.table(results_all, file = results_file,
                  sep = "\t", row.names = FALSE, quote = FALSE)
      print(paste("Combined results saved to:", results_file))
  } else {
      warning("No results were generated to save.")
      results_all <- NULL
  }

  return(results_all)
}

# Run the function
all_results <- plot_survival_metrics(cancer_types, batch_size = 3)
