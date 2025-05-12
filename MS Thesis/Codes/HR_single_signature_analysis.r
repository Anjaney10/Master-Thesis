# Author: Anjaney J Pandey
# Code for single signature analysis and plotting forest plot for different survival metric

library(dplyr)
library(survival)
library(ggplot2)
library(rlang)
library(stringr)
library(purrr)
library(tidyr)
library(readr)
library(scales)

cancer_types <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")


base_dir <- "./Transcriptomics/TCGA_DATA"
ssgsea_base_dir <- file.path(base_dir, "ssgsea_tpm_results")
survival_base_dir <- "./Transcriptomics/tcga/survival_data"
output_base_dir <- file.path(base_dir, "Single_Signature_Survival_Analysis")
dir.create(output_base_dir, showWarnings = FALSE, recursive = TRUE)

signatures_level_cols <- c(
    "EMT_level", "Glycolysis_level", "OXPHOS_level",
    "FAO_level", "Hypoxia_level", "Ferroptosis_level"
)
signature_short_forms <- c("EMT", "Gly", "OXP", "FAO", "Hyp", "Fer")
names(signature_short_forms) <- signatures_level_cols

surv_time_col <- "OS.time"
surv_status_col <- "OS"
surv_metric_name <- "OS"
p_value_filter <- 0.05
hr_filter_upper <- 5

get_significance_stars <- function(p_values) {
  case_when( is.na(p_values) ~ "", p_values < 0.001 ~ "***", p_values < 0.01 ~ "**", p_values < 0.05 ~ "*", TRUE ~ "" ) }

# Function for Pan-Cancer Forest Plot for Single Signatures
create_pan_cancer_single_sig_plot <- function(plot_data, plot_title) {
    plot_data <- plot_data %>% mutate(cancer_type = factor(cancer_type))

    if(nrow(plot_data) == 0) { warning("No data for plot: ", plot_title); return(NULL) }
    max_ci_high <- max(plot_data$conf.high, na.rm = TRUE)
    min_ci_low <- min(plot_data$conf.low, na.rm = TRUE)
    x_limits = c(max(0, min_ci_low - 0.2), max_ci_high + 0.2)
    x_axis_expansion <- expansion(mult = c(0.05, 0.15))

    p_pan <- ggplot(plot_data, aes(x = HR, y = cancer_type)) +
        geom_vline(xintercept = 1, linetype = "dashed", color = "gray60", linewidth=0.5) +
        geom_pointrange(
            aes(xmin = conf.low, xmax = conf.high),
            color = "black",
            size = 0.5,
            linewidth = 0.7
        ) +
        geom_text(
            aes(label = p_label_text, x = conf.high),
            hjust = -0.1,
            vjust = 0.5,
            size = 2.5
        ) +
        scale_x_continuous(expand = x_axis_expansion) +
        theme_bw(base_size = 10) +
        theme(
            plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8),
            axis.title = element_text(size = 10),
            legend.position = "none",
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(),
            plot.margin = margin(t = 10, r = 15, b = 10, l = 10, unit = "pt")
        ) +
        labs(
            title = plot_title,
            x = "Hazard Ratio (HR) [+ vs -]",
            y = "Cancer Type"
        )
    return(p_pan)
}


all_single_sig_results <- list()

for (cancer in cancer_types) {
  print(paste("Processing Cancer:", cancer))

  ssgsea_file <- file.path(ssgsea_base_dir, paste0("TCGA-", cancer, "_ssgsea_results.tsv"))
  survival_file <- file.path(survival_base_dir, paste0(cancer, "_survival.txt"))
  if (!file.exists(ssgsea_file)) { warning(paste("ssGSEA not found:", ssgsea_file)); next }
  if (!file.exists(survival_file)) { warning(paste("Survival not found:", survival_file)); next }
  data_merged <- tryCatch({
     ssgsea_data <- read_tsv(ssgsea_file, show_col_types = FALSE); sample_id_col <- colnames(ssgsea_data)[1]
     colnames(ssgsea_data)[1] <- "sample_id_ssgsea"
     ssgsea_data <- ssgsea_data %>% select(sample_id_ssgsea, all_of(signatures_level_cols)) %>%
        mutate(across(all_of(signatures_level_cols), as.character))
     survival_data <- read_tsv(survival_file, show_col_types = FALSE)
     survival_data$sample <- str_sub(survival_data$sample, 1, 15)
     ssgsea_data$sample_id_ssgsea <- str_sub(ssgsea_data$sample_id_ssgsea, 1, 15)
     survival_data <- survival_data %>% distinct(sample, .keep_all = TRUE); ssgsea_data <- ssgsea_data %>% distinct(sample_id_ssgsea, .keep_all = TRUE)
     merged <- inner_join(survival_data, ssgsea_data, by = c("sample" = "sample_id_ssgsea"))
     merged <- merged %>% filter(!is.na(!!sym(surv_time_col)) & !is.na(!!sym(surv_status_col)) & !!sym(surv_time_col) > 0)
     if(nrow(merged) == 0) { stop("No rows after merge/filter.") }; merged
  }, error = function(e) { warning(paste("Error loading/merging", cancer,":", conditionMessage(e))); return(NULL) })
  if (is.null(data_merged)) next

  for (current_sig_col in signatures_level_cols) {
      current_sig_short <- signature_short_forms[current_sig_col]
      print(paste(".... Signature:", current_sig_short))

      results_sig <- tryCatch({
          if (!current_sig_col %in% colnames(data_merged)) {
              stop("Signature column '", current_sig_col, "' not found in merged data.")
          }
          sig_levels <- unique(data_merged[[current_sig_col]])
          if (length(sig_levels) < 2 || !all(c(paste0(current_sig_short,"+"), paste0(current_sig_short,"-")) %in% sig_levels)) {
               stop("Signature '", current_sig_short, "' does not have both '+' and '-' levels present.")
          }
          data_model <- data_merged %>%
              select(all_of(c(surv_time_col, surv_status_col, current_sig_col))) %>%
              mutate(!!current_sig_col := factor(!!sym(current_sig_col),
                                                 levels = c(paste0(current_sig_short,"-"), paste0(current_sig_short,"+"))))

          surv_obj_str <- paste0("Surv(", surv_time_col, ", ", surv_status_col, ")")
          formula_str <- paste0(surv_obj_str, " ~ `", current_sig_col, "`")

          # Fit model
          mod <- coxph(as.formula(formula_str), data = data_model)
          s <- summary(mod)
          if (is.null(s$coefficients) || nrow(s$coefficients) == 0) { stop("Cox model summary empty.") }

          term_name <- rownames(s$coefficients)[1]

          results_df <- data.frame(
              signature = current_sig_short,
              cancer_type = cancer,
              comparison = paste0(current_sig_short, "+ vs ", current_sig_short, "-"),
              term = term_name,
              HR = signif(s$coefficients[1, "exp(coef)"], 3),
              conf.low = signif(s$conf.int[1, "lower .95"], 3),
              conf.high = signif(s$conf.int[1, "upper .95"], 3),
              p.value = signif(s$coefficients[1, "Pr(>|z|)"], 3),
              n = s$n,
              stringsAsFactors = FALSE
          )
          results_df <- results_df %>%
              mutate(
                  significance = get_significance_stars(p.value),
                  p_label_text = paste0("p=", scales::pvalue(p.value, accuracy = 0.001, add_p = FALSE), significance)
              )
          results_df

      }, error = function(e) {
          warning(paste("Error processing signature", current_sig_short, "for cancer", cancer, ":", conditionMessage(e), ". Skipping signature."))
          return(NULL)
      })
      if (!is.null(results_sig) && nrow(results_sig) > 0) {
          all_single_sig_results[[length(all_single_sig_results) + 1]] <- results_sig
      }

  }
}

print("--- Combining and Saving All Single Signature Results ---")

if (length(all_single_sig_results) > 0) {
    combined_single_sig_df <- bind_rows(all_single_sig_results)

    # Save combined TSV
    combined_tsv_filename <- file.path(output_base_dir, paste0("PanCancer_Single_Signature_", surv_metric_name, "_HR.tsv"))
    tryCatch({
        write.table(combined_single_sig_df, file = combined_tsv_filename, sep = "\t", row.names = FALSE, quote = FALSE)
        print(paste("Combined single signature results saved:", basename(combined_tsv_filename)))
    }, error = function(e) { warning(paste("Failed to save combined TSV:", conditionMessage(e))) })


    print("--- Generating Pan-Cancer Forest Plots for Each Signature ---")
    unique_signatures <- unique(combined_single_sig_df$signature)

    for (current_sig_short in unique_signatures) {
        print(paste("...... Generating plot for signature:", current_sig_short))

        plot_data_sig <- combined_single_sig_df %>%
            filter(signature == current_sig_short, HR <= hr_filter_upper)

        if (nrow(plot_data_sig) == 0) {
            warning("No data remaining for signature '", current_sig_short, "' after filtering HR <= ", hr_filter_upper)
            next
        }
        plot_title_pan <- paste0("Pan-Cancer HR: ", current_sig_short, "+ vs ", current_sig_short, "- (", surv_metric_name, ")")
        pan_cancer_plot <- create_pan_cancer_single_sig_plot(plot_data_sig, plot_title_pan)

        # Save the plot
        if (!is.null(pan_cancer_plot)) {
            plot_filename <- file.path(output_base_dir, paste0("PanCancer_", current_sig_short, "_", surv_metric_name, "_Plot_Filtered.png"))
            n_cancers_in_plot <- n_distinct(plot_data_sig$cancer_type)
            plot_height <- max(5, 2 + n_cancers_in_plot * 0.25)
             tryCatch({
                ggsave(filename = plot_filename, plot = pan_cancer_plot, width = 7, height = plot_height, dpi = 300, bg = "white", limitsize = FALSE)
                print(paste("......... Plot saved:", basename(plot_filename)))
             }, error = function(e) { warning(paste("Failed to save plot:", conditionMessage(e))) })
        } else {
            warning("Plot creation failed for signature: ", current_sig_short)
        }
    }

} else {
    print("No single signature results were generated.")
}

print("--- Processing Complete ---")
