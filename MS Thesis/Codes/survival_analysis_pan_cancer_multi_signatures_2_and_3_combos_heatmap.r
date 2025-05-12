# Author: Anjaney J Pandey
# FINAL CODE FOR Metabolic Reprogramming Signatures and PAN CANCER SURVIVAL METRICS
# Here we are going to generate 2 and 3 signature combinations and calculate the hazard ratios for the combinations
# We will then plot the heatmap for such combinations

library(dplyr)
library(survival)
library(ggplot2)
library(rlang)
library(stringr)
library(tidyr)
library(scales)
library(readr)

cancer_types <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
base_dir <- "./Transcriptomics/TCGA_DATA"
ssgsea_base_dir <- file.path("./Transcriptomics/TCGA_DATA/MR_Combination_Survival_Analysis/ssgsea_tpm_results")
survival_base_dir <- "./Transcriptomics/tcga/survival_data"

output_base_dir <- file.path(base_dir, "Combinatorial_Survival_AllNegRef")
pan_cancer_output_dir <- file.path(output_base_dir, "pan_cancer_combinations")

signatures_of_interest <- c("EMT_level", "Glycolysis_level", "OXPHOS_level", "FAO_level", "Hypoxia_level", "Ferroptosis_level")

signature_info <- list(
  EMT_level = list(short_form = "EMT"),
  Glycolysis_level = list(short_form = "Gly"),
  OXPHOS_level = list(short_form = "OXP"),
  FAO_level = list(short_form = "FAO"),
  Hypoxia_level = list(short_form = "Hyp"),
  Ferroptosis_level = list(short_form = "Fer")
)
surv_time_col <- "OS.time"
surv_status_col <- "OS"
surv_metric_name <- "OS"

# Significance Stars
get_significance_stars <- function(p_values) {
  case_when( is.na(p_values) ~ "", p_values < 0.001 ~ "***", p_values < 0.01 ~ "**", p_values < 0.05 ~ "*", TRUE ~ "" ) }

get_all_negative_reference <- function(combo_signatures_cols, sig_info) {
    ref_parts <- sapply(combo_signatures_cols, function(sig_col_name) {
        info <- sig_info[[sig_col_name]]
        if (is.null(info)) { stop("Info not found: ", sig_col_name) }
        sig_prefix <- info$short_form
        paste0(sig_prefix, "-")
    })
    paste(ref_parts, collapse = "/")
}

#################################################################################################
create_pan_cancer_heatmap_all_neg_ref <- function(plot_data, plot_title, conf_high_limit) {
    plot_data_filtered <- plot_data %>%
        filter(is.finite(HR) & is.finite(conf.high) & conf.high <= conf_high_limit)

    if(nrow(plot_data_filtered) == 0) {
        warning("No data remaining after filtering (conf.high <= ", conf_high_limit, ") for heatmap: ", plot_title)
        return(NULL)
    }
    plot_data_ready <- plot_data_filtered %>%
        mutate(
            plot_labels = gsub("^Combination_Level", "", term),
            plot_labels = ifelse(plot_labels == "", term, plot_labels),
            plot_labels = factor(plot_labels, levels = rev(sort(unique(plot_labels)))),
            cancer_type = factor(cancer_type),
            significance = get_significance_stars(p.value)
        ) %>%
        filter(!is.na(plot_labels) & plot_labels != "")

     if(nrow(plot_data_ready) == 0) {
        warning("No data remaining for '", plot_title, "' after cleaning labels. Skipping heatmap.")
        return(NULL)
    }
    max_abs_dev <- max(abs(plot_data_ready$HR - 1), 0.1, na.rm = TRUE)
    color_limits <- c(max(0, 1 - max_abs_dev), 1 + max_abs_dev)

    # Create Heatmap
    heatmap_plot <- ggplot(plot_data_ready, aes(x = cancer_type, y = plot_labels, fill = HR)) +
      geom_tile(color = "grey90", size=0.2) +
      geom_text(aes(label = significance), color = "black", size = 3.5) +
      scale_fill_gradient2(
            name = "Hazard Ratio (HR)", low = muted("blue"), mid = "white", high = muted("red"),
            midpoint = 1, limit = color_limits, oob = scales::squish ) +
      theme_minimal(base_size = 10) +
      theme(
        axis.text.x = element_text(angle = 60, hjust = 1, size = 9), axis.text.y = element_text(size = 9),
        axis.title = element_blank(), panel.grid = element_blank(),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5), legend.position = "right",
        legend.title = element_text(size=9), legend.text = element_text(size=8),
        legend.key.height = unit(1.2, "cm") ) +
      labs( title = plot_title )

    return(heatmap_plot)
}


###########################################################
###########################################################

dir.create(output_base_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(pan_cancer_output_dir, showWarnings = FALSE, recursive = TRUE)

pan_cancer_results_by_size <- list("2" = list(), "3" = list())

for (cancer in cancer_types) {
  print(paste("Processing Cancer:", cancer))
  ssgsea_file <- file.path(ssgsea_base_dir, paste0("TCGA-", cancer, "_ssgsea_results.tsv"))
  survival_file <- file.path(survival_base_dir, paste0(cancer, "_survival.txt"))
  if (!file.exists(ssgsea_file)) { warning(paste("ssGSEA not found:", ssgsea_file)); next }
  if (!file.exists(survival_file)) { warning(paste("Survival not found:", survival_file)); next }
  data_merged <- tryCatch({
     ssgsea_data <- read_tsv(ssgsea_file, show_col_types = FALSE); sample_id_col <- colnames(ssgsea_data)[1]
     colnames(ssgsea_data)[1] <- "sample_id_ssgsea"
     ssgsea_data <- ssgsea_data %>% select(sample_id_ssgsea, all_of(signatures_of_interest)) %>%
        mutate(across(all_of(signatures_of_interest), as.character))
     survival_data <- read_tsv(survival_file, show_col_types = FALSE)
     survival_data$sample <- str_sub(survival_data$sample, 1, 15); ssgsea_data$sample_id_ssgsea <- str_sub(ssgsea_data$sample_id_ssgsea, 1, 15)
     survival_data <- survival_data %>% distinct(sample, .keep_all = TRUE); ssgsea_data <- ssgsea_data %>% distinct(sample_id_ssgsea, .keep_all = TRUE)
     merged <- inner_join(survival_data, ssgsea_data, by = c("sample" = "sample_id_ssgsea"))
     merged <- merged %>% filter(!is.na(!!sym(surv_time_col)) & !is.na(!!sym(surv_status_col)) & !!sym(surv_time_col) > 0)
     if(nrow(merged) == 0) { stop("No rows after merge/filter.") }; merged
  }, error = function(e) { warning(paste("Error loading/merging", cancer,":", conditionMessage(e))); return(NULL) })
  if (is.null(data_merged)) next

########################################################################

  for (combo_size in 2:3) {
    print(paste(".. Processing combinations of size:", combo_size))

    signature_combinations <- combn(signatures_of_interest, combo_size, simplify = FALSE)

    for (current_combo_signatures in signature_combinations) {
      combo_prefix <- paste(sapply(current_combo_signatures, function(x) signature_info[[x]]$short_form), collapse = "_")
      print(paste(".... Combination:", combo_prefix))

      results_combo <- tryCatch({
          interaction_col_name <- "Combination_Level"
          data_model <- data_merged %>%
              select(all_of(c(surv_time_col, surv_status_col, current_combo_signatures))) %>%
              unite(!!interaction_col_name, all_of(current_combo_signatures), sep = "/", remove = FALSE, na.rm = TRUE) %>%
              mutate(!!interaction_col_name := factor(!!sym(interaction_col_name)))

          if (nlevels(data_model[[interaction_col_name]]) < 2) { stop("Less than 2 levels.") }

          reference_level <- get_all_negative_reference(current_combo_signatures, signature_info)
          if (!reference_level %in% levels(data_model[[interaction_col_name]])) {
              stop(paste("Reference level '", reference_level, "' not found in data.")) }
          data_model[[interaction_col_name]] <- relevel(data_model[[interaction_col_name]], ref = reference_level)
          print(paste("...... Reference Level set to:", reference_level))

          # Run Cox Model
          surv_obj_str <- paste0("Surv(", surv_time_col, ", ", surv_status_col, ")")
          formula_str <- paste0(surv_obj_str, " ~ ", interaction_col_name)
          mod <- coxph(as.formula(formula_str), data = data_model)
          s <- summary(mod)

          # Extract Results
          if (is.null(s$coefficients) || nrow(s$coefficients) == 0) { stop("Cox model summary empty.") }
          results_df <- data.frame( term = rownames(s$coefficients), HR = signif(s$coefficients[, "exp(coef)"], 3),
              conf.low = signif(s$conf.int[, "lower .95"], 3), conf.high = signif(s$conf.int[, "upper .95"], 3),
              p.value = signif(s$coefficients[, "Pr(>|z|)"], 3), Log_Rank_pvalue = signif(s$logtest[['pvalue']], 3),
              n = s$n, stringsAsFactors = FALSE )

          # Add metadata
          results_df$cancer_type <- cancer; results_df$combination <- combo_prefix;
          results_df$reference_level <- reference_level; results_df$survival_metric <- surv_metric_name
          results_df
        }, error = function(e) { warning(paste("Error combo", combo_prefix,"for",cancer,":", conditionMessage(e))); return(NULL) })

      # Collect Results for Pan-Cancer
      if (!is.null(results_combo) && nrow(results_combo) > 0) {
          size_char <- as.character(combo_size)
          pan_cancer_results_by_size[[size_char]] <- c(pan_cancer_results_by_size[[size_char]], list(results_combo))
      }

    }
  }
}


# Save Combined Results & Generate Pan-Cancer Heatmaps
print("--- Starting Pan-Cancer Results Saving and Heatmap Generation")

for (combo_size_char in c("2", "3")) {
    combo_size_num <- as.numeric(combo_size_char)
    print(paste(".. Processing Pan-Cancer size:", combo_size_num))

    results_list_current_size <- pan_cancer_results_by_size[[combo_size_char]]
    if (length(results_list_current_size) == 0) { print(paste(".... No results for size", combo_size_num)); next }

    combined_df_current_size <- bind_rows(results_list_current_size)

    # Save combined TSV
    combo_size_name <- case_when(combo_size_num == 2 ~ "Two", combo_size_num == 3 ~ "Three", TRUE~"Other")
    combined_tsv_filename <- file.path(pan_cancer_output_dir, paste0("PanCancer_", combo_size_name, "_Signature_Combinations_AllNegRef_HR.tsv"))
    tryCatch({
        write.table(combined_df_current_size, file = combined_tsv_filename, sep = "\t", row.names = FALSE, quote = FALSE)
        print(paste(".... Combined TSV saved:", basename(combined_tsv_filename)))
    }, error = function(e) { warning(paste("Failed to save combined TSV:", conditionMessage(e))) })

    # Generate Heatmap
    unique_combinations <- unique(combined_df_current_size$combination)
    conf_filter_limit <- 7

    for (current_combo_name in unique_combinations) {
        print(paste("...... Generating heatmap for combination:", current_combo_name))
        plot_data_combo <- combined_df_current_size %>% filter(combination == current_combo_name)
        plot_title_heatmap <- paste("Pan-Cancer HR:", current_combo_name, "(Ref: All Negative)")

        # Create the heatmap
        pan_cancer_heatmap <- create_pan_cancer_heatmap_all_neg_ref(
            plot_data = plot_data_combo,
            plot_title = plot_title_heatmap,
            conf_high_limit = conf_filter_limit
            )

        # Save the heatmap
        if (!is.null(pan_cancer_heatmap)) {
            heatmap_filename <- file.path(pan_cancer_output_dir, paste0("PanCancer_", current_combo_name, "_", surv_metric_name, "_AllNegRef_Heatmap_Filtered.png"))
            n_cancers_plot <- n_distinct(plot_data_combo$cancer_type)
            n_labels_plot <- n_distinct(plot_data_combo$term)
            plot_width <- 12
            plot_height <- 5
             tryCatch({
                ggsave(filename = heatmap_filename, plot = pan_cancer_heatmap, width = plot_width, height = plot_height, dpi = 300, bg = "white", limitsize = FALSE)
                print(paste("......... Heatmap saved:", basename(heatmap_filename)))
             }, error = function(e) { warning(paste("Failed save heatmap:", conditionMessage(e))) })
        }
    }
}

print("--- Processing Complete")
