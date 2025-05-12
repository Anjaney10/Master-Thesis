# Author: Anjaney J Pandey
# FINAL CODE FOR PAN CANCER SURVIVAL METRICS - Forest Plot

library(dplyr)
library(survival)
library(ggplot2)
library(rlang)
library(stringr)
library(purrr)
library(tidyr)
library(scales)

cancer_types <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
base_dir <- "./Transcriptomics/TCGA_DATA/MR_Combination_Survival_Analysis"
ssgsea_base_dir <- file.path(base_dir, "ssgsea_tpm_results")
survival_base_dir <- "./Transcriptomics/tcga/survival_data"
output_base_dir <- file.path(base_dir, "Combinatorial_Survival_Analysis_New")

signatures_of_interest <- c("EMT_level", "Glycolysis_level", "OXPHOS_level", "FAO_level", "Hypoxia_level", "Ferroptosis_level")
signature_info <- list(
  EMT_level = list(type = "pro", normal_level = "EMT-", short_form = "EMT"),
  Glycolysis_level = list(type = "pro", normal_level = "Gly-", short_form = "Gly"),
  OXPHOS_level = list(type = "pro", normal_level = "OXP-", short_form = "OXP"),
  FAO_level = list(type = "pro", normal_level = "FAO-", short_form = "FAO"),
  Hypoxia_level = list(type = "pro", normal_level = "Hyp-", short_form = "Hyp"),
  Ferroptosis_level = list(type = "pro", normal_level = "Fer-", short_form = "Fer")
)
surv_time_col <- "OS.time"
surv_status_col <- "OS"
surv_metric_name <- "OS"

# Helper Functions

get_significance_stars <- function(p_values) {
  case_when( is.na(p_values) ~ "", p_values < 0.001 ~ "***", p_values < 0.01 ~ "**", p_values < 0.05 ~ "*", TRUE ~ "" ) }

get_reference_level <- function(combo_signatures_cols, sig_info) {
  ref_parts <- sapply(combo_signatures_cols, function(sig_col_name) {
    info <- sig_info[[sig_col_name]]; if (is.null(info)) { stop("Info not found: ", sig_col_name) }
    sig_prefix <- info$short_form; normal_suffix <- ifelse(info$type == "pro", "-", "+")
    paste0(sig_prefix, normal_suffix)
  }); paste(ref_parts, collapse = "/") }

create_forest_plot <- function(results_df, plot_title) {
    results_df <- results_df %>% mutate(
            cleaned_label_temp = gsub("^Combination_Level", "", term),
            plot_labels = factor(cleaned_label_temp, levels = rev(unique(cleaned_label_temp))) )
    results_df <- results_df %>% filter(!is.na(plot_labels))
    if(nrow(results_df) == 0 || nlevels(results_df$plot_labels) == 0) { warning("Plotting failed: No valid levels."); return(NULL) }
    y_axis_expansion <- expansion(mult = c(0.05, 0.20))
    p <- ggplot(results_df, aes(x = plot_labels, y = HR)) +
        geom_hline(yintercept = 1, linetype = "dashed", color = "gray60", linewidth=0.5) +
        geom_pointrange(aes(ymin = conf.low, ymax = conf.high), color = "black", size = 0.6, linewidth = 0.8) +
        geom_text(aes(label = p_label_text, y = conf.high), hjust = -0.1, vjust = 0.5, size = 2.8) +
        scale_y_continuous(expand = y_axis_expansion) +
        coord_flip() + theme_bw(base_size = 11) +
        theme( legend.position = "none", plot.title = element_text(size = 12, face = "bold", hjust = 0.5, lineheight=1.1),
            axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
            axis.title.x = element_text(size = 11), axis.title.y = element_blank(),
            plot.margin = margin(t = 10, r = 15, b = 10, l = 10, unit = "pt") ) +
        labs( title = plot_title, x = NULL, y = "Hazard Ratio (HR)" )
    return(p) }


# Main Loop

pan_cancer_results_by_size <- list(
  "2" = list(),
  "3" = list()
)

for (cancer in cancer_types) {
  print(paste("Processing Cancer:", cancer))
  ssgsea_file <- file.path(ssgsea_base_dir, paste0("TCGA-", cancer, "_ssgsea_results.tsv"))
  survival_file <- file.path(survival_base_dir, paste0(cancer, "_survival.txt"))
  if (!file.exists(ssgsea_file)) { warning(paste("ssGSEA not found:", ssgsea_file)); next }
  if (!file.exists(survival_file)) { warning(paste("Survival not found:", survival_file)); next }
  data_merged <- tryCatch({
     ssgsea_data <- read.delim(ssgsea_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
     sample_id_col <- colnames(ssgsea_data)[1]; colnames(ssgsea_data)[1] <- "sample_id_ssgsea"
     ssgsea_data <- ssgsea_data %>% select(sample_id_ssgsea, all_of(signatures_of_interest)) %>%
        mutate(across(all_of(signatures_of_interest), as.character))
     survival_data <- read.delim(survival_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
     survival_data$sample <- str_sub(survival_data$sample, 1, 15)
     ssgsea_data$sample_id_ssgsea <- str_sub(ssgsea_data$sample_id_ssgsea, 1, 15)
     survival_data <- survival_data %>% distinct(sample, .keep_all = TRUE)
     ssgsea_data <- ssgsea_data %>% distinct(sample_id_ssgsea, .keep_all = TRUE)
     merged <- inner_join(survival_data, ssgsea_data, by = c("sample" = "sample_id_ssgsea"))
     merged <- merged %>% filter(!is.na(!!sym(surv_time_col)) & !is.na(!!sym(surv_status_col)) & !!sym(surv_time_col) > 0)
     if(nrow(merged) == 0) { stop("No rows after merge/filter.") }; merged
  }, error = function(e) { warning(paste("Error loading/merging", cancer,":", conditionMessage(e))); return(NULL) })
  if (is.null(data_merged)) next

#############################################################
  for (combo_size in 2:3) {
    print(paste(".. Processing combinations of size:", combo_size))

    combo_size_name <- case_when(
        combo_size == 2 ~ "Two",
        combo_size == 3 ~ "Three",
        TRUE ~ paste0("Size_", combo_size)
    )
    output_subdir_name <- paste0(combo_size_name, "_Signature_Combination")
    output_dir_combo <- file.path(output_base_dir, output_subdir_name)
    dir.create(output_dir_combo, showWarnings = FALSE, recursive = TRUE)

    # Generate combinations
    signature_combinations <- combn(signatures_of_interest, combo_size, simplify = FALSE)

###############################################################
    for (current_combo_signatures in signature_combinations) {
      combo_prefix <- paste(sapply(current_combo_signatures, function(x) signature_info[[x]]$short_form), collapse = "_")
      print(paste(".... Combination:", combo_prefix))
      results_combo <- tryCatch({
          interaction_col_name <- "Combination_Level"
          data_model <- data_merged %>%
              select(all_of(c(surv_time_col, surv_status_col, current_combo_signatures))) %>%
              unite(!!interaction_col_name, all_of(current_combo_signatures), sep = "/", remove = FALSE, na.rm = TRUE) %>%
              mutate(!!interaction_col_name := factor(!!sym(interaction_col_name)))

          if (nlevels(data_model[[interaction_col_name]]) < 2) { stop("Less than 2 combination levels found.") }

          reference_level <- get_reference_level(current_combo_signatures, signature_info)
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
          results_df <- results_df %>% mutate(
                  significance = get_significance_stars(p.value),
                  p_label_text = paste0("p=", scales::pvalue(p.value, accuracy = 0.001, add_p = FALSE), significance) )
          results_df$cancer_type <- cancer; results_df$combination <- combo_prefix;
          results_df$reference_level <- reference_level; results_df$survival_metric <- surv_metric_name
          results_df
        }, error = function(e) {
          warning(paste("Error processing combo", combo_prefix, "for", cancer, ":", conditionMessage(e))); return(NULL) })

        # Save Results and Plot
        if (!is.null(results_combo) && nrow(results_combo) > 0) {
        results_filename <- file.path(output_dir_combo, paste0(cancer, "_", combo_prefix, "_", surv_metric_name, "_HR.tsv"))
        tryCatch({ write.table(results_combo, file = results_filename, sep = "\t", row.names = FALSE, quote = FALSE)
          }, error = function(e) { warning(paste("Failed save results:", conditionMessage(e))) })

        plot_logrank_p <- results_combo$Log_Rank_pvalue[1]
        plot_title <- paste0(cancer, ": ", combo_prefix, " (vs ", results_combo$reference_level[1], ")\n", surv_metric_name, ", LogRank p = ", signif(plot_logrank_p, 2))
        forest_plot <- create_forest_plot(results_combo, plot_title)
        if (!is.null(forest_plot)) {
            plot_filename <- file.path(output_dir_combo, paste0(cancer, "_", combo_prefix, "_", surv_metric_name, "_Plot.png"))
            plot_dynamic_height = max(4, 2 + nrow(results_combo)*0.35)
             tryCatch({ ggsave(filename = plot_filename, plot = forest_plot, width = 7, height = plot_dynamic_height, dpi = 300, bg = "white", limitsize = FALSE)
             }, error = function(e) { warning(paste("Failed save plot:", conditionMessage(e))) })
        } else { warning("Plot creation failed for ", combo_prefix, " in ", cancer) }

        if (combo_size %in% c(2, 3)) {
             size_char <- as.character(combo_size)
             pan_cancer_results_by_size[[size_char]] <- c(pan_cancer_results_by_size[[size_char]], list(results_combo))
        }

      }
   }
  }
}


# Pan-Cancer Analysis Starts Here for 2 and 3 combination

print("Starting Pan-Cancer Analysis for Combinations of Size 2 & 3...")

pan_cancer_output_dir <- file.path(output_base_dir, "pan_cancer_combinations")
dir.create(pan_cancer_output_dir, showWarnings = FALSE, recursive = TRUE)

create_pan_cancer_plot <- function(plot_data, plot_title, conf_high_limit = 7) {
    plot_data <- plot_data %>%
        filter(is.finite(HR) & is.finite(conf.low) & is.finite(conf.high) & conf.high <= conf_high_limit) %>%
        mutate(
            plot_labels = gsub("^Combination_Level", "", term),
            plot_labels = ifelse(plot_labels == "", term, plot_labels),
            plot_labels = factor(plot_labels),
            cancer_type = factor(cancer_type)
        ) %>%
        filter(!is.na(plot_labels) & plot_labels != "")

    if(nrow(plot_data) == 0) {
        warning("No data remaining after filtering (conf.high <= ", conf_high_limit, ", finite HR/CI, non-empty labels) for pan-cancer plot: ", plot_title)
        return(NULL)
    }
    max_hr_filtered <- max(plot_data$conf.high, na.rm = TRUE)
    min_hr_filtered <- min(plot_data$conf.low, na.rm = TRUE)
    x_axis_expansion <- expansion(mult = c(0.1, 0.15))
    n_colors_needed <- n_distinct(plot_data$plot_labels)
    color_palette_name <- "Set1"

    # Create the plot
    p_pan <- ggplot(plot_data, aes(x = HR, y = cancer_type)) +
        geom_vline(xintercept = 1, linetype = "dashed", color = "gray60", linewidth=0.5) +
        geom_pointrange(
            aes(xmin = conf.low, xmax = conf.high, color = plot_labels),
            size = 0.5, linewidth = 0.7, position = position_dodge(width = 0.6)
        ) +
        geom_text(
            aes(label = p_label_text, x = conf.high, color = plot_labels),
            hjust = -0.1, vjust = 0.5, size = 2.2,
            position = position_dodge(width = 0.6), show.legend = FALSE
        ) +
        scale_x_continuous(expand = x_axis_expansion) +
        scale_color_brewer(palette = color_palette_name, name = "Comparison vs Ref.") +
        theme_bw(base_size = 10) +
        theme(
            plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8),
            axis.title = element_text(size = 10),
            legend.position = "right",
            legend.title = element_text(size=9),
            legend.text = element_text(size=8),
            plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")
        ) +
        labs(
            title = plot_title,
            x = "Hazard Ratio (HR)",
            y = "Cancer Type"
        )
    return(p_pan)
}
# Process collected results
for (combo_size_char in c("2", "3")) {
    combo_size_num <- as.numeric(combo_size_char)
    print(paste(".. Processing Pan-Cancer results for combinations of size:", combo_size_num))

    results_list_current_size <- pan_cancer_results_by_size[[combo_size_char]]
    if (length(results_list_current_size) == 0) {
        print(paste(".... No results collected for size", combo_size_num))
        next
    }
    combined_df_current_size <- bind_rows(results_list_current_size)

    # Save combined TSV
    combo_size_name <- case_when(combo_size_num == 2 ~ "Two", combo_size_num == 3 ~ "Three", TRUE~"Other")
    combined_tsv_filename <- file.path(pan_cancer_output_dir, paste0("PanCancer_", combo_size_name, "_Signature_Combinations_HR.tsv"))
    tryCatch({
        write.table(combined_df_current_size, file = combined_tsv_filename, sep = "\t", row.names = FALSE, quote = FALSE)
        print(paste(".... Combined TSV saved:", basename(combined_tsv_filename)))
    }, error = function(e) { warning(paste("Failed to save combined TSV for size", combo_size_num, ":", conditionMessage(e))) })

    # Generate Pan-Cancer Plot
    unique_combinations <- unique(combined_df_current_size$combination)

    for (current_combo_name in unique_combinations) {
        print(paste("...... Generating Pan-Cancer plot for combination:", current_combo_name))

        # Filter data
        plot_data_combo <- combined_df_current_size %>%
            filter(combination == current_combo_name)
        plot_title_pan <- paste0("Pan-Cancer Analysis: ", current_combo_name, "\n", surv_metric_name)
        pan_cancer_plot <- create_pan_cancer_plot(plot_data_combo, plot_title_pan)

        # Save the plot
        if (!is.null(pan_cancer_plot)) {
            pan_plot_filename <- file.path(pan_cancer_output_dir, paste0("PanCancer_", current_combo_name, "_", surv_metric_name, "_Plot.png"))
            n_cancers_in_plot <- n_distinct(plot_data_combo$cancer_type)
            pan_plot_height <- max(5, 2 + n_cancers_in_plot * 0.3)
             tryCatch({
                ggsave(filename = pan_plot_filename, plot = pan_cancer_plot, width = 9, height = pan_plot_height, dpi = 300, bg = "white", limitsize = FALSE)
                print(paste("......... Pan-Cancer Plot saved:", basename(pan_plot_filename)))
             }, error = function(e) { warning(paste("Failed to save pan-cancer plot:", conditionMessage(e))) })
        } else {
            warning("Pan-cancer plot creation failed for combination: ", current_combo_name)
        }
    }
}
