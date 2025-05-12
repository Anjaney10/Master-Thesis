# Author: Anjaney J Pandey
# Code for subtype wise KM plot for TCGA BRCA data

####################### SUBTYPE WISE KM PLOT FOR TCGA DATA #######################

library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(stringr)
library(patchwork)
library(grid)

status_col <- "OS"
time_col <- "OS.time"
subtype_col <- "PAM50Call_RNAseq"
subtypes_to_plot <- c("LumA", "LumB", "Her2", "Basal")
signature_bases <- c("EMT", "Glycolysis", "OXPHOS", "FAO", "Hypoxia", "Ferroptosis", "RPMS", "BPMS")
signature_level_cols <- paste0(signature_bases, "_level")
level_colors <- c("red", "blue")
legend_labels_generic <- c("High Level", "Low Level")
output_dir <- "KM_Plots_Subtype_HR_Pval_FINAL"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

format_pval <- function(pv) {
  if (is.na(pv)) return("NA")
  if (pv < 0.001) {
    return("p < 0.001")
  } else {
    return(paste("p =", format(round(pv, 3), nsmall = 3)))
  }
}

tcga_final_df_filtered <- tcga_final_df %>%
  filter(.data[[subtype_col]] %in% subtypes_to_plot) %>%
  mutate(!!subtype_col := factor(.data[[subtype_col]], levels = subtypes_to_plot)) %>%
  filter(!is.na(.data[[time_col]]), !is.na(.data[[status_col]]))

if (nrow(tcga_final_df_filtered) == 0) {
  stop("No data remaining after filtering for specified subtypes.")
}

for (i in seq_along(signature_bases)) {
  sig_base <- signature_bases[i]
  sig_level_col <- signature_level_cols[i]

  cat("\nProcessing Signature:", sig_base, "\n")

  # Data Prep
  if (!sig_level_col %in% colnames(tcga_final_df_filtered)) { cat(" Column not found. Skipping.\n"); next }
  df_signature_filtered <- tcga_final_df_filtered %>%
      filter(!is.na(.data[[sig_level_col]]) & .data[[sig_level_col]] != "")
  if (nrow(df_signature_filtered) < 2) { cat(" Insufficient data after NA filter. Skipping.\n"); next }
  actual_levels <- unique(df_signature_filtered[[sig_level_col]])
  if (length(actual_levels) != 2) { cat(" Not exactly 2 levels found. Skipping.\n"); next }
  positive_level <- actual_levels[str_detect(actual_levels, "\\+")]
  negative_level <- actual_levels[str_detect(actual_levels, "\\-")]
  if (length(positive_level) != 1 || length(negative_level) != 1) { cat(" Could not ID +/- levels. Skipping.\n"); next }

  df_signature_filtered <- df_signature_filtered %>%
      mutate(!!sig_level_col := factor(.data[[sig_level_col]], levels = c(positive_level, negative_level)))

  plot_list <- list()
  risk_table_list <- list()
  is_first_plot <- TRUE
  common_legend <- NULL

  for (subtype in subtypes_to_plot) {
    cat("  Processing subtype:", subtype, "\n")
    df_subtype_signature <- df_signature_filtered %>%
      filter(.data[[subtype_col]] == subtype)

    n_levels_subtype <- length(unique(df_subtype_signature[[sig_level_col]]))
    if(nrow(df_subtype_signature) < 2 || n_levels_subtype < 2) {
        cat("    Skipping plot for subtype", subtype, "- insufficient data or levels.\n")
        plot_list[[subtype]] <- ggplot() + theme_void() + ggtitle(paste(subtype, "\n(No Data)")) + theme(plot.title = element_text(hjust=0.5, size=10))
        risk_table_list[[subtype]] <- NULL
        next
    }

    surv_obj_subtype <- Surv(time = df_subtype_signature[[time_col]], event = df_subtype_signature[[status_col]])

    # Calculate HR and P-values
    hr_text <- "HR = NA"
    logrank_p_text <- "Log-rank p = NA"
    stats_text <- "HR = NA\nLog-rank p = NA"

    df_cox <- df_subtype_signature
    df_cox[[sig_level_col]] <- factor(df_cox[[sig_level_col]], levels = c(negative_level, positive_level))

    # Fit Cox model
    cox_model <- tryCatch({
        coxph(surv_obj_subtype ~ get(sig_level_col), data = df_cox)
    }, error = function(e) {
        cat("    Cox model failed for", subtype, ":", conditionMessage(e), "\n")
        NULL
    })

    if (!is.null(cox_model)) {
        cox_summary <- summary(cox_model)
        hr <- round(cox_summary$conf.int[1, "exp(coef)"], 2)
        hr_lci <- round(cox_summary$conf.int[1, "lower .95"], 2)
        hr_uci <- round(cox_summary$conf.int[1, "upper .95"], 2)
        cox_p <- cox_summary$logtest["pvalue"]
        hr_text <- sprintf("HR = %.2f (%.2f-%.2f)", hr, hr_lci, hr_uci)
    }
    logrank_test <- tryCatch({
       survdiff(surv_obj_subtype ~ get(sig_level_col), data = df_subtype_signature)
    }, error = function(e) {
       cat("    survdiff failed for", subtype, ":", conditionMessage(e), "\n")
       NULL
    })

    if(!is.null(logrank_test)){
        logrank_p_val <- 1 - pchisq(logrank_test$chisq, length(logrank_test$n) - 1)
        logrank_p_text <- format_pval(logrank_p_val)
    }

    # Combine HR and Log-rank p for annotation
    stats_text <- paste(hr_text, logrank_p_text, sep = "\n")

    # Fit for KM plot
    fit_subtype <- survfit(surv_obj_subtype ~ get(sig_level_col), data = df_subtype_signature)

    # Create the plot
    ggsurv_obj <- tryCatch({
        ggsurvplot(
            fit_subtype,
            data = df_subtype_signature,
            title = subtype,
            pval=FALSE,
            conf.int = FALSE,
            risk.table = TRUE,
            risk.table.col = "strata",
            risk.table.y.text = is_first_plot,
            risk.table.y.text.col = TRUE,
            tables.y.lab = "Number at risk",
            palette = level_colors,
            legend = "none",
            xlab = "Time (Days)",
            ylab = "Survival Probability",
            ggtheme = theme_minimal(base_size = 11) +
                      theme(
                         plot.title = element_text(hjust = 0.5, face="bold", size=rel(1.3)),
                         axis.text = element_text(face="bold", size=rel(1.1)),
                         axis.title.x = element_text(face="bold", size=rel(1.1)),
                         axis.title.y = element_text(face="bold", size=rel(0.7))
                      )
        )
      }, error = function(e) {
          cat("    Error generating ggsurvplot for subtype", subtype, ":", conditionMessage(e), "\n")
          return(NULL)
      }
    )

    if (!is.null(ggsurv_obj)) {
         max_time <- max(df_subtype_signature[[time_col]], na.rm=TRUE) * 0.95
         y_pos = 0.95

         plot_annotated <- ggsurv_obj$plot +
               annotate("text", x = max_time, y = y_pos, label = stats_text,
                        hjust = 1, vjust = 1, size = 5, fontface="bold")
         plot_list[[subtype]] <- plot_annotated
         risk_table_list[[subtype]] <- ggsurv_obj$table +
                                          theme(axis.title.y = element_blank(),
                                                axis.text.y = element_blank(),
                                                axis.ticks.y = element_blank(),
                                                plot.title = element_blank(),
                                                legend.position = "none")
         if (is.null(common_legend)) {
             cat("    Attempting to extract legend from:", subtype, "\n")
             ggsurv_temp <- ggsurvplot(fit_subtype, data = df_subtype_signature,
                                         palette = level_colors,
                                         legend.title = sig_base,
                                         legend.labs = legend_labels_generic,
                                         font.legend = list(size=17, face="bold"))
             common_legend <- get_legend(ggsurv_temp)
             if (!is.null(common_legend)) cat("    Legend extracted.\n") else cat("    Failed legend extraction.\n")
         }
         is_first_plot <- FALSE
    } else {
         plot_list[[subtype]] <- ggplot() + theme_void() + ggtitle(paste(subtype, "\n(Plot Error)")) + theme(plot.title = element_text(hjust=0.5, size=10))
         risk_table_list[[subtype]] <- NULL
    }

  }
  successful_plots <- plot_list[!sapply(plot_list, is.null)]
  if(length(successful_plots) > 0) {

      plots_to_arrange <- plot_list[subtypes_to_plot]
      plot_layout <- wrap_plots(plots_to_arrange, ncol = 2)
      table_placeholders <- lapply(subtypes_to_plot, function(st) if(is.null(risk_table_list[[st]])) {plot_spacer()} else {risk_table_list[[st]]})
      table_layout <- wrap_plots(table_placeholders, ncol = 2)
      combined_plot_table <- plot_layout / table_layout +
                                plot_layout(heights = c(3.5, 1))

      final_plot_intermediate <- combined_plot_table +
          plot_annotation(title = paste("Overall Survival for", sig_base)) +
          plot_layout(axes = 'collect') &
          theme(plot.title = element_text(hjust = 0.5, face="bold", size=16))

     final_plot_labeled <- final_plot_intermediate +
                    labs(x = "Time (Days)", y = "Survival Probability") +
                    theme(axis.title = element_text(size=12, face="bold"))

      final_plot_with_legend <- final_plot_labeled
      if (!is.null(common_legend) && inherits(common_legend, "grob")) {
            legend_patch <- wrap_elements(common_legend) + theme(plot.background = element_blank())
            combined_layout <- legend_patch / final_plot_labeled
            final_plot_with_legend <- combined_layout +
                                       plot_layout(heights = c(1, 15))

      } else {
          cat("  Warning: Common legend is NULL or not a grob. Skipping legend addition.\n")
      }

      filename <- file.path(output_dir, paste0("KM_Plot_", sig_base, "_Final_HR_Pval.png"))
      ggsave(filename, plot = final_plot_with_legend, width = 13, height = 11, dpi = 300, units = "in")
      cat("  Saved combined plot:", filename, "\n")

  } else {
      cat("  Skipping combination/saving for", sig_base, "- no successful subplots generated.\n")
  }

}
