# Author: Anjaney J Pandey
# Code for TCGA BRCA subtype KM plots signature wise

################################### KM PLOT FOR TCGA BRCA SIGNATURE WISE ######################################

library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(stringr)
library(RColorBrewer)

tcga_final_df <- tcga_final_df %>%
  filter(!is.na(PAM50Call_RNAseq) &
         PAM50Call_RNAseq != "" &
         PAM50Call_RNAseq != "Normal" &
         !is.na(OS) &
         !is.na(OS.time))

status_col <- "OS"
time_col <- "OS.time"
subtype_col <- "PAM50Call_RNAseq"

signature_bases <- c("EMT", "Glycolysis", "OXPHOS", "FAO", "Hypoxia", "Ferroptosis", "RPMS", "BPMS")
signature_level_cols <- paste0(signature_bases, "_level")

subtype_levels <- sort(unique(tcga_final_df[[subtype_col]]))
num_subtypes <- length(subtype_levels)

if (num_subtypes <= 12 && num_subtypes > 0) {
  palette_name <- "Paired"
  n_colors_to_request <- max(3, num_subtypes)
  if(palette_name %in% rownames(RColorBrewer::brewer.pal.info) && RColorBrewer::brewer.pal.info[palette_name, "maxcolors"] >= n_colors_to_request) {
      subtype_colors <- RColorBrewer::brewer.pal(n_colors_to_request, palette_name)[1:num_subtypes]
  } else {
      palette_name <- "Set1"
      n_colors_to_request <- max(3, num_subtypes)
       if(palette_name %in% rownames(RColorBrewer::brewer.pal.info) && RColorBrewer::brewer.pal.info[palette_name, "maxcolors"] >= n_colors_to_request) {
          subtype_colors <- RColorBrewer::brewer.pal(n_colors_to_request, palette_name)[1:num_subtypes]
       } else {
          subtype_colors <- scales::hue_pal()(num_subtypes)
       }
  }

} else if (num_subtypes > 0) {
  subtype_colors <- scales::hue_pal()(num_subtypes)
} else {
    stop("No valid subtypes found after filtering. Cannot generate palette.")
}

names(subtype_colors) <- subtype_levels
print("Subtype Colors Assigned:")
print(subtype_colors)

output_dir <- "KM_Plots_Signature_Stratified_Subtypes"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Plotting Loop

for (i in seq_along(signature_bases)) {
  sig_base <- signature_bases[i]
  sig_level_col <- signature_level_cols[i]

  if (!sig_level_col %in% colnames(tcga_final_df)) {
    cat("\nSkipping Signature:", sig_base, "- Column", sig_level_col, "not found.\n")
    next
  }
  present_sig_levels <- unique(na.omit(tcga_final_df[[sig_level_col]]))
  positive_level <- present_sig_levels[str_detect(present_sig_levels, "\\+")]
  negative_level <- present_sig_levels[str_detect(present_sig_levels, "\\-")]

  if (length(positive_level) != 1 || length(negative_level) != 1) {
       cat("\nSkipping Signature:", sig_base, "- could not reliably identify one '+' and one '-' level. Found:",
           paste(present_sig_levels, collapse=", "), "\n")
       next
  }


  cat("\nProcessing Signature:", sig_base, "\n")

  tcga_final_df_high <- tcga_final_df %>%
    filter(!is.na(.data[[sig_level_col]]), .data[[sig_level_col]] == positive_level)
  plot_high <- NULL

  if (nrow(tcga_final_df_high) < 2 || length(unique(tcga_final_df_high[[subtype_col]])) < 1) {
      cat("  Skipping HIGH level plot for", sig_base, "- insufficient data or groups.\n")
  } else {
      tcga_final_df_high[[subtype_col]] <- factor(tcga_final_df_high[[subtype_col]], levels = subtype_levels)

      surv_obj_high <- Surv(time = tcga_final_df_high[[time_col]], event = tcga_final_df_high[[status_col]])
      fit_high <- survfit(surv_obj_high ~ get(subtype_col), data = tcga_final_df_high)
      present_subtypes_high <- levels(droplevels(tcga_final_df_high[[subtype_col]]))
      legend_labs_high <- present_subtypes_high
      palette_high <- subtype_colors[present_subtypes_high]

      plot_high <- ggsurvplot(
        fit_high,
        data = tcga_final_df_high,
        title = paste("High", sig_base, "Level"),
        pval = FALSE,
        conf.int = FALSE,
        risk.table = TRUE,
        risk.table.col = "strata",
        risk.table.y.text = FALSE,
        legend.title = "Subtype",
        legend.labs = legend_labs_high,
        palette = palette_high,
        xlab = "Time (Days)",
        ylab = "Survival Probability",
        ggtheme = theme_minimal(base_size = 10) +
                  theme(plot.title = element_text(hjust = 0.5, size=rel(1.2)),
                  axis.text = element_text(face = "bold", size = rel(1.1)),
                        legend.text = element_text(size=rel(1.3), face="bold"),
                        legend.title = element_text(size=rel(1.0), face="bold")),
        break.time.by = if(max(tcga_final_df_high[[time_col]], na.rm=TRUE) > 2000) 1000 else 500
      )
      plot_high$table <- plot_high$table + theme(axis.line = element_blank(),
                                                 axis.ticks = element_blank(),
                                                 axis.text.y = element_blank(),
                                                 axis.title.y = element_blank(),
                                                 plot.title = element_blank(),
                                                 legend.position="none")
  }

  tcga_final_df_low <- tcga_final_df %>%
    filter(!is.na(.data[[sig_level_col]]), .data[[sig_level_col]] == negative_level)

  plot_low <- NULL

  if (nrow(tcga_final_df_low) < 2 || length(unique(tcga_final_df_low[[subtype_col]])) < 1) {
      cat("  Skipping LOW level plot for", sig_base, "- insufficient data or groups.\n")
  } else {
      tcga_final_df_low[[subtype_col]] <- factor(tcga_final_df_low[[subtype_col]], levels = subtype_levels)

      surv_obj_low <- Surv(time = tcga_final_df_low[[time_col]], event = tcga_final_df_low[[status_col]])
      fit_low <- survfit(surv_obj_low ~ get(subtype_col), data = tcga_final_df_low)

      present_subtypes_low <- levels(droplevels(tcga_final_df_low[[subtype_col]]))
      legend_labs_low <- present_subtypes_low
      palette_low <- subtype_colors[present_subtypes_low]

      plot_low <- ggsurvplot(
        fit_low,
        data = tcga_final_df_low,
        title = paste("Low", sig_base, "Level"),
        # *** Removed p-value display ***
        pval = FALSE,
        conf.int = FALSE,
        risk.table = TRUE,
        risk.table.col = "strata",
        risk.table.y.text = FALSE,
        # *** Keep legend info for collection ***
        legend.title = "Subtype",
        legend.labs = legend_labs_low,
        palette = palette_low,
        xlab = "Time (Days)",
        ylab = "Survival Probability",
        ggtheme = theme_minimal(base_size = 10) +
                   theme(plot.title = element_text(hjust = 0.5, size=rel(1.2)),
                   axis.text = element_text(face="bold", size=rel(1.1)),
                        legend.text = element_text(size=rel(1.3), face="bold"),
                        legend.title = element_text(size=rel(1.0), face="bold")),
        break.time.by = if(max(tcga_final_df_low[[time_col]], na.rm=TRUE) > 2000) 1000 else 500
      )
       plot_low$table <- plot_low$table + theme(axis.line = element_blank(),
                                                 axis.ticks = element_blank(),
                                                 axis.text.y = element_blank(),
                                                 axis.title.y = element_blank(),
                                                 plot.title = element_blank(),
                                                 legend.position="none")
  }

  plot_list_to_arrange <- list()
  if (!is.null(plot_high)) plot_list_to_arrange$high <- plot_high
  if (!is.null(plot_low)) plot_list_to_arrange$low <- plot_low

  if (length(plot_list_to_arrange) > 0) {

    arranged_plots_object <- arrange_ggsurvplots(
                                plot_list_to_arrange,
                                print = FALSE,
                                ncol = 2,
                                nrow = 1,
                                title = paste("Overall Survival Stratified by", sig_base, "Level"),
                                legend = "top",
                                risk.table.height = 0.25
                               )

    filename <- file.path(output_dir, paste0("KM_Plot_", sig_base, "_Stratified_by_Subtype.png"))

    plot_width <- if(length(plot_list_to_arrange) == 2) 14 else 7
    ggsave(filename, plot = arranged_plots_object, width = 13, height = 9, dpi = 300, units = "in")

    cat("  Saved combined plot:", filename, "\n")

  } else {
      cat("  No plots generated for", sig_base, "due to insufficient data in both groups.\n")
  }

}
