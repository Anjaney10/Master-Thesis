# Author: Anjaney J Pandey
# Code for breast cancer survival analysis for different sub types. Different scripts are written for different plots for subtypes

### For TCGA BRCA samples
library(survival)
library(broom)
library(dplyr)
library(ggplot2)

# Download TCGA BRCA clinical data
tcga_brca_clinical <- read.delim("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.BRCA.sampleMap%2FBRCA_clinicalMatrix")

ssgsea_tcga_brca <- read.table("./Transcriptomics/TCGA_DATA/ssgsea_tpm_results/TCGA-BRCA_ssgsea_results.tsv", header = TRUE, sep = "\t")
# Name the first column as Name
names(ssgsea_tcga_brca)[1] <- "Name"

survival_tcga_brca <- read.table("./Transcriptomics/TCGA_DATA/survival_data/BRCA_survival.txt", header = TRUE, sep = "\t")

tcga_final_df <- ssgsea_tcga_brca %>%
  left_join(tcga_brca_clinical %>% select(sampleID, PAM50Call_RNAseq),
            by = c("Name" = "sampleID")) %>%
  left_join(survival_tcga_brca %>% select(-X_PATIENT, -Redaction),
            by = c("Name" = "sample"))

# Unique PAM50 subtypes
unique(tcga_final_df$PAM50Call_RNAseq)
# "Normal" "LumA"   "LumB"   "Basal"  "Her2"

tcga_final_df <- tcga_final_df %>%
 filter(!is.na(PAM50Call_RNAseq) & PAM50Call_RNAseq != "" &  PAM50Call_RNAseq != "Normal" &
         !is.na(OS) & !is.na(OS.time) & OS.time > 0) %>% 
  mutate(Type = factor(Type, levels = c("RPMS-/BPMS-", "RPMS+/BPMS+",
                                          "RPMS+/BPMS-", "RPMS-/BPMS+")))

cox_list <- tcga_final_df %>%
  group_by(PAM50Call_RNAseq) %>% 
  do({
    mod <- coxph(Surv(OS.time, OS) ~ Type, data = .)
    s <- summary(fit.coxph)

    results_metric <- data.frame(
                  metric = metric_id, # Store the metric ID (e.g., "OS", "PFI")
                  HR = signif(s$coefficients[, "exp(coef)"], 3),
                  HR.confint.lower = signif(s$conf.int[, "lower .95"], 3),
                  HR.confint.upper = signif(s$conf.int[, "upper .95"], 3),
                  p.value = signif(s$coefficients[, "Pr(>|z|)"], 3),
                  labels = rownames(s$coefficients),
                  Log_Rank_pvalue = signif(s$logtest[['pvalue']], 3), # Overall model p-value
                  n = s$n,
                  cancer_type = cancer # Add cancer type
              )
              results_metric$significance <- sapply(results_metric$p.value, get_significance_stars)

    tidy_mod <- tidy(mod, exponentiate = TRUE, conf.int = TRUE) %>%
                  filter(term != "(Intercept)")
    tidy_mod
  }) %>% 
  ungroup() %>%
  mutate(sig = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01 ~ "**",
    p.value < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Save cox_list to file
write.table(cox_list, 
            file = "./Transcriptomics/TCGA_DATA/TCGA_BRCA_Subtype_Survival_Analysis/TCGA_BRCA_OS_HR.tsv",
            sep = "\t", 
            row.names = FALSE,
            quote = FALSE)

ggplot(cox_list, aes(x = estimate, y = term)) +
  geom_point() +
  scale_y_discrete(labels = function(x) gsub("Type", "", x)) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0) +
  geom_text(aes(label = paste0("p=", round(p.value, 3), " ", sig),
                x = conf.high),
            position = position_nudge(y = 0.1),
            hjust = 1.5, size = 3) +
  #facet_wrap(~ PAM50Call_RNAseq, scales = "free_x") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  scale_x_log10(expand = expansion(mult = c(0.05, 0.2))) + 
  coord_cartesian(clip = "off") +
  labs(x = "Hazard Ratio", y = "Type",
       title = "TCGA BRCA\nOverall Survival") +
  theme_bw() +
  theme(plot.margin = margin(10, 40, 10, 10))

# Save the plot
ggsave(filename = "./Transcriptomics/TCGA_DATA/TCGA_BRCA_Subtype_Survival_Analysis/TCGA_BRCA_OS.png",
       plot = last_plot(),
       width = 10, 
       height = 7, 
       dpi = 300,      
       type = "cairo",
       bg = "white")

################################################################################################

library(dplyr)
library(survival)
library(ggplot2)
library(rlang)
library(stringr)

metrics_config <- list(
  OS = list(time_col = "OS.time", status_col = "OS", pretty_name = "Overall Survival", short_name = "OS"),
  PFI = list(time_col = "PFI.time", status_col = "PFI", pretty_name = "Progression Free Interval", short_name = "PFI"),
  tcga_final_dfI = list(time_col = "tcga_final_dfI.time", status_col = "tcga_final_dfI", pretty_name = "Disease Free Interval", short_name = "tcga_final_dfI"),
  DSS = list(time_col = "DSS.time", status_col = "DSS", pretty_name = "Disease Specific Survival", short_name = "DSS")
)
base_output_dir <- "./Transcriptomics/TCGA_DATA/TCGA_BRCA_Subtype_Survival_Analysis"
dir.create(base_output_dir, showWarnings = FALSE, recursive = TRUE)

tcga_filtered_base <- tcga_final_df %>%
  filter(!is.na(PAM50Call_RNAseq) & PAM50Call_RNAseq != "" & PAM50Call_RNAseq != "Normal") %>%
  mutate(Type = factor(Type, levels = c("RPMS-/BPMS-", "RPMS+/BPMS+", "RPMS+/BPMS-", "RPMS-/BPMS+")))

if (nrow(tcga_filtered_base) == 0) {
    stop("Initial filtering removed all data. Check PAM50Call_RNAseq and Type columns.")
}

all_metrics_results <- list()

for (metric_id in names(metrics_config)) {
  metric_info <- metrics_config[[metric_id]]
  time_col <- metric_info$time_col
  status_col <- metric_info$status_col
  pretty_name <- metric_info$pretty_name
  short_name <- metric_info$short_name

  print(paste("Processing Metric:", pretty_name))

  if (!all(c(time_col, status_col) %in% colnames(tcga_filtered_base))) {
      warning(paste("Columns", time_col, "or", status_col, "not found. Skipping metric:", pretty_name))
      next
  }
  tcga_filtered_metric <- tcga_filtered_base %>%
    filter(!is.na(!!sym(time_col)) & !is.na(!!sym(status_col)) & !!sym(time_col) > 0)

  if (nrow(tcga_filtered_metric) == 0) {
      warning(paste("No data remaining after filtering for metric:", pretty_name))
      next
  }
  print(paste("Rows after filtering for", short_name, ":", nrow(tcga_filtered_metric)))

  cox_results_metric <- tcga_filtered_metric %>%
    group_by(PAM50Call_RNAseq) %>%
    do({
      group_data <- .
      results_group <- tryCatch({
          surv_obj_str <- paste0("Surv(", time_col, ", ", status_col, ")")
          formula_str <- paste0(surv_obj_str, " ~ Type")
          mod <- coxph(as.formula(formula_str), data = group_data)

          s <- summary(mod)

          if(is.null(s$coefficients) || nrow(s$coefficients) == 0) {
             warning(paste("Cox summary empty for subtype", unique(group_data$PAM50Call_RNAseq), "metric", short_name))
             return(NULL)
          }

          tcga_final_df_out <- data.frame(
                    term = rownames(s$coefficients),
                    estimate = signif(s$coefficients[, "coef"], 3),
                    HR = signif(s$coefficients[, "exp(coef)"], 3),
                    conf.low = signif(s$conf.int[, "lower .95"], 3),
                    conf.high = signif(s$conf.int[, "upper .95"], 3),
                    p.value = signif(s$coefficients[, "Pr(>|z|)"], 3),
                    Log_Rank_pvalue = signif(s$logtest[['pvalue']], 3),
                    n = s$n,
                    stringsAsFactors = FALSE
                  )

          tcga_final_df_out
        },
        error = function(e) {
          warning(paste("Cox model failed for subtype", unique(group_data$PAM50Call_RNAseq), "and metric", short_name, ":", conditionMessage(e)))
          return(NULL)
        })
      results_group
    }) %>%
    ungroup() %>%
    filter(!is.null(.))

  cox_results_metric <- cox_results_metric %>%
      mutate(sig = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01 ~ "**",
        p.value < 0.05 ~ "*",
        TRUE ~ ""
      ),
      metric = short_name
    )

  all_metrics_results[[metric_id]] <- cox_results_metric

  results_filename <- file.path(base_output_dir, paste0("TCGA_BRCA_", short_name, "_HR.tsv"))
  write.table(cox_results_metric,
              file = results_filename,
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)
  print(paste("Results table saved to:", results_filename))

#############################################################################################################
library(dplyr)
library(survival)
library(ggplot2)
library(rlang)
library(stringr)
library(scales)

metrics_config <- list(
  OS = list(time_col = "OS.time", status_col = "OS", pretty_name = "Overall Survival", short_name = "OS"),
  PFI = list(time_col = "PFI.time", status_col = "PFI", pretty_name = "Progression Free Interval", short_name = "PFI"),
  tcga_final_dfI = list(time_col = "tcga_final_dfI.time", status_col = "tcga_final_dfI", pretty_name = "Disease Free Interval", short_name = "tcga_final_dfI"),
  DSS = list(time_col = "DSS.time", status_col = "DSS", pretty_name = "Disease Specific Survival", short_name = "DSS")
)
base_output_dir <- "./Transcriptomics/TCGA_DATA/TCGA_BRCA_Subtype_Survival_Analysis"
dir.create(base_output_dir, showWarnings = FALSE, recursive = TRUE)


# Significance Function
get_significance_stars <- function(p_values) {
  case_when(
    is.na(p_values) ~ "",
    p_values < 0.001 ~ "***",
    p_values < 0.01 ~ "**",
    p_values < 0.05 ~ "*",
    TRUE ~ ""
  )
}

all_metrics_results <- list()

for (metric_id in names(metrics_config)) {
  metric_info <- metrics_config[[metric_id]]
  time_col <- metric_info$time_col
  status_col <- metric_info$status_col
  pretty_name <- metric_info$pretty_name
  short_name <- metric_info$short_name

  if (!all(c(time_col, status_col) %in% colnames(tcga_filtered_base))) {
      warning(paste("Columns", time_col, "or", status_col, "not found. Skipping metric:", pretty_name)); next }
  tcga_filtered_metric <- tcga_filtered_base %>%
    filter(!is.na(!!sym(time_col)) & !is.na(!!sym(status_col)) & !!sym(time_col) > 0)
  if (nrow(tcga_filtered_metric) == 0) {
      warning(paste("No data remaining after filtering for metric:", pretty_name)); next }
  print(paste("Rows after filtering for", short_name, ":", nrow(tcga_filtered_metric)))

  cox_results_metric <- tcga_filtered_metric %>%
    filter(PAM50Call_RNAseq %in% c("LumA", "LumB", "Her2", "Basal")) %>%
    group_by(PAM50Call_RNAseq) %>%
    do({
      group_data <- .
      if(n_distinct(group_data$Type) < 2 || nrow(group_data) < 5) {
         warning(paste("Skipping Cox model for subtype", unique(group_data$PAM50Call_RNAseq), "- insufficient data/levels for metric", short_name))
         return(NULL)
      }
      results_group <- tryCatch({
          surv_obj_str <- paste0("Surv(", time_col, ", ", status_col, ")")
          formula_str <- paste0(surv_obj_str, " ~ Type")
          mod <- coxph(as.formula(formula_str), data = group_data)
          s <- summary(mod)
          if(is.null(s$coefficients) || nrow(s$coefficients) == 0) {
             warning(paste("Cox summary empty for subtype", unique(group_data$PAM50Call_RNAseq), "metric", short_name)); return(NULL) }
          tcga_final_df_out <- data.frame( term = rownames(s$coefficients), estimate = signif(s$coefficients[, "coef"], 3),
                    HR = signif(s$coefficients[, "exp(coef)"], 3), conf.low = signif(s$conf.int[, "lower .95"], 3),
                    conf.high = signif(s$conf.int[, "upper .95"], 3), p.value = signif(s$coefficients[, "Pr(>|z|)"], 3),
                    Log_Rank_pvalue = signif(s$logtest[['pvalue']], 3), n = s$n, stringsAsFactors = FALSE )
          tcga_final_df_out
        }, error = function(e) { warning(paste("Cox model failed for subtype", unique(group_data$PAM50Call_RNAseq), "metric", short_name,":", conditionMessage(e))); return(NULL) })
      results_group
    }) %>% ungroup() %>% filter(!is.null(.))

  if (nrow(cox_results_metric) == 0) {
      warning(paste("Cox model failed or yielded no results for any included subtype for metric:", pretty_name)); next }

  cox_results_metric <- cox_results_metric %>%
      mutate(
          significance = get_significance_stars(p.value),
          metric = short_name,
          plot_labels = gsub("Type", "", term),
          plot_labels = factor(plot_labels, levels = rev(c("RPMS+/BPMS+", "RPMS+/BPMS-", "RPMS-/BPMS+"))),
          logrank_text = paste("LogRank p =", signif(Log_Rank_pvalue, 2))
       )

  # Store results
  all_metrics_results[[metric_id]] <- cox_results_metric
  results_filename <- file.path(base_output_dir, paste0("TCGA_BRCA_", short_name, "_HR_by_Subtype.tsv"))
  write.table(cox_results_metric, file = results_filename, sep = "\t", row.names = FALSE, quote = FALSE)
  print(paste("Results table saved to:", results_filename))


# Generate Combined Plot
cox_results_metric <- cox_results_metric %>%
    mutate(
        plot_labels = gsub("Type", "", term),
        plot_labels = factor(plot_labels, levels = rev(c("RPMS+/BPMS+", "RPMS+/BPMS-", "RPMS-/BPMS+"))),
        facet_title = paste0(PAM50Call_RNAseq, "\n", "LogRank p = ", signif(Log_Rank_pvalue, 2))
     )
if(any(is.na(cox_results_metric$plot_labels))) {
    warning("NA values generated in plot labels ('plot_labels') for metric ", short_name) }

n_subtypes_plot <- n_distinct(cox_results_metric$PAM50Call_RNAseq)
plot_ncol <- if (n_subtypes_plot <= 4) 2 else ceiling(sqrt(n_subtypes_plot))

faceted_forest_plot <- ggplot(cox_results_metric, aes(x = plot_labels, y = HR)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray60", linewidth=0.5) +
    geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                    color = "black", size = 0.6, linewidth = 0.8) +
    geom_text(aes(label = paste0("p=", round(p.value, 3), " ", significance),
                  y = conf.high),
              color = "black",
              vjust = -0.7,
              hjust = 1.5,
              size = 4.5) +
    facet_wrap(~ facet_title,
               ncol = plot_ncol,
               scales = "free_x") +
    coord_flip() +
    theme_bw(base_size = 12) +
    theme(
        legend.position = "none",
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        strip.background = element_rect(fill="grey90", color=NA),
        strip.text = element_text(face="bold", size=10, lineheight=1.0, margin = margin(t = 2, b = 2)),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        plot.margin = margin(t = 15, r = 10, b = 10, l = 10, unit = "pt")
    ) +
    labs(
        title = paste0("TCGA BRCA Subtype Analysis\n", pretty_name),
        x = NULL,
        y = "Hazard Ratio (HR)"
    )

  # Save the Plot
  plot_filename <- file.path(base_output_dir, paste0("TCGA_BRCA_", short_name, "_Subtype_Comparison_Plot.png"))


  tryCatch({
      ggsave(filename = plot_filename,
             plot = faceted_forest_plot,
             width = 10,
             height = 7,
             dpi = 300,
             bg = "white")
       print(paste("Faceted plot saved to:", plot_filename))
     }, error = function(e){
          warning(paste("Failed to save faceted plot for", metric_id, ":", conditionMessage(e)))
     })

}



###################################################  For Metabric BRCA samples  ###########################################################

metabric_data <- read.table("./Transcriptomics/brca_metabric/data_mrna_illumina_microarray_zscores_ref_diploid_samples_final_unique.txt",
header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)


RPMS_genes <- scan(text = "PEBP1\tIGF2BP2\tGOLT1B\tNRAS\tRDX\tNF2\tDMD\tBACH1\tMAP4K4\tCOL5A2\tCDV3\tCOL3A1\tHMGA2\tCOL1A2\tTGFBR1\tPAPPA\tIGF2BP3\tARID3B\tLRIG2\tIRS2\tHIC2\tIL13\tGNPTAB\tKRT16\tTNFRSF12A\tENO1\tEN1\tTAGLN2\tBAG2\tHMGA1\tDSC2\tGPR56\tLAMC1\tS100A10\tPKP3\tNDRG2\tEMP3\tEIF4G1\tPSMD4\tPIM1\tGTPBP2\tFADS3\tOSBP\tEPHB2\tPSME4\tDYSF\tBRD2\tMAST2\tMAP2\tUBE2E3\tSQSTM1\tRGS2\tSYNGR1\tELK3\tUBE2C\tDIAPH1\tABCD1\tCOL7A1\tCAPN6\tMMP3\tLAMB3\tPPP1R15A\tGNAI1\tRPS6KA4\tABCA2\tF2RL2\tLAMC2\tTRAPPC3\tOR2F1\tAKAP1\tFAP\tGIT1\tKRT8\tRAB3D\tATP1B1\tELAVL2\tPHLDA2\tPCDH17\tENO2\tBDKRB1\tCA7\tGKN1\tADRA1A\tPSMD7\tBMP2\tGPR3\tBTK\tSNCB\tLPP\tNR1I3\tHOXB3\tCDKN1A\tTIAL1\tAPOBEC1\tDLG4\tIL1RN\tMPV17\tHCN3\tSHC3\tRPL23A\tKCNA2\tDYRK1A\tFOXN1\tDCTN2\tPADI4\tMARK1\tRAB30\tGFAP\tAKT3\tEEF1A1\tLMNA\tLTB4R2\tSCRN1\tHMGA2\tMMP1\tSPP1\tCXCR4", what = "character", sep = "\t", quiet = TRUE)
BPMS_genes <- scan(text = "BMPER\tDYM\tFBXO42\tFRMPD4\tHERC3\tHS3ST3B1\tIL1RAP\tIL7\tMAGEC1\tMYCT1\tPDE1C\tPRDM1\tRCAN3", what = "character", sep = "\t", quiet = TRUE)

library(GSVA)
library(tibble)
library(dplyr)

# ssGSEA analysis

genesets <- list(RPMS_genes = RPMS_genes, BPMS_genes = BPMS_genes)

metabric_data_matrix <- as.matrix(metabric_data)
param <- ssgseaParam(exprData = metabric_data_matrix, geneSets = genesets, minSize = 5, use = "na.rm")
ssgsea_results <- gsva(param, verbose = TRUE)

# Transpose the results
ssgsea_results <- t(ssgsea_results)
ssgsea_results <- as.data.frame(ssgsea_results)

# Add a Type column to the results based on RPMS and BPMS genes expression levels compared to median expression
rpms_median <- median(ssgsea_results$RPMS_genes)
bpms_median <- median(ssgsea_results$BPMS_genes)

ssgsea_results$Type <- ifelse(ssgsea_results$RPMS_genes > rpms_median & ssgsea_results$BPMS_genes > bpms_median, "RPMS+/BPMS+",
                 ifelse(ssgsea_results$RPMS_genes > rpms_median & ssgsea_results$BPMS_genes <= bpms_median, "RPMS+/BPMS-",
                 ifelse(ssgsea_results$RPMS_genes <= rpms_median & ssgsea_results$BPMS_genes > bpms_median, "RPMS-/BPMS+",
                       "RPMS-/BPMS-")))

# Save the results
write.table(ssgsea_results, file = "./Transcriptomics/brca_metabric/ssgsea_brca_metabric.txt", sep = "\t", quote = FALSE)

metabric_clincial <- read.table("./Transcriptomics/brca_metabric/data_clinical_patient.txt", header = TRUE, sep = "\t", skip = 1)

metabric_ssgsea <- read.table("./Transcriptomics/brca_metabric/METABRIC_all_signatures_ssgsea_results.tsv", header = TRUE, sep = "\t")


metabric_final_df <- metabric_ssgsea %>% 
  left_join(metabric_clincial %>% select(PATIENT_ID, CLAUDIN_SUBTYPE, OS_MONTHS, OS_STATUS, RFS_MONTHS, RFS_STATUS),
            by = c("Name" = "PATIENT_ID"))

# Keep only the number before : from OS_STATUS and RFS_STATUS columns
metabric_final_df$OS_STATUS <- as.numeric(gsub(":.*", "", metabric_final_df$OS_STATUS))
metabric_final_df$RFS_STATUS <- as.numeric(gsub(":.*", "", metabric_final_df$RFS_STATUS))

# Unique PAM50 subtypes
unique(metabric_final_df$CLAUDIN_SUBTYPE)
# "LumA"        "Her2"        "LumB"        "claudin-low" "Basal"  "Normal"      "NC"  

library(survival)
library(broom)
library(dplyr)
library(ggplot2)

# Set the levels of Type column
metabric_final_df <- metabric_final_df %>% filter(!is.na(CLAUDIN_SUBTYPE) & CLAUDIN_SUBTYPE != "NC" &
                                                  !is.na(OS_STATUS) & !is.na(OS_MONTHS) & OS_MONTHS > 0 &
                                                  !is.na(RFS_STATUS) & !is.na(RFS_MONTHS) & RFS_MONTHS > 0) %>% 
                                                  mutate(Type = factor(Type, levels = c("RPMS-/BPMS-", "RPMS+/BPMS+",
                                          "RPMS+/BPMS-", "RPMS-/BPMS+")))

cox_list <- metabric_final_df %>%
  filter(!is.na(CLAUDIN_SUBTYPE) & CLAUDIN_SUBTYPE != "") %>%
  group_by(CLAUDIN_SUBTYPE) %>% 
  do({
    mod <- coxph(Surv(RFS_MONTHS, RFS_STATUS) ~ Type, data = .)
    tidy_mod <- tidy(mod, exponentiate = TRUE, conf.int = TRUE) %>%
                  filter(term != "(Intercept)")
    #tidy_mod$CLAUDIN_SUBTYPE <- unique(.$CLAUDIN_SUBTYPE)
    tidy_mod
  }) %>%
  ungroup() %>%
  mutate(sig = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01 ~ "**",
    p.value < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Save cox_list to file
write.table(cox_list, 
            file = "./Transcriptomics/TCGA_DATA/METABRIC_Subtype_Survival_Analysis/METABRIC_Subtype_OS_HR.tsv",
            sep = "\t", 
            row.names = FALSE,
            quote = FALSE)

ggplot(cox_list, aes(x = estimate, y = term)) +
  geom_point() +
  scale_y_discrete(labels = function(x) gsub("Type", "", x)) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0) +
  geom_text(aes(label = paste0("p=", round(p.value, 3), " ", sig),
                x = conf.high),
            position = position_nudge(y = 0.1),
            hjust = 1.5, size = 3) +
 # facet_wrap(~ CLAUDIN_SUBTYPE, scales = "free_x") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  scale_x_log10(expand = expansion(mult = c(0.05, 0.2))) + 
  coord_cartesian(clip = "off") +
  labs(x = "Hazard Ratio", y = "Type",
       title = "METABRIC BRCA Overall Survival") +
  theme_bw() +
  theme(plot.margin = margin(10, 40, 10, 10))

# Save the plot
ggsave(filename = "./Transcriptomics/TCGA_DATA/METABRIC_Subtype_Survival_Analysis/Metabric_Subtype_OS.png",
       plot = last_plot(),
       width = 10, 
       height = 8, 
       dpi = 600,      
       type = "cairo",
       bg = "white")



#################################### KAPLAN MEIER PLOT FOR TCGA BRCA SUBTYPE WISE ######################################
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(tidyr)
library(broom)
library(tibble)

# Filter and prepare data
km_data <- tcga_final_df %>%
  filter(!is.na(PAM50Call_RNAseq) & 
         PAM50Call_RNAseq != "" &
         !is.na(OS) & 
         !is.na(OS.time) & 
         OS.time > 0) %>%
  mutate(PAM50Call_RNAseq = factor(PAM50Call_RNAseq))

# Create Kaplan-Meier plot
km_plot <- ggsurvplot(
  fit = survfit(Surv(OS.time, OS) ~ PAM50Call_RNAseq, data = km_data),
  data = km_data,
  pval = FALSE,
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.height = 0.25,
  ggtheme = theme_bw(),
  palette = c("LumA" = "#00BFC4", "LumB" = "#7CAE00", 
              "Basal" = "#F8766D", "Her2" = "#C77CFF", "Normal" = "black"),
  linetype = c("LumA" = "twodash", "LumB" = "dashed",
              "Basal" = "dotted", "Her2" = "longdash", "Normal" = "solid"),
  xlab = "Time (Days)",
  ylab = "Overall Survival Probability",
  title = "TCGA BRCA Survival by PAM50 Subtype",
  legend.title = "PAM50 Subtype",
  font.main = c(14, "bold"),
  font.x = c(12),
  font.y = c(12),
  font.legend = c(10),
  risk.table.fontsize = 3
)

# Save the plot
ggsave(filename = "./Transcriptomics/TCGA_DATA/TCGA_BRCA_KM_OS_Subtype.png",
       plot = km_plot$plot,
       width = 8,
       height = 10,
       dpi = 300,
       bg = "white")
