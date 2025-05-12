# Author: Anjaney J Pandey
# The goal is to visualize the hazard ratio for different breast cancer subtypes as dictated by the different signatures
# The code can be re-used for TCGA-BRCA by replacing the file names accordingly

library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(stringr)
library(rlang)
library(tidyr)

input_df_name <- "metabric_final_df"

# Clinical data columns
status_col <- "OS_STATUS"
time_col <- "OS_MONTHS"
time_units <- "Months"

# Subtype columns
subtype_col <- "CLAUDIN_SUBTYPE"
subtypes_to_plot <- c("LumA", "LumB", "Her2", "Basal")

# Signatures
signature_bases <- c("EMT", "Glycolysis", "OXPHOS", "FAO", "Hypoxia", "Ferroptosis", "RPMS", "BPMS")
signature_level_cols <- paste0(signature_bases, "_level")

# Output
output_dir <- "Heatmaps"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
output_filename <- file.path(output_dir, "METABRIC_HR_Heatmap_metabric_Subtypes_Signatures.png")

# Significance stars
get_signif_stars <- function(p_value) {
  if (is.na(p_value)) { return("") }
  if (p_value < 0.001) { return("***") }
  if (p_value < 0.01)  { return("**") }
  if (p_value < 0.05)  { return("*") }
  return("")
}

# Load data
metabric_df_base <- get(input_df_name)

# Verify Status column
status_values <- unique(na.omit(metabric_df_base[[status_col]]))
numeric_status_col <- status_col

# Filter data
metabric_subtype_filtered <- metabric_df_base %>%
  filter(!is.na(!!sym(subtype_col)) &
         !!sym(subtype_col) %in% subtypes_to_plot) %>%
  mutate(!!sym(subtype_col) := factor(!!sym(subtype_col), levels = subtypes_to_plot)) %>%
  filter(!is.na(!!sym(time_col)), !is.na(!!sym(numeric_status_col)))

if (nrow(metabric_subtype_filtered) == 0) { stop("No data remaining after filtering.") }


# Calculate HR and P-values for each combination
results_list <- list()

for (i in seq_along(signature_bases)) {
  sig_base <- signature_bases[i]
  sig_level_col <- signature_level_cols[i]
  cat("Processing Signature:", sig_base, "\n")

  if (!sig_level_col %in% colnames(metabric_subtype_filtered)) {
    cat("  Signature column", sig_level_col, "not found. Skipping.\n")
    next
  }
  df_signature_filtered <- metabric_subtype_filtered %>%
      filter(!is.na(!!sym(sig_level_col)) & !!sym(sig_level_col) != "")

  actual_levels <- unique(df_signature_filtered[[sig_level_col]])
  if (length(actual_levels) != 2) { cat("  Not exactly 2 levels found. Skipping subtype loop.\n"); next }
  positive_level <- actual_levels[str_detect(actual_levels, "\\+")]
  negative_level <- actual_levels[str_detect(actual_levels, "\\-")]
  if (length(positive_level) != 1 || length(negative_level) != 1) { cat("  Could not ID +/- levels. Skipping subtype loop.\n"); next }

  for (subtype in subtypes_to_plot) {
    cat("  Processing subtype:", subtype, "...\n")

    df_subtype_signature <- df_signature_filtered %>%
      filter(!!sym(subtype_col) == subtype)

    n_samples <- nrow(df_subtype_signature)
    n_levels_subtype <- length(unique(df_subtype_signature[[sig_level_col]]))

    hr_val <- NA_real_
    p_val <- NA_real_

    {
        # Create Surv object
        surv_obj_cell <- Surv(time = df_subtype_signature[[time_col]], event = df_subtype_signature[[numeric_status_col]])

        df_cox_cell <- df_subtype_signature %>%
            mutate(!!sym(sig_level_col) := factor(!!sym(sig_level_col), levels = c(negative_level, positive_level)))

        # Fit Cox model
        cox_formula <- as.formula(paste("surv_obj_cell ~", rlang::sym(sig_level_col)))
        cox_model <- tryCatch({
            coxph(cox_formula, data = df_cox_cell)
        }, error = function(e) {
            cat("    Cox model failed:", conditionMessage(e), "\n")
            NULL
        })

        if (!is.null(cox_model)) {
            cox_summary <- summary(cox_model)
            # Extract HR and P-value
            if ("coefficients" %in% names(cox_summary) && !is.null(cox_summary$conf.int) && nrow(cox_summary$conf.int) > 0) {
                hr_val <- cox_summary$conf.int[1, "exp(coef)"]
                p_val <- cox_summary$sctest["pvalue"]
            } else {
                 cat("    Cox summary missing components.\n")
            }
        }
    } else {
      cat("    Skipping Cox model due to insufficient data (n=", n_samples, ", levels=", n_levels_subtype, ").\n")
    }

    results_list[[length(results_list) + 1]] <- data.frame(
      Subtype = subtype,
      Signature = sig_base,
      HR = hr_val,
      PValue = p_val,
      N = n_samples
    )

  }
}

# Combine results
results_df <- bind_rows(results_list)

# Prepare Data for Plotting
results_df <- results_df %>%
  mutate(
    SignifStars = sapply(PValue, get_signif_stars),
    HR_Formatted = ifelse(is.na(HR), "NA", sprintf("%.2f", HR)),
    Subtype = factor(Subtype, levels = rev(subtypes_to_plot)),
    Signature = factor(Signature, levels = signature_bases)
  )

hr_range <- results_df$HR[!is.na(results_df$HR) & is.finite(results_df$HR)]
if (length(hr_range) > 0) {
    max_abs_log_hr <- max(abs(log2(hr_range[hr_range > 0])), na.rm = TRUE)
    color_limit <- max(1.1, ceiling(max(hr_range)))
} else {
    color_limit <- 2
}


# Create Heatmap
p <- ggplot(results_df, aes(x = Signature, y = Subtype, fill = HR)) +
    geom_tile(color = "grey80") +
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 1,
      limit = c(0, color_limit),
      na.value = "grey90",
      name = "Hazard Ratio (HR)"
    ) +
    geom_text(aes(label = SignifStars), color = "black", size = 5, vjust = 0.75) +
    labs(
      title = "Hazard Ratios of Signatures within Breast Cancer Subtypes (METABRIC)",
      x = "Gene Signature",
      y = "Breast Cancer Subtype"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )

# Save Plot
ggsave(output_filename, plot = p, width = 10, height = 7, dpi = 300)
