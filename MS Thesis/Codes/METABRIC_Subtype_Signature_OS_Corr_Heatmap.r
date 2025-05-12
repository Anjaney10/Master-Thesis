# Author: Anjaney J Pandey
# Code for METABRIC_Subtype_Signature_OS_Corr_Heatmap

########################## CORRELATION FOR METABRIC DATA ##########################
library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(forcats)
library(scales)
library(stats)

get_significance_stars <- function(p_value) {
  if (is.na(p_value)) {
    return("")
  }
  if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else {
    return("")
  }
}

signature_score_columns <- c(
    "EMT", "Glycolysis", "OXPHOS", "FAO",
    "Hypoxia", "Ferroptosis", "RPMS", "BPMS"
)
surv_time_col <- "OS_MONTHS"
subtype_col <- "CLAUDIN_SUBTYPE"
required_cols <- c(subtype_col, surv_time_col, signature_score_columns)
output_dir <- "./Final_Plots_n_Tables"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
output_heatmap_file <- file.path(output_dir, "METABRIC_Subtype_Signature_OS_Corr_Heatmap.png")

missing_req_cols <- setdiff(required_cols, colnames(metabric_final_df))
if (length(missing_req_cols) > 0) {
    stop("Error: The following required columns are missing in 'metabric_final_df': ",
         paste(missing_req_cols, collapse=", "))
}

data_for_corr <- metabric_final_df %>%
  filter(
    !is.na(.data[[subtype_col]]) & .data[[subtype_col]] != "Normal" & .data[[subtype_col]] != "NC" & .data[[subtype_col]] != "claudin-low",
    !is.na(.data[[surv_time_col]]) & .data[[surv_time_col]] > 0
  ) %>%
  select(all_of(required_cols)) %>%
  mutate(across(all_of(signature_score_columns), ~ suppressWarnings(as.numeric(as.character(.))))) %>%
  tidyr::drop_na(all_of(signature_score_columns)) %>%
  tidyr::drop_na(all_of(surv_time_col))

min_samples_per_subtype <- 5
valid_subtypes <- subtype_counts %>%
    filter(n >= min_samples_per_subtype) %>%
    pull(.data[[subtype_col]])

if (length(valid_subtypes) == 0) {
  stop(paste("No subtypes have enough data (>=", min_samples_per_subtype, "samples) after filtering NAs."))
}

data_for_corr <- data_for_corr %>% filter(.data[[subtype_col]] %in% valid_subtypes)

# Calculate Spearman Correlations and P-values per Subtype

correlation_results_long <- data_for_corr %>%
  group_by(.data[[subtype_col]]) %>%
  group_modify(~ {
    subtype_data <- .x
    current_subtype_name <- unique(subtype_data[[subtype_col]])
    results <- list()

    for (sig_col in signature_score_columns) {
      valid_x <- subtype_data[[surv_time_col]]
      valid_y <- subtype_data[[sig_col]]
      complete_idx <- complete.cases(valid_x, valid_y)
      n_complete <- sum(complete_idx)

      if (n_complete < 3) {
         results[[sig_col]] <- data.frame(
             Signature = sig_col, Correlation = NA_real_, p.value = NA_real_, n_pairs = n_complete
         )
         next
      }
      var_x <- stats::sd(valid_x[complete_idx], na.rm = TRUE)
      var_y <- stats::sd(valid_y[complete_idx], na.rm = TRUE)

      if (is.na(var_x) || is.na(var_y) || var_x == 0 || var_y == 0) {
         warning(paste("Zero variance for OS time or", sig_col, "in subtype", current_subtype_name, "- skipping correlation for this pair."))
         results[[sig_col]] <- data.frame(
             Signature = sig_col, Correlation = NA_real_, p.value = NA_real_, n_pairs = n_complete
         )
         next
      }
      test_result <- tryCatch(
          stats::cor.test(valid_x, valid_y,
                   method = "spearman", exact = FALSE, use = "pairwise.complete.obs"),
          error = function(e) {
              warning(paste("cor.test failed for", sig_col, "in subtype", current_subtype_name, ":", conditionMessage(e)))
              return(NULL)
          }
      )

      if (!is.null(test_result)) {
        results[[sig_col]] <- data.frame(
          Signature = sig_col,
          Correlation = test_result$estimate,
          p.value = test_result$p.value,
          n_pairs = n_complete
        )
      } else {
         results[[sig_col]] <- data.frame(
             Signature = sig_col, Correlation = NA_real_, p.value = NA_real_, n_pairs = n_complete
         )
      }
    }
    bind_rows(results)
  }, .keep = TRUE) %>%
  ungroup()

if (nrow(correlation_results_long) == 0) {
  stop("Correlation calculation failed or yielded no results for any subtype/signature.")
}

correlation_heatmap_data <- correlation_results_long %>%
  rowwise() %>%
  mutate(significance = get_significance_stars(p.value)) %>%
  ungroup() %>%
  mutate(
    Signature = factor(Signature, levels = rev(sort(unique(Signature)))),
    !!subtype_col := factor(.data[[subtype_col]], levels = sort(valid_subtypes))
  )

print("Correlation results prepared for plotting (showing NAs and stars):")
print(head(correlation_heatmap_data))
print(tail(correlation_heatmap_data))

max_abs_corr <- max(abs(correlation_heatmap_data$Correlation), 0.1, na.rm = TRUE)
corr_limits <- c(-max_abs_corr, max_abs_corr)
corr_limits[1] <- max(-1, corr_limits[1])
corr_limits[2] <- min(1, corr_limits[2])


heatmap_plot <- ggplot(correlation_heatmap_data, aes(y = .data[[subtype_col]], x = Signature, fill = Correlation)) +
  geom_tile(color = "grey85", linewidth=0.3) +
  geom_text(aes(label = significance), color = "black", size = 4, na.rm = TRUE) +
  scale_fill_gradient2(
    name = "Correlation",
    low = muted("blue"), mid = "white", high = muted("red"),
    midpoint = 0,
    limit = corr_limits,
    oob = scales::squish,
    breaks = scales::breaks_pretty(n=5),
    na.value = "grey80"
  ) +
  scale_x_discrete(limits=rev(levels(correlation_heatmap_data$Signature))) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
    axis.text.y = element_text(size = 15),
    axis.title = element_blank(),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, margin=margin(b=10)),
    legend.position = "right",
    legend.title = element_text(size=9, face="bold"),
    legend.text = element_text(size=8),
    legend.key.height = unit(2, "cm")
  ) +
  labs(
    title = "METABRIC Correlations: Signature ssGSEA vs OS Time by BRCA Subtype",
  )

tryCatch({
    ggsave(
        filename = output_heatmap_file,
        plot = heatmap_plot,
        width = 11,
        height = 7,
        units = "in",
        dpi = 300,
        bg = "white"
    )
    print(paste("Heatmap saved to:", output_heatmap_file))
    }, error = function(e) {
         warning(paste("Failed to save heatmap:", conditionMessage(e)))
    })
