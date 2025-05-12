# Author: Anjaney J Pandey
# Code for correlation between RPMS and BPMS across cancer types

library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(stringr)
library(stats)

data_dir <- "./Transcriptomics/TCGA_DATA/ssgsea_tpm_results"
output_dir <- "./Transcriptomics/TCGA_DATA/"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

pathways <- c("EMT", "Glycolysis", "OXPHOS", "FAO", "Hypoxia", "Ferroptosis")
comparison_vars <- c("RPMS", "BPMS")

correlation_pairs <- list()
for (pathway in pathways) {
  for (comp_var in comparison_vars) {
    correlation_pairs[[length(correlation_pairs) + 1]] <- c(pathway, comp_var)
  }
}
correlation_pairs[[length(correlation_pairs) + 1]] <- c("RPMS", "BPMS")

all_files <- list.files(path = data_dir, pattern = "^TCGA-[A-Z0-9]+_ssgsea_results\\.tsv$", full.names = TRUE, ignore.case = TRUE)

if (length(all_files) == 0) {
  stop("No files matching 'TCGA-XXX_ssgsea_results.tsv' found in: ", data_dir)
}

print(paste("Found", length(all_files), "TCGA ssgsea result files."))

correlation_results_list <- list()

for (file_path in all_files) {
  file_name <- basename(file_path)
  cancer_type <- stringr::str_extract(file_name, "^TCGA-[A-Z0-9]+")

  if (is.na(cancer_type)) {
      warning(paste("Could not extract TCGA cancer type from:", file_name, "- Skipping."))
      next
  }

  print(paste("Processing:", cancer_type))
  data <- NULL
  tryCatch({
      data <- read.delim(file_path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, row.names=NULL)
  }, error = function(e) {
      warning(paste("Could not read file:", file_name, "as TSV. Skipping. Error:", e$message))
      data <<- NULL
  })
  required_cols <- unique(unlist(correlation_pairs))
  if (!all(required_cols %in% colnames(data))) {
    missing_cols <- setdiff(required_cols, colnames(data))
    warning(paste("Missing required columns in", file_name, ":", paste(missing_cols, collapse=", "), ". Skipping correlations for this file."))
    next
  }
  for (pair in correlation_pairs) {
    var1 <- pair[1]
    var2 <- pair[2]

    cor_test_result <- tryCatch({
        valid_data <- data[, c(var1, var2)]
        valid_data <- valid_data[complete.cases(valid_data), ]
        cor.test(valid_data[[var1]], valid_data[[var2]], method = "spearman", exact = FALSE)
        }, error = function(e) {
        NULL
    })

    if (!is.null(cor_test_result)) {
      correlation_results_list[[length(correlation_results_list) + 1]] <- data.frame(
        CancerType = cancer_type,
        Variable1 = var1,
        Variable2 = var2,
        Pair = paste(var1, "vs", var2),
        SpearmanRho = cor_test_result$estimate,
        Pvalue = cor_test_result$p.value,
        stringsAsFactors = FALSE
      )
    } else {
       correlation_results_list[[length(correlation_results_list) + 1]] <- data.frame(
        CancerType = cancer_type,
        Variable1 = var1,
        Variable2 = var2,
        Pair = paste(var1, "vs", var2),
        SpearmanRho = NA,
        Pvalue = NA,
        stringsAsFactors = FALSE
      )
    }
  }
}

if (length(correlation_results_list) == 0) {
  stop("No correlations could be calculated. Check input files and column names.")
}
correlation_results_df <- bind_rows(correlation_results_list)

output_table_file <- file.path(output_dir, "TCGA_Pathway_RPMS_BPMS_SpearmanCorrelations.tsv")
write.table(correlation_results_df, file = output_table_file, sep = "\t", row.names = FALSE, quote = FALSE)
print(paste("Correlation results saved to:", output_table_file))

cor_matrix <- correlation_results_df %>%
  select(CancerType, Pair, SpearmanRho) %>%
  pivot_wider(names_from = CancerType, values_from = SpearmanRho) %>%
  tibble::column_to_rownames(var = "Pair")

pval_matrix <- correlation_results_df %>%
  select(CancerType, Pair, Pvalue) %>%
  pivot_wider(names_from = CancerType, values_from = Pvalue) %>%
  tibble::column_to_rownames(var = "Pair")

common_pairs <- intersect(rownames(cor_matrix), rownames(pval_matrix))
common_cancers <- intersect(colnames(cor_matrix), colnames(pval_matrix))

cor_matrix <- cor_matrix[common_pairs, common_cancers, drop = FALSE]
pval_matrix <- pval_matrix[common_pairs, common_cancers, drop = FALSE]

colnames(cor_matrix) <- gsub("^TCGA-", "", colnames(cor_matrix))
colnames(pval_matrix) <- gsub("^TCGA-", "", colnames(pval_matrix))

significance_matrix <- matrix(ifelse(pval_matrix < 0.001 & !is.na(pval_matrix), "***",
                              ifelse(pval_matrix < 0.01 & !is.na(pval_matrix), "**",
                              ifelse(pval_matrix < 0.05 & !is.na(pval_matrix), "*", ""))),
                              nrow = nrow(pval_matrix))
rownames(significance_matrix) <- rownames(pval_matrix)
colnames(significance_matrix) <- colnames(pval_matrix)

rpms_bpms <- "RPMS vs BPMS"
oxphos_pairs <- grep("^OXPHOS", rownames(cor_matrix), value = TRUE)
fao_pairs <- grep("^FAO", rownames(cor_matrix), value = TRUE)
all_other_pairs <- setdiff(rownames(cor_matrix), c(rpms_bpms, oxphos_pairs, fao_pairs))

desired_row_order <- c(rpms_bpms, sort(all_other_pairs), sort(oxphos_pairs), sort(fao_pairs))
desired_row_order <- desired_row_order[desired_row_order %in% rownames(cor_matrix)]

cor_matrix <- cor_matrix[desired_row_order, , drop = FALSE]
significance_matrix <- significance_matrix[desired_row_order, , drop = FALSE]
sorted_cancer_types <- sort(colnames(cor_matrix))
cor_matrix <- cor_matrix[, sorted_cancer_types, drop = FALSE]
significance_matrix <- significance_matrix[, sorted_cancer_types, drop = FALSE]


# Generate Heatmap
color_palette <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)
breaks <- seq(-1, 1, length.out = 101)
output_plot_file_png <- file.path(output_dir, "TCGA_Pathway_RPMS_BPMS_SpearmanCorrelation_Heatmap.png")

# Save as PNG
png(output_plot_file_png, width = 7, height = 9, units = "in", res = 300)
pheatmap(
  mat = cor_matrix,
  color = color_palette,
  breaks = breaks
  border_color = "grey60",
  cellwidth = NA,
  cellheight = NA,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = significance_matrix,
  fontsize_number = 8,
  fontsize_row = 10,
  fontsize_col = 10,
  angle_col = 45,
  legend = TRUE,
  legend_breaks = c(-1, -0.5, 0, 0.5, 1),
  legend_labels = c("-1.0", "-0.5", "0.0", "0.5", "1.0"),
  main = "Spearman Correlation: Pathways vs RPMS/BPMS across TCGA Cancers",
  na_col = "grey80"
)

dev.off()
