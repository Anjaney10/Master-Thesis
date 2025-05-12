library(TCGAbiolinks)
library(SummarizedExperiment)
library(ggplot2)

# List of TCGA cancer projects
projects <- TCGAbiolinks:::getGDCprojects()$project_id
elf3_expr_list <- list()

projects <- projects[grepl('^TCGA',projects,perl=T)]

# Loop through each project and download ELF3 expression data
for (project in projects) {
  query <- GDCquery(
    project = projects,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )

  tcga_data <- GDCprepare(query)

  elf3_expr <- assay(tcga_data)["ELF3", ]
  elf3_expr_list[[project]] <- data.frame(SampleID = colnames(tcga_data), ELF3_expr = elf3_expr, CancerType = projects)
}

# Combine expression data from all cancer types
elf3_expr_data <- do.call(rbind, elf3_expr_list)

correlation_data <- elf3_expr_data %>%
  group_by(CancerType) %>%
  summarize(
    corr_epithelial = cor(ELF3_expr, epithelial_scores, method = "pearson"),
    corr_mesenchymal = cor(ELF3_expr, mesenchymal_scores, method = "pearson")
  )

# Plot the correlations
ggplot(correlation_data, aes(x = corr_epithelial, y = corr_mesenchymal)) +
  geom_point(size = 3, color = "blue") +
  geom_text_repel(aes(label = CancerType), size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Correlations of ELF3 with Epithelial and Mesenchymal Scores Across Cancer Types",
       x = "Correlation with Epithelial Scores",
       y = "Correlation with Mesenchymal Scores") +
  theme_minimal()

# Boxplot for ELF3 expression
ggplot(elf3_expr_data, aes(x = reorder(CancerType, epithelial_scores), y = ELF3_expr)) +
  geom_boxplot() +
  coord_flip() +
  labs(title = "ELF3 Expression Across Cancer Types Ordered by Median Epithelial Scores",
       x = "Cancer Type",
       y = "ELF3 Expression") +
  theme_minimal()
