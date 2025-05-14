library(Seurat)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tidyr)
library(Matrix)

# Function to read uncompressed 10X data
Read10X_uncompressed <- function(data.dir) {
  barcode.path <- file.path(data.dir, "count_matrix_barcodes.tsv")
  features.path <- file.path(data.dir, "count_matrix_genes.tsv")
  matrix.path <- file.path(data.dir, "count_matrix_sparse.mtx")
  
  barcodes <- read.table(barcode.path, header = FALSE, stringsAsFactors = FALSE)
  features <- read.table(features.path, header = FALSE, stringsAsFactors = FALSE)
  matrix <- Matrix::readMM(matrix.path)
  
  colnames(matrix) <- barcodes$V1
  rownames(matrix) <- features$V1
  
  matrix <- as(matrix, "CsparseMatrix")
  
  return(matrix)
}

# Create Seurat object for GSE176078
data_dir_gse176078 <- "./Transcriptomics/Sarthak_2024/GSE176078_RAW/"
subfolders_gse176078 <- list.dirs(data_dir_gse176078, full.names = TRUE, recursive = FALSE)

seurat_list_gse176078 <- lapply(subfolders_gse176078, function(folder) {
  data <- Read10X_uncompressed(data.dir = folder)
  CreateSeuratObject(counts = data)
})

seurat_gse176078 <- merge(seurat_list_gse176078[[1]], y = seurat_list_gse176078[-1], add.cell.ids = basename(subfolders_gse176078))

# Create Seurat object for GSE173634
data_dir_gse173634 <- "./Transcriptomics/Sarthak_2024/GSE173634_RAW/"
data_gse173634 <- Read10X_uncompressed(data.dir = data_dir_gse173634)
seurat_gse173634 <- CreateSeuratObject(counts = data_gse173634)

# Normalize and scale the data for both datasets
# seurat_gse176078 <- NormalizeData(seurat_gse176078)
# seurat_gse176078 <- ScaleData(seurat_gse176078)

# seurat_gse173634 <- NormalizeData(seurat_gse173634)
# seurat_gse173634 <- ScaleData(seurat_gse173634)

####################
library(Seurat)
library(AnnotationHub)
library(ensembldb)
library(dplyr)

# Extract Ensembl IDs
ensembl_ids <- rownames(seurat_gse173634[["RNA"]])

# Create an AnnotationHub object
ah <- AnnotationHub()

# Query for the appropriate ensembledb package based on the species (human in this case)
edb <- query(ah, pattern = c("EnsDb", "Homo sapiens", "105"))[[1]]

# Map ensembl IDs to gene symbols using the columns function on the EnsDb object
genes <- ensembldb::select(edb, keys = ensembl_ids,
                      keytype = "GENEID", columns = c("GENENAME","GENEID")) %>%
        dplyr::rename(ensembl_id = GENEID, gene_name = GENENAME)

# Join with the original object's feature name
df <- data.frame(ensembl_id = ensembl_ids)
df <- dplyr::left_join(df, genes, by="ensembl_id")

#replace row names with gene names
rownames(seurat_gse173634) <- df$gene_name


# Compute Scores
luminal_genes <- c("ERa66", "GATA3", "PGR", "FOXA1")
basal_genes <- c("DeltaNP63", "SLUG")
epithelial_genes <- c("CDH1", "miR-200")
mesenchymal_genes <- c("ZEB1", "SLUG")

seurat_gse173634 <- AddModuleScore(seurat_gse173634, features = list(luminal_genes), name = "Luminal")
seurat_gse173634 <- AddModuleScore(seurat_gse173634, features = list(basal_genes), name = "Basal")
seurat_gse173634 <- AddModuleScore(seurat_gse173634, features = list(epithelial_genes), name = "Epithelial")
seurat_gse173634 <- AddModuleScore(seurat_gse173634, features = list(mesenchymal_genes), name = "Mesenchymal")


# Panel A
## i) Kernel Density Plot
data_density <- data.frame(
  VIM_CD = seurat_gse173634$Mesenchymal1 - seurat_gse173634$Epithelial1
)

ggplot(data_density, aes(x = VIM_CD)) +
  geom_density(color = "black") +
  geom_vline(xintercept = c(-1, 1), color = "red") +
  theme_minimal() +
  labs(title = "Kernel Density of VIM-CDH1", x = "VIM-CDH1", y = "Density")

## ii & iii) Bar Plot for Subtype Composition
seurat_gse173634$State <- cut(
  seurat_gse173634$Mesenchymal1 - seurat_gse173634$Epithelial1,
  breaks = c(-Inf, -1, 1, Inf),
  labels = c("Epithelial", "Hybrid", "Mesenchymal")
)

cell_line_composition <- seurat_gse173634@meta.data %>%
  group_by(State, CellLine = orig.ident) %>%
  summarize(Fraction = n()) %>%
  group_by(CellLine) %>%
  mutate(Fraction = Fraction / sum(Fraction))

ggplot(cell_line_composition, aes(x = CellLine, y = Fraction, fill = State)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(title = "Subtype Composition by Cell Line", x = "Cell Line", y = "Fraction", fill = "State") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Panel B
## Gene-Gene Pairwise Correlation Heatmap
tf_genes <- rownames(seurat_gse173634)[1:20]
score_matrix <- seurat_gse173634@assays$RNA@scale.data[tf_genes, ]
cor_matrix <- cor(score_matrix)
pheatmap(cor_matrix, main = "TF Pairwise Correlation", clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean")

# Panel D
## i) Scatterplot for GSE173634
scatter_data_173634 <- data.frame(
  Mesenchymal_Epithelial = seurat_gse173634$Mesenchymal1 - seurat_gse173634$Epithelial1,
  Basal_Luminal = seurat_gse173634$Basal1 - seurat_gse173634$Luminal1,
  Subtype = seurat_gse173634$orig.ident
)

ggplot(scatter_data_173634, aes(x = Mesenchymal_Epithelial, y = Basal_Luminal, color = Subtype)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Scatterplot of EMT and Luminal-Basal Scores (GSE173634)", x = "Epithelial-Mesenchymal Score", y = "Basal-Luminal Score", color = "Subtype")

## ii) Scatterplot for GSE176078
seurat_gse176078 <- AddModuleScore(seurat_gse176078, features = list(luminal_genes), name = "Luminal")
seurat_gse176078 <- AddModuleScore(seurat_gse176078, features = list(basal_genes), name = "Basal")
seurat_gse176078 <- AddModuleScore(seurat_gse176078, features = list(epithelial_genes), name = "Epithelial")
seurat_gse176078 <- AddModuleScore(seurat_gse176078, features = list(mesenchymal_genes), name = "Mesenchymal")

scatter_data_176078 <- data.frame(
  Mesenchymal_Epithelial = seurat_gse176078$Mesenchymal1 - seurat_gse176078$Epithelial1,
  Basal_Luminal = seurat_gse176078$Basal1 - seurat_gse176078$Luminal1,
  PatientType = ifelse(seurat_gse176078$orig.ident == "TNBC", "TNBC", "ER+")
)

ggplot(scatter_data_176078, aes(x = Mesenchymal_Epithelial, y = Basal_Luminal, color = PatientType)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Scatterplot of EMT and Luminal-Basal Scores (GSE176078)", x = "Epithelial-Mesenchymal Score", y = "Basal-Luminal Score", color = "Patient Type")
