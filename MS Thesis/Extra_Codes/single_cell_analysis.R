library(BiocManager)
library(usethis)
library(devtools)
library(DrImpute)
library(Seurat)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tidyr)
library(Matrix)
library(tidyverse)
library(GEOquery)
library(AUCell)
library(GSEABase)

######################################################## 
Read10X_uncompressed <- function(file_path) {
  barcode.path <- file.path(file_path, "barcodes.tsv")
  features.path <- file.path(file_path, "genes.tsv")
  matrix.path <- file.path(file_path, "matrix.mtx")
  
  barcodes <- read.table(barcode.path, header = FALSE, stringsAsFactors = FALSE)
  features <- read.table(features.path, header = FALSE, stringsAsFactors = FALSE)
  matrix <- readMM(matrix.path)  

  colnames(matrix) <- barcodes$V1
  rownames(matrix) <- features$V1
  
  matrix <- as.matrix(matrix)
  
  imputed_matrix <- DrImpute::DrImpute(matrix)
  colnames(imputed_matrix) <- colnames(matrix)
  rownames(imputed_matrix) <- rownames(matrix)
 
  matrix <- as.matrix(imputed_matrix)
  print(paste0("The shape of exp_matrix:", dim(matrix)))
  return(matrix)
}

list_sample <-  list("BT474", "MDAMB361")
path_sample <- "./Transcriptomics/Sarthak_2024/GSE173634_RAW"

list_sample <- paste0(path_sample, list_sample)

list_exp_mat <- list()
list_seurat_obj <- list()

for (sample in list_sample) {
  sample_name <- basename(sample)
  exp_mat <- Read10X_uncompressed(sample)
  list_exp_mat[[sample_name]] <- exp_mat
  list_seurat_obj[[sample_name]] <- CreateSeuratObject(counts = exp_mat, project = sample_name, assay = "RNA", min.features = 200)
}

length(list_exp_mat)
length(list_seurat_obj)

all(rownames(list_exp_mat[[1]]) == rownames(list_exp_mat[[2]]) & rownames(list_exp_mat[1]) == rownames(list_exp_mat[3]))

total_cells <- sum(ncol((list_exp_mat[[1]])), ncol((list_exp_mat[[2]])))
total_genes <- sum(nrow((list_exp_mat[[1]])))

merged_seurat_object <- merge(list_seurat_obj[[1]], y = list_seurat_obj[2:length(list_seurat_obj)], add.cell.ids = names(list_seurat_obj))
glimpse(merged_seurat_object)

view(merged_seurat_object@meta.data)

options(future.globals.maxSize = 2 * 1024^3)

merged_seurat_object <- SCTransform(merged_seurat_object, verbose = FALSE)

Assays(merged_seurat_object)
Layers(merged_seurat_object)

################################################################## Plotting ####################################################################
exp_data <- GetAssayData(object = seurat_integrated, assay = "RNA", layer = "data")
exp_data[1:5, 1:5]

class(exp_data)
exp_data <- as.matrix(exp_data)

################################################################## FigA1 ###################################################################### 
data_VIM <- exp_data["VIM", ]
data_CDH1 <- exp_data["CDH1", ]

data_req = data_VIM - data_CDH1

plot_figA_i <- ggplot(data.frame(x = data_req), aes(x)) +
  geom_density(fill = "blue", alpha = 0.5) +
  geom_vline(xintercept = -0.9, color = "red", linetype = "dashed", linewidth = 1) +
  annotate("text", x = -0.5, y = 0.3, label = "x = -0.5", color = "black", size = 5, vjust = -1) + 
  geom_vline(xintercept = 1.4, color = "red", linetype = "dashed", linewidth = 1) +
  annotate("text", x = 1.4, y = 0.3, label = "x = 1.4", color = "black", size = 5, vjust = -1) + 
  labs(title = "KDE Plot", x = "VIM - CDH1", y = "Density") +
  theme_minimal()

ggsave("figA_i.png", plot = plot_figA_i, width = 6, height = 4, dpi = 300)

################################################################## plot A2 #####################################################################
data_req<- data.frame(data_req)
exp_data$cell_type <- ifelse(data_req < -0.5, 
                             "Epithelial", 
                             ifelse(data_req > 1.4, 
                                    "Mesenchymal", 
                                    "Hybrid"))

data_req$sample_name <- sapply(strsplit(row.names(data_req), "_"), function(x) x[1])

HER2_pos <- list("AU565", "EVSAT", "HCC1954", "JIMT1", "MDAMB453")
Basal_like <- list("MCF12A")
luminal_A <- list("BT483", "CAMA1", "EFM19", "HCC1500", "KPL1", "MCF7", "MDAMB415", "T47D", "ZR751")
luminal_B <- list("BT474", "MDAMB361")
TNBC_A <- list("BT20", "DU4475", "HCC1143", "HCC1187", "HCC1937", "HCC70", "MDAMB436", "MDAMB468")
TNBC_B <- list("BT549", "CAL51", "CAL851", "HCC38", "HDQP1", "HS578T", "MX1")

data_req$cancer_type <- sapply(data_req$sample_name, function(sample) {
  if (sample %in% HER2_pos) {
    return("HER2_pos")
  } else if (sample %in% Basal_like) {
    return("Basal_like")
  } else if (sample %in% luminal_A) {
    return("luminal_A")
  } else if (sample %in% luminal_B) {
    return("luminal_B")
  } else if (sample %in% TNBC_A) {
    return("TNBC_A")
  } else if (sample %in% TNBC_B) {
    return("TNBC_B")
  } else {
    return(NA)
  }
})

data_req$cell_type <- factor(data_req$cell_type, levels = c("Epithelial", "Hybrid", "Mesenchymal"))

view(data_req)

data_luminal_A <- subset(data_req, cancer_type == "luminal_A")
data_luminal_B <- subset(data_req, cancer_type == "luminal_B")
data_TNBC_A <- subset(data_req, cancer_type == "TNBC_A")
data_TNBC_B <- subset(data_req, cancer_type == "TNBC_B")
data_HER2_pos <- subset(data_req, cancer_type == "HER2_pos")
data_Basal_like <- subset(data_req, cancer_type == "Basal_like")

#calculate the percentage of each cell type in each cancer type
prcnt_data_luminal_A <- data_luminal_A %>%
  group_by(cell_type) %>%
  summarize(count = n()) %>%
  mutate(percentage = count / sum(count) * 100) 

prcnt_data_luminal_B <- data_luminal_B %>%
  group_by(cell_type) %>%
  summarize(count = n()) %>%
  mutate(percentage = count / sum(count) * 100) 

prcnt_data_TNBC_A <- data_TNBC_A %>%
  group_by(cell_type) %>%
  summarize(count = n()) %>%
  mutate(percentage = count / sum(count) * 100) 

prcnt_data_TNBC_B <- data_TNBC_B %>%
  group_by(cell_type) %>%
  summarize(count = n()) %>%
  mutate(percentage = count / sum(count) * 100) 

prcnt_data_HER2_pos <- data_HER2_pos %>%
  group_by(cell_type) %>%
  summarize(count = n()) %>%
  mutate(percentage = count / sum(count) * 100) 

prcnt_data_Basal_like <- data_Basal_like %>%
  group_by(cell_type) %>%
  summarize(count = n()) %>%
  mutate(percentage = count / sum(count) * 100) 


combined_prcnt_data <- bind_rows(
  prcnt_data_luminal_A %>% mutate(cancer_type = "luminal_A"),
  prcnt_data_luminal_B %>% mutate(cancer_type = "luminal_B"),
  prcnt_data_TNBC_A %>% mutate(cancer_type = "TNBC_A"),
  prcnt_data_TNBC_B %>% mutate(cancer_type = "TNBC_B"),
  prcnt_data_HER2_pos %>% mutate(cancer_type = "HER2_pos"),
  prcnt_data_Basal_like %>% mutate(cancer_type = "Basal_like")
)

ggplot(combined_prcnt_data, aes(x = cancer_type, y = percentage, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Percentage Composition of Cell Types in Different Cancer Types", x = "Cancer Type", y = "Percentage") +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  
ggsave("figA_ii.png", width = 6, height = 4, dpi = 300)

############################################################## plot D1 ###########################################################################


gene_set_path <- "./data/gene_sets/gene_sets_related_RPMS_BPMS.gmt"
gene_sets <- load_AUCell_genesets(gene_set_path)
gene_set_file <- read_tsv(gene_set_path, col_names = FALSE)
view(gene_set_file)

gmtFile <- "./data/gene_sets/gene_sets_related_RPMS_BPMS.gmt"
geneSets <- getGmt(gmtFile)
class(geneSets)


list_gene_sets <- list()
for (i in 1:length(geneSets)) {
  gene_set <- geneSets[[i]]
  
  print(gene_set@setName)
  gs <- gene_set@geneIds
  gs <- gs[gs != "" & !is.na(gs)]
  list_gene_sets[[gene_set@setName]] <- gs
  print(length(gs))
}

class(list_gene_sets)
str(list_gene_sets)
names(list_gene_sets)

# Calculate the AUCell scores for each cell in the Seurat object
gene_sets <- list_gene_sets

exp_mat <- GetAssayData(object = merged_seurat_object, assay = "SCT", layer = "data")
exp_mat <- as(exp_mat, "dgCMatrix")

cells_ranking <- AUCell_buildRankings(exp_mat,  plotStats = FALSE)
cells_AUC <- AUCell_calcAUC(gene_sets, cells_ranking)

class(cells_AUC)
glimpse(cells_AUC)
str(cells_AUC)

AUC_score_table <- cells_AUC@assays@data$AUC
AUC_score_table <- data.frame(AUC_score_table)
#sort the data frame alphabetically by the cell names
AUC_score_table <- AUC_score_table[, order(colnames(AUC_score_table))]
view(AUC_score_table)

exp_mat <- exp_mat[, order(colnames(exp_mat))]
RKIP_data <- as.numeric(exp_mat["PEBP1", ])
view(RKIP_data)
BACH1_data <- as.numeric(exp_mat["BACH1", ])
view(BACH1_data)

write.csv(AUC_score_table, "AUC_score_table.csv")
write.csv(RKIP_data, "RKIP_data.csv")
write.csv(BACH1_data, "BACH1_data.csv")

corr_list <- c("FAO", "HALLMARK_GLYCOLYSIS", "HALLMARK_OxPhos", "HIF1_metagene", 
             "HIF1", "HALLMARK_HYPOXIA", "Epi_Tumor", "Mes_Tumor", "luminal", "Basal", "pEMT", 
             "EMT_down", "EMT_partial", "EMT_up", "RPMS", "BPMS", 
             "HALLMARK_EMT")

r_values <- data.frame(matrix(ncol = 17, nrow = 4))
p_values <- data.frame(matrix(ncol = 17, nrow = 4))
colnames(r_values) <- corr_list
rownames(r_values) <- c("RKIP_data", "BACH1_data", "RPMS", "BPMS")
colnames(p_values) <- corr_list
rownames(p_values) <- c("RKIP_data", "BACH1_data", "RPMS", "BPMS")

for (col in colnames(exp_mat)) {
  print(paste0("Calculating correlation for: ", col))
  
  RKIP <- cor.test(RKIP_data, as.numeric(AUC_score_table[col, ]), method = "spearman")
  r_values["RKIP_data", col] <- RKIP$estimate
  p_values["RKIP_data", col] <- RKIP$p.value
  print("RKIP_done")
  
  # Calculate Spearman correlation for BACH1_data
  BACH1 <- cor.test(BACH1_data, as.numeric(AUC_score_table[col, ]), method = "spearman")
  r_values["BACH1_data", col] <- BACH1$estimate
  p_values["BACH1_data", col] <- BACH1$p.value
  print("BACH1_done")
  
  # Calculate Spearman correlation for RPMS
  RPMS <- cor.test(AUC_score_table["RPMS", ], as.numeric(AUC_score_table[col, ]), method = "spearman")
  r_values["RPMS", col] <- RPMS$estimate
  p_values["RPMS", col] <- RPMS$p.value
  print("RPMS_done")
  
  # Calculate Spearman correlation for BPMS
  BPMS <- cor.test(AUC_score_table["BPMS", ], as.numeric(AUC_score_table[col, ]), method = "spearman")
  r_values["BPMS", col] <- BPMS$estimate
  p_values["BPMS", col] <- BPMS$p.value
  print("BPMS_done")
}

view(as.numeric(AUC_score_table["FAO", ]))
view(RKIP_data)
