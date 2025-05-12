library(Seurat)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tidyr)
library(Matrix)

# Function to read uncompressed 10X data
Read10X_uncompressed <- function(data.dir) {
  barcode.path <- file.path(data.dir, "barcodes.tsv")
  features.path <- file.path(data.dir, "genes.tsv")
  matrix.path <- file.path(data.dir, "matrix.mtx")

  barcodes <- read.table(barcode.path, header = FALSE, stringsAsFactors = FALSE)
  features <- read.table(features.path, header = FALSE, stringsAsFactors = FALSE)
  matrix <- Matrix::readMM(matrix.path)

  colnames(matrix) <- barcodes$V1
  rownames(matrix) <- features$V1

  matrix <- as(matrix, "CsparseMatrix")

  return(matrix)
}

# # Create Seurat object for GSE176078
# data_dir_gse176078 <- "/home/csb/Anjaney/Transcriptomics/Sarthak_2024/GSE176078_RAW/"
# subfolders_gse176078 <- list.dirs(data_dir_gse176078, full.names = TRUE, recursive = FALSE)

# seurat_list_gse176078 <- lapply(subfolders_gse176078, function(folder) {
#   data <- Read10X_uncompressed(data.dir = folder)
#   CreateSeuratObject(counts = data)
# })

# seurat_gse176078 <- merge(seurat_list_gse176078[[1]], y = seurat_list_gse176078[-1], add.cell.ids = basename(subfolders_gse176078))

# # Create Seurat object for GSE173634
# data_dir_gse173634 <- "/home/csb/Anjaney/Transcriptomics/Sarthak_2024/GSE173634_RAW/"
# subfolders_gse173634 <- list.dirs(data_dir_gse173634, full.names = TRUE, recursive = FALSE)

# seurat_list_gse173634 <- lapply(subfolders_gse173634, function(folder) {
#   data <- Read10X_uncompressed(data.dir = folder)
#   CreateSeuratObject(counts = data)
# })

# seurat_gse173634 <- merge(seurat_list_gse173634[[1]], y = seurat_list_gse173634[-1], add.cell.ids = basename(subfolders_gse173634))

# Path to the directory containing your sample subfolders
data_dir_gse173634 <- "/home/csb/Anjaney/Transcriptomics/Sarthak_2024/GSE173634_RAW/"
sample_subfolders <- list.dirs(data_dir_gse173634, full.names = TRUE, recursive = FALSE)

# Define the metadata mapping (you can also read this from a file)
metadata_mapping <- data.frame(
  sample = c("AU565", "BT20", "BT474", "BT483", "BT549", "CAL51", "CAL851", "CAMA1", "DU4475", "EFM19", "EVSAT", "HCC1143", "HCC1187", "HCC1500", "HCC1937", "HCC1954", "HCC38", "HCC70", "HDQP1", "HS578T", "JIMT1", "KPL1", "MCF12A", "MCF7_1", "MDAMB361", "MDAMB415", "MDAMB436_1", "MDAMB453", "MDAMB468", "MX1", "T47D", "ZR751", "MCF7_2", "MDAMB436_2"),
  cell_line = c("AU565", "BT20", "BT474", "BT483", "BT549", "CAL51", "CAL851", "CAMA1", "DU4475", "EFM19", "EVSAT", "HCC1143", "HCC1187", "HCC1500", "HCC1937", "HCC1954", "HCC38", "HCC70", "HDQP1", "HS578T", "JIMT1", "KPL1", "MCF12A", "MCF7", "MDAMB361", "MDAMB415", "MDAMB436", "MDAMB453", "MDAMB468", "MX1", "T47D", "ZR751", "MCF7", "MDAMB436"),
  subtype = c("H", "TNA", "LB", "LA", "TNB", "TNB", "TNB", "LA", "TNA", "LA", "H", "TNA", "TNA", "LA", "TNA", "H", "TNB", "TNA", "TNB", "TNB", "H", "LA", "Basal-like", "LA", "LB", "LA", "TNA", "H", "TNA", "TNB", "LA", "LA", "LA", "TNA")
)
# Create a list to store Seurat objects
seurat_objects_list <- list()

# Loop through each sample subfolder
for (folder in sample_subfolders) {
  # Read in the count matrix
  data <- Read10X_uncompressed(data.dir = folder)

  # Extract sample name from folder name
  sample_name <- basename(folder)

  # Create Seurat object without metadata
  seurat_obj <- CreateSeuratObject(counts = data, project = sample_name)

  # Store the Seurat object in the list
  seurat_objects_list[[sample_name]] <- seurat_obj
}

# Merge the Seurat objects
seurat_gse173634 <- merge(seurat_objects_list[[1]], 
                          y = seurat_objects_list[-1],
                          add.cell.ids = names(seurat_objects_list),
                          project = "GSE173634")

# Rename layers (assuming only "counts" and "data" exist before merging)
#current.layers <- names(seurat_gse173634@assays[[DefaultAssay(seurat_gse173634)]]@layers)
#counts.layers <- grep("^counts", current.layers, value = TRUE)

# Rename only if the layer names follow the pattern identified
#if (length(counts.layers) > 0) {
#  names(seurat_gse173634@assays[[DefaultAssay(seurat_gse173634)]]@layers)[match(counts.layers, names(seurat_gse173634@assays[[DefaultAssay(seurat_gse173634)]]@layers))] <- "counts"
#}

# Add metadata after merging
# Extract sample names from barcodes in the merged object
barcode_names <- colnames(seurat_gse173634)
sample_names_from_barcodes <- sub("_.+", "", barcode_names)

# Create metadata for the merged object
merged_metadata <- data.frame(
  row.names = barcode_names, # Set row names as barcodes
  barcode = barcode_names,
  sample = sample_names_from_barcodes
)

# Merge with the metadata mapping
merged_metadata <- merge(merged_metadata, metadata_mapping, by = "sample", all.x = TRUE)

# Ensure row names are preserved
rownames(merged_metadata) <- merged_metadata$barcode

# Reorder metadata rows to match the order of cells in the Seurat object
merged_metadata <- merged_metadata[colnames(seurat_gse173634), ]

# Add metadata to the merged Seurat object
seurat_gse173634 <- AddMetaData(seurat_gse173634, metadata = merged_metadata)

# Normalize and scale the data for both datasets
# seurat_gse176078 <- NormalizeData(seurat_gse176078)
# seurat_gse176078 <- ScaleData(seurat_gse176078)

seurat_gse173634 <- NormalizeData(seurat_gse173634)
seurat_gse173634 <- ScaleData(object = seurat_gse173634, features = rownames(seurat_gse173634), 
                              vars.to.regress = "nCount_RNA", assay = "RNA",
                              block.size = 100, verbose =  TRUE)
####################
####################
# Compute Scores (Corrected for Seurat v5)
luminal_genes <- c("ESR1", "GATA3", "PGR", "FOXA1")
basal_genes <- c("TP63", "SNAI1")
epithelial_genes <- c("CDH1", "MIR200C")
mesenchymal_genes <- c("ZEB1", "SNAI1")

# Function to calculate scores (Corrected for Seurat v5)
calculate_scores_manual <- function(seurat_obj, assay = "RNA") {
  luminal_genes <- c("ESR1", "GATA3", "PGR", "FOXA1")
  basal_genes <- c("TP63", "SNAI1")
  epithelial_genes <- c("CDH1", "MIR200C")
  mesenchymal_genes <- c("ZEB1", "SNAI1")
  
  # Get the normalized data (data layer)
  normalized_data <- LayerData(seurat_obj, assay = assay, layer = "data")
  
  # Calculate module scores manually
  seurat_obj$Luminal1 <- colMeans(normalized_data[rownames(normalized_data) %in% luminal_genes, , drop = FALSE])
  seurat_obj$Basal1 <- colMeans(normalized_data[rownames(normalized_data) %in% basal_genes, , drop = FALSE])
  seurat_obj$Epithelial1 <- colMeans(normalized_data[rownames(normalized_data) %in% epithelial_genes, , drop = FALSE])
  seurat_obj$Mesenchymal1 <- colMeans(normalized_data[rownames(normalized_data) %in% mesenchymal_genes, , drop = FALSE])
  
  return(seurat_obj)
}

# Function to create panel A (i): Density plot of VIM - CDH1 (Corrected for Seurat v5)
plot_panel_A_i <- function(seurat_obj, assay = "RNA") {
  if(!("VIM" %in% rownames(seurat_obj[[assay]]))) {
    print("VIM gene not found in the dataset, please make sure the gene name is correct")
    return(NULL)
  }
  if(!("CDH1" %in% rownames(seurat_obj[[assay]]))) {
    print("CDH1 gene not found in the dataset, please make sure the gene name is correct")
    return(NULL)
  }
  
  vim_cdh1_diff <- LayerData(seurat_obj, assay = assay, layer = "scale.data")["VIM",] -
    LayerData(seurat_obj, assay = assay, layer = "scale.data")["CDH1",]
  
  df <- data.frame(VIM_CDH1 = vim_cdh1_diff)
  gc()
  p <- ggplot(df, aes(x = VIM_CDH1)) +
    geom_density(fill = "gray", alpha = 0.7) +
    geom_vline(xintercept = c(0, 5), color = "red") +
    annotate("text", x = -2, y = 0.3, label = "Epithelial") +
    annotate("text", x = 2.5, y = 0.15, label = "Hybrid") +
    annotate("text", x = 7, y = 0.3, label = "Mesenchymal") +
    xlab("VIM - CDH1") +
    ylab("Density") +
    theme_minimal()
  return(p)
}

# Function to classify cells into epithelial, hybrid, and mesenchymal states (Corrected for Seurat v5)
classify_cells <- function(seurat_obj, assay = "RNA") {
  if(!("VIM" %in% rownames(seurat_obj[[assay]]))) {
    print("VIM gene not found in the dataset, please make sure the gene name is correct")
    return(NULL)
  }
  if(!("CDH1" %in% rownames(seurat_obj[[assay]]))) {
    print("CDH1 gene not found in the dataset, please make sure the gene name is correct")
    return(NULL)
  }
  
  vim_cdh1_diff <- LayerData(seurat_obj, assay = assay, layer = "scale.data")["VIM",] -
    LayerData(seurat_obj, assay = assay, layer = "scale.data")["CDH1",]
  
  cell_states <- cut(vim_cdh1_diff, breaks = c(-Inf, 0, 5, Inf), labels = c("Epithelial", "Hybrid", "Mesenchymal"))
  
  seurat_obj$cell_state <- cell_states
  gc()
  return(seurat_obj)
}

# Function to create panel A (ii) and (iii): Stacked bar plots (Corrected for Seurat v5)
plot_panel_A_ii_iii <- function(seurat_obj, grouping_var, cell_state_var = "cell_state") {
  if(is.null(seurat_obj[[grouping_var]])) {
    print("Grouping variable not found in the object metadata, please make sure it is correct")
    return(NULL)
  }
  if(is.null(seurat_obj[[cell_state_var]])) {
    print("Cell state variable not found in the object metadata, please classify cells first")
    return(NULL)
  }
  
  df <- seurat_obj@meta.data %>%
    group_by(!!sym(grouping_var), !!sym(cell_state_var)) %>%
    summarize(count = n(), .groups = 'drop') %>%
    group_by(!!sym(grouping_var)) %>%
    mutate(fraction = count / sum(count))
  
  plot <- ggplot(df, aes(x = !!sym(grouping_var), y = fraction, fill = !!sym(cell_state_var))) +
    geom_bar(stat = "identity", position = "stack") +
    ylab("Fraction of cases") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = c("Epithelial" = "gray", "Hybrid" = "orange", "Mesenchymal" = "red"))
  gc()
  return(plot)
}

# Function to create panel B: Gene-gene correlation heatmap (Corrected for Seurat v5)
plot_panel_B <- function(seurat_obj, top_n = 20, assay = "RNA") {
  # Define the gene sets based on figure
  genes_luminal <- c("FOXA1", "SPDEF", "XBP1", "DLX3", "GATA3")
  genes_basal <- c("ZNF552", "TRPS1", "SOX13", "MYB", "ZNF652", "GRHL2", "BAZ2B")
  genes_epithelial <- c("POU2F3", "BATF", "DLX6", "MAFF", "ELK3", "ARNTL2")
  genes_mesenchymal <- c("FOSL1", "ETS1", "ZEB2", "YBX1", "ZNF280C", "PBX4", "HMGA1", "ETS2", "MYBL2", "E2F2", "MXD3", "MYBL1", "E2F7", "TFDP1", "FOXM1")
  all_genes <- c(genes_luminal, genes_basal, genes_epithelial, genes_mesenchymal)
  
  # Filter data to only include genes of interest
  if(!all(all_genes %in% rownames(seurat_obj[[assay]]))) {
    missing_genes <- all_genes[!all_genes %in% rownames(seurat_obj[[assay]])]
    print(paste("These genes are missing in the dataset:", paste(missing_genes, collapse = ","), ". Please check gene names and try again."))
    return(NULL)
  }
  
  filtered_data <- LayerData(seurat_obj, assay = assay, layer = "scale.data")[all_genes, , drop = FALSE]
  
  # Check for genes with zero variance
  gene_variances <- apply(filtered_data, 1, var)
  if (any(gene_variances == 0, na.rm = TRUE)) {
    zero_var_genes <- names(gene_variances)[gene_variances == 0]
    print(paste("Warning: Genes with zero variance found and removed:", paste(zero_var_genes, collapse = ",")))
    filtered_data <- filtered_data[gene_variances > 0, , drop = FALSE] # Ensure this remains a matrix
  }
  
  # Check for NA values and handle them
  if (any(is.na(filtered_data))) {
    print("Warning: NA values found in filtered_data. Replacing with 0.")
    filtered_data[is.na(filtered_data)] <- 0 # Replace NA with 0 or handle as appropriate
  }
  
  # Check if there are enough genes for correlation calculation
  if(nrow(filtered_data) < 2) {
    print("Error: Not enough genes with non-zero variance for correlation analysis.")
    return(NULL)
  }
  
  # Calculate correlation matrix
  correlation_matrix <- cor(t(filtered_data))
  
  # Check if correlation matrix is valid for pheatmap
  if(nrow(correlation_matrix) < 2 || ncol(correlation_matrix) < 2) {
    print("Error: Correlation matrix has dimensions less than 2x2, which is not suitable for pheatmap.")
    return(NULL)
  }
  
  # Define groups for annotation
  group_labels <- factor(c(rep("Luminal", length(genes_luminal)),
                           rep("Basal", length(genes_basal)),
                           rep("Epithelial", length(genes_epithelial)),
                           rep("Mesenchymal", length(genes_mesenchymal))),
                         levels = c("Luminal", "Basal", "Epithelial", "Mesenchymal"))
  
  # Create annotation data
  annotation_col <- data.frame(Group = group_labels)
  rownames(annotation_col) <- colnames(correlation_matrix)
  
  # Generate heatmap
  pheatmap(correlation_matrix,
           annotation_col = annotation_col,
           color = colorRampPalette(c("navyblue", "white", "firebrick3"))(100),
           border_color = NA,
           labRow = rownames(correlation_matrix),
           labCol = rownames(correlation_matrix),
           margins = c(5,5))
  
  gc()
  return(NULL)
}

# Function to create panel D: Scatter plot of cells on 2D epithelial-mesenchymal and luminal-basal space (Corrected for Seurat v5)
plot_panel_D <- function(seurat_obj, grouping_var, assay = "RNA") {
  if(is.null(seurat_obj[[grouping_var]])) {
    print("Grouping variable not found in the object metadata, please make sure it is correct")
    return(NULL)
  }
  if(!all(c("Luminal1", "Basal1", "Epithelial1", "Mesenchymal1") %in% colnames(seurat_obj@meta.data))) {
    print("Not all score columns (Luminal1, Basal1, Epithelial1, Mesenchymal1) are present in the metadata. Skipping Panel D.")
    return(NULL)
  }
  
  # Calculate differences
  seurat_obj$EM_diff <- seurat_obj$Mesenchymal1 - seurat_obj$Epithelial1
  seurat_obj$LB_diff <- seurat_obj$Basal1 - seurat_obj$Luminal1
  
  # Create plot
  tryCatch({
    p <- ggplot(seurat_obj@meta.data, aes(x = EM_diff, y = LB_diff, color = !!sym(grouping_var))) +
      geom_point(size = 1, alpha = 0.7) +
      xlab("Mesenchymal - Epithelial") +
      ylab("Basal - Luminal") +
      theme_minimal() +
      theme(aspect.ratio = 1)
    gc()
    return(p)
  }, error = function(e){
    print(paste("Error in Panel D plotting:", e$message))
    return(NULL)
  })
}

# Main function to orchestrate the analysis and plotting (Corrected for Seurat v5)
create_plots <- function(seurat_obj, dataset_name, grouping_var, assay = "RNA") {
  # Make sure that grouping_var is a string
  grouping_var <- as.character(grouping_var)
  
  # Ensure necessary columns are present before proceeding
  if (!all(c("Mesenchymal1", "Epithelial1", "Basal1", "Luminal1") %in% colnames(seurat_obj@meta.data))) {
    print("Required module scores not found in metadata. Running calculate_scores to add them.")
    seurat_obj <- calculate_scores(seurat_obj, assay = assay)
  }
  
  # Classify cells - ensure this happens after score calculation
  seurat_obj <- classify_cells(seurat_obj, assay = assay)
  
  # Panel A (i)
  plot_A_i <- plot_panel_A_i(seurat_obj, assay = assay)
  if (!is.null(plot_A_i)) {
    png(filename = paste0(dataset_name, "_panel_A_i.png"), width = 6, height = 4, units = "in", res = 300)
    print(plot_A_i)
    dev.off()
  } else {
    print("Panel A (i) could not be generated due to missing data.")
  }
  
  # Panel A (ii and iii)
  plot_A_ii_iii <- plot_panel_A_ii_iii(seurat_obj, grouping_var)
  if (!is.null(plot_A_ii_iii)) {
    png(filename = paste0(dataset_name, "_panel_A_ii_iii.png"), width = 6, height = 4, units = "in", res = 300)
    print(plot_A_ii_iii)
    dev.off()
  } else {
    print("Panel A (ii and iii) could not be generated due to missing data.")
  }
  
  # Panel B
  # plot_panel_B(seurat_obj, assay = assay)
  # dev.copy(png, filename = paste0(dataset_name, "_panel_B.png"), width = 6, height = 6, units = "in", res = 300)
  # dev.off()
  
  # Panel D
  plot_D <- plot_panel_D(seurat_obj, grouping_var, assay = assay)
  if (!is.null(plot_D)) {
    png(filename = paste0(dataset_name, "_panel_D.png"), width = 6, height = 6, units = "in", res = 300)
    print(plot_D)
    dev.off()
  } else {
    print("Panel D could not be generated due to missing data.")
  }
  
  gc()
  return(list(plot_A_i = plot_A_i, plot_A_ii_iii = plot_A_ii_iii, plot_D = plot_D))
}

# Check if required metadata columns exist, if not, calculate scores
if (!all(c("Mesenchymal1", "Epithelial1", "Basal1", "Luminal1") %in% colnames(seurat_gse173634@meta.data))) {
  seurat_gse173634 <- calculate_scores_manual(seurat_gse173634)
}

plots_173634 <- create_plots(seurat_obj = seurat_gse173634, dataset_name = "GSE173634", grouping_var = "orig.ident", assay = "RNA")
