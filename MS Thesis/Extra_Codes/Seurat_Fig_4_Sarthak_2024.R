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
seurat_gse176078 <- NormalizeData(seurat_gse176078)
seurat_gse176078 <- ScaleData(seurat_gse176078)

seurat_gse173634 <- NormalizeData(seurat_gse173634)
seurat_gse173634 <- ScaleData(seurat_gse173634)

####################
library(Seurat)
library(AnnotationHub)
library(ensembldb)
library(dplyr)

# Extract Ensembl IDs
ensembl_ids <- rownames(seurat_gse173634[["RNA"]])

# Create an AnnotationHub object
ah <- AnnotationHub()

# Query for the appropriate ensembledb package based on the species
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

# Function to create panel A (i): Density plot of VIM - CDH1
plot_panel_A_i <- function(seurat_obj, assay = "RNA") {
    
    if(!("VIM" %in% rownames(seurat_obj[[assay]]))) {
        print("VIM gene not found in the dataset, please make sure the gene name is correct")
        return(NULL)
    }
    if(!("CDH1" %in% rownames(seurat_obj[[assay]]))) {
        print("CDH1 gene not found in the dataset, please make sure the gene name is correct")
        return(NULL)
    }
    
    vim_cdh1_diff <- GetAssayData(seurat_obj, assay = assay, layer = "data")["VIM",] -
        GetAssayData(seurat_obj, assay = assay, layer = "data")["CDH1",]
    
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

# Function to classify cells into epithelial, hybrid, and mesenchymal states
classify_cells <- function(seurat_obj, assay = "RNA") {
    if(!("VIM" %in% rownames(seurat_obj[[assay]]))) {
        print("VIM gene not found in the dataset, please make sure the gene name is correct")
        return(NULL)
    }
    if(!("CDH1" %in% rownames(seurat_obj[[assay]]))) {
        print("CDH1 gene not found in the dataset, please make sure the gene name is correct")
        return(NULL)
    }
    
    vim_cdh1_diff <- GetAssayData(seurat_obj, assay = assay, layer="data")["VIM",] -
        GetAssayData(seurat_obj, assay = assay, layer="data")["CDH1",]
    
    cell_states <- cut(vim_cdh1_diff, breaks = c(-Inf, 0, 5, Inf), labels = c("Epithelial", "Hybrid", "Mesenchymal"))
    
    seurat_obj$cell_state <- cell_states
      gc()
    return(seurat_obj)
}

# Function to create panel A (ii) and (iii): Stacked bar plots
plot_panel_A_ii_iii <- function(seurat_obj, grouping_var, cell_state_var="cell_state") {
    
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
        summarize(count = n()) %>%
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

# Function to create panel B: Gene-gene correlation heatmap
plot_panel_B <- function(seurat_obj, top_n = 20, assay = "RNA") {
    
    # Define the gene sets
    genes_luminal <- c("FOXA1", "SPDEF", "XBP1", "DLX3","GATA3")
    genes_basal <- c("ZNF552", "TRPS1", "SOX13", "MYB","ZNF652","GRHL2", "BAZ2B")
    genes_epithelial <- c("POU2F3","BATF","DLX6","MAFF","ELK3","ARNTL2")
    genes_mesenchymal <- c("FOSL1", "ETS1", "ZEB2", "YBX1","ZNF280C","PBX4","HMGA1","ETS2","MYBL2","E2F2","MXD3","MYBL1","E2F7","TFDP1","FOXM1")
    all_genes <- c(genes_luminal, genes_basal, genes_epithelial, genes_mesenchymal)
    
    # Filter data to only include genes of interest
    if(!all(all_genes %in% rownames(seurat_obj[[assay]]))) {
        missing_genes <- all_genes[!all_genes %in% rownames(seurat_obj[[assay]])]
        print(paste("These genes are missing in the dataset:", paste(missing_genes, collapse = ","), ". Please check gene names and try again."))
        return(NULL)
    }
    
    filtered_data <- GetAssayData(seurat_obj, assay = assay, layer = "data")[all_genes,]
    
    # Calculate correlation matrix
     gc()
    correlation_matrix <- cor(t(as.matrix(filtered_data)))
    
  # Define groups for annotation
    group_labels <- factor(c(rep("Luminal", length(genes_luminal)),
                             rep("Basal", length(genes_basal)),
                             rep("Epithelial", length(genes_epithelial)),
                             rep("Mesenchymal", length(genes_mesenchymal))),
                           levels = c("Luminal", "Basal", "Epithelial", "Mesenchymal"))
    
    # Create annotation data
    annotation_col <- data.frame(Group = group_labels)
    rownames(annotation_col) <- colnames(correlation_matrix)
    
    
    pheatmap(correlation_matrix,
           annotation_col = annotation_col,
           color = colorRampPalette(c("navyblue", "white", "firebrick3"))(100),
           border_color = NA,
            labRow = rownames(correlation_matrix), labCol = rownames(correlation_matrix),
           margins = c(5,5))
    
     gc()
    return(NULL)
}

# Function to create panel D: Scatter plot of cells on 2D epithelial-mesenchymal and luminal-basal space
plot_panel_D <- function(seurat_obj, grouping_var, assay = "RNA") {
    
    if(is.null(seurat_obj[[grouping_var]])) {
        print("Grouping variable not found in the object metadata, please make sure it is correct")
        return(NULL)
    }
    if(is.null(seurat_obj$EpithelialScore1) | is.null(seurat_obj$MesenchymalScore1) | is.null(seurat_obj$LuminalScore1) | is.null(seurat_obj$BasalScore1)) {
        print("Please calculate scores first")
        return(NULL)
    }
    
    # Calculate differences
    seurat_obj$EM_diff <- seurat_obj$MesenchymalScore1 - seurat_obj$EpithelialScore1
    seurat_obj$LB_diff <- seurat_obj$BasalScore1 - seurat_obj$LuminalScore1
    
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

# Main function to orchestrate the analysis and plotting
create_plots <- function(seurat_obj, dataset_name, grouping_var, assay = "RNA") {
  # make sure that grouping_var is a string
  grouping_var <- as.character(grouping_var)

  # 0. Make row names unique and remove empty strings
  gene_names <- rownames(seurat_obj[["RNA"]])
  gene_names <- gene_names[gene_names != ""]
  unique_gene_names <- make.unique(gene_names)
  rownames(seurat_obj[["RNA"]]) <- unique_gene_names
  
  # 1. Calculate scores
   gc()
  seurat_obj <- calculate_scores(seurat_obj, assay=assay)

    # 2. Classify cells
  seurat_obj <- classify_cells(seurat_obj, assay = assay)
    

  # 3. Panel A (i)
  plot_A_i <- plot_panel_A_i(seurat_obj, assay = assay)
  
   # 4. Save panel A(i)
   png(filename = paste0(dataset_name, "_panel_A_i.png"), width = 6, height = 4, units = "in", res = 300)
   print(plot_A_i)
   dev.off()
  
  # 5. Panel A (ii and iii)
  plot_A_ii_iii <- plot_panel_A_ii_iii(seurat_obj, grouping_var)
  
   # 6. Save panel A (ii and iii)
    png(filename = paste0(dataset_name, "_panel_A_ii_iii.png"), width = 6, height = 4, units = "in", res = 300)
    print(plot_A_ii_iii)
    dev.off()
  # 7. Panel B
  plot_panel_B(seurat_obj, assay = assay)
    
  # 8. Save panel B
   dev.copy(png,filename=paste0(dataset_name, "_panel_B.png"), width=6, height=6, units="in", res = 300)
   dev.off()

  # 9. Panel D
  plot_D <- plot_panel_D(seurat_obj, grouping_var, assay = assay)
  
   # 10. Save panel D
    png(filename = paste0(dataset_name, "_panel_D.png"), width = 6, height = 6, units = "in", res = 300)
    print(plot_D)
    dev.off()
    
  gc()
  return(list(plot_A_i = plot_A_i, plot_A_ii_iii = plot_A_ii_iii, plot_D = plot_D))
}

# Make row names unique and remove empty names from the original object
gene_names <- rownames(seurat_gse173634[["RNA"]])
gene_names <- gene_names[gene_names != ""]
unique_gene_names <- make.unique(gene_names)
rownames(seurat_gse173634[["RNA"]]) <- unique_gene_names
seurat_gse173634 <- seurat_gse173634[rownames(seurat_gse173634) != "", ]


# Process and plot for GSE173634
plots_173634 <- create_plots(seurat_obj = seurat_gse173634, dataset_name = "GSE173634", grouping_var = "orig.ident", assay = "RNA")
