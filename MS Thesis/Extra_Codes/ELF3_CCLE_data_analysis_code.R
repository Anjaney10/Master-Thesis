# Load necessary libraries

library(R.utils)
# library(GEOquery)
library(GSVA)
library(ggplot2)
library(dplyr)
library(readxl)
library(biomaRt)
#################################################

# Read the RNA-seq file
ccle_data <- read.table("CCLE_RNAseq_rsem_genes_tpm_20180929.txt", header = TRUE, sep = "\t")

# Remove the version suffix from Ensembl IDs in ccle_data
ccle_data$gene_id <- sub("\\..*", "", ccle_data$gene_id)

# Set up the ensembl dataset
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Get the mapping from Ensembl IDs to gene symbols
mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ccle_data$gene_id,
  mart = ensembl
)

# Merge the mapping with the CCLE data
ccle_data <- merge(ccle_data, mapping, by.x = "gene_id", by.y = "ensembl_gene_id")

# Drop rows where there is no corresponding gene symbol
ccle_data <- ccle_data[ccle_data$hgnc_symbol != "", ]

# Remove rows with duplicated gene symbols
ccle_data_unique <- ccle_data[!duplicated(ccle_data$hgnc_symbol), ]

# Set the unique gene symbols as row names
rownames(ccle_data_unique) <- ccle_data_unique$hgnc_symbol

# Now remove the 'hgnc_symbol' column
ccle_data_unique <- ccle_data_unique[, -which(colnames(ccle_data_unique) == "hgnc_symbol")]

# Removing gene_id, transcript_ids columns

ccle_expr_data <- ccle_data_unique[, -c(1, 2)]


# Ensure your expression data is in matrix format with genes as rows
ccle_expr_matrix <- as.matrix(ccle_expr_data)

# Load the gene signature file
gene_signature <- read_excel("EM_gene_signature_tumor_KS.xlsx")

# Split genes based on Epi/Mes category
epithelial_genes <- gene_signature$Gene[gene_signature$Type == "Epi"]
mesenchymal_genes <- gene_signature$Gene[gene_signature$Type == "Mes"]

# Create a gene set list for ssGSEA
gene_sets <- list(
  epithelial = epithelial_genes,
  mesenchymal = mesenchymal_genes
)


# Perform ssGSEA
ssgsea_param <- ssgseaParam(
  ccle_expr_matrix, 
  gene_sets
)

ssgsea_scores <- gsva(ssgsea_param)

# Transpose the expression matrix for correlation calculation
ccle_expr_t <- t(ccle_expr_data)

# Calculate correlations with epithelial and mesenchymal scores
cor_results <- data.frame(
  gene = rownames(ccle_expr_data),
  cor_epithelial = apply(ccle_expr_t, 2, function(x) cor(x, ssgsea_scores["epithelial", ], method = "pearson")),
  cor_mesenchymal = apply(ccle_expr_t, 2, function(x) cor(x, ssgsea_scores["mesenchymal", ], method = "pearson"))
)

# List of genes to highlight
highlight_genes <- c("VIM", "ZEB1", "SNAI1", "SNAI2", "GRHL2", "OVOL2", "KLF4", "CDH1", "ELF3")

# Assign colors for the specific gene categories
cor_results <- cor_results %>%
  mutate(
    category = case_when(
      gene %in% mesenchymal_genes ~ "Mesenchymal",
      gene %in% epithelial_genes ~ "Epithelial",
      gene == "ELF3" ~ "ELF3",
      TRUE ~ "Other"
    ),
    color = case_when(
      category == "Mesenchymal" ~ "blue",
      category == "Epithelial" ~ "orange",
      category == "ELF3" ~ "red",
      TRUE ~ "grey"
    )
  )


# Create the scatter plot
library(ggplot2)
library(ggrepel)

# Create the plot
plot <- ggplot(cor_results, aes(x = cor_epithelial, y = cor_mesenchymal)) +
  geom_point(aes(color = ifelse(gene %in% mesenchymal_genes, "Mesenchymal",
                                ifelse(gene %in% epithelial_genes, "Epithelial",
                                       ifelse(gene == "ELF3", "ELF3", "Other")))),
             shape = 19, size = 1) +
  scale_color_manual(values = c("Mesenchymal" = "blue",
                                "Epithelial" = "orange",
                                "ELF3" = "red",
                                "Other" = "gray")) +
  xlim(-1, 1) +  # Set x-axis limits
  ylim(-1, 1) +  # Set y-axis limits
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Horizontal line
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # Vertical line
  labs(title = "Epithelial vs. Mesenchymal Scores",
       x = "Epithelial Score",
       y = "Mesenchymal Score") +
  theme_minimal() +
  theme(legend.position = "none") +  # Remove the legend
  geom_text_repel(aes(label = ifelse(gene %in% c("CDH1", "GRHL2", "VIM", "ZEB1", "SNAI1", "SNAI2", "ELF3", "OVOL2", "KLF4"), gene, "")),
                  box.padding = 0.5,   # Adjust padding
                  point.padding = 0.5,  # Padding around points
                  max.overlaps = Inf)   # Allow for all overlaps

# Display the plot
print(plot)

######################################################

# List of project directories (each containing the extracted files)
project_dirs <- list.dirs("~/Anjaney/Transcriptomics/TCGA_Data", recursive = FALSE)

# Initialize an empty dataframe to store combined data
combined_data <- data.frame()

# Load necessary libraries
library(dplyr)

# Initialize an empty dataframe to store combined data
combined_data <- data.frame()

# Loop through each project directory to read and combine data
for (proj_dir in project_dirs) {
  
  # Path to the expression file for this project
  expr_file <- file.path(proj_dir, "data_mrna_seq_tpm.txt")
  
  # Path to the meta study file for this project
  meta_file <- file.path(proj_dir, "meta_study.txt")
  
  # Check if both files exist
  if (file.exists(expr_file) & file.exists(meta_file)) {
    
    # Read the expression file
    expr_data <- read.delim(expr_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    # Read the meta study file
    meta_data <- readLines(meta_file)
    
    # Extract the short_name from meta_data
    short_name_line <- meta_data[grepl("short_name:", meta_data)]
    short_name <- gsub("short_name: ", "", short_name_line)  # Extract short_name by removing the label
    
    # Create a row with the short_name as the first column (and NA for all sample columns)
    short_name_row <- c(short_name, rep(NA, ncol(expr_data) - 1))
    
    # Convert the short_name_row to a data frame with the same structure as expr_data
    short_name_df <- as.data.frame(t(short_name_row), stringsAsFactors = FALSE)
    colnames(short_name_df) <- colnames(expr_data)
    
    # Combine the short_name row with the expression data
    expr_data_with_short_name <- rbind(short_name_df, expr_data)
    
    # Combine the data: assuming the first column is gene names and others are samples
    if (nrow(combined_data) == 0) {
      combined_data <- expr_data_with_short_name  # If first time, just assign it
    } else {
      combined_data <- bind_rows(combined_data, expr_data_with_short_name)
    }
  }
}

# Now `combined_data` contains the combined expression data with an additional row containing 'short_name' before each project's samples.

# View combined data structure
#str(combined_data)

############################################

# Filter for ELF3 gene
elf3_expression <- combined_data %>% filter(Entrez_Gene_Id == "1999") #corresponding to ELF3 gene

################################################

library(data.table)
library(fst)

# Set up an empty file to store results
combined_file <- "combined_expression_data.fst"

# Initialize variables
first_run <- TRUE  # To track if it's the first file
output_file <- tempfile("combined_data_chunk_", fileext = ".fst")  # Temp file for combined data

# Helper function to read a file in chunks
read_data_in_chunks <- function(file_path, chunk_size, header = TRUE) {
  con <- file(file_path, "r")
  # Read the header
  if (header) {
    header_line <- readLines(con, n = 1)
    col_names <- strsplit(header_line, "\t")[[1]]
  }
  
  data_chunks <- list()
  
  repeat {
    chunk <- read.table(con, sep = "\t", nrows = chunk_size, header = FALSE, fill = TRUE, stringsAsFactors = FALSE)
    
    if (nrow(chunk) == 0) {
      break
    }
    
    data_chunks[[length(data_chunks) + 1]] <- chunk
  }
  
  close(con)
  
  # Combine chunks and assign column names
  combined_chunk <- do.call(rbind, data_chunks)
  if (header) {
    setnames(combined_chunk, col_names)
  }
  
  return(combined_chunk)
}

# Loop through each project directory to read and combine data
for (proj_dir in project_dirs) {
  
  # Path to the expression file for this project
  expr_file <- file.path(proj_dir, "data_mrna_seq_tpm.txt")
  
  # Path to the meta study file for this project
  meta_file <- file.path(proj_dir, "meta_study.txt")
  
  # Check if both files exist
  if (file.exists(expr_file) & file.exists(meta_file)) {
    
    # Read the meta study file to get the short_name
    meta_data <- readLines(meta_file)
    short_name_line <- meta_data[grepl("short_name:", meta_data)]
    short_name <- gsub("short_name: ", "", short_name_line)  # Extract short_name
    
    # Read the expression data in chunks manually
    chunk_size <- 1000  # Set a chunk size (adjust based on memory limits)
    expr_data <- read_data_in_chunks(expr_file, chunk_size)
    
    # Add the short_name as a new column to the expression data
    expr_data$Study_Short_Name <- short_name
    
    # If it's the first run, create the output file with the first chunk
    if (first_run) {
      write_fst(expr_data, combined_file)  # Create the fst file with the first chunk
      first_run <- FALSE  # Mark that first run is done
    } else {
      # Append subsequent chunks to the combined data file
      append_fst(expr_data, combined_file)
    }
  }
}

# Load the final dataset into memory (if needed)
final_combined_data <- read_fst(combined_file, as.data.table = TRUE)


###############################

# Install necessary packages

library(data.table)
library(fst)

# Set up an empty file to store results
combined_file <- "combined_expression_data.fst"

# Initialize variables
first_run <- TRUE  # To track if it's the first file
output_file <- tempfile("combined_data_chunk_", fileext = ".fst")  # Temp file for combined data

# Helper function to read a file in chunks
read_data_in_chunks <- function(file_path, chunk_size, header = TRUE) {
  con <- file(file_path, "r")
  # Read the header
  if (header) {
    header_line <- readLines(con, n = 1)
    col_names <- strsplit(header_line, "\t")[[1]]
  }
  
  data_chunks <- list()
  
  repeat {
    chunk <- read.table(con, sep = "\t", nrows = chunk_size, header = FALSE, fill = TRUE, stringsAsFactors = FALSE)
    
    if (nrow(chunk) == 0) {
      break
    }
    
    data_chunks[[length(data_chunks) + 1]] <- chunk
  }
  
  close(con)
  
  # Combine chunks and assign column names
  combined_chunk <- do.call(rbind, data_chunks)
  if (header) {
    setnames(combined_chunk, col_names)
  }
  
  return(combined_chunk)
}

# Loop through each project directory to read and combine data
for (proj_dir in project_dirs) {
  
  # Extract the project identifier from the parent folder (e.g., 'ACC' from 'acc_tcga_gdc')
  project <- basename(proj_dir)
  project_name <- toupper(strsplit(project, "_")[[1]][1])  # Convert to uppercase
  
  # Path to the expression file for this project
  expr_file <- file.path(proj_dir, "data_mrna_seq_tpm.txt")
  
  # Check if the expression file exists
  if (file.exists(expr_file)) {
    
    # Read the expression data in chunks manually
    chunk_size <- 1000  # Set a chunk size (adjust based on memory limits)
    expr_data <- read_data_in_chunks(expr_file, chunk_size)
    
    # Add the project name as a new column to the expression data
    expr_data$Study_Short_Name <- project_name  # Insert project name (e.g., ACC, BRCA)
    
    # If it's the first run, create the output file with the first chunk
    if (first_run) {
      write_fst(expr_data, combined_file)  # Create the fst file with the first chunk
      first_run <- FALSE  # Mark that first run is done
    } else {
      # Append subsequent chunks to the combined data file
      append_fst(expr_data, combined_file)
    }
  }
}

# Load the final dataset into memory (if needed)
final_combined_data <- read_fst(combined_file, as.data.table = TRUE)

