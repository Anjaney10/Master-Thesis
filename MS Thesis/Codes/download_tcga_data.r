# Author: Anjaney J Pandey
# Code to download TCGA data for analysis

########### Downloading the TCGA data ###########
library("httr")
library("R.utils")

cancer_types <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")

# Download the files using regex for the source links

for (cancer in cancer_types) {
    url <- paste0("https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-", cancer, ".star_tpm.tsv.gz")
    destfile <- paste0("./Transcriptomics/TCGA_DATA/TCGA-", cancer, ".star_tpm.tsv.gz")
    GET(url, write_disk(destfile, overwrite = TRUE))
}

for (cancer in cancer_types) {
    url <- paste0("https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-", cancer, ".survival.tsv.gz")
    destfile <- paste0("./Transcriptomics/TCGA_DATA/survival_data/TCGA-", cancer, ".survival.tsv.gz")
    GET(url, write_disk(destfile, overwrite = TRUE))
}

# Unzip the TPM and survival data

for (cancer in cancer_types) {
    tpm_file <- paste0("./Transcriptomics/TCGA_DATA/TCGA-", cancer, ".star_tpm.tsv.gz")
    survival_file <- paste0("./Transcriptomics/TCGA_DATA/survival_data/TCGA-", cancer, ".survival.tsv.gz")

    gunzip(tpm_file, remove = TRUE)
    gunzip(survival_file, remove = TRUE)
}

# Download the probemap for the gene names
probemap_url <- "https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v36.annotation.gtf.gene.probemap"
probemap_destfile <- "./Transcriptomics/TCGA_DATA/gencode.v36.annotation.gtf.gene.probemap"
GET(probemap_url, write_disk(probemap_destfile, overwrite = TRUE))

probemap <- read.table(probemap_destfile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Replace the Ensembl_ID with the Gene Symbol so that it can be used for downstream analysis

replace_ensembl_with_gene_symbol <- function(tpm_file, probemap) {
    tpm_data <- read.table(tpm_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    tpm_data$Gene <- probemap$gene[match(tpm_data$Ensembl_ID, probemap$id)]
    tpm_data <- tpm_data[, c("Gene", colnames(tpm_data)[-1])]
    write.table(tpm_data, file = tpm_file, sep = "\t", row.names = FALSE, quote = FALSE)
}

for (cancer in cancer_types) {
    tpm_file <- paste0("./Transcriptomics/TCGA_DATA/TCGA-", cancer, ".star_tpm.tsv")
    replace_ensembl_with_gene_symbol(tpm_file, probemap)
}

# These are old files that were downloaded for analysis but later on we switched to TPM files. So this can be neglected.

# for (cancer in cancer_types) {
#     url <- paste0("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.", cancer, ".sampleMap%2FHiSeqV2.gz")
#     destfile <- paste0("./Transcriptomics/tcga/TCGA-", cancer, ".HiSeqV2.gz")
#     GET(url, write_disk(destfile, overwrite = TRUE))
# }
#
# for (cancer in cancer_types){
#     url <- paste0("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/survival%2F", cancer, "_survival.txt")
#     destfile <- paste0("./Transcriptomics/tcga/survival_data/", cancer, "_survival.txt")
#     GET(url, write_disk(destfile, overwrite = TRUE))
# }
