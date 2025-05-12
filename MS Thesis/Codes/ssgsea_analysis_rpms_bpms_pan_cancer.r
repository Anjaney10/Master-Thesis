# Author: Anjaney J Pandey
# Code to do ssGSEA analysis for RPMS and BPMS using GSVA package and storing the results for future analysis

################### ssGSEA analysis ####################
library("GSVA")

RPMS_genes <- scan(text = "PEBP1\tIGF2BP2\tGOLT1B\tNRAS\tRDX\tNF2\tDMD\tBACH1\tMAP4K4\tCOL5A2\tCDV3\tCOL3A1\tHMGA2\tCOL1A2\tTGFBR1\tPAPPA\tIGF2BP3\tARID3B\tLRIG2\tIRS2\tHIC2\tIL13\tGNPTAB\tKRT16\tTNFRSF12A\tENO1\tEN1\tTAGLN2\tBAG2\tHMGA1\tDSC2\tGPR56\tLAMC1\tS100A10\tPKP3\tNDRG2\tEMP3\tEIF4G1\tPSMD4\tPIM1\tGTPBP2\tFADS3\tOSBP\tEPHB2\tPSME4\tDYSF\tBRD2\tMAST2\tMAP2\tUBE2E3\tSQSTM1\tRGS2\tSYNGR1\tELK3\tUBE2C\tDIAPH1\tABCD1\tCOL7A1\tCAPN6\tMMP3\tLAMB3\tPPP1R15A\tGNAI1\tRPS6KA4\tABCA2\tF2RL2\tLAMC2\tTRAPPC3\tOR2F1\tAKAP1\tFAP\tGIT1\tKRT8\tRAB3D\tATP1B1\tELAVL2\tPHLDA2\tPCDH17\tENO2\tBDKRB1\tCA7\tGKN1\tADRA1A\tPSMD7\tBMP2\tGPR3\tBTK\tSNCB\tLPP\tNR1I3\tHOXB3\tCDKN1A\tTIAL1\tAPOBEC1\tDLG4\tIL1RN\tMPV17\tHCN3\tSHC3\tRPL23A\tKCNA2\tDYRK1A\tFOXN1\tDCTN2\tPADI4\tMARK1\tRAB30\tGFAP\tAKT3\tEEF1A1\tLMNA\tLTB4R2\tSCRN1\tHMGA2\tMMP1\tSPP1\tCXCR4", what = "character", sep = "\t", quiet = TRUE)
BPMS_genes <- scan(text = "BMPER\tDYM\tFBXO42\tFRMPD4\tHERC3\tHS3ST3B1\tIL1RAP\tIL7\tMAGEC1\tMYCT1\tPDE1C\tPRDM1\tRCAN3", what = "character", sep = "\t", quiet = TRUE)

# ssGSEA analysis starts here
perform_ssgsea <- function(tpm_file, genesets) {
    tpm_data <- read.table(tpm_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

    # Handle duplicate gene names
    if ("Gene" %in% colnames(tpm_data)) {
        tpm_data <- tpm_data[!duplicated(tpm_data$Gene), ]
        rownames(tpm_data) <- tpm_data$Gene
        tpm_data <- tpm_data[, -which(colnames(tpm_data) == "Gene")]
    }

    tpm_matrix <- as.matrix(tpm_data)

    # Perform ssGSEA
    param <- ssgseaParam(exprData = tpm_matrix, geneSets = genesets, minSize = 5, use = "na.rm")
    ssgsea_results <- gsva(param, verbose = TRUE)

    return(ssgsea_results)
}

genesets <- list(RPMS_genes = RPMS_genes, BPMS_genes = BPMS_genes)

for (cancer in cancer_types) {
    tpm_file <- paste0("./Transcriptomics/TCGA_DATA/tpm_expression_data/TCGA-", cancer, ".star_tpm.tsv")
    ssgsea_results <- perform_ssgsea(tpm_file, genesets)

    ssgsea_results_file <- paste0("./Transcriptomics/TCGA_DATA/tpm_ssgsea_results/TCGA-", cancer, "_ssgsea_results.tsv")
    write.table(ssgsea_results, file = ssgsea_results_file, sep = "\t", col.names = NA, quote = FALSE)
}

# Since the results are not in the layout as expected. We do additional pre-procesing to get it in the intended format
transpose_and_replace <- function(file) {
    data <- read.table(file, header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
    transposed_data <- t(data)
    rownames(transposed_data) <- gsub("\\.", "-", rownames(transposed_data))
    write.table(transposed_data, file = file, sep = "\t", col.names = NA, quote = FALSE)
}

for (cancer in cancer_types) {
    ssgsea_results_file <- paste0("./Transcriptomics/TCGA_DATA/tpm_ssgsea_results/TCGA-", cancer, "_ssgsea_results.tsv")
    transpose_and_replace(ssgsea_results_file)
}

# We add the Type column which will be useful for doing survival analysis later on
# Currently we are doing only for RPMS and BPMS, but later we add more signatures

add_type_column <- function(file) {
    data <- read.table(file, header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

    rpms_median <- median(data$RPMS_genes)
    bpms_median <- median(data$BPMS_genes)

    data$Type <- ifelse(data$RPMS_genes > rpms_median & data$BPMS_genes > bpms_median, "RPMS+/BPMS+",
                 ifelse(data$RPMS_genes > rpms_median & data$BPMS_genes <= bpms_median, "RPMS+/BPMS-",
                 ifelse(data$RPMS_genes <= rpms_median & data$BPMS_genes > bpms_median, "RPMS-/BPMS+",
                       "RPMS-/BPMS-")))

    write.table(data, file = file, sep = "\t", col.names = NA, quote = FALSE)
}

for (cancer in cancer_types) {
    ssgsea_results_file <- paste0("./Transcriptomics/TCGA_DATA/tpm_ssgsea_results/TCGA-", cancer, "_ssgsea_results.tsv")
    add_type_column(ssgsea_results_file)
}
