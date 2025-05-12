# Author: Anjaney J Pandey
# Code for reproducing biomolecules 2022 by Srinath. At the end, there is code for METABRIC ssgsea for all signatures

library(survival)
library(broom)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(GSVA)

cancer_types <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")

# Gene sets for metabolic signatures

EMT <- c("ABI3BP","ACTA2","ADAM12","ANPEP","APLP1","AREG","BASP1","BDNF","BGN","BMP1","CADM1","CALD1","CALU","CAP2","CAPG","CCN1","CCN2","CD44","CD59","CDH11","CDH2","CDH6","COL11A1","COL12A1","COL16A1","COL1A1","COL1A2","COL3A1","COL4A1","COL4A2","COL5A1","COL5A2","COL5A3","COL6A2","COL6A3","COL7A1","COL8A2","COLGALT1","COMP","COPA","CRLF1","CTHRC1","CXCL1","CXCL12","CXCL6","CXCL8","DAB2","DCN","DKK1","DPYSL3","DST","ECM1","ECM2","EDIL3","EFEMP2","ELN","EMP3","ENO2","FAP","FAS","FBLN1","FBLN2","FBLN5","FBN1","FBN2","FERMT2","FGF2","FLNA","FMOD","FN1","FOXC2","FSTL1","FSTL3","FUCA1","FZD8","GADD45A","GADD45B","GAS1","GEM","GJA1","GLIPR1","GPC1","GPX7","GREM1","HTRA1","ID2","IGFBP2","IGFBP3","IGFBP4","IL15","IL32","IL6","INHBA","ITGA2","ITGA5","ITGAV","ITGB1","ITGB3","ITGB5","JUN","LAMA1","LAMA2","LAMA3","LAMC1","LAMC2","LGALS1","LOX","LOXL1","LOXL2","LRP1","LRRC15","LUM","MAGEE1","MATN2","MATN3","MCM7","MEST","MFAP5","MGP","MMP1","MMP14","MMP2","MMP3","MSX1","MXRA5","MYL9","MYLK","NID2","NNMT","NOTCH2","NT5E","NTM","OXTR","P3H1","PCOLCE","PCOLCE2","PDGFRB","PDLIM4","PFN2","PLAUR","PLOD1","PLOD2","PLOD3","PMEPA1","PMP22","POSTN","PPIB","PRRX1","PRSS2","PTHLH","PTX3","PVR","QSOX1","RGS4","RHOB","SAT1","SCG2","SDC1","SDC4","SERPINE1","SERPINE2","SERPINH1","SFRP1","SFRP4","SGCB","SGCD","SGCG","SLC6A8","SLIT2","SLIT3","SNAI2","SNTB1","SPARC","SPOCK1","SPP1","TAGLN","TFPI2","TGFB1","TGFBI","TGFBR3","TGM2","THBS1","THBS2","THY1","TIMP1","TIMP3","TNC","TNFAIP3","TNFRSF11B","TNFRSF12A","TPM1","TPM2","TPM4","VCAM1","VCAN","VEGFA","VEGFC","VIM","WIPF1","WNT5A")

Glycolysis <- c("ABCB6","ADORA2B","AGL","AGRN","AK3","AK4","AKR1A1","ALDH7A1","ALDH9A1","ALDOA","ALDOB","ALG1","ANG","ANGPTL4","ANKZF1","ARPP19","ARTN","AURKA","B3GALT6","B3GAT1","B3GAT3","B3GNT3","B4GALT1","B4GALT2","B4GALT4","B4GALT7","BIK","BPNT1","CACNA1H","CAPN5","CASP6","CD44","CDK1","CENPA","CHPF","CHPF2","CHST1","CHST12","CHST2","CHST4","CHST6","CITED2","CLDN3","CLDN9","CLN6","COG2","COL5A1","COPB2","CTH","CXCR4","CYB5A","DCN","DDIT4","DEPDC1","DLD","DPYSL4","DSC2","ECD","EFNA3","EGFR","EGLN3","ELF3","ENO1","ENO2","ERO1A","EXT1","EXT2","FAM162A","FBP2","FKBP4","FUT8","G6PD","GAL3ST1","GALE","GALK1","GALK2","GAPDHS","GCLC","GFPT1","GLCE","GLRX","GMPPA","GMPPB","GNE","GNPDA1","GOT1","GOT2","GPC1","GPC3","GPC4","GPR87","GUSB","GYS1","GYS2","HAX1","HDLBP","HK2","HMMR","HOMER1","HS2ST1","HS6ST2","HSPA5","IDH1","IDUA","IER3","IGFBP3","IL13RA1","IRS2","ISG20","KDELR3","KIF20A","KIF2A","LCT","LDHA","LDHC","LHPP","LHX9","MDH1","MDH2","ME1","ME2","MED24","MERTK","MET","MIF","MIOX","MPI","MXI1","NANP","NASP","NDST3","NDUFV3","NOL3","NSDHL","NT5E","P4HA1","P4HA2","PAM","PAXIP1","PC","PDK3","PFKFB1","PFKP","PGAM1","PGAM2","PGK1","PGLS","PGM2","PHKA2","PKM","PKP2","PLOD1","PLOD2","PMM2","POLR3K","PPFIA4","PPIA","PPP2CB","PRPS1","PSMC4","PYGB","PYGL","QSOX1","RARS1","RBCK1","RPE","RRAGD","SAP30","SDC1","SDC2","SDC3","SDHC","SLC16A3","SLC25A10","SLC25A13","SLC35A3","SLC37A4","SOD1","SOX9","SPAG4","SRD5A3","STC1","STC2","STMN1","TALDO1","TFF3","TGFA","TGFBI","TKTL1","TPBG","TPI1","TPST1","TSTA3","TXN","UGP2","VCAN","VEGFA","VLDLR","XYLT2","ZNF292")

OXPHOS <- c("ABCB7","ACAA1","ACAA2","ACADM","ACADSB","ACADVL","ACAT1","ACO2","AFG3L2","AIFM1","ALAS1","ALDH6A1","ATP1B1","ATP5F1A","ATP5F1B","ATP5F1C","ATP5F1D","ATP5F1E","ATP5MC1","ATP5MC2","ATP5MC3","ATP5ME","ATP5MF","ATP5MG","ATP5PB","ATP5PD","ATP5PF","ATP5PO","ATP6AP1","ATP6V0B","ATP6V0C","ATP6V0E1","ATP6V1C1","ATP6V1D","ATP6V1E1","ATP6V1F","ATP6V1G1","ATP6V1H","BAX","BCKDHA","BDH2","CASP7","COX10","COX11","COX15","COX17","COX4I1","COX5A","COX5B","COX6A1","COX6B1","COX6C","COX7A2","COX7A2L","COX7B","COX7C","COX8A","CPT1A","CS","CYB5A","CYB5R3","CYC1","CYCS","DECR1","DLAT","DLD","DLST","ECH1","ECHS1","ECI1","ETFA","ETFB","ETFDH","FDX1","FH","FXN","GLUD1","GOT2","GPI","GPX4","GRPEL1","HADHA","HADHB","HCCS","HSD17B10","HSPA9","HTRA2","IDH1","IDH2","IDH3A","IDH3B","IDH3G","IMMT","ISCA1","ISCU","LDHA","LDHB","LRPPRC","MAOB","MDH1","MDH2","MFN2","MGST3","MPC1","MRPL11","MRPL15","MRPL34","MRPL35","MRPS11","MRPS12","MRPS15","MRPS22","MRPS30","MTRF1","MTRR","MTX2","NDUFA1","NDUFA2","NDUFA3","NDUFA4","NDUFA5","NDUFA6","NDUFA7","NDUFA8","NDUFA9","NDUFAB1","NDUFB1","NDUFB2","NDUFB3","NDUFB4","NDUFB5","NDUFB6","NDUFB7","NDUFB8","NDUFC1","NDUFC2","NDUFS1","NDUFS2","NDUFS3","NDUFS4","NDUFS6","NDUFS7","NDUFS8","NDUFV1","NDUFV2","NNT","NQO2","OAT","OGDH","OPA1","OXA1L","PDHA1","PDHB","PDHX","PDK4","PDP1","PHB2","PHYH","PMPCA","POLR2F","POR","PRDX3","RETSAT","RHOT1","RHOT2","SDHA","SDHB","SDHC","SDHD","SLC25A11","SLC25A12","SLC25A20","SLC25A3","SLC25A4","SLC25A5","SLC25A6","SUCLA2","SUCLG1","SUPV3L1","SURF1","TCIRG1","TIMM10","TIMM13","TIMM17A","TIMM50","TIMM8B","TIMM9","TOMM22","TOMM70","UQCR10","UQCR11","UQCRB","UQCRC1","UQCRC2","UQCRFS1","UQCRH","UQCRQ","VDAC1","VDAC2","VDAC3")

FAO <- c("ACAD8", "HADH", "ACAA1", "CD36", "ECHS1", "ACADVL", "ACAA2", "ACADSB", "SLC25A20", "ACAD9", "IVD", "ACADS", "GCDH")

Hypoxia <- c("ACKR3", "ADM", "ADORA2B", "AK4", "AKAP12", "ALDOA", "ALDOB", "ALDOC", "AMPD3", "ANGPTL4", "ANKZF1", "ANXA2", "ATF3", "ATP7A", "B3GALT6", "B4GALNT2", "BCAN", "BCL2", "BGN", "BHLHE40", "BNIP3L", "BRS3", "BTG1", "CA12", "CASP6", "CAV1", "CAVIN1", "CAVIN3", "CCN1", "CCN2", "CCN5", "CCNG2", "CDKN1A", "CDKN1B", "CDKN1C", "CHST2", "CHST3", "CITED2", "COL5A1", "CP", "CSRP2", "CXCR4", "DCN", "DDIT3", "DDIT4", "DPYSL4", "DTNA", "DUSP1", "EDN2", "EFNA1", "EFNA3", "EGFR", "ENO1", "ENO2", "ENO3", "ERO1A", "ERRFI1", "ETS1", "EXT1", "F3", "FAM162A", "FBP1", "FOS", "FOSL2", "FOXO3", "GAA", "GALK1", "GAPDH", "GAPDHS", "GBE1", "GCK", "GCNT2", "GLRX", "GPC1", "GPC3", "GPC4", "GPI", "GRHPR", "GYS1", "HAS1", "HDLBP", "HEXA", "HK1", "HK2", "HMOX1", "HOXB9", "HS3ST1", "HSPA5", "IDS", "IER3", "IGFBP1", "IGFBP3", "IL6", "ILVBL", "INHA", "IRS2", "ISG20", "JMJD6", "JUN", "KDELR3", "KDM3A", "KIF5A", "KLF6", "KLF7", "KLHL24", "LALBA", "LARGE1", "LDHA", "LDHC", "LOX", "LXN", "MAFF", "MAP3K1", "MIF", "MT1E", "MT2A", "MXI1", "MYH9", "NAGK", "NCAN", "NDRG1", "NDST1", "NDST2", "NEDD4L", "NFIL3", "NOCT", "NR3C1", "P4HA1", "P4HA2", "PAM", "PCK1", "PDGFB", "PDK1", "PDK3", "PFKFB3", "PFKL", "PFKP", "PGAM2", "PGF", "PGK1", "PGM1", "PGM2", "PHKG1", "PIM1", "PKLR", "PKP1", "PLAC8", "PLAUR", "PLIN2", "PNRC1", "PPARGC1A", "PPFIA4", "PPP1R15A", "PPP1R3C", "PRDX5", "PRKCA", "PYGM", "RBPJ", "RORA", "RRAGD", "S100A4", "SAP30", "SCARB1", "SDC2", "SDC3", "SDC4", "SELENBP1", "SERPINE1", "SIAH2", "SLC25A1", "SLC2A1", "SLC2A3", "SLC2A5", "SLC37A4", "SLC6A6", "SRPX", "STBD1", "STC1", "STC2", "SULT2B1", "TES", "TGFB3", "TGFBI", "TGM2", "TIPARP", "TKTL1", "TMEM45A", "TNFAIP3", "TPBG", "TPD52", "TPI1", "TPST2", "UGP2", "VEGFA", "VHL", "VLDLR", "WSB1", "XPNPEP1", "ZFP36", "ZNF292")

Ferroptosis <- c("GPX4", "PANX1", "BID", "SLC3A2", "CD44", "ALOX5", "LPIN1", "ATM", "GCLC", "HMOX1", "SESN2", "OTUB1", "ZFP36", "TF", "CHAC1", "SCP2", "BRD4", "MAPK3", "ENPP2", "TP63", "CA9", "ALOX12", "PGD", "DPP4", "LAMP2", "AKR1C1", "CISD2", "MUC1", "CAV1", "SOCS1", "SLC7A11", "ALOX15B", "CDKN2A", "PRKAA2", "FADS2", "MT1G", "CBS", "TFRC", "TLR4", "PTGS2", "ACSL4", "HELLS", "ZEB1", "TNFAIP3")

ROS <- c("ABCC1", "ATOX1", "CAT", "CDKN2D", "EGLN2", "ERCC2", "FES", "FTL", "G6PD", "GCLC", "GCLM", "GLRX", "GLRX2", "GPX3", "GPX4", "GSR", "HHEX", "HMOX2", "IPCEF1", "JUNB", "LAMTOR5", "LSP1", "MBP", "MGST1", "MPO", "MSRA", "NDUFA6", "NDUFB4", "NDUFS2", "NQO1", "OXSR1", "PDLIM1", "PFKP", "PRDX1", "PRDX2", "PRDX4", "PRDX6", "PRNP", "PTPA", "SBNO2", "SCAF4", "SELENOS", "SOD1", "SOD2", "SRXN1", "STK25", "TXN", "TXNRD1", "TXNRD2")

RPMS <- c("PEBP1", "IGF2BP2", "GOLT1B", "NRAS", "RDX", "NF2", "DMD", "BACH1", "MAP4K4", "COL5A2", "CDV3", "COL3A1", "HMGA2", "COL1A2", "TGFBR1", "PAPPA", "IGF2BP3", "ARID3B", "LRIG2", "IRS2", "HIC2", "IL13", "GNPTAB", "KRT16", "TNFRSF12A", "ENO1", "EN1", "TAGLN2", "BAG2", "HMGA1", "DSC2", "GPR56", "LAMC1", "S100A10", "PKP3", "NDRG2", "EMP3", "EIF4G1", "PSMD4", "PIM1", "GTPBP2", "FADS3", "OSBP", "EPHB2", "PSME4", "DYSF", "BRD2", "MAST2", "MAP2", "UBE2E3", "SQSTM1", "RGS2", "SYNGR1", "ELK3", "UBE2C", "DIAPH1", "ABCD1", "COL7A1", "CAPN6", "MMP3", "LAMB3", "PPP1R15A", "GNAI1", "RPS6KA4", "ABCA2", "F2RL2", "LAMC2", "TRAPPC3", "OR2F1", "AKAP1", "FAP", "GIT1", "KRT8", "RAB3D", "ATP1B1", "ELAVL2", "PHLDA2", "PCDH17", "ENO2", "BDKRB1", "CA7", "GKN1", "ADRA1A", "PSMD7", "BMP2", "GPR3", "BTK", "SNCB", "LPP", "NR1I3", "HOXB3", "CDKN1A", "TIAL1", "APOBEC1", "DLG4", "IL1RN", "MPV17", "HCN3", "SHC3", "RPL23A", "KCNA2", "DYRK1A", "FOXN1", "DCTN2", "PADI4", "MARK1", "RAB30", "GFAP", "AKT3", "EEF1A1", "LMNA", "LTB4R2", "SCRN1", "HMGA2", "MMP1", "SPP1", "CXCR4")

BPMS <- c("BMPER", "DYM", "FBXO42", "FRMPD4", "HERC3", "HS3ST3B1", "IL1RAP", "IL7", "MAGEC1", "MYCT1", "PDE1C", "PRDM1", "RCAN3")

# ssGSEA analysis
genesets <- list(EMT = EMT, Glycolysis = Glycolysis, OXPHOS = OXPHOS, FAO = FAO, Hypoxia = Hypoxia, Ferroptosis = Ferroptosis, ROS = ROS, RPMS = RPMS, BPMS = BPMS)

perform_ssgsea <- function(expression_file, genesets) {
    # Read data with Gene as rownames
    expression_data <- read.table(expression_file, header = TRUE, sep = "\t", 
                                stringsAsFactors = FALSE)
    
    # Handle duplicate genes by keeping the first occurrence
    expression_data <- expression_data[!duplicated(expression_data$Gene), ]
    
    # Set Gene column as rownames and remove it from the data
    rownames(expression_data) <- expression_data$Gene
    expression_data <- expression_data[, -1]
    
    # Convert to matrix
    expression_matrix <- as.matrix(expression_data)
    
    # Perform ssGSEA
    param <- ssgseaParam(exprData = expression_matrix, 
                        geneSets = genesets, 
                        normalize = TRUE, 
                        minSize = 5, 
                        checkNA = "yes", 
                        use = "na.rm")
    ssgsea_results <- gsva(param, verbose = TRUE)
    
    return(ssgsea_results)
}

for (cancer in cancer_types) {
    expression_file <- paste0("./Transcriptomics/TCGA_DATA/TCGA-", cancer, ".star_tpm.tsv")
    ssgsea_results <- perform_ssgsea(expression_file, genesets)
    
    ssgsea_results_file <- paste0("./Transcriptomics/srinath_2022/ssgsea_tpm_results/TCGA-", cancer, "_ssgsea_results.tsv")
    write.table(ssgsea_results, file = ssgsea_results_file, sep = "\t", col.names = NA, quote = FALSE)
}

transpose_and_replace <- function(file) {
    data <- read.table(file, header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
    transposed_data <- t(data)
    rownames(transposed_data) <- gsub("\\.", "-", rownames(transposed_data))
    write.table(transposed_data, file = file, sep = "\t", col.names = NA, quote = FALSE)
}

for (cancer in cancer_types) {
    ssgsea_results_file <- paste0("./Transcriptomics/srinath_2022/ssgsea_tpm_results/TCGA-", cancer, "_ssgsea_results.tsv")
    transpose_and_replace(ssgsea_results_file)
}

add_type_column <- function(file) {
  data <- read.table(file, header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
  
  data$EMT_level <- ifelse(data$EMT > median(data$EMT), "EMT+", "EMT-")
  data$Glycolysis_level <- ifelse(data$Glycolysis > median(data$Glycolysis), "Gly+", "Gly-")
  data$OXPHOS_level <- ifelse(data$OXPHOS > median(data$OXPHOS), "OXP+", "OXP-")
  data$FAO_level <- ifelse(data$FAO > median(data$FAO), "FAO+", "FAO-")
  data$Hypoxia_level <- ifelse(data$Hypoxia > median(data$Hypoxia), "Hyp+", "Hyp-")
  data$Ferroptosis_level <- ifelse(data$Ferroptosis > median(data$Ferroptosis), "Fer+", "Fer-")
  data$ROS_level <- ifelse(data$ROS > median(data$ROS), "ROS+", "ROS-")
  data$RPMS_level <- ifelse(data$RPMS > median(data$RPMS), "RPMS+", "RPMS-")
  data$BPMS_level <- ifelse(data$BPMS > median(data$BPMS), "BPMS+", "BPMS-")
  write.table(data, file = file, sep = "\t", col.names = NA, quote = FALSE)
}

for (cancer in cancer_types) {
    ssgsea_results_file <- paste0("./Transcriptomics/srinath_2022/ssgsea_tpm_results/TCGA-", cancer, "_ssgsea_results.tsv")
    add_type_column(ssgsea_results_file)
}

# Function to read and merge survival and ssGSEA data
merge_survival_ssgsea <- function(cancer) {
    # Read survival data
    survival_file <- paste0("./Transcriptomics/TCGA_DATA/survival_data/TCGA-", cancer, ".survival.tsv")
    survival_data <- read.table(survival_file, header = TRUE, sep = "\t")
    
    # Read ssGSEA results
    ssgsea_file <- paste0("./Transcriptomics/srinath_2022/ssgsea_tpm_results/TCGA-", cancer, "_ssgsea_results.tsv")
    ssgsea_data <- read.table(ssgsea_file, header = TRUE, sep = "\t", row.names = 1)
    
    # Add sample IDs as column
    ssgsea_data$sample <- rownames(ssgsea_data)
    
    # Merge data
    merged_data <- merge(survival_data, ssgsea_data, by = "sample")
    merged_data$cancer_type <- cancer
    
    return(merged_data)
}

# Combine data for all cancer types
tcga_data <- do.call(rbind, lapply(cancer_types, merge_survival_ssgsea))

# Function to clean data
clean_data <- function(data) {
  data <- data %>%
    filter(!is.na(OS.time) & !is.na(OS)) %>%
    filter(OS.time > 0)
    data <- data %>%
      select(-matches("X_PATIENT|Redaction"))
  return(data)
}

# Clean the data
tcga_data <- clean_data(tcga_data)

# Set up factors with explicit reference levels
tcga_data <- tcga_data %>%
  mutate(
    EMT_level = factor(EMT_level, levels = c("EMT-", "EMT+")),      
    Glycolysis_level = factor(Glycolysis_level, levels = c("Gly-", "Gly+")),  
    OXPHOS_level = factor(OXPHOS_level, levels = c("OXP-", "OXP+")),
    FAO_level = factor(FAO_level, levels = c("FAO-", "FAO+")),
    Hypoxia_level = factor(Hypoxia_level, levels = c("Hyp-", "Hyp+")),
    Ferroptosis_level = factor(Ferroptosis_level, levels = c("Fer-", "Fer+")),
    ROS_level = factor(ROS_level, levels = c("ROS-", "ROS+")),
    RPMS_level = factor(RPMS_level, levels = c("RPMS-", "RPMS+")),
    BPMS_level = factor(BPMS_level, levels = c("BPMS-", "BPMS+"))
    )

# Function to calculate Cox model for each cancer type and signature
calculate_cox <- function(data, signature) {
  data %>%
    group_by(cancer_type) %>%
    do({
      
      # Fit the Cox model with signature
      mod <- coxph(Surv(OS.time, OS) ~ (get(signature)), data = .)
      # Tidy the model output
      tidy_mod <- tidy(mod, exponentiate = TRUE, conf.int = TRUE)
      tidy_mod$signature <- signature
      tidy_mod$cancer_type <- unique(.$cancer_type)
      tidy_mod
    }) %>%
    ungroup()
}

# Calculate Cox models for each signature
cox_emt <- calculate_cox(tcga_data, "EMT_level")
cox_glycolysis <- calculate_cox(tcga_data, "Glycolysis_level")
cox_oxphos <- calculate_cox(tcga_data, "OXPHOS_level")
cox_fao <- calculate_cox(tcga_data, "FAO_level")
cox_hypoxia <- calculate_cox(tcga_data, "Hypoxia_level")
cox_ferroptosis <- calculate_cox(tcga_data, "Ferroptosis_level")
cox_ros <- calculate_cox(tcga_data, "ROS_level")
cox_rpms <- calculate_cox(tcga_data, "RPMS_level")
cox_bpms <- calculate_cox(tcga_data, "BPMS_level")

# Combine results
cox_results <- bind_rows(cox_emt, cox_glycolysis, cox_oxphos, cox_fao, cox_hypoxia, cox_ferroptosis, cox_ros, cox_rpms, cox_bpms)

# Add significance stars
cox_results <- cox_results %>%
  mutate(sig = case_when(
    p.value < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Filter the cox_results data frame
filtered_cox_results_1 <- cox_results %>% 
  filter(signature %in% c("OXPHOS_level", "FAO_level", "Glycolysis_level"))

filtered_cox_results_2 <- cox_results %>%
  filter(signature %in% c("EMT_level", "RPMS_level", "BPMS_level"))

filtered_cox_results_3 <- cox_results %>%
  filter(signature %in% c("Hypoxia_level", "Ferroptosis_level", "ROS_level")) %>% 
  filter(!is.na(conf.high) & is.finite(conf.high))

# Create forest plot with vertical panels sharing x-axis
ggplot(filtered_cox_results_3, aes(x = cancer_type, y = estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0) +
  geom_text(aes(label = paste0(" ", sig),
                y = conf.high),
            position = position_nudge(x = 0.1),
            size = 15) +
  facet_wrap(~ signature, ncol = 1, scales = "free_y") +
  geom_hline(yintercept = 1, linetype = "solid", color = "red") +
  scale_y_log10(expand = expansion(mult = c(0.05, 0.2))) +
  coord_cartesian(clip = "off") +
  labs(y = "Hazard Ratio", x = "Cancer Type", size = 15) +
  theme_bw() +
  theme(plot.margin = margin(10, 40, 10, 10),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 15, face = "bold"),
        axis.text.y = element_text(size = 18, face = "bold"),
        panel.spacing = unit(1, "cm"),
        strip.text = element_text(size = 23, face = "bold"))


###############################################################################################
# Create forest plot with vertical panels sharing x-axis using log2 scale
ggplot(cox_results, aes(x = cancer_type, y = log2(estimate))) +
  geom_point() +
  geom_errorbar(aes(ymin = log2(conf.low), ymax = log2(conf.high)), width = 0) +
  geom_text(aes(label = paste0(" ", sig),
                y = log2(conf.high)),
            position = position_nudge(x = 0.1),
            hjust = -0.2, size = 10) +
  facet_wrap(~ signature, ncol = 1, scales = "free_y") +
  geom_hline(yintercept = 0, linetype = "solid", color = "red") +
  coord_cartesian(clip = "off") +
  labs(y = "log2(Hazard Ratio)", x = "Cancer Type",
       title = "Log2 Hazard Ratios") +
  theme_bw() +
  theme(plot.margin = margin(10, 40, 10, 10),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        panel.spacing = unit(1, "cm"))

# Save the plot with adjusted dimensions for vertical layout
ggsave(filename = "./Transcriptomics/srinath_2022/survival_plots/TCGA_Metabolic_Signatures.png",
       plot = last_plot(),
       width = 12,
       height = 18,
       dpi = 300,      
       type = "cairo",
       bg = "white")
################################################################################################

# Define the orders
cancer_order <- c("PCPG", "DLBC", "PRAD", "KICH", "THCA", "READ", "UVM", "UCEC", 
                 "ESCA", "LUSC", "ACC", "OV", "SARC", "BRCA", "SKCM", "KIRC", 
                 "COAD", "STAD", "LGG", "MESO", "PAAD", "BLCA", "CESC", "LUAD", 
                 "GBM", "LIHC", "CHOL", "KIRP", "HNSC", "UCS", "TGCT", "THYM")

type_order <- c("E+/G+/O-", "E-/G+/O+", "E-/G-/O+", "E+/G-/O-", 
                "E+/G-/O+", "E+/G+/O+", "E-/G+/O-")

# Calculate hazard ratios and p-values
hr_results <- filtered_tcga_data %>%
  group_by(cancer_type) %>%
  do({
    .$Type <- factor(.$Type, levels = c("E-/G-/O-", type_order))
    
    # Fit Cox model
    mod <- tryCatch({
      coxph(Surv(OS.time, OS) ~ Type, data = .)
    }, error = function(e) NULL)
    
    if (!is.null(mod)) {
      # Extract coefficients and p-values
      summary_mod <- summary(mod)
      coef_matrix <- summary_mod$coefficients
      hr <- exp(coef_matrix[, "coef"])
      pvals <- coef_matrix[, "Pr(>|z|)"]
      
      data.frame(
        type = names(hr),
        hr = hr,
        p.value = pvals,
        cancer_type = unique(.$cancer_type)
      )
    } else {
      tibble()(
        type = character(),
        hr = numeric(),
        p.value = numeric(),
        cancer_type = unique(.$cancer_type)
      )
    }
  }) %>%
  ungroup()

# Reshape data for heatmap
heatmap_data <- hr_results %>%
  mutate(
    display_value = ifelse(p.value < 0.05, "*", "")
  ) %>%
  select(cancer_type, type, hr, p.value, display_value) %>%
  pivot_wider(names_from = type, 
              values_from = c(hr, p.value, display_value))

# Create matrices for values and significance indicators
value_matrix <- as.matrix((select(heatmap_data, starts_with("hr_"))))
sig_matrix <- as.matrix(select(heatmap_data, starts_with("display_value_")))

rownames(value_matrix) <- heatmap_data$cancer_type
rownames(sig_matrix) <- heatmap_data$cancer_type

# Create heatmap
final_heatmap <- pheatmap(-log2(value_matrix),
     color = colorRampPalette(rev(brewer.pal(8, "RdBu")))(25),
     cluster_rows = TRUE,
     cluster_cols = TRUE,
     display_numbers = sig_matrix,
     fontsize = 20,
     fontsize_number = 40,
     angle_col = 45,
     main = "-log2 (Hazard Ratio)\n(Reference: E-/G-/O-)",
     labels_row = cancer_order,
     labels_col = type_order,
     legend = TRUE,
     legend_breaks = seq(-5, 5, by = 1),
     legend_labels = sprintf("%.0f" , seq(-5, 5, by = 1)),
     fontsize_row = 22,
     fontsize_col = 22)

# Save the heatmap
ggsave("./Transcriptomics/srinath_2022/survival_plots/Final_Heatmap.png",
       plot = final_heatmap,
       width = 12, height = 18, dpi = 300, bg = "white")


# Calculate hazard ratios and p-values from tcga_data
hr_results <- tcga_data %>%
  group_by(cancer_type) %>%
  do({
    # Use Type column with reference level
    .$Type <- factor(.$Type, levels = c("E-/G-/O-",
                                      "E+/G+/O-", "E-/G+/O+", "E-/G-/O+", 
                                      "E+/G-/O-", "E+/G-/O+", "E+/G+/O+", 
                                      "E-/G+/O-"))
    
    # Fit Cox model
    mod <- tryCatch({
      coxph(Surv(OS.time, OS) ~ Type, data = .)
    }, error = function(e) NULL)
    
    if (!is.null(mod)) {
      # Extract coefficients and p-values
      summary_mod <- summary(mod)
      coef_matrix <- summary_mod$coefficients
      hr <- exp(coef_matrix[, "coef"])
      pvals <- coef_matrix[, "Pr(>|z|)"]
      
      data.frame(
        type = names(hr),
        hr = hr,
        p.value = pvals,
        cancer_type = unique(.$cancer_type)
      )
    } else {
      data.frame(
        type = character(),
        hr = numeric(),
        p.value = numeric(),
        cancer_type = unique(.$cancer_type)
      )
    }
  }) %>%
  ungroup()

# Create matrix for heatmap
hr_matrix <- hr_results %>%
  mutate(
    display_value = ifelse(p.value < 0.05, 
                          sprintf("%.2f*", hr),
                          sprintf("%.2f", hr))
  ) %>%
  select(cancer_type, type, hr, display_value) %>%
  pivot_wider(names_from = type,
              values_from = c(hr, display_value))

# Prepare matrices for values and display
value_matrix <- as.matrix(select(hr_matrix, starts_with("hr_")))
display_matrix <- as.matrix(select(hr_matrix, starts_with("display_value_")))
rownames(value_matrix) <- hr_matrix$cancer_type
rownames(display_matrix) <- hr_matrix$cancer_type

# Create heatmap with clustering
pheatmap(value_matrix,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(0, 3, length.out = 101),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize = 17,
         fontsize_number = 8,
         angle_col = 45,
         main = "Hazard Ratios\n(Reference: E-/G-/O-)")


##############################################################################

library(RColorBrewer)
library(pheatmap)
library(tidyr)
library(tibble)

# Calculate hazard ratios and p-values
hr_results <- filtered_tcga_data %>%
  group_by(cancer_type) %>%
  do({
    .$Type <- factor(.$Type, levels = c("E-/G-/O-", type_order))
    
    # Fit Cox model
    mod <- tryCatch({
      coxph(Surv(OS.time, OS) ~ Type, data = .)
    }, error = function(e) NULL)
    
    if (!is.null(mod)) {
      # Extract coefficients and p-values
      summary_mod <- summary(mod)
      coef_matrix <- summary_mod$coefficients
      hr <- exp(coef_matrix[, "coef"])
      pvals <- coef_matrix[, "Pr(>|z|)"]
      
      data.frame(
        type = names(hr),
        hr = hr,
        p.value = pvals,
        cancer_type = unique(.$cancer_type)
      )
    } else {
      data.frame(
        type = character(),
        hr = numeric(),
        p.value = numeric(),
        cancer_type = unique(.$cancer_type)
      )
    }
  }) %>%
  ungroup()

# Reshape data for heatmap
heatmap_data <- hr_results %>%
  mutate(
    display_value = ifelse(p.value < 0.05, "*", "")
  ) %>%
  select(cancer_type, type, hr, p.value, display_value) %>%
  pivot_wider(names_from = type, 
              values_from = c(hr, p.value, display_value))

# Create matrices for values and significance indicators
value_matrix <- as.matrix(select(heatmap_data, starts_with("hr_")))
sig_matrix <- as.matrix(select(heatmap_data, starts_with("display_value_")))

rownames(value_matrix) <- heatmap_data$cancer_type
rownames(sig_matrix) <- heatmap_data$cancer_type

# Create heatmap
heatmap_plot <- pheatmap(value_matrix,
         color = colorRampPalette(rev(brewer.pal(8, "RdBu")))(25),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = sig_matrix,
         fontsize = 20,
         fontsize_number = 40,
         angle_col = 45,
         main = "Reference: E-/G-/O-",
         labels_row = cancer_order,
         labels_col = type_order)

# Save the heatmap
png("./Transcriptomics/srinath_2022/survival_plots/Heatmap.png",
    width = 12*300,  # multiply by dpi
    height = 16*300,
    res = 300)
print(heatmap_plot)
dev.off()



########################### METABRIC ssGSEA for all signatures ############################
metabric_data <- paste0("./Transcriptomics/brca_metabric/data_mrna_illumina_microarray_zscores_ref_diploid_samples_final_unique.txt")
metabric_data_df <- read.table(metabric_data, header = TRUE, sep = "\t", row.names = 1)

perform_ssgsea <- function(expression_file, genesets) {
    expression_data <- read.table(expression_file, header = TRUE, sep = "\t", 
                                  stringsAsFactors = FALSE)
    rownames(expression_data) <- expression_data[, 1]
    expression_data <- expression_data[, -1]
    expression_data <- expression_data[!duplicated(rownames(expression_data)), ]
    
    # Convert to matrix
    expression_matrix <- as.matrix(expression_data)
    
    # Perform ssGSEA
    param <- ssgseaParam(exprData = expression_matrix, 
                         geneSets = genesets, 
                         normalize = TRUE, 
                         minSize = 5, 
                         checkNA = "yes", 
                         use = "na.rm")
    ssgsea_results <- gsva(param, verbose = TRUE)
    
    return(ssgsea_results)
}

# run ssGSEA using all signatures
metabric_ssgsea <- perform_ssgsea(metabric_data, genesets)

# Save the ssGSEA results
ssgsea_results <- as.data.frame(metabric_ssgsea)
ssgsea_results <- t(ssgsea_results)

ssgsea_results_file <- paste0("./Transcriptomics/brca_metabric/METABRIC_all_signatures_ssgsea_results.tsv")
write.table(ssgsea_results, file = ssgsea_results_file, sep = "\t", col.names = NA, quote = FALSE)


ssgsea_results <- transpose_and_replace(ssgsea_results_file)
ssgsea_results <- add_type_column(ssgsea_results_file)
