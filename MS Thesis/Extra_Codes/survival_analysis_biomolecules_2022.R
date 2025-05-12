library(survival)
library(broom)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(GSVA)

# TCGA Cancer Types

cancer_types <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")

# Gene sets

EMT <- c("ABI3BP","ACTA2","ADAM12","ANPEP","APLP1","AREG","BASP1","BDNF","BGN","BMP1","CADM1","CALD1","CALU","CAP2","CAPG","CCN1","CCN2","CD44","CD59","CDH11","CDH2","CDH6","COL11A1","COL12A1","COL16A1","COL1A1","COL1A2","COL3A1","COL4A1","COL4A2","COL5A1","COL5A2","COL5A3","COL6A2","COL6A3","COL7A1","COL8A2","COLGALT1","COMP","COPA","CRLF1","CTHRC1","CXCL1","CXCL12","CXCL6","CXCL8","DAB2","DCN","DKK1","DPYSL3","DST","ECM1","ECM2","EDIL3","EFEMP2","ELN","EMP3","ENO2","FAP","FAS","FBLN1","FBLN2","FBLN5","FBN1","FBN2","FERMT2","FGF2","FLNA","FMOD","FN1","FOXC2","FSTL1","FSTL3","FUCA1","FZD8","GADD45A","GADD45B","GAS1","GEM","GJA1","GLIPR1","GPC1","GPX7","GREM1","HTRA1","ID2","IGFBP2","IGFBP3","IGFBP4","IL15","IL32","IL6","INHBA","ITGA2","ITGA5","ITGAV","ITGB1","ITGB3","ITGB5","JUN","LAMA1","LAMA2","LAMA3","LAMC1","LAMC2","LGALS1","LOX","LOXL1","LOXL2","LRP1","LRRC15","LUM","MAGEE1","MATN2","MATN3","MCM7","MEST","MFAP5","MGP","MMP1","MMP14","MMP2","MMP3","MSX1","MXRA5","MYL9","MYLK","NID2","NNMT","NOTCH2","NT5E","NTM","OXTR","P3H1","PCOLCE","PCOLCE2","PDGFRB","PDLIM4","PFN2","PLAUR","PLOD1","PLOD2","PLOD3","PMEPA1","PMP22","POSTN","PPIB","PRRX1","PRSS2","PTHLH","PTX3","PVR","QSOX1","RGS4","RHOB","SAT1","SCG2","SDC1","SDC4","SERPINE1","SERPINE2","SERPINH1","SFRP1","SFRP4","SGCB","SGCD","SGCG","SLC6A8","SLIT2","SLIT3","SNAI2","SNTB1","SPARC","SPOCK1","SPP1","TAGLN","TFPI2","TGFB1","TGFBI","TGFBR3","TGM2","THBS1","THBS2","THY1","TIMP1","TIMP3","TNC","TNFAIP3","TNFRSF11B","TNFRSF12A","TPM1","TPM2","TPM4","VCAM1","VCAN","VEGFA","VEGFC","VIM","WIPF1","WNT5A")

Glycolysis <- c("ABCB6","ADORA2B","AGL","AGRN","AK3","AK4","AKR1A1","ALDH7A1","ALDH9A1","ALDOA","ALDOB","ALG1","ANG","ANGPTL4","ANKZF1","ARPP19","ARTN","AURKA","B3GALT6","B3GAT1","B3GAT3","B3GNT3","B4GALT1","B4GALT2","B4GALT4","B4GALT7","BIK","BPNT1","CACNA1H","CAPN5","CASP6","CD44","CDK1","CENPA","CHPF","CHPF2","CHST1","CHST12","CHST2","CHST4","CHST6","CITED2","CLDN3","CLDN9","CLN6","COG2","COL5A1","COPB2","CTH","CXCR4","CYB5A","DCN","DDIT4","DEPDC1","DLD","DPYSL4","DSC2","ECD","EFNA3","EGFR","EGLN3","ELF3","ENO1","ENO2","ERO1A","EXT1","EXT2","FAM162A","FBP2","FKBP4","FUT8","G6PD","GAL3ST1","GALE","GALK1","GALK2","GAPDHS","GCLC","GFPT1","GLCE","GLRX","GMPPA","GMPPB","GNE","GNPDA1","GOT1","GOT2","GPC1","GPC3","GPC4","GPR87","GUSB","GYS1","GYS2","HAX1","HDLBP","HK2","HMMR","HOMER1","HS2ST1","HS6ST2","HSPA5","IDH1","IDUA","IER3","IGFBP3","IL13RA1","IRS2","ISG20","KDELR3","KIF20A","KIF2A","LCT","LDHA","LDHC","LHPP","LHX9","MDH1","MDH2","ME1","ME2","MED24","MERTK","MET","MIF","MIOX","MPI","MXI1","NANP","NASP","NDST3","NDUFV3","NOL3","NSDHL","NT5E","P4HA1","P4HA2","PAM","PAXIP1","PC","PDK3","PFKFB1","PFKP","PGAM1","PGAM2","PGK1","PGLS","PGM2","PHKA2","PKM","PKP2","PLOD1","PLOD2","PMM2","POLR3K","PPFIA4","PPIA","PPP2CB","PRPS1","PSMC4","PYGB","PYGL","QSOX1","RARS1","RBCK1","RPE","RRAGD","SAP30","SDC1","SDC2","SDC3","SDHC","SLC16A3","SLC25A10","SLC25A13","SLC35A3","SLC37A4","SOD1","SOX9","SPAG4","SRD5A3","STC1","STC2","STMN1","TALDO1","TFF3","TGFA","TGFBI","TKTL1","TPBG","TPI1","TPST1","TSTA3","TXN","UGP2","VCAN","VEGFA","VLDLR","XYLT2","ZNF292")

OXPHOS <- c("ABCB7","ACAA1","ACAA2","ACADM","ACADSB","ACADVL","ACAT1","ACO2","AFG3L2","AIFM1","ALAS1","ALDH6A1","ATP1B1","ATP5F1A","ATP5F1B","ATP5F1C","ATP5F1D","ATP5F1E","ATP5MC1","ATP5MC2","ATP5MC3","ATP5ME","ATP5MF","ATP5MG","ATP5PB","ATP5PD","ATP5PF","ATP5PO","ATP6AP1","ATP6V0B","ATP6V0C","ATP6V0E1","ATP6V1C1","ATP6V1D","ATP6V1E1","ATP6V1F","ATP6V1G1","ATP6V1H","BAX","BCKDHA","BDH2","CASP7","COX10","COX11","COX15","COX17","COX4I1","COX5A","COX5B","COX6A1","COX6B1","COX6C","COX7A2","COX7A2L","COX7B","COX7C","COX8A","CPT1A","CS","CYB5A","CYB5R3","CYC1","CYCS","DECR1","DLAT","DLD","DLST","ECH1","ECHS1","ECI1","ETFA","ETFB","ETFDH","FDX1","FH","FXN","GLUD1","GOT2","GPI","GPX4","GRPEL1","HADHA","HADHB","HCCS","HSD17B10","HSPA9","HTRA2","IDH1","IDH2","IDH3A","IDH3B","IDH3G","IMMT","ISCA1","ISCU","LDHA","LDHB","LRPPRC","MAOB","MDH1","MDH2","MFN2","MGST3","MPC1","MRPL11","MRPL15","MRPL34","MRPL35","MRPS11","MRPS12","MRPS15","MRPS22","MRPS30","MTRF1","MTRR","MTX2","NDUFA1","NDUFA2","NDUFA3","NDUFA4","NDUFA5","NDUFA6","NDUFA7","NDUFA8","NDUFA9","NDUFAB1","NDUFB1","NDUFB2","NDUFB3","NDUFB4","NDUFB5","NDUFB6","NDUFB7","NDUFB8","NDUFC1","NDUFC2","NDUFS1","NDUFS2","NDUFS3","NDUFS4","NDUFS6","NDUFS7","NDUFS8","NDUFV1","NDUFV2","NNT","NQO2","OAT","OGDH","OPA1","OXA1L","PDHA1","PDHB","PDHX","PDK4","PDP1","PHB2","PHYH","PMPCA","POLR2F","POR","PRDX3","RETSAT","RHOT1","RHOT2","SDHA","SDHB","SDHC","SDHD","SLC25A11","SLC25A12","SLC25A20","SLC25A3","SLC25A4","SLC25A5","SLC25A6","SUCLA2","SUCLG1","SUPV3L1","SURF1","TCIRG1","TIMM10","TIMM13","TIMM17A","TIMM50","TIMM8B","TIMM9","TOMM22","TOMM70","UQCR10","UQCR11","UQCRB","UQCRC1","UQCRC2","UQCRFS1","UQCRH","UQCRQ","VDAC1","VDAC2","VDAC3")


# ssGSEA analysis
genesets <- list(
    EMT = EMT,
    Glycolysis = Glycolysis,
    OXPHOS = OXPHOS
)

perform_ssgsea <- function(expression_file, genesets) {
    expression_data <- read.table(expression_file, header = TRUE, sep = "\t", 
                                stringsAsFactors = FALSE)
    
    # Handle duplicate genes
    expression_data <- expression_data[!duplicated(expression_data$Gene), ]
    
    # Set Gene column as rownames
    rownames(expression_data) <- expression_data$Gene
    expression_data <- expression_data[, -1]

    expression_matrix <- as.matrix(expression_data)
    
    # Run ssGSEA
    param <- ssgseaParam(exprData = expression_matrix, 
                        geneSets = genesets, 
                       # normalize = TRUE, 
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

############################################################################################

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

############################################################################################

add_type_column <- function(file) {
    data <- read.table(file, header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
    
    EMT_median <- median(data$EMT)
    Glycolysis_median <- median(data$Glycolysis)
    OXPHOS_median <- median(data$OXPHOS)
    
    data$Type <- ifelse(data$EMT > EMT_median & 
                       data$Glycolysis > Glycolysis_median & 
                       data$OXPHOS > OXPHOS_median, "E+/G+/O+",
                 ifelse(data$EMT > EMT_median & 
                       data$Glycolysis > Glycolysis_median & 
                       data$OXPHOS <= OXPHOS_median, "E+/G+/O-",
                 ifelse(data$EMT > EMT_median & 
                       data$Glycolysis <= Glycolysis_median & 
                       data$OXPHOS > OXPHOS_median, "E+/G-/O+",
                 ifelse(data$EMT <= EMT_median & 
                       data$Glycolysis > Glycolysis_median & 
                       data$OXPHOS > OXPHOS_median, "E-/G+/O+",
                 ifelse(data$EMT > EMT_median & 
                       data$Glycolysis <= Glycolysis_median & 
                       data$OXPHOS <= OXPHOS_median, "E+/G-/O-",
                 ifelse(data$EMT <= EMT_median & 
                       data$Glycolysis > Glycolysis_median & 
                       data$OXPHOS <= OXPHOS_median, "E-/G+/O-",
                 ifelse(data$EMT <= EMT_median & 
                       data$Glycolysis <= Glycolysis_median & 
                       data$OXPHOS > OXPHOS_median, "E-/G-/O+",
                       "E-/G-/O-")))))))
    
    write.table(data, file = file, sep = "\t", col.names = NA, quote = FALSE)
}

for (cancer in cancer_types) {
    ssgsea_results_file <- paste0("./Transcriptomics/srinath_2022/ssgsea_tpm_results/TCGA-", cancer, "_ssgsea_results.tsv")
    add_type_column(ssgsea_results_file)
}

############################################################################################

# Function to merge survival and ssGSEA data
merge_survival_ssgsea <- function(cancer) {
    survival_file <- paste0("./Transcriptomics/TCGA_DATA/survival_data/TCGA-", cancer, ".survival.tsv")
    survival_data <- read.table(survival_file, header = TRUE, sep = "\t")
    
    ssgsea_file <- paste0("./Transcriptomics/srinath_2022/ssgsea_tpm_results/TCGA-", cancer, "_ssgsea_results.tsv")
    ssgsea_data <- read.table(ssgsea_file, header = TRUE, sep = "\t", row.names = 1)
    
    # Add sample IDs
    ssgsea_data$sample <- rownames(ssgsea_data)
    
    # Merge data
    merged_data <- merge(survival_data, ssgsea_data, by = "sample")
    merged_data$cancer_type <- cancer
    
    return(merged_data)
}

# Combine data for all cancer types
tcga_data <- do.call(rbind, lapply(cancer_types, merge_survival_ssgsea))

clean_data <- function(data) {
  data <- data %>%
    filter(!is.na(OS.time) & !is.na(OS) & !is.na(EMT) & !is.na(Glycolysis) & !is.na(OXPHOS)) %>%
    filter(OS.time > 0)
    data <- data %>%
      select(-matches("X_PATIENT|Redaction"))
  return(data)
}

# Clean the data
tcga_data <- clean_data(tcga_data)

########################################################################################

# Function to calculate Cox model for each cancer type and signature
calculate_cox <- function(data, signature) {
  data %>%
    group_by(cancer_type) %>%
    do({
      mod <- coxph(Surv(OS.time, OS) ~ get(signature), data = .)
      tidy_mod <- tidy(mod, exponentiate = TRUE, conf.int = TRUE)
      tidy_mod$signature <- signature
      tidy_mod$cancer_type <- unique(.$cancer_type)
      tidy_mod
    }) %>%
    ungroup()
}

# Calculate Cox models for each signature
cox_emt <- calculate_cox(tcga_data, "EMT")
cox_glycolysis <- calculate_cox(tcga_data, "Glycolysis")
cox_oxphos <- calculate_cox(tcga_data, "OXPHOS")

# Combine results
cox_results <- bind_rows(cox_emt, cox_glycolysis, cox_oxphos)

# Add significance stars
cox_results <- cox_results %>%
  mutate(sig = case_when(
    p.value < 0.05 ~ "*",
    TRUE ~ ""
  ))

selected_cancer_types <- c("BLCA", "CESC", "CHOL", "COAD", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "MESO", "PAAD", "UVM")

# Subset the cox_results
filtered_cox_results <- cox_results %>% filter(cancer_type %in% selected_cancer_types)

#################################################################################################

ggplot(filtered_cox_results, aes(x = cancer_type, y = estimate)) +
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
  labs(y = "Hazard Ratio", x = "Cancer Type") +
  theme_bw() +
  theme(plot.margin = margin(10, 40, 10, 10),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 18, face = "bold"),
        axis.text.y = element_text(size = 18, face = "bold"),
        panel.spacing = unit(1, "cm")) 


#################################################################################################

# Heatmap

cancer_order <- c("PCPG", "DLBC", "PRAD", "KICH", "THCA", "READ", "UVM", "UCEC", 
                 "ESCA", "LUSC", "ACC", "OV", "SARC", "BRCA", "SKCM", "KIRC", 
                 "COAD", "STAD", "LGG", "MESO", "PAAD", "BLCA", "CESC", "LUAD", 
                 "GBM", "LIHC", "CHOL", "KIRP", "HNSC", "UCS", "TGCT", "THYM")

type_order <- c("E+/G+/O-", "E-/G+/O+", "E-/G-/O+", "E+/G-/O-", 
                "E+/G-/O+", "E+/G+/O+", "E-/G+/O-")


# Calculate HR and p-values
hr_results <- tcga_data %>%
  group_by(cancer_type) %>%
  do({
    .$Type <- factor(.$Type, levels = c("E-/G-/O-", type_order))
    
    # Fit Cox model
    mod <- tryCatch({
      coxph(Surv(OS.time, OS) ~ Type, data = .)
    }, error = function(e) NULL)
    
    if (!is.null(mod)) {
      # Model Summary
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

heatmap_data <- hr_results %>%
  mutate(
    display_value = ifelse(p.value < 0.05, "*", "")
  ) %>%
  select(cancer_type, type, hr, p.value, display_value) %>%
  pivot_wider(names_from = type, 
              values_from = c(hr, p.value, display_value))

value_matrix <- as.matrix(select(heatmap_data, starts_with("hr_")))
sig_matrix <- as.matrix(select(heatmap_data, starts_with("display_value_")))

rownames(value_matrix) <- heatmap_data$cancer_type
rownames(sig_matrix) <- heatmap_data$cancer_type

# Heatmap
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

png("./Transcriptomics/srinath_2022/survival_plots/Heatmap.png",
    width = 12*300,
    height = 16*300,
    res = 300)
print(heatmap_plot)
dev.off()
