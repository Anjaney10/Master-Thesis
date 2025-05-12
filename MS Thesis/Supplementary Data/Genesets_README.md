# FAO_gene_set.gmt
- Source:   
    - Paper "Drug-Tolerant Idling Melanoma Cells Exhibit Theory-Predicted Metabolic Low-Low Phenotype", Frontiers in Oncology   
    - Dongya Jia 1,†, B Bishal Paudel 2,3,4,†, Corey E Hayford 2,3,5, Keisha N Hardeman 2,3, Herbert Levine 1,6,7,*,†, José N Onuchic 1,8,9,10,*,†, Vito Quaranta 2,3,*,†
    - DOI: https://doi.org/10.3389/fonc.2020.01426
- Number of gene-sets: 1
- Description: 
- Used for:

# Hypoxia_AR.gmt
- Source: Soundharya 
- Number of gene-sets: 18 
    - (HIF1_metagene	HIF1	GLYCOLYSIS	PARTIAL_EMT	EPI	MES	HALLMARK_EMT HALLMARK_HYPOXIA	GOBP_ANOIKIS	GOBP_REGULATION_OF_ANOIKIS	GOBP_NEGATIVE_REGULATION_OF_ANOIKIS	(COAD) DEGs|Genes realted Anoikis	126 ANRGs from the Genecards and Harmonizome portals	BRCA 16 ANRGs associated with prognosis	Negative Regulation Of Anoikis - HARMONIZOME	clusterA_t	clusterB_t clusterC_t)
- Description: 
- Used for:

# OXPHOS_geneset.gmt
- Source: 
    - Paper "Quantifying the Patterns of Metabolic Plasticity and Heterogeneity along the Epithelial–Hybrid–Mesenchymal Spectrum in Cancer" 
    - Srinath Muralidharan 1, Sarthak Sahoo 2, Aryamaan Saha 1,†, Sanjay Chandran 1,†, Sauma Suvra Majumdar 3,†, Susmita Mandal 2, Herbert Levine 4,*, Mohit Kumar Jolly 2,*
    - DOI: https://doi.org/10.3390/biom12020297 
- Number of genes-sets: 1
- Description:
- Used for: 

# Glycolysis_geneset.gmt
- Source: 
    - Paper "Quantifying the Patterns of Metabolic Plasticity and Heterogeneity along the Epithelial–Hybrid–Mesenchymal Spectrum in Cancer" 
    - Srinath Muralidharan 1, Sarthak Sahoo 2, Aryamaan Saha 1,†, Sanjay Chandran 1,†, Sauma Suvra Majumdar 3,†, Susmita Mandal 2, Herbert Levine 4,*, Mohit Kumar Jolly 2,*
    - DOI: https://doi.org/10.3390/biom12020297  
- Number of genes-sets: 1
- Description: 
- Used for: 

# Basal_luminal_gene_sets.gmt
- Source: 
    - Citations(26, 27, 28, 29) of paper "Increased prevalence of hybrid epithelial/ mesenchymal state and enhanced phenotypic heterogeneity in basal breast cancer"
    - Sarthak Sahoo, Soundharya Ramu, Madhumathy G. Nair, ..., Jyothi S. Prabhu, Jason A. Somarelli, Mohit Kumar Jolly
    - DOI: https://doi.org/10.1016/j.isci.2024.110116
- Number of genes-sets: 13
    - (Epi_Tumor	Epi_Cellline	Mes_Tumor	Mes_Cellline	luminal	Basal	pEMT	EMT_partial_celline	EMT_up_celline	EMT_down_celline	EMT_down	EMT_partial	EMT_up)
- Description:
- Used for: /home/csb/Ritesh/iscience_analysis/codes/meta_analysis_basal_luminal_gmt.ipynb

# Mammary_sp_sig2.gmt
- Source: Soundharya 
- Number of genes-sets: 17
    - (Epi_Tumor	Epi_Cellline	Mes_Tumor	Mes_Cellline	basal	luminal HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION	pEMT	EMT_partial_celline	EMT_up_celline	EMT_down_celline	EMT_down	EMT_partial	EMT_up	clusterA_mb	clusterB_mb	clusterC_mb)
- Description: 
- Used for: 

# Metabolic_gene_sets.gmt
- Source: subset of gene sets in Hypoxia_AR.gmt, FAO_gene_set.gmt, Glycolysis_geneset.gmt, OXPHOS_geneset.gmt
- NUmber of gene sets: 6
    - (FAO, HALLMARK_GLYCOLYSIS, HALLMARK_OXIDATIVE_PHOSPHORYLATION, HIF1_metagene, HIF1, HALLMARK_HYPOXIA)
- Description: 
- Used for: /home/csb/Ritesh/iscience_analysis/codes/meta_analysis_metabolic_gmt.ipynb

# PDH_related_gene_stes.gmt
- Source: MSigDB
- Number of gene sets: 4
- Description: 
- Used for: 

# Heme_related_gene_stes.gmt
- Source: MSigDB
- Number of gene sets: 5
- Description: 
- Used for: 

# Glutathione_related_gene_stes.gmt
- Source: MSigDB
- Number of gene sets: 5
- Description: 
- Used for: 

# GenSig_Ferroptosis_Immune_45.gmt
- Source: 
    - Paper: "Systematic profiling of ferroptosis gene signatures predicts prognostic factors in esophageal squamous cell carcinoma"
    - (Tong Lu 1, Ran Xu 1, Qi Li 2, Jia-ying Zhao 1, Bo Peng 1, Han Zhang 1 , Ji-da Guo 1 , Sheng-qiang Zhang 1 , Hua-wei Li 1 , Jun Wang 1 , Lin-you Zhang 1)
    - DOI: https://doi.org/10.1016/j.omto.2021.02.011 
- Number of gene sets: 1
- Description: This geneset was used in RKIP-BACH1 paper as ferroptosis gene signature
- Used for: 

# RPMS_BPMS.gmt
- Source: Sai, RKIP-BACH1 paper's resources
- Number of gene sets: 2
- Description: 
- Used for: 

# HALLMARK_ESTROGEN_RESPONSE_EARLY_LATE.gmt
- Source: MSigDB
- Number of gene sets: 2
- Description:
- Used for: 

# gene_sets_corellated_with_RPMS_BPMS.gmt
- Source: compilation of saperate gene sets
- Number of gene sets: 18
- Description: analysing corellation of RPMS/ BPMS with other related genesets
- Used for: "/home/csb/Ritesh/iscience_analysis/codes/meta_analysis_RPMS_BPMS.ipynb"