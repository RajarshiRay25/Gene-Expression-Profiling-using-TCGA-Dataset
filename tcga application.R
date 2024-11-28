## SCRIPT FOR PERFORMING TRANSCRIPTOMIC PROFILING FOR CANCER SAMPLES ##

## DESIGN TYPE - TYPE 1 (LOCALISED) VS TYPE 4 (MALIGNANT)

## PLEASE ENSURE INSTALLATION OF RTOOLS IN THE R FOLDER BEFORE RUNNING THIS SCRIPT

## FROM THE RAW COHORT DATA OPEN IN EXCEL AND USE CONCAT FUNCTION TO ADD "", TO IDS

# get packages


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

remotes::install_github("r-lib/httr2")

BiocManager::install("TCGAbiolinks",force = TRUE)
BiocManager::install("EDASeq",force = TRUE)
BiocManager::install("edgeR",force = TRUE)
BiocManager::install("clusterProfiler",force = TRUE)
BiocManager::install("pathview",force = TRUE)

# Define library

library("TCGAbiolinks")

# Obtain Samples from TCGA with Type 1 cancer condition

bar1 <- c("TCGA-A2-A0YF",
          "TCGA-BH-A0H6",
          "TCGA-A8-A0AD",
          "TCGA-Z7-A8R6",
          "TCGA-GM-A2DD",
          "TCGA-BH-A0BP",
          "TCGA-BH-A209",
          "TCGA-GM-A2DL",
          "TCGA-BH-A18P",
          "TCGA-BH-A1FG",
          "TCGA-AR-A255",
          "TCGA-AR-A1AX",
          "TCGA-E2-A1IH",
          "TCGA-BH-A0W7",
          "TCGA-E2-A156",
          "TCGA-BH-A0BQ",
          "TCGA-E2-A152",
          "TCGA-BH-A1EU",
          "TCGA-GM-A2DK",
          "TCGA-A1-A0SE",
          "TCGA-BH-A0DO",
          "TCGA-E2-A1LH",
          "TCGA-AR-A24S",
          "TCGA-A2-A0YI",
          "TCGA-BH-A0BO",
          "TCGA-GM-A2DI",
          "TCGA-A8-A06O",
          "TCGA-BH-A0C3",
          "TCGA-GM-A2D9",
          "TCGA-BH-A0DX",
          "TCGA-AR-A1AP",
          "TCGA-AR-A1AK",
          "TCGA-AR-A1AJ",
          "TCGA-OL-A66J",
          "TCGA-B6-A0X0",
          "TCGA-BH-A18K",
          "TCGA-AR-A1AY",
          "TCGA-EW-A1IY",
          "TCGA-BH-A0BR",
          "TCGA-A8-A08A",
          "TCGA-A2-A3XZ",
          "TCGA-BH-A0HY",
          "TCGA-B6-A40B",
          "TCGA-BH-A1ET",
          "TCGA-AR-A24N",
          "TCGA-E2-A154",
          "TCGA-BH-A0H5",
          "TCGA-EW-A1J6",
          "TCGA-BH-A0B6",
          "TCGA-AR-A2LR",
          "TCGA-A7-A0CD",
          "TCGA-GM-A2DO",
          "TCGA-A2-A0EP",
          "TCGA-BH-A1FD",
          "TCGA-BH-A0HA",
          "TCGA-AR-A2LE",
          "TCGA-AR-A252",
          "TCGA-A2-A259",
          "TCGA-E2-A1IO",
          "TCGA-B6-A1KI",
          "TCGA-A1-A0SB",
          "TCGA-BH-A0B0",
          "TCGA-A8-A095",
          "TCGA-BH-A0BW",
          "TCGA-E2-A15O",
          "TCGA-E2-A15F",
          "TCGA-E2-A1IN",
          "TCGA-E2-A14S",
          "TCGA-E2-A1IJ",
          "TCGA-AO-A03M",
          "TCGA-AR-A24P",
          "TCGA-BH-A0B8",
          "TCGA-BH-A0BL",
          "TCGA-E2-A14Z",
          "TCGA-E2-A15C",
          "TCGA-E2-A1IF",
          "TCGA-AO-A03V",
          "TCGA-B6-A402",
          "TCGA-BH-A0WA",
          "TCGA-BH-A0AV",
          "TCGA-A7-A3IY",
          "TCGA-BH-A0H3",
          "TCGA-E2-A1II",
          "TCGA-BH-A18S",
          "TCGA-GM-A2DH",
          "TCGA-AO-A03U",
          "TCGA-S3-AA14",
          "TCGA-E2-A15J",
          "TCGA-E2-A14U",
          "TCGA-BH-A0BG"
)


# Create and call the query object to fetch Type 1 Cancer dataset

query_brca_type1 <- GDCquery(project = "TCGA-BRCA",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       sample.type = c("Primary Tumor"),
                       barcode = bar1)


# Download the fetched dataset

GDCdownload(query_brca_type1)

# Obtain Samples from TCGA with Type 4 cancer condition

bar2 <- c("TCGA-A2-A0T2",
          "TCGA-B6-A0I9",
          "TCGA-UU-A93S",
          "TCGA-A8-A07W",
          "TCGA-A8-A08O",
          "TCGA-A2-A0CS",
          "TCGA-B6-A0IB",
          "TCGA-A2-A0SW",
          "TCGA-PL-A8LX",
          "TCGA-A8-A08T",
          "TCGA-BH-A18J",
          "TCGA-5L-AAT1",
          "TCGA-AO-A0J5",
          "TCGA-AC-A62V",
          "TCGA-A2-A0SV",
          "TCGA-A8-A08J",
          "TCGA-AN-A0FJ",
          "TCGA-LL-A73Z",
          "TCGA-B6-A3ZX",
          "TCGA-BH-A1FH"
)

# Create and call the query object to fetch Type 4 Cancer dataset

query_brca_type4 <- GDCquery(project = "TCGA-BRCA",
                             data.category = "Transcriptome Profiling",
                             data.type = "Gene Expression Quantification",
                             sample.type = c("Primary Tumor"),
                             barcode = bar2)

# Download the fetched dataset

GDCdownload(query_brca_type4)


# Prepare Files of GDC dataset to work with R environment 

query_brca_type1_prepare <- GDCprepare(query = query_brca_type1,
                                       save = TRUE,
                                       summarizedExperiment = TRUE,
                                       save.filename = "type1seq.rda")

query_brca_type4_prepare <- GDCprepare(query = query_brca_type4,
                                       save = TRUE,
                                       summarizedExperiment = TRUE,
                                       save.filename = "type4seq.rda")



# Data processing to perform correlation analysis to get unique genes

query_brca_type1_process <- TCGAanalyze_Preprocessing(object = query_brca_type1_prepare,
                                                      cor.cut = 0.7,
                                                      datatype = "unstranded",
                                                      filename = "type1.png")


query_brca_type4_process <- TCGAanalyze_Preprocessing(object = query_brca_type4_prepare,
                                                      cor.cut = 0.7,
                                                      datatype = "unstranded",
                                                      filename = "type4.png")


# Normalisation  of Gene counts data to maintain stability

dataNormalisation <- TCGAanalyze_Normalization(tabDF = cbind(query_brca_type1_process,query_brca_type4_process),
                                               geneInfo = geneInfoHT)

# Filter data with removing genes with no differences

dataFiltering <- TCGAanalyze_Filtering(tabDF = dataNormalisation,
                                       method = "quantile",
                                       qnt.cut = 0.05)


# Save normalised file

save(dataFiltering,file = paste0("full dataset normalised.rda"))


# Split datasets according to data groups

data_filter_type1 <- subset(dataFiltering, select = substr(colnames(dataFiltering),1,12) %in% bar1)

data_filter_type4 <- subset(dataFiltering, select = substr(colnames(dataFiltering),1,12) %in% bar2)

# Perform differential gene exp

diff_exp_genes <- TCGAanalyze_DEA(mat1 = data_filter_type1,
                                  mat2 = data_filter_type4,
                                  Cond1type = "Type 1",
                                  Cond2type = "Type 4",
                                  fdr.cut = 0.01,
                                  logFC.cut = 1,
                                  method = "glmLRT")
# Check number of DEGS

nrow(diff_exp_genes)

# Gene set enrichment analysis

gsea <- TCGAanalyze_EAcomplete(TFname = "DEA Analysis",
                               RegulonList = diff_exp_genes$gene_name)


# Display results in graphical plots

TCGAvisualize_EAbarplot(
  tf = rownames(gsea$ResBP),
  GOBPTab = gsea$ResBP,
  GOCCTab = gsea$ResCC,
  GOMFTab = gsea$ResMF,
  PathTab = gsea$ResPat,
  nRGTab = diff_exp_genes,
  nBar = 10
)

TCGAvisualize_EAbarplot(
  tf = rownames(gsea$ResBP),
  GOBPTab = gsea$ResBP,
  filename = NULL,
  nRGTab = rownames(diff_exp_genes),
  nBar = 10
)

TCGAvisualize_EAbarplot(
  tf = rownames(gsea$ResBP),
  GOCCTab = gsea$ResCC,
  filename = NULL,
  nRGTab = rownames(diff_exp_genes),
  nBar = 10
)

TCGAvisualize_EAbarplot(
  tf = rownames(gsea$ResBP),
  GOMFTab = gsea$ResMF,
  filename = NULL,
  nRGTab = rownames(diff_exp_genes),
  nBar = 10
)

TCGAvisualize_EAbarplot(
  tf = rownames(gsea$ResBP),
  PathTab = gsea$ResPat,
  filename = NULL,
  nRGTab = rownames(diff_exp_genes),
  nBar = 10
)


# Pathway analysis

library(clusterProfiler)
library(pathview)
library(SummarizedExperiment)

# DEG toptable

dataDEGsFiltLevel <- TCGAanalyze_LevelTab(FC_FDR_table_mRNA = diff_exp_genes,
                                          typeCond1 = "Type 1",
                                          typeCond2 = "Type 2",
                                          TableCond1 = dataFiltering[,colnames(data_filter_type1)],
                                          TableCond2 = dataFiltering[,colnames(data_filter_type4)])



dataDEGsFiltLevel$GeneID <- 0

# Converting gene symbol to Gene ID

eg = as.data.frame(
  bitr(dataDEGsFiltLevel$mRNA,
       fromType = "ENSEMBL",
       toType = c("ENTREZID","SYMBOL"),
       OrgDb = "org.Hs.eg.db")
)


eg <- eg[!duplicated(eg$SYMBOL),]
eg <- eg[order(eg$SYMBOL,decreasing=FALSE),]

dataDEGsFiltLevel <- dataDEGsFiltLevel[dataDEGsFiltLevel$mRNA %in% eg$ENSEMBL,]
dataDEGsFiltLevel <- dataDEGsFiltLevel[eg$ENSEMBL,]

rownames(dataDEGsFiltLevel) <- eg$SYMBOL

all(eg$SYMBOL == rownames(dataDEGsFiltLevel))

dataDEGsFiltLevel$GeneID <- eg$ENTREZID

dataDEGsFiltLevel_sub <- subset(dataDEGsFiltLevel,select = c("GeneID","logFC"))
genelistDEGs <- as.numeric(dataDEGsFiltLevel_sub$logFC)


names(genelistDEGs) <- dataDEGsFiltLevel_sub$GeneID



# Export gene ids

test <- names(genelistDEGs)
write.csv(test,"gene_list.csv")


# Extract the KEGG pathway using pathview which we got from david

hsa04080 <- pathview::pathview(
  gene.data = genelistDEGs,
  pathway.id = "hsa04080",
  species = "hsa",
  limit = list(gene=as.integer(max(abs(genelistDEGs))))
)
