# import required libraries
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(PCAtools)
library(dplyr)

# ------------ read in the counts and separate to HEP ------------- #
raw_counts <- read.csv('cleanedEstrogenCounts.csv', row.names = 1)

# split the dataframe into HEP and BIL
bec_cols <- grep("BEC", colnames(raw_counts), value = TRUE)
hep_cols <- grep("HEP", colnames(raw_counts), value = TRUE)

# create new dataframes
bec_counts <- raw_counts[, bec_cols]
hep_counts <- raw_counts[, hep_cols]

# ------------ DESEQ2 data set up for HEPs ------------- #

# get the metadata information from the col names
colnames_hep <- colnames(hep_counts)

# split by . and extract condition and replicate
coldata_hep <- data.frame(
  condition = sapply(strsplit(colnames_hep, "\\."), `[`, 1),
  replicate = sapply(strsplit(colnames_hep, "\\."), `[`, 3)
)
rownames(coldata_hep) <- colnames_hep

# reorder to match counts col names
coldata_hep <- coldata_hep[colnames(hep_counts), ]
all(rownames(coldata_hep) == colnames(hep_counts))  # assert

# set condition and replicate as a factor
coldata_hep$condition <- as.factor(coldata_hep$condition)
coldata_hep$replicate <- as.factor(coldata_hep$replicate)

# ----------- running DESEQ2 on HEPs ----------------- #

# create the deseq object
dds_hep <- DESeqDataSetFromMatrix(countData = hep_counts,
                                  colData = coldata_hep,
                                  design = ~ replicate + condition)

# pre-filter
smallestGroupSize <- 4
keep <- rowSums(counts(dds_hep) >= 10) >= smallestGroupSize
dds_hep <- dds_hep[keep,]

# set the reference level
dds_hep$condition <- relevel(dds_hep$condition, ref = "EtOH")

# run DE
dds_hep <- DESeq(dds_hep)

# shrinkage
resLFChep <- lfcShrink(dds_hep, coef="condition_E2_vs_EtOH", type="apeglm")
resLFChep_df <- as.data.frame(resLFChep)

# extract significant results
sig_hep <- resLFChep_df[
  resLFChep_df$padj < 0.05 &
  !is.na(resLFChep_df$padj) &
  abs(resLFChep_df$log2FoldChange) > 1,
]

# VST counts & PCA plot
vst <- assay(vst(dds_hep))
p <- pca(vst, metadata = colData(dds_hep), removeVar = 0.1)
biplot(p)

# PCA plot
vsd <- vst(dds_hep)
plotPCA(vsd, intgroup=c("condition", "replicate"))

# remove 'batch' effects - only run if doing downstream analyses
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$replicate)

# cooks distance
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds_hep)[["cooks"]]), range=0, las=2)        

# sample distance heatmap
sampleDists <- dist(t(assay(vst(dds_hep))))
sampleDistMatrix <- as.matrix(sampleDists)
pheatmap::pheatmap(sampleDistMatrix, labels_row=colnames(dds_hep))

# ----------- running DESEQ2 on HEPs (excluding rep 1 samples) ----------------- #

# exclude rep 1 from counts and coldata
hep_counts_norep1 <- hep_counts[, !grepl("\\.1$", colnames(hep_counts))]
coldata_hep_norep1 <- coldata_hep[coldata_hep$replicate != 1, ]

# ensure coldata are factors
coldata_hep_norep1$condition <- as.factor(coldata_hep_norep1$condition)
coldata_hep_norep1$replicate <- as.factor(coldata_hep_norep1$replicate)

# create the deseq object
dds_hep_norep1 <- DESeqDataSetFromMatrix(countData = hep_counts_norep1,
                                  colData = coldata_hep_norep1,
                                  design = ~ replicate + condition)

# pre-filter
smallestGroupSize <- 3
keep <- rowSums(counts(dds_hep_norep1) >= 10) >= smallestGroupSize
dds_hep_norep1 <- dds_hep_norep1[keep,]

# set the reference level
dds_hep_norep1$condition <- relevel(dds_hep_norep1$condition, ref = "EtOH")

# run DE
dds_hep_norep1 <- DESeq(dds_hep_norep1)

# shrinkage
resLFChep_norep1 <- lfcShrink(dds_hep_norep1, coef="condition_E2_vs_EtOH", type="apeglm")
resLFChep_df_norep1 <- as.data.frame(resLFChep_norep1)

# extract significant results
sig_hep_norep1 <- resLFChep_df_norep1[
  resLFChep_df_norep1$padj < 0.05 &
    !is.na(resLFChep_df_norep1$padj) &
    abs(resLFChep_df_norep1$log2FoldChange) > 1,
]

# vst
vst <- assay(vst(dds_hep_norep1))
write.csv(vst, file = "/Users/sm2949/Desktop/vstCountsHEP.csv", row.names = TRUE)

# ------------------------ add symbol and human orthologs for HEPs ------------------------ #
# read in conversion file
conversion <- read.csv('/Users/sm2949/Desktop/mart_export-2.txt', sep='\t')
# filter for unique
conversion_unique <- conversion[!duplicated(conversion$Gene.stable.ID), ]

# add ensembl as a column
resLFChep_df$Gene.stable.ID <- rownames(resLFChep_df)
# merge with conversion info
resLFChep_final <- merge(resLFChep_df, conversion_unique, by = "Gene.stable.ID", all.x = TRUE)

# add ensembl as a column
resLFChep_df_norep1$Gene.stable.ID <- rownames(resLFChep_df_norep1)
# merge with conversion info
resLFChep_final_norep1 <- left_join(resLFChep_df_norep1, conversion_unique, by = "Gene.stable.ID")
# set ensembl id as row names
rownames(resLFChep_final_norep1) <- resLFChep_final_norep1[[6]]   
resLFChep_final_norep1 <- resLFChep_final_norep1[, -6]
# save
write.csv(resLFChep_final_norep1, file = "/Users/sm2949/Desktop/deseqResultsHEP.csv", row.names = TRUE)



# ------------------------ volcano plots ------------------------ #
# Define thresholds for labeling
resLFChep_final$label <- ifelse(resLFChep_final$padj < 0.05 & abs(resLFChep_final$log2FoldChange) > 1.5,
                             resLFChep_final$Gene.name,
                             NA)

# plot
ggplot(resLFChep_final, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05 & abs(log2FoldChange) > 1), alpha = 0.7) +
  scale_color_manual(values = c("grey", "red")) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 15) +
  theme_minimal() +
  labs(title = "Volcano Plot",
       x = "log2 Fold Change",
       y = "-log10 Adjusted p-value",
       color = "Significant") +
  theme(legend.position = "bottom")

# Define thresholds for labeling
resLFChep_final_norep1$label <- ifelse(resLFChep_final_norep1$padj < 0.05 & abs(resLFChep_final_norep1$log2FoldChange) > 1.5,
                                    resLFChep_final_norep1$Gene.name,
                                    NA)

# plot
EnhancedVolcano(resLFChep_final_norep1,
                lab = resLFChep_final_norep1$label,  
                x = 'log2FoldChange',
                y = 'padj',
                title = 'BEC E2 vs EtOH',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2,
                labSize = 3.5,
                col = c('grey70', 'grey70', '#0072B2', '#D73027'),  
                legendLabels = c('NS', 'Log2FC', 'padj', 'padj & Log2FC'),
                legendPosition = 'bottom',
                legendLabSize = 10,
                legendIconSize = 3.5,
                drawConnectors = TRUE,
                widthConnectors = 0.4,
                colConnectors = 'grey50',
                max.overlaps = 20
)

