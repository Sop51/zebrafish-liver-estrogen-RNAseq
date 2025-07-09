# import required libraries
library(DESeq2)
library(ggplot2)
library(ggrepel)  
library(biomaRt)
library(EnhancedVolcano)

# ------------ clean the data from featureCounts output ------------- #

# read in the counts file
raw_counts <- read.table("/Users/sm2949/Desktop/estrogenCounts.txt", sep = "\t", header = TRUE, comment.char = "#")

# clean the column names
colnames(raw_counts) <- gsub("^X\\.data\\.mashed\\.liver\\.sm2949\\.estrogenRNAseqResults\\.aligned_bams\\.|\\.Aligned\\.sortedByCoord\\.out\\.bam$", "", colnames(raw_counts))

# remove transcript versions
raw_counts$Geneid <- sub("\\.\\d+$", "", raw_counts$Geneid)

# remove unwanted columns
raw_counts <- raw_counts[,c(1,7:22)]

# set Geneid as row names
rownames(raw_counts) <- raw_counts[[1]]   
raw_counts <- raw_counts[, -1]       

# save the cleaned count matrix as a csv
write.csv(raw_counts, file = "/Users/sm2949/Desktop/cleanedEstrogenCounts.csv", row.names = TRUE)

# split the dataframe into HEP and BIL
bec_cols <- grep("BEC", colnames(raw_counts), value = TRUE)
hep_cols <- grep("HEP", colnames(raw_counts), value = TRUE)

# create new dataframes
bec_counts <- raw_counts[, bec_cols]
hep_counts <- raw_counts[, hep_cols]

# save the cleaned counts matricies as csvs
write.csv(bec_counts, file = "/Users/sm2949/Desktop/rawCountsBEC.csv", row.names = TRUE)
write.csv(hep_counts, file = "/Users/sm2949/Desktop/rawCountsHEP.csv", row.names = TRUE)

# ------------ DESEQ2 data set up for BECs ------------- #

# get the metadata information from the col names
colnames_bec <- colnames(bec_counts)

# split by . and extract condition and replicate
coldata_bec <- data.frame(
  condition = sapply(strsplit(colnames_bec, "\\."), `[`, 1),
  replicate = sapply(strsplit(colnames_bec, "\\."), `[`, 3)
)
rownames(coldata_bec) <- colnames_bec

# reorder to match counts col names
coldata_bec <- coldata_bec[colnames(bec_counts), ]
all(rownames(coldata_bec) == colnames(bec_counts))  # assert

# set condition and replicate as a factor
coldata_bec$condition <- as.factor(coldata_bec$condition)
coldata_bec$replicate <- as.factor(coldata_bec$replicate)

# ----------- running DESEQ2 on BECs ----------------- #

# create the deseq object
dds_bec <- DESeqDataSetFromMatrix(countData = bec_counts,
                              colData = coldata_bec,
                              design = ~ replicate + condition)

# pre-filter
smallestGroupSize <- 4
keep <- rowSums(counts(dds_bec) >= 10) >= smallestGroupSize
dds_bec <- dds_bec[keep,]

# set the reference level
dds_bec$condition <- relevel(dds_bec$condition, ref = "EtOH")

# run DE
dds_bec <- DESeq(dds_bec)

# shrinkage
resLFC <- lfcShrink(dds_bec, coef="condition_E2_vs_EtOH", type="apeglm")
resLFC_df <- as.data.frame(resLFC)

# filter for sig results
sig_bec <- resLFC_df[
  resLFC_df$padj < 0.05 &
    abs(resLFC_df$log2FoldChange) > 1,
]

# VST counts & PCA plot
vst <- assay(vst(dds_bec))
write.csv(vst, file = "/Users/sm2949/Desktop/patrice/estrogenRNAseq/vstCountsBEC.csv", row.names = TRUE)
p <- pca(vst, metadata = colData(dds_bec), removeVar = 0.1)
biplot(p)

# PCA plot
vsd <- vst(dds_bec) 
plotPCA(vsd, intgroup=c("condition", "replicate"))

# cooks distance
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds_bec)[["cooks"]]), range=0, las=2)        

# sample distance heatmap
sampleDists <- dist(t(assay(vst(dds_bec))))
sampleDistMatrix <- as.matrix(sampleDists)
pheatmap::pheatmap(sampleDistMatrix, labels_row=colnames(dds_bec))

# ----------- running DESEQ2 on BECs (excluding rep 1 samples) ----------------- #

# remove sample one from count data and coldata
bec_counts_norep1 <- bec_counts[, !grepl("\\.1$", colnames(bec_counts))]
coldata_bec_norep1 <- coldata_bec[coldata_bec$replicate != 1, ]

# create the deseq object
dds_bec_norep1 <- DESeqDataSetFromMatrix(countData = bec_counts_norep1,
                                  colData = coldata_bec_norep1,
                                  design = ~ replicate + condition)

# pre-filter
smallestGroupSize <- 3
keep <- rowSums(counts(dds_bec_norep1) >= 10) >= smallestGroupSize
dds_bec_norep1 <- dds_bec_norep1[keep,]

# set the reference level
dds_bec_norep1$condition <- relevel(dds_bec_norep1$condition, ref = "EtOH")

# run DE
dds_bec_norep1 <- DESeq(dds_bec_norep1)

# shrinkage
resLFC_norep1 <- lfcShrink(dds_bec_norep1, coef="condition_E2_vs_EtOH", type="apeglm")
resLFC_df_norep1 <- as.data.frame(resLFC_norep1)

# filter for significant results
sig_bec_norep1 <- resLFC_df_norep1[
  resLFC_df_norep1$padj < 0.05 &
    abs(resLFC_df_norep1$log2FoldChange) > 1,
]

# VST counts
vsd <- vst(dds_bec_norep1) 

# PCA plot
plotPCA(vsd, intgroup=c("condition", "replicate"))

# remove batch effects - only run if doing downstream analyses
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$replicate)

# cooks distance
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds_bec_norep1)[["cooks"]]), range=0, las=2)        

# sample distance heatmap
sampleDists <- dist(t(assay(vst(dds_bec_norep1))))
sampleDistMatrix <- as.matrix(sampleDists)
pheatmap::pheatmap(sampleDistMatrix, labels_row=colnames(dds_bec_norep1))

# ------------------------ add symbol and human orthologs for BECs ------------------------ #
# read in conversion file
conversion <- read.csv('/Users/sm2949/Desktop/mart_export-2.txt', sep='\t')
# filter for unique
conversion_unique <- conversion[!duplicated(conversion$Gene.stable.ID), ]

# add ensembl as a column
resLFC_df$Gene.stable.ID <- rownames(resLFC_df)
# merge with conversion info
resLFC_final <- left_join(resLFC_df, conversion_unique, by = "Gene.stable.ID")
# set ensembl id as row names
rownames(resLFC_final) <- resLFC_final[[6]]   
resLFC_final <- resLFC_final[, -6]
# save
write.csv(resLFC_final, file = "/Users/sm2949/Desktop/deseqResultsBEC.csv", row.names = TRUE)

# add ensembl as a column
resLFC_df_norep1$Gene.stable.ID <- rownames(resLFC_df_norep1)
# merge with conversion info
resLFC_final_norep1 <- merge(resLFC_df_norep1, conversion_unique, by = "Gene.stable.ID", all.x = TRUE)

# ------------------------ volcano plots ------------------------ #
# define thresholds for labeling
resLFC_final$label <- ifelse(resLFC_final$padj < 0.05 & abs(resLFC_final$log2FoldChange) > 1.5,
                             resLFC_final$Gene.name,
                            NA)

# plot
EnhancedVolcano(resLFC_final,
                lab = resLFC_final$label,  
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
                max.overlaps = 40
)



# define thresholds for labeling
resLFC_final_norep1$label <- ifelse(resLFC_final_norep1$padj < 0.05 & abs(resLFC_final_norep1$log2FoldChange) > 1.5,
                             resLFC_final_norep1$Gene.name,
                             NA)

# plot
ggplot(resLFC_final_norep1, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05 & abs(log2FoldChange) > 1), alpha = 0.7) +
  scale_color_manual(values = c("grey", "red")) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 15) +
  theme_minimal() +
  labs(title = "Volcano Plot",
       x = "log2 Fold Change",
       y = "-log10 Adjusted p-value",
       color = "Significant") +
  theme(legend.position = "bottom")
