# import required libraries
library(DESeq2)

# ------------ clean the data from featureCounts output ------------- #

# read in the counts file
raw_counts <- read.table("/Users/sophiemarcotte/Desktop/estrogenCounts.txt", sep = "\t", header = TRUE, comment.char = "#")

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
write.csv(raw_counts, file = "cleanedEstrogenCounts.csv", row.names = TRUE)

# split the dataframe into HEP and BIL
bec_cols <- grep("BEC", colnames(raw_counts), value = TRUE)
hep_cols <- grep("HEP", colnames(raw_counts), value = TRUE)

# create new dataframes
bec_counts <- raw_counts[, bec_cols]
hep_counts <- raw_counts[, hep_cols]

# ------------ DESEQ2 data set up for BECs ------------- #

# get the metadata information from the col names
colnames_bec <- colnames(bec_counts)

# split by '.' and extract condition and replicate
coldata_bec <- data.frame(
  condition = sapply(strsplit(colnames_bec, "\\."), `[`, 1),
  replicate = sapply(strsplit(colnames_bec, "\\."), `[`, 3)
)
rownames(coldata_bec) <- colnames_bec

# reorder to match counts col names
coldata_bec <- coldata_bec[colnames(bec_counts), ]
all(rownames(coldata_bec) == colnames(bec_counts))  # assert

# set condition as a factor
coldata_bec$condition <- as.factor(coldata_bec$condition)

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

# VST counts
vsd <- vst(dds_bec) 

# PCA plot
plotPCA(vsd, intgroup=c("condition", "replicate"))

# remove 'batch' effects
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$replicate)

plotPCA(vsd, intgroup=c("condition", "replicate"))

par(mar=c(8,5,2,2))
boxplot(log10(assays(dds_bec)[["cooks"]]), range=0, las=2)        
