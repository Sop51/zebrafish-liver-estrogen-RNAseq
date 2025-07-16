library(tidyverse)
library(DESeq2)

# read in the normalized counts for each dataset
hep_norm_counts <- read.csv('/Users/sm2949/Desktop/patrice/estrogenRNAseq/vstCountsHEP.csv', row.names = 1)
bec_norm_counts <- read.csv('/Users/sm2949/Desktop/patrice/estrogenRNAseq/vstCountsBEC.csv', row.names = 1)

# read in the raw counts for each dataset
hep_raw_counts <- read.csv('/Users/sm2949/Desktop/patrice/estrogenRNAseq/rawCountsHEP.csv', row.names = 1)
bec_raw_counts <- read.csv('/Users/sm2949/Desktop/patrice/estrogenRNAseq/rawCountsBEC.csv', row.names = 1)

# prepare the count matrix for deseq ------------------------
# pull out only control samples
etoh_bec_raw <- bec_raw_counts[, grepl("EtOH", colnames(bec_raw_counts))]
etoh_hep_raw <- hep_raw_counts[, grepl("EtOH", colnames(hep_raw_counts))]

# merge
etoh_combined_raw <- merge(etoh_hep_raw, etoh_bec_raw, 
                           by = "row.names", all = TRUE)

# set gene ID as rownames again
rownames(etoh_combined_raw) <- etoh_combined_raw$Row.names
etoh_combined_raw$Row.names <- NULL

# create the col data for deseq ----------------------------
colnames_combined <- colnames(etoh_combined_raw)

# split by . and extract condition and replicate
coldata_combined <- data.frame(
  cell_type = sapply(strsplit(colnames_combined, "\\."), `[`, 2),
  replicate = sapply(strsplit(colnames_combined, "\\."), `[`, 3)
)
rownames(coldata_combined) <- colnames_combined

# reorder to match counts col names
coldata_combined <- coldata_combined[colnames(etoh_combined_raw), ]
all(rownames(coldata_combined) == colnames(etoh_combined_raw))  # assert

# set cell type and replicate as a factor
coldata_combined$cell_type <- as.factor(coldata_combined$cell_type)
coldata_combined$replicate <- as.factor(coldata_combined$replicate)

# running deseq ----------------------------------------
# create the deseq object
dds_comb <- DESeqDataSetFromMatrix(countData = etoh_combined_raw,
                                  colData = coldata_combined,
                                  design = ~ replicate + cell_type)

# pre-filter
smallestGroupSize <- 3
keep <- rowSums(counts(dds_comb) >= 10) >= smallestGroupSize
dds_comb <- dds_comb[keep,]

# set the reference level
dds_comb$cell_type <- relevel(dds_comb$cell_type, ref = "HEP")

# run DE
dds_comb <- DESeq(dds_comb)

# shrinkage
resLFCcomb <- lfcShrink(dds_comb, coef="cell_type_BEC_vs_HEP", type="apeglm")
resLFCcomb_df <- as.data.frame(resLFCcomb)

# extract significant results
sig_comb <- resLFCcomb_df[
  resLFCcomb_df$padj < 0.05 &
    !is.na(resLFCcomb_df$padj) &
    abs(resLFCcomb_df$log2FoldChange) > 1,
]

# read in conversion file to add 
conversion <- read.csv('/Users/sm2949/Desktop/mart_export-2.txt', sep='\t')
# filter for unique
conversion_unique <- conversion[!duplicated(conversion$Gene.stable.ID), ]

# add ensembl as a column
resLFCcomb_df$Gene.stable.ID <- rownames(resLFCcomb_df)
# merge with conversion info
resLFCcomb_final <- merge(resLFCcomb_df, conversion_unique, by = "Gene.stable.ID", all.x = TRUE)

# ------------------------------- create plots ------------------------------------

# hepatocyte genes to plot
hep_ensembl_markers <- c("ENSDARG00000042780", "ENSDARG00000015866", 
                       "ENSDARG00000008969", "ENSDARG00000090850")

hep_gene_markers <- c("apoba", "apoa2", "fgb", "serpina1l")

# bec genes to plot
bec_ensembl_markers <- c("ENSDARG00000036456", "ENSDARG00000043923",
                         "ENSDARG00000028618", "ENSDARG00000026531")


bec_gene_markers <- c("anxa4", "sox9b", "krt18b", "alcama")

# pull out only control samples
etoh_bec_samples <- bec_norm_counts[, grepl("EtOH", colnames(bec_norm_counts))]
etoh_hep_samples <- hep_norm_counts[, grepl("EtOH", colnames(hep_norm_counts))]

# subset to marker genes
hep_bec_genes <- union(hep_ensembl_markers, bec_ensembl_markers)

# create a lookup for gene labels
gene_lookup <- data.frame(
  ensembl = c(hep_ensembl_markers, bec_ensembl_markers),
  gene = c(hep_gene_markers, bec_gene_markers),
  stringsAsFactors = FALSE
) %>% mutate(label = gene)

# define function to prep long data for either gene set 
prep_plot_data <- function(gene_ids, gene_lookup) {
  hep_sub <- etoh_hep_samples[gene_ids, ] %>%
    as.data.frame() %>%
    mutate(ensembl = rownames(.), Identity = "Hepatocyte")
  
  bec_sub <- etoh_bec_samples[gene_ids, ] %>%
    as.data.frame() %>%
    mutate(ensembl = rownames(.), Identity = "BEC")
  
  combined <- bind_rows(hep_sub, bec_sub)
  
  # reshape
  long_df <- combined %>%
    pivot_longer(cols = -c(ensembl, Identity), names_to = "sample", values_to = "vst") %>%
    left_join(gene_lookup, by = "ensembl")
  
  return(long_df)
}

# hepatocyte marker plot
hep_data_long <- prep_plot_data(hep_ensembl_markers, gene_lookup)

ggplot(hep_data_long, aes(x = label, y = vst, fill = Identity)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(width = 0.7)) +
  geom_jitter(aes(color = Identity), size = 1, position = position_dodge(width = 0.7)) +
  labs(title = "Hepatocyte Marker Gene Expression",
       x = "Gene", y = "Normalized Expression Level (VST)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Hepatocyte" = "#413C58", "BEC" = "#79B86C")) +
  scale_color_manual(values = c("Hepatocyte" = "#413C58", "BEC" = "#79B86C"))


# BEC marker plot 
bec_data_long <- prep_plot_data(bec_ensembl_markers, gene_lookup)

ggplot(bec_data_long, aes(x = label, y = vst, fill = Identity)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(width = 0.7)) +
  geom_jitter(aes(color = Identity), size = 1, position = position_dodge(width = 0.7)) +
  labs(title = "BEC Marker Gene Expression",
       x = "Gene", y = "Normalized Expression Level (VST)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Hepatocyte" = "#413C58", "BEC" = "#79B86C")) +
  scale_color_manual(values = c("Hepatocyte" = "#413C58", "BEC" = "#79B86C"))




