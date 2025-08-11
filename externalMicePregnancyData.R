library(DESeq2)
library(biomaRt)
library(dplyr)
library(PCAtools)
library(org.Mm.eg.db)
library(clusterProfiler)

# comparing to mice bec
# read in the counts file
mice_bec_cts <- read.csv('/Users/sm2949/Desktop/micePregnancyRawCounts.txt', sep='\t')

# remove the dot and number at the end of geneid
mice_bec_cts$Geneid <- gsub("\\.\\d+$", "", mice_bec_cts$Geneid)

# subset to only include wanted columns
mice_bec_cts <- mice_bec_cts[, c(1, 7:12)]

# set geneid to row names
mice_bec_cts <- na.omit(mice_bec_cts)
rownames(mice_bec_cts) <- mice_bec_cts$Geneid
mice_bec_cts <- mice_bec_cts[ , -1]

# create the col data for deseq
coldata <- data.frame(
  replicate = sub(".*_(rep[0-9]+)", "\\1", colnames(mice_bec_cts)),
  status = sub("_rep[0-9]+", "", colnames(mice_bec_cts)),
  row.names = colnames(mice_bec_cts)
)

# set status and replicate as a factor
coldata$status <- as.factor(coldata$status)
coldata$replicate <- as.factor(coldata$replicate)

# create the deseq object
dds_bec_mice <- DESeqDataSetFromMatrix(countData = mice_bec_cts,
                                  colData = coldata,
                                  design = ~ status)

# pre-filter
smallestGroupSize <- 3
keep <- rowSums(counts(dds_bec_mice) >= 10) >= smallestGroupSize
dds_bec_mice <- dds_bec_mice[keep,]

# set the reference level
dds_bec_mice$status <- relevel(dds_bec_mice$status, ref = "ctrl")

# run DE
dds_bec_mice <- DESeq(dds_bec_mice)

# shrinkage
resLFC <- lfcShrink(dds_bec_mice, coef="status_preg_vs_ctrl", type="apeglm")
resLFC_df <- as.data.frame(resLFC)

# filter for sig results
sig_bec_mice <- resLFC_df %>%
  filter(!is.na(padj) & padj < 0.05)

# VST counts & PCA plot
vst <- assay(vst(dds_bec_mice))
p <- pca(vst, metadata = colData(dds_bec_mice), removeVar = 0.1)
biplot(p)

# Use Ensembl marts for mouse and zebrafish
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
zebrafish = useMart("ensembl", dataset = "drerio_gene_ensembl")

orthologs <- getBM(
  attributes = c(
    "external_gene_name",          # mouse gene symbol
    "drerio_homolog_associated_gene_name" # zebrafish ortholog gene symbol
  ),
  filters = "external_gene_name",
  values = rownames(sig_bec_mice),
  mart = mouse
)

# clean 
orthologs <- orthologs %>% filter(!is.na(drerio_homolog_associated_gene_name) & drerio_homolog_associated_gene_name != "")
sig_bec_mice$external_gene_name <- rownames(sig_bec_mice)
orthologs <- left_join(orthologs, sig_bec_mice, by = 'external_gene_name')

# merge with previous BEC results
becDeseqResults <- read.csv('/Users/sm2949/Desktop/patrice/estrogenRNAseq/deseqResultsBEC.csv', row.names=1)
becDeseqResults$Gene.name <- tolower(becDeseqResults$Gene.name)

# pull out rows that are present in other de
shared_becDeseq <- inner_join(orthologs, becDeseqResults, by = c('drerio_homolog_associated_gene_name' = 'Gene.name'))

# ----------------------------- kegg and go ----------------------------- #
# pull out upregulated gene ids
gene_bec_names <- rownames(sig_bec_mice)

# convert to entrez IDs
entrez_ids_bec <- mapIds(org.Mm.eg.db,
                            keys = gene_bec_names,
                            column = "ENTREZID",
                            keytype = "SYMBOL",
                            multiVals = "first")

# run go
go_results_bec <- enrichGO(gene = entrez_ids_bec,
                              OrgDb = org.Mm.eg.db,
                              keyType = "ENTREZID",
                              ont = "BP", 
                              pAdjustMethod = "BH",
                              qvalueCutoff = 0.05,
                              readable = TRUE)
# plot
dotplot(go_results_bec, showCategory = 15) + ggtitle("E2 Upregulated GO Biological Process")

# KEGG enrichment
kegg_bec <- enrichKEGG(gene = entrez_ids_bec,
                          organism = "mmu",
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05)

# plot
dotplot(kegg_bec, showCategory = 15) + ggtitle("E2 Upregulated KEGG")



