library(DESeq2)
library(biomaRt)
library(dplyr)
library(PCAtools)
library(org.Mm.eg.db)
library(clusterProfiler)

# comparing to mice bec
# read in the counts file
mice_bec_cts <- read.csv('/Users/sophiemarcotte/Desktop/fragmentCountsGencode.txt', sep='\t')

# remove the dot and number at the end of geneid
mice_bec_cts$Geneid <- gsub("\\.\\d+$", "", mice_bec_cts$Geneid)

# subset to only include wanted columns
mice_bec_cts <- mice_bec_cts[, c(1, 7:12)]

# set geneid to row names
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
    "ensembl_gene_id",
    "external_gene_name",          # mouse gene symbol
    "drerio_homolog_ensembl_gene" # zebrafish ortholog gene symbol
  ),
  filters = "ensembl_gene_id",
  values = rownames(sig_bec_mice),
  mart = mouse
)

# clean 
orthologs <- orthologs %>% filter(!is.na(drerio_homolog_ensembl_gene) & drerio_homolog_ensembl_gene != "")
sig_bec_mice$ensembl_gene_id <- rownames(sig_bec_mice)
orthologs <- inner_join(orthologs, sig_bec_mice, by = 'ensembl_gene_id')

# merge with previous BEC results
becDeseqResults <- read.csv('/Users/sophiemarcotte/Desktop/patrice/estrogenRNAseq/deseqResultsBEC.csv', row.names=1)
becDeseqResults$ensembl_gene_id <- rownames(becDeseqResults)

# pull out rows that are present in other de
shared_becDeseq <- inner_join(orthologs, becDeseqResults, 
                              by = c('drerio_homolog_ensembl_gene' = 'ensembl_gene_id'),
                              suffix = c("_mice", "_zebrafish")) %>%
  dplyr::select(-2) %>%  # remove the second column
  rename(zebrafish.gene.name = Gene.name) %>%
  filter(padj_zebrafish < 0.05)

write.csv(shared_becDeseq, "/Users/sophiemarcotte/Desktop/patrice/estrogenRNAseq/shared_DE_bec_pregMice_and_E2zebrafish.csv", row.names = FALSE)

# ----------------------------- kegg and go ----------------------------- #
# pull out upregulated gene ids
gene_bec_names <- rownames(sig_bec_mice)

# convert to entrez IDs
entrez_ids_bec <- mapIds(org.Mm.eg.db,
                            keys = gene_bec_names,
                            column = "ENTREZID",
                            keytype = "ENSEMBL",
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
dotplot(go_results_bec, showCategory = 15) + ggtitle("Mice Pregnancy GO for Sig DE Genes")

# KEGG enrichment
kegg_bec <- enrichKEGG(gene = entrez_ids_bec,
                          organism = "mmu",
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05)

# plot
dotplot(kegg_bec, showCategory = 15) + ggtitle("Mice Pregnancy KEGG for Sig DE Genes")

# order results
kegg_results <- kegg_bec@result[order(kegg_bec@result$p.adjust), ]
go_results <- go_results_bec@result[order(go_results_bec@result$p.adjust), ]

write.csv(kegg_results, "/Users/sophiemarcotte/Desktop/patrice/estrogenRNAseq/kegg_pregMiceBEC.csv", row.names = TRUE)
write.csv(go_results, "/Users/sophiemarcotte/Desktop/patrice/estrogenRNAseq/go_pregMiceBEC.csv", row.names = TRUE)

# ------------------ create a plot showing the correlation between the two DE gene sets -----------------------#
shared_becDeseq <- read.csv("/Users/sm2949/Desktop/shared_DE_bec_pregMice_and_E2zebrafish.csv") %>% 
  mutate(zebrafish.gene.name = tolower(zebrafish.gene.name))

# create a lm between the two lfc 
model <- lm(log2FoldChange_mice ~ log2FoldChange_zebrafish, data = shared_becDeseq)
r2 <- summary(model)$r.squared

# plot
ggplot(shared_becDeseq, aes(x = log2FoldChange_mice, y = log2FoldChange_zebrafish)) +
  geom_point(alpha = 0.6) +
  geom_point(
    data = subset(shared_becDeseq, zebrafish.gene.name == "sox17"),
    color = "red",
    size = 2
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +  # vertical line
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +  # horizontal line
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  geom_text_repel(
    data = subset(shared_becDeseq, zebrafish.gene.name == "sox17"),
    aes(label = zebrafish.gene.name),
    color = "red",
    size = 4
  ) +
  labs(
    x = "log2FoldChange (Mice)",
    y = "log2FoldChange (Zebrafish)",
    title = "log2FC of Shared DE Genes: E2-Treated Zebrafish BECs and Pregnant Mice BECs"
  ) +
  annotate("text", x = Inf, y = -Inf, 
           label = paste0("R² = ", round(r2, 3)), 
           hjust = 1.1, vjust = -0.5, size = 5) +
  theme_minimal()

