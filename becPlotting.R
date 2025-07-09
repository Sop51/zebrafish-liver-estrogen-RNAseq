library(pheatmap)
library(clusterProfiler)
library(org.Dr.eg.db)

# ------------------ set up data for analyses --------------------- #
# read in the normalized counts & deseq results
becNormalizedCounts <- read.csv('/Users/sm2949/Desktop/patrice/estrogenRNAseq/vstCountsBEC.csv', row.names=1)
becDeseqResults <- read.csv('/Users/sm2949/Desktop/patrice/estrogenRNAseq/deseqResultsBEC.csv', row.names=1)

# if no gene symbol, replace with ensembl
becDeseqResults$Gene.name <- ifelse(
  is.na(becDeseqResults$Gene.name) | becDeseqResults$Gene.name == "",
  rownames(becDeseqResults),
  becDeseqResults$Gene.name
)

# add gene symbols to normalized counts
becNormalizedCounts$Gene.name <- becDeseqResults$Gene.name[match(rownames(becNormalizedCounts), rownames(becDeseqResults))]

# pull out only significant genes
sig_bec_genes <- becDeseqResults[!is.na(becDeseqResults$padj) & becDeseqResults$padj < 0.05, ]

# ---------------------- plotting a heatmap of up & downreg genes ------------------- #

# sort into up and down regulated
up_bec_genes <- sig_bec_genes[order(-sig_bec_genes$log2FoldChange), ][1:25, ]
down_bec_genes <- sig_bec_genes[order(sig_bec_genes$log2FoldChange), ][1:25, ]

# combine and get gene names
top_bec_de_genes <- unique(c(rownames(up_bec_genes), rownames(down_bec_genes)))

# subset the normalized matrix
heatmap_bec <- becNormalizedCounts[rownames(becNormalizedCounts) %in% top_bec_de_genes, ]

# set gene symbols as row names
rownames(heatmap_bec) <- heatmap_bec[[9]]   
heatmap_bec <- heatmap_bec[, -9]

# z-score normalize
heatmap_bec <- t(scale(t(heatmap_bec)))

# plot
pheatmap(heatmap_bec,
         show_rownames = TRUE,
         show_colnames = TRUE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_row = 8,
         fontsize_col = 10,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

# ----------------------------- kegg and go ----------------------------- #
# pull out upregulated gene ids
up_bec_genes <- sig_bec_genes[sig_bec_genes$log2FoldChange > 1,]
up_gene_bec_names <- rownames(up_bec_genes)

# convert to entrez IDs
entrez_ids_bec_up <- mapIds(org.Dr.eg.db,
                     keys = up_gene_bec_names,
                     column = "ENTREZID",
                     keytype = "ENSEMBL",
                     multiVals = "first")

# run go
go_results_bec_up <- enrichGO(gene = entrez_ids_bec_up,
                       OrgDb = org.Dr.eg.db,
                       keyType = "ENTREZID",
                       ont = "MF", 
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable = TRUE)
# plot
dotplot(go_results_bec_up, showCategory = 15) + ggtitle("E2 Upregulated GO Biological Process")

# KEGG enrichment
kegg_bec_up <- enrichKEGG(gene = entrez_ids_bec_up,
                      organism = "dre",
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05)

# plot
dotplot(kegg_bec_up, showCategory = 15) + ggtitle("E2 Upregulated KEGG")


# pull out downregulated gene ids
down_bec_genes <- sig_bec_genes[sig_bec_genes$log2FoldChange < 1,]
down_gene_bec_names <- rownames(down_bec_genes)

# convert to entrez IDs
entrez_ids_bec_down <- mapIds(org.Dr.eg.db,
                            keys = down_gene_bec_names,
                            column = "ENTREZID",
                            keytype = "ENSEMBL",
                            multiVals = "first")

# run go
down_results_bec_up <- enrichGO(gene = entrez_ids_bec_down,
                              OrgDb = org.Dr.eg.db,
                              keyType = "ENTREZID",
                              ont = "MF", 
                              pAdjustMethod = "BH",
                              qvalueCutoff = 0.05,
                              readable = TRUE)
# plot
dotplot(down_results_bec_up, showCategory = 15) + 
  ggtitle("E2 Downregulated GO Biological Process") +
  theme(axis.text.y = element_text(size = 8, face = "plain"))

# KEGG enrichment
kegg_bec_down <- enrichKEGG(gene = entrez_ids_bec_down,
                          organism = "dre",
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05)

# plot
dotplot(kegg_bec_down, showCategory = 15) + ggtitle("E2 Downregulated KEGG")

# ------------------ bar plot of up, down, ns gene count ---------------- #

# categorize genes
becDeseqResults$regulation <- ifelse(
  becDeseqResults$padj < 0.05 & becDeseqResults$log2FoldChange > 1, "Upregulated",
  ifelse(becDeseqResults$padj < 0.05 & becDeseqResults$log2FoldChange < 1, "Downregulated", "Not Significant")
)

# count how many genes in each category
reg_counts <- table(becDeseqResults$regulation)
reg_df <- as.data.frame(reg_counts)
colnames(reg_df) <- c("Regulation", "Count")

# plot
ggplot(reg_df, aes(x = Regulation, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = Count), vjust = -0.3, size = 5) +
  scale_fill_manual(values = c("Upregulated" = "#D73027", 
                               "Downregulated" = "#4575B4", 
                               "Not Significant" = "grey70")) +
  theme_minimal(base_size = 14) +
  labs(title = "BEC DE Results Summary",
       x = NULL, y = "Number of Genes") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))
