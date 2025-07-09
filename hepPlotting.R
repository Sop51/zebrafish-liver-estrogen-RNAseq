library(pheatmap)
library(clusterProfiler)
library(org.Dr.eg.db)

# ------------------ set up data for analyses --------------------- #
# read in the normalized counts & deseq results
hepNormalizedCounts <- read.csv('/Users/sm2949/Desktop/patrice/estrogenRNAseq/vstCountsHEP.csv', row.names=1)
hepDeseqResults <- read.csv('/Users/sm2949/Desktop/patrice/estrogenRNAseq/deseqResultsHEP.csv', row.names=1)

# if no gene symbol, replace with ensembl
hepDeseqResults$Gene.name <- ifelse(
  is.na(hepDeseqResults$Gene.name) | hepDeseqResults$Gene.name == "",
  rownames(hepDeseqResults),
  hepDeseqResults$Gene.name
)

# add gene symbols to normalized counts
hepNormalizedCounts$Gene.name <- hepDeseqResults$Gene.name[match(rownames(hepNormalizedCounts), rownames(hepDeseqResults))]

# pull out only significant genes
sig_hep_genes <- hepDeseqResults[!is.na(hepDeseqResults$padj) & hepDeseqResults$padj < 0.05, ]

# ---------------------- plotting a heatmap of up & downreg genes ------------------- #

# sort into up and down regulated
up_hep_genes <- sig_hep_genes[order(-sig_hep_genes$log2FoldChange), ][1:25, ]
down_hep_genes <- sig_hep_genes[order(sig_hep_genes$log2FoldChange), ][1:25, ]

# combine and get gene names
top_hep_de_genes <- unique(c(rownames(up_hep_genes), rownames(down_hep_genes)))

# subset the normalized matrix
heatmap_hep <- hepNormalizedCounts[rownames(hepNormalizedCounts) %in% top_hep_de_genes, ]

# set gene symbols as row names
rownames(heatmap_hep) <- heatmap_hep[[7]]   
heatmap_hep <- heatmap_hep[, -7]

# z-score normalize
heatmap_hep <- t(scale(t(heatmap_hep)))

# plot
pheatmap(heatmap_hep,
         show_rownames = TRUE,
         show_colnames = TRUE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_row = 8,
         fontsize_col = 10,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

# ----------------------------- kegg and go ----------------------------- #
# pull out upregulated gene ids
up_hep_genes <- sig_hep_genes[sig_hep_genes$log2FoldChange > 1,]
up_gene_hep_names <- rownames(up_hep_genes)

# convert to entrez IDs
entrez_ids_hep_up <- mapIds(org.Dr.eg.db,
                            keys = up_gene_hep_names,
                            column = "ENTREZID",
                            keytype = "ENSEMBL",
                            multiVals = "first")

# run go
go_results_hep_up <- enrichGO(gene = entrez_ids_hep_up,
                              OrgDb = org.Dr.eg.db,
                              keyType = "ENTREZID",
                              ont = "MF", 
                              pAdjustMethod = "BH",
                              qvalueCutoff = 0.05,
                              readable = TRUE)
# plot
dotplot(go_results_hep_up, showCategory = 15) + ggtitle("E2 Upregulated GO Biological Process")

# KEGG enrichment
kegg_hep_up <- enrichKEGG(gene = entrez_ids_hep_up,
                          organism = "dre",
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05)

# plot
dotplot(kegg_hep_up, showCategory = 15) + ggtitle("E2 Upregulated KEGG")


# pull out downregulated gene ids
down_hep_genes <- sig_hep_genes[sig_hep_genes$log2FoldChange < 1,]
down_gene_hep_names <- rownames(down_hep_genes)

# convert to entrez IDs
entrez_ids_hep_down <- mapIds(org.Dr.eg.db,
                              keys = down_gene_hep_names,
                              column = "ENTREZID",
                              keytype = "ENSEMBL",
                              multiVals = "first")

# run go
down_results_hep_up <- enrichGO(gene = entrez_ids_hep_down,
                                OrgDb = org.Dr.eg.db,
                                keyType = "ENTREZID",
                                ont = "MF", 
                                pAdjustMethod = "BH",
                                qvalueCutoff = 0.05,
                                readable = TRUE)
# plot
dotplot(down_results_hep_up, showCategory = 15) + 
  ggtitle("E2 Downregulated GO Biological Process") 

# KEGG enrichment
kegg_hep_down <- enrichKEGG(gene = entrez_ids_hep_down,
                            organism = "dre",
                            pAdjustMethod = "BH",
                            qvalueCutoff = 0.05)

# plot
dotplot(kegg_hep_down, showCategory = 15) + ggtitle("E2 Downregulated KEGG")

# ------------------ bar plot of up, down, ns gene count ---------------- #

# categorize genes
hepDeseqResults$regulation <- ifelse(
  hepDeseqResults$padj < 0.05 & hepDeseqResults$log2FoldChange > 1, "Upregulated",
  ifelse(hepDeseqResults$padj < 0.05 & hepDeseqResults$log2FoldChange < 1, "Downregulated", "Not Significant")
)

# count how many genes in each category
reg_counts <- table(hepDeseqResults$regulation)
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
  labs(title = "Hepatocyte DE Results Summary",
       x = NULL, y = "Number of Genes") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))
