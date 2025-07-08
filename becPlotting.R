library(pheatmap)
library(biomaRt)

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

# sort into up and down regulated
up_bec_genes <- sig_bec_genes[order(-sig_bec_genes$log2FoldChange), ][1:50, ]
down_bec_genes <- sig_bec_genes[order(sig_bec_genes$log2FoldChange), ][1:50, ]

# combine and get gene names
top_bec_de_genes <- unique(c(rownames(up_bec_genes), rownames(down_bec_genes)))

# subset the normalized matrix
heatmap_bec <- becNormalizedCounts[rownames(becNormalizedCounts) %in% top_bec_de_genes, ]

# set gene symbols as row names
rownames(heatmap_bec) <- heatmap_bec[[9]]   
heatmap_bec <- heatmap_bec[, -9]

# z-score normalize
heatmap_bec <- t(scale(t(heatmap_bec)))

pheatmap(heatmap_bec,
         show_rownames = TRUE,
         show_colnames = TRUE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_row = 6,
         fontsize_col = 8,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
