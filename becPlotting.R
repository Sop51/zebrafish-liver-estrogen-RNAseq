library(pheatmap)
library(clusterProfiler)
library(org.Dr.eg.db)
library(ggplot2)
library(dplyr)
library(genekitr)

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
rownames(heatmap_bec) <- heatmap_bec[[7]]   
heatmap_bec <- heatmap_bec[, -7]

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
                       ont = "BP", 
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
go_results_bec_down <- enrichGO(gene = entrez_ids_bec_down,
                              OrgDb = org.Dr.eg.db,
                              keyType = "ENTREZID",
                              ont = "MF", 
                              pAdjustMethod = "BH",
                              qvalueCutoff = 0.05,
                              readable = TRUE)
# plot
dotplot(go_results_bec_down, showCategory = 15) + 
  ggtitle("E2 Downregulated GO Biological Process") +
  theme(axis.text.y = element_text(size = 8, face = "plain"))

# KEGG enrichment
kegg_bec_down <- enrichKEGG(gene = entrez_ids_bec_down,
                          organism = "dre",
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05)

# plot
dotplot(kegg_bec_down, showCategory = 15) + ggtitle("E2 Downregulated KEGG")

# ------------------ simplifying GO and KEGG results ------------------- #

# simplify each GO result ----
go_bec_up_simplified <- simplify(go_results_bec_up, cutoff = 0.7, by = "p.adjust", select_fun = min)
go_bec_down_simplified <- simplify(go_results_bec_down, cutoff = 0.7, by = "p.adjust", select_fun = min)

# take the top 10 from each
up_go <- head(go_bec_up_simplified@result, 10)
down_go <- head(go_bec_down_simplified@result, 10)

# try plotting
# add regulation labels
up_go$Regulation <- "Up"
down_go$Regulation <- "Down"

# combine
go_combined <- bind_rows(up_go, down_go) %>%
  mutate(logp = -log10(p.adjust),
         signed_logp = ifelse(regulation == "Down", -logp, logp)) %>%
  arrange(signed_logp) %>%
  mutate(term = factor(Description, levels = unique(Description)))

# plot
reg_colors <- c("Down" = "#523095",  
                "Up"   = "#EB712A")

ggplot(go_combined, aes(x = signed_logp, y = term)) +
  geom_segment(aes(x = 0, xend = signed_logp, y = term, yend = term), color = "black") +
  geom_point(aes(color = Regulation), size = 3) +
  scale_color_manual(values = reg_colors) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 10),
    labels = function(x) abs(x)  
  ) +
  labs(x = "-log10 adjusted p-value", y = NULL,
       title = "BEC GO Term Enrichment: Up and Down Regulation") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

# black background plot
ggplot(go_combined, aes(x = signed_logp, y = term)) +
  geom_segment(aes(x = 0, xend = signed_logp, y = term, yend = term), color = "white") +
  geom_point(aes(color = Regulation), size = 3) +
  scale_color_manual(values = reg_colors) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 10),
    labels = function(x) abs(x)
  ) +
  labs(
    x = "-log10 adjusted p-value", 
    y = NULL,
    title = "BEC GO Term Enrichment: Up and Down Regulation"
  ) +
  theme_minimal(base_family = "Arial") +
  theme(
    plot.background = element_rect(fill = "black", color = NA),
    panel.background = element_rect(fill = "black", color = NA),
    panel.grid.major = element_line(color = "gray30"),
    panel.grid.minor = element_line(color = "gray20"),
    axis.text = element_text(color = "white", size = 8),
    axis.title = element_text(color = "white"),
    plot.title = element_text(color = "white", face = "bold", size = 14),
    legend.background = element_rect(fill = "black", color = NA),
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white")
  )

# take the top 10 from each kegg results ----
up_kegg <- head(kegg_bec_up@result, 10)
down_kegg <- head(kegg_bec_down@result, 10)

# try plotting
# add regulation labels
up_kegg$Regulation <- "Up"
down_kegg$Regulation <- "Down"

# combine
kegg_combined <- bind_rows(up_kegg, down_kegg) %>%
  mutate(logp = -log10(p.adjust),
         signed_logp = ifelse(Regulation == "Down", -logp, logp)) %>%
  arrange(signed_logp) %>%
  mutate(term = factor(Description, levels = unique(Description)))

# plot
reg_colors <- c("Down" = "#523095",  
                "Up"   = "#EB712A")

ggplot(kegg_combined, aes(x = signed_logp, y = term)) +
  geom_segment(aes(x = 0, xend = signed_logp, y = term, yend = term), color = "black") +
  geom_point(aes(color = Regulation), size = 3) +
  scale_color_manual(values = reg_colors) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 10),
    labels = function(x) abs(x)  
  ) +
  labs(x = "-log10 adjusted p-value", y = NULL,
       title = "BEC KEGG Pathway Enrichment: Up and Down Regulation") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

# black background plot
ggplot(kegg_combined, aes(x = signed_logp, y = term)) +
  geom_segment(aes(x = 0, xend = signed_logp, y = term, yend = term), color = "white") +
  geom_point(aes(color = Regulation), size = 3) +
  scale_color_manual(values = reg_colors) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 10),
    labels = function(x) abs(x)
  ) +
  labs(
    x = "-log10 adjusted p-value", 
    y = NULL,
    title = "BEC KEGG Pathway Enrichment: Up and Down Regulation"
  ) +
  theme_minimal(base_family = "Arial") +
  theme(
    plot.background = element_rect(fill = "black", color = NA),
    panel.background = element_rect(fill = "black", color = NA),
    panel.grid.major = element_line(color = "gray30"),
    panel.grid.minor = element_line(color = "gray20"),
    axis.text = element_text(color = "white", size = 8),
    axis.title = element_text(color = "white"),
    plot.title = element_text(color = "white", face = "bold", size = 14),
    legend.background = element_rect(fill = "black", color = NA),
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white")
  )

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

# ----------------- plotting BEC genes in EtOH ------------ #

# pull out only the normal samples
etoh_bec_samples <- becNormalizedCounts[, grepl("EtOH", colnames(becNormalizedCounts))]

# name genes
genes_of_interest <- c("ENSDARG00000036456", "ENSDARG00000038153", 
                       "ENSDARG00000040747", "ENSDARG00000100844", "ENSDARG00000005565", 
                       "ENSDARG00000075163")

gene_names <- c("anxa4", "lgals2b", "tm4sf4", "cldn15lb", "entpd8", "cxcl20")

# map ensembl to name
gene_mapping <- data.frame(
  ensembl_id = genes_of_interest,
  gene_name = gene_names
)

# extract the genes and add ensembl id col
gene_expression_data <- etoh_bec_samples[genes_of_interest, ]
gene_expression_data$ensembl_id <- rownames(gene_expression_data)

# long format 
long_data <- melt(gene_expression_data, id.vars = "ensembl_id", variable.name = "sample", value.name = "expression")

# merge with gene names
long_data <- merge(long_data, gene_mapping, by = "ensembl_id")

# plot
ggplot(long_data, aes(x = gene_name, y = expression)) +
  geom_boxplot(fill = "#EE9F4A", color = "black", outlier.size = 1) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5, color = "darkgray") +
  
  theme_classic() +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black", face = "italic"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 14, color = "black", margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, color = "black", margin = margin(r = 10)),
    plot.title = element_text(size = 16, color = "black", hjust = 0.5, face = "bold"),
    
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    
    panel.grid = element_blank(),
    plot.margin = margin(20, 20, 20, 20),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  
  labs(
    title = "Expression of BEC Markers in EtOH Samples",
    x = "Gene",
    y = "VST Normalized Counts"
  ) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) # put y at 0

# plot with black background
ggplot(long_data, aes(x = gene_name, y = expression)) +
  geom_boxplot(fill = "#EE9F4A", color = "white", outlier.size = 1) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5, color = "gray80") +
  
  theme_classic() +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "white", face = "italic"),
    axis.text.y = element_text(size = 12, color = "white"),
    axis.title.x = element_text(size = 14, color = "white", margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, color = "white", margin = margin(r = 10)),
    plot.title = element_text(size = 16, color = "white", hjust = 0.5, face = "bold"),
    
    axis.line = element_line(color = "white", size = 0.5),
    axis.ticks = element_line(color = "white", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    
    panel.grid = element_blank(),
    plot.margin = margin(20, 20, 20, 20),
    
    panel.background = element_rect(fill = "black", color = NA),
    plot.background = element_rect(fill = "black", color = NA)
  ) +
  
  labs(
    title = "Expression of BEC Markers in EtOH Samples",
    x = "Gene",
    y = "VST Normalized Counts"
  ) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))

