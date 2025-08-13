library(pheatmap)
library(clusterProfiler)
library(org.Dr.eg.db)
library(ggplot2)
library(tidyr)
library(ggrepel)
library(reshape2)

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

# ------------------ simplifying GO and KEGG results ------------------- #

# simplify each GO result ----
go_hep_up_simplified <- simplify(go_results_hep_up, cutoff = 0.7, by = "p.adjust", select_fun = min)

# take the top 10 from each
up_go <- head(go_hep_up_simplified@result, 10)

# negative log pval
up_go <- up_go %>%
  mutate(
    signed_logp = -log10(p.adjust),  
  )

# plot
ggplot(up_go, aes(x = signed_logp, y = reorder(Description, signed_logp))) +
  geom_segment(aes(x = 0, xend = signed_logp, y = Description, yend = Description), color = "black") +
  geom_point(color = "#EB712A", size = 3) +  
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 5)
  ) +
  labs(
    x = "-log10 adjusted p-value",
    y = NULL,
    title = "Hepatocyte GO Term Enrichment: Upregulated"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9))

# black background plot
ggplot(up_go, aes(x = signed_logp, y = reorder(Description, signed_logp))) +
  geom_segment(aes(x = 0, xend = signed_logp, y = Description, yend = Description), color = "white") +
  geom_point(color = "#EB712A", size = 3) +  # orange/red for upregulated
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 5)
  ) +
  labs(
    x = "-log10 adjusted p-value",
    y = NULL,
    title = "Hepatocyte GO Term Enrichment: Upregulated"
  ) +
  theme_minimal(base_family = "Arial") +
  theme(
    plot.background = element_rect(fill = "black", color = NA),
    panel.background = element_rect(fill = "black", color = NA),
    panel.grid.major = element_line(color = "gray30"),
    panel.grid.minor = element_line(color = "gray20"),
    axis.text = element_text(color = "white", size = 9),
    axis.title = element_text(color = "white"),
    plot.title = element_text(color = "white", face = "bold", size = 14)
  )


# take the top 10 from kegg results ----
up_kegg <- head(kegg_hep_up@result, 10)

# negative log pval
up_kegg <- up_kegg %>%
  mutate(
    signed_logp = -log10(p.adjust),  
  )

# plot
ggplot(up_kegg, aes(x = signed_logp, y = reorder(Description, signed_logp))) +
  geom_segment(aes(x = 0, xend = signed_logp, y = Description, yend = Description), color = "black") +
  geom_point(color = "#EB712A", size = 3) +  
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 5)
  ) +
  labs(
    x = "-log10 adjusted p-value",
    y = NULL,
    title = "Hepatocyte KEGG Pathway Enrichment: Upregulated"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9))

# black background plot
ggplot(up_kegg, aes(x = signed_logp, y = reorder(Description, signed_logp))) +
  geom_segment(aes(x = 0, xend = signed_logp, y = Description, yend = Description), color = "white") +
  geom_point(color = "#EB712A", size = 3) +  # orange/red for upregulated
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 5)
  ) +
  labs(
    x = "-log10 adjusted p-value",
    y = NULL,
    title = "Hepatocyte KEGG Pathway Enrichment: Upregulated"
  ) +
  theme_minimal(base_family = "Arial") +
  theme(
    plot.background = element_rect(fill = "black", color = NA),
    panel.background = element_rect(fill = "black", color = NA),
    panel.grid.major = element_line(color = "gray30"),
    panel.grid.minor = element_line(color = "gray20"),
    axis.text = element_text(color = "white", size = 9),
    axis.title = element_text(color = "white"),
    plot.title = element_text(color = "white", face = "bold", size = 14)
  )

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

# ----------------- plotting hepatocyte genes in EtOH ------------ #

# pull out only the normal samples
etoh_hep_samples <- hepNormalizedCounts[, grepl("EtOH", colnames(hepNormalizedCounts))]

# name genes
genes_of_interest <- c("ENSDARG00000038439", "ENSDARG00000015866", 
                       "ENSDARG00000037281", "ENSDARG00000008969", "ENSDARG00000090850", 
                       "ENSDARG00000090286")

gene_names <- c("fabp10a", "apoa2", "fgg", "fgb", "serpina1l", "serpina1")

# map ensembl to name
gene_mapping <- data.frame(
  ensembl_id = genes_of_interest,
  gene_name = gene_names
)

# extract the genes and add ensembl id col
gene_expression_data <- etoh_hep_samples[genes_of_interest, ]
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
    title = "Expression of Hepatocyte Markers in EtOH Samples",
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
    title = "Expression of Hepatocyte Markers in EtOH Samples",
    x = "Gene",
    y = "VST Normalized Counts"
  ) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))

# ----------------- plotting hepatocyte protein processing in endoplasmic reticulum KEGG genes  ------------ #

# name the ER genes
er_genes_to_plot <- c('hyou1', 'ero1a', 'calr3b', 'pdia6', 'pdia4', 'hsp90b1', 'hspa5', 'calr')

# ensembl names
ensembl_er_genes_to_plot <- c('ENSDARG00000013670', 'ENSDARG00000015228', 
                              'ENSDARG00000102808', 'ENSDARG00000009001', 
                              'ENSDARG00000018491', 'ENSDARG00000003570', 
                              'ENSDARG00000103846', 'ENSDARG00000076290')

# subset the count df
er_subset <- hepNormalizedCounts[rownames(hepNormalizedCounts) %in% ensembl_er_genes_to_plot, ]

# set Gene.name as rownames
rownames(er_subset) <- er_subset$Gene.name

# remove the Gene.name column before plotting
er_subset$Gene.name <- NULL

scaled_er_subset <- t(scale(t(er_subset)))

# plot heatmap
pheatmap(scaled_er_subset,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = colorRampPalette(c("#523095", "white", "#EB712A"))(100),
         main = "KEGG Protein Processing in Endoplasmic Reticulum Genes")

# ---------------------- create a plot of normalized counts in a scatter plot ----------------- #
# only subset to controls
hep_norm_counts_etoh <- hepNormalizedCounts[, grep("EtOH", colnames(hepNormalizedCounts))]

# calculate avg across samples
gene_mean_exp_hep <- data.frame(
  gene = rownames(hep_norm_counts_etoh),
  x = rowMeans(hep_norm_counts_etoh, na.rm = TRUE),
  y = rowMeans(hep_norm_counts_etoh, na.rm = TRUE)
)

# merge with gene symbol information, biomart
mart <- useEnsembl(biomart = "genes", dataset = "drerio_gene_ensembl")
mapping <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = gene_mean_exp_hep$gene,
  mart = mart
)

# merge
gene_mean_exp_hep <- gene_mean_exp_hep %>%
  left_join(mapping, by = c("gene" = "ensembl_gene_id")) %>% 
  mutate(external_gene_name = tolower(external_gene_name))

# edit er gene names
gene_mean_exp_hep[gene_mean_exp_hep$gene == "ENSDARG00000034181", "external_gene_name"] <- "esr2b"
gene_mean_exp_hep[gene_mean_exp_hep$gene == "ENSDARG00000016454", "external_gene_name"] <- "esr2a"

# create labels for genes
genes_to_label <- c("apoa2", "serpina1l", "krt8", "anxa4","esr1", "esr2a", "esr2b", "gper1")
gene_mean_exp_hep$label <- ifelse(
  gene_mean_exp_hep$external_gene_name %in% genes_to_label,
  gene_mean_exp_hep$external_gene_name,
  NA
)

# column to indicate if labeled or not
gene_mean_exp_hep <- gene_mean_exp_hep %>%
  mutate(is_labeled = !is.na(label))

# create a color column based on gene name
gene_mean_exp_hep$gene_color <- case_when(
  gene_mean_exp_hep$external_gene_name %in% c("anxa4", "krt8") ~ "#43CCB8",
  gene_mean_exp_hep$external_gene_name %in% c("serpina1l", "apoa2") ~ "magenta",
  gene_mean_exp_hep$external_gene_name %in% c("esr2b", "esr2a", "esr1", "gper1") ~ "black",
  TRUE ~ NA_character_
)

# plot
ggplot(gene_mean_exp_hep, aes(x = x, y = y)) +
  geom_point(
    data = subset(gene_mean_exp_hep, !is_labeled),
    shape = 21,
    color = "#DADADA",
    fill = NA,
    size = 2
  ) +
  geom_point(
    data = subset(gene_mean_exp_hep, is_labeled),
    aes(fill = gene_color),
    shape = 21,
    color = "black",  # outline
    size = 3,
    show.legend = FALSE
  ) +
  geom_text_repel(
    data = subset(gene_mean_exp_hep, is_labeled),
    aes(label = label),
    size = 3.5,
    max.overlaps = 20,
    arrow = arrow(length = unit(0.02, "npc")),
    box.padding = 0.4,
    point.padding = 0.5,
    segment.color = 'black'
  ) +
  scale_fill_identity() +
  theme_minimal() +
  labs(
    title = "Gene Expression Across EtOH Samples in Hepatocytes",
    x = "Avg VST Normalized Expression",
    y = "Avg VST Normalized Expression"
  )
