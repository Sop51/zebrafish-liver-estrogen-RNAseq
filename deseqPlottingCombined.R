library(tidyverse)
library(DESeq2)
library(pheatmap)
library(grid)
library(clusterProfiler)
library(org.Dr.eg.db)
library(VennDiagram)

# read in the normalized counts for each dataset
hep_norm_counts <- read.csv('/Users/sophiemarcotte/Desktop/patrice/estrogenRNAseq/vstCountsHEP.csv', row.names = 1)
bec_norm_counts <- read.csv('/Users/sophiemarcotte/Desktop/patrice/estrogenRNAseq/vstCountsBEC.csv', row.names = 1)

# read in the raw counts for each dataset
hep_raw_counts <- read.csv('/Users/sophiemarcotte/Desktop/patrice/estrogenRNAseq/rawCountsHEP.csv', row.names = 1)
bec_raw_counts <- read.csv('/Users/sophiemarcotte/Desktop/patrice/estrogenRNAseq/rawCountsBEC.csv', row.names = 1)

# read in the DE results
hepDeseqResults <- read.csv('/Users/sm2949/Desktop/patrice/estrogenRNAseq/deseqResultsHEP.csv', row.names=1)
becDeseqResults <- read.csv('/Users/sm2949/Desktop/patrice/estrogenRNAseq/deseqResultsBEC.csv', row.names=1)

# filter for sig results
hepDeseqResults_sig <- hepDeseqResults %>%
  filter(abs(log2FoldChange) > 1, padj < 0.05)
becDeseqResults_sig <- becDeseqResults %>%
  filter(abs(log2FoldChange) > 1, padj < 0.05)

# -------------------- prepare the count matrix for deseq ONLY ETOH ------------------------ #
# pull out only control samples
etoh_bec_raw <- bec_raw_counts[, grepl("EtOH", colnames(bec_raw_counts))]
etoh_hep_raw <- hep_raw_counts[, grepl("EtOH", colnames(hep_raw_counts))]

# merge
etoh_combined_raw <- merge(etoh_hep_raw, etoh_bec_raw, 
                           by = "row.names", all = TRUE)

# set gene ID as rownames again
rownames(etoh_combined_raw) <- etoh_combined_raw$Row.names
etoh_combined_raw$Row.names <- NULL

# ------------------------- create the col data for deseq ONLY ETOH ----------------------------
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

# ----------------------- running deseq ETOH ONLY ----------------------------------------
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

# ------------------------------- create plots for marker genes in etoh ------------------------------------

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

# combine samples
combined_counts <- cbind(etoh_hep_samples, etoh_bec_samples)

# subset to marker genes
hep_bec_genes <- union(hep_ensembl_markers, bec_ensembl_markers)
subset_counts <- combined_counts[hep_bec_genes, ]

# custom color assignment for annotations
annotation_colors <- list(
  Cell_Type = c(
    "Hepatocyte" = "#FF0090",  
    "BEC" = "#97E997"   
  )
)

# create annotation for sample cell type
sample_annotations <- data.frame(
  Cell_Type = ifelse(grepl("EtOH.*HEP", colnames(subset_counts)), "Hepatocyte", "BEC")
)
rownames(sample_annotations) <- colnames(subset_counts)

# add gene symbols for row names over ensembl
gene_lookup <- data.frame(
  ensembl = c(hep_ensembl_markers, bec_ensembl_markers),
  gene = c(hep_gene_markers, bec_gene_markers),
  Marker_For = c(rep("Hepatocyte", length(hep_gene_markers)),
                  rep("BEC", length(bec_gene_markers)))
)

# rename rownames to gene symbols
rownames(subset_counts) <- gene_lookup$gene[match(rownames(subset_counts), gene_lookup$ensembl)]
rownames(gene_annotations) <- rownames(subset_counts)

# plot heatmap
p <- pheatmap(subset_counts,
         annotation_col = sample_annotations,
         annotation_colors = annotation_colors,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_colnames = TRUE,
         display_numbers = TRUE,
         show_rownames = TRUE,
         fontsize_row = 10,
         fontsize_col = 8,
         scale = "row", # scale 
         color = colorRampPalette(c("#963489", "white", "#E6C67B"))(100),
         main = "Marker Gene Expression (EtOH Samples)")

# give a black background
grid.draw(rectGrob(gp=gpar(fill="black", lwd=0)))
grid.draw(p)
grid.gedit("layout", gp = gpar(col = "white", text = ""))

# -------------------- prepare the count matrix for deseq COMBINED ------------------------ #
# merge
combined_raw <- merge(hep_raw_counts, bec_raw_counts, 
                           by = "row.names", all = TRUE)

# set gene ID as rownames again
rownames(combined_raw) <- combined_raw$Row.names
combined_raw$Row.names <- NULL

# ------------------------- create the col data for deseq COMBINED ----------------------------
colnames_combined <- colnames(combined_raw)

# split by . and extract condition and replicate
coldata_combined <- data.frame(
  cell_type = sapply(strsplit(colnames_combined, "\\."), `[`, 2),
  replicate = sapply(strsplit(colnames_combined, "\\."), `[`, 3)
)
rownames(coldata_combined) <- colnames_combined

# reorder to match counts col names
coldata_combined <- coldata_combined[colnames(combined_raw), ]
all(rownames(coldata_combined) == colnames(combined_raw))  # assert

# set cell type and replicate as a factor
coldata_combined$cell_type <- as.factor(coldata_combined$cell_type)
coldata_combined$replicate <- as.factor(coldata_combined$replicate)

# ----------------------- running deseq ETOH ONLY ----------------------------------------
# create the deseq object
dds_comb <- DESeqDataSetFromMatrix(countData = combined_raw,
                                   colData = coldata_combined,
                                   design = ~ 1)

# run DE
dds_comb <- DESeq(dds_comb)

# vst
vst <- assay(vst(dds_comb))
write.csv(vst, file = "/Users/sophiemarcotte/Desktop/patrice/estrogenRNAseq/vstCountsCOMBINED.csv", row.names = TRUE)


# ------------------- plotting estrogen receptor genes ------------------------- #

# read in the combined normalized counts
comb_normalized_counts <- read.csv("/Users/sophiemarcotte/Desktop/patrice/estrogenRNAseq/vstCountsCOMBINED.csv", row.names = 1)

# name the receptor
r_genes_to_plot <- c('gper1', 'esr1', 'esr2a', 'esr2b')

# ensembl names
ensembl_r_genes_to_plot <- c('ENSDARG00000074661','ENSDARG00000004111',
                             'ENSDARG00000016454','ENSDARG00000034181')

# subset the count df
r_subset <- comb_normalized_counts[rownames(comb_normalized_counts) %in% ensembl_r_genes_to_plot, ]

# match order of gene symbols to rows in r_subset
matched_gene_symbols <- r_genes_to_plot[match(rownames(r_subset), ensembl_r_genes_to_plot)]

# set gene symbols as rownames
rownames(r_subset) <- matched_gene_symbols

# custom color assignment for annotations
annotation_colors <- list(
  Cell_Type = c(
    "Hepatocyte" = "#FF0090",  
    "BEC" = "#97E997"   
  )
)

# create annotation for sample cell type
sample_annotations <- data.frame(
  Cell_Type = ifelse(grepl("HEP", colnames(r_subset)), "Hepatocyte", "BEC")
)
rownames(sample_annotations) <- colnames(r_subset)

scaled_r_subset <- t(scale(t(r_subset)))

# plot heatmap
p <- pheatmap(scaled_r_subset,
         annotation_col = sample_annotations,
         annotation_colors = annotation_colors,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_colnames = TRUE,
         display_numbers = TRUE,
         show_rownames = TRUE,
         fontsize_row = 10,
         fontsize_col = 8,
         scale = "row", # scale 
         color = colorRampPalette(c("#963489", "white", "#E6C67B"))(100),
         main = "Estrogen Receptor Genes")

# give a black background
grid.draw(rectGrob(gp=gpar(fill="black", lwd=0)))
grid.draw(p)
grid.gedit("layout", gp = gpar(col = "white", text = ""))

# -------- find overlapping and non - overlapping DE genes between cell types -----------
# get genes that are common
common_de_genes <- intersect(rownames(hepDeseqResults_sig), rownames(becDeseqResults_sig))

# pull out the rows of the common genes
hepOverlap <- hepDeseqResults_sig[common_de_genes,]
becOverlap <- becDeseqResults_sig[common_de_genes,]

# get de genes unique to hep
hep_unique_genes <- setdiff(rownames(hepDeseqResults_sig), common_de_genes)
hepUnique <- hepDeseqResults_sig[hep_unique_genes, ]

# get de genes unique to bec
bec_unique_genes <- setdiff(rownames(becDeseqResults_sig), common_de_genes)
becUnique <- becDeseqResults_sig[bec_unique_genes, ]

# create a venn diagram ----
# clear the page
grid.newpage()

# plot
venn.plot <- draw.pairwise.venn(
  area1 = length(rownames(hepDeseqResults_sig)),
  area2 = length(rownames(becDeseqResults_sig)),
  cross.area = length(common_de_genes),
  category = c("Hep", "BEC"),
  fill = c("#FF0090", "#97E997"),
  scaled = TRUE,
  ext.text = FALSE,
  cex = 2,
  cat.cex = 2,
  cat.pos = c(0, 180), 
  cat.dist = c(0.05, 0.025)
)

# common genes into kegg / go ----
# convert to entrez id
entrez_ids_common <- mapIds(org.Dr.eg.db,
                            keys = common_de_genes,
                            column = "ENTREZID",
                            keytype = "ENSEMBL",
                            multiVals = "first")

# run go
go_results_common <- enrichGO(gene = entrez_ids_common,
                              OrgDb = org.Dr.eg.db,
                              keyType = "ENTREZID",
                              ont = "BP", 
                              pAdjustMethod = "BH",
                              qvalueCutoff = 0.05,
                              readable = TRUE)
# plot
dotplot(go_results_common, showCategory = 15) + ggtitle("Common DE Genes GO Biological Process")

# KEGG enrichment
kegg_bec_common <- enrichKEGG(gene = entrez_ids_common,
                          organism = "dre",
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05)

# plot
dotplot(kegg_bec_common, showCategory = 15) + ggtitle("Common DE Genes KEGG")

# plotting some of the DE common genes -----

er_stress_ensembl <- c('ENSDARG00000102808', 'ENSDARG00000009001', 'ENSDARG00000018491',
                       'ENSDARG00000003570', 'ENSDARG00000103846', 'ENSDARG00000076290')

er_stress_gene_name <- c('calr3b', 'pdia6', 'pdia4', 'hsp90b1', 'hspa5', 'calr')
  
dna_rep_ensembl <- c('ENSDARG00000002304', 'ENSDARG00000024204', 'ENSDARG00000011404',
                     'ENSDARG00000040041', 'ENSDARG00000019507', 'ENSDARG00000057683')
  
dna_rep_gene_name <- c('gins2', 'mcm3', 'fen1', 'mcm4', 'mcm5', 'mcm6')

# create a df for annotations and gene name mapping
er_dna_genes <- data.frame(
  ensembl = c(er_stress_ensembl, dna_rep_ensembl),
  gene_name = c(er_stress_gene_name, dna_rep_gene_name),
  category = c(rep("Response to ER Stress", length(er_stress_ensembl)),
               rep("DNA Replication", length(dna_rep_ensembl)))
)

# extract rows from deseq results
hep_lfc <- hepDeseqResults[er_dna_genes$ensembl, "log2FoldChange", drop = FALSE]
bec_lfc <- becDeseqResults[er_dna_genes$ensembl, "log2FoldChange", drop = FALSE]

# combine and put gene symbols as row names
log2fc_mat <- cbind(hep_lfc, bec_lfc)
colnames(log2fc_mat) <- c("EtOH Hep vs E2 Hep", "EtOH BEC vs E2 BEC")
rownames(log2fc_mat) <- er_dna_genes$gene_name

# create row annotation df 
row_anno <- data.frame(GO_BP = er_dna_genes$category)
rownames(row_anno) <- er_dna_genes$gene_name

# create col annotation
col_anno <- data.frame(Cell_Type = c("Hepatocyte", "BEC"))
rownames(col_anno) <- colnames(log2fc_mat)

# set annotation colors
anno_colors <- list(
  GO_BP = c("Response to ER Stress" = "#66c2a5", "DNA Replication" = "#fc8d62"),
  Cell_Type = c("Hepatocyte" = "#FF0090", "BEC" = "#97E997")
)

# plot
p <- pheatmap(log2fc_mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_row = row_anno,
         annotation_col = col_anno,
         annotation_colors = anno_colors,
         color = colorRampPalette(c(
           "#ffba08ff",
           "#faa307ff",
           "#f48c06ff",
           "#e85d04ff",
           "#dc2f02ff",
           "#d00000ff",
           "#9d0208ff",
           "#6a040fff",
           "#370617ff",
           "#03071eff"
         ))(100),
         main = "LFC of ER Stress & DNA Replication Genes")

# give a black background
grid.draw(rectGrob(gp=gpar(fill="black", lwd=0)))
grid.draw(p)
grid.gedit("layout", gp = gpar(col = "white", text = ""))

# ------------------- plotting bile pump genes ------------------------- #

# read in the combined normalized counts
comb_normalized_counts <- read.csv("/Users/sm2949/Desktop/patrice/estrogenRNAseq/vstCountsCOMBINED.csv", row.names = 1)

# name the bile bump genes 
bp_genes_to_plot <- c('abcb11a', 'abcb11b', 'abcb4', 'abcc3', 'abcc4', 'slc10a1')

# ensembl names
ensembl_bp_genes_to_plot <- c('ENSDARG00000011573','ENSDARG00000070078',
                             'ENSDARG00000010936','ENSDARG00000096662',
                             'ENSDARG00000058953', 'ENSDARG00000030588')

# subset the count df
bp_subset <- comb_normalized_counts[rownames(comb_normalized_counts) %in% ensembl_bp_genes_to_plot, ]

# match order of gene symbols to rows in bp_subset
matched_gene_symbols <- bp_genes_to_plot[match(rownames(bp_subset), ensembl_bp_genes_to_plot)]

# set gene symbols as rownames
rownames(bp_subset) <- matched_gene_symbols

# custom color assignment for annotations
annotation_colors <- list(
  Cell_Type = c(
    "Hepatocyte" = "#FF0090",  
    "BEC" = "#97E997"   
  )
)

# create annotation for sample cell type
sample_annotations <- data.frame(
  Cell_Type = ifelse(grepl("HEP", colnames(bp_subset)), "Hepatocyte", "BEC")
)
rownames(sample_annotations) <- colnames(bp_subset)

scaled_bp_subset <- t(scale(t(bp_subset)))

# plot heatmap
p <- pheatmap(scaled_bp_subset,
              annotation_col = sample_annotations,
              annotation_colors = annotation_colors,
              cluster_rows = TRUE,
              cluster_cols = TRUE,
              show_colnames = TRUE,
              display_numbers = TRUE,
              show_rownames = TRUE,
              fontsize_row = 10,
              fontsize_col = 8,
              scale = "row", # scale 
              color = colorRampPalette(c("#963489", "white", "#E6C67B"))(100),
              main = "Bile Pump Genes")

# give a black background
grid.draw(rectGrob(gp=gpar(fill="black", lwd=0)))
grid.draw(p)
grid.gedit("layout", gp = gpar(col = "white", text = ""))
