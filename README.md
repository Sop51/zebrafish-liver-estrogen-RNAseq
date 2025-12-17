# zebrafish-liver-estrogen-RNAseq

This repository contains bulk RNA-seq analysis code and resources for studying estrogen (E2) responses in zebrafish liver cell types.

## Overview

Bulk RNA-seq was performed on zebrafish liver **hepatocytes** and **biliary epithelial cells (BECs)** under two conditions:

- **EtOH**
- **E2 treatment**

Each condition includes **four biological replicates per cell type**.

### Pre-processing

- **`pre-processing-linux-command-parameters`**  
  Linux command-line parameters used for:
  - FASTQ quality control
  - Read alignment
  - Generation of count matrices

### Differential Expression Analysis

- **`deseqHEP.R`**  
  Runs DESeq2 analysis on hepatocyte samples and generates initial result plots.

- **`deseqBEC.R`**  
  Runs DESeq2 analysis on biliary epithelial cell (BEC) samples and generates initial result plots.

### Plotting and Visualization

- **`hepPlotting.R`**  
  Generates plots using hepatocyte count data and differential expression results.

- **`becPlotting.R`**  
  Generates plots using BEC count data and differential expression results.

- **`deseqPlottingCombined.R`**  
  Creates combined plots comparing hepatocytes and BECs.

### Cross-species Comparison

- **`externalMicePregnancyData.R`**  
  Compares differential expression results from:
  - E2-treated zebrafish BECs  
  - BECs from pregnant mice
