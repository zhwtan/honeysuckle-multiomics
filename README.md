# Honeysuckle transcriptome-metabolome analysis

This repository contains the minimum processed data and R code required to reproduce the transcriptome, metabolome, and gene-metabolite correlation analyses reported in the manuscript.

## Files

- `data/transcriptome_matrix_part1.tsv.gz` and `data/transcriptome_matrix_part2.tsv.gz`: two row-wise parts of one transcriptome matrix containing raw count values, FPKM values, upstream DESeq2 results, gene annotations, and KEGG annotations for 37,608 genes and 36 samples. The script combines them automatically.
- `data/metabolome_matrix.tsv.gz`: metabolite annotations and relative LC-MS peak-area values for 1,326 metabolites and 36 samples. The six non-final differential-metabolite flag columns from the source table were excluded.
- `Gemh_analysis.R`: portable analysis script using relative paths only.

Raw sequencing reads, raw mass-spectrometry files, vendor reports, manuscript files, figures, and temporary analysis objects are not included.

## Requirements

```r
install.packages(c("data.table", "mixOmics", "ggplot2"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("clusterProfiler")
```

## Run

Run from the repository root:

```bash
Rscript Gemh_analysis.R
```

The script performs the following analyses:

1. extraction of upstream DESeq2 results and identification of core differentially expressed genes across six stage-matched comparisons (`BH-FDR < 0.05` and `|log2FC| > 1`);
2. PCA of the top 1,000 variable genes based on FPKM values;
3. stage-wise PLS-DA and Welch's t-tests for metabolites, followed by identification of core differential metabolites (`VIP > 1` and `P < 0.05`);
4. KEGG over-representation analysis;
5. Pearson correlations with BH-FDR correction and partial Pearson correlations controlling for cultivar and developmental stage.

Successful execution reproduces 1,248 core differentially expressed genes, 259 core differential metabolites, and 323,232 gene-metabolite correlation tests. It generates Tables S13-S17 and the analysis-derived panels for Figures 3B-D, 4A-D, S3, S4, S5A-B/D, and S6, together with the PLS-DA and volcano plots. Photographs, genomic-variation schematics, and experimental-validation figures are outside the scope of these expression matrices. Generated files are written to `results/` and are not stored in the repository.

