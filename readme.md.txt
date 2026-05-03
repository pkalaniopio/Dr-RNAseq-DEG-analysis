# RNA-seq Differential Expression Analysis

This repository contains an RNA-seq differential expression workflow for zebrafish samples analyzed with DESeq2.

## Overview

The analysis compares each treatment condition against vehicle control (`Veh`) while accounting for batch effects.

The workflow includes:

- DESeq2 differential expression analysis
- Batch-adjusted design: `~ Batch + Condition`
- Log2 fold-change shrinkage using `apeglm`
- Zebrafish gene annotation using Entrez IDs, ZFIN IDs, gene symbols, and gene descriptions
- PCA plot for sample-level QC
- Volcano plots for each treatment vs vehicle
- Treatment-collapsed heatmaps of top DEGs
- Excel workbooks containing complete DEG tables and top DEG summaries

## Experimental Design

Samples are described in `data/metadata.csv`.

Required metadata columns:

| Column | Description |
|---|---|
| `Sample` | Sample ID matching count matrix column names |
| `Condition` | Treatment group, with `Veh` as reference |
| `Batch` | Batch identifier |

Raw counts are provided in `data/readcount.csv`.

Required count matrix format:

| geneID | Sample1 | Sample2 | Sample3 |
|---|---:|---:|---:|
| 100000000 | 120 | 98 | 300 |

Rows are genes and columns are samples.

## Analysis Pipeline

The main analysis is performed in:

```text
analysis.Rmd