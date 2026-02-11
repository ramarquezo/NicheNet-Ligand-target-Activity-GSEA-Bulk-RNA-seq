# NicheNet Ligand Activity + Target-GSEA (Bulk RNA-seq, Two-Dataset Workflow)

Run a **bulk NicheNet** ligand activity analysis between a **Sender** and **Receiver** dataset, visualize the **top ligands with sender DE evidence**, and then perform **GSEA on predicted ligand target sets** in the receiver to identify ligands whose downstream targets are coordinately up/downregulated.

This workflow is designed for **bioinformatics** use cases where you have two differential expression (DE) tables representing communicating cell populations (e.g., astrocytes → tumor, T cells → macrophages, stroma → cancer).

---

## What this script does

Given two DE tables:

- **Sender** (candidate ligands; e.g., secreted factors expressed/induced in sender)
- **Receiver** (responding program; DE signature used for NicheNet)

The script:

1. Loads sender + receiver DE tables (Excel).
2. Builds a receiver DE **gene signature** and **background**.
3. Runs **NicheNet bulk ligand activity inference** (sender → receiver).
4. Plots:
   - ligand activity scatter (rank vs AUPR corrected)
   - top ligands barplot (AUPR corrected)
   - sender log2FC evidence for the same top ligands
   - optional composite panel (if `{patchwork}` is installed)
5. Builds ligand-specific **target gene sets** (top *N* targets per ligand by regulatory potential).
6. Runs **GSEA** on the **receiver ranked list** using those ligand target gene sets.
7. Exports tables (TSV), plots (PNG/PDF), and an optional `.rds` with all objects for reproducibility.

---

## Requirements

R ≥ 4.1 recommended

### Install packages

```r
install.packages(c(
  "dplyr", "tidyr", "tibble", "readxl", "readr", "ggplot2"
))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")

# NicheNet
install.packages("remotes")
remotes::install_github("saeyslab/nichenetr")

# Optional (for composite plots)
install.packages("patchwork")
