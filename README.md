# IsoDecipher 🧬

**IsoDecipher** is a targeted isoform quantification tool designed for 3' single-cell RNA-seq (scRNA-seq) data. It recovers hidden biological signals—such as Alternative Polyadenylation (APA) and B-cell antibody switching—by re-interpreting 3' read patterns that standard pipelines typically collapse into gene-level counts.

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)  
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue)]()  
[![Compatible with Cell Ranger](https://img.shields.io/badge/compatible-CellRanger%20BAMs-green)]()

---
## 🔥 Why IsoDecipher? (Core Advantages)
### 1. Solving the Noise Problem: True Signal Extraction
* **The Problem**: Existing de-novo tools (e.g., scAPA, Sierra) often fail to distinguish True Biological Signals from Internal Priming artifacts. Research indicates that 30-50% of raw 3' clusters in scRNA-seq can be technical noise (A-rich internal sequences), leading to high False Discovery Rates (FDR) and unintelligible feature names like chr14:105741338-105741400.

* **The IsoDecipher Solution**: Inspired by the modeling principles of MAAPER (Li et al., 2021), we implement Biologically Anchored Mapping. By guiding read assignment with a curated GTF-based "Targeted Panel," we filter out intronic noise and internal priming at the source, ensuring every count represents a validated transcript isoform.

### 2. Zero-Ambiguity Indexing: Canonical Ordering
* **The Problem**: Standard tools output raw genomic coordinates. When analyzing negative-strand genes, researchers must manually flip their logic (where a larger coordinate actually means a more proximal site), which is highly error-prone in large-scale datasets.

* **The IsoDecipher Solution**: We introduce Strand-Aware Grouping.

* **Group 0**  is ALWAYS the Proximal site (closest to the TSS/Shortest UTR), regardless of the strand.

* This standardized indexing allows for direct, automated calculation of the **PUI (Proximal Usage Index)** and consistent cross-sample comparisons.


```text
Forward Strand (+):  TSS ---> [Exon] ---(Group 0: Proximal)---(Group 1: Distal)--->
Reverse Strand (-):  <---(Group 1: Distal)---(Group 0: Proximal)--- [Exon] <--- TSS

Index Logic: Group 0 is ALWAYS the closest to TSS (Shortest UTR).
```

### 3. Semantic Labeling: User-Defined Context
* **Industry-First Customization**: Unlike other tools that provide abstract IDs, IsoDecipher allows users to inject Semantic Labels during the panel building stage.

* **Practical Value**: You can pre-label specific clusters (e.g., Secreted vs. Membrane for IGs). When the assignment script runs, the resulting count matrix automatically carries these biological names. This makes the data "Analysis-Ready" the moment it is generated, skipping hours of manual coordinate cross-referencing.

---

### Summary Comparison Table

| Feature | De-novo Tools (scAPA/Sierra) | **IsoDecipher** |
| :--- | :--- | :--- |
| **Feature Naming** | Random coordinates | **Gene Symbol + Group Index** |
| **Biological Context** | None (Post-hoc mapping needed) | **Built-in (Group 0 = Proximal)** |
| **Noise Handling** | High (Includes internal priming) | **Low (GTF-guided filtering)** |
| **Specialized Logic** | General purpose | **Curated IG (Secreted/Membrane) labels** |
| **Analysis Readiness** | Requires heavy cleaning | **Scanpy/VAE-ready output** |


## 🎯 Key Features

* **Distance-based PolyA Clustering**: Automatically groups transcript ends into discrete `polyA_groups` using a strand-aware algorithm.
* **Immune Cell Specialization**: Built-in logic to distinguish **Secreted** vs. **Membrane-bound** isoforms for Immunoglobulin genes (e.g., `IGHM`, `IGHG1-4`).
* **Customizable Precision**: Support for gene-specific clustering tolerances to handle complex genomic loci.
* **scRNA-seq Optimized**: Designed to work directly with Cell Ranger output (BAM files with `CB` and `UB` tags).

---

## 🚀 Installation

We recommend using **Miniforge (Mamba)** for fast and stable dependency management.

```zsh
# Create a dedicated environment
mamba create -n iso_decipher python=3.10 -y
mamba activate iso_decipher

# Install dependencies
pip install gffutils pysam pandas numpy scanpy anndata matplotlib seaborn scipy
```

---

## 🛠️ Workflow

### Step 1: Build Feature Panel
Convert a standard GTF annotation into a quantifiable isoform feature table.

```zsh
python IsoDecipher/scripts/build_panel_features.py \
    --gtf data/Homo_sapiens.GRCh38.115.gtf \
    --genes data/gene_list.txt \
    --out results/panel_features.csv
```

### Step 2: Read Assignment (In Development)
Assign reads from BAM files to the specific `polyA_groups` defined in Step 1.

### Step 3: VAE Analysis (Planned)
Utilize a Variational Autoencoder (VAE) to cluster cells based on their isoform usage profiles.

---

## 📊 Output Schema (`panel_features.csv`)

| Column | Description |
| :--- | :--- |
| `gene` | Gene symbol (e.g., CD44, IGHG1) |
| `polyA_group` | Index of the site (0 = most proximal to TSS) |
| `rep_coord` | Representative genomic coordinate (cluster mean) |
| `user_label` | Custom biological context (e.g., Secreted, Membrane, Distal_Exon, Brain_Specific) |
| `avg_spliced_utr` | Mean length of the 3' UTR across transcripts in the group |

---

## 📚 Technical Logic

### Strand-Aware Ordering
IsoDecipher ensures that `polyA_group 0` is always the **proximal** site (closest to the start of the gene), regardless of whether the gene is on the positive or negative strand. This allows for consistent analysis of "UTR lengthening" or "shortening" across different samples.

### Automated Biological Labeling
For curated Immunoglobulin genes, the tool automatically maps the first cluster (proximal) as the Secreted isoform and distal clusters as Membrane-bound. For other genes, users can manually define labels in the user_label column to facilitate downstream analysis.

---

## 📄 License
MIT License © 2026 Rene Yu-Hong Cheng
