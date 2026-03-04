# IsoDecipher 🧬

**IsoDecipher** is a targeted isoform quantification tool designed for 3' single-cell RNA-seq (scRNA-seq) data. It recovers hidden biological signals—such as Alternative Polyadenylation (APA) and B-cell antibody switching—by re-interpreting 3' read patterns that standard pipelines typically collapse into gene-level counts.

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)  
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue)]()  
[![Compatible with Cell Ranger](https://img.shields.io/badge/compatible-CellRanger%20BAMs-green)]()

---

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
| `ig_label` | Antibody class label (`Secreted`, `Membrane`, or `N/A`) |
| `avg_spliced_utr` | Mean length of the 3' UTR across transcripts in the group |

---

## 📚 Technical Logic

### Strand-Aware Ordering
IsoDecipher ensures that `polyA_group 0` is always the **proximal** site (closest to the start of the gene), regardless of whether the gene is on the positive or negative strand. This allows for consistent analysis of "UTR lengthening" or "shortening" across different samples.

### IG Labeling Heuristics
For curated Immunoglobulin genes, the tool automatically maps the first cluster (proximal) as the **Secreted** isoform and distal clusters as **Membrane-bound**, reflecting the biological mechanism of antibody production in Plasma cells.

---

## 📄 License
MIT License © 2026 Rene Yu-Hong Cheng
