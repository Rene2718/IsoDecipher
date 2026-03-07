
# IsoDecipher 

**IsoDecipher** is a targeted isoform quantification tool designed for **3' single-cell RNA-seq (scRNA-seq)** data.

It recovers hidden biological signals—such as **Alternative Polyadenylation (APA)** and **B-cell antibody isoform switching**—by re-interpreting **3' read distributions** that standard pipelines collapse into gene-level counts.

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)  
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue)]()  
[![Compatible with Cell Ranger](https://img.shields.io/badge/compatible-CellRanger%20BAMs-green)]()

---

### Development Status

⚠️ IsoDecipher is currently under active development.

Current progress:

- Panel builder ✔
- Read assignment ✔
- Downstream analysis ⏳

---

## Concept

Standard scRNA-seq pipelines collapse transcript isoforms into **gene-level counts**.

```
scRNA-seq reads
       │
       ▼
   Cell Ranger
   (gene counts)

Cell        IGHM
--------------------
cell_1        12
cell_2        18
cell_3         9
```

However, **isoform information is lost**.

IsoDecipher instead **reinterprets the genomic positions of 3' reads** to recover polyA site usage.

Example: **IGHM gene**

```
TSS → Exon ────────┬───────────────>
                   │
              polyA site 1
              (Secreted)

TSS → Exon ────────┬────────────────────────>
                   │
              polyA site 2
              (Membrane)
```

Recovered counts:

```
Cell      IGHM_G0   IGHM_G1
-----------------------------
cell_1       10        2
cell_2        4       14
cell_3        8        1
```

This enables detection of:

- Isoform switching
- Alternative Polyadenylation (APA)
- Secreted vs membrane immunoglobulin transcripts

at **single-cell resolution**.

---

## Why IsoDecipher?

### 1. Solving the Noise Problem: True Signal Extraction

### The Problem

Existing de-novo APA tools (e.g. **scAPA**, **Sierra**) identify polyA sites through read clustering **without transcript annotation**.

While powerful for discovery, these approaches may introduce:

- internal priming artifacts in A-rich regions
- false positive clusters
- unintelligible feature names such as

```
chr14:105741338-105741400
```

### IsoDecipher Solution

Inspired by **MAAPER (Li et al., 2021)**, IsoDecipher uses a **Biologically Anchored Mapping strategy**.

Reads are assigned to **GTF-derived polyA groups**, which:

- filters intronic noise
- removes internal priming artifacts
- ensures every count corresponds to a **validated transcript isoform**

---

### 2. Zero-Ambiguity Indexing

Most tools output **raw genomic coordinates**, which becomes confusing for **negative strand genes**.

Example problem:

```
larger coordinate ≠ distal site
```

IsoDecipher introduces **strand-aware canonical ordering**.

```
Forward Strand (+)

TSS ---> [Exon] ---(Group 0)---(Group 1)--->

Reverse Strand (-)

<---(Group 1)---(Group 0)--- [Exon] <--- TSS
```

Rule:

```
Group 0 = ALWAYS the proximal polyA site
```

This enables consistent calculation of the **Proximal Usage Index (PUI)**.

---

### 3. Semantic Labeling

Most APA tools output **abstract site IDs**.

IsoDecipher supports **user-defined semantic labels**.

Example:

| gene | group | user_label |
|-----|------|------|
| IGHM | 0 | Secreted |
| IGHM | 1 | Membrane |

Resulting features:

```
IGHM_G0_secreted
IGHM_G1_membrane
```

The count matrix becomes **biologically interpretable immediately**.

---

## Comparison

| Feature | De-novo Tools | IsoDecipher |
|---|---|---|
Feature naming | Genomic coordinates | Gene + Group index |
Biological context | None | Built-in |
Noise handling | High | GTF-guided |
Immunoglobulin logic | None | Secreted/Membrane |
Analysis readiness | Heavy cleaning | Scanpy ready |

---

##  Key Features

- Distance-based polyA clustering
- Strand-aware site ordering
- Immunoglobulin isoform detection
- Custom gene panels
- Cell Ranger BAM compatibility
- Scanpy / AnnData friendly outputs

---

##  Installation

Recommended environment manager: **Miniforge / Mamba**

```
mamba create -n iso_decipher python=3.10 -y
mamba activate iso_decipher

pip install gffutils pysam pandas numpy scanpy anndata matplotlib seaborn scipy
```

---

##  Quick Start

### Step 1: Build isoform feature panel

Convert a standard GTF annotation into a quantifiable isoform feature table.

```
python IsoDecipher/scripts/build_panel_features.py  --gtf data/Homo_sapiens.GRCh38.115.gtf  --genes data/gene_list.txt  --out results/panel_features.csv
```

Output:

```
panel_features.csv
```

---

### Step 2: Read Assignment (In Development)

Reads from **Cell Ranger BAM files** will be assigned to the predefined `polyA_groups`.

Input:

```
Cell Ranger BAM
(CB + UB tags)
```

Output:

```
Isoform Count Matrix
cells × polyA_groups
```

---

### Step 3: Downstream Analysis (Planned)

Isoform usage matrices can be directly used with:

- Scanpy
- AnnData
- Variational Autoencoders (VAE)
- machine learning pipelines

---

## Workflow Overview

```
GTF annotation + gene panel
            │
            ▼
     Panel Builder
(build_panel_features.py)
            │
            ▼
   panel_features.csv
            │
            ▼
Cell Ranger BAM
(CB + UB tags)
            │
            ▼
    IsoDecipher
   Read Assignment
            │
            ▼
Isoform Count Matrix
(cells × polyA_groups)
            │
            ▼
Scanpy / VAE / ML analysis
```

---

##  Output Schema

`panel_features.csv`

| Column | Description |
|---|---|
| gene | Gene symbol |
| polyA_group | Index (0 = proximal) |
| rep_coord | Representative coordinate |
| user_label | Biological label |
| avg_spliced_utr | Average UTR length |

---

##  Technical Notes

### Strand-Aware Ordering

IsoDecipher guarantees:

```
polyA_group 0 = proximal site
polyA_group 1 = distal site
```

independent of genomic strand orientation.

---

### Immunoglobulin Isoforms

For curated immunoglobulin loci:

```
Group 0 → Secreted isoform
Group 1 → Membrane isoform
```

This enables detection of **antibody secretion switching** directly from 3' scRNA-seq data.

---

##  License

MIT License  
© 2026 Rene Yu-Hong Cheng
