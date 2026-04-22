# IsoDecipher

**IsoDecipher** is a targeted isoform quantification tool designed for **3' single-cell RNA-seq (scRNA-seq)** data.

It recovers hidden biological signals—such as **Alternative Polyadenylation (APA)** and **B-cell antibody isoform switching**—by re-interpreting **3' read distributions** that standard pipelines collapse into gene-level counts.

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue)]()
[![Compatible with Cell Ranger](https://img.shields.io/badge/compatible-CellRanger%20BAMs-green)]()

---

## Development Status

⚠️ IsoDecipher is currently under active development.

| Step | Status |
|------|--------|
| Panel builder | ✅ Complete |
| Annotation filter | ✅ Complete |
| Read assignment | ✅ Complete |
| Downstream analysis (PUI/PSI/Entropy) | ⏳ In progress |

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

**Example: IGHM gene**

```
TSS → Exon ────────┬─────>
                   │
              polyA site 1 (Secreted)

TSS → Exon ────────┬────────────────────────>
                   │
              polyA site 2 (Membrane)
```

Recovered counts:

```
Cell      IGHM_G0_Secreted   IGHM_G1_Membrane
----------------------------------------------
cell_1          10                  2
cell_2           4                 14
cell_3           8                  1
```

This enables detection of isoform switching, Alternative Polyadenylation (APA), and secreted vs membrane immunoglobulin transcripts at **single-cell resolution**.

---

## Why IsoDecipher?

### 1. Biologically Anchored Noise Suppression

Existing de-novo APA tools (e.g. **scAPA**, **Sierra**) identify polyA sites through read clustering **without transcript annotation**, introducing internal priming artifacts and unintelligible feature names (`chr14:105741338-105741400`).

Inspired by **MAAPER (Li et al., 2021)**, IsoDecipher uses **GTF-derived polyA groups** to:
- Filter intronic noise and internal priming artifacts
- Ensure every count corresponds to a **validated transcript isoform**
- Produce biologically interpretable feature names

### 2. Annotation Quality Filter

IsoDecipher automatically filters poorly annotated polyA groups during panel construction:

- Skips **retained introns**, **nonsense-mediated decay**, and **CDS-not-defined** transcripts at the GTF parsing stage
- Removes **zero-UTR singleton groups** (single transcript with no UTR annotation) from the final panel
- Always preserves the **dominant group** (most supporting transcripts) per gene

This reduces panel noise from ~47% zero-UTR features to under 6%, with no loss of biologically meaningful sites.

Use `--no-filter` to disable filtering and retain all groups.

### 3. Zero-Ambiguity Strand-Aware Indexing

IsoDecipher introduces **strand-aware canonical ordering** so that Group 0 always refers to the proximal polyA site, regardless of genomic strand:

```
Forward (+):  TSS ---> [Exon] ---(G0)---(G1)--->
Reverse (-):  <---(G1)---(G0)--- [Exon] <--- TSS
```

This enables consistent calculation of the **Proximal Usage Index (PUI)**.

### 4. Adaptive Downstream Metrics

IsoDecipher automatically selects the appropriate quantification metric based on the number of polyA groups per gene:

| Groups | Method | Description |
|--------|--------|-------------|
| 2 | **PUI** (Proximal Usage Index) | `G0 / (G0 + G1)` — directional, simple |
| 3–5 | **PSI** (Proximal Shift Index) | Weighted average of group index |
| 6+ | **Entropy** | Shannon entropy of polyA site distribution |

### 5. Semantic Labeling

IsoDecipher supports **user-defined biological labels**. For immunoglobulin loci, labels are assigned automatically:

| gene | group | user_label |
|------|-------|------------|
| IGHM | G0 | Secreted |
| IGHM | G1 | Membrane |

---

## Comparison

| Feature | De-novo Tools | IsoDecipher |
|---------|--------------|-------------|
| Feature naming | Genomic coordinates | Gene + Group index |
| Biological context | None | Built-in |
| Noise handling | High artifact rate | GTF-guided, filtered |
| Annotation filter | None | Biotype + UTR quality filter |
| Immunoglobulin logic | None | Secreted/Membrane auto-label |
| Analysis readiness | Heavy cleaning required | Scanpy ready |

---

## Scientific Scope

> **"Is APA switching a bifurcation event, and does this distinguish normal differentiation from oncogenic transformation?"**

IsoDecipher is designed to address this question by combining:

- **APA trajectory analysis** — PUI/PSI/Entropy along pseudotime
- **Bifurcation detection** — bimodality coefficient along differentiation
- **Waddington landscape reconstruction** — potential landscape from polyA site distributions
- **Information-theoretic analysis** — Shannon entropy, JSD, mutual information, transfer entropy
- **Cell stability scoring** — combining Markov-based DEn with APA variance

### Technical Limitation

IsoDecipher quantifies **3' end polyA site usage** from genomic coordinates. It does not resolve internal exon skipping or full-length transcript isoforms, which require complementary long-read platforms (PacBio MAS-seq, Oxford Nanopore). Within its scope, IsoDecipher provides GTF-anchored, noise-suppressed APA quantification at single-cell resolution.

---

## Installation

Recommended: **Miniforge / Mamba**

```bash
mamba create -n iso_decipher python=3.10 -y
mamba activate iso_decipher

pip install gffutils pysam pandas numpy scanpy anndata \
            matplotlib seaborn scipy scikit-learn
```

For pipeline orchestration:

```bash
mamba install -c conda-forge -c bioconda snakemake -y
```

---

## Quick Start

### Step 1: Build isoform feature panel

```bash
python IsoDecipher/scripts/build_panel_features.py \
    --gtf data/Homo_sapiens.GRCh38.115.gtf \
    --genes data/gene_list.txt \
    --out results/panel_features.csv
```

Options:

```
--custom_params   Custom clustering tolerance per gene (TSV)
--no-filter       Disable zero-UTR singleton filter
```

Output summary example:

```
[FILTER] Removed 407 zero-UTR singleton groups
[SUMMARY] Group count distribution:
  1 group  (no analysis): 54 genes
  2 groups (PUI):         75 genes
  3-5 groups (PSI):      106 genes
  6+ groups (Entropy):    44 genes
```

### Step 2: Assign reads

```bash
python IsoDecipher/scripts/assign_reads.py \
    --bam data/samples/exp93/possorted_genome_bam.bam \
    --panel results/panel_features.csv \
    --out results/counts/exp93_isoform_count.csv
```

### Step 3: Run full pipeline with Snakemake

```bash
snakemake --cores 4
```

The Snakemake pipeline automatically:
- Builds the panel once
- Assigns reads for all samples in parallel
- Logs errors per sample to `logs/`

---

## Workflow Overview

```
GTF annotation + gene panel (289 genes)
            │
            ▼
     Panel Builder
(build_panel_features.py)
    + Annotation Filter
            │
            ▼
   panel_features.csv
   (960 features, 279 genes)
            │
            ▼
  Cell Ranger BAM × N samples
      (CB + UB tags)
            │
            ▼
    IsoDecipher Read Assignment
    (strand-aware, UMI-deduplicated)
            │
            ▼
  Isoform Count Matrix
  (cells × polyA_groups)
            │
            ▼
  PUI / PSI / Entropy Analysis
  Pseudotime + Bifurcation Detection
  Waddington Landscape
```

---

## Output Schema

### panel_features.csv

| Column | Description |
|--------|-------------|
| gene | Gene symbol |
| polyA_group | Index (0 = proximal) |
| rep_coord | Representative coordinate |
| strand | Genomic strand |
| chrom | Chromosome |
| avg_spliced_utr | Average spliced UTR length |
| avg_genomic_utr | Average genomic UTR length |
| num_transcripts | Supporting transcripts per group |
| transcript_ids | Ensembl transcript IDs |
| user_label | Biological label (e.g. Secreted/Membrane) |

### isoform_count.csv

Sparse matrix: `cells × polyA_group_features`

Feature naming convention: `{GENE}_G{index}_{label}`

Example: `IGHM_G0_Secreted`, `IGHM_G1_Membrane`

---

## Gene Panel

IsoDecipher ships with a curated panel of **289 genes** across 19 biological categories:

- Housekeeping controls (ACTB, GAPDH, B2M)
- Immunoglobulins (IGHM, IGHG1-4, IGHA1-2, IGHE)
- Plasma cell & secretory pathway (XBP1, MZB1, HSPA5, FKBP11...)
- Plasma cell transcription factors (IRF4, PRDM1, PAX5, BCL6...)
- APA & RNA processing machinery (CSTF1-3, NUDT21, CPSF1-4, FIP1L1...)
- BCR signaling & surface markers (CD79A/B, CD19, CD22, BTK, SYK...)
- T cell markers (CD4, CD8A/B, FOXP3, TOX...)
- Immune checkpoints (CD274/PD-L1, PDCD1, CTLA4, LAG3, TIGIT...)
- Cancer metabolism & proliferation (PKM, MYC, TP53, EGFR, KRAS...)
- Adhesion & EMT (CD44, ITGA/B family, ICAM family, VIM, FN1...)
- EMT transcription factors (ZEB1/2, TWIST1/2, SNAI1/2...)
- S100 family (S100A1-16, S100B, S100P, S100Z)
- SEC-seq discovered markers (CD59, NEAT1, SDC1, TXNIP, MIF...)

---

## License

MIT License
© 2026 Rene Yu-Hong Cheng
