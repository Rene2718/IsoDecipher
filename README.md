# IsoDecipher 🧬

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Python 3.10](https://img.shields.io/badge/python-3.10-blue)]()
[![Build with Cell Ranger / STARsolo](https://img.shields.io/badge/build-CellRanger%2FSTARsolo-orange)]()

**IsoDecipher: Revealing immune and cancer cell states through membrane vs secreted isoform usage and 3′UTR dynamics in 3′ scRNA-seq.**

---

## 🚀 Overview
Most single-cell RNA-seq pipelines (Cell Ranger, STARsolo) collapse all isoforms into a single *gene count*, discarding transcript-level information.  
**IsoDecipher** is a lightweight framework that works *on top of existing Cell Ranger outputs* to reconstruct isoform usage patterns by leveraging:

- **Alternative last exons (ALEs)**  
- **Alternative polyadenylation (APA)**  
- **Membrane vs secreted isoforms** of immune and cancer-related genes  

This enables per-cell isoform biology without re-aligning FASTQs.

---

## 🔹 Why two “panel builders”?
IsoDecipher separates **annotation** from **quantification** so you can choose the right tool for your analysis:

| Script | Purpose | Output | Why it matters |
|--------|---------|--------|----------------|
| `build_panel.py` | Annotation of transcripts | Table of isoforms, transcript IDs, **UTR lengths** | Human-readable reference — lets you inspect which isoforms exist, how long their 3′UTRs are, and which are biologically interesting. Useful for figures and reports. |
| `build_panel_features.py` | Features for isoform assignment | Collapsed **polyA groups** (`GENE_polyA1`, `GENE_polyA2`) with last exon + polyA windows | Machine-readable features — tells IsoDecipher how to assign UMIs to isoform groups when parsing BAM files. Scales from simple genes (IGHM: secreted vs membrane) to complex ones (CD44, S100 with dozens of isoforms). |

**In short:**  
- Use **`build_panel.py`** when you want to *understand the isoforms themselves*.  
- Use **`build_panel_features.py`** when you want to *actually count isoforms per cell* from BAM.  

---


## 🔹 Key Features
- 🧾 **PolyA grouping**: collapse transcripts with the same (or nearby) polyA ends into interpretable groups (`GENE_polyA1`, `GENE_polyA2`).  
- 🧬 **Isoform Quantification**: parse Cell Ranger BAMs (CB/UB tags) and assign UMIs to isoform groups.  
- 📏 **3′UTR Metrics**: compute weighted average UTR length per cell, detect shortening/lengthening events.  
- 📊 **Visualization Ready**: outputs simple matrices/CSV that integrate directly into Scanpy or Seurat.  
- 🌍 **Genome-wide APA option**: extend panel to all genes, or integrate with DaPars2 / scDaPars2 for discovery.  

---

## 🔹 Important Notes
- **Use Cell Ranger with the *standard* GTF** (Ensembl/GENCODE). Do **not** give Cell Ranger the custom isoform GTF.  
- IsoDecipher re-mines the BAM (`possorted_genome_bam.bam`), which contains **all aligned reads with CB/UB tags**, even those not assigned to any gene in the standard matrix.  
- The custom panels you build (`panel_features.csv`) are used only in the IsoDecipher quantifier, not in Cell Ranger.  
- `gene_list.txt` can include comments starting with `#` — these lines will be 

---

## Diagram
```
FASTQ
   │
   ▼
Cell Ranger count  (with standard GTF)
   │
   ├── Gene count matrix (standard)
   │
   └── possorted_genome_bam.bam
            │
            ▼
      IsoDecipher
        ├── build_panel.py         → UTR lengths (annotation)
        ├── build_panel_features.py → polyA groups (short/long for IGH)
        └── quantify_isoforms_from_bam.py
               │
               ▼
        Isoform-aware counts + UTR usage
```


## 📂 Repository Structure
```
IsoDecipher/
├── data/                             # Input/output reference files
│   ├── gene_list.txt                 # Curated immune/cancer panel
│   ├── Homo_sapiens.GRCh38.115.gtf   # Annotation (not stored in repo)
│   └── example_counts.tsv            # Example counts
│
├── scripts/                          # Python utilities
│   ├── build_panel.py                # Isoform + UTR annotation table
│   ├── build_panel_features.py       # PolyA group + IGH short/long feature builder
│   └── quantify_isoforms_from_bam.py # Assign UMIs to isoform groups
│
├── notebooks/                        # Example Jupyter workflows
│   └── demo_analysis.ipynb
│
├── docs/                             # Documentation, diagrams
├── README.md                         # You are here
├── LICENSE
└── requirements.txt
```

---

## 🔧 Installation
```bash
git clone https://github.com/rene2718/IsoDecipher.git
cd IsoDecipher

conda create -n isodecipher python=3.10
conda activate isodecipher

pip install -r requirements.txt
```

---

## 🧪 Usage

### 1. Annotate isoforms and UTR lengths
```bash
python scripts/build_panel.py \
  --gtf Homo_sapiens.GRCh38.115.gtf \
  --genes data/gene_list.txt \
  --out data/isoform_panel.csv
```
Output: table of transcripts with UTR lengths and IDs (for inspection/plots).

### 2. Build isoform groups for quantification
```bash
python scripts/build_panel_features.py \
  --gtf Homo_sapiens.GRCh38.115.gtf \
  --genes data/gene_list.txt \
  --out data/panel_features.csv
```
Output: feature table of polyA groups (`GENE_polyA1`, `GENE_polyA2`), with last-exon and polyA-window intervals.  
For IGH constant region genes (`IGHM`, `IGHG1–4`, `IGHA1–2`, `IGHE`), IsoDecipher automatically labels groups as **`_short`** vs **`_long`**.

---

### 🔹 Demo gene list
To quickly test IsoDecipher with IGH short/long labeling, you can start with a minimal gene list:

```text
# Demo gene list for IsoDecipher
# Includes IGH genes (short/long labeling) and a control gene (CD44)
IGHM
IGHG1
CD44
```

Save this as `demo_gene_list.txt` and run:

```bash
python scripts/build_panel_features.py   --gtf Homo_sapiens.GRCh38.115.gtf   --genes demo_gene_list.txt   --out demo_panel_features.csv
```

**Expected output:**  
- IGHM is split into `IGHM_short` and `IGHM_long`.  
- IGHG1 is split into `IGHG1_short` and `IGHG1_long` (with the `IGHG1-203` exception labeled `short`).  
- CD44 remains labeled as `CD44_polyA1`, `CD44_polyA2`, etc.

---

### 3. Quantify isoform usage from Cell Ranger BAM
```bash
python scripts/quantify_isoforms_from_bam.py \
  --bam cellranger/outs/possorted_genome_bam.bam \
  --panel data/panel_features.csv \
  --genes data/gene_list.txt \
  --out_prefix results/igh
```
---

## 📊 Outputs

| File | Description |
|------|-------------|
| `*_cell_x_polyA_counts.csv` | UMI counts per isoform group per cell |
| `*_cell_x_gene_isoform_fraction.csv` | Isoform fractions per gene per cell |
| `*_isoform_qc.tsv` | QC: UMIs assigned by tier (unique last exon vs polyA window vs ambiguous) |

---

## 📊 Example Outputs

**`build_panel.py` (annotation table)**  
```
gene,transcript_id,transcript_name,utr_length
IGHM,ENST00000390655,IGHM-201,1127
IGHM,ENST00000610629,IGHM-202,310
```

**`build_panel_features.py` (polyA groups for quantification)**  
```
gene,polyA_group,feature_type,feature_id,chrom,start,end,strand,transcripts,transcript_names,avg_utr_length,min_utr_length,max_utr_length,evidence_tier
IGHM,IGHM_short,polyA_window,IGHM_short_window,chr14,105122000,105122600,+,ENST00000390655;ENST00000610629,IGHM-201;IGHM-202,620,310,930,tier3_polyA
IGHM,IGHM_short,last_exon_group,IGHM_short_exon_ENST00000390655,chr14,105121800,105122000,+,ENST00000390655,IGHM-201,620,310,930,tier1_group
```

---

## 🔹 IsoDecipher vs. Genome-wide UTR Tools

| Tool | Scope | Input | Output | Strengths | Limitations |
|------|-------|-------|--------|------------|-------------|
| **IsoDecipher** (this repo) | *Panel-based* (immune/cancer gene list, e.g. IGHM, CD44) | Cell Ranger BAM + panel_features.csv | - `*_cell_x_polyA_counts.csv` → UMI counts per isoform group per cell<br>- `*_cell_x_gene_isoform_fraction.csv` → isoform usage fractions<br>- `*_isoform_qc.tsv` → QC by evidence tier<br>- Weighted **UTR length per cell** (from avg UTRs) | - Very interpretable (e.g. “IGHM_short vs IGHM_long”)<br>- Keeps CB/UB context<br>- Works directly on Cell Ranger BAMs you already have | Not genome-wide by default (limited to curated panel unless expanded) |
| **scUTRquant** | *Genome-wide* | FASTQs or BAM + GTF | - Cell × isoform count matrix (similar to Cell Ranger output)<br>- Isoform fractions<br>- Integrates directly into Seurat/Scanpy | - Lightweight (kallisto|bustools pseudoalignment)<br>- Designed for scRNA-seq 3′ end capture<br>- Gives genome-wide isoform quantification | Requires rerun from FASTQ (not just BAM)<br>- Less flexible naming/grouping (IDs not curated like IsoDecipher) |
| **scDaPars2** | *Genome-wide* | BAM + GTF | - **PDUI matrix**: Percentage Distal polyA site Usage Index per gene per cell/cluster<br>- Genome-wide APA calls | - Extends DaPars (well-known APA tool) to single-cell<br>- Produces PDUI which is standard APA metric | Heavier compute<br>- More statistical modeling, less direct UMI counting<br>- Output is PDUI (relative usage) rather than raw counts |
| **scAPA / scAPAtrap** | *Genome-wide* | BAM | - Detected polyA sites<br>- Counts per polyA site per cell<br>- APA scores (per gene/cell) | - Call novel polyA sites directly from scRNA-seq<br>- Genome-wide | Often noisier polyA detection<br>- More complex workflows |

---

## 📊 Example Applications
- **Plasma B cells**: IGHM/IGHG isoform switch (membrane → secreted) + UTR shortening.  
- **T cells**: Soluble receptor isoforms (IL7R, TNFRSF1A) distinguish activation states.  
- **Tumor subsets**: Global 3′UTR shortening in CD44, VEGFA, S100 family genes.  

## 💻 System Requirements
IsoDecipher is designed to be lightweight and run comfortably on a **local workstation or laptop**:

- **OS**: Linux / macOS (Windows with WSL)  
- **CPU**: 2–4 cores recommended  
- **RAM**: 8–16 GB is sufficient for panels of 50–200 genes  
- **Disk**: <5 GB for typical outputs (depends on BAM size)  

Why so light?  
- IsoDecipher does **not** re-align FASTQs.  
- It works directly from the `possorted_genome_bam.bam` output of Cell Ranger.  
- Only the **genes in your panel** are processed (instead of genome-wide).  

Typical runtime on a laptop: **minutes to an hour** depending on panel size.  

By contrast, genome-wide APA tools (DaPars2, scDaPars2, scUTRquant) may require HPC or cloud environments because they re-process all reads across the transcriptome.

---

## 🗺️ Roadmap
- [ ] Refine last-exon uniqueness (clip overlaps).  
- [ ] Expand gene panels beyond B cells (CD44, IL7R, VEGFA).  
- [ ] Genome-wide APA discovery integration (DaPars2, scDaPars2).  
- [ ] Build visualization notebooks (UMAP, violin, heatmap).  
- [ ] Release preprint + example dataset.  

---

## 📜 License
MIT License © 2025 Rene Yu-Hong Cheng
