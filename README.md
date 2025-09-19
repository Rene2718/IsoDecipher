# IsoDecipher 🧬

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Python 3.10](https://img.shields.io/badge/python-3.10-blue)]()
[![Build with Cell Ranger / STARsolo](https://img.shields.io/badge/build-CellRanger%2FSTARsolo-orange)]()

**IsoDecipher: Revealing immune and cancer cell states through membrane vs secreted isoform usage and 3′UTR dynamics in 3′ scRNA-seq.**

---

## 🚀 Overview
IsoDecipher is a pipeline for **deciphering hidden isoform usage patterns** from **3′ single-cell RNA-seq**.  
By focusing on **alternative last exons (ALEs)** and **alternative polyadenylation (APA)**, IsoDecipher uncovers **immune and cancer cell states** that standard pipelines miss.

IsoDecipher combines:
- 🧾 **Isoform Curation Module** → custom GTF with membrane vs secreted isoforms (e.g., IGHM-201 vs IGHM-202).  
- 📏 **UTR Module** → automated quantification of 3′UTR length dynamics using transcriptome annotations and APA-aware tools.  

---

## 🔹 Features
- 🧾 **Custom GTF Builder**: curate immune/cancer isoforms (membrane vs secreted, soluble vs receptor).  
- 🧬 **Isoform Quantification**: per-cell transcript counts using Cell Ranger / STARsolo / Kallisto / Alevin-fry.  
- 📏 **3′UTR Dynamics**: compute weighted average UTR length per cell and detect shortening/lengthening events.  
- 🌐 **Immune & Cancer Focus**: IGH family, CD79A, TNFSF10, S100 family, IL7R, CD44, VEGFA, and more.  
- 📊 **Visualization Tools**: UMAPs, violin/dot plots, isoform usage heatmaps, pseudotime trajectories.  

---

## 📂 Repository Structure
```
IsoDecipher/
├── data/              # Example custom GTF/FASTA, UTR length tables
├── scripts/           # build_gtf.py, compute_utr_length.py, analysis utils
├── notebooks/         # Demo Jupyter notebooks (Scanpy, visualization)
├── docs/              # Documentation
├── README.md          # You are here
├── LICENSE
└── requirements.txt
```

---

## 🔧 Installation
```bash
git clone https://github.com/<your-username>/IsoDecipher.git
cd IsoDecipher

conda create -n isodecipher python=3.10
conda activate isodecipher

pip install -r requirements.txt
```

---

## 🧪 Usage

### 1. Isoform Curation (Membrane vs Secreted)
```bash
python scripts/build_gtf.py --genes immune_panel.csv --gtf gencode.v44.gtf --out custom.gtf
cellranger mkref --genome=custom_ref --fasta custom.fa --genes custom.gtf
```

### 2. 3′UTR Length Module
```bash
python scripts/compute_utr_length.py --gtf gencode.v44.gtf --genes utr_panel.csv --counts isoform_counts.tsv --out utr_length_table.csv
```

Outputs per-cell average UTR length and isoform ratios.

---

## 📊 Example Applications
- **B cells**: IGHM/IGHG isoform switch (membrane → secreted) + UTR shortening.  
- **T cells**: Soluble receptor isoforms (IL7R, TNFRSF1A) distinguish activation states.  
- **Tumor subsets**: Global 3′UTR shortening in CD44, VEGFA, S100 family genes.  

---

## 🗺️ Roadmap
- [ ] Finalize curated isoform panel (IGH, TNFSF, CD79A, etc.)  
- [ ] Implement automated UTR extraction from GTF (compute_utr_length.py)  
- [ ] Build visualization notebooks (UMAP, violin, heatmap)  
- [ ] Benchmark against PolyASite / scDaPars2  
- [ ] Release preprint + example dataset  

---

## 📜 License
MIT License © 2025 Rene Yu-Hong Cheng
