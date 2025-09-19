# IsoDecipher 🧬

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Python 3.10](https://img.shields.io/badge/python-3.10-blue)]()
[![Build with Cell Ranger / STARsolo](https://img.shields.io/badge/build-CellRanger%2FSTARsolo-orange)]()

**IsoDecipher: Revealing immune and cancer cell states through membrane vs secreted isoform usage and 3′UTR dynamics in 3′ scRNA-seq.**

---

## 🚀 Overview
IsoDecipher is a pipeline for **deciphering hidden isoform usage patterns** from **3′ single-cell RNA-seq**.  
By focusing on **alternative last exons (ALEs)** and **alternative polyadenylation (APA)**, IsoDecipher uncovers **immune and cancer cell states** that standard pipelines miss.

IsoDecipher provides two complementary modes:
- 🧾 **Panel Mode (annotation-based)** → custom GTF + curated gene list; computes isoform usage and 3′UTR lengths for selected immune/cancer genes.  
- 🌍 **Genome-wide APA Mode (discovery)** → optional integration with tools like DaPars2 / scDaPars2 to infer polyA usage across *all genes*.  

---

## 🔹 Features
- 🧾 **Isoform Curation Module**: curate immune/cancer isoforms (membrane vs secreted, soluble vs receptor).  
- 🧬 **Isoform Quantification**: per-cell transcript counts using Cell Ranger / STARsolo / Kallisto / Alevin-fry.  
- 📏 **3′UTR Dynamics**: compute weighted average UTR length per cell and detect shortening/lengthening events.  
- 🌐 **Immune & Cancer Focus**: IGH family, CD79A, TNFSF10, S100 family, IL7R, CD44, VEGFA, and more.  
- 📊 **Visualization Tools**: UMAPs, violin/dot plots, isoform usage heatmaps, pseudotime trajectories.  
- 🌍 **APA Discovery Mode**: genome-wide scan for APA shifts using external packages (DaPars2, scDaPars2, scAPA).  

---

## 📂 Repository Structure
```
IsoDecipher/
├── data/                     # Input/output reference files
│   ├── gene_list.txt         # Curated immune/cancer panel
│   ├── gencode.v44.gtf       # Annotation (not stored in repo)
│   └── example_counts.tsv    # Example counts
│
├── scripts/                  # Python utilities
│   ├── build_panel.py        # Auto-generate isoform panel CSV
│   ├── build_gtf.py          # Create custom GTF (membrane vs secreted)
│   └── compute_utr_length.py # Calculate UTR lengths from isoform counts
│
├── notebooks/                # Example Jupyter workflows
│   └── demo_analysis.ipynb
│
├── docs/                     # Documentation, diagrams
├── README.md                 # You are here
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

### 1. Panel Mode (curated genes)
```bash
python scripts/build_panel.py --gtf data/gencode.v44.gtf --genes data/gene_list.txt --out data/isoform_panel.csv
```

- Builds an isoform panel (with UTR lengths) for your selected genes.  
- Run quantification with Cell Ranger / Kallisto / Alevin-fry.  
- Analyze isoform ratios + UTR shortening in clusters.  

### 2. Genome-wide APA Mode (optional)
```bash
# Example with DaPars2/scDaPars2 (not bundled in repo)
python run_dapars2.py --bam sample.bam --gtf gencode.v44.gtf --out genomewide_APA_results.txt
```

- Computes **Percentage Distal Usage Index (PDUI)** or equivalent metrics across *all genes*.  
- Allows discovery of novel APA-regulated genes to expand your panel.  

---

## 📊 Example Applications
- **Plasma B cells**: IGHM/IGHG isoform switch (membrane → secreted) + UTR shortening.  
- **T cells**: Soluble receptor isoforms (IL7R, TNFRSF1A) distinguish activation states.  
- **Tumor subsets**: Global 3′UTR shortening in CD44, VEGFA, S100 family genes.  

---

## 🗺️ Roadmap
- [ ] Finalize curated isoform panel (IGH, TNFSF, CD79A, CD44, etc.)  
- [ ] Implement automated UTR extraction from GTF (compute_utr_length.py)  
- [ ] Add integration hooks for DaPars2/scDaPars2 (APA discovery mode)  
- [ ] Build visualization notebooks (UMAP, violin, heatmap)  
- [ ] Release preprint + example dataset  

---

## 📜 License
MIT License © 2025 Rene Yu-Hong Cheng
