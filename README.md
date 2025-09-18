# IsoDecipher 🧬

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Python 3.10](https://img.shields.io/badge/python-3.10-blue)]()
[![Build with Cell Ranger / STARsolo](https://img.shields.io/badge/build-CellRanger%2FSTARsolo-orange)]()

**IsoDecipher: Revealing immune and cancer cell states through membrane vs secreted isoform usage in 3′ scRNA-seq.**

---

## 🚀 Overview
IsoDecipher is a pipeline for **deciphering isoform usage** (membrane vs secreted, soluble vs transmembrane) directly from **3′ single-cell RNA-seq**.  
By focusing on alternative last exons (ALEs) and alternative polyadenylation (APA), IsoDecipher uncovers **hidden immune and cancer cell states** that standard pipelines miss.

---

## 🔹 Features
- 🧾 **Custom GTF Builder**: add isoform entries (e.g., `IGHM-201` = secreted, `IGHM-202` = membrane).  
- 🧬 **Isoform Quantification**: count matrices at transcript/isoform resolution.  
- 🔍 **Isoform Usage Analysis**: distinguish functional states (e.g., antibody-secreting vs membrane-bound B cells).  
- 🌐 **Immune & Cancer Focus**: curated panel (IGH isoforms, CD79A, TNFSF10, S100 family, IL7R, CD44, VEGFA, etc.).  
- 📊 **Visualization**: UMAPs, violin/dot plots of isoform usage.  

---

## 📂 Repository Structure
- in progress
---

## 🔧 Installation
```bash
# clone the repo
git clone https://github.com/<your-username>/IsoDecipher.git
cd IsoDecipher

# create environment
conda create -n isodecipher python=3.10
conda activate isodecipher

# install dependencies
pip install -r requirements.txt
