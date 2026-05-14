# IsoDecipher — Master Scientific Framework

## The Central Question

> **"Is APA switching a bifurcation event, and does this distinguish normal differentiation from oncogenic transformation?"**

---

## Chapter 1: What is IsoDecipher?

A targeted isoform quantification tool for **3' scRNA-seq data** that recovers:
- Alternative Polyadenylation (APA) — which polyA site does each cell use?
- 3' UTR length dynamics — short vs long UTR per cell
- Proximal Usage Index (PUI) — per cell, per gene

### Technical Limitation (be honest in paper):
- ✅ polyA site usage, 3' UTR length, APA
- ❌ internal exon skipping, cassette exons, full-length isoforms
- ⚠️ PUI is based on genomic coordinate UTR, not complete transcript

**This is also a feature:** GTF-anchored assignment is more accurate and less noisy than de novo APA tools (scAPA, Sierra).

> Paper framing: "IsoDecipher quantifies APA and 3' UTR usage from genomic coordinates. Internal splicing events require complementary long-read approaches."

---

## Chapter 2: Biological Observations

### IGHM — Gradual shift (NOT bifurcation)
```
Naive B cell:     G1 (membrane) dominant
Activated B cell: G0 + G1 coexist
Plasma cell:      G0 (secreted) dominant, G1 still present
```
→ Gradual, reversible. Plasma cells still need membrane IgM for BCR signaling.

### IGHG — Bifurcation-like (commitment event)
```
Memory B cell:  G1 (membrane) only
Plasma cell:    G0 (secreted) completely dominant, G1 disappears
```
→ Irreversible. True commitment event.

### CD44 — Complex, cancer-relevant
- Normal activation: gradual PUI changes with cell cycle
- Cancer/EMT: PUI shifts correlate with invasive state
- May show bifurcation in EMT context

### CD274 (PD-L1) — Cancer-specific
- Normal: activation-regulated, reversible
- Cancer: persistent expression → possible bifurcation lock-in
- Key immunotherapy target — APA affects mRNA stability

---

## Chapter 3: Theoretical Framework

### 3.1 Waddington Landscape

Cell state x evolves according to:
```
dx/dt = -∇U(x) + η(t)
```
- U(x) = potential landscape ("hills and valleys")
- -∇U(x) = deterministic drift toward attractors
- η(t) = stochastic noise (gene expression noise)

**Key insight:** RNA velocity estimates dx/dt directly → we can reconstruct ∇U(x).

| Cell State | Landscape | PUI | Stability |
|-----------|-----------|-----|-----------|
| Naive/stem | Shallow valley | Variable | Medium |
| Transitioning | Ridge/saddle | Bimodal | Lowest |
| Terminal | Deep valley | Low variance | Highest |
| Cancer EMT | Between attractors | Fluctuating | Very low |
| Proliferating | Outside attractor | High noise | Low |

### 3.2 Bifurcation Theory Applied to APA

A bifurcation = qualitative change in system behavior from small parameter change.

**PUI bifurcation signature:**
```
Before: unimodal PUI distribution
At:     bimodal PUI distribution (two stable states)  ← bifurcation point
After:  unimodal PUI distribution (new committed state)
```

**Detection:**
- Bimodality Coefficient (BC) > 0.555
- Shannon entropy peaks at transition
- PUI variance peaks then drops sharply

### 3.3 Markov Chain & Differentiation Entropy (DEn)

```
DEn(i) = -Σ_j T(i,j) * log(T(i,j))
```

| DEn | APA variance | Interpretation |
|-----|-------------|----------------|
| Low | Low | TRUE terminal cell (deep attractor) |
| High | High | TRUE transitioning cell (saddle point) |
| Low | High | Fate committed, APA still resolving (lag!) |
| High | Low | APA committed first, fate still plastic |

**The last two quadrants reveal temporal ordering — a testable original hypothesis.**

---

## Chapter 4: The Central Hypothesis

### Gene-specific APA bifurcation classification:

| Class | Example Genes | Biology | Reversible? |
|-------|-------------|---------|-------------|
| Bifurcation | IGHG, CD274 in cancer | True commitment | No |
| Gradual | IGHM, CD44 normal | Continuous regulation | Yes |
| Bimodal stable | Some surface markers | Two coexisting states | Partially |
| Noise | Most genes | No APA regulation | N/A |

### The Cancer Hypothesis:
> "Oncogenic APA switching is irreversible (bifurcation-like), while normal differentiation APA switching can be gradual or reversible."

**Clinical implications:**
- Bifurcation point = therapeutic window (intervene before commitment)
- APA bifurcation signature = early epigenetic cancer biomarker

---

## Chapter 5: Information Theory Layer

| Method | Measures | IsoDecipher Application |
|--------|----------|------------------------|
| Shannon Entropy | Diversity/uncertainty | Cell state heterogeneity in APA |
| JSD | Distribution distance | APA divergence between cell types |
| Mutual Information | Statistical dependency | Which TFs drive APA switching? |
| Transfer Entropy | Causal direction in time | Does APA precede or follow fate decision? |
| Total Correlation | Multi-variable shared info | Transcriptome + secretome + APA co-regulation |

**Unique to your data (SEC-seq + IsoDecipher):**
Total Correlation across all three modalities simultaneously — no other dataset or tool can do this.

### Implementation Priority:
1. **JSD** — easiest, most interpretable, direct paper figure
2. **Shannon Entropy** — one function, immediate biological insight
3. **Mutual Information** — find APA machinery coupling (CSTF2, NUDT21)
4. **Total Correlation** — unique to SEC-seq + IsoDecipher
5. **Transfer Entropy** — most novel, most publication-worthy

---

## Chapter 6: Dimensionality Reduction Reference

| Method | Linear? | Global structure | Pseudotime? | Use for |
|--------|---------|-----------------|-------------|---------|
| PCA | ✅ Yes | ✅ Good | Barely | Initial exploration |
| Diffusion Map | ❌ No | ✅ Good | ✅ Designed for this | Pseudotime calculation |
| UMAP | ❌ No | ⚠️ Partial | ⚠️ Not reliable | Visualization |
| t-SNE | ❌ No | ❌ Poor | ❌ No | Rare cell populations |

**Diffusion Map math:** Random walk on KNN graph. Pseudotime = diffusion distance from root cell. Non-linear because distance follows data manifold, not straight line.

**For your data:**
- No velocity data → use DPT (Diffusion Pseudotime)
- Have velocity → use CellRank + scVelo latent time

---

## Chapter 7: Data Strategy

### Your own data:
- Ex vivo B cell culture → plasma cell differentiation
- Best for: IGHM, IGHG APA dynamics
- Limitation: ex vivo ≠ in vivo

### Public cancer datasets:
- **GSE139555** — Tumor infiltrating B cells, multiple cancers
- **GSE181061** — Breast cancer single-cell
- **10x Genomics public** — processed h5ad, ready to use
- Key genes: CD274, CD44, MYC, ITGA series

### BAM files:
- Check CB + UB tags → run velocyto → scVelo
- No tags → use DPT

---

## Chapter 8: Adaptive Metric Strategy

Not all genes behave the same way. IsoDecipher automatically selects the appropriate metric based on the number of valid polyA groups per gene.

### Panel composition (after annotation filter):
- 54 genes: 1 group — no APA analysis possible
- **75 genes: 2 groups → PUI**
- **106 genes: 3-5 groups → PSI**
- **44 genes: 6+ groups → Entropy**

Total: **225 genes with meaningful APA signal**

---

### Method 1: PUI — Proximal Usage Index (2 groups)

```python
def compute_pui(adata, gene, pseudocount=1e-10):
    g0 = adata[:, f'{gene}_G0'].X.toarray().flatten()
    g1 = adata[:, f'{gene}_G1'].X.toarray().flatten()
    return g0 / (g0 + g1 + pseudocount)
# PUI → 1: proximal dominant (short UTR)
# PUI → 0: distal dominant (long UTR)
```

Best for: IGHM, IGHG (secreted/membrane switching)

---

### Method 2: PSI — Proximal Shift Index (3-5 groups)

```python
def compute_psi(adata, gene, groups):
    """Weighted average of group index, normalized to 0-1."""
    cols = [f'{gene}_G{i}' for i in range(groups)]
    counts = np.array([adata[:, c].X.toarray().flatten() for c in cols]).T
    total = counts.sum(axis=1, keepdims=True) + 1e-10
    probs = counts / total
    weights = np.arange(groups) / (groups - 1)  # [0, ..., 1]
    return (probs * weights).sum(axis=1)
# PSI → 0: proximal dominant
# PSI → 1: distal dominant
```

Best for: CD274 (3-5 groups), most checkpoint genes

---

### Method 3: Entropy (6+ groups)

```python
def compute_apa_entropy(adata, gene, groups):
    """Shannon entropy of polyA site distribution per cell."""
    cols = [f'{gene}_G{i}' for i in range(groups)]
    counts = np.array([adata[:, c].X.toarray().flatten() for c in cols]).T
    total = counts.sum(axis=1, keepdims=True) + 1e-10
    probs = counts / total
    return -np.sum(probs * np.log(probs + 1e-10), axis=1)
# High entropy = heterogeneous APA = transitioning/plastic cell
# Low entropy = committed to one site = stable cell state
```

Best for: CD44 (~12 groups), MEF2C (~17 groups), HNRNPC (~14 groups)

---

### Gene-specific strategy table:

| Gene | Groups | Method | Biological relevance |
|------|--------|--------|---------------------|
| IGHM | 2 | PUI | Secreted/Membrane switching |
| IGHG1-4 | 2 | PUI | Class switch commitment |
| CD274 | 3-5 | PSI | PD-L1 checkpoint dynamics |
| EGFR | ~6 | PSI/Entropy | Cancer proliferation |
| CD44 | ~12 | Entropy | EMT plasticity |
| MEF2C | ~17 | Entropy | TF regulatory complexity |
| HNRNPC | ~14 | Entropy | APA machinery self-regulation |

---

### Bifurcation detection by method:

| Method | Bifurcation signature |
|--------|----------------------|
| PUI | Bimodality Coefficient > 0.555 along pseudotime |
| PSI | Variance peak + bimodal PSI distribution |
| Entropy | Entropy peak then sharp drop (commitment) |

---

## Chapter 8b: Key Code Modules

### Bimodality Coefficient (Bifurcation Detection)
```python
def bimodality_coefficient(x):
    from scipy.stats import skew, kurtosis
    n, s, k = len(x), skew(x), kurtosis(x)
    return (s**2 + 1) / (k + 3*(n-1)**2 / ((n-2)*(n-3)))
# BC > 0.555 → bimodal → bifurcation point
```

### Potential Landscape Reconstruction
```python
def reconstruct_potential(pui_values, n_bins=50):
    hist, edges = np.histogram(pui_values, bins=n_bins, range=(0,1), density=True)
    U = -np.log(hist + 1e-10)
    U = U - U.min()
    return (edges[:-1] + edges[1:]) / 2, U  # bin_centers, potential
```

### Combined Stability Score
```python
# DEn (Markov) + APA variance (IsoDecipher)
# Low DEn + Low APA var  = terminal (both agree)
# High DEn + High APA var = transitioning (both agree)
# Low DEn + High APA var  = fate committed, APA resolving (lag)
# High DEn + Low APA var  = APA committed, fate still plastic
```

---

## Chapter 9: Paper Story

### Title options:
1. *"APA switching is a gene-specific bifurcation event: single-cell evidence from IsoDecipher"*
2. *"Reconstructing the Waddington landscape from APA dynamics in single-cell RNA-seq"*
3. *"IsoDecipher: information-theoretic dissection of APA dynamics distinguishes normal differentiation from oncogenic transformation"*

### Figure roadmap:
1. **Figure 1** — IsoDecipher tool overview, IGHM example
2. **Figure 2** — APA trajectory along pseudotime (B cell differentiation)
3. **Figure 3** — Bifurcation detection: IGHM (gradual) vs IGHG (bifurcation)
4. **Figure 4** — Waddington landscape reconstruction from PUI
5. **Figure 5** — Combined stability score quadrant (DEn + APA variance)
6. **Figure 6** — Cancer application: CD274, CD44 APA in tumor microenvironment
7. **Figure 7** — Transfer entropy: temporal ordering of APA vs fate commitment

---

## Chapter 10: Implementation Roadmap

### Current status (Week of Apr 22):
- ✅ build_panel_features.py — complete with biotype + UTR filter
- ✅ assign_reads.py — complete with UMI deduplication
- ✅ Snakemake pipeline — build_panel + assign_reads rules
- ✅ gene_list.txt — 289 genes, 19 categories
- ✅ panel_features.csv — 960 features, 279 genes, 5.6% zero UTR
- ⏳ GCS upload — BAM files uploading to gs://isodecipher-bam/samples/bcell/
- ⏳ assign_reads — pending GCS setup

### Week 1 — Core downstream analysis
| Day | Task |
|-----|------|
| Wednesday | PUI/PSI/Entropy notebook — basic per-cell metrics |
| Thursday | APA trajectory — pseudotime + smoothing + JSD |
| Friday | Bifurcation detection — bimodality coefficient |

### Week 2 — Advanced analysis
| Day | Task |
|-----|------|
| Monday | Waddington landscape visualization |
| Tuesday | Gene classification (bifurcation vs gradual vs noise) |
| Wednesday | CellRank integration + DEn computation |
| Thursday | Combined stability quadrant (DEn + APA variance) |
| Friday | Public cancer dataset download + IsoDecipher application |

### Week 3 — Publication preparation
| Day | Task |
|-----|------|
| Monday | Information theory metrics (MI, Transfer Entropy) |
| Tuesday | Total Correlation (transcriptome + secretome + APA) |
| Wednesday | Package IsoDecipher (pyproject.toml, pip installable) |
| Thursday | GitHub tutorial notebooks |
| Friday | bioRxiv preprint draft outline |

---

## Chapter 11: Why This Is Competitive

### Your unique advantages:
1. **SEC-seq inventor** — only you have transcriptome + secretome simultaneously
2. **IsoDecipher** — APA at single-cell resolution
3. **ODE modeling background** — attractor theory is natural to you (Nat Commun 2022)
4. **Flow cytometry weekly** — can validate computational predictions experimentally
5. **B cell biology expertise** — know exactly what IGHM/IGHG switching means
6. **10x Genomics collaboration** — hands-on platform experience

### What no one else can do:
Total Correlation across transcriptome + secretome + APA:
- Requires SEC-seq (yours) + IsoDecipher (yours)
- This analysis is impossible without both tools

### Novel algorithms:
- Bifurcation detection in APA = **novel**
- Gene-specific APA classification = **novel**
- Combined DEn + APA stability score = **novel**
- Transfer entropy for temporal ordering of APA vs fate = **novel**

None of these exist in current tools.
