#!/usr/bin/env python3
"""
IsoDecipher: Quantify isoforms from Cell Ranger BAM
---------------------------------------------------
Assign UMIs to isoform groups defined in panel_features.csv
"""

import argparse
import pandas as pd
import pysam
from collections import defaultdict

def load_features(panel_file):
    df = pd.read_csv(panel_file)
    features = defaultdict(list)
    for _, row in df.iterrows():
        features[row["gene"]].append(row)
    return features

def quantify_isoforms(bam_file, panel_file, out_prefix):
    features = load_features(panel_file)
    bam = pysam.AlignmentFile(bam_file, "rb")

    umi_dict = defaultdict(set)  # (cell, gene, group) -> set of UMIs
    evidence_counts = defaultdict(int)

    for read in bam.fetch(until_eof=True):
        if not read.has_tag("CB") or not read.has_tag("UB"):
            continue
        cb = read.get_tag("CB")
        ub = read.get_tag("UB")
        pos = read.reference_end if not read.is_reverse else read.reference_start
        chrom = bam.get_reference_name(read.reference_id)
        strand = "-" if read.is_reverse else "+"

        for gene, feats in features.items():
            for f in feats:
                if f["chrom"] == chrom and f["strand"] == strand:
                    if f["start"] <= pos <= f["end"]:
                        group = f["polyA_group"]
                        tier = f["evidence_tier"]
                        umi_dict[(cb, gene, group)].add(ub)
                        evidence_counts[tier] += 1

    # Build cell x group count matrix
    rows = []
    for (cb, gene, group), umis in umi_dict.items():
        rows.append([cb, gene, group, len(umis)])
    counts = pd.DataFrame(rows, columns=["cell", "gene", "polyA_group", "UMIs"])
    counts.to_csv(f"{out_prefix}_cell_x_polyA_counts.csv", index=False)

    # Fractions per gene
    frac_rows = []
    for cell, df_cell in counts.groupby("cell"):
        for gene, df_gene in df_cell.groupby("gene"):
            total = df_gene["UMIs"].sum()
            for _, r in df_gene.iterrows():
                frac = r["UMIs"] / total if total > 0 else 0
                frac_rows.append([cell, gene, r["polyA_group"], frac])
    fractions = pd.DataFrame(frac_rows, columns=["cell", "gene", "polyA_group", "fraction"])
    fractions.to_csv(f"{out_prefix}_cell_x_gene_isoform_fraction.csv", index=False)

    # QC
    qc = pd.DataFrame(list(evidence_counts.items()), columns=["evidence_tier", "count"])
    qc.to_csv(f"{out_prefix}_isoform_qc.tsv", sep="\t", index=False)

    print(f"[IsoDecipher] Written outputs with prefix {out_prefix}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Quantify isoform usage from Cell Ranger BAM")
    ap.add_argument("--bam", required=True, help="Cell Ranger BAM (possorted_genome_bam.bam)")
    ap.add_argument("--panel", required=True, help="Panel features CSV")
    ap.add_argument("--out_prefix", required=True, help="Output file prefix")
    args = ap.parse_args()

    quantify_isoforms(args.bam, args.panel, args.out_prefix)
