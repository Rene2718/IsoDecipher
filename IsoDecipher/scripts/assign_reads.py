import pysam
import pandas as pd
from collections import defaultdict
import argparse
import os

# 1. Setup Argument Parsing
parser = argparse.ArgumentParser(description="Assign reads to isoforms per sample")
parser.add_argument("--bam", required=True, help="Path to input BAM file")
parser.add_argument("--panel", default="results/panel_features.csv", help="Path to panel features CSV")
parser.add_argument("--out", required=True, help="Path to save output CSV")
args = parser.parse_args()

panel = pd.read_csv(args.panel)
targets = defaultdict(list)

for row in panel.itertuples(index=False):
    targets[row.chrom].append({
        'gene': row.gene,
        'pos': row.rep_coord,
        'strand': row.strand,
        'group': row.polyA_group,
        'label': row.user_label,
        'spliced_utr': row.avg_spliced_utr,
        'genomic_utr': row.avg_genomic_utr,
    })


bam = pysam.AlignmentFile(args.bam, "rb")

counts = defaultdict(set)
window = 200
umi_best_match ={}

for chrom, sites in targets.items():
    current_chrom = chrom
    if current_chrom not in bam.references:
        if f"chr{current_chrom}" in bam.references:
            current_chrom = f"chr{current_chrom}"
        else:
            print(f"⚠️ Warning: Contig {chrom} not found in BAM. Skipping...")
            continue
    
    print(f"Current Chromosome: {current_chrom} | Total unique UMIs captured so far: {len(umi_best_match)}")
    
    for site in sites:
        for read in bam.fetch(current_chrom, site['pos']-window, site['pos']+window):
            if not (read.has_tag("CB") and read.has_tag("UB")):
                continue
            if (site["strand"] == "+" and read.is_reverse) or (site["strand"] == "-" and not read.is_reverse):
                continue
            read_3_prime = read.reference_end if site['strand'] == '+' else read.reference_start
            dist = abs(read_3_prime - site['pos'])


            if dist <= window:
                cb = read.get_tag("CB")
                ub = read.get_tag("UB")
                umi_key = (cb, ub)
                feature = f"{site['gene']}_G{site['group']}_{site['label']}"
                
                if umi_key not in umi_best_match or dist < umi_best_match[umi_key][1]:
                    umi_best_match[umi_key] = (feature, dist)


for (cb, ub), (feature, dist) in umi_best_match.items():
    counts[(cb, feature)].add(ub)


results = []
for (cb, feature), umis in counts.items():
    results.append({
        'cell_barcode': cb, 
        'feature': feature,
        'count' : len(umis)
    })

if results:
    df = pd.DataFrame(results)
    matrix = df.pivot(index='cell_barcode',columns='feature',values = 'count').fillna(0)
    matrix.to_csv(args.out)
    print(f"✅ Success: {args.out} saved.")
else:
    print(f"⚠️ No reads assigned for {args.bam}")

