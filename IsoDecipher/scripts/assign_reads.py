import pysam
import pandas as pd
from collections import defaultdict

panel = pd.read_csv("results/panel_features.csv")
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


bam_path = "data/possorted_genome_bam.bam"
bam = pysam.AlignmentFile(bam_path, "rb")

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

df = pd.DataFrame(results)
matrix = df.pivot(index='cell_barcode',columns='feature',values = 'count').fillna(0)
matrix.to_csv("results/isoform_count.csv")

var_info = panel.copy()
var_info['feature'] = var_info.apply(lambda r: f"{r.gene}_G{r.polyA_group}_{r.user_label}", axis=1)
var_info.set_index('feature').to_csv("results/feature_metadata.csv")