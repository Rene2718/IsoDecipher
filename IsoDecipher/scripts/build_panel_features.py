#!/usr/bin/env python3
"""
IsoDecipher: Build feature panel for isoform quantification
-----------------------------------------------------------
From an Ensemble GTF + list of genes, collapse isoforms into
polyA groups (shared 3'ends) and emit features for transcript 
assignment from Cell Ranger BAMs.

Features:
 - polyA_group: merge window around transcript ends
 - last_exon_group: last exon coordinate for transcripts in group
 - avg_utr_length: mean 3'UTR length across transcripts in group
 - min_utr_length: max_utr_length: spread of UTR lengths

Usage:
python IsoDecipher/scripts/build_panel_features.py \
    --gtf data/Homo_sapiens.GRCh38.115.gtf \
    --genes data/gene_list.txt \
    --out results/panel_features.csv \
    [--custom_params custom.tsv]
    [--no-filter]

By default, zero-UTR singleton groups are filtered out.
Use --no-filter to keep all groups.
"""

import argparse
import gffutils
import pandas as pd
from collections import defaultdict
import os

# Biotypes to skip during transcript collection
SKIP_BIOTYPES = {
    "retained_intron",
    "nonsense_mediated_decay",
    "misc_RNA",
    # lncRNA and processed_transcript removed:
    # genes like NEAT1 are lncRNA with meaningful APA signal
}

# Biotypes that don't have CDS by definition — don't penalize them
NO_CDS_BIOTYPES = {
    "lncRNA",
    "processed_transcript",
    "non_stop_decay",
    "sense_intronic",
    "sense_overlapping",
    "antisense",
}

IG_WHITELIST = {"IGHM", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHE"}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build cleavage-centered polyA panel (distance-based model)"
    )
    parser.add_argument("--gtf", required=True, help="Input GTF file")
    parser.add_argument("--genes", required=True, help="Gene list file (one gene per line)")
    parser.add_argument("--db", help="Optional: specific path for gffutils DB")
    parser.add_argument("--out", required=True, help="Output CSV")
    parser.add_argument("--custom_params", help="Custom tolerance TSV/CSV")
    parser.add_argument("--no-filter", action="store_true",
                        help="Disable zero-UTR singleton filter (keep all groups)")
    return parser.parse_args()


def load_custom_parameters(file):
    """
    Load user-defined clustering tolerance per gene.
    Return dict: {gene: end_tolerance}
    """
    if not file:
        return {}
    df = pd.read_csv(file, sep=None, engine="python")
    if not {"gene", "end_tolerance"} <= set(df.columns):
        raise ValueError("Custom parameter must contain columns: gene, end_tolerance")
    param_dict = df.set_index("gene")["end_tolerance"].to_dict()
    print(f"[CUSTOM] Loaded {len(param_dict)} gene-specific tolerances")
    return param_dict  # fixed: was missing return


def load_or_build_db(gtf):
    db_path = gtf + ".db"
    build = False

    if os.path.exists(db_path):
        if os.path.getmtime(gtf) > os.path.getmtime(db_path):
            print("[IsoDecipher] GTF newer than DB. Rebuilding...")
            build = True
        else:
            print("[IsoDecipher] Using existing DB")
    else:
        print("[IsoDecipher] No DB found. Building...")
    if build:
        gffutils.create_db(
            gtf,
            dbfn=db_path,
            force=True,
            keep_order=True,
            merge_strategy="merge",
            disable_infer_genes=True,
            disable_infer_transcripts=True
        )
    return gffutils.FeatureDB(db_path, keep_order=True)


def collect_transcript_end(db, gene_list):
    """
    Collect strand-aware transcript 3' boundaries for selected genes.
    Skips non-coding biotypes (retained intron, NMD, etc.)
    and transcripts with no CDS (CDS not defined).

    Return:
        gene_data = {
            gene_name: [
                {
                    coord,
                    transcript_id,
                    transcript_name,
                    utr_length,
                    chrom,
                    strand
                }
            ]
        }
    """
    gene_data = defaultdict(list)
    all_genes = {
        g.attributes.get("gene_name", [g.id])[0].upper(): g
        for g in db.features_of_type("gene")
    }
    missing_genes = []
    skipped_biotype = 0
    skipped_no_cds = 0

    for gene_name in gene_list:
        gene = all_genes.get(gene_name)

        if gene is None:
            missing_genes.append(gene_name)
            continue

        chrom = gene.seqid
        transcripts = db.children(gene, featuretype="transcript")

        for tx in transcripts:

            tx_id = tx.attributes.get("transcript_id", [""])[0]
            tx_name = tx.attributes.get("transcript_name", [""])[0]
            strand = tx.strand

            # --- Biotype filter ---
            biotype = tx.attributes.get("transcript_biotype", [""])[0]
            if biotype in SKIP_BIOTYPES:
                skipped_biotype += 1
                continue

            # --- CDS check ---
            # Skip CDS-not-defined for protein coding transcripts
            # but allow lncRNA/processed_transcript (they have no CDS by definition)
            cds = list(db.children(tx, featuretype="CDS"))
            if not cds and biotype not in NO_CDS_BIOTYPES:
                skipped_no_cds += 1
                continue

            exons = list(db.children(tx, featuretype="exon"))

            # strand-aware 3' end coordinate
            if exons:
                if strand == "+":
                    last_exon = max(exons, key=lambda e: e.end)
                    coord = last_exon.end
                else:
                    last_exon = min(exons, key=lambda e: e.start)
                    coord = last_exon.start
            else:
                coord = tx.end if tx.strand == "+" else tx.start

            # Compute genomic UTR length (only for protein coding)
            genomic_utr_dist = None
            if cds:
                if strand == "+":
                    cds_end = max(c.end for c in cds)
                    genomic_utr_dist = coord - cds_end
                else:
                    cds_start = min(c.start for c in cds)
                    genomic_utr_dist = cds_start - coord

            # Compute spliced UTR length
            utr_features = list(db.children(tx, featuretype="three_prime_UTR"))
            if utr_features:
                spliced_utr_len = sum(u.end - u.start + 1 for u in utr_features)
            elif genomic_utr_dist is not None:
                spliced_utr_len = genomic_utr_dist
            else:
                spliced_utr_len = None  # lncRNA — no UTR by definition

            gene_data[gene_name].append({
                "coord": coord,
                "transcript_id": tx_id,
                "transcript_name": tx_name,
                "genomic_utr_length": genomic_utr_dist,
                "spliced_utr_length": spliced_utr_len,
                "chrom": chrom,
                "strand": strand
            })

        # Sanity check: all transcripts same strand
        strands = set(t['strand'] for t in gene_data[gene_name])
        if len(strands) > 1:
            print(f"[ERROR] Gene {gene_name} has transcripts on multiple strands: {strands}")

    if missing_genes:
        print(f"\n[WARN] The following {len(missing_genes)} genes were not found in the GTF:")
        print(f"       {', '.join(missing_genes)}")
        print(f"       Please check your gene list spelling or GTF version.\n")

    print(f"[FILTER] Skipped {skipped_biotype} transcripts by biotype "
          f"(retained_intron, NMD, etc.)")
    print(f"[FILTER] Skipped {skipped_no_cds} transcripts with no CDS annotation")

    return gene_data


def cluster_transcript_ends(transcript_data, gene_name, param_dict, default_tolerance=50):
    """
    Groups transcript ends based on a gene-specific or default tolerance.

    Returns:
        list of groups (each group is list of transcript dicts)
    """
    tolerance = param_dict.get(gene_name.upper(), default_tolerance)
    sorted_tx = sorted(transcript_data, key=lambda x: x['coord'])

    if not sorted_tx:
        return []

    clusters = []
    current_cluster = [sorted_tx[0]]

    for i in range(1, len(sorted_tx)):
        dist = sorted_tx[i]["coord"] - sorted_tx[i-1]["coord"]
        if dist <= tolerance:
            current_cluster.append(sorted_tx[i])
        else:
            clusters.append(current_cluster)
            current_cluster = [sorted_tx[i]]
    clusters.append(current_cluster)

    if gene_name.upper() in param_dict:
        print(f"  [CLUSTER] {gene_name}: Applying custom tolerance of {tolerance}bp")

    return clusters


def filter_panel_features(df_panel):
    """
    Remove poorly annotated polyA groups.

    Rules:
    1. Drop groups where avg_spliced_utr == 0 AND num_transcripts == 1
       (likely retained intron, CDS incomplete, or CDS not defined)
    2. Always keep the group with the most transcripts per gene
       (even if it has zero UTR — this is the dominant/canonical site)
    3. Keep groups with multiple transcripts even if zero UTR
       (multiple supporting transcripts = more likely real site)

    Returns filtered DataFrame with re-indexed polyA_group per gene.
    """
    filtered_rows = []
    removed = 0

    for gene, gdf in df_panel.groupby('gene'):
        # Find the dominant group (most transcripts) — always keep
        dominant_idx = gdf['num_transcirpts'].idxmax()

        for idx, row in gdf.iterrows():
            is_zero_utr = row['avg_spliced_utr'] == 0
            is_singleton = row['num_transcirpts'] == 1
            is_dominant = (idx == dominant_idx)

            if is_dominant:
                filtered_rows.append(row)
            elif is_zero_utr and is_singleton:
                removed += 1  # drop
            else:
                filtered_rows.append(row)

    df_filtered = pd.DataFrame(filtered_rows).reset_index(drop=True)

    # Re-index polyA_group per gene (0, 1, 2... after filtering)
    df_filtered['polyA_group'] = df_filtered.groupby('gene').cumcount()

    # Re-assign IG labels after re-indexing
    ig_mask = df_filtered['gene'].str.upper().isin(IG_WHITELIST)
    df_filtered.loc[ig_mask, 'user_label'] = df_filtered.loc[ig_mask, 'polyA_group'].apply(
        lambda i: "Secreted" if i == 0 else "Membrane"
    )

    print(f"\n[FILTER] Removed {removed} zero-UTR singleton groups")
    print(f"[FILTER] Kept {len(df_filtered)} groups ({len(df_panel)} before filter)")

    return df_filtered


def print_panel_summary(df):
    """Print group count distribution for downstream method selection."""
    groups_per_gene = df.groupby('gene')['polyA_group'].count()
    print(f"\n[SUMMARY] Group count distribution:")
    print(f"  1 group  (no analysis): {(groups_per_gene == 1).sum()} genes")
    print(f"  2 groups (PUI):         {(groups_per_gene == 2).sum()} genes")
    print(f"  3+ groups (Entropy):    {(groups_per_gene >= 3).sum()} genes")


def main():
    args = parse_args()

    # Load DB and gene list
    db = load_or_build_db(args.gtf)
    if not os.path.exists(args.genes):
        raise FileNotFoundError(f"Gene list file not found: {args.genes}")
    with open(args.genes, "r") as f:
        gene_list = [line.strip().upper()
                     for line in f
                     if line.strip() and not line.strip().startswith("#")]

    # Load custom gene tolerance
    custom_params = load_custom_parameters(args.custom_params)

    print(f"[IsoDecipher] Collecting transcript ends for {len(gene_list)} genes...")
    raw_gene_data = collect_transcript_end(db, gene_list)

    panel_rows = []

    for gene_name in gene_list:
        transcripts = raw_gene_data.get(gene_name, [])
        if not transcripts:
            continue

        # Annotation-level clustering
        clusters = cluster_transcript_ends(transcripts, gene_name, custom_params)

        # Strand-aware sorting: group_0 = always proximal
        strand = transcripts[0]['strand']
        if strand == "+":
            clusters.sort(key=lambda c: min(tx['coord'] for tx in c))
        else:
            clusters.sort(key=lambda c: max(tx['coord'] for tx in c), reverse=True)

        # Build feature rows
        for i, cluster in enumerate(clusters):
            coords = [tx['coord'] for tx in cluster]
            rep_coord = int(sum(coords) / len(coords))

            spliced_utr_lengths = [tx['spliced_utr_length'] for tx in cluster
                                   if tx['spliced_utr_length'] is not None]
            avg_spliced_utr = (sum(spliced_utr_lengths) / len(spliced_utr_lengths)
                               if spliced_utr_lengths else 0)

            genomic_utr_lengths = [tx['genomic_utr_length'] for tx in cluster
                                   if tx['genomic_utr_length'] is not None]
            avg_genomic_utr = (sum(genomic_utr_lengths) / len(genomic_utr_lengths)
                               if genomic_utr_lengths else 0)

            row = {
                "gene": gene_name,
                "polyA_group": i,
                "rep_coord": rep_coord,
                "strand": strand,
                "chrom": cluster[0]['chrom'],
                "avg_spliced_utr": round(avg_spliced_utr, 2),
                "avg_genomic_utr": round(avg_genomic_utr, 2),
                "num_transcirpts": len(cluster),
                "transcript_ids": ";".join([tx['transcript_id'] for tx in cluster]),
                "transcript_names": ";".join([tx['transcript_name'] for tx in cluster])
            }

            if gene_name.upper() in IG_WHITELIST:
                row["user_label"] = "Secreted" if i == 0 else "Membrane"
            else:
                row["user_label"] = "N/A"

            panel_rows.append(row)

    df_panel = pd.DataFrame(panel_rows)

    # Apply filter unless --no-filter flag is set
    if not args.no_filter:
        df_panel = filter_panel_features(df_panel)
    else:
        print("\n[FILTER] Skipping filter (--no-filter flag set)")

    # Print summary
    print_panel_summary(df_panel)

    # Save
    df_panel.to_csv(args.out, index=False)

    print(f"\n[SUCCESS] IsoDecipher Panel complete!")
    print(f" - Output file: {args.out}")
    print(f" - Total features: {len(df_panel)}")
    print(f" - Genes processed: {df_panel['gene'].nunique()}")


if __name__ == "__main__":
    main()
