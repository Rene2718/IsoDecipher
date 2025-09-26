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
 - ave_utr_length: mean 3'UTR length across transcripts in group
 - min_utr_length: max_utr_length: spread of UTR lengths

Usage:
python build_panel_features.py \\\\
    --gtf Homo_sapiens.GRCh38.115.gtf \\\\
    --genes data/gene_list.txt \\\\
    --out panel_features.csv

"""
import argparse
import gffutils
import pandas as pd

IG_WHITELIST = {"IGHM", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHE"}

def get_ig_label(gene, tx_name):
    """Return short/long labels for immunoglobulin transcripts."""
    if gene not in IG_WHITELIST:
        return None
    if tx_name.endswith("-201"):
        return "short"
    elif tx_name.endswith("-202"):
        return "long"
    elif tx_name == "IGHG1-203":
        return "short"
    return None



def build_panel_features(gtf, gene_list_file, out_csv, polyA_window=300, end_tolerance=20):
   
    #read gene list
    with open(gene_list_file) as f:
        gene_list= [line.strip() for line in f if line.strip() and not line.startswith("#")]

    #load GTF into gffutils database
    db = gffutils.create_db(
        gtf, dbfn=":memory:", force=True, keep_order=True, merge_strategy="merge",disable_infer_genes=True,disable_infer_transcripts=True)

    rows = []

    for gene in gene_list:
        try:
            transcripts = list(db.children(gene, featuretype='transcript'))
            if not transcripts:
                print (f"[WARN] No transcripts found for {gene}")
                continue
            # Collect polyA coordinate + UTR lenghths 
            polyA_coords = []
            for tx in transcripts:
                tx_id=tx.attributes.get("transcript_id", [""])[0]
                tx_name=tx.attributes.get("transcript_name",[""])[0]
                # Get CDS -> compute UTR length
                cds = list(db.children(tx,featuretype="CDS"))
                if cds:
                    if tx.strand == "+":
                        cds_end = max(c.end for c in cds)
                        utr_length = tx.end - cds_end
                        coord = tx.end
                    else:
                        cds_start = min(c.start for c in cds)
                        utr_length = cds_start - tx.start
                        coord = tx.start
                else:
                    utr_length = pd.NA
                    coord = tx.end if tx.strand =='+' else tx.start
                polyA_coords.append((coord, tx, tx_id,tx_name, utr_length))
            
            # Group transcripts by polyA ends within tolerance
            polyA_coords.sort(key=lambda x: x[0])
            groups, current_group = [], []
            for coord, tx,tx_id, tx_name, utr in polyA_coords:
                #First item starts a new cluster.
                if not current_group:
                    current_group = [(coord, tx, tx_id, tx_name,utr)]
                else:
                    if abs(coord - current_group[-1][0]) <= end_tolerance:
                        current_group.append((coord, tx, tx_id, tx_name, utr))
                    else:
                        groups.append(current_group)
            if current_group:
                groups.append(current_group)

            # Emit feature per group
            for i, grp in enumerate(groups,1):
                coords = [c for c, *_ in grp]
                txs = [tx for _, tx, *_ in grp]
                tx_ids = [tid for _, _, tid, *_ in grp]
                tx_names = [tname for _,_,_, tname,_ in grp]
                utrs = [ u for *_, u in grp if pd.notna(u)]
                
                # Default group ID
                group_id = f"{gene}_polyA{i}"
                
                # Try immunoglobulin-specific naming (short/long)
                labels = set (get_ig_label(gene, name) for name in tx_names if get_ig_label(gene,name))
                if labels: 
                    if len(labels) == 1:
                        label = labels.pop()
                        group_id= f"{gene}_{label}"
                    else:
                        print(f"[WARN] Mixed IG labels for {gene}:{tx_names} -> using default poly{i}")



                mean_end = int(sum(coords)/ len(coords))

                # Compute UTR stats
                if utrs:
                    avg_utr = float(pd.Series(utrs).mean())
                    min_utr = int(min(utrs))
                    max_utr = int(max(utrs))
                else:
                    avg_utr = min_utr = max_utr = pd.NA      
                
                # polyA_window feature
                pa_start, pa_end = mean_end - polyA_window, mean_end + polyA_window
                rows.append({
                    "gene": gene,
                    "polyA_group": group_id,
                    "feature_type": "polyA_window",
                    "feature_id": f"{group_id}_window",
                    "chrom": txs[0].seqid,
                    "start": pa_start,
                    "end": pa_end,
                    "strand": txs[0].strand,
                    "transcripts": ";".join(tx_ids),
                    "transcripts_names": ";".join(tx_names),
                    "avg_utr_length": avg_utr, 
                    "min_utr_length": min_utr,
                    "max_utr_length": max_utr,
                    "evidence_tier":"tier3_polyA"

                })
                # last_exon features for group (union of last exons)
                for tx, tx_id, tx_name in zip(txs, tx_ids, tx_names):
                    exons = list(db.children(tx, featuretype='exon'))
                    if not exons:
                        continue
                    last_exon = (max(exons, key= lambda e: e.end) 
                                if tx.strand == "+" 
                                else min(exons, key=lambda e: e.start)
                    )
                    rows.append({
                        "gene": gene,
                        "polyA_group": group_id,
                        "feature_type": "last_exon_group",
                        "feature_id": f"{group_id}_exon_{tx_id}",
                        "chrom": last_exon.seqid,
                        "start": last_exon.start,
                        "end": last_exon.end,
                        "strand": last_exon.strand,
                        "transcripts": tx_id,
                        "transcript_names": tx_name,
                        "avg_utr_length": avg_utr,
                        "min_utr_length": min_utr,
                        "max_utr_length": max_utr,
                        "evidence_tier": "tier1_group"

                    })
        except Exception as e:
            print (f"[WARN] Failed on {gene}:{e}")

    df= pd.DataFrame(rows)
    df.to_csv(out_csv, index=False)
    print(f"[IsoDecipher] Feature panel written to {out_csv}") 
    return df    

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Build isoform feature panel with polyA grouping and UTR lengths")
    ap.add_argument("--gtf", required=True, help="Ensembl GTF file")
    ap.add_argument("--genes", required=True, help="Gene list file (symbols)")
    ap.add_argument("--out", required=True, help="Output CSV")
    ap.add_argument("--polyA_window", type=int, default=300,
                    help="Window size (bp) around transcript end for polyA features")
    ap.add_argument("--end_tolerance", type=int, default=20,
                    help="Collapse transcripts with ends within this bp distance")
    args = ap.parse_args()

    build_panel_features(args.gtf, args.genes, args.out, args.polyA_window, args.end_tolerance)  
