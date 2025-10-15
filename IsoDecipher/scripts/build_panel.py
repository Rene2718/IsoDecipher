#!/usr/bin/env python3
"""
IsoDecipher: Build isoform annotation table
-------------------------------------------
From an Ensemble GTF + list of genes, extract transcript-level
annotation including last exon, UTR lengths, and IG labels.

Features:
 - per-transcript polyA (last exon end/start depending on strand)
 - UTR length relative to CDS
 - immunoglobulin short/long labels
 - concise summary (#transcripts per gene, ave UTR length)

Usage:
python IsoDecipher/sscripts/build_panel.py \
    --gtf data/Homo_sapiens.GRCh38.115.gtf \
    --genes data/gene_list.txt \
    --out results/isoform_panel.csv \
   [--debug]
"""
import gffutils
import pandas as pd
import argparse 
import time
import os 
import numpy as np

IG_WHITELIST = {"IGHM", "IGHG1", "IGHG2", "IGHG3",
                "IGHG4", "IGHA1", "IGHA2", "IGHE"}

def load_or_build_db(gtf):
    """
    Load a gffutils DB, rebuilding if the GTF file is newer or if DB doesn't exist.
    """
    db_path = gtf + ".db"   # e.g. Homo_sapiens.GRCh38.115.gtf.db
    #Check if DB exists and is up-to-date
    if os.path.exists(db_path):
        gtf_mtime = os.path.getmtime(gtf)
        db_mtime = os.path.getmtime(db_path)
        if gtf_mtime > db_mtime:
            print(f"[IsoDecipher] GTF {gtf} is newer than DB {db_path}. Rebuiding...")
            build = True
        else:
            print(f"[IsoDecipher] Re-using existing DB at {db_path}")
            build = False
    else:
        print(f"[IsoDecipher] No DB found at {db_path}. Building new one…")
        build = True
    if build:
        start = time.time()
        gffutils.create_db(
            gtf,
            dbfn=db_path,
            force=True,
            keep_order=True,
            merge_strategy="merge",
            disable_infer_genes=True,
            disable_infer_transcripts=True
        )
        print(f"[IsoDecipher] Finished building DB in {time.time()-start:.1f}s")
    return gffutils.FeatureDB(db_path, keep_order=True)

def get_ig_label(gene, tx_name):
    """Return short/long labels for immunoglobulin transcripts."""
    if gene not in IG_WHITELIST:
        return ""
    if not tx_name:
        return ""

    if tx_name:
        if tx_name.endswith("-201"):
            return "short"
        elif tx_name.endswith("-202"):
            return "long"
        elif tx_name == "IGHG1-203":
            return 'short'
        return "short"
    
    print (f"[WARN] unknown transcript name pattern: {gene} {tx_name}")
    return ""


def build_panel(gtf, gene_list_file, out_csv, debug=False):
    # Read the list 
    with open(gene_list_file) as f:
        gene_list=[line.strip() for line in f if line.strip() and not line.startswith("#")]
    gene_list = sorted(set(gene_list))

    print(f"[IsoDecipher] Loading GTF {gtf} into gffutils DB...")
  
    # Build or load persistent gffutils DB
    
    db = load_or_build_db(gtf)

    # Map gene symbol -> gene_id
    name2id = {
        g.attributes.get("gene_name",[""])[0]: g.id 
        for g in db.features_of_type("gene")
        }
    rows = []

    print(f"[IsoDecipher] Processing {len(gene_list)} genes")
    for gene in gene_list:
        gene_id = name2id.get(gene)
        if gene_id is None:
            print(f"[WARN] No transcripts for {gene}")
            continue
        transcripts = list(db.children(gene_id, featuretype='transcript'))
        if not transcripts:
            print(f"[WARN] No transcripts for {gene}")
            continue
        for tx in transcripts:
            tx_id = tx.attributes.get("transcript_id",[""])[0]
            tx_name = tx.attributes.get("transcript_name",[""])[0]
            exons = list(db.children(tx, featuretype='exon'))

            if exons:
                if tx.strand == "+":
                    last_exon = max(exons, key=lambda e: e.end)
                else:
                    last_exon = min(exons, key=lambda e: e.start)
                last_start, last_end = last_exon.start, last_exon.end
            else:
                last_start, last_end = tx.start, tx.end

            #Compute UTR length
            cds = list(db.children(tx, featuretype='CDS'))
            if cds:
                if tx.strand == "+":
                    cds_end = max(c.end for c in cds)
                    utr_len = last_end - cds_end
                else:
                    cds_start = min(c.start for c in cds)
                    utr_len = cds_start - last_start
                # safety: discard invalid/negative or 0-length UTRs
                if utr_len is None or utr_len <= 0:
                    utr_len = pd.NA
                    utr_status = "invalid_or_negative"
                else:
                    utr_status = "valid"
            else:
                utr_len = pd.NA
                utr_status = "missing_CDS"
            
            rows.append({
                "gene": gene,
                "transcript_id": tx_id,
                "transcript_name": tx_name,
                "strand": tx.strand,
                "last_exon_start": last_start,
                "last_exon_end": last_end,
                "utr_length": utr_len,
                "utr_status":utr_status,
                "ig_label": get_ig_label(gene, tx_name)
            })
            if debug:
                print(f"[DEBUG] {gene} {tx_name} ({tx_id}) "
                    f"strand={tx.strand} "
                    f"last_exon=({last_start}-{last_end}) "
                    f"utr={utr_len} status={utr_status}")
    # --- Write per-transcript table ---
    df = pd.DataFrame(rows)
    df.to_csv(out_csv, index= False)
    # --- Compute per-gene summary ---
    def safe_stat(series, func):
        return func(series.dropna()) if series.dropna().size >0 else np.nan
    

    
    summary= (
        df.groupby("gene")
        .agg(
            num_tx=("transcript_id","count"),
            avg_utr=("utr_length", lambda x: safe_stat(x, np.mean)),
            min_utr=("utr_length", lambda x: safe_stat(x, np.min)),
            max_utr=("utr_length", lambda x: safe_stat(x, np.max)),
            utr_range=("utr_length", lambda x: safe_stat(x, lambda v: np.max(v)-np.min(v))),
            utr_std=('utr_length', lambda x: safe_stat(x, np.std)),
        )
        .reset_index()
    )
    summary = summary.round({
        "avg_utr": 1,
        "min_utr": 0,
        "max_utr": 0,
        "utr_range": 0,
        "utr_std": 1
    })


    summary_path = out_csv.replace(".csv", "_summary.csv")
    summary.to_csv(summary_path, index=False)

    for _, row in summary.iterrows():
        if pd.notna(row["utr_range"]) and row["utr_range"] > 1000:
            print(f"⚠️  {row['gene']}: large UTR diversity ({row['utr_range']:.0f} bp range)") 
    print("\n[SUMMARY] Transcript stat per gene:")
    print(summary.to_string(index=False))
    print(f"[IsoDecipher] Annotation table written to {out_csv}")
    print(f"[IsoDecipher] Summary written to {summary_path}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Build isoform annotation panel")
    ap.add_argument("--gtf", required=True,help="Ensembl GTF file")
    ap.add_argument("--genes", required=True, help="Gene list file")
    ap.add_argument("--out", required=True, help="Out CSV")
    ap.add_argument("--debug", action="store_true", help="Print per-transcript debug info")
    args= ap.parse_args()

    build_panel(args.gtf, args.genes, args.out)

