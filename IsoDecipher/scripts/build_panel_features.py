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
    --gtf Homo_sapiens.GRCh38.115.gtf \
    --genes data/gene_list.txt \
    --out result/panel_features.csv \
    [--custom_params custom.tsv] \
    [--no-skip_singleton] \
    [--no-skip_collapsed]

By default, single-transcript genes and collapsed groups are skipped.
Use --no-skip_singleton and/or --no-skip_collapsed to include them.

"""

import argparse
import gffutils
import pandas as pd
import time
import os
import matplotlib.pyplot as plt

def load_custom_parameters(file):
    """
    Load user-defined parameters from TSV/CSV
    Return dict:{gene: (polyA_window, end_tolerance)}
    """
    if not file:
        return {}
    df = pd.read_csv(file, sep=None, engine='python')
    if not {"gene", "polyA_window", "end_tolerance"} <= set(df.columns):
        raise ValueError("Custom parameter must contain columns: gene, polyA_window, end_tolerance")

    param_dict={
        row["gene"]:(int(row["polyA_window"]),int(row["end_tolerance"]))
        for _, row in df.iterrows()
    }
    print(f"[CUSTOM] Load {len(param_dict)} gene-specific parameters")
    return param_dict

def get_dynamic_parameters(gene, num_transcripts, custom_params=None, base_window=200, base_tolerance=40, strategy="balanced"):
    """
    Dynamically adjust polyA_window and end_tolerance per gene. 
    Priority:
    1. User custom list (highest priority)
    2. Smart adjustment based on user's preference
    Args:
        strategy: 
            "precise"   
            "balanced"  
            "sensitive" 
    """
    if custom_params and gene in custom_params:
        return custom_params[gene]

    STRATEGIES = {
        "precise":{
            "complex_multiplier": 1.0,
            "standard_multiplier": 0.5 ,
            "description": "High precision - avoids false merging, best for known functional isoforms"
        },
        "balanced": {  
            "complex_multiplier": 1.5,
            "standard_multiplier": 1.0 ,
            "description": "Balanced approach - good for most use cases"
        },
        "sensitive": { 
            "complex_multiplier": 2.0,
            "standard_multiplier": 1.5 ,
            "description": "High sensitivity - detects more signals, good for exploratory analysis"
        }
    }
    if strategy not in STRATEGIES:
        raise ValueError(f"Unknown strategy: {strategy}.Choose from {list(STRATEGIES.keys())}")
    params = STRATEGIES[strategy]   

    if num_transcripts >15: #complex gene, many isoform
        return base_window, int(base_tolerance*params["complex_multiplier"])

    else:
        return base_window, int(base_tolerance*params["standard_multiplier"])
    

IG_WHITELIST = {"IGHM", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHE"}

def get_ig_label(gene, tx_name):
    """Return short/long labels for immunoglobulin transcripts."""
    if gene not in IG_WHITELIST:
        return None
    if tx_name:
        if tx_name.endswith("-201"):
            return "short"
        elif tx_name.endswith("-202"):
            return "long"
        elif tx_name == "IGHG1-203":
            return "short"
        # fallback on transcript_id if no usable name

    print(f"[WARN] Unknown transcript name pattern: {gene} {tx_name}")
    return None

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

        


def build_panel_features(gtf, gene_list_file, out_csv, polyA_window=200, end_tolerance=40, custom_params=None, strategy="balanced",skip_singleton=True, skip_collapsed=True):
   
    #read gene list
    with open(gene_list_file) as f:
        gene_list= [line.strip() for line in f if line.strip() and not line.startswith("#")]
    gene_list = sorted(set(gene_list))

    print(f"[IsoDecipher] Loading GTF {gtf} into gffutils DB (this may take several minutes)…")
    
    db = load_or_build_db(gtf)
    
    name2id = {
        g.attributes.get("gene_name",[""])[0]:g.id
        for g in db.features_of_type("gene")
    }
    
    rows = []
    summary_records = []
    
    print(f"[IsoDecipher] Processing {len(gene_list)} genes from {gene_list_file}")

    for gene in gene_list:
        gene_id = name2id.get(gene)
        if gene_id is None:
            print(f"[WARN] Gene {gene} not found in GTF")
            summary_records.append({"gene": gene, "num_groups": 0, "num_transcripts":0, "status": "not_found"})
            continue
        try:
            transcripts = list(db.children(gene_id, featuretype='transcript'))
            num_transcripts = len(transcripts)
            if not transcripts:
                print (f"[WARN] No transcripts found for {gene}")
                summary_records.append({"gene": gene, "num_groups": 0, "num_transcripts": 0, "status": "no_transcripts"})
                continue
            if skip_singleton and num_transcripts <= 1:
                print(f"[SKIP] {gene}: only one transcript")
                summary_records.append({"gene":gene, "num_groups": 0, "num_transcripts": num_transcripts, "status":"single_transcript"})
                continue


            # set dynamic parmas for this gene
            gene_specific_window, gene_specific_tolerance = get_dynamic_parameters(
               gene= gene,
               num_transcripts=len(transcripts),
               custom_params=custom_params,
               base_window=polyA_window, 
               base_tolerance=end_tolerance,
               strategy=strategy
            )
            print (f"[PARAMS] {gene}: window={gene_specific_window}, tolerance={gene_specific_tolerance}")
            
            # Collect polyA coordinate + UTR lengths 
            polyA_coords = []
            for tx in transcripts:
                tx_id=tx.attributes.get("transcript_id", [""])[0]
                tx_name=tx.attributes.get("transcript_name",[""])[0]
                # --- Get last exon boundary ---
                exons = list(db.children(tx,featuretype='exon'))
                if exons: 
                    if tx.strand == "+":
                        last_exon = max(exons, key=lambda e: e.end)
                        coord = last_exon.end
                    else:
                        last_exon = min(exons, key=lambda e: e.start)
                        coord = last_exon.start
                else:
                    coord = tx.end if tx.strand == "+" else tx.start


                # Get CDS -> compute UTR length
                cds = list(db.children(tx,featuretype="CDS"))
                if cds:
                    if tx.strand == "+":
                        cds_end = max(c.end for c in cds)
                        utr_length = coord - cds_end
                        cds_pos = cds_end
                    else:
                        cds_start = min(c.start for c in cds)
                        utr_length = cds_start - coord
                        cds_pos = cds_start
                else:
                    utr_length = pd.NA
                    cds_pos = "NA"

                polyA_coords.append((coord, tx, tx_id,tx_name, utr_length))
            
                print(f"[DEBUG] {gene} {tx_name} strand={tx.strand} polyA={coord} cds_pos={cds_pos} utr={utr_length}")

            # Group transcripts by polyA ends within tolerance
            polyA_coords.sort(key=lambda x: x[0])
            groups, current_group = [], []
            for coord, tx,tx_id, tx_name, utr in polyA_coords:
                #First item starts a new cluster.
                if not current_group:
                    current_group = [(coord, tx, tx_id, tx_name,utr)]
                else:
                    if abs(coord - current_group[-1][0]) <= gene_specific_tolerance:
                        current_group.append((coord, tx, tx_id, tx_name, utr))
                    else:
                        groups.append(current_group)
                        current_group = [(coord, tx, tx_id, tx_name, utr)] 

            if current_group:
                groups.append(current_group)

            # 🔍 Debug: show all transcripts and their chosen coords
            print(f"[DEBUG] {gene} transcripts coordinates: {[(tx_name, coord) for coord, tx, tx_id, tx_name, utr in polyA_coords]}")
            print(f"[DEBUG] {gene} formed {len(groups)} groups")

            
            
            # Skip non-informative
            num_groups= len(groups)
            if skip_collapsed and num_groups <= 1:
                print(f"[SKIP] {gene}: collapsed into one group")
                summary_records.append({"gene": gene, "num_groups": num_groups, "num_transcripts": num_transcripts, "status": "collapsed"})
                continue
            
            summary_records.append({"gene": gene, "num_groups": num_groups, "num_transcripts": num_transcripts, "status": "informative"})
            

            # Emit feature per group
            for i, grp in enumerate(groups,1):
                coords = [c for c, *_ in grp]
                txs = [tx for _, tx, *_ in grp]
                tx_ids = [tid for _, _, tid, *_ in grp]
                tx_names = [tname for _,_,_, tname,_ in grp]
                utrs = [ u for *_, u in grp if pd.notna(u)]
                
                # Default group ID
                group_id = f"{gene}_polyA{i}"
                mean_end = int(sum(coords)/ len(coords))

                # Compute UTR stats
                if utrs:
                    avg_utr = float(pd.Series(utrs).mean())
                    min_utr = int(min(utrs))
                    max_utr = int(max(utrs))
                else:
                    avg_utr = min_utr = max_utr = pd.NA      
                
                # polyA_window feature
                pa_start, pa_end = mean_end - gene_specific_window, mean_end + gene_specific_window
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
                    "transcript_names": ";".join(tx_names),
                    "avg_utr_length": avg_utr, 
                    "min_utr_length": min_utr,
                    "max_utr_length": max_utr,
                    "evidence_tier":"tier3_polyA"

                })

                # last_exon features for group (union of last exons)
                for tx, tx_id, tx_name in zip(txs, tx_ids, tx_names):
                    exons = list(db.children(tx, featuretype='exon'))
                    
                    if exons:
                        if tx.strand == '+':
                            last_exon = max(exons,key=lambda e:e.end)
                        else:
                            last_exon = min(exons,key=lambda e:e.start)
                        chrom = last_exon.seqid
                        start, end, strand = last_exon.start, last_exon.end, last_exon.strand
                    else:
                        chrom = tx.seqid
                        start, end, strand = tx.start, tx.end, tx.strand

                    rows.append({
                        "gene": gene,
                        "polyA_group": group_id,
                        "feature_type": "last_exon_group",
                        "feature_id": f"{group_id}_exon_{tx_id}",
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "transcripts": tx_id,
                        "transcript_names": tx_name,
                        "avg_utr_length": avg_utr,
                        "min_utr_length": min_utr,
                        "max_utr_length": max_utr,
                        "evidence_tier": "tier1_group"

                    })
        except Exception as e:
            print (f"[WARN] Failed on {gene}:{e}")
            summary_records.append({"gene": gene, "num_groups": 0, "num_transcripts": 0, "status": f"error:{e}"})


    df= pd.DataFrame(rows)

    # Relabel IG genes to short/long
    def relabel_ig_groups(df):
        out_rows = []
        for _, row in df.iterrows():
            gene = row["gene"]
            if gene in IG_WHITELIST:
                tx_names = row["transcript_names"].split(";")
                labels = []
                for tx_name in tx_names:
                    label = get_ig_label(gene, tx_name)
                    if label is not None:
                        labels.append(label)
                if labels:
                    if len(set(labels)) == 1:
                        label = labels[0]
                        row["polyA_group"] = f"{gene}_{label}"
                    else:
                        print(f"[WARN] Mixed IG labels for {gene}:{tx_names}")
            out_rows.append(row)
        return pd.DataFrame(out_rows)

    df = relabel_ig_groups(df)
    df.to_csv(out_csv, index=False)

    # Save summary
    summary_df = pd.DataFrame(summary_records)
    summary_path = out_csv.replace(".csv","_summary.csv")
    summary_df.to_csv(summary_path, index=False)


    print("\n[SUMMARY] PolyA group per gene:")
    print(summary_df.to_string(index=False))
    print(f"[IsoDecipher] Feature panel written to {out_csv}") 
    print(f"[IsoDecipher] Summary table written to {summary_path}")


  

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Build isoform feature panel with polyA grouping and UTR lengths")
    ap.add_argument("--gtf", required=True, help="Ensembl GTF file")
    ap.add_argument("--genes", required=True, help="Gene list file (symbols)")
    ap.add_argument("--out", required=True, help="Output CSV")
    ap.add_argument("--polyA_window", type=int, default=200,
                    help="Window size (bp) around transcript end for polyA features")
    ap.add_argument("--end_tolerance", type=int, default=40,
                    help="Collapse transcripts with ends within this bp distance")
    ap.add_argument("--strategy", choices=["precise","balanced","sensitive"], default='balanced', help="Analysis strategy: precise, balanced, or sensitive")
    ap.add_argument("--custom_params",help="Optional TSV/CSV with custom per-gene parameters (gene, polyA_window, end_tolerance)",default=None)
    # Skip flags with both enable/disable options
    ap.add_argument("--skip_singleton", dest="skip_singleton", action="store_true",
                    help="Skip genes with only one transcript (default: on)")
    ap.add_argument("--no-skip_singleton", dest="skip_singleton", action="store_false",
                    help="Include single-transcript genes")
    ap.set_defaults(skip_singleton=True)

    ap.add_argument("--skip_collapsed", dest="skip_collapsed", action="store_true",
                    help="Skip genes collapsed into one group (default: on)")
    ap.add_argument("--no-skip_collapsed", dest="skip_collapsed", action="store_false",
                    help="Include collapsed genes")
    ap.set_defaults(skip_collapsed=True)

    args = ap.parse_args()
    custom_dict = load_custom_parameters(args.custom_params)
    # Ensure results/ exists and adjust output path
    out_path = args.out
    if not ("/" in out_path or "\\" in out_path):  # user gave just filename
        os.makedirs("results", exist_ok=True)
        out_path = os.path.join("results", out_path)
    else:
        os.makedirs(os.path.dirname(out_path), exist_ok=True)


    build_panel_features(
        args.gtf, 
        args.genes, 
        out_path, 
        args.polyA_window, 
        args.end_tolerance, 
        custom_params=custom_dict,
        strategy=args.strategy,
        skip_singleton=args.skip_singleton,
        skip_collapsed=args.skip_collapsed
        )  
