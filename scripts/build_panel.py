#!/usr/bin/env python3
"""
IsoDecipher: Build an isoform panel fron GTF annotation
This script parses a GTF file and a list of genes,
then outputs a CSV panel with all isoforms, their transcript IDs,
and estimated 3′UTR lengths.

Usage:
    python build_panel.py --gtf data.gtf --genes gene_list.txt --out panel.csv
"""
import gffutils
import pandas as pd

def build_isoform_panel(gtf, gene_list_file, out_csv):
    """
    Parse a GTF file for transcripts of given genes
    and calculate 3′UTR lengths.

    Args:
        gtf (str): Path to GTF annotation file
        gene_list_file (str): Path to text file with gene symbols
        out_csv (str): Path to save output CSV panel

    Returns:
        pd.DataFrame: Table of isoforms and their properties
    """
    #read gene list
    with open(gene_list_file) as f:
        gene_list=[line.strip() for line in f if line.strip()]

    #load GTF into gffutils database
    db = gffutils.create_db(gtf, dbfn=":memory:", force=True, keep_order=True, merge_strategy="merge")

    results = []

    for gene in gene_list:
        try:
            #loop through transcripts of this gene
            for tx in db.children(gene, featuretype='transcript'):
                tx_id = tx.attributes.get("transcript_id",[""])[0]
                tx_name = tx.attributes.get("transcript_name",[""])[0]

                #get coding end -> calculat 3'UTR length
                cds=list(db.children(tx, featuretype='CDS'))
                if cds:
                    if tx.strand =="+":
                        cds_end = max(c.end for c in cds)
                        utr_length = tx.end -cds_end
                    else:
                        cds_start = min(c.start for c in cds)
                        utr_length = cds_start - tx.start
                else:
                    utr_length = "NA"
                results.append({
                    "gene": gene,
                    "transcript_id": tx_id,
                    "transcript_name": tx_name,
                    "utr_length":utr_length,
                    "isoform_type":"APA/ALE candidate",
                    "detectability": "3' scRNA-seq (likely)"

                })
        except Exception as e:
            print(f"Warning: gene {gene} not found in GTF({e})")
    #save results 
    df = pd.DataFrame(results)
    df.to_csv(out_csv, index=False)
    print(f"[IsoDecipher] Wrote panel to {out_csv}")
    return df
    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Build isoform panel from GTF")
    parser.add_argument("--gtf", required=True, help="GENCODE/Ensembl GTF file")
    parser.add_argument("--genes", required=True, help="Gene list (txt)")
    parser.add_argument("--out", required=True, help="Output CSV panel")
    args = parser.parse_args()

    build_isoform_panel(args.gtf, args.genes, args.out)  

        
