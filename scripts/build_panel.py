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
    df = gffutils.create_db(gtf,dbfn=":memory:", force=True. keep_order=True, merge_strategy="merge")
    
    result = []

    for gene in gene_list:
        try:
            #loop through transcripts of this gene
            for tx in db.children(gene, featuretype='transcript'):
                tx_id = tx.attribute.get("transcript_id",[""])[0]
                tc_name = tx.attribute.get("transcript_name",[""])[0]

                #get coding end -> calculat 3'UTR length
        
