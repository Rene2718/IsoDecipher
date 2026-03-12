#!/usr/bin/env python3
"""
Usage:
python IsoDecipher/scripts/integrate_samples.py \
    --exp_list data/samples.txt \
    --data_dir data \
    --iso_dir results \
    --out results/master_mosaic_combined.h5ad
"""
import scanpy as sc
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
import argparse
import gc
import os

def main():
    parser = argparse.ArgumentParser(description="Systematically integrate Genes, Isoforms, and CITE-seq")
    parser.add_argument("--exp_list", required=True, help="Text file with one experiment ID per line")
    parser.add_argument("--data_dir", default="data", help="Directory containing experiment folders")
    parser.add_argument("--iso_dir", default="results", help="Directory containing isoform CSVs")
    parser.add_argument("--out", default="results/master_mosaic_combined.h5ad", help="Path for output H5AD")
    args = parser.parse_args()

    with open(args.exp_list, 'r') as f:
        exps = [line.strip() for line in f if line.strip()]

    adatas = []

    for exp in exps:
        print(f"--- Processing {exp} ---")
        
        # 1. Load H5
        h5_path = os.path.join(args.data_dir, exp, "filtered_feature_bc_matrix.h5")
        if not os.path.exists(h5_path):
            continue
            
        adata_h5 = sc.read_10x_h5(h5_path, gex_only=False)
        adata_h5.var_names_make_unique()
        
        # Split features
        genes = adata_h5[:, adata_h5.var['feature_types'] == 'Gene Expression'].copy()
        proteins = adata_h5[:, adata_h5.var['feature_types'] == 'Antibody Capture'].copy()
        
        # 2. Load Isoform Counts
        iso_path = os.path.join(args.iso_dir, f"{exp}_isoform_counts.csv")
        if os.path.exists(iso_path):
            df_iso = pd.read_csv(iso_path, index_col=0)
            
            # Create AnnData with metadata preserved
            adata_iso = ad.AnnData(X=csr_matrix(df_iso.values))
            adata_iso.obs_names = df_iso.index
            adata_iso.var_names = df_iso.columns
            adata_iso.var['feature_types'] = 'Isoform'
            
            # Ensure barcodes match (stripping -1 if one has it and the other doesn't)
            genes.obs_names = genes.obs_names.str.replace("-1", "")
            adata_iso.obs_names = adata_iso.obs_names.str.replace("-1", "")
            
            common_cells = genes.obs_names.intersection(adata_iso.obs_names)
            genes = genes[common_cells].copy()
            adata_iso = adata_iso[common_cells].copy()
            
            # Combine GEX + Iso
            combined_sample = ad.concat([genes, adata_iso], axis=1, merge="first")
            
            # Add ADT if exists
            if proteins.n_vars > 0:
                proteins.obs_names = proteins.obs_names.str.replace("-1", "")
                common_pro = common_cells.intersection(proteins.obs_names)
                proteins = proteins[common_pro].copy()
                proteins.var_names = [f"prot_{n}" for n in proteins.var_names]
                proteins.var['feature_types'] = 'ADT'
                combined_sample = ad.concat([combined_sample, proteins], axis=1, merge="first")
        else:
            combined_sample = genes
        
        # 3. Finalize sample metadata
        combined_sample.obs['batch'] = exp
        combined_sample.obs_names = [f"{exp}_{bc}" for bc in combined_sample.obs_names]
        
        if not isinstance(combined_sample.X, csr_matrix):
            combined_sample.X = csr_matrix(combined_sample.X)
            
        adatas.append(combined_sample)
        gc.collect()

    print("\nMerging all samples...")
    # join='outer' is necessary, but we manually re-tag feature_types after
    adata_final = ad.concat(adatas, join='outer', fill_value=0)
    
    # SYSTEMATIC FIX: Re-apply feature_types to the final object
    adata_final.var['feature_types'] = 'Gene Expression'
    adata_final.var.loc[adata_final.var_names.str.contains(r'_G\d+', regex=True), 'feature_types'] = 'Isoform'
    adata_final.var.loc[adata_final.var_names.str.startswith('prot_'), 'feature_types'] = 'ADT'

    adata_final.write(args.out)
    print(f"✅ Saved to: {args.out} ({adata_final.n_obs} cells)")

if __name__ == "__main__":
    main()