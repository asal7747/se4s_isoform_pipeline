#!/usr/bin/env python
"""Cluster QC'd scRNA cells using PCA + Leiden/Louvain."""
import sys, scanpy as sc, pandas as pd

def main():
    if len(sys.argv) < 3:
        print("Usage: cluster_cells.py <qc_h5ad> <out_csv>"); sys.exit(2)
    h5ad, out_csv = sys.argv[1], sys.argv[2]
    
    adata = sc.read_h5ad(h5ad)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.leiden(adata, resolution=0.5)
    
    clusters = pd.DataFrame({
        "cell": adata.obs_names,
        "cluster": adata.obs["leiden"].astype(str)
    })
    clusters.to_csv(out_csv, index=False)
    print(f"Clustered {len(clusters)} cells into {len(set(clusters['cluster']))} clusters -> {out_csv}")

if __name__ == "__main__":
    main()
