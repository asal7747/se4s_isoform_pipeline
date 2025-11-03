#!/usr/bin/env python
"""Cluster QC'd scRNA cells using PCA + Leiden/Louvain and emit a static UMAP plot.

Produces `outputs/umap_clusters.png` after clustering completes.
"""

import os
import sys

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc


def main():
    if len(sys.argv) < 3:
        print("Usage: cluster_cells.py <qc_h5ad> <out_csv>")
        sys.exit(2)
    h5ad, out_csv = sys.argv[1], sys.argv[2]

    adata = sc.read_h5ad(h5ad)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.leiden(adata, resolution=0.5)

    clusters = pd.DataFrame(
        {"cell": adata.obs_names, "cluster": adata.obs["leiden"].astype(str)}
    )
    clusters.to_csv(out_csv, index=False)
    print(
        f"Clustered {len(clusters)} cells into {len(set(clusters['cluster']))} clusters -> {out_csv}"
    )

    # Create outputs directory if it doesn't exist and save UMAP plot
    os.makedirs("outputs", exist_ok=True)
    try:
        sc.tl.umap(adata)
        # return_fig=True gives a matplotlib Figure we can save; show=False avoids GUI
        fig = sc.pl.umap(
            adata, color="leiden", legend_loc="on data", return_fig=True, show=False
        )
        # sc.pl.umap may return a single Figure or a list; handle both
        if isinstance(fig, (list, tuple)) and len(fig) > 0:
            fig = fig[0]
        plt.figure(fig.number)
        plt.savefig("outputs/umap_clusters.png", dpi=200, bbox_inches="tight")
        plt.close(fig)
        print("Saved UMAP plot -> outputs/umap_clusters.png")
    except Exception as e:
        # Don't crash the script if plotting fails; report the error for reviewers
        print(f"UMAP plotting failed: {e}")


if __name__ == "__main__":
    main()
