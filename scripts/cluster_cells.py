#!/usr/bin/env python
"""Cluster QC'd scRNA cells using PCA + Leiden/Louvain and emit a static UMAP plot.

Produces `outputs/umap_clusters.png` after clustering completes.
"""

import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc


def main():
    if len(sys.argv) < 3:
        print("Usage: cluster_cells.py <qc_h5ad> <out_csv>", file=sys.stderr)
        sys.exit(2)
    h5ad, out_csv = sys.argv[1], sys.argv[2]

    h5ad_path = Path(h5ad)
    if not h5ad_path.is_file():
        print(f"ERROR: QC h5ad not found: {h5ad_path}", file=sys.stderr)
        sys.exit(1)

    adata = sc.read_h5ad(h5ad_path)
    if adata.n_obs == 0 or adata.n_vars == 0:
        print(
            f"ERROR: AnnData loaded from {h5ad_path} is empty "
            f"({adata.n_obs} cells × {adata.n_vars} genes).",
            file=sys.stderr,
        )
        sys.exit(1)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

    try:
        sc.tl.leiden(adata, resolution=0.5)
    except ImportError as e:
        print(
            "ERROR: Leiden clustering requires `python-igraph` and `leidenalg` "
            "to be installed in the environment.",
            file=sys.stderr,
        )
        print(f"Underlying error: {e}", file=sys.stderr)
        sys.exit(1)

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
        fig = sc.pl.umap(
            adata, color="leiden", legend_loc="on data", return_fig=True, show=False
        )
        if isinstance(fig, (list, tuple)) and len(fig) > 0:
            fig = fig[0]
        plt.figure(fig.number)
        plt.savefig("outputs/umap_clusters.png", dpi=200, bbox_inches="tight")
        plt.close(fig)
        print("Saved UMAP plot -> outputs/umap_clusters.png")
    except Exception as e:
        # Don't crash the script if plotting fails; report the error for reviewers
        print(f"UMAP plotting failed: {e}", file=sys.stderr)


if __name__ == "__main__":
    main()
