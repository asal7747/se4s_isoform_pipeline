#!/usr/bin/env python
"""Cluster QC'd scRNA cells using PCA + Leiden/Louvain and generate optional dimensionality reduction plots.

Default behavior produces `outputs/umap_clusters.png` after clustering completes.

Can be run as:
    python cluster_cells.py <qc_h5ad> <out_csv> [--pca] [--tsne] [--umap] [--color-by VARIABLE]
"""

import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

import argparse


def normalize_and_transform_adata(
    adata: sc.AnnData,
    target_sum: float = 1e4,
    min_mean: float = 0.0125,
    max_mean: float = 3,
    min_disp: float = 0.5,
    output_dir: str = "outputs/anndata",
    show_plots: bool = True
) -> sc.AnnData:
    """Normalize, log transform, and identify highly variable genes in an AnnData object.
    
    Performs standard scanpy normalization pipeline: total count normalization,
    log1p transformation, and identification of highly variable genes. The raw
    data is stored in adata.raw before any filtering.
    
    Args:
        adata: The AnnData object to process.
        target_sum: Target sum for normalization.
        min_mean: Minimum mean for highly variable genes.
        max_mean: Maximum mean for highly variable genes.
        min_disp: Minimum dispersion for highly variable genes.
        output_path: Path to save the processed AnnData object.
        show_plots: Whether to show plots.
    
    Returns:
        The normalized and transformed AnnData object.
    """
    # Normalize (i.e. library-size correct) the data matrix
    # so that counts become comparable among cells
    sc.pp.normalize_total(adata, target_sum=target_sum)
    print(adata.X[:20, :20])

    # Log transform the data
    sc.pp.log1p(adata)
    print(adata.X[:20, :20])

    # Identify highly variable genes
    sc.pp.highly_variable_genes(adata, min_mean=min_mean, max_mean=max_mean, min_disp=min_disp)
    
    if show_plots:
        sc.pl.highly_variable_genes(adata)

    # This simply freezes the state of the AnnData object
    adata.raw = adata
    print(adata)

    # Save normalized and transformed anndata object
    if output_dir:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # Use dataset_name from uns if available, otherwise use a default name
        dataset_name = adata.uns.get('dataset_name', 'combined_dataset')
        norm_h5ad_path = output_path / f"{dataset_name}_norm_log_transformed.h5ad"
        adata.write(norm_h5ad_path)
        print(f"Normalized and transformed AnnData saved to: {norm_h5ad_path}")
    
    return adata


def run_pca_analysis(adata, svd_solver='arpack', color_by='SampleType', n_pcs=50, 
                     annotate_var_explained=True, show_plots=True):
    """
    Run PCA analysis and create visualization plots.
    
    Parameters:
    -----------
    adata : AnnData
        The AnnData object to perform PCA on
    svd_solver : str
        SVD solver to use (default: 'arpack')
    color_by : str
        Variable to color the PCA plot by (default: 'SampleType')
    n_pcs : int
        Number of principal components to show in variance plot (default: 50)
    annotate_var_explained : bool
        Whether to annotate variance explained (default: True)
    show_plots : bool
        Whether to show plots (default: True)
    
    Returns:
    --------
    AnnData
        The AnnData object with PCA results added
    """
    # run PCA
    sc.tl.pca(adata, svd_solver=svd_solver)
    
    if show_plots:
        # plot first 2 principal components
        sc.pl.pca(adata, annotate_var_explained=annotate_var_explained, color=color_by)
        
        # inspect the contribution of single PCs to the total variance in the data
        sc.pl.pca_variance_ratio(adata, log=True, n_pcs=n_pcs)
    
    return adata

def compute_tsne(
    adata: sc.AnnData,
    perplexity: float = 20,
    learning_rate: float = 1000,
    random_state: int = 42,
    n_pcs: int = 40,
    color_by: list[str] | str | None = None,
    show_plots: bool = True
) -> sc.AnnData:
    """Compute t-SNE dimensionality reduction and generate plots.
    
    Args:
        adata: The AnnData object to compute t-SNE on.
        perplexity: The perplexity parameter for t-SNE.
        learning_rate: The learning rate for t-SNE optimization.
        random_state: Random seed for reproducibility.
        n_pcs: Number of principal components to use.
        color_by: Variable(s) to color the plot by (e.g., 'SampleType', 'CellType').
        show_plots: Whether to display plots.
    
    Returns:
        The AnnData object with t-SNE coordinates added to obsm.
    """
    # Compute t-SNE
    sc.tl.tsne(adata, perplexity=perplexity, learning_rate=learning_rate,
               random_state=random_state, n_pcs=n_pcs)
    
    # Generate plots if requested
    if show_plots and color_by:
        if isinstance(color_by, str):
            color_by = [color_by]
        for color in color_by:
            sc.pl.tsne(adata, color=color)
    
    return adata


def compute_umap(
    adata: sc.AnnData,
    n_neighbors: int = 10,
    n_pcs: int = 50,
    color_by: list[str] | str | None = None,
    show_plots: bool = True
) -> sc.AnnData:
    """Compute UMAP dimensionality reduction and generate plots.
    
    Args:
        adata: The AnnData object to compute UMAP on.
        n_neighbors: Number of neighbors for UMAP computation.
        n_pcs: Number of principal components to use.
        color_by: Variable(s) to color the plot by (e.g., 'SampleType', 'CellType').
        show_plots: Whether to display plots.
    
    Returns:
        The AnnData object with UMAP coordinates added to obsm.
    """
    # Compute neighbors and UMAP
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.umap(adata)
    
    # Generate plots if requested
    if show_plots and color_by:
        if isinstance(color_by, str):
            color_by = [color_by]
        for color in color_by:
            sc.pl.umap(adata, color=color)
    
    return adata

def main():
    # Parse arguments using argparse
    parser = argparse.ArgumentParser(
        description='Cluster QC\'d scRNA cells using PCA + Leiden/Louvain and generate optional dimensionality reduction plots.',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Positional arguments
    parser.add_argument('qc_h5ad', type=str, help='Path to QC\'d h5ad file')
    parser.add_argument('out_csv', type=str, help='Path to output CSV file for cluster assignments')
    
    # Optional flags
    parser.add_argument('--pca', action='store_true', help='Generate PCA plots')
    parser.add_argument('--tsne', action='store_true', help='Generate t-SNE plots')
    parser.add_argument('--umap', action='store_true', help='Generate UMAP plots')
    parser.add_argument('--color-by', type=str, default='leiden', 
                        help='Variable to color plots by (default: leiden)')
    
    args = parser.parse_args()

    h5ad_path = Path(args.qc_h5ad)
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

    # Normalize and transform
    adata = normalize_and_transform_adata(adata, show_plots=False)
    
    # Run PCA
    if args.pca:
        adata = run_pca_analysis(adata, color_by=args.color_by, show_plots=True)
    else:
        # Still need PCA for downstream analysis
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
    clusters.to_csv(args.out_csv, index=False)
    print(
        f"Clustered {len(clusters)} cells into {len(set(clusters['cluster']))} clusters -> {args.out_csv}"
    )

    # Create outputs directory if it doesn't exist
    os.makedirs("outputs", exist_ok=True)
    
    # Generate optional dimensionality reduction plots
    if args.tsne:
        print("Generating t-SNE plot...")
        try:
            compute_tsne(adata, color_by=args.color_by, show_plots=True)
        except Exception as e:
            print(f"t-SNE plotting failed: {e}", file=sys.stderr)
    
    if args.umap:
        print("Generating UMAP plot...")
        try:
            compute_umap(adata, color_by=args.color_by, show_plots=True)
        except Exception as e:
            print(f"UMAP plotting failed: {e}", file=sys.stderr)
    
    # Always save a static UMAP plot for the default pipeline output
    try:
        if 'X_umap' not in adata.obsm:
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
