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
    show_plots: bool = True,
    skip_if_normalized: bool = True
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
        skip_if_normalized: Skip normalization if data appears already normalized.
    
    Returns:
        The normalized and transformed AnnData object.
    """
    # Check if data is already normalized (log-transformed values are typically < 10)
    # and has highly variable genes marked
    is_normalized = False
    if skip_if_normalized and adata.raw is not None:
        print("Data appears already normalized (adata.raw exists). Skipping normalization.")
        is_normalized = True
    elif skip_if_normalized and 'highly_variable' in adata.var.columns:
        print("Data appears already processed (highly_variable genes marked). Skipping normalization.")
        is_normalized = True
    
    if not is_normalized:
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
    else:
        print(f"Using existing normalization. Data shape: {adata.shape}")

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
                     annotate_var_explained=True, show_plots=True, save_figures=False,
                     output_dir='outputs/clustering', legend_loc=None):
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
    save_figures : bool
        Whether to save figures to disk (default: False)
    output_dir : str
        Directory to save figures (default: 'outputs/clustering')
    legend_loc : str, optional
        Location of the legend. Use 'on data' to show labels on the plot itself
    
    Returns:
    --------
    AnnData
        The AnnData object with PCA results added
    """
    # run PCA if not already present
    if 'X_pca' not in adata.obsm:
        sc.tl.pca(adata, svd_solver=svd_solver)
    else:
        print("PCA already computed, using existing results")
    
    if show_plots or save_figures:
        if save_figures:
            os.makedirs(output_dir, exist_ok=True)
        
        # Normalize column name to handle case variations
        actual_color_by = color_by
        if color_by not in adata.obs.columns:
            for col in adata.obs.columns:
                if col.lower() == color_by.lower():
                    actual_color_by = col
                    break
        
        # plot first 2 principal components
        sc.pl.pca(adata, annotate_var_explained=annotate_var_explained, 
                 color=actual_color_by, show=False if save_figures else True)
        if save_figures:
            plt.savefig(f"{output_dir}/pca_{color_by}.png", dpi=200, bbox_inches="tight")
            plt.close()
            print(f"Saved PCA plot -> {output_dir}/pca_{color_by}.png")
        
        # inspect the contribution of single PCs to the total variance in the data
        sc.pl.pca_variance_ratio(adata, log=True, n_pcs=n_pcs, show=False if save_figures else True)
        if save_figures:
            plt.savefig(f"{output_dir}/pca_variance_ratio.png", dpi=200, bbox_inches="tight")
            plt.close()
            print(f"Saved PCA variance ratio plot -> {output_dir}/pca_variance_ratio.png")
    
    return adata

def compute_tsne(
    adata: sc.AnnData,
    perplexity: float = 20,
    learning_rate: float = 1000,
    random_state: int = 42,
    n_pcs: int = 40,
    color_by: list[str] | str | None = None,
    show_plots: bool = True,
    save_figures: bool = False,
    output_dir: str = "outputs/clustering"
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
        save_figures: Whether to save figures to disk.
        output_dir: Directory to save figures.
    
    Returns:
        The AnnData object with t-SNE coordinates added to obsm.
    """
    # Compute t-SNE if not already present
    if 'X_tsne' not in adata.obsm:
        sc.tl.tsne(adata, perplexity=perplexity, learning_rate=learning_rate,
                   random_state=random_state, n_pcs=n_pcs)
    else:
        print("t-SNE already computed, using existing results")
    
    # Generate plots if requested
    if (show_plots or save_figures) and color_by:
        if save_figures:
            os.makedirs(output_dir, exist_ok=True)
        
        if isinstance(color_by, str):
            color_by = [color_by]
        
        for color in color_by:
            # Normalize column name to handle case variations
            actual_color = color
            if color not in adata.obs.columns:
                for col in adata.obs.columns:
                    if col.lower() == color.lower():
                        actual_color = col
                        break
            
            sc.pl.tsne(adata, color=actual_color, show=False if save_figures else True)
            if save_figures:
                plt.savefig(f"{output_dir}/tsne_{color}.png", dpi=200, bbox_inches="tight")
                plt.close()
                print(f"Saved t-SNE plot -> {output_dir}/tsne_{color}.png")
    
    return adata


def compute_umap(
    adata: sc.AnnData,
    n_neighbors: int = 10,
    n_pcs: int = 50,
    color_by: list[str] | str | None = None,
    show_plots: bool = True,
    save_figures: bool = False,
    output_dir: str = "outputs/clustering",
    legend_loc: str | None = None
) -> sc.AnnData:
    """Compute UMAP dimensionality reduction and generate plots.
    
    Args:
        adata: The AnnData object to compute UMAP on.
        n_neighbors: Number of neighbors for UMAP computation.
        n_pcs: Number of principal components to use.
        color_by: Variable(s) to color the plot by (e.g., 'SampleType', 'CellType').
        show_plots: Whether to display plots.
        save_figures: Whether to save figures to disk.
        output_dir: Directory to save figures.
        legend_loc: Location of the legend. Use 'on data' to show labels on the plot.
    
    Returns:
        The AnnData object with UMAP coordinates added to obsm.
    """
    # Compute neighbors and UMAP if not already present
    if 'X_umap' not in adata.obsm:
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
        sc.tl.umap(adata)
    else:
        print("UMAP already computed, using existing results")
    
    # Generate plots if requested
    if (show_plots or save_figures) and color_by:
        if save_figures:
            os.makedirs(output_dir, exist_ok=True)
        
        if isinstance(color_by, str):
            color_by = [color_by]
        
        for color in color_by:
            # Normalize column name to handle case variations
            actual_color = color
            if color not in adata.obs.columns:
                for col in adata.obs.columns:
                    if col.lower() == color.lower():
                        actual_color = col
                        break
            
            sc.pl.umap(adata, color=actual_color, legend_loc=legend_loc, show=False if save_figures else True)
            if save_figures:
                plt.savefig(f"{output_dir}/umap_{color}.png", dpi=200, bbox_inches="tight")
                plt.close()
                print(f"Saved UMAP plot -> {output_dir}/umap_{color}.png")
    
    return adata


def normalize_column_name(adata: sc.AnnData, column_name: str) -> str | None:
    """Find the actual column name in adata.obs, handling case variations.
    
    Args:
        adata: The AnnData object to search.
        column_name: The column name to find (case-insensitive).
    
    Returns:
        The actual column name as it exists in adata.obs, or None if not found.
    """
    # First check exact match
    if column_name in adata.obs.columns:
        return column_name
    
    # Check case-insensitive match
    for col in adata.obs.columns:
        if col.lower() == column_name.lower():
            return col
    
    # Not found
    return None


def main():
    # Parse arguments using argparse
    parser = argparse.ArgumentParser(
        description='Cluster QC\'d scRNA cells using PCA + Leiden/Louvain and generate optional dimensionality reduction plots.',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Positional arguments
    parser.add_argument('qc_h5ad', type=str, help='Path to QC\'d h5ad file')
    parser.add_argument('out_csv', type=str, help='Path to output CSV file for cluster assignments')
    
    # Optional plotting flags
    parser.add_argument('--pca', action='store_true', help='Generate PCA plots')
    parser.add_argument('--tsne', action='store_true', help='Generate t-SNE plots')
    parser.add_argument('--umap', action='store_true', help='Generate UMAP plots')
    parser.add_argument('--color-by', type=str, default='leiden', 
                        help='Variable to color plots by (default: leiden)')
    parser.add_argument('--save-figures', action='store_true',
                        help='Save generated plots to outputs/clustering/')
    parser.add_argument('--output-dir', type=str, default='outputs/clustering',
                        help='Directory to save figures (default: outputs/clustering)')
    
    # Clustering parameters
    parser.add_argument('--n-neighbors', type=int, default=10,
                        help='Number of neighbors for neighborhood graph (default: 10). '
                             'Higher values produce more connected graphs.')
    parser.add_argument('--n-pcs', type=int, default=40,
                        help='Number of principal components to use (default: 40). '
                             'Should be less than min(n_cells, n_genes).')
    parser.add_argument('--resolution', type=float, default=0.5,
                        help='Leiden clustering resolution (default: 0.5). '
                             'Higher values produce more clusters.')
    parser.add_argument('--force-recluster', action='store_true',
                        help='Force re-clustering even if leiden results already exist in the dataset')
    
    args = parser.parse_args()
    
    # Validate parameter ranges
    if args.n_neighbors < 2:
        print("ERROR: --n-neighbors must be at least 2", file=sys.stderr)
        sys.exit(1)
    if args.n_pcs < 1:
        print("ERROR: --n-pcs must be at least 1", file=sys.stderr)
        sys.exit(1)
    if args.resolution <= 0:
        print("ERROR: --resolution must be positive", file=sys.stderr)
        sys.exit(1)

    h5ad_path = Path(args.qc_h5ad)
    if not h5ad_path.is_file():
        print(f"ERROR: QC h5ad not found: {h5ad_path}", file=sys.stderr)
        sys.exit(1)

    # Load and validate h5ad file
    print(f"\nLoading dataset: {h5ad_path}")
    try:
        adata = sc.read_h5ad(h5ad_path)
        print(f"  ✓ Loaded {adata.n_obs} cells × {adata.n_vars} genes")
    except Exception as e:
        print(f"ERROR: Failed to read h5ad file '{h5ad_path}': {e}", file=sys.stderr)
        sys.exit(1)
    
    # Validate basic structure
    if adata.n_obs == 0 or adata.n_vars == 0:
        print(
            f"ERROR: AnnData loaded from {h5ad_path} is empty "
            f"({adata.n_obs} cells × {adata.n_vars} genes).",
            file=sys.stderr,
        )
        sys.exit(1)
    
    if adata.X is None:
        print(f"ERROR: AnnData object has no expression matrix (X is None)", file=sys.stderr)
        sys.exit(1)
    
    # Check for minimum cell/gene counts for meaningful clustering
    MIN_CELLS = 50
    MIN_GENES = 100
    if adata.n_obs < MIN_CELLS:
        print(
            f"WARNING: Dataset has only {adata.n_obs} cells (recommended minimum: {MIN_CELLS}). "
            "Clustering results may not be reliable.",
            file=sys.stderr
        )
    if adata.n_vars < MIN_GENES:
        print(
            f"WARNING: Dataset has only {adata.n_vars} genes (recommended minimum: {MIN_GENES}). "
            "Clustering results may not be reliable.",
            file=sys.stderr
        )
    
    # Validate output directory if saving figures
    if args.save_figures:
        output_dir = Path(args.output_dir)
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
        except (PermissionError, OSError) as e:
            print(
                f"ERROR: Cannot create output directory '{output_dir}': {e}",
                file=sys.stderr
            )
            sys.exit(1)
        
        # Test write permissions
        test_file = output_dir / ".write_test"
        try:
            test_file.touch()
            test_file.unlink()
        except (PermissionError, OSError) as e:
            print(
                f"ERROR: No write permission for output directory '{output_dir}': {e}",
                file=sys.stderr
            )
            sys.exit(1)

    # Display clustering parameters
    print(f"\nClustering parameters:")
    print(f"  Neighbors: {args.n_neighbors}")
    print(f"  PCs: {args.n_pcs}")
    print(f"  Resolution: {args.resolution}")
    
    # Normalize and transform
    print("\nNormalizing and transforming data...")
    try:
        adata = normalize_and_transform_adata(adata, show_plots=False)
        print("  ✓ Normalization complete")
    except Exception as e:
        print(
            f"ERROR: Normalization and transformation failed: {e}\n"
            "This may indicate corrupted data or incompatible format.",
            file=sys.stderr
        )
        sys.exit(1)
    
    # Check if leiden clustering already exists (case-insensitive) - do this early
    leiden_col = None
    if 'leiden' in adata.obs.columns:
        leiden_col = 'leiden'
    elif 'Leiden' in adata.obs.columns:
        leiden_col = 'Leiden'
    
    # Validate n_pcs against dataset dimensions
    max_pcs = min(adata.n_obs - 1, adata.n_vars - 1)
    if args.n_pcs > max_pcs:
        print(
            f"WARNING: --n-pcs ({args.n_pcs}) exceeds dataset dimensions. "
            f"Using maximum possible: {max_pcs}",
            file=sys.stderr
        )
        n_pcs_to_use = max_pcs
    else:
        n_pcs_to_use = args.n_pcs
    
    # Compute PCA (but defer plotting until after clustering if needed)
    print("\nComputing PCA...")
    if 'X_pca' not in adata.obsm:
        sc.tl.pca(adata, svd_solver="arpack")
        print("  ✓ PCA computed")
    else:
        print("  ✓ Using existing PCA from input file")
    
    # Build neighborhood graph (required for clustering)
    print("\nBuilding neighborhood graph...")
    try:
        sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, n_pcs=n_pcs_to_use)
        print(f"  ✓ Computed neighbors (k={args.n_neighbors}, n_pcs={n_pcs_to_use})")
    except Exception as e:
        print(
            f"ERROR: Failed to compute neighborhood graph: {e}\n"
            "This typically indicates issues with PCA results or data structure.",
            file=sys.stderr
        )
        sys.exit(1)
    
    # Determine if we need to run clustering
    should_cluster = leiden_col is None or args.force_recluster
    
    if should_cluster:
        if args.force_recluster and leiden_col:
            print(f"\nRe-clustering (--force-recluster enabled, overwriting existing '{leiden_col}' column)...")
        else:
            print("\nPerforming Leiden clustering...")
        
        # Run leiden clustering
        try:
            sc.tl.leiden(adata, resolution=args.resolution)
            leiden_col = 'leiden'
            n_clusters = len(adata.obs['leiden'].unique())
            print(f"  ✓ Identified {n_clusters} clusters (resolution={args.resolution})")
        except ImportError as e:
            print(
                "ERROR: Leiden clustering requires `python-igraph` and `leidenalg` "
                "to be installed in the environment.",
                file=sys.stderr,
            )
            print(f"Install with: pip install python-igraph leidenalg", file=sys.stderr)
            print(f"Underlying error: {e}", file=sys.stderr)
            sys.exit(1)
        except Exception as e:
            print(
                f"ERROR: Leiden clustering failed: {e}\n"
                "This may indicate insufficient connectivity in the neighborhood graph.\n"
                "Try adjusting n_neighbors or ensure you have enough cells.",
                file=sys.stderr
            )
            sys.exit(1)
    else:
        print(f"Using existing {leiden_col} clustering from input file")
    
    # Now validate the color_by column (after clustering is done)
    color_by = normalize_column_name(adata, args.color_by)
    if color_by is None:
        # Column doesn't exist at all
        available_cols = sorted(adata.obs.columns.tolist())
        print(
            f"ERROR: Column '{args.color_by}' not found in dataset.\n"
            f"Available columns: {', '.join(available_cols[:10])}"
            f"{' ...' if len(available_cols) > 10 else ''}",
            file=sys.stderr
        )
        sys.exit(1)
    
    if color_by != args.color_by:
        print(f"Note: Using '{color_by}' column (matched from '{args.color_by}')")

    # Save clustering results to CSV
    try:
        clusters = pd.DataFrame(
            {"cell": adata.obs_names, "cluster": adata.obs[leiden_col].astype(str)}
        )
        clusters.to_csv(args.out_csv, index=False)
        print(
            f"Clustered {len(clusters)} cells into {len(set(clusters['cluster']))} clusters -> {args.out_csv}"
        )
    except (PermissionError, OSError) as e:
        print(
            f"ERROR: Failed to write clustering results to '{args.out_csv}': {e}",
            file=sys.stderr
        )
        sys.exit(1)
    except Exception as e:
        print(
            f"ERROR: Failed to create clustering dataframe: {e}",
            file=sys.stderr
        )
        sys.exit(1)

    # Create outputs directory if it doesn't exist
    os.makedirs("outputs", exist_ok=True)
    
    # Generate optional dimensionality reduction plots
    # Use legend_loc="on data" only when coloring by leiden clusters
    legend_loc = "on data" if (leiden_col and color_by.lower() == leiden_col.lower()) else None
    
    if args.pca:
        print("Generating PCA plot...")
        try:
            adata = run_pca_analysis(adata, color_by=color_by, show_plots=True,
                                    save_figures=args.save_figures, output_dir=args.output_dir,
                                    legend_loc=legend_loc)
        except KeyError as e:
            print(
                f"ERROR: PCA plotting failed - column not found: {e}\n"
                f"Ensure '{color_by}' exists in the dataset.",
                file=sys.stderr
            )
        except Exception as e:
            print(
                f"ERROR: PCA plotting failed: {e}",
                file=sys.stderr
            )
    
    if args.tsne:
        print("Generating t-SNE plot...")
        try:
            compute_tsne(adata, color_by=color_by, show_plots=True,
                        save_figures=args.save_figures, output_dir=args.output_dir,
                        legend_loc=legend_loc)
        except KeyError as e:
            print(
                f"ERROR: t-SNE plotting failed - column not found: {e}\n"
                f"Ensure '{color_by}' exists in the dataset.",
                file=sys.stderr
            )
        except Exception as e:
            print(
                f"ERROR: t-SNE computation or plotting failed: {e}\n"
                "This may indicate issues with PCA results or insufficient memory.",
                file=sys.stderr
            )
    
    if args.umap:
        print("Generating UMAP plot...")
        try:
            compute_umap(adata, color_by=color_by, show_plots=True,
                        save_figures=args.save_figures, output_dir=args.output_dir,
                        legend_loc=legend_loc)
        except KeyError as e:
            print(
                f"ERROR: UMAP plotting failed - column not found: {e}\n"
                f"Ensure '{color_by}' exists in the dataset.",
                file=sys.stderr
            )
        except Exception as e:
            print(
                f"ERROR: UMAP computation or plotting failed: {e}\n"
                "This may indicate issues with neighborhood graph or insufficient memory.",
                file=sys.stderr
            )
    
    # Always save a static UMAP plot colored by leiden for the default pipeline output
    # Skip if user already requested UMAP with leiden coloring (case-insensitive)
    if not (args.umap and color_by.lower() == leiden_col.lower() and args.save_figures):
        try:
            if 'X_umap' not in adata.obsm:
                sc.tl.umap(adata)
            fig = sc.pl.umap(
                adata, color=leiden_col, legend_loc="on data", return_fig=True, show=False
            )
            if isinstance(fig, (list, tuple)) and len(fig) > 0:
                fig = fig[0]
            plt.figure(fig.number)
            plt.savefig("outputs/umap_clusters.png", dpi=200, bbox_inches="tight")
            plt.close(fig)
            print("Saved UMAP plot -> outputs/umap_clusters.png")
        except Exception as e:
            # Non-fatal: clustering results are still valid even if plotting fails
            print(
                f"WARNING: Default UMAP plot generation failed: {e}\n"
                "Clustering results were still saved successfully.",
                file=sys.stderr
            )
    
    # Final summary
    print("\n" + "="*60)
    print("CLUSTERING COMPLETE")
    print("="*60)
    print(f"Results saved to: {args.out_csv}")
    if args.save_figures:
        print(f"Figures saved to: {args.output_dir}")
    print("="*60)


if __name__ == "__main__":
    main()
