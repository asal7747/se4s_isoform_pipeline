#!/usr/bin/env python
"""Cluster QC'd scRNA cells using PCA + Leiden/Louvain and generate optional dimensionality reduction plots.

Default behavior produces `outputs/umap_clusters.png` after clustering completes.

Can be run as:
    python cluster_cells.py <qc_h5ad> <out_csv> [--pca] [--tsne] [--umap] [--color-by VARIABLE]

Another example:
    python3 scripts/cluster_cells.py outputs/anndata/combined_short_read_qc.h5ad clusters2.csv --resolution 0.8 --n-pcs 15 --force-recluster --save-figures --pca --tsne --umap

    This would create:
    - outputs/clustering/pca_<color_by>.png
    - outputs/clustering/tsne_<color_by>.png
    - outputs/clustering/umap_<color_by>.png
"""

import os
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc


def normalize_and_transform_adata(
    adata: sc.AnnData,
    target_sum: float = 1e4,
    min_mean: float = 0.0125,
    max_mean: float = 3,
    min_disp: float = 0.5,
    output_dir: str = "outputs/anndata",
    show_plots: bool = True,
    skip_if_normalized: bool = True,
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
        print(
            "Data appears already normalized (adata.raw exists). Skipping normalization."
        )
        is_normalized = True
    elif skip_if_normalized and "highly_variable" in adata.var.columns:
        print(
            "Data appears already processed (highly_variable genes marked). Skipping normalization."
        )
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
        sc.pp.highly_variable_genes(
            adata, min_mean=min_mean, max_mean=max_mean, min_disp=min_disp
        )

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
        dataset_name = adata.uns.get("dataset_name", "combined_dataset")
        norm_h5ad_path = output_path / f"{dataset_name}_norm_log_transformed.h5ad"
        adata.write(norm_h5ad_path)
        print(f"Normalized and transformed AnnData saved to: {norm_h5ad_path}")

    return adata


def run_pca_analysis(
    adata,
    svd_solver="arpack",
    color_by="SampleType",
    n_pcs=50,
    annotate_var_explained=True,
    show_plots=True,
    save_figures=False,
    output_dir="outputs/clustering",
    legend_loc=None,
):
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
    if "X_pca" not in adata.obsm:
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
        sc.pl.pca(
            adata,
            annotate_var_explained=annotate_var_explained,
            color=actual_color_by,
            show=False if save_figures else True,
        )
        if save_figures:
            plt.savefig(
                f"{output_dir}/pca_{color_by}.png", dpi=200, bbox_inches="tight"
            )
            plt.close()
            print(f"Saved PCA plot -> {output_dir}/pca_{color_by}.png")

        # inspect the contribution of single PCs to the total variance in the data
        # Make sure we don't request more PCs than were computed
        n_pcs_available = adata.obsm['X_pca'].shape[1]
        n_pcs_to_plot = min(n_pcs, n_pcs_available)
        
        if n_pcs_to_plot > 0:
            sc.pl.pca_variance_ratio(
                adata, log=True, n_pcs=n_pcs_to_plot, show=False if save_figures else True
            )
            if save_figures:
                plt.savefig(
                    f"{output_dir}/pca_variance_ratio.png", dpi=200, bbox_inches="tight"
                )
                plt.close()
                print(
                    f"Saved PCA variance ratio plot -> {output_dir}/pca_variance_ratio.png"
                )

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
    output_dir: str = "outputs/clustering",
    legend_loc: str | None = None,
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
    if "X_tsne" not in adata.obsm:
        sc.tl.tsne(
            adata,
            perplexity=perplexity,
            learning_rate=learning_rate,
            random_state=random_state,
            n_pcs=n_pcs,
        )
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

            sc.pl.tsne(adata, color=actual_color, legend_loc=legend_loc, 
                      show=False if save_figures else True)
            if save_figures:
                plt.savefig(
                    f"{output_dir}/tsne_{color}.png", dpi=200, bbox_inches="tight"
                )
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
    legend_loc: str | None = None,
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
    if "X_umap" not in adata.obsm:
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

            sc.pl.umap(
                adata,
                color=actual_color,
                legend_loc=legend_loc,
                show=False if save_figures else True,
            )
            if save_figures:
                plt.savefig(
                    f"{output_dir}/umap_{color}.png", dpi=200, bbox_inches="tight"
                )
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


def transfer_umap_coordinates(source_adata, target_adata, plot_output_dir=None, color_by=None):
    """Transfer UMAP coordinates from source to target AnnData based on matching cell barcodes.
    
    Args:
        source_adata: AnnData with UMAP coordinates (e.g., short-read data)
        target_adata: AnnData to transfer coordinates to (e.g., long-read data)
        plot_output_dir: Optional directory to save comparison plots
        color_by: Optional variable to color cells by (must exist in both .obs)
        
    Returns:
        target_adata with transferred UMAP coordinates in .obsm['X_umap_transferred']
    """
    import numpy as np
    
    # Check that source has UMAP
    if 'X_umap' not in source_adata.obsm:
        raise ValueError("Source AnnData does not have UMAP coordinates. Run sc.tl.umap() first.")
    
    # Find matching barcodes
    source_barcodes = set(source_adata.obs_names)
    target_barcodes = set(target_adata.obs_names)
    common_barcodes = source_barcodes & target_barcodes
    
    print(f"\nTransferring UMAP coordinates...")
    print(f"  Source cells: {len(source_barcodes)}")
    print(f"  Target cells: {len(target_barcodes)}")
    print(f"  Common cells: {len(common_barcodes)}")
    
    if len(common_barcodes) == 0:
        raise ValueError("No common barcodes between source and target datasets!")
    
    # Create a mapping of barcode to UMAP coordinates
    umap_coords = pd.DataFrame(
        source_adata.obsm['X_umap'],
        index=source_adata.obs_names,
        columns=['UMAP1', 'UMAP2']
    )
    
    # Transfer coordinates to target (only for matching cells)
    target_umap = np.full((target_adata.n_obs, 2), np.nan)
    for i, barcode in enumerate(target_adata.obs_names):
        if barcode in common_barcodes:
            target_umap[i] = umap_coords.loc[barcode].values
    
    # Store in target
    target_adata.obsm['X_umap_transferred'] = target_umap
    
    # Debug: Check if color_by column exists in target
    if color_by:
        if color_by in target_adata.obs.columns:
            print(f"  Color by '{color_by}' found in target data:")
            print(f"    {target_adata.obs[color_by].value_counts().to_dict()}")
        else:
            print(f"  WARNING: '{color_by}' not found in target data columns")
            print(f"  Available columns: {list(target_adata.obs.columns)}")
    
    # Create comparison plot if requested
    if plot_output_dir is not None:
        plot_output_dir = Path(plot_output_dir)
        plot_output_dir.mkdir(parents=True, exist_ok=True)
        
        # Create side-by-side comparison
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        
        # Plot source data
        ax1 = axes[0]
        if color_by and color_by in source_adata.obs.columns:
            # Color by variable
            categories = source_adata.obs[color_by].unique()
            colors = plt.cm.tab20(np.linspace(0, 1, len(categories)))
            for i, cat in enumerate(categories):
                mask = source_adata.obs[color_by] == cat
                ax1.scatter(
                    source_adata.obsm['X_umap'][mask, 0],
                    source_adata.obsm['X_umap'][mask, 1],
                    c=[colors[i]], label=cat, s=10, alpha=0.7
                )
            ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
        else:
            ax1.scatter(
                source_adata.obsm['X_umap'][:, 0],
                source_adata.obsm['X_umap'][:, 1],
                c='blue', s=10, alpha=0.5
            )
        ax1.set_xlabel('UMAP1', fontsize=12)
        ax1.set_ylabel('UMAP2', fontsize=12)
        ax1.set_title(f'Source Data (n={source_adata.n_obs})', fontsize=14, fontweight='bold')
        ax1.grid(alpha=0.3)
        
        # Plot target data with transferred coordinates
        ax2 = axes[1]
        mask = ~np.isnan(target_umap[:, 0])  # Only plot cells with transferred coords
        if color_by and color_by in target_adata.obs.columns:
            # Get categories, filtering out NaN and 'Unknown' values for cleaner visualization
            categories = target_adata.obs.loc[target_adata.obs_names[mask], color_by].unique()
            # Filter out None/NaN values
            categories = [cat for cat in categories if pd.notna(cat)]
            
            if len(categories) > 0:
                colors = plt.cm.tab20(np.linspace(0, 1, len(categories)))
                for i, cat in enumerate(categories):
                    cat_mask = (target_adata.obs[color_by] == cat) & mask
                    if cat_mask.sum() > 0:  # Only plot if there are cells in this category
                        ax2.scatter(
                            target_umap[cat_mask, 0],
                            target_umap[cat_mask, 1],
                            c=[colors[i]], label=cat, s=10, alpha=0.7
                        )
                ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
            else:
                # Fallback if no valid categories
                ax2.scatter(
                    target_umap[mask, 0],
                    target_umap[mask, 1],
                    c='red', s=10, alpha=0.5
                )
        else:
            ax2.scatter(
                target_umap[mask, 0],
                target_umap[mask, 1],
                c='red', s=10, alpha=0.5
            )
        ax2.set_xlabel('UMAP1', fontsize=12)
        ax2.set_ylabel('UMAP2', fontsize=12)
        ax2.set_title(f'Target Data with Transferred UMAP (n={mask.sum()})', fontsize=14, fontweight='bold')
        ax2.grid(alpha=0.3)
        
        plt.tight_layout()
        output_path = plot_output_dir / "umap_transfer_comparison.png"
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {output_path}")
    
    return target_adata


def main():
    """Main entry point for the cluster_cells module.
    
    This function is intentionally simple and serves as a placeholder.
    The actual CLI functionality is in run_cluster_cells.py.
    """
    print("cluster_cells.py: This module provides clustering functions.")
    print("To run clustering from the command line, use: python run_cluster_cells.py")


if __name__ == "__main__":
    main()
