#!/usr/bin/env python3
"""
Run the cluster_cells.py script to cluster QC'd scRNA cells.

Can be used to cluster cells from a QC'd h5ad file using the functions
in cluster_cells.py. Options allow generating PCA, t-SNE, and UMAP plots.

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

import argparse
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

from . import cluster_cells as cc


def main():
    # Parse arguments using argparse
    parser = argparse.ArgumentParser(
        description="Cluster QC'd scRNA cells using PCA + Leiden/Louvain and generate optional dimensionality reduction plots.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Positional arguments
    parser.add_argument("qc_h5ad", type=str, help="Path to QC'd h5ad file")
    parser.add_argument(
        "out_csv", type=str, help="Path to output CSV file for cluster assignments"
    )

    # Optional plotting flags
    parser.add_argument("--pca", action="store_true", help="Generate PCA plots")
    parser.add_argument("--tsne", action="store_true", help="Generate t-SNE plots")
    parser.add_argument("--umap", action="store_true", help="Generate UMAP plots")
    parser.add_argument(
        "--color-by",
        type=str,
        default="leiden",
        help="Variable to color plots by (default: leiden)",
    )
    parser.add_argument(
        "--save-figures",
        action="store_true",
        help="Save generated plots to outputs/clustering/",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="outputs/clustering",
        help="Directory to save figures (default: outputs/clustering)",
    )

    # Clustering parameters
    parser.add_argument(
        "--n-neighbors",
        type=int,
        default=10,
        help="Number of neighbors for neighborhood graph (default: 10). "
        "Higher values produce more connected graphs.",
    )
    parser.add_argument(
        "--n-pcs",
        type=int,
        default=40,
        help="Number of principal components to use (default: 40). "
        "Should be less than min(n_cells, n_genes).",
    )
    parser.add_argument(
        "--resolution",
        type=float,
        default=0.5,
        help="Leiden clustering resolution (default: 0.5). "
        "Higher values produce more clusters.",
    )
    parser.add_argument(
        "--force-recluster",
        action="store_true",
        help="Force re-clustering even if leiden results already exist in the dataset",
    )

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
        print(
            "ERROR: AnnData object has no expression matrix (X is None)",
            file=sys.stderr,
        )
        sys.exit(1)

    # Check for minimum cell/gene counts for meaningful clustering
    MIN_CELLS = 50
    MIN_GENES = 100
    if adata.n_obs < MIN_CELLS:
        print(
            f"WARNING: Dataset has only {adata.n_obs} cells (recommended minimum: {MIN_CELLS}). "
            "Clustering results may not be reliable.",
            file=sys.stderr,
        )
    if adata.n_vars < MIN_GENES:
        print(
            f"WARNING: Dataset has only {adata.n_vars} genes (recommended minimum: {MIN_GENES}). "
            "Clustering results may not be reliable.",
            file=sys.stderr,
        )

    # Validate output directory if saving figures
    if args.save_figures:
        output_dir = Path(args.output_dir)
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
        except (PermissionError, OSError) as e:
            print(
                f"ERROR: Cannot create output directory '{output_dir}': {e}",
                file=sys.stderr,
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
                file=sys.stderr,
            )
            sys.exit(1)

    # Display clustering parameters
    print("\nClustering parameters:")
    print(f"  Neighbors: {args.n_neighbors}")
    print(f"  PCs: {args.n_pcs}")
    print(f"  Resolution: {args.resolution}")

    # Normalize and transform
    print("\nNormalizing and transforming data...")
    try:
        adata = cc.normalize_and_transform_adata(adata, show_plots=False)
        print("  ✓ Normalization complete")
    except Exception as e:
        print(
            f"ERROR: Normalization and transformation failed: {e}\n"
            "This may indicate corrupted data or incompatible format.",
            file=sys.stderr,
        )
        sys.exit(1)

    # Check if leiden clustering already exists (case-insensitive) - do this early
    leiden_col = None
    if "leiden" in adata.obs.columns:
        leiden_col = "leiden"
    elif "Leiden" in adata.obs.columns:
        leiden_col = "Leiden"

    # Validate n_pcs against dataset dimensions
    max_pcs = min(adata.n_obs - 1, adata.n_vars - 1)
    if max_pcs < 1:
        max_pcs = 1  # safeguard for extremely tiny datasets

    if args.n_pcs > max_pcs:
        print(
            f"WARNING: --n-pcs ({args.n_pcs}) exceeds dataset dimensions. "
            f"Using maximum possible: {max_pcs}",
            file=sys.stderr,
        )
        n_pcs_to_use = max_pcs
    else:
        n_pcs_to_use = args.n_pcs

    # Decide whether we can/should run PCA
    use_pca = False

    # Estimate how many features are available for PCA
    if "highly_variable" in adata.var:
        n_vars_for_pca = int(adata.var["highly_variable"].sum())
    else:
        n_vars_for_pca = adata.n_vars

    if n_vars_for_pca <= 0:
        print(
            "WARNING: No usable features available for PCA (0 highly variable genes). "
            "Skipping PCA and using raw expression for neighbors/clustering.",
            file=sys.stderr,
        )
        n_pcs_to_use = 0
    elif n_pcs_to_use <= 0:
        print(
            "WARNING: n_pcs_to_use <= 0. Skipping PCA and using raw expression for "
            "neighbors/clustering.",
            file=sys.stderr,
        )
        n_pcs_to_use = 0
    else:
        # Compute PCA (but defer plotting until after clustering if needed)
        print("\nComputing PCA...")
        if "X_pca" not in adata.obsm:
            sc.tl.pca(adata, n_comps=n_pcs_to_use, svd_solver="arpack")
            print("  ✓ PCA computed")
        else:
            print("  ✓ Using existing PCA from input file")
        use_pca = True

    # Build neighborhood graph (required for clustering)
    print("\nBuilding neighborhood graph...")
    try:
        if use_pca and n_pcs_to_use > 0:
            sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, n_pcs=n_pcs_to_use)
        else:
            # no PCA / no PCs to use — compute neighbors directly on the (possibly small) X
            sc.pp.neighbors(adata, n_neighbors=args.n_neighbors)
        print(f"  ✓ Computed neighbors (k={args.n_neighbors}, n_pcs={n_pcs_to_use})")
    except Exception as e:
        print(
            f"ERROR: Failed to compute neighborhood graph: {e}\n"
            "This typically indicates issues with PCA results or data structure.",
            file=sys.stderr,
        )
        sys.exit(1)

    # Determine if we need to run clustering
    should_cluster = leiden_col is None or args.force_recluster

    if should_cluster:
        if args.force_recluster and leiden_col:
            print(
                f"\nRe-clustering (--force-recluster enabled, overwriting existing '{leiden_col}' column)..."
            )
        else:
            print("\nPerforming Leiden clustering...")

        # Run leiden clustering
        try:
            sc.tl.leiden(adata, resolution=args.resolution)
            leiden_col = "leiden"
            n_clusters = len(adata.obs["leiden"].unique())
            print(
                f"  ✓ Identified {n_clusters} clusters (resolution={args.resolution})"
            )
        except ImportError as e:
            print(
                "ERROR: Leiden clustering requires `python-igraph` and `leidenalg` "
                "to be installed in the environment.",
                file=sys.stderr,
            )
            print("Install with: pip install python-igraph leidenalg", file=sys.stderr)
            print(f"Underlying error: {e}", file=sys.stderr)
            sys.exit(1)
        except Exception as e:
            print(
                f"ERROR: Leiden clustering failed: {e}\n"
                "This may indicate insufficient connectivity in the neighborhood graph.\n"
                "Try adjusting n_neighbors or ensure you have enough cells.",
                file=sys.stderr,
            )
            sys.exit(1)
    else:
        print(f"Using existing {leiden_col} clustering from input file")

    # Now validate the color_by column (after clustering is done)
    color_by = cc.normalize_column_name(adata, args.color_by)
    if color_by is None:
        # Column doesn't exist at all
        available_cols = sorted(adata.obs.columns.tolist())
        print(
            f"ERROR: Column '{args.color_by}' not found in dataset.\n"
            f"Available columns: {', '.join(available_cols[:10])}"
            f"{' ...' if len(available_cols) > 10 else ''}",
            file=sys.stderr,
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
            file=sys.stderr,
        )
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Failed to create clustering dataframe: {e}", file=sys.stderr)
        sys.exit(1)

    # Create outputs directory if it doesn't exist
    os.makedirs("outputs", exist_ok=True)

    # Generate optional dimensionality reduction plots
    # Use legend_loc="on data" only when coloring by leiden clusters
    legend_loc = (
        "on data" if (leiden_col and color_by.lower() == leiden_col.lower()) else None
    )

    if args.pca:
        print("Generating PCA plot...")
        try:
            adata = cc.run_pca_analysis(
                adata,
                color_by=color_by,
                show_plots=True,
                save_figures=args.save_figures,
                output_dir=args.output_dir,
                legend_loc=legend_loc,
            )
        except KeyError as e:
            print(
                f"ERROR: PCA plotting failed - column not found: {e}\n"
                f"Ensure '{color_by}' exists in the dataset.",
                file=sys.stderr,
            )
        except Exception as e:
            print(f"ERROR: PCA plotting failed: {e}", file=sys.stderr)

    if args.tsne:
        print("Generating t-SNE plot...")
        try:
            cc.compute_tsne(
                adata,
                color_by=color_by,
                show_plots=True,
                save_figures=args.save_figures,
                output_dir=args.output_dir,
                legend_loc=legend_loc,
                n_pcs=n_pcs_to_use,
            )
        except KeyError as e:
            print(
                f"ERROR: t-SNE plotting failed - column not found: {e}\n"
                f"Ensure '{color_by}' exists in the dataset.",
                file=sys.stderr,
            )
        except Exception as e:
            print(
                f"ERROR: t-SNE computation or plotting failed: {e}\n"
                "This may indicate issues with PCA results or insufficient memory.",
                file=sys.stderr,
            )

    if args.umap:
        print("Generating UMAP plot...")
        try:
            cc.compute_umap(
                adata,
                color_by=color_by,
                show_plots=True,
                save_figures=args.save_figures,
                output_dir=args.output_dir,
                legend_loc=legend_loc,
                n_pcs=n_pcs_to_use,
            )
        except KeyError as e:
            print(
                f"ERROR: UMAP plotting failed - column not found: {e}\n"
                f"Ensure '{color_by}' exists in the dataset.",
                file=sys.stderr,
            )
        except Exception as e:
            print(
                f"ERROR: UMAP computation or plotting failed: {e}\n"
                "This may indicate issues with neighborhood graph or insufficient memory.",
                file=sys.stderr,
            )

    # Always save a static UMAP plot colored by leiden for the default pipeline output
    # Skip if user already requested UMAP with leiden coloring (case-insensitive)
    if not (args.umap and color_by.lower() == leiden_col.lower() and args.save_figures):
        try:
            if "X_umap" not in adata.obsm:
                sc.tl.umap(adata)
            fig = sc.pl.umap(
                adata,
                color=leiden_col,
                legend_loc="on data",
                return_fig=True,
                show=False,
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
                file=sys.stderr,
            )

    # Final summary
    print("\n" + "=" * 60)
    print("CLUSTERING COMPLETE")
    print("=" * 60)
    print(f"Results saved to: {args.out_csv}")
    if args.save_figures:
        print(f"Figures saved to: {args.output_dir}")
    print("=" * 60)


if __name__ == "__main__":
    main()
