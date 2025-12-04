#!/usr/bin/env python
"""
Single-cell RNA-seq analysis pipeline.
This is meant to be an example of how to use utils.py and cluster_cells.py.

This script performs the complete analysis workflow:
1. Load and combine short-read datasets
2. Load and combine long-read datasets
3. QC and filter short-read data
4. Normalize and transform data
5. Run dimensionality reduction (PCA, t-SNE, UMAP)
6. Transfer UMAP coordinates from short-read to long-read data

Only writes combined long-read datasets to outputs/anndata, as well as
combined short-read dataset after QC/filtering.
See utils.py and cluster_cells.py for helper functions.
"""

import os

import anndata
import cluster_cells as cc
import utils


def main():
    """Run the complete single-cell analysis pipeline."""

    # Load and combine short read datasets
    print("Loading and combining short-read datasets...")
    combined_short_read_adata, combined_short_read_df = (
        utils.load_and_combine_short_reads()
    )

    # Load long read anndata
    print("\nLoading long-read datasets...")
    long_read_datasets = utils.load_long_read_datasets(data_dir="data")
    for name, adata in long_read_datasets.items():
        print(f"Dataset: {name}")
        print(adata)
        print("\n")

    # Instantiate individual long read anndata objects
    long_gene = long_read_datasets["long_gene"]
    long_transcript = long_read_datasets["long_transcript"]
    long_myotube_gene = long_read_datasets["long_myotube_gene"]
    long_myotube_transcript = long_read_datasets["long_myotube_transcript"]
    long_nuc_gene = long_read_datasets["long_nuc_gene"]
    long_nuc_transcript = long_read_datasets["long_nuc_transcript"]

    # Datasets ending in _transcript provide transcript-level quantification
    # Datasets ending in _gene provide gene-level quantification
    # We want to examine both separately,
    # so we will make 2 different combined anndata objects

    # Combine long read anndata objects ending in "transcript" or "gene" into 1
    # Use merge='same' to keep var metadata that is the same across all datasets
    print("\nCombining long-read datasets...")
    combined_long_read_adata_transcript = anndata.concat(
        [long_transcript, long_myotube_transcript, long_nuc_transcript],
        axis=0,
        merge="same",
    )
    combined_long_read_adata_gene = anndata.concat(
        [long_gene, long_myotube_gene, long_nuc_gene],
        axis=0,
        merge="same",
    )
    print(
        f"Combined long read transcript anndata shape: "
        f"{combined_long_read_adata_transcript.shape}"
    )
    print(
        f"Combined long read gene anndata shape: {combined_long_read_adata_gene.shape}"
    )

    # Add combined dataset name to uns for later use
    combined_long_read_adata_transcript.uns["dataset_name"] = (
        "combined_long_read_transcript"
    )
    combined_long_read_adata_gene.uns["dataset_name"] = "combined_long_read_gene"

    # Do these all have matching short reads? This is the 464 cell library
    # Check if barcode IDs in long read anndata are present
    # in combined short read anndata
    print("\nChecking for matching barcodes between long and short reads...")
    combined_long_read_adata_gene_names = set(combined_long_read_adata_gene.obs_names)
    combined_short_read_barcodes = set(combined_short_read_adata.obs_names)
    matching_barcodes = combined_long_read_adata_gene_names.intersection(
        combined_short_read_barcodes
    )
    print(
        f"Number of barcodes in long read gene: "
        f"{len(combined_long_read_adata_gene_names)}"
    )
    print(
        f"Number of barcodes in combined_short_read_adata: "
        f"{len(combined_short_read_barcodes)}"
    )
    print(f"Number of matching barcodes: "
          f"{len(matching_barcodes)}")

    # Transfer SampleType metadata from short-read
    # to long-read based on matching barcodes
    print("\nTransferring SampleType metadata from short-read to long-read data...")
    for barcode in combined_long_read_adata_gene.obs_names:
        if barcode in combined_short_read_adata.obs_names:
            # Get SampleType from short-read data
            sample_type = combined_short_read_adata.obs.loc[barcode, "SampleType"]
            combined_long_read_adata_gene.obs.loc[barcode, "SampleType"] = sample_type
        else:
            # Mark as unknown if not found in short-read data
            combined_long_read_adata_gene.obs.loc[barcode, "SampleType"] = "Unknown"

    # Do the same for transcript-level data
    for barcode in combined_long_read_adata_transcript.obs_names:
        if barcode in combined_short_read_adata.obs_names:
            sample_type = combined_short_read_adata.obs.loc[barcode, "SampleType"]
            combined_long_read_adata_transcript.obs.loc[barcode, "SampleType"] = (
                sample_type
            )
        else:
            combined_long_read_adata_transcript.obs.loc[barcode, "SampleType"] = (
                "Unknown"
            )

    # Print summary of transferred SampleType
    print("Long-read gene SampleType distribution:")
    print(combined_long_read_adata_gene.obs["SampleType"].value_counts())
    print("\nLong-read transcript SampleType distribution:")
    print(combined_long_read_adata_transcript.obs["SampleType"].value_counts())

    # Save combined long-read datasets WITH SampleType metadata
    print("\nSaving combined long-read datasets with SampleType metadata...")
    os.makedirs("outputs/anndata", exist_ok=True)
    combined_long_read_adata_transcript.write_h5ad(
        "outputs/anndata/combined_long_read_transcript.h5ad"
    )
    combined_long_read_adata_gene.write_h5ad(
        "outputs/anndata/combined_long_read_gene.h5ad"
    )
    print("  Saved: outputs/anndata/combined_long_read_transcript.h5ad")
    print("  Saved: outputs/anndata/combined_long_read_gene.h5ad")

    # Test QC and filter function on combined short read anndata
    print("\nRunning QC and filtering on short-read data...")
    qc_combined_short_read_adata = utils.qc_and_filter(
        combined_short_read_adata, output_dir="outputs/anndata"
    )

    # Will not QC/filter long read data since it is already very sparse

    # Normalize and log transform combined short read anndata
    print("\nNormalizing and transforming data...")
    norm_combined_short_read_adata = cc.normalize_and_transform_adata(
        qc_combined_short_read_adata, output_dir="outputs/anndata"
    )

    # Run PCA analysis
    print("\nRunning PCA analysis...")
    norm_combined_short_read_adata = cc.run_pca_analysis(
        norm_combined_short_read_adata, color_by="SampleType"
    )

    # Make t-SNE plot
    print("\nGenerating t-SNE plots...")
    norm_combined_short_read_adata = cc.compute_tsne(
        norm_combined_short_read_adata, color_by=["SampleType", "CellType"]
    )

    # Try UMAP plot
    print("\nGenerating UMAP plots...")
    norm_combined_short_read_adata = cc.compute_umap(
        norm_combined_short_read_adata, color_by="SampleType"
    )

    # Transfer UMAP coordinates to long-read data and create comparison plots
    print("\nTransferring UMAP coordinates from short-read to long-read data...")
    combined_long_read_adata_gene = cc.transfer_umap_coordinates(
        source_adata=norm_combined_short_read_adata,
        target_adata=combined_long_read_adata_gene,
        plot_output_dir="outputs/anndata",
        color_by="SampleType",
    )

    # Save updated long-read data with transferred UMAP
    print("\nSaving long-read data with transferred UMAP coordinates...")
    combined_long_read_adata_gene.write_h5ad(
        "outputs/anndata/combined_long_read_gene_with_umap.h5ad"
    )
    print("  Saved: outputs/anndata/combined_long_read_gene_with_umap.h5ad")

    print("\nAnalysis complete!")


if __name__ == "__main__":
    main()
