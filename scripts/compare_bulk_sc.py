#!/usr/bin/env python3
"""
Compare transcript diversity between bulk and single-cell long-read data.

This script analyzes:
1. Gene/transcript overlap between bulk TALON and single-cell h5ad datasets
2. Isoform diversity per gene in bulk data
3. Transcript novelty distribution in bulk data
4. Expression patterns in single-cell data

Usage:
    python compare_bulk_sc.py \
        --bulk-talon tests/bulk_sc_talon_read_annot.tsv \
        --output outputs/comparison/

With custom paths:
    python compare_bulk_sc.py \
        --bulk-talon tests/bulk_sc_talon_read_annot.tsv \
        --sc-gene-data outputs/anndata/combined_long_read_gene.h5ad \
        --sc-transcript-data outputs/anndata/
        combined_long_read_transcript.h5ad \
        --output outputs/comparison/

Optional:
    --min-reads 5
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Import from local utils
try:
    from utils import get_anndata, load_long_read_datasets
except ImportError:
    # Try relative import if running as a script
    import sys

    try:
        script_dir = Path(__file__).parent
    except NameError:
        # __file__ not defined in interactive mode, use cwd
        script_dir = Path.cwd() / "scripts"
    sys.path.insert(0, str(script_dir))
    from utils import get_anndata


def load_talon_data(talon_path, min_reads=5):
    """Load TALON TSV and extract gene/isoform information."""
    print(f"Loading TALON data from {talon_path}...")
    df = pd.read_csv(talon_path, sep="\t")

    # Filter by minimum reads if specified
    if min_reads > 0 and "read_count" in df.columns:
        df = df[df["read_count"] >= min_reads]

    print(f"  Loaded {len(df)} transcript records")

    # Extract key columns (adjust based on your TALON output format)
    required_cols = ["annot_gene_name", "annot_transcript_id"]
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        print(f"WARNING: Missing columns in TALON data: {missing}")
        print(f"Available columns: {list(df.columns)}")

    return df


def extract_gene_set_from_talon(
    talon_df, dataset_label="dataset", exclude_talon_ids=True
):
    """Extract unique genes, transcripts, and transcript names from TALON df.

    Args:
        talon_df: TALON dataframe
        dataset_label: Label for print statements
        exclude_talon_ids: If True,
        exclude gene/transcript names starting with 'TALON' (default: True)
    """
    # Extract genes, optionally filtering TALON-generated gene names
    all_genes = set(talon_df["annot_gene_name"].dropna().unique())
    if exclude_talon_ids:
        genes = {gene for gene in all_genes if not str(gene).startswith("TALON")}
        excluded_genes = len(all_genes) - len(genes)
    else:
        genes = all_genes
        excluded_genes = 0

    transcripts = set(talon_df["annot_transcript_id"].dropna().unique())

    # Check for transcript name column
    transcript_names = None
    for col in ["annot_transcript_name", "transcript_name", "transcript"]:
        if col in talon_df.columns:
            all_names = set(talon_df[col].dropna().unique())

            if exclude_talon_ids:
                # Filter out names that start with 'TALON'
                transcript_names = {
                    name for name in all_names if not str(name).startswith("TALON")
                }
                excluded_count = len(all_names) - len(transcript_names)
                print(
                    f"{dataset_label}: "
                    f"{len(genes)} genes ({excluded_genes} TALON genes "
                    f"excluded), {len(transcripts)} transcripts"
                )
                print(
                    f"{len(all_names)} total transcript names "
                    f"({excluded_count} TALON IDs excluded)"
                )
                print(f"{len(transcript_names)} annotated transcript names (non-TALON)")
            else:
                transcript_names = all_names
                print(
                    f"{dataset_label}: {len(genes)} genes, "
                    f"{len(transcripts)} transcripts, "
                    f"{len(transcript_names)} transcript names"
                )

            return genes, transcripts, transcript_names

    print(f"{dataset_label}: {len(genes)} genes, {len(transcripts)} transcripts")
    return genes, transcripts, None


def extract_gene_set_from_h5ad(adata, dataset_label="dataset"):
    """Extract unique genes from single-cell h5ad AnnData object."""
    # Try multiple possible locations for gene names
    if "gene_name" in adata.var.columns:
        genes = set(adata.var["gene_name"].dropna().unique())
    elif "gene_id" in adata.var.columns:
        genes = set(adata.var["gene_id"].dropna().unique())
    else:
        # Fall back to var_names (index)
        genes = set(adata.var_names)

    print(f"  {dataset_label}: {len(genes)} genes from {adata.n_obs} cells")
    return genes


def extract_transcript_set_from_h5ad(adata, dataset_label="dataset"):
    """Extract unique transcripts and their genes from
    transcript-level h5ad.

    Note: Uses var_names (index) as the primary transcript identifier for
    comparisons, since this typically contains gene-based names
    (e.g., Mrpl15-206) that match TALON's annot_transcript_name format.
    """
    # For transcript-level data, use var_names (index) as primary identifier
    # This is typically in gene-based format like "Mrpl15-206" which matches
    # TALON annot_transcript_name
    transcripts = set(adata.var_names)

    # Try to get associated gene names
    if "gene_name" in adata.var.columns:
        genes = set(adata.var["gene_name"].dropna().unique())
    elif "gene_id" in adata.var.columns:
        genes = set(adata.var["gene_id"].dropna().unique())
    else:
        # Extract gene names from transcript IDs ("Mrpl15-206" -> "Mrpl15")
        # This is a common format in single-cell data where transcript IDs are
        # gene-isoform
        genes = set()
        for tid in transcripts:
            # Split on hyphen and take the gene name part
            if "-" in str(tid):
                gene_name = str(tid).rsplit("-", 1)[0]
                genes.add(gene_name)

        if genes:
            print(
                "NOTE: Extracted gene names from transcript IDs "
                "(format: GeneName-IsoformNumber)"
            )

    print(f"{dataset_label}: {len(transcripts)} transcripts from {adata.n_obs} cells")
    if genes:
        print(f"Associated with {len(genes)} genes")

    # Check for novelty information
    novelty_col = None
    for col in ["transcript_novelty", "ISM_subtype", "novelty", "transcript_status"]:
        if col in adata.var.columns:
            novelty_col = col
            print(f"    ✓ Found novelty column: '{col}'")
            # Show distribution if available
            novelty_counts = adata.var[col].value_counts()
            print(f"      Categories: {dict(novelty_counts.head(5))}")
            break

    if novelty_col is None:
        print("    ✗ No transcript novelty information found")
        print(f"      Available columns: {list(adata.var.columns)[:10]}")

    return transcripts, genes, adata.var


def compare_gene_overlap(bulk_genes, sc_genes, output_dir):
    """Compare gene overlap between bulk and single-cell datasets."""
    print("\nComparing gene overlap...")

    overlap = bulk_genes & sc_genes
    bulk_only = bulk_genes - sc_genes
    sc_only = sc_genes - bulk_genes

    print(f"  Genes in both datasets: {len(overlap)}")
    print(f"  Bulk only: {len(bulk_only)}")
    print(f"  Single-cell only: {len(sc_only)}")
    print(f"  Jaccard index: {len(overlap) / len(bulk_genes | sc_genes):.3f}")

    # Create comparison plot
    fig, ax = plt.subplots(figsize=(8, 6))

    categories = ["Both\ndatasets", "Bulk\nonly", "Single-cell\nonly"]
    counts = [len(overlap), len(bulk_only), len(sc_only)]
    colors = ["#2ecc71", "#e74c3c", "#3498db"]

    ax.bar(categories, counts, color=colors, alpha=0.7, edgecolor="black")
    ax.set_ylabel("Number of genes", fontsize=12)
    ax.set_title(
        "Gene Overlap: Bulk vs. Single-Cell Long-read", fontsize=14, fontweight="bold"
    )
    ax.grid(axis="y", alpha=0.3)

    # Add count labels on bars
    for i, (cat, count) in enumerate(zip(categories, counts)):
        ax.text(
            i,
            count + max(counts) * 0.02,
            str(count),
            ha="center",
            va="bottom",
            fontsize=11,
            fontweight="bold",
        )

    plt.tight_layout()
    plt.savefig(output_dir / "gene_overlap.png", dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {output_dir / 'gene_overlap.png'}")

    return overlap


def compare_transcript_overlap(bulk_transcript_names, sc_transcript_ids, output_dir):
    """Compare transcript name overlap between bulk and single-cell datasets.

    Args:
        bulk_transcript_names: Set of transcript names from bulk TALON
        sc_transcript_ids: Set of transcript IDs from SC
        output_dir: Directory to save plots

    Returns:
        Set of overlapping transcript names
    """
    print("\nComparing transcript name overlap...")

    overlap = bulk_transcript_names & sc_transcript_ids
    bulk_only = bulk_transcript_names - sc_transcript_ids
    sc_only = sc_transcript_ids - bulk_transcript_names

    print(f"  Transcripts in both datasets: {len(overlap)}")
    print(f"  Bulk only: {len(bulk_only)}")
    print(f"  Single-cell only: {len(sc_only)}")
    if len(bulk_transcript_names | sc_transcript_ids) > 0:
        print(
            f"  Jaccard index: {
                len(overlap) / len(bulk_transcript_names | sc_transcript_ids):.3f}"
        )

    # Create comparison plot
    fig, ax = plt.subplots(figsize=(8, 6))

    categories = ["Both\ndatasets", "Bulk\nonly", "Single-cell\nonly"]
    counts = [len(overlap), len(bulk_only), len(sc_only)]
    colors = ["#9b59b6", "#e74c3c", "#3498db"]

    ax.bar(categories, counts, color=colors, alpha=0.7, edgecolor="black")
    ax.set_ylabel("Number of transcripts", fontsize=12)
    ax.set_title(
        "Transcript Name Overlap: Bulk vs. Single-Cell Long-read",
        fontsize=14,
        fontweight="bold",
    )
    ax.grid(axis="y", alpha=0.3)

    # Add count labels on bars
    for i, (cat, count) in enumerate(zip(categories, counts)):
        ax.text(
            i,
            count + max(counts) * 0.02 if max(counts) > 0 else 1,
            str(count),
            ha="center",
            va="bottom",
            fontsize=11,
            fontweight="bold",
        )

    # Add percentage labels
    total = len(bulk_transcript_names | sc_transcript_ids)
    if total > 0:
        for i, count in enumerate(counts):
            pct = count / total * 100
            ax.text(
                i,
                count / 2,
                f"{pct:.1f}%",
                ha="center",
                va="center",
                fontsize=10,
                fontweight="bold",
                color="white",
            )

    plt.tight_layout()
    plt.savefig(output_dir / "transcript_overlap.png", dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {output_dir / 'transcript_overlap.png'}")

    return overlap


def compare_transcript_names(bulk_transcript_names, sc_transcript_ids, bulk_df=None):
    """Compare transcript names between bulk and SC, trying to find matches.

    Args:
        bulk_transcript_names: Set of transcript names from bulk TALON
        sc_transcript_ids: Set of transcript IDs from SC
            (format: GeneName-IsoformNumber)
        bulk_df: Optional TALON dataframe to show matched examples

    Returns:
        Tuple of (exact_matches, gene_matches)
    """
    print("\nComparing transcript names...")

    if bulk_transcript_names is None:
        print("  No transcript names available in bulk TALON data")
        return set(), set()

    # Debug: show types and sizes
    print(
        f"  Bulk transcript names: {len(bulk_transcript_names)} items, "
        f"type: {type(bulk_transcript_names)}"
    )
    print(
        f"  SC transcript IDs: {len(sc_transcript_ids)} items, type: "
        f"{type(sc_transcript_ids)}"
    )

    # Try exact matches
    exact_matches = bulk_transcript_names & sc_transcript_ids
    print(f"  Exact transcript name matches: {len(exact_matches)}")
    if len(exact_matches) > 0:
        print(f"    Examples: {list(exact_matches)[:5]}")
    else:
        # Show a few examples to see format differences
        print(f"    Bulk examples: {list(bulk_transcript_names)[:3]}")
        print(f"    SC examples: {list(sc_transcript_ids)[:3]}")

    # Try gene-level matching (extract gene from SC transcript ID)
    bulk_genes_from_names = set()
    for name in bulk_transcript_names:
        # Try to extract gene name (common formats:
        # GeneName-IsoformNum, GeneName.IsoformNum)
        if "-" in str(name):
            gene = str(name).rsplit("-", 1)[0]
            bulk_genes_from_names.add(gene)
        elif "." in str(name):
            gene = str(name).rsplit(".", 1)[0]
            bulk_genes_from_names.add(gene)

    sc_genes_from_ids = set()
    for tid in sc_transcript_ids:
        if "-" in str(tid):
            gene = str(tid).rsplit("-", 1)[0]
            sc_genes_from_ids.add(gene)

    gene_overlap = bulk_genes_from_names & sc_genes_from_ids
    print(f"  Genes with matching transcript names: {len(gene_overlap)}")
    if len(gene_overlap) > 0:
        print(f"    Examples: {list(gene_overlap)[:10]}")

    # Show examples of bulk transcript names
    print(f"\n  Sample bulk transcript names: {list(bulk_transcript_names)[:10]}")
    print(f"  Sample SC transcript IDs: {list(sc_transcript_ids)[:10]}")

    return exact_matches, gene_overlap


def compare_isoform_diversity_h5ad(
    bulk_df, sc_transcript_ids, common_genes, output_dir, exclude_talon_ids=True
):
    """Compare isoform diversity between bulk TALON and single-cell data.

    Args:
        bulk_df: TALON df with annot_gene_name and annot_transcript_id columns
        sc_transcript_ids: Set of transcript IDs from SC data
            (format: GeneName-IsoformNumber)
        common_genes: Set of genes present in both datasets
        output_dir: Directory to save plots
        exclude_talon_ids: If True, only count annotated transcripts
            (exclude TALON IDs) for bulk
    """
    print("\nComparing isoform diversity...")

    # Filter bulk to common genes
    bulk_common = bulk_df[bulk_df["annot_gene_name"].isin(common_genes)]

    # Count isoforms per gene in bulk
    if exclude_talon_ids and "annot_transcript_name" in bulk_df.columns:
        # Filter to only annotated transcripts and genes
        # (not starting with TALON)
        bulk_annotated = bulk_common[
            (~bulk_common["annot_transcript_name"].str.startswith("TALON", na=False))
            & (~bulk_common["annot_gene_name"].str.startswith("TALON", na=False))
        ]
        bulk_isoforms = bulk_annotated.groupby("annot_gene_name")[
            "annot_transcript_name"
        ].nunique()
        print(
            "  Note: Using annotated genes and transcript names only "
            "(excluding TALON IDs)"
        )
        print(f"    Bulk transcripts before filtering: {len(bulk_common)}")
        print(f"    Bulk transcripts after filtering: {len(bulk_annotated)}")
    else:
        bulk_isoforms = bulk_common.groupby("annot_gene_name")[
            "annot_transcript_id"
        ].nunique()

    # Count isoforms per gene in SC by parsing transcript IDs
    sc_gene_isoforms = {}
    for tid in sc_transcript_ids:
        if "-" in str(tid):
            gene_name = str(tid).rsplit("-", 1)[0]
            if gene_name in common_genes:
                sc_gene_isoforms[gene_name] = sc_gene_isoforms.get(gene_name, 0) + 1

    sc_isoforms = pd.Series(sc_gene_isoforms)

    print("\n  Bulk dataset:")
    print(f"    Mean isoforms per gene: {bulk_isoforms.mean():.2f}")
    print(f"    Median isoforms per gene: {bulk_isoforms.median():.0f}")
    print(f"    Max isoforms per gene: {bulk_isoforms.max()}")

    print("\n  Single-cell dataset:")
    print(f"    Mean isoforms per gene: {sc_isoforms.mean():.2f}")
    print(f"    Median isoforms per gene: {sc_isoforms.median():.0f}")
    print(f"    Max isoforms per gene: {sc_isoforms.max()}")

    # Create comparison DataFrame
    comparison = pd.DataFrame(
        {
            "gene": list(set(bulk_isoforms.index) | set(sc_isoforms.index)),
        }
    )
    comparison["bulk_isoforms"] = (
        comparison["gene"].map(bulk_isoforms).fillna(0).astype(int)
    )
    comparison["sc_isoforms"] = (
        comparison["gene"].map(sc_isoforms).fillna(0).astype(int)
    )
    comparison["difference"] = comparison["bulk_isoforms"] - comparison["sc_isoforms"]
    comparison = comparison.sort_values("bulk_isoforms", ascending=False)

    # Plot comparison
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # 1. Histogram comparison
    ax1 = axes[0, 0]
    ax1.hist(
        bulk_isoforms,
        bins=30,
        alpha=0.6,
        label="Bulk",
        color="#e74c3c",
        edgecolor="black",
    )
    ax1.hist(
        sc_isoforms,
        bins=30,
        alpha=0.6,
        label="Single-cell",
        color="#3498db",
        edgecolor="black",
    )
    ax1.set_xlabel("Number of isoforms per gene", fontsize=11)
    ax1.set_ylabel("Number of genes", fontsize=11)
    ax1.set_title("Isoform Diversity Distribution", fontsize=12, fontweight="bold")
    ax1.legend()
    ax1.grid(alpha=0.3)

    # 2. Scatter plot
    ax2 = axes[0, 1]
    genes_both = set(bulk_isoforms.index) & set(sc_isoforms.index)
    bulk_vals = [bulk_isoforms.get(g, 0) for g in genes_both]
    sc_vals = [sc_isoforms.get(g, 0) for g in genes_both]

    ax2.scatter(bulk_vals, sc_vals, alpha=0.5, s=30, color="#9b59b6")
    max_val = max(max(bulk_vals) if bulk_vals else 1, max(sc_vals) if sc_vals else 1)
    ax2.plot([0, max_val], [0, max_val], "k--", linewidth=1, alpha=0.5, label="y=x")
    ax2.set_xlabel("Bulk isoforms per gene", fontsize=11)
    ax2.set_ylabel("Single-cell isoforms per gene", fontsize=11)
    ax2.set_title("Isoform Count Correlation", fontsize=12, fontweight="bold")
    ax2.legend()
    ax2.grid(alpha=0.3)

    # 3. Top genes in bulk
    ax3 = axes[1, 0]
    top_bulk = comparison.nlargest(15, "bulk_isoforms")
    x_pos = range(len(top_bulk))
    ax3.barh(x_pos, top_bulk["bulk_isoforms"], alpha=0.6, label="Bulk", color="#e74c3c")
    ax3.barh(
        x_pos, top_bulk["sc_isoforms"], alpha=0.6, label="Single-cell", color="#3498db"
    )
    ax3.set_yticks(x_pos)
    ax3.set_yticklabels(top_bulk["gene"], fontsize=9)
    ax3.set_xlabel("Number of isoforms", fontsize=11)
    ax3.set_title("Top 15 Genes by Bulk Isoform Count", fontsize=12, fontweight="bold")
    ax3.legend()
    ax3.grid(axis="x", alpha=0.3)
    ax3.invert_yaxis()

    # 4. Difference distribution
    ax4 = axes[1, 1]
    ax4.hist(
        comparison["difference"], bins=30, color="#34495e", alpha=0.7, edgecolor="black"
    )
    ax4.axvline(0, color="red", linestyle="--", linewidth=2, label="Equal diversity")
    ax4.set_xlabel("Isoform count difference (Bulk - SC)", fontsize=11)
    ax4.set_ylabel("Number of genes", fontsize=11)
    ax4.set_title("Isoform Diversity Difference", fontsize=12, fontweight="bold")
    ax4.legend()
    ax4.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(
        output_dir / "isoform_diversity_comparison.png", dpi=300, bbox_inches="tight"
    )
    plt.close()
    print(f"  Saved: {output_dir / 'isoform_diversity_comparison.png'}")

    # Save comparison table
    comparison.to_csv(output_dir / "isoform_diversity_per_gene.csv", index=False)
    print(f"  Saved: {output_dir / 'isoform_diversity_per_gene.csv'}")

    return comparison


def compare_isoform_diversity(bulk_df, sc_df, common_genes, output_dir):
    """Compare isoform diversity per gene between bulk and single-cell."""
    print("\nComparing isoform diversity...")

    # Filter to common genes
    bulk_common = bulk_df[bulk_df["annot_gene_name"].isin(common_genes)]
    sc_common = sc_df[sc_df["annot_gene_name"].isin(common_genes)]

    # Count isoforms per gene in each dataset
    bulk_isoforms = bulk_common.groupby("annot_gene_name")[
        "annot_transcript_id"
    ].nunique()
    sc_isoforms = sc_common.groupby("annot_gene_name")["annot_transcript_id"].nunique()

    print("\n  Bulk dataset:")
    print(f"    Mean isoforms per gene: {bulk_isoforms.mean():.2f}")
    print(f"    Median isoforms per gene: {bulk_isoforms.median():.0f}")
    print(f"    Max isoforms per gene: {bulk_isoforms.max()}")

    print("\n  Single-cell dataset:")
    print(f"    Mean isoforms per gene: {sc_isoforms.mean():.2f}")
    print(f"    Median isoforms per gene: {sc_isoforms.median():.0f}")
    print(f"    Max isoforms per gene: {sc_isoforms.max()}")

    # Create comparison DataFrame
    comparison = pd.DataFrame(
        {
            "gene": list(set(bulk_isoforms.index) | set(sc_isoforms.index)),
        }
    )
    comparison["bulk_isoforms"] = (
        comparison["gene"].map(bulk_isoforms).fillna(0).astype(int)
    )
    comparison["sc_isoforms"] = (
        comparison["gene"].map(sc_isoforms).fillna(0).astype(int)
    )
    comparison["difference"] = comparison["bulk_isoforms"] - comparison["sc_isoforms"]
    comparison = comparison.sort_values("bulk_isoforms", ascending=False)

    # Plot comparison
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # 1. Histogram comparison
    ax1 = axes[0, 0]
    ax1.hist(
        bulk_isoforms,
        bins=30,
        alpha=0.6,
        label="Bulk",
        color="#e74c3c",
        edgecolor="black",
    )
    ax1.hist(
        sc_isoforms,
        bins=30,
        alpha=0.6,
        label="Single-cell",
        color="#3498db",
        edgecolor="black",
    )
    ax1.set_xlabel("Number of isoforms per gene", fontsize=11)
    ax1.set_ylabel("Number of genes", fontsize=11)
    ax1.set_title("Isoform Diversity Distribution", fontsize=12, fontweight="bold")
    ax1.legend()
    ax1.grid(alpha=0.3)

    # 2. Scatter plot
    ax2 = axes[0, 1]
    genes_both = set(bulk_isoforms.index) & set(sc_isoforms.index)
    bulk_vals = [bulk_isoforms.get(g, 0) for g in genes_both]
    sc_vals = [sc_isoforms.get(g, 0) for g in genes_both]

    ax2.scatter(bulk_vals, sc_vals, alpha=0.5, s=30, color="#9b59b6")
    max_val = max(max(bulk_vals), max(sc_vals))
    ax2.plot([0, max_val], [0, max_val], "k--", linewidth=1, alpha=0.5, label="y=x")
    ax2.set_xlabel("Bulk isoforms per gene", fontsize=11)
    ax2.set_ylabel("Single-cell isoforms per gene", fontsize=11)
    ax2.set_title("Isoform Count Correlation", fontsize=12, fontweight="bold")
    ax2.legend()
    ax2.grid(alpha=0.3)

    # 3. Top genes in bulk
    ax3 = axes[1, 0]
    top_bulk = comparison.nlargest(15, "bulk_isoforms")
    y_pos = range(len(top_bulk))
    ax3.barh(
        y_pos,
        top_bulk["bulk_isoforms"],
        alpha=0.7,
        color="#e74c3c",
        edgecolor="black",
        label="Bulk",
    )
    ax3.barh(
        y_pos,
        top_bulk["sc_isoforms"],
        alpha=0.7,
        color="#3498db",
        edgecolor="black",
        label="Single-cell",
    )
    ax3.set_yticks(y_pos)
    ax3.set_yticklabels(top_bulk["gene"], fontsize=9)
    ax3.set_xlabel("Number of isoforms", fontsize=11)
    ax3.set_title("Top 15 Genes by Bulk Isoform Count", fontsize=12, fontweight="bold")
    ax3.legend()
    ax3.grid(axis="x", alpha=0.3)
    ax3.invert_yaxis()

    # 4. Boxplot comparison
    ax4 = axes[1, 1]
    box_data = [bulk_isoforms.values, sc_isoforms.values]
    bp = ax4.boxplot(box_data, labels=["Bulk", "Single-cell"], patch_artist=True)
    bp["boxes"][0].set_facecolor("#e74c3c")
    bp["boxes"][1].set_facecolor("#3498db")
    ax4.set_ylabel("Isoforms per gene", fontsize=11)
    ax4.set_title("Isoform Diversity Summary", fontsize=12, fontweight="bold")
    ax4.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    plt.savefig(
        output_dir / "isoform_diversity_comparison.png", dpi=300, bbox_inches="tight"
    )
    plt.close()
    print(f"  Saved: {output_dir / 'isoform_diversity_comparison.png'}")

    # Save comparison table
    comparison.to_csv(output_dir / "isoforms_comparison.csv", index=False)
    print(f"  Saved: {output_dir / 'isoforms_comparison.csv'}")

    return bulk_isoforms, sc_isoforms


def compare_transcript_novelty(bulk_df, sc_df, output_dir):
    """Compare transcript novelty categories between bulk and single-cell.

    Args:
        bulk_df: TALON dataframe with novelty column
        sc_df: SC var dataframe with 'transcript_novelty' column
               (standardized name)
        output_dir: Directory to save plots
    """
    print("\nComparing transcript novelty...")

    # Check if transcript_novelty column exists in bulk
    novelty_col = None
    for col in ["transcript_novelty", "ISM_subtype", "novelty"]:
        if col in bulk_df.columns:
            novelty_col = col
            break

    if novelty_col is None:
        print(
            "  WARNING: No novelty column found in bulk data."
            " Skipping novelty analysis."
        )
        return None, None

    # Use standardized 'transcript_novelty' column name for SC
    if "transcript_novelty" not in sc_df.columns:
        print(
            "  WARNING: No 'transcript_novelty' column in SC data."
            " Skipping novelty analysis."
        )
        return None, None

    bulk_novelty = bulk_df[novelty_col].value_counts()
    sc_novelty = sc_df["transcript_novelty"].value_counts()

    print("\n  Bulk transcript categories:")
    for category, count in bulk_novelty.items():
        print(f"    {category}: {count} ({count / len(bulk_df) * 100:.1f}%)")

    print("\n  Single-cell transcript categories:")
    for category, count in sc_novelty.items():
        print(f"    {category}: {count} ({count / len(sc_df) * 100:.1f}%)")

    # Create comparison plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    # Get all categories
    all_categories = sorted(set(bulk_novelty.index) | set(sc_novelty.index))

    # Bulk novelty
    bulk_counts = [bulk_novelty.get(cat, 0) for cat in all_categories]
    bulk_pcts = [count / len(bulk_df) * 100 for count in bulk_counts]
    colors1 = plt.cm.Set3(range(len(all_categories)))

    bars1 = ax1.bar(
        range(len(all_categories)),
        bulk_counts,
        color=colors1,
        edgecolor="black",
        alpha=0.8,
    )
    ax1.set_xticks(range(len(all_categories)))
    ax1.set_xticklabels(all_categories, rotation=45, ha="right", fontsize=10)
    ax1.set_ylabel("Number of transcripts", fontsize=12)
    ax1.set_title("Bulk Long-read Novelty Categories", fontsize=13, fontweight="bold")
    ax1.grid(axis="y", alpha=0.3)

    for i, (count, pct) in enumerate(zip(bulk_counts, bulk_pcts)):
        if count > 0:
            ax1.text(
                i,
                count + max(bulk_counts) * 0.02,
                f"{pct:.1f}%",
                ha="center",
                va="bottom",
                fontsize=9,
            )

    # Single-cell novelty
    sc_counts = [sc_novelty.get(cat, 0) for cat in all_categories]
    sc_pcts = [count / len(sc_df) * 100 for count in sc_counts]

    bars2 = ax2.bar(
        range(len(all_categories)),
        sc_counts,
        color=colors1,
        edgecolor="black",
        alpha=0.8,
    )
    ax2.set_xticks(range(len(all_categories)))
    ax2.set_xticklabels(all_categories, rotation=45, ha="right", fontsize=10)
    ax2.set_ylabel("Number of transcripts", fontsize=12)
    ax2.set_title(
        "Single-Cell Long-read Novelty Categories", fontsize=13, fontweight="bold"
    )
    ax2.grid(axis="y", alpha=0.3)

    for i, (count, pct) in enumerate(zip(sc_counts, sc_pcts)):
        if count > 0:
            ax2.text(
                i,
                count + max(sc_counts) * 0.02,
                f"{pct:.1f}%",
                ha="center",
                va="bottom",
                fontsize=9,
            )

    plt.tight_layout()
    plt.savefig(
        output_dir / "transcript_novelty_comparison.png", dpi=300, bbox_inches="tight"
    )
    plt.close()
    print(f"  Saved: {output_dir / 'transcript_novelty_comparison.png'}")

    return bulk_novelty, sc_novelty


def find_cluster_specific_isoforms(
    sc_transcript_adata,
    cluster_col="SampleType",
    output_dir=None,
    top_n=5,
    min_expression=5,
):
    """Find genes with cluster-specific isoform usage patterns.

    Args:
        sc_transcript_adata: AnnData object with transcript-level data
        cluster_col: Column name for cluster annotations
        output_dir: Directory to save plots
        top_n: Number of top genes to show in plots
        min_expression: Minimum total expression to consider a transcript

    Returns:
        DataFrame with cluster-specific isoform information
    """
    print("\nFinding cluster-specific isoform patterns...")

    if cluster_col not in sc_transcript_adata.obs.columns:
        print(f"  WARNING: '{cluster_col}' not found. Skipping.")
        return None

    import numpy as np

    # Get clusters in specific order: MB_cells, MB_nuclei, MT_nuclei
    all_clusters = sc_transcript_adata.obs[cluster_col].unique()
    desired_order = ["MB_cells", "MB_nuclei", "MT_nuclei"]
    clusters = [c for c in desired_order if c in all_clusters]
    # Add any remaining clusters not in desired order
    clusters.extend([c for c in all_clusters if c not in desired_order])

    # For each gene, find dominant isoform in each cluster
    gene_isoform_patterns = []

    # Get all genes
    all_genes = set()
    for tid in sc_transcript_adata.var_names:
        if "-" in str(tid):
            gene = str(tid).rsplit("-", 1)[0]
            all_genes.add(gene)

    print(f"  Analyzing {len(all_genes)} genes across {len(clusters)} clusters...")

    # Pre-group isoforms by gene for efficient lookup
    gene_to_isoforms = {}
    for tid in sc_transcript_adata.var_names:
        if "-" in str(tid):
            gene = str(tid).rsplit("-", 1)[0]
            if gene not in gene_to_isoforms:
                gene_to_isoforms[gene] = []
            gene_to_isoforms[gene].append(tid)

    # Pre-compute cluster masks and subset data once per cluster
    cluster_masks = {
        cluster: sc_transcript_adata.obs[cluster_col] == cluster for cluster in clusters
    }

    for gene in all_genes:
        # Get all isoforms of this gene
        gene_isoforms = gene_to_isoforms.get(gene, [])

        if len(gene_isoforms) < 2:
            continue  # Skip genes with only one isoform

        # Get expression per cluster per isoform
        cluster_isoform_expr = {}

        for cluster in clusters:
            cluster_mask = cluster_masks[cluster]
            cluster_data = sc_transcript_adata[cluster_mask, gene_isoforms]

            # Sum expression for each isoform in this cluster
            if hasattr(cluster_data.X, "toarray"):
                expr_sums = np.asarray(cluster_data.X.toarray()).sum(axis=0)
            else:
                expr_sums = np.asarray(cluster_data.X.sum(axis=0)).flatten()

            isoform_exprs = dict(zip(gene_isoforms, expr_sums))

            # Find dominant isoform for this cluster
            if sum(isoform_exprs.values()) >= min_expression:
                dominant = max(isoform_exprs, key=isoform_exprs.get)
                cluster_isoform_expr[cluster] = {
                    "dominant_isoform": dominant,
                    "expression": isoform_exprs[dominant],
                    "total_expression": sum(isoform_exprs.values()),
                }

        # Check if different clusters have different dominant isoforms
        if len(cluster_isoform_expr) >= 2:
            dominant_isoforms = [
                info["dominant_isoform"] for info in cluster_isoform_expr.values()
            ]

            if len(set(dominant_isoforms)) > 1:  # Isoform switching detected!
                gene_isoform_patterns.append(
                    {
                        "gene": gene,
                        "num_isoforms": len(gene_isoforms),
                        **{
                            f"{cluster}_isoform": cluster_isoform_expr.get(
                                cluster, {}
                            ).get("dominant_isoform", "None")
                            for cluster in clusters
                        },
                        **{
                            f"{cluster}_expr": cluster_isoform_expr.get(
                                cluster, {}
                            ).get("total_expression", 0)
                            for cluster in clusters
                        },
                    }
                )

    if len(gene_isoform_patterns) == 0:
        print("  No cluster-specific isoform switching found.")
        return None

    switching_df = pd.DataFrame(gene_isoform_patterns)
    switching_df = switching_df.sort_values(
        [f"{clusters[0]}_expr"], ascending=False
    ).head(top_n)

    print(
        f"  Found {len(gene_isoform_patterns)} genes with cluster-specific "
        "isoform usage"
    )
    print(f"  Top {top_n} shown by expression")

    # Create visualization
    if output_dir:
        output_dir = Path(output_dir)

        fig, ax = plt.subplots(figsize=(10, max(6, top_n * 0.8)))

        # Collect all unique isoform numbers and assign colors
        all_isoform_nums = set()
        for i, (idx, row) in enumerate(switching_df.iterrows()):
            for cluster in clusters:
                isoform = row[f"{cluster}_isoform"]
                if isoform != "None":
                    isoform_num = isoform.split("-")[-1] if "-" in isoform else isoform
                    all_isoform_nums.add(isoform_num)

        # Create color map for isoform numbers
        isoform_nums_sorted = sorted(all_isoform_nums)
        isoform_colors = {}
        cmap = plt.cm.get_cmap("tab20")
        for idx, iso_num in enumerate(isoform_nums_sorted):
            isoform_colors[iso_num] = cmap(idx % 20)

        # For each gene, show which isoform is dominant in each cluster
        for i, (idx, row) in enumerate(switching_df.iterrows()):
            gene = row["gene"]

            for j, cluster in enumerate(clusters):
                isoform = row[f"{cluster}_isoform"]
                if isoform != "None":
                    # Extract isoform number
                    isoform_num = isoform.split("-")[-1] if "-" in isoform else isoform
                    # Use color based on isoform number
                    ax.text(
                        j,
                        i,
                        isoform_num,
                        ha="center",
                        va="center",
                        fontsize=10,
                        fontweight="bold",
                        bbox=dict(
                            boxstyle="round,pad=0.5",
                            facecolor=isoform_colors[isoform_num],
                            edgecolor="black",
                            linewidth=1.5,
                            alpha=0.8,
                        ),
                    )

        ax.set_xlim(-0.5, len(clusters) - 0.5)
        ax.set_ylim(-0.5, len(switching_df) - 0.5)
        ax.set_xticks(range(len(clusters)))
        ax.set_xticklabels(clusters, fontsize=11, fontweight="bold")
        ax.set_yticks(range(len(switching_df)))
        ax.set_yticklabels(switching_df["gene"], fontsize=10)
        ax.set_xlabel("Cluster", fontsize=12, fontweight="bold")
        ax.set_ylabel("Gene", fontsize=12, fontweight="bold")
        ax.set_title(
            f"Top {top_n} Genes with Cluster-Specific Isoform Switching\n"
            f"(Colors indicate same isoform across clusters)",
            fontsize=13,
            fontweight="bold",
        )
        ax.invert_yaxis()
        ax.grid(True, alpha=0.3, linestyle="--")

        plt.tight_layout()
        plt.savefig(
            output_dir / "cluster_isoform_switching.png", dpi=300, bbox_inches="tight"
        )
        plt.close()
        print(f"  Saved: {output_dir / 'cluster_isoform_switching.png'}")

        # Save table
        switching_df.to_csv(output_dir / "cluster_isoform_switching.csv", index=False)
        print(f"  Saved: {output_dir / 'cluster_isoform_switching.csv'}")

    return switching_df


def compare_bulk_vs_cluster_isoforms(
    bulk_df,
    sc_transcript_adata,
    cluster_col="SampleType",
    output_dir=None,
    exclude_talon_ids=True,
):
    """Compare bulk isoform diversity to cluster-specific diversity.

    Args:
        bulk_df: TALON dataframe
        sc_transcript_adata: AnnData object with transcript-level data
        cluster_col: Column name for cluster annotations
        output_dir: Directory to save plots
        exclude_talon_ids: Whether to exclude TALON IDs from bulk

    Returns:
        Comparison DataFrame
    """
    print("\nComparing bulk vs cluster-specific isoform diversity...")

    # Check if cluster column exists
    if cluster_col not in sc_transcript_adata.obs.columns:
        print(f"  WARNING: '{cluster_col}' not found. Skipping.")
        return None

    # Compute cluster diversity
    clusters = sc_transcript_adata.obs[cluster_col].unique()
    cluster_diversity = {}

    print(f"  Computing isoform diversity for {len(clusters)} clusters...")
    for cluster in clusters:
        # Get cells in this cluster
        cluster_mask = sc_transcript_adata.obs[cluster_col] == cluster
        cluster_data = sc_transcript_adata[cluster_mask, :]

        # Get expressed transcripts
        if hasattr(cluster_data.X, "toarray"):
            expr_sums = np.asarray(cluster_data.X.sum(axis=0)).flatten()
        else:
            expr_sums = cluster_data.X.sum(axis=0)

        expressed_mask = expr_sums > 0
        expressed_transcripts = set(cluster_data.var_names[expressed_mask])

        # Count isoforms per gene
        gene_isoforms = {}
        for tid in expressed_transcripts:
            if "-" in str(tid):
                gene_name = str(tid).rsplit("-", 1)[0]
                gene_isoforms[gene_name] = gene_isoforms.get(gene_name, 0) + 1

        cluster_diversity[cluster] = pd.Series(gene_isoforms)
        print(
            f"    {cluster}: Mean isoforms/gene = "
            f"{pd.Series(gene_isoforms.values()).mean():.2f}"
        )

    # Get bulk isoform counts per gene
    if exclude_talon_ids and "annot_transcript_name" in bulk_df.columns:
        bulk_filtered = bulk_df[
            (~bulk_df["annot_transcript_name"].str.startswith("TALON", na=False))
            & (~bulk_df["annot_gene_name"].str.startswith("TALON", na=False))
        ]
    else:
        bulk_filtered = bulk_df

    bulk_isoforms = bulk_filtered.groupby("annot_gene_name")[
        "annot_transcript_name"
    ].nunique()

    # Create comparison plot
    if output_dir:
        output_dir = Path(output_dir)
        fig, ax = plt.subplots(figsize=(10, 6))

        # Prepare data for violin plot
        all_data = []
        all_labels = []

        # Add bulk data
        all_data.append(bulk_isoforms.values)
        all_labels.append("Bulk")

        # Add each cluster
        for cluster_name, diversity in cluster_diversity.items():
            all_data.append(diversity.values)
            all_labels.append(cluster_name)

        # Create violin plot
        parts = ax.violinplot(
            all_data, positions=range(len(all_data)), showmeans=True, showmedians=True
        )

        # Color the violins
        colors = ["#e74c3c"] + (list(plt.cm.Set3(range(len(cluster_diversity)))))
        for i, pc in enumerate(parts["bodies"]):
            pc.set_facecolor(colors[i])
            pc.set_alpha(0.7)

        ax.set_xticks(range(len(all_labels)))
        ax.set_xticklabels(all_labels, rotation=45, ha="right")
        ax.set_ylabel("Isoforms per gene", fontsize=12)
        ax.set_title(
            "Isoform Diversity: Bulk vs Clusters", fontsize=13, fontweight="bold"
        )
        ax.grid(axis="y", alpha=0.3)

        # Add mean values as text
        for i, (data, label) in enumerate(zip(all_data, all_labels)):
            mean_val = np.mean(data)
            ax.text(
                i,
                ax.get_ylim()[1] * 0.95,
                f"μ={mean_val:.2f}",
                ha="center",
                va="top",
                fontsize=9,
                fontweight="bold",
            )

        plt.tight_layout()
        plt.savefig(
            output_dir / "bulk_vs_cluster_isoform_diversity.png",
            dpi=300,
            bbox_inches="tight",
        )
        plt.close()
        print(f"  Saved: {output_dir / 'bulk_vs_cluster_isoform_diversity.png'}")

        # Create summary statistics table
        summary = pd.DataFrame(
            {
                "Dataset": all_labels,
                "Mean_isoforms_per_gene": [np.mean(data) for data in all_data],
                "Median_isoforms_per_gene": [np.median(data) for data in all_data],
                "Max_isoforms": [np.max(data) for data in all_data],
                "Genes_with_multiple_isoforms": [(data > 1).sum() for data in all_data],
            }
        )

        summary.to_csv(
            output_dir / "bulk_vs_cluster_diversity_summary.csv", index=False
        )
        print(f"  Saved: {output_dir / 'bulk_vs_cluster_diversity_summary.csv'}")

        return summary

    return None


def main():
    parser = argparse.ArgumentParser(
        description="Compare transcript diversity between bulk TALON and "
        "single-cell long-read data"
    )
    parser.add_argument(
        "--bulk-talon", required=True, help="Path to bulk TALON TSV file"
    )
    parser.add_argument(
        "--sc-gene-data",
        default="outputs/anndata/combined_long_read_gene.h5ad",
        help="Path to single-cell gene-level h5ad file "
        "(default: outputs/anndata/combined_long_read_gene"
        ".h5ad)",
    )
    parser.add_argument(
        "--sc-transcript-data",
        default="outputs/anndata/combined_long_read_transcript.h5ad",
        help="Path to single-cell transcript-level h5ad file. "
        "Use data/long_transcript.h5ad for original data "
        "or outputs/anndata/"
        "combined_long_read_transcript.h5ad"
        " (default) for combined datasets.",
    )
    parser.add_argument(
        "--output", required=True, help="Output directory for plots and tables"
    )
    parser.add_argument(
        "--min-reads",
        type=int,
        default=5,
        help="Minimum read count for TALON transcripts (default: 5)",
    )
    parser.add_argument(
        "--include-talon-ids",
        action="store_true",
        help="Include TALON IDs in transcript comparisons (default: exclude them)",
    )

    args = parser.parse_args()

    # Determine whether to exclude TALON IDs (opposite of include flag)
    exclude_talon_ids = not args.include_talon_ids

    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {output_dir}")

    # Load data
    print("\n" + "=" * 60)
    print("LOADING DATA")
    print("=" * 60)

    # Load bulk TALON data
    bulk_df = load_talon_data(args.bulk_talon, args.min_reads)
    bulk_genes, bulk_transcripts, bulk_transcript_names = extract_gene_set_from_talon(
        bulk_df, "Bulk TALON", exclude_talon_ids=exclude_talon_ids
    )

    # Load single-cell gene-level data
    print(f"Loading single-cell gene data from {args.sc_gene_data}...")
    sc_gene_adata = get_anndata(args.sc_gene_data)
    sc_genes = extract_gene_set_from_h5ad(sc_gene_adata, "Single-cell (gene-level)")

    # Load transcript-level data
    print(f"Loading sc transcript data from {args.sc_transcript_data}...")
    sc_transcript_adata = get_anndata(args.sc_transcript_data)
    sc_transcripts, sc_transcript_genes, sc_var = extract_transcript_set_from_h5ad(
        sc_transcript_adata, "Single-cell (transcript-level)"
    )

    # Compare gene overlap
    common_genes = compare_gene_overlap(bulk_genes, sc_genes, output_dir)

    # Compare transcript overlap by ID
    print("\nComparing transcript overlap by ID...")
    transcript_overlap = bulk_transcripts & sc_transcripts
    bulk_only_tx = bulk_transcripts - sc_transcripts
    sc_only_tx = sc_transcripts - bulk_transcripts
    print(f"  Transcripts in both datasets: {len(transcript_overlap)}")
    print(f"  Bulk only: {len(bulk_only_tx)}")
    print(f"  Single-cell only: {len(sc_only_tx)}")
    if len(bulk_transcripts | sc_transcripts) > 0:
        total_transcripts = len(bulk_transcripts | sc_transcripts)
        jaccard = len(transcript_overlap) / total_transcripts
        print(f"  Jaccard index: {jaccard:.3f}")

    # Try comparing by transcript names
    print(f"\n  NOTE: SC transcript data loaded from: {args.sc_transcript_data}")
    print(f"  Sample SC transcript IDs: {list(sc_transcripts)[:5]}")
    if bulk_transcript_names:
        print(f"  Sample bulk transcript names: {list(bulk_transcript_names)[:5]}")

    exact_name_matches, gene_name_matches = compare_transcript_names(
        bulk_transcript_names, sc_transcripts, bulk_df
    )

    # Create visualization for transcript name overlap
    if bulk_transcript_names:
        common_transcripts = compare_transcript_overlap(
            bulk_transcript_names, sc_transcripts, output_dir
        )

    # Compare isoform diversity between bulk and SC transcript data
    # Note: SC transcript IDs are in format "GeneName-IsoformNumber"
    # (e.g., "Mrpl15-206")
    # while bulk TALON uses TALON IDs or Ensembl IDs, so we compare at the gene
    # level
    # Analyze transcript novelty
    novelty_col = None
    for col in ["transcript_novelty", "ISM_subtype", "novelty"]:
        if col in bulk_df.columns:
            novelty_col = col
            break

    # Check if SC data also has novelty info
    sc_has_novelty = False
    sc_novelty_col = None
    if len(sc_var.columns) > 0:
        for col in ["transcript_novelty", "ISM_subtype", "novelty"]:
            if col in sc_var.columns:
                sc_novelty_col = col
                sc_has_novelty = True
                break

    # Initialize for summary report
    novelty_counts = None

    if novelty_col and sc_has_novelty:
        print("\n✓ Both datasets have novelty information!")
        print(f"  Bulk: {novelty_col}, SC: {sc_novelty_col}")

        # Create SC dataframe with novelty info for comparison
        sc_novelty_df = sc_var[[sc_novelty_col]].copy()
        # Standardize column name
        sc_novelty_df.columns = ["transcript_novelty"]

        # Filter bulk to common transcripts if needed
        bulk_novelty, sc_novelty = compare_transcript_novelty(
            bulk_df, sc_novelty_df, output_dir
        )
        novelty_counts = bulk_novelty
    elif novelty_col:
        print("\nBulk transcript novelty categories (SC data has no novelty info):")
        novelty_counts = bulk_df[novelty_col].value_counts()
        for category, count in novelty_counts.items():
            print(f"  {category}: {count} ({count / len(bulk_df) * 100:.1f}%)")
    else:
        print("\n✗ No novelty information available in either dataset")

    # Cluster-specific isoform analysis
    print("\n" + "=" * 60)
    print("CLUSTER-SPECIFIC ISOFORM ANALYSIS")
    print("=" * 60)

    # Find cluster-specific isoform switching
    switching_df = find_cluster_specific_isoforms(
        sc_transcript_adata,
        cluster_col="SampleType",
        output_dir=output_dir,
        top_n=5,
        min_expression=5,
    )

    # Summary report
    print("\n" + "=" * 60)
    print("COMPARISON SUMMARY")
    print("=" * 60)
    print("Bulk TALON:")
    print(f"  Genes: {len(bulk_genes)}")
    print(f"  Transcripts: {len(bulk_transcripts)}")
    if novelty_col and novelty_counts is not None:
        print(f"  Novelty categories: {len(novelty_counts)}")
    print("\nSingle-cell:")
    print(f"  Genes detected: {len(sc_genes)}")
    print(f"  Transcripts detected: {len(sc_transcripts)}")
    print("\nOverlap:")
    print(
        f"  Common genes: {len(common_genes)} "
        f"({len(common_genes) / len(bulk_genes) * 100:.1f}% of bulk)"
    )
    print(
        f"  Common transcripts: {len(transcript_overlap)} "
        f"({len(transcript_overlap) / len(bulk_transcripts) * 100:.1f}% of "
        f"bulk)"
    )
    print(f"\nOutput directory: {output_dir}")
    print("=" * 60)


if __name__ == "__main__":
    main()
