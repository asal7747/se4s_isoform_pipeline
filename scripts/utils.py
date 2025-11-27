#!/usr/bin/env python
"""
Minimal utilities for loading and QC'ing ENCODE scRNA h5ad files.

Preserves the essential steps from the exploratory script:
- read h5ad
- standardize obs_names (if Unnamed: 0 exists) and var_names (if column 'x' exists)
- compute QC metrics (including mitochondrial content)
- filter cells and genes
- normalize total counts and log1p
- write a QC'd h5ad next to the input (with _qc suffix)

Usage:
  python scripts/utils.py <input.h5ad> <output_dir> [min_counts=1000] [min_genes=750]
"""

import sys
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

import scanpy as sc
import anndata
from anndata import AnnData

import cluster_cells as cc # for PCA, UMAP, t-SNE functions


def get_anndata(in_h5ad: str) -> sc.AnnData:
    in_path = Path(in_h5ad)
    if not in_path.is_file():
        raise FileNotFoundError(f"Input .h5ad not found: {in_path}")
    adata = sc.read(in_path)
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError(
            f"Input AnnData from {in_path} is empty "
            f"({adata.n_obs} cells × {adata.n_vars} genes)."
        )
    return adata


def load_short_read_datasets(data_dir: str = "outputs/anndata") -> dict[str, sc.AnnData]:
    """
    Load all short-read h5ad files from a directory.
    
    Searches for files starting with 'short' and ending with '.h5ad',
    loads each using get_anndata, and returns them in a dictionary
    keyed by the filename (without extension).
    
    Args:
        data_dir: Directory containing .h5ad files (default: "outputs/anndata")
        
    Returns:
        Dictionary mapping dataset names to AnnData objects
        
    Raises:
        FileNotFoundError: If data_dir doesn't exist
    """
    data_path = Path(data_dir)
    if not data_path.is_dir():
        raise FileNotFoundError(f"Data directory not found: {data_path}")
    
    datasets = {}
    for fpath in data_path.glob("short*.h5ad"):
        name = fpath.stem  # filename without extension
        adata = get_anndata(str(fpath))
        datasets[name] = adata
    
    return datasets


def short_read_adata_to_df(
    adata: sc.AnnData,
    make_unique: bool = True
) -> pd.DataFrame:
    """
    Convert a short-read AnnData object to a pandas DataFrame.
    
    Uses obs_names as the index (cells/observations) and var_names 
    as columns (genes/variables). Converts sparse matrices to dense
    automatically. Optionally makes names unique to avoid pandas errors.
    
    Args:
        adata: Input AnnData object
        make_unique: If True, make obs_names and var_names unique
                     by appending suffixes (default: True)
        
    Returns:
        DataFrame with cells as rows and genes as columns
        
    Raises:
        ValueError: If obs_names or var_names length doesn't match matrix shape
    """
    # Convert sparse to dense if needed
    data = adata.X.todense() if hasattr(adata.X, "todense") else adata.X

    n_obs, n_vars = data.shape
    if len(adata.obs_names) != n_obs:
        raise ValueError(
            f"adata.obs_names length ({len(adata.obs_names)}) != n_obs ({n_obs})"
        )
    if len(adata.var_names) != n_vars:
        raise ValueError(
            f"adata.var_names length ({len(adata.var_names)}) != n_vars ({n_vars})"
        )

    # Make names unique if requested
    if make_unique:
        if not adata.var_names.is_unique:
            adata.var_names_make_unique(join="-")
        if not adata.obs_names.is_unique:
            adata.obs_names_make_unique(join="-")

    return pd.DataFrame(data, index=adata.obs_names, columns=adata.var_names)


def short_read_adata_index(adata: sc.AnnData) -> sc.AnnData:
    """
    Standardize obs_names and var_names for short-read AnnData.
    
    Renames 'Unnamed: 0' column to 'cell_id' and sets it as obs_names.
    Renames 'x' column to 'transcript_id' and sets it as var_names.
    
    Args:
        adata: Input AnnData object with 'Unnamed: 0' in obs and 'x' in var
        
    Returns:
        AnnData with standardized obs_names and var_names
        
    Raises:
        KeyError: If 'Unnamed: 0' or 'x' columns are missing
    """
    adata.obs.rename({"Unnamed: 0": "cell_id"}, axis=1, inplace=True)
    adata.var.rename({"x": "transcript_id"}, axis=1, inplace=True)
    adata.obs_names = adata.obs["cell_id"].values
    adata.var_names = adata.var["transcript_id"].values
    return adata


def load_long_read_datasets(data_dir: str = "outputs/anndata") -> dict[str, sc.AnnData]:
    """
    Load all long-read h5ad files from a directory.
    
    Searches for files starting with 'long' and ending with '.h5ad',
    loads each using get_anndata, and returns them in a dictionary
    keyed by the filename (without extension).
    
    Args:
        data_dir: Directory containing .h5ad files (default: "outputs/anndata")
        
    Returns:
        Dictionary mapping dataset names to AnnData objects
        
    Raises:
        FileNotFoundError: If data_dir doesn't exist
    """
    data_path = Path(data_dir)
    if not data_path.is_dir():
        raise FileNotFoundError(f"Data directory not found: {data_path}")
    
    datasets = {}
    for fpath in data_path.glob("long*.h5ad"):
        name = fpath.stem  # filename without extension
        adata = get_anndata(str(fpath))
        datasets[name] = adata
    
    return datasets


def qc_and_filter(
    adata: sc.AnnData,
    min_counts: int = 1000,
    min_genes: int = 750,
    mt_prefix: str = "mt-",
    max_genes: int = 200000,
    max_mt_pct: float = 20.0,
    output_dir: str = "outputs/anndata"
) -> sc.AnnData:
    """
    Perform quality control, filtering, and normalization on an AnnData object.
    
    Steps:
    1. Calculate QC metrics (including mitochondrial content)
    2. Filter cells by min_counts and min_genes
    3. Filter genes by min_cells and min_counts
    4. Normalize total counts to 10,000 per cell
    5. Log1p transform
    
    Args:
        adata: Input AnnData object
        min_counts: Minimum total counts per cell
        min_genes: Minimum genes detected per cell
        mt_prefix: Prefix for mitochondrial genes (case-insensitive)
        
    Returns:
        Filtered and normalized AnnData object
        
    Raises:
        ValueError: If min_counts or min_genes are invalid
        RuntimeError: If all cells/genes are filtered out
    """
    if min_counts <= 0 or min_genes <= 0:
        raise ValueError(
            f"min_counts and min_genes must be positive; "
            f"got min_counts={min_counts}, min_genes={min_genes}"
        )
    
    # print highly expressed genes
    print(sc.pl.highest_expr_genes(adata, n_top=20))

    # QC metrics with mitochondrial flag
    adata.var["mt"] = adata.var_names.str.lower().str.startswith(mt_prefix)
    qc = sc.pp.calculate_qc_metrics(adata)
    cell_qc_dataframe = qc[0]
    gene_qc_dataframe = qc[1]

    # plot histograms of QC metrics to validate parameters are reasonable
    # add red line showing min_counts and min_genes thresholds
    plt.figure(figsize=(12, 8))
    # total counts per cell
    plt.subplot(2, 1, 1)
    plt.hist(cell_qc_dataframe['total_counts'], bins=1000, color='skyblue')
    plt.xlim(0,50000)
    plt.ylim(0,2000)
    plt.axvline(min_counts, color='red', linestyle='dashed', linewidth=1)
    plt.title('Total Counts per Cell')
    plt.xlabel('Total Counts')
    plt.ylabel('Number of Cells')
    # number of genes detected per cell
    plt.subplot(2, 1, 2)
    plt.hist(cell_qc_dataframe['n_genes_by_counts'], bins=1000, color='salmon')
    plt.axvline(min_genes, color='red', linestyle='dashed', linewidth=1)
    plt.title('Number of Genes Detected per Cell')
    plt.xlabel('Number of Genes')
    plt.ylabel('Number of Cells')
    plt.tight_layout()

    # Filters
    # remove cells with fewer than N number of reads
    print('Before cell filters: \n', adata)
    sc.pp.filter_cells(adata, min_counts=min_counts)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    print('After cell filters: \n', adata)

    # next, remove genes whose expression level is considered "undetectable"
    # we define a gene as detectable if at least two cells contain more than 5 reads from the gene
    # but the threshold strongly depends on the sequencing depth
    # filter genes
    # e.g. set the threshold to be at least 2 cells with at least 5 reads

    print('Before gene filters: \n', adata)
    sc.pp.filter_genes(adata, min_cells=2)
    sc.pp.filter_genes(adata, min_counts=5)
    print('After gene filters: \n', adata)

    # find transcript with MT prefix
    mt_transcripts = adata.var_names.str.lower().str.contains(mt_prefix)
    # make mt- column in var
    adata.var['mt'] = mt_transcripts

    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # make a violin plot of some of the computed quality measures:
    # the number of genes expressed in the count matrix
    # the total counts per cell
    # the percentage of counts in mitochondrial genes

    sc.pl.violin(adata,
                 ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                 jitter=0.4, multi_panel=True)

    # or, visualize mitochondrial counts as scatterplots 
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')

    # additional hard filters to remove potential outliers (these thresholds can be adjusted)
    adata = adata[adata.obs.n_genes_by_counts < max_genes, :]
    adata = adata[adata.obs.pct_counts_mt < max_mt_pct, :]

    # print final shape after filtering
    print('After additional hard filters: \n', adata)

    # save the QC'ed anndata object
    ica_path = Path(output_dir)
    ica_path.mkdir(parents=True, exist_ok=True)
    
    # Use dataset_name from uns if available, otherwise use a default name
    dataset_name = adata.uns.get('dataset_name', 'combined_dataset')
    qc_h5ad_path = ica_path / f"{dataset_name}_qc.h5ad"
    adata.write(qc_h5ad_path)
    print(f"QC'ed AnnData saved to: {qc_h5ad_path}")

    # Guard against everything being filtered out and warn on tiny outputs
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise RuntimeError(
            "All cells or genes were filtered out during QC "
            f"(n_obs={adata.n_obs}, n_vars={adata.n_vars}). "
            "Consider relaxing min_counts/min_genes or gene filters."
        )

    if adata.n_obs < 200 or adata.n_vars < 200:
        print(
            f"WARNING: QC left a very small dataset "
            f"({adata.n_obs} cells × {adata.n_vars} genes); "
            "downstream clustering may be unstable.",
            file=sys.stderr,
        )

    # Normalize and log1p
    # sc.pp.normalize_total(adata, target_sum=1e4)
    # sc.pp.log1p(adata)
    
    return adata


def load_and_combine_short_reads(
    data_dir: str = "data",
    verbose: bool = True
) -> tuple[sc.AnnData, pd.DataFrame]:
    """
    Load, standardize, and combine all short-read datasets.
    
    This function:
    1. Loads all short-read h5ad files from the directory
    2. Standardizes obs_names and var_names for each dataset
    3. Concatenates all datasets into a single AnnData
    4. Converts to a pandas DataFrame
    
    Args:
        data_dir: Directory containing short-read .h5ad files
        verbose: If True, print dataset information and shapes
        
    Returns:
        Tuple of (combined_adata, combined_dataframe)
    """
    # Load all short read datasets
    short_read_datasets = load_short_read_datasets(data_dir)
    
    if verbose:
        for name, adata in short_read_datasets.items():
            print(f'Dataset: {name}')
            print(adata)
            print('\n')
    
    # Standardize indexes for all datasets
    for name, adata in short_read_datasets.items():
        short_read_datasets[name] = short_read_adata_index(adata)
    
    # Concatenate into single AnnData
    combined_adata = anndata.concat(list(short_read_datasets.values()), axis=0)
    
    # Add dataset name to uns for later use (e.g., in QC output filenames)
    combined_adata.uns['dataset_name'] = 'combined_short_read'
    
    if verbose:
        print(f'Combined short read anndata shape: {combined_adata.shape}')
    
    # Convert to DataFrame
    combined_df = short_read_adata_to_df(combined_adata)
    if verbose:
        print(f'Combined short read dataframe shape: {combined_df.shape}')
    
    return combined_adata, combined_df


def main():
    if len(sys.argv) < 3:
        print(
            "Usage: utils.py <input.h5ad> <output_dir> [min_counts=1000] [min_genes=750]",
            file=sys.stderr,
        )
        sys.exit(2)

    in_h5ad = sys.argv[1]
    out_dir = sys.argv[2]
    try:
        min_counts = int(sys.argv[3]) if len(sys.argv) > 3 else 1000
        min_genes = int(sys.argv[4]) if len(sys.argv) > 4 else 750
    except ValueError:
        print("min_counts and min_genes must be integers.", file=sys.stderr)
        sys.exit(2)
    # try loading short read data and doing qc and filtering
    try:
        combined_short_read_adata, combined_short_read_df = load_and_combine_short_reads()
        qc_combined_short_read_adata = qc_and_filter(
            combined_short_read_adata,
            min_counts=min_counts,
            min_genes=min_genes,
            output_dir=out_dir,
        )
        print("QC and filtering completed successfully.")
        print(f"QC'ed AnnData shape: {qc_combined_short_read_adata.shape}")

    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
