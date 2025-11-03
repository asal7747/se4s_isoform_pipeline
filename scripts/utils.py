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
import scanpy as sc

def load_and_qc_h5ad(
    in_h5ad: str,
    out_dir: str,
    min_counts: int = 1000,
    min_genes: int = 750,
    mt_prefix: str = "mt-",
) -> str:
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    adata = sc.read(in_h5ad)

    # Standardize obs_names if present
    if "Unnamed: 0" in adata.obs.columns:
        adata.obs_names = adata.obs["Unnamed: 0"].astype(str)
        adata.obs.drop(columns=["Unnamed: 0"], inplace=True)

    # Standardize var_names if present
    if "x" in adata.var.columns:
        adata.var_names = adata.var["x"].astype(str)
        adata.var.rename(columns={"x": "transcript_id"}, inplace=True)

    # QC metrics with mitochondrial flag
    adata.var["mt"] = adata.var_names.str.lower().str.startswith(mt_prefix)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    # Filters
    sc.pp.filter_cells(adata, min_counts=min_counts)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=2)
    sc.pp.filter_genes(adata, min_counts=5)

    # Normalize and log1p
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    out_file = out_path / f"{Path(in_h5ad).stem}_qc.h5ad"
    adata.write(out_file)
    print(f"Wrote QC'd h5ad: {out_file} ({adata.n_obs} cells × {adata.n_vars} genes)")
    return str(out_file)

def main():
    if len(sys.argv) < 3:
        print("Usage: utils.py <input.h5ad> <output_dir> [min_counts=1000] [min_genes=750]")
        sys.exit(2)
    in_h5ad = sys.argv[1]
    out_dir = sys.argv[2]
    min_counts = int(sys.argv[3]) if len(sys.argv) > 3 else 1000
    min_genes  = int(sys.argv[4]) if len(sys.argv) > 4 else 750
    load_and_qc_h5ad(in_h5ad, out_dir, min_counts=min_counts, min_genes=min_genes)

if __name__ == "__main__":
    main()
