#!/usr/bin/env python
"""
Run the full isoform proxy pipeline on a TALON TSV.

Usage:
  python scripts/run_full_pipeline.py <talon_tsv> <outdir>

This will:
  0) Download single-cell data files if not present in data/ directory
  1) Combine short+long read datasets from data/ directory
  2) QC the scRNA data
  3) Cluster cells and save clusters + UMAP
  4) Map TALON genes to scRNA symbols
  5) Assign cluster-specific isoform proxies
  6) Compare bulk vs single-cell isoforms
"""

import os
import sys
from pathlib import Path

import subprocess


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent
SCRIPTS = REPO_ROOT / "scripts"


def run(cmd: list[str]) -> None:
    print(f"[run_full_pipeline] Running: {' '.join(cmd)}")
    result = subprocess.run(cmd)
    if result.returncode != 0:
        print(
            f"[run_full_pipeline] ERROR: command failed with "
            f"code {result.returncode}"
        )
        sys.exit(result.returncode)


def main():
    if len(sys.argv) < 3:
        print(
            "Usage: run_full_pipeline.py <talon_tsv> <outdir>",
            file=sys.stderr,
        )
        sys.exit(2)

    talon_tsv = Path(sys.argv[1]).resolve()
    outdir = Path(sys.argv[2]).resolve()

    # Basic input checks
    if not talon_tsv.is_file():
        print(f"ERROR: TALON TSV not found: {talon_tsv}", file=sys.stderr)
        sys.exit(1)

    # Create output structure
    (outdir / "tables").mkdir(parents=True, exist_ok=True)
    (outdir / "plots").mkdir(parents=True, exist_ok=True)

    # 0) Download single-cell data if needed
    data_dir = REPO_ROOT / "data"
    required_files = [
        "short_shallow.h5ad",
        "short_shallow_nuc.h5ad",
        "short_deep.h5ad",
        "short_deep_nuc.h5ad",
        "long_transcript.h5ad",
        "long_gene.h5ad",
        "long_nuc_transcript.h5ad",
        "long_nuc_gene.h5ad",
        "short_shallow_myotube.h5ad",
        "short_deep_myotube.h5ad",
        "long_myotube_transcript.h5ad",
        "long_myotube_gene.h5ad",
    ]
    
    missing_files = [f for f in required_files 
                     if not (data_dir / f).is_file()]
    
    if missing_files:
        print(
            f"[run_full_pipeline] Missing {len(missing_files)} data files. "
            "Downloading..."
        )
        download_script = SCRIPTS / "download_single_cell_data.sh"
        run(["bash", str(download_script), str(data_dir)])
    else:
        print("[run_full_pipeline] All data files present, skipping download.")

    # 1) Combine short+long read datasets
    print("[run_full_pipeline] Combining short+long read datasets...")
    run(
        [
            sys.executable,
            str(SCRIPTS / "single_cell_analysis.py"),
        ]
    )
    # This creates outputs/anndata/combined_short_read_qc.h5ad
    qc_h5ad = REPO_ROOT / "outputs" / "anndata" / "combined_short_read_qc.h5ad"
    if not qc_h5ad.is_file():
        print(
            f"ERROR: Combined QC file not created: {qc_h5ad}",
            file=sys.stderr
        )
        sys.exit(1)    
    
    # 2) Cluster cells
    clusters_csv = outdir / "tables" / "cell_clusters.csv"
    plots_dir = outdir / "plots"
    run(
        [
            sys.executable,
            str(SCRIPTS / "run_cluster_cells.py"),
            str(qc_h5ad),
            str(clusters_csv),
            "--output-dir",
            str(plots_dir),
        ]
    )
    # run_cluster_cells.py will write umap_clusters.png to the plots directory

    # 3) Map IDs by symbol
    symbol_map_csv = outdir / "tables" / "talon_scrna_symbol_map.csv"
    run(
        [
            sys.executable,
            str(SCRIPTS / "map_ids_by_symbol.py"),
            str(talon_tsv),
            str(qc_h5ad),
            str(symbol_map_csv),
        ]
    )

    # 4) Assign isoform proxies
    proxies_csv = outdir / "tables" / "isoform_proxies.csv"
    run(
        [
            sys.executable,
            str(SCRIPTS / "assign_isoform_proxies.py"),
            str(talon_tsv),
            str(qc_h5ad),
            str(clusters_csv),
            str(symbol_map_csv),
            str(proxies_csv),
        ]
    )

    # 5) Long-read data already combined by step 1
    # (single_cell_analysis.py created combined_long_read_transcript.h5ad)
    print("[run_full_pipeline] Long-read datasets already combined.")

    # 6) Compare bulk and single-cell isoforms
    comparison_dir = outdir / "comparison"
    comparison_dir.mkdir(parents=True, exist_ok=True)
    run(
        [
            sys.executable,
            str(SCRIPTS / "compare_bulk_sc.py"),
            "--bulk-talon",
            str(talon_tsv),
            "--output",
            str(comparison_dir),
        ]
    )
    # This creates comparison plots including cluster_isoform_switching.png

    print(f"[run_full_pipeline] Done. Outputs under {outdir}")


if __name__ == "__main__":
    main()
