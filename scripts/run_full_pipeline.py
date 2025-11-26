#!/usr/bin/env python
"""
Run the full isoform proxy pipeline on a TALON TSV + scRNA h5ad.

Usage:
  python scripts/run_full_pipeline.py <talon_tsv> <scrna_h5ad> <outdir>

This will:
  1) QC the scRNA data -> <outdir>/scrna_qc.h5ad
  2) Cluster cells and save clusters + UMAP
  3) Map TALON genes to scRNA symbols
  4) Assign cluster-specific isoform proxies
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
        print(f"[run_full_pipeline] ERROR: command failed with code {result.returncode}")
        sys.exit(result.returncode)


def main():
    if len(sys.argv) < 4:
        print(
            "Usage: run_full_pipeline.py <talon_tsv> <scrna_h5ad> <outdir>",
            file=sys.stderr,
        )
        sys.exit(2)

    talon_tsv = Path(sys.argv[1]).resolve()
    scrna_h5ad = Path(sys.argv[2]).resolve()
    outdir = Path(sys.argv[3]).resolve()

    # Basic input checks
    if not talon_tsv.is_file():
        print(f"ERROR: TALON TSV not found: {talon_tsv}", file=sys.stderr)
        sys.exit(1)
    if not scrna_h5ad.is_file():
        print(f"ERROR: scRNA h5ad not found: {scrna_h5ad}", file=sys.stderr)
        sys.exit(1)

    # Create output structure
    (outdir / "tables").mkdir(parents=True, exist_ok=True)
    (outdir / "plots").mkdir(parents=True, exist_ok=True)

    # 1) QC scRNA
    qc_h5ad = outdir / "scrna_qc.h5ad"
    run(
        [
            sys.executable,
            str(SCRIPTS / "utils.py"),
            str(scrna_h5ad),
            str(outdir),
        ]
    )
    # utils.py writes <stem>_qc.h5ad; infer that path
    qc_candidates = list(outdir.glob(f"{scrna_h5ad.stem}_qc.h5ad"))
    if qc_candidates:
        qc_h5ad = qc_candidates[0]
    if not qc_h5ad.is_file():
        print(f"ERROR: QC output not found in {outdir}", file=sys.stderr)
        sys.exit(1)

    # 2) Cluster cells
    clusters_csv = outdir / "tables" / "cell_clusters.csv"
    run(
        [
            sys.executable,
            str(SCRIPTS / "cluster_cells.py"),
            str(qc_h5ad),
            str(clusters_csv),
        ]
    )
    # cluster_cells.py writes outputs/umap_clusters.png relative to repo root;
    # optionally copy it into outdir/plots for convenience.
    umap_src = REPO_ROOT / "outputs" / "umap_clusters.png"
    umap_dst = outdir / "plots" / "umap_clusters.png"
    if umap_src.is_file():
        umap_dst.write_bytes(umap_src.read_bytes())

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

    print(f"[run_full_pipeline] Done. Outputs under {outdir}")


if __name__ == "__main__":
    main()
