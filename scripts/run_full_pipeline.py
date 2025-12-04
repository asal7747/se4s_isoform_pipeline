#!/usr/bin/env python
"""
Run the full isoform proxy pipeline on a TALON TSV.

Usage:
  python scripts/run_full_pipeline.py <talon_tsv> <outdir>
  python scripts/run_full_pipeline.py <talon_tsv> <scrna_h5ad> <outdir>

With 2 arguments:
  - Downloads ENCODE .h5ad files if not present
  - Combines short+long read datasets automatically
  - Runs full pipeline including bulk vs single-cell comparison

With 3 arguments:
  - Uses your provided scRNA .h5ad file directly
  - Runs QC on it, then clustering, mapping, and proxy assignment
  - Skips the download/combine and comparison steps
"""

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent
SCRIPTS = REPO_ROOT / "scripts"


def run(cmd: list[str]) -> None:
    print(f"[run_full_pipeline] Running: {' '.join(cmd)}")
    result = subprocess.run(cmd)
    if result.returncode != 0:
        print(
            f"[run_full_pipeline] ERROR: command failed with code {result.returncode}"
        )
        sys.exit(result.returncode)


def main():
    if len(sys.argv) < 3:
        print(
            "Usage:\n"
            "  run_full_pipeline.py <talon_tsv> <outdir>\n"
            "  run_full_pipeline.py <talon_tsv> <scrna_h5ad> <outdir>",
            file=sys.stderr,
        )
        sys.exit(2)

    talon_tsv = Path(sys.argv[1]).resolve()

    # Determine mode: 2 args (auto) or 3 args (user-provided h5ad)
    if len(sys.argv) == 3:
        # 2-arg mode: auto-download and combine
        scrna_h5ad = None
        outdir = Path(sys.argv[2]).resolve()
        auto_mode = True
    else:
        # 3-arg mode: user provides h5ad
        scrna_h5ad = Path(sys.argv[2]).resolve()
        outdir = Path(sys.argv[3]).resolve()
        auto_mode = False

    # Basic input checks
    if not talon_tsv.is_file():
        print(f"ERROR: TALON TSV not found: {talon_tsv}", file=sys.stderr)
        sys.exit(1)
    if scrna_h5ad is not None and not scrna_h5ad.is_file():
        print(f"ERROR: scRNA h5ad not found: {scrna_h5ad}", file=sys.stderr)
        sys.exit(1)

    # Create output structure
    (outdir / "tables").mkdir(parents=True, exist_ok=True)
    (outdir / "plots").mkdir(parents=True, exist_ok=True)

    if auto_mode:
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
        missing_files = [f for f in required_files if not (data_dir / f).is_file()]

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
            print(f"ERROR: Combined QC file not created: {qc_h5ad}", file=sys.stderr)
            sys.exit(1)
    else:
        # 3-arg mode: QC the user-provided h5ad
        print(f"[run_full_pipeline] Running QC on {scrna_h5ad}...")
        run(
            [
                sys.executable,
                str(SCRIPTS / "utils.py"),
                str(scrna_h5ad),
                str(outdir),
            ]
        )
        # utils.py writes <stem>_qc.h5ad; find it
        qc_candidates = list(outdir.glob(f"{scrna_h5ad.stem}_qc.h5ad"))
        if qc_candidates:
            qc_h5ad = qc_candidates[0]
        else:
            qc_h5ad = outdir / f"{scrna_h5ad.stem}_qc.h5ad"
        if not qc_h5ad.is_file():
            print(f"ERROR: QC output not found: {qc_h5ad}", file=sys.stderr)
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

    # 5) Compare bulk and single-cell isoforms (only in auto mode)
    if auto_mode:
        print("[run_full_pipeline] Long-read datasets already combined.")
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
