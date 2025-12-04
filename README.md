# se4s_isoform_pipeline

A reproducible Python pipeline to bridge short-read single-cell RNA-seq with long-read isoform resolution, demonstrated on ENCODE C2C12 muscle cells. The pipeline runs QC and clustering on ENCODE scRNA data, maps genes to TALON long-read annotations, and assigns cluster-level "isoform proxy" transcripts. It can be run on any dataset where you have a TALON long-read catalog and a compatible scRNA .h5ad.

---

## What this repo contains

- Long-read utilities built on TALON 5.0 (`se4s` CLI and `se4s_isoform` library) for validation and counts.
- Short-read Scanpy pipeline (QC, clustering, UMAP) implemented as reusable functions (`cluster_cells.py`) and a CLI wrapper (`run_cluster_cells.py`).
- Integration scripts to map TALON gene symbols and assign isoform proxies per gene and cluster, plus an optional isoform-diversity benchmark.
- A Snakemake workflow (`Snakefile`) and Python driver (`run_full_pipeline.py`) that encode the end-to-end DAG: QC -> clustering -> symbol mapping -> proxies.
- Lightweight pytest plus functional shell tests for QC, clustering, mapping, proxies, and the full pipeline.

---

## Requirements

- **Python 3.10, 3.11, 3.12, or 3.13** (Python 3.14+ is not yet supported due to numba/scanpy compatibility)

Check your Python version:

```
python3 --version
```

If you have Python 3.14+, you'll need to create the virtual environment with an older Python. Tools like [pyenv](https://github.com/pyenv/pyenv) (macOS/Linux) or the [Python installer](https://www.python.org/downloads/) (Windows) can install multiple Python versions side-by-side. Then create the venv with the specific version:

```
python3.13 -m venv .venv && source .venv/bin/activate
```

---

## Install

Create and activate a virtual environment, then install the package and its dependencies:

```
python3 -m venv .venv && source .venv/bin/activate
pip install -e .
pip install scanpy anndata numpy pandas matplotlib scikit-learn pytest python-igraph leidenalg
```

Optionally, install Snakemake if you want to use the Snakefile workflow:

```
pip install snakemake
```

---

## Core workflows

### 1. TALON validation and summaries (long-read)

Validate a TALON TSV and QC log:

```
se4s validate --tsv /path/to/bulk_sc_talon_read_annot.tsv --qc /path/to/bulk_run_local_QC.log
```

Summarize top isoforms/genes:

```
se4s counts --tsv /path/to/bulk_sc_talon_read_annot.tsv \
  --out outputs/tables \
  --dataset ENCFF003OWX \
  --top 50
```

Optional isoform QC table:

```
from se4s_isoform.talon_to_counts import write_isoform_qc_table
write_isoform_qc_table("/path/to/bulk_sc_talon_read_annot.tsv", "outputs/tables")
```

### 2. Short-read QC, clustering, mapping, proxies

These individual scripts are useful for custom workflows. For most users, we recommend using the **one-shot pipeline** below instead.

Download ENCODE .h5ad files (once):

```
bash scripts/download_single_cell_data.sh data
```

The remaining steps require combined data with proper gene symbols. See the one-shot pipeline section for the recommended approach.

---

## One-shot pipeline (Python driver)

**This is the recommended way to run the pipeline.**

### Option 1: Auto-download and combine (2 arguments)

To run the full pipeline with automatic data download:

```
python scripts/run_full_pipeline.py /path/to/bulk_sc_talon_read_annot.tsv outputs/
```

This will:
- Download ENCODE .h5ad files to `data/` if not present
- Combine short+long read datasets
- Run QC, clustering, mapping, and proxy assignment
- Run bulk vs single-cell comparison

### Option 2: Use your own h5ad (3 arguments)

To run with a specific scRNA .h5ad file you already have:

```
python scripts/run_full_pipeline.py /path/to/bulk_sc_talon_read_annot.tsv /path/to/your_scrna.h5ad outputs/
```

This will:
- Run QC on your provided h5ad
- Run clustering, mapping, and proxy assignment
- Skip download/combine and comparison steps

Replace the paths with your actual file locations.

Outputs go under `outputs/`:

- `outputs/anndata/combined_short_read_qc.h5ad` (QC'd combined short-read data)
- `outputs/tables/cell_clusters.csv`
- `outputs/tables/talon_scrna_symbol_map.csv`
- `outputs/tables/isoform_proxies.csv`
- `outputs/plots/umap_clusters.png`
- `outputs/comparison/` (bulk vs single-cell comparison results)

---

## Snakemake workflow (optional)

The same DAG is encoded in `Snakefile`. **Unlike the Python driver, Snakemake does NOT auto-download data.** You must first download ENCODE files:

```
bash scripts/download_single_cell_data.sh data
```

Then run with Snakemake installed:

```
snakemake -j 1 --config talon_tsv=/path/to/bulk_sc_talon_read_annot.tsv
```

Replace `/path/to/bulk_sc_talon_read_annot.tsv` with the actual path to your TALON TSV file.

The workflow will:
1. Combine short+long read datasets via `single_cell_analysis.py`
2. Cluster cells and generate UMAP
3. Map TALON gene symbols to scRNA
4. Assign isoform proxies
5. Compare bulk vs single-cell isoforms

Outputs are the same as `run_full_pipeline.py`.

---

## Tests

Run the unit tests:

```
pytest -v
```

Functional shell tests (using `ssshtest`) live under `tests/func/` and exercise:

- `scripts/run_cluster_cells.py` on a small test h5ad.
- `scripts/run_full_pipeline.py` on small TALON + scRNA test data.

Large TALON-dependent smoke tests are skipped automatically if the big TSV/log are not present.

---

## Troubleshooting

**Error: "Cannot install on Python version 3.14.0"**
- Numba (required by scanpy) doesn't support Python 3.14 yet. Create the virtual environment with Python 3.13 or earlier (see Requirements section).

**Error: "Leiden clustering requires `python-igraph`"**
- Install missing clustering dependencies:
  ```
  pip install python-igraph leidenalg
  ```

---

## Notes

- Large files and generated outputs under `outputs/` are not tracked by Git; re-run the pipeline or Snakemake to regenerate them.
- On full ENCODE plus large TALON TSV, some steps may currently stop with explicit, data-format-related errors (for example, no overlapping gene symbols); this is by design rather than failing silently.
- See `notebooks/scRNA_QC_EDA.ipynb` for exploratory QC and clustering decisions.
- TALON itself and ENCODE inputs remain under their original licenses and usage policies.
