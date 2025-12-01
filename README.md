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

## Install

Install the core Python package (TALON utilities) in a virtual environment:

```
python3 -m venv .venv && source .venv/bin/activate
pip install -e .
```

Create or activate a conda/mamba environment for the short-read pipeline (example):

- Python 3.12
- scanpy, anndata, numpy, pandas, matplotlib, scikit-learn, pytest
- snakemake (optional, for the Snakefile workflow)

```
mamba activate se4s_isoform_env
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

Download ENCODE .h5ad files (once):

```
bash scripts/download_single_cell_data.sh outputs/anndata
```

QC a short-read dataset:

```
python scripts/utils.py outputs/anndata/short_shallow.h5ad outputs
# -> outputs/short_shallow_qc.h5ad
```

Cluster cells (PCA + neighbors + Leiden) and write a UMAP:

```
python scripts/run_cluster_cells.py outputs/short_shallow_qc.h5ad \
  outputs/tables/cell_clusters.csv --output-dir outputs/plots
# -> outputs/plots/umap_clusters.png
```

Map TALON gene symbols to scRNA:

```
python scripts/map_ids_by_symbol.py /path/to/bulk_sc_talon_read_annot.tsv \
  outputs/short_shallow_qc.h5ad outputs/tables/talon_scrna_symbol_map.csv
```

Assign isoform proxies (cluster, gene -> top TALON transcript):

```
python scripts/assign_isoform_proxies.py /path/to/bulk_sc_talon_read_annot.tsv \
  outputs/short_shallow_qc.h5ad outputs/tables/cell_clusters.csv \
  outputs/tables/talon_scrna_symbol_map.csv outputs/tables/isoform_proxies.csv
```

Optional: benchmark isoform diversity (TALON vs proxy isoforms per gene):

```
python scripts/benchmark_isoform_diversity.py /path/to/bulk_sc_talon_read_annot.tsv \
  outputs/tables/isoform_proxies.csv outputs/tables/isoform_diversity_benchmark.csv
```

---

## One-shot pipeline (Python driver)

To run QC -> clustering -> mapping -> isoform proxies in one command:

```
python scripts/run_full_pipeline.py /path/to/bulk_sc_talon_read_annot.tsv \
  outputs/anndata/short_shallow.h5ad outputs/
```

Outputs go under `outputs/`:

- `outputs/short_shallow_qc.h5ad`
- `outputs/tables/` (clusters, symbol map, proxies, summaries)
- `outputs/plots/umap_clusters.png`

---

## Snakemake workflow (optional)

The same DAG is encoded in `Snakefile`. With Snakemake installed in your environment:

```
snakemake -j 1
```

By default, this will:

- Take `outputs/anndata/short_shallow.h5ad` and `/path/to/bulk_sc_talon_read_annot.tsv` as inputs.
- Run:
  - `scripts/utils.py` (QC)
  - `scripts/run_cluster_cells.py` (clustering + UMAP)
  - `scripts/map_ids_by_symbol.py` (symbol mapping)
  - `scripts/assign_isoform_proxies.py` (proxy assignment)
- Produce the same outputs under `outputs/` as `run_full_pipeline.py`.

You can edit the paths at the top of `Snakefile` to point to different TALON/scRNA inputs.

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

## Notes

- Large files and generated outputs under `outputs/` are not tracked by Git; re-run the pipeline or Snakemake to regenerate them.
- On full ENCODE plus large TALON TSV, some steps may currently stop with explicit, data-format-related errors (for example, no overlapping gene symbols); this is by design rather than failing silently.
- See `notebooks/scRNA_QC_EDA.ipynb` for exploratory QC and clustering decisions.
- TALON itself and ENCODE inputs remain under their original licenses and usage policies.
