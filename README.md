# se4s_isoform_pipeline

A reproducible Python pipeline to bridge short‑read single‑cell RNA‑seq with long‑read isoform resolution in C2C12 muscle cells. It runs QC and clustering on ENCODE scRNA data, maps genes to TALON long‑read annotations, and assigns cluster‑specific “isoform proxy” transcripts.

---

## What this repo contains

- Long‑read utilities built on **TALON 5.0** (`se4s` CLI and `se4s_isoform` library).
- Short‑read **Scanpy** pipeline (QC, clustering, UMAP).
- Integration scripts to map TALON gene symbols and assign isoform proxies per gene and cluster.
- Lightweight **pytest** tests for QC, clustering, mapping, and proxy assignment.

---

## Install

Install the Python package (for the TALON utilities) in a virtual env:

```
python3 -m venv .venv && source .venv/bin/activate
pip install -e .
```

Create or activate a conda/mamba env for the short‑read side:

- Python 3.12
- scanpy, anndata, numpy, pandas, matplotlib, scikit‑learn, pytest

Example:

```
mamba activate se4s_isoform_env
```

---

## Core workflows

### 1. TALON validation and summaries (long‑read)

Validate a TALON TSV and QC log:

```
se4s validate --tsv outputs/tables/bulk_sc_talon_read_annot.tsv --qc outputs/bulk_run_local_QC.log
```

Summarize top isoforms/genes:

```
se4s counts \
  --tsv outputs/tables/bulk_sc_talon_read_annot.tsv \
  --out outputs/tables \
  --dataset ENCFF003OWX \
  --top 50
```

Optional isoform QC table:

```
python - <<'PY'
from se4s_isoform.talon_to_counts import write_isoform_qc_table
write_isoform_qc_table('outputs/tables/bulk_sc_talon_read_annot.tsv','outputs/tables')
PY
```

### 2. Short‑read QC, clustering, mapping, proxies

Download ENCODE `.h5ad` files (once):

```
bash scripts/download_single_cell_data.sh outputs/anndata
```

QC a short‑read dataset:

```
python scripts/utils.py \
  outputs/anndata/short_shallow.h5ad \
  outputs/anndata
# -> outputs/anndata/short_shallow_qc.h5ad
```

Cluster cells (PCA + neighbors + Leiden) and write a UMAP:

```
python scripts/cluster_cells.py \
  outputs/anndata/short_shallow_qc.h5ad \
  outputs/tables/cell_clusters.csv
# -> outputs/umap_clusters.png
```

Map TALON gene symbols to scRNA:

```
python scripts/map_ids_by_symbol.py \
  outputs/tables/bulk_sc_talon_read_annot.tsv \
  outputs/anndata/short_shallow_qc.h5ad \
  outputs/tables/talon_scrna_symbol_map.csv
```

Assign isoform proxies (cluster, gene → top TALON transcript):

```
python scripts/assign_isoform_proxies.py \
  outputs/tables/bulk_sc_talon_read_annot.tsv \
  outputs/anndata/short_shallow_qc.h5ad \
  outputs/tables/cell_clusters.csv \
  outputs/tables/talon_scrna_symbol_map.csv \
  outputs/tables/isoform_proxies.csv
```

Optional: benchmark isoform diversity (TALON vs proxy isoforms per gene):

```
python scripts/benchmark_isoform_diversity.py \
  outputs/tables/bulk_sc_talon_read_annot.tsv \
  outputs/tables/isoform_proxies.csv \
  outputs/tables/isoform_diversity_benchmark.csv
```

---

## One‑shot pipeline

To run QC → clustering → mapping → isoform proxies in one command:

```
python scripts/run_full_pipeline.py \
  /path/to/bulk_sc_talon_read_annot.tsv \
  outputs/anndata/short_shallow_qc.h5ad \
  outputs/
```

Outputs go under `outputs/anndata/`, `outputs/tables/`, and `outputs/umap_clusters.png`.

---

## Tests

Run the unit tests:

```
pytest -v
```

Large TALON‑dependent smoke tests are skipped automatically if the big TSV/log are not present.

---

## Notes

- Large files and generated outputs under `outputs/` are not tracked by Git; re‑run the pipeline to regenerate them.
- See `notebooks/scRNA_QC_EDA.ipynb` for exploratory QC and clustering decisions.
- TALON itself and ENCODE inputs remain under their original licenses and usage policies.
```
