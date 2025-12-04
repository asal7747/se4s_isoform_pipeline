from pathlib import Path
import os

# --------------------------------------------------------------------
# CONFIG
# --------------------------------------------------------------------

# Input files can be provided via Snakemake `config` (e.g. --config talon_tsv=...) or
# via environment variables (TALON_TSV, SCRNA_H5AD). If neither is provided we
# fall back to sensible repository-relative defaults so the Snakefile is
# functional out-of-the-box.
TALON_TSV = (
    config.get("talon_tsv")
    if "config" in globals() and "talon_tsv" in config
    else os.environ.get(
        "TALON_TSV", "data/bulk_sc_talon_read_annot.tsv"
    )
)

SCRNA_H5AD = (
    config.get("scrna_h5ad")
    if "config" in globals() and "scrna_h5ad" in config
    else os.environ.get("SCRNA_H5AD", "outputs/anndata/short_shallow.h5ad")
)

OUTDIR = "outputs"
TABLES = f"{OUTDIR}/tables"
PLOTS = f"{OUTDIR}/plots"

# Derive the stem from the short-read filename so rules that write
# `<stem>_qc.h5ad` are consistent with the provided input file name.
SCRNA_STEM = (
    config.get("scrna_stem")
    if "config" in globals() and "scrna_stem" in config
    else Path(SCRNA_H5AD).stem
)

# --------------------------------------------------------------------
# FINAL TARGETS
# --------------------------------------------------------------------

rule all:
    input:
        f"{OUTDIR}/{SCRNA_STEM}_qc.h5ad",
        f"{TABLES}/cell_clusters.csv",
        f"{TABLES}/talon_scrna_symbol_map.csv",
        f"{TABLES}/isoform_proxies.csv",
        f"{PLOTS}/umap_clusters.png"

# --------------------------------------------------------------------
# RULES
# --------------------------------------------------------------------

# 1) QC scRNA (wraps scripts/utils.py)
rule qc_scrna:
    input:
        SCRNA_H5AD
    output:
        f"{OUTDIR}/{SCRNA_STEM}_qc.h5ad"
    shell:
        r"""
        mkdir -p {OUTDIR}
        python scripts/utils.py {input} {OUTDIR}
        """

# 2) Cluster cells (wraps scripts/run_cluster_cells.py)
rule cluster_cells:
    input:
        qc_h5ad=f"{OUTDIR}/{SCRNA_STEM}_qc.h5ad"
    output:
        clusters=f"{TABLES}/cell_clusters.csv",
        umap=f"{PLOTS}/umap_clusters.png"
    shell:
        r"""
        mkdir -p {TABLES} {PLOTS}
        python scripts/run_cluster_cells.py \
          {input.qc_h5ad} \
          {output.clusters} \
          --output-dir {PLOTS}
        """

# 3) Map TALON gene symbols to scRNA
rule symbol_map:
    input:
        talon=TALON_TSV,
        qc_h5ad=f"{OUTDIR}/{SCRNA_STEM}_qc.h5ad"
    output:
        f"{TABLES}/talon_scrna_symbol_map.csv"
    shell:
        r"""
        mkdir -p {TABLES}
        python scripts/map_ids_by_symbol.py \
          {input.talon} \
          {input.qc_h5ad} \
          {output}
        """

# 4) Assign isoform proxies
rule isoform_proxies:
    input:
        talon=TALON_TSV,
        qc_h5ad=f"{OUTDIR}/{SCRNA_STEM}_qc.h5ad",
        clusters=f"{TABLES}/cell_clusters.csv",
        symbol_map=f"{TABLES}/talon_scrna_symbol_map.csv"
    output:
        f"{TABLES}/isoform_proxies.csv"
    shell:
        r"""
        python scripts/assign_isoform_proxies.py \
          {input.talon} \
          {input.qc_h5ad} \
          {input.clusters} \
          {input.symbol_map} \
          {output}
        """
