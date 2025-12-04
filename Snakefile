from pathlib import Path
import os

# --------------------------------------------------------------------
# CONFIG
# --------------------------------------------------------------------

# Input files can be provided via Snakemake `config` (e.g. --config talon_tsv=...) or
# via environment variables (TALON_TSV). If neither is provided we fall back to
# sensible repository-relative defaults so the Snakefile is functional out-of-the-box.
TALON_TSV = (
    config.get("talon_tsv")
    if "config" in globals() and "talon_tsv" in config
    else os.environ.get(
        "TALON_TSV", "data/bulk_sc_talon_read_annot.tsv"
    )
)

OUTDIR = "outputs"
ANNDATA = f"{OUTDIR}/anndata"
TABLES = f"{OUTDIR}/tables"
PLOTS = f"{OUTDIR}/plots"
COMPARISON = f"{OUTDIR}/comparison"

# The combined QC'd short-read file produced by single_cell_analysis.py
COMBINED_QC_H5AD = f"{ANNDATA}/combined_short_read_qc.h5ad"

# --------------------------------------------------------------------
# FINAL TARGETS
# --------------------------------------------------------------------

rule all:
    input:
        COMBINED_QC_H5AD,
        f"{ANNDATA}/combined_long_read_gene.h5ad",
        f"{ANNDATA}/combined_long_read_transcript.h5ad",
        f"{TABLES}/cell_clusters.csv",
        f"{TABLES}/talon_scrna_symbol_map.csv",
        f"{TABLES}/isoform_proxies.csv",
        f"{PLOTS}/umap_clusters.png",
        f"{COMPARISON}/cluster_isoform_switching.csv"

# --------------------------------------------------------------------
# RULES
# --------------------------------------------------------------------

# 0) Combine and QC all datasets (wraps scripts/single_cell_analysis.py)
#    This downloads data if needed, combines short+long read datasets,
#    runs QC, and produces the combined h5ad files.
rule combine_and_qc:
    output:
        qc_h5ad=COMBINED_QC_H5AD,
        long_gene=f"{ANNDATA}/combined_long_read_gene.h5ad",
        long_transcript=f"{ANNDATA}/combined_long_read_transcript.h5ad"
    shell:
        r"""
        mkdir -p {ANNDATA}
        python scripts/single_cell_analysis.py
        """

# 1) Cluster cells (wraps scripts/run_cluster_cells.py)
rule cluster_cells:
    input:
        qc_h5ad=COMBINED_QC_H5AD
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

# 2) Map TALON gene symbols to scRNA
rule symbol_map:
    input:
        talon=TALON_TSV,
        qc_h5ad=COMBINED_QC_H5AD
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

# 3) Assign isoform proxies
rule isoform_proxies:
    input:
        talon=TALON_TSV,
        qc_h5ad=COMBINED_QC_H5AD,
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

# 4) Compare bulk and single-cell isoforms
rule compare_bulk_sc:
    input:
        talon=TALON_TSV,
        long_gene=f"{ANNDATA}/combined_long_read_gene.h5ad",
        long_transcript=f"{ANNDATA}/combined_long_read_transcript.h5ad"
    output:
        f"{COMPARISON}/cluster_isoform_switching.csv"
    shell:
        r"""
        mkdir -p {COMPARISON}
        python scripts/compare_bulk_sc.py \
          --bulk-talon {input.talon} \
          --output {COMPARISON}
        """
