from pathlib import Path

# --------------------------------------------------------------------
# CONFIG
# --------------------------------------------------------------------

# Real input files
TALON_TSV = "/Users/asheralbrecht/Desktop/bulk_sc_talon_read_annot.tsv"
SCRNA_H5AD = "outputs/anndata/short_shallow.h5ad"  # adjust if your raw ENCODE file has a different name

OUTDIR = "outputs"
TABLES = f"{OUTDIR}/tables"
PLOTS = f"{OUTDIR}/plots"

# Stem of the raw scRNA file (short_shallow.h5ad -> short_shallow_qc.h5ad)
SCRNA_STEM = "short_shallow"

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
