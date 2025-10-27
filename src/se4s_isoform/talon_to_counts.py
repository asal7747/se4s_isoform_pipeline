from pathlib import Path

import pandas as pd


def isoform_qc_table(tsv_path, dataset=None):
    """Aggregate isoform novelty counts per gene.

    Returns a DataFrame indexed by `annot_gene_id` with columns for each
    novelty class (Known, NIC, NNC, ISM, Antisense, Intergenic, Genomic).

    Parameters
    - tsv_path: path to TALON read-wise TSV
    - dataset: optional dataset filter (string)
    - chunksize: rows per chunk when streaming large TSVs
    """
    cols = ["dataset", "annot_gene_id", "gene_novelty"]
    counts = {}  # (gene, novelty) -> int
    # stream with a sensible default chunk size
    chunk_size = 500_000
    for chunk in pd.read_csv(tsv_path, sep="\t", usecols=cols, chunksize=chunk_size):
        if dataset:
            chunk = chunk[chunk["dataset"] == dataset]
        vc = chunk.groupby(["annot_gene_id", "gene_novelty"]).size()
        for (gene, novelty), v in vc.items():
            counts[(gene, novelty)] = counts.get((gene, novelty), 0) + int(v)
    # Build DataFrame
    if not counts:
        return pd.DataFrame()
    rows = []
    for (gene, novelty), v in counts.items():
        rows.append({"annot_gene_id": gene, "gene_novelty": novelty, "count": v})
    df = pd.DataFrame(rows)
    pivot = df.pivot_table(
        index="annot_gene_id", columns="gene_novelty", values="count", fill_value=0
    )
    # ensure column order is reasonable
    preferred = ["Known", "NIC", "NNC", "ISM", "Antisense", "Intergenic", "Genomic"]
    cols = [c for c in preferred if c in pivot.columns] + [
        c for c in pivot.columns if c not in preferred
    ]
    return pivot[cols]


def write_isoform_qc_table(tsv_path, out_dir, dataset=None):
    """Write isoform QC pivot table CSV and return the file path as a string.

    The filename is deterministic so callers can glob predictably:
    `isoform_qc_table.csv` or `isoform_qc_table_<dataset>.csv`.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    df = isoform_qc_table(tsv_path, dataset)
    name = f"isoform_qc_table{'_' + dataset if dataset else ''}.csv"
    out_path = out / name
    df.to_csv(out_path)
    return str(out_path)


def talon_counts_by_transcript(tsv_path, dataset=None, chunksize=500_000):
    cols = ["dataset", "transcript_ID"]
    counts = {}
    for chunk in pd.read_csv(tsv_path, sep="\t", usecols=cols, chunksize=chunksize):
        if dataset:
            chunk = chunk[chunk["dataset"] == dataset]
        vc = chunk["transcript_ID"].value_counts()
        for k, v in vc.items():
            counts[k] = counts.get(k, 0) + int(v)
    return pd.Series(counts, name="count").sort_values(ascending=False)


def talon_counts_by_gene(tsv_path, dataset=None, chunksize=500_000):
    cols = ["dataset", "annot_gene_id"]
    counts = {}
    for chunk in pd.read_csv(tsv_path, sep="\t", usecols=cols, chunksize=chunksize):
        if dataset:
            chunk = chunk[chunk["dataset"] == dataset]
        vc = chunk["annot_gene_id"].value_counts()
        for k, v in vc.items():
            counts[k] = counts.get(k, 0) + int(v)
    return pd.Series(counts, name="count").sort_values(ascending=False)


def write_top_tables(tsv_path, out_dir, dataset=None, top=100):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    s_tx = talon_counts_by_transcript(tsv_path, dataset)
    s_g = talon_counts_by_gene(tsv_path, dataset)
    s_tx.head(top).to_csv(
        out / f"top_transcripts{('_' + dataset) if dataset else ''}.csv"
    )
    s_g.head(top).to_csv(out / f"top_genes{('_' + dataset) if dataset else ''}.csv")
    return (
        out / f"top_transcripts{('_' + dataset) if dataset else ''}.csv",
        out / f"top_genes{('_' + dataset) if dataset else ''}.csv",
    )
