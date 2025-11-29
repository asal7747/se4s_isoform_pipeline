#!/usr/bin/env python
"""
Simple benchmark of isoform diversity:
- Compares number of TALON isoforms per gene (bulk long-read)
  to number of proxy isoforms per gene across clusters (scRNA).

Usage:
  benchmark_isoform_diversity.py <talon_tsv> <isoform_proxies_csv> <out_csv>
"""

import sys
from pathlib import Path

import pandas as pd


def main():
    if len(sys.argv) < 4:
        print(
            "Usage: benchmark_isoform_diversity.py <talon_tsv> <isoform_proxies_csv> <out_csv>",
            file=sys.stderr,
        )
        sys.exit(2)
    talon_tsv, proxies_csv, out_csv = sys.argv[1:4]

    talon_path = Path(talon_tsv)
    proxies_path = Path(proxies_csv)

    if not talon_path.is_file():
        print(f"ERROR: TALON TSV not found: {talon_path}", file=sys.stderr)
        sys.exit(1)
    if not proxies_path.is_file():
        print(f"ERROR: isoform_proxies CSV not found: {proxies_path}", file=sys.stderr)
        sys.exit(1)

    # TALON: number of distinct isoforms per gene
    talon = pd.read_csv(
        talon_path,
        sep="\t",
        usecols=["annot_gene_name", "annot_transcript_id"],
    ).dropna()
    talon["annot_gene_name"] = talon["annot_gene_name"].astype(str)
    talon["annot_transcript_id"] = talon["annot_transcript_id"].astype(str)
    talon_isoforms = (
        talon.groupby("annot_gene_name")["annot_transcript_id"]
        .nunique()
        .rename("talon_n_isoforms")
    )

    # Proxies: number of distinct proxy transcripts per gene across clusters
    proxies = pd.read_csv(proxies_path)
    if not {"gene", "proxy_transcript"}.issubset(proxies.columns):
        print(
            "ERROR: isoform_proxies CSV must have columns 'gene' and 'proxy_transcript'.",
            file=sys.stderr,
        )
        sys.exit(1)

    proxy_isoforms = (
        proxies.groupby("gene")["proxy_transcript"]
        .nunique()
        .rename("proxy_n_isoforms")
    )

    # Join on gene symbol
    df = pd.concat([talon_isoforms, proxy_isoforms], axis=1).fillna(0)
    df["talon_n_isoforms"] = df["talon_n_isoforms"].astype(int)
    df["proxy_n_isoforms"] = df["proxy_n_isoforms"].astype(int)

    df.to_csv(out_csv, index_label="gene")
    print(f"Wrote benchmark isoform diversity table -> {out_csv}")


if __name__ == "__main__":
    main()
