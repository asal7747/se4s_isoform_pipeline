#!/usr/bin/env python
"""
Assign cluster-specific isoform proxies from TALON annotations, using a symbol map.

Usage:
  assign_isoform_proxies.py <talon_tsv> <qc_h5ad> luster_csv> <symbol_map_csv> <out_csv>

Notes:
- <symbol_map_csv> is outputs/tables/talon_scrna_symbol_map.csv from map_ids_by_symbol.py with columns:
  match_type,gene   (gene = biological symbol present in both TALON and AnnData)
"""
import sys
from pathlib import Path

import anndata as ad
import pandas as pd


def derive_symbols(adata: ad.AnnData) -> pd.Series:
    # Highest priority: explicit gene_name column
    for col in ["gene_name", "gene_names", "name", "symbol"]:
        if col in adata.var.columns:
            s = adata.var[col].astype(str)
            if s.notna().any():
                return s

    # Next: transcript_id may embed "-Symbol" suffix
    if "transcript_id" in adata.var.columns:
        t = adata.var["transcript_id"].astype(str)
        parts = t.str.rsplit("-", n=1, expand=True)
        if isinstance(parts, pd.DataFrame) and parts.shape[1] == 2:
            return parts[1].astype(str)

    # Fallback: var_names may embed "-Symbol" suffix
    vn = pd.Index(adata.var_names.astype(str))
    parts = vn.str.rsplit("-", n=1, expand=True)
    if isinstance(parts, pd.DataFrame) and parts.shape[1] == 2:
        return parts[1].astype(str)

    # Last resort: use var_names directly (will likely fail to match many)
    return pd.Series(vn, index=adata.var.index, dtype=str)


def main():
    if len(sys.argv) < 6:
        print(
            "Usage: assign_isoform_proxies.py <talon_tsv> <qc_h5ad> "
            "luster_csv> <symbol_map_csv> <out_csv>",
            file=sys.stderr,
        )
        sys.exit(2)
    talon_tsv, h5ad, cluster_csv, symbol_map_csv, out_csv = sys.argv[1:6]

    talon_path = Path(talon_tsv)
    h5ad_path = Path(h5ad)
    cluster_path = Path(cluster_csv)
    symmap_path = Path(symbol_map_csv)

    # Basic input checks
    for p, label in [
        (talon_path, "TALON TSV"),
        (h5ad_path, "QC h5ad"),
        (cluster_path, "cluster CSV"),
        (symmap_path, "symbol map CSV"),
    ]:
        if not p.is_file():
            print(f"ERROR: {label} not found: {p}", file=sys.stderr)
            sys.exit(1)

    # 1) TALON transcript counts grouped by annot_gene_name (assumed biological symbol in symbol map)
    header = pd.read_csv(talon_path, sep="\t", nrows=1)
    missing_cols = [
        c for c in ["annot_gene_name", "annot_transcript_id"] if c not in header.columns
    ]
    if missing_cols:
        print(
            f"ERROR: TALON TSV missing required columns: {', '.join(missing_cols)}",
            file=sys.stderr,
        )
        sys.exit(1)

    talon = (
        pd.read_csv(
            talon_path,
            sep="\t",
            usecols=["annot_gene_name", "annot_transcript_id"],
        )
        .dropna()
    )
    talon["annot_gene_name"] = talon["annot_gene_name"].astype(str)
    talon["annot_transcript_id"] = talon["annot_transcript_id"].astype(str)
    tx_counts = (
        talon.groupby(["annot_gene_name", "annot_transcript_id"])
        .size()
        .reset_index(name="read_count")
    )

    # 2) Shared symbols
    symmap = pd.read_csv(symmap_path)
    if "gene" not in symmap.columns:
        print(
            "ERROR: symbol map CSV missing required column 'gene'.",
            file=sys.stderr,
        )
        sys.exit(1)
    shared_symbols = set(symmap["gene"].astype(str))
    if not shared_symbols:
        print(
            "WARNING: symbol map contains no genes; no proxies will be assigned.",
            file=sys.stderr,
        )

    # 3) Load AnnData and clusters
    adata = ad.read_h5ad(h5ad_path)

    # Derive symbols robustly
    symbols = derive_symbols(adata)
    adata.var = adata.var.copy()
    adata.var["symbol"] = symbols.values

    clusters = pd.read_csv(cluster_path)
    if not {"cell", "cluster"}.issubset(clusters.columns):
        print(
            "ERROR: cluster CSV must have columns 'cell' and 'cluster'.",
            file=sys.stderr,
        )
        sys.exit(1)

    cluster_map = clusters.set_index("cell")["cluster"].astype(str)
    adata = adata[adata.obs_names.isin(cluster_map.index)].copy()
    if adata.n_obs == 0:
        print(
            "ERROR: No overlap between AnnData cells and cluster CSV 'cell' column.",
            file=sys.stderr,
        )
        sys.exit(1)
    adata.obs["cluster"] = adata.obs_names.map(cluster_map)

    # Build expression matrix at symbol level (sum duplicated symbols)
    expr = adata.to_df()
    expr.columns = adata.var["symbol"].values
    expr = expr.groupby(axis=1, level=0).sum()

    # Restrict to shared symbols
    expr = expr.loc[:, expr.columns.isin(shared_symbols)]

    # Compute mean per cluster
    cluster_expr = expr.groupby(adata.obs["cluster"]).mean()

    # 4) Assign top TALON transcript by read_count for expressed genes per cluster
    tx_counts_sym = tx_counts[tx_counts["annot_gene_name"].isin(shared_symbols)]
    proxies = []
    for cluster_id in cluster_expr.index:
        row = cluster_expr.loc[cluster_id]
        expressed = row[row > 0]
        if expressed.empty:
            continue
        for gene, mean_expr in expressed.items():
            gene_txs = tx_counts_sym[tx_counts_sym["annot_gene_name"] == gene]
            if gene_txs.empty:
                continue
            top_tx = gene_txs.sort_values("read_count", ascending=False).iloc[0]
            proxies.append(
                {
                    "cluster": cluster_id,
                    "gene": gene,
                    "mean_expr": float(mean_expr),
                    "proxy_transcript": top_tx["annot_transcript_id"],
                    "proxy_read_count": int(top_tx["read_count"]),
                }
            )

    out_df = pd.DataFrame(
        proxies,
        columns=["cluster", "gene", "mean_expr", "proxy_transcript", "proxy_read_count"],
    )
    out_df.to_csv(out_csv, index=False)

    print(f"Assigned {len(out_df)} cluster-gene proxies -> {out_csv}")
    if out_df.empty:
        print(
            "WARNING: No proxies were assigned. Check symbol overlap and expression thresholds.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
