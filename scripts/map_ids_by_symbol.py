#!/usr/bin/env python
import sys, pandas as pd, anndata as ad

def main():
    if len(sys.argv) < 4:
        print("Usage: map_ids_by_symbol.py <talon_tsv> <h5ad> <out_csv>")
        sys.exit(2)
    tsv, h5ad, out_csv = sys.argv[1], sys.argv[2], sys.argv[3]

    # Load TALON symbols
    talon = pd.read_csv(tsv, sep="\t", usecols=lambda c: c in {"annot_gene_name"}).dropna().drop_duplicates()
    talon_syms = set(talon["annot_gene_name"].astype(str))

    # Load AnnData
    adata = ad.read_h5ad(h5ad)

    # Candidate symbol sources in AnnData
    candidates = []

    # 1) Parse from var_names that may look like "ENSMUSG...-Symbol"
    vn = pd.Index(adata.var_names.astype(str))
    parts = vn.str.rsplit("-", n=1, expand=True)
    if isinstance(parts, pd.DataFrame) and parts.shape[1] == 2:
        candidates.append(parts[1].astype(str))

    # 2) Common symbol columns if present
    for col in ["gene_name", "gene_names", "name", "symbol"]:
        if col in adata.var.columns:
            candidates.append(adata.var[col].astype(str))

    # 3) Fallback: if transcript_id contains a trailing "-Symbol", parse it
    if "transcript_id" in adata.var.columns:
        tparts = adata.var["transcript_id"].astype(str).str.rsplit("-", n=1, expand=True)
        if isinstance(tparts, pd.DataFrame) and tparts.shape[1] == 2:
            candidates.append(tparts[1].astype(str))

    if not candidates:
        print("ERROR: Could not derive gene symbols from AnnData (no '-' split or symbol columns found)")
        sys.exit(1)

    # Union of all candidate symbols
    ad_syms = set()
    for s in candidates:
        ad_syms |= set(s.dropna().astype(str))

    matched_syms = sorted(ad_syms & talon_syms)
    out = pd.DataFrame({"match_type": "symbol", "gene": matched_syms})
    out.to_csv(out_csv, index=False)
    print(f"Matched {len(matched_syms)} by symbol -> {out_csv}")

if __name__ == "__main__":
    main()
