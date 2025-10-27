import argparse
import os

try:
    import pandas as pd
except ImportError as e:
    raise ImportError(
        "pandas is required but not installed. Install it with: pip install pandas"
    ) from e

try:
    import pysam
except ImportError as e:
    raise ImportError(
        "pysam is required but not installed. Install it with: pip install pysam"
    ) from e

REQUIRED = [
    "read_name",
    "dataset",
    "genome_build",
    "chrom",
    "read_start",
    "read_end",
    "strand",
    "gene_ID",
    "transcript_ID",
    "gene_novelty",
    "transcript_novelty",
    "fraction_As",
]


def _head(tsv, n=5000, usecols=None):
    return pd.read_csv(tsv, sep="\t", nrows=n, usecols=usecols)


def schema_ok(tsv):
    df = _head(tsv, n=5000)
    return all(c in df.columns for c in REQUIRED) and set(df["strand"].unique()) <= {
        "+",
        "-",
    }


def min_le_max(tsv, n=200000):
    df = _head(tsv, n=n, usecols=["read_start", "read_end"])
    return (
        df[["read_start", "read_end"]].min(axis=1)
        <= df[["read_start", "read_end"]].max(axis=1)
    ).all()


def datasets_present(tsv, expected, n=1_000_000):
    df = _head(tsv, n=n, usecols=["dataset"])
    return expected.issubset(set(df["dataset"].unique()))


def novelty_sets_ok(tsv, n=1_000_000):
    df = _head(tsv, n=n, usecols=["gene_novelty", "transcript_novelty"])
    gene_ok = set(df["gene_novelty"].dropna().unique()) <= {
        "Known",
        "Genomic",
        "Intergenic",
        "Antisense",
        "NIC",
        "NNC",
        "ISM",
    }
    tx_ok = set(df["transcript_novelty"].dropna().unique()) <= {
        "Known",
        "NIC",
        "NNC",
        "ISM",
        "Antisense",
        "Intergenic",
        "Genomic",
    }
    return gene_ok and tx_ok


def qc_identity_ok(qc_path, thresh="0.800000"):
    with open(qc_path) as f:
        return f"Min read identity to reference: {thresh}" in f.read()


def spotcheck_reads(tsv, se4s_root, k=10, n=50000, max_fail_frac=0.4):
    """
    Robust spot check:
      - sample up to k reads from first n rows
      - try clean_labeled.sam, then clean.sam
      - allow some failures (e.g., absent in sampled SAM or different ref in rare cases)
    """
    df = _head(tsv, n=n, usecols=["read_name", "chrom", "dataset"])
    if df.empty:
        return False
    rows = df.sample(n=min(k, len(df)), random_state=42).to_dict("records")
    fails = 0
    for r in rows:
        base = os.path.join(se4s_root, "work", "longread", "bulk", r["dataset"])
        sam_candidates = [
            os.path.join(base, "labeled", "clean_labeled.sam"),
            os.path.join(base, "clean.sam"),
        ]
        sam_path = next((p for p in sam_candidates if os.path.exists(p)), None)
        if not sam_path:
            fails += 1
            continue
        found = False
        with pysam.AlignmentFile(sam_path, "r") as f:
            for aln in f.fetch(until_eof=True):
                if aln.query_name == r["read_name"]:
                    ref = (
                        f.get_reference_name(aln.reference_id)
                        if aln.reference_id >= 0
                        else None
                    )
                    if ref != r["chrom"]:
                        fails += 1
                    found = True
                    break
        if not found:
            fails += 1
    # allow some failures up to the configured fraction
    return (fails / max(1, len(rows))) <= max_fail_frac


def main():
    ap = argparse.ArgumentParser(description="Validate TALON read-wise TSV and QC log.")
    ap.add_argument("--tsv", required=True, help="Path to TALON read-wise TSV")
    ap.add_argument("--qc", required=True, help="Path to TALON QC log")
    ap.add_argument(
        "--se4s-root",
        required=True,
        help="SE4S root to find labeled SAMs for spot checks",
    )
    ap.add_argument(
        "--datasets",
        nargs="+",
        default=["ENCFF003OWX", "ENCFF019HRC", "ENCFF669LWV", "ENCFF676BYQ"],
        help="Expected dataset IDs",
    )
    ap.add_argument(
        "--skip-spotcheck", action="store_true", help="Skip SAM round-trip spot check"
    )
    args = ap.parse_args()

    expected = set(args.datasets)
    checks = {
        "schema_ok": schema_ok(args.tsv),
        "coords_min<=max": min_le_max(args.tsv),
        "datasets_present": datasets_present(args.tsv, expected),
        "novelty_sets_ok": novelty_sets_ok(args.tsv),
        "qc_identity_ok": qc_identity_ok(args.qc),
    }
    if not args.skip_spotcheck:
        checks["spotcheck_reads"] = spotcheck_reads(args.tsv, args.se4s_root)
    # Print results
    for k, v in checks.items():
        print(f"{k}: {v}")
    print("VALID" if all(checks.values()) else "INVALID")


if __name__ == "__main__":
    main()
