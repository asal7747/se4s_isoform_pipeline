# tests/test_assign_isoform_proxies.py

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

# Ensure project root on path
PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts import assign_isoform_proxies


@pytest.fixture
def small_talon_tsv(tmp_path):
    """
    Tiny TALON-like TSV with two genes and multiple transcripts per gene,
    with differing implicit read counts (via repeated rows).
    """
    base = pd.DataFrame(
        {
            "annot_gene_name": ["GeneA", "GeneA", "GeneB", "GeneB"],
            "annot_transcript_id": ["TX1", "TX2", "TX3", "TX4"],
        }
    )

    rows = []
    # TX1 x3
    rows.append(base.iloc[[0]].copy())
    rows.append(base.iloc[[0]].copy())
    rows.append(base.iloc[[0]].copy())
    # TX2 x1
    rows.append(base.iloc[[1]].copy())
    # TX3 x2
    rows.append(base.iloc[[2]].copy())
    rows.append(base.iloc[[2]].copy())
    # TX4 x1
    rows.append(base.iloc[[3]].copy())

    df = pd.concat(rows, ignore_index=True)

    tsv_path = tmp_path / "talon_small.tsv"
    df.to_csv(tsv_path, sep="\t", index=False)
    return tsv_path


@pytest.fixture
def small_qc_h5ad(tmp_path):
    """
    Tiny AnnData with two genes (GeneA, GeneB) and two clusters.
    """
    X = np.array(
        [
            [10, 0],  # cell1: GeneA high, GeneB 0
            [5,  2],  # cell2: both expressed
        ],
        dtype=float,
    )
    obs = pd.DataFrame(index=["cell1", "cell2"])
    var = pd.DataFrame(index=["GeneA", "GeneB"])
    adata = ad.AnnData(X=X, obs=obs, var=var)

    h5ad_path = tmp_path / "scrna_qc_small.h5ad"
    adata.write_h5ad(h5ad_path)
    return h5ad_path


@pytest.fixture
def small_cluster_csv(tmp_path):
    """
    Simple cluster assignment: both cells in cluster '0'.
    """
    df = pd.DataFrame({"cell": ["cell1", "cell2"], "cluster": ["0", "0"]})
    csv_path = tmp_path / "cell_clusters_small.csv"
    df.to_csv(csv_path, index=False)
    return csv_path


@pytest.fixture
def small_symbol_map_csv(tmp_path):
    """
    Symbol map listing GeneA and GeneB as shared symbols.
    """
    df = pd.DataFrame({"match_type": ["symbol", "symbol"], "gene": ["GeneA", "GeneB"]})
    csv_path = tmp_path / "talon_scrna_symbol_map_small.csv"
    df.to_csv(csv_path, index=False)
    return csv_path


def test_assign_isoform_proxies_basic(tmp_path, small_talon_tsv, small_qc_h5ad, small_cluster_csv, small_symbol_map_csv, monkeypatch):
    """
    assign_isoform_proxies.py should:
    - Produce at least one proxy row.
    - Choose the highest-read-count transcript per gene (TX1 for GeneA, TX3 for GeneB in this setup).
    """
    out_csv = tmp_path / "isoform_proxies_small.csv"

    argv = [
        "assign_isoform_proxies.py",
        str(small_talon_tsv),
        str(small_qc_h5ad),
        str(small_cluster_csv),
        str(small_symbol_map_csv),
        str(out_csv),
    ]
    monkeypatch.setattr(sys, "argv", argv)
    assign_isoform_proxies.main()

    assert out_csv.is_file()
    df = pd.read_csv(out_csv)

    # Should have at least GeneA and GeneB proxies for cluster "0"
    assert not df.empty
    assert set(df["cluster"].astype(str)) == {"0"}

    proxies_by_gene = df.set_index("gene")["proxy_transcript"].to_dict()
    # Expect highest-read-count transcripts chosen by our synthetic tallies
    assert proxies_by_gene.get("GeneA") == "TX1"
    assert proxies_by_gene.get("GeneB") == "TX3"
