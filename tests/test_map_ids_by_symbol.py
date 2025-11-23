# tests/test_map_ids_by_symbol.py

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

from scripts import map_ids_by_symbol


@pytest.fixture
def small_talon_tsv(tmp_path):
    """
    Tiny TALON-like TSV with annot_gene_name column.
    """
    df = pd.DataFrame(
        {
            "annot_gene_name": ["GeneA", "GeneB", "GeneC", "GeneA"],
            "annot_transcript_id": ["TX1", "TX2", "TX3", "TX4"],
        }
    )
    tsv_path = tmp_path / "talon_small.tsv"
    df.to_csv(tsv_path, sep="\t", index=False)
    return tsv_path


@pytest.fixture
def small_scrna_h5ad(tmp_path):
    """
    Tiny AnnData where gene symbols are stored in a 'gene_name' column.
    """
    X = np.array(
        [
            [1, 0, 2],
            [0, 3, 1],
        ],
        dtype=float,
    )
    obs = pd.DataFrame(index=["cell1", "cell2"])
    var = pd.DataFrame(
        {
            "gene_name": ["GeneA", "GeneB", "GeneX"],
        },
        index=["id1", "id2", "id3"],
    )
    adata = ad.AnnData(X=X, obs=obs, var=var)

    h5ad_path = tmp_path / "scrna_small.h5ad"
    adata.write_h5ad(h5ad_path)
    return h5ad_path



def test_map_ids_by_symbol_basic_overlap(tmp_path, small_talon_tsv, small_scrna_h5ad, monkeypatch):
    """
    map_ids_by_symbol.py should find overlapping symbols between TALON and AnnData.
    In this synthetic example, GeneA and GeneB should overlap.
    """
    out_csv = tmp_path / "talon_scrna_symbol_map.csv"

    argv = [
        "map_ids_by_symbol.py",
        str(small_talon_tsv),
        str(small_scrna_h5ad),
        str(out_csv),
    ]
    monkeypatch.setattr(sys, "argv", argv)
    map_ids_by_symbol.main()

    assert out_csv.is_file()
    df = pd.read_csv(out_csv)
    # Expect at least GeneA and GeneB
    genes = set(df["gene"].astype(str))
    assert {"GeneA", "GeneB"}.issubset(genes)
    # Match type should be 'symbol'
    assert set(df["match_type"]) == {"symbol"}
