# tests/test_cluster_cells.py

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scanpy as sc

# Make sure project root is on sys.path so `scripts` is importable
PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts import cluster_cells


@pytest.fixture
def small_qc_h5ad(tmp_path):
    """
    Create a tiny QC'd AnnData object suitable for clustering.
    We keep it simple but non-degenerate.
    """
    X = np.array(
        [
            [5, 1, 0],
            [3, 0, 2],
            [4, 1, 1],
            [2, 0, 3],
        ],
        dtype=float,
    )
    obs = pd.DataFrame(index=[f"cell{i}" for i in range(X.shape[0])])
    var = pd.DataFrame(index=["GeneA", "GeneB", "GeneC"])
    adata = ad.AnnData(X=X, obs=obs, var=var)

    h5ad_path = tmp_path / "small_qc.h5ad"
    adata.write_h5ad(h5ad_path)
    return h5ad_path


def test_cluster_cells_creates_outputs(tmp_path, small_qc_h5ad, monkeypatch):
    """
    Smoke test for cluster_cells.main:
    - Runs on a small QC h5ad.
    - Creates a clusters CSV with expected columns.
    - Creates a UMAP PNG in the outputs/ directory.
    """

    out_csv = tmp_path / "cell_clusters.csv"

    # Run the script's main() with fake argv
    argv = ["cluster_cells.py", str(small_qc_h5ad), str(out_csv)]
    monkeypatch.setattr(sys, "argv", argv)
    cluster_cells.main()

    # Check clusters CSV exists and has the expected structure
    assert out_csv.is_file()
    df = pd.read_csv(out_csv)
    assert set(df.columns) == {"cell", "cluster"}
    assert len(df) == 4
    assert df["cluster"].nunique() >= 1

    # Check that a UMAP plot was written under outputs/
    umap_png = PROJECT_ROOT / "outputs" / "umap_clusters.png"
    assert umap_png.is_file()
