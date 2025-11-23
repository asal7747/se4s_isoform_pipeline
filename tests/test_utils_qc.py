import sys
from pathlib import Path

# Ensure project root is on sys.path so `scripts` is importable
PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scanpy as sc  # only used to reload the QC output for sanity checks

from scripts.utils import load_and_qc_h5ad


@pytest.fixture
def small_h5ad(tmp_path):
    """
    Create a tiny AnnData with a few cells/genes and save as .h5ad.

    This is intentionally small so we can:
    - Run through the QC pipeline quickly.
    - Trigger filtering of at least one low-coverage cell.
    """
    X = np.array(
        [
            [10, 0, 3],
            [5, 2, 0],
            [0, 0, 0],  # low-coverage cell that may be filtered
        ],
        dtype=float,
    )
    obs = pd.DataFrame(index=[f"cell{i}" for i in range(X.shape[0])])
    var = pd.DataFrame(index=["mt-Nd1", "GeneA", "GeneB"])  # tests mt_prefix handling

    adata = ad.AnnData(X=X, obs=obs, var=var)

    in_path = tmp_path / "small.h5ad"
    adata.write_h5ad(in_path)
    return in_path


def test_load_and_qc_h5ad_happy_path(small_h5ad, tmp_path):
    """
    End-to-end check that load_and_qc_h5ad:
    - Runs without error on a small input.
    - Writes a *_qc.h5ad file.
    - Leaves at least one cell and gene after QC.
    """
    out_dir = tmp_path / "out"

    out_path = load_and_qc_h5ad(
        in_h5ad=str(small_h5ad),
        out_dir=str(out_dir),
        min_counts=1,
        min_genes=1,
    )

    out_file = Path(out_path)
    assert out_file.is_file()

    adata_qc = sc.read(out_file)
    assert adata_qc.n_obs > 0
    assert adata_qc.n_vars > 0


def test_load_and_qc_h5ad_missing_file(tmp_path):
    """
    Calling load_and_qc_h5ad on a non-existent file should raise FileNotFoundError.
    """
    missing = tmp_path / "does_not_exist.h5ad"
    out_dir = tmp_path / "out"

    with pytest.raises(FileNotFoundError):
        load_and_qc_h5ad(str(missing), str(out_dir))


def test_load_and_qc_h5ad_bad_thresholds(small_h5ad, tmp_path):
    """
    Non-positive QC thresholds should raise ValueError.
    """
    out_dir = tmp_path / "out"

    with pytest.raises(ValueError):
        load_and_qc_h5ad(
            in_h5ad=str(small_h5ad),
            out_dir=str(out_dir),
            min_counts=0,
            min_genes=10,
        )

    with pytest.raises(ValueError):
        load_and_qc_h5ad(
            in_h5ad=str(small_h5ad),
            out_dir=str(out_dir),
            min_counts=10,
            min_genes=-1,
        )


def test_load_and_qc_h5ad_all_filtered_raises(small_h5ad, tmp_path):
    """
    If QC thresholds are so strict that all cells/genes are filtered out,
    the function should raise RuntimeError instead of silently writing
    an empty AnnData.
    """
    out_dir = tmp_path / "out"

    with pytest.raises(RuntimeError):
        load_and_qc_h5ad(
            in_h5ad=str(small_h5ad),
            out_dir=str(out_dir),
            min_counts=1000,  # high enough to drop all cells in this tiny dataset
            min_genes=1000,
        )
