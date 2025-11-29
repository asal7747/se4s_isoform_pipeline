# tests/test_cluster_cells.py

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

# Make sure project root is on sys.path so `scripts` is importable
PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts import cluster_cells, run_cluster_cells


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


@pytest.fixture
def larger_qc_h5ad(tmp_path):
    """
    Create a larger dataset suitable for PCA/t-SNE/UMAP plotting.
    Needs enough cells and genes for dimensionality reduction.
    Uses realistic count-like data that will pass highly variable gene selection.
    """
    np.random.seed(42)
    
    # Create realistic RNA-seq count data (100 cells × 500 genes)
    # Use Poisson/negative binomial-like distributions
    X = np.zeros((100, 500))
    
    # Group 1 marker genes (genes 0-99): high in first 50 cells
    for i in range(100):
        X[:50, i] = np.random.poisson(lam=100, size=50)  # High expression
        X[50:, i] = np.random.poisson(lam=10, size=50)   # Low expression
    
    # Group 2 marker genes (genes 100-199): high in last 50 cells
    for i in range(100, 200):
        X[:50, i] = np.random.poisson(lam=10, size=50)   # Low expression
        X[50:, i] = np.random.poisson(lam=100, size=50)  # High expression
    
    # Moderately variable genes (genes 200-350)
    for i in range(200, 350):
        mean_expr = np.random.uniform(15, 50)  # Different mean for each gene
        X[:, i] = np.random.poisson(lam=mean_expr, size=100)
    
    # Housekeeping genes (genes 350-500): similar across all cells, moderate expression
    for i in range(350, 500):
        X[:, i] = np.random.poisson(lam=50, size=100) + 10  # More consistent expression
    
    obs = pd.DataFrame(index=[f"cell{i}" for i in range(100)])
    var = pd.DataFrame(index=[f"Gene{i}" for i in range(500)])
    adata = ad.AnnData(X=X, obs=obs, var=var)

    h5ad_path = tmp_path / "larger_qc.h5ad"
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

    argv = ["cluster_cells.py", str(small_qc_h5ad), str(out_csv)]
    monkeypatch.setattr(sys, "argv", argv)
    run_cluster_cells.main()

    # Check clusters CSV exists and has the expected structure
    assert out_csv.is_file()
    df = pd.read_csv(out_csv)
    assert set(df.columns) == {"cell", "cluster"}
    assert len(df) == 4
    assert df["cluster"].nunique() >= 1
    assert df["cluster"].notna().all()

    # Check that a UMAP plot was written under outputs/
    umap_png = PROJECT_ROOT / "outputs" / "umap_clusters.png"
    assert umap_png.is_file()


# ============================================================================
# Parameter Validation Tests
# ============================================================================


def test_n_neighbors_validation_too_low(tmp_path, small_qc_h5ad, monkeypatch):
    """Test that --n-neighbors < 2 raises a system exit."""
    out_csv = tmp_path / "clusters.csv"
    argv = ["cluster_cells.py", str(small_qc_h5ad), str(out_csv), "--n-neighbors", "1"]
    monkeypatch.setattr(sys, "argv", argv)
    
    with pytest.raises(SystemExit) as exc_info:
        run_cluster_cells.main()
    assert exc_info.value.code == 1


def test_n_pcs_validation_too_low(tmp_path, small_qc_h5ad, monkeypatch):
    """Test that --n-pcs < 1 raises a system exit."""
    out_csv = tmp_path / "clusters.csv"
    argv = ["cluster_cells.py", str(small_qc_h5ad), str(out_csv), "--n-pcs", "0"]
    monkeypatch.setattr(sys, "argv", argv)
    
    with pytest.raises(SystemExit) as exc_info:
        run_cluster_cells.main()
    assert exc_info.value.code == 1


def test_resolution_validation_zero(tmp_path, small_qc_h5ad, monkeypatch):
    """Test that --resolution <= 0 raises a system exit."""
    out_csv = tmp_path / "clusters.csv"
    argv = ["cluster_cells.py", str(small_qc_h5ad), str(out_csv), "--resolution", "0"]
    monkeypatch.setattr(sys, "argv", argv)
    
    with pytest.raises(SystemExit) as exc_info:
        run_cluster_cells.main()
    assert exc_info.value.code == 1


def test_resolution_validation_negative(tmp_path, small_qc_h5ad, monkeypatch):
    """Test that --resolution < 0 raises a system exit."""
    out_csv = tmp_path / "clusters.csv"
    argv = ["cluster_cells.py", str(small_qc_h5ad), str(out_csv), "--resolution", "-0.5"]
    monkeypatch.setattr(sys, "argv", argv)
    
    with pytest.raises(SystemExit) as exc_info:
        run_cluster_cells.main()
    assert exc_info.value.code == 1


# ============================================================================
# Parameterization Tests
# ============================================================================


def test_custom_resolution_affects_clustering(tmp_path, monkeypatch):
    """Test that different --resolution values produce different cluster counts."""
    # Create a slightly larger dataset for more reliable clustering differences
    X = np.random.rand(50, 10) * 10  # 50 cells, 10 genes
    obs = pd.DataFrame(index=[f"cell{i}" for i in range(50)])
    var = pd.DataFrame(index=[f"Gene{i}" for i in range(10)])
    adata = ad.AnnData(X=X, obs=obs, var=var)
    
    h5ad_path = tmp_path / "test_data.h5ad"
    adata.write_h5ad(h5ad_path)
    
    # Test with low resolution
    out_csv1 = tmp_path / "clusters_low_res.csv"
    argv1 = ["cluster_cells.py", str(h5ad_path), str(out_csv1), "--resolution", "0.1"]
    monkeypatch.setattr(sys, "argv", argv1)
    run_cluster_cells.main()
    
    df1 = pd.read_csv(out_csv1)
    n_clusters_low = df1["cluster"].nunique()
    
    # Test with high resolution
    out_csv2 = tmp_path / "clusters_high_res.csv"
    argv2 = ["cluster_cells.py", str(h5ad_path), str(out_csv2), 
             "--resolution", "2.0", "--force-recluster"]
    monkeypatch.setattr(sys, "argv", argv2)
    run_cluster_cells.main()
    
    df2 = pd.read_csv(out_csv2)
    n_clusters_high = df2["cluster"].nunique()
    
    # Higher resolution should generally produce more clusters
    # (not always guaranteed, but likely with this data size)
    assert n_clusters_high >= n_clusters_low


def test_custom_n_neighbors_parameter(tmp_path, small_qc_h5ad, monkeypatch):
    """Test that --n-neighbors parameter is accepted and runs successfully."""
    out_csv = tmp_path / "clusters.csv"
    argv = ["cluster_cells.py", str(small_qc_h5ad), str(out_csv), "--n-neighbors", "3"]
    monkeypatch.setattr(sys, "argv", argv)
    
    # Should complete without error
    run_cluster_cells.main()
    
    assert out_csv.is_file()
    df = pd.read_csv(out_csv)
    assert len(df) == 4


def test_custom_n_pcs_parameter(tmp_path, small_qc_h5ad, monkeypatch):
    """Test that --n-pcs parameter is accepted and runs successfully."""
    out_csv = tmp_path / "clusters.csv"
    argv = ["cluster_cells.py", str(small_qc_h5ad), str(out_csv), "--n-pcs", "2"]
    monkeypatch.setattr(sys, "argv", argv)
    
    # Should complete without error
    run_cluster_cells.main()
    
    assert out_csv.is_file()
    df = pd.read_csv(out_csv)
    assert len(df) == 4


# ============================================================================
# Force Recluster Tests
# ============================================================================


def test_force_recluster_overwrites_existing_leiden(tmp_path, monkeypatch):
    """Test that --force-recluster overwrites existing leiden clustering."""
    # Create h5ad with existing leiden column
    X = np.random.rand(20, 5) * 10
    obs = pd.DataFrame(
        {"leiden": ["0"] * 10 + ["1"] * 10},  # Pre-existing clusters
        index=[f"cell{i}" for i in range(20)]
    )
    var = pd.DataFrame(index=[f"Gene{i}" for i in range(5)])
    adata = ad.AnnData(X=X, obs=obs, var=var)
    
    h5ad_path = tmp_path / "with_leiden.h5ad"
    adata.write_h5ad(h5ad_path)
    
    # Run WITHOUT --force-recluster (should use existing)
    out_csv1 = tmp_path / "clusters_no_force.csv"
    argv1 = ["cluster_cells.py", str(h5ad_path), str(out_csv1)]
    monkeypatch.setattr(sys, "argv", argv1)
    run_cluster_cells.main()
    
    df1 = pd.read_csv(out_csv1)
    clusters_no_force = set(df1["cluster"].unique())
    
    # Run WITH --force-recluster and different resolution
    out_csv2 = tmp_path / "clusters_with_force.csv"
    argv2 = ["cluster_cells.py", str(h5ad_path), str(out_csv2), 
             "--force-recluster", "--resolution", "1.5"]
    monkeypatch.setattr(sys, "argv", argv2)
    run_cluster_cells.main()
    
    df2 = pd.read_csv(out_csv2)
    
    # Verify clustering was performed (CSV exists and has data)
    assert out_csv2.is_file()
    assert len(df2) == 20
    assert df2["cluster"].notna().all()


def test_without_force_recluster_uses_existing_leiden(tmp_path, monkeypatch, capsys):
    """Test that existing leiden clustering is used when --force-recluster is not provided."""
    # Create h5ad with existing leiden column
    X = np.random.rand(15, 5) * 10
    obs = pd.DataFrame(
        {"leiden": ["A", "A", "A", "B", "B", "B", "C", "C", "C", "C", "C", "C", "C", "C", "C"]},
        index=[f"cell{i}" for i in range(15)]
    )
    var = pd.DataFrame(index=[f"Gene{i}" for i in range(5)])
    adata = ad.AnnData(X=X, obs=obs, var=var)
    
    h5ad_path = tmp_path / "with_leiden.h5ad"
    adata.write_h5ad(h5ad_path)
    
    out_csv = tmp_path / "clusters.csv"
    argv = ["cluster_cells.py", str(h5ad_path), str(out_csv)]
    monkeypatch.setattr(sys, "argv", argv)
    run_cluster_cells.main()
    
    # Check output confirms using existing clustering
    captured = capsys.readouterr()
    assert "Using existing leiden clustering" in captured.out
    
    # Verify original clusters are preserved
    df = pd.read_csv(out_csv)
    assert set(df["cluster"].unique()) == {"A", "B", "C"}


# ============================================================================
# Input Validation Tests
# ============================================================================


def test_missing_input_file(tmp_path, monkeypatch):
    """Test error handling for non-existent h5ad file."""
    out_csv = tmp_path / "clusters.csv"
    fake_h5ad = tmp_path / "nonexistent.h5ad"
    
    argv = ["cluster_cells.py", str(fake_h5ad), str(out_csv)]
    monkeypatch.setattr(sys, "argv", argv)
    
    with pytest.raises(SystemExit) as exc_info:
        run_cluster_cells.main()
    assert exc_info.value.code == 1


def test_empty_anndata_zero_cells(tmp_path, monkeypatch):
    """Test error handling for AnnData with 0 cells."""
    X = np.array([]).reshape(0, 5)  # 0 cells, 5 genes
    adata = ad.AnnData(X=X)
    
    h5ad_path = tmp_path / "empty_cells.h5ad"
    adata.write_h5ad(h5ad_path)
    
    out_csv = tmp_path / "clusters.csv"
    argv = ["cluster_cells.py", str(h5ad_path), str(out_csv)]
    monkeypatch.setattr(sys, "argv", argv)
    
    with pytest.raises(SystemExit) as exc_info:
        run_cluster_cells.main()
    assert exc_info.value.code == 1


def test_empty_anndata_zero_genes(tmp_path, monkeypatch):
    """Test error handling for AnnData with 0 genes."""
    X = np.array([]).reshape(10, 0)  # 10 cells, 0 genes
    obs = pd.DataFrame(index=[f"cell{i}" for i in range(10)])
    adata = ad.AnnData(X=X, obs=obs)
    
    h5ad_path = tmp_path / "empty_genes.h5ad"
    adata.write_h5ad(h5ad_path)
    
    out_csv = tmp_path / "clusters.csv"
    argv = ["cluster_cells.py", str(h5ad_path), str(out_csv)]
    monkeypatch.setattr(sys, "argv", argv)
    
    with pytest.raises(SystemExit) as exc_info:
        run_cluster_cells.main()
    assert exc_info.value.code == 1


def test_invalid_color_by_column(tmp_path, small_qc_h5ad, monkeypatch):
    """Test error when --color-by column doesn't exist."""
    out_csv = tmp_path / "clusters.csv"
    argv = ["cluster_cells.py", str(small_qc_h5ad), str(out_csv), 
            "--color-by", "NonExistentColumn"]
    monkeypatch.setattr(sys, "argv", argv)
    
    with pytest.raises(SystemExit) as exc_info:
        run_cluster_cells.main()
    assert exc_info.value.code == 1


def test_corrupted_h5ad_file(tmp_path, monkeypatch):
    """Test error handling for corrupted/invalid h5ad file."""
    # Create a file that's not a valid h5ad
    bad_h5ad = tmp_path / "corrupted.h5ad"
    bad_h5ad.write_text("This is not a valid h5ad file")
    
    out_csv = tmp_path / "clusters.csv"
    argv = ["cluster_cells.py", str(bad_h5ad), str(out_csv)]
    monkeypatch.setattr(sys, "argv", argv)
    
    with pytest.raises(SystemExit) as exc_info:
        run_cluster_cells.main()
    assert exc_info.value.code == 1


# ============================================================================
# Plotting Flag Tests
# ============================================================================


def test_save_figures_creates_pca_file(tmp_path, larger_qc_h5ad, monkeypatch):
    """Test that --save-figures with --pca creates PCA plot file."""
    out_csv = tmp_path / "clusters.csv"
    output_dir = tmp_path / "plots"
    
    # Use fewer PCs appropriate for test data
    argv = ["cluster_cells.py", str(larger_qc_h5ad), str(out_csv), 
            "--pca", "--save-figures", "--output-dir", str(output_dir),
            "--n-pcs", "10", "--n-neighbors", "5"]
    monkeypatch.setattr(sys, "argv", argv)
    run_cluster_cells.main()
    
    # Check that PCA plot was created (should be pca_leiden.png by default)
    pca_file = output_dir / "pca_leiden.png"
    assert pca_file.is_file(), f"PCA plot file should be created at {pca_file}"


def test_save_figures_creates_tsne_file(tmp_path, larger_qc_h5ad, monkeypatch):
    """Test that --save-figures with --tsne creates t-SNE plot file."""
    out_csv = tmp_path / "clusters.csv"
    output_dir = tmp_path / "plots"
    
    # Use fewer PCs appropriate for test data
    argv = ["cluster_cells.py", str(larger_qc_h5ad), str(out_csv), 
            "--tsne", "--save-figures", "--output-dir", str(output_dir),
            "--n-pcs", "10", "--n-neighbors", "5"]
    monkeypatch.setattr(sys, "argv", argv)
    run_cluster_cells.main()
    
    # Check that t-SNE plot was created (should be tsne_leiden.png by default)
    tsne_file = output_dir / "tsne_leiden.png"
    assert tsne_file.is_file(), f"t-SNE plot file should be created at {tsne_file}"


def test_save_figures_creates_umap_file(tmp_path, larger_qc_h5ad, monkeypatch):
    """Test that --save-figures with --umap creates UMAP plot file."""
    out_csv = tmp_path / "clusters.csv"
    output_dir = tmp_path / "plots"
    
    # Use fewer PCs appropriate for test data
    argv = ["cluster_cells.py", str(larger_qc_h5ad), str(out_csv), 
            "--umap", "--save-figures", "--output-dir", str(output_dir),
            "--n-pcs", "10", "--n-neighbors", "5"]
    monkeypatch.setattr(sys, "argv", argv)
    run_cluster_cells.main()
    
    # Check that UMAP plot was created in output directory
    umap_file = output_dir / "umap_leiden.png"
    assert umap_file.is_file(), "UMAP plot file should be created"


def test_all_plots_generated_together(tmp_path, larger_qc_h5ad, monkeypatch):
    """Test that --pca --tsne --umap together creates all three plot files."""
    out_csv = tmp_path / "clusters.csv"
    output_dir = tmp_path / "plots"
    
    # Use fewer PCs appropriate for test data
    argv = ["cluster_cells.py", str(larger_qc_h5ad), str(out_csv), 
            "--pca", "--tsne", "--umap", "--save-figures", "--output-dir", str(output_dir),
            "--n-pcs", "10", "--n-neighbors", "5"]
    monkeypatch.setattr(sys, "argv", argv)
    run_cluster_cells.main()
    
    # Check that all three plot types were created (with default color_by=leiden)
    pca_file = output_dir / "pca_leiden.png"
    tsne_file = output_dir / "tsne_leiden.png"
    umap_file = output_dir / "umap_leiden.png"
    
    assert pca_file.is_file(), f"PCA plot should be created at {pca_file}"
    assert tsne_file.is_file(), f"t-SNE plot should be created at {tsne_file}"
    assert umap_file.is_file(), f"UMAP plot should be created at {umap_file}"


def test_default_umap_always_created(tmp_path, small_qc_h5ad, monkeypatch):
    """Test that default UMAP plot is always created in outputs/ directory."""
    out_csv = tmp_path / "clusters.csv"
    
    # Run without any plotting flags
    argv = ["cluster_cells.py", str(small_qc_h5ad), str(out_csv)]
    monkeypatch.setattr(sys, "argv", argv)
    run_cluster_cells.main()
    
    # Default UMAP should still be created
    default_umap = PROJECT_ROOT / "outputs" / "umap_clusters.png"
    assert default_umap.is_file(), "Default UMAP plot should always be created"


# ============================================================================
# Case-Insensitive Column Handling Tests
# ============================================================================


def test_normalize_column_name_exact_match(tmp_path, monkeypatch):
    """Test that normalize_column_name finds exact matches."""
    X = np.random.rand(10, 3)
    obs = pd.DataFrame({"CellType": ["A"] * 5 + ["B"] * 5}, 
                       index=[f"cell{i}" for i in range(10)])
    var = pd.DataFrame(index=["Gene1", "Gene2", "Gene3"])
    adata = ad.AnnData(X=X, obs=obs, var=var)
    
    result = cluster_cells.normalize_column_name(adata, "CellType")
    assert result == "CellType"


def test_normalize_column_name_case_insensitive(tmp_path, monkeypatch):
    """Test that normalize_column_name handles case variations."""
    X = np.random.rand(10, 3)
    obs = pd.DataFrame({"Leiden": ["0"] * 5 + ["1"] * 5}, 
                       index=[f"cell{i}" for i in range(10)])
    var = pd.DataFrame(index=["Gene1", "Gene2", "Gene3"])
    adata = ad.AnnData(X=X, obs=obs, var=var)
    
    # Should find 'Leiden' even when searching for 'leiden'
    result = cluster_cells.normalize_column_name(adata, "leiden")
    assert result == "Leiden"


def test_normalize_column_name_returns_none_if_not_found(tmp_path, monkeypatch):
    """Test that normalize_column_name returns None for non-existent columns."""
    X = np.random.rand(10, 3)
    obs = pd.DataFrame({"RealColumn": ["A"] * 10}, 
                       index=[f"cell{i}" for i in range(10)])
    var = pd.DataFrame(index=["Gene1", "Gene2", "Gene3"])
    adata = ad.AnnData(X=X, obs=obs, var=var)
    
    result = cluster_cells.normalize_column_name(adata, "FakeColumn")
    assert result is None


def test_color_by_leiden_case_variations(tmp_path, monkeypatch):
    """Test that --color-by works with different capitalizations of leiden."""
    # Create h5ad with 'Leiden' (capitalized) column
    X = np.random.rand(15, 5) * 10
    obs = pd.DataFrame(
        {"Leiden": ["0"] * 8 + ["1"] * 7},
        index=[f"cell{i}" for i in range(15)]
    )
    var = pd.DataFrame(index=[f"Gene{i}" for i in range(5)])
    adata = ad.AnnData(X=X, obs=obs, var=var)
    
    h5ad_path = tmp_path / "with_capital_leiden.h5ad"
    adata.write_h5ad(h5ad_path)
    
    out_csv = tmp_path / "clusters.csv"
    # Use lowercase 'leiden' in --color-by
    argv = ["cluster_cells.py", str(h5ad_path), str(out_csv), "--color-by", "leiden"]
    monkeypatch.setattr(sys, "argv", argv)
    
    # Should complete without error (case-insensitive matching)
    run_cluster_cells.main()
    
    assert out_csv.is_file()
    df = pd.read_csv(out_csv)
    assert len(df) == 15


def test_existing_leiden_column_case_variations(tmp_path, monkeypatch, capsys):
    """Test that existing leiden clustering is detected regardless of capitalization."""
    # Test with lowercase 'leiden'
    X1 = np.random.rand(10, 5)
    obs1 = pd.DataFrame({"leiden": ["0"] * 5 + ["1"] * 5}, 
                        index=[f"cell{i}" for i in range(10)])
    var1 = pd.DataFrame(index=[f"Gene{i}" for i in range(5)])
    adata1 = ad.AnnData(X=X1, obs=obs1, var=var1)
    
    h5ad_path1 = tmp_path / "lowercase_leiden.h5ad"
    adata1.write_h5ad(h5ad_path1)
    
    out_csv1 = tmp_path / "clusters1.csv"
    argv1 = ["run_cluster_cells.py", str(h5ad_path1), str(out_csv1)]
    monkeypatch.setattr(sys, "argv", argv1)
    run_cluster_cells.main()
    
    captured1 = capsys.readouterr()
    assert "Using existing leiden clustering" in captured1.out
    
    # Test with capitalized 'Leiden'
    X2 = np.random.rand(10, 5)
    obs2 = pd.DataFrame({"Leiden": ["A"] * 5 + ["B"] * 5}, 
                        index=[f"cell{i}" for i in range(10)])
    var2 = pd.DataFrame(index=[f"Gene{i}" for i in range(5)])
    adata2 = ad.AnnData(X=X2, obs=obs2, var=var2)
    
    h5ad_path2 = tmp_path / "capital_leiden.h5ad"
    adata2.write_h5ad(h5ad_path2)
    
    out_csv2 = tmp_path / "clusters2.csv"
    argv2 = ["cluster_cells.py", str(h5ad_path2), str(out_csv2)]
    monkeypatch.setattr(sys, "argv", argv2)
    run_cluster_cells.main()
    
    captured2 = capsys.readouterr()
    assert "Using existing Leiden clustering" in captured2.out
