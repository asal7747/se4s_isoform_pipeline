# test run full pipeline to make sure it makes expected output files
test -e ssshtest || curl -s -O https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest
. ssshtest

# Input and output paths
talon_tsv="tests/bulk_sc_talon_read_annot.tsv"
scrna_h5ad="tests/test_data_cluster_cells.h5ad"
output_dir="tests/func/full_pipeline_output"

# Skip test if required input files don't exist
if [ ! -f "$talon_tsv" ]; then
    echo "SKIP: TALON TSV not found at $talon_tsv"
    exit 0
fi

if [ ! -f "$scrna_h5ad" ]; then
    echo "SKIP: scRNA h5ad not found at $scrna_h5ad"
    exit 0
fi

# Expected output files
qc_h5ad="$output_dir/test_data_cluster_cells_qc.h5ad"
clusters_csv="$output_dir/tables/cell_clusters.csv"
symbol_map_csv="$output_dir/tables/talon_scrna_symbol_map.csv"
proxies_csv="$output_dir/tables/isoform_proxies.csv"
umap_plot="$output_dir/plots/umap_clusters.png"

# Clean up any existing output files before test
rm -rf "$output_dir"

# test 1: run full pipeline and check all expected outputs
run test_run_full_pipeline python3 scripts/run_full_pipeline.py \
"$talon_tsv" \
"$scrna_h5ad" \
"$output_dir"

# Check exit code
assert_exit_code 0

# Check that QC'd h5ad was created
assert_equal "$qc_h5ad" "$( ls $qc_h5ad )"

# Check that cluster assignments CSV was created
assert_equal "$clusters_csv" "$( ls $clusters_csv )"

# Check that symbol map CSV was created
assert_equal "$symbol_map_csv" "$( ls $symbol_map_csv )"

# Check that isoform proxies CSV was created
assert_equal "$proxies_csv" "$( ls $proxies_csv )"

# Check that UMAP plot was created
assert_equal "$umap_plot" "$( ls $umap_plot )"

# Clean up
rm -rf "$output_dir"