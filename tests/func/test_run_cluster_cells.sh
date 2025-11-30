# test run_cluster_cells to make sure it makes expected default output file

test -e ssshtest || curl -s -O https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest
. ssshtest

# Input and output paths
qc_h5ad="tests/test_data_cluster_cells.h5ad"
out_csv="tests/func/cluster_cells_output.csv"
test_output_dir="tests/func"
default_umap="$test_output_dir/umap_clusters.png"

# Skip tests if required input file doesn't exist
if [ ! -f "$qc_h5ad" ]; then
    echo "SKIP: Test data not found at $qc_h5ad"
    exit 0
fi

# Clean up any existing output files before test
rm -f "$out_csv" "$default_umap"

# test 1: test run_cluster_cells creates CSV and default UMAP
run test_run_cluster_cells python3 scripts/run_cluster_cells.py \
"$qc_h5ad" \
"$out_csv" \
--output-dir "$test_output_dir"
assert_exit_code 0
assert_equal "$out_csv" "$( ls $out_csv )"
assert_equal "$default_umap" "$( ls $default_umap )"

# test 2: test run_cluster_cells with t-SNE plotting
out_csv_2="tests/func/cluster_cells_output_tsne.csv"
tsne_plot="$test_output_dir/tsne_leiden.png"

# Clean up any existing output files before test 2
rm -f "$out_csv_2" "$tsne_plot"

run test_run_cluster_cells_tsne python3 scripts/run_cluster_cells.py \
"$qc_h5ad" \
"$out_csv_2" \
--output-dir "$test_output_dir" \
--tsne \
--save-figures

assert_exit_code 0
assert_equal "$out_csv_2" "$( ls $out_csv_2 )"
assert_equal "$tsne_plot" "$( ls $tsne_plot )"

# Clean up all test outputs
rm -f "$out_csv" "$default_umap" "$out_csv_2" "$tsne_plot"