import pandas as pd

from se4s_isoform.talon_validate import (
    datasets_present,
    min_le_max,
    novelty_sets_ok,
    qc_identity_ok,
    schema_ok,
)


def test_smoke_validate(tmp_path):
    # Create a small TALON-like TSV satisfying required schema and values
    df = pd.DataFrame(
        {
            "read_name": ["r1", "r2", "r3"],
            "dataset": ["ds1", "ds2", "ds1"],
            "genome_build": ["mm10", "mm10", "mm10"],
            "chrom": ["chr1", "chr2", "chr1"],
            "read_start": [100, 200, 150],
            "read_end": [200, 300, 250],
            "strand": ["+", "-", "+"],
            "gene_ID": ["G1", "G2", "G1"],
            "transcript_ID": ["T1", "T2", "T3"],
            "gene_novelty": ["Known", "NIC", "NNC"],
            "transcript_novelty": ["Known", "NIC", "NNC"],
            "fraction_As": [0.0, 0.1, 0.0],
        }
    )

    tsv_path = tmp_path / "talon_test.tsv"
    df.to_csv(tsv_path, sep="\t", index=False)

    # Create a QC log file containing the required identity string
    qc_path = tmp_path / "bulk_run_test_QC.log"
    qc_path.write_text(
        "Some header\nMin read identity to reference: 0.800000\nMore info\n"
    )

    assert schema_ok(str(tsv_path))
    assert min_le_max(str(tsv_path))
    assert novelty_sets_ok(str(tsv_path))
    assert qc_identity_ok(str(qc_path))
    assert datasets_present(str(tsv_path), {"ds1", "ds2"})
