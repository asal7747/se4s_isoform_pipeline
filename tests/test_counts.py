import pandas as pd
import pytest

from se4s_isoform.talon_to_counts import (
    talon_counts_by_gene,
    talon_counts_by_transcript,
)


@pytest.fixture
def small_talon_tsv(tmp_path):
    """
    Create a tiny TALON-like TSV suitable for unit tests. Columns include
    `dataset`, `transcript_ID`, and `annot_gene_id` which are consumed by
    the utility functions under test.
    """
    df = pd.DataFrame(
        {
            "dataset": ["ds1", "ds1", "ds2", "ds1", "ds2"],
            "transcript_ID": ["TX1", "TX1", "TX2", "TX3", "TX2"],
            "annot_gene_id": ["G1", "G1", "G2", "G3", "G2"],
        }
    )

    tsv_path = tmp_path / "talon_small.tsv"
    df.to_csv(tsv_path, sep="\t", index=False)
    return str(tsv_path)


def test_counts_tx_nonempty(small_talon_tsv):
    s = talon_counts_by_transcript(small_talon_tsv, dataset="ds1")
    # TX1 appears twice in ds1, so it should be present and integer-typed
    assert "TX1" in s.index
    assert s.dtype.kind in "iu"


def test_counts_gene_nonempty(small_talon_tsv):
    s = talon_counts_by_gene(small_talon_tsv, dataset="ds2")
    # For ds2, G2 has two observations in our mini-TSV
    assert "G2" in s.index
    assert s.loc["G2"] >= 1
    assert s.dtype.kind in "iu"
