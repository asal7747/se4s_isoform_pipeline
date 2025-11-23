import os
from pathlib import Path

import pytest

from se4s_isoform.talon_to_counts import (
    talon_counts_by_gene,
    talon_counts_by_transcript,
)

REPO = os.path.expanduser(
    "~/Library/CloudStorage/GoogleDrive-asal7747@colorado.edu/My Drive/Fall Semester 2025/SE4S/se4s_isoform_pipeline"
)
TSV = os.path.join(REPO, "outputs", "tables", "bulk_sc_talon_read_annot.tsv")


pytestmark = pytest.mark.skipif(
    not Path(TSV).is_file(),
    reason="Requires local TALON TSV (bulk_sc_talon_read_annot.tsv) not in repo",
)


def test_counts_tx_nonempty():
    s = talon_counts_by_transcript(TSV, dataset="ENCFF003OWX")
    assert len(s) > 0 and s.dtype.kind in "iu"


def test_counts_gene_nonempty():
    s = talon_counts_by_gene(TSV, dataset="ENCFF019HRC")
    assert len(s) > 0 and s.dtype.kind in "iu"
