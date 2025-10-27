import os

from se4s_isoform.talon_validate import (
    datasets_present,
    min_le_max,
    novelty_sets_ok,
    qc_identity_ok,
    schema_ok,
)

REPO = os.path.expanduser(
    "~/Library/CloudStorage/GoogleDrive-asal7747@colorado.edu/My Drive/Fall Semester 2025/SE4S/se4s_isoform_pipeline"
)
TSV = os.path.join(REPO, "outputs", "tables", "bulk_sc_talon_read_annot.tsv")
QC = os.path.join(REPO, "outputs", "bulk_run_local_QC.log")


def test_smoke_validate():
    assert schema_ok(TSV)
    assert min_le_max(TSV)
    assert novelty_sets_ok(TSV)
    assert qc_identity_ok(QC)
    assert datasets_present(
        TSV, {"ENCFF003OWX", "ENCFF019HRC", "ENCFF669LWV", "ENCFF676BYQ"}
    )
