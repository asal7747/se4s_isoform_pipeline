"""se4s_isoform package: small utilities to summarize TALON read-wise TSVs.

This package exposes validator and counting helpers used by the `se4s` CLI.
"""

from .talon_to_counts import (
    isoform_qc_table,
    talon_counts_by_gene,
    talon_counts_by_transcript,
    write_isoform_qc_table,
    write_top_tables,
)
from .talon_validate import (
    datasets_present,
    min_le_max,
    novelty_sets_ok,
    qc_identity_ok,
    schema_ok,
    spotcheck_reads,
)

__all__ = [
    "talon_counts_by_transcript",
    "talon_counts_by_gene",
    "write_top_tables",
    "isoform_qc_table",
    "write_isoform_qc_table",
    "schema_ok",
    "min_le_max",
    "datasets_present",
    "novelty_sets_ok",
    "qc_identity_ok",
    "spotcheck_reads",
]
