
# Changelog

All notable changes to this project are documented in this file. The format
follows the "Keep a Changelog" convention and is intended to be human
readable.

## [0.1.0] - 2025-10-26

### Added
- TALON validator (`se4s validate`) with schema, coordinate, novelty-set,
    and QC string checks. A SAM round‑trip spot‑check is provided as an
    opt-in flag (`--enable-spotcheck`) and requires `--se4s-root` to locate
    labeled SAMs.
- Counting utilities and CLI (`se4s counts`) to aggregate long‑read support
    by `transcript_ID` and `annot_gene_id`. `write_top_tables` emits top‑N CSVs
    to `outputs/tables` (e.g., `top_transcripts_<DATASET>.csv`,
    `top_genes_<DATASET>.csv`).
- Isoform QC utilities: `isoform_qc_table` and
    `write_isoform_qc_table` generate `outputs/tables/isoform_qc_table.csv`, a
    per‑gene pivot of novelty category counts (columns include Known, NIC,
    NNC, ISM, Antisense, Intergenic, Genomic when present).
- Basic pytest tests covering validator invariants and counting outputs.
- Editable install and console script entry point `se4s` via
    `pyproject.toml`.

### Changed
- Repository reorganized as a minimal Python package under `src/se4s_isoform`
    for clean imports, CLI exposure, and packaging.

### Removed
- None for this initial release.

### Fixed
- Not applicable for this initial tagged release.

### Notes
- v0.1.0 intentionally focuses on reliable validation and simple summaries
    for long‑read outputs. Future work (v0.2) will expand short‑read ingestion,
    mapping utilities, and downstream integration with single‑cell analyses.


    v0.1.0 focuses on reliable validation and simple summaries for long‑read outputs and is intentionally compact to enable rapid iteration toward scRNA ingestion and integration in v0.2.

​