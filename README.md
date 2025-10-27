# se4s_isoform_pipeline

A reproducible Python pipeline to bridge short‑read single‑cell RNA‑seq with long‑read isoform resolution in C2C12 muscle cells. Includes QC, clustering, isoform‑aware proxies, and a TALON-based read‑wise annotation step.

---

## Features
- **Long‑read read‑wise annotation** built on TALON 5.0 with a config‑driven workflow and mm10 vM21 gene models.
- **Containerized labeling option** and local x86_64 conda environment for robust Apple Silicon support.
- **Clear separation of inputs, scripts, and outputs** to enable reruns and downstream analyses.

# se4s_isoform_pipeline

A reproducible Python pipeline to bridge short‑read single‑cell RNA‑seq with
long‑read isoform resolution in C2C12 muscle cells. Includes QC, clustering,
isoform‑aware proxies, and a TALON-based read‑wise annotation step.

---

## Features
- **Long‑read read‑wise annotation** built on TALON 5.0 with a config‑driven
	workflow and mm10 vM21 gene models.
- **Containerized labeling option** and local x86_64 conda environment for
	robust Apple Silicon support.
- **Clear separation of inputs, scripts, and outputs** to enable reruns and
	downstream analyses.

## External Inputs
- **ENCODE long‑read datasets** (PacBio Iso‑Seq; mouse C2C12):
	- ENCFF003OWX, ENCFF019HRC, ENCFF669LWV, ENCFF676BYQ (stored under
		`work/longread/bulk/<dataset>/`).
- **Mouse genome (mm10/GRCm38) FASTA** from UCSC bigZips (e.g., `mm10.fa.gz`)
	with samtools index (`mm10.fa.fai`), placed under `refs/mm10.fa(.gz)`.
- **GENCODE M21 annotation** for mm10 with UCSC names
	(`gencode.vM21.primary_assembly.annotation_UCSC_names.gtf.gz`) used to
	initialize and interpret the TALON database.
- **TALON database** for mm10 vM21 (`mm10_vM21.db`) stored under
	`work/longread/talon/` as the target DB for annotation and exporter steps.

## Requirements
- Docker Desktop on macOS (optional, for containerized labeling and
	consistent runtimes)
- Bioconda TALON 5.0 under an x86_64 (Rosetta) conda environment with Python
	3.7
- TALON’s config‑driven interface and a four‑column CSV listing dataset,
	description, platform, and labeled SAM paths

## Quickstart: TALON Read‑wise TSV
1. **Prepare labeled SAMs** for each ENCODE dataset (`clean_labeled.sam`
	 under `work/longread/bulk/<dataset>/labeled/`) using `talon_label_reads`, or
	 verify they already exist.
2. **Activate the x86_64 conda env** and run TALON 5.0 with minimal supported
	 flags: `--f`, `--db`, `--build`, `--threads`, `--cov`, `--identity`,
	 `--o`.
3. **Write outprefix and QC** to a local, non‑synced scratch directory for
	 performance, then copy the final TSV and log into `outputs/`.

### Example
- **Labeling (Docker; optional):** Mount project at `/data` and run
	`talon_label_reads` against `mm10.fa` with temporary space in dataset‑scoped
	tmp.
- **Annotation (conda; TALON 5.0):** Build `config.csv` with absolute host
	paths to labeled SAMs, run `talon` to produce `OUTPREFIX_talon_read_annot.tsv`
	# se4s_isoform_pipeline

	A small, reproducible pipeline that validates TALON long‑read annotations, summarizes isoforms, and ships a simple CLI for quick exploration; built for C2C12 with room to plug in scRNA‑seq next.

	## Why it exists

	Show what long reads add beyond a standard short‑read baseline by validating TALON outputs and producing isoform summaries that are easy to reuse.

	## Install

	Python environment:

	```bash
	python3 -m venv .venv && source .venv/bin/activate
	pip install -e .
	```

	## Quickstart

	Validate TALON outputs:

	```bash
	se4s validate --tsv outputs/tables/bulk_sc_talon_read_annot.tsv --qc outputs/bulk_run_local_QC.log
	```

	Optional SAM spot‑check (when SAMs exist):

	```bash
	se4s validate --tsv ... --qc ... --enable-spotcheck --se4s-root "/path/to/SE4S"
	```

	Summarize counts:

	```bash
	se4s counts --tsv outputs/tables/bulk_sc_talon_read_annot.tsv --out outputs/tables --dataset ENCFF003OWX --top 50
	```

	Isoform QC table (per‑gene novelty pivot):

	```bash
	python - <<'PY'
	from se4s_isoform.talon_to_counts import write_isoform_qc_table
	write_isoform_qc_table('outputs/tables/bulk_sc_talon_read_annot.tsv','outputs/tables')
	PY
	```

	## Outputs

	`outputs/bulk_run_local_QC.log`
	: TALON run parameters and QC summary.

	`outputs/tables/top_transcripts_<DATASET>.csv` and `outputs/tables/top_genes_<DATASET>.csv`
	: top isoforms/genes by read support.

	`outputs/tables/isoform_qc_table.csv`
	: per‑gene counts of novelty classes (e.g., Known, NIC, NNC, ISM).

	## What’s included (v0.1)

	- CLI: `se4s validate` (schema/coords/novelty/QC; optional spot‑check), `se4s counts` (top lists, CSVs).
	- Library: `talon_to_counts` (counts by transcript/gene, isoform_qc_table), `talon_validate` (core checks).
	- Tests: lightweight `pytest` smoke tests for validator and counts.

	## Minimal requirements

	- TALON‑produced TSV and QC log for validation/summaries.
	- Optional: labeled SAMs only if enabling spot‑check.

	## Notes / next steps

	- Keep large artifacts (TSV/H5/H5AD) out of Git; small CSVs are fine to track.
	- CI: run `pytest` on push/PR (recommended).
	- See `configs/example.yml` for path examples you can adapt.

	## License

	Project license is included in `LICENSE` and covers the pipeline code and repository assets; external inputs remain under their respective source licenses and data use policies.


