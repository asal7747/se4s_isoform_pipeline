<<<<<<< HEAD

# se4s_isoform_pipeline

A reproducible Python pipeline to bridge short‑read single‑cell RNA‑seq with long‑read isoform resolution in C2C12 muscle cells. Includes QC, clustering, isoform‑aware proxies, and a TALON-based read‑wise annotation step.

---

## Features
- **Long‑read read‑wise annotation** built on TALON 5.0 with a config‑driven workflow and mm10 vM21 gene models.
- **Containerized labeling option** and local x86_64 conda environment for robust Apple Silicon support.
- **Clear separation of inputs, scripts, and outputs** to enable reruns and downstream analyses.

## External Inputs
- **ENCODE long‑read datasets** (PacBio Iso‑Seq; mouse C2C12):
	- ENCFF003OWX, ENCFF019HRC, ENCFF669LWV, ENCFF676BYQ (stored under `work/longread/bulk/<dataset>/`).
- **Mouse genome (mm10/GRCm38) FASTA** from UCSC bigZips (e.g., `mm10.fa.gz`) with samtools index (`mm10.fa.fai`), placed under `refs/mm10.fa(.gz)`.
- **GENCODE M21 annotation** for mm10 with UCSC names (`gencode.vM21.primary_assembly.annotation_UCSC_names.gtf.gz`) used to initialize and interpret the TALON database.
- **TALON database** for mm10 vM21 (`mm10_vM21.db`) stored under `work/longread/talon/` as the target DB for annotation and exporter steps.

## Requirements
- Docker Desktop on macOS (optional, for containerized labeling and consistent runtimes)
- Bioconda TALON 5.0 under an x86_64 (Rosetta) conda environment with Python 3.7
- TALON’s config‑driven interface and a four‑column CSV listing dataset, description, platform, and labeled SAM paths

## Quickstart: TALON Read‑wise TSV
1. **Prepare labeled SAMs** for each ENCODE dataset (`clean_labeled.sam` under `work/longread/bulk/<dataset>/labeled/`) using `talon_label_reads`, or verify they already exist.
2. **Activate the x86_64 conda env** and run TALON 5.0 with minimal supported flags: `--f`, `--db`, `--build`, `--threads`, `--cov`, `--identity`, `--o`.
3. **Write outprefix and QC** to a local, non‑synced scratch directory for performance, then copy the final TSV and log into `outputs/`.

### Example
- **Labeling (Docker; optional):** Mount project at `/data` and run `talon_label_reads` against `mm10.fa` with temporary space in dataset‑scoped tmp.
- **Annotation (conda; TALON 5.0):** Build `config.csv` with absolute host paths to labeled SAMs, run `talon` to produce `OUTPREFIX_talon_read_annot.tsv` and QC log, and validate dataset and novelty counts.

## Reproducibility
- Use `scripts/run_talon_end_to_end.sh` to verify inputs, optionally label reads in Docker, run TALON annotate locally, validate outputs, and copy artifacts into the repo.
- For Apple Silicon, launch Terminal “Open using Rosetta,” pin `CONDA_SUBDIR=osx‑64`, install TALON 5.0 from Bioconda, and execute the script to reproduce the TSV.

## Repository Layout
- `scripts/`: Runnable pipeline and utility scripts, including `run_talon_end_to_end.sh` for the TALON step.
- `docs/`: `Inputs.md` listing all externally sourced inputs and `ReproducibleSteps.md` with exact commands and environment details.
- `outputs/`: QC logs and tables; `outputs/tables/bulk_sc_talon_read_annot.tsv` is the unified read‑wise annotation for downstream analysis.

## Notes
- TALON’s config‑driven approach standardizes multi‑dataset runs and supports per‑read novelty calls aligned to mm10 GENCODE M21.
- UCSC bigZips provides mm10 FASTA and related indexes, which should be placed under `refs/` to ensure consistent labeling and genomic mapping.

## License
Project license is included in `LICENSE` and covers the pipeline code and repository assets; external inputs remain under their respective source licenses and data use policies.
=======
# se4s_isoform_pipeline
Reproducible Python pipeline bridging short‑read single‑cell RNA‑seq and long‑read isoform resolution in C2C12 muscle cells, with QC, clustering, isoform‑aware proxies, and benchmarks to guide when long reads are needed.
>>>>>>> origin/main
