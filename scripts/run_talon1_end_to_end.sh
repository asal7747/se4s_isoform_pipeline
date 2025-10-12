#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# TALON end-to-end: label reads (optional), annotate, validate, and persist
# Tested on Apple Silicon with Terminal “Open using Rosetta” for x86_64 conda
# Outputs: bulk_sc_talon_read_annot.tsv and QC log into repo/outputs
###############################################################################

# -----------------------------
# 0) Project paths and inputs
# -----------------------------
# Root of project (Google Drive path used for persistent inputs/DB)
PROJECT_ROOT="$HOME/Library/CloudStorage/GoogleDrive-asal7747@colorado.edu/My Drive/Fall Semester 2025/SE4S"  # [ref] project tree on host [web:147]
# Repo working copy
REPO_DIR="$PROJECT_ROOT/se4s_isoform_pipeline"  # [ref] repo root [web:147]
# Reference genome for labeling (already used previously)
REF_FA="$PROJECT_ROOT/refs/mm10.fa"  # [ref] mm10 fasta used by talon_label_reads [web:5]
# TALON database (mm10 vM21)
TALON_DB="$PROJECT_ROOT/work/longread/talon/mm10_vM21.db"  # [ref] TALON DB path [web:5]
# Datasets (ENCODE-style basenames)
DATASETS=("ENCFF003OWX" "ENCFF019HRC" "ENCFF669LWV" "ENCFF676BYQ")  # [ref] four bulk datasets [web:5]
# Labeled SAMs (produced earlier; this script will verify presence and can rebuild if needed)
declare -A LABELED_SAMS  # [ref] map dataset -> labeled SAM path [web:5]
for BN in "${DATASETS[@]}"; do
  LABELED_SAMS["$BN"]="$PROJECT_ROOT/work/longread/bulk/$BN/labeled/clean_labeled.sam"  # [ref] expected labeled SAM outputs [web:5]
done

# Local scratch for heavy I/O (non-synced; avoids Drive stalls)
SCRATCH="$HOME/Desktop/talon_scratch"  # [ref] use local disk for temp/out [web:147]
mkdir -p "$SCRATCH" "$SCRATCH/tmp"  # [ref] ensure scratch exists [web:147]

# Repo outputs
mkdir -p "$REPO_DIR/outputs/tables"  # [ref] ensure repo outputs directory [web:147]

# -----------------------------
# 1) Print and verify inputs
# -----------------------------
echo "== Inputs =="
echo "PROJECT_ROOT: $PROJECT_ROOT"  # [ref] project tree summary [web:147]
echo "REPO_DIR:     $REPO_DIR"  # [ref] repo location [web:147]
echo "REF_FA:       $REF_FA"  # [ref] reference genome used in labeling [web:5]
echo "TALON_DB:     $TALON_DB"  # [ref] TALON database [web:5]
for BN in "${DATASETS[@]}"; do
  echo "DATASET: $BN"  # [ref] dataset list [web:5]
  echo "  LABELED_SAM: ${LABELED_SAMS[$BN]}"  # [ref] expected labeled SAM [web:5]
done

# Fail early if DB or reference is missing
[ -f "$TALON_DB" ] || { echo "Missing TALON DB: $TALON_DB" >&2; exit 1; }  # [ref] guard presence [web:5]
[ -f "$REF_FA" ] || echo "WARN: Reference FASTA not found (only needed if re-labeling): $REF_FA"  # [ref] optional if relabeling [web:5]

# Check labeled SAMs; collect missing
MISSING=()
for BN in "${DATASETS[@]}"; do
  SAM="${LABELED_SAMS[$BN]}"
  if [ ! -f "$SAM" ]; then
    echo "MISSING labeled SAM: $SAM"  # [ref] will trigger optional labeling [web:5]
    MISSING+=("$BN")
  else
    ls -lh "$SAM"  # [ref] report size for context [web:147]
  fi
done

# -----------------------------
# 2) Optional: build labeled SAMs via Docker (only if missing)
# -----------------------------
# Requires Docker Desktop installed; runs talon_label_reads inside container and writes back to /data [web:147]
if [ "${#MISSING[@]}" -gt 0 ]; then
  echo "== Labeling missing SAMs in Docker =="
  # Each dataset must have a base SAM to label (clean.sam); rebuild from BAM if needed
  for BN in "${MISSING[@]}"; do
    D="$PROJECT_ROOT/work/longread/bulk/$BN"  # [ref] dataset directory [web:147]
    BASE_SAM="$D/clean.sam"  # [ref] base SAM expected by talon_label_reads [web:5]
    if [ ! -f "$BASE_SAM" ]; then
      # Convert an aligned BAM to SAM if present
      IN="$D/clean.bam"; [ -f "$IN" ] || IN="$D/aligned.bam"  # [ref] choose available BAM [web:5]
      [ -f "$IN" ] || { echo "Missing BAM for $BN at $D"; exit 1; }  # [ref] fail if no BAM [web:5]
      echo "Converting BAM->SAM for $BN"  # [ref] prepare input for labeling [web:5]
      samtools view -h "$IN" > "$BASE_SAM"  # [ref] create SAM with headers [web:5]
    fi
  done

  # talon_label_reads in container; pin platform if needed on Apple Silicon [web:147]
  docker run --rm --platform linux/amd64 -v "$PROJECT_ROOT":/data kwellswrasman/talon bash -lc '
FA="/data/refs/mm10.fa";
for BN in '"$(printf "'%s' " "${MISSING[@]}")"'; do
  D="/data/work/longread/bulk/$BN";
  SAM="$D/clean.sam"; LBL="$D/labeled"; mkdir -p "$LBL";
  talon_label_reads --f "$SAM" --g "$FA" --t 4 --tmpDir "$LBL/tmp" --deleteTmp --o "$LBL/clean";
done'  # [ref] label reads per official guidance [web:5][web:147]
fi

# -----------------------------
# 3) Intel (Rosetta) TALON environment (local x86_64 conda)
# -----------------------------
# This run uses TALON 5.0 from Bioconda inside an x86_64 Python 3.7 environment [web:5]
# Precondition: Terminal launched with “Open using Rosetta” and CONDA_SUBDIR=osx-64 [web:150]
if ! command -v talon >/dev/null 2>&1; then
  echo "ERROR: talon not on PATH; ensure Intel shell + conda env active (see README)" >&2  # [ref] env guard [web:150]
  exit 1
fi

# -----------------------------
# 4) Build TALON config for annotation (host paths)
# -----------------------------
CFG_LOCAL="$SCRATCH/config_local.csv"  # [ref] four-column CSV required by TALON [web:5]
: > "$CFG_LOCAL"  # [ref] truncate/create [web:147]
for BN in "${DATASETS[@]}"; do
  echo "$BN,C2C12_bulk,PacBio-Sequel2,${LABELED_SAMS[$BN]}" >> "$CFG_LOCAL"  # [ref] dataset,desc,platform,path [web:5]
done
echo "Built config: $CFG_LOCAL"  # [ref] report config path [web:147]

# -----------------------------
# 5) TALON annotate to local scratch
# -----------------------------
OUT_PREFIX="$SCRATCH/bulk_run_local"  # [ref] outprefix for TSV and QC [web:5]
LOG="$SCRATCH/bulk_run_local.talon.log"  # [ref] capture log for provenance [web:5]
echo "== Running TALON annotate =="  # [ref] annotate step [web:5]
talon \
  --f "$CFG_LOCAL" \
  --db "$TALON_DB" \
  --build mm10 \
  --threads 1 \
  --cov 0.9 \
  --identity 0.8 \
  --o "$OUT_PREFIX" | tee "$LOG"  # [ref] minimal supported flags for this build [web:5]

# Expected outputs (may be named by tool as *_talon_read_annot.tsv)
TSV_CANDIDATES=(
  "${OUT_PREFIX}_read_annot.tsv"
  "${OUT_PREFIX}_talon_read_annot.tsv"
)  # [ref] handle variant filenames from exporter [web:5]
QC_FILE="${OUT_PREFIX}_QC.log"  # [ref] QC log path [web:5]

TSV_FOUND=""
for f in "${TSV_CANDIDATES[@]}"; do
  if [ -f "$f" ]; then TSV_FOUND="$f"; break; fi  # [ref] pick existing TSV [web:5]
done
if [ -z "${TSV_FOUND:-}" ]; then
  echo "ERROR: read_annot TSV not found; inspect log: $LOG" >&2  # [ref] failure guard [web:5]
  exit 1
fi

# -----------------------------
# 6) Validation summaries (lightweight)
# -----------------------------
echo "== Validation summaries =="# [ref] basic QC on TSV [web:5]
TSV="$TSV_FOUND"  # [ref] resolved TSV path [web:5]
echo "TSV: $TSV"  # [ref] show final TSV path [web:5]
echo "QC:  $QC_FILE"  # [ref] show QC log path [web:5]
wc -l "$TSV"  # [ref] row count (includes header) [web:5]
# Dataset counts
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i} NR>1{c[$(h["dataset"])]++} END{for(k in c) printf "%s\t%d\n", k, c[k]}' "$TSV" | sort  # [ref] per-dataset totals [web:5]
# Novelty distributions
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i} NR>1 && $(h["gene_novelty"])!=""{g[$(h["gene_novelty"])]++} END{for(k in g) printf "gene_novelty\t%s\t%d\n", k, g[k]}' "$TSV" | sort  # [ref] gene novelty [web:5]
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i} NR>1 && $(h["transcript_novelty"])!=""{t[$(h["transcript_novelty"])]++} END{for(k in t) printf "transcript_novelty\t%s\t%d\n", k, t[k]}' "$TSV" | sort  # [ref] transcript novelty [web:5]
# Coordinate sanity (min<=max)
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i} NR>1{ s=$(h["read_start"])+0; e=$(h["read_end"])+0; if((s<e? s:e) <= (s<e? e:s)) ok++ } END{print "min<=max rows:",ok}' "$TSV"  # [ref] orientation-agnostic check [web:5]

# -----------------------------
# 7) Persist into repo
# -----------------------------
FINAL_TSV="$REPO_DIR/outputs/tables/bulk_sc_talon_read_annot.tsv"  # [ref] repo location for TSV [web:147]
cp -p "$TSV" "$FINAL_TSV"  # [ref] copy final TSV [web:147]
cp -p "$QC_FILE" "$REPO_DIR/outputs/bulk_run_local_QC.log"  # [ref] copy QC log [web:147]
echo "Wrote:"  # [ref] report copied artifacts [web:147]
ls -lh "$FINAL_TSV" "$REPO_DIR/outputs/bulk_run_local_QC.log"  # [ref] verify sizes [web:147]

# -----------------------------
# 8) (Optional) Git commit
# -----------------------------
if [ -d "$REPO_DIR/.git" ]; then
  cd "$REPO_DIR"  # [ref] enter repo [web:147]
  git add outputs/tables/bulk_sc_talon_read_annot.tsv outputs/bulk_run_local_QC.log || true  # [ref] stage files [web:147]
  git commit -m "Add TALON read-wise TSV and QC log (reproducible run script included)" || true  # [ref] commit artifacts [web:147]
  echo "Committed changes (if any) to repo"  # [ref] notify commit status [web:147]
fi

echo "DONE: TALON read-wise TSV produced and saved to repo."  # [ref] completion message [web:5][web:147]
