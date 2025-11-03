#!/usr/bin/env bash
set -euo pipefail

# Download ENCODE single-cell h5ad files into a target directory
OUT_DIR="${1:-outputs/anndata}"
mkdir -p "$OUT_DIR"

echo "Downloading ENCODE scRNA-seq h5ad files to $OUT_DIR ..."

# Short-read shallow (9k cells)
curl -L --fail --retry 3 --create-dirs \
  -o "$OUT_DIR/short_shallow.h5ad" \
  "https://www.encodeproject.org/files/ENCFF914DEE/@@download/ENCFF914DEE.h5ad"

# Short-read deep (1k cells)
curl -L --fail --retry 3 \
  -o "$OUT_DIR/short_deep.h5ad" \
  "https://www.encodeproject.org/files/ENCFF418TAK/@@download/ENCFF418TAK.h5ad"

# Long-read transcript counts
curl -L --fail --retry 3 \
  -o "$OUT_DIR/long_transcript.h5ad" \
  "https://www.encodeproject.org/files/ENCFF513MSS/@@download/ENCFF513MSS.h5ad"

# Long-read gene counts
curl -L --fail --retry 3 \
  -o "$OUT_DIR/long_gene.h5ad" \
  "https://www.encodeproject.org/files/ENCFF508RNZ/@@download/ENCFF508RNZ.h5ad"

echo "Downloaded 4 files to $OUT_DIR"
ls -lh "$OUT_DIR"
