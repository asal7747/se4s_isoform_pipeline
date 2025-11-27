#!/usr/bin/env bash
set -euo pipefail

# Download ENCODE single-cell h5ad files into a target directory
OUT_DIR="${1:-outputs/anndata}"
mkdir -p "$OUT_DIR"

echo "Downloading ENCODE scRNA-seq h5ad files to $OUT_DIR ..."

# short-read shallow Split-seq single-cell myoblast (9k cells)
curl -L --fail --retry 3 --create-dirs \
  -o "$OUT_DIR/short_shallow.h5ad" \
  "https://www.encodeproject.org/files/ENCFF914DEE/@@download/ENCFF914DEE.h5ad"
# short-read shallow Split-seq single-nucleus myoblast (9k cells)
curl -L --fail --retry 3 \
  -o "$OUT_DIR/short_shallow_nuc.h5ad" \
  "https://www.encodeproject.org/files/ENCFF129LKB/@@download/ENCFF129LKB.h5ad"

# Short-read deep (1k cells)
curl -L --fail --retry 3 \
  -o "$OUT_DIR/short_deep.h5ad" \
  "https://www.encodeproject.org/files/ENCFF418TAK/@@download/ENCFF418TAK.h5ad"
# Short-read deep single-nucleus (1k cells)
curl -L --fail --retry 3 \
  -o "$OUT_DIR/short_deep_nuc.h5ad" \
  "https://www.encodeproject.org/files/ENCFF755KGW/@@download/ENCFF755KGW.h5ad"

# Long-read Split-seq single-cell myoblast data
# Long-read transcript counts
curl -L --fail --retry 3 \
  -o "$OUT_DIR/long_transcript.h5ad" \
  "https://www.encodeproject.org/files/ENCFF513MSS/@@download/ENCFF513MSS.h5ad"

# Long-read gene counts
curl -L --fail --retry 3 \
  -o "$OUT_DIR/long_gene.h5ad" \
  "https://www.encodeproject.org/files/ENCFF508RNZ/@@download/ENCFF508RNZ.h5ad"

# long-read Split-seq single-nucleus myoblast data
# Long-read nucleus transcript counts
curl -L --fail --retry 3 \
  -o "$OUT_DIR/long_nuc_transcript.h5ad" \
  "https://www.encodeproject.org/files/ENCFF033NVO/@@download/ENCFF033NVO.h5ad"
# Long-read nucleus gene counts
curl -L --fail --retry 3 \
  -o "$OUT_DIR/long_nuc_gene.h5ad" \
  "https://www.encodeproject.org/files/ENCFF301DZH/@@download/ENCFF301DZH.h5ad"

# short-read shallow Split-seq single-cell myotube (9k cells)
curl -L --fail --retry 3 \
  -o "$OUT_DIR/short_shallow_myotube.h5ad" \
  "https://www.encodeproject.org/files/ENCFF545EUV/@@download/ENCFF545EUV.h5ad"
# short-read deep Split-seq single-cell myotube (1k cells)
curl -L --fail --retry 3 \
  -o "$OUT_DIR/short_deep_myotube.h5ad" \
  "https://www.encodeproject.org/files/ENCFF588HWS/@@download/ENCFF588HWS.h5ad"

# long-read Split-seq single-cell myotube data
# Long-read myotube transcript counts
curl -L --fail --retry 3 \
  -o "$OUT_DIR/long_myotube_transcript.h5ad" \
  "https://www.encodeproject.org/files/ENCFF362NGZ/@@download/ENCFF362NGZ.h5ad"
# Long-read myotube gene counts
curl -L --fail --retry 3 \
  -o "$OUT_DIR/long_myotube_gene.h5ad" \
  "https://www.encodeproject.org/files/ENCFF049LST/@@download/ENCFF049LST.h5ad"

echo "Downloaded 12 files to $OUT_DIR"
ls -lh "$OUT_DIR"
