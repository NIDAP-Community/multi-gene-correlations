#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
    cat <<'USAGE'
Usage: run_clinical_plots.sh <path-to-sif> [entrypoint args]

Environment variables:
  DATA_DIR   Directory containing counts/metadata files (default: $PWD/data)
  OUTPUT_DIR Directory for analysis outputs (default: $PWD/output)
  EXTRA_BINDS Optional additional Singularity bind mounts (comma-separated)

Example:
  DATA_DIR=$PWD/test_data OUTPUT_DIR=$PWD/test_output ./run_clinical_plots.sh \
      ~/images/multi-gene-correlations.sif \
      --counts /data/counts.tsv \
      --metadata /data/metadata.csv \
      --gene_column Gene --sample_column SampleID \
      --genes "GATM,GAMT,CKM,SLC6A8" \
      --category_column Tissue --categories "Normal,Tumor" \
      --signature_name Creatine --signature_genes "PTEN,SPRY2" \
      --clinical_columns "Stage,PSA" --output_dir /output
USAGE
    exit 1
fi

IMAGE_PATH="$1"
shift

if [[ ! -f "$IMAGE_PATH" ]]; then
    echo "Error: SIF image not found at $IMAGE_PATH" >&2
    exit 2
fi

DATA_DIR="${DATA_DIR:-$PWD/data}"
OUTPUT_DIR="${OUTPUT_DIR:-$PWD/output}"
EXTRA_BINDS="${EXTRA_BINDS:-}"

mkdir -p "$OUTPUT_DIR"

BIND_ARGS=("${DATA_DIR}:/data"
           "${OUTPUT_DIR}:/output")

if [[ -n "$EXTRA_BINDS" ]]; then
    IFS=',' read -r -a extra_array <<< "$EXTRA_BINDS"
    for item in "${extra_array[@]}"; do
        BIND_ARGS+=("$item")
    done
fi

BIND_FLAGS=()
for bind in "${BIND_ARGS[@]}"; do
    BIND_FLAGS+=("--bind" "$bind")
done

singularity run \
    --pwd /output \
    "${BIND_FLAGS[@]}" \
    --env R_LIBS_USER=/usr/local/lib/R/site-library \
    "$IMAGE_PATH" \
    "$@"
