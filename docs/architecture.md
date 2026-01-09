# Architecture

## Overview
Containerized R CLI that computes gene-to-signature correlations and produces plots/tables. Canonical entrypoint is `entrypoint.R`, and core analysis lives in `Multi_Gene_Correlations_to_Signature.R`.

## Execution Flow
1) `entrypoint.R` parses CLI flags, loads counts/metadata, resolves columns, and calls the core function.
2) `Multi_Gene_Correlations_to_Signature_IODC_Beta` performs correlation analysis, builds plots (bar, scatter, heatmap), and returns data frames + ggplot/ComplexHeatmap objects.
3) `entrypoint.R` writes outputs under the user-specified `output_dir` (tables/, barplots/, scatter/, heatmap/, metadata/).

## Container Layout
- `container/Dockerfile`: builds runtime image (rocker/r-ver:4.5.2) with required R packages; copies root scripts into `/app` and sets `entrypoint.R` as entrypoint.
- `container/README.md`: documents CAF model, build targets (`linux/amd64`), and runtime guidance.

## Inputs and Outputs
- Inputs: normalized counts (CSV/TSV), metadata (CSV/TSV), gene lists, signature genes, category labels.
- Outputs: `tables/pathways.csv`, `metadata/run_parameters.json`, `barplots/*.png`, `scatter/*.png`, `heatmap/heatmap.png`.

## Platform Notes
- Canonical build/target: `linux/amd64` (Biowulf-compatible).
- Apple Silicon: build with `--platform linux/amd64`.
# Architecture Overview

## Purpose
This repository implements a containerized analytical tool for computing
gene correlation–based signatures from expression data.

The goal is reproducible, non-interactive execution suitable for local,
server, and HPC environments.

---

## Execution Model
This tool follows the **Container-as-a-Function (CAF)** model:

outputs = container(inputs, configuration)

Properties:
- single canonical entrypoint (`entrypoint.R`)
- non-interactive execution
- explicit inputs via CLI arguments
- explicit outputs written to a user-specified directory
- no hidden state

---

## Key Components

- `entrypoint.R`
  - Canonical entrypoint
  - Parses CLI arguments
  - Orchestrates execution
  - Performs minimal validation
  - Delegates computation to core functions

- `Multi_Gene_Correlations_to_Signature.R`
  - Core computational logic
  - No side effects outside provided inputs/outputs

- `container/`
  - Container specification and runtime wrappers
  - Docker image is the canonical runtime artifact
  - Singularity/Apptainer supported via Docker image conversion

- `docs/`
  - `decision_log.md`: architectural and workflow decisions (append-only)
  - `session_notes.md`: development progress notes (append-only)

---

## Data Flow

1. Inputs provided via CLI arguments
2. `entrypoint.R` validates inputs and configuration
3. Core computation performed in `Multi_Gene_Correlations_to_Signature.R`
4. Results written to structured subdirectories under `output_dir`
   (e.g., tables/, plots/, metadata/)

---

## Assumptions and Constraints

- Target runtime architecture: linux/amd64 (x86_64)
- Primary HPC target: Biowulf
- Container binaries (.sif) are not committed to git
- Docker image published to GHCR is the canonical artifact
