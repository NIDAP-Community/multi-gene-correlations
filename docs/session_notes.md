# Session Notes (append-only)

## 2026-01-09
- Added container/README.md clarifying canonical entrypoint and container distribution.
- Updated README to remove ambiguity about Singularity availability.
- No functional code changes.
- Relocated Dockerfile into container/ and rebuilt image (linux/amd64). Verified container run with TCGA PRAD example after move.
- Added docs/architecture.md outlining components, flow, I/O, and platform notes.
- Renamed core function to `Multi_Gene_Correlations_to_Signature` and updated entrypoint/docs; retest pending after rename.

## 2026-01-12
- Extended `Multi_Gene_Correlations_to_Signature` and `entrypoint.R` with optional clinical correlation support (numeric + categorical metadata), including new CLI flags and outputs for tables/plots.
- Updated README to document the clinical workflow and note the additional artifacts.
- Added a local/core script locator in `entrypoint.R` so `Rscript entrypoint.R` works outside the container; generated synthetic `test_data/` inputs and verified the new clinical correlations by running the entrypoint into `test_output/clinical_demo`.

