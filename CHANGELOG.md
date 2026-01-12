# Changelog

All notable changes to this project will be documented in this file.

## Unreleased
- Initial repository scaffolding and documentation standards.
- Move Dockerfile to container/ and add container README; update build instructions to point to new path.

## 2026-01-12
- Added optional clinical correlation workflow (numeric + categorical metadata) with new CLI flags and outputs for signature/gene associations.
- Enabled local `Rscript entrypoint.R` runs by auto-locating the core script; added synthetic `test_data/` fixtures and verified both local and Docker executions.
- Built and published `ghcr.io/nidap-community/multi-gene-correlations:latest` (linux/amd64) from container/Dockerfile to capture the new features.
- Added clinical signature- and gene-level bar plots so clinical outputs now mirror the categorical bar plots available for gene vs signature analyses.
- Hardened the container build: preinstalled tidyverse dependencies, set `R_LIBS_USER=/usr/local/lib/R/site-library`, bundled the synthetic smoke-test data, and wired `scripts/generate_clinical_correlation_plots.R` into the Docker build.
- Added `container/run_clinical_plots.sh` helper plus README updates documenting the standard `/data` and `/output` binds for Singularity users.
- Clinical numeric scatter plots now show **gene expression vs clinical measurements** (with Pearson r and p in facet labels) instead of signature-only scatters, matching the intended visualization requirements.
- Per-category signature scatter plots now annotate each gene facet with the corresponding Pearson r/p so reviewers can tie the bar plots back to the underlying point clouds.
