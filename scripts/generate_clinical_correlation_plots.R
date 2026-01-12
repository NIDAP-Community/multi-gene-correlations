#!/usr/bin/env Rscript

# Smoke test for the clinical correlation workflow. This script runs the
# canonical entrypoint against bundled synthetic data to ensure all runtime
# dependencies are installed before publishing the container image.

suppressPackageStartupMessages({
    library(jsonlite)
})

counts_path <- Sys.getenv("SMOKE_COUNTS", "/app/test_data/counts.tsv")
metadata_path <- Sys.getenv("SMOKE_METADATA", "/app/test_data/metadata.csv")
output_dir <- Sys.getenv("SMOKE_OUTPUT", "/tmp/smoke_clinical")

if (!file.exists(counts_path)) {
    stop(sprintf("Counts file missing at %s", counts_path))
}
if (!file.exists(metadata_path)) {
    stop(sprintf("Metadata file missing at %s", metadata_path))
}
if (dir.exists(output_dir)) {
    unlink(output_dir, recursive = TRUE, force = TRUE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

entry_args <- c(
    "/app/entrypoint.R",
    "--counts", counts_path,
    "--metadata", metadata_path,
    "--gene_column", "Gene",
    "--sample_column", "SampleID",
    "--genes", "GATM,GAMT,CKM,SLC6A8",
    "--category_column", "Tissue",
    "--categories", "Normal,Tumor",
    "--signature_name", "Creatine",
    "--signature_genes", "PTEN,SPRY2",
    "--clinical_columns", "Stage,PSA",
    "--clinical_use_signature", "true",
    "--clinical_use_genes", "true",
    "--output_dir", output_dir
)

status <- system2("Rscript", args = entry_args)
if (!identical(status, 0L)) {
    stop("Smoke test failed when running entrypoint.")
}

required_files <- c(
    file.path(output_dir, "tables", "pathways.csv"),
    file.path(output_dir, "tables", "clinical_signature_correlations.csv"),
    file.path(output_dir, "tables", "clinical_gene_correlations.csv"),
    file.path(output_dir, "clinical", "clinical_bar_signature_vs_clinical.png")
)

missing <- required_files[!file.exists(required_files)]
if (length(missing) > 0) {
    stop(sprintf("Smoke test missing expected outputs: %s", paste(missing, collapse = ", ")))
}

cat("Smoke test completed successfully. Outputs stored in", output_dir, "\n")
