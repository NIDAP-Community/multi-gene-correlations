#!/usr/bin/env Rscript
suppressMessages({
    lib_candidates <- c("/usr/lib/R/site-library", "/usr/lib/R/library")
    lib_candidates <- lib_candidates[dir.exists(lib_candidates)]
    if (length(lib_candidates) > 0) {
        .libPaths(unique(c(.libPaths(), lib_candidates)))
    }
})

suppressPackageStartupMessages({
    # jsonlite optional; fallback to RDS when unavailable
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
        message("jsonlite not available; will save run_parameters as RDS.")
    }
    library(ggplot2)
    library(ComplexHeatmap)
})

# Helpers
read_table_auto <- function(path) {
    ext <- tolower(sub(".*\\.", "", path))
    if (ext %in% c("tsv", "txt")) {
        utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
    } else {
        utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    }
}

parse_list <- function(x) {
    if (is.null(x) || length(x) == 0 || x == "") {
        return(character())
    }
    # If x is a readable file, load lines
    if (file.exists(x)) {
        vals <- readLines(x, warn = FALSE)
        return(trimws(vals[vals != ""]))
    }
    trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

parse_bool <- function(x, default = FALSE) {
    if (is.null(x) || x == "") {
        return(default)
    }
    t <- tolower(x)
    if (t %in% c("true", "t", "1", "yes", "y")) {
        return(TRUE)
    }
    if (t %in% c("false", "f", , "no", "n")) {
        return(FALSE)
    }
    default
}

ensure_dir <- function(path) {
    if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

safe_name <- function(x) {
    gsub("[^-A-Za-z0-9_]+", "_", x)
}

resolve_column <- function(cols, desired) {
    if (is.null(desired) || !nzchar(desired)) {
        return(desired)
    }
    hit <- which(tolower(cols) == tolower(desired))
    if (length(hit) > 0) {
        return(cols[hit[1]])
    }
    desired
}

parse_kv_args <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    out <- list()
    i <- 1
    while (i <= length(args)) {
        a <- args[[i]]
        if (startsWith(a, "--")) {
            a2 <- sub("^--", "", a)
            if (grepl("=", a2, fixed = TRUE)) {
                parts <- strsplit(a2, "=", fixed = TRUE)[[1]]
                key <- parts[1]
                val <- paste(parts[-1], collapse = "=")
                out[[key]] <- val
            } else {
                key <- a2
                val <- if (i + 1 <= length(args) && !startsWith(args[[i + 1]], "--")) args[[i + 1]] else ""
                out[[key]] <- val
                if (val != "") i <- i + 1
            }
        }
        i <- i + 1
    }
    out
}

# Defaults and parse
opt <- modifyList(list(
    counts = "",
    metadata = "",
    gene_column = "",
    sample_column = "",
    samples = "",
    genes = "",
    category_column = "",
    categories = "",
    signature_name = "Signature",
    signature_genes = "",
    output_dir = "/output",
    barplots_dir = "",
    heatmap_dir = "",
    scatter_dir = "",
    tables_dir = "",
    margin_adjust = "-11",
    show_row_dendrogram = "false",
    heatmap_width = "6",
    add_text_to_heatmap = "true",
    rename_categories = "",
    sum_duplicates = "false"
), parse_kv_args())

# Validate required inputs
required <- c("counts", "metadata", "gene_column", "sample_column", "genes", "category_column", "categories", "signature_genes")
missing <- required[vapply(required, function(k) !nzchar(opt[[k]] %||% ""), logical(1))]
if (length(missing) > 0) {
    cat("Missing required options:\n")
    cat(paste0("  --", missing, "\n"), sep = "")
    quit(status = 2)
}

# Resolve lists and flags
samples_to_include <- parse_list(opt$samples)
genes <- parse_list(opt$genes)
categories <- parse_list(opt$categories)
signature_genes <- parse_list(opt$signature_genes)
rename_categories <- parse_list(opt$rename_categories)
show_row_dendrogram <- parse_bool(opt$show_row_dendrogram, FALSE)
add_text_to_heatmap <- parse_bool(opt$add_text_to_heatmap, TRUE)
sum_duplicates <- parse_bool(opt$sum_duplicates, FALSE)

# Prepare output folders
base_out <- opt$output_dir
bar_out <- if (nzchar(opt$barplots_dir)) opt$barplots_dir else file.path(base_out, "barplots")
heat_out <- if (nzchar(opt$heatmap_dir)) opt$heatmap_dir else file.path(base_out, "heatmap")
scatter_out <- if (nzchar(opt$scatter_dir)) opt$scatter_dir else file.path(base_out, "scatter")
table_out <- if (nzchar(opt$tables_dir)) opt$tables_dir else file.path(base_out, "tables")
meta_out <- file.path(base_out, "metadata")

lapply(c(base_out, bar_out, heat_out, scatter_out, table_out, meta_out), ensure_dir)

# Load data
counts <- read_table_auto(opt$counts)
meta <- read_table_auto(opt$metadata)

# Resolve column names case-insensitively to actual headers
gene_col <- resolve_column(colnames(counts), opt$gene_column)
sample_col <- resolve_column(colnames(meta), opt$sample_column)
category_col <- resolve_column(colnames(meta), opt$category_column)

# Load function
source("/app/Multi_Gene_Correlations_to_Signature.R")

# Run
res <- Multi_Gene_Correlations_to_Signature_IODC_Beta(
    normalized_counts = counts,
    sample_metadata = meta,
    gene_column = gene_col,
    sample_column = sample_col,
    samples_to_include = if (length(samples_to_include) == 0) NULL else samples_to_include,
    genes = genes,
    category_column = category_col,
    categories = categories,
    signature_name = opt$signature_name,
    signature_genes = signature_genes,
    margin_adjust = as.numeric(opt$margin_adjust),
    show_row_dendrogram = show_row_dendrogram,
    heatmap_width = as.numeric(opt$heatmap_width),
    add_text_to_heatmap = add_text_to_heatmap,
    rename_categories = rename_categories,
    sum_duplicates = sum_duplicates
)

# Save tables
pathways_path <- file.path(table_out, "pathways.csv")
write.csv(res$pathways, pathways_path, row.names = FALSE)

# Save run parameters (JSON if jsonlite present, else RDS)
params_json <- file.path(meta_out, "run_parameters.json")
params_rds <- file.path(meta_out, "run_parameters.rds")
if (requireNamespace("jsonlite", quietly = TRUE)) {
    write(jsonlite::toJSON(res$run_parameters, pretty = TRUE, auto_unbox = TRUE), params_json)
} else {
    saveRDS(res$run_parameters, params_rds)
}

# Save barplots
if (!is.null(res$plots$barplots) && length(res$plots$barplots) > 0) {
    for (nm in names(res$plots$barplots)) {
        plt <- res$plots$barplots[[nm]]
        if (inherits(plt, "ggplot")) {
            fn <- file.path(bar_out, paste0("barplot_", safe_name(nm), ".png"))
            ggplot2::ggsave(fn, plt, width = 10, height = 8, dpi = 150)
        }
    }
}

# Save heatmap
if (!is.null(res$plots$heatmap)) {
    fn <- file.path(heat_out, "heatmap.png")
    png(fn, width = 1600, height = 2000, res = 150)
    try(
        {
            ComplexHeatmap::draw(res$plots$heatmap)
        },
        silent = TRUE
    )
    dev.off()
}

# Save scatter plots
if (!is.null(res$plots$scatter) && length(res$plots$scatter) > 0) {
    for (nm in names(res$plots$scatter)) {
        plt <- res$plots$scatter[[nm]]
        if (inherits(plt, "ggplot")) {
            fn <- file.path(scatter_out, paste0("scatter_", safe_name(nm), ".png"))
            ggplot2::ggsave(fn, plt, width = 12, height = 8, dpi = 150)
        }
    }
}

cat("\nCompleted. Outputs:\n")
cat(paste0("  Tables: ", pathways_path, "\n"))
if (requireNamespace("jsonlite", quietly = TRUE)) {
    cat(paste0("  Params: ", params_json, "\n"))
} else {
    cat(paste0("  Params (RDS): ", params_rds, "\n"))
}
cat(paste0("  Barplots: ", bar_out, "\n"))
cat(paste0("  Heatmap: ", heat_out, "\n"))
cat(paste0("  Scatter: ", scatter_out, "\n\n"))
