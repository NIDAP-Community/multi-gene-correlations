# Multi-Gene Correlations to Signature [IODC] [Beta] - Version 14
# Extracted from NIDAP template: ri.vector.main.template.363cef88-16c4-433b-a709-58aa605f0958
# Description: Correlate a gene panel against a signature score across user-selected categories and
# visualise both individual bar plots and an aggregate heatmap.

utils::globalVariables(c(
    "annotation",
    "corr",
    "gene",
    "p_value",
    "percorder"
))

#' Multi-Gene Correlations to Signature [IODC] [Beta] (v14)
#'
#' Reimplementation of the NIDAP "Multi-Gene Correlations to Signature" template (version 14).
#' The helper computes Pearson correlations between a user-supplied gene panel and a signature score,
#' generates per-category correlation bar plots, and renders a combined heatmap. Results are returned
#' alongside the generated plots and run metadata so downstream notebooks can reproduce the template output.
#'
#' @param normalized_counts Data frame of normalized counts (genes in rows, samples in columns).
#' @param sample_metadata Data frame of sample metadata containing category assignments.
#' @param gene_column Character string naming the gene identifier column in `normalized_counts`.
#' @param sample_column Character string naming the sample identifier column in `sample_metadata`.
#' @param samples_to_include Character vector of sample identifiers to analyse. Defaults to the
#'   intersection of samples shared by counts and metadata when `NULL`.
#' @param genes Character vector of genes to correlate against the signature.
#' @param category_column Character string naming the metadata column that defines categories.
#' @param categories Character vector of category values to include.
#' @param signature_name Character string used in plot titles for the signature.
#' @param signature_genes Character vector of genes defining the signature score (averaged per sample).
#' @param margin_adjust Numeric margin adjustment (cm) for shifting row labels in the heatmap.
#' @param show_row_dendrogram Logical controlling whether to draw the row dendrogram.
#' @param heatmap_width Numeric width (cm) applied to each heatmap column.
#' @param add_text_to_heatmap Logical; when `TRUE`, prints correlation values inside the heatmap cells.
#' @param rename_categories Optional character vector supplying display labels for each category in order.
#' @param sum_duplicates Logical; when `TRUE`, duplicate gene rows are summed, otherwise the maximum row is kept.
#'
#' @return List with components `pathways` (merged correlation table), `plots` (per-category bar plots and the
#'   heatmap object), and `run_parameters` (metadata describing the run setup).
#' @export
#'
Multi_Gene_Correlations_to_Signature <- function(
  normalized_counts,
  sample_metadata,
  gene_column,
  sample_column,
  samples_to_include = NULL,
  genes = character(),
  category_column,
  categories = character(),
  signature_name = "Signature",
  signature_genes = character(),
  margin_adjust = 4,
  show_row_dendrogram = FALSE,
  heatmap_width = 6,
  add_text_to_heatmap = TRUE,
  rename_categories = character(),
  sum_duplicates = FALSE
) {
    suppressPackageStartupMessages({
        library(dplyr)
        library(tidyr)
        library(stringr)
        library(ggplot2)
        library(ComplexHeatmap)
        library(circlize)
        library(RColorBrewer)
        library(grid)
    })

    if (!requireNamespace("l2psupp", quietly = TRUE)) {
        stop("Package 'l2psupp' is required for gene symbol updates.")
    }

    if (!is.data.frame(normalized_counts)) {
        stop("normalized_counts must be a data frame.")
    }
    if (!is.data.frame(sample_metadata)) {
        stop("sample_metadata must be a data frame.")
    }
    if (!gene_column %in% colnames(normalized_counts)) {
        stop("gene_column not found in normalized_counts.")
    }
    if (!sample_column %in% colnames(sample_metadata)) {
        stop("sample_column not found in sample_metadata.")
    }
    if (!category_column %in% colnames(sample_metadata)) {
        stop("category_column not found in sample_metadata.")
    }

    if (length(genes) == 0) {
        stop("genes vector must not be empty.")
    }
    if (length(signature_genes) == 0) {
        stop("signature_genes vector must not be empty.")
    }
    if (length(categories) == 0) {
        stop("categories vector must not be empty.")
    }

    categories <- unique(categories)
    if (length(rename_categories) > 0 && length(rename_categories) != length(categories)) {
        stop("rename_categories must be empty or match the length of categories.")
    }

    meta_filtered <- sample_metadata %>%
        filter(.data[[category_column]] %in% categories)
    if (nrow(meta_filtered) == 0) {
        stop("No samples remain after filtering to the requested categories.")
    }

    available_from_counts <- setdiff(colnames(normalized_counts), gene_column)
    if (is.null(samples_to_include) || length(samples_to_include) == 0) {
        samples_to_include <- intersect(meta_filtered[[sample_column]], available_from_counts)
    } else {
        samples_to_include <- intersect(unique(samples_to_include), available_from_counts)
        samples_to_include <- intersect(samples_to_include, meta_filtered[[sample_column]])
    }
    if (length(samples_to_include) == 0) {
        stop("No overlapping samples between counts and metadata after applying selections.")
    }

    missing_samples <- setdiff(meta_filtered[[sample_column]], available_from_counts)
    if (length(missing_samples) > 0) {
        warning(sprintf("%d metadata samples are missing in the counts matrix and will be ignored.", length(missing_samples)))
    }

    meta_filtered <- meta_filtered %>%
        filter(.data[[sample_column]] %in% samples_to_include)

    normalized_counts <- normalized_counts %>%
        select(dplyr::all_of(c(gene_column, samples_to_include)))

    numeric_columns <- setdiff(colnames(normalized_counts), gene_column)
    numeric_columns <- numeric_columns[sapply(normalized_counts[numeric_columns], is.numeric)]
    if (length(numeric_columns) == 0) {
        stop("No numeric sample columns detected in normalized_counts after filtering.")
    }

    update_gene_symbols <- function(gene_vec) {
        if (length(gene_vec) == 0) {
            mapping <- data.frame(old = character(), newname = character(), stringsAsFactors = FALSE)
            return(list(updated = character(), mapping = mapping))
        }
        upd <- l2psupp::updategenes(gene_vec, trust = 1)
        if (is.data.frame(upd) && "newname" %in% colnames(upd)) {
            mapping <- upd
            new_names <- upd$newname
        } else if (is.list(upd) && !is.null(upd$newname)) {
            mapping <- as.data.frame(upd, stringsAsFactors = FALSE)
            new_names <- mapping$newname
        } else {
            new_names <- as.character(upd)
            mapping <- data.frame(
                old = gene_vec,
                newname = new_names,
                stringsAsFactors = FALSE
            )
        }
        if (!"old" %in% colnames(mapping)) {
            mapping$old <- gene_vec
        }
        if (!"newname" %in% colnames(mapping)) {
            mapping$newname <- unname(new_names)
        }
        list(updated = unname(new_names), mapping = mapping)
    }

    extract_mapping_column <- function(mapping, column) {
        if (is.null(mapping) || !column %in% colnames(mapping)) {
            return(character())
        }
        mapping[[column]]
    }

    genes_info <- update_gene_symbols(genes)
    genes_requested_original <- extract_mapping_column(genes_info$mapping, "old")
    genes_requested_new <- extract_mapping_column(genes_info$mapping, "newname")
    genes_requested_original <- genes_requested_original[!is.na(genes_requested_original)]
    genes_requested_new <- genes_requested_new[!is.na(genes_requested_new)]
    genes <- genes_info$updated

    signature_info <- update_gene_symbols(signature_genes)
    signature_requested_original <- extract_mapping_column(signature_info$mapping, "old")
    signature_requested_new <- extract_mapping_column(signature_info$mapping, "newname")
    signature_requested_original <- signature_requested_original[!is.na(signature_requested_original)]
    signature_requested_new <- signature_requested_new[!is.na(signature_requested_new)]
    signature_genes <- signature_info$updated

    if (nrow(genes_info$mapping) > 0) {
        message("Gene name updates (requested → used):")
        apply(genes_info$mapping, 1, function(row) message(sprintf("  %s → %s", row[["old"]], row[["newname"]])))
    }
    if (nrow(signature_info$mapping) > 0) {
        message("Signature gene updates (requested → used):")
        apply(signature_info$mapping, 1, function(row) message(sprintf("  %s → %s", row[["old"]], row[["newname"]])))
    }

    normcounts_filt <- normalized_counts %>%
        filter(.data[[gene_column]] != "---") %>%
        select(dplyr::all_of(c(gene_column, numeric_columns)))

    if (sum_duplicates) {
        normcounts_filt <- normcounts_filt %>%
            group_by(.data[[gene_column]]) %>%
            summarise(across(everything(), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
    } else {
        normcounts_filt <- normcounts_filt %>%
            mutate(.max_row = rowSums(across(-dplyr::all_of(gene_column), ~ replace(.x, is.na(.x), 0)))) %>%
            group_by(.data[[gene_column]]) %>%
            filter(.max_row == max(.max_row, na.rm = TRUE)) %>%
            select(-.max_row) %>%
            ungroup()
    }

    normcounts_matrix <- normcounts_filt %>%
        select(-dplyr::all_of(gene_column)) %>%
        as.matrix()
    rownames(normcounts_matrix) <- normcounts_filt[[gene_column]]

    genes_present <- intersect(genes, rownames(normcounts_matrix))
    if (length(genes_present) == 0) {
        stop("None of the requested genes are present in the normalized counts matrix.")
    }

    signature_present <- intersect(signature_genes, rownames(normcounts_matrix))
    if (length(signature_present) == 0) {
        stop("None of the signature genes are present in the normalized counts matrix.")
    }

    gene_expression <- normcounts_matrix[genes_present, samples_to_include, drop = FALSE]
    signature_counts <- normcounts_matrix[signature_present, samples_to_include, drop = FALSE]
    signature_vector <- colMeans(signature_counts, na.rm = TRUE)

    signature_used_display <- if (length(signature_present) > 0) signature_present else signature_genes
    signature_label <- paste(signature_name, "Signature")
    signature_gene_label <- paste(signature_used_display, collapse = ", ")
    full_signature_title <- sprintf("%s (%s)", signature_label, signature_gene_label)

    rd_bu_colors <- RColorBrewer::brewer.pal(11, "RdBu")

    calc_correlation <- function(df, vec) {
        apply(df, 1, function(row) {
            ok <- complete.cases(row, vec)
            if (sum(ok) >= 3) {
                test_result <- suppressWarnings(stats::cor.test(row[ok], vec[ok], method = "pearson"))
                ci <- test_result$conf.int
                if (is.null(ci) || length(ci) < 2) {
                    ci_vals <- c(NA_real_, NA_real_)
                } else {
                    ci_vals <- ci[1:2]
                }
                list(
                    corr = round(unname(test_result$estimate), 3),
                    pval = test_result$p.value,
                    conf_low = ci_vals[1],
                    conf_high = ci_vals[2]
                )
            } else {
                list(corr = NA_real_, pval = NA_real_, conf_low = NA_real_, conf_high = NA_real_)
            }
        })
    }

    calculate_correlation <- function(filtcounts, signature_values, category_label) {
        results <- calc_correlation(filtcounts, signature_values)
        correlations <- vapply(results, function(x) x$corr, numeric(1))
        p_values <- vapply(results, function(x) x$pval, numeric(1))
        conf_low <- vapply(results, function(x) x$conf_low, numeric(1))
        conf_high <- vapply(results, function(x) x$conf_high, numeric(1))

        corr_table <- data.frame(
            gene = names(correlations),
            corr = correlations,
            p_value = p_values,
            conf_low = conf_low,
            conf_high = conf_high,
            stringsAsFactors = FALSE
        )
        corr_table$p_adj <- stats::p.adjust(corr_table$p_value, method = "BH")
        corr_table$abs_corr <- abs(corr_table$corr)
        corr_table <- corr_table %>%
            mutate(
                annotation = dplyr::case_when(
                    is.na(.data$p_adj) ~ "",
                    .data$p_adj < 0.001 ~ "***",
                    .data$p_adj < 0.01 ~ "**",
                    .data$p_adj < 0.05 ~ "*",
                    TRUE ~ ""
                ),
                gene = reorder(.data$gene, .data$corr, na.rm = TRUE),
                text_hjust = ifelse(.data$corr >= 0, -0.4, 1.4)
            )

        signature_df <- data.frame(
            sample = names(signature_values),
            signature = as.numeric(signature_values),
            stringsAsFactors = FALSE
        )
        expression_wide <- as.data.frame(t(filtcounts), stringsAsFactors = FALSE)
        expression_wide$sample <- rownames(expression_wide)
        scatter_df <- expression_wide %>%
            tidyr::pivot_longer(
                cols = -sample,
                names_to = "gene",
                values_to = "expression"
            ) %>%
            left_join(signature_df, by = "sample") %>%
            mutate(
                category = category_label
            ) %>%
            filter(!is.na(.data$expression), !is.na(.data$signature))

        title_text <- sprintf(
            "Correlation to %s in %s (n=%d)",
            full_signature_title,
            category_label,
            length(signature_values)
        )
        wrapped_title <- stringr::str_wrap(title_text, width = ifelse(nchar(title_text) > 80, 60, 80))
        title_size <- ifelse(nchar(title_text) > 80, 12, 15)

        max_abs_plot <- max(abs(corr_table$corr), na.rm = TRUE)
        if (!is.finite(max_abs_plot) || max_abs_plot == 0) {
            max_abs_plot <- 1
        }
        rounding_step_plot <- max(0.1, signif(max_abs_plot / 5, 1))
        max_abs_plot <- ceiling(max_abs_plot / rounding_step_plot) * rounding_step_plot
        mid_plot <- max_abs_plot / 2
        fill_breaks <- c(-max_abs_plot, -mid_plot, 0, mid_plot, max_abs_plot)
        fill_labels <- c(
            sprintf("-%0.2f", max_abs_plot),
            sprintf("-%0.2f", mid_plot),
            "0.00",
            sprintf("%0.2f", mid_plot),
            sprintf("%0.2f", max_abs_plot)
        )

        gp <- ggplot(corr_table, aes(x = gene, y = corr, fill = corr)) +
            geom_bar(stat = "identity") +
            theme_bw() +
            coord_flip() +
            geom_text(
                aes(label = annotation, hjust = text_hjust),
                color = "black"
            ) +
            scale_fill_gradientn(
                colors = rev(rd_bu_colors),
                na.value = "grey80",
                guide = "colourbar",
                limits = c(-max_abs_plot, max_abs_plot),
                breaks = fill_breaks,
                labels = fill_labels
            ) +
            scale_y_continuous(limits = c(-1, 1)) +
            theme(
                axis.text.y = element_text(hjust = 1, size = 12),
                axis.text.x = element_text(size = 12),
                axis.title.y = element_text(size = 14, margin = margin(t = 10, b = 10)),
                axis.title.x = element_text(size = 14, margin = margin(t = 10, b = 10)),
                plot.title = element_text(hjust = 0.5, size = title_size)
            ) +
            labs(
                title = wrapped_title,
                x = "Gene",
                y = "Correlation Coefficient"
            )
        print(gp)

        corr_table <- corr_table %>% arrange(.data$corr)
        tidy_table <- corr_table %>%
            as.data.frame() %>%
            select(gene, corr, p_value, p_adj, conf_low, conf_high)
        colnames(tidy_table)[colnames(tidy_table) == "corr"] <- paste0(category_label, "_corr")
        colnames(tidy_table)[colnames(tidy_table) == "p_value"] <- paste0(category_label, "_p_value")
        colnames(tidy_table)[colnames(tidy_table) == "p_adj"] <- paste0(category_label, "_p_adj")
        colnames(tidy_table)[colnames(tidy_table) == "conf_low"] <- paste0(category_label, "_conf_low")
        colnames(tidy_table)[colnames(tidy_table) == "conf_high"] <- paste0(category_label, "_conf_high")
        colnames(tidy_table) <- gsub(" ", "", colnames(tidy_table), fixed = TRUE)

        list(table = tidy_table, plot = gp, scatter = scatter_df)
    }

    category_labels <- if (length(rename_categories) > 0) rename_categories else categories
    names(category_labels) <- categories

    gene_tables <- list()
    bar_plots <- list()
    samples_by_category <- list()
    scatter_records <- list()

    for (cat in categories) {
        label <- category_labels[[cat]]
        cat_meta <- meta_filtered %>% filter(.data[[category_column]] == cat)
        current_samples <- intersect(cat_meta[[sample_column]], colnames(gene_expression))
        samples_by_category[[label]] <- length(current_samples)
        if (length(current_samples) == 0) {
            warning(sprintf("Category '%s' has no matching samples after filtering and will be skipped.", cat))
            next
        }
        expr_subset <- gene_expression[, current_samples, drop = FALSE]
        sig_subset <- signature_vector[current_samples]
        cat_result <- calculate_correlation(expr_subset, sig_subset, label)
        gene_tables[[cat]] <- cat_result$table
        bar_plots[[label]] <- cat_result$plot
        scatter_records[[label]] <- cat_result$scatter
    }

    gene_tables <- gene_tables[!vapply(gene_tables, is.null, logical(1))]
    if (length(gene_tables) == 0) {
        stop("No correlation results could be generated for the requested categories.")
    }

    all_result <- calculate_correlation(gene_expression, signature_vector, "all")
    gene_tables[["all"]] <- all_result$table
    bar_plots[["all"]] <- all_result$plot
    scatter_records[["all"]] <- all_result$scatter

    correlation_vectors <- lapply(names(gene_tables), function(name) {
        tbl <- gene_tables[[name]]
        corr_col <- grep("_corr$", colnames(tbl), value = TRUE)[1]
        stats::setNames(tbl[[corr_col]], tbl$gene)
    })
    names(correlation_vectors) <- names(gene_tables)
    all_genes <- unique(unlist(lapply(correlation_vectors, names)))
    correlation_matrix <- do.call(cbind, lapply(correlation_vectors, function(vec) vec[all_genes]))
    rownames(correlation_matrix) <- all_genes

    heatmap_object <- NULL
    if (!is.null(correlation_matrix) && ncol(correlation_matrix) > 0) {
        display_names <- names(correlation_vectors)
        display_names <- ifelse(display_names %in% names(category_labels), category_labels[display_names], display_names)
        colnames(correlation_matrix) <- display_names

        heatmap_title_text <- sprintf("Gene correlations to %s", full_signature_title)
        heatmap_title <- stringr::str_wrap(heatmap_title_text, width = 40)

        max_abs_corr <- max(abs(correlation_matrix), na.rm = TRUE)
        if (!is.finite(max_abs_corr) || max_abs_corr == 0) {
            max_abs_corr <- 1
        }
        rounding_step <- max(0.1, signif(max_abs_corr / 5, 1))
        max_abs_corr <- ceiling(max_abs_corr / rounding_step) * rounding_step
        mid_corr <- max_abs_corr / 2
        color_breaks <- c(-max_abs_corr, -mid_corr, 0, mid_corr, max_abs_corr)
        legend_breaks <- color_breaks
        max_label <- formatC(max_abs_corr, format = "f", digits = 2)
        mid_label <- formatC(mid_corr, format = "f", digits = 2)
        legend_labels <- c(
            paste0("-", max_label),
            paste0("-", mid_label),
            "0.00",
            mid_label,
            max_label
        )
        col_fun <- circlize::colorRamp2(color_breaks, c("darkblue", "lightblue", "white", "pink", "darkred"))
        heatmap_object <- ComplexHeatmap::Heatmap(
            correlation_matrix,
            name = "Correlation",
            column_title = heatmap_title,
            col = col_fun,
            show_row_dend = show_row_dendrogram,
            width = grid::unit(rep(heatmap_width, ncol(correlation_matrix)), "cm"),
            row_names_side = "left",
            row_names_max_width = grid::unit(6, "cm"),
            heatmap_legend_param = list(
                title = "Correlation",
                at = legend_breaks,
                labels = legend_labels
            )
        )
        if (isTRUE(add_text_to_heatmap)) {
            heatmap_object@matrix_param$cell_fun <- function(j, i, x, y, width, height, fill) {
                grid::grid.text(sprintf("%.1f", correlation_matrix[i, j]), x, y, gp = grid::gpar(fontsize = 10))
            }
        }
        # print(heatmap_object)
    }

    merged_df <- Reduce(function(acc, tbl) dplyr::full_join(acc, tbl, by = "gene"), gene_tables)

    scatter_plots <- lapply(names(scatter_records), function(name) {
        scatter_df <- scatter_records[[name]]
        if (is.null(scatter_df) || nrow(scatter_df) == 0) {
            return(NULL)
        }
        ggplot(scatter_df, aes(x = signature, y = expression)) +
            geom_point(alpha = 0.6, size = 1.5, color = "#336699") +
            geom_smooth(method = "lm", se = FALSE, color = "#FF8C00") +
            facet_wrap(~gene, scales = "free") +
            theme_bw() +
            labs(
                title = sprintf("%s vs Gene Expression in %s", full_signature_title, name),
                x = sprintf("%s Score", full_signature_title),
                y = "Gene Expression"
            )
    })
    names(scatter_plots) <- names(scatter_records)

    run_parameters <- list(
        template_id = "ri.vector.main.template.363cef88-16c4-433b-a709-58aa605f0958",
        version = 14,
        signature_name = signature_name,
        signature_title = full_signature_title,
        genes_requested = length(genes_requested_original),
        genes_used = length(genes_present),
        genes_missing = setdiff(genes_requested_new, genes_present),
        signature_genes_requested = length(signature_requested_original),
        signature_genes_used = length(signature_present),
        signature_genes_missing = setdiff(signature_requested_new, signature_present),
        signature_genes_display = signature_used_display,
        categories_requested = categories,
        samples_total = length(samples_to_include),
        samples_by_category = samples_by_category,
        sum_duplicates = sum_duplicates,
        heatmap_width = heatmap_width,
        margin_adjust = margin_adjust,
        show_row_dendrogram = show_row_dendrogram,
        add_text_to_heatmap = add_text_to_heatmap,
        fdr_correction = "BH"
    )

    plots <- list(
        barplots = bar_plots,
        heatmap = heatmap_object,
        scatter = scatter_plots
    )

    list(pathways = merged_df, plots = plots, run_parameters = run_parameters)
}
