#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, digits = 17)

parse_args <- function(args) {
    expected <- c("--samplesheet", "--output-dir")
    if (length(args) != 4L || !all(expected %in% args[c(1L, 3L)])) {
        stop("usage: run_deseq2.R --samplesheet FILE --output-dir DIR")
    }
    values <- setNames(args[c(2L, 4L)], args[c(1L, 3L)])
    list(samplesheet = values[["--samplesheet"]], output_dir = values[["--output-dir"]])
}

ensure_packages <- function() {
    required <- c(DESeq2 = "1.46.0", airway = "1.26.0")
    installed <- vapply(
        names(required),
        function(package) requireNamespace(package, quietly = TRUE),
        logical(1)
    )
    versions_match <- installed
    versions_match[installed] <- vapply(
        names(required)[installed],
        function(package) as.character(packageVersion(package)) == required[[package]],
        logical(1)
    )
    if (!all(versions_match)) {
        stop(
            "Pinned packages are unavailable. Build containers/Dockerfile or install ",
            "DESeq2 1.46.0 and airway 1.26.0 before running."
        )
    }
}

write_png <- function(path, expression, width = 1400L, height = 900L) {
    grDevices::png(path, width = width, height = height, res = 150)
    on.exit(grDevices::dev.off(), add = TRUE)
    force(expression)
}

json_escape <- function(value) {
    value <- gsub("\\\\", "\\\\\\\\", value)
    gsub('"', '\\\\"', value)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
ensure_packages()
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(airway))

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
samples <- read.csv(args$samplesheet, check.names = FALSE)
stopifnot(nrow(samples) == 8L, length(unique(samples$cell)) == 4L)

data("airway", package = "airway")
stopifnot(setequal(samples$sample_id, colnames(airway)))
airway <- airway[, samples$sample_id]

airway_cells <- as.character(SummarizedExperiment::colData(airway)$cell)
airway_dex <- as.character(SummarizedExperiment::colData(airway)$dex)
stopifnot(identical(samples$cell, airway_cells), identical(samples$dex, airway_dex))

SummarizedExperiment::colData(airway)$cell <- factor(airway_cells)
SummarizedExperiment::colData(airway)$dex <- factor(airway_dex, levels = c("untrt", "trt"))
stopifnot(levels(SummarizedExperiment::colData(airway)$dex)[1L] == "untrt")

dds <- DESeqDataSet(airway, design = ~ cell + dex)
genes_before <- nrow(dds)
keep <- rowSums(counts(dds) >= 10L) >= 4L
dds <- dds[keep, ]
genes_after <- nrow(dds)
stopifnot(genes_after > 0L, genes_after < genes_before)

raw_counts <- counts(dds)
library_sizes <- colSums(raw_counts)
write_png(file.path(args$output_dir, "library_sizes.png"), {
    barplot(
        library_sizes,
        las = 2,
        col = ifelse(samples$dex == "trt", "#C44E52", "#4C72B0"),
        ylab = "Library size (processed gene counts)",
        main = "GSE52778 processed-count library sizes"
    )
})

write_png(file.path(args$output_dir, "count_distributions.png"), {
    boxplot(
        log10(raw_counts + 1),
        names = samples$sample_id,
        las = 2,
        col = ifelse(samples$dex == "trt", "#C44E52", "#4C72B0"),
        ylab = "log10(count + 1)",
        main = "Filtered processed-count distributions"
    )
})

dds <- DESeq(dds, quiet = TRUE)
res <- results(
    dds,
    contrast = c("dex", "trt", "untrt"),
    alpha = 0.05,
    independentFiltering = TRUE
)

feature_data <- as.data.frame(SummarizedExperiment::rowData(airway))
symbol_column <- intersect(c("symbol", "gene_name"), colnames(feature_data))[1L]
if (is.na(symbol_column)) {
    stop("airway row metadata does not contain symbol or gene_name")
}
symbols <- as.character(feature_data[rownames(res), symbol_column])

full <- data.frame(
    gene_id = rownames(res),
    gene_symbol = symbols,
    baseMean = res$baseMean,
    log2FoldChange = res$log2FoldChange,
    lfcSE = res$lfcSE,
    stat = res$stat,
    pvalue = res$pvalue,
    padj = res$padj,
    check.names = FALSE
)
full <- full[order(full$gene_id), , drop = FALSE]
significant <- full[!is.na(full$padj) & full$padj <= 0.05, , drop = FALSE]
significant <- significant[order(significant$padj, significant$gene_id), , drop = FALSE]

crispld2 <- full[full$gene_symbol == "CRISPLD2", , drop = FALSE]
if (nrow(crispld2) < 1L || !any(crispld2$log2FoldChange > 0, na.rm = TRUE)) {
    stop("CRISPLD2 sanity check failed for trt versus untrt")
}

write.csv(
    full,
    file.path(args$output_dir, "full_results.csv"),
    row.names = FALSE,
    na = "",
    quote = FALSE
)
write.csv(
    significant,
    file.path(args$output_dir, "fdr_0.05_results.csv"),
    row.names = FALSE,
    na = "",
    quote = FALSE
)

transformed <- vst(dds, blind = FALSE)
pca_data <- plotPCA(transformed, intgroup = c("dex", "cell"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))
write_png(file.path(args$output_dir, "pca.png"), {
    plot(
        pca_data$PC1,
        pca_data$PC2,
        pch = 19,
        col = ifelse(pca_data$dex == "trt", "#C44E52", "#4C72B0"),
        xlab = sprintf("PC1: %s%% variance", percent_var[1L]),
        ylab = sprintf("PC2: %s%% variance", percent_var[2L]),
        main = "Variance-stabilised PCA"
    )
    text(pca_data$PC1, pca_data$PC2, labels = rownames(pca_data), pos = 3, cex = 0.65)
    legend("topright", legend = c("untreated", "treated"), col = c("#4C72B0", "#C44E52"), pch = 19)
})

sample_distances <- as.matrix(dist(t(assay(transformed))))
write_png(file.path(args$output_dir, "sample_distance_heatmap.png"), {
    heatmap(
        sample_distances,
        symm = TRUE,
        margins = c(9, 9),
        col = hcl.colors(50, "Blues 3", rev = TRUE),
        main = "Sample distances after variance stabilisation"
    )
})

metadata <- c(
    "{",
    '  "dataset_accession": "GSE52778",',
    '  "input_type": "processed_counts",',
    '  "design": "~ cell + dex",',
    '  "reference_level": "untrt",',
    '  "contrast": {',
    '    "factor": "dex",',
    '    "numerator": "trt",',
    '    "denominator": "untrt"',
    '  },',
    '  "low_count_filter": {',
    '    "rule": "rowSums(counts >= 10) >= 4",',
    '    "applied": true,',
    sprintf('    "genes_before": %d,', genes_before),
    sprintf('    "genes_after": %d', genes_after),
    '  },',
    '  "alpha": 0.05,',
    sprintf('  "crispld2_log2_fold_change": %.17g', max(crispld2$log2FoldChange, na.rm = TRUE)),
    "}"
)
writeLines(metadata, file.path(args$output_dir, "analysis_metadata.json"), useBytes = TRUE)

session_lines <- capture.output(sessionInfo())
writeLines(session_lines, file.path(args$output_dir, "sessionInfo.txt"), useBytes = TRUE)
