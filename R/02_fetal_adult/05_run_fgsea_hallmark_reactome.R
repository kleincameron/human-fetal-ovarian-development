# fgsea analysis of fetal-vs-adult pseudobulk edgeR results.
#
# This script reads standardized pseudobulk edgeR DEG tables, builds ranked
# gene statistics per cell type, and runs fgsea against local MSigDB Hallmark
# and Reactome GMT files.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(fgsea)
  library(ggplot2)
  library(tibble)
})

set.seed(42)
options(bitmapType = "cairo")

project_root <- "/home/liyan/liyan/Final/github_code_for_publication"
results_base <- "/home/liyan/liyan/Final/github_code_for_publication_results"

source(file.path(project_root, "config", "plotting.R"))
source(file.path(project_root, "config", "labels_colors.R"))

paths_local <- file.path(project_root, "config", "paths_local.R")
if (file.exists(paths_local)) {
  source(paths_local)
}

deg_results_root <- file.path(results_base, "fetal_adult_pseudobulk_edgeR")
deg_table_dir <- file.path(deg_results_root, "tables")

qc_summary_csv <- file.path(
  deg_table_dir,
  "PSEUDOBULK_edgeR_QC_summary.csv"
)

default_gmt_hallmark <- "/home/liyan/liyan/references/msigdb/h.all.v2023.2.Hs.symbols.gmt"
default_gmt_reactome <- "/home/liyan/liyan/references/msigdb/c2.cp.reactome.v2023.2.Hs.symbols.gmt"

gmt_hallmark <- if (exists("msigdb_hallmark_gmt", inherits = FALSE)) {
  msigdb_hallmark_gmt
} else {
  default_gmt_hallmark
}

gmt_reactome <- if (exists("msigdb_reactome_gmt", inherits = FALSE)) {
  msigdb_reactome_gmt
} else {
  default_gmt_reactome
}

results_root <- file.path(results_base, "fetal_adult_fgsea_hallmark_reactome")

out_table_dir <- file.path(results_root, "tables")
out_plot_dir <- file.path(results_root, "figures")
out_log_dir <- file.path(results_root, "logs")

dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_log_dir, recursive = TRUE, showWarnings = FALSE)

minSize <- 15
maxSize <- 500
nPermSimple <- 10000
top_n_plot <- 20
padj_plot_cut <- 0.25
skip_celltypes <- c("germ")
plot_celltype_colors <- c(
  celltype_colors,
  degenerated = "#CFCFCF",
  epithelial = "#BFE6E2",
  theca = "#FB8072",
  "glial cell" = "#F04E37",
  "plasma cell" = "#00E41A",
  "smooth muscle cell" = "#1F2BFF"
)

sanitize_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

pretty_pathway <- function(pathway) {
  pathway |>
    str_replace_all("^HALLMARK_", "") |>
    str_replace_all("^REACTOME_", "") |>
    str_replace_all("_", " ") |>
    str_to_lower() |>
    str_to_sentence()
}

save_publication_plot_pdf_png <- function(plot, pdf_file, png_file, width, height, dpi = 600) {
  ggsave(
    pdf_file,
    plot,
    width = width,
    height = height,
    device = cairo_pdf,
    bg = "white",
    limitsize = FALSE
  )

  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(
      png_file,
      plot,
      width = width,
      height = height,
      dpi = dpi,
      device = ragg::agg_png,
      background = "white",
      limitsize = FALSE
    )
  } else {
    ggsave(
      png_file,
      plot,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white",
      limitsize = FALSE
    )
  }
}

make_rank <- function(df) {
  required_cols <- c("gene", "logFC")
  missing_cols <- setdiff(required_cols, colnames(df))

  if (length(missing_cols) > 0) {
    stop("DE table missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  df <- df |>
    mutate(
      gene = as.character(gene),
      logFC = as.numeric(logFC)
    ) |>
    filter(
      !is.na(gene),
      gene != "",
      !is.na(logFC)
    )
  if ("F" %in% colnames(df) && any(is.finite(df$F), na.rm = TRUE)) {
    statistic <- sign(df$logFC) * as.numeric(df$F)
    statistic_name <- "sign_logFC_times_edgeR_QLF"
  } else if ("PValue" %in% colnames(df) && any(is.finite(df$PValue), na.rm = TRUE)) {
    statistic <- sign(df$logFC) * (-log10(pmax(as.numeric(df$PValue), 1e-300)))
    statistic_name <- "sign_logFC_times_neglog10_PValue"
  } else {
    statistic <- df$logFC
    statistic_name <- "logFC"
  }

  ranked_table <- df |>
    mutate(rank_statistic = statistic) |>
    filter(is.finite(rank_statistic)) |>
    group_by(gene) |>
    slice_max(order_by = abs(rank_statistic), n = 1, with_ties = FALSE) |>
    ungroup() |>
    arrange(desc(rank_statistic))

  ranks <- ranked_table$rank_statistic
  names(ranks) <- ranked_table$gene
  ranks <- sort(ranks, decreasing = TRUE)

  list(
    ranks = ranks,
    ranked_table = ranked_table,
    statistic_name = statistic_name
  )
}

run_fgsea_one <- function(ranks, pathways, collection_name) {
  res <- fgseaMultilevel(
    pathways = pathways,
    stats = ranks,
    minSize = minSize,
    maxSize = maxSize,
    eps = 0.0,
    nPermSimple = 10000
  ) |>
    as_tibble() |>
    arrange(padj, desc(abs(NES))) |>
    mutate(
      collection = collection_name,
      direction = if_else(NES > 0, "Higher in fetal", "Higher in adult"),
      leadingEdge = vapply(leadingEdge, paste, collapse = ";", character(1))
    ) |>
    relocate(collection, direction, pathway, NES, padj, pval, size, leadingEdge)

  res
}

select_pathways_for_plot <- function(res) {
  res_ranked <- res |>
    filter(!is.na(padj), !is.na(NES)) |>
    arrange(padj, desc(abs(NES)))

  res_preferred <- res_ranked |>
    filter(padj <= padj_plot_cut)

  if (nrow(res_preferred) == 0) {
    res_preferred <- res_ranked
  }

  bind_rows(
    res_preferred |>
      filter(NES > 0) |>
      slice_head(n = ceiling(top_n_plot / 2)),
    res_preferred |>
      filter(NES < 0) |>
      slice_head(n = floor(top_n_plot / 2))
  ) |>
    distinct(pathway, .keep_all = TRUE) |>
    mutate(
      pathway_pretty = pretty_pathway(pathway),
      is_significant = padj < 0.05
    ) |>
    arrange(NES) |>
    mutate(pathway_pretty = factor(pathway_pretty, levels = pathway_pretty))
}
make_fgsea_barplot <- function(plot_df, celltype) {
  base_color <- plot_celltype_colors[[celltype]]
  if (is.null(base_color) || is.na(base_color)) {
    base_color <- "#4D4D4D"
  }

  ggplot(plot_df, aes(x = NES, y = pathway_pretty)) +
    geom_vline(
      xintercept = 0,
      color = "grey55",
      linewidth = 0.25
    ) +
    geom_col(
      fill = base_color,
      width = 0.72
    ) +
    geom_point(
      data = plot_df |> filter(is_significant),
      aes(x = NES, y = pathway_pretty),
      inherit.aes = FALSE,
      shape = 21,
      size = 1.2,
      stroke = 0.25,
      fill = "white",
      color = "black"
    ) +
    labs(
      x = "Normalized enrichment score",
      y = NULL
    ) +
    theme_publication()
}

stopifnot(dir.exists(deg_table_dir))
stopifnot(file.exists(qc_summary_csv))
stopifnot(file.exists(gmt_hallmark))
stopifnot(file.exists(gmt_reactome))

message("Loading GMT files.")
pathways_hallmark <- fgsea::gmtPathways(gmt_hallmark)
pathways_reactome <- fgsea::gmtPathways(gmt_reactome)

collections <- list(
  HALLMARK = pathways_hallmark,
  REACTOME = pathways_reactome
)

qc_summary <- read_csv(qc_summary_csv, show_col_types = FALSE)

required_qc_cols <- c("celltype", "status", "output_file")
missing_qc_cols <- setdiff(required_qc_cols, colnames(qc_summary))

if (length(missing_qc_cols) > 0) {
  stop("Pseudobulk QC summary missing required columns: ", paste(missing_qc_cols, collapse = ", "))
}

analysis_qc <- qc_summary |>
  mutate(
    celltype = as.character(celltype),
    status = as.character(status),
    output_file = as.character(output_file)
  ) |>
  filter(
    status == "ok_saved",
    !tolower(celltype) %in% tolower(skip_celltypes)
  )

if (nrow(analysis_qc) == 0) {
  stop("No successful pseudobulk DE results available for fgsea.")
}

missing_deg_files <- analysis_qc$output_file[!file.exists(analysis_qc$output_file)]
if (length(missing_deg_files) > 0) {
  stop("Missing pseudobulk DEG files:\n", paste(missing_deg_files, collapse = "\n"))
}

write_csv(
  analysis_qc,
  file.path(out_table_dir, "fgsea_input_pseudobulk_tables.csv")
)
summary_rows <- list()
selected_plot_rows <- list()

for (i in seq_len(nrow(analysis_qc))) {
  celltype <- analysis_qc$celltype[[i]]
  deg_file <- analysis_qc$output_file[[i]]

  message("\n=============================")
  message("Cell type: ", celltype)

  deg_table <- read_csv(deg_file, show_col_types = FALSE)

  rank_object <- make_rank(deg_table)
  ranks <- rank_object$ranks

  message("Ranked genes: ", length(ranks))
  message("Rank statistic: ", rank_object$statistic_name)

  write_csv(
    rank_object$ranked_table,
    file.path(
      out_table_dir,
      paste0("fgsea_ranked_genes__", sanitize_filename(celltype), ".csv")
    )
  )

  for (collection_name in names(collections)) {
    message("Running fgsea: ", collection_name)

    res <- run_fgsea_one(
      ranks = ranks,
      pathways = collections[[collection_name]],
      collection_name = collection_name
    ) |>
      mutate(
        celltype = celltype,
        rank_statistic = rank_object$statistic_name
      ) |>
      relocate(celltype, collection, direction)

    out_csv <- file.path(
      out_table_dir,
      paste0("fgsea_", collection_name, "__", sanitize_filename(celltype), ".csv")
    )

    out_sig_csv <- file.path(
      out_table_dir,
      paste0("fgsea_", collection_name, "__", sanitize_filename(celltype), "__padj0.05.csv")
    )

    write_csv(res, out_csv)
    write_csv(
      res |> filter(padj < 0.05) |> arrange(padj, desc(abs(NES))),
      out_sig_csv
    )

    plot_df <- select_pathways_for_plot(res)

    if (nrow(plot_df) > 0) {
      file_prefix <- paste0(
        "fgsea_",
        collection_name,
        "_barplot_fetal_vs_adult__",
        sanitize_filename(celltype)
      )

      pdf_file <- file.path(out_plot_dir, paste0(file_prefix, ".pdf"))
      png_file <- file.path(out_plot_dir, paste0(file_prefix, ".png"))

      plot_height <- max(2.6, 0.20 * nrow(plot_df) + 0.8)

      save_publication_plot_pdf_png(
        plot = make_fgsea_barplot(plot_df, celltype),
        pdf_file = pdf_file,
        png_file = png_file,
        width = 3.8,
        height = plot_height,
        dpi = 600
      )
      selected_plot_rows[[length(selected_plot_rows) + 1]] <- plot_df |>
        mutate(
          celltype = celltype,
          collection = collection_name
        ) |>
        select(
          celltype,
          collection,
          pathway,
          pathway_pretty,
          NES,
          padj,
          direction,
          is_significant
        )
    } else {
      pdf_file <- NA_character_
      png_file <- NA_character_
    }

    summary_rows[[length(summary_rows) + 1]] <- tibble(
      celltype = celltype,
      collection = collection_name,
      rank_statistic = rank_object$statistic_name,
      n_ranked_genes = length(ranks),
      n_pathways_tested = nrow(res),
      n_padj_lt_0_05 = sum(res$padj < 0.05, na.rm = TRUE),
      n_padj_le_plot_cut = sum(res$padj <= padj_plot_cut, na.rm = TRUE),
      table_file = out_csv,
      significant_table_file = out_sig_csv,
      pdf_file = pdf_file,
      png_file = png_file
    )

    message(collection_name, ": saved tables and plot.")
  }
}

fgsea_summary <- bind_rows(summary_rows)

write_csv(
  fgsea_summary,
  file.path(out_table_dir, "fgsea_summary_counts.csv")
)

if (length(selected_plot_rows) > 0) {
  write_csv(
    bind_rows(selected_plot_rows),
    file.path(out_table_dir, "fgsea_selected_pathways_for_barplots.csv")
  )
}

write_csv(
  tibble(
    setting = c(
      "minSize",
      "maxSize",
      "nPermSimple",
      "top_n_plot",
      "padj_plot_cut",
      "skip_celltypes",
      "gmt_hallmark",
      "gmt_reactome"
    ),
    value = c(
      as.character(minSize),
      as.character(maxSize),
      as.character(nPermSimple),
      as.character(top_n_plot),
      as.character(padj_plot_cut),
      paste(skip_celltypes, collapse = "; "),
      gmt_hallmark,
      gmt_reactome
    )
  ),
  file.path(out_table_dir, "fgsea_settings.csv")
)

summary_lines <- c(
  paste("Input pseudobulk table directory:", deg_table_dir),
  paste("Input pseudobulk QC summary:", qc_summary_csv),
  paste("Hallmark GMT:", gmt_hallmark),
  paste("Reactome GMT:", gmt_reactome),
  paste("Hallmark pathways:", length(pathways_hallmark)),
  paste("Reactome pathways:", length(pathways_reactome)),
  paste("Output directory:", results_root),
  paste("Cell types analyzed:", paste(unique(fgsea_summary$celltype), collapse = ", ")),
  paste("Collections:", paste(names(collections), collapse = ", ")),
  paste("Minimum pathway size:", minSize),
  paste("Maximum pathway size:", maxSize),
  paste("fgsea nPermSimple:", nPermSimple),
  paste("Skipped cell types:", paste(skip_celltypes, collapse = ", "))
)

writeLines(
  summary_lines,
  file.path(out_log_dir, "fgsea_hallmark_reactome_summary.txt")
)

sink(file.path(out_log_dir, "fgsea_hallmark_reactome_sessionInfo.txt"))
print(sessionInfo())
sink()

cat("\nDONE.\n")
cat("fgsea outputs written to:\n  ", results_root, "\n", sep = "")
cat("Tables:\n  ", out_table_dir, "\n", sep = "")
cat("Plots:\n  ", out_plot_dir, "\n", sep = "")
