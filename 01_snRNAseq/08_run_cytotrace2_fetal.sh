#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-/home/liyan/liyan/Final/github_code_for_publication}"
RESULTS_BASE="${RESULTS_BASE:-/home/liyan/liyan/Final/github_code_for_publication_results}"

ANNOTATED_RDS="${ANNOTATED_RDS:-${RESULTS_BASE}/snRNAseq_annotated_object/objects/fetal_ovary_snRNAseq_canonical_annotated.rds}"

OUT_ROOT="${OUT_ROOT:-${RESULTS_BASE}/cytotrace2_fetal}"
IN_DIR="${OUT_ROOT}/inputs"
CT2_DIR="${OUT_ROOT}/cytotrace2"
PLOT_DIR="${OUT_ROOT}/figures"
TABLE_DIR="${OUT_ROOT}/tables"
LOG_DIR="${OUT_ROOT}/logs"

CYTOTRACE2_BIN="${CYTOTRACE2_BIN:-cytotrace2}"
CYTOTRACE2_SPECIES="${CYTOTRACE2_SPECIES:-human}"
CYTOTRACE2_BATCH_SIZE="${CYTOTRACE2_BATCH_SIZE:-10000}"
CYTOTRACE2_SUB_BATCH_SIZE="${CYTOTRACE2_SUB_BATCH_SIZE:-1000}"
CYTOTRACE2_THREADS="${CYTOTRACE2_THREADS:-32}"

mkdir -p "${IN_DIR}" "${CT2_DIR}" "${PLOT_DIR}" "${TABLE_DIR}" "${LOG_DIR}"

exec > >(tee "${LOG_DIR}/cytotrace2_fetal_pipeline.log") 2>&1

EXPR_TSV_GZ="${IN_DIR}/fetal_expression_raw_counts.tsv.gz"
EXPR_TSV_GZ_TMP="${EXPR_TSV_GZ}.tmp"
ANNO_TSV="${IN_DIR}/fetal_cytotrace2_annotation.tsv"
META_TSV="${IN_DIR}/fetal_metadata_for_plotting.tsv"

echo "[INFO] Project root: ${PROJECT_ROOT}"
echo "[INFO] Results base: ${RESULTS_BASE}"
echo "[INFO] Input annotated RDS: ${ANNOTATED_RDS}"
echo "[INFO] Output root: ${OUT_ROOT}"

if [[ ! -f "${ANNOTATED_RDS}" ]]; then
  echo "[ERROR] Annotated RDS not found: ${ANNOTATED_RDS}" >&2
  exit 1
fi

if [[ -s "${EXPR_TSV_GZ}" && -s "${ANNO_TSV}" && -s "${META_TSV}" ]] && gzip -t "${EXPR_TSV_GZ}" >/dev/null 2>&1; then
  echo "[INFO] Existing valid CytoTRACE2 input files found; skipping export."
else
  rm -f "${EXPR_TSV_GZ_TMP}"
  echo "[INFO] Exporting fetal metadata and raw counts for CytoTRACE2."

Rscript - <<RSCRIPT
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(readr)
  library(tibble)
})

project_root <- "${PROJECT_ROOT}"
rds <- "${ANNOTATED_RDS}"
expr_out <- "${EXPR_TSV_GZ_TMP}"
anno_out <- "${ANNO_TSV}"
meta_out <- "${META_TSV}"
table_dir <- "${TABLE_DIR}"

source(file.path(project_root, "config", "labels_colors.R"))
message("[R] Reading: ", rds)
seu <- readRDS(rds)
stopifnot(inherits(seu, "Seurat"))

metadata_raw <- seu@meta.data

if ("dataset" %in% colnames(metadata_raw) && any(as.character(metadata_raw\$dataset) == "fetal", na.rm = TRUE)) {
  cells <- rownames(metadata_raw)[as.character(metadata_raw\$dataset) == "fetal"]
} else {
  cells <- colnames(seu)
}

message("[R] Cells exported: ", length(cells))

required_cols <- c("gestational_week", "fine_subcluster", "final_annotation")
missing_cols <- setdiff(required_cols, colnames(metadata_raw))
if (length(missing_cols) > 0) {
  stop("[R] Missing required metadata columns: ", paste(missing_cols, collapse = ", "))
}

if ("major_cell_type" %in% colnames(metadata_raw)) {
  major_col <- "major_cell_type"
} else if ("cell_type" %in% colnames(metadata_raw)) {
  major_col <- "cell_type"
} else {
  stop("[R] Missing major cell-type column. Expected major_cell_type or cell_type.")
}

metadata <- metadata_raw[cells, , drop = FALSE] |>
  rownames_to_column("cell_id") |>
  transmute(
    cell_id = cell_id,
    gestational_week = as.numeric(gestational_week),
    major_cell_type = as.character(.data[[major_col]]),
    fine_subcluster = as.character(fine_subcluster),
    final_annotation = as.character(final_annotation)
  ) |>
  mutate(
    major_cell_type = recode(major_cell_type, quiescent = "degenerated")
  ) |>
  filter(
    !is.na(cell_id),
    cell_id != "",
    !is.na(major_cell_type),
    major_cell_type != "",
    !is.na(fine_subcluster),
    fine_subcluster != "",
    !is.na(final_annotation),
    final_annotation != ""
  )

extra_cell_types <- setdiff(unique(metadata\$major_cell_type), celltype_order)
if (length(extra_cell_types) > 0) {
  stop("[R] Unexpected major cell types: ", paste(extra_cell_types, collapse = ", "))
}

write_tsv(metadata, meta_out)
message("[R] Wrote metadata: ", meta_out)

annotation <- metadata |>
  transmute(cell_id = cell_id, label = fine_subcluster)
write_tsv(annotation, anno_out)
message("[R] Wrote CytoTRACE2 annotation: ", anno_out)

DefaultAssay(seu) <- "RNA"

rna_layers <- tryCatch(Layers(seu[["RNA"]]), error = function(e) character())
preferred_layers <- c("counts", "counts.SeuratProject")
counts_layer <- preferred_layers[preferred_layers %in% rna_layers][1]

if (is.na(counts_layer)) {
  counts_layers <- grep("^counts", rna_layers, value = TRUE)
  if (length(counts_layers) == 1) {
    counts_layer <- counts_layers
  }
}

if (is.na(counts_layer) || length(counts_layer) == 0) {
  stop("[R] Could not identify a raw-count RNA layer. Available RNA layers: ", paste(rna_layers, collapse = ", "))
}

message("[R] Using RNA count layer: ", counts_layer)
counts <- GetAssayData(seu, assay = "RNA", layer = counts_layer)

common_cells <- intersect(metadata\$cell_id, colnames(counts))
if (length(common_cells) != nrow(metadata)) {
  stop("[R] Metadata/count matrix cell mismatch. Metadata cells: ", nrow(metadata),
       "; count-layer cells matched: ", length(common_cells))
}

metadata <- metadata[match(common_cells, metadata\$cell_id), , drop = FALSE]
counts <- counts[, metadata\$cell_id, drop = FALSE]

if (!inherits(counts, "dgCMatrix")) {
  counts <- as(counts, "dgCMatrix")
}

message("[R] Count matrix dimensions, genes x cells: ", nrow(counts), " x ", ncol(counts))

write_tsv(
  tibble(
    input_rds = rds,
    assay = "RNA",
    counts_layer = counts_layer,
    genes = nrow(counts),
    cells = ncol(counts),
    metadata_file = meta_out,
    annotation_file = anno_out,
    expression_file = expr_out
  ),
  file.path(table_dir, "cytotrace2_fetal_input_summary.tsv")
)

message("[R] Writing gzipped dense expression table for CytoTRACE2.")
message("[R] This is expected to be large and may take a while.")

genes <- rownames(counts)
cells <- colnames(counts)

con <- gzfile(expr_out, open = "wt")
on.exit(close(con), add = TRUE)

writeLines(paste(c("gene", cells), collapse = "\t"), con)

chunk_size <- 200L
ngenes <- nrow(counts)
starts <- seq.int(1L, ngenes, by = chunk_size)

for (start in starts) {
  end <- min(start + chunk_size - 1L, ngenes)
  block <- as.matrix(counts[start:end, , drop = FALSE])
  rownames(block) <- genes[start:end]

  out_block <- cbind(gene = rownames(block), block)
  write.table(
    out_block,
    con,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  if ((end %% 2000L) == 0L || end == ngenes) {
    message("[R] Wrote genes: ", end, "/", ngenes)
  }
}

close(con)
con <- NULL

message("[R] Export complete.")
RSCRIPT

echo "[INFO] Validating exported gzip expression file."
gzip -t "${EXPR_TSV_GZ_TMP}"
mv -f "${EXPR_TSV_GZ_TMP}" "${EXPR_TSV_GZ}"
echo "[INFO] Expression file validated and moved to: ${EXPR_TSV_GZ}"
fi

echo "[INFO] Running CytoTRACE2."
CT2_RESULTS_FILE="${CT2_DIR}/cytotrace2_results.txt"

if [[ -s "${CT2_RESULTS_FILE}" ]]; then
  echo "[INFO] Existing CytoTRACE2 results found; skipping CytoTRACE2 run."
else
"${CYTOTRACE2_BIN}" \
    -f "${EXPR_TSV_GZ}" \
    -a "${ANNO_TSV}" \
    -sp "${CYTOTRACE2_SPECIES}" \
    -o "${CT2_DIR}" \
    -bs "${CYTOTRACE2_BATCH_SIZE}" \
    -sbs "${CYTOTRACE2_SUB_BATCH_SIZE}" \
    -mc "${CYTOTRACE2_THREADS}" \
    -dpl \
    -dpa
fi

echo "[INFO] CytoTRACE2 output files:"
find "${CT2_DIR}" -maxdepth 2 -type f | sort

echo "[INFO] Generating manuscript-ready CytoTRACE2 plots."

Rscript - <<RSCRIPT
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  library(scales)
  library(cowplot)
})

set.seed(42)

project_root <- "${PROJECT_ROOT}"
ct2_results_file <- file.path("${CT2_DIR}", "cytotrace2_results.txt")
metadata_file <- "${META_TSV}"
figure_dir <- "${PLOT_DIR}"
table_dir <- "${TABLE_DIR}"
log_dir <- "${LOG_DIR}"

source(file.path(project_root, "config", "plotting.R"))
source(file.path(project_root, "config", "labels_colors.R"))

fine_subcluster_colors <- c(
  germ_1 = "#F8766D",
  germ_2 = "#C49A00",
  germ_3 = "#53B400",
  germ_4 = "#00C094",
  germ_5 = "#00B6EB",
  germ_7 = "#A58AFF",
  germ_8 = "#FB61D7",
  germ_0 = "#6B6B6B",
  germ_6 = "#B0B0B0",

  granulosa_0 = "#F8766D",
  granulosa_1 = "#D39200",
  granulosa_2 = "#93AA00",
  granulosa_3 = "#00BA38",
  granulosa_4 = "#00C19F",
  granulosa_5 = "#00B9E3",
  granulosa_6 = "#619CFF",
  granulosa_7 = "#DB72FB",
  granulosa_8 = "#FF61C3",

  stroma_0 = "#F8766D",
  stroma_1 = "#A3A500",
  stroma_2 = "#00BF7D",
  stroma_3 = "#00B0F6",
  stroma_4 = "#E76BF3",

  endothelial_0 = "#E7298A",

  mural_0 = "#66A61E",
  mural_1 = "#B2DF8A",

  immune_0 = "#E6AB02",
  immune_1 = "#FDBF6F",

  erythroid_0 = "#A6761D",
  erythroid_1 = "#D95F0E"
)
merged_cell_types <- c("endothelial", "mural", "erythroid", "immune")

plot_definitions <- list(
  list(
    plot_id = "germ_including_degenerated",
    output_prefix = "cytotrace2_violin_germ_including_degenerated",
    include_cell_types = c("germ", "degenerated"),
    title = "Germ and degenerating follicle-cell states",
    block_by_cell_type = FALSE
  ),
  list(
    plot_id = "granulosa",
    output_prefix = "cytotrace2_violin_granulosa",
    include_cell_types = c("granulosa"),
    title = "Granulosa states",
    block_by_cell_type = FALSE
  ),
  list(
    plot_id = "stroma",
    output_prefix = "cytotrace2_violin_stroma",
    include_cell_types = c("stroma"),
    title = "Stromal states",
    block_by_cell_type = FALSE
  ),
  list(
    plot_id = "merged_endothelial_mural_immune_erythroid",
    output_prefix = "cytotrace2_violin_merged_endothelial_mural_immune_erythroid",
    include_cell_types = merged_cell_types,
    title = "Endothelial, mural, immune, and erythroid states",
    block_by_cell_type = TRUE
  )
)

display_annotation <- function(fine_subcluster, final_annotation) {
  case_when(
    fine_subcluster == "germ_0" ~ "Atresia-stressed",
    fine_subcluster == "germ_6" ~ "Clearance-associated",
    TRUE ~ final_annotation
  )
}

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

make_fallback_colors <- function(missing_subclusters) {
  if (length(missing_subclusters) == 0) {
    return(character(0))
  }

  fallback <- grDevices::hcl.colors(length(missing_subclusters), palette = "Dynamic")
  names(fallback) <- missing_subclusters
  fallback
}

read_cytotrace2_results <- function(file) {
  stopifnot(file.exists(file))
  raw <- read_tsv(file, show_col_types = FALSE)
  lower_names <- setNames(colnames(raw), tolower(colnames(raw)))

  cell_col <- lower_names[intersect(
    c("cell_id", "cell", "cellid", "barcode"),
    names(lower_names)
  )][1]

  if (is.na(cell_col)) {
    cell_col <- colnames(raw)[1]
  }

  score_col <- lower_names[intersect(
    c(
      "ct2_score",
      "cytotrace2_score",
      "cytotrace_score",
      "score",
      "potency",
      "developmental_potential"
    ),
    names(lower_names)
  )][1]

  if (is.na(score_col)) {
    candidate_score_cols <- colnames(raw)[
      grepl("score|poten|cytotrace", colnames(raw), ignore.case = TRUE)
    ]
    score_col <- candidate_score_cols[1]
  }

  if (is.na(score_col) || length(score_col) == 0) {
    stop("Could not identify CytoTRACE2 score column. Columns: ", paste(colnames(raw), collapse = ", "))
  }

  raw |>
    transmute(
      cell_id = as.character(.data[[cell_col]]),
      cytotrace2_score = as.numeric(.data[[score_col]])
    ) |>
    filter(!is.na(cell_id), cell_id != "", !is.na(cytotrace2_score))
}

make_plot_table <- function(df, plot_definition) {
  plot_data <- df |>
    filter(major_cell_type %in% plot_definition\$include_cell_types) |>
    mutate(
      display_annotation = display_annotation(fine_subcluster, final_annotation)
    )

  if (nrow(plot_data) == 0) {
    return(NULL)
  }

  order_table <- plot_data |>
    group_by(major_cell_type, fine_subcluster, display_annotation) |>
    summarise(
      n_cells = n(),
      mean_cytotrace2_score = mean(cytotrace2_score, na.rm = TRUE),
      .groups = "drop"
    )

  if (isTRUE(plot_definition\$block_by_cell_type)) {
    order_table <- order_table |>
      mutate(
        major_cell_type = factor(as.character(major_cell_type), levels = plot_definition\$include_cell_types)
      ) |>
      arrange(major_cell_type, desc(mean_cytotrace2_score), fine_subcluster)
  } else {
    order_table <- order_table |>
      arrange(desc(mean_cytotrace2_score), fine_subcluster)
  }

  order_table <- order_table |>
    mutate(
      group_id = fine_subcluster,
      axis_label = paste0(display_annotation, "\n(n=", n_cells, ")"),
      axis_label = factor(axis_label, levels = axis_label)
    )

  plot_data |>
    left_join(
      order_table |>
        select(
          major_cell_type,
          fine_subcluster,
          display_annotation,
          n_cells,
          mean_cytotrace2_score,
          group_id,
          axis_label
        ),
      by = c("major_cell_type", "fine_subcluster", "display_annotation")
    ) |>
    mutate(
      axis_label = factor(axis_label, levels = levels(order_table\$axis_label)),
      plot_id = plot_definition\$plot_id
    )
}
make_violin_plot <- function(plot_data, plot_definition) {
  n_groups <- length(unique(plot_data\$axis_label))

  panel_width <- max(3.8, 0.42 * n_groups + 1.0)
  panel_height <- 2.65
  label_height <- 0.85
  figure_height <- panel_height + label_height

  plot_subclusters <- unique(as.character(plot_data\$fine_subcluster))
  missing_colors <- setdiff(plot_subclusters, names(fine_subcluster_colors))
  plot_subcluster_colors <- c(
    fine_subcluster_colors,
    make_fallback_colors(missing_colors)
  )
  plot_subcluster_colors <- plot_subcluster_colors[plot_subclusters]

  axis_levels <- levels(plot_data\$axis_label)
  label_data <- tibble(
    axis_label = factor(axis_levels, levels = axis_levels),
    x = seq_along(axis_levels),
    y = 1
  )

  panel_plot <- ggplot(
    plot_data,
    aes(
      x = axis_label,
      y = cytotrace2_score,
      fill = fine_subcluster
    )
  ) +
    geom_violin(
      scale = "width",
      trim = TRUE,
      color = "black",
      linewidth = 0.25
    ) +
    stat_summary(
      fun = mean,
      geom = "point",
      size = 0.8,
      color = "black"
    ) +
    scale_fill_manual(values = plot_subcluster_colors, drop = FALSE, guide = "none") +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.25),
      expand = expansion(mult = c(0.02, 0.04))
    ) +
    labs(
      x = NULL,
      y = "CytoTRACE2 developmental potential"
    ) +
    theme_publication() +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_line(linewidth = 0.25),
      axis.text.y = element_text(
        family = publication_font_family,
        size = publication_base_size
      ),
      legend.position = "none",
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.25),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.margin = margin(4, 4, 0, 4)
    )

  label_plot <- ggplot(label_data, aes(x = axis_label, y = y, label = axis_label)) +
    geom_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      family = publication_font_family,
      size = publication_base_size / ggplot2::.pt
    ) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void(base_family = publication_font_family, base_size = publication_base_size) +
    theme(
      plot.margin = margin(0, 4, 4, 4)
    )

  combined_plot <- cowplot::plot_grid(
    panel_plot,
    label_plot,
    ncol = 1,
    rel_heights = c(panel_height, label_height),
    align = "v",
    axis = "lr"
  )

  save_publication_plot(
    combined_plot,
    file.path(figure_dir, paste0(plot_definition\$output_prefix, ".png")),
    width = panel_width,
    height = figure_height,
    dpi = 600
  )

  combined_plot
}

stopifnot(file.exists(metadata_file))

ct2 <- read_cytotrace2_results(ct2_results_file)
metadata <- read_tsv(metadata_file, show_col_types = FALSE)

required_metadata_cols <- c(
  "cell_id",
  "gestational_week",
  "major_cell_type",
  "fine_subcluster",
  "final_annotation"
)

missing_metadata_cols <- setdiff(required_metadata_cols, colnames(metadata))
if (length(missing_metadata_cols) > 0) {
  stop("Missing metadata columns: ", paste(missing_metadata_cols, collapse = ", "))
}

plot_metadata <- metadata |>
  mutate(
    cell_id = as.character(cell_id),
    major_cell_type = as.character(major_cell_type),
    major_cell_type = recode(major_cell_type, quiescent = "degenerated"),
    fine_subcluster = as.character(fine_subcluster),
    final_annotation = as.character(final_annotation)
  ) |>
  inner_join(ct2, by = "cell_id")

if (nrow(plot_metadata) == 0) {
  stop("No cells overlapped between metadata and CytoTRACE2 results.")
}
write_tsv(
  plot_metadata,
  file.path(table_dir, "cytotrace2_fetal_results_with_metadata.tsv")
)

plot_tables <- lapply(
  plot_definitions,
  function(plot_definition) {
    plot_data <- make_plot_table(plot_metadata, plot_definition)

    if (is.null(plot_data)) {
      message("Skipping ", plot_definition\$plot_id, ": no cells found.")
      return(NULL)
    }

    message("Plotting: ", plot_definition\$plot_id)
    make_violin_plot(plot_data, plot_definition)

    write_tsv(
      plot_data,
      file.path(table_dir, paste0(plot_definition\$output_prefix, "_plot_data.tsv"))
    )

    plot_data
  }
)

plot_tables <- plot_tables[!vapply(plot_tables, is.null, logical(1))]
all_plot_data <- bind_rows(plot_tables)

write_tsv(
  all_plot_data,
  file.path(table_dir, "cytotrace2_fetal_all_plot_data.tsv")
)

plot_summary <- all_plot_data |>
  group_by(plot_id, major_cell_type, fine_subcluster, final_annotation) |>
  summarise(
    n_cells = n(),
    mean_cytotrace2_score = mean(cytotrace2_score, na.rm = TRUE),
    median_cytotrace2_score = median(cytotrace2_score, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(plot_id, major_cell_type, desc(mean_cytotrace2_score))

write_tsv(
  plot_summary,
  file.path(table_dir, "cytotrace2_fetal_plot_summary.tsv")
)

summary_lines <- c(
  paste("CytoTRACE2 results:", ct2_results_file),
  paste("Metadata:", metadata_file),
  paste("Cells with CytoTRACE2 scores:", nrow(plot_metadata)),
  paste("Plots generated:", length(plot_tables)),
  paste("Plot IDs:", paste(vapply(plot_definitions, function(x) x[["plot_id"]], character(1)), collapse = ", ")),
  "Major cell type labels use degenerated, not quiescent.",
  "Figure-only abbreviated labels: germ_0 = Atresia-stressed; germ_6 = Clearance-associated.",
  paste("Output directory:", "${OUT_ROOT}")
)

writeLines(
  summary_lines,
  file.path(log_dir, "cytotrace2_fetal_plotting_summary.txt")
)

sink(file.path(log_dir, "cytotrace2_fetal_plotting_sessionInfo.txt"))
print(sessionInfo())
sink()

cat("[R] CytoTRACE2 plotting complete.\n")
RSCRIPT

echo "[INFO] DONE."
echo "[INFO] Outputs written to: ${OUT_ROOT}"
