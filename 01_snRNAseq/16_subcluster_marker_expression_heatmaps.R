#!/usr/bin/env Rscript

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "dplyr",
  "readr",
  "tibble",
  "stringr",
  "pheatmap",
  "grid"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(pheatmap)
  library(grid)
})

set.seed(42)
options(bitmapType = "cairo")

project_root <- normalizePath(
  Sys.getenv("PROJECT_ROOT", unset = getwd()),
  mustWork = TRUE
)

if (!file.exists(file.path(project_root, "config", "plotting.R"))) {
  stop("PROJECT_ROOT does not point to the repository root. Run from the repo root or set PROJECT_ROOT.")
}

paths_local <- file.path(project_root, "config", "paths_local.R")
if (file.exists(paths_local)) {
  source(paths_local)
}

results_base <- if (exists("results_root", inherits = FALSE)) {
  results_root
} else {
  Sys.getenv(
    "FETAL_OVARY_RESULTS_ROOT",
    unset = file.path(dirname(project_root), paste0(basename(project_root), "_results"))
  )
}
results_base <- normalizePath(results_base, mustWork = FALSE)

source(file.path(project_root, "config", "plotting.R"))
source(file.path(project_root, "config", "labels_colors.R"))

if (!exists("publication_font_family", inherits = TRUE)) {
  publication_font_family <- "Helvetica"
}
if (!exists("publication_base_size", inherits = TRUE)) {
  publication_base_size <- 8
}

annotated_rds <- if (exists("fetal_annotated_rds", inherits = FALSE)) {
  fetal_annotated_rds
} else {
  file.path(
    results_base,
    "snRNAseq_annotated_object",
    "objects",
    "fetal_ovary_snRNAseq_canonical_annotated.rds"
  )
}

deg_results_root <- if (exists("snrna_deg_results_root", inherits = FALSE)) {
  snrna_deg_results_root
} else {
  file.path(results_base, "snRNAseq_DEGs_celltype_subcluster")
}

out_root <- if (exists("snrna_marker_heatmap_results_root", inherits = FALSE)) {
  snrna_marker_heatmap_results_root
} else {
  file.path(results_base, "snRNAseq_subcluster_marker_expression_heatmaps")
}

out_figure_dir <- file.path(out_root, "figures")
out_table_dir <- file.path(out_root, "tables")
out_log_dir <- file.path(out_root, "logs")

dir.create(out_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_log_dir, recursive = TRUE, showWarnings = FALSE)

assay_use <- "RNA"

target_cell_types <- trimws(unlist(strsplit(
  Sys.getenv("MARKER_HEATMAP_CELL_TYPES", unset = "germ,granulosa,stroma"),
  ","
)))
target_cell_types <- target_cell_types[target_cell_types != ""]

top_genes_per_subcluster <- as.integer(Sys.getenv("MARKER_HEATMAP_TOP_GENES_PER_SUBCLUSTER", unset = "10"))
max_genes_per_heatmap <- as.integer(Sys.getenv("MARKER_HEATMAP_MAX_GENES", unset = "75"))

z_score_cap <- as.numeric(Sys.getenv("MARKER_HEATMAP_Z_CAP", unset = "2.5"))
heatmap_width <- as.numeric(Sys.getenv("MARKER_HEATMAP_WIDTH", unset = "8.5"))
heatmap_min_height <- as.numeric(Sys.getenv("MARKER_HEATMAP_MIN_HEIGHT", unset = "5.5"))
heatmap_row_height <- as.numeric(Sys.getenv("MARKER_HEATMAP_ROW_HEIGHT", unset = "0.13"))
heatmap_dpi <- as.integer(Sys.getenv("MARKER_HEATMAP_DPI", unset = "600"))

cluster_rows <- tolower(Sys.getenv("MARKER_HEATMAP_CLUSTER_ROWS", unset = "true")) %in% c("true", "t", "1", "yes", "y")
cluster_cols <- tolower(Sys.getenv("MARKER_HEATMAP_CLUSTER_COLUMNS", unset = "false")) %in% c("true", "t", "1", "yes", "y")

label_mode <- Sys.getenv("MARKER_HEATMAP_LABEL_MODE", unset = "annotation")
label_wrap_width <- as.integer(Sys.getenv("MARKER_HEATMAP_LABEL_WRAP_WIDTH", unset = "18"))

subcluster_top_table <- file.path(
  deg_results_root,
  "tables",
  "subcluster_top50",
  "DEG_subcluster__top50_adaptive_thresholds_combined.csv"
)

stopifnot(file.exists(annotated_rds))
if (!file.exists(subcluster_top_table)) {
  stop(
    "Adaptive subcluster top DEG table not found: ", subcluster_top_table,
    "\nRun 01_snRNAseq/14_generate_deg_tables_celltype_subcluster.R first."
  )
}

safe_file_component <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

mode_string <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) {
    return(NA_character_)
  }
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

derive_cell_type <- function(meta) {
  if ("major_cell_type" %in% colnames(meta)) {
    ct <- as.character(meta$major_cell_type)
  } else if ("cell_type_quiescent" %in% colnames(meta)) {
    ct <- as.character(meta$cell_type_quiescent)
  } else if ("cell_type" %in% colnames(meta)) {
    ct <- as.character(meta$cell_type)
  } else {
    fs <- as.character(meta$fine_subcluster)
    ct <- sub("_.*$", "", fs)
  }

  ct[ct == "quiescent"] <- "degenerated"
  ct
}

natural_subcluster_order <- function(x) {
  x <- unique(as.character(x))
  prefix <- sub("_[0-9]+$", "", x)
  number <- suppressWarnings(as.integer(sub("^.*_([0-9]+)$", "\\1", x)))

  tibble(x = x, prefix = prefix, number = number) |>
    arrange(prefix, number, x) |>
    pull(x)
}
order_subclusters <- function(cell_type, subclusters) {
  subclusters <- unique(as.character(subclusters))

  candidate_order <- character(0)

  if (exists("fine_subcluster_order", inherits = TRUE)) {
    candidate_order <- get("fine_subcluster_order", inherits = TRUE)
  } else if (exists("subcluster_order", inherits = TRUE)) {
    candidate_order <- get("subcluster_order", inherits = TRUE)
  }

  candidate_order <- as.character(candidate_order)
  candidate_order <- candidate_order[startsWith(candidate_order, paste0(cell_type, "_"))]

  if (length(candidate_order) > 0) {
    c(intersect(candidate_order, subclusters), setdiff(natural_subcluster_order(subclusters), candidate_order))
  } else {
    natural_subcluster_order(subclusters)
  }
}

shorten_label <- function(x) {
  dplyr::recode(
    as.character(x),
    "Atresia-stressed degenerating follicle cells" = "Atresia-stressed",
    "Clearance-associated degenerating follicle cells" = "Clearance-associated",
    .default = as.character(x)
  )
}

make_column_labels <- function(label_df, subcluster_order) {
  label_df <- label_df |>
    filter(fine_subcluster %in% subcluster_order) |>
    mutate(fine_subcluster = factor(fine_subcluster, levels = subcluster_order)) |>
    arrange(fine_subcluster)

  if (label_mode == "id") {
    labels <- as.character(label_df$fine_subcluster)
  } else if (label_mode == "both") {
    labels <- paste0(label_df$fine_subcluster, ": ", label_df$display_label)
  } else {
    labels <- label_df$display_label
  }

  labels <- stringr::str_wrap(labels, width = label_wrap_width)
  names(labels) <- as.character(label_df$fine_subcluster)
  labels
}

read_subcluster_top_table <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE)

  required <- c(
    "cell_type",
    "fine_subcluster",
    "final_annotation",
    "gene",
    "avg_log2FC",
    "pct.1",
    "pct.2",
    "p_val_adj",
    "score"
  )

  missing_required <- setdiff(required, colnames(df))
  if (length(missing_required) > 0) {
    stop(
      "Adaptive subcluster top DEG table is missing required columns: ",
      paste(missing_required, collapse = ", ")
    )
  }

  df |>
    mutate(
      top_table_row = row_number(),
      gene = as.character(gene),
      cell_type = as.character(cell_type),
      fine_subcluster = as.character(fine_subcluster),
      final_annotation = as.character(final_annotation),
      avg_log2FC = as.numeric(avg_log2FC),
      pct.1 = as.numeric(pct.1),
      pct.2 = as.numeric(pct.2),
      p_val_adj = as.numeric(p_val_adj),
      score = as.numeric(score)
    )
}
read_celltype_deg_tables <- function(cell_type) {
  out <- subcluster_top_combined |>
    filter(cell_type == !!cell_type)

  if (nrow(out) == 0) {
    warning("Skipping ", cell_type, ": no rows found in adaptive subcluster top DEG table.")
  }

  out
}

select_marker_genes <- function(deg_tbl, cell_type) {
  subclusters <- order_subclusters(cell_type, unique(deg_tbl$fine_subcluster))

  selected <- deg_tbl |>
    filter(
      !is.na(gene),
      !is.na(avg_log2FC),
      !is.na(pct.1),
      !is.na(pct.2),
      !is.na(p_val_adj),
      !is.na(score)
    ) |>
    mutate(selection_subcluster = fine_subcluster) |>
    group_by(fine_subcluster) |>
    arrange(top_table_row, .by_group = TRUE) |>
    slice_head(n = top_genes_per_subcluster) |>
    ungroup()

  if (nrow(selected) == 0) {
    return(tibble())
  }

  selected |>
    group_by(gene) |>
    arrange(top_table_row, .by_group = TRUE) |>
    slice_head(n = 1) |>
    ungroup() |>
    arrange(
      factor(selection_subcluster, levels = subclusters),
      top_table_row,
      gene
    ) |>
    mutate(selection_rank = row_number()) |>
    slice_head(n = max_genes_per_heatmap)
}


get_data_matrix <- function(obj, assay = "RNA") {
  if (inherits(obj[[assay]], "Assay5")) {
    obj <- tryCatch(
      JoinLayers(obj, assay = assay),
      error = function(e) obj
    )
  }

  layers_available <- tryCatch(
    Layers(obj[[assay]]),
    error = function(e) character(0)
  )

  if (!any(grepl("^data", layers_available))) {
    obj <- NormalizeData(obj, assay = assay, verbose = FALSE)
  }

  data_mat <- tryCatch(
    GetAssayData(obj, assay = assay, layer = "data"),
    error = function(e) GetAssayData(obj, assay = assay, slot = "data")
  )

  data_mat
}
average_expression_by_subcluster <- function(data_mat, meta, genes, subclusters) {
  genes <- intersect(genes, rownames(data_mat))
  subclusters <- unique(as.character(subclusters))

  if (length(genes) == 0) {
    stop("No selected genes were found in the RNA assay.")
  }

  avg_list <- lapply(subclusters, function(sc) {
    cells <- meta$cell_id[meta$fine_subcluster == sc]
    cells <- intersect(cells, colnames(data_mat))

    if (length(cells) == 0) {
      rep(NA_real_, length(genes))
    } else {
      Matrix::rowMeans(data_mat[genes, cells, drop = FALSE])
    }
  })

  avg_mat <- do.call(cbind, avg_list)
  rownames(avg_mat) <- genes
  colnames(avg_mat) <- subclusters

  avg_mat
}

row_zscore <- function(mat) {
  z <- t(scale(t(mat)))
  z[is.na(z)] <- 0
  z <- pmax(pmin(z, z_score_cap), -z_score_cap)
  z
}

save_pheatmap_dual <- function(ph, pdf_out, png_out, width, height, dpi) {
  pdf(
    pdf_out,
    width = width,
    height = height,
    family = publication_font_family,
    useDingbats = FALSE
  )
  grid::grid.newpage()
  grid::grid.draw(ph$gtable)
  dev.off()

  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(
      png_out,
      width = width,
      height = height,
      units = "in",
      res = dpi,
      background = "white"
    )
  } else {
    png(
      png_out,
      width = width,
      height = height,
      units = "in",
      res = dpi,
      bg = "white"
    )
  }

  grid::grid.newpage()
  grid::grid.draw(ph$gtable)
  dev.off()
}

plot_celltype_heatmap <- function(z_mat, column_labels, cell_type) {
  heatmap_colors <- grDevices::colorRampPalette(
    c("#2166AC", "#F7F7F7", "#B2182B")
  )(101)

  heatmap_breaks <- seq(
    -z_score_cap,
    z_score_cap,
    length.out = length(heatmap_colors) + 1
  )

  figure_height <- max(
    heatmap_min_height,
    2.2 + nrow(z_mat) * heatmap_row_height
  )

  title_text <- paste0(
    "Subcluster marker expression - ", cell_type,
    "\nrow-scaled average normalized expression"
  )

  ph <- pheatmap::pheatmap(
    z_mat,
    color = heatmap_colors,
    breaks = heatmap_breaks,
    cluster_rows = cluster_rows && nrow(z_mat) > 1,
    cluster_cols = cluster_cols && ncol(z_mat) > 1,
    labels_col = column_labels[colnames(z_mat)],
    border_color = NA,
    fontsize = publication_base_size,
    fontsize_row = max(4, publication_base_size - 2),
    fontsize_col = max(5, publication_base_size - 1),
    angle_col = 90,
    main = title_text,
    silent = TRUE
  )

  pdf_out <- file.path(
    out_figure_dir,
    paste0("Subcluster_marker_expression_heatmap_", safe_file_component(cell_type), ".pdf")
  )

  png_out <- file.path(
    out_figure_dir,
    paste0("Subcluster_marker_expression_heatmap_", safe_file_component(cell_type), ".png")
  )

  save_pheatmap_dual(
    ph = ph,
    pdf_out = pdf_out,
    png_out = png_out,
    width = heatmap_width,
    height = figure_height,
    dpi = heatmap_dpi
  )
}

cat("Loading annotated snRNA-seq object:\n  ", annotated_rds, "\n", sep = "")
obj <- readRDS(annotated_rds)
stopifnot(inherits(obj, "Seurat"))
stopifnot("fine_subcluster" %in% colnames(obj@meta.data))

DefaultAssay(obj) <- assay_use

meta <- obj@meta.data |>
  rownames_to_column("cell_id")

meta$fine_subcluster <- as.character(meta$fine_subcluster)
meta$cell_type_heatmap <- derive_cell_type(meta)

if ("final_annotation" %in% colnames(meta)) {
  meta$final_annotation <- as.character(meta$final_annotation)
} else {
  meta$final_annotation <- as.character(meta$fine_subcluster)
}

label_map <- meta |>
  group_by(cell_type_heatmap, fine_subcluster) |>
  summarise(
    display_label = mode_string(final_annotation),
    n_cells = n(),
    .groups = "drop"
  ) |>
  mutate(display_label = shorten_label(display_label))

write_csv(
  label_map,
  file.path(out_table_dir, "subcluster_marker_heatmap_label_map.csv")
)
subcluster_top_combined <- read_subcluster_top_table(subcluster_top_table)

write_csv(
  subcluster_top_combined |>
    count(cell_type, fine_subcluster, final_annotation, name = "n_top_table_rows") |>
    arrange(cell_type, fine_subcluster),
  file.path(out_table_dir, "QC_adaptive_subcluster_top_table_rows.csv")
)

data_mat <- get_data_matrix(obj, assay = assay_use)

summary_rows <- list()
selected_gene_tables <- list()

for (cell_type in target_cell_types) {
  cat("\n============================================================\n")
  cat("Generating marker-expression heatmap: ", cell_type, "\n", sep = "")
  cat("============================================================\n")

  deg_tbl <- read_celltype_deg_tables(cell_type)

  if (nrow(deg_tbl) == 0) {
    warning("No DEG rows available for cell type: ", cell_type)
    next
  }

  markers <- select_marker_genes(deg_tbl, cell_type)

  if (nrow(markers) == 0) {
    warning("No marker genes passed filters for cell type: ", cell_type)
    next
  }

  subclusters <- order_subclusters(cell_type, unique(markers$fine_subcluster))

  label_df <- label_map |>
    filter(cell_type_heatmap == cell_type, fine_subcluster %in% subclusters)

  column_labels <- make_column_labels(label_df, subclusters)

  genes <- markers$gene
  genes <- genes[genes %in% rownames(data_mat)]

  avg_mat <- average_expression_by_subcluster(
    data_mat = data_mat,
    meta = meta,
    genes = genes,
    subclusters = subclusters
  )

  avg_mat <- avg_mat[genes, subclusters, drop = FALSE]
  z_mat <- row_zscore(avg_mat)

  plot_celltype_heatmap(
    z_mat = z_mat,
    column_labels = column_labels,
    cell_type = cell_type
  )

  avg_out <- as.data.frame(avg_mat) |>
    rownames_to_column("gene")

  z_out <- as.data.frame(z_mat) |>
    rownames_to_column("gene")

  write_csv(
    avg_out,
    file.path(out_table_dir, paste0("average_expression_matrix_", safe_file_component(cell_type), ".csv"))
  )

  write_csv(
    z_out,
    file.path(out_table_dir, paste0("row_zscore_matrix_", safe_file_component(cell_type), ".csv"))
  )

  markers_out <- markers |>
    mutate(cell_type_heatmap = cell_type) |>
    relocate(cell_type_heatmap)

  write_csv(
    markers_out,
    file.path(out_table_dir, paste0("selected_marker_genes_", safe_file_component(cell_type), ".csv"))
  )

  selected_gene_tables[[cell_type]] <- markers_out

  summary_rows[[cell_type]] <- tibble(
    cell_type = cell_type,
    n_subclusters = length(subclusters),
    n_selected_genes = length(genes),
    n_adaptive_top_table_rows_available = nrow(deg_tbl),
    top_genes_per_subcluster = top_genes_per_subcluster,
    max_genes_per_heatmap = max_genes_per_heatmap
  )
}

if (length(selected_gene_tables) > 0) {
  write_csv(
    bind_rows(selected_gene_tables),
    file.path(out_table_dir, "COMBINED_selected_marker_genes_for_heatmaps.csv")
  )
}

if (length(summary_rows) > 0) {
  write_csv(
    bind_rows(summary_rows),
    file.path(out_table_dir, "subcluster_marker_heatmap_summary.csv")
  )
}

writeLines(
  c(
    paste("Input annotated snRNA-seq object:", annotated_rds),
    paste("Input adaptive subcluster top DEG table:", subcluster_top_table),
    paste("Output root:", out_root),
    paste("Target cell types:", paste(target_cell_types, collapse = ", ")),
    paste("Assay:", assay_use),
    paste("Top genes per subcluster:", top_genes_per_subcluster),
    paste("Maximum genes per heatmap:", max_genes_per_heatmap),
    paste("Z-score cap:", z_score_cap),
    paste("Column label mode:", label_mode),
    "Heatmap values are row-scaled average normalized expression by fine subcluster.",
    "Marker genes are selected directly from the adaptive subcluster top DEG table generated by script 14.",
    "This preserves the script-14 marker-selection logic, including score ranking, adaptive thresholds, priority high-prevalence/high-score genes, and top-table exclusion rules."
  ),
  file.path(out_log_dir, "subcluster_marker_expression_heatmaps_notes.txt")
)

sink(file.path(out_log_dir, "sessionInfo_subcluster_marker_expression_heatmaps.txt"))
print(sessionInfo())
sink()

cat("\nDone. Marker-expression heatmaps written to:\n", out_root, "\n")
