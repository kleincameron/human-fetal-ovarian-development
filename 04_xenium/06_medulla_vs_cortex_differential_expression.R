#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(future)
})

set.seed(42)
options(bitmapType = "cairo")

future::plan("sequential")
options(future.globals.maxSize = 100 * 1024^3)

project_root <- normalizePath(
  Sys.getenv("PROJECT_ROOT", unset = getwd()),
  mustWork = TRUE
)

if (!file.exists(file.path(project_root, "config", "labels_colors.R"))) {
  stop("PROJECT_ROOT does not point to the repository root. Run from the repo root or set PROJECT_ROOT.")
}

source(file.path(project_root, "config", "labels_colors.R"))

paths_local <- file.path(project_root, "config", "paths_local.R")
if (file.exists(paths_local)) {
  source(paths_local)
}

data_root <- if (exists("data_root", inherits = FALSE)) {
  data_root
} else {
  Sys.getenv(
    "FETAL_OVARY_DATA_ROOT",
    unset = file.path(dirname(project_root), "github_code_for_publication_controlled_data", "HRA019091")
  )
}
data_root <- normalizePath(data_root, mustWork = FALSE)

results_base <- if (exists("results_root", inherits = FALSE)) {
  results_root
} else {
  Sys.getenv(
    "FETAL_OVARY_RESULTS_ROOT",
    unset = file.path(dirname(project_root), paste0(basename(project_root), "_results"))
  )
}
results_base <- normalizePath(results_base, mustWork = FALSE)

manifest_file <- if (exists("xenium_manifest", inherits = FALSE)) {
  xenium_manifest
} else {
  Sys.getenv(
    "XENIUM_SAMPLE_MANIFEST",
    unset = file.path(project_root, "config", "xenium_samples.csv")
  )
}

xenium_rds <- if (exists("xenium_annotated_rds", inherits = FALSE)) {
  xenium_annotated_rds
} else {
  file.path(
    results_base,
    "xenium_annotated_object",
    "objects",
    "fetal_ovary_xenium_annotated.rds"
  )
}

snrna_annotation_metadata <- if (exists("snrna_annotation_metadata", inherits = FALSE)) {
  snrna_annotation_metadata
} else {
  file.path(
    project_root,
    "metadata",
    "snRNAseq_cell_annotations.csv"
  )
}

medulla_roi_csv <- if (exists("xenium_medulla_roi_csv", inherits = FALSE)) {
  xenium_medulla_roi_csv
} else {
  file.path(
    project_root,
    "metadata",
    "xenium_medulla_roi_coordinates.csv"
  )
}

results_root <- if (exists("xenium_medulla_cortex_de_results_root", inherits = FALSE)) {
  xenium_medulla_cortex_de_results_root
} else {
  file.path(
    results_base,
    "xenium_medulla_cortex_differential_expression"
  )
}

if (dir.exists(results_root)) {
  unlink(results_root, recursive = TRUE, force = TRUE)
}

out_table_dir <- file.path(results_root, "tables")
out_complete_dir <- file.path(out_table_dir, "complete_deg_tables")
out_log_dir <- file.path(results_root, "logs")

dir.create(out_complete_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_log_dir, recursive = TRUE, showWarnings = FALSE)

fallback_celltype_order <- c(
  "germ",
  "degenerated",
  "granulosa",
  "stroma",
  "endothelial",
  "mural",
  "immune",
  "erythroid"
)

if (!exists("celltype_order", inherits = FALSE)) {
  celltype_order <- fallback_celltype_order
}

celltype_order <- unique(c(celltype_order, fallback_celltype_order))

min_cells_per_region <- 100L
top_genes_per_region <- 50L
fdr_cutoff <- 0.05

resolve_data_path <- function(path, data_root) {
  path <- path.expand(as.character(path))

  if (is.na(path) || path == "") {
    return(path)
  }

  is_absolute <- grepl("^/", path) || grepl("^[A-Za-z]:[\\/]", path)

  if (is_absolute) {
    normalizePath(path, mustWork = FALSE)
  } else {
    normalizePath(file.path(data_root, path), mustWork = FALSE)
  }
}

normalize_include <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

safe_file_label <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

pick_column <- function(df, candidates, required = TRUE) {
  hit <- candidates[candidates %in% colnames(df)]

  if (length(hit) > 0) {
    return(hit[[1]])
  }

  if (isTRUE(required)) {
    stop("None of these columns were found: ", paste(candidates, collapse = ", "))
  }

  NA_character_
}

first_existing_column <- function(cols, candidates, label) {
  hit <- candidates[candidates %in% cols]

  if (length(hit) == 0) {
    stop("Could not find ", label, " column. Tried: ", paste(candidates, collapse = ", "))
  }

  hit[[1]]
}

read_polygon_coordinates <- function(path, label) {
  stopifnot(file.exists(path))

  for (skip_n in c(0, 1, 2, 3, 4)) {
    dat <- suppressWarnings(
      tryCatch(
        read_csv(path, skip = skip_n, show_col_types = FALSE, progress = FALSE),
        error = function(e) NULL
      )
    )

    if (is.null(dat) || ncol(dat) == 0) {
      next
    }

    x_col <- pick_column(
      dat,
      c("x", "X", "x_um", "X_um", "x_centroid", "x_centroid_um"),
      required = FALSE
    )
    y_col <- pick_column(
      dat,
      c("y", "Y", "y_um", "Y_um", "y_centroid", "y_centroid_um"),
      required = FALSE
    )

    if (!is.na(x_col) && !is.na(y_col)) {
      out <- tibble(
        x = suppressWarnings(as.numeric(dat[[x_col]])),
        y = suppressWarnings(as.numeric(dat[[y_col]]))
      ) |>
        filter(is.finite(x), is.finite(y))

      if (nrow(out) >= 3) {
        return(out)
      }
    }
  }

  stop("Could not read ", label, " polygon coordinates from: ", path)
}

point_in_polygon <- function(px, py, poly_x, poly_y) {
  n <- length(poly_x)
  inside <- rep(FALSE, length(px))
  j <- n

  for (i in seq_len(n)) {
    xi <- poly_x[i]
    yi <- poly_y[i]
    xj <- poly_x[j]
    yj <- poly_y[j]

    crosses <- ((yi > py) != (yj > py)) &
      (px < (xj - xi) * (py - yi) / (yj - yi + .Machine$double.eps) + xi)

    inside <- xor(inside, crosses)
    j <- i
  }

  inside
}

empty_de_result <- function(
  comparison,
  comparison_label,
  cell_type,
  cells_medulla,
  cells_cortex,
  status,
  reason
) {
  list(
    de = tibble(),
    summary = tibble(
      comparison = comparison,
      comparison_label = comparison_label,
      cell_type = cell_type,
      cells_medulla = cells_medulla,
      cells_cortex = cells_cortex,
      genes_tested = 0L,
      significant_genes_FDR_0_05 = 0L,
      medulla_enriched_genes_FDR_0_05 = 0L,
      cortex_enriched_genes_FDR_0_05 = 0L,
      complete_deg_table = NA_character_,
      status = status,
      reason = reason
    )
  )
}

run_de <- function(obj, comparison, comparison_label, cell_type = "all_cells") {
  cat("\n==============================\n")
  cat("DE comparison: ", comparison_label, "\n", sep = "")
  cat("==============================\n")

  group_vec <- obj$compartment
  keep_cells <- colnames(obj)[!is.na(group_vec)]

  if (length(keep_cells) == 0) {
    return(empty_de_result(
      comparison,
      comparison_label,
      cell_type,
      cells_medulla = 0L,
      cells_cortex = 0L,
      status = "skipped",
      reason = "No cells with compartment labels."
    ))
  }

  obj <- subset(obj, cells = keep_cells)
  group_vec <- obj$compartment
  n_tab <- table(group_vec)

  cells_medulla <- if ("Medulla" %in% names(n_tab)) as.integer(n_tab[["Medulla"]]) else 0L
  cells_cortex <- if ("Cortex" %in% names(n_tab)) as.integer(n_tab[["Cortex"]]) else 0L

  if (cells_medulla < min_cells_per_region || cells_cortex < min_cells_per_region) {
    return(empty_de_result(
      comparison,
      comparison_label,
      cell_type,
      cells_medulla = cells_medulla,
      cells_cortex = cells_cortex,
      status = "skipped",
      reason = paste0(
        "Fewer than ",
        min_cells_per_region,
        " cells in cortex or medulla."
      )
    ))
  }

  Idents(obj) <- factor(obj$compartment, levels = c("Cortex", "Medulla"))

  de <- FindMarkers(
    obj,
    ident.1 = "Medulla",
    ident.2 = "Cortex",
    assay = DefaultAssay(obj),
    slot = "data",
    test.use = "wilcox",
    min.pct = 0,
    logfc.threshold = 0,
    only.pos = FALSE,
    verbose = FALSE
  )

  if (is.null(de) || nrow(de) == 0) {
    return(empty_de_result(
      comparison,
      comparison_label,
      cell_type,
      cells_medulla = cells_medulla,
      cells_cortex = cells_cortex,
      status = "skipped",
      reason = "FindMarkers returned no genes."
    ))
  }

  de <- de |>
    rownames_to_column("gene") |>
    as_tibble()

  if (!"avg_log2FC" %in% colnames(de) && "avg_logFC" %in% colnames(de)) {
    de <- de |>
      rename(avg_log2FC = avg_logFC)
  }

  if (!"pct.1" %in% colnames(de)) {
    de$pct.1 <- NA_real_
  }

  if (!"pct.2" %in% colnames(de)) {
    de$pct.2 <- NA_real_
  }

  if (!"p_val" %in% colnames(de)) {
    de$p_val <- NA_real_
  }

  if (!"p_val_adj" %in% colnames(de)) {
    de$p_val_adj <- NA_real_
  }

  de <- de |>
    transmute(
      comparison = comparison,
      comparison_label = comparison_label,
      cell_type = cell_type,
      ident_1 = "Medulla",
      ident_2 = "Cortex",
      cells_medulla = cells_medulla,
      cells_cortex = cells_cortex,
      gene = gene,
      avg_log2FC_medulla_vs_cortex = avg_log2FC,
      pct_medulla = pct.1,
      pct_cortex = pct.2,
      p_val = p_val,
      p_val_adj = p_val_adj,
      enriched_region = case_when(
        avg_log2FC_medulla_vs_cortex > 0 ~ "Medulla",
        avg_log2FC_medulla_vs_cortex < 0 ~ "Cortex",
        TRUE ~ "No direction"
      ),
      significant_FDR_0_05 = is.finite(p_val_adj) & p_val_adj < fdr_cutoff
    ) |>
    arrange(
      p_val_adj,
      desc(abs(avg_log2FC_medulla_vs_cortex)),
      gene
    ) |>
    mutate(rank_within_complete_deg_table = row_number())

  out_file <- file.path(
    out_complete_dir,
    paste0("xenium_medulla_cortex_de_", safe_file_label(comparison), ".csv")
  )

  write_csv(de, out_file)
  summary <- tibble(
    comparison = comparison,
    comparison_label = comparison_label,
    cell_type = cell_type,
    cells_medulla = cells_medulla,
    cells_cortex = cells_cortex,
    genes_tested = nrow(de),
    significant_genes_FDR_0_05 = sum(de$significant_FDR_0_05, na.rm = TRUE),
    medulla_enriched_genes_FDR_0_05 = sum(
      de$significant_FDR_0_05 & de$enriched_region == "Medulla",
      na.rm = TRUE
    ),
    cortex_enriched_genes_FDR_0_05 = sum(
      de$significant_FDR_0_05 & de$enriched_region == "Cortex",
      na.rm = TRUE
    ),
    complete_deg_table = out_file,
    status = "completed",
    reason = NA_character_
  )

  list(de = de, summary = summary)
}

# -------------------------------------------------------------------------
# Inputs.
# -------------------------------------------------------------------------
stopifnot(file.exists(manifest_file))
stopifnot(file.exists(xenium_rds))
stopifnot(file.exists(snrna_annotation_metadata))
stopifnot(file.exists(medulla_roi_csv))

manifest <- read_csv(manifest_file, show_col_types = FALSE)

required_manifest_cols <- c("sample_id", "gestational_week", "xenium_dir", "include")
missing_manifest_cols <- setdiff(required_manifest_cols, colnames(manifest))

if (length(missing_manifest_cols) > 0) {
  stop("Manifest missing columns: ", paste(missing_manifest_cols, collapse = ", "))
}

manifest <- manifest |>
  mutate(
    include = normalize_include(include),
    sample_id = as.character(sample_id),
    gestational_week = as.character(gestational_week),
    xenium_dir = as.character(xenium_dir)
  ) |>
  filter(include) |>
  mutate(
    xenium_dir = vapply(
      xenium_dir,
      resolve_data_path,
      data_root = data_root,
      FUN.VALUE = character(1)
    )
  )

if (nrow(manifest) != 1) {
  stop("Expected exactly one included Xenium sample; found: ", nrow(manifest))
}

sample_id <- manifest$sample_id[[1]]
gestational_week <- manifest$gestational_week[[1]]
cells_csv <- file.path(manifest$xenium_dir[[1]], "cells.csv.gz")

stopifnot(file.exists(cells_csv))

cat("Loading Xenium object:\n  ", xenium_rds, "\n", sep = "")
xenium <- readRDS(xenium_rds)
stopifnot(inherits(xenium, "Seurat"))

if (length(xenium@images) > 0) {
  xenium@images <- list()
}

xen_assay <- if ("Xenium" %in% Assays(xenium)) "Xenium" else DefaultAssay(xenium)
DefaultAssay(xenium) <- xen_assay

if (xen_assay != "Xenium") {
  stop("Expected Xenium assay, found: ", xen_assay)
}

if (exists("JoinLayers", mode = "function")) {
  xenium <- tryCatch(
    JoinLayers(xenium, assay = xen_assay),
    error = function(e) xenium
  )
  DefaultAssay(xenium) <- xen_assay
}

metadata_cols <- colnames(xenium@meta.data)

fine_col <- first_existing_column(
  metadata_cols,
  c(
    "pred_fine_subcluster_hybrid_filtered",
    "pred_fine_subcluster_hybrid",
    "pred_fine_subcluster_2stage_filtered"
  ),
  "fine-subcluster"
)

annotation_col <- first_existing_column(
  metadata_cols,
  c(
    "pred_final_annotation_hybrid_filtered",
    "pred_final_annotation_hybrid",
    "pred_final_annotation_2stage_filtered"
  ),
  "final-annotation"
)

cat("Using fine-subcluster column: ", fine_col, "\n", sep = "")
cat("Using final-annotation column: ", annotation_col, "\n", sep = "")

cat("Loading snRNA-seq annotation metadata:\n  ", snrna_annotation_metadata, "\n", sep = "")
snrna_md <- read_csv(snrna_annotation_metadata, show_col_types = FALSE)

required_snrna_cols <- c("fine_subcluster", "major_cell_type", "final_annotation")
missing_snrna_cols <- setdiff(required_snrna_cols, colnames(snrna_md))

if (length(missing_snrna_cols) > 0) {
  stop("snRNA-seq metadata missing columns: ", paste(missing_snrna_cols, collapse = ", "))
}

subcluster_map <- snrna_md |>
  transmute(
    fine_subcluster = as.character(fine_subcluster),
    major_cell_type = recode(as.character(major_cell_type), quiescent = "degenerated"),
    final_annotation_reference = as.character(final_annotation)
  ) |>
  distinct(fine_subcluster, .keep_all = TRUE)

fine_to_celltype <- setNames(subcluster_map$major_cell_type, subcluster_map$fine_subcluster)
fine_to_annotation <- setNames(
  subcluster_map$final_annotation_reference,
  subcluster_map$fine_subcluster
)

cat("Loading Xenium cell coordinates:\n  ", cells_csv, "\n", sep = "")
cells_raw <- read_csv(cells_csv, show_col_types = FALSE)

cell_id_col <- pick_column(cells_raw, c("cell_id", "CellID", "cell", "barcode"))
x_col <- pick_column(
  cells_raw,
  c("x_centroid", "x_center", "x_location", "xglobal_px", "X", "x_centroid_px", "x_centroid_um")
)
y_col <- pick_column(
  cells_raw,
  c("y_centroid", "y_center", "y_location", "yglobal_px", "Y", "y_centroid_px", "y_centroid_um")
)

cells <- cells_raw |>
  transmute(
    cell_id = as.character(.data[[cell_id_col]]),
    x = suppressWarnings(as.numeric(.data[[x_col]])),
    y = suppressWarnings(as.numeric(.data[[y_col]]))
  ) |>
  filter(!is.na(cell_id), is.finite(x), is.finite(y))

cat("Loading medulla ROI coordinates:\n  ", medulla_roi_csv, "\n", sep = "")
medulla <- read_polygon_coordinates(medulla_roi_csv, "medulla ROI")

metadata <- xenium@meta.data |>
  rownames_to_column("cell_id") |>
  transmute(
    cell_id = cell_id,
    fine_subcluster = as.character(.data[[fine_col]]),
    final_annotation = as.character(.data[[annotation_col]])
  ) |>
  mutate(
    final_annotation_from_reference = unname(fine_to_annotation[fine_subcluster]),
    final_annotation = ifelse(
      is.na(final_annotation) | final_annotation == "" | final_annotation == fine_subcluster,
      final_annotation_from_reference,
      final_annotation
    ),
    cell_type = unname(fine_to_celltype[fine_subcluster]),
    cell_type = recode(as.character(cell_type), quiescent = "degenerated"),
    is_ambiguous = fine_subcluster %in% c("Ambiguous", "Ambiguous_subcluster") |
      final_annotation %in% c("Ambiguous", "Ambiguous_subcluster") |
      is.na(fine_subcluster) |
      is.na(final_annotation) |
      is.na(cell_type)
  ) |>
  filter(!is_ambiguous, !is.na(cell_type), cell_type %in% celltype_order) |>
  select(cell_id, fine_subcluster, final_annotation, cell_type)
cell_join_rate <- mean(metadata$cell_id %in% cells$cell_id)

if (!is.finite(cell_join_rate) || cell_join_rate < 0.90) {
  stop(
    "Cell coordinate join rate is too low: ",
    round(100 * cell_join_rate, 2),
    "%. Check that the Xenium object, metadata, and cells.csv.gz are from the same dataset."
  )
}

cat("Assigning medulla/cortex compartment labels.\n")

analysis_md <- metadata |>
  inner_join(cells, by = "cell_id") |>
  mutate(
    compartment = ifelse(
      point_in_polygon(x, y, medulla$x, medulla$y),
      "Medulla",
      "Cortex"
    ),
    compartment = factor(compartment, levels = c("Cortex", "Medulla")),
    cell_type = factor(cell_type, levels = celltype_order)
  )

assignment_summary <- tibble(
  sample_id = sample_id,
  gestational_week = gestational_week,
  total_nonambiguous_cells_with_coordinates = nrow(analysis_md),
  coordinate_join_rate = cell_join_rate,
  cells_cortex = sum(analysis_md$compartment == "Cortex"),
  cells_medulla = sum(analysis_md$compartment == "Medulla"),
  medulla_roi_vertices = nrow(medulla),
  min_cells_per_region = min_cells_per_region,
  top_genes_per_region = top_genes_per_region
)

write_csv(
  assignment_summary,
  file.path(out_table_dir, "xenium_medulla_cortex_de_assignment_summary.csv")
)

compartment_counts <- analysis_md |>
  count(cell_type, compartment, name = "cells") |>
  complete(
    cell_type = factor(celltype_order, levels = celltype_order),
    compartment = factor(c("Cortex", "Medulla"), levels = c("Cortex", "Medulla")),
    fill = list(cells = 0L)
  ) |>
  group_by(cell_type) |>
  mutate(
    cells_in_celltype = sum(cells),
    region_threshold_pass = all(cells >= min_cells_per_region),
    percent_within_celltype = ifelse(cells_in_celltype > 0, 100 * cells / cells_in_celltype, 0)
  ) |>
  ungroup() |>
  arrange(cell_type, compartment)

write_csv(
  compartment_counts,
  file.path(out_table_dir, "xenium_medulla_cortex_de_compartment_counts.csv")
)

common_cells <- intersect(colnames(xenium), analysis_md$cell_id)

analysis_md <- analysis_md |>
  filter(cell_id %in% common_cells) |>
  arrange(match(cell_id, common_cells))

xenium <- subset(xenium, cells = common_cells)

analysis_md <- analysis_md |>
  arrange(match(cell_id, colnames(xenium)))

if (!identical(analysis_md$cell_id, colnames(xenium))) {
  stop("Internal metadata ordering failure after subsetting Xenium object.")
}

xenium$compartment <- analysis_md$compartment
xenium$fine_subcluster <- analysis_md$fine_subcluster
xenium$final_annotation <- analysis_md$final_annotation
xenium$cell_type <- analysis_md$cell_type

cat("Normalizing Xenium expression data.\n")
xenium <- NormalizeData(
  xenium,
  assay = xen_assay,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

DefaultAssay(xenium) <- xen_assay

de_results <- list()

de_results[["all_cells"]] <- run_de(
  obj = xenium,
  comparison = "all_cells",
  comparison_label = "All annotated cells",
  cell_type = "all_cells"
)

celltypes_to_run <- celltype_order[
  celltype_order %in% as.character(unique(analysis_md$cell_type))
]

# Germ cells are retained in the all-cells medulla/cortex comparison but are
# excluded from the cell-type-specific DE comparisons.
celltypes_to_run <- setdiff(celltypes_to_run, "germ")

for (ct in celltypes_to_run) {
  cells_ct <- colnames(xenium)[xenium$cell_type == ct & !is.na(xenium$cell_type)]

  obj_ct <- subset(xenium, cells = cells_ct)

  de_results[[paste0("celltype_", ct)]] <- run_de(
    obj = obj_ct,
    comparison = paste0("celltype_", ct),
    comparison_label = paste("Cell type:", ct),
    cell_type = ct
  )

  rm(obj_ct)
  gc()
}

de_complete_combined <- bind_rows(lapply(de_results, `[[`, "de"))
run_summary <- bind_rows(lapply(de_results, `[[`, "summary"))

if (nrow(de_complete_combined) == 0) {
  warning("No completed DE comparisons produced genes.")
  de_top50 <- tibble()
} else {
  de_top50 <- de_complete_combined |>
    filter(
      significant_FDR_0_05,
      enriched_region %in% c("Medulla", "Cortex")
    ) |>
    arrange(
      comparison,
      enriched_region,
      p_val_adj,
      desc(abs(avg_log2FC_medulla_vs_cortex)),
      gene
    ) |>
    group_by(comparison, comparison_label, cell_type, enriched_region) |>
    mutate(rank_within_comparison_and_region = row_number()) |>
    slice_head(n = top_genes_per_region) |>
    ungroup()
}

write_csv(
  de_top50,
  file.path(out_table_dir, "xenium_medulla_cortex_de_top50_by_comparison_and_region.csv")
)

write_csv(
  run_summary,
  file.path(out_table_dir, "xenium_medulla_cortex_de_run_summary.csv")
)

summary_lines <- c(
  paste("Input Xenium object:", xenium_rds),
  paste("Input Xenium manifest:", manifest_file),
  paste("Input cells.csv.gz:", cells_csv),
  paste("Input snRNA-seq annotation metadata:", snrna_annotation_metadata),
  paste("Input medulla ROI:", medulla_roi_csv),
  paste("Sample ID:", sample_id),
  paste("Gestational week:", gestational_week),
  paste("Fine-subcluster column:", fine_col),
  paste("Final-annotation column:", annotation_col),
  paste("Coordinate join rate:", round(100 * cell_join_rate, 3), "%"),
  paste("Nonambiguous annotated cells with coordinates:", nrow(analysis_md)),
  paste("Cortex cells:", sum(analysis_md$compartment == "Cortex")),
  paste("Medulla cells:", sum(analysis_md$compartment == "Medulla")),
  paste("Minimum cells per region for each DE comparison:", min_cells_per_region),
  paste("Top genes per comparison and enriched region:", top_genes_per_region),
  "Germ cells are included in the all-cells comparison but excluded from cell-type-specific comparisons.",
  "",
  "DE design:",
  "FindMarkers ident.1 = Medulla, ident.2 = Cortex.",
  "Positive avg_log2FC_medulla_vs_cortex indicates medulla-enriched expression.",
  "Negative avg_log2FC_medulla_vs_cortex indicates cortex-enriched expression.",
  "",
  "Generated complete DEG tables:",
  paste(
    "  ",
    list.files(out_complete_dir, pattern = "\\.csv$", full.names = TRUE),
    collapse = "\n"
  ),
  "",
  "Generated summary tables:",
  paste("  ", file.path(out_table_dir, "xenium_medulla_cortex_de_assignment_summary.csv")),
  paste("  ", file.path(out_table_dir, "xenium_medulla_cortex_de_compartment_counts.csv")),
  paste("  ", file.path(out_table_dir, "xenium_medulla_cortex_de_top50_by_comparison_and_region.csv")),
  paste("  ", file.path(out_table_dir, "xenium_medulla_cortex_de_run_summary.csv")),
  "",
  "Assignment summary:",
  paste(capture.output(print(as.data.frame(assignment_summary))), collapse = "\n"),
  "",
  "DE run summary:",
  paste(capture.output(print(as.data.frame(run_summary))), collapse = "\n")
)

writeLines(
  summary_lines,
  file.path(out_log_dir, "xenium_medulla_cortex_differential_expression_summary.txt")
)
sink(file.path(out_log_dir, "xenium_medulla_cortex_differential_expression_sessionInfo.txt"))
print(sessionInfo())
sink()

cat("\nDone.\n")
cat("Xenium medulla/cortex differential-expression outputs written to:\n  ", results_root, "\n", sep = "")
