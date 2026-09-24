#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
})

set.seed(42)
options(bitmapType = "cairo")

# ============================================================
# Paths
# ============================================================

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

out_root <- file.path(results_base, "snRNAseq_DEGs_celltype_subcluster")

out_celltype_full_dir <- file.path(out_root, "tables", "cell_type_full")
out_celltype_top_dir <- file.path(out_root, "tables", "cell_type_top50")
out_subcluster_full_root <- file.path(out_root, "tables", "subcluster_full_by_cell_type")
out_subcluster_top_dir <- file.path(out_root, "tables", "subcluster_top50")
out_qc_dir <- file.path(out_root, "tables", "qc")
out_log_dir <- file.path(out_root, "logs")

dir.create(out_celltype_full_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_celltype_top_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_subcluster_full_root, recursive = TRUE, showWarnings = FALSE)
dir.create(out_subcluster_top_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_log_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Parameters
# ============================================================

assay_use <- "RNA"
test_use <- Sys.getenv("DEG_TEST_USE", unset = "wilcox")

full_logfc_threshold <- as.numeric(Sys.getenv("DEG_FULL_LOGFC_THRESHOLD", unset = "0"))
full_min_pct <- as.numeric(Sys.getenv("DEG_FULL_MIN_PCT", unset = "0"))

top_n <- as.integer(Sys.getenv("DEG_TOP_N", unset = "50"))
top_logfc_threshold <- as.numeric(Sys.getenv("DEG_TOP_LOGFC_THRESHOLD", unset = "0.25"))

# Priority marker rule for top tables:
# genes passing these stricter expression/score criteria are selected first,
# then remaining top-N slots are filled by the normal adaptive selection.
priority_pct1_threshold <- as.numeric(Sys.getenv("DEG_TOP_PRIORITY_PCT1_THRESHOLD", unset = "0.75"))
priority_score_threshold <- as.numeric(Sys.getenv("DEG_TOP_PRIORITY_SCORE_THRESHOLD", unset = "0.5"))

min_cells_per_group <- as.integer(Sys.getenv("DEG_MIN_CELLS_PER_GROUP", unset = "3"))

threshold_schedule <- tibble(
  pct1_threshold = c(0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10),
  padj_threshold = c(0.05, 0.10, 0.15, 0.25, 0.40, 0.60, 0.80, 0.90, 1.00)
)

write_csv(
  threshold_schedule,
  file.path(out_qc_dir, "DEG_top50_adaptive_threshold_schedule.csv")
)

# ============================================================
# Helpers
# ============================================================

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

is_mito_or_ribo_gene <- function(gene) {
  g <- toupper(gene)

  grepl("^MT-", g) |
    grepl("^RPL", g) |
    grepl("^RPS", g) |
    grepl("^MRPL", g) |
    grepl("^MRPS", g)
}

is_excluded_from_top_marker_table <- function(gene) {
  g <- toupper(gene)

  grepl("^AC[0-9]+(\\.[0-9]+)?$", g) |
    grepl("^AL[0-9]+(\\.[0-9]+)?$", g) |
    grepl("^LINC[0-9]+", g)
}

derive_cell_type <- function(meta) {
  if ("major_cell_type" %in% colnames(meta)) {
    ct <- as.character(meta$major_cell_type)
  } else if ("cell_type_quiescent" %in% colnames(meta)) {
    ct <- as.character(meta$cell_type_quiescent)
  } else if ("cell_type" %in% colnames(meta)) {
    ct <- as.character(meta$cell_type)
  } else if ("fine_subcluster" %in% colnames(meta)) {
    fs <- as.character(meta$fine_subcluster)
    ct <- sub("_.*$", "", fs)
    ct[fs %in% c("germ_0", "germ_6")] <- "degenerated"
  } else {
    stop("Could not derive cell type. Need major_cell_type, cell_type_quiescent, cell_type, or fine_subcluster.")
  }

  ct[ct == "quiescent"] <- "degenerated"
  ct
}

normalize_deg_table <- function(deg) {
  deg <- deg |>
    rownames_to_column("gene")

  if (!"avg_log2FC" %in% colnames(deg)) {
    if ("avg_logFC" %in% colnames(deg)) {
      deg <- deg |>
        rename(avg_log2FC = avg_logFC)
    } else {
      stop("DEG table does not contain avg_log2FC or avg_logFC.")
    }
  }

  required <- c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")
  missing_required <- setdiff(required, colnames(deg))

  if (length(missing_required) > 0) {
    stop("DEG table is missing required columns: ", paste(missing_required, collapse = ", "))
  }

  deg |>
    filter(!is_mito_or_ribo_gene(gene)) |>
    mutate(
      gene = as.character(gene),
      p_val = as.numeric(p_val),
      avg_log2FC = as.numeric(avg_log2FC),
      pct.1 = as.numeric(pct.1),
      pct.2 = as.numeric(pct.2),
      p_val_adj = as.numeric(p_val_adj),
      score = (pct.1 - pct.2) * avg_log2FC
    ) |>
    arrange(desc(score), desc(avg_log2FC), p_val_adj, desc(pct.1), gene)
}
select_top_adaptive <- function(deg, top_n, logfc_threshold, threshold_schedule) {
  selected <- NULL
  chosen_pct <- NA_real_
  chosen_padj <- NA_real_
  chosen_n_candidates <- 0L
  chosen_n_priority_candidates <- 0L
  chosen_n_priority_selected <- 0L
  chosen_n_fill_selected <- 0L
  reached_top_n <- FALSE

  for (i in seq_len(nrow(threshold_schedule))) {
    pct_cut <- threshold_schedule$pct1_threshold[[i]]
    padj_cut <- threshold_schedule$padj_threshold[[i]]

    candidates <- deg |>
      filter(
        !is_excluded_from_top_marker_table(gene),
        avg_log2FC > logfc_threshold,
        pct.1 > pct_cut,
        p_val_adj < padj_cut
      ) |>
      arrange(desc(score), desc(avg_log2FC), p_val_adj, desc(pct.1), gene)

    priority_candidates <- candidates |>
      filter(
        pct.1 > priority_pct1_threshold,
        score > priority_score_threshold
      ) |>
      arrange(desc(score), desc(avg_log2FC), p_val_adj, desc(pct.1), gene)

    priority_selected <- priority_candidates |>
      slice_head(n = top_n)

    n_remaining <- max(top_n - nrow(priority_selected), 0L)

    fill_selected <- candidates |>
      filter(!gene %in% priority_selected$gene) |>
      arrange(desc(score), desc(avg_log2FC), p_val_adj, desc(pct.1), gene) |>
      slice_head(n = n_remaining)

    selected_this_round <- bind_rows(priority_selected, fill_selected) |>
      mutate(
        top50_selection_class = ifelse(
          gene %in% priority_selected$gene,
          "priority_pct1_gt_0.75_score_gt_0.5",
          "standard_adaptive_fill"
        )
      ) |>
      arrange(
        factor(
          top50_selection_class,
          levels = c("priority_pct1_gt_0.75_score_gt_0.5", "standard_adaptive_fill")
        ),
        desc(score),
        desc(avg_log2FC),
        p_val_adj,
        desc(pct.1),
        gene
      )

    if (nrow(selected_this_round) >= top_n || i == nrow(threshold_schedule)) {
      selected <- selected_this_round |>
        slice_head(n = top_n)

      chosen_pct <- pct_cut
      chosen_padj <- padj_cut
      chosen_n_candidates <- nrow(candidates)
      chosen_n_priority_candidates <- nrow(priority_candidates)
      chosen_n_priority_selected <- sum(selected$top50_selection_class == "priority_pct1_gt_0.75_score_gt_0.5")
      chosen_n_fill_selected <- sum(selected$top50_selection_class == "standard_adaptive_fill")
      reached_top_n <- nrow(selected_this_round) >= top_n
      break
    }
  }

  selected |>
    mutate(
      adaptive_pct1_threshold = chosen_pct,
      adaptive_padj_threshold = chosen_padj,
      adaptive_n_candidates = chosen_n_candidates,
      adaptive_n_priority_candidates = chosen_n_priority_candidates,
      adaptive_n_priority_selected = chosen_n_priority_selected,
      adaptive_n_standard_fill_selected = chosen_n_fill_selected,
      adaptive_reached_top_n = reached_top_n,
      top_priority_pct1_threshold = priority_pct1_threshold,
      top_priority_score_threshold = priority_score_threshold,
      top_logfc_threshold = logfc_threshold,
      top_n_requested = top_n
    )
}

run_findmarkers <- function(obj, ident_1, ident_2 = NULL, features_use) {
  FindMarkers(
    object = obj,
    ident.1 = ident_1,
    ident.2 = ident_2,
    assay = assay_use,
    slot = "data",
    features = features_use,
    logfc.threshold = full_logfc_threshold,
    min.pct = full_min_pct,
    test.use = test_use,
    only.pos = FALSE,
    verbose = FALSE
  )
}

# ============================================================
# Load object
# ============================================================

cat("Loading annotated snRNA-seq object:\n  ", annotated_rds, "\n", sep = "")
stopifnot(file.exists(annotated_rds))

obj <- readRDS(annotated_rds)
stopifnot(inherits(obj, "Seurat"))

DefaultAssay(obj) <- assay_use
stopifnot(assay_use %in% names(obj@assays))

if (!"fine_subcluster" %in% colnames(obj@meta.data)) {
  stop("Input object must contain fine_subcluster metadata.")
}

if (!"final_annotation" %in% colnames(obj@meta.data)) {
  warning("final_annotation metadata not found; subcluster display labels will use fine_subcluster.")
}

# Join layers if needed for Seurat v5 DE.
if (inherits(obj[[assay_use]], "Assay5")) {
  layers_before <- Layers(obj[[assay_use]])

  if (any(grepl("^data\\.", layers_before)) || any(grepl("^counts\\.", layers_before))) {
    cat("Joining RNA layers before DEG analysis.\n")
    obj <- tryCatch(
      JoinLayers(obj, assay = assay_use),
      error = function(e) {
        warning("JoinLayers failed; continuing with existing layers. Error: ", e$message)
        obj
      }
    )
  }
}

if (inherits(obj[[assay_use]], "Assay5")) {
  layers_after <- Layers(obj[[assay_use]])
  if (!any(grepl("^data", layers_after))) {
    cat("No normalized RNA data layer detected. Running NormalizeData().\n")
    obj <- NormalizeData(obj, assay = assay_use, verbose = FALSE)
  }
}
meta <- obj@meta.data |>
  rownames_to_column("cell_id")

meta$fine_subcluster <- as.character(meta$fine_subcluster)

if ("final_annotation" %in% colnames(meta)) {
  meta$final_annotation <- as.character(meta$final_annotation)
} else {
  meta$final_annotation <- as.character(meta$fine_subcluster)
}

meta$cell_type_deg <- derive_cell_type(meta)

bad_genes <- rownames(obj[[assay_use]])[is_mito_or_ribo_gene(rownames(obj[[assay_use]]))]
features_use <- setdiff(rownames(obj[[assay_use]]), bad_genes)

if (length(features_use) == 0) {
  stop("No genes remain after mitochondrial/ribosomal gene filtering.")
}

write_csv(
  tibble(gene = bad_genes, excluded_reason = "mitochondrial_or_ribosomal"),
  file.path(out_qc_dir, "DEG_excluded_mitochondrial_ribosomal_genes.csv")
)

write_csv(
  meta |>
    count(cell_type_deg, name = "n_cells") |>
    arrange(cell_type_deg),
  file.path(out_qc_dir, "DEG_QC_cell_counts_by_cell_type.csv")
)

subcluster_label_map <- meta |>
  group_by(cell_type_deg, fine_subcluster) |>
  summarise(
    final_annotation = mode_string(final_annotation),
    n_cells = n(),
    .groups = "drop"
  ) |>
  arrange(cell_type_deg, fine_subcluster)

write_csv(
  subcluster_label_map,
  file.path(out_qc_dir, "DEG_QC_subcluster_label_map.csv")
)

write_csv(
  subcluster_label_map |>
    arrange(cell_type_deg, fine_subcluster),
  file.path(out_qc_dir, "DEG_QC_cell_counts_by_subcluster.csv")
)

# Store metadata back into object.
obj$cell_type_deg <- meta$cell_type_deg[match(colnames(obj), meta$cell_id)]
obj$fine_subcluster_deg <- meta$fine_subcluster[match(colnames(obj), meta$cell_id)]

# ============================================================
# Cell-type DEGs: one cell type vs all other cells
# ============================================================

cell_types <- sort(unique(obj$cell_type_deg))

if (exists("celltype_order", inherits = TRUE)) {
  cell_types <- c(
    intersect(celltype_order, cell_types),
    setdiff(cell_types, celltype_order)
  )
}

cat("\nCell types for DEG analysis:\n")
print(cell_types)

Idents(obj) <- obj$cell_type_deg

celltype_full_list <- list()
celltype_top_list <- list()
celltype_summary <- list()

for (ct in cell_types) {
  cat("\n============================================================\n")
  cat("Cell-type DEG: ", ct, " vs all other cells\n", sep = "")

  n_ident <- sum(obj$cell_type_deg == ct)
  n_other <- sum(obj$cell_type_deg != ct)

  if (n_ident < min_cells_per_group || n_other < min_cells_per_group) {
    warning("Skipping cell type ", ct, ": insufficient cells.")
    celltype_summary[[ct]] <- tibble(
      comparison_level = "cell_type",
      group = ct,
      cell_type = ct,
      n_ident = n_ident,
      n_reference = n_other,
      n_full_deg_rows = NA_integer_,
      n_top_rows = NA_integer_,
      status = "skipped_insufficient_cells"
    )
    next
  }

  deg <- run_findmarkers(
    obj = obj,
    ident_1 = ct,
    ident_2 = NULL,
    features_use = features_use
  ) |>
    normalize_deg_table() |>
    mutate(
      comparison_level = "cell_type",
      comparison = paste0(ct, "_vs_all_other_cells"),
      cell_type = ct,
      group = ct,
      n_ident = n_ident,
      n_reference = n_other
    ) |>
    relocate(
      comparison_level,
      comparison,
      cell_type,
      group,
      n_ident,
      n_reference,
      gene
    )

  full_file <- file.path(
    out_celltype_full_dir,
    paste0("DEG_cell_type__", safe_file_component(ct), "__vs_all_other_cells.csv")
  )
  write_csv(deg, full_file)

  top <- select_top_adaptive(
    deg = deg,
    top_n = top_n,
    logfc_threshold = top_logfc_threshold,
    threshold_schedule = threshold_schedule
  ) |>
    mutate(full_deg_file = full_file) |>
    relocate(full_deg_file, .after = n_reference)

  celltype_full_list[[ct]] <- deg
  celltype_top_list[[ct]] <- top

  celltype_summary[[ct]] <- tibble(
    comparison_level = "cell_type",
    group = ct,
    cell_type = ct,
    n_ident = n_ident,
    n_reference = n_other,
    n_full_deg_rows = nrow(deg),
    n_top_rows = nrow(top),
    adaptive_pct1_threshold = unique(top$adaptive_pct1_threshold),
    adaptive_padj_threshold = unique(top$adaptive_padj_threshold),
    adaptive_n_candidates = unique(top$adaptive_n_candidates),
    adaptive_reached_top_n = unique(top$adaptive_reached_top_n),
    status = "ok"
  )
}

celltype_full_combined <- bind_rows(celltype_full_list)
celltype_top_combined <- bind_rows(celltype_top_list)
celltype_summary_tbl <- bind_rows(celltype_summary)

write_csv(
  celltype_full_combined,
  file.path(out_celltype_full_dir, "DEG_cell_type__ALL_cell_types_vs_all_other_cells_combined.csv")
)

write_csv(
  celltype_top_combined,
  file.path(out_celltype_top_dir, "DEG_cell_type__top50_adaptive_thresholds_combined.csv")
)

write_csv(
  celltype_summary_tbl,
  file.path(out_qc_dir, "DEG_QC_cell_type_summary.csv")
)

# ============================================================
# Subcluster DEGs: one subcluster vs other subclusters within same cell type
# ============================================================

subcluster_top_list <- list()
subcluster_summary <- list()

for (ct in cell_types) {
  cat("\n============================================================\n")
  cat("Subcluster DEGs within cell type: ", ct, "\n", sep = "")

  ct_cells <- colnames(obj)[obj$cell_type_deg == ct]
  subclusters_ct <- sort(unique(obj$fine_subcluster_deg[ct_cells]))

  if (length(subclusters_ct) < 2) {
    cat("Skipping ", ct, ": fewer than two subclusters.\n", sep = "")

    subcluster_summary[[paste0(ct, "__skipped")]] <- tibble(
      comparison_level = "subcluster_within_cell_type",
      cell_type = ct,
      fine_subcluster = NA_character_,
      final_annotation = NA_character_,
      n_ident = NA_integer_,
      n_reference = NA_integer_,
      n_full_deg_rows = NA_integer_,
      n_top_rows = NA_integer_,
      status = "skipped_fewer_than_two_subclusters"
    )
    next
  }

  ct_dir <- file.path(out_subcluster_full_root, safe_file_component(ct))
  dir.create(ct_dir, recursive = TRUE, showWarnings = FALSE)

  obj_ct <- subset(obj, cells = ct_cells)
  Idents(obj_ct) <- obj_ct$fine_subcluster_deg

  for (sc in subclusters_ct) {
    cat("Subcluster DEG: ", sc, " within ", ct, "\n", sep = "")

    n_ident <- sum(obj_ct$fine_subcluster_deg == sc)
    n_reference <- sum(obj_ct$fine_subcluster_deg != sc)

    label_sc <- subcluster_label_map |>
      filter(cell_type_deg == ct, fine_subcluster == sc) |>
      pull(final_annotation)

    if (length(label_sc) == 0) {
      label_sc <- sc
    }

    if (n_ident < min_cells_per_group || n_reference < min_cells_per_group) {
      warning("Skipping subcluster ", sc, " within ", ct, ": insufficient cells.")

      subcluster_summary[[paste0(ct, "__", sc)]] <- tibble(
        comparison_level = "subcluster_within_cell_type",
        cell_type = ct,
        fine_subcluster = sc,
        final_annotation = label_sc[[1]],
        n_ident = n_ident,
        n_reference = n_reference,
        n_full_deg_rows = NA_integer_,
        n_top_rows = NA_integer_,
        status = "skipped_insufficient_cells"
      )
      next
    }
    other_subclusters <- setdiff(subclusters_ct, sc)

    deg <- run_findmarkers(
      obj = obj_ct,
      ident_1 = sc,
      ident_2 = other_subclusters,
      features_use = features_use
    ) |>
      normalize_deg_table() |>
      mutate(
        comparison_level = "subcluster_within_cell_type",
        comparison = paste0(sc, "_vs_other_", ct, "_subclusters"),
        cell_type = ct,
        fine_subcluster = sc,
        final_annotation = label_sc[[1]],
        n_ident = n_ident,
        n_reference = n_reference
      ) |>
      relocate(
        comparison_level,
        comparison,
        cell_type,
        fine_subcluster,
        final_annotation,
        n_ident,
        n_reference,
        gene
      )

    full_file <- file.path(
      ct_dir,
      paste0(
        "DEG_subcluster__",
        safe_file_component(sc),
        "__within_",
        safe_file_component(ct),
        ".csv"
      )
    )

    write_csv(deg, full_file)

    top <- select_top_adaptive(
      deg = deg,
      top_n = top_n,
      logfc_threshold = top_logfc_threshold,
      threshold_schedule = threshold_schedule
    ) |>
      mutate(full_deg_file = full_file) |>
      relocate(full_deg_file, .after = n_reference)

    subcluster_top_list[[paste0(ct, "__", sc)]] <- top

    subcluster_summary[[paste0(ct, "__", sc)]] <- tibble(
      comparison_level = "subcluster_within_cell_type",
      cell_type = ct,
      fine_subcluster = sc,
      final_annotation = label_sc[[1]],
      n_ident = n_ident,
      n_reference = n_reference,
      n_full_deg_rows = nrow(deg),
      n_top_rows = nrow(top),
      adaptive_pct1_threshold = unique(top$adaptive_pct1_threshold),
      adaptive_padj_threshold = unique(top$adaptive_padj_threshold),
      adaptive_n_candidates = unique(top$adaptive_n_candidates),
      adaptive_reached_top_n = unique(top$adaptive_reached_top_n),
      full_deg_file = full_file,
      status = "ok"
    )
  }
}

subcluster_top_combined <- bind_rows(subcluster_top_list)
subcluster_summary_tbl <- bind_rows(subcluster_summary)

write_csv(
  subcluster_top_combined,
  file.path(out_subcluster_top_dir, "DEG_subcluster__top50_adaptive_thresholds_combined.csv")
)

write_csv(
  subcluster_summary_tbl,
  file.path(out_qc_dir, "DEG_QC_subcluster_summary.csv")
)
# ============================================================
# Notes and session info
# ============================================================

notes <- c(
  paste("Input annotated snRNA-seq object:", annotated_rds),
  paste("Output root:", out_root),
  "This script generates publication-ready DEG tables from the GitHub-generated annotated fetal ovary snRNA-seq object.",
  "Cell-type DEGs compare each broad cell type against all other cells.",
  "Subcluster DEGs compare each fine subcluster against all other subclusters within the same broad cell type.",
  "Cell types with fewer than two subclusters are skipped for subcluster-level DEG testing.",
  "Mitochondrial and ribosomal genes are excluded before DEG testing.",
  paste("Excluded gene patterns:", "^MT-, ^RPL, ^RPS, ^MRPL, ^MRPS"),
  "AC*, AL*, and LINC* genes are retained in full DEG tables if tested, but excluded from adaptive top50 marker tables.",
  "All DEG tables include score = (pct.1 - pct.2) * avg_log2FC and are ordered by descending score.",
  paste("DEG test:", test_use),
  paste("Full-table logfc.threshold:", full_logfc_threshold),
  paste("Full-table min.pct:", full_min_pct),
  paste("Minimum cells per group:", min_cells_per_group),
  paste("Top N requested:", top_n),
  paste("Top-table logFC threshold:", top_logfc_threshold),
  "Top-table adaptive filter begins at avg_log2FC > 0.25, pct.1 > 0.5, and p_val_adj < 0.05.",
  "Within each adaptive threshold round, genes with pct.1 > 0.75 and score > 0.5 are selected first, then remaining top-50 slots are filled by the normal score-ranked adaptive selection.",
  paste("Top-table priority pct.1 threshold:", priority_pct1_threshold),
  paste("Top-table priority score threshold:", priority_score_threshold),
  "If fewer than 50 genes pass, pct.1 is decreased and p_val_adj is increased according to tables/qc/DEG_top50_adaptive_threshold_schedule.csv.",
  "The final fallback threshold is pct.1 > 0.1 and p_val_adj < 1.",
  "",
  "Main outputs:",
  paste("Cell-type full combined table:", file.path(out_celltype_full_dir, "DEG_cell_type__ALL_cell_types_vs_all_other_cells_combined.csv")),
  paste("Cell-type top50 combined table:", file.path(out_celltype_top_dir, "DEG_cell_type__top50_adaptive_thresholds_combined.csv")),
  paste("Subcluster full tables root:", out_subcluster_full_root),
  paste("Subcluster top50 combined table:", file.path(out_subcluster_top_dir, "DEG_subcluster__top50_adaptive_thresholds_combined.csv"))
)

writeLines(
  notes,
  file.path(out_log_dir, "DEG_celltype_subcluster_notes.txt")
)

sink(file.path(out_log_dir, "DEG_celltype_subcluster_sessionInfo.txt"))
print(sessionInfo())
sink()

cat("\nDone.\n")
cat("DEG outputs written to:\n  ", out_root, "\n", sep = "")
cat("Cell-type top50 table:\n  ", file.path(out_celltype_top_dir, "DEG_cell_type__top50_adaptive_thresholds_combined.csv"), "\n", sep = "")
cat("Subcluster top50 table:\n  ", file.path(out_subcluster_top_dir, "DEG_subcluster__top50_adaptive_thresholds_combined.csv"), "\n", sep = "")
