#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(monocle3)
  library(readr)
  library(tibble)
  library(igraph)
})

analysis_seed <- 42L
set.seed(analysis_seed)
options(bitmapType = "cairo")

command_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", command_args, value = TRUE)

if (length(file_arg) != 1L) {
  stop("Run this script with Rscript.")
}

script_path <- normalizePath(
  sub("^--file=", "", file_arg[[1]]),
  winslash = "/",
  mustWork = TRUE
)

project_root <- normalizePath(
  file.path(dirname(script_path), "..", ".."),
  winslash = "/",
  mustWork = TRUE
)

default_results_base <- file.path(
  dirname(project_root),
  paste0(basename(project_root), "_results")
)

results_base <- normalizePath(
  Sys.getenv("FETAL_OVARY_RESULTS_ROOT", unset = default_results_base),
  winslash = "/",
  mustWork = FALSE
)

source(file.path(project_root, "config", "plotting.R"))
source(file.path(project_root, "config", "labels_colors.R"))

required_plot_config <- c(
  "publication_font_family",
  "publication_base_size",
  "theme_publication",
  "save_publication_plot"
)

missing_plot_config <- required_plot_config[
  !vapply(required_plot_config, exists, logical(1), inherits = TRUE)
]

if (length(missing_plot_config) > 0L) {
  stop("Missing plotting configuration: ", paste(missing_plot_config, collapse = ", "))
}
annotated_rds <- file.path(
  results_base,
  "snRNAseq_annotated_object",
  "objects",
  "fetal_ovary_snRNAseq_canonical_annotated.rds"
)

results_root <- file.path(results_base, "snRNAseq_monocle3_pseudotime")
out_figure_dir <- file.path(results_root, "figures")
out_table_dir <- file.path(results_root, "tables")
out_object_dir <- file.path(results_root, "objects")
out_log_dir <- file.path(results_root, "logs")

for (path in c(out_figure_dir, out_table_dir, out_object_dir, out_log_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

assay_use <- "RNA"
germ_all_clusters <- paste0("germ_", 0:8)
germ_exclude_from_pseudotime <- c("germ_0", "germ_6")
germ_root_cluster <- "germ_3"

analysis_plan <- tribble(
  ~celltype,      ~workflow,          ~resolution, ~expected_clusters, ~pca_npcs, ~show_trajectory_graph,
  "germ",        "rna",                    0.14,                  9L,      50L, TRUE,
  "granulosa",   "sct_split_layers",       0.20,                  9L,      30L, FALSE,
  "stroma",      "sct_split_layers",       0.08,                  5L,      30L, FALSE,
  "endothelial", "sct_split_layers",       0.001,                 1L,      30L, FALSE,
  "mural",       "sct_split_layers",       0.40,                  2L,      30L, FALSE,
  "immune",      "sct_split_layers",       0.26,                  2L,      30L, FALSE,
  "erythroid",   "sct_split_layers",       0.30,                  2L,      15L, FALSE
)

if (!exists("fine_subcluster_colors", inherits = TRUE)) {
  fine_subcluster_colors <- c(
    germ_0 = "#6B6B6B",
    germ_1 = "#F8766D",
    germ_2 = "#C49A00",
    germ_3 = "#53B400",
    germ_4 = "#00C094",
    germ_5 = "#00B6EB",
    germ_6 = "#B0B0B0",
    germ_7 = "#A58AFF",
    germ_8 = "#FB61D7",
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
}

sanitize_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

display_label <- function(x) {
  x <- as.character(x)

  if (exists("subcluster_label_map", inherits = TRUE)) {
    label_map <- get("subcluster_label_map", inherits = TRUE)
    matched <- x %in% names(label_map)
    x[matched] <- unname(label_map[x[matched]])
  }

  x
}
sample_column <- function(obj) {
  candidates <- c("sampleID", "sample_id", "orig.ident")
  hit <- candidates[candidates %in% colnames(obj@meta.data)]

  if (length(hit) == 0L) {
    return(NA_character_)
  }

  hit[[1]]
}

get_counts_layer <- function(obj, assay = "RNA", layer = "counts") {
  if (!assay %in% names(obj@assays)) {
    stop("Assay not found: ", assay)
  }

  matrix <- LayerData(obj, assay = assay, layer = layer)

  if (!inherits(matrix, "dgCMatrix")) {
    matrix <- as(matrix, "dgCMatrix")
  }

  matrix
}

prepare_metadata <- function(metadata, sample_col) {
  if (!sample_col %in% colnames(metadata)) {
    stop("Sample column not found: ", sample_col)
  }

  week_source <- if ("gestational_week" %in% colnames(metadata)) {
    metadata$gestational_week
  } else if ("Week" %in% colnames(metadata)) {
    metadata$Week
  } else {
    stop("No gestational week column found.")
  }

  week_numeric <- suppressWarnings(readr::parse_number(as.character(week_source)))

  if (all(is.na(week_numeric))) {
    stop("Gestational week values could not be parsed.")
  }

  week_levels_numeric <- sort(unique(week_numeric[!is.na(week_numeric)]))

  format_week <- function(x) {
    formatC(x, format = "fg", digits = 4, drop0trailing = TRUE)
  }

  metadata |>
    mutate(
      major_cell_type = recode(as.character(major_cell_type), quiescent = "degenerated"),
      fine_subcluster = as.character(fine_subcluster),
      final_annotation = as.character(final_annotation),
      sample_for_split = as.character(.data[[sample_col]]),
      week_numeric = week_numeric,
      week_for_plot = factor(
        format_week(week_numeric),
        levels = format_week(week_levels_numeric)
      ),
      cluster_for_pseudotime = fine_subcluster,
      cluster_display = display_label(fine_subcluster)
    )
}

make_germ_rna_object <- function(counts, metadata) {
  cells <- rownames(metadata)[metadata$fine_subcluster %in% germ_all_clusters]

  if (length(cells) == 0L) {
    stop("No germ cells found.")
  }

  CreateSeuratObject(
    counts = counts[, cells, drop = FALSE],
    meta.data = metadata[cells, , drop = FALSE],
    assay = assay_use,
    project = "fetal_germ_rna"
  )
}
make_sample_split_object <- function(counts, metadata) {
  obj <- CreateSeuratObject(
    counts = counts,
    meta.data = metadata[colnames(counts), , drop = FALSE],
    assay = assay_use,
    project = "fetal_sample_split"
  )

  if (anyNA(obj$sample_for_split) || any(obj$sample_for_split == "")) {
    stop("Sample identifiers contain missing or empty values.")
  }

  DefaultAssay(obj) <- assay_use

  obj[[assay_use]] <- split(
    obj[[assay_use]],
    f = factor(obj$sample_for_split)
  )

  obj
}

make_non_germ_sct_subset <- function(obj, celltype) {
  cells <- rownames(obj@meta.data)[obj$major_cell_type == celltype]

  if (length(cells) == 0L) {
    stop("No cells found for cell type: ", celltype)
  }

  subset(obj, cells = cells)
}

run_germ_rna_workflow <- function(obj, resolution, pca_npcs) {
  DefaultAssay(obj) <- assay_use

  obj <- NormalizeData(obj, verbose = FALSE)

  obj <- FindVariableFeatures(
    obj,
    selection.method = "vst",
    nfeatures = 3000,
    verbose = FALSE
  )

  obj <- ScaleData(
    obj,
    features = VariableFeatures(obj),
    verbose = FALSE
  )

  obj <- RunPCA(
    obj,
    features = VariableFeatures(obj),
    npcs = pca_npcs,
    verbose = FALSE,
    seed.use = analysis_seed
  )

  dims_use <- seq_len(min(30L, ncol(Embeddings(obj, "pca"))))

  obj <- FindNeighbors(
    obj,
    reduction = "pca",
    dims = dims_use,
    verbose = FALSE
  )

  obj <- FindClusters(
    obj,
    resolution = resolution,
    algorithm = 1,
    random.seed = analysis_seed,
    verbose = FALSE
  )

  obj <- RunUMAP(
    obj,
    reduction = "pca",
    dims = dims_use,
    umap.method = "uwot",
    metric = "cosine",
    seed.use = analysis_seed,
    verbose = FALSE
  )

  obj$recluster_resolution <- resolution
  obj$recluster_workflow <- "RNA_NormalizeData_FindVariableFeatures_ScaleData_RunPCA_FindNeighbors_FindClusters_RunUMAP"
  obj$n_pcs_used <- length(dims_use)

  obj
}
run_non_germ_sct_workflow <- function(obj, resolution, pca_npcs) {
  DefaultAssay(obj) <- assay_use

  message("  RNA layers entering SCTransform: ", paste(Layers(obj[[assay_use]]), collapse = ", "))

  obj <- SCTransform(obj, verbose = FALSE)

  obj <- RunPCA(
    obj,
    npcs = pca_npcs,
    verbose = FALSE,
    seed.use = analysis_seed
  )

  dims_use <- seq_len(min(30L, ncol(Embeddings(obj, "pca"))))

  obj <- RunUMAP(
    obj,
    reduction = "pca",
    dims = dims_use,
    umap.method = "uwot",
    metric = "cosine",
    seed.use = analysis_seed,
    verbose = FALSE
  )

  obj <- FindNeighbors(
    obj,
    reduction = "pca",
    dims = dims_use,
    verbose = FALSE
  )

  obj <- FindClusters(
    obj,
    resolution = resolution,
    algorithm = 1,
    random.seed = analysis_seed,
    verbose = FALSE
  )

  obj$recluster_resolution <- resolution
  obj$recluster_workflow <- "SampleSplitRNA_SCTransform_RunPCA_RunUMAP_FindNeighbors_FindClusters"
  obj$n_pcs_used <- length(dims_use)

  obj
}

validate_cluster_count <- function(obj, celltype, expected_clusters) {
  observed_clusters <- length(unique(as.character(obj$seurat_clusters)))

  if (observed_clusters != expected_clusters) {
    stop("[", celltype, "] expected ", expected_clusters, " Seurat clusters but observed ", observed_clusters, ".")
  }

  observed_clusters
}

filter_pseudotime_cells <- function(obj, celltype) {
  if (celltype != "germ") {
    return(obj)
  }

  keep <- colnames(obj)[
    !as.character(obj$fine_subcluster) %in% germ_exclude_from_pseudotime
  ]

  subset(obj, cells = keep)
}

choose_root_cluster <- function(obj, celltype) {
  if (celltype == "germ") {
    return(germ_root_cluster)
  }

  root_table <- obj@meta.data |>
    transmute(
      cluster = as.character(cluster_for_pseudotime),
      gestational_week = as.numeric(week_numeric)
    ) |>
    filter(
      !is.na(cluster),
      cluster != "",
      is.finite(gestational_week)
    ) |>
    group_by(cluster) |>
    summarise(
      median_week = median(gestational_week),
      n_cells = n(),
      .groups = "drop"
    ) |>
    arrange(median_week, desc(n_cells), cluster)

  if (nrow(root_table) == 0L) {
    stop("[", celltype, "] could not determine root cluster.")
  }

  root_table$cluster[[1]]
}
seurat_to_cds <- function(obj, counts_source) {
  counts <- counts_source[, colnames(obj), drop = FALSE]

  cell_metadata <- obj@meta.data[colnames(obj), , drop = FALSE]

  gene_metadata <- data.frame(
    gene_short_name = rownames(counts),
    row.names = rownames(counts),
    stringsAsFactors = FALSE
  )

  new_cell_data_set(
    expression_data = counts,
    cell_metadata = cell_metadata,
    gene_metadata = gene_metadata
  )
}

run_monocle_on_existing_umap <- function(obj, counts_source, celltype, root_cluster) {
  cds <- seurat_to_cds(obj, counts_source)

  num_dim <- min(50L, ncol(obj) - 1L, nrow(counts_source) - 1L)

  if (num_dim < 2L) {
    stop("[", celltype, "] too few cells or genes for Monocle3.")
  }

  cds <- suppressMessages(
    preprocess_cds(
      cds,
      num_dim = num_dim
    )
  )

  existing_umap <- Embeddings(obj, "umap")
  existing_umap <- existing_umap[colnames(cds), , drop = FALSE]
  reducedDims(cds)$UMAP <- existing_umap

  cds <- cluster_cells(cds, reduction_method = "UMAP")
  cds <- learn_graph(cds, use_partition = FALSE)

  cluster_values <- as.character(colData(cds)$cluster_for_pseudotime)
  names(cluster_values) <- colnames(cds)

  if (!root_cluster %in% unique(cluster_values)) {
    stop("[", celltype, "] root cluster is absent: ", root_cluster)
  }

  root_cells <- names(cluster_values)[cluster_values == root_cluster]

  if (length(root_cells) == 0L) {
    stop("[", celltype, "] no root cells found.")
  }

  cds <- order_cells(cds, root_cells = root_cells)

  pseudotime_values <- as.numeric(monocle3::pseudotime(cds))
  names(pseudotime_values) <- colnames(cds)
  colData(cds)$pseudotime <- pseudotime_values

  list(
    cds = cds,
    pseudotime = pseudotime_values,
    root_cluster = root_cluster,
    root_cells = root_cells
  )
}

principal_graph_segments <- function(cds) {
  graph <- principal_graph(cds)[["UMAP"]]
  aux <- cds@principal_graph_aux[["UMAP"]]

  if (is.null(graph) || is.null(aux$dp_mst)) {
    return(tibble())
  }

  coords <- as.matrix(aux$dp_mst)

  if (nrow(coords) == 2L) {
    coords <- t(coords)
  }

  coords <- as.data.frame(coords[, 1:2, drop = FALSE])
  colnames(coords) <- c("UMAP_1", "UMAP_2")

  if (length(igraph::V(graph)$name) == nrow(coords)) {
    rownames(coords) <- igraph::V(graph)$name
  }

  edges <- igraph::ends(graph, igraph::E(graph), names = TRUE)
  tibble(
    from = edges[, 1],
    to = edges[, 2],
    x = coords[from, "UMAP_1"],
    y = coords[from, "UMAP_2"],
    xend = coords[to, "UMAP_1"],
    yend = coords[to, "UMAP_2"]
  ) |>
    filter(
      is.finite(x),
      is.finite(y),
      is.finite(xend),
      is.finite(yend)
    )
}

make_plot_data <- function(obj, pseudotime, celltype, root_cluster) {
  umap <- Embeddings(obj, "umap") |>
    as.data.frame() |>
    rownames_to_column("cell_id")

  colnames(umap)[2:3] <- c("UMAP_1", "UMAP_2")

  obj@meta.data |>
    rownames_to_column("cell_id") |>
    transmute(
      cell_id,
      celltype = .env$celltype,
      major_cell_type = as.character(major_cell_type),
      fine_subcluster = as.character(fine_subcluster),
      final_annotation = as.character(final_annotation),
      cluster_for_pseudotime = as.character(cluster_for_pseudotime),
      cluster_display = as.character(cluster_display),
      seurat_cluster = as.character(seurat_clusters),
      gestational_week = as.numeric(week_numeric),
      sample_id = as.character(sample_for_split),
      root_cluster = .env$root_cluster,
      pseudotime = pseudotime[cell_id]
    ) |>
    left_join(umap, by = "cell_id")
}

resolve_subcluster_colors <- function(plot_data) {
  key <- plot_data |>
    distinct(cluster_for_pseudotime, cluster_display) |>
    arrange(cluster_for_pseudotime)

  colors <- fine_subcluster_colors[key$cluster_for_pseudotime]
  missing <- is.na(colors)

  if (any(missing)) {
    colors[missing] <- grDevices::hcl.colors(
      sum(missing),
      palette = "Dark 3"
    )
  }

  names(colors) <- key$cluster_display
  colors
}

resolve_week_colors <- function(obj) {
  levels_present <- levels(droplevels(obj$week_for_plot))

  configured <- NULL

  for (candidate in c("week_colors", "gestational_week_colors")) {
    if (exists(candidate, inherits = TRUE)) {
      configured <- get(candidate, inherits = TRUE)
      break
    }
  }
  if (!is.null(configured) && all(levels_present %in% names(configured))) {
    return(configured[levels_present])
  }

  colors <- setNames(
    grDevices::hcl.colors(length(levels_present), palette = "Dark 3"),
    levels_present
  )

  colors
}

point_size_for_cells <- function(n_cells, celltype = NULL) {
  if (!is.null(celltype) && celltype == "germ") return(0.60)
  if (!is.null(celltype) && celltype %in% c("granulosa", "stroma")) return(0.65)

  if (n_cells >= 20000L) return(0.20)
  if (n_cells >= 10000L) return(0.25)
  if (n_cells >= 5000L) return(0.32)
  if (n_cells >= 1000L) return(0.48)
  if (n_cells >= 300L) return(0.80)
  1.15
}

publication_umap_theme <- function() {
  theme_publication() +
    theme(
      plot.title = element_blank(),
      axis.title = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(
        color = "black",
        linewidth = 0.35
      ),
      axis.ticks = element_line(
        color = "black",
        linewidth = 0.25
      ),
      axis.text = element_text(
        family = publication_font_family,
        size = publication_base_size
      ),
      aspect.ratio = 1,
      legend.title = element_text(
        family = publication_font_family,
        size = publication_base_size
      ),
      legend.text = element_text(
        family = publication_font_family,
        size = publication_base_size
      ),
      legend.key.height = grid::unit(0.30, "cm"),
      legend.key.width = grid::unit(0.30, "cm")
    )
}

make_two_panel_figure <- function(cds, plot_data, show_trajectory_graph, cell_size) {
  segments <- if (isTRUE(show_trajectory_graph)) {
    principal_graph_segments(cds)
  } else {
    tibble()
  }

  graph_layer <- if (nrow(segments) > 0L) {
    geom_segment(
      data = segments,
      aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE,
      color = "grey15",
      linewidth = 0.45,
      lineend = "round"
    )
  } else {
    NULL
  }

  p_cluster <- ggplot(
    plot_data,
    aes(x = UMAP_1, y = UMAP_2, color = cluster_display)
  ) +
    geom_point(size = cell_size, stroke = 0, alpha = 1, shape = 16) +
    graph_layer +
    scale_color_manual(
      values = resolve_subcluster_colors(plot_data),
      drop = FALSE,
      name = NULL
    ) +
    labs(x = NULL, y = NULL, color = NULL) +
    guides(
      color = guide_legend(
        title = NULL,
        override.aes = list(size = 2.4, alpha = 1)
      )
    ) +
    publication_umap_theme()

  p_pseudotime <- ggplot(
    plot_data,
    aes(x = UMAP_1, y = UMAP_2, color = pseudotime)
  ) +
    geom_point(size = cell_size, stroke = 0, alpha = 1, shape = 16) +
    graph_layer +
    scale_color_viridis_c(
      option = "C",
      direction = 1,
      na.value = "grey85",
      name = "Pseudotime"
    ) +
    labs(x = NULL, y = NULL) +
    publication_umap_theme()

  p_cluster + p_pseudotime + plot_layout(ncol = 2, widths = c(1, 1))
}

save_plot_pair <- function(plot, filename_stem, width, height) {
  save_publication_plot(
    plot,
    paste0(filename_stem, ".png"),
    width = width,
    height = height,
    dpi = 600
  )
  save_publication_plot(
    plot,
    paste0(filename_stem, ".pdf"),
    width = width,
    height = height
  )
}

save_qc_umaps <- function(obj, celltype) {
  safe <- sanitize_filename(celltype)
  point_size <- point_size_for_cells(ncol(obj), celltype)

  umap <- Embeddings(obj, "umap") |>
    as.data.frame() |>
    rownames_to_column("cell_id")

  colnames(umap)[2:3] <- c("UMAP_1", "UMAP_2")

  plot_data <- obj@meta.data |>
    rownames_to_column("cell_id") |>
    transmute(
      cell_id,
      week_for_plot = week_for_plot
    ) |>
    left_join(umap, by = "cell_id")

  p_week <- ggplot(
    plot_data,
    aes(x = UMAP_1, y = UMAP_2, color = week_for_plot)
  ) +
    geom_point(size = point_size, stroke = 0, alpha = 1, shape = 16) +
    scale_color_manual(values = resolve_week_colors(obj), drop = FALSE) +
    labs(x = NULL, y = NULL, color = "Gestational week") +
    guides(
      color = guide_legend(
        override.aes = list(size = 2.4, alpha = 1)
      )
    ) +
    publication_umap_theme()

  save_plot_pair(
    p_week,
    file.path(out_figure_dir, paste0("qc_umap_by_week__", safe)),
    width = 5.2,
    height = 4.8
  )
}

message("Loading annotated snRNA-seq object:")
message("  ", annotated_rds)

if (!file.exists(annotated_rds)) {
  stop("Annotated object not found: ", annotated_rds)
}

seu <- readRDS(annotated_rds)

if (!inherits(seu, "Seurat")) {
  stop("Input is not a Seurat object: ", annotated_rds)
}

required_metadata <- c(
  "major_cell_type",
  "fine_subcluster",
  "final_annotation"
)

missing_metadata <- setdiff(required_metadata, colnames(seu@meta.data))

if (length(missing_metadata) > 0L) {
  stop("Missing metadata columns: ", paste(missing_metadata, collapse = ", "))
}

DefaultAssay(seu) <- assay_use

sample_col <- sample_column(seu)

if (is.na(sample_col)) {
  stop("No sample identifier column found.")
}

counts <- get_counts_layer(seu, assay = assay_use, layer = "counts")

metadata <- prepare_metadata(
  seu@meta.data[colnames(counts), , drop = FALSE],
  sample_col = sample_col
)

if (!identical(colnames(counts), rownames(metadata))) {
  stop("Count matrix and metadata are not identically ordered.")
}
message("Repository root: ", project_root)
message("Results directory: ", results_root)
message("Sample column: ", sample_col)
message("Input cells: ", ncol(counts))
message("Input genes: ", nrow(counts))

sample_split_object <- make_sample_split_object(counts, metadata)
split_layers <- Layers(sample_split_object[[assay_use]])

write_csv(
  tibble(RNA_layer = split_layers),
  file.path(out_table_dir, "sample_split_rna_layers_used_for_sct.csv")
)

write_csv(
  analysis_plan,
  file.path(out_table_dir, "analysis_parameters.csv")
)

summary_rows <- list()

for (i in seq_len(nrow(analysis_plan))) {
  entry <- analysis_plan[i, ]

  celltype <- entry$celltype[[1]]
  workflow <- entry$workflow[[1]]
  resolution <- entry$resolution[[1]]
  expected_clusters <- entry$expected_clusters[[1]]
  pca_npcs <- entry$pca_npcs[[1]]
  show_trajectory_graph <- entry$show_trajectory_graph[[1]]

  message(strrep("=", 72))
  message("Cell type: ", celltype)
  message("Workflow: ", workflow)
  message("Fixed resolution: ", resolution)

  if (workflow == "rna") {
    obj_all <- make_germ_rna_object(counts, metadata)
    obj_all <- run_germ_rna_workflow(obj_all, resolution, pca_npcs)
  } else if (workflow == "sct_split_layers") {
    obj_all <- make_non_germ_sct_subset(sample_split_object, celltype)
    obj_all <- run_non_germ_sct_workflow(obj_all, resolution, pca_npcs)
  } else {
    stop("Unsupported workflow: ", workflow)
  }

  observed_clusters <- validate_cluster_count(
    obj_all,
    celltype,
    expected_clusters
  )

  obj_pt <- filter_pseudotime_cells(obj_all, celltype)
  root_cluster <- choose_root_cluster(obj_pt, celltype)

  message("Cells used for UMAP: ", ncol(obj_all))
  message("Cells used for pseudotime: ", ncol(obj_pt))
  message("Observed Seurat clusters: ", observed_clusters)
  message("Root cluster: ", root_cluster)

  pseudotime_result <- run_monocle_on_existing_umap(
    obj = obj_pt,
    counts_source = counts,
    celltype = celltype,
    root_cluster = root_cluster
  )

  plot_data <- make_plot_data(
    obj = obj_pt,
    pseudotime = pseudotime_result$pseudotime,
    celltype = celltype,
    root_cluster = pseudotime_result$root_cluster
  )

  safe <- sanitize_filename(celltype)
  subset_all_rds <- file.path(
    out_object_dir,
    paste0("reclustered_subset_all_cells__", safe, ".rds")
  )

  subset_pt_rds <- file.path(
    out_object_dir,
    paste0("reclustered_subset_pseudotime_cells__", safe, ".rds")
  )

  pseudotime_csv <- file.path(
    out_table_dir,
    paste0("monocle3_pseudotime_per_cell__", safe, ".csv")
  )

  saveRDS(obj_all, subset_all_rds)
  saveRDS(obj_pt, subset_pt_rds)
  write_csv(plot_data, pseudotime_csv)

  save_qc_umaps(obj_all, celltype)

  final_figure <- make_two_panel_figure(
    cds = pseudotime_result$cds,
    plot_data = plot_data,
    show_trajectory_graph = show_trajectory_graph,
    cell_size = point_size_for_cells(ncol(obj_pt), celltype)
  )

  final_figure_stem <- file.path(
    out_figure_dir,
    paste0("monocle3_umap_subcluster_pseudotime__", safe)
  )

  save_plot_pair(
    final_figure,
    final_figure_stem,
    width = 10.5,
    height = 5.2
  )

  summary_rows[[length(summary_rows) + 1L]] <- tibble(
    celltype = celltype,
    workflow = workflow,
    input_object = annotated_rds,
    uses_old_subset_rds = FALSE,
    uses_recreated_sample_split_rna_layers = workflow == "sct_split_layers",
    cells_for_umap = ncol(obj_all),
    cells_for_pseudotime = ncol(obj_pt),
    expected_seurat_clusters = expected_clusters,
    observed_seurat_clusters = observed_clusters,
    recluster_resolution = resolution,
    n_pcs_used = unique(obj_all$n_pcs_used)[[1]],
    root_cluster = pseudotime_result$root_cluster,
    root_cell_count = length(pseudotime_result$root_cells),
    trajectory_graph_shown = show_trajectory_graph,
    subset_all_rds = subset_all_rds,
    subset_pseudotime_rds = subset_pt_rds,
    pseudotime_csv = pseudotime_csv,
    final_figure_png = paste0(final_figure_stem, ".png"),
    final_figure_pdf = paste0(final_figure_stem, ".pdf")
  )
}

summary_table <- bind_rows(summary_rows)

write_csv(
  summary_table,
  file.path(out_table_dir, "monocle3_pseudotime_summary.csv")
)
writeLines(
  c(
    paste("Input annotated object:", annotated_rds),
    paste("Repository root:", project_root),
    paste("Results directory:", results_root),
    paste("Random seed:", analysis_seed),
    paste("Sample metadata column:", sample_col),
    "No old or non-GitHub subset objects are read by this workflow.",
    "Germ uses only the validated RNA workflow.",
    "Non-germ cell types use recreated sample-specific RNA count layers before SCTransform.",
    "Fixed resolutions are supplied by a separate exploratory resolution scan.",
    "Germ_0 and germ_6 are used for germ UMAP construction and excluded before pseudotime inference.",
    "Monocle3 learns trajectories on the validated Seurat cell-type UMAP coordinates.",
    "Trajectory graph lines are displayed only for germ.",
    "Exactly one QC figure type is generated per cell type: gestational week.",
    "All figures are written as PNG and PDF using the centralized plotting configuration."
  ),
  file.path(out_log_dir, "monocle3_pseudotime_run_summary.txt")
)

sink(file.path(out_log_dir, "sessionInfo.txt"))
print(sessionInfo())
sink()

message(strrep("=", 72))
message("Completed cell-type-specific UMAP and Monocle3 analysis.")
message("Outputs written to: ", results_root)
