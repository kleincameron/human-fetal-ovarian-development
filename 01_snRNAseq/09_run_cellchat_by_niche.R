#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(CellChat)
  library(igraph)
  library(dplyr)
  library(readr)
  library(tibble)
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

results_base <- Sys.getenv(
  "RESULTS_BASE",
  unset = file.path(dirname(project_root), paste0(basename(project_root), "_results"))
)

source(file.path(project_root, "config", "plotting.R"))
source(file.path(project_root, "config", "labels_colors.R"))

paths_local <- file.path(project_root, "config", "paths_local.R")
if (file.exists(paths_local)) {
  source(paths_local)
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

cellchat_results_root <- if (exists("snrna_cellchat_results_root", inherits = FALSE)) {
  snrna_cellchat_results_root
} else {
  file.path(results_base, "snRNAseq_cellchat_4_niches")
}

out_figure_dir <- file.path(cellchat_results_root, "figures")
out_table_dir <- file.path(cellchat_results_root, "tables")
out_log_dir <- file.path(cellchat_results_root, "logs")
out_object_dir <- file.path(cellchat_results_root, "objects")

dir.create(out_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_log_dir, recursive = TRUE, showWarnings = FALSE)

assay_use <- "RNA"

max_cells_per_group <- as.numeric(Sys.getenv("CELLCHAT_MAX_CELLS_PER_GROUP", unset = "Inf"))
communication_mean_type <- Sys.getenv("CELLCHAT_MEAN_TYPE", unset = "truncatedMean")
communication_trim <- as.numeric(Sys.getenv("CELLCHAT_TRIM", unset = "0.1"))
min_cells_per_group <- as.integer(Sys.getenv("CELLCHAT_MIN_CELLS_PER_GROUP", unset = "10"))
top_n_interactions <- as.integer(Sys.getenv("CELLCHAT_TOP_N_INTERACTIONS", unset = "200"))
save_cellchat_objects <- tolower(Sys.getenv("CELLCHAT_SAVE_OBJECTS", unset = "false")) %in% c("true", "t", "1", "yes", "y")

# CellChat circle-plot display scaling.
# These affect figure readability only; exported tables retain raw CellChat values.
circle_plot_width <- as.numeric(Sys.getenv("CELLCHAT_CIRCLE_WIDTH", unset = "12.5"))
circle_plot_height <- as.numeric(Sys.getenv("CELLCHAT_CIRCLE_HEIGHT", unset = "12.5"))

# Node display scale.
circle_vertex_size_max <- as.numeric(Sys.getenv("CELLCHAT_CIRCLE_VERTEX_SIZE_MAX", unset = "34"))
circle_vertex_min_scaled <- as.numeric(Sys.getenv("CELLCHAT_CIRCLE_VERTEX_MIN_SCALED", unset = "0.30"))

# Edge display scale. Positive edges are linearly rescaled but not allowed to disappear.
circle_edge_width_max <- as.numeric(Sys.getenv("CELLCHAT_CIRCLE_EDGE_WIDTH_MAX", unset = "16"))
circle_edge_min_scaled <- as.numeric(Sys.getenv("CELLCHAT_CIRCLE_EDGE_MIN_SCALED", unset = "0.24"))
circle_edge_alpha <- as.numeric(Sys.getenv("CELLCHAT_CIRCLE_EDGE_ALPHA", unset = "0.68"))
circle_edge_curved <- as.numeric(Sys.getenv("CELLCHAT_CIRCLE_EDGE_CURVED", unset = "0.22"))

# Arrows are scaled with edge display width.
circle_arrow_size_min <- as.numeric(Sys.getenv("CELLCHAT_CIRCLE_ARROW_SIZE_MIN", unset = "0.55"))
circle_arrow_size_max <- as.numeric(Sys.getenv("CELLCHAT_CIRCLE_ARROW_SIZE_MAX", unset = "1.25"))

# Labels are drawn manually outside the graph.
circle_vertex_label_cex <- as.numeric(Sys.getenv("CELLCHAT_CIRCLE_VERTEX_LABEL_CEX", unset = "1.25"))
circle_label_radius <- as.numeric(Sys.getenv("CELLCHAT_CIRCLE_LABEL_RADIUS", unset = "1.48"))
circle_plot_xlim <- as.numeric(Sys.getenv("CELLCHAT_CIRCLE_PLOT_XLIM", unset = "2.45"))

circle_margin <- as.numeric(Sys.getenv("CELLCHAT_CIRCLE_MARGIN", unset = "0.03"))
circle_base_par_cex <- as.numeric(Sys.getenv("CELLCHAT_CIRCLE_BASE_PAR_CEX", unset = "1.0"))

if (!communication_mean_type %in% c("truncatedMean", "triMean")) {
  stop("CELLCHAT_MEAN_TYPE must be 'truncatedMean' or 'triMean'.")
}

niche_list <- list(
  early_outer_cortex = c(
    "Meiotic-entry germ cells",
    "Mitotic oogonia",
    "Leptotene/zygotene germ cells",
    "Supportive pre-granulosa",
    "Morphogenetic granulosa RELN+",
    "Cortical stroma",
    "Tissue macrophages",
    "NK T cells"
  ),

  outer_cortex = c(
    "Pachytene/diplotene germ cells",
    "Leptotene/zygotene germ cells",
    "Stalled meiotic germ cells",
    "Signaling granulosa",
    "Supportive pre-granulosa",
    "Cortical stroma",
    "Signaling stroma",
    "Tissue macrophages",
    "NK T cells"
  ),
  follicle_related = c(
    "Pachytene/diplotene germ cells",
    "Primordial follicle oocytes",
    "Primordial follicle granulosa",
    "Cortical stroma",
    "Tissue macrophages",
    "NK T cells"
  ),

  degeneration_related = c(
    "Stalled meiotic germ cells",
    "Degenerating germ cells",
    "Signaling granulosa",
    "Supportive pre-granulosa",
    "Cortical stroma",
    "Signaling stroma",
    "Tissue macrophages",
    "NK T cells",
    "Atresia-stressed degenerating follicle cells",
    "Clearance-associated degenerating follicle cells"
  )
)

# ============================================================
# XENIUM-STANDARDIZED SUBCLUSTER COLORS
# ============================================================
# Names are the final_annotation labels used by the snRNA-seq object.
# Colors are matched to the established Xenium subcluster palette.
# If a displayed group is missing from this palette, it falls back to
# the established broad cell-type color.

subcluster_colors <- c(
  "Cortical stroma" = "#9C27B0",
  "Supportive pre-granulosa" = "#E69F00",
  "Primordial follicle granulosa" = "#FDB462",
  "Tissue macrophages" = "#FFD92F",
  "Morphogenetic granulosa RELN+" = "#FFA54F",
  "Signaling granulosa" = "#D55E00",
  "Primordial follicle oocytes" = "#2C7FB8",
  "Ambiguous_subcluster" = "#999999",
  "Stalled meiotic germ cells" = "#4DBBD5",
  "Angiogenic endothelia" = "#FF4FA3",
  "NK T cells" = "#E6AB02",
  "Pachytene/diplotene germ cells" = "#00C1A2",
  "Leptotene/zygotene germ cells" = "#00A087",
  "Proliferative stromal progenitors" = "#8E63CE",
  "Pericytes" = "#E41A1C",
  "Epithelial-like granulosa" = "#F28E2B",
  "Clearance-associated degenerating follicle cells" = "#B0B0B0",
  "Medullary stroma" = "#7B61FF",
  "Proliferative granulosa progenitors" = "#FF7F00",
  "Contractile VSMC" = "#B22222",
  "Signaling stroma" = "#6A3D9A",
  "Atresia-stressed degenerating follicle cells" = "#6B6B6B",
  "Late erythroid" = "#6F3F1F",
  "Perineural stroma" = "#C77CFF",
  "Matrix-remodeling granulosa" = "#C77C2E",
  "Early erythroid" = "#A66A3F",
  "Degenerating germ cells" = "#1B9E77"
)

get_celltype_color <- function(cell_type) {
  if (
    exists("celltype_colors", inherits = TRUE) &&
      cell_type %in% names(celltype_colors)
  ) {
    return(unname(celltype_colors[[cell_type]]))
  }

  fallback <- c(
    germ = "#1B9E77",
    degenerated = "#CFCFCF",
    granulosa = "#D95F02",
    stroma = "#7570B3",
    endothelial = "#E7298A",
    mural = "#66A61E",
    immune = "#E6AB02",
    erythroid = "#A6761D",
    ambiguous = "#999999",
    unknown = "#999999"
  )

  if (cell_type %in% names(fallback)) {
    unname(fallback[[cell_type]])
  } else {
    unname(fallback[["unknown"]])
  }
}

annotation_to_celltype <- function(x) {
  dplyr::case_when(
    x == "Ambiguous_subcluster" ~ "ambiguous",
    grepl("Atresia-stressed|Clearance-associated", x, ignore.case = TRUE) ~ "degenerated",
    grepl("granulosa", x, ignore.case = TRUE) ~ "granulosa",
    grepl("stroma", x, ignore.case = TRUE) ~ "stroma",
    grepl("macrophage|NK T|immune", x, ignore.case = TRUE) ~ "immune",
    grepl("endothel", x, ignore.case = TRUE) ~ "endothelial",
    grepl("pericyte|VSMC|mural", x, ignore.case = TRUE) ~ "mural",
    grepl("erythroid", x, ignore.case = TRUE) ~ "erythroid",
    grepl("germ|oogonia|oocyte|meiotic|pachytene|diplotene|leptotene|zygotene", x, ignore.case = TRUE) ~ "germ",
    TRUE ~ "unknown"
  )
}

resolve_subcluster_color_table <- function(group_names, context = "CellChat plot") {
  direct_color <- unname(subcluster_colors[group_names])
  has_direct_color <- !is.na(direct_color)

  inferred_cell_type <- annotation_to_celltype(group_names)
  fallback_color <- vapply(inferred_cell_type, get_celltype_color, character(1))

  color <- ifelse(has_direct_color, direct_color, fallback_color)

  color_source <- ifelse(
    has_direct_color,
    "xenium_subcluster_palette",
    paste0("broad_celltype_fallback:", inferred_cell_type)
  )

  missing_direct <- group_names[!has_direct_color]

  if (length(missing_direct) > 0) {
    message(
      "Using broad cell-type fallback colors for groups in ",
      context,
      ":
  ",
      paste(missing_direct, collapse = "
  ")
    )
  }

  tibble(
    group = group_names,
    color = unname(color),
    color_source = color_source,
    inferred_cell_type = inferred_cell_type
  )
}

resolve_subcluster_colors <- function(group_names, context = "CellChat plot") {
  color_table <- resolve_subcluster_color_table(group_names, context = context)
  out <- color_table$color
  names(out) <- color_table$group
  out
}


format_niche_label <- function(x) {
  tools::toTitleCase(gsub("_", " ", x))
}

extract_data_matrix <- function(seu, cells, assay = "RNA") {
  DefaultAssay(seu) <- assay

  if (!inherits(seu[[assay]], "Assay5")) {
    mat <- GetAssayData(seu, assay = assay, slot = "data")
    mat <- mat[, cells, drop = FALSE]
    if (!inherits(mat, "dgCMatrix")) mat <- as(mat, "dgCMatrix")
    return(mat)
  }

  layers <- Layers(seu[[assay]])
  data_layers <- layers[grepl("^data", layers)]

  if (length(data_layers) == 0) {
    stop("No normalized RNA data layers were found.")
  }

  preferred_layers <- c("data", "data.SeuratProject")
  preferred_layers <- preferred_layers[preferred_layers %in% data_layers]

  for (layer in preferred_layers) {
    mat <- LayerData(seu, assay = assay, layer = layer)

    if (all(cells %in% colnames(mat))) {
      mat <- mat[, cells, drop = FALSE]
      if (!inherits(mat, "dgCMatrix")) mat <- as(mat, "dgCMatrix")
      return(mat)
    }
  }

  mats <- lapply(data_layers, function(layer) {
    mat <- LayerData(seu, assay = assay, layer = layer)
    overlap <- intersect(cells, colnames(mat))

    if (length(overlap) == 0) {
      return(NULL)
    }

    mat[, overlap, drop = FALSE]
  })

  mats <- mats[!vapply(mats, is.null, logical(1))]

  if (length(mats) == 0) {
    stop("No normalized RNA data layer overlaps requested cells.")
  }
  all_genes <- sort(unique(unlist(lapply(mats, rownames))))

  mats <- lapply(mats, function(mat) {
    missing_genes <- setdiff(all_genes, rownames(mat))

    if (length(missing_genes) > 0) {
      zeros <- Matrix::Matrix(
        0,
        nrow = length(missing_genes),
        ncol = ncol(mat),
        sparse = TRUE
      )
      rownames(zeros) <- missing_genes
      colnames(zeros) <- colnames(mat)
      mat <- rbind(mat, zeros)
    }

    mat[all_genes, , drop = FALSE]
  })

  mat <- if (length(mats) == 1) {
    mats[[1]]
  } else {
    Reduce(Matrix::cbind2, mats)
  }

  mat <- mat[, !duplicated(colnames(mat)), drop = FALSE]

  missing_cells <- setdiff(cells, colnames(mat))
  if (length(missing_cells) > 0) {
    stop("Combined data layers are missing requested cells: ", length(missing_cells))
  }

  mat <- mat[, cells, drop = FALSE]

  if (!inherits(mat, "dgCMatrix")) {
    mat <- as(mat, "dgCMatrix")
  }

  mat
}

save_base_dual <- function(plot_fun, pdf_out, png_out, width = 9, height = 9, res = 600) {
  pdf(
    pdf_out,
    width = width,
    height = height,
    family = publication_font_family,
    useDingbats = FALSE
  )
  par(
    family = publication_font_family,
    cex = circle_base_par_cex,
    mar = c(0.3, 0.3, 0.3, 0.3),
    xpd = NA
  )
  plot_fun()
  dev.off()

  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(
      png_out,
      width = width,
      height = height,
      units = "in",
      res = res,
      background = "white"
    )
  } else {
    png(
      png_out,
      width = width,
      height = height,
      units = "in",
      res = res,
      bg = "white"
    )
  }

  par(
    family = publication_font_family,
    cex = circle_base_par_cex,
    mar = c(0.3, 0.3, 0.3, 0.3),
    xpd = NA
  )
  plot_fun()
  dev.off()
}


rescale_positive_values_for_circle <- function(x, min_scaled = 0.1, max_scaled = 1) {
  out <- as.numeric(x)
  names(out) <- names(x)

  positive <- out > 0 & !is.na(out)
  if (!any(positive)) return(out)

  rng <- range(out[positive], na.rm = TRUE)
  if (!all(is.finite(rng))) return(out)

  if (diff(rng) == 0) {
    out[positive] <- max_scaled
    return(out)
  }

  out[positive] <- min_scaled + (out[positive] - rng[1]) / diff(rng) * (max_scaled - min_scaled)
  out
}

rescale_positive_matrix_for_circle <- function(mat, min_scaled = 0.1, max_scaled = 1) {
  out <- mat
  positive <- out > 0 & !is.na(out)
  if (!any(positive)) return(out)

  rng <- range(out[positive], na.rm = TRUE)
  if (!all(is.finite(rng))) return(out)

  if (diff(rng) == 0) {
    out[positive] <- max_scaled
    return(out)
  }

  out[positive] <- min_scaled + (out[positive] - rng[1]) / diff(rng) * (max_scaled - min_scaled)
  out
}

abbreviate_circle_group_labels <- function(x) {
  dplyr::recode(
    x,
    "Atresia-stressed degenerating follicle cells" = "Atresia-stressed",
    "Clearance-associated degenerating follicle cells" = "Clearance-associated",
    .default = x
  )
}

make_circle_layout <- function(group_names) {
  n <- length(group_names)
  theta <- seq(from = pi / 2, to = pi / 2 - 2 * pi, length.out = n + 1)[seq_len(n)]

  layout <- cbind(
    x = cos(theta),
    y = sin(theta)
  )

  rownames(layout) <- group_names

  list(
    layout = layout,
    theta = theta
  )
}

plot_cellchat_strength_circle <- function(weight_mat, group_size, color_use, display_labels) {
  stopifnot(length(display_labels) == nrow(weight_mat))
  stopifnot(length(display_labels) == ncol(weight_mat))

  plot_weight_matrix <- rescale_positive_matrix_for_circle(
    weight_mat,
    min_scaled = circle_edge_min_scaled,
    max_scaled = 1
  )

  rownames(plot_weight_matrix) <- display_labels
  colnames(plot_weight_matrix) <- display_labels

  plot_group_size <- rescale_positive_values_for_circle(
    group_size,
    min_scaled = circle_vertex_min_scaled,
    max_scaled = 1
  )
  names(plot_group_size) <- display_labels

  color_use_display <- color_use
  names(color_use_display) <- display_labels

  circle_layout <- make_circle_layout(display_labels)
  layout_mat <- circle_layout$layout
  theta <- circle_layout$theta
  names(theta) <- display_labels

  g <- igraph::graph_from_adjacency_matrix(
    as.matrix(plot_weight_matrix),
    mode = "directed",
    weighted = TRUE,
    diag = TRUE
  )

  vertex_names <- igraph::V(g)$name

  vertex_size <- plot_group_size[vertex_names] * circle_vertex_size_max
  vertex_color <- unname(color_use_display[vertex_names])

  if (length(igraph::E(g)) > 0) {
    edge_weight <- igraph::E(g)$weight
    edge_width <- edge_weight * circle_edge_width_max

    edge_arrow_size <- circle_arrow_size_min +
      edge_weight * (circle_arrow_size_max - circle_arrow_size_min)

    edge_ends <- igraph::ends(g, igraph::E(g), names = TRUE)
    edge_source <- edge_ends[, 1]
    edge_color <- grDevices::adjustcolor(
      color_use_display[edge_source],
      alpha.f = circle_edge_alpha
    )
    edge_loop_angle <- theta[edge_source]
  } else {
    edge_width <- numeric(0)
    edge_arrow_size <- numeric(0)
    edge_color <- character(0)
    edge_loop_angle <- numeric(0)
  }

  graphics::plot(
    g,
    layout = layout_mat[vertex_names, , drop = FALSE],
    rescale = FALSE,
    xlim = c(-circle_plot_xlim, circle_plot_xlim),
    ylim = c(-circle_plot_xlim, circle_plot_xlim),
    asp = 1,
    margin = rep(circle_margin, 4),
    vertex.color = vertex_color,
    vertex.frame.color = "white",
    vertex.frame.width = 0.4,
    vertex.size = vertex_size,
    vertex.label = NA,
    edge.color = edge_color,
    edge.width = edge_width,
    edge.arrow.size = edge_arrow_size,
    edge.curved = circle_edge_curved,
    edge.loop.angle = edge_loop_angle
  )

  label_xy <- layout_mat[display_labels, , drop = FALSE] * circle_label_radius

  for (i in seq_along(display_labels)) {
    x <- label_xy[i, 1]
    y <- label_xy[i, 2]

    adj_x <- if (x > 0.10) {
      0
    } else if (x < -0.10) {
      1
    } else {
      0.5
    }

    graphics::text(
      x = x,
      y = y,
      labels = display_labels[i],
      adj = c(adj_x, 0.5),
      cex = circle_vertex_label_cex,
      family = publication_font_family
    )
  }
}


make_group_pair_table <- function(cellchat, niche_name) {
  count_mat <- cellchat@net$count
  weight_mat <- cellchat@net$weight

  pair_grid <- expand.grid(
    source = rownames(count_mat),
    target = colnames(count_mat),
    stringsAsFactors = FALSE
  )

  pair_grid |>
    mutate(
      niche = niche_name,
      n_interactions = as.numeric(count_mat[cbind(source, target)]),
      interaction_weight = as.numeric(weight_mat[cbind(source, target)])
    ) |>
    select(niche, source, target, n_interactions, interaction_weight) |>
    arrange(niche, desc(interaction_weight), desc(n_interactions), source, target)
}

run_cellchat_for_niche <- function(seu, niche_name, niche_labels) {
  cat("\n============================================================\n")
  cat("Running CellChat niche: ", niche_name, "\n", sep = "")
  cat("============================================================\n")

  missing_niche_labels <- setdiff(niche_labels, unique(seu$final_annotation))

  if (length(missing_niche_labels) > 0) {
    stop(
      "Niche '", niche_name, "' contains labels not found in final_annotation:\n",
      paste(missing_niche_labels, collapse = "\n")
    )
  }
  invisible(resolve_subcluster_colors(
    niche_labels,
    context = paste0("niche definition: ", niche_name)
  ))

  cells_initial <- rownames(seu@meta.data)[seu$final_annotation %in% niche_labels]

  meta_use <- seu@meta.data[cells_initial, , drop = FALSE] |>
    mutate(
      cell = rownames(seu@meta.data)[seu$final_annotation %in% niche_labels],
      cellchat_group = factor(as.character(final_annotation), levels = niche_labels)
    )

  qc_counts <- meta_use |>
    count(cellchat_group, name = "n_cells") |>
    mutate(niche = niche_name, downsampled = FALSE) |>
    relocate(niche, downsampled) |>
    arrange(cellchat_group)

  too_small <- qc_counts |>
    filter(n_cells < min_cells_per_group)

  if (nrow(too_small) > 0) {
    warning(
      "Some groups in niche '", niche_name, "' have fewer than ",
      min_cells_per_group, " cells:\n",
      paste(too_small$cellchat_group, too_small$n_cells, sep = "=", collapse = ", ")
    )
  }

  if (is.finite(max_cells_per_group)) {
    set.seed(1234)

    meta_use <- meta_use |>
      group_by(cellchat_group) |>
      group_modify(~ {
        n_take <- min(max_cells_per_group, nrow(.x))
        slice_sample(.x, n = n_take)
      }) |>
      ungroup()

    meta_use$cellchat_group <- factor(
      as.character(meta_use$cellchat_group),
      levels = niche_labels
    )

    qc_counts <- meta_use |>
      count(cellchat_group, name = "n_cells") |>
      mutate(niche = niche_name, downsampled = TRUE) |>
      relocate(niche, downsampled) |>
      arrange(cellchat_group)
  }

  cells_use <- meta_use$cell
  cat("Cells used: ", length(cells_use), "\n", sep = "")

  data_input <- extract_data_matrix(
    seu = seu,
    cells = cells_use,
    assay = assay_use
  )

  meta <- data.frame(
    group = as.character(meta_use$cellchat_group[match(colnames(data_input), meta_use$cell)]),
    row.names = colnames(data_input),
    stringsAsFactors = FALSE
  )

  meta$group <- factor(meta$group, levels = niche_labels)

  stopifnot(identical(colnames(data_input), rownames(meta)))

  cellchat <- createCellChat(
    object = data_input,
    meta = meta,
    group.by = "group"
  )

  cellchat <- addMeta(cellchat, meta = meta)
  cellchat <- setIdent(cellchat, ident.use = "group")

  cellchat@DB <- CellChatDB.human

  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)

  if (communication_mean_type == "truncatedMean") {
    cellchat <- computeCommunProb(
      cellchat,
      type = "truncatedMean",
      trim = communication_trim,
      population.size = TRUE
    )
  } else {
    cellchat <- computeCommunProb(
      cellchat,
      type = "triMean",
      population.size = TRUE
    )
  }

  cellchat <- filterCommunication(
    cellchat,
    min.cells = min_cells_per_group
  )
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)

  if (save_cellchat_objects) {
    dir.create(out_object_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(
      cellchat,
      file.path(out_object_dir, paste0("cellchat_", niche_name, ".rds"))
    )
  }

  comm_all <- subsetCommunication(cellchat)

  top_interactions <- if (nrow(comm_all) > 0 && "prob" %in% colnames(comm_all)) {
    comm_all |>
      arrange(desc(prob)) |>
      slice_head(n = top_n_interactions) |>
      mutate(niche = niche_name) |>
      relocate(niche)
  } else {
    tibble(niche = niche_name)
  }

  comm_pathway <- tryCatch(
    subsetCommunication(cellchat, slot.name = "netP"),
    error = function(e) tibble()
  )

  pathway_interactions <- if (nrow(comm_pathway) > 0) {
    comm_pathway |>
      mutate(niche = niche_name) |>
      relocate(niche)
  } else {
    tibble(niche = niche_name)
  }

  aggregate_network <- make_group_pair_table(cellchat, niche_name)

  group_names <- rownames(cellchat@net$count)

  color_use <- resolve_subcluster_colors(
    group_names,
    context = paste0("CellChat circle plot: ", niche_name)
  )

  group_size <- as.numeric(table(cellchat@idents))
  names(group_size) <- names(table(cellchat@idents))
  group_size <- group_size[group_names]

  if (any(is.na(group_size))) {
    stop("Could not match group sizes to CellChat network groups for niche: ", niche_name)
  }

  display_group_names <- abbreviate_circle_group_labels(group_names)

  save_base_dual(
    plot_fun = function() {
      plot_cellchat_strength_circle(
        weight_mat = cellchat@net$weight,
        group_size = group_size,
        color_use = color_use,
        display_labels = display_group_names
      )
    },
    pdf_out = file.path(out_figure_dir, paste0("Circle_", niche_name, "_interaction_strength.pdf")),
    png_out = file.path(out_figure_dir, paste0("Circle_", niche_name, "_interaction_strength.png")),
    width = circle_plot_width,
    height = circle_plot_height
  )

  summary_row <- tibble(
    niche = niche_name,
    cells_used = length(cells_use),
    n_groups = length(niche_labels),
    groups = paste(niche_labels, collapse = ";"),
    n_lr_interactions_inferred = nrow(comm_all),
    n_pathway_interactions_inferred = nrow(comm_pathway),
    max_interaction_probability = ifelse(
      nrow(comm_all) > 0 && "prob" %in% colnames(comm_all),
      max(comm_all$prob, na.rm = TRUE),
      NA_real_
    ),
    total_interaction_weight = sum(cellchat@net$weight, na.rm = TRUE),
    n_nonzero_group_edges = sum(cellchat@net$count > 0, na.rm = TRUE)
  )
  list(
    qc_counts = qc_counts,
    top_interactions = top_interactions,
    pathway_interactions = pathway_interactions,
    aggregate_network = aggregate_network,
    summary_row = summary_row
  )
}

cat("Loading annotated snRNA-seq object:\n", annotated_rds, "\n")

stopifnot(file.exists(annotated_rds))

seu <- readRDS(annotated_rds)
stopifnot(inherits(seu, "Seurat"))

DefaultAssay(seu) <- assay_use

required_metadata <- c("fine_subcluster", "final_annotation")
missing_metadata <- setdiff(required_metadata, colnames(seu@meta.data))

if (length(missing_metadata) > 0) {
  stop("Annotated object is missing required metadata columns: ", paste(missing_metadata, collapse = ", "))
}

seu@meta.data <- seu@meta.data |>
  mutate(
    fine_subcluster = as.character(fine_subcluster),
    final_annotation = as.character(final_annotation)
  )

niche_labels_all <- unique(unlist(niche_list))
missing_labels <- setdiff(niche_labels_all, unique(seu$final_annotation))

if (length(missing_labels) > 0) {
  stop(
    "The annotated object does not contain all expected niche labels:\n",
    paste(missing_labels, collapse = "\n")
  )
}

color_map_all <- bind_rows(lapply(names(niche_list), function(niche_name) {
  resolve_subcluster_color_table(
    niche_list[[niche_name]],
    context = paste0("niche definition: ", niche_name)
  ) |>
    mutate(niche = niche_name) |>
    relocate(niche)
}))

write_csv(
  color_map_all,
  file.path(out_table_dir, "CellChat_snRNAseq_4_niches_color_map.csv")
)
results <- lapply(names(niche_list), function(niche_name) {
  run_cellchat_for_niche(
    seu = seu,
    niche_name = niche_name,
    niche_labels = niche_list[[niche_name]]
  )
})

names(results) <- names(niche_list)

qc_all <- bind_rows(lapply(results, `[[`, "qc_counts"))
summary_table <- bind_rows(lapply(results, `[[`, "summary_row"))
top200_all <- bind_rows(lapply(results, `[[`, "top_interactions"))
pathway_all <- bind_rows(lapply(results, `[[`, "pathway_interactions"))
aggregate_network_all <- bind_rows(lapply(results, `[[`, "aggregate_network"))

write_csv(
  qc_all,
  file.path(out_table_dir, "QC_cellchat_group_counts_by_niche.csv")
)

write_csv(
  summary_table,
  file.path(out_table_dir, "CellChat_snRNAseq_4_niches_summary.csv")
)

write_csv(
  top200_all,
  file.path(out_table_dir, "COMBINED_all_niches_top200_by_probability.csv")
)

write_csv(
  pathway_all,
  file.path(out_table_dir, "COMBINED_all_niches_pathway_level_interactions.csv")
)

write_csv(
  aggregate_network_all,
  file.path(out_table_dir, "COMBINED_all_niches_aggregate_network_by_group_pair.csv")
)
writeLines(
  c(
    paste("Input annotated snRNA-seq object:", annotated_rds),
    paste("Output root:", cellchat_results_root),
    "This script runs CellChat for four biologically defined snRNA-seq niches.",
    "It starts from the GitHub-generated annotated fetal snRNA-seq object.",
    "Minimal retained outputs only.",
    "Final retained figures: one enlarged custom-rendered interaction-strength CellChat circle plot per niche.",
    "Final retained tables: combined top LR interactions per niche, pathway-level interactions, aggregate group-pair network summary, CellChat niche summary, group-count QC, and the CellChat color map.",
    "Circle plots are rendered from CellChat aggregate interaction weights using custom igraph plotting for larger nodes, thicker minimum edges, larger arrows, outside labels, shortened degenerating-state labels, and no plot titles. This scaling affects figures only, not exported interaction tables.",
    "No full per-niche LR tables, per-niche pathway tables, raw count/weight matrices, bubble plots, count-circle plots, signaling heatmaps, or per-niche sessionInfo files are saved.",
    paste("Assay:", assay_use),
    paste("Niches:", paste(names(niche_list), collapse = ", ")),
    paste("Communication mean type:", communication_mean_type),
    paste("Communication trim:", communication_trim),
    paste("Minimum cells per group:", min_cells_per_group),
    paste("Maximum cells per group:", max_cells_per_group),
    paste("Top interactions retained per niche:", top_n_interactions),
    paste("Save full CellChat objects:", save_cellchat_objects),
    paste("Circle plot width:", circle_plot_width),
    paste("Circle plot height:", circle_plot_height),
    paste("Circle vertex size max:", circle_vertex_size_max),
    paste("Circle vertex minimum scaled display weight:", circle_vertex_min_scaled),
    paste("Circle edge width max:", circle_edge_width_max),
    paste("Circle edge minimum scaled display weight:", circle_edge_min_scaled),
    paste("Circle vertex label cex:", circle_vertex_label_cex),
    paste("Circle label radius:", circle_label_radius),
    paste("Circle arrow size min:", circle_arrow_size_min),
    paste("Circle arrow size max:", circle_arrow_size_max),
    "CellChat circle plots use Xenium-standardized subcluster colors with broad cell-type fallback for missing subcluster-specific colors."
  ),
  file.path(out_log_dir, "CellChat_snRNAseq_4_niches_notes.txt")
)

sink(file.path(out_log_dir, "sessionInfo_all_niches.txt"))
print(sessionInfo())
sink()

cat("\nDone. Minimal CellChat outputs written to:\n", cellchat_results_root, "\n")
