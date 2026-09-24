#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tibble)
  library(ggtext)
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

out_root <- file.path(results_base, "snRNAseq_ligand_receptor_grouped_dotplot")
out_figure_dir <- file.path(out_root, "figures")
out_table_dir <- file.path(out_root, "tables")
out_log_dir <- file.path(out_root, "logs")

dir.create(out_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_log_dir, recursive = TRUE, showWarnings = FALSE)

pdf_out <- file.path(out_figure_dir, "DotPlot_ligand_receptor_grouped_selected_subclusters.pdf")
png_out <- file.path(out_figure_dir, "DotPlot_ligand_receptor_grouped_selected_subclusters.png")

gene_map_out <- file.path(out_table_dir, "ligand_receptor_grouped_genes_found.csv")
block_out <- file.path(out_table_dir, "ligand_receptor_grouped_block_spans.csv")
group_out <- file.path(out_table_dir, "ligand_receptor_grouped_pair_boundaries.csv")
source_data_out <- file.path(out_table_dir, "ligand_receptor_grouped_dotplot_source_data.csv")
cell_counts_out <- file.path(out_table_dir, "QC_selected_subcluster_cell_counts.csv")

assay_use <- "RNA"

match_gene <- function(gene, all_genes) {
  hit <- all_genes[toupper(all_genes) == toupper(gene)]
  if (length(hit) == 0) {
    return(NA_character_)
  }
  hit[1]
}
get_hline_positions <- function(top_down_labels, type_map) {
  n <- length(top_down_labels)
  if (n <= 1) {
    return(numeric(0))
  }

  pos <- c()

  for (i in seq_len(n - 1)) {
    t1 <- type_map[top_down_labels[i]]
    t2 <- type_map[top_down_labels[i + 1]]

    if (!is.na(t1) && !is.na(t2) && t1 != t2) {
      pos <- c(pos, (n - i) + 0.5)
    }
  }

  pos
}

extract_data_matrix <- function(seu, cells, genes, assay = "RNA") {
  DefaultAssay(seu) <- assay

  if (!inherits(seu[[assay]], "Assay5")) {
    mat <- GetAssayData(seu, assay = assay, slot = "data")
    mat <- mat[genes, cells, drop = FALSE]
    if (!inherits(mat, "dgCMatrix")) {
      mat <- as(mat, "dgCMatrix")
    }
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

    if (all(cells %in% colnames(mat)) && all(genes %in% rownames(mat))) {
      mat <- mat[genes, cells, drop = FALSE]
      if (!inherits(mat, "dgCMatrix")) {
        mat <- as(mat, "dgCMatrix")
      }
      return(mat)
    }
  }

  mats <- lapply(data_layers, function(layer) {
    mat <- LayerData(seu, assay = assay, layer = layer)
    overlap_cells <- intersect(cells, colnames(mat))

    if (length(overlap_cells) == 0) {
      return(NULL)
    }

    missing_genes <- setdiff(genes, rownames(mat))

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

    mat[genes, overlap_cells, drop = FALSE]
  })
  mats <- mats[!vapply(mats, is.null, logical(1))]

  if (length(mats) == 0) {
    stop("No normalized RNA data layer overlaps the requested cells.")
  }

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

  mat <- mat[genes, cells, drop = FALSE]

  if (!inherits(mat, "dgCMatrix")) {
    mat <- as(mat, "dgCMatrix")
  }

  mat
}

make_dotplot_source_data <- function(expr_mat, group_vector, group_levels) {
  rows <- lapply(group_levels, function(group_name) {
    cells <- names(group_vector)[group_vector == group_name]

    if (length(cells) == 0) {
      return(tibble())
    }

    mat <- expr_mat[, cells, drop = FALSE]

    tibble(
      final_annotation = group_name,
      gene = rownames(mat),
      avg_exp = as.numeric(Matrix::rowMeans(mat)),
      pct_exp = as.numeric(Matrix::rowMeans(mat > 0) * 100)
    )
  })

  out <- bind_rows(rows)

  out <- out |>
    group_by(gene) |>
    mutate(
      avg_exp_scaled = as.numeric(scale(avg_exp)),
      avg_exp_scaled = ifelse(is.na(avg_exp_scaled), 0, avg_exp_scaled),
      avg_exp_scaled = pmax(pmin(avg_exp_scaled, 2.5), -2.5)
    ) |>
    ungroup()

  out
}

save_plot_dual <- function(p, pdf_out, png_out, width = 16, height = 10, dpi = 600) {
  ggsave(
    pdf_out,
    plot = p,
    width = width,
    height = height,
    device = "pdf",
    family = publication_font_family,
    useDingbats = FALSE,
    limitsize = FALSE
  )

  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(
      png_out,
      plot = p,
      width = width,
      height = height,
      dpi = dpi,
      device = ragg::agg_png,
      background = "white",
      limitsize = FALSE
    )
  } else {
    ggsave(
      png_out,
      plot = p,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white",
      limitsize = FALSE
    )
  }
}

subcluster_order <- c(
  "Mitotic oogonia",
  "Meiotic-entry germ cells",
  "Leptotene/zygotene germ cells",
  "Pachytene/diplotene germ cells",
  "Primordial follicle oocytes",
  "Stalled meiotic germ cells",
  "Degenerating germ cells",
  "Atresia-stressed degenerating follicle cells",
  "Clearance-associated degenerating follicle cells",
  "Supportive pre-granulosa",
  "Signaling granulosa",
  "Primordial follicle granulosa",
  "Morphogenetic granulosa RELN+",
  "Cortical stroma",
  "Signaling stroma",
  "Tissue macrophages",
  "NK T cells"
)
subcluster_type_map <- c(
  "Mitotic oogonia" = "germ",
  "Meiotic-entry germ cells" = "germ",
  "Leptotene/zygotene germ cells" = "germ",
  "Pachytene/diplotene germ cells" = "germ",
  "Primordial follicle oocytes" = "germ",
  "Stalled meiotic germ cells" = "germ",
  "Degenerating germ cells" = "germ",
  "Atresia-stressed degenerating follicle cells" = "degenerated",
  "Clearance-associated degenerating follicle cells" = "degenerated",
  "Supportive pre-granulosa" = "granulosa",
  "Signaling granulosa" = "granulosa",
  "Primordial follicle granulosa" = "granulosa",
  "Morphogenetic granulosa RELN+" = "granulosa",
  "Cortical stroma" = "stroma",
  "Signaling stroma" = "stroma",
  "Tissue macrophages" = "immune",
  "NK T cells" = "immune"
)

celltype_label_colors <- c(
  germ = if ("germ" %in% names(celltype_colors)) celltype_colors[["germ"]] else "#1B9E77",
  degenerated = "#6B6B6B",
  granulosa = if ("granulosa" %in% names(celltype_colors)) celltype_colors[["granulosa"]] else "#D95F02",
  stroma = if ("stroma" %in% names(celltype_colors)) celltype_colors[["stroma"]] else "#7570B3",
  immune = if ("immune" %in% names(celltype_colors)) celltype_colors[["immune"]] else "#E6AB02"
)

pair_defs <- list(
  list(
    group_id = "pair_01",
    group_label = "NTN1 -> UNC5D",
    ligands = c("NTN1"),
    receptors = c("UNC5D")
  ),
  list(
    group_id = "pair_02",
    group_label = "SLIT2/SLIT3 -> ROBO1/ROBO2",
    ligands = c("SLIT2", "SLIT3"),
    receptors = c("ROBO1", "ROBO2")
  ),
  list(
    group_id = "pair_03",
    group_label = "TENM3/TENM4/FLRT3 -> ADGRL2/ADGRL3",
    ligands = c("TENM3", "TENM4", "FLRT3"),
    receptors = c("ADGRL2", "ADGRL3")
  ),
  list(
    group_id = "pair_04",
    group_label = "NRXN1/NRXN3 -> NLGN1/LRRTM3",
    ligands = c("NRXN1", "NRXN3"),
    receptors = c("NLGN1", "LRRTM3")
  ),
  list(
    group_id = "pair_05",
    group_label = "CNTN4 -> PTPRG",
    ligands = c("CNTN4"),
    receptors = c("PTPRG")
  ),
  list(
    group_id = "pair_06",
    group_label = "NRG3/NRG4 -> ERBB4",
    ligands = c("NRG3", "NRG4"),
    receptors = c("ERBB4")
  ),
  list(
    group_id = "pair_07",
    group_label = "RELN -> VLDLR",
    ligands = c("RELN"),
    receptors = c("VLDLR")
  ),
  list(
    group_id = "pair_08",
    group_label = "LAMA2 -> ADGRG6",
    ligands = c("LAMA2"),
    receptors = c("ADGRG6")
  ),
  list(
    group_id = "pair_09",
    group_label = "LRFN5 -> PTPRD",
    ligands = c("LRFN5"),
    receptors = c("PTPRD")
  ),
  list(
    group_id = "pair_10",
    group_label = "LRRC4C -> PTPRF",
    ligands = c("LRRC4C"),
    receptors = c("PTPRF")
  )
)

homophilic_genes <- c(
  "CADM1", "CADM2",
  "NECTIN3",
  "PTPRM", "PTPRK", "PTPRT", "PTPRU",
  "NCAM1", "NCAM2",
  "NEGR1"
)

requested_tbl <- bind_rows(
  lapply(pair_defs, function(x) {
    bind_rows(
      data.frame(
        overall_group_id = x$group_id,
        overall_group_label = x$group_label,
        subgroup_id = paste0(x$group_id, "_ligand"),
        subgroup_label = paste0(x$group_label, " [ligand]"),
        gene_class = "Ligand",
        requested_gene = x$ligands,
        stringsAsFactors = FALSE
      ),
      data.frame(
        overall_group_id = x$group_id,
        overall_group_label = x$group_label,
        subgroup_id = paste0(x$group_id, "_receptor"),
        subgroup_label = paste0(x$group_label, " [receptor]"),
        gene_class = "Receptor",
        requested_gene = x$receptors,
        stringsAsFactors = FALSE
      )
    )
  }),
  data.frame(
    overall_group_id = "homophilic",
    overall_group_label = "Homophilic",
    subgroup_id = "homophilic",
    subgroup_label = "Homophilic",
    gene_class = "Homophilic",
    requested_gene = homophilic_genes,
    stringsAsFactors = FALSE
  )
)
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

missing_requested_subclusters <- setdiff(
  subcluster_order,
  unique(seu$final_annotation)
)

if (length(missing_requested_subclusters) > 0) {
  warning(
    "These requested annotations were not found in the object and will be omitted:\n",
    paste(missing_requested_subclusters, collapse = ", ")
  )
}

subcluster_order_present <- subcluster_order[
  subcluster_order %in% unique(seu$final_annotation)
]

if (length(subcluster_order_present) == 0) {
  stop("None of the requested annotations were found in the object.")
}

missing_type_map <- setdiff(subcluster_order_present, names(subcluster_type_map))
if (length(missing_type_map) > 0) {
  stop("Missing subcluster_type_map entries for: ", paste(missing_type_map, collapse = ", "))
}

subcluster_factor_levels <- rev(subcluster_order_present)

cells_use <- rownames(seu@meta.data)[seu$final_annotation %in% subcluster_order_present]

cell_counts <- seu@meta.data[cells_use, , drop = FALSE] |>
  mutate(final_annotation = factor(final_annotation, levels = subcluster_order_present)) |>
  count(final_annotation, name = "n_cells") |>
  arrange(final_annotation)

write_csv(cell_counts, cell_counts_out)

all_genes <- rownames(seu[[assay_use]])

gene_map_tbl <- requested_tbl |>
  mutate(
    matched_gene = vapply(requested_gene, match_gene, character(1), all_genes = all_genes),
    found = !is.na(matched_gene)
  )

write_csv(gene_map_tbl, gene_map_out)

genes_plot <- gene_map_tbl |>
  filter(found) |>
  distinct(matched_gene, .keep_all = TRUE) |>
  pull(matched_gene)

if (length(genes_plot) == 0) {
  stop("None of the requested genes were found in the object.")
}

gene_positions <- data.frame(
  matched_gene = genes_plot,
  x_pos = seq_along(genes_plot),
  stringsAsFactors = FALSE
)

block_spans <- gene_map_tbl |>
  filter(found) |>
  distinct(
    overall_group_id,
    overall_group_label,
    subgroup_id,
    subgroup_label,
    gene_class,
    matched_gene
  ) |>
  left_join(gene_positions, by = "matched_gene") |>
  group_by(
    overall_group_id,
    overall_group_label,
    subgroup_id,
    subgroup_label,
    gene_class
  ) |>
  summarise(
    start = min(x_pos),
    end = max(x_pos),
    n_genes = n(),
    .groups = "drop"
  ) |>
  arrange(start)
pair_boundaries <- gene_map_tbl |>
  filter(found) |>
  distinct(overall_group_id, overall_group_label, matched_gene) |>
  left_join(gene_positions, by = "matched_gene") |>
  group_by(overall_group_id, overall_group_label) |>
  summarise(
    start = min(x_pos),
    end = max(x_pos),
    .groups = "drop"
  ) |>
  arrange(start)

write_csv(block_spans, block_out)
write_csv(pair_boundaries, group_out)

expr_mat <- extract_data_matrix(
  seu = seu,
  cells = cells_use,
  genes = genes_plot,
  assay = assay_use
)

group_vector <- seu@meta.data[colnames(expr_mat), "final_annotation", drop = TRUE]
group_vector <- factor(group_vector, levels = subcluster_factor_levels)
names(group_vector) <- colnames(expr_mat)

dp_df <- make_dotplot_source_data(
  expr_mat = expr_mat,
  group_vector = group_vector,
  group_levels = subcluster_factor_levels
)

gene_annotation_tbl <- gene_map_tbl |>
  filter(found) |>
  distinct(
    matched_gene,
    overall_group_id,
    overall_group_label,
    subgroup_id,
    subgroup_label,
    gene_class
  )

dp_df <- dp_df |>
  left_join(gene_annotation_tbl, by = c("gene" = "matched_gene")) |>
  mutate(
    gene = factor(gene, levels = genes_plot),
    final_annotation = factor(final_annotation, levels = subcluster_factor_levels)
  )

write_csv(
  dp_df |>
    mutate(
      gene = as.character(gene),
      final_annotation = as.character(final_annotation)
    ),
  source_data_out
)

y_axis_label_map <- setNames(
  paste0(
    "<span style='color:",
    celltype_label_colors[subcluster_type_map[subcluster_order_present]],
    ";'>",
    subcluster_order_present,
    "</span>"
  ),
  subcluster_order_present
)

hline_pos <- get_hline_positions(
  top_down_labels = subcluster_order_present,
  type_map = subcluster_type_map
)

bg_fill_values <- c(
  Ligand = "grey96",
  Receptor = "grey88",
  Homophilic = "grey80"
)

vline_pos <- pair_boundaries$end + 0.5
if (length(vline_pos) > 0) {
  vline_pos <- vline_pos[-length(vline_pos)]
}

plot_width <- max(16, min(28, 0.42 * length(genes_plot) + 6))
plot_height <- max(8, 0.38 * length(subcluster_order_present) + 3)
p <- ggplot(dp_df, aes(x = gene, y = final_annotation)) +
  geom_rect(
    data = block_spans,
    inherit.aes = FALSE,
    aes(
      xmin = start - 0.5,
      xmax = end + 0.5,
      ymin = -Inf,
      ymax = Inf,
      fill = gene_class
    ),
    color = NA
  ) +
  scale_fill_manual(
    name = "Gene class",
    values = bg_fill_values,
    breaks = c("Ligand", "Receptor", "Homophilic")
  ) +
  geom_vline(
    xintercept = vline_pos,
    linewidth = 0.5,
    color = "grey35"
  ) +
  geom_hline(
    yintercept = hline_pos,
    linewidth = 0.5,
    color = "grey35"
  ) +
  geom_point(
    aes(size = pct_exp, color = avg_exp_scaled)
  ) +
  scale_size(
    name = "% expressing",
    range = c(0, 6),
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100)
  ) +
  scale_color_gradientn(
    name = "Scaled average expression",
    colors = c("grey85", "#9ECAE1", "#3182BD", "#08519C"),
    limits = c(-2.5, 2.5)
  ) +
  scale_y_discrete(labels = y_axis_label_map) +
  theme_publication() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.4),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
    axis.text.y = ggtext::element_markdown(size = 8),
    axis.title = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 8),
    plot.margin = margin(10, 12, 10, 10),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title = "Grouped ligand-receptor and homophilic gene expression by fetal ovary subcluster"
  )

save_plot_dual(
  p,
  pdf_out = pdf_out,
  png_out = png_out,
  width = plot_width,
  height = plot_height,
  dpi = 600
)
writeLines(
  c(
    paste("Input annotated snRNA-seq object:", annotated_rds),
    paste("Output root:", out_root),
    "This script generates a grouped ligand-receptor and homophilic gene-expression dot plot for selected fetal ovary snRNA-seq subclusters.",
    "It starts from the GitHub-generated annotated fetal snRNA-seq object and uses final_annotation directly.",
    "No old manuscript-object path is used.",
    paste("Assay:", assay_use),
    paste("Selected annotations retained:", paste(subcluster_order_present, collapse = "; ")),
    paste("Genes requested:", nrow(gene_map_tbl)),
    paste("Genes found:", sum(gene_map_tbl$found)),
    paste("Genes plotted:", length(genes_plot)),
    "Dot size encodes percentage of cells expressing each gene.",
    "Dot color encodes per-gene scaled average normalized expression across selected annotations.",
    "Retained outputs: final PDF/PNG figure, figure source data, gene matching table, group/block coordinate tables, selected subcluster cell-count QC, notes, and sessionInfo."
  ),
  file.path(out_log_dir, "ligand_receptor_grouped_dotplot_notes.txt")
)

sink(file.path(out_log_dir, "sessionInfo_ligand_receptor_grouped_dotplot.txt"))
print(sessionInfo())
sink()

cat("\nDone. Outputs written to:\n", out_root, "\n")
cat("Dotplot PDF:\n", pdf_out, "\n")
cat("Dotplot PNG:\n", png_out, "\n")
cat("Source data:\n", source_data_out, "\n")
