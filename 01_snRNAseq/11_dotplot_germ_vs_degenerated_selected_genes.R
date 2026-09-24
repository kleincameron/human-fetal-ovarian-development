#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(grid)
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

out_root <- file.path(results_base, "snRNAseq_germ_vs_degenerated_selected_genes")
out_figure_dir <- file.path(out_root, "figures")
out_table_dir <- file.path(out_root, "tables")
out_log_dir <- file.path(out_root, "logs")

dir.create(out_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_log_dir, recursive = TRUE, showWarnings = FALSE)

assay_use <- "RNA"

# ============================================================
# Figure / test settings
# ============================================================

genes_to_plot <- c(
  "ANXA1", "CCL2", "CCL3", "CCL4", "CXCL2", "CXCL8",
  "FOS", "JUN", "MKI67", "POU5F1"
)

divider_after <- 6L
fdr_cutoff <- 0.05

figure_width <- 8.0
figure_height <- 3.8
figure_dpi <- 600

point_size_min <- 1.0
point_size_max <- 9.0

expr_low <- "#9FD3FF"
expr_mid <- "#4C9FE6"
expr_high <- "#0B2E59"
# ============================================================
# Helpers
# ============================================================

match_gene <- function(gene, all_genes) {
  hit <- all_genes[toupper(all_genes) == toupper(gene)]
  if (length(hit) == 0) {
    return(NA_character_)
  }
  hit[1]
}

save_plot_dual <- function(plot, filename_prefix,
                           width = figure_width,
                           height = figure_height,
                           dpi = figure_dpi) {
  pdf_file <- paste0(filename_prefix, ".pdf")
  png_file <- paste0(filename_prefix, ".png")

  ggsave(
    filename = pdf_file,
    plot = plot,
    width = width,
    height = height,
    device = "pdf",
    family = publication_font_family,
    useDingbats = FALSE,
    bg = "white",
    limitsize = FALSE
  )

  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(
      filename = png_file,
      plot = plot,
      width = width,
      height = height,
      dpi = dpi,
      device = ragg::agg_png,
      background = "white",
      limitsize = FALSE
    )
  } else {
    ggsave(
      filename = png_file,
      plot = plot,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white",
      limitsize = FALSE
    )
  }

  invisible(c(pdf = pdf_file, png = png_file))
}

build_germ_degenerated_group <- function(meta) {
  if ("major_cell_type" %in% colnames(meta)) {
    grp <- as.character(meta$major_cell_type)
    grp[grp == "quiescent"] <- "degenerated"
    grp[!(grp %in% c("germ", "degenerated"))] <- NA_character_
    return(grp)
  }

  if ("cell_type_quiescent" %in% colnames(meta)) {
    grp <- as.character(meta$cell_type_quiescent)
    grp[grp == "quiescent"] <- "degenerated"
    grp[!(grp %in% c("germ", "degenerated"))] <- NA_character_
    return(grp)
  }

  if ("fine_subcluster" %in% colnames(meta)) {
    fs <- as.character(meta$fine_subcluster)

    grp <- rep(NA_character_, length(fs))
    grp[grepl("^germ_", fs)] <- "germ"
    grp[fs %in% c("germ_0", "germ_6")] <- "degenerated"

    return(grp)
  }

  stop(
    "Could not construct germ/degenerated grouping. ",
    "Expected one of: major_cell_type, cell_type_quiescent, or fine_subcluster."
  )
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

calculate_dotplot_statistics <- function(expr_mat, groups, genes_use) {
  stats_list <- lapply(genes_use, function(gene) {
    x_germ <- as.numeric(expr_mat[gene, groups == "germ"])
    x_deg <- as.numeric(expr_mat[gene, groups == "degenerated"])

    tibble(
      gene = gene,
      group = c("germ", "degenerated"),
      avg_expr = c(mean(x_germ), mean(x_deg)),
      pct_expressing = c(mean(x_germ > 0) * 100, mean(x_deg > 0) * 100),
      n_cells = c(sum(groups == "germ"), sum(groups == "degenerated"))
    )
  })

  stats_df <- bind_rows(stats_list) |>
    group_by(gene) |>
    mutate(
      avg_expr_scaled = as.numeric(scale(avg_expr)),
      avg_expr_scaled = ifelse(is.na(avg_expr_scaled), 0, avg_expr_scaled),
      avg_expr_scaled = pmax(pmin(avg_expr_scaled, 2.5), -2.5)
    ) |>
    ungroup()

  stats_df
}

run_significance_tests <- function(expr_mat, groups, genes_use, fdr_cutoff = 0.05) {
  sig_df <- lapply(genes_use, function(gene) {
    x_germ <- as.numeric(expr_mat[gene, groups == "germ"])
    x_deg <- as.numeric(expr_mat[gene, groups == "degenerated"])

    wt <- suppressWarnings(
      wilcox.test(x = x_germ, y = x_deg, exact = FALSE)
    )

    mean_germ <- mean(x_germ)
    mean_deg <- mean(x_deg)

    tibble(
      gene = gene,
      p_value = wt$p.value,
      mean_expr_germ = mean_germ,
      mean_expr_degenerated = mean_deg,
      higher_group = ifelse(mean_germ >= mean_deg, "germ", "degenerated")
    )
  }) |>
    bind_rows() |>
    mutate(
      FDR = p.adjust(p_value, method = "BH"),
      significant = FDR < fdr_cutoff,
      label = ifelse(significant, "*", "")
    )

  sig_df
}
# ============================================================
# Load object
# ============================================================

cat("Loading annotated snRNA-seq object:\n  ", annotated_rds, "\n", sep = "")
stopifnot(file.exists(annotated_rds))

seu <- readRDS(annotated_rds)
stopifnot(inherits(seu, "Seurat"))

DefaultAssay(seu) <- assay_use

stopifnot(assay_use %in% names(seu@assays))

group_vec <- build_germ_degenerated_group(seu@meta.data)
names(group_vec) <- colnames(seu)

keep_cells <- names(group_vec)[!is.na(group_vec)]
group_vec <- group_vec[keep_cells]

if (!all(c("germ", "degenerated") %in% unique(group_vec))) {
  stop(
    "Both germ and degenerated groups must be present. Found groups: ",
    paste(sort(unique(group_vec)), collapse = ", ")
  )
}

all_genes <- rownames(seu[[assay_use]])

gene_map <- tibble(
  requested_gene = genes_to_plot,
  matched_gene = vapply(genes_to_plot, match_gene, character(1), all_genes = all_genes),
  found = !is.na(matched_gene)
)

write_csv(
  gene_map,
  file.path(out_table_dir, "DotPlot_selected_genes_germ_vs_degenerated_gene_map.csv")
)

genes_present <- gene_map |>
  filter(found) |>
  pull(matched_gene)

genes_missing <- gene_map |>
  filter(!found) |>
  pull(requested_gene)

if (length(genes_present) == 0) {
  stop("None of the requested genes were found in the RNA assay.")
}

if (length(genes_missing) > 0) {
  warning(
    "The following requested genes were not found and will be skipped:\n",
    paste(genes_missing, collapse = ", ")
  )
}

expr_mat <- extract_data_matrix(
  seu = seu,
  cells = keep_cells,
  genes = genes_present,
  assay = assay_use
)

group_vec <- group_vec[colnames(expr_mat)]

cat("Cells used:\n")
print(table(group_vec))

cat("Genes plotted:\n")
print(genes_present)

# ============================================================
# Compute stats
# ============================================================

dot_stats <- calculate_dotplot_statistics(
  expr_mat = expr_mat,
  groups = group_vec,
  genes_use = genes_present
)

sig_stats <- run_significance_tests(
  expr_mat = expr_mat,
  groups = group_vec,
  genes_use = genes_present,
  fdr_cutoff = fdr_cutoff
)

plot_df <- dot_stats |>
  left_join(
    sig_stats |>
      select(gene, p_value, FDR, significant, higher_group, label),
    by = "gene"
  ) |>
  mutate(
    gene = factor(gene, levels = genes_present),
    group = factor(group, levels = c("germ", "degenerated")),
    y_num = ifelse(group == "germ", 1, 2)
  )

asterisk_df <- sig_stats |>
  filter(significant) |>
  transmute(
    gene = factor(gene, levels = genes_present),
    y_num = ifelse(higher_group == "germ", 1.18, 2.18),
    label = "*"
  )

expr_abs_max <- max(abs(plot_df$avg_expr_scaled), na.rm = TRUE)
expr_limits <- c(-expr_abs_max, expr_abs_max)
# ============================================================
# Figure
# ============================================================

p <- ggplot(plot_df, aes(x = gene, y = y_num)) +
  geom_point(
    aes(size = pct_expressing, color = avg_expr_scaled)
  ) +
  geom_text(
    data = asterisk_df,
    aes(x = gene, y = y_num, label = label),
    inherit.aes = FALSE,
    size = 4.5,
    family = publication_font_family,
    fontface = "bold"
  ) +
  geom_vline(
    xintercept = divider_after + 0.5,
    linewidth = 0.3,
    color = "grey60"
  ) +
  scale_y_continuous(
    breaks = c(1, 2),
    labels = c("germ", "degenerated"),
    limits = c(0.5, 2.5),
    expand = c(0, 0)
  ) +
  scale_size_continuous(
    name = "% expressing",
    range = c(point_size_min, point_size_max),
    breaks = c(0, 3, 6, 9, 12)
  ) +
  scale_color_gradient2(
    name = "Avg expr\n(scaled)",
    low = expr_low,
    mid = expr_mid,
    high = expr_high,
    midpoint = 0,
    limits = expr_limits
  ) +
  guides(
    color = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = grid::unit(1.00, "in"),
      barwidth = grid::unit(0.18, "in")
    ),
    size = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      override.aes = list(color = "black")
    )
  ) +
  labs(
    title = "Selected genes in germ vs degenerated groups",
    x = "Genes",
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(
    base_size = publication_base_size,
    base_family = publication_font_family
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0, size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 9),
    legend.title = element_text(size = 8, lineheight = 0.9),
    legend.text = element_text(size = 7),
    legend.box = "vertical",
    legend.box.just = "top",
    legend.box.margin = margin(t = 14, r = 8, b = 4, l = 8),
    legend.margin = margin(t = 8, r = 4, b = 4, l = 4),
    legend.spacing.y = grid::unit(0.20, "in"),
    plot.margin = margin(t = 12, r = 28, b = 8, l = 6)
  )

figure_prefix <- file.path(
  out_figure_dir,
  "DotPlot_selected_genes_germ_vs_degenerated"
)

save_plot_dual(
  plot = p,
  filename_prefix = figure_prefix,
  width = figure_width,
  height = figure_height,
  dpi = figure_dpi
)

# ============================================================
# Outputs
# ============================================================

write_csv(
  plot_df |>
    mutate(
      gene = as.character(gene),
      group = as.character(group)
    ) |>
    select(
      gene,
      group,
      avg_expr,
      avg_expr_scaled,
      pct_expressing,
      n_cells,
      p_value,
      FDR,
      significant,
      higher_group
    ),
  file.path(out_table_dir, "DotPlot_selected_genes_germ_vs_degenerated_source_data.csv")
)

write_csv(
  sig_stats,
  file.path(out_table_dir, "DotPlot_selected_genes_germ_vs_degenerated_significance.csv")
)

write_csv(
  tibble(group = names(table(group_vec)), n_cells = as.integer(table(group_vec))),
  file.path(out_table_dir, "QC_germ_vs_degenerated_cell_counts.csv")
)

writeLines(
  c(
    paste("Input annotated snRNA-seq object:", annotated_rds),
    paste("Output root:", out_root),
    "This script generates a selected-gene dot plot comparing germ and degenerated fetal ovary snRNA-seq groups.",
    "It starts from the GitHub-generated annotated fetal snRNA-seq object.",
    "No old manuscript-object path is used.",
    paste("Assay:", assay_use),
    paste("Genes requested:", paste(genes_to_plot, collapse = ", ")),
    paste("Genes plotted:", paste(genes_present, collapse = ", ")),
    paste("Missing genes:", ifelse(length(genes_missing) == 0, "none", paste(genes_missing, collapse = ", "))),
    paste("Germ cells:", sum(group_vec == "germ")),
    paste("Degenerated cells:", sum(group_vec == "degenerated")),
    paste("FDR cutoff:", fdr_cutoff),
    "Dot size encodes percentage of cells expressing each gene.",
    "Dot color encodes per-gene scaled average normalized expression across the two groups.",
    "A single '*' is shown when the Wilcoxon rank-sum test passes BH-adjusted FDR cutoff.",
    "The asterisk is placed above the group with higher mean expression.",
    paste("Figure PDF:", paste0(figure_prefix, ".pdf")),
    paste("Figure PNG:", paste0(figure_prefix, ".png"))
  ),
  file.path(out_log_dir, "DotPlot_selected_genes_germ_vs_degenerated_notes.txt")
)

sink(file.path(out_log_dir, "sessionInfo_DotPlot_selected_genes_germ_vs_degenerated.txt"))
print(sessionInfo())
sink()

cat("\nDone.\n")
cat("Outputs written to:\n  ", out_root, "\n", sep = "")
cat("Figure PDF:\n  ", paste0(figure_prefix, ".pdf"), "\n", sep = "")
cat("Figure PNG:\n  ", paste0(figure_prefix, ".png"), "\n", sep = "")
