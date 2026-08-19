#!/usr/bin/env Rscript

required_packages <- c(
  "dplyr",
  "readr",
  "tibble",
  "pheatmap",
  "grid"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
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

deg_results_root <- if (exists("snrna_deg_results_root", inherits = FALSE)) {
  snrna_deg_results_root
} else {
  file.path(results_base, "snRNAseq_DEGs_celltype_subcluster")
}

out_root <- if (exists("snrna_hallmark_ora_results_root", inherits = FALSE)) {
  snrna_hallmark_ora_results_root
} else {
  file.path(results_base, "snRNAseq_Hallmark_ORA_by_celltype")
}

subcluster_deg_root <- file.path(
  deg_results_root,
  "tables",
  "subcluster_full_by_cell_type"
)

if (!dir.exists(subcluster_deg_root)) {
  stop(
    "Subcluster DEG directory not found: ", subcluster_deg_root,
    "\nRun 01_snRNAseq/14_generate_deg_tables_celltype_subcluster.R first."
  )
}

out_figure_dir <- file.path(out_root, "figures")
out_table_dir <- file.path(out_root, "tables")
out_log_dir <- file.path(out_root, "logs")

dir.create(out_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_log_dir, recursive = TRUE, showWarnings = FALSE)

deg_padj_cut <- as.numeric(Sys.getenv("HALLMARK_ORA_DEG_PADJ_CUT", unset = "0.05"))
deg_logfc_min <- as.numeric(Sys.getenv("HALLMARK_ORA_DEG_LOGFC_MIN", unset = "0"))
deg_pct1_min <- as.numeric(Sys.getenv("HALLMARK_ORA_DEG_PCT1_MIN", unset = "0.25"))
universe_pct_min <- as.numeric(Sys.getenv("HALLMARK_ORA_UNIVERSE_PCT_MIN", unset = "0.10"))

ora_padj_cut <- as.numeric(Sys.getenv("HALLMARK_ORA_PADJ_CUT", unset = "0.05"))
min_overlap <- as.integer(Sys.getenv("HALLMARK_ORA_MIN_OVERLAP", unset = "2"))
top_terms_per_subcluster <- as.integer(Sys.getenv("HALLMARK_ORA_TOP_TERMS_PER_SUBCLUSTER", unset = "6"))
max_terms_per_celltype <- as.integer(Sys.getenv("HALLMARK_ORA_MAX_TERMS_PER_CELLTYPE", unset = "30"))
max_significant_fraction <- as.numeric(Sys.getenv("HALLMARK_ORA_MAX_SIGNIFICANT_FRACTION", unset = "0.80"))
min_specificity_range <- as.numeric(Sys.getenv("HALLMARK_ORA_MIN_SPECIFICITY_RANGE", unset = "0.50"))

heatmap_value_cap <- as.numeric(Sys.getenv("HALLMARK_ORA_HEATMAP_VALUE_CAP", unset = "4"))
heatmap_width <- as.numeric(Sys.getenv("HALLMARK_ORA_HEATMAP_WIDTH", unset = "9"))
heatmap_height <- as.numeric(Sys.getenv("HALLMARK_ORA_HEATMAP_HEIGHT", unset = "7"))
heatmap_dpi <- as.integer(Sys.getenv("HALLMARK_ORA_HEATMAP_DPI", unset = "600"))
cluster_columns <- tolower(Sys.getenv("HALLMARK_ORA_CLUSTER_COLUMNS", unset = "true")) %in% c("true", "t", "1", "yes", "y")

hallmark_gmt <- if (exists("hallmark_gmt", inherits = FALSE)) {
  hallmark_gmt
} else {
  Sys.getenv("MSIGDB_HALLMARK_GMT", unset = "")
}

available_cell_types <- basename(list.dirs(subcluster_deg_root, recursive = FALSE))

cell_types_env <- Sys.getenv("HALLMARK_ORA_CELL_TYPES", unset = "")
if (nzchar(cell_types_env)) {
  target_cell_types <- trimws(unlist(strsplit(cell_types_env, ",")))
} else {
  target_cell_types <- c(
    celltype_order[celltype_order %in% available_cell_types],
    setdiff(sort(available_cell_types), celltype_order)
  )
}
target_cell_types <- setdiff(target_cell_types, "endothelial")

if (length(target_cell_types) == 0) {
  stop("No target cell types with subcluster DEG directories were found.")
}

read_hallmark_sets_from_gmt <- function(gmt_path) {
  if (!file.exists(gmt_path)) {
    stop("Hallmark GMT file not found: ", gmt_path)
  }

  lines <- readLines(gmt_path, warn = FALSE)
  lines <- lines[nzchar(lines)]

  out <- lapply(lines, function(line) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(fields) < 3) return(NULL)

    tibble(
      pathway = fields[[1]],
      gene = unique(fields[-c(1, 2)])
    )
  })

  bind_rows(out) |>
    filter(!is.na(pathway), pathway != "", !is.na(gene), gene != "") |>
    mutate(
      pathway = toupper(pathway),
      gene = toupper(gene)
    ) |>
    distinct(pathway, gene)
}

read_hallmark_sets_from_msigdbr <- function() {
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop(
      "No Hallmark GMT file was provided and package 'msigdbr' is not installed.\n",
      "Set MSIGDB_HALLMARK_GMT to a local Hallmark GMT file, define hallmark_gmt in config/paths_local.R, or install msigdbr.",
      call. = FALSE
    )
  }

  msig_fun <- getExportedValue("msigdbr", "msigdbr")
  args <- names(formals(msig_fun))

  msig <- if ("collection" %in% args) {
    msig_fun(species = "Homo sapiens", collection = "H")
  } else {
    msig_fun(species = "Homo sapiens", category = "H")
  }

  pathway_col <- if ("gs_name" %in% colnames(msig)) "gs_name" else "gs_exact_source"
  gene_col <- if ("gene_symbol" %in% colnames(msig)) "gene_symbol" else "human_gene_symbol"

  tibble(
    pathway = toupper(as.character(msig[[pathway_col]])),
    gene = toupper(as.character(msig[[gene_col]]))
  ) |>
    filter(!is.na(pathway), pathway != "", !is.na(gene), gene != "") |>
    distinct(pathway, gene)
}

read_hallmark_sets <- function() {
  if (nzchar(hallmark_gmt)) {
    read_hallmark_sets_from_gmt(hallmark_gmt)
  } else {
    read_hallmark_sets_from_msigdbr()
  }
}

first_existing_col <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) NA_character_ else hit[[1]]
}

first_matching_col <- function(df, pattern) {
  hit <- grep(pattern, colnames(df), value = TRUE, ignore.case = TRUE)
  if (length(hit) == 0) NA_character_ else hit[[1]]
}

infer_subcluster_from_filename <- function(path, cell_type) {
  filename <- tools::file_path_sans_ext(basename(path))
  cell_type_pattern <- paste0(cell_type, "_[0-9]+")
  hit <- regmatches(filename, regexpr(cell_type_pattern, filename))
  if (length(hit) == 1 && nzchar(hit)) return(hit)
  hit <- regmatches(filename, regexpr("[A-Za-z]+_[0-9]+", filename))
  if (length(hit) == 1 && nzchar(hit)) return(hit)
  sub("^DEG_", "", sub("_vs.*$", "", filename))
}

label_subcluster <- function(x) {
  out <- as.character(x)
  if (exists("subcluster_label_map", inherits = TRUE)) {
    mapped <- unname(subcluster_label_map[out])
    out[!is.na(mapped)] <- mapped[!is.na(mapped)]
  }
  out
}

order_subclusters <- function(cell_type, subclusters) {
  if (!exists("subcluster_label_map", inherits = TRUE)) {
    return(sort(subclusters))
  }

  id_order <- names(subcluster_label_map)
  label_order <- unname(subcluster_label_map)

  if (cell_type == "degenerated") {
    preferred_ids <- c("germ_0", "germ_6")
  } else {
    preferred_ids <- id_order[startsWith(id_order, paste0(cell_type, "_"))]
  }

  preferred <- c(preferred_ids, label_order[match(preferred_ids, id_order)])
  c(intersect(preferred, subclusters), setdiff(sort(subclusters), preferred))
}

standardize_deg_table <- function(path, cell_type) {
  df <- readr::read_csv(path, show_col_types = FALSE)

  if (nrow(df) == 0) return(tibble())

  gene_col <- first_existing_col(df, c("gene", "Gene", "symbol", "SYMBOL", "features", "feature"))
  if (is.na(gene_col)) gene_col <- colnames(df)[[1]]

  padj_col <- first_existing_col(df, c("p_val_adj", "padj", "p_adj", "adj_p_val", "FDR", "q_value"))
  logfc_col <- first_matching_col(df, "^avg_log2FC$|^avg_logFC$|log2FC|logFC")
  pct1_col <- first_existing_col(df, c("pct.1", "pct_1", "pct1", "pct_in_cluster", "pct_in"))
  pct2_col <- first_existing_col(df, c("pct.2", "pct_2", "pct2", "pct_out_cluster", "pct_out"))

  if (any(is.na(c(padj_col, logfc_col, pct1_col)))) {
    stop("Could not identify required DEG columns in: ", path)
  }

  subcluster_col <- first_existing_col(df, c(
    "subcluster",
    "fine_subcluster",
    "cluster",
    "ident.1",
    "group",
    "comparison_group"
  ))
  subcluster_value <- if (is.na(subcluster_col)) {
    infer_subcluster_from_filename(path, cell_type)
  } else {
    as.character(df[[subcluster_col]])
  }

  tibble(
    gene = toupper(as.character(df[[gene_col]])),
    cell_type = cell_type,
    subcluster = as.character(subcluster_value),
    subcluster_label = label_subcluster(as.character(subcluster_value)),
    p_val_adj = suppressWarnings(as.numeric(df[[padj_col]])),
    avg_log2FC = suppressWarnings(as.numeric(df[[logfc_col]])),
    pct.1 = suppressWarnings(as.numeric(df[[pct1_col]])),
    pct.2 = if (!is.na(pct2_col)) suppressWarnings(as.numeric(df[[pct2_col]])) else NA_real_
  ) |>
    filter(!is.na(gene), gene != "", !is.na(subcluster), subcluster != "") |>
    distinct(cell_type, subcluster, gene, .keep_all = TRUE)
}

read_celltype_degs <- function(cell_type) {
  deg_dir <- file.path(subcluster_deg_root, cell_type)

  if (!dir.exists(deg_dir)) {
    warning("Skipping missing DEG directory: ", deg_dir)
    return(tibble())
  }

  deg_files <- list.files(deg_dir, pattern = "\\.csv$", full.names = TRUE)
  deg_files <- deg_files[
    !grepl("top50|summary|qc|combined", basename(deg_files), ignore.case = TRUE)
  ]

  if (length(deg_files) == 0) {
    warning("No full DEG CSV files found for cell type: ", cell_type)
    return(tibble())
  }

  bind_rows(lapply(deg_files, standardize_deg_table, cell_type = cell_type))
}

format_hallmark_label <- function(pathway) {
  x <- gsub("^HALLMARK_", "", pathway)
  toupper(gsub("[^A-Za-z0-9]+", "_", x))
}

run_hallmark_ora <- function(deg_tbl, hallmark_tbl, cell_type) {
  subclusters <- order_subclusters(cell_type, unique(deg_tbl$subcluster))
  hallmark_tbl <- hallmark_tbl |>
    filter(gene %in% deg_tbl$gene) |>
    group_by(pathway) |>
    filter(n_distinct(gene) >= min_overlap) |>
    ungroup()

  if (nrow(hallmark_tbl) == 0) return(tibble())

  term_genes_by_pathway <- split(hallmark_tbl$gene, hallmark_tbl$pathway)
  pathways <- names(term_genes_by_pathway)

  records <- vector("list", length(subclusters) * length(pathways))
  idx <- 0L

  for (subcluster_i in subclusters) {
    deg_i <- deg_tbl |>
      filter(subcluster == subcluster_i)

    universe_i <- deg_i |>
      filter(
        !is.na(pct.1),
        pct.1 >= universe_pct_min | (!is.na(pct.2) & pct.2 >= universe_pct_min)
      ) |>
      pull(gene) |>
      unique() |>
      sort()

    candidate_i <- deg_i |>
      filter(
        !is.na(p_val_adj),
        !is.na(avg_log2FC),
        !is.na(pct.1),
        p_val_adj <= deg_padj_cut,
        avg_log2FC > deg_logfc_min,
        pct.1 >= deg_pct1_min
      ) |>
      pull(gene) |>
      unique() |>
      intersect(universe_i)

    universe_n <- length(universe_i)
    candidate_n <- length(candidate_i)

    for (pathway_i in pathways) {
      term_genes <- intersect(term_genes_by_pathway[[pathway_i]], universe_i)

      a <- length(intersect(candidate_i, term_genes))
      b <- candidate_n - a
      c <- length(term_genes) - a
      d <- universe_n - a - b - c

      odds_ratio <- ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))

      p_value <- if (candidate_n == 0 || length(term_genes) == 0 || universe_n == 0) {
        1
      } else {
        suppressWarnings(
          fisher.test(
            matrix(c(a, b, c, d), nrow = 2, byrow = TRUE),
            alternative = "greater"
          )$p.value
        )
      }
      idx <- idx + 1L
      records[[idx]] <- tibble(
        cell_type = cell_type,
        subcluster = subcluster_i,
        subcluster_label = label_subcluster(subcluster_i),
        pathway = pathway_i,
        pathway_label = format_hallmark_label(pathway_i),
        universe_n = universe_n,
        candidate_gene_n = candidate_n,
        pathway_size_in_universe = length(term_genes),
        overlap_n = a,
        overlap_genes = paste(sort(intersect(candidate_i, term_genes)), collapse = ";"),
        p_value = p_value,
        odds_ratio = odds_ratio,
        log2_odds_ratio = log2(odds_ratio)
      )
    }
  }

  bind_rows(records) |>
    group_by(subcluster) |>
    mutate(p_adj = p.adjust(p_value, method = "BH")) |>
    ungroup() |>
    arrange(cell_type, subcluster, p_adj, desc(log2_odds_ratio), pathway)
}

select_heatmap_terms <- function(ora_tbl) {
  n_subclusters <- length(unique(ora_tbl$subcluster))

  term_stats <- ora_tbl |>
    group_by(pathway, pathway_label) |>
    summarise(
      best_p_adj = min(p_adj, na.rm = TRUE),
      max_log2_odds_ratio = max(log2_odds_ratio, na.rm = TRUE),
      min_log2_odds_ratio = min(log2_odds_ratio, na.rm = TRUE),
      range_log2_odds_ratio = max_log2_odds_ratio - min_log2_odds_ratio,
      n_significant_positive = sum(
        p_adj <= ora_padj_cut &
          log2_odds_ratio > 0 &
          overlap_n >= min_overlap,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) |>
    filter(
      n_significant_positive > 0,
      n_significant_positive <= ceiling(n_subclusters * max_significant_fraction),
      range_log2_odds_ratio >= min_specificity_range
    )

  selected_by_subcluster <- ora_tbl |>
    semi_join(term_stats, by = c("pathway", "pathway_label")) |>
    filter(
      p_adj <= ora_padj_cut,
      log2_odds_ratio > 0,
      overlap_n >= min_overlap
    ) |>
    group_by(subcluster) |>
    arrange(p_adj, desc(log2_odds_ratio), desc(overlap_n), .by_group = TRUE) |>
    slice_head(n = top_terms_per_subcluster) |>
    ungroup() |>
    distinct(pathway, pathway_label)

  term_stats |>
    semi_join(selected_by_subcluster, by = c("pathway", "pathway_label")) |>
    arrange(best_p_adj, desc(range_log2_odds_ratio), desc(max_log2_odds_ratio), pathway_label) |>
    slice_head(n = max_terms_per_celltype)
}

make_heatmap_matrix <- function(ora_tbl, selected_terms, subclusters) {
  term_order <- selected_terms$pathway
  display_labels <- label_subcluster(subclusters)

  mat <- matrix(
    0,
    nrow = length(term_order),
    ncol = length(subclusters),
    dimnames = list(term_order, subclusters)
  )

  star_mat <- matrix(
    "",
    nrow = length(term_order),
    ncol = length(subclusters),
    dimnames = list(term_order, subclusters)
  )

  for (i in seq_len(nrow(ora_tbl))) {
    pathway_i <- ora_tbl$pathway[[i]]
    subcluster_i <- ora_tbl$subcluster[[i]]

    if (pathway_i %in% term_order && subcluster_i %in% subclusters) {
      value <- ora_tbl$log2_odds_ratio[[i]]
      mat[pathway_i, subcluster_i] <- value

      if (
        !is.na(ora_tbl$p_adj[[i]]) &&
          ora_tbl$p_adj[[i]] <= ora_padj_cut &&
          value > 0 &&
          ora_tbl$overlap_n[[i]] >= min_overlap
      ) {
        star_mat[pathway_i, subcluster_i] <- "*"
      }
    }
  }

  mat <- pmax(pmin(mat, heatmap_value_cap), -heatmap_value_cap)

  row_labels <- selected_terms$pathway_label[match(rownames(mat), selected_terms$pathway)]
  row_labels <- make.unique(row_labels, sep = "_")

  rownames(mat) <- row_labels
  rownames(star_mat) <- row_labels
  colnames(mat) <- make.unique(display_labels, sep = "_")
  colnames(star_mat) <- colnames(mat)

  list(matrix = mat, stars = star_mat)
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

plot_celltype_heatmap <- function(ora_tbl, selected_terms, cell_type) {
  subclusters <- order_subclusters(cell_type, unique(ora_tbl$subcluster))
  hm <- make_heatmap_matrix(ora_tbl, selected_terms, subclusters)

  heatmap_colors <- grDevices::colorRampPalette(
    c("#313695", "#FFFFFF", "#A50026")
  )(101)

  heatmap_breaks <- seq(
    -heatmap_value_cap,
    heatmap_value_cap,
    length.out = length(heatmap_colors) + 1
  )
  title_text <- paste0(
    "Hallmark ORA - ", cell_type,
    "\nvalue = signed log2(odds ratio); '*' ORA padj < ", ora_padj_cut,
    "\nDEG filter: padj <= ", deg_padj_cut,
    ", log2FC > ", deg_logfc_min,
    ", pct.1 >= ", deg_pct1_min
  )

  ph <- pheatmap::pheatmap(
    hm$matrix,
    color = heatmap_colors,
    breaks = heatmap_breaks,
    cluster_rows = nrow(hm$matrix) > 1,
    cluster_cols = cluster_columns && ncol(hm$matrix) > 1,
    display_numbers = hm$stars,
    number_color = "black",
    fontsize = publication_base_size,
    fontsize_row = max(5, publication_base_size - 1),
    fontsize_col = max(5, publication_base_size - 1),
    angle_col = 90,
    border_color = NA,
    main = title_text,
    silent = TRUE
  )

  pdf_out <- file.path(out_figure_dir, paste0("Hallmark_ORA_heatmap_", cell_type, ".pdf"))
  png_out <- file.path(out_figure_dir, paste0("Hallmark_ORA_heatmap_", cell_type, ".png"))

  save_pheatmap_dual(
    ph = ph,
    pdf_out = pdf_out,
    png_out = png_out,
    width = heatmap_width,
    height = heatmap_height,
    dpi = heatmap_dpi
  )

  readr::write_csv(
    as.data.frame(hm$matrix) |>
      rownames_to_column("pathway_label"),
    file.path(out_table_dir, paste0("Hallmark_ORA_heatmap_matrix_", cell_type, ".csv"))
  )
}

hallmark_tbl <- read_hallmark_sets()

readr::write_csv(
  hallmark_tbl |>
    count(pathway, name = "n_genes") |>
    arrange(pathway),
  file.path(out_table_dir, "Hallmark_gene_set_sizes.csv")
)

all_ora <- list()
all_selected <- list()
summary_rows <- list()

for (cell_type in target_cell_types) {
  cat("\n============================================================\n")
  cat("Running Hallmark ORA for cell type: ", cell_type, "\n", sep = "")
  cat("============================================================\n")

  deg_tbl <- read_celltype_degs(cell_type)

  if (nrow(deg_tbl) == 0) {
    warning("No DEG rows available for cell type: ", cell_type)
    next
  }

  ora_tbl <- run_hallmark_ora(deg_tbl, hallmark_tbl, cell_type)

  if (nrow(ora_tbl) == 0) {
    warning("No Hallmark ORA results available for cell type: ", cell_type)
    next
  }

  selected_terms <- select_heatmap_terms(ora_tbl)

  readr::write_csv(
    ora_tbl,
    file.path(out_table_dir, paste0("Hallmark_ORA_all_terms_", cell_type, ".csv"))
  )

  readr::write_csv(
    selected_terms,
    file.path(out_table_dir, paste0("Hallmark_ORA_selected_heatmap_terms_", cell_type, ".csv"))
  )

  if (nrow(selected_terms) > 0) {
    plot_celltype_heatmap(ora_tbl, selected_terms, cell_type)
  } else {
    warning("No specific significant Hallmark terms selected for heatmap: ", cell_type)
  }

  all_ora[[cell_type]] <- ora_tbl
  all_selected[[cell_type]] <- selected_terms |>
    mutate(cell_type = cell_type) |>
    relocate(cell_type)

  summary_rows[[cell_type]] <- tibble(
    cell_type = cell_type,
    n_subclusters = length(unique(deg_tbl$subcluster)),
    n_deg_rows = nrow(deg_tbl),
    n_candidate_gene_subcluster_pairs = deg_tbl |>
      filter(
        !is.na(p_val_adj),
        !is.na(avg_log2FC),
        !is.na(pct.1),
        p_val_adj <= deg_padj_cut,
        avg_log2FC > deg_logfc_min,
        pct.1 >= deg_pct1_min
      ) |>
      distinct(subcluster, gene) |>
      nrow(),
    n_ora_rows = nrow(ora_tbl),
    n_selected_heatmap_terms = nrow(selected_terms)
  )
}

if (length(all_ora) > 0) {
  readr::write_csv(
    bind_rows(all_ora),
    file.path(out_table_dir, "COMBINED_Hallmark_ORA_all_terms_by_celltype.csv")
  )
}

if (length(all_selected) > 0) {
  readr::write_csv(
    bind_rows(all_selected),
    file.path(out_table_dir, "COMBINED_Hallmark_ORA_selected_heatmap_terms_by_celltype.csv")
  )
}

if (length(summary_rows) > 0) {
  readr::write_csv(
    bind_rows(summary_rows),
    file.path(out_table_dir, "Hallmark_ORA_summary_by_celltype.csv")
  )
}

writeLines(
  c(
    paste("Input DEG root:", deg_results_root),
    paste("Output root:", out_root),
    paste("Hallmark source:", ifelse(nzchar(hallmark_gmt), hallmark_gmt, "msigdbr")),
    paste("Cell types:", paste(target_cell_types, collapse = ", ")),
    paste("Excluded cell types:", "endothelial"),
    paste("DEG adjusted p-value cutoff:", deg_padj_cut),
    paste("DEG log2FC minimum:", deg_logfc_min),
    paste("DEG pct.1 minimum:", deg_pct1_min),
    paste("Universe pct minimum:", universe_pct_min),
    paste("ORA adjusted p-value cutoff:", ora_padj_cut),
    paste("Minimum overlap:", min_overlap),
    paste("Maximum significant subcluster fraction:", max_significant_fraction),
    paste("Minimum specificity range:", min_specificity_range),
    paste("Top terms per subcluster:", top_terms_per_subcluster),
    paste("Maximum heatmap terms per cell type:", max_terms_per_celltype),
    paste("Heatmap value cap:", heatmap_value_cap),
    "ORA uses Fisher's exact test for positive subcluster DEGs against each subcluster's tested gene universe.",
    "Heatmap rows are selected from significant, subcluster-specific Hallmark terms.",
    "Heatmap columns use final publication annotation labels."
  ),
  file.path(out_log_dir, "Hallmark_ORA_heatmaps_by_celltype_notes.txt")
)

sink(file.path(out_log_dir, "sessionInfo_Hallmark_ORA_heatmaps_by_celltype.txt"))
print(sessionInfo())
sink()

cat("\nDone. Hallmark ORA outputs written to:\n", out_root, "\n")
