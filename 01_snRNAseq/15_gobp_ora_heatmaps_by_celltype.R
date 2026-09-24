#!/usr/bin/env Rscript

required_packages <- c(
  "AnnotationDbi",
  "org.Hs.eg.db",
  "GO.db",
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
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(GO.db)
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

deg_results_root <- if (exists("snrna_deg_results_root", inherits = FALSE)) {
  snrna_deg_results_root
} else {
  file.path(results_base, "snRNAseq_DEGs_celltype_subcluster")
}

out_root <- if (exists("snrna_gobp_ora_results_root", inherits = FALSE)) {
  snrna_gobp_ora_results_root
} else {
  file.path(results_base, "snRNAseq_GOBP_ORA_by_celltype")
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

deg_padj_cut <- as.numeric(Sys.getenv("GOBP_ORA_DEG_PADJ_CUT", unset = "0.05"))
deg_logfc_min <- as.numeric(Sys.getenv("GOBP_ORA_DEG_LOGFC_MIN", unset = "0"))
deg_pct1_min <- as.numeric(Sys.getenv("GOBP_ORA_DEG_PCT1_MIN", unset = "0.25"))
universe_pct_min <- as.numeric(Sys.getenv("GOBP_ORA_UNIVERSE_PCT_MIN", unset = "0.10"))

ora_padj_cut <- as.numeric(Sys.getenv("GOBP_ORA_PADJ_CUT", unset = "0.05"))
min_go_size <- as.integer(Sys.getenv("GOBP_ORA_MIN_GO_SIZE", unset = "10"))
max_go_size <- as.integer(Sys.getenv("GOBP_ORA_MAX_GO_SIZE", unset = "500"))
min_overlap <- as.integer(Sys.getenv("GOBP_ORA_MIN_OVERLAP", unset = "2"))

top_terms_per_subcluster <- as.integer(Sys.getenv("GOBP_ORA_TOP_TERMS_PER_SUBCLUSTER", unset = "6"))
max_terms_per_celltype <- as.integer(Sys.getenv("GOBP_ORA_MAX_TERMS_PER_CELLTYPE", unset = "28"))
max_significant_fraction <- as.numeric(Sys.getenv("GOBP_ORA_MAX_SIGNIFICANT_FRACTION", unset = "0.70"))
min_specificity_range <- as.numeric(Sys.getenv("GOBP_ORA_MIN_SPECIFICITY_RANGE", unset = "1.00"))

heatmap_value_cap <- as.numeric(Sys.getenv("GOBP_ORA_HEATMAP_VALUE_CAP", unset = "4"))
heatmap_width <- as.numeric(Sys.getenv("GOBP_ORA_HEATMAP_WIDTH", unset = "10"))
heatmap_height <- as.numeric(Sys.getenv("GOBP_ORA_HEATMAP_HEIGHT", unset = "8"))
heatmap_dpi <- as.integer(Sys.getenv("GOBP_ORA_HEATMAP_DPI", unset = "600"))
cluster_columns <- tolower(Sys.getenv("GOBP_ORA_CLUSTER_COLUMNS", unset = "true")) %in% c("true", "t", "1", "yes", "y")

available_cell_types <- basename(list.dirs(subcluster_deg_root, recursive = FALSE))

cell_types_env <- Sys.getenv("GOBP_ORA_CELL_TYPES", unset = "")
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
    gene = as.character(df[[gene_col]]),
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

get_gobp_gene_map <- function(universe_genes) {
  go_raw <- AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = universe_genes,
    keytype = "SYMBOL",
    columns = c("SYMBOL", "GO", "ONTOLOGY")
  ) |>
    as_tibble()

  required_cols <- c("SYMBOL", "GO", "ONTOLOGY")
  missing_cols <- setdiff(required_cols, colnames(go_raw))
  if (length(missing_cols) > 0) {
    stop("org.Hs.eg.db GO mapping did not return expected columns: ", paste(missing_cols, collapse = ", "))
  }

  go_map <- tibble(
    gene = as.character(go_raw[["SYMBOL"]]),
    go_id = as.character(go_raw[["GO"]]),
    ontology = as.character(go_raw[["ONTOLOGY"]])
  ) |>
    filter(!is.na(gene), !is.na(go_id), ontology == "BP") |>
    distinct(gene, go_id)

  if (nrow(go_map) == 0) return(tibble())

  term_sizes <- go_map |>
    count(go_id, name = "term_size") |>
    filter(term_size >= min_go_size, term_size <= max_go_size)

  go_map |>
    semi_join(term_sizes, by = "go_id") |>
    left_join(term_sizes, by = "go_id")
}

get_go_term_names <- function(go_ids) {
  vapply(go_ids, function(go_id) {
    term <- GO.db::GOTERM[[go_id]]
    if (is.null(term)) NA_character_ else AnnotationDbi::Term(term)
  }, character(1))
}

format_go_label <- function(term_name) {
  x <- toupper(gsub("[^A-Za-z0-9]+", "_", term_name))
  gsub("^_+|_+$", "", x)
}

run_gobp_ora <- function(deg_tbl, cell_type) {
  subclusters <- order_subclusters(cell_type, unique(deg_tbl$subcluster))
  universe_genes <- sort(unique(deg_tbl$gene))

  gene_go <- get_gobp_gene_map(universe_genes)
  if (nrow(gene_go) == 0) {
    warning("No GOBP mappings found for cell type: ", cell_type)
    return(tibble())
  }

  term_genes_by_go <- split(gene_go$gene, gene_go$go_id)
  term_genes_by_go <- lapply(term_genes_by_go, unique)
  term_ids <- names(term_genes_by_go)
  term_names <- get_go_term_names(term_ids)

  records <- vector("list", length(subclusters) * length(term_ids))
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

    for (go_id in term_ids) {
      term_genes <- intersect(term_genes_by_go[[go_id]], universe_i)

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
        go_id = go_id,
        term_name = unname(term_names[[go_id]]),
        term_label = format_go_label(unname(term_names[[go_id]])),
        universe_n = universe_n,
        candidate_gene_n = candidate_n,
        term_size_in_universe = length(term_genes),
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
    arrange(cell_type, subcluster, p_adj, desc(log2_odds_ratio), go_id)
}

select_heatmap_terms <- function(ora_tbl) {
  n_subclusters <- length(unique(ora_tbl$subcluster))

  term_stats <- ora_tbl |>
    group_by(go_id, term_name, term_label) |>
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
    semi_join(term_stats, by = c("go_id", "term_name", "term_label")) |>
    filter(
      p_adj <= ora_padj_cut,
      log2_odds_ratio > 0,
      overlap_n >= min_overlap
    ) |>
    group_by(subcluster) |>
    arrange(p_adj, desc(log2_odds_ratio), desc(overlap_n), .by_group = TRUE) |>
    slice_head(n = top_terms_per_subcluster) |>
    ungroup() |>
    distinct(go_id, term_name, term_label)

  term_stats |>
    semi_join(selected_by_subcluster, by = c("go_id", "term_name", "term_label")) |>
    arrange(best_p_adj, desc(range_log2_odds_ratio), desc(max_log2_odds_ratio), term_label) |>
    slice_head(n = max_terms_per_celltype)
}

make_heatmap_matrix <- function(ora_tbl, selected_terms, subclusters) {
  term_order <- selected_terms$go_id
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
    go_id <- ora_tbl$go_id[[i]]
    subcluster_i <- ora_tbl$subcluster[[i]]

    if (go_id %in% term_order && subcluster_i %in% subclusters) {
      value <- ora_tbl$log2_odds_ratio[[i]]
      mat[go_id, subcluster_i] <- value

      if (
        !is.na(ora_tbl$p_adj[[i]]) &&
          ora_tbl$p_adj[[i]] <= ora_padj_cut &&
          value > 0 &&
          ora_tbl$overlap_n[[i]] >= min_overlap
      ) {
        star_mat[go_id, subcluster_i] <- "*"
      }
    }
  }

  mat <- pmax(pmin(mat, heatmap_value_cap), -heatmap_value_cap)

  row_labels <- selected_terms$term_label[match(rownames(mat), selected_terms$go_id)]
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
    "GOBP ORA - ", cell_type,
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
    fontsize_row = max(4, publication_base_size - 2),
    fontsize_col = max(5, publication_base_size - 1),
    angle_col = 90,
    border_color = NA,
    main = title_text,
    silent = TRUE
  )

  pdf_out <- file.path(out_figure_dir, paste0("GOBP_ORA_heatmap_", cell_type, ".pdf"))
  png_out <- file.path(out_figure_dir, paste0("GOBP_ORA_heatmap_", cell_type, ".png"))

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
      rownames_to_column("term_label"),
    file.path(out_table_dir, paste0("GOBP_ORA_heatmap_matrix_", cell_type, ".csv"))
  )
}

all_ora <- list()
all_selected <- list()
summary_rows <- list()

for (cell_type in target_cell_types) {
  cat("\n============================================================\n")
  cat("Running GOBP ORA for cell type: ", cell_type, "\n", sep = "")
  cat("============================================================\n")

  deg_tbl <- read_celltype_degs(cell_type)

  if (nrow(deg_tbl) == 0) {
    warning("No DEG rows available for cell type: ", cell_type)
    next
  }

  ora_tbl <- run_gobp_ora(deg_tbl, cell_type)

  if (nrow(ora_tbl) == 0) {
    warning("No ORA results available for cell type: ", cell_type)
    next
  }

  selected_terms <- select_heatmap_terms(ora_tbl)

  readr::write_csv(
    ora_tbl,
    file.path(out_table_dir, paste0("GOBP_ORA_all_terms_", cell_type, ".csv"))
  )

  readr::write_csv(
    selected_terms,
    file.path(out_table_dir, paste0("GOBP_ORA_selected_heatmap_terms_", cell_type, ".csv"))
  )

  if (nrow(selected_terms) > 0) {
    plot_celltype_heatmap(ora_tbl, selected_terms, cell_type)
  } else {
    warning("No specific significant GOBP terms selected for heatmap: ", cell_type)
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
    file.path(out_table_dir, "COMBINED_GOBP_ORA_all_terms_by_celltype.csv")
  )
}

if (length(all_selected) > 0) {
  readr::write_csv(
    bind_rows(all_selected),
    file.path(out_table_dir, "COMBINED_GOBP_ORA_selected_heatmap_terms_by_celltype.csv")
  )
}

if (length(summary_rows) > 0) {
  readr::write_csv(
    bind_rows(summary_rows),
    file.path(out_table_dir, "GOBP_ORA_summary_by_celltype.csv")
  )
}

writeLines(
  c(
    paste("Input DEG root:", deg_results_root),
    paste("Output root:", out_root),
    paste("Cell types:", paste(target_cell_types, collapse = ", ")),
    paste("Excluded cell types:", "endothelial"),
    paste("DEG adjusted p-value cutoff:", deg_padj_cut),
    paste("DEG log2FC minimum:", deg_logfc_min),
    paste("DEG pct.1 minimum:", deg_pct1_min),
    paste("Universe pct minimum:", universe_pct_min),
    paste("ORA adjusted p-value cutoff:", ora_padj_cut),
    paste("GO term size range:", paste(min_go_size, max_go_size, sep = "-")),
    paste("Minimum overlap:", min_overlap),
    paste("Maximum significant subcluster fraction:", max_significant_fraction),
    paste("Minimum specificity range:", min_specificity_range),
    paste("Top terms per subcluster:", top_terms_per_subcluster),
    paste("Maximum heatmap terms per cell type:", max_terms_per_celltype),
    paste("Heatmap value cap:", heatmap_value_cap),
    "ORA uses Fisher's exact test for positive subcluster DEGs against each subcluster's tested gene universe.",
    "Heatmap rows are selected from significant, subcluster-specific GOBP terms.",
    "Heatmap columns use final publication annotation labels."
  ),
  file.path(out_log_dir, "GOBP_ORA_heatmaps_by_celltype_notes.txt")
)

sink(file.path(out_log_dir, "sessionInfo_GOBP_ORA_heatmaps_by_celltype.txt"))
print(sessionInfo())
sink()

cat("\nDone. GOBP ORA outputs written to:\n", out_root, "\n")
