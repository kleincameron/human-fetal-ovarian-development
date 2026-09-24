project_root <- normalizePath(
  Sys.getenv("PROJECT_ROOT", unset = getwd()),
  mustWork = TRUE
)

data_root <- normalizePath(
  Sys.getenv(
    "FETAL_OVARY_DATA_ROOT",
    unset = file.path(
      dirname(project_root),
      paste0(basename(project_root), "_controlled_data"),
      "HRA019091"
    )
  ),
  mustWork = FALSE
)

results_root <- normalizePath(
  Sys.getenv(
    "FETAL_OVARY_RESULTS_ROOT",
    unset = file.path(
      dirname(project_root),
      paste0(basename(project_root), "_results")
    )
  ),
  mustWork = FALSE
)

object_dir <- file.path(results_root, "objects")
figure_dir <- file.path(results_root, "figures")
table_dir <- file.path(results_root, "tables")
log_dir <- file.path(results_root, "logs")

snrna_fastq_manifest <- file.path(project_root, "metadata", "input_manifest_snRNAseq_fastq.csv")
snrna_sample_metadata <- file.path(project_root, "metadata", "sample_metadata_snRNAseq.csv")

snrna_object <- file.path(object_dir, "snRNAseq", "fetal_ovary_snRNAseq_canonical.rds")

dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
