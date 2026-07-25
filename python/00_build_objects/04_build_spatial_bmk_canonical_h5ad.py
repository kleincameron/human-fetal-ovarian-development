#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import os
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def default_results_root(project_root: Path) -> Path:
    return Path(
        os.environ.get(
            "FETAL_OVARY_RESULTS_ROOT",
            str(project_root.parent / "github_code_for_publication_results"),
        )
    ).resolve()


def default_spatial_data_root(project_root: Path) -> Path:
    return Path(
        os.environ.get(
            "FETAL_OVARY_SPATIAL_BMK_ROOT",
            str(project_root / "data" / "processed_spatial"),
        )
    ).resolve()


def parse_bool(value: str) -> bool:
    return str(value).strip().lower() in {"true", "t", "1", "yes", "y"}


def resolve_input_path(path_value: str, data_root: Path) -> Path:
    p = Path(str(path_value))
    if p.is_absolute():
        return p
    return (data_root / p).resolve()


def validate_10x_dir(mtx_dir: Path) -> None:
    if not mtx_dir.exists():
        raise FileNotFoundError(f"Matrix directory does not exist: {mtx_dir}")

    present = {p.name for p in mtx_dir.iterdir()}

    if not ({"matrix.mtx", "matrix.mtx.gz"} & present):
        raise FileNotFoundError(f"matrix.mtx(.gz) not found in {mtx_dir}")

    if not ({"barcodes.tsv", "barcodes.tsv.gz"} & present):
        raise FileNotFoundError(f"barcodes.tsv(.gz) not found in {mtx_dir}")

    if not (
        {"features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz"} & present
    ):
        raise FileNotFoundError(
            f"features.tsv(.gz) or genes.tsv(.gz) not found in {mtx_dir}"
        )
def sparse_fractional_report(x: sparse.spmatrix) -> dict:
    x = x.tocsr()
    data = x.data

    if data.size == 0:
        return {
            "nnz": 0,
            "max_fractional_deviation": 0.0,
            "fraction_non_integer": 0.0,
            "min_value": 0.0,
            "max_value": 0.0,
        }

    rounded = np.rint(data)
    dev = np.abs(data - rounded)

    return {
        "nnz": int(data.size),
        "max_fractional_deviation": float(dev.max()),
        "fraction_non_integer": float(np.mean(dev > 1e-6)),
        "min_value": float(data.min()),
        "max_value": float(data.max()),
    }


def ensure_integer_like_counts(x, tolerance: float = 1e-6) -> sparse.csr_matrix:
    if sparse.issparse(x):
        x = x.tocsr()
        report = sparse_fractional_report(x)
        if report["max_fractional_deviation"] > tolerance:
            raise ValueError(
                "Input matrix is not integer-like. "
                f"max_fractional_deviation={report['max_fractional_deviation']}"
            )
        x.data = np.rint(x.data).astype(np.float32, copy=False)
        x.data[x.data < 0] = 0.0
        return x.astype(np.float32, copy=False)

    arr = np.asarray(x)
    dev = np.abs(arr - np.rint(arr))
    if float(dev.max()) > tolerance:
        raise ValueError(
            "Input matrix is not integer-like. "
            f"max_fractional_deviation={float(dev.max())}"
        )
    arr = np.rint(arr).astype(np.float32)
    arr[arr < 0] = 0.0
    return sparse.csr_matrix(arr)


def read_coords(coords_path: Path) -> pd.DataFrame:
    if not coords_path.exists():
        raise FileNotFoundError(f"Coordinate file does not exist: {coords_path}")

    coords = pd.read_csv(coords_path, sep="\t", header=None)

    if coords.shape[1] < 3:
        raise ValueError(
            f"Coordinate file has {coords.shape[1]} columns; expected at least 3."
        )

    coords = coords.iloc[:, :3].copy()
    coords.columns = ["barcode", "x", "y"]
    coords["barcode"] = coords["barcode"].astype(str)
    coords["x"] = pd.to_numeric(coords["x"], errors="raise")
    coords["y"] = pd.to_numeric(coords["y"], errors="raise")
    coords = coords.drop_duplicates(subset=["barcode"], keep="first")
    coords = coords.set_index("barcode")

    return coords
def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w") as handle:
        json.dump(payload, handle, indent=2)
    tmp.replace(path)


def build_one_sample(
    row: dict,
    data_root: Path,
    out_object_dir: Path,
    out_table_dir: Path,
    min_coordinate_fraction: float,
) -> dict:
    sample_id = row["sample_id"]
    gestational_week = row["gestational_week"]
    platform_name = row.get("platform", "BMK_ST")

    mtx_dir = resolve_input_path(row["mtx_dir"], data_root)
    coords_path = resolve_input_path(row["coords_path"], data_root)

    print("=" * 80)
    print(f"Sample: {sample_id}")
    print(f"  gestational_week: {gestational_week}")
    print(f"  platform: {platform_name}")
    print(f"  matrix directory: {mtx_dir.name}")
    print(f"  coordinate file: {coords_path.name}")

    validate_10x_dir(mtx_dir)

    adata = sc.read_10x_mtx(
        str(mtx_dir),
        var_names="gene_symbols",
        cache=False,
    )

    adata.var_names = adata.var_names.astype(str)
    adata.obs_names = adata.obs_names.astype(str)

    n_locations_input = int(adata.n_obs)
    n_genes_input = int(adata.n_vars)

    count_report_before = (
        sparse_fractional_report(adata.X)
        if sparse.issparse(adata.X)
        else {
            "nnz": int(np.count_nonzero(adata.X)),
            "max_fractional_deviation": float(
                np.abs(np.asarray(adata.X) - np.rint(np.asarray(adata.X))).max()
            ),
            "fraction_non_integer": float(
                np.mean(np.abs(np.asarray(adata.X) - np.rint(np.asarray(adata.X))) > 1e-6)
            ),
            "min_value": float(np.asarray(adata.X).min()),
            "max_value": float(np.asarray(adata.X).max()),
        }
    )

    adata.X = ensure_integer_like_counts(adata.X)

    libsize = np.asarray(adata.X.sum(axis=1)).ravel().astype(np.float64)
    keep_nonzero = libsize > 0

    n_zero_locations = int(np.sum(~keep_nonzero))
    print(f"  input locations: {n_locations_input}")
    print(f"  zero-count locations removed: {n_zero_locations}")

    adata = adata[keep_nonzero, :].copy()
    libsize = libsize[keep_nonzero]

    coords = read_coords(coords_path)
    common = adata.obs_names.intersection(coords.index)
    coordinate_fraction = len(common) / max(adata.n_obs, 1)

    if coordinate_fraction < min_coordinate_fraction:
        raise ValueError(
            f"Only {len(common)} / {adata.n_obs} locations have coordinates "
            f"({coordinate_fraction:.3f}); expected at least {min_coordinate_fraction:.3f}."
        )

    if len(common) < adata.n_obs:
        print(
            f"  subsetting to locations with coordinates: "
            f"{len(common)} / {adata.n_obs}"
        )
        adata = adata[common, :].copy()

    coords = coords.loc[adata.obs_names]

    adata.obs["sample_id"] = sample_id
    adata.obs["gestational_week"] = str(gestational_week)
    adata.obs["platform"] = platform_name
    adata.obs["source_type"] = "vendor_processed_spatial_matrix"
    adata.obs["barcode"] = adata.obs_names.astype(str)
    adata.obs["x"] = coords["x"].to_numpy(dtype=np.float32)
    adata.obs["y"] = coords["y"].to_numpy(dtype=np.float32)

    adata.obsm["spatial"] = coords[["x", "y"]].to_numpy(dtype=np.float32)

    adata.uns["spatial_input"] = {
        "sample_id": sample_id,
        "gestational_week": str(gestational_week),
        "platform": platform_name,
        "input_type": "vendor_processed_10x_style_matrix_plus_barcode_coordinates",
        "source_matrix_directory_name": mtx_dir.name,
        "source_coordinate_file_name": coords_path.name,
        "note": "Input paths are resolved from a user-provided manifest and are not stored in the canonical object.",
    }
    output_h5ad = out_object_dir / f"{sample_id}.spatial_bmk_canonical.h5ad"
    output_h5ad.parent.mkdir(parents=True, exist_ok=True)

    print(f"  writing: {output_h5ad}")
    adata.write_h5ad(output_h5ad)

    summary = {
        "sample_id": sample_id,
        "gestational_week": gestational_week,
        "platform": platform_name,
        "source_matrix_directory_name": mtx_dir.name,
        "source_coordinate_file_name": coords_path.name,
        "output_h5ad_name": output_h5ad.name,
        "n_locations_input": n_locations_input,
        "n_genes_input": n_genes_input,
        "n_zero_count_locations_removed": n_zero_locations,
        "n_locations_final": int(adata.n_obs),
        "n_genes_final": int(adata.n_vars),
        "coordinate_fraction_after_zero_filter": float(coordinate_fraction),
        "nnz_final": int(adata.X.nnz),
        "min_total_counts": float(libsize[libsize > 0].min()) if np.any(libsize > 0) else 0.0,
        "median_total_counts": float(np.median(libsize)) if libsize.size else 0.0,
        "max_total_counts": float(libsize.max()) if libsize.size else 0.0,
        "count_max_fractional_deviation_input": count_report_before["max_fractional_deviation"],
        "count_fraction_non_integer_input": count_report_before["fraction_non_integer"],
        "count_min_value_input": count_report_before["min_value"],
        "count_max_value_input": count_report_before["max_value"],
    }

    write_json(
        out_table_dir / f"{sample_id}.spatial_bmk_canonical_summary.json",
        summary,
    )

    return summary


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build canonical BMK-ST spatial AnnData objects from processed matrix and coordinate files."
    )
    parser.add_argument(
        "--manifest",
        default=None,
        help="CSV manifest. Defaults to metadata/spatial_bmk_samples.csv.",
    )
    parser.add_argument(
        "--results-root",
        default=None,
        help="Output root. Defaults to FETAL_OVARY_RESULTS_ROOT or ../github_code_for_publication_results.",
    )
    parser.add_argument(
        "--spatial-data-root",
        default=None,
        help="Root used to resolve relative paths in the manifest. Defaults to FETAL_OVARY_SPATIAL_BMK_ROOT or data/processed_spatial.",
    )
    parser.add_argument(
        "--min-coordinate-fraction",
        type=float,
        default=0.80,
        help="Minimum fraction of nonzero locations that must have coordinates.",
    )

    args = parser.parse_args()
    project_root = repo_root()
    manifest = Path(args.manifest).resolve() if args.manifest else project_root / "metadata" / "spatial_bmk_samples.csv"
    results_root = Path(args.results_root).resolve() if args.results_root else default_results_root(project_root)
    spatial_data_root = Path(args.spatial_data_root).resolve() if args.spatial_data_root else default_spatial_data_root(project_root)

    out_root = results_root / "spatial_bmk_canonical_objects"
    out_object_dir = out_root / "objects"
    out_table_dir = out_root / "tables"
    out_log_dir = out_root / "logs"

    out_object_dir.mkdir(parents=True, exist_ok=True)
    out_table_dir.mkdir(parents=True, exist_ok=True)
    out_log_dir.mkdir(parents=True, exist_ok=True)

    if not manifest.exists():
        raise FileNotFoundError(f"Manifest not found: {manifest}")

    print("Building canonical BMK-ST spatial objects")
    print(f"  project_root: {project_root}")
    print(f"  manifest file: {manifest.name}")
    print("  input paths are resolved from the user-provided manifest")
    print(f"  results_root: {results_root}")

    rows = []
    with manifest.open() as handle:
        reader = csv.DictReader(handle)
        required = {"sample_id", "gestational_week", "platform", "mtx_dir", "coords_path", "include"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Manifest is missing columns: {sorted(missing)}")

        for row in reader:
            if parse_bool(row["include"]):
                rows.append(row)

    if not rows:
        raise ValueError("No manifest rows with include=TRUE.")

    summaries = []
    for row in rows:
        summaries.append(
            build_one_sample(
                row=row,
                data_root=spatial_data_root,
                out_object_dir=out_object_dir,
                out_table_dir=out_table_dir,
                min_coordinate_fraction=args.min_coordinate_fraction,
            )
        )
    summary_df = pd.DataFrame(summaries)
    summary_csv = out_table_dir / "spatial_bmk_canonical_summary.csv"
    summary_df.to_csv(summary_csv, index=False)

    package_versions = {
        "python": sys.version,
        "platform": platform.platform(),
        "scanpy": sc.__version__,
        "numpy": np.__version__,
        "pandas": pd.__version__,
    }

    write_json(out_log_dir / "spatial_bmk_canonical_package_versions.json", package_versions)

    with (out_log_dir / "spatial_bmk_canonical_summary.txt").open("w") as handle:
        handle.write(f"Manifest file: {manifest.name}\n")
        handle.write("Input paths are resolved from the user-provided manifest.\n")
        handle.write("Local input paths are not stored in canonical H5AD metadata or public summary tables.\n")
        handle.write(f"Samples: {', '.join(summary_df['sample_id'].astype(str))}\n")
        handle.write(f"Summary table name: {summary_csv.name}\n")

    print("\nDone.")
    print(f"Summary table: {summary_csv}")
    print(f"Output directory: {out_root}")


if __name__ == "__main__":
    main()
