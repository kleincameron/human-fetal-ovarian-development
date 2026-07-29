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

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


SAMPLE_LABELS = {
    "24w_bmk_cellbin": "24w",
    "25w_bmk_cellbin": "25w",
    "26p4w_bmk_cellbin": "26.5w",
}

SAMPLE_ORDER = ["24w_bmk_cellbin", "25w_bmk_cellbin", "26p4w_bmk_cellbin"]


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def default_results_root(project_root: Path) -> Path:
    return Path(
        os.environ.get(
            "FETAL_OVARY_RESULTS_ROOT",
            str(project_root.parent / "github_code_for_publication_results"),
        )
    ).resolve()


def parse_bool(value: str) -> bool:
    return str(value).strip().lower() in {"true", "t", "1", "yes", "y"}


def setup_plot_style() -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 8,
        "axes.titlesize": 8,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.dpi": 600,
        "figure.dpi": 120,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
    })
def safe_write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w") as handle:
        json.dump(payload, handle, indent=2)
    tmp.replace(path)


def save_figure(fig: plt.Figure, out_prefix: Path) -> None:
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_prefix.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out_prefix.with_suffix(".png"), bbox_inches="tight", dpi=600)
    plt.close(fig)


def load_manifest(manifest: Path) -> list[dict]:
    if not manifest.exists():
        raise FileNotFoundError(f"Manifest not found: {manifest}")

    rows = []

    with manifest.open() as handle:
        reader = csv.DictReader(handle)

        required = {"sample_id", "gestational_week", "platform", "include"}
        missing = required - set(reader.fieldnames or [])

        if missing:
            raise ValueError(f"Manifest is missing columns: {sorted(missing)}")

        for row in reader:
            if parse_bool(row["include"]):
                rows.append(row)

    if not rows:
        raise ValueError("No manifest rows with include=TRUE.")

    return rows


def quantiles(values: np.ndarray, prefix: str) -> dict:
    values = np.asarray(values, dtype=np.float64)

    if values.size == 0:
        return {
            f"{prefix}_min": 0.0,
            f"{prefix}_q25": 0.0,
            f"{prefix}_median": 0.0,
            f"{prefix}_q75": 0.0,
            f"{prefix}_q90": 0.0,
            f"{prefix}_q99": 0.0,
            f"{prefix}_max": 0.0,
        }

    q = np.quantile(values, [0, 0.25, 0.5, 0.75, 0.9, 0.99, 1.0])

    return {
        f"{prefix}_min": float(q[0]),
        f"{prefix}_q25": float(q[1]),
        f"{prefix}_median": float(q[2]),
        f"{prefix}_q75": float(q[3]),
        f"{prefix}_q90": float(q[4]),
        f"{prefix}_q99": float(q[5]),
        f"{prefix}_max": float(q[6]),
    }
def gene_category_masks(var_names: pd.Index) -> dict[str, np.ndarray]:
    genes = pd.Index(var_names.astype(str))

    mt = np.asarray(genes.str.startswith("MT-"), dtype=bool)
    ribo = (
        np.asarray(genes.str.startswith("RPL"), dtype=bool) |
        np.asarray(genes.str.startswith("RPS"), dtype=bool) |
        np.asarray(genes.str.startswith("MRPL"), dtype=bool) |
        np.asarray(genes.str.startswith("MRPS"), dtype=bool)
    )

    return {
        "mt": mt,
        "ribo": ribo,
    }


def compute_location_qc(adata: sc.AnnData) -> pd.DataFrame:
    x = adata.X.tocsr() if sparse.issparse(adata.X) else sparse.csr_matrix(np.asarray(adata.X))
    masks = gene_category_masks(pd.Index(adata.var_names))

    total_counts = np.asarray(x.sum(axis=1)).ravel().astype(np.float64)
    detected_genes = np.asarray((x > 0).sum(axis=1)).ravel().astype(np.float64)

    mt_counts = (
        np.asarray(x[:, masks["mt"]].sum(axis=1)).ravel().astype(np.float64)
        if np.any(masks["mt"])
        else np.zeros(adata.n_obs)
    )
    ribo_counts = (
        np.asarray(x[:, masks["ribo"]].sum(axis=1)).ravel().astype(np.float64)
        if np.any(masks["ribo"])
        else np.zeros(adata.n_obs)
    )

    total_safe = np.where(total_counts > 0, total_counts, np.nan)
    pct_mt = np.nan_to_num(100.0 * mt_counts / total_safe, nan=0.0, posinf=0.0, neginf=0.0)
    pct_ribo = np.nan_to_num(100.0 * ribo_counts / total_safe, nan=0.0, posinf=0.0, neginf=0.0)

    obs = adata.obs.copy()

    if "x" in obs.columns and "y" in obs.columns:
        x_coord = obs["x"].to_numpy(dtype=np.float64)
        y_coord = obs["y"].to_numpy(dtype=np.float64)
    elif "spatial" in adata.obsm:
        x_coord = adata.obsm["spatial"][:, 0].astype(np.float64)
        y_coord = adata.obsm["spatial"][:, 1].astype(np.float64)
    else:
        x_coord = np.full(adata.n_obs, np.nan)
        y_coord = np.full(adata.n_obs, np.nan)

    sample_id = str(obs["sample_id"].iloc[0]) if "sample_id" in obs.columns and len(obs) else ""
    gestational_week = str(obs["gestational_week"].iloc[0]) if "gestational_week" in obs.columns and len(obs) else ""
    platform_name = str(obs["platform"].iloc[0]) if "platform" in obs.columns and len(obs) else ""
    barcode = obs["barcode"].astype(str).values if "barcode" in obs.columns else adata.obs_names.astype(str)

    return pd.DataFrame({
        "sample_id": sample_id,
        "gestational_week": gestational_week,
        "platform": platform_name,
        "barcode": barcode,
        "x": x_coord,
        "y": y_coord,
        "total_counts": total_counts,
        "detected_genes": detected_genes,
        "mt_counts": mt_counts,
        "pct_mt": pct_mt,
        "ribo_counts": ribo_counts,
        "pct_ribo": pct_ribo,
    })
def summarize_sample(
    sample_id: str,
    gestational_week: str,
    platform_name: str,
    adata: sc.AnnData,
    location_qc: pd.DataFrame,
) -> dict:
    x = adata.X.tocsr() if sparse.issparse(adata.X) else sparse.csr_matrix(np.asarray(adata.X))
    masks = gene_category_masks(pd.Index(adata.var_names))

    coordinate_complete = (
        location_qc["x"].notna().all() and
        location_qc["y"].notna().all() and
        "spatial" in adata.obsm and
        adata.obsm["spatial"].shape[0] == adata.n_obs
    )

    summary = {
        "sample_id": sample_id,
        "gestational_week": gestational_week,
        "platform": platform_name,
        "canonical_h5ad_name": f"{sample_id}.spatial_bmk_canonical.h5ad",
        "n_locations": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "nnz": int(x.nnz),
        "matrix_density": float(x.nnz / (adata.n_obs * adata.n_vars)),
        "n_mt_genes": int(np.sum(masks["mt"])),
        "n_ribo_genes": int(np.sum(masks["ribo"])),
        "has_spatial_obsm": bool("spatial" in adata.obsm),
        "coordinate_complete": bool(coordinate_complete),
        "x_min": float(np.nanmin(location_qc["x"])),
        "x_max": float(np.nanmax(location_qc["x"])),
        "y_min": float(np.nanmin(location_qc["y"])),
        "y_max": float(np.nanmax(location_qc["y"])),
    }

    for col in ["total_counts", "detected_genes", "pct_mt", "pct_ribo"]:
        summary.update(quantiles(location_qc[col].to_numpy(), col))

    return summary


def plot_metric_distributions(location_df: pd.DataFrame, figure_dir: Path) -> None:
    metrics = [
        ("total_counts", "log10(total counts + 1)", True),
        ("detected_genes", "Detected genes", False),
        ("pct_mt", "Mitochondrial reads (%)", False),
        ("pct_ribo", "Ribosomal reads (%)", False),
    ]

    fig, axes = plt.subplots(1, 4, figsize=(7.4, 2.15), constrained_layout=True)

    for ax, (col, ylabel, log_transform) in zip(axes, metrics):
        data = []
        labels = []

        for sample_id in SAMPLE_ORDER:
            vals = location_df.loc[location_df["sample_id"] == sample_id, col].astype(float).to_numpy()

            if log_transform:
                vals = np.log10(vals + 1.0)

            data.append(vals)
            labels.append(SAMPLE_LABELS.get(sample_id, sample_id))

        parts = ax.violinplot(
            data,
            showmeans=False,
            showmedians=True,
            showextrema=False,
        )
        for body in parts["bodies"]:
            body.set_facecolor("#BDBDBD")
            body.set_edgecolor("black")
            body.set_alpha(1.0)
            body.set_linewidth(0.4)

        if "cmedians" in parts:
            parts["cmedians"].set_color("black")
            parts["cmedians"].set_linewidth(0.8)

        ax.set_ylabel(ylabel)
        ax.set_xticks(range(1, len(labels) + 1))
        ax.set_xticklabels(labels)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    save_figure(fig, figure_dir / "spatial_bmk_qc_metric_distributions")


def plot_spatial_maps(location_qc: pd.DataFrame, sample_id: str, figure_dir: Path, max_points: int | None) -> None:
    sample = location_qc.loc[location_qc["sample_id"] == sample_id].copy()

    if max_points is not None and len(sample) > max_points:
        sample = sample.sample(n=max_points, random_state=1)

    metrics = [
        ("total_counts", "log10 counts", True),
        ("detected_genes", "Detected genes", False),
        ("pct_mt", "Mitochondrial %", False),
        ("pct_ribo", "Ribosomal %", False),
    ]

    fig, axes = plt.subplots(1, 4, figsize=(7.4, 2.25), constrained_layout=True)

    for ax, (col, title, log_transform) in zip(axes, metrics):
        values = sample[col].astype(float).to_numpy()

        if log_transform:
            values = np.log10(values + 1.0)

        lo, hi = np.quantile(values, [0.01, 0.99])

        scatter = ax.scatter(
            sample["x"],
            sample["y"],
            c=values,
            s=0.15,
            linewidths=0,
            cmap="viridis",
            vmin=lo,
            vmax=hi,
            rasterized=True,
        )

        ax.set_title(title)
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_xticks([])
        ax.set_yticks([])

        for spine in ax.spines.values():
            spine.set_visible(False)

        cbar = fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.02)
        cbar.ax.tick_params(labelsize=6, width=0.4)

    fig.suptitle(SAMPLE_LABELS.get(sample_id, sample_id), y=1.02, fontsize=8)

    save_figure(fig, figure_dir / f"{sample_id}.spatial_bmk_qc_maps")


def write_index(figure_dir: Path) -> None:
    lines = [
        "BMK-ST spatial sequencing QC figures",
        "",
        "These figures use canonical BMK-ST spatial H5AD objects.",
        "",
        "Generated figures:",
        "- spatial_bmk_qc_metric_distributions.pdf/png",
        "- <sample>.spatial_bmk_qc_maps.pdf/png",
    ]
    (figure_dir / "spatial_bmk_qc_figure_index.txt").write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="QC for canonical BMK-ST spatial sequencing objects."
    )
    parser.add_argument("--manifest", default=None, help="CSV manifest. Defaults to metadata/spatial_bmk_samples.csv.")
    parser.add_argument("--results-root", default=None, help="Results root.")
    parser.add_argument(
        "--max-points-per-sample",
        type=int,
        default=None,
        help="Optional downsampling for spatial maps. Default uses all locations.",
    )

    args = parser.parse_args()

    project_root = repo_root()
    results_root = Path(args.results_root).resolve() if args.results_root else default_results_root(project_root)
    manifest = Path(args.manifest).resolve() if args.manifest else project_root / "metadata" / "spatial_bmk_samples.csv"

    canonical_dir = results_root / "spatial_bmk_canonical_objects" / "objects"
    out_root = results_root / "spatial_bmk_canonical_qc"
    table_dir = out_root / "tables"
    figure_dir = out_root / "figures"
    log_dir = out_root / "logs"

    table_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    setup_plot_style()

    print("Running BMK-ST spatial sequencing QC")
    print(f"  manifest file: {manifest.name}")
    print(f"  canonical object directory: {canonical_dir}")
    print(f"  output directory: {out_root}")

    rows = load_manifest(manifest)

    all_location_qc = []
    summaries = []

    for row in rows:
        sample_id = row["sample_id"]
        gestational_week = row.get("gestational_week", "")
        platform_name = row.get("platform", "")

        h5ad_path = canonical_dir / f"{sample_id}.spatial_bmk_canonical.h5ad"

        if not h5ad_path.exists():
            raise FileNotFoundError(f"Canonical H5AD not found: {h5ad_path}")

        print("=" * 80)
        print(f"Sample: {sample_id}")
        print(f"Canonical H5AD: {h5ad_path.name}")

        adata = sc.read_h5ad(h5ad_path)
        adata.X = adata.X.tocsr() if sparse.issparse(adata.X) else sparse.csr_matrix(np.asarray(adata.X))

        location_qc = compute_location_qc(adata)
        location_qc.to_csv(table_dir / f"{sample_id}.spatial_bmk_location_qc.csv", index=False)

        summary = summarize_sample(
            sample_id=sample_id,
            gestational_week=gestational_week,
            platform_name=platform_name,
            adata=adata,
            location_qc=location_qc,
        )
        summaries.append(summary)
        all_location_qc.append(location_qc)

        print(f"  locations: {summary['n_locations']}")
        print(f"  genes: {summary['n_genes']}")
        print(f"  median counts: {summary['total_counts_median']:.3f}")
        print(f"  median detected genes: {summary['detected_genes_median']:.3f}")
        print(f"  median mitochondrial %: {summary['pct_mt_median']:.3f}")
        print(f"  median ribosomal %: {summary['pct_ribo_median']:.3f}")

    summary_df = pd.DataFrame(summaries)
    summary_df.to_csv(table_dir / "spatial_bmk_qc_summary.csv", index=False)

    location_df = pd.concat(all_location_qc, ignore_index=True)

    plot_metric_distributions(location_df, figure_dir)

    for sample_id in SAMPLE_ORDER:
        if sample_id in set(location_df["sample_id"]):
            plot_spatial_maps(
                location_qc=location_df,
                sample_id=sample_id,
                figure_dir=figure_dir,
                max_points=args.max_points_per_sample,
            )

    write_index(figure_dir)

    package_versions = {
        "python": sys.version,
        "platform": platform.platform(),
        "scanpy": sc.__version__,
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "matplotlib": matplotlib.__version__,
    }

    safe_write_json(log_dir / "spatial_bmk_qc_package_versions.json", package_versions)

    with (log_dir / "spatial_bmk_qc_summary.txt").open("w") as handle:
        handle.write("BMK-ST spatial sequencing QC\n")
        handle.write("Input: canonical BMK-ST spatial H5AD objects\n")
        handle.write(f"Manifest file: {manifest.name}\n")
        handle.write(f"Samples: {', '.join(summary_df['sample_id'].astype(str))}\n")
        handle.write("Summary table: spatial_bmk_qc_summary.csv\n")
        handle.write("Location-level tables: <sample>.spatial_bmk_location_qc.csv\n")

    print("\nDone.")
    print(f"QC tables: {table_dir}")
    print(f"QC figures: {figure_dir}")


if __name__ == "__main__":
    main()
