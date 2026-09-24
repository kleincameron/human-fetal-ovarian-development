#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import os
import platform
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


SAMPLE_LABELS = {
    "24w_bmk_cellbin": "24w",
    "25w_bmk_cellbin": "25w",
    "26p4w_bmk_cellbin": "26.5w",
}

SAMPLE_ORDER = ["24w_bmk_cellbin", "25w_bmk_cellbin", "26p4w_bmk_cellbin"]

# This intentionally matches the old interactive cell-type logic.
# Do not add "degenerated" here if the goal is to reproduce that interactive map.
MAJOR_ORDER = [
    "germ",
    "granulosa",
    "stroma",
    "endothelial",
    "mural",
    "immune",
    "erythroid",
    "Unassigned",
]

PLOT_CELL_TYPES = [
    "germ",
    "granulosa",
    "stroma",
    "endothelial",
    "mural",
    "immune",
    "erythroid",
]

MAJOR_COLORS = {
    "germ": "#1B9E77",
    "granulosa": "#D95F02",
    "stroma": "#7570B3",
    "endothelial": "#E7298A",
    "mural": "#66A61E",
    "immune": "#E6AB02",
    "erythroid": "#A6761D",
    "Unassigned": "#D9D9D9",
}

LINEAGE_RELATIVE_THRESHOLD = {
    "germ": 0.40,
    "granulosa": 0.30,
    "stroma": 0.10,
    "endothelial": 0.10,
    "mural": 0.10,
    "immune": 0.10,
    "erythroid": 0.10,
}

BACKGROUND_COLOR = "#BFBFBF"
BACKGROUND_ALPHA = 0.36
FOREGROUND_ALPHA = 0.95
BACKGROUND_POINT_SIZE = 0.22
FOREGROUND_POINT_SIZE = 0.26
FIGURE_DPI = 600


def repo_root() -> Path:
    configured_root = os.environ.get("PROJECT_ROOT")
    if configured_root:
        return Path(configured_root).expanduser().resolve()

    return Path(__file__).resolve().parents[1]


def default_results_root(project_root: Path) -> Path:
    return Path(
        os.environ.get(
            "FETAL_OVARY_RESULTS_ROOT",
            str(project_root.parent / f"{project_root.name}_results"),
        )
    ).expanduser().resolve()
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
        "savefig.dpi": FIGURE_DPI,
        "figure.dpi": 120,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
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
    fig.savefig(out_prefix.with_suffix(".png"), bbox_inches="tight", dpi=FIGURE_DPI)
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


def strip_factor_prefix(name: str) -> str:
    s = str(name)

    if "mu_fg_" in s:
        s = s.split("mu_fg_", 1)[1]

    m = re.search(r"([A-Za-z]+_\d+)$", s)
    if m:
        return m.group(1)

    return s


def relabel_for_interactive_compatibility(fine: str) -> str:
    fine = str(fine)
    # Preserve original interactive behavior:
    # germ_0 and germ_6 were relabeled internally as quiescent states,
    # then counted as germ at the coarse cell-type level.
    if fine == "germ_0":
        return "quiescent_0"
    if fine == "germ_6":
        return "quiescent_1"

    return fine


def lineage_from_plot_label(plot_label: str) -> str:
    s = str(plot_label)

    if s.startswith("quiescent_"):
        return "germ"

    if s == "Unassigned":
        return "Unassigned"

    return s.split("_", 1)[0] if "_" in s else s


def get_abundance_and_factors(adata: sc.AnnData, abundance_key: str) -> tuple[np.ndarray, list[str]]:
    if abundance_key not in adata.obsm:
        raise KeyError(f"Missing adata.obsm['{abundance_key}']. Available obsm keys: {list(adata.obsm.keys())}")

    W = adata.obsm[abundance_key]

    if isinstance(W, pd.DataFrame):
        factor_names = [strip_factor_prefix(x) for x in W.columns.astype(str)]
        W = W.to_numpy(dtype=np.float64)
    else:
        W = np.asarray(W, dtype=np.float64)
        factor_names = adata.uns.get("mod", {}).get("factor_names", None)

        if factor_names is None:
            raise KeyError(
                "Abundance matrix is not a DataFrame and adata.uns['mod']['factor_names'] is missing."
            )

        factor_names = [strip_factor_prefix(x) for x in factor_names]

    if W.ndim != 2:
        raise ValueError(f"Abundance matrix is not 2D: shape={W.shape}")

    if W.shape[1] != len(factor_names):
        raise ValueError(
            f"Factor name mismatch: abundance has {W.shape[1]} columns, but {len(factor_names)} names were found."
        )

    return W, list(map(str, factor_names))


def get_spatial_xy(adata: sc.AnnData) -> tuple[np.ndarray, np.ndarray]:
    if "spatial" in adata.obsm:
        xy = np.asarray(adata.obsm["spatial"])
        if xy.ndim == 2 and xy.shape[1] >= 2:
            return xy[:, 0].astype(float), xy[:, 1].astype(float)

    if {"x", "y"}.issubset(adata.obs.columns):
        x = pd.to_numeric(adata.obs["x"], errors="coerce").to_numpy(dtype=float)
        y = pd.to_numeric(adata.obs["y"], errors="coerce").to_numpy(dtype=float)
        ok = np.isfinite(x) & np.isfinite(y)

        if ok.mean() > 0.99:
            return x, y

    raise KeyError("No spatial coordinates found. Expected adata.obsm['spatial'] or obs columns x/y.")


def relative_abundance(W: np.ndarray) -> np.ndarray:
    row_sum = W.sum(axis=1, keepdims=True)
    row_sum = np.maximum(row_sum, 1e-12)
    return W / row_sum


def compute_lineage_sums(R_common: np.ndarray, common_plot_labels: list[str]) -> dict[str, np.ndarray]:
    lineage_sums = {
        lin: np.zeros(R_common.shape[0], dtype=float)
        for lin in PLOT_CELL_TYPES
    }

    for j, label in enumerate(common_plot_labels):
        lineage = lineage_from_plot_label(label)
        if lineage in lineage_sums:
            lineage_sums[lineage] += R_common[:, j]

    return lineage_sums


def thresholded_dominant_celltype(
    R_common: np.ndarray,
    common_plot_labels: list[str],
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    n = R_common.shape[0]
    lineage_sums = compute_lineage_sums(R_common, common_plot_labels)

    label_arr = np.full(n, "Unassigned", dtype=object)

    for i in range(n):
        eligible = []
        for lineage in PLOT_CELL_TYPES:
            threshold = LINEAGE_RELATIVE_THRESHOLD[lineage]
            if lineage_sums[lineage][i] >= threshold:
                eligible.append(lineage)

        if len(eligible) == 0:
            continue

        label_arr[i] = max(eligible, key=lambda lineage: lineage_sums[lineage][i])

    return label_arr, lineage_sums


def topk_fine_states(R_common: np.ndarray, common_plot_labels: list[str], k: int = 3) -> np.ndarray:
    k = min(k, R_common.shape[1])

    idx = np.argpartition(-R_common, kth=k - 1, axis=1)[:, :k]
    vals = np.take_along_axis(R_common, idx, axis=1)
    order = np.argsort(-vals, axis=1)

    idx_sorted = np.take_along_axis(idx, order, axis=1)
    vals_sorted = np.take_along_axis(vals, order, axis=1)

    names = np.array(common_plot_labels, dtype=object)
    top_names = names[idx_sorted]

    out = []

    for i in range(R_common.shape[0]):
        parts = [f"{top_names[i, j]}:{vals_sorted[i, j]:.3f}" for j in range(k)]
        out.append(", ".join(parts))

    return np.array(out, dtype=object)


def scan_common_plot_labels(sample_ids: list[str], mapped_dir: Path, abundance_key: str) -> list[str]:
    factor_sets = []

    print("Scanning factor names for common fine-subcluster labels.")

    for sample_id in sample_ids:
        mapped_path = mapped_dir / f"{sample_id}.cell2location_mapped.h5ad"

        if not mapped_path.exists():
            raise FileNotFoundError(f"Mapped H5AD not found: {mapped_path}")

        adata = sc.read_h5ad(mapped_path)
        _, factors = get_abundance_and_factors(adata, abundance_key)

        plot_labels = [relabel_for_interactive_compatibility(x) for x in factors]
        factor_sets.append(set(plot_labels))

        print(f"  {sample_id}: {len(plot_labels)} factors")

    common = sorted(set.intersection(*factor_sets))

    if len(common) == 0:
        raise RuntimeError("No common fine-subcluster labels across selected samples.")

    print(f"Common fine-subcluster labels: {len(common)}")

    return common


def load_sample_payload(
    sample_id: str,
    mapped_path: Path,
    abundance_key: str,
    common_plot_labels: list[str],
) -> tuple[pd.DataFrame, dict]:
    if not mapped_path.exists():
        raise FileNotFoundError(f"Mapped H5AD not found: {mapped_path}")

    adata = sc.read_h5ad(mapped_path)

    x, y = get_spatial_xy(adata)
    W, factor_names_raw = get_abundance_and_factors(adata, abundance_key)

    factor_names_plot = [relabel_for_interactive_compatibility(x) for x in factor_names_raw]

    R = relative_abundance(W)

    idx_map = {factor: i for i, factor in enumerate(factor_names_plot)}
    missing = [factor for factor in common_plot_labels if factor not in idx_map]

    if missing:
        raise RuntimeError(f"{sample_id}: missing common factors, example: {missing[:10]}")

    cols = [idx_map[factor] for factor in common_plot_labels]
    R_common = R[:, cols]

    dominant_major, lineage_sums = thresholded_dominant_celltype(
        R_common=R_common,
        common_plot_labels=common_plot_labels,
    )
    top3 = topk_fine_states(R_common, common_plot_labels, k=3)

    plot_df = pd.DataFrame({
        "sample_id": sample_id,
        "sample_label": SAMPLE_LABELS.get(sample_id, sample_id),
        "barcode": adata.obs["barcode"].astype(str).values if "barcode" in adata.obs.columns else adata.obs_names.astype(str),
        "x": x,
        "y": y,
        "dominant_major_cell_type": dominant_major,
        "top3_fine_states": top3,
    })

    for lineage in PLOT_CELL_TYPES:
        plot_df[f"relative_abundance_{lineage}"] = lineage_sums[lineage]

    counts = (
        plot_df["dominant_major_cell_type"]
        .value_counts(dropna=False)
        .rename_axis("major_cell_type")
        .reset_index(name="n_locations")
    )
    counts["sample_id"] = sample_id
    counts["sample_label"] = SAMPLE_LABELS.get(sample_id, sample_id)
    counts["fraction_locations"] = counts["n_locations"] / plot_df.shape[0]

    summary = {
        "sample_id": sample_id,
        "sample_label": SAMPLE_LABELS.get(sample_id, sample_id),
        "mapped_h5ad_name": mapped_path.name,
        "n_locations": int(plot_df.shape[0]),
        "n_genes": int(adata.n_vars),
        "n_factors_raw": int(W.shape[1]),
        "n_common_factors_after_relabeling": int(len(common_plot_labels)),
        "abundance_key": abundance_key,
        "assigned_fraction": float(np.mean(dominant_major != "Unassigned")),
        "counts": counts,
    }

    return plot_df, summary


def plot_one_celltype(
    cell_type: str,
    payloads: dict[str, pd.DataFrame],
    figure_dir: Path,
) -> None:
    color = MAJOR_COLORS.get(cell_type, "#4D4D4D")

    fig, axes = plt.subplots(
        1,
        len(SAMPLE_ORDER),
        figsize=(7.4, 2.6),
        constrained_layout=True,
    )

    if len(SAMPLE_ORDER) == 1:
        axes = [axes]

    for ax, sample_id in zip(axes, SAMPLE_ORDER):
        df = payloads[sample_id]

        x = df["x"].to_numpy(dtype=float)
        y = df["y"].to_numpy(dtype=float)
        mask = df["dominant_major_cell_type"].astype(str).to_numpy() == cell_type

        ax.scatter(
            x,
            y,
            c=BACKGROUND_COLOR,
            s=BACKGROUND_POINT_SIZE,
            alpha=BACKGROUND_ALPHA,
            linewidths=0,
            rasterized=True,
        )

        if mask.sum() > 0:
            ax.scatter(
                x[mask],
                y[mask],
                c=color,
                s=FOREGROUND_POINT_SIZE,
                alpha=FOREGROUND_ALPHA,
                linewidths=0,
                rasterized=True,
            )
        ax.set_title(SAMPLE_LABELS.get(sample_id, sample_id))
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_xticks([])
        ax.set_yticks([])

        for spine in ax.spines.values():
            spine.set_visible(False)

        pct = 100.0 * float(mask.mean()) if len(mask) > 0 else 0.0
        ax.text(
            0.02,
            0.02,
            f"{int(mask.sum()):,} locations\n{pct:.1f}%",
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=7,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.80, pad=1.5),
        )

    fig.suptitle(cell_type.capitalize(), fontweight="bold", y=1.03)

    out_prefix = figure_dir / f"spatial_bmk_celltype_highlight_{cell_type}"
    save_figure(fig, out_prefix)


def write_figure_index(figure_dir: Path, cell_types: list[str]) -> None:
    lines = [
        "BMK-ST cell-type highlight spatial figures",
        "",
        "Assignment logic matches the original interactive thresholded dominant CELLTYPE script.",
        "Grey points are all spatial locations.",
        "Colored points are locations assigned to the target major cell type.",
        "",
        "Important assignment rule:",
        "- germ_0 and germ_6 are internally relabeled as quiescent_0 and quiescent_1, then counted as germ.",
        "- No separate degenerated major-cell-type category is used in this script.",
        "",
        "Cell types:",
    ]

    lines.extend([f"- {x}" for x in cell_types])

    (figure_dir / "spatial_bmk_celltype_highlight_figure_index.txt").write_text(
        "\n".join(lines) + "\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate static BMK-ST major-cell-type highlight maps from cell2location outputs."
    )
    parser.add_argument("--manifest", default=None, help="CSV manifest. Defaults to metadata/spatial_bmk_samples.csv.")
    parser.add_argument("--results-root", default=None, help="Results root.")
    parser.add_argument("--abundance-key", default="means_cell_abundance_w_sf")
    parser.add_argument(
        "--cell-types",
        default=",".join(PLOT_CELL_TYPES),
        help="Comma-separated major cell types to plot.",
    )

    args = parser.parse_args()

    project_root = repo_root()
    results_root = Path(args.results_root).resolve() if args.results_root else default_results_root(project_root)
    manifest = Path(args.manifest).resolve() if args.manifest else project_root / "metadata" / "spatial_bmk_samples.csv"

    mapped_dir = results_root / "spatial_bmk_cell2location" / "objects"

    out_root = results_root / "spatial_bmk_celltype_highlight_figures"
    figure_dir = out_root / "figures"
    table_dir = out_root / "tables"
    log_dir = out_root / "logs"

    figure_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    setup_plot_style()
    print("Generating BMK-ST cell-type highlight spatial figures")
    print(f"  manifest file: {manifest.name}")
    print(f"  mapped object directory: {mapped_dir}")
    print(f"  output directory: {out_root}")
    print("  assignment logic: original interactive-compatible thresholded dominant cell type")
    print("  germ_0/germ_6 handling: counted as germ, not as separate degenerated")

    rows = load_manifest(manifest)

    selected_samples = [row["sample_id"] for row in rows]
    missing_order = [x for x in SAMPLE_ORDER if x not in selected_samples]
    if missing_order:
        raise ValueError(f"Expected samples are missing from manifest/include=TRUE: {missing_order}")

    common_plot_labels = scan_common_plot_labels(
        sample_ids=SAMPLE_ORDER,
        mapped_dir=mapped_dir,
        abundance_key=args.abundance_key,
    )

    payloads = {}
    summaries = []
    count_tables = []

    for sample_id in SAMPLE_ORDER:
        mapped_path = mapped_dir / f"{sample_id}.cell2location_mapped.h5ad"

        print("=" * 80)
        print(f"Sample: {sample_id}")
        print(f"Mapped H5AD: {mapped_path.name}")

        plot_df, summary = load_sample_payload(
            sample_id=sample_id,
            mapped_path=mapped_path,
            abundance_key=args.abundance_key,
            common_plot_labels=common_plot_labels,
        )

        payloads[sample_id] = plot_df
        summaries.append({k: v for k, v in summary.items() if k != "counts"})
        count_tables.append(summary["counts"])

        plot_df.to_csv(table_dir / f"{sample_id}.spatial_bmk_celltype_highlight_assignments.csv", index=False)

        print(f"  locations: {summary['n_locations']}")
        print(f"  genes: {summary['n_genes']}")
        print(f"  raw factors: {summary['n_factors_raw']}")
        print(f"  common factors after relabeling: {summary['n_common_factors_after_relabeling']}")
        print(f"  assigned fraction: {summary['assigned_fraction']:.3f}")
        print("  assigned counts:")
        print(summary["counts"].sort_values(["major_cell_type"]).to_string(index=False))

    requested_cell_types = [x.strip() for x in args.cell_types.split(",") if x.strip()]
    cell_types = [x for x in requested_cell_types if x in PLOT_CELL_TYPES]

    if not cell_types:
        raise ValueError("No valid cell types selected for plotting.")

    for cell_type in cell_types:
        plot_one_celltype(
            cell_type=cell_type,
            payloads=payloads,
            figure_dir=figure_dir,
        )
        print(f"  wrote figure for: {cell_type}")
    summary_df = pd.DataFrame(summaries)
    counts_df = pd.concat(count_tables, ignore_index=True)

    summary_df.to_csv(table_dir / "spatial_bmk_celltype_highlight_sample_summary.csv", index=False)
    counts_df.to_csv(table_dir / "spatial_bmk_celltype_highlight_location_counts.csv", index=False)

    pd.DataFrame({"common_plot_label": common_plot_labels}).to_csv(
        table_dir / "spatial_bmk_celltype_highlight_common_factors.csv",
        index=False,
    )

    write_figure_index(figure_dir, cell_types)

    safe_write_json(
        log_dir / "spatial_bmk_celltype_highlight_package_versions.json",
        {
            "python": sys.version,
            "platform": platform.platform(),
            "scanpy": sc.__version__,
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "matplotlib": matplotlib.__version__,
        },
    )

    with (log_dir / "spatial_bmk_celltype_highlight_summary.txt").open("w") as handle:
        handle.write("BMK-ST cell-type highlight spatial figures\n")
        handle.write("Input: cell2location-mapped BMK-ST H5AD objects\n")
        handle.write("Assignment logic: original interactive-compatible thresholded dominant cell type\n")
        handle.write("germ_0 and germ_6 are counted as germ through quiescent_0/quiescent_1 internal relabeling.\n")
        handle.write("No separate degenerated major-cell-type category is used.\n")
        handle.write(f"Manifest file: {manifest.name}\n")
        handle.write(f"Mapped object directory: {mapped_dir.name}\n")
        handle.write(f"Cell types plotted: {', '.join(cell_types)}\n")
        handle.write("Output figures: spatial_bmk_celltype_highlight_<cell_type>.pdf/png\n")

    print("\nDone.")
    print(f"Figures: {figure_dir}")
    print(f"Tables: {table_dir}")
    print(f"Logs: {log_dir}")


if __name__ == "__main__":
    main()
