#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import inspect
import json
import os
import platform
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
import torch

from cell2location.models import Cell2location
from lightning.pytorch.callbacks import Callback, ModelCheckpoint


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


def safe_write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w") as handle:
        json.dump(payload, handle, indent=2)
    tmp.replace(path)


def load_reference_signatures(path: Path) -> pd.DataFrame:
    ref = pd.read_csv(path, index_col=0)
    ref = ref[~ref.index.duplicated(keep="first")]
    ref.index = ref.index.astype(str)
    ref.columns = ref.columns.astype(str)
    return ref
def load_fine_to_major(path: Path) -> dict[str, str]:
    if not path.exists():
        raise FileNotFoundError(f"Fine-to-major mapping not found: {path}")

    df = pd.read_csv(path)

    required = {"fine_subcluster", "major_cell_type"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Fine-to-major mapping missing columns: {sorted(missing)}")

    return dict(zip(df["fine_subcluster"].astype(str), df["major_cell_type"].astype(str)))


def strip_factor_prefix(factor_name: str) -> str:
    s = str(factor_name)

    if "mu_fg_" in s:
        return s.split("mu_fg_", 1)[1]

    if s.startswith("means_cell_abundance_w_sf_"):
        return s.replace("means_cell_abundance_w_sf_", "", 1)

    parts = s.split("_")
    if len(parts) >= 2 and parts[-1].isdigit():
        return "_".join(parts[-2:])

    return s


def plan_kwargs_with_lr(lr: float) -> dict:
    try:
        from scvi.train._trainingplans import PyroTrainingPlan as plan_cls
    except Exception as exc:
        print(f"Could not import PyroTrainingPlan: {exc}")
        return {}

    try:
        sig = inspect.signature(plan_cls.__init__)
        params = set(sig.parameters.keys())

        candidates = [
            ("lr", {"lr": lr}),
            ("optimizer_kwargs", {"optimizer_kwargs": {"lr": lr}}),
            ("optim_kwargs", {"optim_kwargs": {"lr": lr}}),
        ]

        for name, payload in candidates:
            if name in params:
                print(f"Setting learning rate through plan_kwargs['{name}']")
                return payload

        print("PyroTrainingPlan does not expose an LR kwarg; using default optimizer learning rate.")
        return {}

    except Exception as exc:
        print(f"Could not inspect PyroTrainingPlan signature ({exc}); using optim_kwargs fallback.")
        return {"optim_kwargs": {"lr": lr}}


class StopOnNonFiniteLoss(Callback):
    def on_train_batch_end(self, trainer, pl_module, outputs, batch, batch_idx):
        metrics = trainer.callback_metrics

        for key in ["loss", "elbo_train", "train_loss"]:
            if key not in metrics:
                continue

            value = metrics[key]

            try:
                val = float(value.detach().cpu().item())
            except Exception:
                try:
                    val = float(value)
                except Exception:
                    continue

            if not np.isfinite(val):
                raise RuntimeError(f"Non-finite training metric detected: {key}={val}")
def ensure_counts_csr(x) -> sparse.csr_matrix:
    if sparse.issparse(x):
        x = x.tocsr()
        if x.nnz > 0:
            x.data = np.rint(x.data).astype(np.float32, copy=False)
            x.data[x.data < 0] = 0.0
        return x.astype(np.float32, copy=False)

    arr = np.asarray(x)
    arr = np.rint(arr).astype(np.float32)
    arr[arr < 0] = 0.0
    return sparse.csr_matrix(arr)


def describe_libsizes(x_csr: sparse.csr_matrix) -> dict:
    lib = np.asarray(x_csr.sum(axis=1)).ravel().astype(np.float64)

    if lib.size == 0:
        return {
            "min": 0.0,
            "q25": 0.0,
            "median": 0.0,
            "q75": 0.0,
            "q90": 0.0,
            "q99": 0.0,
            "max": 0.0,
        }

    q = np.quantile(lib, [0, 0.25, 0.5, 0.75, 0.9, 0.99, 1.0])

    return {
        "min": float(q[0]),
        "q25": float(q[1]),
        "median": float(q[2]),
        "q75": float(q[3]),
        "q90": float(q[4]),
        "q99": float(q[5]),
        "max": float(q[6]),
    }


def subset_to_reference_genes(adata: sc.AnnData, ref_sig: pd.DataFrame) -> tuple[sc.AnnData, pd.DataFrame]:
    spatial_genes = pd.Index(adata.var_names.astype(str))
    ref_genes = pd.Index(ref_sig.index.astype(str))

    genes_use = ref_genes.intersection(spatial_genes)

    if len(genes_use) < 500:
        raise RuntimeError(f"Too few overlapping genes between spatial object and reference: {len(genes_use)}")

    adata = adata[:, genes_use].copy()
    adata.X = ensure_counts_csr(adata.X)

    ref_sub = ref_sig.loc[genes_use, :].copy()

    lib = np.asarray(adata.X.sum(axis=1)).ravel().astype(np.float64)
    keep = lib > 0

    if int(keep.sum()) < 1000:
        raise RuntimeError("Too few locations remain after filtering zero counts over reference genes.")

    if int(keep.sum()) < adata.n_obs:
        print(f"Filtering zero-count locations after gene intersection: keeping {int(keep.sum())} / {adata.n_obs}")
        adata = adata[keep, :].copy()

    return adata, ref_sub


def get_factor_names(adata: sc.AnnData, ref_sub: pd.DataFrame) -> list[str]:
    factor_names = adata.uns.get("mod", {}).get("factor_names", None)

    if factor_names is not None:
        factor_names = [strip_factor_prefix(x) for x in list(factor_names)]

    if factor_names is None or len(factor_names) != ref_sub.shape[1]:
        factor_names = [str(x) for x in ref_sub.columns]

    return factor_names
def add_hard_labels(adata: sc.AnnData, ref_sub: pd.DataFrame, fine_to_major: dict[str, str]) -> pd.DataFrame:
    if "means_cell_abundance_w_sf" not in adata.obsm:
        raise KeyError("means_cell_abundance_w_sf not found in adata.obsm after posterior export.")

    factor_names = get_factor_names(adata, ref_sub)
    abundance = np.asarray(adata.obsm["means_cell_abundance_w_sf"])

    if abundance.shape[1] != len(factor_names):
        raise ValueError(
            f"Abundance matrix has {abundance.shape[1]} columns but factor_names has {len(factor_names)} entries."
        )

    imax = np.argmax(abundance, axis=1)
    fine_labels = [factor_names[i] for i in imax]
    major_labels = [fine_to_major.get(x, "unmapped") for x in fine_labels]

    adata.obs["snRNA_fine_subcluster_cell2loc"] = pd.Categorical(fine_labels)
    adata.obs["snRNA_cell_type_cell2loc"] = pd.Categorical(major_labels)
    adata.obs["cell2loc_max_abundance"] = abundance[np.arange(abundance.shape[0]), imax].astype(np.float32)

    return pd.DataFrame(abundance, index=adata.obs_names, columns=factor_names)


def run_one_sample(
    sample_id: str,
    canonical_h5ad: Path,
    ref_sig: pd.DataFrame,
    fine_to_major: dict[str, str],
    out_root: Path,
    max_epochs: int,
    posterior_samples: int,
    batch_size: int,
    dl_num_workers: int,
    log_every_n_steps: int,
    learning_rate: float,
    precision: str,
    n_cells_per_location: int,
    detection_alpha: float,
    save_checkpoints: bool,
) -> None:
    print("=" * 80)
    print(f"Sample: {sample_id}")
    print(f"Canonical H5AD: {canonical_h5ad.name}")

    if not canonical_h5ad.exists():
        raise FileNotFoundError(f"Canonical spatial H5AD not found: {canonical_h5ad}")

    object_dir = out_root / "objects"
    table_dir = out_root / "tables"
    log_dir = out_root / "logs"
    checkpoint_dir = out_root / "checkpoints" / sample_id

    object_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()

    adata = sc.read_h5ad(canonical_h5ad)
    adata.var_names = adata.var_names.astype(str)
    adata.obs_names = adata.obs_names.astype(str)
    adata.X = ensure_counts_csr(adata.X)

    print(f"Loaded spatial object: {adata.shape}")

    adata, ref_sub = subset_to_reference_genes(adata, ref_sig)
    lib_report = describe_libsizes(adata.X)

    print(f"Overlap genes with reference: {adata.n_vars}")
    print(f"Locations after reference-gene filtering: {adata.n_obs}")
    print(f"Total-count median after gene filtering: {lib_report['median']}")

    Cell2location.setup_anndata(adata)

    model = Cell2location(
        adata=adata,
        cell_state_df=ref_sub,
        N_cells_per_location=n_cells_per_location,
        detection_alpha=detection_alpha,
    )
    try:
        import scvi

        scvi.settings.dl_num_workers = dl_num_workers

        if hasattr(scvi.settings, "dl_pin_memory"):
            scvi.settings.dl_pin_memory = True

        if hasattr(scvi.settings, "dl_persistent_workers"):
            scvi.settings.dl_persistent_workers = True

    except Exception as exc:
        print(f"scvi dataloader settings were not fully configurable: {exc}")

    callbacks = [StopOnNonFiniteLoss()]

    if save_checkpoints:
        callbacks.append(
            ModelCheckpoint(
                dirpath=str(checkpoint_dir),
                filename=f"{sample_id}" + "__{epoch:04d}",
                save_top_k=-1,
                every_n_epochs=250,
                save_last=True,
                monitor=None,
            )
        )

    print(
        f"Training cell2location: epochs={max_epochs}, batch_size={batch_size}, "
        f"lr={learning_rate}, precision={precision}"
    )

    model.train(
        max_epochs=max_epochs,
        batch_size=batch_size,
        train_size=1.0,
        accelerator="gpu" if torch.cuda.is_available() else "cpu",
        log_every_n_steps=log_every_n_steps,
        enable_checkpointing=save_checkpoints,
        default_root_dir=str(checkpoint_dir),
        callbacks=callbacks,
        precision=precision,
        plan_kwargs=plan_kwargs_with_lr(learning_rate),
        enable_progress_bar=True,
    )

    print(f"Exporting posterior: samples={posterior_samples}")

    model.export_posterior(
        adata,
        sample_kwargs={
            "num_samples": posterior_samples,
            "batch_size": batch_size,
        },
    )

    abundance_df = add_hard_labels(adata, ref_sub, fine_to_major)

    mapped_h5ad = object_dir / f"{sample_id}.cell2location_mapped.h5ad"
    abundance_csv = table_dir / f"{sample_id}.cell2loc_abundance_means.csv"
    label_csv = table_dir / f"{sample_id}.cell2loc_location_labels.csv"
    summary_json = log_dir / f"{sample_id}.cell2location_summary.json"

    print(f"Writing mapped H5AD: {mapped_h5ad.name}")
    adata.write_h5ad(mapped_h5ad)

    abundance_df.to_csv(abundance_csv)

    label_cols = [
        "sample_id",
        "gestational_week",
        "platform",
        "barcode",
        "x",
        "y",
        "snRNA_fine_subcluster_cell2loc",
        "snRNA_cell_type_cell2loc",
        "cell2loc_max_abundance",
    ]
    label_cols = [x for x in label_cols if x in adata.obs.columns]

    adata.obs[label_cols].to_csv(label_csv)

    elapsed_min = (time.time() - t0) / 60.0
    summary = {
        "sample_id": sample_id,
        "canonical_h5ad_name": canonical_h5ad.name,
        "mapped_h5ad_name": mapped_h5ad.name,
        "reference_signature_name": "reference_signatures_fine_subcluster.csv",
        "n_locations": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_factors": int(ref_sub.shape[1]),
        "nnz": int(adata.X.nnz),
        "total_counts_after_gene_filter": lib_report,
        "max_epochs": int(max_epochs),
        "posterior_samples": int(posterior_samples),
        "batch_size": int(batch_size),
        "learning_rate": float(learning_rate),
        "precision": precision,
        "n_cells_per_location": int(n_cells_per_location),
        "detection_alpha": float(detection_alpha),
        "torch_version": torch.__version__,
        "cuda_available": bool(torch.cuda.is_available()),
        "cuda_device": torch.cuda.get_device_name(0) if torch.cuda.is_available() else None,
        "elapsed_minutes": float(elapsed_min),
    }

    safe_write_json(summary_json, summary)

    print(f"Done sample {sample_id}. Elapsed minutes: {elapsed_min:.2f}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run cell2location on canonical BMK-ST spatial objects.")
    parser.add_argument("--manifest", default=None, help="CSV manifest. Defaults to metadata/spatial_bmk_samples.csv.")
    parser.add_argument("--results-root", default=None, help="Results root.")
    parser.add_argument("--sample-id", default=None, help="Optional single sample_id to run.")
    parser.add_argument("--epochs", type=int, default=5000)
    parser.add_argument("--posterior-samples", type=int, default=500)
    parser.add_argument("--batch-size", type=int, default=2048)
    parser.add_argument("--dl-num-workers", type=int, default=8)
    parser.add_argument("--log-every-n-steps", type=int, default=50)
    parser.add_argument("--learning-rate", type=float, default=1e-3)
    parser.add_argument("--precision", default="32-true")
    parser.add_argument("--torch-matmul-precision", default="high")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--n-cells-per-location", type=int, default=30)
    parser.add_argument("--detection-alpha", type=float, default=200.0)
    parser.add_argument("--no-checkpoints", action="store_true")

    args = parser.parse_args()

    project_root = repo_root()
    results_root = Path(args.results_root).resolve() if args.results_root else default_results_root(project_root)
    manifest = Path(args.manifest).resolve() if args.manifest else project_root / "metadata" / "spatial_bmk_samples.csv"

    canonical_dir = results_root / "spatial_bmk_canonical_objects" / "objects"
    reference_dir = results_root / "cell2location_reference"
    out_root = results_root / "spatial_bmk_cell2location"

    ref_path = reference_dir / "reference_signatures_fine_subcluster.csv"
    fine_to_major_path = reference_dir / "tables" / "fine_subcluster_to_major_cell_type.csv"

    if not manifest.exists():
        raise FileNotFoundError(f"Manifest not found: {manifest}")
    if not ref_path.exists():
        raise FileNotFoundError(f"Reference signatures not found: {ref_path}")
    if not canonical_dir.exists():
        raise FileNotFoundError(f"Canonical spatial object directory not found: {canonical_dir}")

    np.random.seed(args.seed)
    torch.manual_seed(args.seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(args.seed)
        torch.set_float32_matmul_precision(args.torch_matmul_precision)

    print("Running BMK-ST cell2location")
    print(f"  manifest file: {manifest.name}")
    print(f"  canonical object directory: {canonical_dir}")
    print(f"  reference signature file: {ref_path.name}")
    print(f"  output directory: {out_root}")
    print(f"  torch: {torch.__version__}")
    print(f"  cuda available: {torch.cuda.is_available()}")

    if torch.cuda.is_available():
        print(f"  cuda device: {torch.cuda.get_device_name(0)}")

    ref_sig = load_reference_signatures(ref_path)
    fine_to_major = load_fine_to_major(fine_to_major_path)

    rows = []
    with manifest.open() as handle:
        reader = csv.DictReader(handle)
        required = {"sample_id", "include"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Manifest is missing columns: {sorted(missing)}")

        for row in reader:
            if not parse_bool(row["include"]):
                continue
            if args.sample_id is not None and row["sample_id"] != args.sample_id:
                continue
            rows.append(row)

    if not rows:
        raise ValueError("No samples selected for cell2location.")

    out_root.mkdir(parents=True, exist_ok=True)

    package_versions = {
        "python": sys.version,
        "platform": platform.platform(),
        "scanpy": sc.__version__,
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "torch": torch.__version__,
    }
    safe_write_json(out_root / "logs" / "cell2location_package_versions.json", package_versions)

    for row in rows:
        sample_id = row["sample_id"]
        canonical_h5ad = canonical_dir / f"{sample_id}.spatial_bmk_canonical.h5ad"

        run_one_sample(
            sample_id=sample_id,
            canonical_h5ad=canonical_h5ad,
            ref_sig=ref_sig,
            fine_to_major=fine_to_major,
            out_root=out_root,
            max_epochs=args.epochs,
            posterior_samples=args.posterior_samples,
            batch_size=args.batch_size,
            dl_num_workers=args.dl_num_workers,
            log_every_n_steps=args.log_every_n_steps,
            learning_rate=args.learning_rate,
            precision=args.precision,
            n_cells_per_location=args.n_cells_per_location,
            detection_alpha=args.detection_alpha,
            save_checkpoints=not args.no_checkpoints,
        )

    print("\nAll selected samples complete.")


if __name__ == "__main__":
    main()
