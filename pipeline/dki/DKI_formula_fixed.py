#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Formula-fixed DKI/cNODE script.

Core fix:
    wrong old formula:
        dy = y * (out - sum(y))

    fixed formula:
        dy = y * (out - sum(y * out))

Recommended stable full-run setting:
    --time-mode endpoint
    --ode-method rk4
    --ode-step-size 0.01
"""

import argparse
import copy
import csv
import json
import os
import random
import sys
from typing import Dict, List, Tuple

import numpy as np
import torch
from torch.optim.lr_scheduler import ReduceLROnPlateau
from torchdiffeq import odeint


def parse_bool(value):
    if isinstance(value, bool):
        return value
    value = str(value).strip().lower()
    if value in {"1", "true", "yes", "y"}:
        return True
    if value in {"0", "false", "no", "n"}:
        return False
    raise argparse.ArgumentTypeError(f"Cannot parse boolean value: {value}")


def parse_args():
    parser = argparse.ArgumentParser(description="Formula-fixed cNODE DKI runner")

    parser.add_argument("--train-file", required=True, help="Absolute path to Ptrain.csv")
    parser.add_argument("--test-file", required=True, help="Absolute path to Ptest.csv")
    parser.add_argument("--output-dir", required=True, help="Absolute output directory")
    parser.add_argument("--project", default=None, help="Project ID for metrics")

    parser.add_argument("--max-epochs", type=int, default=3000)
    parser.add_argument("--lr", type=float, default=0.005)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--early-stop-patience", type=int, default=200)
    parser.add_argument("--weight-decay", type=float, default=1e-4)
    parser.add_argument("--seed", type=int, default=20260629)

    parser.add_argument("--device", default="auto", choices=["auto", "cpu", "cuda"])

    parser.add_argument(
        "--activation",
        default="relu",
        choices=["relu", "tanh", "gelu", "elu", "leaky_relu", "none"],
        help="Default relu keeps the current engineering version. Use none to be closer to original DKI.",
    )

    parser.add_argument(
        "--bottleneck",
        default="auto",
        choices=["auto", "true", "false"],
        help="auto uses bottleneck when N > 500.",
    )

    parser.add_argument("--normalize-z", type=parse_bool, default=True)
    parser.add_argument("--chunk-size", type=int, default=2000)

    parser.add_argument(
        "--time-mode",
        default="endpoint",
        choices=["endpoint", "old100"],
        help="endpoint uses [0, 0.99]. old100 uses 100 output time points from 0 to 0.99.",
    )

    parser.add_argument(
        "--ode-method",
        default="rk4",
        choices=["rk4", "dopri5", "euler", "midpoint"],
        help="Use rk4 for stable fixed-step integration. dopri5 may fail with underflow on some projects.",
    )

    parser.add_argument(
        "--ode-step-size",
        type=float,
        default=0.01,
        help="Fixed step size for rk4/euler/midpoint.",
    )

    parser.add_argument("--rtol", type=float, default=1e-3)
    parser.add_argument("--atol", type=float, default=1e-4)

    return parser.parse_args()


def set_seed(seed: int):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def resolve_device(device_arg: str) -> str:
    if device_arg == "auto":
        return "cuda" if torch.cuda.is_available() else "cpu"
    if device_arg == "cuda" and not torch.cuda.is_available():
        print("⚠️ CUDA requested but unavailable. Falling back to CPU.")
        return "cpu"
    return device_arg


def load_matrix(path: str, expected_rows: int = None) -> np.ndarray:
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}")

    x = np.loadtxt(path, delimiter=",")
    x = np.asarray(x, dtype=np.float32)

    if x.ndim == 1:
        if expected_rows is not None and x.size == expected_rows:
            x = x.reshape((expected_rows, 1))
        else:
            x = x.reshape((-1, 1))

    if x.ndim != 2:
        raise ValueError(f"Expected 2D matrix, got shape {x.shape}: {path}")

    x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)
    x[x < 0] = 0.0
    return x


def normalize_columns(matrix: np.ndarray) -> np.ndarray:
    x = np.asarray(matrix, dtype=np.float32).copy()
    x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)
    x[x < 0] = 0.0
    sums = x.sum(axis=0, keepdims=True)
    return np.divide(x, sums, out=np.zeros_like(x, dtype=np.float32), where=sums > 0)


def normalize_prediction_rows(raw: np.ndarray) -> Tuple[np.ndarray, Dict[str, float]]:
    x = np.asarray(raw, dtype=np.float64).copy()
    x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)
    x[x < 0] = 0.0

    row_sums_before = x.sum(axis=1, keepdims=True)
    zero_rows = int(np.sum(row_sums_before.ravel() == 0))

    norm = np.divide(x, row_sums_before, out=np.zeros_like(x), where=row_sums_before > 0)
    row_sums_after = norm.sum(axis=1)

    metrics = {
        "row_sum_min": float(np.min(row_sums_after)) if row_sums_after.size else None,
        "row_sum_median": float(np.median(row_sums_after)) if row_sums_after.size else None,
        "row_sum_max": float(np.max(row_sums_after)) if row_sums_after.size else None,
        "min_value": float(np.min(norm)) if norm.size else None,
        "max_value": float(np.max(norm)) if norm.size else None,
        "zero_rows_after_normalization": zero_rows,
    }

    return norm.astype(np.float32), metrics


def raw_row_sum_summary(raw: np.ndarray) -> Dict[str, float]:
    x = np.asarray(raw, dtype=np.float64)
    x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)
    row_sums = x.sum(axis=1)
    return {
        "min": float(np.min(row_sums)) if row_sums.size else None,
        "median": float(np.median(row_sums)) if row_sums.size else None,
        "max": float(np.max(row_sums)) if row_sums.size else None,
    }


def process_data(P_species_by_samples: np.ndarray, normalize_z: bool) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    Input:
        P_species_by_samples: species x samples

    Output:
        P: samples x species relative abundance composition
        Z: samples x species binary assemblage vector, optionally normalized by number of present species
    """
    P_raw = np.asarray(P_species_by_samples, dtype=np.float32).copy()
    P_raw = np.nan_to_num(P_raw, nan=0.0, posinf=0.0, neginf=0.0)
    P_raw[P_raw < 0] = 0.0

    Z = (P_raw > 0).astype(np.float32)
    P = normalize_columns(P_raw).astype(np.float32)

    if normalize_z:
        col_sum_Z = Z.sum(axis=0, keepdims=True)
        Z = np.divide(Z, col_sum_Z, out=np.zeros_like(Z, dtype=np.float32), where=col_sum_Z > 0)

    return torch.from_numpy(P.T), torch.from_numpy(Z.T)


def get_activation(name: str):
    name = name.lower()
    if name == "none":
        return lambda x: x
    if name == "relu":
        return torch.relu
    if name == "tanh":
        return torch.tanh
    if name == "gelu":
        return torch.nn.functional.gelu
    if name == "elu":
        return torch.nn.functional.elu
    if name == "leaky_relu":
        return torch.nn.functional.leaky_relu
    raise ValueError(f"Unsupported activation: {name}")


class ODEFunc(torch.nn.Module):
    def __init__(self, N: int, activation: str, bottleneck: str):
        super().__init__()

        self.N = N
        self.activation_name = activation
        self.activation = get_activation(activation)

        if bottleneck == "auto":
            self.use_bottleneck = N > 500
        elif bottleneck == "true":
            self.use_bottleneck = True
        else:
            self.use_bottleneck = False

        if self.use_bottleneck:
            hidden_dim = max(2, N // 4)
            self.fcc1 = torch.nn.Linear(N, hidden_dim)
            self.fcc2 = torch.nn.Linear(hidden_dim, hidden_dim)
            self.fcc3 = torch.nn.Linear(hidden_dim, N)

            print(f"✅ Using bottleneck: {N} -> {hidden_dim} -> {hidden_dim} -> {N}")

            torch.nn.init.xavier_uniform_(self.fcc1.weight)
            torch.nn.init.xavier_uniform_(self.fcc2.weight)
            torch.nn.init.xavier_uniform_(self.fcc3.weight)
        else:
            self.fcc1 = torch.nn.Linear(N, N)
            self.fcc2 = torch.nn.Linear(N, N)

            print(f"✅ Using full architecture: {N} -> {N} -> {N}")

            torch.nn.init.xavier_uniform_(self.fcc1.weight)
            torch.nn.init.xavier_uniform_(self.fcc2.weight)

    def forward(self, t, y):
        if self.use_bottleneck:
            out = self.activation(self.fcc1(y))
            out = self.activation(self.fcc2(out))
            out = self.fcc3(out)
        else:
            out = self.activation(self.fcc1(y))
            out = self.fcc2(out)

        # ============================================================
        # Core formula fix:
        # wrong old formula: dy = y * (out - sum(y))
        # fixed formula:     dy = y * (out - sum(y * out))
        # ============================================================
        f = torch.sum(y * out, dim=1, keepdim=True)
        return y * (out - f)


def make_time_points(time_mode: str, device: str) -> torch.Tensor:
    if time_mode == "endpoint":
        return torch.tensor([0.0, 0.99], dtype=torch.float32, device=device)

    full_t = torch.arange(0.0, 100.0, 0.01, dtype=torch.float32)
    return full_t[:100].to(device)


def ode_solve(func, y0, time_points, args):
    if args.ode_method in {"rk4", "euler", "midpoint"}:
        return odeint(
            func,
            y0,
            time_points,
            method=args.ode_method,
            options={"step_size": args.ode_step_size},
        )

    return odeint(
        func,
        y0,
        time_points,
        method="dopri5",
        rtol=args.rtol,
        atol=args.atol,
    )


def get_batch(ztrn, ptrn, mb_size: int, device: str):
    num_samples = ptrn.size(dim=0)
    actual_mb = min(mb_size, num_samples)

    if actual_mb == num_samples:
        sample_ids = torch.from_numpy(np.arange(num_samples, dtype=np.int64))
    else:
        sample_ids = torch.from_numpy(
            np.random.choice(np.arange(num_samples, dtype=np.int64), actual_mb, replace=False)
        )

    return ztrn[sample_ids, :].to(device), ptrn[sample_ids, :].to(device)


def bray_curtis_loss_batch(pred: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
    num = torch.sum(torch.abs(pred - target), dim=1)
    den = torch.sum(torch.abs(pred + target), dim=1) + 1e-8
    return torch.mean(num / den)


def write_loss_tsv(path: str, rows: List[Tuple[int, float]], value_name: str):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["epoch", value_name])
        writer.writerows(rows)


def predict_in_chunks(func, data_z, time_points, chunk_size, device, args):
    preds = []
    num_samples = data_z.size(0)

    t_infer = torch.tensor(
        [float(time_points[0].item()), float(time_points[-1].item())],
        dtype=torch.float32,
        device=device,
    )

    func.eval()
    with torch.no_grad():
        for i in range(0, num_samples, chunk_size):
            chunk_z = data_z[i:i + chunk_size].to(device)
            pred_full = ode_solve(func, chunk_z, t_infer, args)
            pred = pred_full[-1, :, :].detach().cpu().numpy()
            preds.append(pred)

            del chunk_z, pred_full
            if torch.cuda.is_available():
                torch.cuda.empty_cache()

    return np.vstack(preds)


def train_model(ztrn, ptrn, zval, pval, ztst, zall, N, args, device):
    time_points = make_time_points(args.time_mode, device)

    func = ODEFunc(N=N, activation=args.activation, bottleneck=args.bottleneck).to(device)

    optimizer = torch.optim.Adam(func.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    scheduler = ReduceLROnPlateau(optimizer, mode="min", factor=0.5, patience=50)

    best_validation_loss = float("inf")
    best_state = None
    patience_counter = 0
    early_stop_triggered = False
    epoch_stopped = -1

    train_losses = []
    val_losses = []

    print(f"🚀 Device: {device}")
    print(f"🚀 Start training: N={N}, parameters={sum(p.numel() for p in func.parameters())}")
    print(f"🚀 time_mode={args.time_mode}, ode_method={args.ode_method}, step_size={args.ode_step_size}")
    print(f"🚀 activation={args.activation}, bottleneck={args.bottleneck}")

    for epoch in range(args.max_epochs):
        epoch_stopped = epoch

        func.train()
        optimizer.zero_grad(set_to_none=True)

        batch_z, batch_p = get_batch(ztrn, ptrn, args.batch_size, device)
        pred_full = ode_solve(func, batch_z, time_points, args)
        pred = pred_full[-1, :, :]

        train_loss = bray_curtis_loss_batch(pred, batch_p)
        train_loss.backward()
        optimizer.step()

        train_loss_value = float(train_loss.item())
        train_losses.append((epoch, train_loss_value))

        func.eval()
        val_total = 0.0
        val_n = zval.size(0)
        val_batch_size = 16

        with torch.no_grad():
            for i in range(0, val_n, val_batch_size):
                batch_z_val = zval[i:i + val_batch_size].to(device)
                batch_p_val = pval[i:i + val_batch_size].to(device)

                pred_val_full = ode_solve(func, batch_z_val, time_points, args)
                pred_val = pred_val_full[-1, :, :]

                batch_loss = bray_curtis_loss_batch(pred_val, batch_p_val)
                val_total += float(batch_loss.item()) * batch_z_val.size(0)

        validation_loss = val_total / max(1, val_n)
        val_losses.append((epoch, validation_loss))
        scheduler.step(validation_loss)

        if epoch % 100 == 0 or epoch == args.max_epochs - 1:
            print(
                f"Epoch {epoch}/{args.max_epochs}: "
                f"Train Loss={train_loss_value:.6f}, Validation Loss={validation_loss:.6f}",
                flush=True,
            )

        if validation_loss <= best_validation_loss:
            best_validation_loss = validation_loss
            best_state = copy.deepcopy(func.state_dict())
            patience_counter = 0
        else:
            patience_counter += 1

        if patience_counter >= args.early_stop_patience:
            early_stop_triggered = True
            print(f"🛑 Early stopping at epoch {epoch}: no validation improvement for {args.early_stop_patience} epochs.")
            break

    if best_state is None:
        raise RuntimeError("No best model state was saved.")

    func.load_state_dict(best_state)
    func.eval()

    model_path = os.path.join(args.output_dir, "best_model.pth")
    torch.save(best_state, model_path)

    write_loss_tsv(os.path.join(args.output_dir, "training_loss.tsv"), train_losses, "train_loss")
    write_loss_tsv(os.path.join(args.output_dir, "validation_loss.tsv"), val_losses, "validation_loss")

    print("🔮 Predicting final matrices...")
    qtst_raw = predict_in_chunks(func, ztst, time_points, args.chunk_size, device, args)
    qtrn_raw = predict_in_chunks(func, zall, time_points, args.chunk_size, device, args)

    return {
        "qtst_raw": qtst_raw,
        "qtrn_raw": qtrn_raw,
        "best_validation_loss": float(best_validation_loss),
        "final_train_loss": float(train_losses[-1][1]) if train_losses else None,
        "epoch_stopped": int(epoch_stopped),
        "early_stop_triggered": bool(early_stop_triggered),
        "used_bottleneck": bool(func.use_bottleneck),
        "model_path": model_path,
    }


def infer_project(train_file: str, output_dir: str, explicit_project: str = None) -> str:
    if explicit_project:
        return explicit_project

    train_abs = os.path.abspath(train_file)
    parent = os.path.basename(os.path.dirname(os.path.dirname(train_abs)))
    if parent:
        return parent

    out_parent = os.path.basename(os.path.dirname(os.path.abspath(output_dir)))
    return out_parent or "unknown"


def split_train_validation(P: np.ndarray, seed: int):
    n_cols = P.shape[1]
    if n_cols < 2:
        raise ValueError("Ptrain.csv must contain at least 2 samples.")

    val_count = max(1, int(round(0.2 * n_cols)))
    val_count = min(val_count, n_cols - 1)

    rng = np.random.default_rng(seed)
    val_idx = np.sort(rng.choice(np.arange(n_cols), size=val_count, replace=False))

    train_mask = np.ones(n_cols, dtype=bool)
    train_mask[val_idx] = False

    return P[:, train_mask], P[:, val_idx]


def main():
    args = parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    set_seed(args.seed)
    device = resolve_device(args.device)

    project = infer_project(args.train_file, args.output_dir, args.project)

    print(f"Project: {project}")
    print(f"Train file: {args.train_file}")
    print(f"Test file: {args.test_file}")
    print(f"Output dir: {args.output_dir}")

    P_all = load_matrix(args.train_file)
    P_test = load_matrix(args.test_file, expected_rows=P_all.shape[0])

    if P_all.shape[0] != P_test.shape[0]:
        raise ValueError(
            f"Train/test species dimension mismatch: "
            f"Ptrain rows={P_all.shape[0]}, Ptest rows={P_test.shape[0]}"
        )

    P_train, P_val = split_train_validation(P_all, args.seed)

    ptrn, ztrn = process_data(P_train, args.normalize_z)
    pval, zval = process_data(P_val, args.normalize_z)
    pall, zall = process_data(P_all, args.normalize_z)
    ptst, ztst = process_data(P_test, args.normalize_z)

    _, N = ptrn.shape

    print(
        f"Dimensions: N_species={N}, "
        f"M_train_samples={P_train.shape[1]}, "
        f"M_val_samples={P_val.shape[1]}, "
        f"M_test_perturbations={P_test.shape[1]}",
        flush=True,
    )

    result = train_model(
        ztrn=ztrn,
        ptrn=ptrn,
        zval=zval,
        pval=pval,
        ztst=ztst,
        zall=zall,
        N=N,
        args=args,
        device=device,
    )

    qtst_raw = result["qtst_raw"]
    qtrn_raw = result["qtrn_raw"]

    qtst, qtst_metrics = normalize_prediction_rows(qtst_raw)
    qtrn, qtrn_metrics = normalize_prediction_rows(qtrn_raw)

    qtst_raw_sums = raw_row_sum_summary(qtst_raw)
    qtrn_raw_sums = raw_row_sum_summary(qtrn_raw)

    np.savetxt(os.path.join(args.output_dir, "qtst_raw.csv"), qtst_raw, delimiter=",", fmt="%.10g")
    np.savetxt(os.path.join(args.output_dir, "qtrn_raw.csv"), qtrn_raw, delimiter=",", fmt="%.10g")
    np.savetxt(os.path.join(args.output_dir, "qtst.csv"), qtst, delimiter=",", fmt="%.10g")
    np.savetxt(os.path.join(args.output_dir, "qtrn.csv"), qtrn, delimiter=",", fmt="%.10g")

    metrics = {
        "project": project,
        "N_species": int(N),
        "M_train_samples": int(P_train.shape[1]),
        "M_val_samples": int(P_val.shape[1]),
        "M_test_samples": int(P_test.shape[1]),
        "best_validation_loss": result["best_validation_loss"],
        "final_train_loss": result["final_train_loss"],
        "epoch_stopped": result["epoch_stopped"],
        "early_stop_triggered": result["early_stop_triggered"],
        "device": device,
        "activation": args.activation,
        "bottleneck": args.bottleneck,
        "used_bottleneck": result["used_bottleneck"],
        "normalize_z": bool(args.normalize_z),
        "seed": int(args.seed),
        "time_mode": args.time_mode,
        "ode_method": args.ode_method,
        "ode_step_size": args.ode_step_size,
        "rtol": args.rtol,
        "atol": args.atol,
        "chunk_size": int(args.chunk_size),
        "model_path": result["model_path"],

        "qtst_raw_row_sum_min": qtst_raw_sums["min"],
        "qtst_raw_row_sum_median": qtst_raw_sums["median"],
        "qtst_raw_row_sum_max": qtst_raw_sums["max"],
        "qtrn_raw_row_sum_min": qtrn_raw_sums["min"],
        "qtrn_raw_row_sum_median": qtrn_raw_sums["median"],
        "qtrn_raw_row_sum_max": qtrn_raw_sums["max"],

        "qtst_row_sum_min": qtst_metrics["row_sum_min"],
        "qtst_row_sum_median": qtst_metrics["row_sum_median"],
        "qtst_row_sum_max": qtst_metrics["row_sum_max"],
        "qtrn_row_sum_min": qtrn_metrics["row_sum_min"],
        "qtrn_row_sum_median": qtrn_metrics["row_sum_median"],
        "qtrn_row_sum_max": qtrn_metrics["row_sum_max"],

        "qtst_min_value": qtst_metrics["min_value"],
        "qtrn_min_value": qtrn_metrics["min_value"],
        "qtst_zero_rows_after_normalization": qtst_metrics["zero_rows_after_normalization"],
        "qtrn_zero_rows_after_normalization": qtrn_metrics["zero_rows_after_normalization"],
    }

    with open(os.path.join(args.output_dir, "model_metrics.json"), "w") as handle:
        json.dump(metrics, handle, indent=2, sort_keys=True, allow_nan=False)

    print(f"✅ Saved outputs to: {args.output_dir}")
    print("✅ Formula fixed: dy = y * (out - sum(y * out))")
    print(f"✅ ODE solver: {args.ode_method}")
    print("✅ qtst.csv and qtrn.csv are non-negative and row-normalized.")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
