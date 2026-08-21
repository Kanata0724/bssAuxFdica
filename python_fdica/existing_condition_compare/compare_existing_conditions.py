"""Compare MATLAB/Python AuxFDICA using existing main.m/main.py conditions."""

from __future__ import annotations

import csv
import json
import subprocess
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.io
import scipy.signal
import torch

from python_fdica.bssAuxFdica import local_calcFdicaCost, local_projectionBack
from python_fdica.permSolverCor import permSolverCor
from python_fdica.permSolverIps import permSolverIps
from python_fdica.stft import dgt_istft, dgt_stft


REPO = Path(__file__).resolve().parents[2]
TOOL_DIR = REPO / "python_fdica" / "existing_condition_compare"
OUT_DIR = REPO / "output" / "fdica_existing_conditions_compare"


def run_matlab(script: Path) -> None:
    subprocess.run(["matlab", "-batch", f"run('{script.as_posix()}')"], cwd=REPO, check=True)


def to_python(value: Any) -> Any:
    if isinstance(value, np.void) and value.dtype.names:
        return {name: to_python(value[name]) for name in value.dtype.names}
    if isinstance(value, np.ndarray) and value.dtype.names:
        if value.size == 1:
            item = value.item()
            if isinstance(item, tuple):
                return {name: to_python(field) for name, field in zip(value.dtype.names, item, strict=True)}
            return to_python(item)
        return {name: to_python(value[name]) for name in value.dtype.names}
    if isinstance(value, np.ndarray) and value.dtype.kind in {"U", "S"}:
        return "".join(value.reshape(-1).astype(str).tolist())
    if isinstance(value, np.ndarray) and value.size == 1:
        item = value.item()
        return item.decode() if isinstance(item, bytes) else item
    return value


def load_params(mat: dict[str, Any]) -> dict[str, Any]:
    raw = to_python(mat["params"])
    return {
        "dataNo": int(raw["dataNo"]),
        "seed": int(raw["seed"]),
        "fftSize": int(raw["fftSize"]),
        "shiftSize": int(raw["shiftSize"]),
        "nSrc": int(raw["nSrc"]),
        "nIter": int(raw["nIter"]),
        "isWhiten": bool(raw["isWhiten"]),
        "srcModel": str(raw["srcModel"]),
        "refMic": int(raw["refMic"]),
        "permSolver": str(raw["permSolver"]),
        "isDraw": bool(raw["isDraw"]),
        "isFilt": bool(raw["isFilt"]),
        "isPowRatio": bool(raw["isPowRatio"]),
        "typeCor": str(raw["typeCor"]),
        "deltaFreq": int(raw["deltaFreq"]),
        "ratioFreq": int(raw["ratioFreq"]),
    }


def array_stats(a: np.ndarray, b: np.ndarray) -> dict[str, float]:
    diff = np.asarray(a) - np.asarray(b)
    abs_diff = np.abs(diff)
    denom = np.maximum(np.abs(a), np.finfo(float).eps)
    return {
        "max_abs_error": float(np.max(abs_diff)),
        "mean_abs_error": float(np.mean(abs_diff)),
        "sum_abs_error": float(np.sum(abs_diff)),
        "sum_squared_abs_error": float(np.sum(abs_diff**2)),
        "rmse": float(np.sqrt(np.mean(abs_diff**2))),
        "max_relative_error": float(np.max(abs_diff / denom)),
        "mean_relative_error": float(np.mean(abs_diff / denom)),
    }


def complex_metrics(a: np.ndarray, b: np.ndarray, *, align_scale: bool = False) -> dict[str, Any]:
    a_arr = np.asarray(a)
    b_arr = np.asarray(b)
    if align_scale:
        denom = np.vdot(b_arr.reshape(-1), b_arr.reshape(-1))
        gain = np.vdot(b_arr.reshape(-1), a_arr.reshape(-1)) / denom if abs(denom) > 0 else 1.0
        b_cmp = gain * b_arr
    else:
        gain = 1.0 + 0.0j
        b_cmp = b_arr
    diff = a_arr - b_cmp
    abs_diff = np.abs(diff)
    a_norm2 = float(np.sum(np.abs(a_arr) ** 2))
    b_norm2 = float(np.sum(np.abs(b_cmp) ** 2))
    corr_denom = np.sqrt(a_norm2 * b_norm2)
    corr = np.vdot(a_arr.reshape(-1), b_cmp.reshape(-1)) / corr_denom if corr_denom > 0 else 0.0
    return {
        "sum_squared_abs_error": float(np.sum(abs_diff**2)),
        "rmse": float(np.sqrt(np.mean(abs_diff**2))),
        "mae": float(np.mean(abs_diff)),
        "max_abs_error": float(np.max(abs_diff)),
        "normalized_rmse": float(np.sqrt(np.sum(abs_diff**2) / max(a_norm2, np.finfo(float).eps))),
        "complex_correlation_abs": float(abs(corr)),
        "complex_correlation_real": float(np.real(corr)),
        "complex_correlation_imag": float(np.imag(corr)),
        "alignment_gain_real": float(np.real(gain)),
        "alignment_gain_imag": float(np.imag(gain)),
    }


def sourcewise_complex_metrics(a: np.ndarray, b: np.ndarray, *, align_scale: bool = False) -> list[dict[str, Any]]:
    return [
        complex_metrics(np.asarray(a)[:, :, source], np.asarray(b)[:, :, source], align_scale=align_scale)
        for source in range(np.asarray(a).shape[2])
    ]


def waveform_metrics(a: np.ndarray, b: np.ndarray, *, align_scale: bool = False, align_delay: bool = False) -> list[dict[str, float]]:
    a_arr = np.asarray(a, dtype=np.float64)
    b_arr = np.asarray(b, dtype=np.float64)
    out: list[dict[str, float]] = []
    for source in range(a_arr.shape[1]):
        x = a_arr[:, source]
        y = b_arr[:, source]
        lag = 0
        if align_delay:
            corr_full = scipy.signal.correlate(y, x, mode="full", method="fft")
            lags = scipy.signal.correlation_lags(y.size, x.size, mode="full")
            lag = int(lags[int(np.argmax(np.abs(corr_full)))])
            if lag > 0:
                y_cmp = y[lag:]
                x_cmp = x[: y_cmp.size]
            elif lag < 0:
                x_cmp = x[-lag:]
                y_cmp = y[: x_cmp.size]
            else:
                x_cmp = x
                y_cmp = y
        else:
            x_cmp = x
            y_cmp = y
        if align_scale:
            denom = float(np.dot(y_cmp, y_cmp))
            scale = float(np.dot(y_cmp, x_cmp) / denom) if denom > 0 else 1.0
        else:
            scale = 1.0
        y_scaled = scale * y_cmp
        diff = x_cmp - y_scaled
        x_energy = float(np.sum(x_cmp**2))
        y_center = y_scaled - np.mean(y_scaled)
        x_center = x_cmp - np.mean(x_cmp)
        corr_denom = float(np.linalg.norm(x_center) * np.linalg.norm(y_center))
        pearson = float(np.dot(x_center, y_center) / corr_denom) if corr_denom > 0 else 0.0
        ncc_denom = float(np.linalg.norm(x_cmp) * np.linalg.norm(y_scaled))
        ncc = float(np.dot(x_cmp, y_scaled) / ncc_denom) if ncc_denom > 0 else 0.0
        out.append(
            {
                "sum_squared_error": float(np.sum(diff**2)),
                "rmse": float(np.sqrt(np.mean(diff**2))),
                "mae": float(np.mean(np.abs(diff))),
                "max_abs_error": float(np.max(np.abs(diff))),
                "normalized_rmse": float(np.sqrt(np.sum(diff**2) / max(x_energy, np.finfo(float).eps))),
                "pearson_correlation": pearson,
                "normalized_cross_correlation": ncc,
                "alignment_scale": scale,
                "alignment_delay_samples": float(lag),
            }
        )
    return out


def cost_stats(matlab_cost: np.ndarray, python_cost: np.ndarray) -> dict[str, float]:
    m = matlab_cost.reshape(-1)
    p = python_cost.reshape(-1)
    abs_diff = np.abs(m - p)
    return {
        "matlab_initial": float(m[0]),
        "python_initial": float(p[0]),
        "matlab_final": float(m[-1]),
        "python_final": float(p[-1]),
        "sum_abs_error": float(np.sum(abs_diff)),
        "sum_squared_abs_error": float(np.sum(abs_diff**2)),
        "mean_abs_error": float(np.mean(abs_diff)),
        "rmse": float(np.sqrt(np.mean(abs_diff**2))),
        "max_abs_error": float(np.max(abs_diff)),
        "max_relative_error": float(np.max(abs_diff / np.maximum(np.abs(m), np.finfo(float).eps))),
    }


def first_noticeable_cost_error(matlab_cost: np.ndarray, python_cost: np.ndarray) -> dict[str, Any]:
    m = matlab_cost.reshape(-1)
    p = python_cost.reshape(-1)
    abs_diff = np.abs(m - p)
    checks = {}
    for threshold in (1e-6, 1e-3, 1.0, 10.0, 100.0, 1000.0):
        idx = np.where(abs_diff > threshold)[0]
        checks[str(threshold)] = int(idx[0]) if idx.size else None
    return checks


def local_aux_fdica_current_traced(
    x: torch.Tensor, n_iter: int, src_model: str
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, dict[str, torch.Tensor]]:
    n_freq, n_frame, n_ch = x.shape
    w = torch.eye(n_ch, dtype=x.dtype, device=x.device).repeat(n_freq, 1, 1)
    y = x.clone()
    cost = torch.zeros(n_iter + 1, dtype=x.real.dtype, device=x.device)
    cost[0] = local_calcFdicaCost(y, w, src_model)
    threshold = 10000 * torch.finfo(x.real.dtype).eps
    eye = torch.eye(n_ch, dtype=x.dtype, device=x.device)
    debug: dict[str, torch.Tensor] = {
        "Y0": y.clone(),
        "W0": w.permute(1, 2, 0).contiguous().clone(),
    }
    for i_iter in range(n_iter):
        radius = torch.abs(y) if src_model == "LAP" else torch.abs(y) ** 2
        inv_radius = radius.clamp_min(threshold).reciprocal()
        for n in range(n_ch):
            for f in range(n_freq):
                xf = x[f].mT
                vn = (xf * inv_radius[f, :, n][None, :]) @ xf.mH / n_frame
                mean_power = torch.real(torch.trace(vn)) / n_ch
                ridge = 100 * torch.finfo(x.real.dtype).eps * mean_power.clamp_min(1.0)
                vn_reg = vn + ridge * eye
                system = w[f] @ vn_reg
                try:
                    wn = torch.linalg.solve(system, eye[:, n])
                except RuntimeError as exc:
                    if "singular" not in str(exc).lower():
                        raise
                    wn = torch.linalg.pinv(system) @ eye[:, n]
                norm = torch.sqrt(torch.real(wn.conj() @ vn_reg @ wn).clamp_min(threshold))
                wn = wn / norm
                w[f, n] = wn.conj()
                y[f, :, n] = x[f] @ wn.conj()
        cost[i_iter + 1] = local_calcFdicaCost(y, w, src_model)
        if i_iter == 0:
            debug["Y_after_iter1"] = y.clone()
            debug["W_after_iter1"] = w.permute(1, 2, 0).contiguous().clone()
        if i_iter == 9:
            debug["Y_after_iter10"] = y.clone()
            debug["W_after_iter10"] = w.permute(1, 2, 0).contiguous().clone()
        if i_iter == 24:
            debug["Y_after_iter25"] = y.clone()
            debug["W_after_iter25"] = w.permute(1, 2, 0).contiguous().clone()
    debug["Y_final"] = y.clone()
    debug["W_final"] = w.permute(1, 2, 0).contiguous().clone()
    return y, w.permute(1, 2, 0).contiguous(), cost, debug


def local_aux_fdica_matlab_like(
    x: torch.Tensor, n_iter: int, src_model: str
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, dict[str, torch.Tensor]]:
    n_freq, n_frame, n_ch = x.shape
    w = torch.eye(n_ch, dtype=x.dtype, device=x.device).repeat(n_freq, 1, 1)
    y = x.clone()
    cost = torch.zeros(n_iter + 1, dtype=x.real.dtype, device=x.device)
    cost[0] = local_calcFdicaCost(y, w, src_model)
    threshold = 10000 * torch.finfo(x.real.dtype).eps
    eye = torch.eye(n_ch, dtype=x.dtype, device=x.device)
    debug: dict[str, torch.Tensor] = {
        "Y0": y.clone(),
        "W0": w.permute(1, 2, 0).contiguous().clone(),
    }
    for i_iter in range(n_iter):
        radius = torch.abs(y) if src_model == "LAP" else torch.abs(y) ** 2
        inv_radius = radius.clamp_min(threshold).reciprocal()
        for n in range(n_ch):
            for f in range(n_freq):
                xf = x[f].mT
                vn = (xf * inv_radius[f, :, n][None, :]) @ xf.mH / n_frame
                wn = torch.linalg.solve(w[f] @ vn, eye[:, n])
                norm = torch.sqrt(wn.conj() @ vn @ wn)
                wn = wn / norm
                if i_iter == 0 and n == 0 and f == 0:
                    debug["V_first"] = vn.clone()
                    debug["wn_first"] = wn[:, None].clone()
                    debug["norm_first"] = norm.reshape(1, 1).clone()
                w[f, n] = wn.conj()
                y[f, :, n] = x[f] @ wn.conj()
        cost[i_iter + 1] = local_calcFdicaCost(y, w, src_model)
        if i_iter == 0:
            debug["Y_after_iter1"] = y.clone()
            debug["W_after_iter1"] = w.permute(1, 2, 0).contiguous().clone()
        if i_iter == 9:
            debug["Y_after_iter10"] = y.clone()
            debug["W_after_iter10"] = w.permute(1, 2, 0).contiguous().clone()
        if i_iter == 24:
            debug["Y_after_iter25"] = y.clone()
            debug["W_after_iter25"] = w.permute(1, 2, 0).contiguous().clone()
    debug["Y_final"] = y.clone()
    debug["W_final"] = w.permute(1, 2, 0).contiguous().clone()
    return y, w.permute(1, 2, 0).contiguous(), cost, debug


def run_python_fdica(
    obs_sig: torch.Tensor,
    obs_spec: torch.Tensor,
    ref_spec: torch.Tensor,
    src_spec: torch.Tensor,
    params: dict[str, Any],
    *,
    matlab_like: bool,
) -> dict[str, torch.Tensor]:
    if params["isWhiten"]:
        raise RuntimeError("This comparison tool expects existing main conditions isWhiten=false.")
    obs_spec_input = obs_spec[:, :, : params["nSrc"]]
    if matlab_like:
        y, w, cost, debug = local_aux_fdica_matlab_like(obs_spec_input, params["nIter"], params["srcModel"])
    else:
        y, w, cost, debug = local_aux_fdica_current_traced(obs_spec_input, params["nIter"], params["srcModel"])
    fixed, demix_fixed = local_projectionBack(y, ref_spec, w)
    est_none = dgt_istft(fixed, params["fftSize"], params["shiftSize"], length=obs_sig.shape[0])
    est_cor_spec, perm_cor = permSolverCor(
        fixed, params["isPowRatio"], params["typeCor"], params["deltaFreq"], params["ratioFreq"]
    )
    est_cor = dgt_istft(est_cor_spec, params["fftSize"], params["shiftSize"], length=obs_sig.shape[0])
    est_ips_spec, perm_ips = permSolverIps(fixed, src_spec)
    est_ips = dgt_istft(est_ips_spec, params["fftSize"], params["shiftSize"], length=obs_sig.shape[0])
    return {
        "obsSpecInput": obs_spec_input,
        "Y": y,
        "W": w,
        "cost": cost,
        "fixed": fixed,
        "demixFixed": demix_fixed,
        "estNone": est_none,
        "estCor": est_cor,
        "estIps": est_ips,
        "estCorSpec": est_cor_spec,
        "estIpsSpec": est_ips_spec,
        "permCor": perm_cor,
        "permIps": perm_ips,
        **debug,
    }


def extract_metric_table(metrics: Any) -> dict[str, dict[str, list[float]]]:
    raw = to_python(metrics)
    table: dict[str, dict[str, list[float]]] = {}
    for method, values in raw.items():
        if isinstance(values, dict):
            table[method] = {}
            for key in ("sdr", "sir", "sar", "perm"):
                if key in values:
                    table[method][key] = np.asarray(values[key]).reshape(-1).astype(float).tolist()
    return table


def permutation_summary(matlab_perm_one_based: np.ndarray, python_perm_zero_based: np.ndarray) -> dict[str, Any]:
    matlab_perm = np.asarray(matlab_perm_one_based, dtype=int) - 1
    python_perm = np.asarray(python_perm_zero_based, dtype=int)
    diff_bins = np.where(np.any(matlab_perm != python_perm, axis=1))[0]
    return {
        "n_frequency_bins": int(matlab_perm.shape[0]),
        "n_mismatch_bins": int(len(diff_bins)),
        "mismatch_bins_zero_based": diff_bins.astype(int).tolist(),
        "matlab_perm_zero_based_at_mismatch": matlab_perm[diff_bins].astype(int).tolist(),
        "python_perm_zero_based_at_mismatch": python_perm[diff_bins].astype(int).tolist(),
    }


def cor_mismatch_energy(
    fixed_spec: np.ndarray,
    matlab_perm_one_based: np.ndarray,
    python_perm_zero_based: np.ndarray,
    fs: float,
    fft_size: int,
) -> dict[str, Any]:
    matlab_perm = np.asarray(matlab_perm_one_based, dtype=int) - 1
    python_perm = np.asarray(python_perm_zero_based, dtype=int)
    diff_bins = np.where(np.any(matlab_perm != python_perm, axis=1))[0]
    bin_energy = np.sum(np.abs(fixed_spec) ** 2, axis=(1, 2))
    total_energy = float(np.sum(bin_energy))
    mismatch_energy = float(np.sum(bin_energy[diff_bins]))
    rows = [
        {
            "bin_zero_based": int(idx),
            "frequency_hz": float(idx * fs / fft_size),
            "matlab_perm_zero_based": matlab_perm[idx].astype(int).tolist(),
            "python_perm_zero_based": python_perm[idx].astype(int).tolist(),
            "bin_energy": float(bin_energy[idx]),
            "bin_energy_ratio": float(bin_energy[idx] / max(total_energy, np.finfo(float).eps)),
        }
        for idx in diff_bins
    ]
    return {
        "total_energy": total_energy,
        "mismatch_energy": mismatch_energy,
        "mismatch_energy_ratio": float(mismatch_energy / max(total_energy, np.finfo(float).eps)),
        "rows": rows,
    }


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    run_matlab(TOOL_DIR / "export_matlab_existing_debug.m")
    mat = scipy.io.loadmat(OUT_DIR / "matlab_existing_debug.mat", squeeze_me=False)
    params = load_params(mat)

    obs_sig_np = np.asarray(mat["obsSig"], dtype=np.float64)
    src_sig_np = np.asarray(mat["srcSig"], dtype=np.float64)
    obs_sig = torch.as_tensor(obs_sig_np, dtype=torch.float64)
    src_sig = torch.as_tensor(src_sig_np, dtype=torch.float64)
    torch.manual_seed(params["seed"])

    obs_spec = dgt_stft(obs_sig, params["fftSize"], params["shiftSize"])
    ref_spec = obs_spec[:, :, params["refMic"] - 1]
    src_spec = dgt_stft(src_sig[:, params["refMic"] - 1, :], params["fftSize"], params["shiftSize"])

    current = run_python_fdica(obs_sig, obs_spec, ref_spec, src_spec, params, matlab_like=False)
    matlab_like = run_python_fdica(obs_sig, obs_spec, ref_spec, src_spec, params, matlab_like=True)

    scipy.io.savemat(
        OUT_DIR / "python_existing_debug.mat",
        {
            "estSigCurrentNone": current["estNone"].numpy(),
            "estSigCurrentCor": current["estCor"].numpy(),
            "estSigCurrentIps": current["estIps"].numpy(),
            "estSigMatlabLikeNone": matlab_like["estNone"].numpy(),
            "estSigMatlabLikeCor": matlab_like["estCor"].numpy(),
            "estSigMatlabLikeIps": matlab_like["estIps"].numpy(),
            "costCurrent": current["cost"].numpy(),
            "costMatlabLike": matlab_like["cost"].numpy(),
        },
    )
    run_matlab(TOOL_DIR / "evaluate_python_existing_outputs.m")
    eval_mat = scipy.io.loadmat(OUT_DIR / "python_existing_bss_eval.mat", squeeze_me=False)

    matlab_cost = np.asarray(mat["cost"]).reshape(-1)
    current_cost = current["cost"].numpy().reshape(-1)
    matlab_like_cost = matlab_like["cost"].numpy().reshape(-1)

    rows = []
    for i, (m, c, ml) in enumerate(zip(matlab_cost, current_cost, matlab_like_cost)):
        rows.append(
            {
                "iteration": i,
                "matlab_cost": float(m),
                "python_current_cost": float(c),
                "python_matlab_like_cost": float(ml),
                "current_abs_error": float(abs(m - c)),
                "current_squared_abs_error": float(abs(m - c) ** 2),
                "current_relative_error": float(abs(m - c) / max(abs(m), np.finfo(float).eps)),
                "matlab_like_abs_error": float(abs(m - ml)),
                "matlab_like_squared_abs_error": float(abs(m - ml) ** 2),
                "matlab_like_relative_error": float(abs(m - ml) / max(abs(m), np.finfo(float).eps)),
            }
        )
    with (OUT_DIR / "existing_condition_cost_comparison.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fs = float(np.asarray(mat["fs"]).reshape(-1)[0])
    cor_energy = cor_mismatch_energy(
        np.asarray(mat["estSpecFixed"]).squeeze(),
        np.asarray(mat["permCor"]),
        current["permCor"].numpy(),
        fs,
        params["fftSize"],
    )
    with (OUT_DIR / "cor_mismatch_bins.csv").open("w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "bin_zero_based",
            "frequency_hz",
            "matlab_perm_zero_based",
            "python_perm_zero_based",
            "bin_energy",
            "bin_energy_ratio",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(cor_energy["rows"])

    plt.figure(figsize=(9, 5))
    x_axis = [row["iteration"] for row in rows]
    plt.plot(x_axis, [row["matlab_cost"] for row in rows], label="MATLAB")
    plt.plot(x_axis, [row["python_current_cost"] for row in rows], "--", label="Python current")
    plt.plot(x_axis, [row["python_matlab_like_cost"] for row in rows], ":", label="Python MATLAB-like")
    plt.xlabel("Iteration")
    plt.ylabel("Cost")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT_DIR / "existing_condition_cost_comparison.png", dpi=160)
    plt.close()

    debug = to_python(mat["debug"])
    matlab_y = np.asarray(mat["estSpecFdica"])
    matlab_fixed = np.asarray(mat["estSpecFixed"]).squeeze()
    matlab_cor_perm = np.asarray(mat["permCor"], dtype=int) - 1
    matlab_ips_perm = np.asarray(mat["permIps"], dtype=int) - 1
    matlab_cor_spec = np.empty_like(matlab_fixed)
    matlab_ips_spec = np.empty_like(matlab_fixed)
    for f in range(matlab_fixed.shape[0]):
        matlab_cor_spec[f] = matlab_fixed[f][:, matlab_cor_perm[f]]
        matlab_ips_spec[f] = matlab_fixed[f][:, matlab_ips_perm[f]]
    matlab_est_none = np.asarray(mat["estSigNone"])
    matlab_est_cor = np.asarray(mat["estSigCor"])
    matlab_est_ips = np.asarray(mat["estSigIps"])
    internal_stats = {
        "stft_obsSpec": array_stats(np.asarray(mat["obsSpec"]), obs_spec.numpy()),
        "obsSpecInput": array_stats(np.asarray(mat["obsSpecInput"]), current["obsSpecInput"].numpy()),
        "matlab_like_Y_after_iter1": array_stats(np.asarray(debug["Y_after_iter1"]), matlab_like["Y_after_iter1"].numpy()),
        "matlab_like_W_after_iter1": array_stats(np.asarray(debug["W_after_iter1"]), matlab_like["W_after_iter1"].numpy()),
        "current_Y_after_fdica": array_stats(np.asarray(mat["estSpecFdica"]), current["Y"].numpy()),
        "matlab_like_Y_after_fdica": array_stats(np.asarray(mat["estSpecFdica"]), matlab_like["Y"].numpy()),
        "current_projection_back": array_stats(np.asarray(mat["estSpecFixed"]).squeeze(), current["fixed"].numpy()),
        "matlab_like_projection_back": array_stats(np.asarray(mat["estSpecFixed"]).squeeze(), matlab_like["fixed"].numpy()),
        "current_estSig_none": array_stats(np.asarray(mat["estSigNone"]), current["estNone"].numpy()),
        "matlab_like_estSig_none": array_stats(np.asarray(mat["estSigNone"]), matlab_like["estNone"].numpy()),
    }

    stage_metrics = {
        "aux_fdica_raw": {
            "current_before_alignment": sourcewise_complex_metrics(matlab_y, current["Y"].numpy(), align_scale=False),
            "current_after_complex_scale": sourcewise_complex_metrics(matlab_y, current["Y"].numpy(), align_scale=True),
            "matlab_like_before_alignment": sourcewise_complex_metrics(matlab_y, matlab_like["Y"].numpy(), align_scale=False),
            "matlab_like_after_complex_scale": sourcewise_complex_metrics(matlab_y, matlab_like["Y"].numpy(), align_scale=True),
        },
        "projection_back_before_permutation": {
            "current_before_alignment": sourcewise_complex_metrics(matlab_fixed, current["fixed"].numpy(), align_scale=False),
            "current_after_complex_scale": sourcewise_complex_metrics(matlab_fixed, current["fixed"].numpy(), align_scale=True),
            "matlab_like_before_alignment": sourcewise_complex_metrics(matlab_fixed, matlab_like["fixed"].numpy(), align_scale=False),
            "matlab_like_after_complex_scale": sourcewise_complex_metrics(matlab_fixed, matlab_like["fixed"].numpy(), align_scale=True),
        },
        "cor_after_permutation": {
            "current_before_alignment": sourcewise_complex_metrics(matlab_cor_spec, current["estCorSpec"].numpy(), align_scale=False),
            "current_after_complex_scale": sourcewise_complex_metrics(matlab_cor_spec, current["estCorSpec"].numpy(), align_scale=True),
            "matlab_like_before_alignment": sourcewise_complex_metrics(matlab_cor_spec, matlab_like["estCorSpec"].numpy(), align_scale=False),
            "matlab_like_after_complex_scale": sourcewise_complex_metrics(matlab_cor_spec, matlab_like["estCorSpec"].numpy(), align_scale=True),
        },
        "ips_after_permutation": {
            "current_before_alignment": sourcewise_complex_metrics(matlab_ips_spec, current["estIpsSpec"].numpy(), align_scale=False),
            "current_after_complex_scale": sourcewise_complex_metrics(matlab_ips_spec, current["estIpsSpec"].numpy(), align_scale=True),
            "matlab_like_before_alignment": sourcewise_complex_metrics(matlab_ips_spec, matlab_like["estIpsSpec"].numpy(), align_scale=False),
            "matlab_like_after_complex_scale": sourcewise_complex_metrics(matlab_ips_spec, matlab_like["estIpsSpec"].numpy(), align_scale=True),
        },
        "waveform_none": {
            "current_before_alignment": waveform_metrics(matlab_est_none, current["estNone"].numpy()),
            "current_after_scale_alignment": waveform_metrics(matlab_est_none, current["estNone"].numpy(), align_scale=True),
            "current_after_scale_delay_alignment": waveform_metrics(
                matlab_est_none, current["estNone"].numpy(), align_scale=True, align_delay=True
            ),
            "matlab_like_before_alignment": waveform_metrics(matlab_est_none, matlab_like["estNone"].numpy()),
            "matlab_like_after_scale_alignment": waveform_metrics(
                matlab_est_none, matlab_like["estNone"].numpy(), align_scale=True
            ),
        },
        "waveform_cor": {
            "current_before_alignment": waveform_metrics(matlab_est_cor, current["estCor"].numpy()),
            "current_after_scale_alignment": waveform_metrics(matlab_est_cor, current["estCor"].numpy(), align_scale=True),
            "current_after_scale_delay_alignment": waveform_metrics(
                matlab_est_cor, current["estCor"].numpy(), align_scale=True, align_delay=True
            ),
            "matlab_like_before_alignment": waveform_metrics(matlab_est_cor, matlab_like["estCor"].numpy()),
            "matlab_like_after_scale_alignment": waveform_metrics(
                matlab_est_cor, matlab_like["estCor"].numpy(), align_scale=True
            ),
        },
        "waveform_ips": {
            "current_before_alignment": waveform_metrics(matlab_est_ips, current["estIps"].numpy()),
            "current_after_scale_alignment": waveform_metrics(matlab_est_ips, current["estIps"].numpy(), align_scale=True),
            "current_after_scale_delay_alignment": waveform_metrics(
                matlab_est_ips, current["estIps"].numpy(), align_scale=True, align_delay=True
            ),
            "matlab_like_before_alignment": waveform_metrics(matlab_est_ips, matlab_like["estIps"].numpy()),
            "matlab_like_after_scale_alignment": waveform_metrics(
                matlab_est_ips, matlab_like["estIps"].numpy(), align_scale=True
            ),
        },
    }

    error_growth = {
        "stft_rmse": internal_stats["stft_obsSpec"]["rmse"],
        "iter1_current_y_rmse": complex_metrics(np.asarray(debug["Y_after_iter1"]), current["Y_after_iter1"].numpy())["rmse"],
        "iter1_matlab_like_y_rmse": complex_metrics(np.asarray(debug["Y_after_iter1"]), matlab_like["Y_after_iter1"].numpy())["rmse"],
        "iter10_current_y_rmse": complex_metrics(np.asarray(debug["Y_after_iter10"]), current["Y_after_iter10"].numpy())["rmse"],
        "iter10_matlab_like_y_rmse": complex_metrics(np.asarray(debug["Y_after_iter10"]), matlab_like["Y_after_iter10"].numpy())["rmse"],
        "iter25_current_y_rmse": complex_metrics(np.asarray(debug["Y_after_iter25"]), current["Y_after_iter25"].numpy())["rmse"],
        "iter25_matlab_like_y_rmse": complex_metrics(np.asarray(debug["Y_after_iter25"]), matlab_like["Y_after_iter25"].numpy())["rmse"],
        "iter50_current_y_rmse": complex_metrics(matlab_y, current["Y"].numpy())["rmse"],
        "iter50_matlab_like_y_rmse": complex_metrics(matlab_y, matlab_like["Y"].numpy())["rmse"],
        "projection_back_current_rmse": complex_metrics(matlab_fixed, current["fixed"].numpy())["rmse"],
        "projection_back_matlab_like_rmse": complex_metrics(matlab_fixed, matlab_like["fixed"].numpy())["rmse"],
        "cor_after_perm_current_rmse": complex_metrics(matlab_cor_spec, current["estCorSpec"].numpy())["rmse"],
        "cor_after_perm_matlab_like_rmse": complex_metrics(matlab_cor_spec, matlab_like["estCorSpec"].numpy())["rmse"],
        "ips_after_perm_current_rmse": complex_metrics(matlab_ips_spec, current["estIpsSpec"].numpy())["rmse"],
        "ips_after_perm_matlab_like_rmse": complex_metrics(matlab_ips_spec, matlab_like["estIpsSpec"].numpy())["rmse"],
        "cor_waveform_current_rmse_mean": float(np.mean([m["rmse"] for m in stage_metrics["waveform_cor"]["current_before_alignment"]])),
        "ips_waveform_current_rmse_mean": float(np.mean([m["rmse"] for m in stage_metrics["waveform_ips"]["current_before_alignment"]])),
    }

    summary = {
        "conditions_source": {
            "matlab": "main.m parameter block and getInputFileNames(dataNo=1)",
            "python": "python_fdica/main.py dataset 1 path and bssAuxFdica call",
            "note": "isDraw=true is a plotting/debug flag in normal entry points; this tool records cost without opening GUI figures.",
        },
        "params": params,
        "matlab_metrics": extract_metric_table(mat["metrics"]),
        "python_metrics": extract_metric_table(eval_mat["metrics"]),
        "cost_stats_current": cost_stats(matlab_cost, current_cost),
        "cost_stats_matlab_like": cost_stats(matlab_cost, matlab_like_cost),
        "first_noticeable_cost_error_current": first_noticeable_cost_error(matlab_cost, current_cost),
        "first_noticeable_cost_error_matlab_like": first_noticeable_cost_error(matlab_cost, matlab_like_cost),
        "internal_stats": internal_stats,
        "stage_metrics": stage_metrics,
        "error_growth": error_growth,
        "permutation": {
            "cor_current": permutation_summary(np.asarray(mat["permCor"]), current["permCor"].numpy()),
            "cor_matlab_like": permutation_summary(np.asarray(mat["permCor"]), matlab_like["permCor"].numpy()),
            "ips_current": permutation_summary(np.asarray(mat["permIps"]), current["permIps"].numpy()),
            "ips_matlab_like": permutation_summary(np.asarray(mat["permIps"]), matlab_like["permIps"].numpy()),
        },
        "cor_mismatch_energy": {
            "total_energy": cor_energy["total_energy"],
            "mismatch_energy": cor_energy["mismatch_energy"],
            "mismatch_energy_ratio": cor_energy["mismatch_energy_ratio"],
        },
        "notes": [
            "Current Python local_auxFdica uses a tiny ridge regularization and pinv fallback; MATLAB production code does not.",
            "MATLAB-like Python path removes that regularization only inside this comparison tool.",
        ],
    }
    with (OUT_DIR / "existing_condition_summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
