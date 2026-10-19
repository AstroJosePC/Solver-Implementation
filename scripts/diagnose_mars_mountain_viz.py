#!/usr/bin/env python3
"""
diagnose_mars_mountain_viz.py

Diagnose whether suspicious mountain-sweep outputs are primarily visualization/
post-processing artifacts or likely physics/simulation issues.

The script does not rerun simulations. It reuses existing netCDF snapshots and
compares outcomes under controlled plotting and signal-preprocessing variants.
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import re
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset

from wavelet_tc98 import wavelet


RUN_DIR_PATTERN = re.compile(r"run_U(?P<u>[-+0-9.]+)_sig(?P<sig>[-+0-9.]+)$")


@dataclass
class RunInfo:
    run_name: str
    run_dir: str
    u_bar: float
    sigma_km: float


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Diagnose viz-vs-physics causes in Mars mountain sweep outputs")
    p.add_argument("--sweep-dir", default="output/sweep_mars_mountain", help="Directory containing run_U*_sig* folders")
    p.add_argument("--outdir", default="output/mars_mountain_viz_diagnostics", help="Directory for diagnostics outputs")
    p.add_argument("--run-pattern", default="run_U*_sig*", help="Glob pattern to select runs")
    p.add_argument(
        "--run-list",
        default="",
        help="Optional comma-separated exact run folder names (overrides --run-pattern)",
    )
    p.add_argument("--station-x-rel-km", default="20,40,60,80", help="Comma-separated station x positions relative to mountain center")
    p.add_argument("--station-z-km", type=float, default=20.0, help="Station z position in km")
    p.add_argument("--wavelet-station-index", type=int, default=2, help="Zero-based station index for wavelet diagnostics")
    p.add_argument("--snapshot-z-max-km", type=float, default=20.0, help="Max z (km) for snapshot diagnostics")
    p.add_argument("--snapshot-x-limit-km", type=float, default=150.0, help="Half-width x limit (km) for snapshot diagnostics")
    p.add_argument("--fixed-vmax-w", type=float, default=0.2, help="Fixed symmetric vmax for snapshot comparison")
    p.add_argument("--trim-minutes", type=float, default=30.0, help="Transient trim window for wavelet/time-series diagnostics")
    p.add_argument("--highpass-window-minutes", type=float, default=60.0, help="Moving-average window used in high-pass preprocessing")
    p.add_argument("--wavelet-max-freq-mhz", type=float, default=10.0, help="Max frequency shown in wavelet comparison")
    return p.parse_args()


def parse_station_list(text: str) -> list[float]:
    vals = [x.strip() for x in text.split(",") if x.strip()]
    if not vals:
        raise ValueError("No station offsets provided")
    return [float(v) for v in vals]


def parse_run_info(run_dir: str) -> RunInfo | None:
    run_name = os.path.basename(os.path.normpath(run_dir))
    m = RUN_DIR_PATTERN.match(run_name)
    if not m:
        return None
    return RunInfo(run_name=run_name, run_dir=run_dir, u_bar=float(m.group("u")), sigma_km=float(m.group("sig")))


def discover_runs(sweep_dir: str, run_pattern: str) -> list[RunInfo]:
    dirs = sorted(glob.glob(os.path.join(sweep_dir, run_pattern)))
    out: list[RunInfo] = []
    for d in dirs:
        info = parse_run_info(d)
        if info is not None:
            out.append(info)
    out.sort(key=lambda r: (r.u_bar, r.sigma_km, r.run_name))
    return out


def filter_runs_by_name(runs: list[RunInfo], run_list_text: str) -> list[RunInfo]:
    names = [x.strip() for x in run_list_text.split(",") if x.strip()]
    if not names:
        return runs
    keep = set(names)
    return [r for r in runs if r.run_name in keep]


def get_sorted_files_and_times(run_dir: str) -> tuple[list[str], np.ndarray]:
    files = sorted(glob.glob(os.path.join(run_dir, "*.nc")))
    if not files:
        return [], np.array([], dtype=np.float64)
    times = np.empty(len(files), dtype=np.float64)
    for i, f in enumerate(files):
        with Dataset(f, "r") as ds:
            times[i] = float(ds.getncattr("time_s"))
    order = np.argsort(times)
    return [files[i] for i in order], times[order]


def centered_x_km(x_m: np.ndarray) -> np.ndarray:
    center = 0.5 * (float(np.min(x_m)) + float(np.max(x_m)))
    return (x_m - center) / 1000.0


def detrend_linear(y: np.ndarray) -> np.ndarray:
    x = np.arange(len(y), dtype=np.float64)
    c = np.polyfit(x, y, 1)
    trend = c[0] * x + c[1]
    return y - trend


def highpass_moving_average(y: np.ndarray, window_samples: int) -> np.ndarray:
    if window_samples <= 1:
        return y - np.mean(y)
    window_samples = min(window_samples, len(y))
    kernel = np.ones(window_samples, dtype=np.float64) / float(window_samples)
    smooth = np.convolve(y, kernel, mode="same")
    return y - smooth


def compute_wavelet_power(signal: np.ndarray, dt_s: float) -> tuple[np.ndarray, np.ndarray]:
    dj = 0.05
    s0 = 2.0 * dt_s
    j1 = int(np.fix(np.log2(len(signal) * dt_s / s0) / dj))
    wv, period_s, _, _ = wavelet(signal - np.mean(signal), dt_s, pad=1, dj=dj, s0=s0, J1=j1, mother="MORLET", param=6.0)
    power = np.abs(wv) ** 2
    return power, period_s


def normalize_power_max(power: np.ndarray) -> np.ndarray:
    pmax = float(np.max(power))
    if pmax <= 0.0:
        return power
    return power / pmax


def normalize_power_percentile(power: np.ndarray, p: float = 99.5) -> np.ndarray:
    den = float(np.percentile(power, p))
    if den <= 0.0:
        return power
    out = power / den
    return np.clip(out, 0.0, 1.0)


def dominant_frequency_mhz(power: np.ndarray, period_s: np.ndarray, min_freq_mhz: float = 0.2) -> float:
    freq_mhz = 1000.0 / period_s
    scale_mask = freq_mhz >= min_freq_mhz
    if not np.any(scale_mask):
        return float("nan")
    sub = power[scale_mask, :]
    i_scale, _ = np.unravel_index(np.argmax(sub), sub.shape)
    return float(freq_mhz[scale_mask][i_scale])


def align_power_to_period(
    source_power: np.ndarray,
    source_period_s: np.ndarray,
    target_period_s: np.ndarray,
) -> np.ndarray:
    """Interpolate power along frequency axis to match a target period grid."""
    src_freq = 1000.0 / source_period_s
    tgt_freq = 1000.0 / target_period_s

    # Interpolation requires ascending x; wavelet freq is typically descending.
    src_order = np.argsort(src_freq)
    src_freq_sorted = src_freq[src_order]
    src_power_sorted = source_power[src_order, :]

    out = np.empty((len(tgt_freq), source_power.shape[1]), dtype=np.float64)
    for j in range(source_power.shape[1]):
        out[:, j] = np.interp(
            tgt_freq,
            src_freq_sorted,
            src_power_sorted[:, j],
            left=src_power_sorted[0, j],
            right=src_power_sorted[-1, j],
        )
    return out


def plot_snapshot_comparison(
    outpath: str,
    x_km: np.ndarray,
    z_km: np.ndarray,
    w2d: np.ndarray,
    z_max_km: float,
    x_limit_km: float,
    fixed_vmax: float,
) -> dict[str, float]:
    iz = np.where(z_km <= z_max_km)[0]
    w = w2d[iz, :]
    q995 = float(np.quantile(np.abs(w), 0.995))
    q98 = float(np.quantile(np.abs(w), 0.98))
    maxabs = float(np.max(np.abs(w)))

    vmax_q = max(q995, 1.0e-8)
    vmax_fixed = max(fixed_vmax, 1.0e-8)
    vmax_r = max(q98, 1.0e-8)

    fig, axes = plt.subplots(1, 3, figsize=(17, 4.8), constrained_layout=True)

    configs = [
        ("quantile 99.5%", vmax_q),
        (f"fixed +/-{vmax_fixed:.3g}", vmax_fixed),
        ("robust 98%", vmax_r),
    ]

    for ax, (title, vmax) in zip(axes, configs):
        levels = np.linspace(-vmax, vmax, 51)
        im = ax.contourf(x_km, z_km[iz], w, levels=levels, cmap="RdBu_r", extend="both")
        ax.set_title(title)
        ax.set_xlabel("x rel. mountain center (km)")
        ax.set_ylabel("z (km)")
        ax.set_xlim(-x_limit_km, x_limit_km)
        fig.colorbar(im, ax=ax)

    plt.savefig(outpath, dpi=180, bbox_inches="tight")
    plt.close(fig)

    return {
        "w_max_abs": maxabs,
        "w_q98_abs": q98,
        "w_q995_abs": q995,
        "outlier_ratio_max_to_q995": (maxabs / q995) if q995 > 0 else float("inf"),
    }


def plot_station_preproc_comparison(
    outpath: str,
    times_s: np.ndarray,
    raw: np.ndarray,
    trimmed: np.ndarray,
    trimmed_times_s: np.ndarray,
    detrended: np.ndarray,
    highpass: np.ndarray,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharex=False)
    axes = axes.ravel()

    axes[0].plot(times_s, raw, lw=1.0)
    axes[0].set_title("raw station w(t)")
    axes[0].set_xlabel("time (s)")
    axes[0].set_ylabel("w (m/s)")
    axes[0].grid(True, alpha=0.25)

    axes[1].plot(trimmed_times_s, trimmed, lw=1.0)
    axes[1].set_title("trimmed transient")
    axes[1].set_xlabel("time (s)")
    axes[1].set_ylabel("w (m/s)")
    axes[1].grid(True, alpha=0.25)

    axes[2].plot(trimmed_times_s, detrended, lw=1.0)
    axes[2].set_title("trimmed + linear detrend")
    axes[2].set_xlabel("time (s)")
    axes[2].set_ylabel("w (m/s)")
    axes[2].grid(True, alpha=0.25)

    axes[3].plot(trimmed_times_s, highpass, lw=1.0)
    axes[3].set_title("trimmed + moving-average high-pass")
    axes[3].set_xlabel("time (s)")
    axes[3].set_ylabel("w (m/s)")
    axes[3].grid(True, alpha=0.25)

    plt.savefig(outpath, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_wavelet_comparison(
    outpath: str,
    times_s: np.ndarray,
    period_s: np.ndarray,
    p_raw_max: np.ndarray,
    p_trim_max: np.ndarray,
    p_detr_max: np.ndarray,
    p_detr_pct: np.ndarray,
    max_freq_mhz: float,
) -> None:
    freq_mhz = 1000.0 / period_s
    t_min = times_s / 60.0

    fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharex=True, sharey=True)
    axes = axes.ravel()

    panels = [
        (p_raw_max, "raw + max norm"),
        (p_trim_max, "trimmed + max norm"),
        (p_detr_max, "trimmed+detrended + max norm"),
        (p_detr_pct, "trimmed+detrended + p99.5 norm"),
    ]

    for ax, (pw, ttl) in zip(axes, panels):
        im = ax.pcolormesh(t_min, freq_mhz, pw, cmap="jet", shading="auto", vmin=0.0, vmax=1.0)
        ax.axhline(1.8, color="limegreen", ls="--", lw=2)
        ax.set_ylim(0.0, max_freq_mhz)
        ax.set_title(ttl)
        ax.grid(False)

    axes[2].set_xlabel("time (min)")
    axes[3].set_xlabel("time (min)")
    axes[0].set_ylabel("freq (mHz)")
    axes[2].set_ylabel("freq (mHz)")

    cbar = fig.colorbar(im, ax=axes.tolist(), shrink=0.95)
    cbar.set_label("normalized power")
    plt.savefig(outpath, dpi=180, bbox_inches="tight")
    plt.close(fig)


def extract_station_series(
    files: list[str],
    times_s: np.ndarray,
    station_x_rel_km: list[float],
    station_z_km: float,
    station_idx: int,
) -> tuple[np.ndarray, np.ndarray, float, float, str]:
    with Dataset(files[0], "r") as ds0:
        x_m = np.asarray(ds0.variables["x"][:], dtype=np.float64)
        z_km = np.asarray(ds0.variables["z"][:], dtype=np.float64) / 1000.0

    x_rel = centered_x_km(x_m)
    iz = int(np.argmin(np.abs(z_km - station_z_km)))
    ix_list = [int(np.argmin(np.abs(x_rel - xk))) for xk in station_x_rel_km]

    sid = int(np.clip(station_idx, 0, len(ix_list) - 1))
    ix = ix_list[sid]
    label = f"S{sid+1}"

    y = np.empty(len(files), dtype=np.float64)
    for i, f in enumerate(files):
        with Dataset(f, "r") as ds:
            y[i] = float(ds.variables["w"][iz, ix])

    return np.asarray(times_s, dtype=np.float64), y, float(x_rel[ix]), float(z_km[iz]), label


def classify_verdict(
    transient_ratio: float,
    drift_ratio: float,
    f_raw: float,
    f_detr: float,
    outlier_ratio: float,
) -> tuple[str, list[str]]:
    reasons = []

    if transient_ratio >= 3.0:
        reasons.append("strong startup transient dominates early signal")
    if drift_ratio >= 0.8:
        reasons.append("large low-frequency drift relative to variability")
    if np.isfinite(f_raw) and np.isfinite(f_detr):
        if (f_raw < 0.2 and f_detr > 0.5) or (abs(f_raw - f_detr) > 0.8):
            reasons.append("dominant frequency unstable under detrending")
    if outlier_ratio >= 6.0:
        reasons.append("snapshot dynamic range dominated by outliers")

    if len(reasons) >= 2:
        return "Viz Confirmed", reasons
    if len(reasons) == 1:
        return "Mixed", reasons
    return "Physics Suspected", ["diagnostic controls did not materially change interpretation"]


def analyze_run(run: RunInfo, args: argparse.Namespace) -> dict[str, object]:
    run_out = os.path.join(args.outdir, run.run_name)
    os.makedirs(run_out, exist_ok=True)

    files, times_s = get_sorted_files_and_times(run.run_dir)
    if len(files) < 8:
        return {
            "run": run.run_name,
            "status": "failed",
            "reason": "insufficient snapshots",
        }

    # Snapshot diagnostics from latest state.
    with Dataset(files[-1], "r") as ds_last:
        x_m = np.asarray(ds_last.variables["x"][:], dtype=np.float64)
        z_m = np.asarray(ds_last.variables["z"][:], dtype=np.float64)
        w_latest = np.asarray(ds_last.variables["w"][:, :], dtype=np.float64)

    x_km = centered_x_km(x_m)
    z_km = z_m / 1000.0

    snap_stats = plot_snapshot_comparison(
        outpath=os.path.join(run_out, "snapshot_scaling_diagnostic.png"),
        x_km=x_km,
        z_km=z_km,
        w2d=w_latest,
        z_max_km=args.snapshot_z_max_km,
        x_limit_km=args.snapshot_x_limit_km,
        fixed_vmax=args.fixed_vmax_w,
    )

    # Station series extraction.
    st_targets = parse_station_list(args.station_x_rel_km)
    t, raw, station_x_km, station_z_actual_km, station_label = extract_station_series(
        files=files,
        times_s=times_s,
        station_x_rel_km=st_targets,
        station_z_km=args.station_z_km,
        station_idx=args.wavelet_station_index,
    )

    dt_s = float(np.median(np.diff(t)))
    trim_s = args.trim_minutes * 60.0
    keep = t >= (t[0] + trim_s)
    if np.count_nonzero(keep) < 16:
        keep = np.arange(len(t)) >= int(0.1 * len(t))

    t_trim = t[keep]
    y_trim = raw[keep]
    y_detr = detrend_linear(y_trim)
    hp_window_samples = max(3, int(round((args.highpass_window_minutes * 60.0) / max(dt_s, 1.0e-9))))
    y_hp = highpass_moving_average(y_trim, hp_window_samples)

    plot_station_preproc_comparison(
        outpath=os.path.join(run_out, "station_preprocessing_diagnostic.png"),
        times_s=t,
        raw=raw,
        trimmed=y_trim,
        trimmed_times_s=t_trim,
        detrended=y_detr,
        highpass=y_hp,
    )

    # Wavelet diagnostics (same data, different preprocessing/normalization).
    p_raw, per_raw = compute_wavelet_power(raw, dt_s=dt_s)
    p_trim, per_trim = compute_wavelet_power(y_trim, dt_s=dt_s)
    p_detr, per_detr = compute_wavelet_power(y_detr, dt_s=dt_s)

    # Interpolate panel frequencies by using detrended period as a common axis
    # only when dimensions differ; otherwise they are typically equal.
    period_s = per_detr
    p_raw_aligned = align_power_to_period(p_raw, per_raw, period_s)
    p_trim_aligned = align_power_to_period(p_trim, per_trim, period_s)
    p_detr_aligned = p_detr

    p_raw_n = normalize_power_max(p_raw_aligned)
    p_trim_n = normalize_power_max(p_trim_aligned)
    p_detr_n = normalize_power_max(p_detr_aligned)
    p_detr_p = normalize_power_percentile(p_detr_aligned, p=99.5)

    # Align matrix widths for plotting with trimmed time axis.
    n_t = len(t_trim)
    if p_raw_n.shape[1] != n_t:
        p_raw_n = p_raw_n[:, -n_t:]
    if p_trim_n.shape[1] != n_t:
        p_trim_n = p_trim_n[:, -n_t:]
    if p_detr_n.shape[1] != n_t:
        p_detr_n = p_detr_n[:, -n_t:]
    if p_detr_p.shape[1] != n_t:
        p_detr_p = p_detr_p[:, -n_t:]

    plot_wavelet_comparison(
        outpath=os.path.join(run_out, "wavelet_diagnostic_comparison.png"),
        times_s=t_trim,
        period_s=period_s,
        p_raw_max=p_raw_n,
        p_trim_max=p_trim_n,
        p_detr_max=p_detr_n,
        p_detr_pct=p_detr_p,
        max_freq_mhz=args.wavelet_max_freq_mhz,
    )

    # Quantitative checks for verdict.
    split_idx = max(3, int(0.1 * len(raw)))
    rms_early = float(np.sqrt(np.mean(raw[:split_idx] ** 2)))
    rms_late = float(np.sqrt(np.mean(raw[split_idx:] ** 2)))
    transient_ratio = (rms_early / max(rms_late, 1.0e-12))

    drift_ratio = float(abs(np.mean(y_trim)) / max(np.std(y_trim), 1.0e-12))
    f_raw = dominant_frequency_mhz(normalize_power_max(p_raw), per_raw, min_freq_mhz=0.01)
    f_detr = dominant_frequency_mhz(normalize_power_max(p_detr), per_detr, min_freq_mhz=0.2)

    verdict, reasons = classify_verdict(
        transient_ratio=transient_ratio,
        drift_ratio=drift_ratio,
        f_raw=f_raw,
        f_detr=f_detr,
        outlier_ratio=float(snap_stats["outlier_ratio_max_to_q995"]),
    )

    report = {
        "run": run.run_name,
        "status": "ok",
        "U_bar_mps": run.u_bar,
        "sigma_km": run.sigma_km,
        "station": station_label,
        "station_x_actual_rel_km": station_x_km,
        "station_z_actual_km": station_z_actual_km,
        "dt_s": dt_s,
        "trim_minutes": args.trim_minutes,
        "transient_ratio": transient_ratio,
        "drift_ratio": drift_ratio,
        "dominant_freq_raw_mhz": f_raw,
        "dominant_freq_detrended_mhz": f_detr,
        "snapshot_stats": snap_stats,
        "verdict": verdict,
        "reasons": reasons,
    }

    with open(os.path.join(run_out, "diagnostic_report.json"), "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)

    return report


def write_summary_csv(path: str, rows: list[dict[str, object]]) -> None:
    fields = [
        "run",
        "status",
        "U_bar_mps",
        "sigma_km",
        "station",
        "dt_s",
        "trim_minutes",
        "transient_ratio",
        "drift_ratio",
        "dominant_freq_raw_mhz",
        "dominant_freq_detrended_mhz",
        "verdict",
        "reasons",
    ]
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            row_out = {k: row.get(k, "") for k in fields}
            if isinstance(row_out.get("reasons"), list):
                row_out["reasons"] = " | ".join(row_out["reasons"])
            writer.writerow(row_out)


def main() -> None:
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    runs = discover_runs(args.sweep_dir, args.run_pattern)
    runs = filter_runs_by_name(runs, args.run_list)
    if not runs:
        raise ValueError(f"No runs found in {args.sweep_dir} with pattern {args.run_pattern}")

    print(f"Found {len(runs)} runs")
    reports: list[dict[str, object]] = []
    for i, run in enumerate(runs, start=1):
        print(f"[{i}/{len(runs)}] Diagnosing {run.run_name} ...")
        try:
            rep = analyze_run(run, args)
            reports.append(rep)
            print(f"  verdict: {rep.get('verdict', rep.get('status'))}")
        except Exception as exc:
            rep = {
                "run": run.run_name,
                "status": "failed",
                "verdict": "failed",
                "reasons": [str(exc)],
            }
            reports.append(rep)
            print(f"  failed: {exc}")

    summary_csv = os.path.join(args.outdir, "diagnostic_summary.csv")
    write_summary_csv(summary_csv, reports)
    with open(os.path.join(args.outdir, "diagnostic_summary.json"), "w", encoding="utf-8") as f:
        json.dump(reports, f, indent=2)

    n_ok = sum(1 for r in reports if r.get("status") == "ok")
    print(f"Completed diagnostics. ok={n_ok}, failed={len(reports)-n_ok}")
    print(f"Summary: {summary_csv}")


if __name__ == "__main__":
    main()
