#!/usr/bin/env python3
"""
analyze_mars_mountain_sweep.py

Post-process Mars mountain sweep runs produced by scripts/sweep_mars_mountain.jl.
Per run, this script generates:
1) Latest-time snapshot plots (w and density-scaled w)
2) Station time series for w at mountain-centered stations
3) Wavelet spectrogram (Morlet) for a selected station
4) Per-run metrics JSON
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import re
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset

from wavelet_tc98 import wavelet


# Matches both the original sigma-sweep dirs (run_U<u>_sig<sig>) and the
# Arsia Stage 1 h_max-sweep dirs (run_U<u>_h<h>). Exactly one of the sig/h
# groups is populated per directory name.
RUN_DIR_PATTERN = re.compile(
    r"run_U(?P<u>[-+0-9.]+)_(?:sig(?P<sig>[-+0-9.]+)|h(?P<h>[-+0-9.]+))$"
)

# Arsia Stage 1 holds the Gaussian hill half-width fixed; h_max is the swept
# parameter, so sigma is not encoded in those run names.
SIGMA_FIXED_KM = 100.0


@dataclass
class RunInfo:
    run_name: str
    run_dir: str
    u_bar: float
    sigma_km: float
    h_max_km: float | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Analyze Mars mountain sweep outputs")
    parser.add_argument("--sweep-dir", default="output/sweep_mars_mountain", help="Directory containing run_U*_sig* folders")
    parser.add_argument("--outdir", default="output/mars_mountain_sweep_analysis", help="Directory for analysis outputs")
    parser.add_argument(
        "--run-pattern",
        default="run_U*_sig*",
        help="Glob pattern used to select run folders inside --sweep-dir",
    )
    parser.add_argument(
        "--station-x-rel-km",
        default="20,40,60,80",
        help="Comma-separated station x locations in km relative to mountain center",
    )
    parser.add_argument("--station-z-km", type=float, default=20.0, help="Station z location in km")
    parser.add_argument(
        "--wavelet-station-index",
        type=int,
        default=2,
        help="Zero-based station index used for wavelet (default 2 -> third station)",
    )
    parser.add_argument("--snapshot-z-max-km", type=float, default=20.0, help="Max z for snapshot plots in km")
    parser.add_argument("--snapshot-x-limit-km", type=float, default=150.0, help="Half-width x-limit for snapshot plots in km")  
    parser.add_argument("--wavelet-trim-min", type=float, default=120.0, help="Minutes of transient to trim from wavelet analysis")
    parser.add_argument("--wavelet-dj", type=float, default=0.05, help="Scale spacing for wavelet")
    parser.add_argument("--wavelet-morlet-param", type=float, default=6.0, help="Morlet nondimensional parameter")
    parser.add_argument("--wavelet-min-freq-mhz", type=float, default=0.0, help="Min frequency for spectrogram y-axis")
    parser.add_argument("--wavelet-max-freq-mhz", type=float, default=10.0, help="Max frequency for spectrogram y-axis")
    # --- Sweep-level composite (opt-in) ---
    parser.add_argument("--composite", action="store_true", help="After the per-run loop, build a sweep-level composite figure and band response table")
    parser.add_argument("--composite-only", action="store_true", help="Skip per-run analysis; only (re)build the sweep composite and band response (implies --composite)")
    parser.add_argument("--band-label", default="AMEC", help="Label for the response band in outputs")
    parser.add_argument("--band-z-min-km", type=float, default=45.0, help="Lower edge of the response band (km)")
    parser.add_argument("--band-z-max-km", type=float, default=53.0, help="Upper edge of the response band (km)")
    parser.add_argument("--composite-z-max-km", type=float, default=90.0, help="Max z for composite panels (km)")
    parser.add_argument("--composite-x-limit-km", type=float, default=500.0, help="Half-width x-limit for composite panels (km)")
    parser.add_argument("--sponge-base-km", type=float, default=60.0, help="Sponge base altitude to mark on composite panels (km)")
    return parser.parse_args()


def parse_station_x_list(text: str) -> list[float]:
    vals = [v.strip() for v in text.split(",") if v.strip()]
    if not vals:
        raise ValueError("No station x values were provided")
    return [float(v) for v in vals]


def parse_run_info(run_dir: str) -> RunInfo | None:
    run_name = os.path.basename(os.path.normpath(run_dir))
    m = RUN_DIR_PATTERN.match(run_name)
    if m is None:
        return None
    sig = m.group("sig")
    h = m.group("h")
    return RunInfo(
        run_name=run_name,
        run_dir=run_dir,
        u_bar=float(m.group("u")),
        sigma_km=float(sig) if sig is not None else SIGMA_FIXED_KM,
        h_max_km=float(h) if h is not None else None,
    )


def get_snapshot_files_and_times(run_dir: str) -> tuple[list[str], np.ndarray]:
    files = sorted(glob.glob(os.path.join(run_dir, "*.nc")))
    if not files:
        return [], np.array([], dtype=np.float64)

    times = np.empty(len(files), dtype=np.float64)
    for i, fpath in enumerate(files):
        with Dataset(fpath, "r") as ds:
            times[i] = float(ds.getncattr("time_s"))

    order = np.argsort(times)
    sorted_files = [files[i] for i in order]
    sorted_times = times[order]
    return sorted_files, sorted_times


def build_centered_x_km(x_m: np.ndarray) -> np.ndarray:
    x_center = 0.5 * (float(np.min(x_m)) + float(np.max(x_m)))
    return (x_m - x_center) / 1000.0


def write_failed_metrics(outdir: str, run_info: RunInfo, reason: str) -> None:
    payload = {
        "status": "failed",
        "run_name": run_info.run_name,
        "U_bar_mps": run_info.u_bar,
        "sigma_km": run_info.sigma_km,
        "h_max_km": run_info.h_max_km,
        "reason": reason,
    }
    path = os.path.join(outdir, "metrics.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)


def plot_latest_snapshot(
    run_output_dir: str,
    run_info: RunInfo,
    x_m: np.ndarray,
    z_m: np.ndarray,
    w: np.ndarray,
    rho: np.ndarray,
    time_s: float,
    z_max_km: float,
    x_limit_km: float,
) -> None:
    xx_km = build_centered_x_km(x_m)
    z_km = z_m / 1000.0
    iz = np.where(z_km <= z_max_km)[0]
    if len(iz) == 0:
        raise ValueError(f"No z-levels <= {z_max_km} km found for snapshot")

    rho_surface = float(np.mean(rho[0, :]))
    rho_ratio = rho / rho_surface
    w_scaled = w * np.sqrt(rho_ratio)

    w_plot = w[iz, :]
    ws_plot = w_scaled[iz, :]

    vmax_w = float(np.quantile(np.abs(w_plot), 0.995))
    vmax_ws = float(np.quantile(np.abs(ws_plot), 0.995))
    vmax_w = max(vmax_w, 1.0e-8)
    vmax_ws = max(vmax_ws, 1.0e-8)

    levels_w = np.linspace(-vmax_w, vmax_w, 51)
    levels_ws = np.linspace(-vmax_ws, vmax_ws, 51)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)

    ax = axes[0]
    im0 = ax.contourf(xx_km, z_km[iz], w_plot, levels=levels_w, cmap="RdBu_r", extend="both")
    ax.set_title("w simulated (m/s)")
    ax.set_xlabel("x relative to mountain center (km)")
    ax.set_ylabel("z (km)")
    ax.set_xlim(-x_limit_km, x_limit_km)
    fig.colorbar(im0, ax=ax)

    ax = axes[1]
    im1 = ax.contourf(xx_km, z_km[iz], ws_plot, levels=levels_ws, cmap="RdBu_r", extend="both")
    ax.set_title(r"$w\sqrt{\rho/\rho_0}$ (m/s)")
    ax.set_xlabel("x relative to mountain center (km)")
    ax.set_ylabel("z (km)")
    ax.set_xlim(-x_limit_km, x_limit_km)
    fig.colorbar(im1, ax=ax)

    fig.suptitle(
        f"{run_info.run_name} | U={run_info.u_bar:.1f} m/s, sigma={run_info.sigma_km:.1f} km | t={time_s:.1f} s"
    )

    outpath = os.path.join(run_output_dir, "snapshot_latest.png")
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)


def extract_station_series(
    files: list[str],
    times: np.ndarray,
    station_x_rel_km: list[float],
    station_z_km: float,
) -> tuple[np.ndarray, list[float], float, dict[str, np.ndarray]]:
    if not files:
        raise ValueError("No snapshot files provided")

    with Dataset(files[0], "r") as ds0:
        x_m = np.asarray(ds0.variables["x"][:], dtype=np.float64)
        z_km = np.asarray(ds0.variables["z"][:], dtype=np.float64) / 1000.0

    x_rel_km = build_centered_x_km(x_m)
    iz = int(np.argmin(np.abs(z_km - station_z_km)))

    ix_list = [int(np.argmin(np.abs(x_rel_km - sx))) for sx in station_x_rel_km]
    station_x_actual_km = [float(x_rel_km[ix]) for ix in ix_list]
    station_z_actual_km = float(z_km[iz])

    n_files = len(files)
    station_arrays = {f"S{i+1}": np.empty(n_files, dtype=np.float64) for i in range(len(ix_list))}

    for it, fpath in enumerate(files):
        with Dataset(fpath, "r") as ds:
            w_row = np.asarray(ds.variables["w"][iz, :], dtype=np.float64)
        for i, ix in enumerate(ix_list):
            station_arrays[f"S{i+1}"][it] = float(w_row[ix])

    return np.asarray(times, dtype=np.float64), station_x_actual_km, station_z_actual_km, station_arrays


def plot_station_timeseries(
    run_output_dir: str,
    run_info: RunInfo,
    times_s: np.ndarray,
    station_x_actual_km: list[float],
    station_z_actual_km: float,
    station_arrays: dict[str, np.ndarray],
) -> None:
    n_stations = len(station_arrays)
    fig, axes = plt.subplots(n_stations, 1, figsize=(10, 2.2 * n_stations), sharex=True)
    if n_stations == 1:
        axes = [axes]

    for i, (label, arr) in enumerate(station_arrays.items()):
        ax = axes[i]
        ax.plot(times_s, arr, color="tab:blue", linewidth=1.2)
        ax.set_ylabel("w (m/s)")
        ax.grid(True, alpha=0.25)
        ax.text(
            0.99,
            0.88,
            f"{label}: x={station_x_actual_km[i]:.1f} km, z={station_z_actual_km:.1f} km",
            transform=ax.transAxes,
            ha="right",
            fontsize=9,
        )

    axes[-1].set_xlabel("time (s)")
    fig.suptitle(f"{run_info.run_name} station time series")
    outpath = os.path.join(run_output_dir, "station_timeseries_w.png")
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)


def compute_wavelet(
    signal: np.ndarray,
    dt_s: float,
    dj: float,
    morlet_param: float,
) -> tuple[np.ndarray, np.ndarray]:
    s0 = 2.0 * dt_s
    j1 = int(np.fix(np.log2(len(signal) * dt_s / s0) / dj))
    wave, period, _, _ = wavelet(
        signal - np.mean(signal),
        dt_s,
        pad=1,
        dj=dj,
        s0=s0,
        J1=j1,
        mother="MORLET",
        param=morlet_param,
    )
    power = np.abs(wave) ** 2
    pmax = float(np.max(power))
    if pmax > 0.0:
        power = power / pmax
    return power, period


def plot_wavelet(
    run_output_dir: str,
    run_info: RunInfo,
    times_s: np.ndarray,
    power: np.ndarray,
    period_s: np.ndarray,
    station_label: str,
    min_freq_mhz: float,
    max_freq_mhz: float,
) -> tuple[float, float]:
    freq_mhz = 1000.0 / period_s
    t_min = times_s / 60.0

    # Dominant peak across all time-frequency bins.
    i_scale, _ = np.unravel_index(np.argmax(power), power.shape)
    dominant_freq_mhz = float(freq_mhz[i_scale])
    dominant_period_s = float(period_s[i_scale])

    fig, ax = plt.subplots(figsize=(9, 5.5))
    im = ax.pcolormesh(t_min, freq_mhz, power, cmap="jet", shading="auto")
    ax.axhline(1.8, color="limegreen", linestyle="--", linewidth=2, label="acoustic cutoff ~1.8 mHz")
    ax.set_xlabel("time (min)")
    ax.set_ylabel("frequency (mHz)")
    ax.set_ylim(min_freq_mhz, max_freq_mhz)
    ax.set_xlim(float(t_min[0]), float(t_min[-1]))
    ax.set_title(f"{run_info.run_name} wavelet power ({station_label})")
    ax.legend(loc="upper right", fontsize=9)
    fig.colorbar(im, ax=ax, label="normalized power")

    outpath = os.path.join(run_output_dir, "wavelet_w.png")
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return dominant_freq_mhz, dominant_period_s


def plot_mean_flow_evolution(
    run_output_dir: str,
    run_info: RunInfo,
    files: list[str],
    times_s: np.ndarray,
) -> None:
    if not files:
        return

    with Dataset(files[0], "r") as ds0:
        z_km = np.asarray(ds0.variables["z"][:], dtype=np.float64) / 1000.0

    n_times = len(files)
    n_z = len(z_km)
    u_mean = np.empty((n_times, n_z), dtype=np.float64)

    for i, fpath in enumerate(files):
        with Dataset(fpath, "r") as ds:
            u = np.asarray(ds.variables["u"][:, :], dtype=np.float64)
            u_mean[i, :] = np.mean(u, axis=1)

    u_dev = u_mean - run_info.u_bar
    t_hr = times_s / 3600.0

    fig, ax = plt.subplots(figsize=(10, 5))
    vmax = float(np.quantile(np.abs(u_dev), 0.995))
    vmax = max(vmax, 1.0e-3)
    
    im = ax.pcolormesh(t_hr, z_km, u_dev.T, cmap="RdBu_r", shading="auto", vmin=-vmax, vmax=vmax)
    ax.set_xlabel("time (hours)")
    ax.set_ylabel("z (km)")
    ax.set_title(f"{run_info.run_name} mean flow deviation u' = u_bar(z,t) - {run_info.u_bar:.1f} (m/s)")
    fig.colorbar(im, ax=ax, label="u' (m/s)")

    outpath = os.path.join(run_output_dir, "mean_flow_evolution.png")
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)


def analyze_single_run(run_info: RunInfo, args: argparse.Namespace) -> tuple[bool, str]:
    run_output_dir = os.path.join(args.outdir, run_info.run_name)
    os.makedirs(run_output_dir, exist_ok=True)

    files, times = get_snapshot_files_and_times(run_info.run_dir)
    if len(files) == 0:
        write_failed_metrics(run_output_dir, run_info, "no snapshots found")
        return False, "no snapshots"

    try:
        with Dataset(files[-1], "r") as ds_last:
            x_m = np.asarray(ds_last.variables["x"][:], dtype=np.float64)
            z_m = np.asarray(ds_last.variables["z"][:], dtype=np.float64)
            w = np.asarray(ds_last.variables["w"][:, :], dtype=np.float64)
            rho = np.asarray(ds_last.variables["rho"][:, :], dtype=np.float64)
            t_last = float(ds_last.getncattr("time_s"))

        plot_latest_snapshot(
            run_output_dir=run_output_dir,
            run_info=run_info,
            x_m=x_m,
            z_m=z_m,
            w=w,
            rho=rho,
            time_s=t_last,
            z_max_km=args.snapshot_z_max_km,
            x_limit_km=args.snapshot_x_limit_km,
        )

        station_times_s, station_x_actual_km, station_z_actual_km, station_arrays = extract_station_series(
            files=files,
            times=times,
            station_x_rel_km=parse_station_x_list(args.station_x_rel_km),
            station_z_km=args.station_z_km,
        )

        plot_station_timeseries(
            run_output_dir=run_output_dir,
            run_info=run_info,
            times_s=station_times_s,
            station_x_actual_km=station_x_actual_km,
            station_z_actual_km=station_z_actual_km,
            station_arrays=station_arrays,
        )

        station_keys = sorted(station_arrays.keys())
        wavelet_idx = int(np.clip(args.wavelet_station_index, 0, len(station_keys) - 1))
        wavelet_key = station_keys[wavelet_idx]
        wavelet_signal = station_arrays[wavelet_key]

        if len(station_times_s) < 4:
            raise ValueError("Not enough snapshots for wavelet analysis (need >= 4)")

        # --- Transient Diagnostic ---
        # transient_ratio = rms(first 10%) / rms(remaining 90%)
        split_idx = max(3, int(0.1 * len(wavelet_signal)))
        rms_early = float(np.sqrt(np.mean(wavelet_signal[:split_idx]**2)))
        rms_late = float(np.sqrt(np.mean(wavelet_signal[split_idx:]**2)))
        transient_ratio = rms_early / max(rms_late, 1.0e-12)

        # --- Trimming ---
        trim_s = args.wavelet_trim_min * 60.0
        keep = station_times_s >= (station_times_s[0] + trim_s)
        if np.count_nonzero(keep) < 4:
            # Fallback if trimming removes too much
            keep = np.arange(len(station_times_s)) >= int(0.1 * len(station_times_s))
        
        station_times_trimmed = station_times_s[keep]
        wavelet_signal_trimmed = wavelet_signal[keep]

        # --- Drift Diagnostic (on trimmed signal) ---
        drift_ratio = float(abs(np.mean(wavelet_signal_trimmed)) / max(np.std(wavelet_signal_trimmed), 1.0e-12))

        dt_est = float(np.median(np.diff(station_times_s)))
        if dt_est <= 0.0:
            raise ValueError("Non-positive time step detected")

        power, period_s = compute_wavelet(
            signal=wavelet_signal_trimmed,
            dt_s=dt_est,
            dj=args.wavelet_dj,
            morlet_param=args.wavelet_morlet_param,
        )
        dominant_freq_mhz, dominant_period_s = plot_wavelet(
            run_output_dir=run_output_dir,
            run_info=run_info,
            times_s=station_times_trimmed,
            power=power,
            period_s=period_s,
            station_label=wavelet_key,
            min_freq_mhz=args.wavelet_min_freq_mhz,
            max_freq_mhz=args.wavelet_max_freq_mhz,
        )

        plot_mean_flow_evolution(
            run_output_dir=run_output_dir,
            run_info=run_info,
            files=files,
            times_s=station_times_s,
        )

        metrics = {
            "status": "ok",
            "run_name": run_info.run_name,
            "U_bar_mps": run_info.u_bar,
            "sigma_km": run_info.sigma_km,
            "h_max_km": run_info.h_max_km,
            "n_snapshots": int(len(files)),
            "t_start_s": float(station_times_s[0]),
            "t_end_s": float(station_times_s[-1]),
            "dt_median_s": dt_est,
            "snapshot_time_s": t_last,
            "station_z_target_km": float(args.station_z_km),
            "station_z_actual_km": station_z_actual_km,
            "station_x_target_rel_km": parse_station_x_list(args.station_x_rel_km),
            "station_x_actual_rel_km": station_x_actual_km,
            "wavelet_station": wavelet_key,
            "wavelet_trim_min": args.wavelet_trim_min,
            "transient_ratio": transient_ratio,
            "drift_ratio": drift_ratio,
            "dominant_freq_mhz": dominant_freq_mhz,
            "dominant_period_s": dominant_period_s,
            "acoustic_cutoff_mhz": 1.8,
        }
        with open(os.path.join(run_output_dir, "metrics.json"), "w", encoding="utf-8") as f:
            json.dump(metrics, f, indent=2)

        return True, "ok"
    except Exception as exc:
        write_failed_metrics(run_output_dir, run_info, str(exc))
        return False, str(exc)


# =============================================================================
# Sweep-level composite (opt-in via --composite)
# =============================================================================

def load_sweep_results(sweep_dir: str) -> dict[tuple[float, float], dict[str, str]]:
    """Read sweep_results.csv into a dict keyed by (U_bar, h_max_km).

    Pulls Nh_over_U and the per-run wall-clock duration that the Julia runner
    already recorded, so the analyzer does not re-derive any physics constants.
    """
    path = os.path.join(sweep_dir, "sweep_results.csv")
    table: dict[tuple[float, float], dict[str, str]] = {}
    if not os.path.exists(path):
        return table
    with open(path, "r", encoding="utf-8") as f:
        header = f.readline().strip().split(",")
        idx = {name: i for i, name in enumerate(header)}
        if "U_bar" not in idx or "h_max" not in idx:
            return table
        for line in f:
            parts = line.strip().split(",")
            if len(parts) < len(header):
                continue
            try:
                u = float(parts[idx["U_bar"]])
                h_km = float(parts[idx["h_max"]]) / 1000.0
            except ValueError:
                continue
            rec = {col: parts[idx[col]] for col in ("Nh_over_U", "duration_s",
                                                    "final_t", "steps", "converged")
                   if col in idx}
            table[(round(u, 3), round(h_km, 3))] = rec
    return table


def load_final_snapshot(run_dir: str):
    """Return (x_m, z_m, w, rho, t_last) for the latest snapshot, or None.

    Snapshot filenames are zero-padded by step (snapshot_<step>.nc), so a
    lexicographic sort already orders them by time. Taking the last name avoids
    opening every file just to read its time attribute (which would make the
    sweep composite as slow as the full per-run pass over an external drive).
    """
    files = sorted(glob.glob(os.path.join(run_dir, "*.nc")))
    if not files:
        return None
    with Dataset(files[-1], "r") as ds:
        x_m = np.asarray(ds.variables["x"][:], dtype=np.float64)
        z_m = np.asarray(ds.variables["z"][:], dtype=np.float64)
        w = np.asarray(ds.variables["w"][:, :], dtype=np.float64)
        rho = np.asarray(ds.variables["rho"][:, :], dtype=np.float64)
        t_last = float(ds.getncattr("time_s"))
    return x_m, z_m, w, rho, t_last


def compute_band_response(runs: list[RunInfo], args: argparse.Namespace) -> list[dict]:
    """Per run, peak/mean |w| over the response band at the final time.

    Reports the raw vertical velocity (not density-scaled): the actual updraft is
    what drives the adiabatic cooling relevant to cloud formation at the band.

    The primary peak/mean are computed inside the hill window (|x_rel| <=
    composite_x_limit_km), the same horizontal extent shown in the composite
    panels. The unrestricted global peak |w| and its x-location are reported
    separately so a boundary feature (e.g. a stationary inflow-edge perturbation
    advected far upwind of the hill) is documented rather than silently inflating
    the reported orographic amplitude.

    Writes <band>_response.csv and .md into args.outdir and returns the rows.
    """
    sweep = load_sweep_results(args.sweep_dir)
    z_lo = args.band_z_min_km * 1000.0
    z_hi = args.band_z_max_km * 1000.0
    x_win = args.composite_x_limit_km
    rows: list[dict] = []
    for r in runs:
        data = load_final_snapshot(r.run_dir)
        if data is None:
            continue
        x_m, z_m, w, _, t_last = data
        band = (z_m >= z_lo) & (z_m <= z_hi)
        if not np.any(band):
            continue
        x_rel_km = build_centered_x_km(x_m)
        wb = np.abs(w[band, :])
        in_win = np.abs(x_rel_km) <= x_win
        wb_hill = wb[:, in_win]
        jx_glob = int(np.unravel_index(np.argmax(wb), wb.shape)[1])
        key = (round(r.u_bar, 3), round(r.h_max_km, 3)) if r.h_max_km is not None else None
        rec = sweep.get(key, {}) if key is not None else {}
        rows.append({
            "run_name": r.run_name,
            "U_bar": r.u_bar,
            "h_max_km": r.h_max_km,
            "Nh_over_U": rec.get("Nh_over_U", ""),
            "peak_w_hill_mps": float(np.max(wb_hill)),
            "mean_w_hill_mps": float(np.mean(wb_hill)),
            "peak_w_global_mps": float(np.max(wb)),
            "peak_x_global_km": float(x_rel_km[jx_glob]),
            "t_last_s": t_last,
            "duration_s": rec.get("duration_s", ""),
        })

    rows.sort(key=lambda d: (d["U_bar"], d["h_max_km"] if d["h_max_km"] is not None else 0.0))

    base = f"{args.band_label.lower()}_response"
    csv_path = os.path.join(args.outdir, base + ".csv")
    cols = ["run_name", "U_bar", "h_max_km", "Nh_over_U",
            "peak_w_hill_mps", "mean_w_hill_mps",
            "peak_w_global_mps", "peak_x_global_km", "t_last_s", "duration_s"]
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write(",".join(cols) + "\n")
        for d in rows:
            f.write(",".join(str(d[c]) for c in cols) + "\n")

    md_path = os.path.join(args.outdir, base + ".md")
    with open(md_path, "w", encoding="utf-8") as f:
        f.write(f"Hill window: |x_rel| <= {x_win:.0f} km; "
                f"{args.band_label} band {args.band_z_min_km:.0f}-{args.band_z_max_km:.0f} km.\n\n")
        f.write(f"| U_bar | h_max (km) | Nh/U | hill peak |w| (m/s) | "
                f"hill mean |w| (m/s) | global peak |w| (m/s) | global peak x (km) | "
                f"wall-clock (s) |\n")
        f.write("|---|---|---|---|---|---|---|---|\n")
        for d in rows:
            nh = d["Nh_over_U"]
            nh_s = f"{float(nh):.2f}" if nh not in ("", None) else "n/a"
            dur = d["duration_s"]
            dur_s = f"{float(dur):.0f}" if dur not in ("", None) else "n/a"
            f.write(f"| {d['U_bar']:.0f} | {d['h_max_km']:.1f} | {nh_s} | "
                    f"{d['peak_w_hill_mps']:.2f} | {d['mean_w_hill_mps']:.3f} | "
                    f"{d['peak_w_global_mps']:.2f} | {d['peak_x_global_km']:.0f} | {dur_s} |\n")

    print(f"  band response -> {csv_path}")
    print(f"  band response -> {md_path}")
    return rows


def plot_sweep_composite(runs: list[RunInfo], args: argparse.Namespace) -> None:
    """Grid of final-time w*sqrt(rho/rho0) panels, one per (U_bar, h_max) run.

    Density-scaled w removes the trivial rho^-1/2 growth so panels show true
    wave activity on a shared, directly-comparable color scale. Rows are h_max
    (largest on top), columns are U_bar (ascending).
    """
    u_vals = sorted({r.u_bar for r in runs})
    h_vals = sorted({r.h_max_km for r in runs if r.h_max_km is not None})
    if not u_vals or not h_vals:
        print("  composite skipped: runs are not on a U_bar x h_max grid")
        return

    sweep = load_sweep_results(args.sweep_dir)
    z_max = args.composite_z_max_km
    x_lim = args.composite_x_limit_km

    panel: dict[tuple[float, float], tuple] = {}
    pooled: list[np.ndarray] = []
    for r in runs:
        if r.h_max_km is None:
            continue
        data = load_final_snapshot(r.run_dir)
        if data is None:
            continue
        x_m, z_m, w, rho, t_last = data
        xx_km = build_centered_x_km(x_m)
        z_km = z_m / 1000.0
        rho_surface = float(np.mean(rho[0, :]))
        w_scaled = w * np.sqrt(rho / rho_surface)
        iz = np.where(z_km <= z_max)[0]
        ix = np.where(np.abs(xx_km) <= x_lim)[0]
        sub = w_scaled[np.ix_(iz, ix)]
        panel[(r.u_bar, r.h_max_km)] = (xx_km[ix], z_km[iz], sub, t_last)
        pooled.append(np.abs(sub).ravel())

    if not pooled:
        print("  composite skipped: no snapshots found")
        return

    vmax = float(np.quantile(np.concatenate(pooled), 0.995))
    vmax = max(vmax, 1.0e-8)
    levels = np.linspace(-vmax, vmax, 41)

    nrows, ncols = len(h_vals), len(u_vals)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.2 * ncols, 3.2 * nrows),
                             sharex=True, sharey=True, squeeze=False)
    im = None
    for ri, h in enumerate(reversed(h_vals)):       # largest hill on top
        for ci, u in enumerate(u_vals):
            ax = axes[ri][ci]
            key = (u, h)
            if key not in panel:
                ax.set_visible(False)
                continue
            xx, zz, sub, t_last = panel[key]
            im = ax.contourf(xx, zz, sub, levels=levels, cmap="RdBu_r", extend="both")
            ax.axhspan(args.band_z_min_km, args.band_z_max_km, color="k", alpha=0.10)
            ax.axhline(args.sponge_base_km, color="0.35", ls="--", lw=0.9)
            rec = sweep.get((round(u, 3), round(h, 3)), {})
            nh = rec.get("Nh_over_U", "")
            nh_s = f", Nh/U={float(nh):.2f}" if nh not in ("", None) else ""
            ax.set_title(f"U={u:.0f} m/s, h={h:.1f} km{nh_s}", fontsize=9)
            if ci == 0:
                ax.set_ylabel("z (km)")
            if ri == nrows - 1:
                ax.set_xlabel("x rel. center (km)")

    if im is not None:
        fig.colorbar(im, ax=axes.ravel().tolist(),
                     label=r"$w\,\sqrt{\rho/\rho_0}$ (m/s)", shrink=0.85)
    fig.suptitle(
        f"Arsia Mons Stage 1 sweep | final-time $w\\sqrt{{\\rho/\\rho_0}}$ | "
        f"{args.band_label} band {args.band_z_min_km:.0f}-{args.band_z_max_km:.0f} km shaded, "
        f"sponge base {args.sponge_base_km:.0f} km dashed"
    )
    outpath = os.path.join(args.outdir, "sweep_composite_w.png")
    plt.savefig(outpath, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  composite -> {outpath}")


def discover_runs(sweep_dir: str, run_pattern: str) -> list[RunInfo]:
    dirs = sorted(glob.glob(os.path.join(sweep_dir, run_pattern)))
    out = []
    for d in dirs:
        info = parse_run_info(d)
        if info is not None:
            out.append(info)

    out.sort(key=lambda r: (r.u_bar, r.sigma_km, r.run_name))
    return out


def main() -> None:
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    runs = discover_runs(args.sweep_dir, args.run_pattern)
    if not runs:
        raise ValueError(f"No run directories found in {args.sweep_dir}")

    station_targets = parse_station_x_list(args.station_x_rel_km)
    print(f"Found {len(runs)} runs in {args.sweep_dir}")
    print(f"Station x targets (km, relative to center): {station_targets}")
    print(f"Station z target: {args.station_z_km:.2f} km")

    if args.composite_only:
        print("Composite-only: skipping per-run analysis")
    else:
        n_ok = 0
        n_fail = 0
        for i, run_info in enumerate(runs, start=1):
            print(f"[{i}/{len(runs)}] Processing {run_info.run_name} ...")
            ok, msg = analyze_single_run(run_info, args)
            if ok:
                n_ok += 1
                print(f"  done: {msg}")
            else:
                n_fail += 1
                print(f"  failed: {msg}")
        print(f"Completed. successful={n_ok}, failed={n_fail}")

    if args.composite or args.composite_only:
        print(f"Building sweep-level composite ({args.band_label} band "
              f"{args.band_z_min_km:.0f}-{args.band_z_max_km:.0f} km) ...")
        compute_band_response(runs, args)
        plot_sweep_composite(runs, args)

    print(f"Outputs in: {args.outdir}")


if __name__ == "__main__":
    main()
