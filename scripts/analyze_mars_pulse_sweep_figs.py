#!/usr/bin/env python3
"""
Batch generation of Fig 12-15 style diagnostics for every Mars pulse sweep run.

For each folder output/sweep_mars_pulse/run_sx*_sz*:
- Fig 12: pressure perturbation snapshot
- Fig 13: station pressure time series (S1-S4)
- Fig 14: wavelet spectrogram at S3 (or fallback panel if too sparse)
- Fig 15: arrival-time and attenuation fits
"""

import argparse
import csv
import glob
import os
import re
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from scipy.optimize import curve_fit

from wavelet_tc98 import wavelet


RUN_DIR_PATTERN = re.compile(r"run_sx(?P<sx>[-+0-9.]+)_sz(?P<sz>[-+0-9.]+)_pt(?P<pt>[-+0-9.]+)_z0(?P<z0>[-+0-9.]+)$")

plt.rcParams.update({
    "font.size": 12,
    "axes.labelsize": 14,
    "axes.titlesize": 16,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
})


@dataclass
class RunSummary:
    run_name: str
    sigma_x_m: float
    sigma_z_m: float
    p_tilde_Pa: float
    z0_km: float
    n_snapshots: int
    t_min_s: float
    t_max_s: float
    fig12_snapshot_time_s: float
    fig14_mode: str
    fig15_speed_slope_s_per_km: float
    fig15_speed_intercept_s: float
    fig15_speed_r2: float
    fig15_att_a0: float
    fig15_att_c: float
    fig15_att_r2: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Fig12-15 for every sweep run")
    parser.add_argument("--sweep-dir", default="output/sweep_mars_pulse", help="Directory with run_sx*_sz* folders")
    parser.add_argument("--outdir", default="output/mars_pulse_sweep_figs", help="Output directory for per-run figures")
    parser.add_argument("--source-x-km", type=float, default=100.0, help="Pulse center x location in km")
    parser.add_argument("--target-time-s", type=float, default=95.0, help="Target time for Fig 12")
    parser.add_argument("--snapshot-tol-s", type=float, default=10.0, help="Tolerance around target time for Fig 12")
    return parser.parse_args()


def parse_params(run_dir: str) -> tuple[float, float, float, float]:
    name = os.path.basename(os.path.normpath(run_dir))
    m = RUN_DIR_PATTERN.match(name)
    if m is None:
        raise ValueError(f"Could not parse parameters from {run_dir}")
    return (
        float(m.group("sx")),
        float(m.group("sz")),
        float(m.group("pt")),
        float(m.group("z0"))
    )


def build_file_time_index(indir: str) -> tuple[list[str], np.ndarray]:
    files = sorted(glob.glob(os.path.join(indir, "*.nc")))
    if not files:
        raise ValueError(f"No files found in {indir}")

    times = np.empty(len(files), dtype=np.float64)
    for i, f in enumerate(files):
        with Dataset(f, "r") as ds:
            times[i] = float(ds.getncattr("time_s"))

    order = np.argsort(times)
    return [files[i] for i in order], times[order]


def get_snapshot_for_fig12(files: list[str], times: np.ndarray, target_time_s: float, tolerance_s: float) -> tuple[str, float, bool]:
    idx = int(np.argmin(np.abs(times - target_time_s)))
    snap_t = float(times[idx])
    strict = abs(snap_t - target_time_s) <= tolerance_s
    if (not strict) and idx == 0 and len(times) > 1:
        idx = 1
        snap_t = float(times[idx])
    return files[idx], snap_t, strict


def extract_station_and_front_data(files: list[str], times_index: np.ndarray, source_x_km: float):
    with Dataset(files[0], "r") as ds:
        x = np.asarray(ds.variables["x"][:], dtype=np.float64) / 1000.0
        z = np.asarray(ds.variables["z"][:], dtype=np.float64) / 1000.0

    station_offsets = np.array([20.0, 40.0, 60.0, 80.0], dtype=np.float64)
    station_x = source_x_km + station_offsets

    ix_s = [int(np.argmin(np.abs(x - sx))) for sx in station_x]
    iz_20 = int(np.argmin(np.abs(z - 20.0)))

    ix_start = int(np.argmin(np.abs(x - source_x_km)))
    x_right = x[ix_start:] - source_x_km

    max_amps = np.zeros(len(x_right), dtype=np.float64)
    time_of_max = np.zeros(len(x_right), dtype=np.float64)

    times_out = []
    p_s = [[] for _ in range(4)]

    for i, f in enumerate(files):
        t = float(times_index[i])
        with Dataset(f, "r") as ds:
            p_row = np.asarray(ds.variables["P"][iz_20, :], dtype=np.float64)

        p_bar = np.mean(p_row)
        p_pert = p_row - p_bar

        times_out.append(t)
        for j in range(4):
            p_s[j].append(float(p_pert[ix_s[j]]))

        p_pert_abs_right = np.abs(p_pert[ix_start:])
        mask = p_pert_abs_right > max_amps
        max_amps[mask] = p_pert_abs_right[mask]
        time_of_max[mask] = t

    return (
        np.array(times_out),
        np.array(p_s[0]),
        np.array(p_s[1]),
        np.array(p_s[2]),
        np.array(p_s[3]),
        x_right,
        time_of_max,
        max_amps,
    )


def plot_fig12(files: list[str], times: np.ndarray, outdir: str, source_x_km: float, target_time_s: float, tol_s: float) -> float:
    fname, used_t, is_strict = get_snapshot_for_fig12(files, times, target_time_s, tol_s)

    with Dataset(fname, "r") as ds:
        x = np.asarray(ds.variables["x"][:], dtype=np.float64) / 1000.0
        z = np.asarray(ds.variables["z"][:], dtype=np.float64) / 1000.0
        p = np.asarray(ds.variables["P"][:, :], dtype=np.float64)

    p_bar = np.mean(p, axis=1, keepdims=True)
    p_pert = p - p_bar

    x_shift = x - source_x_km
    X, Z = np.meshgrid(x_shift, z)

    fig, ax = plt.subplots(figsize=(8, 3.5))
    vabs = max(1e-12, float(np.max(np.abs(p_pert))))
    levels = np.linspace(-vabs, vabs, 101)
    im = ax.contourf(X, Z, p_pert, levels=levels, cmap="gray", extend="both")

    stations_x = [20, 40, 60, 80]
    colors = ["tab:blue", "tab:orange", "tab:green", "tab:purple"]
    for i, sx in enumerate(stations_x):
        ax.plot(sx, 20.0, "x", color=colors[i], markersize=10, markeredgewidth=3)
        ax.text(sx - 3, 22.0, f"S{i+1}", fontweight="bold", fontsize=12)

    if is_strict:
        txt = f"T = {used_t:.1f} s"
    else:
        txt = f"T = {used_t:.1f} s (closest to {target_time_s:.1f} s)"
    ax.text(-95, 32, txt, fontsize=11, fontweight="bold")

    ax.set_xlabel("x (km)")
    ax.set_ylabel("z (km)")
    ax.set_xlim(-100, 100)
    ax.set_ylim(0, 40)
    ax.set_xticks([-50, 0, 50])

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Pa", loc="top", rotation=0, labelpad=-20)

    outpath = os.path.join(outdir, "mars_fig12.png")
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return used_t


def plot_fig13(times, p_s1, p_s2, p_s3, p_s4, outdir):
    fig, axes = plt.subplots(4, 1, figsize=(10, 8), sharex=True)

    series = [p_s1, p_s2, p_s3, p_s4]
    labels = [
        "S1 (20 km from source)",
        "S2 (40 km from source)",
        "S3 (60 km from source)",
        "S4 (80 km from source)",
    ]

    x_max = float(np.max(times)) if len(times) > 0 else 1.0
    y_abs = max(1e-12, float(np.max(np.abs(np.concatenate(series)))))

    for i in range(4):
        axes[i].plot(times, series[i], "tab:blue", linewidth=1.2, marker="o")
        axes[i].set_ylim(-1.1 * y_abs, 1.1 * y_abs)
        axes[i].set_xlim(float(np.min(times)), x_max)
        axes[i].set_ylabel("Pa", fontsize=12)
        axes[i].text(x_max * 0.98, 0.65 * y_abs, labels[i], horizontalalignment="right", fontsize=10, fontweight="bold")

    axes[3].set_xlabel("t (s)", fontsize=14)
    plt.subplots_adjust(hspace=0.1)

    outpath = os.path.join(outdir, "mars_fig13.png")
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_fig14(times, p_s3, outdir):
    fig, ax = plt.subplots(figsize=(8, 6))

    if len(times) < 8:
        ax.plot(times / 60.0, p_s3, "o-", color="tab:blue")
        ax.set_xlabel("time (min)")
        ax.set_ylabel("P' (Pa)")
        ax.set_title("Fig 14 fallback: insufficient temporal samples for wavelet")
        ax.text(0.03, 0.93, f"Need >= 8 samples, got {len(times)}", transform=ax.transAxes, va="top")
        fig14_mode = "fallback_sparse"
    else:
        dt = float(np.median(np.diff(times)))
        y = p_s3 - np.mean(p_s3)
        dj = 0.05
        s0 = 2 * dt
        J1 = int(np.fix(np.log2(len(times) * dt / s0) / dj))

        if J1 < 1:
            ax.plot(times / 60.0, p_s3, "o-", color="tab:blue")
            ax.set_xlabel("time (min)")
            ax.set_ylabel("P' (Pa)")
            ax.set_title("Fig 14 fallback: insufficient scale depth for wavelet")
            fig14_mode = "fallback_scales"
        else:
            wave, period, _, _ = wavelet(y, dt, pad=1, dj=dj, s0=s0, J1=J1, mother="MORLET", param=6.0)
            power = np.abs(wave) ** 2
            pmax = float(np.max(power))
            if pmax > 0:
                power = power / pmax

            freq_mHz = (1.0 / period) * 1000.0
            times_min = times / 60.0

            im = ax.pcolormesh(times_min, freq_mHz, power, cmap="jet", shading="auto")
            ax.axhline(1.8, color="limegreen", linestyle="--", lw=3)
            ax.set_ylim(0, 30)
            ax.set_xlim(float(np.min(times_min)), float(np.max(times_min)))
            ax.set_xlabel("time (min)", fontsize=14)
            ax.set_ylabel("frequency (mHz)", fontsize=14)
            fig.colorbar(im, ax=ax)
            fig14_mode = "wavelet"

    outpath = os.path.join(outdir, "mars_fig14.png")
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return fig14_mode


def _r2(y_true: np.ndarray, y_fit: np.ndarray) -> float:
    denom = np.sum((y_true - np.mean(y_true)) ** 2)
    if denom <= 0.0:
        return np.nan
    return float(1.0 - np.sum((y_true - y_fit) ** 2) / denom)


def plot_fig15(x_right, time_of_max, max_amps, outdir):
    valid_idx = np.where(x_right > 1.0)[0]
    xr = x_right[valid_idx]
    tm = time_of_max[valid_idx]
    ma = max_amps[valid_idx]

    def linear_fit(x, m, c):
        return m * x + c

    def geom_atten(x, a0, c):
        return a0 / np.sqrt(x) + c

    popt_v = np.array([np.nan, np.nan], dtype=np.float64)
    r2_v = np.nan
    popt_a = np.array([np.nan, np.nan], dtype=np.float64)
    r2_a = np.nan

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    ax1.plot(xr, tm, "k.", markersize=6, label="Simulation data")
    if len(xr) >= 3 and float(np.std(tm)) > 1e-12:
        try:
            popt_v, _ = curve_fit(linear_fit, xr, tm)
            fit_v = linear_fit(xr, *popt_v)
            r2_v = _r2(tm, fit_v)
            ax1.plot(xr, fit_v, "g-", linewidth=2, label="Curve fit")
            ax1.text(0.03, 0.95, f"y = {popt_v[0]:.3f}x + {popt_v[1]:.1f}\n$R^2$ = {r2_v:.4f}", transform=ax1.transAxes, va="top")
        except Exception:
            ax1.text(0.03, 0.95, "Linear fit failed", transform=ax1.transAxes, va="top")
    else:
        ax1.text(0.03, 0.95, "Insufficient variation for linear fit", transform=ax1.transAxes, va="top")

    ax1.set_xlabel("distance from source (km)")
    ax1.set_ylabel("time of arrival (s)")
    ax1.set_xlim(0, max(100, float(np.max(xr)) if len(xr) else 100))
    ax1.legend(loc="lower right")

    ax2.plot(xr, ma, "k.", markersize=6, label="Simulation data")
    if len(xr) >= 3 and float(np.std(ma)) > 1e-12:
        try:
            popt_a, _ = curve_fit(geom_atten, xr, ma, p0=[1.0, 0.0], maxfev=10000)
            fit_a = geom_atten(xr, *popt_a)
            r2_a = _r2(ma, fit_a)
            ax2.plot(xr, fit_a, "b-", linewidth=2, label="Curve fit")
            ax2.text(0.03, 0.95, f"y = {popt_a[0]:.2f}/$\\sqrt{{x}}$ {popt_a[1]:+.2f}\n$R^2$ = {r2_a:.4f}", transform=ax2.transAxes, va="top")
        except Exception:
            ax2.text(0.03, 0.95, "Attenuation fit failed", transform=ax2.transAxes, va="top")
    else:
        ax2.text(0.03, 0.95, "Insufficient variation for attenuation fit", transform=ax2.transAxes, va="top")

    ax2.set_xlabel("distance from source (km)")
    ax2.set_ylabel("Amplitude (Pa)")
    ax2.set_xlim(0, max(100, float(np.max(xr)) if len(xr) else 100))
    ax2.legend(loc="upper right")

    outpath = os.path.join(outdir, "mars_fig15.png")
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)

    return popt_v, r2_v, popt_a, r2_a


def write_summary(out_csv: str, rows: list[RunSummary]) -> None:
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([
            "run_name", "sigma_x_m", "sigma_z_m", "n_snapshots", "t_min_s", "t_max_s",
            "fig12_snapshot_time_s", "fig14_mode",
            "fig15_speed_slope_s_per_km", "fig15_speed_intercept_s", "fig15_speed_r2",
            "fig15_att_a0", "fig15_att_c", "fig15_att_r2",
        ])
        for r in rows:
            w.writerow([
                r.run_name, r.sigma_x_m, r.sigma_z_m, r.n_snapshots, r.t_min_s, r.t_max_s,
                r.fig12_snapshot_time_s, r.fig14_mode,
                r.fig15_speed_slope_s_per_km, r.fig15_speed_intercept_s, r.fig15_speed_r2,
                r.fig15_att_a0, r.fig15_att_c, r.fig15_att_r2,
            ])


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    run_dirs = sorted(glob.glob(os.path.join(args.sweep_dir, "run_sx*_sz*")))
    if not run_dirs:
        raise FileNotFoundError(f"No run directories found in {args.sweep_dir}")

    summaries: list[RunSummary] = []

    for run_dir in run_dirs:
        run_name = os.path.basename(os.path.normpath(run_dir))
        try:
            sigma_x_m, sigma_z_m, p_tilde_Pa, z0_km = parse_params(run_dir)
        except ValueError as e:
            print(f"Error parsing parameters for {run_dir}: {e}")
            continue
        run_outdir = os.path.join(args.outdir, run_name)
        os.makedirs(run_outdir, exist_ok=True)

        files, times = build_file_time_index(run_dir)

        fig12_t = plot_fig12(
            files=files,
            times=times,
            outdir=run_outdir,
            source_x_km=args.source_x_km,
            target_time_s=args.target_time_s,
            tol_s=args.snapshot_tol_s,
        )

        (
            t_arr,
            p_s1,
            p_s2,
            p_s3,
            p_s4,
            x_right,
            time_of_max,
            max_amps,
        ) = extract_station_and_front_data(files, times, source_x_km=args.source_x_km)

        plot_fig13(t_arr, p_s1, p_s2, p_s3, p_s4, run_outdir)
        fig14_mode = plot_fig14(t_arr, p_s3, run_outdir)
        popt_v, r2_v, popt_a, r2_a = plot_fig15(x_right, time_of_max, max_amps, run_outdir)

        summaries.append(
            RunSummary(
                run_name=run_name,
                sigma_x_m=sigma_x_m,
                sigma_z_m=sigma_z_m,
                p_tilde_Pa=p_tilde_Pa,
                z0_km=z0_km,
                n_snapshots=len(files),
                t_min_s=float(np.min(times)),
                t_max_s=float(np.max(times)),
                fig12_snapshot_time_s=fig12_t,
                fig14_mode=fig14_mode,
                fig15_speed_slope_s_per_km=float(popt_v[0]),
                fig15_speed_intercept_s=float(popt_v[1]),
                fig15_speed_r2=float(r2_v),
                fig15_att_a0=float(popt_a[0]),
                fig15_att_c=float(popt_a[1]),
                fig15_att_r2=float(r2_a),
            )
        )

        print(f"Processed {run_name}: Fig12-15 written to {run_outdir}")

    summary_csv = os.path.join(args.outdir, "sweep_fig12_15_summary.csv")
    write_summary(summary_csv, summaries)
    print(f"Summary written: {summary_csv}")


if __name__ == "__main__":
    main()
