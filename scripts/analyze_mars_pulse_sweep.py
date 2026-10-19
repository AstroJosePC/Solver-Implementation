#!/usr/bin/env python3
"""
analyze_mars_pulse_sweep.py
Analyze sweep outputs from scripts/sweep_mars_pulse.jl.

Produces:
- summary CSV with one row per (sigma_x, sigma_z) run
- per-run station time series CSV
- heatmaps for key metrics across the sweep grid
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


RUN_DIR_PATTERN = re.compile(r"run_sx(?P<sx>[-+0-9.]+)_sz(?P<sz>[-+0-9.]+)$")


@dataclass
class RunMetrics:
    sigma_x_m: float
    sigma_z_m: float
    n_snapshots: int
    t_start_s: float
    t_end_s: float
    dt_mean_s: float
    p_rms_final_pa: float
    p_peak_final_pa: float
    p_rms_mean_pa: float
    p_peak_mean_pa: float
    ke_proxy_mean: float
    lowk_frac_mean: float
    station_std_pa: float
    station_peak_pa: float
    station_trough_pa: float
    station_ptp_pa: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Analyze Mars pulse sweep outputs")
    parser.add_argument("--sweep-dir", default="output/sweep_mars_pulse", help="Directory with run_sx*_sz* folders")
    parser.add_argument("--outdir", default="output/mars_pulse_sweep_analysis", help="Directory for analysis outputs")
    parser.add_argument("--large-scale-lambda-km", type=float, default=30.0,
                        help="Horizontal wavelength threshold (km) for low-k fraction")
    parser.add_argument("--station-x-km", type=float, default=160.0, help="Station x position in km")
    parser.add_argument("--station-z-km", type=float, default=20.0, help="Station z position in km")
    return parser.parse_args()


def get_snapshot_files(run_dir: str) -> tuple[list[str], np.ndarray]:
    files = sorted(glob.glob(os.path.join(run_dir, "*.nc")))
    if not files:
        return [], np.array([], dtype=np.float64)

    times = np.empty(len(files), dtype=np.float64)
    for i, fpath in enumerate(files):
        with Dataset(fpath, "r") as ds:
            times[i] = float(ds.getncattr("time_s"))

    order = np.argsort(times)
    return [files[i] for i in order], times[order]


def parse_sigma_from_run_dir(run_dir: str) -> tuple[float, float] | None:
    name = os.path.basename(os.path.normpath(run_dir))
    match = RUN_DIR_PATTERN.match(name)
    if not match:
        return None
    return float(match.group("sx")), float(match.group("sz"))


def lowk_power_fraction(signal: np.ndarray, dx_m: float, lambda_threshold_m: float) -> float:
    if signal.size < 2:
        return 0.0

    signal = signal - np.mean(signal)
    fft_vals = np.fft.rfft(signal)
    power = np.abs(fft_vals) ** 2
    freqs = np.fft.rfftfreq(signal.size, d=dx_m)

    if power.size <= 1:
        return 0.0

    nonzero = freqs > 0
    if not np.any(nonzero):
        return 0.0

    wavelength = np.zeros_like(freqs)
    wavelength[nonzero] = 1.0 / freqs[nonzero]
    lowk = nonzero & (wavelength >= lambda_threshold_m)

    denom = np.sum(power[nonzero])
    if denom <= 0.0:
        return 0.0

    return float(np.sum(power[lowk]) / denom)


def analyze_single_run(
    run_dir: str,
    sigma_x_m: float,
    sigma_z_m: float,
    lambda_threshold_m: float,
    station_x_km: float,
    station_z_km: float,
) -> tuple[RunMetrics, np.ndarray, np.ndarray]:
    files, times = get_snapshot_files(run_dir)
    if len(files) == 0:
        raise ValueError(f"No snapshots found in {run_dir}")

    p_rms_all = []
    p_peak_all = []
    ke_all = []
    lowk_all = []
    station_series = []

    with Dataset(files[0], "r") as ds0:
        x_km = np.asarray(ds0.variables["x"][:], dtype=np.float64) / 1000.0
        z_km = np.asarray(ds0.variables["z"][:], dtype=np.float64) / 1000.0

    ix_station = int(np.argmin(np.abs(x_km - station_x_km)))
    iz_station = int(np.argmin(np.abs(z_km - station_z_km)))

    for fpath in files:
        with Dataset(fpath, "r") as ds:
            p = np.asarray(ds.variables["P"][:, :], dtype=np.float64)
            u = np.asarray(ds.variables["u"][:, :], dtype=np.float64)
            w = np.asarray(ds.variables["w"][:, :], dtype=np.float64)
            x_m = np.asarray(ds.variables["x"][:], dtype=np.float64)

        p_pert = p - np.mean(p, axis=1, keepdims=True)

        p_rms_all.append(float(np.sqrt(np.mean(p_pert ** 2))))
        p_peak_all.append(float(np.max(np.abs(p_pert))))
        ke_all.append(float(np.mean(0.5 * (u ** 2 + w ** 2))))

        dx_m = float(np.median(np.diff(x_m)))
        row = p_pert[iz_station, :]
        lowk_all.append(lowk_power_fraction(row, dx_m=dx_m, lambda_threshold_m=lambda_threshold_m))
        station_series.append(float(row[ix_station]))

    station_arr = np.array(station_series, dtype=np.float64)
    dt_mean = float(np.mean(np.diff(times))) if len(times) > 1 else np.nan

    metrics = RunMetrics(
        sigma_x_m=sigma_x_m,
        sigma_z_m=sigma_z_m,
        n_snapshots=len(files),
        t_start_s=float(times[0]),
        t_end_s=float(times[-1]),
        dt_mean_s=dt_mean,
        p_rms_final_pa=float(p_rms_all[-1]),
        p_peak_final_pa=float(p_peak_all[-1]),
        p_rms_mean_pa=float(np.mean(p_rms_all)),
        p_peak_mean_pa=float(np.mean(p_peak_all)),
        ke_proxy_mean=float(np.mean(ke_all)),
        lowk_frac_mean=float(np.mean(lowk_all)),
        station_std_pa=float(np.std(station_arr)),
        station_peak_pa=float(np.max(station_arr)),
        station_trough_pa=float(np.min(station_arr)),
        station_ptp_pa=float(np.ptp(station_arr)),
    )
    return metrics, times, station_arr


def write_summary_csv(path: str, rows: list[RunMetrics]) -> None:
    fields = [
        "sigma_x_m", "sigma_z_m", "n_snapshots", "t_start_s", "t_end_s", "dt_mean_s",
        "p_rms_final_pa", "p_peak_final_pa", "p_rms_mean_pa", "p_peak_mean_pa",
        "ke_proxy_mean", "lowk_frac_mean", "station_std_pa", "station_peak_pa",
        "station_trough_pa", "station_ptp_pa",
    ]
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row.__dict__)


def write_station_timeseries_csv(path: str, run_series: list[tuple[float, float, np.ndarray, np.ndarray]]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["sigma_x_m", "sigma_z_m", "time_s", "p_pert_station_pa"])
        for sigma_x_m, sigma_z_m, times, values in run_series:
            for t, v in zip(times, values):
                writer.writerow([sigma_x_m, sigma_z_m, float(t), float(v)])


def plot_heatmaps(rows: list[RunMetrics], outdir: str) -> None:
    sigma_x_vals = sorted({r.sigma_x_m for r in rows})
    sigma_z_vals = sorted({r.sigma_z_m for r in rows})

    metric_specs = [
        ("p_rms_final_pa", "Final P' RMS (Pa)"),
        ("p_peak_final_pa", "Final max |P'| (Pa)"),
        ("lowk_frac_mean", "Mean large-scale fraction"),
        ("station_ptp_pa", "Station peak-to-peak (Pa)"),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    axes = axes.ravel()

    for ax, (field, title) in zip(axes, metric_specs):
        grid = np.full((len(sigma_z_vals), len(sigma_x_vals)), np.nan, dtype=np.float64)
        for r in rows:
            ix = sigma_x_vals.index(r.sigma_x_m)
            iz = sigma_z_vals.index(r.sigma_z_m)
            grid[iz, ix] = getattr(r, field)

        im = ax.imshow(grid, origin="lower", aspect="auto", cmap="viridis")
        ax.set_title(title)
        ax.set_xlabel("sigma_x (m)")
        ax.set_ylabel("sigma_z (m)")
        ax.set_xticks(range(len(sigma_x_vals)), [f"{v:.0f}" for v in sigma_x_vals])
        ax.set_yticks(range(len(sigma_z_vals)), [f"{v:.0f}" for v in sigma_z_vals])

        for iz in range(len(sigma_z_vals)):
            for ix in range(len(sigma_x_vals)):
                val = grid[iz, ix]
                if np.isfinite(val):
                    ax.text(ix, iz, f"{val:.2g}", ha="center", va="center", color="white", fontsize=8)
        fig.colorbar(im, ax=ax, shrink=0.9)

    plt.tight_layout()
    outpath = os.path.join(outdir, "sweep_heatmaps.png")
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_station_overlay(series: list[tuple[float, float, np.ndarray, np.ndarray]], outdir: str) -> None:
    fig, ax = plt.subplots(figsize=(10, 6))
    for sigma_x_m, sigma_z_m, times, values in series:
        label = f"sx={sigma_x_m:.0f}, sz={sigma_z_m:.0f}"
        ax.plot(times, values, linewidth=1.5, label=label)

    ax.set_xlabel("time (s)")
    ax.set_ylabel("P' at station (Pa)")
    ax.set_title("Station pressure perturbation across sweep")
    ax.grid(True, alpha=0.25)
    ax.legend(ncol=2, fontsize=8)

    outpath = os.path.join(outdir, "station_timeseries_overlay.png")
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    run_dirs = sorted(glob.glob(os.path.join(args.sweep_dir, "run_sx*_sz*")))
    if not run_dirs:
        raise FileNotFoundError(f"No run folders found in {args.sweep_dir}")

    rows: list[RunMetrics] = []
    station_series: list[tuple[float, float, np.ndarray, np.ndarray]] = []

    print(f"Found {len(run_dirs)} run directories")
    for run_dir in run_dirs:
        sigmas = parse_sigma_from_run_dir(run_dir)
        if sigmas is None:
            print(f"Skipping unrecognized folder name: {run_dir}")
            continue

        sigma_x_m, sigma_z_m = sigmas
        metrics, times, series = analyze_single_run(
            run_dir=run_dir,
            sigma_x_m=sigma_x_m,
            sigma_z_m=sigma_z_m,
            lambda_threshold_m=args.large_scale_lambda_km * 1000.0,
            station_x_km=args.station_x_km,
            station_z_km=args.station_z_km,
        )
        rows.append(metrics)
        station_series.append((sigma_x_m, sigma_z_m, times, series))
        print(f"Analyzed sx={sigma_x_m:.0f}, sz={sigma_z_m:.0f} | snapshots={metrics.n_snapshots}")

    rows.sort(key=lambda r: (r.sigma_x_m, r.sigma_z_m))

    summary_csv = os.path.join(args.outdir, "mars_pulse_sweep_summary.csv")
    station_csv = os.path.join(args.outdir, "mars_pulse_sweep_station_timeseries.csv")

    write_summary_csv(summary_csv, rows)
    write_station_timeseries_csv(station_csv, station_series)
    plot_heatmaps(rows, args.outdir)
    plot_station_overlay(station_series, args.outdir)

    rank = sorted(rows, key=lambda r: (r.lowk_frac_mean, r.station_ptp_pa, r.p_rms_final_pa), reverse=True)

    print("\nTop 3 runs by (large-scale fraction, station PTP, final RMS):")
    for i, r in enumerate(rank[:3], start=1):
        print(
            f"{i}. sx={r.sigma_x_m:.0f}, sz={r.sigma_z_m:.0f} | "
            f"lowk={r.lowk_frac_mean:.3f}, station_ptp={r.station_ptp_pa:.3e} Pa, "
            f"final_rms={r.p_rms_final_pa:.3e} Pa"
        )

    print(f"\nWrote summary: {summary_csv}")
    print(f"Wrote station data: {station_csv}")
    print(f"Wrote plots in: {args.outdir}")


if __name__ == "__main__":
    main()
