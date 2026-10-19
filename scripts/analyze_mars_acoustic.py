#!/usr/bin/env python3
"""
analyze_mars_acoustic.py
Generates the figures for the Acoustic Wave Propagation on Mars benchmark:
- Fig 12: Pressure perturbation snapshot at t=95s
- Fig 13: Time series of pressure at S1-S4
- Fig 14: Wavelet spectrogram for S3
- Fig 15: Wave speed and attenuation fits
"""

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy.optimize import curve_fit
from wavelet_tc98 import wavelet

plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
})

def build_file_time_index(indir):
    files = sorted(glob.glob(os.path.join(indir, "*.nc")))
    if not files:
        raise ValueError(f"No files found in {indir}")

    times = np.empty(len(files), dtype=np.float64)
    for i, f in enumerate(files):
        with Dataset(f, "r") as ds:
            times[i] = float(ds.getncattr("time_s"))

    return files, times

def get_snapshot_by_time(indir, target_time_s, tolerance=10.0, files=None, times=None):
    if files is None or times is None:
        files, times = build_file_time_index(indir)

    idx = int(np.argmin(np.abs(times - target_time_s)))
    closest_file = files[idx]
    min_diff = abs(float(times[idx]) - target_time_s)

    if min_diff <= tolerance:
        print(f"  Picked {os.path.basename(closest_file)} at t = {times[idx]:.2f} s")
        return closest_file
    else:
        raise ValueError(f"Could not find snapshot near {target_time_s} s (closest: {min_diff} s diff)")

def plot_fig12(indir, outdir, files=None, times=None):
    fname = get_snapshot_by_time(indir, 95.0, files=files, times=times)
    print(f"Fig 12: Loading {os.path.basename(fname)} (t=95s)")

    with Dataset(fname, "r") as ds:
        x = ds.variables["x"][:].astype(np.float64) / 1000.0  # km
        z = ds.variables["z"][:].astype(np.float64) / 1000.0  # km
        P = ds.variables["P"][:, :].astype(np.float64)

    P_bar = np.mean(P, axis=1, keepdims=True)
    P_pert = P - P_bar

    X, Z = np.meshgrid(x - 100.0, z) # shift origin to center for plotting

    fig, ax = plt.subplots(figsize=(8, 3.5))
    levels = np.linspace(-20, 5, 101)
    im = ax.contourf(X, Z, P_pert, levels=levels, cmap="gray", extend="both")

    # Mark stations with crosses as in Fig 12
    stations_x = [20, 40, 60, 80] # km from center
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:purple']
    for i, sx in enumerate(stations_x):
        ax.plot(sx, 20.0, 'x', color=colors[i], markersize=10, markeredgewidth=3)
        ax.text(sx - 3, 22.0, f"S{i+1}", fontweight='bold', fontsize=12)

    # Annotation inside as in paper
    ax.text(-80, 32, "T = 95 s", fontsize=14, fontweight='bold')

    ax.set_xlabel("x (km)")
    ax.set_ylabel("z (km)")
    ax.set_xlim(-100, 100)
    ax.set_ylim(0, 40)
    ax.set_xticks([-50, 0, 50])

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Pa", loc='top', rotation=0, labelpad=-20)
    cbar.set_ticks([-20, -15, -10, -5, 0, 5])

    outpath = os.path.join(outdir, "mars_fig12.png")
    plt.savefig(outpath, dpi=200, bbox_inches='tight')
    print(f"Saved: {outpath}")
    plt.close()

def extract_station_and_front_data(files, times_index):
    if len(files) == 0:
        raise ValueError("No files provided")

    # Read grid
    with Dataset(files[0], "r") as ds:
        x = ds.variables["x"][:].astype(np.float64) / 1000.0  # km
        z = ds.variables["z"][:].astype(np.float64) / 1000.0  # km

    # Station locations (x=120, 140, 160, 180 km; z=20 km)
    ix_s1 = np.argmin(np.abs(x - 120.0))
    ix_s2 = np.argmin(np.abs(x - 140.0))
    ix_s3 = np.argmin(np.abs(x - 160.0))
    ix_s4 = np.argmin(np.abs(x - 180.0))
    iz_20 = np.argmin(np.abs(z - 20.0))

    # We only care about x > source (x > 100)
    ix_start = np.argmin(np.abs(x - 100.0))
    x_right = x[ix_start:] - 100.0  # Distance from source in km

    max_amps = np.zeros(len(x_right), dtype=np.float64)
    time_of_max = np.zeros(len(x_right), dtype=np.float64)

    times_out = []
    p_s1 = []
    p_s2 = []
    p_s3 = []
    p_s4 = []

    print("Extracting time series and arrival metrics (this may take a minute)...")
    for i, f in enumerate(files):
        t = float(times_index[i])
        with Dataset(f, "r") as ds:
            # Only read the specific slice to save memory
            P_row = ds.variables["P"][iz_20, :]

        P_bar = np.mean(P_row) # approximate background
        P_pert = P_row - P_bar

        times_out.append(t)
        p_s1.append(P_pert[ix_s1])
        p_s2.append(P_pert[ix_s2])
        p_s3.append(P_pert[ix_s3])
        p_s4.append(P_pert[ix_s4])

        # Update max amplitudes and time of max for Fig 15
        P_pert_abs_right = np.abs(P_pert[ix_start:])
        mask = P_pert_abs_right > max_amps
        max_amps[mask] = P_pert_abs_right[mask]
        time_of_max[mask] = t

    return (
        np.array(times_out),
        np.array(p_s1),
        np.array(p_s2),
        np.array(p_s3),
        np.array(p_s4),
        x_right,
        time_of_max,
        max_amps,
    )

def plot_fig13(times, p_s1, p_s2, p_s3, p_s4, outdir):
    fig, axes = plt.subplots(4, 1, figsize=(10, 8), sharex=True)

    series = [p_s1, p_s2, p_s3, p_s4]
    labels = ["S1 (20 km from source)", "S2 (40 km from source)",
              "S3 (60 km from source)", "S4 (80 km from source)"]

    for i in range(4):
        axes[i].plot(times, series[i], 'tab:blue', linewidth=1.2)
        axes[i].set_ylim(-1, 1)
        axes[i].set_xlim(0, 2000)
        axes[i].set_ylabel("Pa", fontsize=12)
        # Station label inside upper right
        axes[i].text(1950, 0.5, labels[i], horizontalalignment='right',
                     fontsize=11, fontweight='bold')

    axes[3].set_xlabel("t (s)", fontsize=14)

    plt.subplots_adjust(hspace=0.1)

    outpath = os.path.join(outdir, "mars_fig13.png")
    plt.savefig(outpath, dpi=200, bbox_inches='tight')
    print(f"Saved: {outpath}")
    plt.close()

def plot_fig14(times, p_s3, outdir):
    dt = times[1] - times[0]
    # Wavelet params
    dj = 0.05
    s0 = 2 * dt
    J1 = int(np.fix(np.log2(len(times) * dt / s0) / dj))

    y = p_s3 - np.mean(p_s3)
    wave, period, scale, coi = wavelet(y, dt, pad=1, dj=dj, s0=s0, J1=J1, mother="MORLET", param=6.0)

    power = np.abs(wave)**2
    # Normalize by max
    power /= np.max(power)

    freq_mHz = (1.0 / period) * 1000.0
    times_min = times / 60.0

    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.pcolormesh(times_min, freq_mHz, power, cmap="jet", shading='auto')

    # Acoustic cutoff approx 1.8 mHz
    ax.axhline(1.8, color='limegreen', linestyle='--', lw=3)

    ax.set_ylim(0, 30) # paper shows up to ~30 mHz
    ax.set_xlim(0, 15)
    ax.set_xlabel("time (min)", fontsize=14)
    ax.set_ylabel("frequency (mHz)", fontsize=14)

    fig.colorbar(im, ax=ax)

    outpath = os.path.join(outdir, "mars_fig14.png")
    plt.savefig(outpath, dpi=200, bbox_inches='tight')
    print(f"Saved: {outpath}")
    plt.close()

def plot_fig15(x_right, time_of_max, max_amps, outdir):
    # Fit speed (t = m*x + c)
    valid_idx = np.where(x_right > 1.0)[0]
    xr = x_right[valid_idx]
    tm = time_of_max[valid_idx]
    ma = max_amps[valid_idx]

    def linear_fit(x, m, c):
        return m*x + c

    popt_v, _ = curve_fit(linear_fit, xr, tm)
    r2_v = 1 - np.sum((tm - linear_fit(xr, *popt_v))**2) / np.sum((tm - np.mean(tm))**2)

    # Fit attenuation (A = A0 / sqrt(r) + C)
    def geom_atten(x, A0, C):
        return A0 / np.sqrt(x) + C

    popt_a, _ = curve_fit(geom_atten, xr, ma, p0=[10.0, 0.5])
    r2_a = 1 - np.sum((ma - geom_atten(xr, *popt_a))**2) / np.sum((ma - np.mean(ma))**2)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Panel 1: Speed
    ax1.plot(xr, tm, 'k.', markersize=6, label="Simulation data")
    ax1.plot(xr, linear_fit(xr, *popt_v), 'g-', linewidth=2, label="Curve fit")
    ax1.set_xlabel("distance from source (km)")
    ax1.set_ylabel("time of arrival (s)")
    ax1.set_xlim(0, 100)
    ax1.set_ylim(0, 600)
    ax1.text(5, 450, f"y = {popt_v[0]:.3f}x + {popt_v[1]:.1f}\n$R^2 = {r2_v:.4f}$", fontsize=12)
    ax1.legend(loc='lower right')

    # Panel 2: Attenuation
    ax2.plot(xr, ma, 'k.', markersize=6, label="Simulation data")
    ax2.plot(xr, geom_atten(xr, *popt_a), 'b-', linewidth=2, label="Curve fit")
    ax2.set_xlabel("distance from source (km)")
    ax2.set_ylabel("Amplitude (Pa)")
    ax2.set_xlim(0, 100)
    ax2.set_ylim(0, 4)
    ax2.text(40, 2.5, f"y = {popt_a[0]:.2f}/$\sqrt{{x}}$ {popt_a[1]:+.2f}\n$R^2 = {r2_a:.4f}$", fontsize=12)
    ax2.legend(loc='upper right')

    outpath = os.path.join(outdir, "mars_fig15.png")
    plt.savefig(outpath, dpi=200, bbox_inches='tight')
    print(f"Saved: {outpath}")

def main():
    indir = "output/mars_acoustic_test"
    outdir = "output/mars_acoustic_analysis"
    os.makedirs(outdir, exist_ok=True)

    files, times_index = build_file_time_index(indir)

    plot_fig12(indir, outdir, files=files, times=times_index)

    times, p_s1, p_s2, p_s3, p_s4, x_right, time_of_max, max_amps = extract_station_and_front_data(files, times_index)
    plot_fig13(times, p_s1, p_s2, p_s3, p_s4, outdir)
    plot_fig14(times, p_s3, outdir)
    plot_fig15(x_right, time_of_max, max_amps, outdir)

if __name__ == "__main__":
    main()
