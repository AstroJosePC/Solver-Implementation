#!/usr/bin/env python3
"""
wavelet_fig6b_pywt.py
Reproduce Figure 6b using PyWavelets module instead of custom wavelet_tc98.
Usage:
python wavelet_fig6b_pywt.py --indir=output/dispersion_test_fig6_strict
"""
import argparse
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import pywt  # PyWavelets module

def parse_args():
    p = argparse.ArgumentParser(description="Fig 6b wavelet analysis (PyWavelets)")
    p.add_argument("--indir", default="output/dispersion_test_fig6_strict",
                   help="Directory containing .nc snapshots")
    p.add_argument("--outfile", default="wavelet_fig6b_pywt.png",
                   help="Output figure path")
    p.add_argument("--t-start", type=float, default=1800.0,
                   help="Start time in seconds (default: 1800 = 30 min)")
    p.add_argument("--mother", default="morl",
                   help="Wavelet mother (morl, paul, dgauN for DOG)")
    p.add_argument("--scales", type=int, default=128,
                   help="Number of scales for CWT")
    p.add_argument("--pool", default="max", choices=["max", "mean"],
                   help="Pooling over altitude")
    p.add_argument("--lambda-min", type=float, default=5.0,
                   help="Y-axis min (km)")
    p.add_argument("--lambda-max", type=float, default=160.0,
                   help="Y-axis max (km)")
    return p.parse_args()


def load_snapshots(indir):
    """Load all .nc snapshots, return sorted times, profiles, x, z."""
    files = sorted(glob.glob(os.path.join(indir, "*.nc")))
    if not files:
        raise FileNotFoundError(f"No .nc files found in {indir}")

    # Read grid from first file
    with Dataset(files[0], "r") as ds:
        x = ds.variables["x"][:].astype(np.float64)
        z = ds.variables["z"][:].astype(np.float64)
        ix_mid = np.argmin(np.abs(x - 0.5 * (x.min() + x.max())))
        nz = len(z)

    times = []
    profiles = []
    for f in files:
        with Dataset(f, "r") as ds:
            w = ds.variables["w"][:, :].T.astype(np.float64)
            rho = ds.variables["rho"][:, :].T.astype(np.float64)
            t = float(ds.getncattr("time_s"))
            rho_surface = np.mean(rho[:, 0])
            wscaled = w[ix_mid, :] * np.sqrt(rho[ix_mid, :] / rho_surface)
            times.append(t)
            profiles.append(wscaled)

    # Sort by time
    order = np.argsort(times)
    times = np.array(times)[order]
    profiles = [profiles[i] for i in order]

    print(f"Loaded {len(files)} snapshots from {indir}")
    print(f"  z range: {z[0]/1e3:.1f} -- {z[-1]/1e3:.1f} km, nz={nz}")
    print(f"  t range: {times[0]:.0f} -- {times[-1]:.0f} s")

    return times, profiles, x, z, ix_mid


def compute_wavelet_panel(profiles, times, z, tsel,
                          mother="morl", n_scales=128,
                          pool_mode="max"):
    """
    Compute lambda_z vs time power panel using PyWavelets CWT.

    Parameters
    ----------
    profiles : list of 1D arrays, vertical profiles at each time
    times : 1D array, times in seconds
    z : 1D array, altitude grid in meters
    tsel : index array, which time steps to use
    mother : wavelet name (PyWavelets convention)
    n_scales : number of scales for CWT
    pool_mode : 'max' or 'mean' over altitude

    Returns
    -------
    power_panel : (n_scales, n_times) normalized power
    period_km : (n_scales,) vertical wavelengths in km
    scales : (n_scales,) scale values
    """
    nz = len(z)
    dz = z[1] - z[0]  # meters (sampling period)

    # Create logarithmically spaced scales
    # PyWavelets scales are in units of the sampling period
    min_scale = 2  # minimum scale (similar to T&C98 s0=2*dt)
    max_scale = nz / 4  # maximum scale (avoid edge effects)
    scales = np.logspace(np.log10(min_scale), np.log10(max_scale), n_scales)

    # Get wavelet function to compute central frequency for period conversion
    # wavelet = pywt.Wavelet(mother)

    # For Morlet, central frequency is approximately wavelet.wavefun()[2]
    # We'll use the built-in central_frequency if available
    try:
        fc = pywt.central_frequency(mother)
    except:
        # Fallback estimates for common wavelets
        if mother == "morl":
            fc = 0.8125  # Morlet central frequency
        elif mother.startswith("paul"):
            fc = 0.5
        elif mother.startswith("dgau"):
            fc = 0.5
        else:
            fc = 0.5

    # Convert scales to periods (wavelengths in meters)
    # period = scale * sampling_period / central_frequency
    period_m = scales * dz / fc
    period_km = period_m / 1e3

    nt = len(tsel)
    print(f"  Wavelet: {mother}, {n_scales} scales")
    print(f"  Period range: {period_km.min():.1f} -- {period_km.max():.1f} km")

    power_panel = np.zeros((n_scales, nt))

    for col, it in enumerate(tsel):
        y = profiles[it].copy()
        y -= np.mean(y)  # remove mean

        # Compute Continuous Wavelet Transform
        # pywt.cwt returns (coefficients, frequencies)
        coefficients, frequencies = pywt.cwt(y, scales, mother, sampling_period=dz)

        # Power spectrum: |coefficients|^2
        # coefficients shape: (n_scales, n_data_points)
        Wp = np.abs(coefficients) ** 2

        # Pool over altitude
        if pool_mode == "max":
            power_panel[:, col] = np.max(Wp, axis=1)
        else:
            power_panel[:, col] = np.mean(Wp, axis=1)

    return power_panel, period_km, scales, None  # COI not computed in PyWavelets


def normalize_per_time(panel):
    """Normalize each time column by its own maximum."""
    out = panel.copy()
    for col in range(out.shape[1]):
        pmax = out[:, col].max()
        if pmax > 0:
            out[:, col] /= pmax
    return np.clip(out, 0, 1)


def normalize_global(panel):
    """Normalize by global maximum."""
    pmax = panel.max()
    if pmax > 0:
        return np.clip(panel / pmax, 0, 1)
    return panel


def compute_coi_pywt(scales, dz, n_data, mother="morl"):
    """
    Approximate Cone of Influence for PyWavelets.
    This is a simplified version - PyWavelets doesn't provide COI directly.
    """
    wavelet = pywt.Wavelet(mother)
    # Get wavelet support (approximate)
    try:
        support = wavelet.support
    except:
        support = 2.0  # default estimate

    # COI at each scale: scale * support * sampling_period
    coi_scales = scales * support * dz

    # COI at each time point (edge effects)
    coi = np.zeros(n_data)
    for i in range(n_data):
        dist_to_edge = min(i, n_data - 1 - i) * dz
        coi[i] = dist_to_edge / support / dz

    return coi_scales, coi


def main():
    args = parse_args()

    # Load data
    times, profiles, x, z, ix_mid = load_snapshots(args.indir)

    # Select time range
    tsel = np.where(times >= args.t_start)[0]
    if len(tsel) == 0:
        raise ValueError(f"No snapshots with t >= {args.t_start} s")

    print(f"Using {len(tsel)} snapshots with t >= {args.t_start:.0f} s")

    # Compute wavelet panel
    power_panel, period_km, scales, coi = compute_wavelet_panel(
        profiles, times, z, tsel,
        mother=args.mother, n_scales=args.scales,
        pool_mode=args.pool
    )

    t_min = times[tsel] / 60.0  # minutes

    # ---- Figure: 2-panel (per-time norm + global norm) ----
    power_pertime = normalize_per_time(power_panel)
    power_global = normalize_global(power_panel)

    fig, axes = plt.subplots(2, 1, figsize=(10, 10), constrained_layout=True)

    # Panel 1: per-time normalization
    ax = axes[0]
    im = ax.pcolormesh(t_min, period_km, power_pertime,
                       cmap="jet", vmin=0, vmax=1, shading="auto")
    fig.colorbar(im, ax=ax, label="Normalized power")

    # Dominant wavelength
    idx_dom = np.argmax(power_pertime, axis=0)
    ax.plot(t_min, period_km[idx_dom], "w-", lw=1.5)
    ax.set_ylim(args.lambda_min, args.lambda_max)
    ax.set_xlabel("t (min)")
    ax.set_ylabel(r"$\lambda_z$ (km)")
    ax.set_title(f"Per-time normalized wavelet spectrum "
                 f"(PyWavelets, {args.mother})")

    # Panel 2: global normalization
    ax = axes[1]
    im = ax.pcolormesh(t_min, period_km, power_global,
                       cmap="jet", vmin=0, vmax=1, shading="auto")
    fig.colorbar(im, ax=ax, label="Normalized power")

    idx_dom2 = np.argmax(power_global, axis=0)
    ax.plot(t_min, period_km[idx_dom2], "w-", lw=1.5)
    ax.set_ylim(args.lambda_min, args.lambda_max)
    ax.set_xlabel("t (min)")
    ax.set_ylabel(r"$\lambda_z$ (km)")
    ax.set_title(f"Global-norm wavelet spectrum "
                 f"(PyWavelets, {args.mother})")

    outfile = args.outfile
    fig.savefig(outfile, dpi=200, bbox_inches="tight")
    print(f"Saved: {outfile}")
    plt.close(fig)


if __name__ == "__main__":
    main()