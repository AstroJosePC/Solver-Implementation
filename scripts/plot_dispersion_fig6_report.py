#!/usr/bin/env python3
"""
plot_dispersion_fig6_report.py
Generates Figure 6 panels for the validation report.
Panel A (Left): B&W Wave travel diagram (reproduced from Julia logic)
Panel B (Right): Per-time normalized wavelet spectrum
"""

import argparse
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from wavelet_tc98 import wavelet

# Set global plotting style for report consistency
plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'figure.titlesize': 16
})

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--indir", default="output/dispersion_test_fig6_strict")
    p.add_argument("--outdir", default="output/dispersion_analysis")
    return p.parse_args()

def load_data(indir):
    files = sorted(glob.glob(os.path.join(indir, "*.nc")))
    if not files:
        raise FileNotFoundError(f"No snapshots in {indir}")
    
    with Dataset(files[0], "r") as ds:
        x = ds.variables["x"][:].astype(np.float64)
        z = ds.variables["z"][:].astype(np.float64)
    
    ix_mid = np.argmin(np.abs(x - 0.5 * (x.min() + x.max())))
    
    times = []
    profiles = []
    
    for f in files:
        with Dataset(f, "r") as ds:
            # NetCDF variables are loaded as (nz, nx) in Python for these files
            w = ds.variables["w"][:, :].astype(np.float64) 
            rho = ds.variables["rho"][:, :].astype(np.float64)
            t = float(ds.getncattr("time_s"))
        
        rho_surface = np.mean(rho[0, :])
        # w scaled by sqrt(rho/rho0)
        # Vertical profile at mid-x (nz,)
        w_scaled = w[:, ix_mid] * np.sqrt(rho[:, ix_mid] / rho_surface)
        
        times.append(t)
        profiles.append(w_scaled)
    
    times = np.array(times)
    order = np.argsort(times)
    times = times[order]
    
    # Filter for t >= 1800.0 (30 minutes) to avoid messy start-up phase
    tsel = times >= 1800.0
    times = times[tsel]
    
    # profiles is list of 1D arrays (length nz), nt of them.
    # Convert to 2D (nz, nt)
    profiles_arr = np.zeros((len(z), len(times)))
    # We only take the profiles where tsel is True
    valid_indices = np.where(tsel)[0]
    for i, orig_idx in enumerate(valid_indices):
        idx = order[orig_idx]
        profiles_arr[:, i] = profiles[idx]
    
    return times, profiles_arr, z

def plot_panel_a(times, profiles, z, outpath):
    """B&W Wave Travel Diagram."""
    T_min = times / 60.0
    Z_km = z / 1000.0
    
    # Scale for visualization like Julia script (1e4 * w_scaled)
    data = 1e4 * profiles
    cabs = np.max(np.abs(data))
    cmax = np.ceil(cabs) if cabs > 0 else 1.0
    
    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.pcolormesh(T_min, Z_km, data, cmap="gray", 
                       vmin=-cmax, vmax=cmax, shading='auto')
    
    ax.set_xlabel("Time (min)")
    ax.set_ylabel("Altitude z (km)")
    cb = fig.colorbar(im, ax=ax)
    cb.set_label(r"m/s $\times 10^{-4}$")
    
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved Panel A: {outpath}")

def plot_panel_b(times, profiles, z, outpath):
    """Per-time normalized wavelet spectrum."""
    nz, nt = profiles.shape
    dz = z[1] - z[0]
    
    # Wavelet params matching TC98 script
    dj = 0.05
    s0 = 2 * dz
    J1 = int(np.fix(np.log2(nz * dz / s0) / dj))
    
    # Compute power panel
    power_panel = np.zeros((J1 + 1, nt))
    
    # Get period array from a dummy call
    _, period, _, _ = wavelet(profiles[:, 0], dz, pad=1, dj=dj, s0=s0, J1=J1)
    
    for i in range(nt):
        y = profiles[:, i]
        y = y - np.mean(y)
        wave, _, _, _ = wavelet(y, dz, pad=1, dj=dj, s0=s0, J1=J1)
        # Power: |W|^2
        Wp = np.abs(wave[:, :nz])**2
        # Max pool over altitude (like in Fig 6b)
        power_panel[:, i] = np.max(Wp, axis=1)
    
    # Per-time normalization
    for i in range(nt):
        pmax = np.max(power_panel[:, i])
        if pmax > 0:
            power_panel[:, i] /= pmax
            
    T_min = times / 60.0
    Lz_km = period / 1000.0
    
    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.pcolormesh(T_min, Lz_km, power_panel, cmap="jet", 
                       vmin=0, vmax=1, shading='auto')
    
    # Plot dominant wavelength line
    dom_idx = np.argmax(power_panel, axis=0)
    ax.plot(T_min, Lz_km[dom_idx], 'w-', lw=1.5)
    
    ax.set_xlabel("t (min)")
    ax.set_ylabel(r"$\lambda_z$ (km)")
    ax.set_ylim(5, 160)
    
    cb = fig.colorbar(im, ax=ax)
    cb.set_label("Normalized power")
    
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved Panel B: {outpath}")

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    
    times, profiles, z = load_data(args.indir)
    
    path_a = os.path.join(args.outdir, "fig6a_report.png")
    path_b = os.path.join(args.outdir, "fig6b_report.png")
    
    plot_panel_a(times, profiles, z, path_a)
    plot_panel_b(times, profiles, z, path_b)

if __name__ == "__main__":
    main()
