#!/usr/bin/env python3
"""
Tier 1 Analysis Script for MAGNUSP Mars Simulations
Analyzes outputs from Case V1 (no dissipation), Case V2 (with dissipation, long wave),
and Case V3 (with dissipation, short wave).
"""

import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

# Add current script directory to path to import wavelet_tc98
sys.path.append(os.path.dirname(__file__))
from wavelet_tc98 import wavelet

def load_midx_profiles(indir: str):
    files = sorted(glob.glob(os.path.join(indir, "snapshot_*.nc")))
    if not files:
        raise FileNotFoundError(f"No snapshot_*.nc found in {indir}")

    times = []
    w_scaled_list = []
    w_raw_list = []
    
    x = z = None
    ix_mid = None

    for fp in files:
        with Dataset(fp, "r") as ds:
            if x is None:
                x = np.array(ds.variables["x"][:], dtype=float)
                z = np.array(ds.variables["z"][:], dtype=float)
                ix_mid = int(np.argmin(np.abs(x - 0.5 * (x.min() + x.max()))))

            w = np.array(ds.variables["w"][:, :], dtype=float).T
            rho = np.array(ds.variables["rho"][:, :], dtype=float).T
            t = float(ds.getncattr("time_s"))

            rho_surface = np.mean(rho[:, 0])
            weight = np.sqrt(np.clip(rho[ix_mid, :] / rho_surface, 0.0, None))
            w_mid = w[ix_mid, :]

            times.append(t)
            w_raw_list.append(w_mid)
            w_scaled_list.append(w_mid * weight)

    times = np.array(times, dtype=float)
    order = np.argsort(times)
    
    return times[order], z, np.array(w_raw_list, dtype=float)[order, :], np.array(w_scaled_list, dtype=float)[order, :]

def get_snapshot_2d(indir: str, target_time: float):
    files = sorted(glob.glob(os.path.join(indir, "snapshot_*.nc")))
    best_file = None
    best_diff = float('inf')
    
    for fp in files:
        with Dataset(fp, "r") as ds:
            t = float(ds.getncattr("time_s"))
            if abs(t - target_time) < best_diff:
                best_diff = abs(t - target_time)
                best_file = fp
                
    if not best_file:
        return None, None, None, None
        
    with Dataset(best_file, "r") as ds:
        x = np.array(ds.variables["x"][:], dtype=float)
        z = np.array(ds.variables["z"][:], dtype=float)
        w = np.array(ds.variables["w"][:, :], dtype=float).T
        rho = np.array(ds.variables["rho"][:, :], dtype=float).T
        t_actual = float(ds.getncattr("time_s"))
        
    rho_surface = np.mean(rho[:, 0])
    w_scaled = w * np.sqrt(np.clip(rho / rho_surface, 0.0, None))
    
    return x, z, w_scaled, t_actual

def plot_2d_snapshot(x, z, w_scaled, t, title, outpath, cmap='RdBu_r', tiles=3):
    fig, ax = plt.subplots(figsize=(6, 8))
    
    # Guarantee w_scaled is (nz, nx)
    if w_scaled.shape[0] != len(z):
        w_scaled = w_scaled.T
        
    dx = x[1] - x[0]
    w_tiled = np.tile(w_scaled, (1, tiles))
    x_tiled = np.linspace(0, (x[-1] + dx) * tiles, len(x) * tiles, endpoint=False)
    
    x_km = x_tiled / 1000.0
    z_km = z / 1000.0
    
    mask_z = z_km <= 170.0
    z_km = z_km[mask_z]
    
    w_tiled = w_tiled[mask_z, :]
    
    vm = np.max(np.abs(w_tiled)) * 0.8
    if vm == 0: vm = 1.0
    
    im = ax.pcolormesh(x_km, z_km, w_tiled, cmap=cmap, shading='auto', vmin=-vm, vmax=vm)
    
    ax.set_xlabel("Horizontal Distance (km)")
    ax.set_ylabel("Altitude (km)")
    ax.set_title(f"{title} (t = {t:.0f} s)\n[{tiles} horizontal periods shown]")
    cbar = fig.colorbar(im, ax=ax, orientation='horizontal', pad=0.1)
    cbar.set_label("Density-Scaled w' (m/s)")
    
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close(fig)

def plot_travel_diagram(times, z, w_scaled, title, outpath):
    fig, (ax_main, ax_force) = plt.subplots(2, 1, figsize=(8, 8), gridspec_kw={'height_ratios': [4, 1]}, sharex=True)
    times_min = times / 60.0
    z_km = z / 1000.0
    
    mask_z = z_km <= 170.0
    z_km = z_km[mask_z]
    w_scaled = w_scaled[:, mask_z]
    
    vm = np.max(np.abs(w_scaled)) * 0.8
    if vm == 0: vm = 1.0
    
    im = ax_main.imshow(
        w_scaled.T,
        origin='lower',
        aspect='auto',
        extent=[times_min.min(), times_min.max(), z_km.min(), z_km.max()],
        cmap='RdBu_r',
        vmin=-vm, vmax=vm
    )
    
    ax_main.set_ylabel("Altitude (km)")
    ax_main.set_title(title)
    cbar = fig.colorbar(im, ax=ax_main)
    cbar.set_label("Density-Scaled w' (m/s)")
    
    # Bottom panel: Forcing
    t0_s = 1200.0
    sigma_s = 600.0
    omega = 0.005
    w_tilde = 0.001
    
    # Reconstruct forcing over time
    t_seconds = times_min * 60.0
    w_force = w_tilde * np.cos(omega * (t_seconds - t0_s)) * np.exp(-((t_seconds - t0_s)**2) / (2 * sigma_s**2))
    
    ax_force.plot(times_min, w_force, color='black')
    ax_force.set_xlabel("Time (min)")
    ax_force.set_ylabel("Forcing w' (m/s)")
    ax_force.grid(True, linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close(fig)

def plot_damping_comparison(z, profiles, outpath):
    z_km = z / 1000.0
    mask_z = z_km <= 170.0
    z_km = z_km[mask_z]
    
    fig, ax = plt.subplots(figsize=(6, 8))
    
    colors = ['b', 'r', 'g', 'orange']
    for idx, (label, w_raw) in enumerate(profiles.items()):
        env = np.max(np.abs(w_raw), axis=0)
        # Normalize relative to lowest valid altitude to easily compare decay rates
        env = env / env[0]
        ax.plot(env[mask_z], z_km, label=label, color=colors[idx % len(colors)])
    
    ax.set_ylabel("Altitude (km)")
    ax.set_xlabel("Relative Raw w' Amplitude Growth")
    ax.set_title("Amplitude Growth vs Altitude")
    ax.legend()
    ax.grid(True, linestyle='--', alpha=0.6)
    
    # Use log scale for x-axis since the growth is exponential
    ax.set_xscale('log')
    
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close(fig)

def plot_phase_analysis(indir: str, target_z_km: float, outpath: str):
    files = sorted(glob.glob(os.path.join(indir, "snapshot_*.nc")))
    if not files:
        raise FileNotFoundError(f"No snapshot_*.nc found in {indir}")

    times = []
    data_t = {'rho': [], 'P': [], 'w': [], 'u': []}
    
    x = z = None
    
    for fp in files:
        with Dataset(fp, "r") as ds:
            if x is None:
                x = np.array(ds.variables["x"][:], dtype=float)
                z = np.array(ds.variables["z"][:], dtype=float)
            
            data_t['rho'].append(np.array(ds.variables["rho"][:, :], dtype=float))
            data_t['P'].append(np.array(ds.variables["P"][:, :], dtype=float))
            data_t['w'].append(np.array(ds.variables["w"][:, :], dtype=float))
            data_t['u'].append(np.array(ds.variables["u"][:, :], dtype=float))
            times.append(float(ds.getncattr("time_s")))

    times = np.array(times, dtype=float)
    order = np.argsort(times)
    times = times[order]
    
    for k in data_t:
        data_t[k] = np.array(data_t[k], dtype=float)[order]

    # Constants
    R_gas = 188.9  # Mars CO2
    
    # Target indices
    iz = int(np.argmin(np.abs(z / 1000.0 - target_z_km)))
    
    # Midpoint and +140km point
    mid_x = 0.5 * (x.min() + x.max())
    ix1 = int(np.argmin(np.abs(x - mid_x)))
    
    # Handle periodic wrap-around for +140km
    Lx = x.max() - x.min() + (x[1] - x[0])
    target_x2 = mid_x + 140e3
    target_x2_wrapped = (target_x2 - x.min()) % Lx + x.min()
    ix2 = int(np.argmin(np.abs(x - target_x2_wrapped)))

    # Extract time series at the two points
    ts1 = {k: data_t[k][:, ix1, iz] for k in data_t}
    ts2 = {k: data_t[k][:, ix2, iz] for k in data_t}
    
    # Compute T
    T1 = ts1['P'] / (ts1['rho'] * R_gas)
    T2 = ts2['P'] / (ts2['rho'] * R_gas)
    
    # Compute base states and perturbations
    rho0_1 = np.mean(ts1['rho'])
    T0_1 = np.mean(T1)
    
    rho0_2 = np.mean(ts2['rho'])
    T0_2 = np.mean(T2)
    
    drho_rho_1 = (ts1['rho'] - rho0_1) / rho0_1
    dT_T_1 = (T1 - T0_1) / T0_1
    
    drho_rho_2 = (ts2['rho'] - rho0_2) / rho0_2
    dT_T_2 = (T2 - T0_2) / T0_2

    times_min = times / 60.0

    # 4 rows, 2 columns
    fig, axes = plt.subplots(4, 2, figsize=(10, 12), sharex=True)
    
    # Row 0: Density
    axes[0, 0].plot(times_min, drho_rho_1, color='blue')
    axes[0, 0].set_ylabel(r"$\delta \rho / \rho_0$")
    axes[0, 0].set_title(f"x = {x[ix1]/1000.0:.1f} km (Source Center)")
    
    axes[0, 1].plot(times_min, drho_rho_2, color='blue')
    axes[0, 1].set_title(f"x = {x[ix2]/1000.0:.1f} km (+140 km offset)")
    
    # Row 1: Temperature
    axes[1, 0].plot(times_min, dT_T_1, color='red')
    axes[1, 0].set_ylabel(r"$\delta T / T_0$")
    axes[1, 1].plot(times_min, dT_T_2, color='red')
    
    # Row 2: w
    axes[2, 0].plot(times_min, ts1['w'], color='green')
    axes[2, 0].set_ylabel(r"$w$ (m/s)")
    axes[2, 1].plot(times_min, ts2['w'], color='green')
    
    # Row 3: u
    axes[3, 0].plot(times_min, ts1['u'], color='purple')
    axes[3, 0].set_ylabel(r"$u$ (m/s)")
    axes[3, 0].set_xlabel("Time (min)")
    axes[3, 1].plot(times_min, ts2['u'], color='purple')
    axes[3, 1].set_xlabel("Time (min)")
    
    # Share y-axes across columns for each row to allow direct comparison
    for r in range(4):
        y_min = min(axes[r, 0].get_ylim()[0], axes[r, 1].get_ylim()[0])
        y_max = max(axes[r, 0].get_ylim()[1], axes[r, 1].get_ylim()[1])
        axes[r, 0].set_ylim(y_min, y_max)
        axes[r, 1].set_ylim(y_min, y_max)
        axes[r, 0].grid(True, linestyle='--', alpha=0.5)
        axes[r, 1].grid(True, linestyle='--', alpha=0.5)

    fig.suptitle(f"Phase Analysis at z = {target_z_km} km", fontsize=14)
    plt.tight_layout()
    plt.subplots_adjust(top=0.92) # Make room for suptitle
    plt.savefig(outpath, dpi=150)
    plt.close(fig)

def main():
    dirs = {
        "V1 (No Dissipation, λx=150km)": "output/tier1_v1",
        "V2 (Dissipation, λx=150km)": "output/tier1_v2",
        "V3 (Dissipation, λx=50km)": "output/tier1_v3"
    }
    
    for path in dirs.values():
        if not os.path.exists(path):
            print(f"Error: Output directory {path} not found.")
            return

    data = {}
    for label, d in dirs.items():
        print(f"Loading {label}...")
        t, z, w_raw, w_scaled = load_midx_profiles(d)
        data[label] = (t, z, w_raw, w_scaled)

    # Plot travel diagrams with envelope overlay
    print("Plotting travel diagrams...")
    for label, d in dirs.items():
        t, z, w_raw, w_scaled = data[label]
        slug = d.split('/')[-1]
        plot_travel_diagram(t, z, w_scaled, f"Travel Diagram: {label}", f"{d}/{slug}_travel.png")
    
    # Plot 2D snapshots for V3 at 8000s
    print("Extracting 2D snapshots for V3...")
    x3, z3, w3_2d, t3_actual = get_snapshot_2d("output/tier1_v3", 8000.0)
    if x3 is not None:
        plot_2d_snapshot(x3, z3, w3_2d, t3_actual, "V3 2D Snapshot", "output/tier1_v3/tier1_v3_2d_rdbu.png", cmap="RdBu_r")
        plot_2d_snapshot(x3, z3, w3_2d, t3_actual, "V3 2D Snapshot", "output/tier1_v3/tier1_v3_2d_gray.png", cmap="gray")
    
    # Plot Damping Comparison
    print("Plotting damping comparison...")
    profiles_to_compare = {label: d[2] for label, d in data.items()}
    plot_damping_comparison(data[list(data.keys())[0]][1], profiles_to_compare, "output/tier1_v3/tier1_damping_comparison.png")
    
    # Phase Analysis
    print("Plotting phase analysis for V2...")
    plot_phase_analysis("output/tier1_v2", target_z_km=100.0, outpath="output/tier1_v2/tier1_v2_phase_analysis.png")
    
    print("Analysis complete. Check output/tier1_v*/ directories for updated plots.")

if __name__ == "__main__":
    main()
