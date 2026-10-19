#!/usr/bin/env python3
"""
Tier 2 Analysis Script: Extract DSMC Boundary Conditions
Extracts time series of macroscopic fluid fields at the DSMC handoff altitude (z = 100 km)
from the realistic Mars MAGNUSP simulations, performs harmonic analysis, and outputs
parameters ready for ingestion by HARRAH.
"""

import os
import sys
import glob
import json
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

def load_handoff_series(indir: str, target_z_km: float):
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

    R_gas = 188.9  # Mars CO2
    
    # Target index for altitude
    iz = int(np.argmin(np.abs(z / 1000.0 - target_z_km)))
    
    # We will horizontally average to get the general perturbation amplitude 
    # Or we can track a single point in the center. 
    # For harmonic amplitude extraction, picking a single point x_mid is simpler and preserves phase.
    ix_mid = int(np.argmin(np.abs(x - 0.5 * (x.min() + x.max()))))

    # Extract time series at the single point
    ts = {k: data_t[k][:, ix_mid, iz] for k in data_t}
    
    # Compute T
    ts['T'] = ts['P'] / (ts['rho'] * R_gas)
    
    # Background state
    rho0 = np.mean(ts['rho'])
    T0 = np.mean(ts['T'])
    
    # Fractional perturbations
    ts['drho_frac'] = (ts['rho'] - rho0) / rho0
    ts['dT_frac'] = (ts['T'] - T0) / T0
    
    return times, ts, rho0, T0

def extract_dominant_harmonic(times, signal, target_omega=None):
    """
    Extract amplitude and phase of the dominant frequency.
    """
    dt = times[1] - times[0]
    N = len(times)
    Y = np.fft.rfft(signal - np.mean(signal))
    freq = np.fft.rfftfreq(N, d=dt)
    
    # Convert freq to angular freq (omega)
    omega_grid = 2 * np.pi * freq
    
    if target_omega is not None:
        idx = np.argmin(np.abs(omega_grid - target_omega))
    else:
        idx = np.argmax(np.abs(Y[1:])) + 1
        
    amplitude = 2.0 * np.abs(Y[idx]) / N
    phase = np.angle(Y[idx])
    found_omega = omega_grid[idx]
    
    return amplitude, phase, found_omega

def process_case(case_name, indir, target_z_km, expected_omega):
    print(f"Processing {case_name} from {indir} at z = {target_z_km} km")
    if not os.path.exists(indir):
        print(f"Directory {indir} not found. Skipping.")
        return
        
    times, ts, rho0, T0 = load_handoff_series(indir, target_z_km)
    
    # Use only the second half of the time series to ensure steady state
    half_idx = len(times) // 2
    t_steady = times[half_idx:]
    ts_steady = {k: v[half_idx:] for k, v in ts.items()}
    
    amp_rho, phase_rho, om = extract_dominant_harmonic(t_steady, ts_steady['drho_frac'], expected_omega)
    amp_T, phase_T, _ = extract_dominant_harmonic(t_steady, ts_steady['dT_frac'], expected_omega)
    amp_w, phase_w, _ = extract_dominant_harmonic(t_steady, ts_steady['w'], expected_omega)
    amp_u, phase_u, _ = extract_dominant_harmonic(t_steady, ts_steady['u'], expected_omega)
    
    # Calculate relative phase lag (T relative to rho)
    phase_T_rel = phase_T - phase_rho
    # Wrap to [-pi, pi]
    phase_T_rel = (phase_T_rel + np.pi) % (2 * np.pi) - np.pi
    
    results = {
        "case": case_name,
        "handoff_altitude_km": target_z_km,
        "omega_rad_s": float(om),
        "period_min": float(2 * np.pi / om / 60.0),
        "background_rho": float(rho0),
        "background_T": float(T0),
        "delta_n": float(amp_rho),
        "delta_T": float(amp_T),
        "w_amplitude": float(amp_w),
        "u_amplitude": float(amp_u),
        "phase_rho": float(phase_rho),
        "phase_T": float(phase_T),
        "phase_T_rel_to_rho": float(phase_T_rel)
    }
    
    # Output to JSON
    out_json = os.path.join(indir, f"{case_name}_dsmc_bc.json")
    with open(out_json, 'w') as f:
        json.dump(results, f, indent=4)
    print(f"  Saved DSMC parameters to {out_json}")
    
    # Visualization
    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    t_min = times / 60.0
    
    axes[0].plot(t_min, ts['drho_frac'], label='drho/rho0', color='blue')
    axes[0].plot(t_min, ts['dT_frac'], label='dT/T0', color='red', alpha=0.8)
    axes[0].set_ylabel("Fractional Perturbation")
    axes[0].legend()
    axes[0].set_title(f"{case_name} DSMC Handoff Time Series (z={target_z_km}km)")
    
    axes[1].plot(t_min, ts['w'], label='w (m/s)', color='green')
    axes[1].set_ylabel("Vertical Vel (m/s)")
    axes[1].legend()
    
    axes[2].plot(t_min, ts['u'], label='u (m/s)', color='purple')
    axes[2].set_ylabel("Horizontal Vel (m/s)")
    axes[2].set_xlabel("Time (min)")
    axes[2].legend()
    
    # Mark steady state region used for FFT
    for ax in axes:
        ax.axvline(t_steady[0] / 60.0, color='k', linestyle='--', alpha=0.5, label='Steady State Start')
        ax.grid(True, alpha=0.3)
        
    out_png = os.path.join(indir, f"{case_name}_dsmc_bc_timeseries.png")
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"  Saved summary plot to {out_png}")

def main():
    # R1: omega = 0.005
    process_case("R1", "output/tier2_R1", 100.0, 0.005)
    
    # R2: omega = 0.005
    process_case("R2", "output/tier2_R2", 100.0, 0.005)
    
    # R3: omega = 0.02
    process_case("R3", "output/tier2_R3", 100.0, 0.02)

if __name__ == "__main__":
    main()
