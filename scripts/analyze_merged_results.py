#!/usr/bin/env python3
"""
analyze_merged_results.py
Analyzes compiled NetCDF files from dayside (DD2-like) and nightside (DD6-like)
Mars AGW simulations. Compares wave characteristics (amplitude growth, vertical wavelength,
NGIMS matching, and sensible vs. viscous heating) against MAVEN literature constraints.
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import scipy.interpolate as interp

# Mars Physical Constants
G_MARS = 3.72
R_GAS_MARS = 188.92210861167916
GAMMA = 1.31
CP = (GAMMA * R_GAS_MARS) / (GAMMA - 1.0)  # ~798.3 J/(kg*K)

def compile_snapshots(run_dir, compiled_path):
    """Compile individual snapshot files into a single NetCDF if it doesn't exist."""
    if os.path.exists(compiled_path):
        print(f"Compiled file already exists: {compiled_path}")
        return

    files = sorted(glob.glob(os.path.join(run_dir, "snapshot_*.nc")))
    if not files:
        raise FileNotFoundError(f"No snapshots found in {run_dir}")

    times = []
    for f in files:
        with Dataset(f, "r") as ds:
            times.append(float(ds.getncattr("time_s")))
    times = np.array(times)
    order = np.argsort(times)
    files = [files[i] for i in order]
    times = times[order]

    nt = len(files)
    with Dataset(files[0], "r") as ds0:
        x_m = np.asarray(ds0.variables["x"][:], dtype=np.float64)
        z_m = np.asarray(ds0.variables["z"][:], dtype=np.float64)
        field_names = [name for name, var in ds0.variables.items()
                       if var.ndim == 2 and name not in {"x", "z"}]

    os.makedirs(os.path.dirname(compiled_path), exist_ok=True)
    with Dataset(compiled_path, "w", format="NETCDF4") as out:
        out.createDimension("time", nt)
        out.createDimension("z", len(z_m))
        out.createDimension("x", len(x_m))

        vt = out.createVariable("time", "f8", ("time",))
        vt[:] = times
        vz = out.createVariable("z", "f8", ("z",))
        vz[:] = z_m
        vx = out.createVariable("x", "f8", ("x",))
        vx[:] = x_m

        handles = {}
        for name in field_names:
            handles[name] = out.createVariable(name, "f4", ("time", "z", "x"))

        for i, f in enumerate(files):
            with Dataset(f, "r") as ds:
                for name in field_names:
                    handles[name][i] = ds.variables[name][:]

    print(f"Successfully compiled {nt} snapshots -> {compiled_path}")

def analyze_case(case_name, compiled_path):
    print(f"\n================================================================================")
    print(f"ANALYZING CASE: {case_name}")
    print(f"================================================================================")
    
    with Dataset(compiled_path, "r") as ds:
        times = ds.variables["time"][:]
        z = ds.variables["z"][:]       # height in meters
        x = ds.variables["x"][:]       # distance in meters
        w = ds.variables["w"][:]       # (time, z, x)
        u = ds.variables["u"][:]       # (time, z, x)
        rho = ds.variables["rho"][:]   # (time, z, x)
        T = ds.variables["T"][:]       # (time, z, x)
        P = ds.variables["P"][:]       # (time, z, x)
        
    dx = x[1] - x[0]
    dz = z[1] - z[0]
    nz = len(z)
    nx = len(x)
    
    # We analyze the late-time steady state (last 20% of the simulation)
    late_idx = int(0.8 * len(times))
    print(f"Analyzing time-averaged state from t = {times[late_idx]/3600:.1f} to {times[-1]/3600:.1f} hours.")
    
    # Time-averaged background density, temperature and pressure
    rho_bar = np.mean(rho[late_idx:, :, :], axis=(0, 2))
    T_bar = np.mean(T[late_idx:, :, :], axis=(0, 2))
    P_bar = np.mean(P[late_idx:, :, :], axis=(0, 2))
    
    # 1. Amplitude Growth check: energy flux conservation w_scaled = w * sqrt(rho / rho_surf)
    # Average wave amplitude as a function of height
    w_amp = np.zeros(nz)
    w_scaled_amp = np.zeros(nz)
    rho_surf = rho_bar[0]
    
    for iz in range(nz):
        # perturbation velocity standard deviation
        w_pert = w[late_idx:, iz, :] - np.mean(w[late_idx:, iz, :], axis=1, keepdims=True)
        w_amp[iz] = np.mean(np.std(w_pert, axis=1))
        w_scaled_amp[iz] = w_amp[iz] * np.sqrt(rho_bar[iz] / rho_surf)
        
    # 2. Phase-tilt analysis to measure vertical wavelength lambda_z
    # Forced horizontal wavenumber kx
    lambda_x = 200e3
    kx = 2.0 * np.pi / lambda_x
    # Index of forced mode (k = 1 wave across the 200 km domain)
    k_idx = 1  
    
    phases = []
    for iz in range(nz):
        # time-mean complex Fourier coeff at each altitude
        w_ft_z = []
        for it in range(late_idx, len(times)):
            w_z = w[it, iz, :]
            ft = np.fft.fft(w_z - np.mean(w_z))
            w_ft_z.append(ft[k_idx])
        mean_ft = np.mean(w_ft_z)
        phases.append(np.angle(mean_ft))
        
    phases = np.unwrap(phases)
    
    # Derivative of phase w.r.t z gives vertical wavenumber m
    # Smooth phase first to avoid noise in derivative
    window_len = 11
    if nz > window_len:
        smoothed_phases = pd.Series(phases).rolling(window=window_len, center=True).mean().bfill().ffill().values
    else:
        smoothed_phases = phases
        
    m = np.gradient(smoothed_phases, dz)
    lambda_z = np.zeros(nz)
    for iz in range(nz):
        if abs(m[iz]) > 1e-8:
            lambda_z[iz] = (2.0 * np.pi / abs(m[iz])) / 1000.0  # in km
        else:
            lambda_z[iz] = np.nan
            
    # Compute theoretical lambda_z from dispersion relation
    # N^2 = (gamma - 1) * g^2 / (c_s^2) + (g / T) * dT/dz
    # For Mars, c_s^2 = gamma * R_gas * T
    N2 = np.zeros(nz)
    dTdz = np.gradient(T_bar, dz)
    for iz in range(nz):
        cs2 = GAMMA * R_GAS_MARS * T_bar[iz]
        N2[iz] = (GAMMA - 1.0) * G_MARS**2 / cs2 + (G_MARS / T_bar[iz]) * dTdz[iz]
        
    omega_forced = 0.003  # must match run_tier2_R2_kd_120km_merged.jl
    lambda_z_theory = np.zeros(nz)
    for iz in range(nz):
        # m^2 = (N^2 / omega^2 - 1) * kx^2 - 1/(4H^2)
        # where scale height H = R_gas * T / g
        H = R_GAS_MARS * T_bar[iz] / G_MARS
        m2_th = (N2[iz] / omega_forced**2 - 1.0) * kx**2 - 1.0 / (4.0 * H**2)
        if m2_th > 0:
            lambda_z_theory[iz] = (2.0 * np.pi / np.sqrt(m2_th)) / 1000.0 # km
        else:
            lambda_z_theory[iz] = np.nan

    # 3. NGIMS Amplitude Matching: relative density perturbation delta_rho / rho_0
    rel_rho_amp = np.zeros(nz)
    for iz in range(nz):
        rho_pert = (rho[late_idx:, iz, :] - rho_bar[iz]) / rho_bar[iz]
        rel_rho_amp[iz] = np.mean(np.std(rho_pert, axis=1)) * np.sqrt(2.0) # peak relative amplitude
        
    # 4. Sensible heat flux and viscous heating
    Q_sensible = np.zeros(nz)
    Q_visc = np.zeros(nz)
    
    # Pre-calculate derivatives of velocity perturbations
    dudz_list = []
    dudx_list = []
    dwdz_list = []
    dwdx_list = []
    
    for it in range(late_idx, len(times)):
        dudz_list.append(np.gradient(u[it], dz, axis=0))
        dudx_list.append(np.gradient(u[it], dx, axis=1))
        dwdz_list.append(np.gradient(w[it], dz, axis=0))
        dwdx_list.append(np.gradient(w[it], dx, axis=1))
        
    dudz = np.array(dudz_list)
    dudx = np.array(dudx_list)
    dwdz = np.array(dwdz_list)
    dwdx = np.array(dwdx_list)
    
    # Viscosity mu(T) = 4.2e-7 * T^0.69
    for iz in range(nz):
        T_mean_z = T_bar[iz]
        mu = 4.2e-7 * T_mean_z**0.69
        
        # Viscous dissipation rate epsilon
        term1 = 2.0 * dudx[:, iz, :]**2
        term2 = 2.0 * dwdz[:, iz, :]**2
        term3 = (dudz[:, iz, :] + dwdx[:, iz, :])**2
        term4 = -(2.0/3.0) * (dudx[:, iz, :] + dwdz[:, iz, :])**2
        epsilon = mu * (term1 + term2 + term3 + term4)
        
        # Time and horizontal average of dissipation
        epsilon_bar = np.mean(epsilon)
        Q_visc[iz] = epsilon_bar / (rho_bar[iz] * CP) * 86400.0  # K/sol
        
        # Sensible heat flux: F_sens = rho0 * cp * <w' T'>
        w_pert = w[late_idx:, iz, :] - np.mean(w[late_idx:, iz, :], axis=1, keepdims=True)
        T_pert = T[late_idx:, iz, :] - np.mean(T[late_idx:, iz, :], axis=1, keepdims=True)
        wT_cov = np.mean(w_pert * T_pert)
        Q_sensible[iz] = wT_cov # temporary store covariance
        
    # Compute derivative of sensible heat flux: Q_sens = - 1/rho0 * d(rho0 * <w'T'>)/dz
    flux_profile = rho_bar * Q_sensible
    dflux_dz = np.gradient(flux_profile, dz)
    Q_sensible_rate = - dflux_dz / rho_bar * 86400.0  # K/sol
    
    # 5. Summarize key statistics
    z_km = z / 1000.0
    idx_100 = np.argmin(np.abs(z_km - 100.0))
    idx_120 = np.argmin(np.abs(z_km - 120.0))
    
    # Wave energy flux check: look at scaled amplitude flatness in 20-80 km
    mask_flat = (z_km >= 20.0) & (z_km <= 80.0)
    w_scaled_flat = w_scaled_amp[mask_flat]
    growth_std_rel = np.std(w_scaled_flat) / np.mean(w_scaled_flat)
    
    results = {
        "z_km": z_km,
        "w_amp": w_amp,
        "w_scaled_amp": w_scaled_amp,
        "rel_rho_amp": rel_rho_amp,
        "lambda_z": lambda_z,
        "lambda_z_theory": lambda_z_theory,
        "Q_sensible": Q_sensible_rate,
        "Q_visc": Q_visc,
        "w_amp_100": w_amp[idx_100],
        "rel_rho_100": rel_rho_amp[idx_100],
        "rel_rho_120": rel_rho_amp[idx_120],
        "lambda_z_100": lambda_z[idx_100],
        "lambda_z_th_100": lambda_z_theory[idx_100],
        "growth_flatness": growth_std_rel,
        "Q_sens_110": Q_sensible_rate[np.argmin(np.abs(z_km - 110.0))],
        "Q_visc_110": Q_visc[np.argmin(np.abs(z_km - 110.0))],
        "P_surf": P_bar[0]
    }
    
    # Report findings
    print(f"--- Verification against MAVEN Constraints ---")
    print(f"1. Energy Flux Flatness (20-80 km relative std): {growth_std_rel*100:.2f}% (Target: low, showing conservation)")
    print(f"2. Vertical Wavelength at 100 km: {results['lambda_z_100']:.1f} km (Theory: {results['lambda_z_th_100']:.1f} km)")
    print(f"3. Relative Density Perturbation at 100 km: {results['rel_rho_100']*100:.2f}%")
    print(f"4. Relative Density Perturbation at 120 km: {results['rel_rho_120']*100:.2f}% (Target: Dayside 20-40%, Nightside up to 70%)")
    print(f"5. Thermal effects at 110 km:")
    print(f"   - Sensible heat flux rate: {results['Q_sens_110']:.2f} K/sol")
    print(f"   - Viscous dissipation rate: {results['Q_visc_110']:.2f} K/sol")
    print(f"   - Net wave thermal effect: {results['Q_sens_110'] + results['Q_visc_110']:.2f} K/sol")
    
    return results

def plot_verification(results_day, results_night):
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Panel 1: Vertical velocity amplitude profiles
    ax = axes[0, 0]
    ax.plot(results_day["w_amp"], results_day["z_km"], 'b-', label="Day (raw w)")
    ax.plot(results_day["w_scaled_amp"], results_day["z_km"], 'b--', label="Day (scaled w)")
    ax.plot(results_night["w_amp"], results_night["z_km"], 'r-', label="Night (raw w)")
    ax.plot(results_night["w_scaled_amp"], results_night["z_km"], 'r--', label="Night (scaled w)")
    ax.axhline(100.0, color='gray', linestyle=':', label="Handoff (100 km)")
    ax.set_xlabel("Velocity Amplitude (m/s)")
    ax.set_ylabel("Altitude z (km)")
    ax.set_title("Wave Amplitude Profiles")
    ax.legend()
    ax.grid(True)
    
    # Panel 2: Relative density perturbation profiles (MAVEN calibration)
    ax = axes[0, 1]
    ax.plot(results_day["rel_rho_amp"] * 100.0, results_day["z_km"], 'b-', label="Dayside (DD2-like)")
    ax.plot(results_night["rel_rho_amp"] * 100.0, results_night["z_km"], 'r-', label="Nightside (DD6-like)")
    ax.axvspan(20.0, 40.0, color='blue', alpha=0.1, label="Dayside Target (20-40%)")
    ax.axvspan(40.0, 70.0, color='red', alpha=0.1, label="Nightside Target (up to 70%)")
    ax.axhline(100.0, color='gray', linestyle=':')
    ax.axhline(120.0, color='k', linestyle='--', label="Domain Top (120 km)")
    ax.set_xlabel("Relative Density Perturbation delta_rho/rho_0 (%)")
    ax.set_ylabel("Altitude z (km)")
    ax.set_title("Relative Density Perturbation (Calibration)")
    ax.legend(loc="lower right")
    ax.grid(True)
    
    # Panel 3: Vertical wavelength comparison
    ax = axes[1, 0]
    ax.plot(results_day["lambda_z"], results_day["z_km"], 'b-', label="Day (Measured)")
    ax.plot(results_day["lambda_z_theory"], results_day["z_km"], 'b--', label="Day (Theory)")
    ax.plot(results_night["lambda_z"], results_night["z_km"], 'r-', label="Night (Measured)")
    ax.plot(results_night["lambda_z_theory"], results_night["z_km"], 'r--', label="Night (Theory)")
    ax.axhline(100.0, color='gray', linestyle=':')
    ax.set_xlabel("Vertical Wavelength lambda_z (km)")
    ax.set_ylabel("Altitude z (km)")
    ax.set_title("Vertical Wavelength vs. Dispersion Relation")
    ax.legend()
    ax.grid(True)
    
    # Panel 4: Heating/cooling rates
    ax = axes[1, 1]
    ax.plot(results_day["Q_sensible"], results_day["z_km"], 'b-', label="Day Sensible")
    ax.plot(results_day["Q_visc"], results_day["z_km"], 'b--', label="Day Viscous")
    ax.plot(results_night["Q_sensible"], results_night["z_km"], 'r-', label="Night Sensible")
    ax.plot(results_night["Q_visc"], results_night["z_km"], 'r--', label="Night Viscous")
    ax.axvline(0.0, color='black', linestyle='-')
    ax.axhline(100.0, color='gray', linestyle=':')
    ax.set_xlabel("Heating Rate (K/sol)")
    ax.set_ylabel("Altitude z (km)")
    ax.set_title("Thermal Heating and Cooling Rates")
    ax.legend()
    ax.grid(True)
    
    plt.tight_layout()
    outpath = "output/tier2_merged_verification_summary.png"
    plt.savefig(outpath, dpi=150)
    plt.close()
    print(f"\nSaved combined verification plot to: {outpath}")

def main():
    day_run_dir = "output/tier2_R2_kd_120km_merged_day"
    night_run_dir = "output/tier2_R2_kd_120km_merged_night"
    
    day_compiled = os.path.join(day_run_dir, "compiled.nc")
    night_compiled = os.path.join(night_run_dir, "compiled.nc")
    
    print("Compiling dayside snapshots...")
    compile_snapshots(day_run_dir, day_compiled)
    
    print("\nCompiling nightside snapshots...")
    compile_snapshots(night_run_dir, night_compiled)
    
    results_day = analyze_case("Dayside (DD2-like)", day_compiled)
    results_night = analyze_case("Nightside (DD6-like)", night_compiled)
    
    plot_verification(results_day, results_night)
    
    # Write a clean markdown summary of the results to be loaded by main script
    summary_path = "output/verification_summary.md"
    with open(summary_path, "w") as f:
        f.write("# Verification Summary of Merged Atmosphere Simulation Sweep\n\n")
        f.write("## Overview\n")
        f.write("Gravity wave simulations on a 120 km domain using the Klemp-Durran radiation BC were run and validated against MAVEN literature constraints.\n\n")
        f.write("## Key Results Table\n\n")
        f.write("| Parameter | Dayside Case (DD2-like) | Nightside Case (DD6-like) | Target / Constraints |\n")
        f.write("| --- | --- | --- | --- |\n")
        f.write(f"| Surface Pressure $P_0$ | {results_day['P_surf']:.1f} Pa | {results_night['P_surf']:.1f} Pa | Matches MCD at Equator |\n")
        f.write(f"| Energy Flux conservation (20-80 km relative std) | {results_day['growth_flatness']*100:.2f}% | {results_night['growth_flatness']*100:.2f}% | Low value (<15% preferred) |\n")
        f.write(f"| Vertical wavelength $\\lambda_z$ at 100 km (Measured / Theory) | {results_day['lambda_z_100']:.1f} / {results_day['lambda_z_th_100']:.1f} km | {results_night['lambda_z_100']:.1f} / {results_night['lambda_z_th_100']:.1f} km | Matches dispersion relation |\n")
        f.write(f"| Relative density perturbation $\\delta\\rho/\\rho_0$ at 100 km | {results_day['rel_rho_100']*100:.2f}% | {results_night['rel_rho_100']*100:.2f}% | Model coupling boundary |\n")
        f.write(f"| Relative density perturbation $\\delta\\rho/\\rho_0$ at 120 km | {results_day['rel_rho_120']*100:.2f}% | {results_night['rel_rho_120']*100:.2f}% | Day: 20-40%, Night: up to 70% |\n")
        f.write(f"| Net wave thermal rate at 110 km (Sensible / Viscous) | {results_day['Q_sens_110']:.2f} / {results_day['Q_visc_110']:.2f} K/sol | {results_night['Q_sens_110']:.2f} / {results_night['Q_visc_110']:.2f} K/sol | Net Cooling at Thermosphere |\n")
        
    print(f"\nWritten verification summary markdown to {summary_path}")

if __name__ == "__main__":
    main()
