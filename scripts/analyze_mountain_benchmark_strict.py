#!/usr/bin/env python3
"""
analyze_mountain_benchmark_strict.py
Post-processing for mountain wave benchmark (Srivastava et al. 2022, Sec 3.3).
Produces Figure 10 and Table 11 metrics.
"""

import argparse
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy.interpolate import interp1d

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--indir", default="output/mountain_benchmark_strict_dx1000")
    p.add_argument("--outdir", default="output/mountain_analysis")
    p.add_argument("--gamma", type=float, default=1.4)
    p.add_argument("--g", type=float, default=9.81)
    p.add_argument("--R", type=float, default=287.0)
    p.add_argument("--T0", type=float, default=239.0)
    p.add_argument("--b", type=float, default=10.0)      # hill height [m]
    p.add_argument("--a", type=float, default=20000.0)   # half-width [m]
    p.add_argument("--ubar", type=float, default=18.0)   # mean wind [m/s]
    return p.parse_args()

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # 1. Load latest snapshot
    files = sorted(glob.glob(os.path.join(args.indir, "*.nc")))
    if not files:
        print(f"No files found in {args.indir}")
        return
    
    fname = files[-1]
    print(f"Analyzing {fname}...")
    with Dataset(fname, "r") as ds:
        x = ds.variables["x"][:].astype(np.float64)
        z = ds.variables["z"][:].astype(np.float64)
        rho = ds.variables["rho"][:, :].astype(np.float64) # (nz, nx)
        u = ds.variables["u"][:, :].astype(np.float64)
        w = ds.variables["w"][:, :].astype(np.float64)
        t_sim = ds.getncattr("time_s")

    nx = len(x)
    nz = len(z)
    dx = x[1] - x[0]
    dz = z[1] - z[0]

    # Shift x so mountain center is at 0
    x_center = 0.5 * (x.min() + x.max())
    xx = x - x_center

    # 2. Theory: Queney (1948) Eqs 39-40
    # cs = sqrt(gamma*R*T0), N = sqrt(gamma-1)*g/cs, ks = N/ubar
    cs = np.sqrt(args.gamma * args.R * args.T0)
    N = np.sqrt(args.gamma - 1.0) * args.g / cs
    ks = N / args.ubar
    H = args.R * args.T0 / args.g
    
    # Grid for theory
    X, Z = np.meshgrid(xx, z)
    
    # Eq 39: Streamline displacement zeta
    # zeta = b * sqrt(rho0/rho) * (cos(ks*z) - (x/a)*sin(ks*z)) / (1 + (x/a)^2)
    # sqrt(rho0/rho) = exp(z/2H) for isothermal
    rho_ratio_inv = np.exp(Z / (2.0 * H))
    denom = 1.0 + (X / args.a)**2
    zeta = args.b * rho_ratio_inv * (np.cos(ks * Z) - (X / args.a) * np.sin(ks * Z)) / denom
    
    # Eq 40: w' = u_bar * d(zeta)/dx
    # u' = -u_bar * d(zeta)/dz (so u = u_bar + u' = u_bar(1 - d(zeta)/dz))
    dzeta_dx = np.gradient(zeta, dx, axis=1)
    dzeta_dz = np.gradient(zeta, dz, axis=0)
    
    w_theory = args.ubar * dzeta_dx
    u_theory_pert = -args.ubar * dzeta_dz
    u_theory = args.ubar + u_theory_pert

    # 3. Density-scaled w for Figure 10
    # Paper uses w * sqrt(rho/rho0)
    rho_surface = np.mean(rho[0, :])
    rho_ratio = rho / rho_surface
    w_sim_scaled = w * np.sqrt(rho_ratio)
    
    # Theory scaling
    rho_theory_ratio = np.exp(-Z / H)
    w_theory_scaled = w_theory * np.sqrt(rho_theory_ratio)

    # 4. Scalar metric: u_RMS as % of u_bar (Table 11)
    # Eq 42: RMS error in horizontal wind perturbation
    # Paper uses physical domain excluding sponge (0-40 km)
    u_pert_sim = u - args.ubar
    iz_phys = np.where(z <= 40000.0)[0]
    u_rms = np.sqrt(np.mean((u_pert_sim[iz_phys, :] - u_theory_pert[iz_phys, :])**2))
    u_rms_pct = 100.0 * u_rms / args.ubar
    print(f"u_RMS (0-40km) as % of u_bar = {u_rms_pct:.3f}%")

    # 5. Plotting (Fig 10 style)
    # Plot only z up to 10 km as in Figure 10
    iz_plot = np.where(z <= 10000.0)[0]
    z_km = z[iz_plot] / 1000.0
    x_km = xx / 1000.0
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    
    # clims from paper look like -0.01 to 0.01
    vmax = 0.01
    levels = np.linspace(-vmax, vmax, 21)
    
    # Simulation
    ax = axes[0]
    im = ax.contourf(x_km, z_km, w_sim_scaled[iz_plot, :], levels=levels, cmap="RdBu_r", extend="both")
    ax.set_title(r"$w\sqrt{\rho/\rho_0}$ simulated (m/s)")
    ax.set_xlabel("x (km)")
    ax.set_ylabel("z (km)")
    ax.set_xlim(-100, 100)
    fig.colorbar(im, ax=ax)
    
    # Theory
    ax = axes[1]
    im = ax.contourf(x_km, z_km, w_theory_scaled[iz_plot, :], levels=levels, cmap="RdBu_r", extend="both")
    ax.set_title(r"$w\sqrt{\rho/\rho_0}$ analytical (m/s)")
    ax.set_xlabel("x (km)")
    ax.set_ylabel("z (km)")
    ax.set_xlim(-100, 100)
    fig.colorbar(im, ax=ax)
    
    fig.suptitle(f"Mountain Benchmark (Srivastava Fig 10 reproduction) | dx={dx:.0f}m | u_RMS={u_rms_pct:.3f}%")
    
    out_png = os.path.join(args.outdir, f"mountain_fig10_strict_dx{int(dx)}.png")
    plt.savefig(out_png, dpi=200)
    print(f"Saved plot to {out_png}")

    # Save metrics to CSV
    with open(os.path.join(args.outdir, "mountain_metrics_strict.csv"), "a") as f:
        if f.tell() == 0:
            f.write("dx,t_sim,u_rms_pct\n")
        f.write(f"{dx},{t_sim},{u_rms_pct}\n")

if __name__ == "__main__":
    main()
