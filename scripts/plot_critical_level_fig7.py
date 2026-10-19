#!/usr/bin/env python3
"""
plot_critical_level_fig7.py
Generates the Figure 7 reproduction for the critical level test.
Plots density-scaled vertical velocity at t=50, 80, 120, 150 mins.
"""

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

def get_snapshot_by_time(indir, target_time_min, tolerance=5.0):
    files = sorted(glob.glob(os.path.join(indir, "*.nc")))
    target_s = target_time_min * 60.0
    
    closest_file = None
    min_diff = float('inf')
    
    for f in files:
        with Dataset(f, "r") as ds:
            t = float(ds.getncattr("time_s"))
            diff = abs(t - target_s)
            if diff < min_diff:
                min_diff = diff
                closest_file = f
                
    if min_diff <= tolerance * 60.0:
        return closest_file
    else:
        raise ValueError(f"Could not find snapshot near {target_time_min} mins")

def background_wind(z):
    u_0 = 100.0
    z_0 = 80000.0
    sigma_u = 10000.0
    return u_0 * np.exp(-(z - z_0)**2 / (2.0 * sigma_u**2))

def main():
    indir = "output/critical_layer_test"
    outdir = "output/critical_analysis"
    os.makedirs(outdir, exist_ok=True)
    
    times_to_plot = [50, 80, 120, 150]
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)
    axes = axes.flatten()
    
    for i, t_min in enumerate(times_to_plot):
        fname = get_snapshot_by_time(indir, t_min)
        print(f"Plotting t={t_min} min using {os.path.basename(fname)}")
        
        with Dataset(fname, "r") as ds:
            x = ds.variables["x"][:].astype(np.float64)
            z = ds.variables["z"][:].astype(np.float64)
            w = ds.variables["w"][:, :].astype(np.float64) # (nz, nx) Note: verify if it's (nz, nx) or (nx, nz). The previous fix showed it's (nz, nx) for some, but let's be careful.
            # In Julia, write_snapshot_nc writes (x, z). So in python netCDF4 it reads as (z, x).
            rho = ds.variables["rho"][:, :].astype(np.float64)
            
        rho_surface = np.mean(rho[0, :])
        w_scaled = w * np.sqrt(rho / rho_surface)
        
        # Remove horizontal mean to filter out large-scale acoustic transients
        # and background numerical drift, isolating the gravity wave packet.
        w_scaled -= np.mean(w_scaled, axis=1, keepdims=True)
        
        # Grid in km (shift X by -200km to center at 0)
        X, Z = np.meshgrid((x / 1000.0) - 200.0, z / 1000.0)
        
        ax = axes[i]
        
        # Color map range logic from paper: max amplitude changes over time
        cmax = np.max(np.abs(w_scaled))
        # Ensure symmetric colorbar
        levels = np.linspace(-cmax, cmax, 51)
        
        im = ax.contourf(X, Z, w_scaled, levels=levels, cmap="jet", extend="both")
        
        # Set axis limits to match paper
        ax.set_xlim(-200, 200)
        ax.set_ylim(0, 160)
        ax.set_xlabel("x (km)")
        ax.set_ylabel("z (km)")
        
        # Overlay mean wind profile
        # Use a twin top axis for the wind profile
        ax_top = ax.twiny()
        u_bar = background_wind(z)
        ax_top.plot(u_bar, z / 1000.0, 'k-', lw=2.0)
        ax_top.set_xlim(0, 120)  # wind goes up to 100 m/s
        ax_top.set_xlabel(r"$\bar{u}$ (m/s)")
        
        ax.text(0.8, 0.9, f"$t = {t_min}$ min", transform=ax.transAxes, 
                fontsize=16, fontweight='bold', bbox=dict(facecolor='white', alpha=0.7))
        
        fig.colorbar(im, ax=ax, orientation='vertical')

    plt.suptitle("Critical Layer Evolution (Srivastava et al. Fig 7 Reproduction)", fontsize=18)
    
    outpath = os.path.join(outdir, "critical_layer_fig7.png")
    plt.savefig(outpath, dpi=200)
    print(f"Saved: {outpath}")

if __name__ == "__main__":
    main()
