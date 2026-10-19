#!/usr/bin/env python3
"""
plot_critical_profiles_fig8.py
Generates the Figure 8 reproduction: Scaled vertical velocity profiles
at 100, 120, and 150 min, taken at the horizontal center of the wave packet.
"""

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

# Set global plotting style for report consistency
plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12
})

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
    
    times_to_plot = [100, 120, 150]
    colors = ['#32CD32', '#E32636', 'k'] # approx matching paper: lime green, crimson red, black
    
    fig, ax1 = plt.subplots(figsize=(6, 8), constrained_layout=True)
    
    for i, t_min in enumerate(times_to_plot):
        fname = get_snapshot_by_time(indir, t_min)
        print(f"Loading t={t_min} min from {os.path.basename(fname)}")
        
        with Dataset(fname, "r") as ds:
            x = ds.variables["x"][:].astype(np.float64)
            z = ds.variables["z"][:].astype(np.float64)
            w = ds.variables["w"][:, :].astype(np.float64) # (nz, nx)
            rho = ds.variables["rho"][:, :].astype(np.float64)
            
        rho_surface = np.mean(rho[0, :])
        w_scaled = w * np.sqrt(rho / rho_surface)
        
        # Remove horizontal mean to isolate the wave packet
        w_scaled -= np.mean(w_scaled, axis=1, keepdims=True)
        
        # Find the center of the wave packet. We integrate |w_scaled| over altitude 
        # and find the x coordinate where this is maximum.
        w_int = np.sum(np.abs(w_scaled), axis=0)
        idx_center = np.argmax(w_int)
        
        # Extract the profile at the center
        profile = w_scaled[:, idx_center]
        
        # Scale for plotting (paper x-axis is * 10^-4)
        ax1.plot(profile * 1e4, z / 1000.0, color=colors[i], lw=1.5, label=f"t = {t_min} min")

    # Paper limits for bottom x-axis: -1.5 to 1.5 (* 10^-4)
    ax1.set_xlim(-1.5, 1.5)
    ax1.set_xlabel(r"$w\sqrt{\rho/\rho_0}$ (m/s) $\times 10^{-4}$")
    ax1.set_ylim(0, 80)
    ax1.set_ylabel("z (km)")
    
    # Legend position
    ax1.legend(loc='lower left', framealpha=1.0)
    
    # Second axis for (c_x - u_bar)
    ax2 = ax1.twiny()
    
    # Parameters from paper
    omega = 0.01
    lambda_x = 20e3
    kx = 2.0 * np.pi / lambda_x
    c_x = omega / kx # ~ 31.8 m/s
    
    u_bar = background_wind(z)
    c_diff = c_x - u_bar
    
    ax2.plot(c_diff, z / 1000.0, '--', color='#1f77b4', lw=3)
    ax2.set_xlim(-40, 40)
    ax2.set_xlabel(r"$(c_x - \bar{u})$ (m/s)")
    
    # Draw vertical line at c_diff = 0
    ax2.axvline(0, color='k', lw=1.5)

    outpath = os.path.join(outdir, "critical_layer_fig8.png")
    plt.savefig(outpath, dpi=200)
    print(f"Saved: {outpath}")

if __name__ == "__main__":
    main()
