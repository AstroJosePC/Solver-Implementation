import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import os

def plot_snapshots():
    ds = nc.Dataset("output/tier2_R2/compiled/tier2_R2.nc", "r")
    times = ds.variables["time"][:]
    z_km = ds.variables["z"][:]
    x_km = ds.variables["x"][:]
    w = ds.variables["w"][:]
    rho = ds.variables["rho"][:]
    
    rho_surface = np.mean(rho[0, 0, :])
    
    # Pick a few snapshots
    nt = len(times)
    indices = [nt//4, nt//2, 3*nt//4, nt-1]
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10), sharex=True, sharey=True)
    axes = axes.flatten()
    
    for i, idx in enumerate(indices):
        w_scaled = w[idx, :, :] * np.sqrt(rho[idx, :, :] / rho_surface)
        vmax = float(np.quantile(np.abs(w_scaled), 0.995))
        if vmax == 0: vmax = 1.0
        
        im = axes[i].pcolormesh(x_km, z_km, w_scaled, cmap="RdBu_r", shading="auto", vmin=-vmax, vmax=vmax)
        axes[i].set_title(f"t = {times[idx]/60.0:.1f} min")
        axes[i].axhline(100.0, color='yellow', linestyle='--', alpha=0.5)
        fig.colorbar(im, ax=axes[i], label=r"$w\sqrt{\rho/\rho_0}$ (m/s)")
        
    axes[2].set_xlabel("x (km)")
    axes[3].set_xlabel("x (km)")
    axes[0].set_ylabel("Altitude z (km)")
    axes[2].set_ylabel("Altitude z (km)")
    
    plt.tight_layout()
    plt.savefig("output/tier2_R2/snapshots_w_scaled.png", dpi=150)
    print("Saved snapshots to output/tier2_R2/snapshots_w_scaled.png")

if __name__ == "__main__":
    plot_snapshots()
