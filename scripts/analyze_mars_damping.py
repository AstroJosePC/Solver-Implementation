import os
import glob
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

def analyze_damping():
    outdir_off = "output/mars_mountain_damping/cooling_off"
    outdir_on = "output/mars_mountain_damping/cooling_on"
    
    # Get last snapshot
    files_off = sorted(glob.glob(os.path.join(outdir_off, "snapshot_*.nc")))
    files_on = sorted(glob.glob(os.path.join(outdir_on, "snapshot_*.nc")))
    
    if not files_off or not files_on:
        print("Snapshots not found.")
        return
        
    ds_off = nc.Dataset(files_off[-1], "r")
    ds_on = nc.Dataset(files_on[-1], "r")
    
    z = ds_off.variables['z'][:] / 1e3 # km
    x = ds_off.variables['x'][:] / 1e3 # km
    
    rho_off = ds_off.variables['rho'][:]
    w_off = ds_off.variables['w'][:]
    
    rho_on = ds_on.variables['rho'][:]
    w_on = ds_on.variables['w'][:]
    
    # Extract vertical profile at middle of domain (x = 300 km)
    mid_idx = len(x) // 2
    
    # Get background density (from first snapshot or assume rho_0 approx rho(t=end))
    ds0 = nc.Dataset(files_off[0], "r")
    rho0 = ds0.variables['rho'][:, mid_idx]
    
    # Calculate density-scaled w
    w_scaled_off = w_off[:, mid_idx] * np.sqrt(rho_off[:, mid_idx] / rho0)
    w_scaled_on = w_on[:, mid_idx] * np.sqrt(rho_on[:, mid_idx] / rho0)
    
    # Plot
    plt.figure(figsize=(8, 6))
    plt.plot(w_scaled_off, z, label='Cooling OFF', linestyle='--')
    plt.plot(w_scaled_on, z, label='Cooling ON (Altitude-dependent)')
    plt.xlabel('Density-scaled vertical velocity, $w \\sqrt{\\rho/\\rho_0}$ (m/s)')
    plt.ylabel('Altitude (km)')
    plt.title('Effect of CO2 Radiative Damping on Martian Orographic GW')
    plt.legend()
    plt.grid(True)
    plt.xlim(-10, 10)
    
    out_img = "output/mars_mountain_damping/damping_comparison.png"
    plt.savefig(out_img, dpi=150)
    print(f"Saved comparison plot to {out_img}")

if __name__ == "__main__":
    analyze_damping()