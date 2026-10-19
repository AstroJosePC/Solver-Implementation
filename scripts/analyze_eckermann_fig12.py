import os
import glob
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

def analyze_damping():
    outdir_off = "output/eckermann_fig12/cooling_off"
    outdir_on = "output/eckermann_fig12/cooling_on"
    
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
    
    w_off = ds_off.variables['w'][:]
    w_on = ds_on.variables['w'][:]
    
    # Wave parameters
    c = 6.6
    lambda_x = 100e3
    omega = 2 * np.pi * c / lambda_x
    
    # Calculate peak amplitude at each altitude
    # Since it's a harmonic wave, we can just take the max over x (or max-min/2)
    # Exclude ghost cells if necessary, but we are using the middle anyway.
    w_amp_off = (np.max(w_off, axis=1) - np.min(w_off, axis=1)) / 2.0
    w_amp_on = (np.max(w_on, axis=1) - np.min(w_on, axis=1)) / 2.0
    
    # Convert velocity amplitude to displacement amplitude: zeta = w / omega
    zeta_off = w_amp_off / omega
    zeta_on = w_amp_on / omega
    
    # Plot
    plt.figure(figsize=(6, 8))
    plt.plot(zeta_off, z, label='Viscous Damping Only (Run A)', linestyle='--', color='blue')
    plt.plot(zeta_on, z, label='Viscous + Radiative Damping (Run B)', color='cyan')
    plt.xlabel('Peak Vertical Displacement Amplitude (m)')
    plt.ylabel('Altitude (km)')
    plt.title('Replication of Eckermann Fig 12(a)\nEffect of CO2 Radiative Damping')
    plt.legend()
    plt.grid(True)
    # plt.xlim(0, 1000) # Optional limit
    plt.ylim(0, 200)
    
    out_img = "output/eckermann_fig12/eckermann_fig12a_replication.png"
    plt.savefig(out_img, dpi=150)
    print(f"Saved comparison plot to {out_img}")

if __name__ == "__main__":
    analyze_damping()
