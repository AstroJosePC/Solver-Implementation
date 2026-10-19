import os
import glob
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def harmonic_fit(t, A, B, C):
    c = 6.6
    lambda_x = 100e3
    omega = 2 * np.pi * c / lambda_x
    return A * np.cos(omega * t) + B * np.sin(omega * t) + C

def extract_amplitude(outdir):
    files = sorted(glob.glob(os.path.join(outdir, "snapshot_*.nc")))
    if not files:
        return None, None
    
    # Take the last 20 snapshots (~35,000s to 50,000s, ~1 period is 15,151s)
    files = files[-20:]
    
    times = []
    w_series = []
    z_coords = None
    
    for f in files:
        ds = nc.Dataset(f, "r")
        t = getattr(ds, "time_s", None)
        if t is None:
             t = float(ds.filepath().split('_')[-1].split('.')[0]) * 0.85 # fallback approx
        times.append(t)
        
        w_field = ds.variables['w'][:]
        w_mid = w_field[w_field.shape[0]//2, :]
        w_series.append(w_mid)
        
        if z_coords is None:
            z_coords = ds.variables['z'][:] / 1e3 # km
        
        ds.close()
        
    times = np.array(times)
    w_series = np.array(w_series) # shape: (len(times), nz)
    
    nz = len(z_coords)
    w_amp = np.zeros(nz)
    
    for k in range(nz):
        try:
            popt, _ = curve_fit(harmonic_fit, times, w_series[:, k], p0=[0.0, 0.0, 0.0])
            A, B, C = popt
            w_amp[k] = np.sqrt(A**2 + B**2)
        except:
            w_amp[k] = 0.0
        
    return z_coords, w_amp

def analyze_damping():
    outdir_off = "output/eckermann_fig12/cooling_off"
    outdir_on = "output/eckermann_fig12/cooling_on"
    
    z, w_amp_off = extract_amplitude(outdir_off)
    _, w_amp_on = extract_amplitude(outdir_on)
    
    if z is None:
        print("No data found.")
        return
        
    c = 6.6
    lambda_x = 100e3
    omega = 2 * np.pi * c / lambda_x
    
    zeta_off = w_amp_off / omega
    zeta_on = w_amp_on / omega
    
    plt.figure(figsize=(6, 8))
    plt.plot(zeta_off, z, label='Viscous Damping Only', linestyle='--', color='blue')
    plt.plot(zeta_on, z, label='Viscous + Radiative Damping', color='cyan')
    plt.xlabel('Peak Vertical Displacement Amplitude (m)')
    plt.ylabel('Altitude (km)')
    plt.title('Replication of Eckermann Fig 12(a)\n(Harmonic Fit, Log Scale)')
    plt.xscale('log')
    plt.xlim(0.1, 1000)
    plt.ylim(0, 200)
    plt.legend()
    plt.grid(True)
    
    out_img = "output/eckermann_fig12/eckermann_fig12a_replication_log.png"
    plt.savefig(out_img, dpi=150)
    print(f"Saved comparison plot to {out_img}")

if __name__ == "__main__":
    analyze_damping()