import os
import glob
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

def extract_amplitude_fft(outdir):
    files = sorted(glob.glob(os.path.join(outdir, "snapshot_*.nc")))
    if not files:
        return None, None
    
    # Take the last ~3 periods (1 period is ~5262s).
    # t_end is 65000. Snapshots are every 1000s, so ~65 files.
    # Take the last 20 snapshots (spanning 20,000s, ~3.8 periods)
    files = files[-20:]
    
    times = []
    w_series = []
    z_coords = None
    
    for f in files:
        ds = nc.Dataset(f, "r")
        t = getattr(ds, "time_s", None)
        if t is None:
             t = float(f.split('_')[-1].split('.')[0]) * 1.0 # fallback approx
        times.append(t)
        
        w_field = ds.variables['w'][:]
        # extract horizontal center (shape is nz, nx)
        w_mid = w_field[:, w_field.shape[1]//2]
        w_series.append(w_mid)
        
        if z_coords is None:
            z_coords = ds.variables['z'][:] / 1e3 # km
        
        ds.close()
        
    times = np.array(times)
    w_series = np.array(w_series) # shape: (len(times), nz)
    
    # We assume uniform time steps in the snapshots
    dt = times[1] - times[0]
    n_samples = len(times)
    
    # Target frequency
    c = 19.0
    lambda_x = 100e3
    omega_target = 2 * np.pi * c / lambda_x
    freq_target = omega_target / (2 * np.pi)
    
    # Compute FFT frequencies
    freqs = np.fft.rfftfreq(n_samples, d=dt)
    
    # Find the closest frequency bin
    idx_target = np.argmin(np.abs(freqs - freq_target))
    
    nz = len(z_coords)
    w_amp = np.zeros(nz)
    
    # Compute FFT for each altitude and extract magnitude at target bin
    for k in range(nz):
        # Windowing the signal to minimize spectral leakage
        w_signal = w_series[:, k] * np.hanning(n_samples)
        
        fft_vals = np.fft.rfft(w_signal)
        # Multiply by 2.0/n_samples to get physical amplitude, 
        # and correct for the Hanning window's amplitude reduction (factor of ~2)
        # Using 4.0/n_samples overall factor for Hanning window amplitude correction
        amp = np.abs(fft_vals[idx_target]) * (4.0 / n_samples) 
        
        w_amp[k] = amp
        
    return z_coords, w_amp

def analyze_damping():
    outdir_off = "output/eckermann_fig14/cooling_off"
    outdir_on = "output/eckermann_fig14/cooling_on"
    
    z, w_amp_off = extract_amplitude_fft(outdir_off)
    _, w_amp_on = extract_amplitude_fft(outdir_on)
    
    if z is None:
        print("No data found.")
        return
        
    c = 19.0
    lambda_x = 100e3
    omega = 2 * np.pi * c / lambda_x
    
    # Convert from w' amplitude to vertical displacement amplitude zeta'
    zeta_off = w_amp_off / omega
    zeta_on = w_amp_on / omega
    
    plt.figure(figsize=(6, 8))
    plt.plot(zeta_off, z, label='Viscous Damping Only', linestyle='--', color='blue')
    plt.plot(zeta_on, z, label='Viscous + Radiative Damping', color='cyan')
    plt.xlabel('Peak Vertical Displacement Amplitude (m)')
    plt.ylabel('Altitude (km)')
    plt.title('Replication of Eckermann Fig 14\n(c=19.0 m/s, True FFT, Log Scale)')
    plt.xscale('log')
    # Since we reduced forcing amplitude to 1.0m, set the xlim bounds
    # from 0.01 to 1000 to observe the exponential growth and decay
    plt.xlim(0.01, 1000)
    plt.ylim(0, 200)
    plt.legend()
    plt.grid(True)
    
    out_img = "output/eckermann_fig14/eckermann_fig14_replication_fixed.png"
    plt.savefig(out_img, dpi=150)
    print(f"Saved comparison plot to {out_img}")

if __name__ == "__main__":
    analyze_damping()