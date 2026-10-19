import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

sys.path.insert(0, "scripts")
from wavelet_tc98 import wavelet

def plot_wavelet():
    ds = nc.Dataset("output/tier2_R2/dsmc_handoff_100km.nc", "r")
    times = ds.variables["time"][:]
    w = ds.variables["w"][:]
    x = ds.variables["x"][:]
    
    # Mid-x index
    ix_mid = len(x) // 2
    
    # Extract w(t) at 100 km, mid-x
    w_mid = w[ix_mid, :]
    
    dt_s = float(np.median(np.diff(times)))
    s0 = 2.0 * dt_s
    dj = 0.05
    j1 = int(np.fix(np.log2(len(w_mid) * dt_s / s0) / dj))
    
    wave, period_s, _, _ = wavelet(
        w_mid - np.mean(w_mid), dt_s,
        pad=1, dj=dj, s0=s0, J1=j1,
        mother="MORLET", param=6.0
    )
    
    power = np.abs(wave) ** 2
    power = power / max(float(np.max(power)), 1e-30)
    
    freq_mhz = 1000.0 / period_s
    t_hr = times / 3600.0
    
    fig, ax = plt.subplots(figsize=(10, 5))
    im = ax.pcolormesh(t_hr, freq_mhz, power, cmap="jet", shading="auto")
    
    # Theoretical forced frequency
    omega_forced = 0.005 # rad/s
    freq_forced_mhz = (omega_forced / (2*np.pi)) * 1000.0
    
    ax.axhline(freq_forced_mhz, color="w", linestyle="--", linewidth=2, label=f"Forced ~{freq_forced_mhz:.2f} mHz")
    
    ax.set_ylim(0, 3) # zoom in on low frequencies
    ax.set_xlabel("Time (hours)")
    ax.set_ylabel("Frequency (mHz)")
    ax.set_title("Normalized Wavelet Power Spectrum (z=100km, mid-x)")
    ax.legend()
    fig.colorbar(im, ax=ax, label="Normalized Power")
    
    plt.savefig("output/tier2_R2/wavelet_spectrogram.png", dpi=150)
    print("Saved wavelet spectrogram to output/tier2_R2/wavelet_spectrogram.png")

if __name__ == "__main__":
    plot_wavelet()
