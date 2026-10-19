#!/usr/bin/env python3
"""
plot_critical_spectrum_fig9.py
Generates the Figure 9 reproduction:
(Top) Scaled vertical wavenumber spectrogram over time.
(Bottom) Dominant vertical wavelength evolution with time.
"""

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
})

def main():
    indir = "output/critical_layer_test"
    outdir = "output/critical_analysis"
    os.makedirs(outdir, exist_ok=True)
    
    files = sorted(glob.glob(os.path.join(indir, "*.nc")))
    if not files:
        raise ValueError(f"No files found in {indir}")
        
    times = []
    spectra = []
    dom_wavelengths = []
    
    # Parameters for extracting kx
    lambda_x_target = 20e3
    
    print("Processing snapshots for 2D FFT...")
    for f in files:
        with Dataset(f, "r") as ds:
            x = ds.variables["x"][:].astype(np.float64)
            z = ds.variables["z"][:].astype(np.float64)
            w = ds.variables["w"][:, :].astype(np.float64)
            rho = ds.variables["rho"][:, :].astype(np.float64)
            t = float(ds.getncattr("time_s"))
            
        rho_surface = np.mean(rho[0, :])
        w_scaled = w * np.sqrt(rho / rho_surface)
        w_scaled -= np.mean(w_scaled, axis=1, keepdims=True)
        
        # 2D FFT
        fft2d = np.fft.fft2(w_scaled)
        fft2d_shifted = np.fft.fftshift(fft2d)
        power_spectrum = np.abs(fft2d_shifted)**2
        
        # Frequencies
        dx = x[1] - x[0]
        dz = z[1] - z[0]
        nx = len(x)
        nz = len(z)
        
        kx = np.fft.fftshift(np.fft.fftfreq(nx, d=dx)) * 2 * np.pi
        m = np.fft.fftshift(np.fft.fftfreq(nz, d=dz)) * 2 * np.pi
        
        # Find index for target lambda_x
        # We look for the positive kx that matches 2pi/20km
        k_target = 2 * np.pi / lambda_x_target
        idx_kx = np.argmin(np.abs(kx - k_target))
        
        # Extract m spectrum for this kx
        m_spectrum = power_spectrum[:, idx_kx]
        
        # Normalize by max
        if np.max(m_spectrum) > 0:
            m_spectrum /= np.max(m_spectrum)
            
        # Dominant m (look only at negative m for upward propagation, as stated in paper)
        # Note: in numpy fft, the sign convention might differ. Let's find the peak absolute m.
        # Upward propagation implies phase goes down, but mathematically it depends on the wave eq.
        # We will just find the peak power index.
        idx_m_dom = np.argmax(m_spectrum)
        m_dom = m[idx_m_dom]
        
        if m_dom != 0:
            lambda_z = 2 * np.pi / np.abs(m_dom)
        else:
            lambda_z = np.nan
            
        times.append(t)
        spectra.append(m_spectrum)
        dom_wavelengths.append(lambda_z)
        
    times = np.array(times) / 60.0  # Convert to minutes
    spectra = np.array(spectra).T   # Shape (nz, nt)
    dom_wavelengths = np.array(dom_wavelengths) / 1000.0  # Convert to km
    
    # Generate Plots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), constrained_layout=True)
    
    # Top Panel: Spectrogram
    # Paper m axis is * 10^-4 m^-1, from 0 to -15
    M_scaled = m * 1e4
    
    im = ax1.pcolormesh(times, M_scaled, spectra, cmap="jet", shading='auto')
    ax1.set_ylim(-16, 2)
    ax1.set_xlim(30, 150)
    ax1.set_title("Vertical wavenumber spectrum")
    ax1.set_ylabel(r"$m\ (\times 10^{-4}\ \mathrm{m}^{-1})$")
    fig.colorbar(im, ax=ax1)
    
    # Bottom Panel: Dominant Wavelength
    ax2.plot(times, dom_wavelengths, 'k-', lw=2.5)
    ax2.set_xlim(30, 150)
    ax2.set_ylim(5, 15)
    ax2.set_title("Vertical wavelength")
    ax2.set_ylabel(r"$\lambda_z$ (km)")
    ax2.set_xlabel("t (min)")
    
    outpath = os.path.join(outdir, "critical_layer_fig9.png")
    plt.savefig(outpath, dpi=200)
    print(f"Saved: {outpath}")

if __name__ == "__main__":
    main()
