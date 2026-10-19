import argparse
import glob
import os
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset

def analyze_resonance(run_dir, outdir):
    files = sorted(glob.glob(os.path.join(run_dir, "*.nc")))
    if not files:
        print(f"No files in {run_dir}")
        return

    # Use the latest snapshot
    fpath = files[-1]
    with Dataset(fpath, "r") as ds:
        x = np.asarray(ds.variables["x"][:], dtype=np.float64)
        z = np.asarray(ds.variables["z"][:], dtype=np.float64)
        w = np.asarray(ds.variables["w"][:, :], dtype=np.float64)
        time_s = float(ds.getncattr("time_s"))

    nx = len(x)
    dx = x[1] - x[0]
    lx = x[-1] - x[0]
    
    # Horizontal FFT at each z level
    w_fft = np.fft.rfft(w, axis=1)
    kx = 2 * np.pi * np.fft.rfftfreq(nx, d=dx)
    
    # Power spectrum
    power = np.abs(w_fft)**2
    # Mean power across altitude
    mean_power = np.mean(power, axis=0)
    
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.loglog(kx, mean_power, label="Mean power across z")
    
    # Mark domain-scale frequencies
    fundamental = 2 * np.pi / lx
    for n in range(1, 5):
        ax.axvline(n * fundamental, color="red", linestyle="--", alpha=0.5, label=f"{n}x 2pi/Lx" if n==1 else None)
    
    ax.set_xlabel("kx (rad/m)")
    ax.set_ylabel("Power")
    ax.set_title(f"Horizontal Power Spectrum - {os.path.basename(run_dir)} (t={time_s:.1f}s)")
    ax.legend()
    ax.grid(True, which="both", alpha=0.25)
    
    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, f"resonance_kx_{os.path.basename(run_dir)}.png")
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved resonance plot to {outpath}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--outdir", default="output/resonance_analysis")
    args = parser.parse_args()
    analyze_resonance(args.run_dir, args.outdir)
