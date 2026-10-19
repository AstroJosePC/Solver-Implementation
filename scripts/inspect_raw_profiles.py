import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def read_mcd(filepath):
    z = []
    T = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'): continue
            parts = line.split()
            if len(parts) >= 2:
                z.append(float(parts[0]))
                T.append(float(parts[1]))
    return np.array(z), np.array(T)

# Load raw MCS data
df = pd.read_csv("external/digitized_figure_11/temps.csv", header=None)
T_mcs = df.iloc[:, 0].values
z_mcs = df.iloc[:, 1].values * 1000.0

# Load raw MCD data
z_mcd_day, T_mcd_day = read_mcd("scratch/mcd_cache_day_t.txt")
z_mcd_night, T_mcd_night = read_mcd("scratch/mcd_cache_night_t.txt")

fig, axs = plt.subplots(1, 2, figsize=(12, 6))

# Panel 1: Raw Scatter of all data
axs[0].scatter(T_mcs, z_mcs/1000.0, s=10, label="Raw MCS (temps.csv)", color='black', alpha=0.5)
axs[0].plot(T_mcd_day, z_mcd_day/1000.0, label="MCD Day", color='orange')
axs[0].plot(T_mcd_night, z_mcd_night/1000.0, label="MCD Night", color='darkblue')
axs[0].set_title("Raw Data Profiles")
axs[0].set_ylabel("Altitude (km)")
axs[0].set_xlabel("Temperature (K)")
axs[0].legend()
axs[0].grid(True, alpha=0.3)

# Panel 2: The interpolated MCS profile that gets fed into the solver
idx = np.argsort(z_mcs)
z_mcs_sorted = z_mcs[idx]
T_mcs_sorted = T_mcs[idx]

z_grid = np.arange(0, 80000, 500)
T_mcs_interp = np.interp(z_grid, z_mcs_sorted, T_mcs_sorted)

axs[1].plot(T_mcs_interp, z_grid/1000.0, label="Interpolated MCS", color='red')
axs[1].set_title("Resulting Interpolated MCS")
axs[1].set_xlabel("Temperature (K)")
axs[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(r"C:\Users\Jose\.gemini\antigravity\brain\78110b3d-e698-4e6d-8850-c8a95728718b\raw_data_inspection.png", dpi=150)
print("Saved plot to raw_data_inspection.png")
