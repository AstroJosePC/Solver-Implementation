#!/usr/bin/env python3
"""
plot_merged_backgrounds.py
Quality-control visualization of the merged MCS-MCD background atmospheres
(dayside DD2-like and nightside DD6-like) used for the Tier 2 simulations.

Panels:
  Row 1 (raw CSV columns): Temperature, Density, Pressure, Kinematic viscosity
  Row 2 (derived stability/propagation): Lapse rate vs dry adiabatic, N^2 (with
         omega^2 marker), vertical wavenumber m^2, and a 0-5 km surface zoom that
         overlays the previous (unstable) surface graft to document the fix.

The coupling altitude (100 km) and the sponge layer (110-120 km) are marked.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # headless: write straight to file
import matplotlib.pyplot as plt

# Mars constants (matching src/parameters.jl mars_config)
GAMMA = 1.31
R_GAS = 188.92210861167916
G = 3.72
CP = GAMMA * R_GAS / (GAMMA - 1.0)
ADIAB_LAPSE = G / CP * 1000.0  # dry adiabatic lapse rate [K/km]

# Forced wave parameters (run_tier2_R2_kd_120km_merged.jl)
OMEGA = 0.003
KX = 2.0 * np.pi / 200e3

COUPLE_KM = 100.0
SPONGE_LO_KM = 110.0
SPONGE_HI_KM = 120.0

DAY_C = "#e8772e"    # dayside colour
NIGHT_C = "#1f3b8c"  # nightside colour


def compute_properties(df):
    z = df["z_m"].values
    T = df["T"].values
    rho = df["rho"].values
    P = df["P"].values
    nu = df["nu"].values
    dTdz = np.gradient(T, z)
    N2 = (G / T) * (dTdz + G / CP)
    H = R_GAS * T / G
    m2 = ((N2 - OMEGA**2) / OMEGA**2) * KX**2 - 1.0 / (4.0 * H**2)
    lapse = -dTdz * 1000.0  # K/km
    return dict(z=z, T=T, rho=rho, P=P, nu=nu, N2=N2, m2=m2, lapse=lapse)


def previous_surface_graft(case):
    """
    Reconstruct the bottom of the *previous* (unstable) profile for comparison.
    The old logic grafted the MCD surface temperature at z=0 and ramped linearly
    to the lowest digitized MCS sample at ~728 m. Only z < 728 m was affected.
    """
    # Lowest digitized MCS sample (value, z_km) -> coldest near-surface anchor
    mcs = pd.read_csv("external/digitized_figure_11/temps.csv", header=None)
    z_mcs = mcs.iloc[:, 1].values * 1000.0
    T_mcs = mcs.iloc[:, 0].values
    order = np.argsort(z_mcs)
    z_mcs, T_mcs = z_mcs[order], T_mcs[order]
    z_min, T_min = z_mcs[0], T_mcs[0]

    # MCD surface temperature from the cached profile
    cache = f"scratch/mcd_cache_{case}_t.txt"
    T_mcd_surf = None
    with open(cache) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                T_mcd_surf = float(parts[1])
                break

    z_below = np.array([0.0, 500.0])
    T_below = T_mcd_surf + (z_below / z_min) * (T_min - T_mcd_surf)
    # append the shared anchor point so the line connects to the real profile
    return np.append(z_below, z_min), np.append(T_below, T_min)


def mark_levels(ax, legend=False):
    ax.axhline(COUPLE_KM, color="gray", ls=":", lw=1)
    ax.axhspan(SPONGE_LO_KM, SPONGE_HI_KM, color="gray", alpha=0.12)


day = compute_properties(pd.read_csv("external/digitized_figure_11/mars_mcs_mcd_merged_day.csv"))
night = compute_properties(pd.read_csv("external/digitized_figure_11/mars_mcs_mcd_merged_night.csv"))

fig, axs = plt.subplots(2, 4, figsize=(18, 11))

# --- Row 1: raw profile columns ---
def plot_pair(ax, key, title, xlabel, logx=False):
    ax.plot(day[key], day["z"] / 1000, color=DAY_C, label="Dayside (DD2-like)")
    ax.plot(night[key], night["z"] / 1000, color=NIGHT_C, label="Nightside (DD6-like)")
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    if logx:
        ax.set_xscale("log")
    ax.grid(True, alpha=0.3)
    mark_levels(ax)

plot_pair(axs[0, 0], "T", "Temperature", "T (K)")
axs[0, 0].set_ylabel("Altitude (km)")
axs[0, 0].legend(fontsize=8, loc="upper right")
plot_pair(axs[0, 1], "rho", "Density", r"$\rho$ (kg/m$^3$)", logx=True)
plot_pair(axs[0, 2], "P", "Pressure", "P (Pa)", logx=True)
plot_pair(axs[0, 3], "nu", "Kinematic viscosity", r"$\nu$ (m$^2$/s)", logx=True)

# --- Row 2: derived stability/propagation ---
ax = axs[1, 0]
ax.plot(day["lapse"], day["z"] / 1000, color=DAY_C)
ax.plot(night["lapse"], night["z"] / 1000, color=NIGHT_C)
ax.axvline(ADIAB_LAPSE, color="r", ls="--", lw=1, label=f"dry adiabatic ({ADIAB_LAPSE:.2f} K/km)")
ax.axvline(0, color="k", lw=0.8)
ax.set_title("Lapse rate (superadiabatic if > red)")
ax.set_xlabel("-dT/dz (K/km)")
ax.set_ylabel("Altitude (km)")
ax.legend(fontsize=8, loc="upper right")
ax.grid(True, alpha=0.3)
mark_levels(ax)

ax = axs[1, 1]
ax.plot(day["N2"], day["z"] / 1000, color=DAY_C)
ax.plot(night["N2"], night["z"] / 1000, color=NIGHT_C)
ax.axvline(0, color="k", ls="--", lw=1, label="N²=0 (neutral)")
ax.axvline(OMEGA**2, color="r", ls=":", lw=1, label=r"$\omega^2$")
ax.set_title("Brunt-Vaisala N²")
ax.set_xlabel(r"N² (s$^{-2}$)")
ax.legend(fontsize=8, loc="upper right")
ax.grid(True, alpha=0.3)
mark_levels(ax)

ax = axs[1, 2]
ax.plot(day["m2"], day["z"] / 1000, color=DAY_C)
ax.plot(night["m2"], night["z"] / 1000, color=NIGHT_C)
ax.axvline(0, color="k", ls="--", lw=1, label="m²=0 (evanescent if < 0)")
ax.set_title("Vertical wavenumber m²")
ax.set_xlabel(r"m² (m$^{-2}$)")
ax.set_xlim(-1e-8, 1e-8)
ax.legend(fontsize=8, loc="upper right")
ax.grid(True, alpha=0.3)
mark_levels(ax)

# Surface zoom: current (fixed) vs previous (unstable) graft
ax = axs[1, 3]
zoom = 5.0  # km
for case, props, color in [("day", day, DAY_C), ("night", night, NIGHT_C)]:
    msk = props["z"] / 1000 <= zoom
    ax.plot(props["T"][msk], props["z"][msk] / 1000, color=color, lw=2,
            label=f"{case} (fixed)")
    z_old, T_old = previous_surface_graft(case)
    msk_o = z_old / 1000 <= zoom
    ax.plot(T_old[msk_o], z_old[msk_o] / 1000, color=color, ls="--", lw=1.5,
            marker="x", label=f"{case} (previous, unstable)")
ax.set_title("Surface zoom (0-5 km): fix verification")
ax.set_xlabel("T (K)")
ax.set_ylabel("Altitude (km)")
ax.set_ylim(0, zoom)
ax.legend(fontsize=7, loc="upper right")
ax.grid(True, alpha=0.3)

fig.suptitle("Merged MCS-MCD background atmospheres (corrected surface blend)", fontsize=14)
plt.tight_layout(rect=[0, 0, 1, 0.98])

out_dir = "output/mars_profile_qc"
os.makedirs(out_dir, exist_ok=True)
out_path = os.path.join(out_dir, "merged_backgrounds_qc.png")
plt.savefig(out_path, dpi=150)
print(f"Saved plot to {out_path}")
