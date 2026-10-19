import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import os

def plot_travel():
    ds = nc.Dataset("output/tier2_R3/compiled/tier2_R3.nc", "r")
    times = ds.variables["time"][:] / 60.0  # minutes
    z_km = ds.variables["z"][:]
    x_km = ds.variables["x"][:]
    w = ds.variables["w"][:]
    rho = ds.variables["rho"][:]

    # Mid-x index
    ix_mid = len(x_km) // 2

    # Scale by density
    rho_surface = np.mean(rho[0, 0, :])
    w_scaled = w[:, :, ix_mid] * np.sqrt(rho[:, :, ix_mid] / rho_surface)

    vmax = float(np.quantile(np.abs(w_scaled), 0.995))

    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.pcolormesh(times, z_km, w_scaled.T, cmap="RdBu_r", shading="auto", vmin=-vmax, vmax=vmax)
    ax.set_xlabel("Time (min)")
    ax.set_ylabel("Altitude z (km)")
    ax.set_title("Wave Travel Diagram (Tier 2 - R3) at mid-x")
    cb = fig.colorbar(im, ax=ax)
    cb.set_label(r"$w\sqrt{\rho/\rho_0}$ (m/s)")

    # Draw handoff line
    ax.axhline(100.0, color='yellow', linestyle='--', label='DSMC Handoff (100 km)')
    ax.legend()

    plt.savefig("output/tier2_R3/travel_diagram.png", dpi=150)
    print("Saved travel diagram to output/tier2_R3/travel_diagram.png")

if __name__ == "__main__":
    plot_travel()
