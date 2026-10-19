"""
animate_R3.py

Generates an MP4 animation of the density-scaled vertical velocity field
w_scaled = w * sqrt(rho / rho0) from the Tier 2 R3 compiled NetCDF output.
Each frame is an x-z heatmap matching the style of the frozen snapshots.

Output: output/tier2_R3/wave_animation.mp4
"""

import os
import sys
import numpy as np
import netCDF4 as nc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR   = os.path.dirname(SCRIPT_DIR)
DATA_PATH  = os.path.join(ROOT_DIR, "output", "tier2_R3", "compiled", "tier2_R3.nc")
OUT_PATH   = os.path.join(ROOT_DIR, "output", "tier2_R3", "wave_animation.mp4")

# ---------------------------------------------------------------------------
# Animation parameters
# ---------------------------------------------------------------------------
FPS        = 20          # frames per second in the output video
HANDOFF_KM = 100.0       # DSMC handoff altitude marker [km]
CMAP       = "RdBu_r"
VMAX_PCTILE = 99.5       # percentile used to set symmetric colour limits


def load_data(path: str):
    """Load and return arrays from the compiled NetCDF file."""
    ds = nc.Dataset(path, "r")
    time = ds.variables["time"][:]          # (T,)       seconds
    z    = ds.variables["z"][:]             # (Nz,)      km
    x    = ds.variables["x"][:]             # (Nx,)      km
    rho  = ds.variables["rho"][:]           # (T, Nz, Nx)
    w    = ds.variables["w"][:]             # (T, Nz, Nx)
    ds.close()
    return time, z, x, rho, w


def compute_w_scaled(w: np.ndarray, rho: np.ndarray) -> np.ndarray:
    """Return w * sqrt(rho / rho0) where rho0 is the surface density at t=0."""
    rho0 = float(rho[0, 0, :].mean())      # mean surface density at first frame
    return w * np.sqrt(rho / rho0)


def make_animation(time, z, x, w_scaled):
    """Build and return the FuncAnimation object."""
    z_km = z
    x_km = x
    t_min = time / 60.0

    # Fixed symmetric colour limits from full dataset percentile
    vmax = float(np.nanpercentile(np.abs(w_scaled), VMAX_PCTILE))
    print(f"Colour limit: ±{vmax:.4f} m/s  (from {VMAX_PCTILE}th percentile)")

    fig, ax = plt.subplots(figsize=(7, 6))

    pcm = ax.pcolormesh(
        x_km, z_km, w_scaled[0],
        cmap=CMAP, vmin=-vmax, vmax=vmax,
        shading="auto",
    )

    cb = fig.colorbar(pcm, ax=ax, label=r"$w\sqrt{\rho/\rho_0}$ (m/s)")

    # DSMC handoff line
    ax.axhline(HANDOFF_KM, color="gold", linestyle="--", linewidth=1.5,
               label="DSMC handoff (100 km)")
    ax.legend(loc="upper right", fontsize=8)

    ax.set_xlabel("x (km)")
    ax.set_ylabel("Altitude z (km)")
    title = ax.set_title(f"t = {t_min[0]:.1f} min")

    ax.set_xlim(x_km[0], x_km[-1])
    ax.set_ylim(z_km[0], z_km[-1])

    fig.tight_layout()

    n_frames = len(time)

    def update(frame_idx: int):
        data = w_scaled[frame_idx]
        pcm.set_array(data.ravel())
        title.set_text(f"t = {t_min[frame_idx]:.1f} min")
        return pcm, title

    ani = animation.FuncAnimation(
        fig, update,
        frames=n_frames,
        interval=1000 / FPS,   # ms between frames
        blit=True,
    )
    return ani, fig


def main():
    print(f"Loading data from {DATA_PATH} ...")
    time, z, x, rho, w = load_data(DATA_PATH)
    print(f"  Frames: {len(time)}, z-levels: {len(z)}, x-points: {len(x)}")
    print(f"  Time range: {time[0]/60:.1f} – {time[-1]/60:.1f} min")

    print("Computing density-scaled vertical velocity ...")
    w_scaled = compute_w_scaled(w, rho)

    print("Building animation ...")
    ani, fig = make_animation(time, z, x, w_scaled)

    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    print(f"Saving animation to {OUT_PATH}  (fps={FPS}) ...")
    writer = animation.FFMpegWriter(fps=FPS, bitrate=2000,
                                    metadata={"title": "R3 Acoustic Wave"})
    ani.save(OUT_PATH, writer=writer, dpi=120)
    plt.close(fig)
    print("Done.")


if __name__ == "__main__":
    main()
