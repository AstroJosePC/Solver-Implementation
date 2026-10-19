"""
Visualization suite for the R2 KD-120km simulation.

Generates all standard diagnostics for the Klemp-Durran radiation BC run:
1. Compiles individual snapshots into a single NetCDF
2. snapshots_w_scaled.png — 4 density-scaled w field snapshots
3. wave_breaking_animation.mp4 — animated w*sqrt(rho/rho0) heatmap
4. wavelet_spectrogram_late.png — late-time wavelet power at z=100 km
5. dsmc_handoff_verification.png — w, T, n time series at handoff

Output: output/tier2_R2_kd_120km/
"""

import glob
import os
import sys
import time as timer

import matplotlib
matplotlib.use("Agg")
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset

sys.path.insert(0, "scripts")
from wavelet_tc98 import wavelet

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
RUN_DIR = "output/tier2_R2_kd_120km"
COMPILED_DIR = os.path.join(RUN_DIR, "compiled")
COMPILED_PATH = os.path.join(COMPILED_DIR, "tier2_R2_kd_120km.nc")
HANDOFF_PATH = os.path.join(RUN_DIR, "dsmc_handoff_100km.nc")

HANDOFF_KM = 100.0
OMEGA = 0.005  # rad/s (gravity wave)

# ---------------------------------------------------------------------------
# Step 0: Compile snapshots
# ---------------------------------------------------------------------------

def compile_snapshots():
    """Compile individual snapshot files into a single NetCDF."""
    if os.path.exists(COMPILED_PATH):
        print(f"  Compiled file already exists: {COMPILED_PATH}")
        return

    files = sorted(glob.glob(os.path.join(RUN_DIR, "snapshot_*.nc")))
    if not files:
        raise FileNotFoundError(f"No snapshots in {RUN_DIR}")

    times = []
    for f in files:
        with Dataset(f, "r") as ds:
            times.append(float(ds.getncattr("time_s")))
    times = np.array(times)
    order = np.argsort(times)
    files = [files[i] for i in order]
    times = times[order]

    nt = len(files)
    with Dataset(files[0], "r") as ds0:
        x_m = np.asarray(ds0.variables["x"][:], dtype=np.float64)
        z_m = np.asarray(ds0.variables["z"][:], dtype=np.float64)
        field_names = [name for name, var in ds0.variables.items()
                       if var.ndim == 2 and name not in {"x", "z"}]

    x_km = x_m / 1000.0
    z_km = z_m / 1000.0

    os.makedirs(COMPILED_DIR, exist_ok=True)
    with Dataset(COMPILED_PATH, "w", format="NETCDF4") as out:
        out.createDimension("time", nt)
        out.createDimension("z", len(z_km))
        out.createDimension("x", len(x_km))

        vt = out.createVariable("time", "f8", ("time",))
        vt[:] = times
        vz = out.createVariable("z", "f4", ("z",))
        vz[:] = z_km.astype(np.float32)
        vx = out.createVariable("x", "f4", ("x",))
        vx[:] = x_km.astype(np.float32)

        handles = {}
        for name in field_names:
            handles[name] = out.createVariable(name, "f4", ("time", "z", "x"))

        for i, f in enumerate(files):
            if i % 50 == 0:
                print(f"  Compiling {i}/{nt}...", end="\r")
            with Dataset(f, "r") as ds:
                for name in field_names:
                    handles[name][i] = ds.variables[name][:]

    print(f"  Compiled {nt} snapshots -> {COMPILED_PATH}")


# ---------------------------------------------------------------------------
# Step 1: Snapshots (snapshots_w_scaled.png)
# ---------------------------------------------------------------------------

def plot_snapshots():
    """4-panel density-scaled w snapshots at quartile times."""
    ds = Dataset(COMPILED_PATH, "r")
    times = ds.variables["time"][:]
    z_km = ds.variables["z"][:]
    x_km = ds.variables["x"][:]
    w = ds.variables["w"][:]
    rho = ds.variables["rho"][:]
    ds.close()

    rho_surface = float(np.mean(rho[0, 0, :]))
    nt = len(times)
    indices = [nt // 4, nt // 2, 3 * nt // 4, nt - 1]

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True)
    axes = axes.flatten()

    for i, idx in enumerate(indices):
        w_scaled = w[idx] * np.sqrt(rho[idx] / rho_surface)
        vmax = float(np.quantile(np.abs(w_scaled), 0.995))
        if vmax == 0:
            vmax = 1.0

        im = axes[i].pcolormesh(
            x_km, z_km, w_scaled, cmap="RdBu_r",
            shading="auto", vmin=-vmax, vmax=vmax,
        )
        axes[i].set_title(f"t = {times[idx]/60:.1f} min")
        axes[i].axhline(HANDOFF_KM, color="yellow", linestyle="--",
                        alpha=0.7, label="DSMC handoff")
        fig.colorbar(im, ax=axes[i], label=r"$w\sqrt{\rho/\rho_0}$ (m/s)")

    axes[2].set_xlabel("x (km)")
    axes[3].set_xlabel("x (km)")
    axes[0].set_ylabel("Altitude z (km)")
    axes[2].set_ylabel("Altitude z (km)")
    fig.suptitle("R2 KD-120km: Density-Scaled Vertical Velocity", fontsize=12)

    plt.tight_layout()
    outpath = os.path.join(RUN_DIR, "snapshots_w_scaled.png")
    plt.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"  Saved: {outpath}")


# ---------------------------------------------------------------------------
# Step 2: Animation (wave_breaking_animation.mp4)
# ---------------------------------------------------------------------------

def make_animation():
    """Animated w*sqrt(rho/rho0) heatmap."""
    ds = Dataset(COMPILED_PATH, "r")
    times = ds.variables["time"][:]
    z_km = ds.variables["z"][:]
    x_km = ds.variables["x"][:]
    w = ds.variables["w"][:]
    rho = ds.variables["rho"][:]
    ds.close()

    rho0 = float(np.mean(rho[0, 0, :]))
    w_scaled = w * np.sqrt(rho / rho0)

    t_min = times / 60.0
    vmax = float(np.nanpercentile(np.abs(w_scaled), 99.5))
    print(f"  Animation colour limit: +/-{vmax:.4f} m/s")

    fig, ax = plt.subplots(figsize=(7, 5))
    pcm = ax.pcolormesh(
        x_km, z_km, w_scaled[0],
        cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="auto",
    )
    fig.colorbar(pcm, ax=ax, label=r"$w\sqrt{\rho/\rho_0}$ (m/s)")

    ax.axhline(HANDOFF_KM, color="gold", linestyle="--", linewidth=1.5,
               label="DSMC handoff (100 km)")
    ax.legend(loc="upper right", fontsize=8)
    ax.set_xlabel("x (km)")
    ax.set_ylabel("Altitude z (km)")
    title = ax.set_title(f"t = {t_min[0]:.1f} min")
    fig.tight_layout()

    fps = 20

    def update(frame_idx):
        pcm.set_array(w_scaled[frame_idx].ravel())
        title.set_text(f"t = {t_min[frame_idx]:.1f} min")
        return pcm, title

    ani = animation.FuncAnimation(
        fig, update, frames=len(times),
        interval=1000 / fps, blit=True,
    )

    outpath = os.path.join(RUN_DIR, "wave_breaking_animation.mp4")
    writer = animation.FFMpegWriter(
        fps=fps, bitrate=2000,
        metadata={"title": "R2 KD-120km Wave Animation"},
    )
    ani.save(outpath, writer=writer, dpi=120)
    plt.close(fig)
    print(f"  Saved: {outpath}")


# ---------------------------------------------------------------------------
# Step 3: Wavelet spectrogram (wavelet_spectrogram_late.png)
# ---------------------------------------------------------------------------

def plot_wavelet():
    """Late-time wavelet power spectrum at z=100 km."""
    ds = Dataset(HANDOFF_PATH, "r")
    times = ds.variables["time"][:]
    w = ds.variables["w"][:]
    x = ds.variables["x"][:]
    ds.close()

    ix_mid = len(x) // 2
    w_mid = w[ix_mid, :]

    # Filter for t > 4 hours (post-ramp steady state)
    mask = (times / 3600.0) > 4.0
    times_late = times[mask]
    w_late = w_mid[mask]

    if len(times_late) < 10:
        print("  Not enough data after 4 hrs for wavelet")
        return

    dt_s = float(np.median(np.diff(times_late)))
    s0 = 2.0 * dt_s
    dj = 0.05
    j1 = int(np.fix(np.log2(len(w_late) * dt_s / s0) / dj))

    wave, period_s, _, _ = wavelet(
        w_late - np.mean(w_late), dt_s,
        pad=1, dj=dj, s0=s0, J1=j1,
        mother="MORLET", param=6.0,
    )

    power = np.abs(wave) ** 2
    power = power / max(float(np.max(power)), 1e-30)

    freq_mhz = 1000.0 / period_s
    t_hr = times_late / 3600.0
    freq_forced_mhz = (OMEGA / (2 * np.pi)) * 1000.0

    fig, ax = plt.subplots(figsize=(10, 5))
    im = ax.pcolormesh(t_hr, freq_mhz, power, cmap="jet", shading="auto")
    ax.axhline(freq_forced_mhz, color="w", linestyle="--", linewidth=2,
               label=f"Forced ~{freq_forced_mhz:.2f} mHz")

    ax.set_ylim(0, 3)
    ax.set_xlabel("Time (hrs)")
    ax.set_ylabel("Frequency (mHz)")
    ax.set_title("R2 KD-120km: Wavelet Power (z=100 km, mid-x, t > 4 hrs)")
    ax.legend()
    fig.colorbar(im, ax=ax, label="Normalized Power")

    outpath = os.path.join(RUN_DIR, "wavelet_spectrogram_late.png")
    plt.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"  Saved: {outpath}")


# ---------------------------------------------------------------------------
# Step 4: DSMC handoff verification (dsmc_handoff_verification.png)
# ---------------------------------------------------------------------------

def plot_handoff():
    """w, T, n time series at the handoff altitude."""
    ds = Dataset(HANDOFF_PATH, "r")
    times = ds.variables["time"][:]
    w = ds.variables["w"][:]
    T = ds.variables["T"][:]
    n = ds.variables["n"][:]
    x = ds.variables["x"][:]
    ds.close()

    mid_idx = len(x) // 2
    w_mid = w[mid_idx, :]
    T_mid = T[mid_idx, :]
    n_mid = n[mid_idx, :]

    t_hr = times / 3600.0

    fig, axs = plt.subplots(3, 1, figsize=(10, 8), sharex=True)

    axs[0].plot(t_hr, w_mid, color="blue", linewidth=0.5)
    axs[0].set_ylabel("Vertical Velocity (m/s)")
    axs[0].set_title("R2 KD-120km: DSMC Handoff Data at z = 100 km")
    axs[0].grid(True, alpha=0.3)

    axs[1].plot(t_hr, T_mid, color="red", linewidth=0.5)
    axs[1].set_ylabel("Temperature (K)")
    axs[1].grid(True, alpha=0.3)

    axs[2].plot(t_hr, n_mid, color="green", linewidth=0.5)
    axs[2].set_ylabel(r"Number Density (m$^{-3}$)")
    axs[2].set_xlabel("Time (hrs)")
    axs[2].grid(True, alpha=0.3)

    plt.tight_layout()
    outpath = os.path.join(RUN_DIR, "dsmc_handoff_verification.png")
    plt.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"  Saved: {outpath}")

    print(f"  Max |w| = {np.max(np.abs(w_mid)):.4f} m/s")
    print(f"  Mean T = {np.mean(T_mid):.2f} K, "
          f"Max dT = {np.max(T_mid) - np.mean(T_mid):.2f} K")
    print(f"  Mean n = {np.mean(n_mid):.2e} m^-3, "
          f"dn/n = {(np.max(n_mid) - np.mean(n_mid))/np.mean(n_mid)*100:.2f}%")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    wall_start = timer.time()

    print("=" * 60)
    print("R2 KD-120km Visualization Suite")
    print("=" * 60)

    print("\n[1/5] Compiling snapshots...")
    compile_snapshots()

    print("\n[2/5] Generating snapshots_w_scaled.png...")
    plot_snapshots()

    print("\n[3/5] Generating wave_breaking_animation.mp4...")
    make_animation()

    print("\n[4/5] Generating wavelet_spectrogram_late.png...")
    plot_wavelet()

    print("\n[5/5] Generating dsmc_handoff_verification.png...")
    plot_handoff()

    wall_s = timer.time() - wall_start
    print(f"\nAll done. Wall-clock: {wall_s:.1f} s ({wall_s/60:.1f} min)")


if __name__ == "__main__":
    main()
