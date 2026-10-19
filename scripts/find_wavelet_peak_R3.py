"""
Find the time of peak wavelet power at the forced frequency (3.18 mHz)
in the R3 handoff data, then extract and plot a snapshot at that time.
"""
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.insert(0, "scripts")
from wavelet_tc98 import wavelet

# --- Step 1: Find peak power time from handoff wavelet ---
ds = nc.Dataset("output/tier2_R3/dsmc_handoff_100km.nc", "r")
times = ds.variables["time"][:]
w = ds.variables["w"][:]
x = ds.variables["x"][:]
ix_mid = len(x) // 2
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
freq_mhz = 1000.0 / period_s

# Find the frequency bin closest to 3.18 mHz
omega_forced = 0.02
freq_forced_mhz = (omega_forced / (2 * np.pi)) * 1000.0
i_freq = np.argmin(np.abs(freq_mhz - freq_forced_mhz))

print(f"Forced frequency: {freq_forced_mhz:.2f} mHz")
print(f"Closest wavelet bin: {freq_mhz[i_freq]:.2f} mHz")

# Power at the forced frequency as a function of time
power_at_forced = power[i_freq, :]
i_peak = np.argmax(power_at_forced)
t_peak_s = times[i_peak]
t_peak_hr = t_peak_s / 3600.0

print(f"Peak power time: {t_peak_s:.0f} s = {t_peak_hr:.2f} hrs")
print(f"Power at peak: {power_at_forced[i_peak]:.4e}")
print(f"Power at end: {power_at_forced[-1]:.4e}")
print(f"Ratio (end/peak): {power_at_forced[-1]/power_at_forced[i_peak]:.4f}")

# --- Step 2: Find the closest snapshot time in compiled data ---
ds2 = nc.Dataset("output/tier2_R3/compiled/tier2_R3.nc", "r")
snap_times = ds2.variables["time"][:]
z_km = ds2.variables["z"][:]
x_km = ds2.variables["x"][:]

i_snap = np.argmin(np.abs(snap_times - t_peak_s))
t_snap = snap_times[i_snap]
print(f"\nClosest snapshot: index {i_snap}, t = {t_snap:.0f} s = {t_snap/3600:.2f} hrs")

# --- Step 3: Extract and plot snapshot at peak-power time ---
w_snap = ds2.variables["w"][i_snap, :, :]
rho_snap = ds2.variables["rho"][i_snap, :, :]
rho_surface = np.mean(ds2.variables["rho"][0, 0, :])
w_scaled = w_snap * np.sqrt(rho_snap / rho_surface)

# Also grab an early snapshot (during ramp-up) and a late snapshot for comparison
# Early: ~1 hr (3600 s), Late: last snapshot
i_early = np.argmin(np.abs(snap_times - 3600.0))
i_late = len(snap_times) - 1

fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)

for ax, idx, label in zip(axes,
                           [i_early, i_snap, i_late],
                           ["Early (ramp-up)", "Peak wavelet power", "End of simulation"]):
    w_s = ds2.variables["w"][idx, :, :] * np.sqrt(ds2.variables["rho"][idx, :, :] / rho_surface)
    vmax = float(np.quantile(np.abs(w_s), 0.995))
    if vmax == 0:
        vmax = 1.0

    im = ax.pcolormesh(x_km, z_km, w_s, cmap="RdBu_r", shading="auto", vmin=-vmax, vmax=vmax)
    ax.axhline(100.0, color='yellow', linestyle='--', alpha=0.7, label='DSMC handoff')
    ax.set_title(f"{label}\nt = {snap_times[idx]/3600:.2f} hrs", fontsize=11)
    ax.set_xlabel("x (km)")
    fig.colorbar(im, ax=ax, label=r"$w\sqrt{\rho/\rho_0}$ (m/s)", shrink=0.85)

axes[0].set_ylabel("Altitude z (km)")
axes[0].legend(fontsize=8)

fig.suptitle(f"R3 Acoustic Wave Snapshots (peak wavelet power at t = {t_peak_hr:.2f} hrs)", y=1.02)
plt.tight_layout()
plt.savefig("output/tier2_R3/snapshot_wavelet_peak.png", dpi=150, bbox_inches='tight')
print(f"\nSaved to output/tier2_R3/snapshot_wavelet_peak.png")

# --- Step 4: Also plot the power-at-forced-freq time series ---
fig2, ax2 = plt.subplots(figsize=(10, 4))
t_hr = times / 3600.0
ax2.plot(t_hr, power_at_forced / np.max(power_at_forced), 'k-', linewidth=0.8)
ax2.axvline(t_peak_hr, color='red', linestyle='--', label=f'Peak: t = {t_peak_hr:.2f} hrs')
ax2.set_xlabel("Time (hrs)")
ax2.set_ylabel("Normalized power at 3.18 mHz")
ax2.set_title("Wavelet Power at Forced Frequency vs Time")
ax2.legend()
plt.tight_layout()
plt.savefig("output/tier2_R3/wavelet_power_at_forced_freq.png", dpi=150)
print("Saved to output/tier2_R3/wavelet_power_at_forced_freq.png")
