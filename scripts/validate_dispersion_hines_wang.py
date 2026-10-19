"""
Dispersion relation validation of the KD-120km gravity wave simulation.

Compares simulation results against the Hines (1960) acoustic-gravity wave
dispersion relation and Wang et al. (2024) Mars wave parameters.

Generates:
  Fig 1: omega-kx dispersion diagram with simulation point
  Fig 2: kz-kx diagram with constant-omega contours
  Fig 3: Amplitude growth verification (w vs z, density-scaled)
  Fig 4: 2D FFT of w field to measure kx and kz from simulation data

Output: output/dispersion_validation/
"""

import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from scipy.signal import hilbert

# ---------------------------------------------------------------------------
# Plot style
# ---------------------------------------------------------------------------
plt.rcParams.update({
    "font.size": 12,
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 10,
    "figure.titlesize": 16,
})

# ---------------------------------------------------------------------------
# Mars atmosphere parameters
# ---------------------------------------------------------------------------
GAMMA = 1.31
G = 3.72           # m/s^2
R_GAS = 8.314462 / 0.04401  # J/(kg K), specific gas constant for CO2

# Two parameter sets: isothermal reference and actual MCS profile
# Isothermal reference (T0=210 K, used in config)
T0_ISO = 210.0
CS_ISO = np.sqrt(GAMMA * R_GAS * T0_ISO)
N_ISO = np.sqrt(GAMMA - 1.0) * G / CS_ISO
OMAC_ISO = GAMMA * G / (2.0 * CS_ISO)

# MCS profile bulk values (T~161 K above 40 km where wave propagates)
T0_MCS = 160.6
CS_MCS = np.sqrt(GAMMA * R_GAS * T0_MCS)
N_MCS = np.sqrt(GAMMA - 1.0) * G / CS_MCS
OMAC_MCS = GAMMA * G / (2.0 * CS_MCS)

# Use the MCS values as primary (matches actual simulation atmosphere)
CS = CS_MCS
N_BV = N_MCS
OMEGA_AC = OMAC_MCS
H = R_GAS * T0_MCS / G

# Simulation wave parameters
OMEGA_SIM = 0.005                  # rad/s
LAMBDA_X_SIM = 200e3               # m
KX_SIM = 2.0 * np.pi / LAMBDA_X_SIM

# Wang et al. 2024 GW test case
OMEGA_WANG = 2.0 * np.pi / (39.4 * 60.0)   # rad/s (period 39.4 min)
KX_WANG = 2.0 * np.pi / 210e3              # rad/m (lambda_h = 210 km)

# Data paths
RUN_DIR = "output/tier2_R2_kd_120km"
OUT_DIR = "output/dispersion_validation"


def kz_squared(omega, kx, cs=CS, N=N_BV, omega_ac=OMEGA_AC):
    """Hines (1960) Eq. 21' solved for kz^2.

    kz^2 = (omega^2 - omega_ac^2) / cs^2 + kx^2 * (N^2 / omega^2 - 1)
    """
    return (omega**2 - omega_ac**2) / cs**2 + kx**2 * (N**2 / omega**2 - 1.0)


def omega_branches(kx_arr, cs=CS, N=N_BV, omega_ac=OMEGA_AC):
    """Solve dispersion relation for omega given kx and kz=0 (boundary curves).

    Setting kz=0 gives the bounding curves of the propagation regions:
      omega^4 - omega^2 * (cs^2 * kx^2 + omega_ac^2) + N^2 * cs^2 * kx^2 = 0

    Returns (omega_gravity, omega_acoustic) arrays.
    """
    b = -(cs**2 * kx_arr**2 + omega_ac**2)
    c = N**2 * cs**2 * kx_arr**2

    disc = b**2 - 4.0 * c
    disc = np.maximum(disc, 0.0)

    omega2_low = (-b - np.sqrt(disc)) / 2.0
    omega2_high = (-b + np.sqrt(disc)) / 2.0

    return np.sqrt(np.maximum(omega2_low, 0.0)), np.sqrt(np.maximum(omega2_high, 0.0))


# ---------------------------------------------------------------------------
# Figure 1: omega-kx dispersion diagram
# ---------------------------------------------------------------------------
def plot_omega_kx(outdir):
    """Dispersion diagram in omega vs kx space (Hines 1960 style)."""
    fig, ax = plt.subplots(figsize=(9, 7))

    kx_arr = np.linspace(1e-7, 1.2e-4, 500)
    kx_km = kx_arr * 1e3  # rad/km for plotting

    # Boundary curves (kz=0) for MCS atmosphere
    omega_grav_0, omega_acou_0 = omega_branches(kx_arr)

    # Constant-kz curves
    kz_vals = [1e-5, 2e-5, 4e-5, 8e-5]
    colors_kz = ["#aaaaaa", "#888888", "#666666", "#444444"]

    for kz_val, col in zip(kz_vals, colors_kz):
        b_coeff = -(CS**2 * (kx_arr**2 + kz_val**2) + OMEGA_AC**2)
        c_coeff = N_BV**2 * CS**2 * kx_arr**2
        disc = b_coeff**2 - 4.0 * c_coeff
        disc = np.maximum(disc, 0.0)
        o2_low = (-b_coeff - np.sqrt(disc)) / 2.0
        o2_high = (-b_coeff + np.sqrt(disc)) / 2.0

        o_grav = np.sqrt(np.maximum(o2_low, 0.0))
        o_acou = np.sqrt(np.maximum(o2_high, 0.0))

        lz_km = 2 * np.pi / kz_val / 1e3
        ax.plot(kx_km, o_grav * 1e3, color=col, lw=0.8, ls="--")
        ax.plot(kx_km, o_acou * 1e3, color=col, lw=0.8, ls="--",
                label=r"$\lambda_z$ = " + f"{lz_km:.0f} km")

    # Boundary curves (kz=0)
    ax.plot(kx_km, omega_grav_0 * 1e3, "b-", lw=2.0, label="Gravity branch ($k_z=0$)")
    ax.plot(kx_km, omega_acou_0 * 1e3, "r-", lw=2.0, label="Acoustic branch ($k_z=0$)")

    # Asymptotic frequency limits
    ax.axhline(N_BV * 1e3, color="blue", ls=":", lw=1.0, alpha=0.7)
    ax.axhline(OMEGA_AC * 1e3, color="red", ls=":", lw=1.0, alpha=0.7)
    xr = kx_km[-1]
    ax.text(xr * 1.02, N_BV * 1e3, f"$N$ = {N_BV*1e3:.2f} mrad/s",
            va="center", fontsize=9, color="blue")
    ax.text(xr * 1.02, OMEGA_AC * 1e3, r"$\omega_{ac}$" + f" = {OMEGA_AC*1e3:.2f} mrad/s",
            va="center", fontsize=9, color="red")

    # Evanescent band
    ax.axhspan(N_BV * 1e3, OMEGA_AC * 1e3, alpha=0.15, color="gray",
               label="Evanescent band")

    # Simulation point
    ax.plot(KX_SIM * 1e3, OMEGA_SIM * 1e3, "k*", ms=15, zorder=10,
            label=f"This work ($\\tau$={2*np.pi/OMEGA_SIM/60:.0f} min, "
                  f"$\\lambda_x$={LAMBDA_X_SIM/1e3:.0f} km)")

    # Wang et al. 2024 GW point
    ax.plot(KX_WANG * 1e3, OMEGA_WANG * 1e3, "gD", ms=9, zorder=10,
            label=f"Wang+ 2024 GW ($\\tau$=39.4 min, $\\lambda_h$=210 km)")

    ax.set_xlabel(r"$k_x$ (mrad/m)")
    ax.set_ylabel(r"$\omega$ (mrad/s)")
    ax.set_title("Acoustic-Gravity Wave Dispersion (Mars CO$_2$, MCS profile T$\\approx$161 K)")
    ax.set_xlim(0, xr)
    ax.set_ylim(0, 35)
    ax.legend(loc="upper left", fontsize=8, framealpha=0.9)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig1_omega_kx_dispersion.png"), dpi=200)
    plt.close(fig)
    print("  Saved fig1_omega_kx_dispersion.png")


# ---------------------------------------------------------------------------
# Figure 2: kz-kx diagram with constant omega contours
# ---------------------------------------------------------------------------
def plot_kz_kx(outdir, kz_measured=None):
    """Hines Fig 9 style: contours of constant omega in kz vs kx space."""
    fig, ax = plt.subplots(figsize=(9, 7))

    kx_arr = np.linspace(1e-7, 1.5e-4, 500)

    # Gravity wave contours (omega < N)
    gw_omegas = [0.002, 0.003, 0.005, 0.007, 0.009]
    gw_colors = plt.cm.Blues(np.linspace(0.3, 0.9, len(gw_omegas)))

    for omega_val, col in zip(gw_omegas, gw_colors):
        kz2 = kz_squared(omega_val, kx_arr)
        mask = kz2 > 0
        if mask.any():
            kz_val = np.sqrt(kz2[mask])
            period_min = 2 * np.pi / omega_val / 60.0
            ax.plot(kx_arr[mask] * 1e3, kz_val * 1e3, color=col, lw=1.5,
                    label=f"$\\omega$={omega_val*1e3:.0f} mrad/s"
                          f" ($\\tau$={period_min:.0f} min)")

    # Acoustic wave contours (omega > omega_ac)
    aw_omegas = [0.013, 0.015, 0.020, 0.030]
    aw_colors = plt.cm.Reds(np.linspace(0.3, 0.9, len(aw_omegas)))

    for omega_val, col in zip(aw_omegas, aw_colors):
        kz2 = kz_squared(omega_val, kx_arr)
        mask = kz2 > 0
        if mask.any():
            kz_val = np.sqrt(kz2[mask])
            period_min = 2 * np.pi / omega_val / 60.0
            ax.plot(kx_arr[mask] * 1e3, kz_val * 1e3, color=col, lw=1.5, ls="--",
                    label=f"$\\omega$={omega_val*1e3:.0f} mrad/s"
                          f" ($\\tau$={period_min:.1f} min)")

    # Theoretical point for our simulation
    kz2_sim = kz_squared(OMEGA_SIM, KX_SIM)
    if kz2_sim > 0:
        kz_theory = np.sqrt(kz2_sim)
        ax.plot(KX_SIM * 1e3, kz_theory * 1e3, "k*", ms=15, zorder=10,
                label=f"Theory: $\\lambda_z$={2*np.pi/kz_theory/1e3:.0f} km")

    # Measured kz from 2D FFT (if available)
    if kz_measured is not None and kz_measured > 0:
        ax.plot(KX_SIM * 1e3, kz_measured * 1e3, "rs", ms=10, zorder=10,
                label=f"Measured: $\\lambda_z$={2*np.pi/kz_measured/1e3:.0f} km")

    ax.set_xlabel(r"$k_x$ (mrad/m)")
    ax.set_ylabel(r"$k_z$ (mrad/m)")
    ax.set_title("Wavenumber Diagram (Mars CO$_2$, MCS profile)")
    ax.legend(loc="upper left", fontsize=8, framealpha=0.9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 0.15)
    ax.set_ylim(0, 0.15)

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig2_kz_kx_diagram.png"), dpi=200)
    plt.close(fig)
    print("  Saved fig2_kz_kx_diagram.png")


# ---------------------------------------------------------------------------
# Load snapshot data
# ---------------------------------------------------------------------------
def load_late_snapshot(run_dir):
    """Load a late-time snapshot (near 80% of simulation)."""
    files = sorted(glob.glob(os.path.join(run_dir, "snapshot_*.nc")))
    if not files:
        raise FileNotFoundError(f"No snapshots in {run_dir}")

    idx = int(0.8 * len(files))
    fpath = files[idx]

    with Dataset(fpath, "r") as ds:
        x = ds.variables["x"][:].astype(np.float64)
        z = ds.variables["z"][:].astype(np.float64)
        w = ds.variables["w"][:, :].astype(np.float64)
        rho = ds.variables["rho"][:, :].astype(np.float64)
        t = float(ds.getncattr("time_s"))

    print(f"  Loaded snapshot {os.path.basename(fpath)} at t={t:.0f} s ({t/3600:.2f} hrs)")
    return x, z, w, rho, t


# ---------------------------------------------------------------------------
# Figure 3: Amplitude growth verification
# ---------------------------------------------------------------------------
def plot_amplitude_growth(x, z, w, rho, t, outdir):
    """Compare w amplitude growth against rho^(-1/2) prediction."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 6))

    ix_mid = len(x) // 2
    z_km = z / 1e3

    w_profile = w[:, ix_mid]
    rho_profile = rho[:, ix_mid]

    # Use actual background density (horizontal mean, removes wave perturbation)
    rho_bg = np.mean(rho, axis=1)

    # Density-scaled field using actual background density
    w_scaled = w_profile * np.sqrt(rho_bg / rho_bg[0])

    # Hilbert envelope of raw w
    analytic = hilbert(w_profile)
    envelope = np.abs(analytic)

    # Reference altitude: 20 km (above the temperature transition region)
    iz_ref = np.argmin(np.abs(z - 20e3))
    iz_sponge = np.argmin(np.abs(z - 110e3))
    iz_100 = np.argmin(np.abs(z - 100e3))

    # Theoretical envelope: w ~ rho^(-1/2)
    w_ref = envelope[iz_ref]
    theo_envelope = w_ref * np.sqrt(rho_bg[iz_ref] / rho_bg)

    # Panel A: raw w(z) with envelope and theory
    ax = axes[0]
    ax.plot(w_profile, z_km, "b-", lw=0.6, alpha=0.6, label="$w(z)$")
    ax.plot(envelope, z_km, "r-", lw=1.5, label="Hilbert envelope")
    ax.plot(-envelope, z_km, "r-", lw=1.5)
    mask = slice(iz_ref, iz_sponge)
    ax.plot(theo_envelope[mask], z_km[mask], "k--", lw=1.5,
            label=r"$\propto \rho^{-1/2}$")
    ax.plot(-theo_envelope[mask], z_km[mask], "k--", lw=1.5)
    ax.axhline(110, color="gray", ls=":", alpha=0.5)
    ax.text(0.02, 0.88, "Sponge", transform=ax.transAxes, fontsize=8,
            color="gray", rotation=0)
    ax.set_ylabel("Altitude (km)")
    ax.set_xlabel("$w$ (m/s)")
    ax.set_title("Vertical velocity profile")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(True, alpha=0.3)

    # Panel B: density-scaled w
    ax = axes[1]
    ax.plot(w_scaled, z_km, "b-", lw=0.8)
    ax.set_xlabel(r"$w \cdot \sqrt{\rho/\rho_0}$ (m/s)")
    ax.set_ylabel("Altitude (km)")
    ax.set_title("Density-scaled $w$")
    ax.axhline(110, color="gray", ls=":", alpha=0.5)
    ax.grid(True, alpha=0.3)

    # Compute CV in clean region (20-100 km)
    w_scaled_clean = w_scaled[iz_ref:iz_100]
    # Use Hilbert envelope of scaled field for more robust CV
    analytic_scaled = hilbert(w_scaled)
    env_scaled = np.abs(analytic_scaled)
    env_clean = env_scaled[iz_ref:iz_100]
    cv = np.std(env_clean) / np.mean(env_clean) if np.mean(env_clean) > 0 else 0

    ax.text(0.05, 0.95, f"Envelope CV (20-100 km)\n= {cv:.3f}",
            transform=ax.transAxes, fontsize=9, va="top",
            bbox=dict(boxstyle="round", fc="wheat", alpha=0.8))

    # Panel C: normalized growth factor
    ax = axes[2]
    norm_envelope = envelope / envelope[iz_ref]
    norm_theory = np.sqrt(rho_bg[iz_ref] / rho_bg) / np.sqrt(rho_bg[iz_ref] / rho_bg[iz_ref])
    ax.plot(norm_envelope[iz_ref:iz_sponge], z_km[iz_ref:iz_sponge],
            "b-", lw=1.5, label="Measured envelope")
    ax.plot(norm_theory[iz_ref:iz_sponge], z_km[iz_ref:iz_sponge],
            "k--", lw=1.5, label=r"$\sqrt{\rho_0/\rho(z)}$")
    ax.set_xlabel("Normalized amplitude")
    ax.set_ylabel("Altitude (km)")
    ax.set_title("Amplitude growth factor")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    fig.suptitle(f"Amplitude Growth Verification (t = {t/3600:.1f} hrs)", fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig3_amplitude_growth.png"), dpi=200,
                bbox_inches="tight")
    plt.close(fig)
    print("  Saved fig3_amplitude_growth.png")

    return cv


# ---------------------------------------------------------------------------
# Figure 4: Phase-tilt measurement of kz
# ---------------------------------------------------------------------------
def plot_phase_tilt(x, z, w, rho, t, outdir):
    """Measure kz by tracking the horizontal phase of w(x) at each altitude.

    For a wave w ~ cos(kx*x - kz*z - omega*t), the phase at fixed t is
    phi(z) = kz*z + const. So dφ/dz = kz (upward propagation, sign depends
    on convention). We extract the phase of the dominant Fourier mode at
    each altitude and fit a line to get kz.
    """
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.5))

    dx = x[1] - x[0]
    dz = z[1] - z[0]
    nx_pts = len(x)

    # Theoretical kz
    kz2_theory = kz_squared(OMEGA_SIM, KX_SIM)
    kz_theory = np.sqrt(max(kz2_theory, 0))

    # At each altitude, compute FFT along x and extract phase at kx = KX_SIM
    # kx_fft modes: mode index m corresponds to kx = 2*pi*m/(nx*dx)
    # Our wave has kx = 2*pi/lambda_x, so mode index m = nx*dx/lambda_x = 1
    # (since domain width = lambda_x = 200 km and nx*dx = 200 km)
    target_mode = 1  # fundamental horizontal mode

    phases = np.zeros(len(z))
    amplitudes = np.zeros(len(z))

    for j in range(len(z)):
        fft_row = np.fft.fft(w[j, :])
        phases[j] = np.angle(fft_row[target_mode])
        amplitudes[j] = np.abs(fft_row[target_mode]) * 2.0 / nx_pts

    # Unwrap phase for continuous measurement
    phases_unwrapped = np.unwrap(phases)

    # Fit linear slope in clean region (30-100 km, avoid boundaries and sponge)
    iz_fit_lo = np.argmin(np.abs(z - 30e3))
    iz_fit_hi = np.argmin(np.abs(z - 100e3))
    z_fit = z[iz_fit_lo:iz_fit_hi]
    phi_fit = phases_unwrapped[iz_fit_lo:iz_fit_hi]

    # Linear fit: phi(z) = kz * z + phi_0
    coeffs = np.polyfit(z_fit, phi_fit, 1)
    kz_measured = coeffs[0]  # rad/m (slope = dφ/dz = kz)
    phi_0 = coeffs[1]

    lambda_z_meas = 2 * np.pi / abs(kz_measured) if abs(kz_measured) > 0 else np.inf
    lambda_z_theory = 2 * np.pi / kz_theory if kz_theory > 0 else np.inf

    # Residuals from linear fit
    phi_fitted = np.polyval(coeffs, z_fit)
    residuals = phi_fit - phi_fitted
    rms_residual = np.sqrt(np.mean(residuals**2))

    # Panel A: w(x,z) field with phase lines
    ax = axes[0]
    z_km = z / 1e3
    x_km = x / 1e3
    rho_bg = np.mean(rho, axis=1)
    w_scaled_2d = w * np.sqrt(rho_bg[:, np.newaxis] / rho_bg[0])

    vlim = np.percentile(np.abs(w_scaled_2d[10:, :]), 99)
    im = ax.pcolormesh(x_km, z_km, w_scaled_2d, cmap="RdBu_r",
                       vmin=-vlim, vmax=vlim, shading="auto")
    plt.colorbar(im, ax=ax, label=r"$w\sqrt{\rho/\rho_0}$ (m/s)", shrink=0.8)

    # Overlay constant-phase lines (where phase = n*pi for integer n)
    for n in range(-10, 10):
        # phi(z) = kz*z + kx*x + const => x = (n*pi - kz*z - const) / kx
        # From the fit: phase at mode 1 = kz_meas*z + phi_0
        # A crest is where the total phase is 0 (or 2*n*pi)
        x_crest = (2 * n * np.pi - kz_measured * z - phi_0) / KX_SIM
        x_crest_km = x_crest / 1e3
        valid = (x_crest_km > 0) & (x_crest_km < 200)
        if valid.any():
            ax.plot(x_crest_km[valid], z_km[valid], "k-", lw=0.8, alpha=0.5)

    # Also show theoretical phase lines
    for n in range(-10, 10):
        x_crest_t = (2 * n * np.pi - kz_theory * z - phi_0) / KX_SIM
        x_crest_km_t = x_crest_t / 1e3
        valid = (x_crest_km_t > 0) & (x_crest_km_t < 200)
        if valid.any():
            ax.plot(x_crest_km_t[valid], z_km[valid], "g--", lw=0.8, alpha=0.5)

    ax.set_xlabel("$x$ (km)")
    ax.set_ylabel("Altitude (km)")
    ax.set_title(r"$w\sqrt{\rho/\rho_0}$ with phase lines")
    ax.set_ylim(0, 120)
    # Legend for phase lines
    ax.plot([], [], "k-", lw=0.8, label="Measured phase")
    ax.plot([], [], "g--", lw=0.8, label="Theory phase")
    ax.legend(fontsize=8, loc="upper right")

    # Panel B: phase vs altitude
    ax = axes[1]
    ax.plot(phases_unwrapped, z_km, "b-", lw=1.2, label="Measured phase")
    ax.plot(phi_fitted, z_fit / 1e3, "r--", lw=2.0,
            label=f"Linear fit (30-100 km)")
    # Theoretical slope
    phi_theory_line = kz_theory * z_fit + phi_0
    ax.plot(phi_theory_line, z_fit / 1e3, "k:", lw=1.5,
            label=f"Theory slope")
    ax.axhline(30, color="green", ls=":", alpha=0.4)
    ax.axhline(100, color="green", ls=":", alpha=0.4)
    ax.axhline(110, color="gray", ls=":", alpha=0.4)
    ax.set_xlabel("Phase (rad)")
    ax.set_ylabel("Altitude (km)")
    ax.set_title("Phase of dominant $k_x$ mode vs. altitude")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Annotation box
    ax.text(0.05, 0.95,
            f"$k_z$ (measured) = {kz_measured*1e3:.3f} mrad/m\n"
            f"$\\lambda_z$ (measured) = {lambda_z_meas/1e3:.0f} km\n"
            f"$k_z$ (theory) = {kz_theory*1e3:.3f} mrad/m\n"
            f"$\\lambda_z$ (theory) = {lambda_z_theory/1e3:.0f} km\n"
            f"RMS residual = {rms_residual:.3f} rad",
            transform=ax.transAxes, fontsize=8, va="top",
            bbox=dict(boxstyle="round", fc="lightyellow", alpha=0.9))

    # Panel C: amplitude of dominant mode vs altitude
    ax = axes[2]
    ax.plot(amplitudes, z_km, "b-", lw=1.2, label="FFT amplitude at $k_x$")

    # Theoretical growth: A ~ rho^(-1/2)
    iz_ref = np.argmin(np.abs(z - 30e3))
    a_ref = amplitudes[iz_ref]
    theo_amp = a_ref * np.sqrt(rho_bg[iz_ref] / rho_bg)
    iz_sponge = np.argmin(np.abs(z - 110e3))
    ax.plot(theo_amp[iz_ref:iz_sponge], z_km[iz_ref:iz_sponge],
            "k--", lw=1.5, label=r"$\propto \rho^{-1/2}$")

    ax.axhline(110, color="gray", ls=":", alpha=0.5)
    ax.set_xlabel("Amplitude (m/s)")
    ax.set_ylabel("Altitude (km)")
    ax.set_title("Amplitude of $k_x$ mode vs. altitude")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    fig.suptitle(f"Phase-Tilt Analysis (t = {t/3600:.1f} hrs)", fontsize=14)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig4_phase_tilt.png"), dpi=200)
    plt.close(fig)

    print(f"  Saved fig4_phase_tilt.png")
    print(f"  Phase-tilt kz = {kz_measured:.4e} rad/m (lambda_z = {lambda_z_meas/1e3:.0f} km)")
    print(f"  Theory    kz = {kz_theory:.4e} rad/m (lambda_z = {lambda_z_theory/1e3:.0f} km)")
    print(f"  Ratio measured/theory = {abs(kz_measured)/kz_theory:.3f}")
    print(f"  RMS phase residual = {rms_residual:.3f} rad")

    return abs(kz_measured), kz_theory


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
def print_summary(kz_measured, kz_theory, cv_scaled):
    """Print comparison table."""
    lz_m = 2 * np.pi / kz_measured / 1e3 if kz_measured > 0 else np.inf
    lz_t = 2 * np.pi / kz_theory / 1e3 if kz_theory > 0 else np.inf

    print("\n" + "=" * 72)
    print("DISPERSION RELATION VALIDATION SUMMARY")
    print("=" * 72)
    print(f"\nMars MCS atmosphere (T~{T0_MCS} K above 40 km, gamma={GAMMA}, g={G} m/s^2)")
    print(f"  cs        = {CS:.2f} m/s")
    print(f"  H         = {H/1e3:.2f} km")
    print(f"  N         = {N_BV:.6f} rad/s  (period = {2*np.pi/N_BV/60:.1f} min)")
    print(f"  omega_ac  = {OMEGA_AC:.6f} rad/s  (period = {2*np.pi/OMEGA_AC/60:.1f} min)")
    print(f"  Ordering  : omega_sim ({OMEGA_SIM}) < N ({N_BV:.4f}) < omega_ac ({OMEGA_AC:.4f})")
    print(f"  Wave regime: GRAVITY WAVE")
    print(f"\nSimulation wave:")
    print(f"  omega     = {OMEGA_SIM} rad/s  (period = {2*np.pi/OMEGA_SIM/60:.1f} min)")
    print(f"  kx        = {KX_SIM:.4e} rad/m  (lambda_x = {LAMBDA_X_SIM/1e3:.0f} km)")
    print(f"\n{'Quantity':<35} {'Theory':<20} {'Measured':<20}")
    print("-" * 72)
    print(f"{'kz (rad/m)':<35} {kz_theory:.4e}{'':>8} {kz_measured:.4e}")
    print(f"{'lambda_z (km)':<35} {lz_t:<20.0f} {lz_m:<20.0f}")
    print(f"{'Envelope CV (20-100 km)':<35} {'~0 (constant)':<20} {cv_scaled:<20.3f}")

    cp_x = OMEGA_SIM / KX_SIM
    cp_z_t = OMEGA_SIM / kz_theory if kz_theory > 0 else np.inf
    cp_z_m = OMEGA_SIM / kz_measured if kz_measured > 0 else np.inf
    print(f"{'c_px (m/s)':<35} {cp_x:<20.1f} {'(imposed)':<20}")
    print(f"{'c_pz (m/s)':<35} {cp_z_t:<20.1f} {cp_z_m:<20.1f}")

    # Note about method
    ratio = abs(kz_measured) / kz_theory if kz_theory > 0 else np.inf
    print(f"\nkz ratio (measured/theory) = {ratio:.3f}")
    print(f"Note: kz measured via phase-tilt of dominant kx Fourier mode.")
    print(f"Theoretical lambda_z ({lz_t:.0f} km) exceeds the domain height")
    print(f"(120 km), so this is a partial-wavelength measurement.")
    print("=" * 72)

    return {
        "kz_theory": kz_theory, "kz_measured": kz_measured,
        "lambda_z_theory_km": lz_t, "lambda_z_measured_km": lz_m,
        "cv_scaled": cv_scaled, "cp_x": cp_x,
        "cp_z_theory": cp_z_t, "cp_z_measured": cp_z_m,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    print("Dispersion Relation Validation (Hines 1960 / Wang et al. 2024)")
    print("=" * 60)
    print(f"\nUsing MCS atmosphere parameters (T~{T0_MCS} K above 40 km)")
    print(f"  cs={CS:.1f} m/s, N={N_BV:.5f} rad/s, omega_ac={OMEGA_AC:.5f} rad/s\n")

    print("[1/4] omega-kx dispersion diagram...")
    plot_omega_kx(OUT_DIR)

    print("\n[2/4] Loading late-time snapshot...")
    x, z, w, rho, t = load_late_snapshot(RUN_DIR)

    print("\n[3/4] Amplitude growth verification...")
    cv_scaled = plot_amplitude_growth(x, z, w, rho, t, OUT_DIR)

    print("\n[4/4] Phase-tilt kz measurement...")
    kz_measured, kz_theory = plot_phase_tilt(x, z, w, rho, t, OUT_DIR)

    # Generate kz-kx diagram with measured point
    print("\n[2b] kz-kx diagram with measured kz...")
    plot_kz_kx(OUT_DIR, kz_measured=kz_measured)

    summary = print_summary(kz_measured, kz_theory, cv_scaled)
    return summary


if __name__ == "__main__":
    main()
