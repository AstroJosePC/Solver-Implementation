"""
Cavity interference analysis for R2 (gravity wave) and R3 (acoustic wave) simulations.

Models the amplitude envelope at the DSMC handoff altitude (z = 100 km) as the
superposition of the direct forced wave and N reflected copies bouncing between
the bottom boundary and the sponge layer at 170 km.

E(t) = A0 * tanh(t/t0) * |1 + sum_{n=1}^{N(t)} r^n * exp(i*n*phi_rt)|

Produces four figures saved to output/cavity_analysis/.
"""

import sys
import time as timer
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import minimize
from scipy.signal import hilbert

sys.path.insert(0, "scripts")
from wavelet_tc98 import wavelet

import netCDF4 as nc


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

@dataclass
class CavityParams:
    """Parameters for a cavity interference model."""
    name: str
    omega: float              # angular frequency of forced wave (rad/s)
    kz: float                 # vertical wavenumber (rad/m)
    L_cavity: float           # cavity depth in meters (bottom to sponge)
    r_init: float             # per-pass reflection coefficient (from sponge sweep)
    t0_ramp: float            # ramp timescale in seconds
    data_path: str            # path to dsmc_handoff_100km.nc
    c_gz_theory: float        # theoretical group velocity (m/s)
    # Derived / filled after loading
    c_gz_empirical: Optional[float] = None
    T_roundtrip: Optional[float] = None

    @property
    def lambda_z(self) -> float:
        return 2 * np.pi / self.kz

    @property
    def phi_rt_theory(self) -> float:
        """Theoretical round-trip phase (rad)."""
        return 2 * self.kz * self.L_cavity


R3_PARAMS = CavityParams(
    name="R3 (Acoustic)",
    omega=0.02,
    kz=6.83e-5,
    L_cavity=170e3,
    r_init=0.014,
    t0_ramp=1200.0,
    data_path="output/tier2_R3/dsmc_handoff_100km.nc",
    c_gz_theory=181.0,
)

R2_PARAMS = CavityParams(
    name="R2 (Gravity)",
    omega=0.005,
    kz=5.67e-5,
    L_cavity=170e3,
    r_init=0.014,
    t0_ramp=5000.0,
    data_path="output/tier2_R2/dsmc_handoff_100km.nc",
    c_gz_theory=95.0,
)

OUTPUT_DIR = Path("output/cavity_analysis")


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_handoff_data(path: str) -> tuple[np.ndarray, np.ndarray]:
    """Load handoff data, return (times_s, w_mid) at the midpoint horizontal cell."""
    ds = nc.Dataset(path, "r")
    times = ds.variables["time"][:]
    w = ds.variables["w"][:]
    x = ds.variables["x"][:]
    ix_mid = len(x) // 2
    w_mid = w[ix_mid, :]
    ds.close()
    return times.astype(np.float64), w_mid.astype(np.float64)


# ---------------------------------------------------------------------------
# Empirical group velocity estimation
# ---------------------------------------------------------------------------

def estimate_cgz_empirical(
    times: np.ndarray,
    w: np.ndarray,
    z_handoff: float = 100e3,
    frac_threshold: float = 0.10,
    smooth_window: int = 50,
) -> float:
    """
    Estimate vertical group velocity from the first-arrival time of wave
    energy at the handoff altitude.

    Uses a running-max envelope of |w|, then finds the first time it exceeds
    frac_threshold * peak_envelope. The 10% threshold avoids triggering on
    early numerical noise while capturing the physical wave front.
    """
    from scipy.ndimage import maximum_filter1d

    env = maximum_filter1d(np.abs(w), size=smooth_window)
    peak = np.max(env)

    if peak == 0:
        raise ValueError("Signal is identically zero")

    threshold = frac_threshold * peak
    arrival_mask = env > threshold

    if not np.any(arrival_mask):
        raise ValueError("No signal exceeding threshold found")

    i_arrival = np.argmax(arrival_mask)
    t_arrival = times[i_arrival]

    if t_arrival <= 0:
        raise ValueError(f"Arrival time is non-positive: {t_arrival}")

    c_gz = z_handoff / t_arrival
    return c_gz


# ---------------------------------------------------------------------------
# Hilbert envelope
# ---------------------------------------------------------------------------

def compute_hilbert_envelope(
    w: np.ndarray,
    times: np.ndarray,
    omega: float,
    n_smooth_periods: float = 2.0,
) -> np.ndarray:
    """
    Compute the smoothed amplitude envelope via the Hilbert transform.

    Applies Gaussian smoothing over n_smooth_periods wave periods to suppress
    intra-cycle oscillation in the envelope.
    """
    analytic = hilbert(w)
    envelope = np.abs(analytic)

    # Smooth over n_smooth_periods of the forced wave
    period_s = 2 * np.pi / omega
    dt = np.median(np.diff(times))
    sigma_samples = (n_smooth_periods * period_s) / dt
    envelope_smooth = gaussian_filter1d(envelope, sigma=sigma_samples)

    return envelope_smooth


# ---------------------------------------------------------------------------
# Cavity model
# ---------------------------------------------------------------------------

def compute_cavity_model(
    times: np.ndarray,
    params: CavityParams,
    A0: float,
    phi_rt: float,
    r: float,
    gamma_decay: float = 0.0,
) -> np.ndarray:
    """
    Compute the cavity interference model envelope.

    E(t) = A0 * tanh(t / t0) * exp(-gamma * t) * |1 + sum_{n=1}^{N(t)} r^n * exp(i*n*phi_rt)|

    The exponential decay term captures cumulative energy loss through the
    sponge layer under continuous forcing, beyond what individual reflections
    account for. Each reflection n arrives after n round-trips.
    """
    T_rt = params.T_roundtrip
    t0 = params.t0_ramp

    ramp = np.tanh(times / t0)
    decay = np.exp(-gamma_decay * times)

    # Maximum number of round trips possible at the end of the simulation
    N_max = int(np.floor(times[-1] / T_rt)) + 1

    # Build the complex interference sum incrementally
    interference = np.ones(len(times), dtype=complex)

    for n in range(1, N_max + 1):
        arrival_time = n * T_rt
        mask = times >= arrival_time
        interference[mask] += r**n * np.exp(1j * n * phi_rt)

    envelope = A0 * ramp * decay * np.abs(interference)
    return envelope


# ---------------------------------------------------------------------------
# Fitting
# ---------------------------------------------------------------------------

def fit_cavity_model(
    envelope: np.ndarray,
    times: np.ndarray,
    params: CavityParams,
    free_r: bool = False,
) -> dict:
    """
    Fit the cavity model to the observed Hilbert envelope.

    Free parameters: A0, phi_rt, gamma_decay, t0_eff (and optionally r).
    t0_eff is the effective buildup timescale at the handoff altitude, which
    can differ from the source ramp because amplitude grows as wave cycles
    accumulate over several wave periods.

    Uses differential_evolution for global optimization to avoid local minima.
    """
    from scipy.optimize import differential_evolution

    # Normalize envelope for numerical stability
    env_scale = np.max(envelope)
    if env_scale == 0:
        env_scale = 1.0
    env_norm = envelope / env_scale

    # Bounds: [A0_norm, phi_rt, gamma_decay, t0_eff_s, (r)]
    gamma_max = 5.0 / 3600.0  # max decay rate
    t_end = times[-1]
    bounds = [
        (0.5, 10.0),                     # A0 normalized
        (0.0, 8 * np.pi),                # phi_rt (0 to 4 full cycles)
        (0.0, gamma_max),                # gamma_decay (1/s)
        (params.t0_ramp * 0.5, t_end),   # t0_eff (s)
    ]
    if free_r:
        bounds.append((0.001, 0.5))  # r

    def cost(x):
        A0_n = x[0]
        phi_rt = x[1]
        gamma = x[2]
        t0_eff = x[3]
        r = x[4] if free_r else params.r_init

        # Override the ramp timescale with the fitted t0_eff
        ramp = np.tanh(times / t0_eff)
        decay = np.exp(-gamma * times)

        T_rt = params.T_roundtrip
        N_max = int(np.floor(times[-1] / T_rt)) + 1
        interference = np.ones(len(times), dtype=complex)
        for n in range(1, N_max + 1):
            mask = times >= n * T_rt
            interference[mask] += r**n * np.exp(1j * n * phi_rt)

        model_n = A0_n * ramp * decay * np.abs(interference)

        # Weights: down-weight early ramp (t < t0_eff)
        weights = np.ones_like(times)
        weights[times < t0_eff] = 0.3

        residual = weights * (model_n - env_norm) ** 2
        return np.sum(residual)

    # Global optimization
    result = differential_evolution(
        cost, bounds, seed=42, maxiter=500, tol=1e-8, polish=True,
    )

    A0_fit = result.x[0] * env_scale
    phi_fit = result.x[1]
    gamma_fit = result.x[2]
    t0_fit = result.x[3]
    r_fit = result.x[4] if free_r else params.r_init

    # Recompute final model
    model = compute_cavity_model(times, params, A0_fit, phi_fit, r_fit,
                                 gamma_decay=gamma_fit)
    # Correct for fitted t0: recompute with t0_fit instead of params.t0_ramp
    ramp = np.tanh(times / t0_fit)
    decay = np.exp(-gamma_fit * times)
    T_rt = params.T_roundtrip
    N_max = int(np.floor(times[-1] / T_rt)) + 1
    interference = np.ones(len(times), dtype=complex)
    for n in range(1, N_max + 1):
        mask = times >= n * T_rt
        interference[mask] += r_fit**n * np.exp(1j * n * phi_fit)
    model = A0_fit * ramp * decay * np.abs(interference)

    residual_rms = np.sqrt(np.mean((model - envelope) ** 2))
    residual_norm = residual_rms / np.max(envelope)

    return {
        "A0": A0_fit,
        "phi_rt": phi_fit,
        "gamma_decay": gamma_fit,
        "gamma_decay_hr": gamma_fit * 3600,
        "tau_decay_hr": 1 / gamma_fit / 3600 if gamma_fit > 0 else np.inf,
        "t0_eff": t0_fit,
        "t0_eff_hr": t0_fit / 3600,
        "r": r_fit,
        "residual_rms": residual_rms,
        "residual_norm": residual_norm,
        "model": model,
        "success": result.success,
    }


# ---------------------------------------------------------------------------
# Clean window computation
# ---------------------------------------------------------------------------

def compute_clean_window(params: CavityParams) -> dict:
    """
    Compute the clean window: time after ramp and before first reflection.

    Clean window = [3*t0, T_first_reflection_arrival]
    where T_first_reflection_arrival ~ T_roundtrip (first reflected wave
    returns to the handoff altitude).
    """
    t_ramp_done = 3 * params.t0_ramp
    T_rt = params.T_roundtrip

    # First reflection arrives at handoff altitude after one round trip
    t_first_refl = T_rt

    clean_start = t_ramp_done
    clean_end = t_first_refl

    if clean_end <= clean_start:
        # No clean window exists
        return {
            "exists": False,
            "start_s": clean_start,
            "end_s": clean_end,
            "duration_s": 0.0,
            "note": "First reflection arrives before ramp completes",
        }

    return {
        "exists": True,
        "start_s": clean_start,
        "end_s": clean_end,
        "duration_s": clean_end - clean_start,
    }


# ---------------------------------------------------------------------------
# Wavelet power at forced frequency
# ---------------------------------------------------------------------------

def compute_wavelet_power_at_forced(
    w_mid: np.ndarray,
    times: np.ndarray,
    omega: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute wavelet power at the forced frequency as a function of time."""
    dt_s = float(np.median(np.diff(times)))
    s0 = 2.0 * dt_s
    dj = 0.05
    j1 = int(np.fix(np.log2(len(w_mid) * dt_s / s0) / dj))

    wave, period_s, _, _ = wavelet(
        w_mid - np.mean(w_mid), dt_s,
        pad=1, dj=dj, s0=s0, J1=j1,
        mother="MORLET", param=6.0,
    )

    power = np.abs(wave) ** 2
    freq_mhz = 1000.0 / period_s
    freq_forced_mhz = (omega / (2 * np.pi)) * 1000.0
    i_freq = np.argmin(np.abs(freq_mhz - freq_forced_mhz))

    return times, power[i_freq, :]


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def plot_fig1_R3_envelope(times, w_mid, envelope, fit_result, params):
    """Figure 1: R3 Envelope Analysis (two panels)."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    t_hr = times / 3600.0

    # Top: raw w, Hilbert envelope, cavity model, round-trip markers
    ax1.plot(t_hr, w_mid, color="0.7", linewidth=0.3, label="w(t) raw", rasterized=True)
    ax1.plot(t_hr, envelope, "b-", linewidth=1.2, label="Hilbert envelope")
    ax1.plot(t_hr, -envelope, "b-", linewidth=1.2, alpha=0.5)
    ax1.plot(t_hr, fit_result["model"], "r--", linewidth=1.2,
             label=f'Cavity model (r={fit_result["r"]:.3f})')
    ax1.plot(t_hr, -fit_result["model"], "r--", linewidth=1.2, alpha=0.5)

    # Round-trip markers
    T_rt_hr = params.T_roundtrip / 3600.0
    n_rt = int(t_hr[-1] / T_rt_hr) + 1
    for n in range(1, n_rt + 1):
        ax1.axvline(n * T_rt_hr, color="gray", linestyle=":", alpha=0.4, linewidth=0.7)
    ax1.axvline(T_rt_hr, color="gray", linestyle=":", alpha=0.4, linewidth=0.7,
                label=f"T_rt = {params.T_roundtrip:.0f} s")

    ax1.set_ylabel("w (m/s)")
    ax1.set_title(f"{params.name}: Envelope at z = 100 km")
    ax1.legend(fontsize=8, loc="upper right")

    # Bottom: normalized envelope with clean window shading
    env_peak = np.max(envelope)
    ax2.plot(t_hr, envelope / env_peak, "b-", linewidth=1.2, label="Normalized envelope")
    ax2.plot(t_hr, fit_result["model"] / env_peak, "r--", linewidth=1.2,
             label="Cavity model (normalized)")

    cw = compute_clean_window(params)
    if cw["exists"]:
        ax2.axvspan(cw["start_s"] / 3600, cw["end_s"] / 3600,
                     color="green", alpha=0.15, label="Clean window")
    else:
        ax2.annotate("No clean window\n(reflection arrives before ramp ends)",
                     xy=(0.5, 0.5), xycoords="axes fraction", ha="center",
                     fontsize=9, color="red", alpha=0.7)

    # Annotate fit parameters
    tau_str = (f"{fit_result['tau_decay_hr']:.1f} hrs"
               if np.isfinite(fit_result['tau_decay_hr']) else "inf")
    text = (f"A0 = {fit_result['A0']:.3f} m/s\n"
            f"t0_eff = {fit_result['t0_eff_hr']:.2f} hrs\n"
            f"phi_rt = {fit_result['phi_rt']:.2f} rad "
            f"({fit_result['phi_rt']/(2*np.pi):.2f} cycles)\n"
            f"r = {fit_result['r']:.4f}\n"
            f"tau_decay = {tau_str}\n"
            f"RMS residual = {fit_result['residual_norm']:.1%}")
    ax2.text(0.02, 0.95, text, transform=ax2.transAxes, fontsize=8,
             verticalalignment="top", fontfamily="monospace",
             bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

    ax2.set_xlabel("Time (hrs)")
    ax2.set_ylabel("Normalized amplitude")
    ax2.legend(fontsize=8, loc="upper right")

    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "fig1_R3_envelope.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {OUTPUT_DIR / 'fig1_R3_envelope.png'}")


def plot_fig2_R2_envelope(times, w_mid, envelope, fit_result, params):
    """Figure 2: R2 Envelope Analysis (two panels)."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    t_hr = times / 3600.0

    ax1.plot(t_hr, w_mid, color="0.7", linewidth=0.3, label="w(t) raw", rasterized=True)
    ax1.plot(t_hr, envelope, "b-", linewidth=1.2, label="Hilbert envelope")
    ax1.plot(t_hr, -envelope, "b-", linewidth=1.2, alpha=0.5)
    ax1.plot(t_hr, fit_result["model"], "r--", linewidth=1.2,
             label=f'Cavity model (r={fit_result["r"]:.3f})')
    ax1.plot(t_hr, -fit_result["model"], "r--", linewidth=1.2, alpha=0.5)

    T_rt_hr = params.T_roundtrip / 3600.0
    n_rt = int(t_hr[-1] / T_rt_hr) + 1
    for n in range(1, n_rt + 1):
        ax1.axvline(n * T_rt_hr, color="gray", linestyle=":", alpha=0.4, linewidth=0.7)
    ax1.axvline(T_rt_hr, color="gray", linestyle=":", alpha=0.4, linewidth=0.7,
                label=f"T_rt = {params.T_roundtrip:.0f} s")

    ax1.set_ylabel("w (m/s)")
    ax1.set_title(f"{params.name}: Envelope at z = 100 km")
    ax1.legend(fontsize=8, loc="upper right")

    # Annotate N_total and T_rt
    ax1.annotate(f"N_total = {n_rt - 1}, T_rt = {params.T_roundtrip/60:.0f} min",
                 xy=(0.98, 0.05), xycoords="axes fraction", ha="right",
                 fontsize=9, color="gray")

    env_peak = np.max(envelope)
    ax2.plot(t_hr, envelope / env_peak, "b-", linewidth=1.2, label="Normalized envelope")
    ax2.plot(t_hr, fit_result["model"] / env_peak, "r--", linewidth=1.2,
             label="Cavity model (normalized)")

    cw = compute_clean_window(params)
    if cw["exists"]:
        ax2.axvspan(cw["start_s"] / 3600, cw["end_s"] / 3600,
                     color="green", alpha=0.15, label="Clean window")
    else:
        ax2.annotate("No clean window\n(reflection arrives before ramp ends)",
                     xy=(0.5, 0.5), xycoords="axes fraction", ha="center",
                     fontsize=9, color="red", alpha=0.7)

    tau_str = (f"{fit_result['tau_decay_hr']:.1f} hrs"
               if np.isfinite(fit_result['tau_decay_hr']) else "inf")
    text = (f"A0 = {fit_result['A0']:.3f} m/s\n"
            f"t0_eff = {fit_result['t0_eff_hr']:.2f} hrs\n"
            f"phi_rt = {fit_result['phi_rt']:.2f} rad "
            f"({fit_result['phi_rt']/(2*np.pi):.2f} cycles)\n"
            f"r = {fit_result['r']:.4f}\n"
            f"tau_decay = {tau_str}\n"
            f"RMS residual = {fit_result['residual_norm']:.1%}")
    ax2.text(0.02, 0.95, text, transform=ax2.transAxes, fontsize=8,
             verticalalignment="top", fontfamily="monospace",
             bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

    ax2.set_xlabel("Time (hrs)")
    ax2.set_ylabel("Normalized amplitude")
    ax2.legend(fontsize=8, loc="upper right")

    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "fig2_R2_envelope.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {OUTPUT_DIR / 'fig2_R2_envelope.png'}")


def plot_fig3_comparison(
    times_r3, env_r3, fit_r3, params_r3,
    times_r2, env_r2, fit_r2, params_r2,
):
    """Figure 3: R2 vs R3 comparison (two panels)."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

    # Left: normalized envelopes overlaid
    t_hr_r3 = times_r3 / 3600.0
    t_hr_r2 = times_r2 / 3600.0
    env_r3_norm = env_r3 / np.max(env_r3)
    env_r2_norm = env_r2 / np.max(env_r2)

    ax1.plot(t_hr_r3, env_r3_norm, "b-", linewidth=1.2, label="R3 (Acoustic)")
    ax1.plot(t_hr_r2, env_r2_norm, "r-", linewidth=1.2, label="R2 (Gravity)")
    ax1.set_xlabel("Time (hrs)")
    ax1.set_ylabel("Normalized envelope")
    ax1.set_title("Amplitude Envelopes")
    ax1.legend(fontsize=9)
    ax1.set_ylim(0, 1.1)

    # Right: parameter summary table
    ax2.axis("off")
    peak_end_r3 = np.max(env_r3) / env_r3[-1] if env_r3[-1] > 0 else np.inf
    peak_end_r2 = np.max(env_r2) / env_r2[-1] if env_r2[-1] > 0 else np.inf

    tau_r3 = (f"{fit_r3['tau_decay_hr']:.1f} hrs"
              if np.isfinite(fit_r3['tau_decay_hr']) else "inf")
    tau_r2 = (f"{fit_r2['tau_decay_hr']:.1f} hrs"
              if np.isfinite(fit_r2['tau_decay_hr']) else "inf")

    table_data = [
        ["Parameter", "R3 (Acoustic)", "R2 (Gravity)"],
        ["T_roundtrip", f"{params_r3.T_roundtrip/60:.0f} min",
         f"{params_r2.T_roundtrip/60:.0f} min"],
        ["N_total", f"{int(times_r3[-1] / params_r3.T_roundtrip)}",
         f"{int(times_r2[-1] / params_r2.T_roundtrip)}"],
        ["c_gz (empirical)", f"{params_r3.c_gz_empirical:.1f} m/s",
         f"{params_r2.c_gz_empirical:.1f} m/s"],
        ["phi_rt (fit)", f"{fit_r3['phi_rt']:.2f} rad",
         f"{fit_r2['phi_rt']:.2f} rad"],
        ["phi_rt (theory)", f"{params_r3.phi_rt_theory:.2f} rad",
         f"{params_r2.phi_rt_theory:.2f} rad"],
        ["r (fit)", f"{fit_r3['r']:.4f}", f"{fit_r2['r']:.4f}"],
        ["t0_eff", f"{fit_r3['t0_eff_hr']:.2f} hrs",
         f"{fit_r2['t0_eff_hr']:.2f} hrs"],
        ["tau_decay", tau_r3, tau_r2],
        ["Peak/End ratio", f"{peak_end_r3:.1f}x", f"{peak_end_r2:.1f}x"],
        ["Residual (norm)", f"{fit_r3['residual_norm']:.1%}",
         f"{fit_r2['residual_norm']:.1%}"],
    ]

    table = ax2.table(cellText=table_data, loc="center", cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.0, 1.5)

    # Style header row
    for j in range(3):
        table[0, j].set_facecolor("#4472C4")
        table[0, j].set_text_props(color="white", fontweight="bold")

    ax2.set_title("Cavity Parameters", fontsize=11, pad=15)

    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "fig3_comparison.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {OUTPUT_DIR / 'fig3_comparison.png'}")


def plot_fig4_handoff_viability(
    times_r3, w_mid_r3, env_r3, params_r3,
    times_r2, w_mid_r2, env_r2, params_r2,
):
    """Figure 4: Handoff viability (three panels)."""
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 9), sharex=False)

    # --- Helper for interference region shading ---
    def shade_regions(ax, times, envelope, params):
        t_hr = times / 3600.0
        env_peak = np.max(envelope)
        env_norm = envelope / env_peak

        ax.plot(t_hr, env_norm, "k-", linewidth=1.2)

        # Green: before first reflection (clean)
        t_refl = params.T_roundtrip / 3600.0
        ax.axvspan(0, t_refl, color="green", alpha=0.12, label="Clean (pre-reflection)")

        # Yellow: moderate interference (envelope > 50% of peak)
        # Red: substantial interference (envelope < 50% of peak after peak)
        i_peak = np.argmax(envelope)
        if i_peak < len(times) - 1:
            t_peak_hr = t_hr[i_peak]
            # After peak: find where envelope drops below 50%
            post_peak_norm = env_norm[i_peak:]
            i_50 = np.argmax(post_peak_norm < 0.5)
            if i_50 > 0:
                t_50_hr = t_hr[i_peak + i_50]
                ax.axvspan(t_refl, t_50_hr, color="gold", alpha=0.12,
                           label="Moderate interference")
                ax.axvspan(t_50_hr, t_hr[-1], color="red", alpha=0.10,
                           label="Substantial interference")
            else:
                ax.axvspan(t_refl, t_hr[-1], color="gold", alpha=0.12,
                           label="Moderate interference")

        ax.axhline(0.1, color="gray", linestyle="--", alpha=0.5, label="10% threshold")
        ax.set_ylabel("Normalized envelope")
        ax.legend(fontsize=7, loc="upper right")

    # Top: R3
    shade_regions(ax1, times_r3, env_r3, params_r3)
    ax1.set_title(f"{params_r3.name}: Handoff Viability")
    ax1.set_xlabel("Time (hrs)")

    # Middle: R2
    shade_regions(ax2, times_r2, env_r2, params_r2)
    ax2.set_title(f"{params_r2.name}: Handoff Viability")
    ax2.set_xlabel("Time (hrs)")

    # Bottom: Wavelet power at forced frequency for both
    print("  Computing wavelet for R3...")
    t_w3, pow_r3 = compute_wavelet_power_at_forced(w_mid_r3, times_r3, params_r3.omega)
    print("  Computing wavelet for R2...")
    t_w2, pow_r2 = compute_wavelet_power_at_forced(w_mid_r2, times_r2, params_r2.omega)

    pow_r3_norm = pow_r3 / np.max(pow_r3)
    pow_r2_norm = pow_r2 / np.max(pow_r2)

    ax3.plot(t_w3 / 3600, pow_r3_norm, "b-", linewidth=1.0, label="R3 (Acoustic)")
    ax3.plot(t_w2 / 3600, pow_r2_norm, "r-", linewidth=1.0, label="R2 (Gravity)")
    ax3.axhline(0.1, color="gray", linestyle="--", alpha=0.5, label="10% threshold")
    ax3.set_xlabel("Time (hrs)")
    ax3.set_ylabel("Normalized wavelet power\nat forced frequency")
    ax3.set_title("Wavelet Power at Forced Frequency")
    ax3.legend(fontsize=8)

    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "fig4_handoff_viability.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {OUTPUT_DIR / 'fig4_handoff_viability.png'}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def analyze_case(params: CavityParams) -> dict:
    """Run full cavity analysis for one case."""
    print(f"\n{'='*60}")
    print(f"Analyzing: {params.name}")
    print(f"{'='*60}")

    # Load data
    print("  Loading handoff data...")
    times, w_mid = load_handoff_data(params.data_path)
    dt = np.median(np.diff(times))
    print(f"  dt = {dt:.2f} s, N_t = {len(times)}, "
          f"duration = {times[-1]/3600:.2f} hrs")

    # Empirical group velocity
    print("  Estimating empirical c_gz...")
    try:
        c_gz_emp = estimate_cgz_empirical(times, w_mid)
    except ValueError as e:
        print(f"  WARNING: {e}. Falling back to theoretical c_gz.")
        c_gz_emp = params.c_gz_theory
    params.c_gz_empirical = c_gz_emp
    print(f"  c_gz (empirical) = {c_gz_emp:.1f} m/s "
          f"(theory: {params.c_gz_theory:.1f} m/s)")

    # Round-trip time from empirical c_gz
    params.T_roundtrip = 2 * params.L_cavity / c_gz_emp
    print(f"  T_roundtrip = {params.T_roundtrip:.0f} s "
          f"= {params.T_roundtrip/60:.1f} min")
    print(f"  N_total = {int(times[-1] / params.T_roundtrip)} round trips")

    # Hilbert envelope
    print("  Computing Hilbert envelope...")
    envelope = compute_hilbert_envelope(w_mid, times, params.omega)

    # Fit cavity model (first with r fixed)
    print("  Fitting cavity model (r fixed)...")
    fit_fixed = fit_cavity_model(envelope, times, params, free_r=False)
    tau_str = (f"{fit_fixed['tau_decay_hr']:.1f} hrs"
               if np.isfinite(fit_fixed['tau_decay_hr']) else "inf")
    print(f"    A0 = {fit_fixed['A0']:.4f}, t0_eff = {fit_fixed['t0_eff_hr']:.2f} hrs, "
          f"phi_rt = {fit_fixed['phi_rt']:.3f} rad, "
          f"tau = {tau_str}, residual = {fit_fixed['residual_norm']:.1%}")

    # Try with r free
    print("  Fitting cavity model (r free)...")
    fit_free = fit_cavity_model(envelope, times, params, free_r=True)
    tau_str2 = (f"{fit_free['tau_decay_hr']:.1f} hrs"
                if np.isfinite(fit_free['tau_decay_hr']) else "inf")
    print(f"    A0 = {fit_free['A0']:.4f}, t0_eff = {fit_free['t0_eff_hr']:.2f} hrs, "
          f"phi_rt = {fit_free['phi_rt']:.3f} rad, "
          f"r = {fit_free['r']:.4f}, tau = {tau_str2}, "
          f"residual = {fit_free['residual_norm']:.1%}")

    # Choose best fit
    if fit_free["residual_norm"] < fit_fixed["residual_norm"] * 0.9:
        print("  -> Using free-r fit (substantially better)")
        fit_result = fit_free
    else:
        print("  -> Using fixed-r fit (free-r not substantially better)")
        fit_result = fit_fixed

    # Clean window
    cw = compute_clean_window(params)
    if cw["exists"]:
        print(f"  Clean window: {cw['start_s']/3600:.2f} - {cw['end_s']/3600:.2f} hrs "
              f"({cw['duration_s']/60:.0f} min)")
    else:
        print(f"  Clean window: NONE ({cw['note']})")

    return {
        "times": times,
        "w_mid": w_mid,
        "envelope": envelope,
        "fit_result": fit_result,
        "fit_fixed": fit_fixed,
        "fit_free": fit_free,
        "clean_window": cw,
        "params": params,
    }


def main():
    wall_start = timer.time()

    # Analyze both cases
    r3 = analyze_case(R3_PARAMS)
    r2 = analyze_case(R2_PARAMS)

    # Generate figures
    print(f"\n{'='*60}")
    print("Generating figures...")
    print(f"{'='*60}")

    print("\nFigure 1: R3 Envelope Analysis")
    plot_fig1_R3_envelope(
        r3["times"], r3["w_mid"], r3["envelope"],
        r3["fit_result"], r3["params"],
    )

    print("\nFigure 2: R2 Envelope Analysis")
    plot_fig2_R2_envelope(
        r2["times"], r2["w_mid"], r2["envelope"],
        r2["fit_result"], r2["params"],
    )

    print("\nFigure 3: R2 vs R3 Comparison")
    plot_fig3_comparison(
        r3["times"], r3["envelope"], r3["fit_result"], r3["params"],
        r2["times"], r2["envelope"], r2["fit_result"], r2["params"],
    )

    print("\nFigure 4: Handoff Viability")
    plot_fig4_handoff_viability(
        r3["times"], r3["w_mid"], r3["envelope"], r3["params"],
        r2["times"], r2["w_mid"], r2["envelope"], r2["params"],
    )

    wall_end = timer.time()
    wall_s = wall_end - wall_start
    print(f"\nTotal wall-clock time: {wall_s:.1f} s ({wall_s/60:.1f} min)")

    # Print summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    for case in [r3, r2]:
        p = case["params"]
        f = case["fit_result"]
        cw = case["clean_window"]
        env = case["envelope"]
        print(f"\n{p.name}:")
        print(f"  c_gz empirical: {p.c_gz_empirical:.1f} m/s")
        print(f"  T_roundtrip: {p.T_roundtrip:.0f} s ({p.T_roundtrip/60:.1f} min)")
        print(f"  phi_rt (fit/theory): {f['phi_rt']:.2f} / {p.phi_rt_theory:.2f} rad")
        print(f"  r: {f['r']:.4f}")
        print(f"  t0_eff: {f['t0_eff_hr']:.2f} hrs")
        tau_s = (f"{f['tau_decay_hr']:.1f} hrs"
                 if np.isfinite(f['tau_decay_hr']) else "inf")
        print(f"  tau_decay: {tau_s}")
        print(f"  Peak envelope: {np.max(env):.3f} m/s")
        print(f"  End envelope: {env[-1]:.3f} m/s")
        print(f"  Peak/End ratio: {np.max(env)/env[-1]:.1f}x")
        if cw["exists"]:
            print(f"  Clean window: {cw['duration_s']/60:.0f} min")
        else:
            print(f"  Clean window: NONE")
        print(f"  Fit residual (normalized): {f['residual_norm']:.1%}")


if __name__ == "__main__":
    main()
