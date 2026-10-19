"""
Comparison of Klemp-Durran radiation BC vs Eq.19 extrapolation BC
for R2 (gravity wave) and R3 (acoustic wave) simulations.

Loads handoff data from all 4 runs and produces overlay envelope
comparisons to quantify the improvement (or lack thereof) from the
KD radiation condition.

Produces 3 figures saved to output/kd_comparison/.
"""

import time as timer
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import hilbert

import netCDF4 as nc


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

@dataclass
class CaseConfig:
    """Configuration for a single simulation case."""
    name: str
    label: str
    omega: float
    data_path: str
    color: str
    linestyle: str


CASES = {
    "R2_eq19": CaseConfig(
        name="R2 (Gravity)",
        label="R2 Eq.19",
        omega=0.005,
        data_path="output/tier2_R2/dsmc_handoff_100km.nc",
        color="tab:red",
        linestyle="--",
    ),
    "R2_kd": CaseConfig(
        name="R2 (Gravity)",
        label="R2 KD-200km",
        omega=0.005,
        data_path="output/tier2_R2_kd/dsmc_handoff_100km.nc",
        color="tab:blue",
        linestyle="-.",
    ),
    "R2_kd_120": CaseConfig(
        name="R2 (Gravity)",
        label="R2 KD-120km",
        omega=0.005,
        data_path="output/tier2_R2_kd_120km/dsmc_handoff_100km.nc",
        color="tab:green",
        linestyle="-",
    ),
    "R3_eq19": CaseConfig(
        name="R3 (Acoustic)",
        label="R3 Eq.19",
        omega=0.02,
        data_path="output/tier2_R3/dsmc_handoff_100km.nc",
        color="tab:red",
        linestyle="--",
    ),
    "R3_kd": CaseConfig(
        name="R3 (Acoustic)",
        label="R3 Klemp-Durran",
        omega=0.02,
        data_path="output/tier2_R3_kd/dsmc_handoff_100km.nc",
        color="tab:blue",
        linestyle="-",
    ),
}

OUTPUT_DIR = Path("output/kd_comparison")


# ---------------------------------------------------------------------------
# Data loading (reused from analyze_cavity_interference.py)
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
# Hilbert envelope (reused from analyze_cavity_interference.py)
# ---------------------------------------------------------------------------

def compute_hilbert_envelope(
    w: np.ndarray,
    times: np.ndarray,
    omega: float,
    n_smooth_periods: float = 2.0,
) -> np.ndarray:
    """
    Compute the smoothed amplitude envelope via the Hilbert transform.
    Gaussian smoothing over n_smooth_periods wave periods.
    """
    analytic = hilbert(w)
    envelope = np.abs(analytic)

    period_s = 2 * np.pi / omega
    dt = np.median(np.diff(times))
    sigma_samples = (n_smooth_periods * period_s) / dt
    envelope_smooth = gaussian_filter1d(envelope, sigma=sigma_samples)

    return envelope_smooth


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

def compute_metrics(
    times: np.ndarray,
    envelope: np.ndarray,
    omega: float,
) -> dict:
    """Compute comparison metrics for an envelope."""
    peak = np.max(envelope)
    end_val = envelope[-1]
    peak_end_ratio = peak / end_val if end_val > 0 else np.inf

    # Time to reach 90% of peak
    threshold_90 = 0.9 * peak
    i_90 = np.argmax(envelope >= threshold_90)
    t_90 = times[i_90] / 3600.0 if i_90 > 0 else np.nan

    # Envelope flatness: std/mean in the second half of the simulation
    # (after ramp-up transient has settled)
    n_half = len(envelope) // 2
    env_second_half = envelope[n_half:]
    flatness = np.std(env_second_half) / np.mean(env_second_half) if np.mean(env_second_half) > 0 else np.inf

    # Steady-state amplitude: mean of last 20% of envelope
    n_tail = max(1, len(envelope) // 5)
    steady_state = np.mean(envelope[-n_tail:])

    return {
        "peak": peak,
        "end_val": end_val,
        "peak_end_ratio": peak_end_ratio,
        "t_90_hr": t_90,
        "flatness_cv": flatness,
        "steady_state": steady_state,
    }


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def plot_case_comparison(
    key_eq19: str,
    key_kd: str,
    case_label: str,
    data: dict,
    fig_name: str,
):
    """
    Two-panel comparison plot for one wave type (R2 or R3).

    Top panel: raw w(t) for both BCs overlaid with Hilbert envelopes.
    Bottom panel: normalized envelopes with metrics annotated.
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 8), sharex=True)

    d_old = data[key_eq19]
    d_new = data[key_kd]

    t_hr_old = d_old["times"] / 3600.0
    t_hr_new = d_new["times"] / 3600.0

    # --- Top panel: raw signals and envelopes ---
    ax1.plot(t_hr_old, d_old["w_mid"], color="tab:red", alpha=0.15,
             linewidth=0.3, rasterized=True)
    ax1.plot(t_hr_new, d_new["w_mid"], color="tab:blue", alpha=0.15,
             linewidth=0.3, rasterized=True)

    ax1.plot(t_hr_old, d_old["envelope"], color="tab:red", linewidth=1.5,
             linestyle="--", label="Eq.19 envelope")
    ax1.plot(t_hr_old, -d_old["envelope"], color="tab:red", linewidth=1.0,
             linestyle="--", alpha=0.5)
    ax1.plot(t_hr_new, d_new["envelope"], color="tab:blue", linewidth=1.5,
             linestyle="-", label="KD envelope")
    ax1.plot(t_hr_new, -d_new["envelope"], color="tab:blue", linewidth=1.0,
             linestyle="-", alpha=0.5)

    ax1.set_ylabel("w (m/s)")
    ax1.set_title(f"{case_label}: Eq.19 vs Klemp-Durran at z = 100 km")
    ax1.legend(fontsize=9, loc="upper right")

    # --- Bottom panel: normalized envelopes ---
    env_peak_old = np.max(d_old["envelope"])
    env_peak_new = np.max(d_new["envelope"])
    env_peak_max = max(env_peak_old, env_peak_new)

    ax2.plot(t_hr_old, d_old["envelope"] / env_peak_max, color="tab:red",
             linewidth=1.5, linestyle="--", label="Eq.19 (normalized)")
    ax2.plot(t_hr_new, d_new["envelope"] / env_peak_max, color="tab:blue",
             linewidth=1.5, linestyle="-", label="KD (normalized)")

    # Annotate metrics
    m_old = d_old["metrics"]
    m_new = d_new["metrics"]

    text_old = (
        f"Eq.19:\n"
        f"  Peak/End = {m_old['peak_end_ratio']:.1f}x\n"
        f"  Peak = {m_old['peak']:.3f} m/s\n"
        f"  CV(2nd half) = {m_old['flatness_cv']:.2f}\n"
        f"  t(90%) = {m_old['t_90_hr']:.2f} hrs"
    )
    text_new = (
        f"Klemp-Durran:\n"
        f"  Peak/End = {m_new['peak_end_ratio']:.1f}x\n"
        f"  Peak = {m_new['peak']:.3f} m/s\n"
        f"  CV(2nd half) = {m_new['flatness_cv']:.2f}\n"
        f"  t(90%) = {m_new['t_90_hr']:.2f} hrs"
    )

    ax2.text(0.02, 0.95, text_old, transform=ax2.transAxes, fontsize=8,
             verticalalignment="top", fontfamily="monospace",
             color="tab:red",
             bbox=dict(boxstyle="round", facecolor="mistyrose", alpha=0.7))
    ax2.text(0.02, 0.48, text_new, transform=ax2.transAxes, fontsize=8,
             verticalalignment="top", fontfamily="monospace",
             color="tab:blue",
             bbox=dict(boxstyle="round", facecolor="lavender", alpha=0.7))

    ax2.set_xlabel("Time (hrs)")
    ax2.set_ylabel("Normalized amplitude")
    ax2.legend(fontsize=9, loc="upper right")

    plt.tight_layout()
    outpath = OUTPUT_DIR / fig_name
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {outpath}")


def plot_summary(data: dict):
    """
    Summary figure: all 4 normalized envelopes + comparison table.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5),
                                    gridspec_kw={"width_ratios": [1.2, 1]})

    # --- Left: all normalized envelopes ---
    plot_configs = [
        ("R2_eq19", "R2 Eq.19", "tab:red", "--"),
        ("R2_kd", "R2 KD-200km", "tab:blue", "-."),
        ("R2_kd_120", "R2 KD-120km", "tab:green", "-"),
        ("R3_eq19", "R3 Eq.19", "tab:orange", "--"),
        ("R3_kd", "R3 KD", "tab:purple", "-"),
    ]

    for key, label, color, ls in plot_configs:
        if key not in data:
            continue
        d = data[key]
        env_norm = d["envelope"] / np.max(d["envelope"])
        t_hr = d["times"] / 3600.0
        ax1.plot(t_hr, env_norm, color=color, linestyle=ls,
                 linewidth=1.3, label=label)

    ax1.set_xlabel("Time (hrs)")
    ax1.set_ylabel("Normalized envelope")
    ax1.set_title("Amplitude Envelopes (all cases)")
    ax1.legend(fontsize=9)
    ax1.set_ylim(0, 1.15)

    # --- Right: summary table ---
    ax2.axis("off")

    headers = ["Metric"]
    rows_data = [
        ["Peak (m/s)"],
        ["Peak/End"],
        ["CV (2nd half)"],
        ["t(90%) (hrs)"],
        ["Steady-state (m/s)"],
    ]

    for key, label, _, _ in plot_configs:
        if key not in data:
            continue
        m = data[key]["metrics"]
        headers.append(label)
        rows_data[0].append(f"{m['peak']:.4f}")
        rows_data[1].append(f"{m['peak_end_ratio']:.2f}x")
        rows_data[2].append(f"{m['flatness_cv']:.3f}")
        rows_data[3].append(f"{m['t_90_hr']:.2f}")
        rows_data[4].append(f"{m['steady_state']:.4f}")

    table_data = [headers] + rows_data
    table = ax2.table(cellText=table_data, loc="center", cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.0, 1.6)

    # Style header row
    for j in range(len(headers)):
        table[0, j].set_facecolor("#4472C4")
        table[0, j].set_text_props(color="white", fontweight="bold")

    ax2.set_title("Comparison Metrics", fontsize=11, pad=15)

    plt.tight_layout()
    outpath = OUTPUT_DIR / "fig3_summary.png"
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {outpath}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def load_case(key: str, cfg: CaseConfig) -> dict:
    """Load and process one case."""
    path = Path(cfg.data_path)
    if not path.exists():
        print(f"  SKIP: {cfg.data_path} not found")
        return None

    print(f"  Loading {cfg.label} from {cfg.data_path}...")
    times, w_mid = load_handoff_data(cfg.data_path)
    envelope = compute_hilbert_envelope(w_mid, times, cfg.omega)
    metrics = compute_metrics(times, envelope, cfg.omega)

    print(f"    N_t={len(times)}, duration={times[-1]/3600:.2f} hrs, "
          f"peak={metrics['peak']:.4f} m/s, "
          f"peak/end={metrics['peak_end_ratio']:.2f}x, "
          f"CV={metrics['flatness_cv']:.3f}")

    return {
        "times": times,
        "w_mid": w_mid,
        "envelope": envelope,
        "metrics": metrics,
        "config": cfg,
    }


def main():
    wall_start = timer.time()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load all available cases
    print("Loading handoff data...")
    data = {}
    for key, cfg in CASES.items():
        result = load_case(key, cfg)
        if result is not None:
            data[key] = result

    if len(data) == 0:
        print("ERROR: No data files found. Run the simulations first.")
        return

    # Generate figures
    print(f"\n{'='*60}")
    print("Generating comparison figures...")
    print(f"{'='*60}")

    # Fig 1: R2 comparison (if both available)
    if "R2_eq19" in data and "R2_kd" in data:
        print("\nFigure 1: R2 Envelope Comparison (Eq.19 vs KD-200km)")
        plot_case_comparison(
            "R2_eq19", "R2_kd",
            "R2 (Gravity Wave, w=0.05 m/s) Eq.19 vs KD-200km",
            data, "fig1_R2_comparison.png",
        )

    # Fig 1b: R2 Eq.19 vs KD-120km
    if "R2_eq19" in data and "R2_kd_120" in data:
        print("\nFigure 1b: R2 Envelope Comparison (Eq.19 vs KD-120km)")
        plot_case_comparison(
            "R2_eq19", "R2_kd_120",
            "R2 (Gravity Wave, w=0.05 m/s) Eq.19 vs KD-120km",
            data, "fig1b_R2_kd120_comparison.png",
        )
    else:
        print("\nSkipping Fig 1b: R2 KD-120km data not found")

    # Fig 2: R3 comparison (if both available)
    if "R3_eq19" in data and "R3_kd" in data:
        print("\nFigure 2: R3 Envelope Comparison")
        plot_case_comparison(
            "R3_eq19", "R3_kd",
            "R3 (Acoustic Wave, w=0.01 m/s, omega=0.02 rad/s)",
            data, "fig2_R3_comparison.png",
        )
    elif "R3_kd" in data:
        print("\nSkipping Fig 2: R3 Eq.19 data not found (only KD available)")
    else:
        print("\nSkipping Fig 2: R3 KD data not found (run simulation first)")

    # Fig 3: Summary (with whatever data is available)
    if len(data) >= 2:
        print("\nFigure 3: Summary Comparison")
        plot_summary(data)

    # Print summary
    wall_end = timer.time()
    wall_s = wall_end - wall_start
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")

    for key in ["R2_eq19", "R2_kd", "R2_kd_120", "R3_eq19", "R3_kd"]:
        if key not in data:
            continue
        d = data[key]
        m = d["metrics"]
        print(f"\n{d['config'].label}:")
        print(f"  Peak amplitude: {m['peak']:.4f} m/s")
        print(f"  Peak/End ratio: {m['peak_end_ratio']:.2f}x")
        print(f"  Envelope CV (2nd half): {m['flatness_cv']:.3f}")
        print(f"  t(90% peak): {m['t_90_hr']:.2f} hrs")
        print(f"  Steady-state amplitude: {m['steady_state']:.4f} m/s")

    # Improvement summary
    comparisons = [
        ("R2", "R2_eq19", "R2_kd", "KD-200km"),
        ("R2", "R2_eq19", "R2_kd_120", "KD-120km"),
        ("R3", "R3_eq19", "R3_kd", "KD"),
    ]
    for wave, k_old, k_new, label in comparisons:
        if k_old in data and k_new in data:
            r_old = data[k_old]["metrics"]["peak_end_ratio"]
            r_new = data[k_new]["metrics"]["peak_end_ratio"]
            cv_old = data[k_old]["metrics"]["flatness_cv"]
            cv_new = data[k_new]["metrics"]["flatness_cv"]
            print(f"\n{wave} vs {label}:")
            print(f"  Peak/End ratio: {r_old:.2f}x -> {r_new:.2f}x "
                  f"({'improved' if r_new < r_old else 'no improvement'})")
            print(f"  Envelope CV: {cv_old:.3f} -> {cv_new:.3f} "
                  f"({'improved' if cv_new < cv_old else 'no improvement'})")

    print(f"\nTotal wall-clock time: {wall_s:.1f} s")


if __name__ == "__main__":
    main()
