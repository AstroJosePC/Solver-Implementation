#!/usr/bin/env python3
"""
analyze_arsia_stage1_deepdive.py

Deep-dive diagnostics for a single Arsia Mons Stage 1 mountain-wave run
(uniform eastward steering wind, 90 km domain, sponge 60-90 km). Complements
the sweep-level figures from analyze_mars_mountain_sweep.py with:

1) Data inventory (snapshot count, cadence, variables, on-disk size) -> JSON
2) Full-record wavelet spectrogram at station S3 (no trim), so the startup
   transient and any late-time high-frequency burst appear together
3) High-altitude station time series at z ~ 57.5 km (highest sponge-free
   level band; sponge base is 60 km)
4) Vertical profiles below the sponge base: hill-window max/RMS |w|(z) and
   the x-mean momentum flux MF(z) = <rho u' w'>_x at selected times
5) Mean-flow deceleration profile delta_u(z) = <u>_x(z, t) - U_bar
6) Transient-timing metrics: dominant ring-down period (FFT of the early
   window) and the time at which the station envelope falls below a
   threshold, quantifying when the startup transient has passed
7) Convergence time series at fixed heights (AMEC-band hill amplitude, MF at
   selected levels, peak mean-flow deviation below the sponge base) with a
   per-window linear trend, which is what distinguishes a run that has
   settled from one still evolving at t_end. An optional shorter baseline
   run can be overlaid for a like-for-like comparison.

All figures/JSON go to an analysis directory with a "deepdive_" prefix.
Analysis is restricted to z < 60 km (sponge-free region) per the scoping
decision; the sponge region is excluded from metrics.

Usage:
  python3 scripts/analyze_arsia_stage1_deepdive.py \
      --run-dir output/arsia_mons_stage1/run_U57.0_h3.0
  python3 scripts/analyze_arsia_stage1_deepdive.py \
      --run-dir output/arsia_mons_stage1_ext6h/run_U40.0_h3.0 \
      --compare-run output/arsia_mons_stage1/run_U40.0_h3.0

U_bar, the run label and the outdir are inferred from the run directory name
(run_U<U>_h<h>) unless given explicitly.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import re
import time

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset

from wavelet_tc98 import wavelet

SPONGE_BASE_KM = 60.0
HILL_WINDOW_KM = 500.0        # |x_rel| <= 500 km, same as sweep analysis
STATION_X_REL_KM = [50.0, 150.0, 300.0, 500.0]
STATION_Z_LOW_KM = 50.0       # AMEC-band station used by the sweep analysis
STATION_Z_HIGH_KM = 57.5      # highest sponge-free band (below 60 km)
WAVELET_STATION_IDX = 2       # S3, matching the sweep analysis
BUOYANCY_N_REF = 0.018        # s^-1, the near-surface inversion value used for the
                              # Nh/U label in the sweep scripts. It is NOT a
                              # lower-atmosphere mean: N over 10-55 km in this
                              # column is 0.0078-0.0097 s^-1. Kept only as the
                              # launch-level reference line on the wavelets, since
                              # a surface-launched transient does carry it. Every
                              # frequency compared against a measurement made at
                              # altitude uses background_frequencies() instead.
EARLY_WINDOW_S = 4000.0       # window for the transient period estimate
WAVELET_BG_FRAC = 0.5         # the wavelet background is taken over the last
                              # this-fraction of the record. The transient bound
                              # EARLY_WINDOW_S is too early to serve: in the 12 h
                              # runs the transient tail is still visible at 200 min.
WAVELET_REL_FLOOR_DB = -40.0  # cells this far below the record peak are blanked on
                              # the background-normalised panel. Dividing by a
                              # background that is itself the noise floor produces
                              # tens of dB of contrast out of nothing, so the panel
                              # is only drawn where there is power to normalise.
AMEC_BAND_KM = (45.0, 53.0)   # cloud-formation band
MF_LEVELS_KM = [30.0, 45.0, 52.0]   # levels tracked in the convergence series

# Mars thermodynamics, matching mars_config in src/parameters.jl (gamma = 1.31,
# g = 3.72, R_gas = R_universal / M_CO2). kappa = R/cp = (gamma-1)/gamma.
G_MARS = 3.72
GAMMA = 1.31
R_GAS = 188.9
KAPPA = (GAMMA - 1.0) / GAMMA
CP_MARS = GAMMA * R_GAS / (GAMMA - 1.0)   # ~798 J/(kg K)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Deep dive on one Arsia Stage 1 run")
    p.add_argument("--run-dir", default="output/arsia_mons_stage1/run_U57.0_h3.0")
    p.add_argument("--outdir", default=None,
                   help="default: <run parent>_analysis/<run name>")
    p.add_argument("--u-bar", type=float, default=None,
                   help="steering wind [m/s]; default parsed from the run name")
    p.add_argument("--label", default=None,
                   help="figure label; default is the run directory name")
    p.add_argument("--profile-hours", default=None,
                   help="comma-separated times [h] for the vertical profiles")
    p.add_argument("--compare-run", default=None,
                   help="shorter baseline run overlaid on the convergence series")
    return p.parse_args()


def snapshot_files(run_dir: str) -> list[str]:
    files = sorted(glob.glob(os.path.join(run_dir, "snapshot_*.nc")))
    if not files:
        raise FileNotFoundError(f"No snapshots found in {run_dir}")
    return files


def u_bar_from_name(run_dir: str) -> float:
    m = re.search(r"run_U([0-9.]+)_h", os.path.basename(os.path.normpath(run_dir)))
    if not m:
        raise ValueError(f"Cannot parse U_bar from {run_dir}; pass --u-bar")
    return float(m.group(1))


def load_grid(path: str):
    with Dataset(path) as ds:
        x = np.asarray(ds["x"][:], dtype=np.float64)
        z = np.asarray(ds["z"][:], dtype=np.float64)
        varnames = [v for v in ds.variables if v not in ("x", "z")]
    return x, z, varnames


def load_time(path: str) -> float:
    with Dataset(path) as ds:
        return float(ds.getncattr("time_s"))


def background_frequencies(path):
    """Background N(z) and acoustic cutoff omega_a(z) from the first snapshot.

    Taken from the x-mean temperature of the initial state, so the reference
    frequencies drawn on the wavelets and quoted against the station spectra are
    the ones this run actually carries at the station height rather than a single
    near-surface value. omega_a = gamma g / (2 c), the isothermal cutoff, equal to
    c / (2H) with H = R T / g.

    The two are not interchangeable references. A disturbance at omega < N is an
    internal gravity wave, one at omega > omega_a is a freely propagating acoustic
    wave, and the band between them is evanescent in a stratified column.
    """
    with Dataset(path) as ds:
        z = np.asarray(ds["z"][:], dtype=np.float64)
        T = np.asarray(ds["T"][:], dtype=np.float64).mean(axis=1)   # (z, x) -> (z,)
    n = np.sqrt(np.clip((G_MARS / T) * (np.gradient(T, z) + G_MARS / CP_MARS),
                        1e-12, None))
    c = np.sqrt(GAMMA * R_GAS * T)
    return n, GAMMA * G_MARS / (2.0 * c)


def default_profile_hours(t_end_h: float) -> list[float]:
    """Fixed early times plus the record end, so figures stay comparable.

    The 2% tolerance keeps a round-number hour that the run ends just short of
    (a 3 h run stops near 2.99 h, since t_end lands between snapshots).
    """
    hours = [0.5, 1.0, 2.0, 3.0]
    if t_end_h > 3.5:
        hours += [4.5, round(t_end_h, 2)]
    return [h for h in hours if h <= t_end_h * 1.02]


def trend_per_hour(times: np.ndarray, series: np.ndarray,
                   t0: float, t1: float) -> dict:
    """Least-squares slope of `series` over [t0, t1], absolute and normalized.

    The normalization is by the window mean magnitude, so `frac_per_hour` is a
    scale-free measure: values near zero indicate a settled quantity, whereas
    order-unity values indicate one still evolving on the run's own timescale.
    """
    sel = (times >= t0) & (times <= t1)
    if sel.sum() < 3:
        return {"window_s": [t0, t1], "n": int(sel.sum()), "slope_per_hour": None,
                "mean": None, "frac_per_hour": None, "cov": None}
    tt, ss = times[sel], series[sel]
    slope = float(np.polyfit(tt, ss, 1)[0]) * 3600.0
    mean = float(ss.mean())
    scale = float(np.abs(ss).mean())
    return {
        "window_s": [float(t0), float(t1)],
        "n": int(sel.sum()),
        "slope_per_hour": slope,
        "mean": mean,
        "frac_per_hour": (slope / scale) if scale > 0 else None,
        "cov": float(ss.std() / scale) if scale > 0 else None,
    }


def run_series(run_dir: str, u_bar: float) -> dict:
    """Fixed-height time series used for the convergence assessment.

    Returns AMEC-band hill-window mean |w|, MF at MF_LEVELS_KM, and the peak
    mean-flow deviation below the sponge base, all as functions of time.
    """
    files = snapshot_files(run_dir)
    x, z, _ = load_grid(files[0])
    z_km = z / 1e3
    x_rel_km = (x - 0.5 * (x[0] + x[-1])) / 1e3
    hill = np.abs(x_rel_km) <= HILL_WINDOW_KM
    amec = (z_km >= AMEC_BAND_KM[0]) & (z_km <= AMEC_BAND_KM[1])
    jz_clean = z_km < SPONGE_BASE_KM
    iz_mf = [int(np.argmin(np.abs(z_km - lv))) for lv in MF_LEVELS_KM]

    n_t = len(files)
    out = {
        "times": np.empty(n_t),
        "amec_w": np.empty(n_t),
        "mf": np.empty((len(iz_mf), n_t)),
        "du_max": np.empty(n_t),
        "z_du_max": np.empty(n_t),
    }
    for it, f in enumerate(files):
        with Dataset(f) as ds:
            u = np.asarray(ds["u"][:], dtype=np.float64)
            w = np.asarray(ds["w"][:], dtype=np.float64)
            rho = np.asarray(ds["rho"][:], dtype=np.float64)
            out["times"][it] = float(ds.getncattr("time_s"))
        ub, wb = u.mean(axis=1), w.mean(axis=1)
        mf = (rho * (u - ub[:, None]) * (w - wb[:, None])).mean(axis=1)
        out["amec_w"][it] = np.abs(w[np.ix_(amec, hill)]).mean()
        out["mf"][:, it] = mf[iz_mf]
        du = ub[jz_clean] - u_bar
        i_du = int(np.argmax(np.abs(du)))
        out["du_max"][it] = du[i_du]
        out["z_du_max"][it] = z_km[jz_clean][i_du]
    return out


def main() -> None:
    t_wall0 = time.perf_counter()
    args = parse_args()

    run_dir = os.path.normpath(args.run_dir)
    run_name = os.path.basename(run_dir)
    u_bar = args.u_bar if args.u_bar is not None else u_bar_from_name(run_dir)
    label = args.label or f"{os.path.basename(os.path.dirname(run_dir))}/{run_name}"
    outdir = args.outdir or os.path.join(
        os.path.dirname(run_dir) + "_analysis", run_name)
    os.makedirs(outdir, exist_ok=True)

    files = snapshot_files(run_dir)
    x, z, varnames = load_grid(files[0])
    n_bg, wa_bg = background_frequencies(files[0])
    x_rel_km = (x - 0.5 * (x[0] + x[-1])) / 1e3
    z_km = z / 1e3

    times = np.array([load_time(f) for f in files])
    dt_snap = float(np.median(np.diff(times)))
    size_bytes = sum(os.path.getsize(f) for f in files)
    t_end_h = times[-1] / 3600.0

    jz_clean = z_km < SPONGE_BASE_KM              # sponge-free levels
    hill = np.abs(x_rel_km) <= HILL_WINDOW_KM     # hill window columns
    amec = (z_km >= AMEC_BAND_KM[0]) & (z_km <= AMEC_BAND_KM[1])

    iz_low = int(np.argmin(np.abs(z_km - STATION_Z_LOW_KM)))
    iz_high = int(np.argmin(np.abs(z_km - STATION_Z_HIGH_KM)))
    ix_st = [int(np.argmin(np.abs(x_rel_km - sx))) for sx in STATION_X_REL_KM]
    iz_mf = [int(np.argmin(np.abs(z_km - lv))) for lv in MF_LEVELS_KM]

    n_t, n_st = len(files), len(ix_st)
    w_low = np.empty((n_st, n_t))
    w_high = np.empty((n_st, n_t))
    ubar_zt = np.empty((len(z), n_t))             # <u>_x(z, t)
    mf_zt = np.empty((len(z), n_t))               # <rho u' w'>_x(z, t)
    amec_w_t = np.empty(n_t)                      # hill-window mean |w| in band

    if args.profile_hours:
        profile_hours = [float(v) for v in args.profile_hours.split(",")]
    else:
        profile_hours = default_profile_hours(t_end_h)
    # Requested hours can collide on one snapshot (short record, coarse cadence);
    # keep the first of each so the profile arrays stay aligned with the labels.
    it_profiles, keep = [], []
    for h in profile_hours:
        it = int(np.argmin(np.abs(times - h * 3600.0)))
        if it not in it_profiles:
            it_profiles.append(it)
            keep.append(h)
    profile_hours = keep
    wmax_profiles = np.empty((len(it_profiles), len(z)))
    wrms_profiles = np.empty((len(it_profiles), len(z)))

    # Time-height arrays kept for the convergence figures and the cache. The
    # static-stability terms are a wave-breaking indicator: with
    # theta = T (P_surface/P)^kappa and N^2 = (g/theta) dtheta/dz, a column with
    # N^2 < 0 has overturned, so the fraction of unstable hill columns at each
    # level tracks the onset and spread of breaking.
    unstable_frac_zt = np.empty((len(z), n_t))
    umin_hill_zt = np.empty((len(z), n_t))
    wmax_hill_zt = np.empty((len(z), n_t))
    wrms_hill_zt = np.empty((len(z), n_t))

    for it, f in enumerate(files):
        with Dataset(f) as ds:
            w = np.asarray(ds["w"][:], dtype=np.float64)     # (z, x)
            u = np.asarray(ds["u"][:], dtype=np.float64)
            rho = np.asarray(ds["rho"][:], dtype=np.float64)
            T = np.asarray(ds["T"][:], dtype=np.float64)
            P = np.asarray(ds["P"][:], dtype=np.float64)
        for i, ix in enumerate(ix_st):
            w_low[i, it] = w[iz_low, ix]
            w_high[i, it] = w[iz_high, ix]
        ub = u.mean(axis=1)
        wb = w.mean(axis=1)
        ubar_zt[:, it] = ub
        mf_zt[:, it] = (rho * (u - ub[:, None]) * (w - wb[:, None])).mean(axis=1)
        amec_w_t[it] = np.abs(w[np.ix_(amec, hill)]).mean()

        theta = T * (P[0:1, :] / P) ** KAPPA
        n2 = np.gradient(theta, z, axis=0) * G_MARS / theta
        unstable_frac_zt[:, it] = (n2[:, hill] < 0.0).mean(axis=1)
        umin_hill_zt[:, it] = u[:, hill].min(axis=1)
        wmax_hill_zt[:, it] = np.abs(w[:, hill]).max(axis=1)
        wrms_hill_zt[:, it] = np.sqrt((w[:, hill] ** 2).mean(axis=1))

        if it in it_profiles:
            k = it_profiles.index(it)
            wmax_profiles[k] = wmax_hill_zt[:, it]
            wrms_profiles[k] = wrms_hill_zt[:, it]

    du_zt = ubar_zt[jz_clean, :] - u_bar
    i_du = np.argmax(np.abs(du_zt), axis=0)
    du_max_t = du_zt[i_du, np.arange(n_t)]
    z_du_max_t = z_km[jz_clean][i_du]

    # Column-integrated breaking indicators over the sponge-free region.
    unstable_frac_t = unstable_frac_zt[jz_clean].mean(axis=0)
    umin_hill_t = umin_hill_zt[jz_clean].min(axis=0)
    wmax_hill_t = wmax_hill_zt[jz_clean].max(axis=0)
    umean_min_t = ubar_zt[jz_clean].min(axis=0)

    # Shared cache so downstream figures reuse one verified extraction rather
    # than re-reading the snapshots.
    np.savez_compressed(
        os.path.join(outdir, "deepdive_series.npz"),
        times=times, z_km=z_km, u_bar=np.array(u_bar), label=np.array(label),
        sponge_base_km=np.array(SPONGE_BASE_KM),
        amec_band_km=np.array(AMEC_BAND_KM),
        mf_levels_km=np.array(MF_LEVELS_KM),
        ubar_zt=ubar_zt, mf_zt=mf_zt, unstable_frac_zt=unstable_frac_zt,
        umin_hill_zt=umin_hill_zt, wmax_hill_zt=wmax_hill_zt,
        wrms_hill_zt=wrms_hill_zt,
        amec_w_t=amec_w_t, du_max_t=du_max_t, z_du_max_t=z_du_max_t,
        unstable_frac_t=unstable_frac_t, umin_hill_t=umin_hill_t,
        wmax_hill_t=wmax_hill_t, umean_min_t=umean_min_t,
    )

    # --- 1) Station time series at z ~ 57.5 km (high) vs 50 km (low)
    fig, axes = plt.subplots(n_st, 1, figsize=(10, 2.2 * n_st), sharex=True)
    for i, ax in enumerate(np.atleast_1d(axes)):
        ax.plot(times, w_high[i], lw=1.0, label=f"z={z_km[iz_high]:.2f} km")
        ax.plot(times, w_low[i], lw=0.8, alpha=0.6, label=f"z={z_km[iz_low]:.2f} km")
        ax.text(0.985, 0.9, f"S{i+1}: x_rel={x_rel_km[ix_st[i]]:.0f} km",
                transform=ax.transAxes, ha="right", va="top")
        ax.grid(alpha=0.3)
        ax.set_ylabel("w (m/s)")
        if i == 0:
            ax.legend(loc="upper left", fontsize=8, ncol=2)
    np.atleast_1d(axes)[-1].set_xlabel("time (s)")
    fig.suptitle(f"{label} station w near sponge base vs AMEC level")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "deepdive_station_timeseries_w_z57.png"), dpi=150)
    plt.close(fig)

    # --- 2) Full-record wavelet at S3, both altitudes
    for tag, series, izw in (("z50", w_low[WAVELET_STATION_IDX], iz_low),
                             ("z57", w_high[WAVELET_STATION_IDX], iz_high)):
        wave, period, _, coi = wavelet(series, dt_snap, pad=1, dj=0.05,
                                       mother="MORLET", param=6.0)
        power = np.abs(wave) ** 2
        freq_mhz = 1e3 / period
        tmin = times / 60.0
        # The startup transient is two orders of magnitude above everything that
        # follows, so a linear scale normalised by the global peak renders the
        # settled record at a few per cent of full scale and it disappears. Both
        # panels are therefore logarithmic, and they answer different questions.
        # The top one keeps a single absolute scale, so amplitudes stay comparable
        # across time and between the two altitudes. The bottom one divides each
        # frequency by its own median over the settled part of the record, which
        # shows a weak late feature at its own level but discards amplitude and
        # will manufacture contrast wherever the background is near zero. That is
        # why it is blanked below WAVELET_REL_FLOOR_DB and why the two panels are
        # kept together: the top one is the check on the bottom one.
        floor = 1e-300
        p_abs = 10 * np.log10(np.maximum(power, floor) / power.max())
        bg_start = (1.0 - WAVELET_BG_FRAC) * times[-1]
        bg_mask = times > bg_start
        if bg_mask.sum() < 10:                      # short record, e.g. a smoke run
            bg_mask = times > times[len(times) // 2]
        bg = np.median(power[:, bg_mask], axis=1, keepdims=True)
        p_rel = 10 * np.log10(np.maximum(power, floor) / np.maximum(bg, floor))
        # Blank the relative panel wherever there is no absolute power to speak of,
        # so amplified floor cannot be read as a late feature.
        p_rel = np.where(p_abs > WAVELET_REL_FLOOR_DB, p_rel, np.nan)

        # Local references, from the background column at this station height.
        # The band between them is evanescent, so power there is neither a
        # propagating gravity wave nor a propagating acoustic wave at this level.
        f_n = 1e3 * n_bg[izw] / (2 * np.pi)
        f_a = 1e3 * wa_bg[izw] / (2 * np.pi)
        # Cone of influence. Periods longer than coi(t) are contaminated by the
        # zero padding at the record ends, so the shaded wedge is not evidence.
        f_coi = 1e3 / np.maximum(coi, 1e-9)

        fig, ax = plt.subplots(2, 1, figsize=(9, 9), sharex=True)
        for k, (pp, cmap, lo, hi, cl, ttl) in enumerate((
                (p_abs, "magma", -50, 0, "dB below the record peak",
                 "absolute, one scale for the whole record"),
                (p_rel, "RdBu_r", -20, 20, "dB against the settled background",
                 f"relative to the median over t > {bg_start/3600:.1f} h, per "
                 f"frequency; blank below {WAVELET_REL_FLOOR_DB:.0f} dB absolute"))):
            cm = plt.get_cmap(cmap).copy()
            cm.set_bad("0.85")
            mesh = ax[k].pcolormesh(tmin, freq_mhz, np.ma.masked_invalid(pp),
                                    cmap=cm, vmin=lo, vmax=hi, shading="auto")
            ax[k].fill_between(tmin, 0, f_coi, color="0.5", alpha=0.35,
                               step="mid", lw=0)
            ax[k].plot(tmin, f_coi, color="0.5", lw=0.8)
            ax[k].axhline(f_n, color="k", ls=":", lw=1.1,
                          label=f"local N/2pi = {f_n:.2f} mHz")
            ax[k].axhline(f_a, color="lime", ls="--", lw=1.2,
                          label=f"local acoustic cutoff = {f_a:.2f} mHz")
            ax[k].axhline(1e3 * BUOYANCY_N_REF / (2 * np.pi), color="0.35",
                          ls="-.", lw=0.9, label="surface N/2pi = 2.87 mHz")
            ax[k].set_ylim(0, 10)
            ax[k].set_ylabel("frequency (mHz)")
            ax[k].set_title(ttl, fontsize=10)
            fig.colorbar(mesh, ax=ax[k], label=cl)
        ax[0].legend(loc="upper right", fontsize=8)
        ax[1].set_xlabel("time (min)")
        fig.suptitle(f"{label} full-record wavelet power (S3, {tag}); "
                     "shaded wedge is the cone of influence", fontsize=11)
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, f"deepdive_wavelet_w_full_{tag}.png"), dpi=150)
        plt.close(fig)

    # --- 3) |w|(z) profiles in the hill window at selected times (clean region)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 6), sharey=True)
    for k, h in enumerate(profile_hours):
        ax1.plot(wmax_profiles[k][jz_clean], z_km[jz_clean], label=f"t={h:.1f} h")
        ax2.plot(wrms_profiles[k][jz_clean], z_km[jz_clean], label=f"t={h:.1f} h")
    ax1.set_xlabel("hill-window max |w| (m/s)")
    ax2.set_xlabel("hill-window RMS w (m/s)")
    ax1.set_ylabel("z (km)")
    for ax in (ax1, ax2):
        ax.axhspan(*AMEC_BAND_KM, color="gray", alpha=0.2)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.suptitle(f"{label} vertical wave amplitude, |x_rel|<=500 km, z<60 km\n"
                 "(AMEC band 45-53 km shaded)")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "deepdive_w_profiles.png"), dpi=150)
    plt.close(fig)

    # --- 4) Momentum flux MF(z) and mean-flow change at selected times
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 6), sharey=True)
    for k, h in zip(it_profiles, profile_hours):
        ax1.plot(mf_zt[jz_clean, k] * 1e3, z_km[jz_clean], label=f"t={h:.1f} h")
        ax2.plot(ubar_zt[jz_clean, k] - u_bar, z_km[jz_clean], label=f"t={h:.1f} h")
    ax1.set_xlabel("MF = <rho u'w'>_x (mPa)")
    ax2.set_xlabel("<u>_x - U_bar (m/s)")
    ax1.set_ylabel("z (km)")
    for ax in (ax1, ax2):
        ax.axhspan(*AMEC_BAND_KM, color="gray", alpha=0.2)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.suptitle(f"{label} momentum flux and mean-flow change, z<60 km")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "deepdive_mf_ubar_profiles.png"), dpi=150)
    plt.close(fig)

    # --- 5) Convergence time series at fixed heights (+ optional baseline)
    cmp_series = None
    if args.compare_run:
        cmp_series = run_series(os.path.normpath(args.compare_run), u_bar)

    fig, axes = plt.subplots(3, 1, figsize=(10, 8.5), sharex=True)
    th = times / 3600.0
    axes[0].plot(th, amec_w_t, lw=1.2, color="C0", label=label)
    axes[0].set_ylabel("AMEC hill mean |w| (m/s)")
    for k, lv in enumerate(MF_LEVELS_KM):
        axes[1].plot(th, mf_zt[iz_mf[k]] * 1e3, lw=1.2, color=f"C{k}",
                     label=f"z={z_km[iz_mf[k]]:.1f} km")
    axes[1].set_ylabel("MF = <rho u'w'>_x (mPa)")
    axes[2].plot(th, du_max_t, lw=1.2, color="C0", label=label)
    axes[2].set_ylabel("peak <u>_x - U_bar,\nz<60 km (m/s)")
    axes[2].set_xlabel("time (h)")
    if cmp_series is not None:
        cth = cmp_series["times"] / 3600.0
        cmp_lbl = os.path.basename(os.path.dirname(
            os.path.normpath(args.compare_run)))
        axes[0].plot(cth, cmp_series["amec_w"], lw=1.0, ls="--", color="k",
                     alpha=0.7, label=f"{cmp_lbl} (baseline)")
        for k in range(len(MF_LEVELS_KM)):
            axes[1].plot(cth, cmp_series["mf"][k] * 1e3, lw=1.0, ls="--",
                         color=f"C{k}", alpha=0.6)
        axes[2].plot(cth, cmp_series["du_max"], lw=1.0, ls="--", color="k",
                     alpha=0.7, label=f"{cmp_lbl} (baseline)")
        for ax in axes:
            ax.axvline(cth[-1], color="gray", ls=":", lw=1.0)
    for ax in axes:
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8, loc="best")
    fig.suptitle(f"{label} convergence diagnostics (dotted line: baseline t_end)")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "deepdive_convergence.png"), dpi=150)
    plt.close(fig)

    # --- 6) Transient-timing metrics at S3, z=50 km
    s3 = w_low[WAVELET_STATION_IDX]
    late_mean = s3[times > times[-1] - 3600.0].mean()
    resid = np.abs(s3 - late_mean)
    peak = resid.max()
    below = resid < 0.10 * peak
    i_peak = int(np.argmax(resid))
    i_settle = next((i for i in range(i_peak, n_t) if below[i:].all()), n_t - 1)
    t_settle = float(times[i_settle])

    early = s3[times <= EARLY_WINDOW_S]
    early = early - early.mean()
    spec = np.abs(np.fft.rfft(early)) ** 2
    fft_freq = np.fft.rfftfreq(len(early), d=dt_snap)
    i_dom = 1 + int(np.argmax(spec[1:]))
    transient_period_s = float(1.0 / fft_freq[i_dom])

    # Trend windows: the last hour of the record, and the 2-3 h window that the
    # 3 h sweep would have used as its own final hour.
    windows = {"h2_3": (7200.0, 10800.0),
               "last_hour": (times[-1] - 3600.0, times[-1])}
    tracked = {"amec_hill_mean_abs_w": amec_w_t,
               "peak_mean_flow_dev": du_max_t}
    for k, lv in enumerate(MF_LEVELS_KM):
        tracked[f"mf_{lv:.0f}km"] = mf_zt[iz_mf[k]]
    convergence = {name: {wk: trend_per_hour(times, ser, *wv)
                          for wk, wv in windows.items()}
                   for name, ser in tracked.items()}

    i60 = int(np.argmax(np.abs(ubar_zt[jz_clean, -1] - u_bar)))
    metrics = {
        "run_dir": os.path.abspath(run_dir),
        "label": label,
        "U_bar_mps": u_bar,
        "inventory": {
            "n_snapshots": n_t,
            "variables": varnames,
            "grid": {"nx": int(len(x)), "nz": int(len(z)),
                     "dx_m": float(x[1] - x[0]), "dz_m": float(z[1] - z[0])},
            "t_first_s": float(times[0]),
            "t_last_s": float(times[-1]),
            "snapshot_cadence_median_s": dt_snap,
            "size_on_disk_MB": round(size_bytes / 1e6, 1),
        },
        "stations": {
            "x_rel_km": [float(x_rel_km[i]) for i in ix_st],
            "z_low_km": float(z_km[iz_low]),
            "z_high_km": float(z_km[iz_high]),
        },
        "transient": {
            "station": f"S{WAVELET_STATION_IDX + 1} (z_low)",
            "peak_w_mps": float(s3.max()),
            "min_w_mps": float(s3.min()),
            "dominant_period_s_first4000s": transient_period_s,
            "dominant_freq_mHz": 1e3 / transient_period_s,
            "fft_bin_width_mHz": 1e3 * float(fft_freq[1]),
            "buoyancy_period_2pi_over_Nref_s": 2 * np.pi / BUOYANCY_N_REF,
            "local_N_per_s": float(n_bg[iz_low]),
            "local_buoyancy_period_s": float(2 * np.pi / n_bg[iz_low]),
            "local_N_over_2pi_mHz": float(1e3 * n_bg[iz_low] / (2 * np.pi)),
            "local_acoustic_cutoff_mHz": float(1e3 * wa_bg[iz_low] / (2 * np.pi)),
            "surface_N_over_2pi_mHz": float(1e3 * n_bg[0] / (2 * np.pi)),
            "band_at_station": (
                "gravity" if 1e3 / transient_period_s
                < 1e3 * n_bg[iz_low] / (2 * np.pi) else
                ("acoustic" if 1e3 / transient_period_s
                 > 1e3 * wa_bg[iz_low] / (2 * np.pi) else "evanescent")),
            "settle_time_s_10pct_of_peak": t_settle,
            "late_mean_w_mps": float(late_mean),
        },
        "clean_region_z_lt_60km": {
            "final_hill_max_abs_w_mps": float(wmax_profiles[-1][jz_clean].max()),
            "z_of_final_max_km": float(z_km[jz_clean][int(np.argmax(wmax_profiles[-1][jz_clean]))]),
            "final_MF_min_mPa": float(mf_zt[jz_clean, -1].min() * 1e3),
            "final_max_mean_flow_change_mps": float((ubar_zt[jz_clean, -1] - u_bar)[i60]),
            "z_of_max_mean_flow_change_km": float(z_km[jz_clean][i60]),
            "final_amec_hill_mean_abs_w_mps": float(amec_w_t[-1]),
            "z_of_peak_mean_flow_dev_final_km": float(z_du_max_t[-1]),
        },
        "breaking_indicators_hourly": {
            "hours": [float(times[int(np.argmin(np.abs(times - h * 3600.0)))] / 3600.0)
                      for h in profile_hours],
            "unstable_area_fraction": [float(unstable_frac_t[int(np.argmin(np.abs(times - h * 3600.0)))])
                                       for h in profile_hours],
            "min_hill_u_mps": [float(umin_hill_t[int(np.argmin(np.abs(times - h * 3600.0)))])
                               for h in profile_hours],
            "max_hill_abs_w_mps": [float(wmax_hill_t[int(np.argmin(np.abs(times - h * 3600.0)))])
                                   for h in profile_hours],
            "min_zonal_mean_u_mps": [float(umean_min_t[int(np.argmin(np.abs(times - h * 3600.0)))])
                                     for h in profile_hours],
        },
        "convergence_trends": convergence,
        "compare_run": (os.path.abspath(args.compare_run)
                        if args.compare_run else None),
        "analysis_wall_clock_s": None,  # filled below
    }
    if cmp_series is not None:
        cmp_tracked = {"amec_hill_mean_abs_w": cmp_series["amec_w"],
                       "peak_mean_flow_dev": cmp_series["du_max"]}
        for k, lv in enumerate(MF_LEVELS_KM):
            cmp_tracked[f"mf_{lv:.0f}km"] = cmp_series["mf"][k]
        ct = cmp_series["times"]
        metrics["compare_run_final"] = {
            "t_last_s": float(ct[-1]),
            "values": {n: float(s[-1]) for n, s in cmp_tracked.items()},
            "trends_last_hour": {
                n: trend_per_hour(ct, s, ct[-1] - 3600.0, ct[-1])
                for n, s in cmp_tracked.items()},
        }
    metrics["final_values"] = {n: float(s[-1]) for n, s in tracked.items()}
    metrics["analysis_wall_clock_s"] = round(time.perf_counter() - t_wall0, 1)
    out_json = os.path.join(outdir, "deepdive_metrics.json")
    with open(out_json, "w") as fh:
        json.dump(metrics, fh, indent=2)
    print(json.dumps(metrics, indent=2))
    print(f"\nWrote figures + {out_json} in {metrics['analysis_wall_clock_s']} s")


if __name__ == "__main__":
    main()
