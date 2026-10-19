#!/usr/bin/env python
# =============================================================================
# Arsia Mons Stage 2 analysis (real reversing shear, critical levels)
#
# Stage 2 launches a stationary orographic gravity wave into the REAL reversing
# zonal wind u(z). For a stationary wave (phase speed c = 0), every height where
# u(z) = 0 is a critical level. This analyzer quantifies whether the wave
# deposits momentum near the 45-53 km AMEC cloud band, whose top (~53.3 km)
# coincides with the second u(z) zero-crossing.
#
# Core diagnostic: the vertical flux of horizontal momentum (Eliassen-Palm
# stress) MF(z) = mean_x[ rho * u' * w' ], with u' = u - mean_x(u),
# w' = w - mean_x(w). For a steady, non-dissipative wave MF(z) is constant with
# height; a drop in |MF| marks momentum deposition (critical-level absorption or
# breaking). Lessons carried from Stage 1: the orographic amplitude is reported
# as a hill-windowed mean (not a domain peak), and the periodic-x BC removes the
# inflow-edge artifact, so MF is taken as the full-domain zonal mean.
#
# Usage (datasci env; python is not on PATH on the dev box):
#   <datasci python> scripts/analyze_arsia_stage2.py \
#       --sweep-dir output/arsia_mons_stage2 \
#       --outdir output/arsia_mons_stage2_analysis
#   Add --width-mode for the width-sensitivity sweep
#   (--sweep-dir output/arsia_mons_stage2_width).
# =============================================================================

import argparse
import csv
import glob
import json
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset

# run_U<u>_h<h> or run_LX<lx>_U<u>_h<h>  (U may be negative)
RUN_DIR_PATTERN = re.compile(
    r"run_(?:LX(?P<lx>\d+)_)?U(?P<u>-?\d+(?:\.\d+)?)_h(?P<h>-?\d+(?:\.\d+)?)"
)


# ---------------------------------------------------------------------------
# Background atmosphere and critical levels
# ---------------------------------------------------------------------------
def load_background(atm_file):
    """Return (z_km, u_zonal, critical_levels_km) from the MCD CSV column 6."""
    data = np.genfromtxt(atm_file, delimiter=",", names=True)
    z_km = data["z_m"] / 1000.0
    u = data["u_zonal"]
    crit = []
    for i in range(len(u) - 1):
        a, b = u[i], u[i + 1]
        if a == 0.0:
            crit.append(z_km[i])
        elif (a < 0.0) != (b < 0.0):
            zc = z_km[i] + (0.0 - a) / (b - a) * (z_km[i + 1] - z_km[i])
            crit.append(zc)
    # Drop the trivial surface crossing at z=0 (forcing handles the launch).
    crit = [c for c in crit if c > 0.5]
    return z_km, u, crit


# ---------------------------------------------------------------------------
# Run discovery and snapshot loading
# ---------------------------------------------------------------------------
def find_runs(sweep_dir, pattern):
    runs = []
    for d in sorted(glob.glob(os.path.join(sweep_dir, pattern))):
        if not os.path.isdir(d):
            continue
        m = RUN_DIR_PATTERN.search(os.path.basename(d))
        if m is None:
            continue
        runs.append({
            "dir": d,
            "name": os.path.basename(d),
            "u_bar": float(m.group("u")),
            "h_max_km": float(m.group("h")),
            "lx_km": float(m.group("lx")) if m.group("lx") else None,
        })
    return runs


def _orient_zx(arr, nz, nx):
    """Return arr as (nz, nx) regardless of stored axis order."""
    if arr.shape == (nz, nx):
        return arr
    if arr.shape == (nx, nz):
        return arr.T
    raise ValueError("Unexpected field shape %s (nz=%d, nx=%d)" % (arr.shape, nz, nx))


def load_final_snapshot(run_dir):
    """Open only the lexicographically last snapshot (step-padded == temporal)."""
    files = sorted(glob.glob(os.path.join(run_dir, "*.nc")))
    if not files:
        return None
    with Dataset(files[-1], "r") as ds:
        x_m = np.asarray(ds.variables["x"][:], dtype=float)
        z_m = np.asarray(ds.variables["z"][:], dtype=float)
        nz, nx = len(z_m), len(x_m)
        u = _orient_zx(np.asarray(ds.variables["u"][:], dtype=float), nz, nx)
        w = _orient_zx(np.asarray(ds.variables["w"][:], dtype=float), nz, nx)
        rho = _orient_zx(np.asarray(ds.variables["rho"][:], dtype=float), nz, nx)
        t_last = float(getattr(ds, "time_s", np.nan))
    return {"x_m": x_m, "z_m": z_m, "u": u, "w": w, "rho": rho, "t_last": t_last}


# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------
def momentum_flux(u, w, rho):
    """Vertical flux of horizontal momentum MF(z) = mean_x[rho * u' * w'] [Pa]."""
    up = u - u.mean(axis=1, keepdims=True)
    wp = w - w.mean(axis=1, keepdims=True)
    return (rho * up * wp).mean(axis=1)


def nearest_index(z_km, target_km):
    return int(np.argmin(np.abs(z_km - target_km)))


def final_hour_mf(run_dir, stride=3):
    """MF(z) averaged over the final hour, from every `stride`-th snapshot.

    Returns (z_km, mf_mean, n_snapshots, t_last_s), or None without snapshots.
    """
    files = sorted(glob.glob(os.path.join(run_dir, "snapshot_*.nc")))
    if not files:
        return None
    with Dataset(files[-1], "r") as ds:
        t_last = float(getattr(ds, "time_s"))
    acc, n, z_km = None, 0, None
    for f in files[::-1][::stride]:
        with Dataset(f, "r") as ds:
            if float(getattr(ds, "time_s")) < t_last - 3600.0:
                break
            z_m = np.asarray(ds.variables["z"][:], dtype=float)
            nz, nx = len(z_m), len(ds.variables["x"][:])
            fld = {v: _orient_zx(np.asarray(ds.variables[v][:], dtype=float), nz, nx)
                   for v in ("u", "w", "rho")}
        mf = momentum_flux(fld["u"], fld["w"], fld["rho"])
        acc = mf if acc is None else acc + mf
        n += 1
        z_km = z_m / 1000.0
    return z_km, acc / n, n, t_last


# Levels of the Stage 2 flux table in docs/reflection_convergence_tests.md §7
COMPARE_LEVELS_KM = (1.25, 2.75, 10.25, 30.25, 50.25, 55.25)


def compare_sweeps(new_dir, ref_dir, pattern, outdir, bg, args):
    """Final-hour MF of matching runs in two sweeps (e.g. pbar_mode A/B)."""
    new_runs = {r["name"]: r for r in find_runs(new_dir, pattern)}
    ref_runs = {r["name"]: r for r in find_runs(ref_dir, pattern)}
    names = sorted(set(new_runs) & set(ref_runs))
    if not names:
        print("No matching runs between %s and %s" % (new_dir, ref_dir))
        return
    rows = []
    fig, axes = plt.subplots(1, len(names), figsize=(4.5 * len(names), 6),
                             squeeze=False)
    for ax, name in zip(axes[0], names):
        z_n, mf_n, n_n, t_n = final_hour_mf(new_runs[name]["dir"])
        z_r, mf_r, n_r, t_r = final_hour_mf(ref_runs[name]["dir"])
        lev = {}
        for zk in COMPARE_LEVELS_KM:
            a = float(mf_n[nearest_index(z_n, zk)])
            b = float(mf_r[nearest_index(z_r, zk)])
            lev["%.2f" % zk] = {"new_uPa": a * 1e6, "ref_uPa": b * 1e6,
                                "new_over_ref": a / b if b != 0 else float("nan")}
        absorb = lambda mf, z: abs(mf[nearest_index(z, 10.25)] / mf[nearest_index(z, 1.25)])
        rows.append({"run": name, "t_last_new_s": t_n, "t_last_ref_s": t_r,
                     "snapshots_new": n_n, "snapshots_ref": n_r, "levels": lev,
                     "absorption_ratio_new": float(absorb(mf_n, z_n)),
                     "absorption_ratio_ref": float(absorb(mf_r, z_r))})
        zm = z_n <= args.composite_z_max_km
        ax.plot(mf_r[zm] * 1e6, z_r[zm], color="0.5", ls="--", lw=1.5, label="reference")
        ax.plot(mf_n[zm] * 1e6, z_n[zm], color="tab:red", lw=1.5, label="new")
        ax.set_xscale("symlog", linthresh=1.0)
        ax.set_title(name)
        ax.set_xlabel("final-hour <rho u'w'>_x [uPa]")
        _mark_levels(ax, bg, args, profile=True)
        ax.set_ylim(0, args.composite_z_max_km)
        ax.legend(fontsize=8)
    axes[0][0].set_ylabel("z [km]")
    fig.tight_layout()
    png = os.path.join(outdir, "mf_compare.png")
    fig.savefig(png, dpi=130)
    plt.close(fig)
    with open(os.path.join(outdir, "mf_compare.json"), "w") as f:
        json.dump({"new_dir": os.path.abspath(new_dir),
                   "ref_dir": os.path.abspath(ref_dir), "runs": rows}, f, indent=1)
    print("\nFinal-hour MF [uPa], new / reference (every third snapshot):")
    print("| run | " + " | ".join("%.2f km" % z for z in COMPARE_LEVELS_KM) +
          " | |MF(10.25)/MF(1.25)| new, ref |")
    for r in rows:
        cells = ["%.4g / %.4g" % (v["new_uPa"], v["ref_uPa"]) for v in r["levels"].values()]
        print("| %s | %s | %.3f, %.3f |" % (r["run"], " | ".join(cells),
                                           r["absorption_ratio_new"],
                                           r["absorption_ratio_ref"]))
    print("Comparison -> %s, mf_compare.json" % png)


# ---------------------------------------------------------------------------
# Hourly time series of a long record (docs/reflection_convergence_tests.md §12)
# ---------------------------------------------------------------------------
TS_LEVELS_KM = (1.25, 10.25, 30.25, 50.25, 55.25)


def _crossing_near(z_km, u, target=53.3, lo=40.0, hi=70.0):
    """Height [km] of the u = 0 crossing nearest `target` within [lo, hi]."""
    m = (z_km >= lo) & (z_km <= hi)
    zz, uu = z_km[m], u[m]
    idx = np.where(np.sign(uu[:-1]) != np.sign(uu[1:]))[0]
    if len(idx) == 0:
        return float("nan")
    zc = [zz[i] - uu[i] * (zz[i + 1] - zz[i]) / (uu[i + 1] - uu[i]) for i in idx]
    return float(min(zc, key=lambda v: abs(v - target)))


def run_timeseries(run_dir, stride=1, band=(15.0, 45.0)):
    """Hour-mean MF(z) and background diagnostics of one record.

    Per hour: MF at TS_LEVELS_KM, the peak-to-peak of MF over `band` relative to its
    band mean, the x-mean u = 0 level near 53 km and the largest x-mean density
    departure below 45 km (last snapshot of the hour, against snapshot 0), the
    largest |u'| below 60 km and the fraction of points below 60 km with
    d(theta)/dz < 0.
    """
    files = sorted(glob.glob(os.path.join(run_dir, "snapshot_*.nc")))[::stride]
    if not files:
        return None
    with Dataset(files[0], "r") as ds:
        z_km = np.asarray(ds.variables["z"][:], dtype=float) / 1000.0
        rho0_x = np.asarray(ds.variables["rho"][:], dtype=float).mean(axis=1)
        gamma = float(getattr(ds, "gamma"))
    kappa = (gamma - 1.0) / gamma
    hours = {}
    for f in files:
        with Dataset(f, "r") as ds:
            t = float(getattr(ds, "time_s"))
            u = np.asarray(ds.variables["u"][:], dtype=float)
            w = np.asarray(ds.variables["w"][:], dtype=float)
            rho = np.asarray(ds.variables["rho"][:], dtype=float)
            T = np.asarray(ds.variables["T"][:], dtype=float)
            P = np.asarray(ds.variables["P"][:], dtype=float)
        k = int(t // 3600.0)
        h = hours.setdefault(k, {"mf": 0.0, "n": 0, "umax": 0.0, "unstable": 0.0})
        h["mf"] = h["mf"] + momentum_flux(u, w, rho)
        h["n"] += 1
        low = z_km < 60.0
        ubar = u.mean(axis=1, keepdims=True)
        h["umax"] = max(h["umax"], float(np.abs(u - ubar)[low].max()))
        theta = T * (P[0:1, :] / P) ** kappa
        h["unstable"] = max(h["unstable"],
                            float((np.diff(theta, axis=0)[low[:-1]] < 0).mean()))
        h["t_last"] = t
        h["crit_km"] = _crossing_near(z_km, ubar[:, 0])
        dep = np.abs(rho.mean(axis=1) / rho0_x - 1.0)[z_km < 45.0]
        h["drho_max_pct"] = float(100.0 * dep.max())
    rows = []
    bmask = (z_km >= band[0]) & (z_km <= band[1])
    for k in sorted(hours):
        h = hours[k]
        mf = h["mf"] / h["n"]
        mb = mf[bmask]
        row = {"t_h": k + 0.5, "n": h["n"],
               "ptp_band_pct": float(100.0 * (mb.max() - mb.min()) / abs(mb.mean()))
               if mb.mean() != 0 else float("nan"),
               "crit_km": h["crit_km"], "drho_max_below45_pct": h["drho_max_pct"],
               "max_abs_uprime_below60": h["umax"],
               "frac_dthetadz_neg_below60": h["unstable"]}
        for zk in TS_LEVELS_KM:
            row["MF_%.2f_uPa" % zk] = float(mf[nearest_index(z_km, zk)] * 1e6)
        row["MF30_over_absMF1"] = row["MF_30.25_uPa"] / abs(row["MF_1.25_uPa"])
        row["MF55_over_MF50"] = row["MF_55.25_uPa"] / row["MF_50.25_uPa"]
        rows.append(row)
    return {"run_dir": os.path.abspath(run_dir), "z_km": z_km.tolist(), "hours": rows}


def _window(rows, t_lo, t_hi, key):
    v = [r[key] for r in rows if t_lo < r["t_h"] < t_hi]
    return float(np.mean(v)) if v else float("nan")


def timeseries_report(runs, outdir, stride):
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=min(9, len(runs))) as ex:
        res = list(ex.map(run_timeseries, [r["dir"] for r in runs],
                          [stride] * len(runs)))
    out = []
    print("\nWindows: early 6-9 h, final 18-24 h (hour means)")
    print("| run | ptp 15-45 km early, final [%] | MF30/|MF1| early, final | final/early "
          "| MF55/MF50 final | u=0 level t0 -> end [km] | max drho<45 km [%] "
          "| max |u'| [m/s] | max frac dtheta/dz<0 |")
    for run, r in zip(runs, res):
        if r is None:
            continue
        hs = r["hours"]
        e30 = _window(hs, 6, 9, "MF30_over_absMF1")
        f30 = _window(hs, 18, 24, "MF30_over_absMF1")
        summ = {"run": run["name"], "t_end_h": hs[-1]["t_h"] + 0.5,
                "ptp_early": _window(hs, 6, 9, "ptp_band_pct"),
                "ptp_final": _window(hs, 18, 24, "ptp_band_pct"),
                "MF30_rel_early": e30, "MF30_rel_final": f30,
                "MF30_final_over_early": f30 / e30 if e30 else float("nan"),
                "MF55_over_MF50_final": _window(hs, 18, 24, "MF55_over_MF50"),
                "crit_first_km": hs[0]["crit_km"], "crit_last_km": hs[-1]["crit_km"],
                "drho_max_pct": max(x["drho_max_below45_pct"] for x in hs),
                "max_uprime": max(x["max_abs_uprime_below60"] for x in hs),
                "max_unstable_frac": max(x["frac_dthetadz_neg_below60"] for x in hs)}
        r["summary"] = summ
        out.append(r)
        print("| %s | %.1f, %.1f | %.3g, %.3g | %.2f | %.3f | %.2f -> %.2f | %.3f | %.2f | %.2e |"
              % (summ["run"], summ["ptp_early"], summ["ptp_final"], e30, f30,
                 summ["MF30_final_over_early"], summ["MF55_over_MF50_final"],
                 summ["crit_first_km"], summ["crit_last_km"], summ["drho_max_pct"],
                 summ["max_uprime"], summ["max_unstable_frac"]))
    with open(os.path.join(outdir, "timeseries.json"), "w") as f:
        json.dump({"stride": stride, "runs": out}, f, indent=1)
    fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=True)
    for r in out:
        hs = r["hours"]
        t = [x["t_h"] for x in hs]
        lab = r["summary"]["run"].replace("run_", "")
        axes[0].plot(t, [x["ptp_band_pct"] for x in hs], label=lab)
        axes[1].plot(t, [x["MF30_over_absMF1"] for x in hs], label=lab)
        axes[2].plot(t, [x["crit_km"] for x in hs], label=lab)
    axes[0].set_yscale("log")
    axes[0].axhline(10, color="0.5", ls=":")
    axes[0].set_ylabel("MF peak-to-peak 15-45 km [%]")
    axes[1].set_ylabel("MF(30.25) / |MF(1.25)|")
    axes[1].set_yscale("symlog", linthresh=1e-3)
    axes[2].set_ylabel("x-mean u = 0 level near 53 km [km]")
    axes[2].set_xlabel("t [h]")
    for a in axes:
        a.axvspan(6, 9, color="tab:blue", alpha=0.08)
        a.axvspan(18, 24, color="tab:red", alpha=0.08)
    axes[0].legend(fontsize=7, ncol=3)
    fig.tight_layout()
    png = os.path.join(outdir, "timeseries.png")
    fig.savefig(png, dpi=130)
    plt.close(fig)
    print("Time series -> %s, timeseries.json" % png)


def load_durations(sweep_dir):
    """Read sweep_results.csv into a full key (U, h, LX) map and a (U, h)
    fallback map. The main sweep names dirs without LX (run_U<u>_h<h>) while the
    CSV always records LX_km, so the fallback is needed to join wall-clock."""
    path = os.path.join(sweep_dir, "sweep_results.csv")
    by_full, by_uh = {}, {}
    if not os.path.isfile(path):
        return by_full, by_uh
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            try:
                u = round(float(row["U_bar"]), 3)
                h = round(float(row["h_max"]) / 1000.0, 3)
                lx = round(float(row["LX_km"]), 1) if row.get("LX_km") else None
            except (KeyError, ValueError):
                continue
            dur = float(row.get("duration_s", "nan"))
            by_full[(u, h, lx)] = dur
            by_uh[(u, h)] = dur
    return by_full, by_uh


# ---------------------------------------------------------------------------
# Per-run figure
# ---------------------------------------------------------------------------
def plot_run(run, snap, bg, args, outdir):
    z_km = snap["z_m"] / 1000.0
    x_rel_km = (snap["x_m"] - 0.5 * (snap["x_m"][0] + snap["x_m"][-1])) / 1000.0
    w = snap["w"]
    mf = momentum_flux(snap["u"], snap["w"], snap["rho"])

    zmask = z_km <= args.composite_z_max_km
    xmask = np.abs(x_rel_km) <= args.composite_x_limit_km

    fig, (axw, axp) = plt.subplots(
        1, 2, figsize=(13, 6), gridspec_kw={"width_ratios": [2.2, 1.0]}
    )

    # --- w field ---
    wsub = w[np.ix_(zmask, xmask)]
    vmax = max(np.quantile(np.abs(wsub), 0.995), 1e-6)
    levels = np.linspace(-vmax, vmax, 41)
    im = axw.contourf(
        x_rel_km[xmask], z_km[zmask], wsub,
        levels=levels, cmap="RdBu_r", extend="both",
    )
    _mark_levels(axw, bg, args, x_text=x_rel_km[xmask].min())
    axw.set_xlabel("x - x_hill [km]")
    axw.set_ylabel("altitude [km]")
    axw.set_title("w [m/s]   %s   (t=%.0f s)" % (run["name"], snap["t_last"]))
    fig.colorbar(im, ax=axw, shrink=0.85, label="w [m/s]")

    # --- profiles: u(z) and MF(z) sharing the z-axis ---
    axp.plot(bg["u_interp"](z_km), z_km, color="tab:gray", lw=1.6, label="u_zonal")
    axp.axvline(0.0, color="k", lw=0.6)
    axp.set_xlabel("u_zonal [m/s]", color="tab:gray")
    axp.tick_params(axis="x", labelcolor="tab:gray")
    axp.set_ylim(0, args.composite_z_max_km)
    axp.set_ylabel("altitude [km]")

    axmf = axp.twiny()
    axmf.plot(mf[zmask], z_km[zmask], color="tab:red", lw=1.8, label="MF")
    axmf.set_xlabel("MF = <rho u' w'> [Pa]", color="tab:red")
    axmf.tick_params(axis="x", labelcolor="tab:red")
    _mark_levels(axp, bg, args, x_text=None, profile=True)
    axp.set_title("background wind and momentum flux")

    fig.tight_layout()
    path = os.path.join(outdir, "stage2_%s.png" % run["name"])
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def _mark_levels(ax, bg, args, x_text=None, profile=False):
    ax.axhspan(args.band_z_min_km, args.band_z_max_km, color="tab:green", alpha=0.12)
    ax.axhline(args.sponge_base_km, color="0.4", ls=":", lw=1.0)
    for c in bg["critical_levels"]:
        ax.axhline(c, color="tab:purple", ls="--", lw=1.0)
        if x_text is not None:
            ax.text(x_text, c, " crit %.1f km" % c, color="tab:purple",
                    va="bottom", ha="left", fontsize=7)


# ---------------------------------------------------------------------------
# Sweep-level composite of w fields
# ---------------------------------------------------------------------------
def plot_composite(runs_data, bg, args, outdir):
    u_vals = sorted({rd["run"]["u_bar"] for rd in runs_data})
    h_vals = sorted({rd["run"]["h_max_km"] for rd in runs_data})
    lx_vals = sorted({rd["run"]["lx_km"] for rd in runs_data if rd["run"]["lx_km"]})

    # Width-mode: 1 row of widths. Otherwise rows=h (largest on top), cols=U.
    if len(lx_vals) > 1:
        rows, cols = [None], lx_vals
        row_key, col_key = "h_max_km", "lx_km"
    else:
        rows, cols = list(reversed(h_vals)), u_vals
        row_key, col_key = "h_max_km", "u_bar"

    index = {}
    for rd in runs_data:
        index[(rd["run"][row_key], rd["run"][col_key])] = rd

    # Shared color scale on density-scaled w (removes trivial rho^-1/2 growth).
    pooled = []
    for rd in runs_data:
        pooled.append(np.abs(rd["w_scaled_sub"]).ravel())
    vmax = max(np.quantile(np.concatenate(pooled), 0.995), 1e-9)
    levels = np.linspace(-vmax, vmax, 41)

    nrows, ncols = len(rows), len(cols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.0 * ncols, 3.4 * nrows),
                             squeeze=False, sharex=True, sharey=True)
    im = None
    for r, rk in enumerate(rows):
        for c, ck in enumerate(cols):
            ax = axes[r][c]
            key = (rk if rk is not None else runs_data[0]["run"][row_key], ck)
            rd = index.get(key)
            if rd is None:
                ax.set_axis_off()
                continue
            im = ax.contourf(rd["x_rel_km"], rd["z_km"], rd["w_scaled_sub"],
                             levels=levels, cmap="RdBu_r", extend="both")
            _mark_levels(ax, bg, args, x_text=None)
            ax.set_ylim(0, args.composite_z_max_km)
            if len(lx_vals) > 1:
                ax.set_title("LX=%d km" % int(ck))
            else:
                ax.set_title("U=%g  h=%g km  (Nh/U=%.2f)"
                             % (ck, rk, bg["n_ref"] * rk * 1000.0 / abs(ck)))
            if r == nrows - 1:
                ax.set_xlabel("x - x_hill [km]")
            if c == 0:
                ax.set_ylabel("altitude [km]")

    if im is not None:
        fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.85,
                     label="w * sqrt(rho/rho_surface) [m/s]")
    fig.suptitle("Stage 2 composite: density-scaled w with critical levels", y=1.0)
    path = os.path.join(outdir, "stage2_composite_w.png")
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_mf_overlay(runs_data, bg, args, outdir):
    fig, ax = plt.subplots(figsize=(6.5, 7))
    for rd in runs_data:
        ax.plot(rd["mf"], rd["z_km_full"], lw=1.4, label=rd["run"]["name"])
    ax.axvline(0.0, color="k", lw=0.6)
    _mark_levels(ax, bg, args, x_text=None, profile=True)
    ax.set_ylim(0, args.composite_z_max_km)
    ax.set_xlabel("MF = <rho u' w'> [Pa]")
    ax.set_ylabel("altitude [km]")
    ax.set_title("Momentum-flux profiles (drop = deposition)")
    ax.legend(fontsize=7, loc="upper right")
    fig.tight_layout()
    path = os.path.join(outdir, "stage2_mf_overlay.png")
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# Response table
# ---------------------------------------------------------------------------
def write_response(rows, outdir):
    cols = ["run_name", "U_bar", "h_max_km", "LX_km", "Nh_over_U",
            "peak_w_hill_mps", "mean_w_hill_mps", "peak_w_global_mps",
            "mf_launch_Pa", "mf_amec_Pa", "mf_deposited_frac",
            "t_last_s", "duration_s"]
    csv_path = os.path.join(outdir, "stage2_response.csv")
    with open(csv_path, "w", newline="") as f:
        wtr = csv.writer(f)
        wtr.writerow(cols)
        for r in rows:
            wtr.writerow([r[c] for c in cols])

    md_path = os.path.join(outdir, "stage2_response.md")
    with open(md_path, "w") as f:
        f.write("# Arsia Mons Stage 2 AMEC-band response\n\n")
        f.write("Hill metrics use the |x_rel| <= %d km window; MF is the "
                "full-domain zonal mean. mf_deposited_frac = 1 - "
                "|MF(~52 km)| / |MF(~10 km)| (fraction of upward momentum flux "
                "absorbed below the AMEC band).\n\n" % rows[0]["x_window_km"])
        hdr = ("| run | U | h (km) | LX (km) | Nh/U | hill peak |w| | "
               "hill mean |w| | global peak |w| | MF launch | MF AMEC | "
               "deposited frac | wall (s) |")
        f.write(hdr + "\n")
        f.write("|" + "---|" * 12 + "\n")
        for r in rows:
            f.write("| %s | %g | %g | %s | %.2f | %.3f | %.3f | %.3f | "
                    "%.3e | %.3e | %.2f | %.0f |\n" % (
                        r["run_name"], r["U_bar"], r["h_max_km"],
                        ("%g" % r["LX_km"]) if r["LX_km"] else "-",
                        r["Nh_over_U"], r["peak_w_hill_mps"], r["mean_w_hill_mps"],
                        r["peak_w_global_mps"], r["mf_launch_Pa"], r["mf_amec_Pa"],
                        r["mf_deposited_frac"], r["duration_s"]))
    return csv_path, md_path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description="Arsia Mons Stage 2 analysis")
    p.add_argument("--sweep-dir", default="output/arsia_mons_stage2")
    p.add_argument("--run-pattern", default="run_*U*_h*")
    p.add_argument("--outdir", default="output/arsia_mons_stage2_analysis")
    p.add_argument("--atm-file",
                   default="external/mcd_profiles/mars_arsia_mons_my31.csv")
    p.add_argument("--band-z-min-km", type=float, default=45.0)
    p.add_argument("--band-z-max-km", type=float, default=53.0)
    p.add_argument("--composite-z-max-km", type=float, default=90.0)
    p.add_argument("--composite-x-limit-km", type=float, default=500.0)
    p.add_argument("--x-window-km", type=float, default=500.0,
                   help="hill window half-width for the windowed |w| metric")
    p.add_argument("--sponge-base-km", type=float, default=60.0)
    p.add_argument("--mf-launch-km", type=float, default=10.0)
    p.add_argument("--mf-amec-km", type=float, default=52.0)
    p.add_argument("--n-ref", type=float, default=0.018)
    p.add_argument("--width-mode", action="store_true",
                   help="defaults --sweep-dir to the width-sensitivity sweep")
    p.add_argument("--timeseries", action="store_true",
                   help="hourly MF and background diagnostics of long records "
                        "(reflection_convergence_tests.md §12)")
    p.add_argument("--ts-stride", type=int, default=1,
                   help="use every n-th snapshot in --timeseries")
    p.add_argument("--compare-dir", default=None,
                   help="reference sweep; adds a final-hour MF table of the runs "
                        "common to both (e.g. the pbar_mode A/B)")
    return p.parse_args()


def main():
    args = parse_args()
    if args.width_mode and args.sweep_dir == "output/arsia_mons_stage2":
        args.sweep_dir = "output/arsia_mons_stage2_width"
    os.makedirs(args.outdir, exist_ok=True)

    z_km_bg, u_bg, crit = load_background(args.atm_file)
    bg = {
        "z_km": z_km_bg, "u_zonal": u_bg, "critical_levels": crit,
        "n_ref": args.n_ref,
        "u_interp": lambda zq: np.interp(zq, z_km_bg, u_bg),
    }
    print("Background: U(0)=%.2f m/s, critical levels (km): %s"
          % (u_bg[0], ", ".join("%.1f" % c for c in crit)))

    runs = find_runs(args.sweep_dir, args.run_pattern)
    if not runs:
        print("No runs found under %s matching %s"
              % (args.sweep_dir, args.run_pattern))
        return
    print("Found %d run(s) in %s" % (len(runs), args.sweep_dir))

    durations_full, durations_uh = load_durations(args.sweep_dir)
    runs_data, rows = [], []
    for run in runs:
        snap = load_final_snapshot(run["dir"])
        if snap is None:
            print("  %s: no snapshots, skipping" % run["name"])
            continue
        z_km = snap["z_m"] / 1000.0
        x_rel_km = (snap["x_m"]
                    - 0.5 * (snap["x_m"][0] + snap["x_m"][-1])) / 1000.0
        mf = momentum_flux(snap["u"], snap["w"], snap["rho"])

        # AMEC-band amplitude (hill-windowed, lesson from Stage 1)
        band = (z_km >= args.band_z_min_km) & (z_km <= args.band_z_max_km)
        win = np.abs(x_rel_km) <= args.x_window_km
        wb = snap["w"][band, :]
        peak_w_hill = float(np.abs(wb[:, win]).max())
        mean_w_hill = float(np.abs(wb[:, win]).mean())
        peak_w_global = float(np.abs(wb).max())

        # Momentum deposition below the AMEC band
        mf_launch = float(mf[nearest_index(z_km, args.mf_launch_km)])
        mf_amec = float(mf[nearest_index(z_km, args.mf_amec_km)])
        dep_frac = (1.0 - abs(mf_amec) / abs(mf_launch)
                    if abs(mf_launch) > 1e-30 else float("nan"))

        lx_key = round(run["lx_km"], 1) if run["lx_km"] else None
        uh_key = (round(run["u_bar"], 3), round(run["h_max_km"], 3))
        dur = (durations_full.get(uh_key + (lx_key,), float("nan"))
               if lx_key is not None
               else durations_uh.get(uh_key, float("nan")))

        # density-scaled w for the composite
        rho_surface = float(np.mean(snap["rho"][0, :]))
        w_scaled = snap["w"] * np.sqrt(snap["rho"] / rho_surface)
        zmask = z_km <= args.composite_z_max_km
        xmask = np.abs(x_rel_km) <= args.composite_x_limit_km

        runs_data.append({
            "run": run, "z_km": z_km[zmask], "z_km_full": z_km,
            "x_rel_km": x_rel_km[xmask],
            "w_scaled_sub": w_scaled[np.ix_(zmask, xmask)], "mf": mf,
        })
        rows.append({
            "run_name": run["name"], "U_bar": run["u_bar"],
            "h_max_km": run["h_max_km"], "LX_km": run["lx_km"],
            "Nh_over_U": args.n_ref * run["h_max_km"] * 1000.0 / abs(run["u_bar"]),
            "peak_w_hill_mps": peak_w_hill, "mean_w_hill_mps": mean_w_hill,
            "peak_w_global_mps": peak_w_global, "mf_launch_Pa": mf_launch,
            "mf_amec_Pa": mf_amec, "mf_deposited_frac": dep_frac,
            "t_last_s": snap["t_last"], "duration_s": dur,
            "x_window_km": args.x_window_km,
        })
        plot_run(run, snap, bg, args, args.outdir)
        print("  %s: hill mean|w|=%.3f m/s, MF dep frac=%.2f"
              % (run["name"], mean_w_hill, dep_frac))

    if not runs_data:
        print("No snapshots loaded; nothing to plot.")
        return

    comp = plot_composite(runs_data, bg, args, args.outdir)
    mfo = plot_mf_overlay(runs_data, bg, args, args.outdir)
    csv_path, md_path = write_response(rows, args.outdir)
    print("Composite      -> %s" % comp)
    print("MF overlay     -> %s" % mfo)
    print("Response table -> %s" % csv_path)
    print("               -> %s" % md_path)
    if args.timeseries:
        timeseries_report(runs, args.outdir, args.ts_stride)
    if args.compare_dir:
        compare_sweeps(args.sweep_dir, args.compare_dir, args.run_pattern,
                       args.outdir, bg, args)
    print("Outputs in: %s" % args.outdir)


if __name__ == "__main__":
    main()
