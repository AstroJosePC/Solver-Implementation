#!/usr/bin/env python
# =============================================================================
# Section 14 analysis: time-step convergence of the Stage 2 momentum flux and the
# launch-layer hypothesis (docs/reflection_convergence_tests.md Section 14).
#
# Subcommands (all read existing snapshots; none runs the solver):
#   suite    score the Section 14 runs against the pre-registered readings
#   zustar   Z4: the dissipation-controlled wind band U* from the MCD column
#   zprofile Z1: ratio of window-mean MF(z) between two records at different dt
#   zalt     2 dz alternation of the near-ground profile (Section 14 addendum 2)
#   ztheory  cavity theory level under several N(z) treatments (addendum 2)
#   zprobe   ground-face w of the z-predictor against the forced w (addendum 2)
#   zflin    F_lin of the uniform-wind Gaussian ridge, factor by factor, and F/F_lin (addendum 3)
#   zequil   H = 0 equilibrium: drho_L2 from the log, energy and mass drift (addendum 3)
#   zqueney  Queney case, face against ghost path, and the Srivastava Table 11 metric (addendum 3)
#   zdensity Z2: growth law of the near-ground density departure
#   zh3      Z3: precursors of the h = 3 km collapses (ground versus top, spectra)
#   zoverturn first convective and kinematic overturning by band, Lindzen height (Section 15)
#   zblend   Z5: Richardson-number gate on the ublend background before S14E/S14F
#   bitcheck the patched driver reproduces an archived record bit for bit
#
# Usage (python3 with netCDF4):
#   python scripts/analyze_dt_convergence.py suite \
#       --run S04=output/arsia_mons_stage2_24h_K0_0.005_LZ95_pbarP_ramp2h/run_U40.0_h0.1 \
#       --run S14A=output/arsia_mons_stage2_9h_K0_0.005_LZ95_pbarP_ramp2h_cfl0.2/run_U40.0_h0.1 \
#       --run S14B=... --run S14D=... --run S14E=<..._ub8>/run_U-40.0_h0.1 --run S14F=<..._ub8_cfl0.2>/run_U-40.0_h0.1 \
#       --run C2R=output/arsia_mons_linear18h_LX2000_mode6_K0_0.005_LZ95_pbarP_ramp2h/run_U40.0_h0.005 \
#       --run S14C=... --outdir output/reflection_tests/s14
#   python scripts/analyze_dt_convergence.py zustar --dt 0.659 --dz 500 --lx-km 4000
#   python scripts/analyze_dt_convergence.py zprofile --a <LZ95 run> --b <LZ120 run> --t0 18 --t1 24
#   python scripts/analyze_dt_convergence.py zdensity --run <run> [--run <run> ...]
#   python scripts/analyze_dt_convergence.py zh3 --run <U40 h3.0 LZ95 run> --stride 1
#
# Conventions: snapshots store fields as [z, x] after _orient_zx (asserted against the
# grid); times in seconds unless suffixed; MF in Pa (tables in uPa), a domain x-mean,
# so it scales as 1/LX and is only compared between runs of equal LX.
# =============================================================================

import argparse
import glob
import json
import os
import re
import sys

import numpy as np
from netCDF4 import Dataset

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_arsia_stage2 import (  # noqa: E402
    TS_LEVELS_KM, _orient_zx, _window, momentum_flux, nearest_index, run_timeseries)

GAMMA = 1.31      # solver value (mars_config)
R_GAS = 188.9     # J kg^-1 K^-1, solver value
G_MARS = 3.72     # m s^-2, solver value
STEP_RE = re.compile(r"snapshot_(\d+)\.nc$")


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------
def snapshot_files(run_dir):
    files = sorted(glob.glob(os.path.join(run_dir, "snapshot_*.nc")))
    if not files:
        raise FileNotFoundError("no snapshots in %s" % run_dir)
    return files


def mean_dt(run_dir):
    """Mean time step [s] from the first and last snapshot (step index in the name)."""
    files = snapshot_files(run_dir)
    s0, s1 = (int(STEP_RE.search(f).group(1)) for f in (files[0], files[-1]))
    with Dataset(files[0]) as a, Dataset(files[-1]) as b:
        t0, t1 = float(a.time_s), float(b.time_s)
    return (t1 - t0) / (s1 - s0) if s1 > s0 else float("nan")


def load_fields(path, names=("u", "w", "rho")):
    """Return {z_km [nz], x_m [nx], t_s, gamma, <name>: [nz, nx]} from one snapshot."""
    with Dataset(path) as ds:
        z_km = np.asarray(ds.variables["z"][:], float) / 1e3
        x_m = np.asarray(ds.variables["x"][:], float)
        out = {"z_km": z_km, "x_m": x_m, "t_s": float(ds.time_s),
               "gamma": float(getattr(ds, "gamma", GAMMA))}
        for n in names:
            a = _orient_zx(np.asarray(ds.variables[n][:], float), len(z_km), len(x_m))
            assert a.shape == (len(z_km), len(x_m))
            out[n] = a
    return out


def window_metrics(rows, t_lo, t_hi):
    """Window means of MF at the TS levels [uPa] and of the normalized ratios."""
    m = {"MF_%.2f" % z: _window(rows, t_lo, t_hi, "MF_%.2f_uPa" % z) for z in TS_LEVELS_KM}
    m["R10_30"] = m["MF_10.25"] / m["MF_30.25"]
    m["R50_30"] = m["MF_50.25"] / m["MF_30.25"]
    m["R55_50"] = m["MF_55.25"] / m["MF_50.25"]
    m["n_hours"] = len([r for r in rows if t_lo < r["t_h"] < t_hi])
    return m


def rel(a, b):
    """Relative change of b against a."""
    return (b - a) / abs(a) if a else float("nan")


# ---------------------------------------------------------------------------
# suite: the pre-registered readings of Section 14
# ---------------------------------------------------------------------------
def cmd_suite(args):
    runs = dict(kv.split("=", 1) for kv in args.run)
    os.makedirs(args.outdir, exist_ok=True)
    res = {}
    for lab, d in runs.items():
        ts = run_timeseries(d, stride=args.stride)
        # uniform-wind runs (C2R, S14C, the S14G addendum runs and the S14K, S14X and
        # S14R runs of addendum 2) read over 12-17 h; the MCD runs over 6-9 h
        lo, hi = ((12, 17) if lab in ("C2R", "S14C") or lab.startswith(("S14G", "S14K", "S14X", "S14R"))
                  else (6, 9))
        if args.window:
            lo, hi = args.window
        res[lab] = {"dir": os.path.abspath(d), "mean_dt_s": mean_dt(d),
                    "window_h": [lo, hi], "w": window_metrics(ts["hours"], lo, hi),
                    "hours": ts["hours"]}
        w = res[lab]["w"]
        print("%-5s dt=%.4f s  MF1.25=%9.3f MF10=%8.4f MF30=%8.4f MF50=%8.4f uPa  "
              "R10/30=%.3f R50/30=%.3f R55/50=%.3f (%d h)"
              % (lab, res[lab]["mean_dt_s"], w["MF_1.25"], w["MF_10.25"], w["MF_30.25"],
                 w["MF_50.25"], w["R10_30"], w["R50_30"], w["R55_50"], w["n_hours"]))

    verdict = {}

    def pair(a, b, keys):
        if a in res and b in res:
            return {k: rel(res[a]["w"][k], res[b]["w"][k]) for k in keys}
        return None

    # P1, primary criterion as proposed in Section 13: < 5 % between CFL 0.2 and 0.1.
    p = pair("S14A", "S14B", ["MF_30.25", "MF_1.25"])
    if p:
        verdict["P1_primary_cfl0.2_vs_0.1"] = {
            "changes": p, "pass": all(abs(v) < 0.05 for v in p.values()),
            "predicted_under_H_L": "fail, MF_30.25 change >= +15 %"}
    # P2, normalized propagation between CFL 0.2 and 0.1.
    p = pair("S14A", "S14B", ["R10_30", "R50_30", "R55_50"])
    if p:
        verdict["P2_normalized_cfl0.2_vs_0.1"] = {
            "changes": p,
            "pass": abs(p["R10_30"]) < 0.05 and abs(p["R50_30"]) < 0.05 and abs(p["R55_50"]) < 0.10,
            "predicted_under_H_L": "pass"}
    # P3, order diagnostic L = ln G2 / ln G1 over the dt sequence 0.4 -> 0.2 -> 0.1.
    if all(k in res for k in ("S04", "S14A", "S14B")):
        out = {}
        for key in ("MF_30.25", "MF_1.25", "MF_50.25"):
            a, b, c = (res[k]["w"][key] for k in ("S04", "S14A", "S14B"))
            g1, g2 = b / a, c / b
            out[key] = {"G1": g1, "G2": g2,
                        "L": float(np.log(abs(g2)) / np.log(abs(g1)))
                        if g1 > 0 and g2 > 0 and abs(np.log(g1)) > 1e-3 else float("nan")}
        verdict["P3_order"] = {"values": out, "reading": "L~1: no finite dt->0 limit; "
                               "L~0.5: first order, extrapolable; L~0.25: second order"}
    # P4, uniform-wind corner at two time steps (window 12-17 h).
    p = pair("C2R", "S14C", ["MF_10.25", "MF_30.25", "MF_50.25"])
    if p:
        verdict["P4_C2R_cfl0.4_vs_0.2"] = {
            "changes": p, "pass_H_L": all(abs(v) < 0.01 for v in p.values()),
            "intrinsic_solver_dt_dependence": any(abs(v) > 0.02 for v in p.values()),
            "note": "level and ripple share from compare_sponge_reflection.py, C2R vs S14C"}
    # P5, launch-consistent background at two time steps.
    p = pair("S14E", "S14F", ["MF_1.25", "MF_10.25", "MF_30.25", "MF_50.25"])
    if p:
        verdict["P5_ublend_cfl0.4_vs_0.2"] = {
            "changes": p, "pass": abs(p["MF_30.25"]) < 0.05 and abs(p["MF_10.25"]) < 0.05,
            "predicted_under_H_L": "pass", "predicted_under_H_S": "fail"}
    # P6, resolution at fixed dt (S14A dz 500 against S14D dz 250, both dt ~0.33 s).
    p = pair("S14A", "S14D", ["MF_30.25", "MF_1.25", "R50_30", "R10_30"])
    if p:
        verdict["P6_dz500_vs_dz250_fixed_dt"] = {
            "changes": p,
            "dt_match": rel(res["S14A"]["mean_dt_s"], res["S14D"]["mean_dt_s"]),
            "normalized_pass": abs(p["R50_30"]) < 0.05 and abs(p["R10_30"]) < 0.05,
            "predicted_under_H_L": "absolute changes >= 10 %, normalized pass"}
    # Section 15: any pair A:B, change of B against A and the plain ratio B/A
    for pr in args.pair or []:
        a, b = pr.split(":")
        keys = ["MF_%.2f" % z for z in TS_LEVELS_KM] + ["R10_30", "R50_30", "R55_50"]
        p = pair(a, b, keys)
        if p:
            verdict["pair_%s_%s" % (a, b)] = {
                "window_h": [res[a]["window_h"], res[b]["window_h"]], "changes": p,
                "ratios": {k: res[b]["w"][k] / res[a]["w"][k] for k in keys}}
    for k, v in verdict.items():
        print("\n%s\n  %s" % (k, json.dumps(v, indent=None, default=float)))
    with open(os.path.join(args.outdir, "s14_suite.json"), "w") as f:
        json.dump({"runs": res, "verdict": verdict}, f, indent=1, default=float)
    print("\n-> %s" % os.path.join(args.outdir, "s14_suite.json"))


# ---------------------------------------------------------------------------
# Z4: dissipation-controlled wind band from the MCD column
# ---------------------------------------------------------------------------
def ustar_profile(atm_file, dt, dz, k):
    """U*(z) below which the vertical Lax-Wendroff dissipation absorbs a stationary mode.

    For a hydrostatic stationary mode (m = N/|U|) the two-step LW sweep adds, per unit time,
    -kappa4 d^4/dz^4 on the acoustic pair (w, p'/rho0 c), kappa4 = c^2 dt dz^2 / 8. Weighting
    by the mode's energy fractions (k^2/m^2 in w, U^2/c^2 in p') and dividing by
    |c_gz| = k U^2 / N gives the flux attenuation per metre
        2 gamma / c_gz = dt dz^2 N^5 (1 + c^2 k^2 / N^2) / (8 k U^4),
    and integrating toward a critical level with |U| = |U'| s gives optical depth 1 at
        U* = [dt dz^2 N^5 (1 + c^2 k^2/N^2) / (24 k |U'|)]^(1/3).
    Validity: WKB, which fails within about |U|/|U'| < 1/N of a critical level and wherever
    Ri = N^2/U'^2 is O(1); U* is an order-of-magnitude band, not a sharp edge.

    Parameters
    ----------
    atm_file : str, MCD CSV with columns z_m, rho, P, T, nu, u_zonal (SI units)
    dt : float [s];  dz : float [m];  k : float [m^-1], horizontal wavenumber
    Returns dict of arrays over the CSV levels: z_km, u, N, dudz, ustar (m/s).
    """
    d = np.genfromtxt(atm_file, delimiter=",", names=True)
    z, T, u = d["z_m"], d["T"], d["u_zonal"]
    cp = GAMMA * R_GAS / (GAMMA - 1.0)
    n2 = G_MARS / T * (np.gradient(T, z) + G_MARS / cp)
    N = np.sqrt(np.clip(n2, 1e-12, None))
    c = np.sqrt(GAMMA * R_GAS * T)
    dudz = np.abs(np.gradient(u, z)) + 1e-12
    ustar = (dt * dz ** 2 * N ** 5 * (1 + (c * k / N) ** 2) / (24 * k * dudz)) ** (1 / 3)
    return {"z_km": z / 1e3, "u": u, "N": N, "dudz": dudz, "ustar": ustar}


def cmd_zustar(args):
    out = {}
    for n in args.modes:
        k = 2 * np.pi * n / (args.lx_km * 1e3)
        p = ustar_profile(args.atm_file, args.dt, args.dz, k)
        z, u, us = p["z_km"], p["u"], p["ustar"]
        low = z <= args.z_crit_low
        near53 = (z > 45) & (z < 60)
        band = low & (np.abs(u) < us)
        i53 = np.argmin(np.abs(z - 53.3))
        out[n] = {"k_m-1": k, "max_abs_u_below_%.1f_km" % args.z_crit_low: float(np.abs(u[low]).max()),
                  "ustar_median_below": float(np.median(us[low])),
                  "controlled_fraction_of_launch_layer": float(band.sum() / low.sum()),
                  "ustar_at_53.3_km": float(us[i53]),
                  "controlled_depth_at_53.3_km_m": float(us[i53] / p["dudz"][i53]),
                  "levels_below_ustar_45_60_km": int((near53 & (np.abs(u) < us)).sum())}
        print("n=%2d lambda_x=%5.0f km | launch layer: max|u| %.2f m/s, U* median %.2f m/s, "
              "controlled fraction %.2f | 53.3 km: U* %.2f m/s over %.0f m"
              % (n, args.lx_km / n, out[n]["max_abs_u_below_%.1f_km" % args.z_crit_low],
                 out[n]["ustar_median_below"], out[n]["controlled_fraction_of_launch_layer"],
                 out[n]["ustar_at_53.3_km"], out[n]["controlled_depth_at_53.3_km_m"]))
    if args.out:
        with open(args.out, "w") as f:
            json.dump({"dt_s": args.dt, "dz_m": args.dz, "lx_km": args.lx_km, "modes": out}, f,
                      indent=1)


# ---------------------------------------------------------------------------
# Z1: ratio of window-mean MF(z) between two records
# ---------------------------------------------------------------------------
def window_profile(run_dir, t0_h, t1_h, stride=1):
    """Window-mean MF(z) [Pa] over t0_h < t < t1_h; returns (z_km [nz], mf [nz], n)."""
    acc, n, z_km = None, 0, None
    for f in snapshot_files(run_dir)[::stride]:
        with Dataset(f) as ds:
            t = float(ds.time_s) / 3600.0
        if not (t0_h < t < t1_h):
            continue
        s = load_fields(f)
        mf = momentum_flux(s["u"], s["w"], s["rho"])
        acc = mf if acc is None else acc + mf
        n, z_km = n + 1, s["z_km"]
    if n == 0:
        raise ValueError("no snapshots in window %.1f-%.1f h for %s" % (t0_h, t1_h, run_dir))
    return z_km, acc / n, n


def cmd_zprofile(args):
    za, ma, na = window_profile(args.a, args.t0, args.t1, args.stride)
    zb, mb, nb = window_profile(args.b, args.t0, args.t1, args.stride)
    top = min(za.max(), zb.max())
    zz = za[za <= min(top, 60.0)]
    r = np.interp(zz, za, ma) / np.interp(zz, zb, mb)
    lo, hi = zz < args.z_split, (zz >= args.z_split) & (zz <= 50.5)
    rep = {"window_h": [args.t0, args.t1], "n_snap": [na, nb], "z_split_km": args.z_split,
           "ratio_median_above": float(np.median(r[hi])),
           "ratio_spread_above_pct": float(100 * (r[hi].max() - r[hi].min()) / abs(np.median(r[hi]))),
           "ratio_at_1.25": float(r[nearest_index(zz, 1.25)]),
           "ratio_at_10.25": float(r[nearest_index(zz, 10.25)]),
           "ratio_at_30.25": float(r[nearest_index(zz, 30.25)]),
           "ratio_at_50.25": float(r[nearest_index(zz, 50.25)])}
    print(json.dumps(rep, indent=1))
    print("Reading fixed in Section 14: H_L predicts a spread above %.1f km below 10 %% with the "
          "departure from 1 set below it; path dissipation predicts a trend concentrated at "
          "45-53 km." % args.z_split)
    if args.out:
        with open(args.out, "w") as f:
            json.dump({"report": rep, "z_km": zz.tolist(), "ratio": r.tolist()}, f, indent=1)


# ---------------------------------------------------------------------------
# zalt: the near-ground 2 dz alternation (Section 14, addendum 2 step 1)
# ---------------------------------------------------------------------------
def alternation(f):
    """Grid-scale (2 dz) component of a profile f(z), real or complex.

    a_j = (-1)^j (f_{j-1} - 2 f_j + f_{j+1}) / 4 demodulates the 2 dz part to a smooth
    envelope; the smooth part of f leaks into a_j as an alternating term of relative
    size (m dz)^2 / 4, which the [1, 2, 1] / 4 filter removes. Returns (alt, smooth) on
    levels 2 .. n-3 (index offset 2), with smooth = f minus the alternating part.
    """
    # A 2-D f is treated column by column along its first (z) axis.
    sh = (-1,) + (1,) * (np.ndim(f) - 1)
    j = np.arange(1, len(f) - 1).reshape(sh)
    a = (-1.0) ** j * (f[:-2] - 2 * f[1:-1] + f[2:]) / 4
    env = (a[:-2] + 2 * a[1:-1] + a[2:]) / 4
    jj = np.arange(2, len(f) - 2).reshape(sh)
    return env, f[2:-2] - (-1.0) ** jj * env


def decay_fit(relalt, z_km, lev_lo, lev_hi, dz_km):
    """Log-linear fit of the relative 2 dz amplitude over 1-based levels lev_lo..lev_hi.

    relalt and z_km are on levels 3 .. n-2 (the output of alternation). Returns the decay
    length [km] and the amplitude extrapolated to the first cell (z = dz/2), which assumes
    exponential decay below the fitted range, where the filter has no value.
    """
    sel = slice(lev_lo - 3, lev_hi - 2)
    zz, v = z_km[sel], relalt[sel]
    ok = v > 0
    slope, icpt = np.polyfit(zz[ok], np.log(v[ok]), 1)
    return {"decay_km": float(-1.0 / slope), "amp_first_cell": float(np.exp(icpt + slope * 0.5 * dz_km)),
            "fit_km": [float(zz[0]), float(zz[-1])]}


def cmd_zalt(args):
    """2 dz component of the window-mean w(z), its decay length and the descriptive model.

    --field mode (default) takes the complex mode-n coefficient of w; --field xrms takes
    the window-mean w(z, x) of a ridge run and forms rms_x of the 2 dz part over rms_x of
    the smooth part at each level. The decay length is fitted over --fit-levels (1-based)
    and set beside L = dz^2 / (2 c dt) with c = sqrt(gamma R T) at the ground of the MCD
    column (addendum 3 step 1b, descriptive, not registered).
    """
    runs = dict(kv.split("=", 1) for kv in args.run)
    raw = np.genfromtxt(args.atm_file, delimiter=",", names=True)
    res = {}
    for lab, d in runs.items():
        acc_w, acc_mf, n, z = 0, 0, 0, None
        for f in snapshot_files(d)[::args.stride]:
            with Dataset(f) as ds:
                t = float(ds.time_s) / 3600.0
            if not (args.t0 < t < args.t1):
                continue
            s = load_fields(f)
            nx = s["w"].shape[1]
            if args.field == "xrms":
                acc_w = acc_w + s["w"]
            else:
                acc_w = acc_w + np.fft.rfft(s["w"], axis=1)[:, args.mode] / nx
            acc_mf = acc_mf + momentum_flux(s["u"], s["w"], s["rho"])
            n, z, gam = n + 1, s["z_km"], s["gamma"]
        w6, mf = acc_w / n, acc_mf / n
        dt = mean_dt(d)
        dz_km = float(z[1] - z[0])
        c0 = float(np.sqrt(gam * R_GAS * raw["T"][0]))
        out = {"dt_s": dt, "n_snap": n, "field": args.field, "dz_km": dz_km, "c_ground_m_s": c0,
               "L_model_km": (dz_km * 1e3) ** 2 / (2 * c0 * dt) / 1e3}
        zz = z[2:-2]
        band = (zz >= args.z_lo) & (zz <= args.z_hi)
        if args.field == "xrms":
            alt, smooth = alternation(w6)
            profs = {"w_xrms": np.sqrt((alt ** 2).mean(axis=1)) / np.sqrt((smooth ** 2).mean(axis=1))}
            w_abs = np.sqrt((w6 ** 2).mean(axis=1))
        else:
            profs = {}
            for key, prof in (("w_complex", w6), ("w_abs", np.abs(w6))):
                alt, smooth = alternation(prof)
                profs[key] = np.abs(alt) / np.abs(smooth)
            w_abs = np.abs(w6)
        alt, smooth = alternation(mf)
        profs["MF"] = np.abs(alt) / np.abs(smooth)
        for key, relalt in profs.items():
            out[key] = {"rms_rel_alt": float(np.sqrt(np.mean(relalt[band] ** 2))),
                        "rel_alt_by_level": {"%.3f" % zk: float(v)
                                             for zk, v in zip(zz[:16], relalt[:16])}}
        wkey = "w_xrms" if args.field == "xrms" else "w_complex"
        out[wkey]["fit"] = decay_fit(profs[wkey], zz, args.fit_levels[0], args.fit_levels[1], dz_km)
        out["w_abs_mm_s"] = {"%.3f" % zk: float(1e3 * v) for zk, v in zip(z[:12], w_abs[:12])}
        res[lab] = out
        fit = out[wkey]["fit"]
        print("%-6s dt %.4f s  n %3d  rms rel. 2dz amplitude %.2f-%.2f km: w %.4f, MF %.4f | "
              "decay %.3f km (levels %d-%d), amp(first cell) %.3f, L_model %.3f km, fit/model %.3f"
              % (lab, dt, n, args.z_lo, args.z_hi, out[wkey]["rms_rel_alt"], out["MF"]["rms_rel_alt"],
                 fit["decay_km"], args.fit_levels[0], args.fit_levels[1], fit["amp_first_cell"],
                 out["L_model_km"], fit["decay_km"] / out["L_model_km"]))
    if args.out:
        with open(args.out, "w") as f:
            json.dump(res, f, indent=1)


# ---------------------------------------------------------------------------
# zprobe: ground-face w of the z-predictor against the forced w (addendum 2 step 3)
# ---------------------------------------------------------------------------
def cmd_zprobe(args):
    """Mode-n complex ratios of the face_probe.nc records (SimConfig.face_probe).

    w_ghost/wf: the face w the :ghost path's corrector uses; w_rebuilt/wf: the same
    predictor with the ghost rows rebuilt after the x-sweep (component 2 of the :face
    switch); w_avg/wf: the plain cell average, without the predictor's flux and source
    terms; p_ghost/p_rebuilt: face p' of the two predictors. Each is the ratio of the
    window means of the complex mode-n coefficients, over --t0-min to --t1-min, and over
    the sub-windows of --sub-min for stationarity.
    """
    res = {}
    for kv in args.run:
        lab, d = kv.split("=", 1)
        with Dataset(os.path.join(d, "face_probe.nc")) as ds:
            t = np.asarray(ds["time"][:], float) / 60.0
            V = {k: np.asarray(ds[k][:], float) for k in
                 ("wf", "w_ghost", "w_avg", "w_rebuilt", "p_ghost", "p_rebuilt")}
        c = {k: np.fft.rfft(v, axis=1)[:, args.mode] / v.shape[1] for k, v in V.items()}

        def ratios(sel):
            m = {k: v[sel].mean() for k, v in c.items()}
            return {"w_ghost/wf": m["w_ghost"] / m["wf"], "w_rebuilt/wf": m["w_rebuilt"] / m["wf"],
                    "w_avg/wf": m["w_avg"] / m["wf"], "p_ghost/p_rebuilt": m["p_ghost"] / m["p_rebuilt"]}

        main = ratios((t >= args.t0_min) & (t <= args.t1_min))
        subs = {}
        for a in np.arange(args.t0_min, args.t1_min, args.sub_min):
            subs["%g-%g min" % (a, a + args.sub_min)] = ratios((t >= a) & (t < a + args.sub_min))
        fmt = lambda r: {k: [abs(v), float(np.degrees(np.angle(v)))] for k, v in r.items()}
        res[lab] = {"run_dir": d, "n_records": int(((t >= args.t0_min) & (t <= args.t1_min)).sum()),
                    "window": fmt(main), "sub_windows": {k: fmt(v) for k, v in subs.items()}}
        print("%s (%d records, %g-%g min): " % (lab, res[lab]["n_records"], args.t0_min, args.t1_min) +
              "  ".join("%s %.4f (%.2f deg)" % (k, abs(v), np.degrees(np.angle(v))) for k, v in main.items()))
        for k, v in subs.items():
            print("    %-12s |w_ghost/wf| %.4f  |w_rebuilt/wf| %.4f  |w_avg/wf| %.4f" %
                  (k, abs(v["w_ghost/wf"]), abs(v["w_rebuilt/wf"]), abs(v["w_avg/wf"])))
    if args.out:
        with open(args.out, "w") as f:
            json.dump(res, f, indent=1)


# ---------------------------------------------------------------------------
# ztheory: the cavity theory level with the N(z) of the near-ground inversion
# ---------------------------------------------------------------------------
def cmd_ztheory(args):
    """Linear upward level of the single-mode C2R column under several N(z) treatments.

    V0 repeats compare_sponge_reflection.cavity_closure (grid levels, N from the T
    gradient of the interpolated MCD profile, Phi(0) by linear extrapolation from the
    first cell). V1 takes N^2 = g [(dP0/dz)/(gamma P0) - dln(rho0)/dz] from the solver's
    initial column (snapshot 0, x-mean), which is the buoyancy the discrete background
    carries. V2 integrates on a --dzf grid from z = 0 with the piecewise-linear MCD T
    and imposes Phi(0) at z = 0 exactly; V3 adds the compressible term to V2; V4 is V2
    with N^2 from the piecewise-linear rho and P that the solver interpolates. Every
    variant keeps the run's imposed Psi(0) and reads F_up as the 10-55 km mean on the
    solver levels. The measured levels (12-17 h means from the cavity JSON) are
    re-expressed against each variant's layer-theory flux.
    """
    from compare_sponge_reflection import PROFILE, run_setup  # noqa: E402
    from analyze_linear_flux_constancy import rayleigh_sponge_reflection, CP, G  # noqa: E402
    st = run_setup(args.run)
    n, U = st["mode"], st["u_bar"]
    s0 = load_fields(snapshot_files(args.run)[0], ("rho", "P"))
    zg = s0["z_km"] * 1e3
    lx = s0["x_m"].size * (s0["x_m"][1] - s0["x_m"][0])
    k = 2 * np.pi * n / lx
    raw = np.genfromtxt(PROFILE, delimiter=",", names=True)
    psi0 = np.sqrt(raw["rho"][0]) * U * k * st["h_m"] / 2
    band = (zg >= 10e3) & (zg <= 55e3)
    top, base = st["top_km"] * 1e3, st["base_km"] * 1e3

    def levels(z, N, rho, at_ground, c2=None):
        out = {}
        for K0, key in ((0.0, "open"), (st["K0"], "layer")):
            _, Phi, dPhi = rayleigh_sponge_reflection(k, z, N, rho, U, K0, base, top,
                                                      top="radiation", return_field=True,
                                                      c2=c2)
            c = psi0 / (Phi[0] if at_ground else Phi[0] - z[0] * dPhi[0])
            Phi, dPhi = c * Phi, c * dPhi
            m = np.sqrt(np.clip((N / U) ** 2 - k ** 2, 0, None))
            A = 0.5 * (Phi + dPhi / (1j * np.where(m > 0, m, 1.0)))
            fup = 2 * (m / k) * np.abs(A) ** 2
            idx = [int(np.argmin(np.abs(z - zz))) for zz in zg[band]]
            out[key] = {"F_up": float(fup[idx].mean()),
                        "psi_first_cell": float(np.abs(Phi[int(np.argmin(np.abs(z - zg[0])))]) / psi0)}
        return out

    T_g = np.interp(zg, raw["z_m"], raw["T"])
    rho_g = np.interp(zg, raw["z_m"], raw["rho"])
    N_T = lambda T, z: np.sqrt(np.clip((G / T) * (np.gradient(T, z) + G / CP), 1e-12, None))
    gam = s0["gamma"]
    rho0, P0 = s0["rho"].mean(axis=1), s0["P"].mean(axis=1)
    N_model = np.sqrt(np.clip(G_MARS * (np.gradient(P0, zg) / (gam * P0)
                                        - np.gradient(np.log(rho0), zg)), 1e-12, None))
    zf = np.arange(0.0, top + 0.5 * args.dzf, args.dzf)
    T_f = np.interp(zf, raw["z_m"], raw["T"])
    rho_f = np.interp(zf, raw["z_m"], raw["rho"])
    P_f = np.interp(zf, raw["z_m"], raw["P"])
    N_f = N_T(T_f, zf)
    N_fm = np.sqrt(np.clip(G_MARS * (np.gradient(P_f, zf) / (gam * P_f)
                                     - np.gradient(np.log(rho_f), zf)), 1e-12, None))
    c2_f = gam * P_f / rho_f
    var = {"V0 grid, N(T), extrapolated Phi(0)": levels(zg, N_T(T_g, zg), rho_g, False),
           "V1 grid, N of the solver column": levels(zg, N_model, rho0, False),
           "V2 fine, N(T), Phi(0) at z = 0": levels(zf, N_f, rho_f, True),
           "V3 fine, N(T), compressible": levels(zf, N_f, rho_f, True, c2=c2_f),
           "V4 fine, N of interpolated rho and P": levels(zf, N_fm, rho_f, True)}
    meas = {}
    for e in json.load(open(args.cavity)):
        hrs = [h for h in e["hours"] if args.t0 - 0.5 < h["t_h"] < args.t1 + 0.5]
        lab = os.path.basename(os.path.dirname(os.path.normpath(e["run_dir"])))
        meas[lab] = (float(np.mean([h["level_measured"] for h in hrs])),
                     float(np.mean([h["psi_ground_over_imposed"] for h in hrs])))
    v0 = var["V0 grid, N(T), extrapolated Phi(0)"]
    F_open0 = v0["open"]["F_up"]
    rep = {"N_near_ground_s-1": {"%.2f km" % (zz / 1e3): [float(a), float(b)]
                                 for zz, a, b in zip(zg[:8], N_T(T_g, zg)[:8], N_model[:8])},
           "variants": {}}
    print("N at the first levels, from T and from the solver column [s^-1]:")
    for zk, (a, b) in rep["N_near_ground_s-1"].items():
        print("  %s  %.5f  %.5f" % (zk, a, b))
    for name, v in var.items():
        lev = v["layer"]["F_up"] / v["open"]["F_up"]
        r = {"level_theory": lev, "F_layer_over_V0": v["layer"]["F_up"] / v0["layer"]["F_up"],
             "psi_first_cell_layer": v["layer"]["psi_first_cell"],
             "measured_over_layer_theory": {lab: m[0] * F_open0 / v["layer"]["F_up"]
                                            for lab, m in meas.items()}}
        rep["variants"][name] = r
        print("%-40s level %.4f  F_layer/F_layer(V0) %.4f  Psi(0.25)/imp %.4f" %
              (name, lev, r["F_layer_over_V0"], r["psi_first_cell_layer"]))
        for lab, x in r["measured_over_layer_theory"].items():
            print("    %-70s measured/theory %.4f" % (lab, x))
    if args.out:
        with open(args.out, "w") as f:
            json.dump(rep, f, indent=1)


# ---------------------------------------------------------------------------
# zflin: linear flux of the Gaussian ridge, factor by factor (addendum 3, V1)
# ---------------------------------------------------------------------------
def cmd_zflin(args):
    """F_lin of the uniform-wind Gaussian ridge and the measured F(z)/F_lin of each run.

    The ridge h = H exp(-x^2 / (2 sigma^2)) has the hydrostatic Boussinesq drag
    D = rho0 U N H^2 per unit length (Nappo 2002, eq. 3.121 as cited for addendum 3), so
    the domain-mean flux is F_N = -D / LX. The factors are successive ratios, so their
    product telescopes to the final theory:
      Q0  F_N with rho0, N0 at the first cell
      Q1  discrete periodic spectrum, 2 rho0 U N0 sum_n k_n |c_n|^2       f_disc = Q1/Q0
      Q2  non-hydrostatic, m = (N0^2/U^2 - k^2)^1/2                       f_nh   = Q2/Q1
      Q3  compressible, m^2 + k^2 U^2/c0^2 - 1/(4 H_rho^2)                f_comp = Q3/Q2
      Q4  open column, Taylor-Goldstein with N(z), rho(z), c(z) of the MCD column, radiation
          at the layer base, Phi(0) imposed at z = 0 (as ztheory V3), net flux over --band
                                                                          f_col  = Q4/Q3
      Q5  as Q4 with the sin^2 Rayleigh layer, steady (all returns in place)  f_layer = Q5/Q4
      Qm  Q5 for modes whose round trip to the layer base and back to the band top has
          elapsed by --t0 (after --ramp-delay-h), Q4 for the others.
    Primary (registered in addendum 3): per hour of --t0..--t1, the theory Qm(t) with the
    modes classed at the hour midpoint, the measured hour-mean |MF| over --band divided by
    it, and the mean of those hourly ratios. Classification bounds use delays of 0 and 3 h.
    """
    from analyze_linear_flux_constancy import rayleigh_sponge_reflection, CP, G  # noqa: E402
    raw = np.genfromtxt(args.atm_file, delimiter=",", names=True)
    LX, H, U, sig = args.lx_km * 1e3, args.h_m, args.u, args.sigma_km * 1e3
    nx = int(round(LX / args.dx))
    x = (np.arange(nx) + 0.5) * args.dx
    h = H * np.exp(-((x - 0.5 * LX) ** 2) / (2 * sig ** 2))
    c = np.fft.rfft(h) / nx
    kk = 2 * np.pi * np.fft.rfftfreq(nx, args.dx)
    n_use = np.arange(1, len(kk))
    n_use = n_use[kk[n_use] * U < 0.999 * 0.0176]   # propagating at the ground N
    top, base = args.top_km * 1e3, args.base_km * 1e3
    zg = (np.arange(int(round(top / args.dz))) + 0.5) * args.dz
    T_g = np.interp(zg, raw["z_m"], raw["T"])
    rho_g = np.interp(zg, raw["z_m"], raw["rho"])
    N_T = lambda T, z: np.sqrt(np.clip((G / T) * (np.gradient(T, z) + G / CP), 1e-12, None))
    N1, rho1 = float(N_T(T_g, zg)[0]), float(rho_g[0])
    gam = GAMMA
    c1 = float(np.sqrt(gam * R_GAS * T_g[0]))
    Hrho = float(-1.0 / np.gradient(np.log(rho_g), zg)[0])
    w = kk * np.abs(c) ** 2
    Q0 = rho1 * U * N1 * H ** 2 / LX
    Q1 = 2 * rho1 * U * N1 * w[n_use].sum()
    m_nh = np.sqrt(N1 ** 2 / U ** 2 - kk[n_use] ** 2)
    Q2 = 2 * rho1 * U ** 2 * (w[n_use] * m_nh).sum()
    m_c = np.sqrt(N1 ** 2 / U ** 2 - kk[n_use] ** 2 + kk[n_use] ** 2 * U ** 2 / c1 ** 2
                  - 1.0 / (4 * Hrho ** 2))
    Q3 = 2 * rho1 * U ** 2 * (w[n_use] * m_c).sum()
    zf = np.arange(0.0, top + 0.5 * args.dzf, args.dzf)
    T_f = np.interp(zf, raw["z_m"], raw["T"])
    rho_f = np.interp(zf, raw["z_m"], raw["rho"])
    P_f = np.interp(zf, raw["z_m"], raw["P"])
    N_f, c2_f = N_T(T_f, zf), gam * P_f / rho_f
    bsel = (zf >= args.band[0] * 1e3) & (zf <= args.band[1] * 1e3)
    I = np.concatenate([[0.0], np.cumsum(0.5 * (N_f[1:] + N_f[:-1]) * np.diff(zf))])
    I_at = lambda zz: float(np.interp(zz, zf, I))
    per = []
    for n in n_use:
        k = kk[n]
        psi0 = np.sqrt(raw["rho"][0]) * U * k * np.abs(c[n])
        F = {}
        # open: radiation 250 m below the layer base, so nothing above it enters (midway
        # between the 500 m nodes of the MCD table, where the s' of the piecewise-linear
        # rho has no spike); layer: the sin^2 layer up to the model top, radiation there
        for K0, zb, ztop, key in ((0.0, 0.0, base - 250.0, "open"), (args.k0, base, top, "layer")):
            _, Phi, dPhi = rayleigh_sponge_reflection(k, zf, N_f, rho_f, U, K0, zb, ztop,
                                                      top="radiation", return_field=True, c2=c2_f)
            s = psi0 / Phi[0]
            F[key] = float(np.mean((2.0 / k) * np.imag(np.conj(s * Phi[bsel]) * s * dPhi[bsel])))
        t_up = I_at(args.band[1] * 1e3) / (k * U ** 2) / 3600.0
        t_rt = (2 * I_at(base) - I_at(args.band[1] * 1e3)) / (k * U ** 2) / 3600.0
        per.append({"n": int(n), "lambda_km": float(2 * np.pi / k / 1e3), "share": 0.0,
                    "F_open": F["open"], "F_layer": F["layer"], "t_up_band_top_h": t_up,
                    "t_roundtrip_band_top_h": t_rt,
                    "returned": bool(t_rt + args.ramp_delay_h <= args.t0),
                    "arrived": bool(t_up + args.ramp_delay_h <= args.t0)})
    Q4 = sum(p["F_open"] for p in per)
    Q5 = sum(p["F_layer"] for p in per)
    Qm = sum(p["F_layer"] if p["returned"] else p["F_open"] for p in per)
    for p in per:
        p["share"] = p["F_open"] / Q4
    th = {"Q0_nappo_uPa": Q0 * 1e6, "Q1_discrete_uPa": Q1 * 1e6, "Q2_nonhydrostatic_uPa": Q2 * 1e6,
          "Q3_compressible_uPa": Q3 * 1e6, "Q4_open_column_uPa": Q4 * 1e6,
          "Q5_layer_steady_uPa": Q5 * 1e6, "Qm_primary_uPa": Qm * 1e6,
          "f_disc": Q1 / Q0, "f_nh": Q2 / Q1, "f_comp": Q3 / Q2, "f_col": Q4 / Q3,
          "f_layer": Q5 / Q4, "f_mixed": Qm / Q4,
          "inputs": {"rho1": rho1, "N1": N1, "c1": c1, "H_rho_m": Hrho, "rho_ground": float(raw["rho"][0]),
                     "n_modes": int(len(n_use)), "band_km": args.band, "window_h": [args.t0, args.t1]},
          "share_unarrived": float(sum(p["share"] for p in per if not p["arrived"])),
          "share_unreturned": float(sum(p["share"] for p in per if not p["returned"]))}
    print("F_lin factors (uPa, magnitudes): Q0 Nappo %.4f | f_disc %.4f | f_nh %.5f | f_comp %.5f | "
          "f_col %.4f | f_layer %.4f | mixed %.4f" % (th["Q0_nappo_uPa"], th["f_disc"], th["f_nh"],
                                                      th["f_comp"], th["f_col"], th["f_layer"], th["f_mixed"]))
    print("  Q1 %.4f Q2 %.4f Q3 %.4f Q4 %.4f Q5 %.4f Qm(primary) %.4f uPa; flux share unarrived %.4f, "
          "unreturned %.4f at %.1f h" % (Q1 * 1e6, Q2 * 1e6, Q3 * 1e6, Q4 * 1e6, Q5 * 1e6, Qm * 1e6,
                                          th["share_unarrived"], th["share_unreturned"], args.t0))
    hours = np.arange(np.floor(args.t0), np.ceil(args.t1)) + 0.5

    def q_hourly(delay):
        return np.array([sum(p["F_layer"] if p["t_roundtrip_band_top_h"] + delay <= t else p["F_open"]
                             for p in per) for t in hours])

    qh = {dl: q_hourly(dl) for dl in (0.0, args.ramp_delay_h, 3.0)}
    th["hourly_theory_uPa"] = {"%.1f" % t: float(q * 1e6) for t, q in zip(hours, qh[args.ramp_delay_h])}
    th["hourly_theory_over_Q4_by_delay"] = {"%g h" % dl: (q / Q4).tolist() for dl, q in qh.items()}
    print("  hourly theory / Q4 (delay %g h): %s" % (args.ramp_delay_h, " ".join(
        "%.4f" % v for v in qh[args.ramp_delay_h] / Q4)))
    meas = {}
    for kv in args.run or []:
        lab, d = kv.split("=", 1)
        acc, z_km = {}, None
        for f in snapshot_files(d)[::args.stride]:
            with Dataset(f) as ds:
                t = float(ds.time_s) / 3600.0
            if not (args.t0 < t < args.t1):
                continue
            s = load_fields(f)
            k_h = int(t)
            a = acc.setdefault(k_h, [0.0, 0])
            a[0] = a[0] + momentum_flux(s["u"], s["w"], s["rho"])
            a[1] += 1
            z_km = s["z_km"]
        b = (z_km >= args.band[0]) & (z_km <= args.band[1])
        Fh = np.array([abs((acc[int(t)][0] / acc[int(t)][1])[b].mean()) for t in hours])
        mf = sum(acc[int(t)][0] for t in hours) / sum(acc[int(t)][1] for t in hours)
        ratios = {dl: Fh / q for dl, q in qh.items()}
        prim = float(ratios[args.ramp_delay_h].mean())
        Fm = float(abs(mf[b].mean()))
        meas[lab] = {"dir": os.path.abspath(d), "dt_s": mean_dt(d),
                     "n_snap": int(sum(acc[int(t)][1] for t in hours)), "F_band_uPa": Fm * 1e6,
                     "hourly_F_band_uPa": (Fh * 1e6).tolist(),
                     "hourly_F_over_theory": ratios[args.ramp_delay_h].tolist(),
                     "primary_F_over_Flin": prim,
                     "classification_bounds": [float(ratios[0.0].mean()), float(ratios[3.0].mean())],
                     "F_over_Q5": Fm / Q5, "F_over_Q4": Fm / Q4, "F_over_Q3": Fm / Q3,
                     "F_over_Q1": Fm / Q1, "F_over_Q0": Fm / Q0,
                     "profile_over_Q4": {"%.2f" % zk: float(abs(mf[nearest_index(z_km, zk)]) / Q4)
                                         for zk in (1.25, 4.75, 10.25, 15.25, 20.25, 30.25, 40.25, 50.25)}}
        print("%-6s dt %.4f s  |F| band %.4f uPa  primary F/F_lin %.4f (classification %.4f-%.4f)  "
              "F/Q4 %.4f  F/Q1 %.4f  F/Q0 %.4f" % (lab, meas[lab]["dt_s"], Fm * 1e6, prim,
                                                  *meas[lab]["classification_bounds"], Fm / Q4, Fm / Q1, Fm / Q0))
        print("       hourly F/theory: " + " ".join("%.4f" % v for v in ratios[args.ramp_delay_h]))
        print("       F(z)/Q4: " + "  ".join("%s %.3f" % kv2 for kv2 in meas[lab]["profile_over_Q4"].items()))
    if args.out:
        with open(args.out, "w") as f:
            json.dump({"theory": th, "modes": per, "measured": meas}, f, indent=1)


# ---------------------------------------------------------------------------
# zequil: H = 0 equilibrium, drho_L2 from the solver log and energy from snapshots (V2)
# ---------------------------------------------------------------------------
LOG_RE = re.compile(r"^Step\s+(\d+) \| t=\s*([\d.]+) s .*?\| E=\s*([\d.e+-]+) \| KE=\s*([\d.e+-]+) "
                    r"\| drho_L2=([\d.e+-]+)")


def equil_budget(run_dir, t_end):
    """Mass budget of an H = 0 record from its snapshots to t_end (addendum 3, V2 note).

    All terms relative to the snapshot-0 mass sum(rho) (cells of equal area dx dz), time
    integrals by the trapezoid rule over the stored snapshots:
      sponge = int sum K(z) (rho_bar - rho) dt, with K = K0 sin^2(pi/2 eta) above the
               sponge base (src/sponge.jl, Srivastava et al. 2022 eq. 27) and rho_bar the
               snapshot 0 (valid only when the run starts from the background, H = 0);
      top    = -int sum_x rho w (last cell) dt / dz, the outflow through the top (inflow > 0);
      ghost_ground = int sum_x [w_p (rho_bar_p - rho_bar_g) / 2 - dt g rho'_p / 2] dt / dz,
               the leading terms of the ground-face predictor momentum (eq. 10a) when the
               ghost cells carry w_g = -w_p (impose_ghost_w!, src/forcing.jl) over a
               stratified background; rho_bar_g is extrapolated log-linearly. Zero by
               construction on the face path, where the corrector imposes w_f.
    dm_z, dE_z: change of sum_x rho and of sum_x (P/(gamma-1) + rho (u^2 + w^2)/2) per
    level [nz], relative to the snapshot-0 totals.
    """
    files = snapshot_files(run_dir)
    dt = mean_dt(run_dir)
    s0 = load_fields(files[0], ("rho", "u", "w", "P"))
    with Dataset(files[0]) as ds:
        g, k0, frac = float(ds.g), float(ds.K0_sponge), float(ds.sponge_frac)
    z = s0["z_km"] * 1e3
    dz = z[1] - z[0]
    z_max = z[-1] + 0.5 * dz
    z_sp = z_max * (1.0 - frac)
    K = np.where(z >= z_sp, k0 * np.sin(0.5 * np.pi * np.clip((z - z_sp) / (z_max - z_sp), 0, 1)) ** 2, 0.0)
    rb = s0["rho"]
    M0 = rb.sum()
    rb_g = rb[0] ** 2 / rb[1]
    energy = lambda s: s["P"] / (s["gamma"] - 1.0) + 0.5 * s["rho"] * (s["u"] ** 2 + s["w"] ** 2)
    E0z = energy(s0).sum(axis=1)
    t, sp, top, gw, gg = [], [], [], [], []
    for f in files:
        s = load_fields(f, ("rho", "u", "w", "P"))
        if s["t_s"] > t_end + 1.0:
            break
        r = s["rho"]
        t.append(s["t_s"])
        sp.append((K[:, None] * (rb - r)).sum() / M0)
        top.append(-(r[-1] * s["w"][-1]).sum() / dz / M0)
        gw.append((0.5 * s["w"][0] * (rb[0] - rb_g)).sum() / dz / M0)
        gg.append((-0.5 * dt * g * (r[0] - rb[0])).sum() / dz / M0)
        last = s
    integ = lambda y: float(np.trapezoid(np.array(y), np.array(t)))
    return {"t_s": t, "z_km": s0["z_km"].tolist(), "dt_s": dt,
            "sponge": integ(sp), "top": integ(top), "ghost_ground_w": integ(gw),
            "ghost_ground_g": integ(gg), "ghost_ground": integ(gw) + integ(gg),
            "dm_z": ((last["rho"] - rb).sum(axis=1) / M0).tolist(),
            "dE_z": ((energy(last).sum(axis=1) - E0z) / E0z.sum()).tolist()}


def cmd_zequil(args):
    """drho_L2 at every logged step to --t-end and the energy and mass drift to --t-end.

    LABEL=run_dir=log. drho_L2 is the solver's own diagnostic (every 500 steps); total
    energy sum(P/(gamma-1) + rho (u^2 + w^2)/2) and mass sum(rho) over the stored
    physical domain are compared at full precision between snapshot 0 and the snapshot
    nearest --t-end.
    """
    res = {}
    for kv in args.run:
        lab, d, log = kv.split("=", 2)
        rows = [LOG_RE.match(l) for l in open(log, errors="replace")]
        rows = [(int(m.group(1)), float(m.group(2)), float(m.group(3)), float(m.group(5)))
                for m in rows if m and float(m.group(2)) <= args.t_end + 1.0]
        files = snapshot_files(d)
        ts = []
        for f in files:
            with Dataset(f) as ds:
                ts.append(float(ds.time_s))
        f1 = files[int(np.argmin(np.abs(np.array(ts) - args.t_end)))]

        def totals(f):
            s = load_fields(f, ("rho", "u", "w", "P"))
            E = s["P"] / (s["gamma"] - 1.0) + 0.5 * s["rho"] * (s["u"] ** 2 + s["w"] ** 2)
            return float(E.sum()), float(s["rho"].sum()), s["t_s"]

        E0, M0, _ = totals(files[0])
        E1, M1, t1 = totals(f1)
        dr = np.array([r[3] for r in rows])
        tt = np.array([r[1] for r in rows])
        late = tt >= 0.5 * args.t_end
        res[lab] = {"dir": os.path.abspath(d), "log": os.path.abspath(log),
                    "drho_L2_by_step": {str(r[0]): r[3] for r in rows},
                    "drho_L2_peak": float(dr.max()), "t_peak_s": float(tt[dr.argmax()]),
                    "drho_L2_mean_late": float(dr[late].mean()), "drho_L2_last": float(dr[-1]),
                    "E_log_first_last": [rows[0][2], rows[-1][2]],
                    "t_snapshot_s": t1, "energy_drift": E1 / E0 - 1.0, "mass_drift": M1 / M0 - 1.0}
        r = res[lab]
        print("%-6s drho_L2 peak %.3e (%.0f s) mean %.0f-%.0f s %.3e last %.3e | E log %.6e -> %.6e | "
              "snapshot t %.0f s: energy drift %.3e, mass drift %.3e"
              % (lab, r["drho_L2_peak"], r["t_peak_s"], 0.5 * args.t_end, args.t_end, r["drho_L2_mean_late"],
                 r["drho_L2_last"], *r["E_log_first_last"], t1, r["energy_drift"], r["mass_drift"]))
        if args.budget:
            r["budget"] = equil_budget(d, args.t_end)
            b = r["budget"]
            print("       budget to %.0f s: sponge %.3e, top outflow %.3e, sum %.3e | ghost ground-face "
                  "estimate %.3e (w term %.3e, dt g rho' term %.3e)"
                  % (b["t_s"][-1], b["sponge"], b["top"], b["sponge"] + b["top"], b["ghost_ground"],
                     b["ghost_ground_w"], b["ghost_ground_g"]))
    if args.budget and len(res) == 2:
        (la, ra), (lb, rb) = res.items()
        ba, bb = ra["budget"], rb["budget"]
        z = np.array(ba["z_km"])
        pair = {}
        for q in ("dm_z", "dE_z"):
            dd = np.array(bb[q]) - np.array(ba[q])
            cum = np.cumsum(dd)
            pair[q] = {"diff_z": dd.tolist(), "cum_from_ground": cum.tolist(),
                       "z_half_km": float(z[np.searchsorted(cum / cum[-1], 0.5)]),
                       "z_99_km": float(z[np.searchsorted(cum / cum[-1], 0.99)])}
            print("%s - %s, %s: total %.3e; half reached by %.2f km, 99 %% by %.2f km"
                  % (lb, la, q, cum[-1], pair[q]["z_half_km"], pair[q]["z_99_km"]))
        res["pair"] = {"labels": [la, lb], **pair}
    if args.out:
        with open(args.out, "w") as f:
            json.dump(res, f, indent=1)


# ---------------------------------------------------------------------------
# zqueney: Queney case, face against ghost path (addendum 3, V3)
# ---------------------------------------------------------------------------
def queney_upert(x, z, b=10.0, a=20e3, ubar=18.0, gamma=1.4, g=9.81, R=287.0, T0=239.0):
    """Queney (1948) u' of the isothermal hydrostatic bell, as analyze_mountain_benchmark_strict.py."""
    cs = np.sqrt(gamma * R * T0)
    ks = np.sqrt(gamma - 1.0) * g / cs / ubar
    H = R * T0 / g
    X, Z = np.meshgrid(x - 0.5 * (x.min() + x.max()), z)
    zeta = b * np.exp(Z / (2 * H)) * (np.cos(ks * Z) - (X / a) * np.sin(ks * Z)) / (1 + (X / a) ** 2)
    return -ubar * np.gradient(zeta, z[1] - z[0], axis=0)


def cmd_zqueney(args):
    """Face/ghost amplitude of w over 0..--z-top km, |x - x0| <= --half-width km.

    Primary: ratio of the rms of w (face over ghost), each averaged over the snapshots
    in --t0..--t1 h; with it the regression ratio sum(w_f w_g)/sum(w_g^2) and the hourly
    ratios. Also the Srivastava Table 11 metric, the rms of u' - u'_Queney over 0-40 km as
    a percentage of U, at the last snapshot of each run. --metric-tops adds the metric over
    0 to each listed height [km], split into its x-mean part and its x-deviation (wave)
    part, whose squares add to the square of the total (Section 15, S14V3p).
    """
    out = {}
    fields = {}
    for lab, d in (("ghost", args.ghost), ("face", args.face)):
        files = snapshot_files(d)
        sel, tt = [], []
        for f in files:
            with Dataset(f) as ds:
                t = float(ds.time_s) / 3600.0
            if args.t0 <= t <= args.t1 + 1e-6:
                sel.append(f)
                tt.append(t)
        ws = [load_fields(f, ("w", "u", "rho")) for f in sel]
        s0 = ws[0]
        x, z = s0["x_m"], s0["z_km"] * 1e3
        box = np.ix_(z <= args.z_top * 1e3, np.abs(x - 0.5 * (x.min() + x.max())) <= args.half_width * 1e3)
        last = load_fields(files[-1], ("u",))
        up = last["u"] - args.ubar
        th = queney_upert(x, z)
        phys = z <= 40e3
        fields[lab] = {"t_h": tt, "w": [s["w"][box] for s in ws], "dt": mean_dt(d)}
        out[lab] = {"dir": os.path.abspath(d), "dt_s": mean_dt(d), "n_snap": len(sel),
                    "t_last_h": float(last["t_s"] / 3600.0),
                    "rms_w_mm_s": float(1e3 * np.mean([np.sqrt(np.mean(w ** 2)) for w in fields[lab]["w"]])),
                    "u_rms_err_pct_of_U": float(100 * np.sqrt(np.mean((up[phys] - th[phys]) ** 2)) / args.ubar)}
        for ztop in args.metric_tops or []:
            e = (up - th)[z <= ztop * 1e3]
            em = e.mean(axis=1, keepdims=True)
            out[lab]["metric_0_%g_km" % ztop] = {
                "total_pct": float(100 * np.sqrt(np.mean(e ** 2)) / args.ubar),
                "xmean_pct": float(100 * np.sqrt(np.mean(em ** 2)) / args.ubar),
                "wave_pct": float(100 * np.sqrt(np.mean((e - em) ** 2)) / args.ubar),
                "theory_rms_pct": float(100 * np.sqrt(np.mean(th[z <= ztop * 1e3] ** 2)) / args.ubar),
                "max_abs_xmean_u_minus_U_m_s": float(np.abs(up[z <= ztop * 1e3].mean(axis=1)).max())}
            print("%s 0-%g km: metric %.3f %% (x-mean %.3f, wave %.3f; theory rms %.3f %%), max |<u> - U| %.4f m/s"
                  % (lab, ztop, *[out[lab]["metric_0_%g_km" % ztop][k] for k in
                                  ("total_pct", "xmean_pct", "wave_pct", "theory_rms_pct",
                                   "max_abs_xmean_u_minus_U_m_s")]))
    g, f = fields["ghost"], fields["face"]
    n = min(len(g["w"]), len(f["w"]))
    amp = [np.sqrt(np.mean(f["w"][i] ** 2)) / np.sqrt(np.mean(g["w"][i] ** 2)) for i in range(n)]
    reg = [float(np.sum(f["w"][i] * g["w"][i]) / np.sum(g["w"][i] ** 2)) for i in range(n)]
    out["face_over_ghost"] = {"amplitude_primary": out["face"]["rms_w_mm_s"] / out["ghost"]["rms_w_mm_s"],
                              "regression_mean": float(np.mean(reg)),
                              "by_snapshot": {"%.3f h" % g["t_h"][i]: [float(amp[i]), reg[i]] for i in range(n)},
                              "dt_ratio": f["dt"] / g["dt"]}
    o = out["face_over_ghost"]
    print("ghost dt %.5f s, face dt %.5f s (ratio %.5f); rms w %.4f / %.4f mm/s; face/ghost amplitude %.4f, "
          "regression %.4f" % (g["dt"], f["dt"], o["dt_ratio"], out["ghost"]["rms_w_mm_s"],
                               out["face"]["rms_w_mm_s"], o["amplitude_primary"], o["regression_mean"]))
    print("u'_RMS error (Table 11 metric, last snapshot): ghost %.3f %% at %.2f h, face %.3f %% at %.2f h"
          % (out["ghost"]["u_rms_err_pct_of_U"], out["ghost"]["t_last_h"],
             out["face"]["u_rms_err_pct_of_U"], out["face"]["t_last_h"]))
    for k, v in o["by_snapshot"].items():
        print("   %s amplitude %.4f regression %.4f" % (k, *v))
    if args.out:
        with open(args.out, "w") as fo:
            json.dump(out, fo, indent=1)


# ---------------------------------------------------------------------------
# Z2: growth law of the near-ground density departure
# ---------------------------------------------------------------------------
def cmd_zdensity(args):
    for d in args.run:
        rows = run_timeseries(d, stride=args.stride)["hours"]
        t = np.array([r["t_h"] for r in rows])
        y = np.array([r["drho_max_below45_pct"] for r in rows])
        out = {}
        for lo, hi in ((2, 9), (9, 24)):
            m = (t > lo) & (t < hi) & (y > 0)
            if m.sum() >= 3:
                p = np.polyfit(np.log(t[m]), np.log(y[m]), 1)
                out["%d-%d h" % (lo, hi)] = float(p[0])
        print("%s\n  exponent of max x-mean density departure below 45 km in t: %s"
              % (d, json.dumps(out)))
    print("Reading fixed in Section 14: secular pumping of a calm layer gives ~2 early; a "
          "steady wave-driven mean-flow response gives ~1.")


# ---------------------------------------------------------------------------
# Z3: precursors of the h = 3 km collapses
# ---------------------------------------------------------------------------
def cmd_zh3(args):
    files = snapshot_files(args.run)[::args.stride]
    s0 = load_fields(files[0], ("T", "P", "w", "u", "rho"))
    z = s0["z_km"]
    T0 = s0["T"].mean(axis=1, keepdims=True)
    kappa = (s0["gamma"] - 1.0) / s0["gamma"]
    ground, top, free = z < args.z_ground, z > args.z_top, (z > 20) & (z < 60)
    nx = s0["w"].shape[1]
    rows = []
    for f in files:
        s = load_fields(f, ("T", "P", "w"))
        dT = np.abs(s["T"] / T0 - 1.0)
        wt = s["w"][top]
        var = (wt ** 2).mean()
        spec = np.abs(np.fft.rfft(wt, axis=1)) ** 2
        spec[:, 1:] *= 2.0
        tot = spec.sum()
        theta = s["T"] * (s["P"][0:1, :] / s["P"]) ** kappa
        unst = (np.diff(theta, axis=0)[free[:-1]] < 0).mean()
        rows.append({"t_h": s["t_s"] / 3600.0,
                     "dT_ground_max": float(dT[ground].max()),
                     "w_top_max": float(np.abs(wt).max()),
                     "k0_frac_top": float(spec[:, 0].sum() / tot) if tot > 0 else float("nan"),
                     "ridge_frac_top": float(spec[:, 1:args.n_ridge + 1].sum() / tot)
                     if tot > 0 else float("nan"),
                     "unstable_frac_20_60": float(unst), "w_top_var": float(var)})

    def first(key, thr):
        for r in rows:
            if r[key] > thr:
                return r["t_h"]
        return float("nan")

    t_g, t_t, t_u = first("dT_ground_max", 0.5), first("w_top_max", 5.0), first("unstable_frac_20_60", 0.01)
    rep = {"first_dT_ground_gt_0.5_h": t_g, "first_w_top_gt_5_h": t_t,
           "first_unstable_20_60_gt_1pct_h": t_u, "nx": nx, "n_ridge": args.n_ridge,
           "last": rows[-1]}
    print(json.dumps(rep, indent=1))
    print("Readings fixed in Section 14: H_L if the ground crossing precedes the top crossing "
          "and the top variance is mostly k = 0 or outside the ridge band; top reflection if the "
          "top crosses first with ridge_frac_top > 0.5; free-column breaking if the 20-60 km "
          "unstable fraction crosses first.")
    if args.out:
        with open(args.out, "w") as f:
            json.dump({"report": rep, "rows": rows}, f, indent=1)


# ---------------------------------------------------------------------------
# zoverturn: first overturning and the Lindzen saturation height (Section 15)
# ---------------------------------------------------------------------------
def lindzen_height(z, ubar, n_bv, rho_bar, amp_ref, z_ref, z_top):
    """Lowest z in (z_ref, z_top) where the WKB-grown stationary-wave |u'| reaches |ubar|.

    Hydrostatic stationary wave, constant flux rho |u'|^2 k |ubar| / (2 N) (Lindzen 1981,
    c = 0) gives |u'|(z) = amp_ref (N rho_ref |ubar_ref| / (N_ref rho |ubar|))^(1/2).
    z [m], ubar [m/s], n_bv [s^-1], rho_bar [kg m^-3] on the same [nz] grid; amp_ref [m/s]
    at z_ref [m]. Valid where ubar keeps one sign and N^2 > 0 between z_ref and z_top.
    Returns (z_s [m] or nan, |u'|(z) [nz]).
    """
    j0 = nearest_index(z, z_ref)
    amp = amp_ref * np.sqrt(n_bv * rho_bar[j0] * abs(ubar[j0]) /
                            (n_bv[j0] * rho_bar * np.maximum(np.abs(ubar), 1e-6)))
    sel = (z > z[j0]) & (z < z_top) & (amp >= np.abs(ubar))
    return (float(z[sel][0]) if sel.any() else float("nan")), amp


def cmd_zoverturn(args):
    """First convective and kinematic overturning by height band, and the Lindzen height.

    Convective: d theta/dz < 0 between adjacent levels (theta = T (1/P)^kappa, kappa =
    (gamma - 1)/gamma), the face height assigned. Kinematic: u sign(ubar) <= 0 where
    |ubar(z)| >= --u-floor and the current x-mean u keeps the sign of ubar (a local
    reversal, not the descent of a critical level), ubar the x-mean u of snapshot 0 (the
    background; the forcing ramp starts at zero). Both below --z-max [km]. Lindzen: lindzen_height with amp_ref the
    x-max of u minus its current x-mean at --z-ref [km], z_top the first zero of ubar above
    z_ref. Fields [z, x] via load_fields.
    """
    edges = [float(v) for v in args.bands.split(",")]
    out = {}
    for kv in args.run:
        lab, d = kv.split("=", 1)
        files = snapshot_files(d)[::args.stride]
        s0 = load_fields(files[0], ("u", "P", "T", "rho"))
        z = s0["z_km"] * 1e3
        zf = 0.5 * (z[1:] + z[:-1])
        kappa = (s0["gamma"] - 1.0) / s0["gamma"]
        with Dataset(files[0]) as ds:
            g = float(ds.g)
        ubar = s0["u"].mean(axis=1)
        th0 = (s0["T"] * s0["P"] ** (-kappa)).mean(axis=1)
        n2 = g * np.gradient(np.log(th0), z)
        n_bv = np.sqrt(np.maximum(n2, 1e-12))
        rho_bar = s0["rho"].mean(axis=1)
        j_ref = nearest_index(z, args.z_ref * 1e3)
        above = (z > z[j_ref]) & (np.sign(ubar) != np.sign(ubar[j_ref]))
        z_top = float(z[above][0]) if above.any() else float(z[-1])
        kin_ok = (np.abs(ubar) >= args.u_floor) & (z <= args.z_max * 1e3)
        rows = []
        for f in files:
            s = load_fields(f, ("u", "P", "T"))
            th = s["T"] * s["P"] ** (-kappa)
            unst = np.diff(th, axis=0) < 0                     # [nz-1, nx], faces
            keep = np.sign(s["u"].mean(axis=1)) == np.sign(ubar)    # excludes critical-level descent
            rev = (s["u"] * np.sign(ubar)[:, None] <= 0) & (kin_ok & keep)[:, None]
            row = {"t_h": s["t_s"] / 3600.0}
            for lo, hi in zip(edges[:-1], edges[1:]):
                bf = (zf >= lo * 1e3) & (zf < hi * 1e3)
                bc = (z >= lo * 1e3) & (z < hi * 1e3)
                key = "%g-%g" % (lo, hi)
                uf = unst[bf]
                hit = uf.any(axis=1)
                row["unst_frac_" + key] = float(uf.mean())
                row["unst_zlo_km_" + key] = float(zf[bf][hit][0] / 1e3) if hit.any() else None
                rv = rev[bc].any(axis=1)
                row["rev_n_" + key] = int(rev[bc].sum())
                row["rev_zlo_km_" + key] = float(z[bc][rv][0] / 1e3) if rv.any() else None
            up = s["u"][j_ref] - s["u"][j_ref].mean()
            amp = float(np.abs(up).max())
            zs, _ = lindzen_height(z, ubar, n_bv, rho_bar, amp, z[j_ref], z_top)
            row["amp_ref_m_s"] = amp
            row["z_lindzen_km"] = zs / 1e3
            rows.append(row)
        rep = {"dir": os.path.abspath(d), "stride": args.stride, "u_floor": args.u_floor,
               "z_max_km": args.z_max, "z_ref_km": float(z[j_ref] / 1e3), "z_top_km": z_top / 1e3,
               "unstable_faces_at_t0": int((np.diff(s0["T"] * s0["P"] ** (-kappa), axis=0) < 0).sum())}
        for lo, hi in zip(edges[:-1], edges[1:]):
            key = "%g-%g" % (lo, hi)
            for kind in ("unst", "rev"):
                first = next((r for r in rows if r[kind + "_zlo_km_" + key] is not None), None)
                rep["first_%s_%s" % (kind, key)] = (
                    {"t_h": first["t_h"], "z_km": first[kind + "_zlo_km_" + key],
                     "z_lindzen_km": first["z_lindzen_km"]} if first else None)
        out[lab] = {"report": rep, "rows": rows}
        print("%s: %s" % (lab, json.dumps(rep, indent=1)))
    if args.out:
        with open(args.out, "w") as f:
            json.dump(out, f, indent=1)


# ---------------------------------------------------------------------------
# Z5: gate on the launch-consistent (ublend) background before S14E/S14F launch
# ---------------------------------------------------------------------------
def cmd_zblend(args):
    """Richardson number of the blended wind that `ublend=Z` builds in the driver.

    u(z) = u_L (1 - S) + u_MCD(z) S,  S = (1 - cos(pi z / Z)) / 2 for z < Z, with
    u_L = sign(u_MCD(Z)) |U_bar| (make_blended_wind in scripts/sweep_arsia_mons_stage2.jl).
    Ri = N^2 / (du/dz)^2 with N^2 = (g/T)(dT/dz + g/c_p) from the MCD column. Units SI.
    Gate fixed in Section 14: the pair is withdrawn (Z is not retuned) if min Ri below Z < 1
    or if the blended wind changes sign below Z.
    """
    d = np.genfromtxt(args.atm_file, delimiter=",", names=True)
    z, T, u = d["z_m"], d["T"], d["u_zonal"]
    zf = np.linspace(0.0, args.zb_km * 1e3, 801)
    u_mcd = np.interp(zf, z, u)
    uL = np.sign(np.interp(args.zb_km * 1e3, z, u)) * abs(args.u_bar)
    S = 0.5 * (1 - np.cos(np.pi * zf / (args.zb_km * 1e3)))
    ub = uL * (1 - S) + u_mcd * S
    cp = GAMMA * R_GAS / (GAMMA - 1.0)
    n2 = np.interp(zf, z, G_MARS / T * (np.gradient(T, z) + G_MARS / cp))
    shear = np.gradient(ub, zf)
    ri = n2 / np.maximum(shear ** 2, 1e-14)
    rep = {"u_L_m_s": float(uL), "u_MCD_at_Z_m_s": float(u_mcd[-1]),
           "sign_change_below_Z": bool(np.any(np.sign(ub) != np.sign(uL))),
           "min_Ri_below_Z": float(ri.min()), "z_of_min_Ri_km": float(zf[np.argmin(ri)] / 1e3),
           "max_shear_s-1": float(np.abs(shear).max())}
    rep["gate_pass"] = (not rep["sign_change_below_Z"]) and rep["min_Ri_below_Z"] >= 1.0
    print(json.dumps(rep, indent=1))
    print("Pass u=%.0f to the driver so that the forcing wind has the sign of u_L." % uL)


# ---------------------------------------------------------------------------
# bitcheck: bit-for-bit reproduction of an earlier record by the patched driver
# ---------------------------------------------------------------------------
def cmd_bitcheck(args):
    """Compare every snapshot the new run shares (by step) with the reference record.

    Exit status 1 unless at least two snapshots are shared and every field (rho, u, w,
    P, T) and the time attribute agree exactly.
    """
    new = {STEP_RE.search(f).group(1): f for f in snapshot_files(args.new)}
    ref = {STEP_RE.search(f).group(1): f for f in snapshot_files(args.ref)}
    common = sorted(set(new) & set(ref))
    worst = 0.0
    for s in common:
        with Dataset(new[s]) as a, Dataset(ref[s]) as b:
            d = abs(float(a.time_s) - float(b.time_s))
            for v in ("rho", "u", "w", "P", "T"):
                d = max(d, float(np.max(np.abs(np.asarray(a.variables[v][:], float)
                                              - np.asarray(b.variables[v][:], float)))))
        worst = max(worst, d)
        print("step %s  max|diff| = %.3e" % (s, d))
    ok = len(common) >= 2 and worst == 0.0
    print("%s: %d shared snapshots, max|diff| over all = %r -> %s"
          % (os.path.basename(os.path.normpath(args.new)), len(common), worst,
             "IDENTICAL" if ok else "DIFFERENT"))
    sys.exit(0 if ok else 1)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sp = p.add_subparsers(dest="cmd", required=True)
    a = sp.add_parser("suite")
    a.add_argument("--run", action="append", required=True, help="LABEL=run_dir")
    a.add_argument("--outdir", default="output/reflection_tests/s14")
    a.add_argument("--stride", type=int, default=1)
    a.add_argument("--window", type=float, nargs=2, default=None,
                   help="window [h] for every run, in place of the per-label default")
    a.add_argument("--pair", action="append", default=None, help="A:B, change of B against A")
    a = sp.add_parser("zustar")
    a.add_argument("--atm-file", default="external/mcd_profiles/mars_arsia_mons_my31.csv")
    a.add_argument("--dt", type=float, required=True)
    a.add_argument("--dz", type=float, default=500.0)
    a.add_argument("--lx-km", type=float, default=4000.0)
    a.add_argument("--modes", type=int, nargs="+", default=[1, 2, 4, 6, 10])
    a.add_argument("--z-crit-low", type=float, default=4.5)
    a.add_argument("--out", default=None)
    a = sp.add_parser("zprofile")
    a.add_argument("--a", required=True)
    a.add_argument("--b", required=True)
    a.add_argument("--t0", type=float, default=18.0)
    a.add_argument("--t1", type=float, default=24.0)
    a.add_argument("--z-split", type=float, default=10.25)
    a.add_argument("--stride", type=int, default=1)
    a.add_argument("--out", default=None)
    a = sp.add_parser("zalt")
    a.add_argument("--run", action="append", required=True, help="LABEL=run_dir")
    a.add_argument("--mode", type=int, default=6)
    a.add_argument("--t0", type=float, default=12.0)
    a.add_argument("--t1", type=float, default=17.0)
    a.add_argument("--z-lo", type=float, default=0.75)
    a.add_argument("--z-hi", type=float, default=4.75)
    a.add_argument("--stride", type=int, default=1)
    a.add_argument("--field", choices=("mode", "xrms"), default="mode")
    a.add_argument("--fit-levels", type=int, nargs=2, default=[3, 10],
                   help="1-based levels of the log-linear decay fit")
    a.add_argument("--atm-file", default="external/mcd_profiles/mars_arsia_mons_my31.csv")
    a.add_argument("--out", default=None)
    a = sp.add_parser("zprobe")
    a.add_argument("--run", action="append", required=True, help="LABEL=run_dir")
    a.add_argument("--mode", type=int, default=6)
    a.add_argument("--t0-min", type=float, default=40.0)
    a.add_argument("--t1-min", type=float, default=60.0)
    a.add_argument("--sub-min", type=float, default=10.0)
    a.add_argument("--out", default=None)
    a = sp.add_parser("ztheory")
    a.add_argument("--run", required=True, help="single-mode run (setup and snapshot 0)")
    a.add_argument("--cavity", required=True, help="cavity_closure JSON with measured levels")
    a.add_argument("--dzf", type=float, default=25.0, help="fine-grid spacing [m]")
    a.add_argument("--t0", type=float, default=12.0)
    a.add_argument("--t1", type=float, default=17.0)
    a.add_argument("--out", default=None)
    a = sp.add_parser("zflin")
    a.add_argument("--run", action="append", default=None, help="LABEL=run_dir (optional)")
    a.add_argument("--atm-file", default="external/mcd_profiles/mars_arsia_mons_my31.csv")
    a.add_argument("--lx-km", type=float, default=2000.0)
    a.add_argument("--dx", type=float, default=4000.0)
    a.add_argument("--dz", type=float, default=500.0)
    a.add_argument("--sigma-km", type=float, default=50.0)
    a.add_argument("--h-m", type=float, default=30.0)
    a.add_argument("--u", type=float, default=40.0)
    a.add_argument("--k0", type=float, default=0.005)
    a.add_argument("--base-km", type=float, default=60.0)
    a.add_argument("--top-km", type=float, default=95.0)
    a.add_argument("--dzf", type=float, default=25.0)
    a.add_argument("--band", type=float, nargs=2, default=[5.25, 15.25])
    a.add_argument("--t0", type=float, default=12.0)
    a.add_argument("--t1", type=float, default=17.0)
    a.add_argument("--ramp-delay-h", type=float, default=2.0)
    a.add_argument("--stride", type=int, default=1)
    a.add_argument("--out", default=None)
    a = sp.add_parser("zequil")
    a.add_argument("--run", action="append", required=True, help="LABEL=run_dir=log")
    a.add_argument("--t-end", type=float, default=3600.0)
    a.add_argument("--budget", action="store_true",
                   help="mass budget (sponge, top, ghost ground face) and per-level drift; "
                        "with two runs, their difference cumulated from the ground")
    a.add_argument("--out", default=None)
    a = sp.add_parser("zqueney")
    a.add_argument("--ghost", required=True)
    a.add_argument("--face", required=True)
    a.add_argument("--t0", type=float, default=8.0)
    a.add_argument("--t1", type=float, default=10.0)
    a.add_argument("--z-top", type=float, default=10.0, help="km")
    a.add_argument("--half-width", type=float, default=100.0, help="km about the hill")
    a.add_argument("--ubar", type=float, default=18.0)
    a.add_argument("--metric-tops", type=float, nargs="+", default=None, help="km")
    a.add_argument("--out", default=None)
    a = sp.add_parser("zdensity")
    a.add_argument("--run", action="append", required=True)
    a.add_argument("--stride", type=int, default=1)
    a = sp.add_parser("zh3")
    a.add_argument("--run", required=True)
    a.add_argument("--stride", type=int, default=1)
    a.add_argument("--z-ground", type=float, default=5.0)
    a.add_argument("--z-top", type=float, default=85.0)
    a.add_argument("--n-ridge", type=int, default=12)
    a.add_argument("--out", default=None)
    a = sp.add_parser("zoverturn")
    a.add_argument("--run", action="append", required=True, help="LABEL=run_dir")
    a.add_argument("--stride", type=int, default=1)
    a.add_argument("--bands", default="0,5,20,45,55", help="band edges [km]")
    a.add_argument("--z-max", type=float, default=55.0, help="km")
    a.add_argument("--u-floor", type=float, default=5.0, help="|ubar| floor of the reversal test [m/s]")
    a.add_argument("--z-ref", type=float, default=10.25, help="Lindzen reference level [km]")
    a.add_argument("--out", default=None)
    a = sp.add_parser("zblend")
    a.add_argument("--atm-file", default="external/mcd_profiles/mars_arsia_mons_my31.csv")
    a.add_argument("--zb-km", type=float, default=8.0)
    a.add_argument("--u-bar", type=float, default=40.0)
    a = sp.add_parser("bitcheck")
    a.add_argument("--new", required=True, help="run directory written by the patched driver")
    a.add_argument("--ref", required=True, help="archived run directory it must reproduce")
    args = p.parse_args()
    {"suite": cmd_suite, "zustar": cmd_zustar, "zprofile": cmd_zprofile, "zalt": cmd_zalt,
     "ztheory": cmd_ztheory, "zprobe": cmd_zprobe, "zflin": cmd_zflin, "zequil": cmd_zequil, "zqueney": cmd_zqueney,
     "zdensity": cmd_zdensity, "zh3": cmd_zh3, "zoverturn": cmd_zoverturn, "zblend": cmd_zblend,
     "bitcheck": cmd_bitcheck}[args.cmd](args)


if __name__ == "__main__":
    main()
