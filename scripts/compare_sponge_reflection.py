"""Compare the damping-layer variants of the Arsia linear control corner.

Reads the `linear_flux_constancy.json` that analyze_linear_flux_constancy.py writes
for each run and sets the runs side by side on the three questions of
docs/reflection_convergence_tests.md:

1. Does the reflected share in the modes whose round trip from the layer has
   elapsed follow the linear Rayleigh-layer prediction across the variants? The
   split (updown_split) and the variance ripple (ripple_reflection) are shown
   against the per-mode theory (rayleigh_sponge_reflection).
2. How does the band reflected share evolve in time in each run?
3. Does the incident flux level keep decaying, and does the band drift settle?

Figures and a summary JSON go to --outdir. Usage

    python3 scripts/compare_sponge_reflection.py
    python3 scripts/compare_sponge_reflection.py --runs S0=<analysis dir> S1=<...>

With --cavity it instead runs the single-mode cavity closure test on runs made
with `sweep_arsia_mons_stage2.jl linear mode=N`. The ground condition fixes
Psi(0) = A0 + B0, so with r0 = B0/A0 the upward amplitude is A0 = Psi(0)/(1 + r0)
and the upward flux, relative to the flux the same forcing gives in a column with
no damping layer and a radiating top, should be |1 + r_open|^2 / |1 + r0|^2. The
measured level is set against that prediction with r0 read two ways: from the
split at the lowest level, and from the ripple of |Psi|^2 over 10-55 km carried
down to the ground. The second shares no algebra with the split. The linear
Rayleigh-layer theory of the run's own layer is reported beside both.

    python3 scripts/compare_sponge_reflection.py --cavity <run dir> [<run dir> ...]

With --background it differences the x-mean rho, T, u and N^2 of each run against its
initial column (the state the damping layer relaxes toward) at --background-times.

    python3 scripts/compare_sponge_reflection.py --background C0=<run dir> L1=<run dir> \
        --background-times 0.5 6 18 36

With --theory it evaluates the linear reflection of each single-mode run's own mode on
the initial column and on the adjusted x-mean columns at --theory-times, with and
without the run's layer and with and without the compressible term.

    python3 scripts/compare_sponge_reflection.py --theory <run dir> [<run dir> ...]

With --mass it audits the domain mass: its change, the exchange through the top and
bottom boundaries, and the source of the damping layer's relaxation of rho.

    python3 scripts/compare_sponge_reflection.py --mass C0=<run dir> [...]
"""

import argparse
import glob
import json
import os
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_linear_flux_constancy import (CP, G, rayleigh_sponge_reflection,  # noqa: E402
                                           updown_split)

PROFILE = "external/mcd_profiles/mars_arsia_mons_my31.csv"

A = "output/arsia_mons_linear15h_LX2000_sig50{}_analysis/run_U40.0_h0.03"
DEFAULT_RUNS = [
    ("S0", A.format(""), "60-90 km, K0 0.05"),
    ("S2", A.format("_K0_0.0158"), "60-90 km, K0 0.0158"),
    ("S1", A.format("_K0_0.005"), "60-90 km, K0 0.005"),
    ("T1", A.format("_K0_0.005_LZ120"), "60-120 km, K0 0.005"),
    ("T2", A.format("_LZ100_spb70"), "70-100 km, K0 0.05"),
]
# Fixed categorical order (validated for CVD separation on the light surface);
# the marker shape is the secondary encoding, so identity is never colour alone.
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
MARKERS = ["o", "s", "^", "D", "v"]
RETURNED_MODES = range(6, 13)   # round trip from a 60 or 70 km base elapsed by 11 h


def parse_args():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--runs", nargs="*", default=None,
                   help="label=analysis_dir pairs; default the five variants")
    p.add_argument("--outdir", default="output/reflection_tests")
    p.add_argument("--cavity", nargs="*", default=None,
                   help="single-mode run directories for the cavity closure test")
    p.add_argument("--cavity-window", type=float, default=0.5,
                   help="half-width [h] of the windows the closure is averaged over")
    p.add_argument("--cavity-snaps", type=int, default=7,
                   help="snapshots per window")
    p.add_argument("--background", nargs="*", default=None,
                   help="label=run_dir pairs for the x-mean background audit")
    p.add_argument("--background-times", nargs="*", type=float, default=[0.5, 6.0, 18.0],
                   help="times [h] of the audit profiles (nearest snapshot is used)")
    p.add_argument("--mass", nargs="*", default=None,
                   help="label=run_dir pairs for the domain mass audit")
    p.add_argument("--theory", nargs="*", default=None,
                   help="single-mode run directories: theory on the initial and adjusted columns")
    p.add_argument("--theory-times", nargs="*", type=float, default=[6.0, 18.0],
                   help="times [h] of the adjusted columns")
    p.add_argument("--background-levels", nargs="*", type=float,
                   default=[0.25, 5, 10, 20, 30, 40, 50, 55, 58, 60, 62, 65, 70, 80, 90,
                            100, 110, 119.75],
                   help="heights [km] printed in the audit tables")
    return p.parse_args()


def run_setup(run_dir):
    """Mode, amplitude and layer of a run, from the names the driver gives it."""
    parent = os.path.basename(os.path.dirname(os.path.normpath(run_dir)))
    name = os.path.basename(os.path.normpath(run_dir))
    grab = lambda pat, dflt: float(m.group(1)) if (m := re.search(pat, parent)) else dflt
    return {"mode": int(grab(r"_mode(\d+)", 0)),
            "h_m": 1e3 * float(re.search(r"_h([0-9.]+)$", name).group(1)),
            "u_bar": float(re.search(r"run_U([0-9.]+)_", name).group(1)),
            "K0": grab(r"_K0_([0-9.]+)", 0.05),
            "top_km": grab(r"_LZ(\d+)", 90.0),
            "base_km": grab(r"_spb(\d+)", 60.0)}


def ripple_complex(psi2, z, m, sel, iref=0):
    """Ripple fit of m|Psi|^2 returning the complex B/A at grid level iref.

    Same model as ripple_reflection(): with theta = int_{z_ref}^z m dz,
    m|Psi|^2 = a^2 + b^2 + 2ab cos(2 theta + arg A - arg B), so the fitted phase
    psi of the cos/sin pair is arg(B/A) at z_ref, and the depth gives |B/A|.
    """
    th = np.concatenate([[0.0], np.cumsum(0.5 * (m[1:] + m[:-1]) * np.diff(z))])
    th = (th - th[iref])[sel]
    zz = z[sel]
    y = (psi2 * m)[sel]
    X = np.column_stack([np.ones_like(zz), zz - zz.mean(), np.cos(2 * th), np.sin(2 * th)])
    c = np.linalg.lstsq(X, y, rcond=None)[0]
    d = min(np.hypot(c[2], c[3]) / max(c[0], 1e-300), 1.0)
    q = (1.0 - np.sqrt(1.0 - d * d)) / d if d > 0 else 0.0
    return q * np.exp(1j * np.arctan2(c[3], c[2]))


def lower_reflection(k, z, N, rho, U, iref):
    """Reflection of the column below level iref for a downward wave, A/B there.

    Integrates the undamped structure equation upward from a rigid ground
    (Phi(0) = 0, Phi'(0) = 1) and splits the result at z[iref] with the split's
    wavenumber. For a uniform column this is -1 at every level; the MCD column
    departs from it through the near-surface N(z) structure. It is a property of
    the lower column alone and carries no information about the damping layer.
    """
    zz = np.concatenate([[0.0], z[: iref + 1]])
    rr = np.interp(zz, z, rho, left=rho[0])
    NN = np.interp(zz, z, N, left=N[0])
    sr = -0.5 * np.gradient(np.log(rr), zz)
    L = k ** 2 * (NN ** 2 / (U * k) ** 2 - 1.0) - sr ** 2 + np.gradient(sr, zz)
    Phi, dPhi = 0.0 + 0j, 1.0 + 0j
    for j in range(len(zz) - 1):
        h = zz[j + 1] - zz[j]
        Lm = 0.5 * (L[j] + L[j + 1])
        k1p, k1d = dPhi, -L[j] * Phi
        k2p, k2d = dPhi + 0.5 * h * k1d, -Lm * (Phi + 0.5 * h * k1p)
        k3p, k3d = dPhi + 0.5 * h * k2d, -Lm * (Phi + 0.5 * h * k2p)
        k4p, k4d = dPhi + h * k3d, -L[j + 1] * (Phi + h * k3p)
        Phi = Phi + h / 6 * (k1p + 2 * k2p + 2 * k3p + k4p)
        dPhi = dPhi + h / 6 * (k1d + 2 * k2d + 2 * k3d + k4d)
    m = np.sqrt((N[iref] / U) ** 2 - k ** 2)
    return complex((Phi + dPhi / (1j * m)) / (Phi - dPhi / (1j * m)))


def cavity_closure(run_dir, win_h, nsnap, z_lo=10.0, z_hi=55.0):
    st = run_setup(run_dir)
    n, U = st["mode"], st["u_bar"]
    files = sorted(glob.glob(os.path.join(run_dir, "snapshot_*.nc")))
    t = []
    for f in files:
        with Dataset(f) as ds:
            t.append(float(ds.getncattr("time_s")) / 3600.0)
    t = np.array(t)
    with Dataset(files[0]) as ds:
        x = np.asarray(ds["x"][:], float)
        z = np.asarray(ds["z"][:], float)
    zkm, lx = z / 1e3, x.size * (x[1] - x[0])
    k = 2 * np.pi * n / lx
    sel = (zkm >= z_lo) & (zkm <= z_hi)

    raw = np.genfromtxt(PROFILE, delimiter=",", names=True)
    rho = np.interp(z, raw["z_m"], raw["rho"])
    T = np.interp(z, raw["z_m"], raw["T"])
    N = np.sqrt(np.clip((G / T) * (np.gradient(T, z) + G / CP), 1e-12, None))
    s_rho = -0.5 * np.gradient(np.log(rho), z)
    m_ex = np.sqrt(np.clip((N / U) ** 2 - k ** 2 - s_rho ** 2 + np.gradient(s_rho, z), 0, None))
    m_sp = np.sqrt(np.clip((N / U) ** 2 - k ** 2, 0, None))

    # Imposed lower boundary in the split's variable: w = U k h cos(kx) at z = 0,
    # whose one-sided rfft/nx coefficient is U k h / 2, scaled by sqrt(rho_0(0)).
    psi0 = np.sqrt(np.interp(0.0, raw["z_m"], raw["rho"])) * U * k * st["h_m"] / 2

    def column(K0, top_km, base_km):
        """Upward flux and ground reflection of the linear column, forced as the run."""
        R, Phi, dPhi = rayleigh_sponge_reflection(k, z, N, rho, U, K0, base_km * 1e3,
                                                  top_km * 1e3, top="radiation",
                                                  return_field=True)
        c = psi0 / (Phi[0] - z[0] * dPhi[0])
        Phi, dPhi = c * Phi, c * dPhi
        A = 0.5 * (Phi + dPhi / (1j * m_sp))
        B = 0.5 * (Phi - dPhi / (1j * m_sp))
        fup = 2 * (m_sp / k) * np.abs(A) ** 2
        return {"F_up": float(fup[sel].mean()), "R_band": float(np.nanmean(R[sel])),
                "r0": complex(B[0] / A[0]), "r_ref": complex(B[iref] / A[iref])}

    # Reference level for the second closure: the band base. There the lower
    # column's reflection r_low is taken from theory and the upper one, r_top,
    # from the run; the upward level is then 1/|1 - r_low r_top|^2 relative to a
    # column that returns nothing from above.
    iref = int(np.argmin(np.abs(zkm - z_lo)))
    r_low = lower_reflection(k, z, N, rho, U, iref)

    open_col = column(0.0, st["top_km"], st["base_km"])
    layer_col = column(st["K0"], st["top_km"], st["base_km"])

    rows = []
    for hc in np.arange(1.0, np.floor(t[-1]) + 1e-6, 1.0):
        idx = np.where(np.abs(t - hc) <= win_h)[0]
        if idx.size == 0:
            continue
        idx = idx[np.linspace(0, idx.size - 1, min(nsnap, idx.size)).astype(int)]
        acc = {"A": 0, "B": 0, "fup": 0, "fdn": 0, "psi2": 0, "net": 0}
        for i in idx:
            with Dataset(files[i]) as ds:
                u = np.asarray(ds["u"][:], float)
                w = np.asarray(ds["w"][:], float)
                r = np.asarray(ds["rho"][:], float)
                Ts = np.asarray(ds["T"][:], float).mean(1)
            Nl = np.sqrt(np.clip((G / Ts) * (np.gradient(Ts, z) + G / CP), 1e-12, None))
            sp = updown_split(u, w, r, x, z, Nl, u.mean(1))
            acc["A"] = acc["A"] + sp["A"][:, n]
            acc["B"] = acc["B"] + sp["B"][:, n]
            acc["fup"] = acc["fup"] + np.abs(sp["fup"][:, n])
            acc["fdn"] = acc["fdn"] + np.abs(sp["fdn"][:, n])
            acc["psi2"] = acc["psi2"] + np.abs(sp["A"][:, n] + sp["B"][:, n]) ** 2
            acc["net"] = acc["net"] + sp["mf_spec"]
        g = {q: v / idx.size for q, v in acc.items()}
        r_split = complex(g["B"][0] / g["A"][0])
        r_rip = ripple_complex(g["psi2"], z, m_ex, sel)
        rt_split = complex(g["B"][iref] / g["A"][iref])
        rt_rip = ripple_complex(g["psi2"], z, m_ex, sel, iref)
        F_up = float(g["fup"][sel].mean())
        pred = lambda r0: abs(1 + open_col["r0"]) ** 2 / abs(1 + r0) ** 2
        pred_ref = lambda rt: (abs(1 - r_low * open_col["r_ref"]) ** 2
                               / abs(1 - r_low * rt) ** 2)
        rows.append({"t_h": float(t[idx].mean()),
                     "R_split_band": float(np.mean(g["fdn"][sel] / g["fup"][sel])),
                     "R_split_ground": abs(r_split) ** 2,
                     "phase_split_ground_deg": float(np.degrees(np.angle(r_split))),
                     "R_ripple": abs(r_rip) ** 2,
                     "phase_ripple_ground_deg": float(np.degrees(np.angle(r_rip))),
                     "psi_ground_over_imposed": float(abs(g["A"][0] + g["B"][0]) / psi0),
                     "level_measured": F_up / open_col["F_up"],
                     "level_pred_split": pred(r_split),
                     "level_pred_ripple": pred(r_rip),
                     "R_top_ref": abs(rt_split) ** 2,
                     "phase_top_ref_deg": float(np.degrees(np.angle(rt_split))),
                     "phase_top_ref_ripple_deg": float(np.degrees(np.angle(rt_rip))),
                     "level_pred_ref_split": pred_ref(rt_split),
                     "level_pred_ref_ripple": pred_ref(rt_rip),
                     "net_ptp_pct": float(100 * np.ptp(g["net"][sel]) / abs(np.median(g["net"][sel])))})
    return {"run_dir": run_dir, "setup": st, "k_per_m": k,
            "lambda_x_km": lx / n / 1e3, "psi0_imposed": psi0,
            "open_column": {"F_up": open_col["F_up"], "R_band": open_col["R_band"],
                            "r0": [open_col["r0"].real, open_col["r0"].imag],
                            "r_ref": [open_col["r_ref"].real, open_col["r_ref"].imag]},
            "z_ref_km": float(zkm[iref]), "r_low": [r_low.real, r_low.imag],
            "layer_theory": {"R_band": layer_col["R_band"],
                             "R_ground": abs(layer_col["r0"]) ** 2,
                             "phase_ground_deg": float(np.degrees(np.angle(layer_col["r0"]))),
                             "R_ref": abs(layer_col["r_ref"]) ** 2,
                             "phase_ref_deg": float(np.degrees(np.angle(layer_col["r_ref"]))),
                             "level": layer_col["F_up"] / open_col["F_up"]},
            "hours": rows}


def cavity_main(a):
    os.makedirs(a.outdir, exist_ok=True)
    res = [cavity_closure(d, a.cavity_window, a.cavity_snaps) for d in a.cavity]
    for r in res:
        st, th = r["setup"], r["layer_theory"]
        print(f"\n{r['run_dir']}\n  mode {st['mode']} (lambda_x {r['lambda_x_km']:.0f} km), "
              f"h {st['h_m']:.1f} m, layer {st['base_km']:.0f}-{st['top_km']:.0f} km, "
              f"K0 {st['K0']:g}")
        print(f"  linear theory of this layer: R band {th['R_band']:.4f}, R ground "
              f"{th['R_ground']:.4f}, phase {th['phase_ground_deg']:.0f} deg, level "
              f"{th['level']:.4f};  open column R band {r['open_column']['R_band']:.5f}")
        rl = complex(*r["r_low"])
        print(f"  at z_ref = {r['z_ref_km']:.2f} km: r_low = {abs(rl):.3f} at "
              f"{np.degrees(np.angle(rl)):.0f} deg; theory r_top^2 {th['R_ref']:.4f} at "
              f"{th['phase_ref_deg']:.0f} deg")
        print(f"{'t[h]':>6} {'R_band':>7} {'R_rip':>7} {'R_gnd':>7} {'ph_gnd':>7} "
              f"{'Psi0/imp':>8} {'R_ref':>7} {'ph_ref':>7} {'ph_ref_r':>8} {'level':>7} "
              f"{'gnd_sp':>7} {'ref_sp':>7} {'ref_rip':>7} {'net_ptp%':>8}")
        for h in r["hours"]:
            print(f"{h['t_h']:6.2f} {h['R_split_band']:7.4f} {h['R_ripple']:7.4f} "
                  f"{h['R_split_ground']:7.4f} {h['phase_split_ground_deg']:7.0f} "
                  f"{h['psi_ground_over_imposed']:8.4f} {h['R_top_ref']:7.4f} "
                  f"{h['phase_top_ref_deg']:7.0f} {h['phase_top_ref_ripple_deg']:8.0f} "
                  f"{h['level_measured']:7.4f} {h['level_pred_split']:7.4f} "
                  f"{h['level_pred_ref_split']:7.4f} {h['level_pred_ref_ripple']:7.4f} "
                  f"{h['net_ptp_pct']:8.2f}")
    with open(os.path.join(a.outdir, "cavity_closure.json"), "w") as fh:
        json.dump(res, fh, indent=1, default=float)

    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    for r, c, mk in zip(res, COLORS, MARKERS):
        st = r["setup"]
        lab = f"{st['base_km']:.0f}-{st['top_km']:.0f} km, K0 {st['K0']:g}"
        hh = r["hours"]
        tt = [h["t_h"] for h in hh]
        ax[0].plot(tt, [h["level_measured"] for h in hh], "-", color=c, lw=2,
                   label=f"{lab}, measured")
        ax[0].plot(tt, [h["level_pred_split"] for h in hh], "--", color=c, lw=1.5,
                   label=f"{lab}, from r0 (split)")
        ax[0].plot(tt, [h["level_pred_ref_ripple"] for h in hh], "none", marker=mk,
                   mec=c, mfc="none", ms=7, label=f"{lab}, from r_top at 10 km (ripple)")
        ax[1].plot(tt, [h["R_split_band"] for h in hh], "-", color=c, lw=2, label=lab)
        ax[1].plot(tt, [h["R_ripple"] for h in hh], "none", marker=mk, mec=c, mfc="none",
                   ms=7)
        ax[1].axhline(r["layer_theory"]["R_band"], color=c, ls=":", lw=1.5)
    ax[0].set(xlabel="time [h]", ylabel="upward flux / open-column value, 10-55 km",
              title="Upward level: measured (solid) and cavity predictions")
    ax[1].set(xlabel="time [h]", ylabel=r"$|B/A|^2$", yscale="log",
              title="Reflected share: split (solid), ripple (open), theory (dotted)")
    ax[0].legend(fontsize=8)
    ax[1].legend(fontsize=8)
    for x_ in ax:
        x_.grid(True, color="0.9", lw=0.8)
        for sp_ in ("top", "right"):
            x_.spines[sp_].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(a.outdir, "cavity_closure.png"), dpi=150)
    plt.close(fig)
    print(f"\nFigure and JSON in {a.outdir}")


def snapshot_near(run_dir, t_h):
    """The snapshot whose time is closest to t_h, from the step numbers in the names."""
    files = sorted(glob.glob(os.path.join(run_dir, "snapshot_*.nc")))
    steps = np.array([int(re.search(r"_(\d+)\.nc$", f).group(1)) for f in files])
    with Dataset(files[1]) as d:
        dt_snap = float(d.time_s) / steps[1]
    return files[int(np.argmin(np.abs(steps * dt_snap - 3600 * t_h)))]


def xmean_column(path):
    """x-mean rho, T, u, P and N^2 of one snapshot, with its time [h]."""
    with Dataset(path) as d:
        g, gam, t = float(d.g), float(d.gamma), float(d.time_s) / 3600
        z = d["z"][:].data.astype(float)
        col = {k: d[k][:].data.astype(float).mean(axis=1) for k in ("rho", "T", "u", "P")}
    col["N2"] = g * np.gradient(np.log(col["T"]) - (gam - 1) / gam * np.log(col["P"]), z)
    col["c2"] = gam * col["P"] / col["rho"]
    return t, z, col


def background_audit(run_dir, times_h):
    """x-mean departure of the column from its initial state at the requested times."""
    _, z, c0 = xmean_column(os.path.join(run_dir, "snapshot_0000000.nc"))
    rows = []
    for th in times_h:
        t, _, c = xmean_column(snapshot_near(run_dir, th))
        if rows and t - rows[-1]["t_h"] < 0.5:     # a time past the end of the record
            continue
        zk = z / 1e3
        g5 = (zk > 5) & (zk < 50)
        base = (zk >= 59) & (zk <= 61)
        rel = c["rho"] / c0["rho"] - 1
        rows.append({"t_h": t,
                     "dlnrho_slope_pct_per_km": float(100 * np.polyfit(zk[g5], rel[g5], 1)[0]),
                     "N2_min_ratio_59_61km": float(np.min(c["N2"][base] / c0["N2"][base])),
                     "drho_pct": 100 * rel,
                     "dT_pct": 100 * (c["T"] / c0["T"] - 1),
                     "du_ms": c["u"] - c0["u"],
                     "dN2_pct": 100 * (c["N2"] / c0["N2"] - 1)})
    return {"run_dir": run_dir, "setup": run_setup(run_dir), "z_km": z / 1e3,
            "N2_initial": c0["N2"], "rows": rows}


def background_main(a):
    os.makedirs(a.outdir, exist_ok=True)
    res = []
    for spec in a.background:
        lab, run_dir = spec.split("=", 1)
        r = background_audit(run_dir, a.background_times)
        r["label"] = lab
        res.append(r)
        zk, st = r["z_km"], r["setup"]
        print(f"\n{lab}: {run_dir}  layer {st['base_km']:.0f}-{st['top_km']:.0f} km, "
              f"K0 {st['K0']:g}")
        for row in r["rows"]:
            print(f"  t = {row['t_h']:.2f} h   d(ln rho)/dz change 5-50 km "
                  f"{row['dlnrho_slope_pct_per_km']:+.4f} %/km,   min N2/N2_0 over 59-61 km "
                  f"{row['N2_min_ratio_59_61km']:.3f}")
            print(f"  {'z[km]':>7} {'drho%':>8} {'dT%':>8} {'du m/s':>8} {'dN2%':>8}")
            for h in a.background_levels:
                if h > zk.max():
                    continue
                i = int(np.argmin(np.abs(zk - h)))
                print(f"  {zk[i]:7.2f} {row['drho_pct'][i]:+8.3f} {row['dT_pct'][i]:+8.3f} "
                      f"{row['du_ms'][i]:+8.3f} {row['dN2_pct'][i]:+8.2f}")
    with open(os.path.join(a.outdir, "background_audit.json"), "w") as fh:
        json.dump(res, fh, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist")
                  else float(o))

    # Small multiples: one column per run, one row per quantity. Time is a magnitude,
    # so the curves share one sequential hue that darkens with time.
    qty = [("drho_pct", r"$\Delta\rho/\rho_0$ [%]"), ("dT_pct", r"$\Delta T/T_0$ [%]"),
           ("du_ms", r"$\Delta u$ [m/s]"), ("dN2_pct", r"$\Delta N^2/N^2_0$ [%]")]
    fig, ax = plt.subplots(len(qty), len(res), figsize=(4.2 * len(res), 13),
                           sharey=True, squeeze=False)
    for j, r in enumerate(res):
        shades = plt.cm.Blues(np.linspace(0.4, 0.95, len(r["rows"])))
        for i, (k, lab) in enumerate(qty):
            x_ = ax[i, j]
            for row, c in zip(r["rows"], shades):
                x_.plot(row[k], r["z_km"], color=c, lw=1.6, label=f"{row['t_h']:.1f} h")
            x_.axvline(0, color="0.6", lw=0.8)
            x_.axhspan(r["setup"]["base_km"], r["setup"]["top_km"], color="0.93", zorder=0)
            x_.set_xlabel(lab)
            x_.grid(True, color="0.9", lw=0.8)
            for sp_ in ("top", "right"):
                x_.spines[sp_].set_visible(False)
            if k == "dN2_pct":
                x_.set_xlim(-100, 30)
        ax[0, j].set_title(f"{r['label']}, layer {r['setup']['base_km']:.0f}-"
                           f"{r['setup']['top_km']:.0f} km (shaded), K0 {r['setup']['K0']:g}",
                           fontsize=10)
        ax[0, j].legend(fontsize=8, loc="lower left")
    for i in range(len(qty)):
        ax[i, 0].set_ylabel("z [km]")
    fig.tight_layout()
    fig.savefig(os.path.join(a.outdir, "background_audit.png"), dpi=150)
    plt.close(fig)
    print(f"\nFigure and JSON in {a.outdir}")


def mass_audit(run_dir):
    """Domain mass and its exchange through the top, the bottom and the damping layer.

    Per unit length in y. The top exchange is the x-integrated rho*w at the highest
    level and the layer source is sum K (rho_bar - rho), with K(z) rebuilt from the
    run's layer; both are integrated in time by the trapezoid rule over the snapshot
    cadence. x is periodic and the density floor does not act, so the only other
    path is the bottom face. Its flux is not recoverable from cell-centre snapshots
    (the mirror condition cancels rho*w at the face, but the Lax-Wendroff half step
    adds a term from the vertical flux difference across the face), so it is
    reported by elimination: bottom_out = -(dM + top_out - layer_in).
    """
    st = run_setup(run_dir)
    files = sorted(glob.glob(os.path.join(run_dir, "snapshot_*.nc")))
    with Dataset(files[0]) as d:
        x = np.asarray(d["x"][:], float)
        z = np.asarray(d["z"][:], float)
        rho0 = np.asarray(d["rho"][:], float)
    dx, dz = x[1] - x[0], z[1] - z[0]
    base, top = st["base_km"] * 1e3, st["top_km"] * 1e3
    eta = np.clip((z - base) / (top - base), 0.0, 1.0)
    K = np.where(z >= base, st["K0"] * np.sin(0.5 * np.pi * eta) ** 2, 0.0)
    t, M, Ftop, S = [], [], [], []
    for f in files:
        with Dataset(f) as d:
            rho = np.asarray(d["rho"][:], float)
            w = np.asarray(d["w"][:], float)
            t.append(float(d.time_s))
        M.append(rho.sum() * dx * dz)
        Ftop.append((rho[-1] * w[-1]).sum() * dx)        # outward positive
        S.append((K[:, None] * (rho0 - rho)).sum() * dx * dz)
    t, M, Ftop, S = map(np.array, (t, M, Ftop, S))
    cum = lambda y: np.concatenate([[0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(t))])
    out = {"t_h": t / 3600, "dM_pct": 100 * (M / M[0] - 1),
           "top_out_pct": 100 * cum(Ftop) / M[0], "sponge_in_pct": 100 * cum(S) / M[0]}
    out["bottom_out_pct"] = -(out["dM_pct"] + out["top_out_pct"] - out["sponge_in_pct"])
    return {"run_dir": run_dir, "setup": st, **{k: v.tolist() for k, v in out.items()}}


def mass_main(a):
    os.makedirs(a.outdir, exist_ok=True)
    res = []
    for spec in a.mass:
        lab, run_dir = spec.split("=", 1)
        r = mass_audit(run_dir)
        r["label"] = lab
        res.append(r)
        print(f"\n{lab}: {run_dir}")
        print(f"  {'t[h]':>6} {'dM%':>9} {'top out%':>9} {'layer in%':>10} "
              f"{'bottom out% (by elimination)':>29}")
        t = np.array(r["t_h"])
        for th in [0.1, 0.25, 0.5, 1, 2, 3, 4, 6, 9, 12, 15, 18, 24, 30, 36]:
            if th > t[-1] + 0.05:
                continue
            i = int(np.argmin(np.abs(t - th)))
            print(f"  {t[i]:6.2f} {r['dM_pct'][i]:+9.4f} {r['top_out_pct'][i]:+9.4f} "
                  f"{r['sponge_in_pct'][i]:+10.4f} {r['bottom_out_pct'][i]:+29.4f}")
    with open(os.path.join(a.outdir, "mass_audit.json"), "w") as fh:
        json.dump(res, fh)
    print(f"\nJSON in {a.outdir}")


def adjusted_theory(run_dir, times_h, z_lo=10.0, z_hi=55.0):
    """Linear reflected share of the run's own mode on the initial and adjusted columns.

    For each column the theory is evaluated with no damping (K0 = 0, the column's
    own structure alone) and with the run's layer, each in the anelastic form and
    with the compressible term. The top is radiating in every case.
    """
    st = run_setup(run_dir)
    n, U = st["mode"], st["u_bar"]
    t0, z, c0 = xmean_column(os.path.join(run_dir, "snapshot_0000000.nc"))
    with Dataset(os.path.join(run_dir, "snapshot_0000000.nc")) as d:
        x = np.asarray(d["x"][:], float)
    k = 2 * np.pi * n / (x.size * (x[1] - x[0]))
    zkm = z / 1e3
    sel = (zkm >= z_lo) & (zkm <= z_hi)
    base, top = st["base_km"] * 1e3, st["top_km"] * 1e3

    cols = [("initial", t0, c0)]
    for th in times_h:
        t, _, c = xmean_column(snapshot_near(run_dir, th))
        cols.append((f"adjusted {t:.1f} h", t, c))
    out = []
    for name, t, c in cols:
        N = np.sqrt(np.clip(c["N2"], 1e-12, None))
        row = {"column": name, "t_h": t}
        for K0 in (0.0, st["K0"]):
            for comp in (False, True):
                R = rayleigh_sponge_reflection(k, z, N, c["rho"], U, K0, base, top,
                                               top="radiation", c2=c["c2"] if comp else None)
                key = f"K0={K0:g}{' comp' if comp else ''}"
                row[key] = float(np.nanmean(R[sel]))
                row[key + " R(z)"] = {f"{h:g}": float(R[int(np.argmin(np.abs(zkm - h)))])
                                      for h in (10, 30, 50, 55, 58)}
        out.append(row)
    return {"run_dir": run_dir, "setup": st, "k_per_m": k, "columns": out}


def theory_main(a):
    os.makedirs(a.outdir, exist_ok=True)
    meas = {}
    f = os.path.join(a.outdir, "cavity_closure.json")
    if os.path.isfile(f):
        with open(f) as fh:
            for r in json.load(fh):
                hh = [h for h in r["hours"] if 12 <= h["t_h"] <= 17.5]
                meas[os.path.normpath(r["run_dir"])] = (
                    np.mean([h["R_split_band"] for h in hh]),
                    np.mean([h["R_ripple"] for h in hh]),
                    {round(h["t_h"]): (h["R_split_band"], h["R_ripple"]) for h in r["hours"]})
    res = []
    for d in a.theory:
        r = adjusted_theory(d, a.theory_times)
        st = r["setup"]
        m = meas.get(os.path.normpath(d))
        r["measured_12_17h"] = m[:2] if m else None
        res.append(r)
        print(f"\n{d}\n  mode {st['mode']}, layer {st['base_km']:.0f}-{st['top_km']:.0f} km, "
              f"K0 {st['K0']:g}" + (f";  measured 12-17 h: split {m[0]:.4f}, ripple {m[1]:.4f}"
                                    if m else ""))
        keys = [k_ for k_ in r["columns"][0] if k_.startswith("K0=") and "R(z)" not in k_]
        print(f"  {'column':<16}" + "".join(f"{k_:>16}" for k_ in keys))
        for c in r["columns"]:
            print(f"  {c['column']:<16}" + "".join(f"{c[k_]:16.5f}" for k_ in keys))
        if m:
            own = f"K0={st['K0']:g} comp"
            print(f"  {'t[h]':>6} {'K0=0':>8} {own:>16} {'meas split':>11} {'meas ripple':>12}")
            for c in r["columns"][1:]:
                hr = m[2].get(round(c["t_h"]), (np.nan, np.nan))
                print(f"  {c['t_h']:6.1f} {c['K0=0']:8.4f} {c[own]:16.4f} {hr[0]:11.4f} {hr[1]:12.4f}")
        print("  R(z) at 10/30/50/55/58 km:")
        for c in r["columns"]:
            for k_ in keys:
                v = c[k_ + " R(z)"]
                print(f"  {c['column']:<16} {k_:<14}" + " ".join(f"{v[h]:8.4f}" for h in v))
    with open(os.path.join(a.outdir, "adjusted_theory.json"), "w") as fh:
        json.dump(res, fh, indent=1, default=float)
    print(f"\nJSON in {a.outdir}")


def load(runs):
    got = []
    for lab, path, desc in runs:
        f = os.path.join(path, "linear_flux_constancy.json")
        if not os.path.isfile(f):
            print(f"  {lab}: no {f}, skipped")
            continue
        with open(f) as fh:
            got.append((lab, desc, json.load(fh)))
    return got


def final_hour_modes(d):
    rip = d["updown"]["ripple_permode"]
    tmax = max(e["t_h"] for e in rip)
    return {e["n"]: e for e in rip if e["t_h"] == tmax}, tmax


def trace(d):
    """The half-hourly up/down trace as arrays keyed like the JSON rows."""
    rows = d["updown"]["trace"]
    return {k: np.array([r[k] for r in rows]) for k in rows[0]}


def class_steady(d, cls, z):
    for e in d["transient"]["class_split"]:
        if e["class"] == cls and abs(e["z_km"] - z) < 0.3:
            return e["R_steady"]
    return np.nan


def main():
    a = parse_args()
    if a.cavity:
        return cavity_main(a)
    if a.background:
        return background_main(a)
    if a.theory:
        return theory_main(a)
    if a.mass:
        return mass_main(a)
    runs = DEFAULT_RUNS
    if a.runs:
        runs = [(r.split("=", 1)[0], r.split("=", 1)[1], "") for r in a.runs]
    data = load(runs)
    if not data:
        raise SystemExit("no runs found")
    os.makedirs(a.outdir, exist_ok=True)

    # ---- Table ------------------------------------------------------------
    summary = {}
    print(f"\nFinal-hour band-mean reflected share, modes n = {RETURNED_MODES.start}-"
          f"{RETURNED_MODES.stop - 1} (round trip elapsed)")
    print(f"{'run':>4} {'layer':>20} {'t[h]':>6} {'split':>8} {'ripple':>8} "
          f"{'theory':>8} {'split/th':>9} {'ret@9.75':>9} {'noret@9.75':>11} "
          f"{'band R':>7} {'drift':>8}")
    base = None
    for lab, desc, d in data:
        pm, tmax = final_hour_modes(d)
        sel = [pm[n] for n in RETURNED_MODES if n in pm]
        rs = float(np.mean([e["R_split"] for e in sel]))
        rr = float(np.mean([e["R_ripple"] for e in sel]))
        rt = float(np.mean([e["R_theory"] for e in sel]))
        tr = trace(d)
        cb = d["constancy_by_hour"]
        last = max(cb, key=lambda k: float(k.rstrip("h")))
        s = {"layer": desc, "t_h": tmax, "R_split": rs, "R_ripple": rr, "R_theory": rt,
             "R_split_by_mode": {e["n"]: e["R_split"] for e in sel},
             "R_ripple_by_mode": {e["n"]: e["R_ripple"] for e in sel},
             "R_theory_by_mode": {e["n"]: e["R_theory"] for e in sel},
             "returned_steady_9.75": class_steady(d, "returned", 9.75),
             "no_return_steady_9.75": class_steady(d, "no return", 9.75),
             "band_R_final": float(tr["R"][-1]),
             "drift_final_pct_per_km": cb[last]["drift_pct_per_km"],
             "ptp_final": cb[last]["ptp_rel"],
             "hourly_bins": d["transient"].get("hourly_summary", {}),
             "min_u_over_U": d["linearity"].get("min_u_over_U"),
             "sponge_theory": d["updown"].get("sponge_theory")}
        if base is None:
            base = s
        s["split_over_S0"] = rs / base["R_split"]
        s["theory_over_S0"] = rt / base["R_theory"]
        summary[lab] = s
        print(f"{lab:>4} {desc:>20} {tmax:6.2f} {rs:8.4f} {rr:8.4f} {rt:8.4f} "
              f"{rs/rt:9.3f} {s['returned_steady_9.75']:9.4f} "
              f"{s['no_return_steady_9.75']:11.4f} {s['band_R_final']:7.4f} "
              f"{s['drift_final_pct_per_km']:+8.4f}")
    print("\nRatio to the first run, measured (split) against theory:")
    for lab, s in summary.items():
        print(f"  {lab}: measured {s['split_over_S0']:.3f}   theory {s['theory_over_S0']:.3f}")

    with open(os.path.join(a.outdir, "compare_sponge_reflection.json"), "w") as fh:
        json.dump(summary, fh, indent=1, default=float)

    # ---- Figure 1: per-mode reflected share against theory -----------------
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.2))
    for (lab, desc, d), c, mk in zip(data, COLORS, MARKERS):
        pm, _ = final_hour_modes(d)
        ns = sorted(pm)
        lam = np.array([pm[n]["lambda_x_km"] for n in ns])
        est = np.array([pm[n]["established"] for n in ns])
        ax[0].plot(lam, [pm[n]["R_theory"] for n in ns], "--", color=c, lw=1.5)
        ax[0].plot(lam[est], np.array([pm[n]["R_split"] for n in ns])[est], "-",
                   marker=mk, color=c, lw=2, ms=8, label=f"{lab}  {desc}")
        ax[0].plot(lam[est], np.array([pm[n]["R_ripple"] for n in ns])[est], "none",
                   marker=mk, mec=c, mfc="none", ms=11, mew=1.5)
        s = summary[lab]
        ax[1].plot(s["R_theory"], s["R_split"], marker=mk, color=c, ms=10, ls="none",
                   label=lab)
        ax[1].plot(s["R_theory"], s["R_ripple"], marker=mk, mec=c, mfc="none",
                   ms=13, mew=1.5, ls="none")
        ax[1].annotate(lab, (s["R_theory"], s["R_split"]), xytext=(6, 4),
                       textcoords="offset points", fontsize=10, color="0.25")
    ax[0].set(xscale="log", yscale="log", xlabel=r"$\lambda_x$ [km]",
              ylabel=r"reflected share $|B/A|^2$, 10-55 km mean",
              title="Per mode: split (filled), ripple (open), theory (dashed)")
    ax[0].set_xticks([167, 250, 400, 667, 1000, 2000])
    ax[0].set_xticklabels(["167", "250", "400", "667", "1000", "2000"])
    ax[0].minorticks_off()
    ax[0].legend(fontsize=9, loc="lower left")
    lim = [1e-3, 1.0]
    ax[1].plot(lim, lim, color="0.6", lw=1)
    ax[1].set(xscale="log", yscale="log", xlim=lim, ylim=[1e-3, 0.3],
              xlabel="predicted, linear Rayleigh layer",
              ylabel="measured, modes 6-12",
              title="Mean over the returned modes (1:1 line grey)")
    for x in ax:
        x.grid(True, which="major", color="0.9", lw=0.8)
        for sp in ("top", "right"):
            x.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(a.outdir, "reflection_vs_theory.png"), dpi=150)
    plt.close(fig)

    # ---- Figure 2: band reflected share in time ----------------------------
    fig, ax = plt.subplots(figsize=(8, 4.8))
    for (lab, desc, d), c, mk in zip(data, COLORS, MARKERS):
        tr = trace(d)
        ax.plot(tr["t_h"], tr["R"], color=c, lw=2, label=f"{lab}  {desc}")
        ax.plot(tr["t_h"][-1], tr["R"][-1], marker=mk, color=c, ms=8)
    ax.set(xlabel="time [h]", ylabel="band R = |MF_down| / |MF_up|, 10-55 km",
           title="Reflected share of the band flux")
    ax.legend(fontsize=9)
    ax.grid(True, color="0.9", lw=0.8)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(a.outdir, "reflection_in_time.png"), dpi=150)
    plt.close(fig)

    # ---- Figure 3: incident level and band drift ---------------------------
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    for (lab, desc, d), c, mk in zip(data, COLORS, MARKERS):
        tr = trace(d)
        lin = abs(d["linear_steady_MF_uPa"])
        ax[0].plot(tr["t_h"], np.abs(tr["MF_up_uPa"]) / lin, color=c, lw=2,
                   label=f"{lab}  {desc}")
        cb = d["constancy_by_hour"]
        hs = sorted(cb, key=lambda k: float(k.rstrip("h")))
        ax[1].plot([float(h.rstrip("h")) for h in hs],
                   [cb[h]["drift_pct_per_km"] for h in hs], "-", marker=mk,
                   color=c, lw=2, ms=6, label=lab)
    ax[0].set(xlabel="time [h]", ylabel=r"$|MF_{up}|$ / steady linear value",
              title="Incident flux level, 10-55 km")
    ax[1].axhspan(-0.1, 0.1, color="0.92", zorder=0)
    ax[1].set(xlabel="time [h]", ylabel="band drift [% per km]", ylim=[-1.0, 0.5],
              title="Net-flux drift (grey band: 0.1 %/km threshold)")
    ax[0].legend(fontsize=9)
    for x in ax:
        x.grid(True, color="0.9", lw=0.8)
        for sp in ("top", "right"):
            x.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(a.outdir, "incident_level_and_drift.png"), dpi=150)
    plt.close(fig)
    print(f"\nFigures and summary in {a.outdir}")


if __name__ == "__main__":
    main()
