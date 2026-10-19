#!/usr/bin/env python3
"""
analyze_ext6h_mechanism.py

Mechanism attribution for an Arsia Mons Stage 1 time extension (U=40, h=3 km).

The sampled hours are derived from the series itself, so the same script serves
the 6 h, 12 h or any other extension horizon without edits; pass --npz.

Answers, from the cached deep-dive series plus the run's own background column
and solver constants, four questions about the descending static-instability
layer documented in docs/TODO_tracker.md ("Arsia Mons Stage 1 U=40/h=3 time
extension"):

  1. WKB      - does the wave amplitude grow at the conservative (rho N)^-1/2 rate?
  2. ONSET    - where and when is the stationary-wave overturning criterion
                u_hat >= U (equivalently min_x u <= 0) first met, level by level?
  3. BUDGET   - does the mean-flow deceleration equal the resolved wave-momentum
                flux divergence -(1/rho) d<rho u'w'>_x/dz ?
  4. DESCENT  - does the observed layer-base descent match (a) the Lindzen
                saturation level |MF| = F_sat = rho k U^3 / (2N) and (b) the
                kinematic contour speed -(du/dt)/(du/dz) ?

plus two artifact tests: the sponge-pinned reference wind above 60 km, and the
horizontal/vertical resolution of the flux-carrying scale.

Reads only `deepdive_series.npz` (written by analyze_arsia_stage1_deepdive.py)
and the MCD background CSV, so it runs in ~1 s and never touches the NetCDF.

Usage:
  python3 scripts/analyze_ext6h_mechanism.py
  python3 scripts/analyze_ext6h_mechanism.py --npz <path> --profile <csv> --json <out.json>
"""
import argparse
import json
import os

import numpy as np

# --- solver constants, mirrored from src/parameters.jl mars_config() ---
GAMMA = 1.31
G = 3.72
R_GAS = 8.314462 / 44.01e-3          # ~188.9 J/(kg K)
CP = GAMMA * R_GAS / (GAMMA - 1.0)

# --- run configuration, mirrored from scripts/sweep_arsia_mons_stage1.jl ---
U_BAR = 40.0
H_MAX = 3000.0
SIGMA = 100e3
LX, DX, NX = 2000e3, 4000.0, 500
DZ = 500.0
Z_SPONGE, Z_TOP, K0_SPONGE = 60e3, 90e3, 0.05

DEF_NPZ = "output/arsia_mons_stage1_ext6h_analysis/run_U40.0_h3.0/deepdive_series.npz"
DEF_CSV = "external/mcd_profiles/mars_arsia_mons_my31.csv"


def ridge_spectrum(h_max=H_MAX, sigma=SIGMA, lx=LX, dx=DX, nx=NX):
    """Flux-weighted mean horizontal wavenumber of the Gaussian ridge.

    The hydrostatic momentum-flux spectrum of a 2D ridge weights |h(k)|^2 by |k|,
    so the single wavenumber that best represents the drag is the k-weighted mean
    of that distribution rather than the spectral peak of the terrain itself.

    The geometry defaults to this run's configuration but is parameterized so the
    linear control corner (h_max = 100 m) can reuse it; see
    scripts/analyze_linear_flux_constancy.py. Note that k_eff depends on lx only
    through the discrete spectrum, whereas the returned hk, and hence the flux
    built from it in linear_steady_flux, is a domain mean and scales as 1/lx.
    """
    x = (np.arange(nx) + 0.5) * dx
    h = h_max * np.exp(-((x - 0.5 * lx) ** 2) / (2 * sigma ** 2))
    hk = np.fft.rfft(h) / nx
    kk = 2 * np.pi * np.fft.rfftfreq(nx, dx)
    wgt = kk * np.abs(hk) ** 2
    return float((wgt * kk).sum() / wgt.sum()), hk, kk


def linear_steady_flux(rho, N, hk, kk, z, u_bar=U_BAR, nx=NX):
    """Domain-mean <rho u'w'>_x of the steady linear hydrostatic WKB solution.

    zeta ~ (rho N)^-1/2 (wave-action conservation), u' = -U dzeta/dz, w' = U dzeta/dx.
    Non-dissipative, so the result is height-independent; it is the flux the run
    would carry if the stationary wave were fully established and unsaturated.
    """
    m = N / u_bar
    amp = np.sqrt((rho[0] * N[0]) / (rho * N))
    phase = np.concatenate([[0.0], np.cumsum(0.5 * (m[1:] + m[:-1]) * np.diff(z))])
    out = np.empty(len(z))
    for j in range(len(z)):
        spec = hk * amp[j] * np.exp(1j * phase[j])
        u_p = np.fft.irfft(spec * (-u_bar * 1j * m[j]) * nx, nx)
        w_p = np.fft.irfft(spec * (1j * kk * u_bar) * nx, nx)
        out[j] = rho[j] * np.mean(u_p * w_p)
    return out


def hour_grid(th, step=1.0, start=1.0):
    """Hours actually covered by the series `th`, at `step` spacing.

    Derived from the series rather than hardcoded so one script serves any
    extension horizon. Nothing past th[-1] is returned, which matters because
    `it()` below clamps to the nearest index: an out-of-range request would
    otherwise relabel the final column with an hour the run never reached, for
    instance after an early abort on the dt-collapse guard.
    """
    # One snapshot interval of slack, so a run that stops at 5.991 h still
    # reports a "6h" row; anything further out would be extrapolation.
    tol = float(np.median(np.diff(th))) if len(th) > 1 else 0.0
    t_end = float(th[-1]) + tol
    n = int(np.floor((t_end - start) / step + 1e-9))
    return [round(start + i * step, 6) for i in range(max(0, n + 1))]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", default=DEF_NPZ)
    ap.add_argument("--profile", default=DEF_CSV)
    ap.add_argument("--json", default=None, help="write the summary metrics here")
    args = ap.parse_args()

    d = np.load(args.npz, allow_pickle=True)
    t = d["times"]
    th = t / 3600.0
    zkm = d["z_km"]
    z = zkm * 1e3
    ubar, mf = d["ubar_zt"], d["mf_zt"]
    ufrac, wrms, umin = d["unstable_frac_zt"], d["wrms_hill_zt"], d["umin_hill_zt"]

    raw = np.genfromtxt(args.profile, delimiter=",", names=True)
    rho = np.interp(z, raw["z_m"], raw["rho"])
    T = np.interp(z, raw["z_m"], raw["T"])
    N = np.sqrt(np.clip((G / T) * (np.gradient(T, z) + G / CP), 1e-12, None))
    H = R_GAS * T / G

    iz = lambda v: int(np.argmin(np.abs(zkm - v)))
    it = lambda v: int(np.argmin(np.abs(th - v)))
    hl = lambda h: f"{h:g}h"
    clean = (zkm >= 2.0) & (zkm < 60.0)

    # Sample schedules, all clipped to the series (see hour_grid).
    T_END_H = float(th[-1])
    H_FROM1 = hour_grid(th, 1.0, 1.0)     # WKB, equilibrium height
    H_FROM2 = hour_grid(th, 1.0, 2.0)     # momentum budget
    H_FROM3 = hour_grid(th, 1.0, 3.0)     # deposition, resolution
    H_FROM4 = hour_grid(th, 1.0, 4.0)     # instability-per-slope outlier test
    H_HALF3 = hour_grid(th, 0.5, 3.0)     # descent
    H_REF = min(3.0, T_END_H)             # reference hour for the onset table

    k_eff, hk, kk = ridge_spectrum()
    mf_lin = linear_steady_flux(rho, N, hk, kk, z)
    F_sat = lambda j, i: rho[j] * k_eff * ubar[j, i] ** 3 / (2 * N[j])

    out = {"k_eff_per_m": k_eff, "lambda_x_eff_km": 2 * np.pi / k_eff / 1e3,
           "linear_steady_MF_mPa": float(mf_lin[iz(10)] * 1e3)}

    hdr = lambda s: print("\n" + "=" * 82 + f"\n{s}\n" + "=" * 82)

    print(f"Run: {d['label']}   U_bar={U_BAR:.0f} m/s  h_max={H_MAX/1e3:.0f} km  "
          f"Nh/U={0.018*H_MAX/U_BAR:.2f}")
    print(f"Ridge: sigma={SIGMA/1e3:.0f} km -> flux-weighted k_eff={k_eff:.3e} 1/m "
          f"(lambda_x={2*np.pi/k_eff/1e3:.0f} km)")
    print(f"Steady linear hydrostatic flux: {mf_lin[iz(10)]*1e3:.2f} mPa "
          f"(height-independent to {np.ptp(mf_lin)/abs(mf_lin[0])*100:.2e} %)")

    # ---------------- 1. WKB amplitude growth ----------------
    hdr("1. WKB GROWTH  -  simulated w_rms(z) vs conservative (rho N)^-1/2")
    wkb = 1.0 / np.sqrt(rho * N)
    print(f"{'t[h]':>5} {'band':>14} {'sim growth':>11} {'WKB growth':>11} {'ratio':>7}")
    # Restricted to 1-3 h on purpose: agreement with the conservative rate is
    # only expected before overturning begins in the 45-59 km band.
    for hh in [h for h in H_FROM1 if h <= 3.0]:
        i = it(hh)
        for lo, hi in [(10, 30), (30, 45), (45, 59)]:
            g_sim = wrms[iz(hi), i] / wrms[iz(lo), i]
            g_wkb = wkb[iz(hi)] / wkb[iz(lo)]
            print(f"{th[i]:5.2f} {f'{lo}-{hi} km':>14} {g_sim:11.3f} {g_wkb:11.3f} "
                  f"{g_sim/g_wkb:7.3f}")
    out["wkb_ratio_45_59km"] = {hl(hh): float(
        (wrms[iz(59), it(hh)] / wrms[iz(45), it(hh)]) / (wkb[iz(59)] / wkb[iz(45)]))
        for hh in H_FROM1 if hh <= 3.0}

    # ---------------- 2. Overturning onset ----------------
    hdr("2. ONSET  -  first flow reversal (u_hat >= U) and first N^2 < 0, by level")
    print(f"{'z[km]':>7} {'t_reversal[h]':>14} {'t_unstable[h]':>14} "
          f"{f'u_hat/U @{H_REF:g}h':>12}")
    onset = {}
    for zt in [0.75, 5, 10, 20, 30, 40, 45, 50, 53, 55, 57, 58.75, 59.75]:
        j = iz(zt)
        rev = np.where(umin[j] <= 0)[0]
        uns = np.where(ufrac[j] > 0)[0]
        tr = float(th[rev[0]]) if len(rev) else float("nan")
        tu = float(th[uns[0]]) if len(uns) else float("nan")
        r3 = (ubar[j, it(H_REF)] - umin[j, it(H_REF)]) / ubar[j, it(H_REF)]
        print(f"{zkm[j]:7.2f} {tr:14.2f} {tu:14.2f} {r3:12.2f}")
        onset[f"{zkm[j]:.2f}"] = {"t_reversal_h": tr, "t_unstable_h": tu}
    out["onset_by_level"] = onset

    # ---------------- 3. Momentum budget ----------------
    hdr("3. BUDGET  -  d<u>/dt vs wave drag -(1/rho) d<rho u'w'>/dz   [m/s per hour]")
    print("   (sponge K == 0 identically below 60 km; residual carries the lateral")
    print("    flux through the outflow x-boundary, which is not periodic in Stage 1)")
    sm = lambda a: np.apply_along_axis(
        lambda v: np.convolve(v, np.ones(5) / 5.0, mode="same"), 0, a)
    dudt = sm(np.gradient(ubar, t, axis=1) * 3600.0)
    drag = sm(-(np.gradient(mf, z, axis=0) / rho[:, None]) * 3600.0)
    budget = {}
    for hh in H_FROM2:
        i = it(hh)
        print(f"-- t={th[i]:.2f} h")
        print(f"   {'z[km]':>7} {'d<u>/dt':>10} {'wave drag':>10} {'residual':>10} {'ratio':>7}")
        for zt in [20, 30, 40, 45, 50, 53, 55, 57, 58.75]:
            j = iz(zt)
            a, b = dudt[j, i], drag[j, i]
            r = a / b if abs(b) > 1e-6 else float("nan")
            print(f"   {zkm[j]:7.2f} {a:10.3f} {b:10.3f} {a-b:10.3f} {r:7.2f}")
        # closure quality at the level of peak deceleration
        jp = int(np.argmin(np.where(clean, dudt[:, i], 0.0)))
        budget[hl(hh)] = {"z_peak_decel_km": float(zkm[jp]),
                            "dudt": float(dudt[jp, i]), "drag": float(drag[jp, i]),
                            "ratio": float(dudt[jp, i] / drag[jp, i])}
    out["budget_at_peak_deceleration"] = budget

    # ---------------- 4. Descent ----------------
    hdr("4. DESCENT  -  observed layer base vs two independent predictions")
    base = np.full(len(t), np.nan)
    for i in range(len(t)):
        idx = np.where(clean & ((ubar[:, i] - U_BAR) <= -5.0))[0]
        if len(idx):
            base[i] = zkm[idx[0]]
    dudz = np.gradient(ubar, z, axis=0)
    print(f"{'t[h]':>6} {'obs base':>9} {'z(|MF|=F_sat)':>14} {'obs dz/dt':>10} "
          f"{'-(du/dt)/(du/dz)':>17}")
    desc = {}
    for hh in H_HALF3:
        i = it(hh)
        ratio = np.array([abs(mf[j, i]) / F_sat(j, i) for j in range(len(z))])
        idx = np.where(clean & (ratio >= 1.0))[0]
        zsat = zkm[idx[0]] if len(idx) else float("nan")
        w = slice(max(0, i - 8), min(len(t), i + 9))
        obs = np.polyfit(th[w], base[w], 1)[0] if np.isfinite(base[w]).all() else np.nan
        j = iz(base[i]) if np.isfinite(base[i]) else 0
        kin = -dudt[j, i] / (dudz[j, i] * 1e3) if abs(dudz[j, i]) > 1e-9 else np.nan
        print(f"{th[i]:6.2f} {base[i]:9.2f} {zsat:14.2f} {obs:10.2f} {kin:17.2f}")
        desc[hl(hh)] = {"base_km": float(base[i]), "z_saturation_km": float(zsat),
                          "obs_dzdt_km_per_h": float(obs), "kinematic_km_per_h": float(kin)}
    out["descent"] = desc

    hdr("4b. EQUILIBRIUM breaking height for the flux actually supplied")
    Fs40 = rho * k_eff * U_BAR ** 3 / (2 * N)
    eq = {}
    for hh in H_FROM1:
        F = abs(mf[iz(5), it(hh)])
        j = np.where((Fs40 < F) & (zkm > 2))[0]
        zq = float(zkm[j[0]]) if len(j) else float("nan")
        print(f"  |MF(5 km)| at t={hh:g} h = {F*1e3:6.3f} mPa -> F_sat crossing at z = {zq:5.1f} km")
        eq[hl(hh)] = zq
    Fl = abs(mf_lin[iz(10)])
    j = np.where((Fs40 < Fl) & (zkm > 2))[0]
    eq["linear_steady"] = float(zkm[j[0]]) if len(j) else float("nan")
    print(f"  full steady LINEAR flux  = {Fl*1e3:6.3f} mPa -> F_sat crossing at z = "
          f"{eq['linear_steady']:5.1f} km   <-- asymptotic breaking height")
    out["equilibrium_breaking_height_km"] = eq

    # ---------------- 5. Deposition partition ----------------
    hdr("5. DEPOSITION  -  where the incoming flux is actually absorbed")
    bands = [(25, 40), (40, 45), (45, 53), (53, 56), (56, 60)]
    dep = {}
    for hh in H_FROM3:
        i = it(hh)
        Fin = abs(mf[iz(25), i])
        parts = [(abs(mf[iz(a), i]) - abs(mf[iz(b), i])) / Fin * 100 for a, b in bands]
        print(f"  t={th[i]:.2f} h  F_in(25 km)={Fin*1e3:6.3f} mPa | " +
              "  ".join(f"{a}-{b}: {p:5.1f}%" for (a, b), p in zip(bands, parts)))
        dep[hl(hh)] = {f"{a}-{b}km": float(p) for (a, b), p in zip(bands, parts)}
    out["deposition_percent"] = dep

    # ---------------- 6. Sponge and resolution ----------------
    hdr("6. ARTIFACT TESTS")
    print("6a. Sponge reference wind (src/run.jl:57 sets Q_bar = initial state, u = 40 m/s)")
    h_cols = [0.0, round(T_END_H / 2.0), round(T_END_H)]
    print(f"    {'z[km]':>7} {'K[1/s]':>10} {'tau[s]':>9} " +
          " ".join(f"{f'<u>@{h:g}h':>8}" for h in h_cols))
    for zt in [56.75, 58.75, 59.75, 60.25, 62, 65, 70, 80]:
        j = iz(zt)
        eta = max(0.0, (z[j] - Z_SPONGE) / (Z_TOP - Z_SPONGE))
        K = K0_SPONGE * np.sin(0.5 * np.pi * eta) ** 2
        tau = 1 / K if K > 0 else float("inf")
        print(f"    {zkm[j]:7.2f} {K:10.3e} {tau:9.1f} " +
              " ".join(f"{ubar[j, it(h)]:8.2f}" for h in h_cols))

    print("\n6b. Static instability per unit wave slope  S = w_rms N / (k U^2).")
    print("    A band whose ufrac/S is far above the column value has a non-wave source.")
    print(f"    {'t[h]':>5} {'band':>12} {'<ufrac>':>9} {'<S>':>8} {'ufrac/S':>9}")
    art = {}
    for hh in H_FROM4:
        i = it(hh)
        S = wrms[:, i] * N / (k_eff * ubar[:, i] ** 2)
        for lo, hi, lab in [(45, 57, "45-57 km"), (57.5, 60, "57.5-60 km")]:
            b = (zkm >= lo) & (zkm <= hi)
            r = ufrac[b, i].mean() / S[b].mean()
            print(f"    {th[i]:5.2f} {lab:>12} {ufrac[b, i].mean():9.3f} "
                  f"{S[b].mean():8.2f} {r:9.4f}")
            art[f"{hl(hh)}_{lab}"] = float(r)
    out["ufrac_per_slope"] = art

    print("\n6c. Resolution of the flux-carrying scale")
    print(f"    {'t[h]':>5} {'z[km]':>7} {'U':>6} {'lam_z[km]':>10} {'pts/lam_z':>10} "
          f"{'lam_x(MF=Fsat)[km]':>19} {'pts/lam_x':>10}")
    for hh in H_FROM3:
        i = it(hh)
        for zt in [50, 55]:
            j = iz(zt)
            U = ubar[j, i]
            lz = 2 * np.pi * U / N[j]
            lx = 2 * np.pi / (2 * abs(mf[j, i]) * N[j] / (rho[j] * U ** 3))
            print(f"    {th[i]:5.2f} {zkm[j]:7.2f} {U:6.2f} {lz/1e3:10.2f} {lz/DZ:10.1f} "
                  f"{lx/1e3:19.1f} {lx/DX:10.1f}")
    # Evaluated at the last sample, i.e. the most marginally resolved state.
    i_last = len(th) - 1
    out["min_pts_per_lambda_x"] = float(
        2 * np.pi / (2 * abs(mf[iz(55), i_last]) * N[iz(55)] /
                     (rho[iz(55)] * ubar[iz(55), i_last] ** 3)) / DX)
    out["series_end_h"] = T_END_H

    if args.json:
        with open(args.json, "w") as f:
            json.dump(out, f, indent=2)
        print(f"\nwrote {args.json}")


if __name__ == "__main__":
    main()
