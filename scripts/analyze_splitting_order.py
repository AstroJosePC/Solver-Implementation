#!/usr/bin/env python3
"""
Section 16 readings (docs/reflection_convergence_tests.md): temporal self-convergence of
the dimensional and operator splitting on the face-path Gaussian ridge.

Each run is passed as --run SCHEME:LEVEL=<run dir>, with SCHEME one of L (Lie, :xz),
S (:strang + operator_split :strang), A (:alternate + :strang) or X (:strang, Lie
operators), and LEVEL 1 to 4 for CFL 0.4, 0.2, 0.1 and 0.05. Every run must hold a
snapshot at t_end (driver `final=1`); the last snapshot of each run is read and the times
are checked to agree to 1e-6 s.

Metric of record: w' = w minus its x-mean, filtered in z with [1, 2, 1]/4 (which removes
the near-ground 2 dz alternation of Section 14), over 12 to 30 km, above the reach of the
dt-dependent ground layer (decay length dz^2/(2 c dt), 5.9 km at CFL 0.05), and all x.
  E_ab(l)     = ||w'_a,l - w'_b,l|| / ||w'_b,l||        (schemes a, b at one level)
  q_ab        = log2 E_ab(l)/E_ab(l+1)                 (order of their difference)
  D_s(l, l+1) = ||w'_l - w'_{l+1}|| / ||w'_4||        (scheme s, successive levels)
  p_s         = log2 D_s(1,2)/D_s(2,3), log2 D_s(2,3)/D_s(3,4)
Also reported: the same on the unfiltered w' and on 2 to 8 km, and the momentum flux
MF(z) = rho_bar <u' w'> at 10.25 and 20.25 km with its relative change per level.

    python3 scripts/analyze_splitting_order.py --run L:1=<dir> ... --out s16.json
"""
import argparse
import glob
import json
import os

import netCDF4
import numpy as np

BAND = (12.0, 30.0)
LOW = (2.0, 8.0)
MF_Z = (10.25, 20.25)


def load(run):
    snaps = sorted(glob.glob(os.path.join(run, "snapshot_*.nc")))
    with netCDF4.Dataset(snaps[0]) as d:
        rho0 = d["rho"][:].data
    with netCDF4.Dataset(snaps[-1]) as d:
        z = d["z"][:].data / 1e3
        u, w = d["u"][:].data, d["w"][:].data
        t = float(d.getncattr("time_s"))
    wp = w - w.mean(axis=1, keepdims=True)
    up = u - u.mean(axis=1, keepdims=True)
    wf = wp.copy()
    wf[1:-1] = 0.25 * wp[:-2] + 0.5 * wp[1:-1] + 0.25 * wp[2:]
    rb = rho0.mean(axis=1)
    mf = rb * (up * wp).mean(axis=1)
    return dict(t=t, z=z, wf=wf, wr=wp, mf=mf, snap=os.path.basename(snaps[-1]))


def norm(a, z, band):
    m = (z >= band[0]) & (z < band[1])
    return float(np.sqrt((a[m] ** 2).sum()))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", action="append", required=True, help="SCHEME:LEVEL=<run dir>")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    runs = {}
    for spec in a.run:
        key, path = spec.split("=", 1)
        s, l = key.split(":")
        runs[(s, int(l))] = load(path)
    t0 = next(iter(runs.values()))["t"]
    for k, r in runs.items():
        assert abs(r["t"] - t0) < 1e-6, f"{k} ends at {r['t']} s, not {t0} s"
    z = next(iter(runs.values()))["z"]
    out = {"t_end_s": t0, "band_km": BAND, "low_band_km": LOW, "per_scheme": {}, "between": {}}

    for s in sorted({k[0] for k in runs}):
        lv = sorted(l for (ss, l) in runs if ss == s)
        rec = {}
        for field, band, tag in (("wf", BAND, "filtered"), ("wr", BAND, "raw"), ("wf", LOW, "filtered_2_8km")):
            ref = norm(runs[(s, lv[-1])][field], z, band)
            D = [norm(runs[(s, l)][field] - runs[(s, l + 1)][field], z, band) / ref
                 for l in lv[:-1]]
            p = [float(np.log2(D[i] / D[i + 1])) for i in range(len(D) - 1)]
            rec[tag] = {"D": D, "p": p}
        for zz in MF_Z:
            j = int(np.argmin(abs(z - zz)))
            mf = [float(runs[(s, l)]["mf"][j]) for l in lv]
            rec[f"MF_{zz}"] = {"MF_Pa": mf,
                               "change_pct": [100 * (mf[i + 1] / mf[i] - 1) for i in range(len(mf) - 1)]}
        out["per_scheme"][s] = rec

    for sa, sb in (("L", "S"), ("L", "A"), ("A", "S"), ("X", "S"), ("L", "X")):
        common = sorted(l for (s, l) in runs if s == sa and (sb, l) in runs)
        if not common:
            continue
        E = {}
        for l in common:
            ra, rb = runs[(sa, l)], runs[(sb, l)]
            E[l] = {"E_filtered": norm(ra["wf"] - rb["wf"], z, BAND) / norm(rb["wf"], z, BAND),
                    "E_filtered_2_8km": norm(ra["wf"] - rb["wf"], z, LOW) / norm(rb["wf"], z, LOW)}
            for zz in MF_Z:
                j = int(np.argmin(abs(z - zz)))
                E[l][f"MF{zz}_diff_pct"] = 100 * (ra["mf"][j] / rb["mf"][j] - 1)
        q = [float(np.log2(E[l]["E_filtered"] / E[l + 1]["E_filtered"]))
             for l in common[:-1] if l + 1 in E]
        out["between"][f"{sa}-{sb}"] = {"per_level": E, "q": q}

    with open(a.out, "w") as f:
        json.dump(out, f, indent=1)
    for s, rec in out["per_scheme"].items():
        print(f"{s}: D = {['%.3e' % d for d in rec['filtered']['D']]}  p = {['%.2f' % p for p in rec['filtered']['p']]}"
              f"  | raw p = {['%.2f' % p for p in rec['raw']['p']]}  | 2-8 km p = {['%.2f' % p for p in rec['filtered_2_8km']['p']]}")
        print(f"   MF10 change %/level = {['%.4f' % c for c in rec['MF_10.25']['change_pct']]}"
              f"  MF20 = {['%.4f' % c for c in rec['MF_20.25']['change_pct']]}")
    for pair, rec in out["between"].items():
        E = rec["per_level"]
        print(f"{pair}: E = {['%.2e' % E[l]['E_filtered'] for l in E]}  q = {['%.2f' % v for v in rec['q']]}"
              f"  MF10 diff % = {['%.4f' % E[l]['MF10.25_diff_pct'] for l in E]}"
              f"  MF20 diff % = {['%.4f' % E[l]['MF20.25_diff_pct'] for l in E]}")


if __name__ == "__main__":
    main()
