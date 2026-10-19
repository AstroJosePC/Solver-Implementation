#!/usr/bin/env python3
"""
Section 16 post-hoc checks (docs/reflection_convergence_tests.md, Section 16, Results).
Not readings of record: they explain the 16-b result and the raw-field self-convergence.

  pairs    E_ab and q_ab (filtered w', 12 to 30 km) at snapshots common to all four levels,
           which fall on even step counts (outint 100 * 2^(l-1)), so the alternating path
           has completed its last xz/zx pair. Times agree between schemes to ~4e-5 s
           (adaptive dt), which is far below the differences measured.
  mfprof   change per level of the L-scheme flux rho_bar <u' w'>, raw and with u' and w'
           filtered by [1, 2, 1]/4, and the size of the 2 dz part of w', by height.
  lw1d     1-D Richtmyer Lax-Wendroff, linear advection: two dt/2 steps against one dt
           step at a fixed grid, the x-sweep difference between Strang and alternating.

    python3 scripts/s16_posthoc.py <runs root> pairs 24 48 72 97
    python3 scripts/s16_posthoc.py <runs root> mfprof
    python3 scripts/s16_posthoc.py - lw1d
"""
import glob
import sys

import numpy as np

CASE = "arsia_mons_linear2h_LX2000_sig50_K0_0.005_LZ95_pbarP_ramp2h_bcface"
SCHEME_TAG = dict(L="", S="_strang_opS", A="_alt_opS", X="_strang")
LEVEL_TAG = {1: "", 2: "_cfl0.2", 3: "_cfl0.1", 4: "_cfl0.05"}


def snaps(root, s, l):
    return sorted(glob.glob(f"{root}/{CASE}{SCHEME_TAG[s]}{LEVEL_TAG[l]}/run_U40.0_h0.03/snapshot_*.nc"))


def smooth(a):
    f = a.copy()
    f[1:-1] = 0.25 * a[:-2] + 0.5 * a[1:-1] + 0.25 * a[2:]
    return f


def fields(path, rho0_path):
    import netCDF4
    with netCDF4.Dataset(rho0_path) as d:
        rb = d["rho"][:].data.mean(axis=1)
    with netCDF4.Dataset(path) as d:
        z = d["z"][:].data / 1e3
        u, w = d["u"][:].data, d["w"][:].data
        t = float(d.getncattr("time_s"))
    up = u - u.mean(axis=1, keepdims=True)
    wp = w - w.mean(axis=1, keepdims=True)
    return dict(z=z, t=t, rb=rb, up=up, wp=wp, wf=smooth(wp), uf=smooth(up))


def band_norm(a, z):
    m = (z >= 12.0) & (z < 30.0)
    return float(np.sqrt((a[m] ** 2).sum()))


def pairs(root, idx):
    for k in idx:
        print(f"--- snapshot index {k}")
        for a, b in (("L", "S"), ("A", "S"), ("L", "A")):
            E, t = [], []
            for l in (1, 2, 3, 4):
                sa, sb = snaps(root, a, l), snaps(root, b, l)
                fa, fb = fields(sa[k], sa[0]), fields(sb[k], sb[0])
                assert abs(fa["t"] - fb["t"]) < 1e-3
                E.append(band_norm(fa["wf"] - fb["wf"], fa["z"]) / band_norm(fb["wf"], fb["z"]))
                t.append(fa["t"])
            q = [np.log2(E[i] / E[i + 1]) for i in range(3)]
            print(f"{a}-{b} t = {t[0]:.1f} s  E = {['%.2e' % e for e in E]}  q = {['%.2f' % v for v in q]}")


def mfprof(root):
    R = {}
    for l in (1, 2, 3, 4):
        s = snaps(root, "L", l)
        f = fields(s[-1], s[0])
        R[l] = (f["z"], f["rb"] * (f["up"] * f["wp"]).mean(axis=1),
                f["rb"] * (f["uf"] * f["wf"]).mean(axis=1), np.abs(f["wp"] - f["wf"]).max(axis=1))
    z = R[1][0]
    print("z [km] | raw MF change % per level | filtered MF change % per level | max |2dz part of w'| [m/s], L1..L4")
    for zz in (0.75, 1.75, 3.75, 5.75, 7.75, 10.25, 13.75, 20.25, 29.75):
        j = int(np.argmin(abs(z - zz)))
        cr = [100 * (R[l + 1][1][j] / R[l][1][j] - 1) for l in (1, 2, 3)]
        cf = [100 * (R[l + 1][2][j] / R[l][2][j] - 1) for l in (1, 2, 3)]
        print(f"{z[j]:6.2f} | " + " ".join(f"{c:7.2f}" for c in cr) + " | "
              + " ".join(f"{c:6.2f}" for c in cf) + " | " + " ".join(f"{R[l][3][j]:.1e}" for l in (1, 2, 3, 4)))


def lw1d():
    N, a, T = 400, 1.0, 0.5
    dx = 1.0 / N
    x = (np.arange(N) + 0.5) * dx
    u0 = np.exp(-((x - 0.3) / 0.03) ** 2)

    def step(u, dt):
        up = np.roll(u, -1)
        h = 0.5 * (u + up) - 0.5 * a * dt / dx * (up - u)
        return u - a * dt / dx * (h - np.roll(h, 1))

    def run(dt, half):
        u = u0.copy()
        for _ in range(int(round(T / dt))):
            u = step(step(u, dt / 2), dt / 2) if half else step(u, dt)
        return u

    prev = None
    for nu in (0.4, 0.2, 0.1, 0.05, 0.025):
        dt = nu * dx / a
        d = np.linalg.norm(run(dt, True) - run(dt, False)) / np.linalg.norm(u0)
        print(f"nu = {nu:5.3f}  ||LW(dt/2)^2 - LW(dt)|| = {d:.3e}" + (f"  q = {np.log2(prev / d):.2f}" if prev else ""))
        prev = d


if __name__ == "__main__":
    root, what = sys.argv[1], sys.argv[2]
    if what == "pairs":
        pairs(root, [int(k) for k in sys.argv[3:]] or [24, 48, 72, 97])
    elif what == "mfprof":
        mfprof(root)
    elif what == "lw1d":
        lw1d()
    else:
        sys.exit(__doc__)
