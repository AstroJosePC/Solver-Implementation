#!/usr/bin/env python3
"""
Section 16 figures (docs/reflection_convergence_tests.md, Section 16): the old Lie splitting (L)
against the Strang (S) and alternating (A) paths on the 2 h S14V1a ridge, runs S16{L,S,A}{1..4}
at CFL 0.4, 0.2, 0.1, 0.05. Style: scripts/magnusp.mplstyle and
docs/splitting_order_design_philosophy.md ("Halving Cartography").

  extract       read the NetCDF snapshots once, write NPZ caches to output/analysis/s16/
  diff-fields   fig_s16_ls_difference_fields: 2^(l-1) (w'_L - w'_S), filtered w', four levels
  order-courant fig_s16_order_vs_courant: q_AS and q_LS against the x-sweep Courant number,
                with the 1-D Lax-Wendroff half-step test and the two-term fit of E_AS
  ground-2dz    fig_s16_ground_2dz_layer: the 2 dz part of w' (w' minus [1,2,1]/4 w'), L scheme,
                CFL 0.4 against CFL 0.05
  all           extract (if missing) and all three figures

    ~/.venvs/magnusp/bin/python scripts/plot_splitting_order.py all [--runs ~/MAGNUSP_runs]

Arrays: Julia writes (x, z) column-major, which netCDF4 returns as (z, x); every loaded field is
asserted to have shape (NZ, NX) = (190, 500). Times in s, heights in km, w in m/s unless suffixed.
"""
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from s16_posthoc import CASE, LEVEL_TAG, SCHEME_TAG, band_norm, fields, smooth, snaps  # noqa: E402

NZ, NX = 190, 500
DX_M, DZ_M = 4000.0, 500.0
X0_KM = 1000.0                 # ridge centre, x0 = LX/2 (sweep_arsia_mons_stage2.jl)
SIGMA_KM = 50.0                # ridge half-width, sigma=50
U_BAR = 40.0                   # background wind [m/s]
GAMMA = 1.31
CFL = {1: 0.4, 2: 0.2, 3: 0.1, 4: 0.05}
EVEN_IDX = (24, 48, 72, 97)    # snapshots common to all levels, even step counts
XWIN_KM = (X0_KM - 400.0, X0_KM + 400.0)

SCRIPT = Path(__file__).resolve()
STYLE = SCRIPT.parent / "magnusp.mplstyle"
FONT_DIR = Path.home() / ".claude/skills/canvas-design/canvas-fonts"

# Halving Cartography tokens (docs/splitting_order_design_philosophy.md)
INK, SLATE, GRID, MUTED, SAND = "#181615", "#5A5550", "#E6E1DA", "#8C857B", "#F3EFEA"
EMBER, COBALT, TEAL, OCHRE, CHARCOAL = "#C83E2B", "#1E4E79", "#2E777D", "#D97724", "#3B3834"


def main_checkout():
    common = subprocess.run(["git", "-C", str(SCRIPT.parent), "rev-parse", "--git-common-dir"],
                            capture_output=True, text=True, check=True).stdout.strip()
    return (SCRIPT.parent / common).resolve().parent


def git_commit():
    r = subprocess.run(["git", "-C", str(SCRIPT.parent), "rev-parse", "--short", "HEAD"],
                       capture_output=True, text=True)
    d = subprocess.run(["git", "-C", str(SCRIPT.parent), "status", "--porcelain", "--", ":/scripts", ":/src"],
                       capture_output=True, text=True)
    return r.stdout.strip() + ("-dirty" if d.stdout.strip() else "")


def cache_dir():
    return main_checkout() / "output" / "analysis" / "s16"


def load_zx(path, name):
    """One (z, x) field from a snapshot, asserted against the grid."""
    import netCDF4
    with netCDF4.Dataset(path) as d:
        a = d[name][:].data
        assert d[name].dimensions == ("z", "x"), d[name].dimensions
    assert a.shape == (NZ, NX), a.shape
    return a


def snap_attr(path, name):
    import netCDF4
    with netCDF4.Dataset(path) as d:
        return d.getncattr(name)


def run_id(s, l):
    return f"S16{s}{l}"


# ----------------------------------------------------------------------------------- extract

def lw1d_q(nus, sigma_cells, N=400, T=0.5, shape="dgauss"):
    """
    1-D two-step Richtmyer Lax-Wendroff for u_t + a u_x = 0 (a = 1, periodic, N cells):
    D(nu) = ||LW(dt/2)^2 - LW(dt)||_2 / ||u0||_2 after time T at dt = nu dx / a, and
    q = log2 D(nu_i) / D(nu_{i+1}). shape 'gauss' is s16_posthoc.lw1d (exp(-((x-0.3)/0.03)^2));
    'dgauss' is the ridge forcing shape dh/dx of exp(-x^2 / (2 sigma^2)), sigma in cells.
    """
    dx = 1.0 / N
    x = (np.arange(N) + 0.5) * dx
    if shape == "gauss":
        u0 = np.exp(-((x - 0.3) / 0.03) ** 2)
    else:
        s = sigma_cells * dx
        u0 = -(x - 0.5) / s ** 2 * np.exp(-((x - 0.5) ** 2) / (2 * s ** 2))

    def step(u, dt):
        up = np.roll(u, -1)
        h = 0.5 * (u + up) - 0.5 * dt / dx * (up - u)
        return u - dt / dx * (h - np.roll(h, 1))

    D = []
    for nu in nus:
        dt = nu * dx
        uh, uf = u0.copy(), u0.copy()
        for _ in range(int(round(T / dt))):
            uh = step(step(uh, dt / 2), dt / 2)
            uf = step(uf, dt)
        D.append(np.linalg.norm(uh - uf) / np.linalg.norm(u0))
    D = np.array(D)
    return D, np.log2(D[:-1] / D[1:])


def extract(root):
    out = cache_dir()
    out.mkdir(parents=True, exist_ok=True)
    S = {(s, l): snaps(root, s, l) for s in "LSA" for l in CFL}
    for k, v in S.items():
        assert len(v) == 99, (k, len(v))
    z_km = fields(S["L", 1][-1], S["L", 1][0])["z"]
    import netCDF4
    with netCDF4.Dataset(S["L", 1][0]) as d:
        x_km = d["x"][:].data / 1e3
    assert z_km.shape == (NZ,) and x_km.shape == (NX,)

    # mean dt per level from the final snapshot (t_end / steps)
    t_end = np.array([float(snap_attr(S["L", l][-1], "time_s")) for l in CFL])
    steps = np.array([int(snap_attr(S["L", l][-1], "step")) for l in CFL])
    dt_l = t_end / steps

    # Fig. 1: filtered w' difference L - S at t_end
    dwf, E_LS = [], []
    for l in CFL:
        fL, fS = fields(S["L", l][-1], S["L", l][0]), fields(S["S", l][-1], S["S", l][0])
        assert fL["wf"].shape == (NZ, NX) and abs(fL["t"] - fS["t"]) < 1e-3
        dwf.append(fL["wf"] - fS["wf"])
        E_LS.append(band_norm(fL["wf"] - fS["wf"], z_km) / band_norm(fS["wf"], z_km))
        if l == 1:
            wfS = fS["wf"]
    np.savez(out / "S16_lsdiff_xz.npz", dwf_lzx=np.array(dwf), wfS_zx=wfS, z_km=z_km, x_km=x_km,
             E_LS_l=np.array(E_LS), t_s=t_end, dt_s_l=dt_l, steps_l=steps)

    # Fig. 2: q at the even-step snapshots and at t_end, measured Courant numbers, fits
    def E_pair(a, b, idx):
        E, t, st = [], [], []
        for l in CFL:
            fa, fb = fields(S[a, l][idx], S[a, l][0]), fields(S[b, l][idx], S[b, l][0])
            assert abs(fa["t"] - fb["t"]) < 1e-3
            E.append(band_norm(fa["wf"] - fb["wf"], z_km) / band_norm(fb["wf"], z_km))
            t.append(fa["t"])
            st.append(int(snap_attr(S[a, l][idx], "step")))
        return np.array(E), np.array(t), np.array(st)

    E_AS_tl, E_LS_tl, t_k, steps_kl = [], [], [], []
    for idx in EVEN_IDX + (-1,):
        eas, t, st = E_pair("A", "S", idx)
        els, _, _ = E_pair("L", "S", idx)
        E_AS_tl.append(eas), E_LS_tl.append(els), t_k.append(t[0]), steps_kl.append(st)
    E_AS_tl, E_LS_tl = np.array(E_AS_tl), np.array(E_LS_tl)
    q_AS_tk = np.log2(E_AS_tl[:, :-1] / E_AS_tl[:, 1:])
    q_LS_tk = np.log2(E_LS_tl[:, :-1] / E_LS_tl[:, 1:])

    nu_ac, nu_u, cs_ground = [], [], []
    for i, l in enumerate(CFL):
        p = S["L", l][-1]
        rho, u, w, P = (load_zx(p, n) for n in ("rho", "u", "w", "P"))
        cs = np.sqrt(GAMMA * np.maximum(P, 0.0) / rho)
        nu_ac.append(np.max(np.abs(u) + cs) * dt_l[i] / DX_M)
        nu_u.append(np.max(np.abs(u)) * dt_l[i] / DX_M)
        cs_ground.append(float(cs[0].mean()))
    nu_ac, nu_u, cs_ground = np.array(nu_ac), np.array(nu_u), np.array(cs_ground)

    # two-term fit E_AS = a r + b r^2, r = dt/dt_1, relative least squares, at index 97
    r = dt_l / dt_l[0]
    k97 = EVEN_IDX.index(97)
    A = np.vstack([r, r ** 2]).T / E_AS_tl[k97][:, None]
    a, b = np.linalg.lstsq(A, np.ones(4), rcond=None)[0]

    nus_ref = np.array([0.4, 0.2, 0.1, 0.05, 0.025])
    _, q_ref = lw1d_q(nus_ref, None, shape="gauss")
    nus_1d = 0.4 / 2.0 ** np.arange(9)
    sigma_cells = SIGMA_KM * 1e3 / DX_M
    D_1d, q_1d = lw1d_q(nus_1d, sigma_cells)
    np.savez(out / "S16_order_courant.npz", q_AS_tk=q_AS_tk, q_LS_tk=q_LS_tk, E_AS_tl=E_AS_tl,
             E_LS_tl=E_LS_tl, t_k=np.array(t_k), steps_kl=np.array(steps_kl), nu_ac_l=nu_ac,
             nu_u_l=nu_u, dt_s_l=dt_l, fit_ab=np.array([a, b]), nu_1d=nus_1d, D_1d=D_1d, q_1d=q_1d,
             q_ref_posthoc=q_ref, ktheta=1.0 / sigma_cells)

    # Fig. 3: 2 dz part of w', L scheme, levels 1 and 4
    w2, w2max, Ld = [], [], []
    j10 = int(np.argmin(abs(z_km - 10.25)))
    for i, l in enumerate(CFL):
        f = fields(S["L", l][-1], S["L", l][0])
        part = f["wp"] - f["wf"]
        w2max.append(float(np.abs(part[j10]).max()))
        Ld.append(DZ_M ** 2 / (2.0 * cs_ground[i] * dt_l[i]) / 1e3)
        if l in (1, 4):
            w2.append(part)
    np.savez(out / "S16_ground2dz_xz.npz", w2dz_lzx=np.array(w2), levels=np.array([1, 4]),
             Ld_km_l=np.array(Ld), w2dz_max_10km_l=np.array(w2max), z_km=z_km, x_km=x_km,
             cs_ground_l=cs_ground, dt_s_l=dt_l)

    print(f"t_end {t_end}  steps {steps}  dt {dt_l}")
    print(f"E_LS(t_end) {['%.3e' % e for e in E_LS]}")
    print(f"q_LS even idx {np.round(q_LS_tk[:-1], 2).tolist()}  t_end {np.round(q_LS_tk[-1], 2)}")
    print(f"q_AS even idx {np.round(q_AS_tk[:-1], 2).tolist()}  t_end {np.round(q_AS_tk[-1], 2)}")
    print(f"t_k {np.round(t_k, 1)}  steps at idx 97 {steps_kl[k97]}")
    print(f"nu_x acoustic {nu_ac}  advective {nu_u}  c_ground {cs_ground}")
    print(f"fit a {a:.3e} b {b:.3e}  rel. residuals {(a * r + b * r ** 2) / E_AS_tl[k97] - 1}")
    print(f"lw1d posthoc shape q {np.round(q_ref, 2)}  ridge shape q {np.round(q_1d, 2)}")
    print(f"2dz max at 10.25 km {['%.1e' % v for v in w2max]}  L_d km {np.round(Ld, 2)}")


# ----------------------------------------------------------------------------------- figures

def setup_style(font_dir):
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import font_manager
    import matplotlib.pyplot as plt
    for f in ("Outfit-Regular.ttf", "Outfit-Bold.ttf"):
        p = Path(font_dir) / f
        if p.exists():
            font_manager.fontManager.addfont(str(p))
    plt.style.use(str(STYLE))
    return plt


def save(fig, outdir, stem, sidecar):
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / f"{stem}.pdf")
    fig.savefig(outdir / f"{stem}.png", dpi=300)
    sidecar = dict(figure=stem, script=str(SCRIPT.relative_to(SCRIPT.parent.parent)),
                   git_commit=git_commit(), style="scripts/magnusp.mplstyle",
                   philosophy="docs/splitting_order_design_philosophy.md",
                   section="docs/reflection_convergence_tests.md, Section 16",
                   case=f"{CASE}{{scheme tag}}{{cfl tag}}/run_U40.0_h0.03",
                   forcing="face path (bcface=1), linear, h = 0.03 km, sigma = 50 km, U = 40 m/s, ramp 2 h",
                   dx_m=DX_M, dz_m=DZ_M, **sidecar)
    with open(outdir / f"{stem}.json", "w") as f:
        json.dump(sidecar, f, indent=2)
    print(f"wrote {outdir / stem}.{{pdf,png,json}}")


def panel_label(ax, s):
    ax.text(0.03, 0.97, s, transform=ax.transAxes, ha="left", va="top", color=INK,
            fontsize=8.5, fontweight="bold")


def fig_diff_fields(plt, outdir):
    from matplotlib.colors import TwoSlopeNorm
    c = np.load(cache_dir() / "S16_lsdiff_xz.npz")
    z, x = c["z_km"], c["x_km"]
    mx = (x >= XWIN_KM[0]) & (x <= XWIN_KM[1])
    mz = z <= 35.0
    scaled = np.array([2.0 ** (l - 1) * c["dwf_lzx"][l - 1] for l in CFL]) * 1e6   # um/s
    lim = np.percentile(np.abs(scaled[:, mz][:, :, mx]), 99.5)
    norm = TwoSlopeNorm(0.0, -lim, lim)
    E = c["E_LS_l"]
    q = np.log2(E[:-1] / E[1:])
    fig, axes = plt.subplots(1, 4, figsize=(7.0, 2.55), sharey=True)
    ctx = c["wfS_zx"][mz][:, mx]
    for i, ax in enumerate(axes):
        l = i + 1
        pm = ax.pcolormesh(x[mx] - X0_KM, z[mz], scaled[i][mz][:, mx], cmap="RdBu_r", norm=norm,
                           shading="nearest", rasterized=True)
        ax.contour(x[mx] - X0_KM, z[mz], ctx, levels=[-0.5 * np.abs(ctx).max(), 0.5 * np.abs(ctx).max()],
                   colors=MUTED, linewidths=0.35, linestyles="solid")
        for zb in (12.0, 30.0):
            ax.axhline(zb, color=SLATE, lw=0.5, ls=(0, (3, 2)))
        ax.set_xlim(-400, 400)
        ax.set_ylim(0, 35)
        ax.set_xticks([-300, 0, 300])
        ax.set_xlabel(r"$x - x_0$ [km]")
        panel_label(ax, f"({'abcd'[i]})")
        ax.text(0.97, 0.97, rf"CFL ${CFL[l]:g}$", transform=ax.transAxes, ha="right", va="top",
                color=SLATE, fontsize=8)
        txt = rf"$10^{{4}}\,E_{{\mathrm{{LS}}}} = {E[i] * 1e4:.2f}$"
        if i > 0:
            txt += rf"$,\ \ q = {q[i - 1]:.2f}$"
        ax.text(1.0, 1.02, txt, transform=ax.transAxes, ha="right", va="bottom", color=EMBER,
                fontsize=7.5)
    axes[0].set_ylabel(r"altitude $z$ [km]")
    cb = fig.colorbar(pm, ax=axes, shrink=0.92, pad=0.015, aspect=28)
    cb.set_label(r"$2^{\,l-1}\,(w'_{\mathrm{L}} - w'_{\mathrm{S}})$, filtered [$\mu$m s$^{-1}$]")
    cb.outline.set_linewidth(0.4)
    save(fig, outdir, "fig_s16_ls_difference_fields", dict(
        runs=[run_id(s, l) for s in "LS" for l in CFL], t_s=float(c["t_s"][0]),
        dt_s_l=c["dt_s_l"].tolist(), steps_l=c["steps_l"].tolist(), cfl_l=list(CFL.values()),
        field="w' (w minus its x-mean) filtered in z by [1,2,1]/4, L minus S, times 2^(l-1)",
        window=dict(x_km=XWIN_KM, z_km=[0, 35], band_km=[12, 30]), color_limit_um_s=float(lim),
        context_contour="filtered w'_S at CFL 0.4, +-0.5 of its window maximum",
        E_LS_l=E.tolist(), q_LS=q.tolist(), inputs=[str(cache_dir() / "S16_lsdiff_xz.npz")]))
    plt.close(fig)


def fig_order_courant(plt, outdir):
    c = np.load(cache_dir() / "S16_order_courant.npz")
    nu = c["nu_ac_l"]
    nu_h = np.sqrt(nu[:-1] * nu[1:])                     # each halving at the geometric mean
    qa, ql = c["q_AS_tk"][:-1], c["q_LS_tk"][:-1]        # even-step snapshots
    qa_end = c["q_AS_tk"][-1]
    k97 = EVEN_IDX.index(97)
    a, b = c["fit_ab"]
    fig, ax = plt.subplots(figsize=(3.5, 2.75))
    for qq, lab, dy, va in ((1.0, "first order", -0.04, "top"), (2.0, "second order", 0.03, "bottom")):
        ax.axhline(qq, color=GRID, lw=0.8, zorder=0)
        ax.text(2.1e-3, qq + dy, lab, color=SLATE, fontsize=7, ha="left", va=va)
    # 1-D Richtmyer LW(dt/2)^2 against LW(dt), ridge-shaped profile
    nu1 = np.sqrt(c["nu_1d"][:-1] * c["nu_1d"][1:])
    ax.plot(nu1, c["q_1d"], ls=":", color=CHARCOAL, lw=1.1, marker="o", ms=3, mfc="white", mec=CHARCOAL)
    ax.text(0.105, 1.40, "1-D Lax-Wendroff,\ntwo half steps\nagainst one step",
            color=CHARCOAL, fontsize=7, ha="left", va="top", linespacing=1.15)
    # two-term model of E_AS, as a continuous q(nu)
    rr = np.geomspace(1.0, 1 / 8, 200)
    E_of = lambda r: a * r + b * r ** 2
    q_mod = np.log2(E_of(rr) / E_of(rr / 2))
    ax.plot(nu[0] * rr / np.sqrt(2), q_mod, ls="--", color=TEAL, lw=1.0)
    ax.text(0.0036, 1.25, r"fit $a\,\Delta t + b\,\Delta t^{2}$", color=TEAL, fontsize=7,
            ha="center", va="bottom")
    # 2-D points: median over snapshot times, range bar over times
    for qset, col, mk, lab in ((qa, TEAL, "o", r"A$-$S"), (ql, EMBER, "s", r"L$-$S")):
        ax.vlines(nu_h, qset.min(axis=0), qset.max(axis=0), color=col, lw=1.0)
        ax.plot(nu_h, qset[k97], ls="-", color=col, lw=1.0, marker=mk, ms=4, mfc=col, mec=col)
        ax.text(nu_h[0] * 1.12, qset[k97][0], lab, color=col, fontsize=8, ha="left", va="center")
    ax.plot(nu_h, qa_end, ls="none", marker="o", ms=4, mfc="white", mec=TEAL, mew=0.8)
    ax.text(nu_h[-1] * 1.25, qa_end[-1], r"$t_{\mathrm{end}}$, odd step count", color=SLATE,
            fontsize=7, ha="left", va="center")
    ax.set_xscale("log")
    ax.set_xlim(2.0e-3, 0.6)
    ax.set_ylim(0.6, 2.35)
    ax.set_xlabel(r"x-sweep Courant number $\nu_x = \max(|u|+c)\,\Delta t/\Delta x$")
    ax.set_ylabel(r"order per halving $q$")
    save(fig, outdir, "fig_s16_order_vs_courant", dict(
        runs=[run_id(s, l) for s in "LSA" for l in CFL], dt_s_l=c["dt_s_l"].tolist(),
        nu_acoustic_l=nu.tolist(), nu_advective_l=c["nu_u_l"].tolist(),
        snapshot_times_s=c["t_k"].tolist(), snapshot_steps_kl=c["steps_kl"].tolist(),
        q_AS_tk=c["q_AS_tk"].tolist(), q_LS_tk=c["q_LS_tk"].tolist(),
        markers="filled: t = %.0f s (even steps); bars: range over t = %s s; open: t_end (odd steps at levels 1-3)"
                % (c["t_k"][k97], np.round(c["t_k"][:-1]).tolist()),
        fit="E_AS = a r + b r^2, r = dt/dt_1, a = %.3e, b = %.3e (t = %.0f s)" % (a, b, c["t_k"][k97]),
        lw1d=dict(profile="dh/dx of a Gaussian, sigma = %.1f cells (ridge sigma/dx)" % (1 / c["ktheta"]),
                  N=400, T=0.5, nu=c["nu_1d"].tolist(), q=c["q_1d"].tolist(),
                  q_posthoc_shape=c["q_ref_posthoc"].tolist()),
        metric="E_ab: L2 norm of filtered w'_a - w'_b over 12-30 km, all x, relative to w'_b",
        inputs=[str(cache_dir() / "S16_order_courant.npz")]))
    plt.close(fig)


def fig_ground_2dz(plt, outdir):
    from matplotlib.colors import TwoSlopeNorm
    c = np.load(cache_dir() / "S16_ground2dz_xz.npz")
    z, x = c["z_km"], c["x_km"]
    mx = (x >= XWIN_KM[0]) & (x <= XWIN_KM[1])
    mz = z <= 20.0
    f = c["w2dz_lzx"] * 1e3                                    # mm/s
    lim = np.percentile(np.abs(f[1][mz][:, mx]), 99.5)
    norm = TwoSlopeNorm(0.0, -lim, lim)
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.6), sharey=True)
    for i, (ax, l) in enumerate(zip(axes, c["levels"])):
        pm = ax.pcolormesh(x[mx] - X0_KM, z[mz], f[i][mz][:, mx], cmap="RdBu_r", norm=norm,
                           shading="nearest", rasterized=True)
        Ld = c["Ld_km_l"][l - 1]
        ax.axhline(Ld, color=OCHRE, lw=0.9, ls=(0, (4, 2)))
        ax.text(395, Ld + 0.3, rf"$L_d = \Delta z^2/(2c\,\Delta t) = {Ld:.2f}$ km", color=OCHRE,
                fontsize=7.5, ha="right", va="bottom")
        ax.axhline(10.25, color=SLATE, lw=0.5, ls=(0, (1, 2)))
        ax.text(-395, 10.7, rf"$\max_x|w'_{{2\Delta z}}|(10.25\ \mathrm{{km}}) = "
                rf"{c['w2dz_max_10km_l'][l - 1] * 1e3:.2g}$ mm s$^{{-1}}$",
                color=SLATE, fontsize=7.5, ha="left", va="bottom")
        ax.set_xlim(-400, 400)
        ax.set_ylim(0, 20)
        ax.set_yticks([0, 5, 10, 15, 20])
        ax.set_xlabel(r"$x - x_0$ [km]")
        panel_label(ax, f"({'ab'[i]})")
        ax.text(0.97, 0.97, rf"L scheme, CFL ${CFL[int(l)]:g}$", transform=ax.transAxes,
                ha="right", va="top", color=SLATE, fontsize=8)
    axes[0].set_ylabel(r"altitude $z$ [km]")
    cb = fig.colorbar(pm, ax=axes, shrink=0.92, pad=0.015, aspect=24, extend="both")
    cb.set_label(r"$w' - [1,2,1]\,w'/4$ [mm s$^{-1}$]")
    cb.outline.set_linewidth(0.4)
    save(fig, outdir, "fig_s16_ground_2dz_layer", dict(
        runs=["S16L1", "S16L4"], t_s=7200.0, dt_s_l=c["dt_s_l"].tolist(),
        field="2 dz part of w': w' minus w' filtered in z by [1,2,1]/4, L scheme",
        window=dict(x_km=XWIN_KM, z_km=[0, 20]), color_limit_mm_s=float(lim),
        color_limit_from="99.5th percentile of |field| at CFL 0.05 in the window; one scale for both panels",
        Ld_km_l=c["Ld_km_l"].tolist(), cs_ground_m_s_l=c["cs_ground_l"].tolist(),
        w2dz_max_10km_m_s_l=c["w2dz_max_10km_l"].tolist(),
        inputs=[str(cache_dir() / "S16_ground2dz_xz.npz")]))
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["extract", "diff-fields", "order-courant", "ground-2dz", "all"])
    ap.add_argument("--runs", default=os.path.expanduser("~/MAGNUSP_runs"))
    ap.add_argument("--outdir", default=None, help="default: docs/figures in the main checkout")
    ap.add_argument("--font-dir", default=str(FONT_DIR))
    args = ap.parse_args()
    outdir = Path(args.outdir) if args.outdir else main_checkout() / "docs" / "figures"
    if args.cmd == "extract" or (args.cmd == "all" and not (cache_dir() / "S16_ground2dz_xz.npz").exists()):
        extract(args.runs)
    if args.cmd == "extract":
        return
    plt = setup_style(args.font_dir)
    for name, fn in (("diff-fields", fig_diff_fields), ("order-courant", fig_order_courant),
                     ("ground-2dz", fig_ground_2dz)):
        if args.cmd in (name, "all"):
            fn(plt, outdir)


if __name__ == "__main__":
    main()
