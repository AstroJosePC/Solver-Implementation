"""Two summary figures of the Arsia linear control corner, for presentation.

1. benchmark: how the model top is tested. Momentum-flux profiles against the
   steady linear value, the approach to steady state for each damping layer, and
   the measured returned share against the linear Rayleigh-layer prediction.
2. showcase: the L1 run (36 h, h = 30 m, sigma = 50 km, U = 40 m/s, MCD column,
   60-90 km layer). Density-scaled w at five times and the time-height record of
   the x-mean momentum flux.

Context and talking points: docs/reflection_convergence_tests.md, "Summary figures".

Usage:
    python scripts/meeting_figures.py            # both figures
    python scripts/meeting_figures.py --gif      # also an animation of the showcase run
"""
import argparse
import glob
import json
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm
from netCDF4 import Dataset

OUT = "output/meeting_figures"
A = "output/arsia_mons_linear15h_LX2000_sig50{}_analysis/run_U40.0_h0.03"
L1_RUN = "output/arsia_mons_linear36h_LX2000_sig50/run_U40.0_h0.03"
L1_ANA = "output/arsia_mons_linear36h_LX2000_sig50_analysis/run_U40.0_h0.03"
RUNS = [  # label, analysis dir, layer description, layer base and top [km]
    ("S0", A.format(""), "60-90 km, K0 0.05", 60, 90),
    ("S2", A.format("_K0_0.0158"), "60-90 km, K0 0.0158", 60, 90),
    ("S1", A.format("_K0_0.005"), "60-90 km, K0 0.005", 60, 90),
    ("T1", A.format("_K0_0.005_LZ120"), "60-120 km, K0 0.005", 60, 120),
    ("T2", A.format("_LZ100_spb70"), "70-100 km, K0 0.05", 70, 100),
]
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]   # as compare_sponge_reflection.py
MARKERS = ["o", "s", "^", "D", "v"]
INK, MUTED, GRID = "0.15", "0.45", "0.9"
BAND = (10.0, 55.0)          # evaluation band [km]
RIDGE_X = 100                # vertical exaggeration of the 30 m ridge in the sketch


def mf_linear():
    """Steady linear momentum flux [Pa] used as the normaliser in every run."""
    with open(os.path.join(L1_ANA, "linear_flux_constancy.json")) as f:
        return json.load(f)["linear_steady_MF_uPa"] * 1e-6


def load_series(ana):
    d = np.load(os.path.join(ana, "deepdive_series.npz"))
    return d["times"] / 3600.0, d["z_km"], d["mf_zt"]


def style(ax):
    ax.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(colors=INK)


def benchmark():
    lin = mf_linear()
    fig, ax = plt.subplots(1, 3, figsize=(16, 5.4), gridspec_kw={"width_ratios": [1, 1.25, 1.1]})

    # (a) flux profile in the final hour of each run
    a = ax[0]
    a.axhspan(60, 100, color="0.94", zorder=0)
    a.axhspan(*BAND, color="#2a78d6", alpha=0.06, zorder=0)
    a.text(0.03, 95, "damping layers\n(base 60 km; T2 70 km)", fontsize=9, color=MUTED, va="top")
    a.text(0.03, 11.2, "evaluation band 10-55 km", fontsize=9, color=MUTED)
    a.axvline(1.0, color=MUTED, ls="--", lw=1.2)
    a.text(1.02, 94, "steady linear\ntheory", fontsize=9, color=MUTED, va="top")
    for (lab, ana, desc, *_), c in zip(RUNS, COLORS):
        t, z, mf = load_series(ana)
        sel = t >= t[-1] - 1.0
        a.plot(mf[:, sel].mean(1) / lin, z, color=c, lw=1.8, label=f"{lab} (14-15 h)")
    t, z, mf = load_series(L1_ANA)
    sel = t >= t[-1] - 1.0
    a.plot(mf[:, sel].mean(1) / lin, z, color=INK, lw=2.6, label="L1 (35-36 h)")
    a.set_xlim(0, 1.4)
    a.set_ylim(0, 100)
    a.set_xlabel("momentum flux / steady linear value")
    a.set_ylabel("altitude [km]")
    a.set_title("(a) Flux is carried up and absorbed at the top", loc="left", fontsize=11)
    a.legend(fontsize=8.5, loc="center left", bbox_to_anchor=(0.0, 0.4), frameon=False)
    style(a)

    # (b) band-mean flux against time
    b = ax[1]
    for (lab, ana, desc, *_), c, mk in zip(RUNS, COLORS, MARKERS):
        t, z, mf = load_series(ana)
        zb = (z >= BAND[0]) & (z <= BAND[1])
        lev = mf[zb].mean(0) / lin
        b.plot(t, lev, color=c, lw=1.8, label=f"{lab}  {desc}")
        b.plot(t[-1], lev[-1], marker=mk, color=c, ms=7, mec="white", mew=1.2)
    t, z, mf = load_series(L1_ANA)
    zb = (z >= BAND[0]) & (z <= BAND[1])
    b.plot(t, mf[zb].mean(0) / lin, color=INK, lw=2.2, label="L1  S0 continued to 36 h")
    b.axvline(15, color=MUTED, ls=":", lw=1)
    b.text(15.4, 1.22, "15 h: end of the layer sweep", fontsize=9, color=MUTED)
    b.set_xlim(0, 36)
    b.set_ylim(0, 1.3)
    b.set_xlabel("time [h]")
    b.set_ylabel("band-mean momentum flux / steady linear value")
    b.set_title("(b) Approach to steady state, 10-55 km mean", loc="left", fontsize=11)
    b.legend(fontsize=8.5, loc="lower right", frameon=False)
    style(b)

    # (c) returned share against the linear Rayleigh-layer prediction
    cax = ax[2]
    with open("output/reflection_tests/compare_sponge_reflection.json") as f:
        cmp_ = json.load(f)
    lim = (1e-3, 1.0)
    cax.plot(lim, lim, color=MUTED, lw=1)
    cax.text(0.35, 0.5, "measured = theory", rotation=45, fontsize=9, color=MUTED,
             ha="center", va="center")
    for (lab, *_), c, mk in zip(RUNS, COLORS, MARKERS):
        s = cmp_[lab]
        cax.plot(s["R_theory"], s["R_split"], marker=mk, color=c, ms=10, ls="none",
                 mec="white", mew=1.2)
        note = "  (outside linear range)" if lab == "T2" else ""
        cax.annotate(lab + note, (s["R_theory"], s["R_split"]), xytext=(7, 4),
                     textcoords="offset points", fontsize=9, color=INK)
    with open("output/reflection_tests/cavity_closure.json") as f:
        cav = json.load(f)
    for r, lab in zip(cav, ("C0", "C1")):
        hh = [h for h in r["hours"] if 12 <= h["t_h"] <= 17.5]
        meas = np.mean([h["R_split_band"] for h in hh])
        cax.plot(r["layer_theory"]["R_band"], meas, marker="*", color=COLORS[0 if lab == "C0" else 3],
                 ms=15, ls="none", mec=INK, mew=0.8)
        cax.annotate(f"{lab}, single mode", (r["layer_theory"]["R_band"], meas), xytext=(7, -12),
                     textcoords="offset points", fontsize=9, color=INK)
    cax.set_xscale("log")
    cax.set_yscale("log")
    cax.set_xlim(*lim)
    cax.set_ylim(*lim)
    cax.set_xlabel("predicted returned share, linear Rayleigh layer")
    cax.set_ylabel("measured returned share |B/A|$^2$")
    cax.set_title("(c) Returned share does not follow the layer theory", loc="left", fontsize=11)
    style(cax)
    cax.grid(True, which="major", color=GRID, lw=0.8)

    fig.suptitle("Top-boundary and convergence benchmark: linear Arsia control, U = 40 m/s, "
                 "h = 30 m, sigma = 50 km, MCD MY31 column", fontsize=12, color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    path = os.path.join(OUT, "fig1_benchmark_sponge_convergence.png")
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print("wrote", path)


def snapshot_at(run, t_h, files, times):
    i = int(np.argmin(np.abs(times - t_h * 3600.0)))
    with Dataset(files[i]) as d:
        w = d["w"][:].data
        return float(d.time_s), d["x"][:].data, d["z"][:].data, w - w.mean(1, keepdims=True)


def l1_index():
    files = sorted(glob.glob(os.path.join(L1_RUN, "snapshot_*.nc")))
    t, _, _ = load_series(L1_ANA)       # the deep-dive series has one entry per snapshot
    return files, t * 3600.0


def scaled_w(files):
    with Dataset(files[0]) as d:
        rho0 = d["rho"][:, 0].data.astype(float)
    return np.sqrt(rho0 / rho0[0])[:, None]


def showcase():
    files, times = l1_index()
    if len(files) != len(times):
        times = np.array([Dataset(f).time_s for f in files])
    scale = scaled_w(files)
    t_show = [1, 3, 6, 12, 36]
    snaps = [snapshot_at(L1_RUN, th, files, times) for th in t_show]
    vmax = 0.03
    fig = plt.figure(figsize=(16, 8.6))
    gs = fig.add_gridspec(2, len(t_show), height_ratios=[1.15, 1], hspace=0.38, wspace=0.08)
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    for k, (t, x, z, w) in enumerate(snaps):
        a = fig.add_subplot(gs[0, k])
        xr = (x - 0.5 * (x[0] + x[-1])) / 1e3
        im = a.pcolormesh(xr, z / 1e3, w * scale, cmap="RdBu_r", norm=norm, shading="auto",
                          rasterized=True)
        a.axhspan(60, 90, color="white", alpha=0.35, lw=0)
        a.axhline(60, color=INK, ls="--", lw=0.9)
        if k == 0:
            a.text(-480, 62, "damping layer", fontsize=9, color=INK)
        a.fill_between(xr, 0, 30e-3 * RIDGE_X * np.exp(-(xr / 50.0) ** 2), color=INK, lw=0)
        a.set_xlim(-500, 500)
        a.set_ylim(0, 90)
        a.set_title(f"t = {t / 3600:.0f} h", fontsize=11, color=INK)
        a.set_xlabel("x - x_ridge [km]")
        if k == 0:
            a.set_ylabel("altitude [km]")
        else:
            a.set_yticklabels([])
        a.tick_params(colors=INK)
    cb = fig.colorbar(im, ax=fig.axes[:len(t_show)], fraction=0.015, pad=0.01, extend="both")
    cb.set_label("w' (rho/rho_s)$^{1/2}$ [m/s], x-mean removed")

    # time-height record of the x-mean momentum flux
    b = fig.add_subplot(gs[1, :])
    t, z, mf = load_series(L1_ANA)
    lin = mf_linear()
    pm = b.pcolormesh(t, z, mf / lin, cmap="Blues", vmin=0, vmax=1.1, shading="auto", rasterized=True)
    b.axhline(60, color=INK, ls="--", lw=0.9)
    for th in t_show:
        b.axvline(th, color="white", lw=0.8, alpha=0.8)
    b.set_ylim(0, 90)
    b.set_xlim(0, 36)
    b.set_xlabel("time [h]")
    b.set_ylabel("altitude [km]")
    b.set_title("x-mean momentum flux / steady linear value; the front rises and the layer "
                "absorbs it above 60 km", loc="left", fontsize=11, color=INK)
    cb2 = fig.colorbar(pm, ax=b, fraction=0.015, pad=0.01)
    cb2.set_label("MF / MF_linear")
    fig.suptitle("L1: linear mountain wave on the Arsia MCD column (h = 30 m, sigma = 50 km, "
                 "U = 40 m/s; periodic 2000 km x 90 km domain, +-500 km shown; dx = 4 km, "
                 "dz = 500 m; ridge height drawn x100)",
                 fontsize=12, color=INK, x=0.01, ha="left")
    path = os.path.join(OUT, "fig2_showcase_L1_wave.png")
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote", path)
    return files, times, scale


def animate(files, times, scale, step_min=15):
    from matplotlib.animation import FuncAnimation, PillowWriter
    targets = np.arange(0, times[-1] + 1, step_min * 60.0)
    idx = sorted({int(np.argmin(np.abs(times - tt))) for tt in targets})
    with Dataset(files[0]) as d:
        x, z = d["x"][:].data, d["z"][:].data
    xr = (x - 0.5 * (x[0] + x[-1])) / 1e3
    fig, a = plt.subplots(figsize=(7, 5))
    norm = TwoSlopeNorm(vmin=-0.03, vcenter=0, vmax=0.03)
    im = a.pcolormesh(xr, z / 1e3, np.zeros((len(z), len(x))), cmap="RdBu_r", norm=norm,
                      shading="auto")
    a.axhline(60, color=INK, ls="--", lw=0.9)
    a.text(-480, 62, "damping layer", fontsize=9, color=INK)
    a.set_xlim(-500, 500)
    a.set_xlabel("x - x_ridge [km]")
    a.set_ylabel("altitude [km]")
    fig.colorbar(im, ax=a, extend="both").set_label("w' (rho/rho_s)$^{1/2}$ [m/s], x-mean removed")
    ttl = a.set_title("")

    def frame(i):
        with Dataset(files[i]) as d:
            w = d["w"][:].data
            im.set_array(((w - w.mean(1, keepdims=True)) * scale).ravel())
            ttl.set_text(f"L1, t = {float(d.time_s) / 3600:5.2f} h")
        return im, ttl

    anim = FuncAnimation(fig, frame, frames=idx, blit=False)
    path = os.path.join(OUT, "fig2_showcase_L1_wave.gif")
    anim.save(path, writer=PillowWriter(fps=8), dpi=110)
    plt.close(fig)
    print("wrote", path, f"({len(idx)} frames)")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gif", action="store_true", help="also write an animation of the L1 run")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    benchmark()
    files, times, scale = showcase()
    if a.gif:
        animate(files, times, scale)


if __name__ == "__main__":
    main()
