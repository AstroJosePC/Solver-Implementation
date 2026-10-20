#!/usr/bin/env python3
"""
Section 2 figures (docs/diffusion_solver_performance.md, Section 2): the eight :fft_tridiag twins
S02F1..S02F8 (set d2, diffsolver=fft bcface=1) against their archived :lu off members. Style:
scripts/magnusp.mplstyle and docs/diffusion_solver_design_philosophy.md ("Floor and Fork
Cartography"). Reads saved JSON, TSV and run logs only; no snapshot is opened.

  extract       parse the cost histories and run logs once, write
                output/analysis/s02/figcache/cost_series.npz
  convergence   fig_s02_convergence: rel difference against t for the eight pairs, bands to t_c
  observables   fig_s02_observables: MF(30.25) window-mean changes against the band and the hourly
                range; onset and end-time differences of twin and off member
  cost          fig_s02_cost: wall time per step against simulated time, and per-pair ratios
  all           extract (if missing) and the three figures

    ~/.venvs/magnusp/bin/python scripts/plot_diffusion_solver.py all [--outdir docs/figures]

Inputs: output/analysis/s02/readings_d2.json (series t_h is the off member's snapshot time; pairs
are matched by step), output/d2/logs/{d2.cost.tsv,d2.start,S02F*.log},
output/reflection_tests/logs/{S15,S14A3}.{cost.tsv,start} and the off-member logs. Times in h
unless suffixed, cost in s of wall clock per step.
"""
import argparse
import json
import re
import subprocess
from datetime import datetime
from pathlib import Path

import numpy as np

SCRIPT = Path(__file__).resolve()
REPO = SCRIPT.parent.parent
STYLE = SCRIPT.parent / "magnusp.mplstyle"
FONT_DIR = Path.home() / ".claude/skills/canvas-design/canvas-fonts"
OUT = REPO / "output"
S02 = OUT / "analysis" / "s02"
CACHE = S02 / "figcache" / "cost_series.npz"
READINGS = S02 / "readings_d2.json"
D2LOG = OUT / "d2" / "logs"
RTLOG = OUT / "reflection_tests" / "logs"

# Floor and Fork tokens (docs/diffusion_solver_design_philosophy.md)
GROUND, SAND, INK, SLATE, GRID, MUTED = "#FFFFFF", "#F3EFEA", "#181615", "#5A5550", "#E6E1DA", "#8C857B"
EMBER, COBALT, TEAL, OCHRE, CHARCOAL = "#C83E2B", "#1E4E79", "#2E777D", "#D97724", "#3B3834"

# twin, off member, dz [m], U [m/s], h [km], cost history of the off member
PAIRS = [
    ("S02F1", "S14S2", 250, 40, 0.1, "S14A3"),
    ("S02F2", "PR1-U20-h0.1", 500, 20, 0.1, "S15"),
    ("S02F3", "PR1-U57-h0.1", 500, 57, 0.1, "S15"),
    ("S02F4", "PR1-U20-h1.0", 500, 20, 1.0, "S15"),
    ("S02F5", "S15F1", 500, 40, 1.0, "S15"),
    ("S02F6", "PR1-U57-h1.0", 500, 57, 1.0, "S15"),
    ("S02F7", "PR1-U20-h3.0", 500, 20, 3.0, "S15"),
    ("S02F8", "PR1-U40-h3.0", 500, 40, 3.0, "S15"),
]
CFL = 0.4
RAMP_H = 2.0
BAND_THERMO = 1e-11
BAND_W = {"D-1": 1e-7, "D-2": 1e-8}
BAND_MF = 1e-7
SNAP_S = 132.0                 # one scored snapshot (stride 2), D-2 row 1
CAVEAT = ("The twin/off-member ratios mix three changes and are descriptive of the production gain "
          "only: the diffusion solver, the pack and unpack refactor 0f6ff62 (which also speeds up :lu; "
          "the off members ran before it, s15 on 671bf32 and S14S2 in set a3), and concurrency 8 on an "
          "otherwise idle machine for d2 against 12 (s15) or 11 (a3) at launch for the off members.")


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


def readings():
    with open(READINGS) as f:
        return json.load(f)


# ----------------------------------------------------------------------------------- extract

def read_cost_tsv(path):
    """Cost monitor history: epoch seconds [n] and {column: values [n]} (NA as NaN)."""
    lines = path.read_text().splitlines()
    head = lines[0].split("\t")
    rows = [l.split("\t") for l in lines[1:] if l.strip()]
    t = np.array([datetime.strptime(r[0], "%Y-%m-%d %H:%M:%S").timestamp() for r in rows])
    cols = {h: np.array([float(r[i]) if r[i] not in ("NA", "") else np.nan for r in rows])
            for i, h in enumerate(head[1:], 1)}
    return t, cols


def read_step_map(path):
    """(step, simulated time [s]) from the 'Step N | t= ... s' lines of a run log, with (0, 0)."""
    st, ts = [0], [0.0]
    for l in path.read_text().splitlines():
        m = re.match(r"Step\s+(\d+) \| t=\s*([\d.]+) s", l)
        if m:
            st.append(int(m.group(1)))
            ts.append(float(m.group(2)))
    return np.array(st), np.array(ts)


def read_start(path):
    """Start and end epochs of each run of a launcher set from its .start stamp."""
    out = {}
    for l in path.read_text().splitlines():
        m = re.match(r"(\S+) exit=\d+ start=(\d+) end=(\d+)", l)
        if m:
            out[m.group(1)] = (int(m.group(2)), int(m.group(3)))
    return out


def rolling_cost(t, step, wmin_s, nmin_steps, wmax_s=4 * 3600.0):
    """
    Wall time per step from a cost history: at each row, the inverse slope of a least-squares fit
    of step against wall time over a centred window of at least wmin_s seconds and nmin_steps
    steps (widened up to wmax_s). The logged step is quantized to 500, so short windows alias.
    Returns centre step [n] and s per step [n].
    """
    cs, sp = [], []
    for k in range(len(t)):
        half = wmin_s / 2
        while True:
            m = (t >= t[k] - half) & (t <= t[k] + half)
            span_t, span_s = np.ptp(t[m]), np.ptp(step[m])
            if (span_t >= 0.9 * wmin_s and span_s >= nmin_steps) or 2 * half >= wmax_s:
                break
            half *= 1.25
        if m.sum() < 4 or span_s <= 0 or span_t < 0.9 * wmin_s:
            continue
        p = np.polyfit(t[m] - t[k], step[m], 1)
        cs.append(p[1])
        sp.append(1.0 / p[0])
    return np.array(cs), np.array(sp)


def extract():
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    hist = {"d2": (D2LOG / "d2.cost.tsv", D2LOG / "d2.start"),
            "S15": (RTLOG / "S15.cost.tsv", RTLOG / "S15.start"),
            "S14A3": (RTLOG / "S14A3.cost.tsv", RTLOG / "S14A3.start")}
    loaded = {k: (read_cost_tsv(a), read_start(b)) for k, (a, b) in hist.items()}
    arrays, summary = {}, {}
    for twin, off, dz, U, h, offset in PAIRS:
        for rid, setname, logdir, wmin, nmin in ((twin, "d2", D2LOG, 300.0, 8000),
                                                 (off, offset, RTLOG, 900.0, 4000)):
            (t, cols), start = loaded[setname]
            step = cols[f"{rid}_step"]
            ok = ~np.isnan(step)
            t_r, s_r = t[ok], step[ok]
            iend = int(np.argmax(s_r >= np.nanmax(s_r)))       # first row at the last logged step
            t_r, s_r = t_r[:iend + 1], s_r[:iend + 1]
            st, ts = read_step_map(logdir / f"{rid}.log")
            assert st[-1] == int(np.nanmax(step)), (rid, st[-1], np.nanmax(step))
            cstep, sps = rolling_cost(t_r, s_r, wmin, nmin)
            keep = (cstep > 0) & (cstep < st[-1])
            cstep, sps = cstep[keep], sps[keep]
            tsim_h = np.interp(cstep, st, ts) / 3600.0
            # number of runs of the set active at each centre (wall time interpolated from step)
            twall = np.interp(cstep, s_r, t_r)
            active = np.array([sum(a <= tw <= b for a, b in start.values()) for tw in twall])
            arrays[f"{rid}__tsim_h"] = tsim_h
            arrays[f"{rid}__s_per_step"] = sps
            arrays[f"{rid}__active"] = active
            arrays[f"{rid}__tend_h"] = np.array(ts[-1] / 3600.0)
            done = re.search(r"Simulation complete: (\d+) steps, t=([\d.]+) s", (logdir / f"{rid}.log").read_text())
            arrays[f"{rid}__complete_steps"] = np.array(int(done.group(1)))
            arrays[f"{rid}__complete_t_s"] = np.array(float(done.group(2)))
            summary[rid] = dict(set=setname, rows=int(len(t_r)), window_s=wmin, min_steps=nmin,
                                median_s_per_step=float(np.median(sps)),
                                active_min_max=[int(active.min()), int(active.max())])
            print(f"{rid:>13}: {len(sps)} points, median {np.median(sps):.4f} s/step, "
                  f"range {sps.min():.4f} to {sps.max():.4f}, active {active.min()} to {active.max()}")
        fl = re.search(r"(\d+) factorizations", (D2LOG / f"{twin}.log").read_text())
        arrays[f"{twin}__factorizations"] = np.array(int(fl.group(1)))
    np.savez_compressed(CACHE, **arrays)
    with open(CACHE.with_suffix(".json"), "w") as f:
        json.dump(dict(inputs={k: [str(a), str(b)] for k, (a, b) in hist.items()},
                       method="rolling least-squares slope of logged step against wall time; "
                              "simulated time from the run log step lines", runs=summary), f, indent=2)
    print(f"wrote {CACHE}")


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


def tint(color, frac):
    """Mix a token with GROUND: frac 1 is the token, 0 is the ground."""
    from matplotlib.colors import to_rgb
    c, g = np.array(to_rgb(color)), np.array(to_rgb(GROUND))
    return tuple(frac * c + (1 - frac) * g)


def save(fig, outdir, stem, sidecar):
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / f"{stem}.pdf")
    fig.savefig(outdir / f"{stem}.png", dpi=300)
    sidecar = dict(figure=stem, script="scripts/plot_diffusion_solver.py", git_commit=git_commit(),
                   style="scripts/magnusp.mplstyle",
                   philosophy="docs/diffusion_solver_design_philosophy.md",
                   section="docs/diffusion_solver_performance.md, Section 2 (Results)",
                   set="d2 (twins) against archived :lu records",
                   solver="twins :fft_tridiag (diffsolver=fft), off members :lu (global sparse LU, "
                          "refactorized when dt drifts by more than 1 %)",
                   forcing="face path (bcface=1), pbar=1, ramp tanh(t/2 h), K0 0.005, LZ 95 km, "
                           "LX 4000 km, sigma 100 km, dx 4000 m",
                   cfl=CFL, dz_m={p[0]: p[2] for p in PAIRS},
                   pairs={p[0]: p[1] for p in PAIRS}, **sidecar)
    with open(outdir / f"{stem}.json", "w") as f:
        json.dump(sidecar, f, indent=2)
    print(f"wrote {outdir / stem}.{{pdf,png,json}}")


def thermo_max(s):
    return np.max([s[k] for k in ("rel_rho", "rel_u", "rel_P", "rel_T")], axis=0)


def fig_convergence(plt, outdir):
    from matplotlib.ticker import FixedLocator, LogFormatterMathtext, NullLocator
    R = readings()
    fig, axes = plt.subplots(2, 4, figsize=(7.0, 3.95), sharex=True, sharey=True)
    ylo, yhi = 1e-18, 1e5
    wt = tint(COBALT, 0.42)
    we = tint(EMBER, 0.42)
    side = {}
    for i, (ax, (twin, off, dz, U, h, _)) in enumerate(zip(axes.flat, PAIRS)):
        p = R[twin]
        s = p["series"]
        t = np.array(s["t_h"])
        th, w = thermo_max(s), np.array(s["rel_w"])
        k0 = 1                                              # the t = 0 row is exactly zero
        t, th, w = t[k0:], th[k0:], w[k0:]
        tc = p["t_c_h"]
        bw = BAND_W[p["set"]]
        ax.axvspan(0, RAMP_H, color=SAND, lw=0, zorder=0)
        tb = tc if tc is not None else t[-1]
        ax.hlines([BAND_THERMO, bw], 0, tb, color=MUTED, lw=0.6, zorder=1)
        if tc is None:
            ax.plot(t, w, color=wt, lw=0.6, zorder=2)
            ax.plot(t, th, color=COBALT, lw=1.0, zorder=3)
        else:
            j = int(np.searchsorted(t, tc, side="right"))
            ax.plot(t[:j], w[:j], color=wt, lw=0.6, zorder=2)
            ax.plot(t[j - 1:], w[j - 1:], color=we, lw=0.6, zorder=2)
            ax.plot(t[:j], th[:j], color=COBALT, lw=1.0, zorder=3)
            ax.plot(t[j - 1:], th[j - 1:], color=EMBER, lw=1.0, zorder=3)
            ax.axvline(tc, color=OCHRE, lw=0.7, zorder=1.5)
            ax.text(tc + 0.45, 10 ** -17.55, rf"$t_c = {tc:.2f}$ h", color=OCHRE, fontsize=7,
                    ha="left", va="bottom")
            g = p.get("growth_fit") or {}
            if g.get("efold_h"):
                ax.text(23.6, 10 ** -17.55, rf"$\tau_e = {g['efold_h']:.2f}$ h", color=EMBER,
                        fontsize=7, ha="right", va="bottom")
        # every 1.5 h a small open marker on r_w, so the two series differ in grayscale too
        for tm in np.arange(0.75, t[-1], 1.5):
            k = int(np.argmin(abs(t - tm)))
            ax.plot(t[k], w[k], ls="none", marker="o", ms=2.1, mfc=GROUND, mew=0.55,
                    mec=(we if (tc is not None and t[k] > tc) else wt), zorder=4)
        letter = "abcdefgh"[i]
        extra = {"S02F1": r", $\Delta z$ 250 m", "S02F5": " (U40 h1.0)"}.get(twin, "")
        ax.annotate(f"({letter})", (0, 1), xycoords="axes fraction", xytext=(0, 4),
                    textcoords="offset points", ha="left", va="bottom", fontsize=7.5, color=INK,
                    fontweight="bold")
        ax.annotate(f"{twin} / {off}{extra}", (0, 1), xycoords="axes fraction", xytext=(13, 4),
                    textcoords="offset points", ha="left", va="bottom", fontsize=7.5, color=INK)
        side[twin] = dict(off=off, t_c_h=tc, band_thermo=BAND_THERMO, band_w=bw,
                          max_thermo_to_tc=float(th[t <= (tc if tc else 1e9)].max()),
                          max_w_to_tc=float(w[t <= (tc if tc else 1e9)].max()),
                          efold_h=(p.get("growth_fit") or {}).get("efold_h"),
                          t_last_shared_h=p["t_last_shared_h"], n_shared=p["n_shared"])
    ax = axes.flat[0]
    # key in the empty right half of (a) (S02F1 ends at 9 h)
    kx, ky, dy = 0.43, 0.92, 0.085
    items = [(COBALT, 1.0, None, r"$\max(r_\rho, r_u, r_P, r_T)$"),
             (wt, 0.6, "o", r"$r_w$"),
             (EMBER, 1.0, None, r"after $t_c$"),
             (MUTED, 0.6, None, r"band, to $t_c$")]
    for n, (c, lw, mk, lab) in enumerate(items):
        y = ky - n * dy
        ax.plot([kx, kx + 0.07], [y, y], transform=ax.transAxes, color=c, lw=lw, clip_on=False)
        if mk:
            ax.plot([kx + 0.035], [y], transform=ax.transAxes, marker=mk, ms=2.1, mfc=GROUND,
                    mec=c, mew=0.55, clip_on=False)
        ax.text(kx + 0.10, y, lab, transform=ax.transAxes, color=SLATE, fontsize=6.8, ha="left",
                va="center")
    y = ky - len(items) * dy
    ax.add_patch(plt.Rectangle((kx, y - 0.025), 0.07, 0.05, transform=ax.transAxes, color=SAND, lw=0))
    ax.text(kx + 0.10, y, r"ramp, $0$ to $2$ h", transform=ax.transAxes, color=SLATE,
            fontsize=6.8, ha="left", va="center")
    for a in axes.flat:
        a.set_yscale("log")
        a.set_ylim(ylo, yhi)
        a.set_xlim(0, 24)
        a.set_xticks([0, 6, 12, 18, 24])
        a.yaxis.set_major_locator(FixedLocator([1e-16, 1e-12, 1e-8, 1e-4, 1e0, 1e4]))
        a.yaxis.set_major_formatter(LogFormatterMathtext())
        a.yaxis.set_minor_locator(NullLocator())
    fig.supxlabel(r"simulated time $t$ [h]", fontsize=8.5, color=INK)
    fig.supylabel(r"relative difference $r_f$ [$-$]", fontsize=8.5, color=INK)
    caption = (
        "Relative difference r_f = max|f_twin - f_off| / max|f_off| over all cells at each shared "
        "snapshot (pairs matched by step; t is the off member's snapshot time) between the "
        ":fft_tridiag twins S02F1 to S02F8 (set d2) and their archived :lu off members. Thick: "
        "largest of r_rho, r_u, r_P, r_T; thin with open markers: r_w. Cobalt up to the off member's "
        "first overturning t_c (zoverturn, convective onset in 0 to 5 km at 1.0 km), ember after it. "
        "Hairlines: the registered bands, drawn over the stretch they judge (D-1 to the end of S02F1, "
        "S02F2 and to t_c of S02F3; D-2 to t_c): 1e-11 for rho, u, P, T and 1e-7 (D-1) or 1e-8 (D-2) "
        "for w. tau_e: e-folding time of ln rel_max fitted from t_c to the first time rel_max exceeds "
        "1e-3 (or the last shared snapshot). Face path (bcface=1), pbar=1, CFL 0.4, dz 500 m except "
        "S02F1 (250 m, 9 h). Shaded: forcing ramp. S02F5 ends at 15.35 h (twin) against 23.81 h (S15F1); "
        "its last shared snapshots compare the two runs at equal step, not equal time.")
    save(fig, outdir, "fig_s02_convergence", dict(
        inputs=[str(READINGS)], runs=[p[0] for p in PAIRS] + [p[1] for p in PAIRS],
        windows="whole record; bands judged at t <= t_c", pairs_detail=side,
        y_limits=[ylo, yhi], caption=caption))
    plt.close(fig)


def rel_range(run):
    lo, hi = run["hourly_range"]
    return (hi - lo) / run["MF_30.25_uPa"]


def fig_observables(plt, outdir):
    from matplotlib.ticker import FixedLocator, LogFormatterMathtext, NullLocator
    R = readings()
    MF = R["MF"]
    C = R["cost"]
    cc = np.load(CACHE)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.0, 2.95), sharey=True,
                                 gridspec_kw=dict(width_ratios=[1.0, 1.0]))
    rows = {p[0]: i for i, p in enumerate(PAIRS)}
    ylab = [f"{p[0]} / {p[1]}" for p in PAIRS]
    # (a) window-mean MF(30.25) changes
    wins = [("mf_F1_6_9", "S02F1", "o"), ("mf_F2_6_9", "S02F2", "o"), ("mf_F2_18_24", "S02F2", "D"),
            ("mf_F3_6_9", "S02F3", "o"), ("mf_F3_18_24", "S02F3", "D"),
            ("pre_F4", "S02F4", "s"), ("pre_F5", "S02F5", "s"), ("pre_F6", "S02F6", "s")]
    side_a = {}
    for key, twin, mk in wins:
        m = MF[key]
        off = dict((p[0], p[1]) for p in PAIRS)[twin]
        tc = R[twin]["t_c_h"]
        after = tc is not None and m["window_h"][0] >= tc
        dy = {"o": -0.17, "D": 0.17, "s": 0.0}[mk] if twin in ("S02F2", "S02F3") else 0.0
        y = rows[twin] + dy
        chg = abs(m["changes"]["MF_30.25"])
        rr = rel_range(m["runs"][off])
        col = EMBER if after else COBALT
        a1.plot([chg, rr], [y, y], ls=(0, (1, 1.6)), color=CHARCOAL, lw=0.6, zorder=1)
        a1.plot(rr, y, ls="none", marker=mk, ms=3.6, mfc=GROUND, mec=CHARCOAL, mew=0.7, zorder=3)
        a1.plot(chg, y, ls="none", marker=mk, ms=4.0, mfc=col, mec=col, zorder=4)
        side_a[key] = dict(twin=twin, off=off, window_h=m["window_h"], rel_change_MF30=m["changes"]["MF_30.25"],
                           relative_hourly_range_off=rr, after_t_c=after,
                           MF30_uPa_off=m["runs"][off]["MF_30.25_uPa"], MF30_uPa_twin=m["runs"][twin]["MF_30.25_uPa"])
    for twin in ("S02F7", "S02F8"):
        a1.text(2e-13, rows[twin], r"no flux quoted ($h$ = 3 km)", color=SLATE, fontsize=6.8,
                ha="left", va="center")
    a1.vlines(BAND_MF, -0.95, 2.45, color=MUTED, lw=0.6, zorder=0)
    a1.text(BAND_MF * 1.35, -0.88, r"band $10^{-7}$ (D-1)", color=SLATE, fontsize=6.8, ha="left",
            va="center")
    # marker key for windows, in the empty right part of the h = 3 km rows
    for n, (mk, lab) in enumerate((("o", "window 6 to 9 h"), ("D", "window 18 to 24 h"),
                                   ("s", "3 h before collapse"))):
        yy = 6.08 + 0.5 * n
        a1.plot(3e-5, yy, ls="none", marker=mk, ms=3.6, mfc=SLATE, mec=SLATE)
        a1.text(7e-5, yy, lab, color=SLATE, fontsize=6.8, ha="left", va="center")
    # direct labels on the first row
    a1.text(1.6e-12, -0.42, "twin against off member", color=COBALT, fontsize=6.8, ha="left", va="center")
    a1.text(0.5, -0.42, "hourly range / mean", color=CHARCOAL, fontsize=6.8, ha="right", va="center")
    a1.set_xscale("log")
    a1.set_xlim(1e-13, 1e1)
    a1.xaxis.set_major_locator(FixedLocator([1e-12, 1e-9, 1e-6, 1e-3, 1e0]))
    a1.xaxis.set_major_formatter(LogFormatterMathtext())
    a1.xaxis.set_minor_locator(NullLocator())
    a1.set_xlabel(r"relative change of $\langle\mathrm{MF}(30.25)\rangle$ [$-$]")
    a1.set_yticks(range(len(PAIRS)))
    a1.set_yticklabels(ylab, fontsize=7)
    a1.set_ylim(len(PAIRS) - 0.35, -1.15)
    a1.tick_params(axis="y", length=0)
    a1.text(0.0, 1.03, "(a)", transform=a1.transAxes, ha="left", va="bottom", fontsize=8.5, color=INK, fontweight="bold")

    # (b) event-time differences
    ZERO = 1e-14
    side_b = {}
    for twin, off, *_ in PAIRS:
        p = R[twin]
        y = rows[twin]
        on_o, on_t = p["onsets_off"], p["onsets_twin"]
        tc_key = p["t_c_event"]
        ev = []
        for k, a in on_o.items():
            b = on_t.get(k)
            if a and b and a.get("t_h") is not None and b.get("t_h") is not None:
                ev.append((k, (b["t_h"] - a["t_h"]) * 3600.0))
        c = C[twin]
        off_end_s = float(cc[f"{off}__complete_t_s"])
        assert abs(float(cc[f"{twin}__complete_t_s"]) - c["t_end_s"]) < 1e-3
        dend = c["t_end_s"] - off_end_s
        side_b[twin] = dict(onset_dt_s=dict(ev), t_c_event=tc_key, end_twin_s=c["t_end_s"],
                            end_off_s=off_end_s, end_dt_s=dend)
        xs = lambda v: max(abs(v), 0.0) if v != 0 else ZERO
        later = [xs(v) for k, v in ev if k != tc_key]
        if later:
            a2.plot(later, [y] * len(later), ls="none", marker="o", ms=3.6, mfc=GROUND, mec=EMBER,
                    mew=0.75, zorder=3)
        for k, v in ev:
            if k == tc_key:
                a2.plot(xs(v), y, ls="none", marker="^", ms=4.4, mfc=OCHRE, mec=OCHRE, zorder=4)
        a2.plot(xs(dend), y, ls="none", marker="s", ms=3.8, mfc=EMBER if dend else COBALT,
                mec=EMBER if dend else COBALT, zorder=4)
    a2.axvline(SNAP_S, color=CHARCOAL, lw=0.7, ls=(0, (1, 1.6)), zorder=0)
    a2.text(SNAP_S * 0.7, -0.75, "one scored\nsnapshot, 132 s", color=CHARCOAL, fontsize=6.8,
            ha="right", va="center", linespacing=1.0)
    a2.axvline(3.2e-14, color=GRID, lw=0.8, zorder=0)
    a2.set_xscale("log")
    a2.set_xlim(4e-15, 3e5)
    a2.xaxis.set_major_locator(FixedLocator([ZERO, 1e-12, 1e-8, 1e-4, 1e0, 1e4]))
    fmt = LogFormatterMathtext()
    a2.xaxis.set_major_formatter(lambda v, pos: "0" if abs(v - ZERO) < 1e-20 else fmt(v, pos))
    a2.xaxis.set_minor_locator(NullLocator())
    a2.set_xlabel(r"time difference, twin minus off member, $|\Delta t|$ [s]")
    a2.tick_params(axis="y", length=0)
    a2.text(0.0, 1.03, "(b)", transform=a2.transAxes, ha="left", va="bottom", fontsize=8.5, color=INK, fontweight="bold")
    kx, ky = 0.07, 1.03
    keys = [("^", OCHRE, OCHRE, r"$t_c$ onset"), ("o", GROUND, EMBER, "later onsets"),
            ("s", EMBER, EMBER, "end of run")]
    xpos = [0.10, 0.33, 0.62]
    for (mk, fc, ec, lab), xk in zip(keys, xpos):
        a2.plot([xk], [ky + 0.035], transform=a2.transAxes, ls="none", marker=mk, ms=3.8, mfc=fc,
                mec=ec, mew=0.75, clip_on=False)
        a2.text(xk + 0.03, ky + 0.035, lab, transform=a2.transAxes, color=SLATE, fontsize=6.8,
                ha="left", va="center")
    caption = (
        "(a) Relative change of the window mean of MF(30.25), twin against off member (filled; "
        "cobalt for windows before the first overturning t_c, ember after it), and the relative "
        "hourly range (max - min)/mean of the off member's hourly MF(30.25) in the same window (open, "
        "charcoal), the sampling scatter against which the change can be read. Windows: 6 to 9 h "
        "(circles), 18 to 24 h (diamonds; S02F3 after its t_c at 10.69 h), and the off member's "
        "3 h pre-collapse windows of the h = 1 km cases (squares: 12.34 to 15.34, 10.02 to 13.02, "
        "7.15 to 10.15 h), from s14_suite.json at full precision. Hairline: D-1 band 1e-7 over the "
        "three non-breaking pairs; the pre-collapse windows are descriptive. No h = 3 km flux is "
        "quoted. (b) Twin minus off member time of each zoverturn first onset (triangle: the t_c "
        "event; circles: later onsets in other bands and kinds, same band, kind and height in every "
        "pair) and of the end of the run (squares; cobalt: both reach t_end, ember: dt collapse). "
        "Exact zeros are drawn at the column '0'. Dotted: one scored snapshot (stride 2, about 132 s), "
        "the D-2 tolerance on the t_c event. Face path (bcface=1), pbar=1, CFL 0.4, dz 500 m except "
        "S02F1 (250 m).")
    save(fig, outdir, "fig_s02_observables", dict(
        inputs=[str(READINGS)], runs=[p[0] for p in PAIRS] + [p[1] for p in PAIRS],
        windows_h={k: MF[k]["window_h"] for k, *_ in wins}, mf_changes=side_a, event_times=side_b,
        zero_column_x=ZERO, caption=caption))
    plt.close(fig)


def fig_cost(plt, outdir):
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
    R = readings()
    C = R["cost"]
    c = np.load(CACHE)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.0, 2.95), gridspec_kw=dict(width_ratios=[1.55, 1.0]))
    lab_lu, side = [], {}
    for twin, off, dz, U, h, _ in PAIRS:
        for rid, col, lw in ((off, EMBER, 1.0), (twin, COBALT, 1.0)):
            x, y = c[f"{rid}__tsim_h"], c[f"{rid}__s_per_step"]
            a1.plot(x, y, color=col, lw=lw if dz == 500 else 0.8, zorder=3)
            tend = float(c[f"{rid}__tend_h"])
            side[rid] = dict(median_s_per_step=float(np.median(y)), first=float(y[0]), last=float(y[-1]),
                             t_end_h=tend, active_runs_min_max=[int(c[f"{rid}__active"].min()),
                                                                int(c[f"{rid}__active"].max())])
    # dz 250 m pair drawn with markers so it reads apart from the dz 500 m lines in grayscale
    for rid, col in (("S14S2", EMBER), ("S02F1", COBALT)):
        x, y = c[f"{rid}__tsim_h"], c[f"{rid}__s_per_step"]
        idx = np.linspace(0, len(x) - 1, 10).astype(int)
        a1.plot(x[idx], y[idx], ls="none", marker="o", ms=2.6, mfc=GROUND, mec=col, mew=0.7, zorder=4)
    # end of each run by dt collapse
    for twin, off, *_ in PAIRS:
        for rid, col in ((off, EMBER), (twin, COBALT)):
            tend = float(c[f"{rid}__tend_h"])
            if tend < 23.9 and rid not in ("S14S2", "S02F1"):
                x, y = c[f"{rid}__tsim_h"], c[f"{rid}__s_per_step"]
                a1.plot(x[-1], y[-1], ls="none", marker="x", ms=3.6, mew=0.9, color=col, zorder=5)
    # direct labels
    L = dict(fontsize=6.8, va="center")
    a1.text(24.4, 0.101, "PR1-U20-h0.1, PR1-U57-h0.1,\nPR1-U20-h1.0 (24 h)", color=EMBER, ha="left",
            linespacing=1.0, **L)
    # labels of the dt-collapsing :lu records, placed in clear space next to their end crosses
    for rid, (lx, ly), ha in (("PR1-U40-h3.0", (6.75, 0.40), "right"), ("PR1-U57-h1.0", (12.45, 0.75), "left"),
                              ("PR1-U20-h3.0", (12.45, 0.27), "left"), ("S15F1", (24.4, 0.66), "left")):
        x, y = c[f"{rid}__tsim_h"], c[f"{rid}__s_per_step"]
        a1.text(lx, ly, rid, color=EMBER, ha=ha, **L)
        if rid == "PR1-U20-h3.0":
            a1.plot([x[-1] + 0.15, lx - 0.15], [y[-1] * 0.97, ly * 1.02], color=tint(EMBER, 0.6), lw=0.5)
    a1.text(0.3, 0.255, r"S14S2, $\Delta z$ 250 m", color=EMBER, ha="left", **L)
    x, y = c["S02F1__tsim_h"], c["S02F1__s_per_step"]
    a1.text(x[-1] + 0.5, y[-1], r"S02F1, $\Delta z$ 250 m", color=COBALT, ha="left", **L)
    a1.text(24.4, 0.0178, "S02F2 to S02F8", color=COBALT, ha="left", **L)
    a1.text(0.3, 1.25, r":lu off members (s15, a3)", color=EMBER, ha="left", fontsize=7, va="center")
    a1.text(0.3, 0.0125, r":fft_tridiag twins (d2)", color=COBALT, ha="left", fontsize=7, va="center")
    a1.set_yscale("log")
    a1.set_ylim(0.01, 1.6)
    a1.set_xlim(0, 24)
    a1.set_xticks([0, 6, 12, 18, 24])
    a1.yaxis.set_major_locator(FixedLocator([0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]))
    a1.yaxis.set_major_formatter(FuncFormatter(lambda v, p: f"{v:g}"))
    a1.yaxis.set_minor_locator(NullLocator())
    a1.set_xlabel(r"simulated time $t$ [h]")
    a1.set_ylabel(r"wall time per step [s]")
    a1.text(0.0, 1.03, "(a)", transform=a1.transAxes, ha="left", va="bottom", fontsize=8.5, color=INK, fontweight="bold")

    # (b) whole-run s per step, off member against twin
    side_b = {}
    for i, (twin, off, dz, *_) in enumerate(PAIRS):
        q = C[twin]
        a2.plot([q["s_per_step"], q["off_s_per_step"]], [i, i], color=GRID, lw=1.6, zorder=1,
                solid_capstyle="round")
        a2.plot(q["off_s_per_step"], i, ls="none", marker="s", ms=4.0, mfc=EMBER, mec=EMBER, zorder=3)
        a2.plot(q["s_per_step"], i, ls="none", marker="o", ms=4.2, mfc=COBALT, mec=COBALT, zorder=3)
        a2.text(0.98, i, f"{q['per_step_ratio']:.2f}", transform=a2.get_yaxis_transform(), ha="right",
                va="center", fontsize=7, color=INK)
        a2.text(1.14, i, f"{q['real_ratio']:.2f}", transform=a2.get_yaxis_transform(), ha="right",
                va="center", fontsize=7, color=SLATE)
        side_b[twin] = {k: q[k] for k in ("s_per_step", "off_s_per_step", "per_step_ratio", "real_ratio",
                                          "real_s", "off_real_s", "steps", "off_steps")}
    a2.text(0.98, -1.05, "per step", transform=a2.get_yaxis_transform(), ha="right", va="center",
            fontsize=6.8, color=INK)
    a2.text(1.14, -1.05, "real", transform=a2.get_yaxis_transform(), ha="right", va="center",
            fontsize=6.8, color=SLATE)
    a2.text(1.06, -1.6, "twin / off", transform=a2.get_yaxis_transform(), ha="center", va="center",
            fontsize=6.8, color=SLATE)
    tw = sum(C[p[0]]["real_s"] for p in PAIRS) / 3600
    to = sum(C[p[0]]["off_real_s"] for p in PAIRS) / 3600
    a2.text(0.02, len(PAIRS) + 0.05, rf"$\Sigma$ real: {tw:.1f} h (twins) against {to:.1f} h (off members)",
            transform=a2.get_yaxis_transform(), ha="left", va="top", fontsize=6.8, color=SLATE)
    a2.set_xscale("log")
    a2.set_xlim(0.012, 0.32)
    a2.xaxis.set_major_locator(FixedLocator([0.02, 0.05, 0.1, 0.2]))
    a2.xaxis.set_major_formatter(FuncFormatter(lambda v, p: f"{v:g}"))
    a2.xaxis.set_minor_locator(NullLocator())
    a2.set_yticks(range(len(PAIRS)))
    a2.set_yticklabels([p[0] for p in PAIRS], fontsize=7)
    a2.set_ylim(len(PAIRS) + 0.75, -1.9)
    a2.tick_params(axis="y", length=0)
    a2.set_xlabel(r"whole-run wall time per step [s]")
    a2.text(0.0, 1.03, "(b)", transform=a2.transAxes, ha="left", va="bottom", fontsize=8.5, color=INK, fontweight="bold")
    caption = (
        "(a) Wall time per step against simulated time for the :lu off members (ember; s15 batch, "
        "671bf32, concurrency 12 at launch, and S14S2 in set a3, concurrency 11 at launch) and the "
        ":fft_tridiag twins (cobalt; set d2, df52f0e, concurrency 8, otherwise idle machine). "
        "Computed from the cost monitor histories (one row per about 64 s, logged step quantized to "
        "500) as the inverse slope of a least-squares fit of step against wall time over centred "
        "windows of at least 900 s and 4000 steps (:lu) or 300 s and 8000 steps (twins), with "
        "simulated time from the run log step lines. Crosses: end by dt collapse. Open circles: the "
        "dz 250 m pair (S14S2, S02F1). (b) Whole-run s per step (real / steps): off member (squares) "
        "and twin (circles), with the twin/off ratios of s per step and of real time. " + CAVEAT +
        " The number of active runs also fell during each batch (s15 from 12 to 3, a3 from 11 to 1, "
        "d2 from 8 to 1), which lowers the late cost of the longest runs. Face path (bcface=1), pbar=1, "
        "CFL 0.4, julia -t 1, OPENBLAS_NUM_THREADS=1.")
    save(fig, outdir, "fig_s02_cost", dict(
        inputs=[str(CACHE), str(CACHE.with_suffix(".json")), str(READINGS)],
        runs=[p[0] for p in PAIRS] + [p[1] for p in PAIRS], caveat=CAVEAT,
        cost_series=side, whole_run=side_b,
        factorizations={p[0]: int(c[f"{p[0]}__factorizations"]) for p in PAIRS},
        sum_real_h=dict(twins=tw, off_members=to), caption=caption))
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["extract", "convergence", "observables", "cost", "all"])
    ap.add_argument("--outdir", default=None, help="default: docs/figures in the main checkout")
    ap.add_argument("--font-dir", default=str(FONT_DIR))
    args = ap.parse_args()
    outdir = Path(args.outdir) if args.outdir else main_checkout() / "docs" / "figures"
    if args.cmd == "extract" or (args.cmd in ("all", "cost", "observables") and not CACHE.exists()):
        extract()
    if args.cmd == "extract":
        return
    plt = setup_style(args.font_dir)
    for name, fn in (("convergence", fig_convergence), ("observables", fig_observables), ("cost", fig_cost)):
        if args.cmd in (name, "all"):
            fn(plt, outdir)


if __name__ == "__main__":
    main()
