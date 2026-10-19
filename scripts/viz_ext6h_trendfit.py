#!/usr/bin/env python3
"""
viz_ext6h_trendfit.py

Visual proof of the least-squares trend statistics reported for the
U=40 m/s, h=3 km Arsia Mons Stage 1 extended (6 h) run. A previous analysis
quoted normalized trends ("frac_per_hour") for five diagnostics over two
time windows but only ever plotted the raw series, never the fitted lines
behind those numbers. This script draws both, so the claimed early-transient
decay and late-time flattening/plateau can be checked by eye against the
fits rather than taken on faith.

Trend statistic (applied identically to every row and both windows)
---------------------------------------------------------------------
For a series q(t) restricted to a window [t0, t1]:
    slope_per_hour  = np.polyfit(t, q, 1)[0] * 3600
    scale           = mean(|q|) over the window
    frac_per_hour   = slope_per_hour / scale        (normalized trend)
    cov             = std(q) / scale over the window

Two windows are evaluated:
    Window A, "2-3 h": t0=7200 s,  t1=10800 s
        (the final hour of the original, un-extended 3 h run)
    Window B, "last hour": t0=times[-1]-3600 s, t1=times[-1]
        (the final hour of the 6 h extended run)

Each row of the figure shows the raw 6 h series (thin grey line), the two
windows shaded, the actual fitted line drawn thick over its own window and
extended a little (about 0.7 h) beyond it as a thin dashed line, and a text
box with the slope/frac/cov numbers for both windows. Rows 1 (AMEC hill
mean |w|) and 5 (peak mean-flow deviation) additionally overlay the shorter
3 h baseline run of the same configuration, to show the two runs agree over
their shared interval before the extended run diverges from it (or fails
to).

Usage:
    python3 scripts/viz_ext6h_trendfit.py
    python3 scripts/viz_ext6h_trendfit.py --npz <path> --baseline-npz <path> --out <path>
"""

from __future__ import annotations

import argparse
import os
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

DEFAULT_NPZ = (
    "output/arsia_mons_stage1_ext6h_analysis/run_U40.0_h3.0/deepdive_series.npz"
)
DEFAULT_BASELINE_NPZ = (
    "output/arsia_mons_stage1_analysis/run_U40.0_h3.0/deepdive_series.npz"
)
DEFAULT_OUT = (
    "output/arsia_mons_stage1_ext6h_analysis/run_U40.0_h3.0/deepdive_trendfit.png"
)

MF_LEVELS_KM = [30.0, 45.0, 52.0]
RUN_END_3H_HOURS = 2.99

WINDOW_A_COLOR = "tab:blue"
WINDOW_B_COLOR = "tab:red"
WINDOW_A_SPAN_COLOR = "#c6dbef"
WINDOW_B_SPAN_COLOR = "#fcbba1"


def trend_stats(times: np.ndarray, series: np.ndarray, t0: float, t1: float) -> dict:
    """Least-squares slope of `series` over [t0, t1], absolute and normalized.

    See module docstring for the exact definition. Returns None fields if the
    window contains fewer than 3 samples or the scale is zero.
    """
    sel = (times >= t0) & (times <= t1)
    if sel.sum() < 3:
        return {"n": int(sel.sum()), "slope_per_hour": None, "frac_per_hour": None,
                "cov": None, "poly": None}
    tt, ss = times[sel], series[sel]
    poly = np.polyfit(tt, ss, 1)
    slope_per_hour = float(poly[0]) * 3600.0
    scale = float(np.abs(ss).mean())
    cov = float(ss.std() / scale) if scale > 0 else None
    frac_per_hour = (slope_per_hour / scale) if scale > 0 else None
    return {"n": int(sel.sum()), "slope_per_hour": slope_per_hour,
            "frac_per_hour": frac_per_hour, "cov": cov, "poly": poly}


def _fit_line_hours(poly: np.ndarray, t0: float, t1: float) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate a np.polyfit line (in seconds) at the window endpoints, hours out."""
    tt_s = np.array([t0, t1])
    yy = np.polyval(poly, tt_s)
    return tt_s / 3600.0, yy


def draw_trendfit_row(
    ax,
    times_s: np.ndarray,
    series: np.ndarray,
    windows: list[tuple[float, float, str]],
    unit_label: str,
    fmt_slope: str = "{:+.2e}",
    baseline_times_s: np.ndarray | None = None,
    baseline_series: np.ndarray | None = None,
    xlim_hours: tuple[float, float] = (0.0, 6.0),
    extend_hours: float = 0.7,
) -> dict:
    """Draw one row of the trend-fit figure and return the trend stats.

    `windows` is a list of (t0_s, t1_s, color) tuples; window shading and
    labelling of the spans is handled by the caller (once, on the top row)
    so the legend is not duplicated across rows.

    The y axis is widened beyond the true data range to reserve a blank
    strip along the top, reserved purely for the annotation box, so the box
    never sits on top of the raw series, a fitted line, or its extension.
    The strip is computed from the actual plotted values (not by hard-coding
    a corner), so it adapts row by row. Window shading and the 3 h marker
    are clipped to the data sub-range so they do not bleed into that strip.
    """
    t_hours = times_s / 3600.0
    ax.plot(t_hours, series, color="0.5", lw=1.0, zorder=2)

    # Track every y-value that will actually be drawn as a line, so the
    # reserved annotation strip is guaranteed clear of all of them.
    plotted_values = [series]
    if baseline_series is not None:
        plotted_values.append(baseline_series)

    stats_by_window = {}
    text_lines = []
    for (t0, t1, color), label in zip(windows, ["A", "B"]):
        stats = trend_stats(times_s, series, t0, t1)
        stats_by_window[label] = stats
        if stats["poly"] is None:
            continue
        poly = stats["poly"]

        # Solid fit line over the window itself.
        xh, yy = _fit_line_hours(poly, t0, t1)
        ax.plot(xh, yy, color=color, lw=2.5, zorder=4)

        # Thin dashed extension beyond the window, clipped to the axes.
        t0_ext = max(t0 - extend_hours * 3600.0, xlim_hours[0] * 3600.0)
        t1_ext = min(t1 + extend_hours * 3600.0, xlim_hours[1] * 3600.0)
        xh_pre, yy_pre = _fit_line_hours(poly, t0_ext, t0)
        ax.plot(xh_pre, yy_pre, color=color, lw=1.0, ls="--", alpha=0.7, zorder=3)
        xh_post, yy_post = _fit_line_hours(poly, t1, t1_ext)
        ax.plot(xh_post, yy_post, color=color, lw=1.0, ls="--", alpha=0.7, zorder=3)

        # The fit is linear, so its extremes over [t0_ext, t1_ext] are at
        # the four endpoints already computed above.
        plotted_values.append(np.array([yy_pre[0], yy[0], yy[1], yy_post[-1]]))

        win_hours = f"{t0/3600.0:.2f}-{t1/3600.0:.2f} h"
        text_lines.append(
            f"{win_hours}:  {fmt_slope.format(stats['slope_per_hour'])} {unit_label}/h   "
            f"frac {stats['frac_per_hour']:+.2f}   cov {stats['cov']:.2f}"
        )

    if baseline_times_s is not None and baseline_series is not None:
        ax.plot(baseline_times_s / 3600.0, baseline_series, color="black", lw=1.2,
                ls="--", alpha=0.35, zorder=2)

    # Reserve a blank strip above the true data range for the annotation,
    # and clip the window shading / run-end marker to the data sub-range so
    # neither bleeds into that strip.
    all_values = np.concatenate([np.asarray(v).ravel() for v in plotted_values])
    data_min, data_max = float(np.nanmin(all_values)), float(np.nanmax(all_values))
    data_range = data_max - data_min
    if data_range <= 0:
        data_range = max(abs(data_max), 1.0)
    bottom_pad = 0.06 * data_range
    top_pad = 0.32 * data_range
    y0 = data_min - bottom_pad
    y1 = data_max + top_pad
    data_top_frac = (data_max - y0) / (y1 - y0)

    for t0, t1, color in windows:
        ax.axvspan(t0 / 3600.0, t1 / 3600.0, ymin=0.0, ymax=data_top_frac,
                   color=color, alpha=0.15, zorder=1)
    ax.axvline(RUN_END_3H_HOURS, ymin=0.0, ymax=data_top_frac,
               color="0.2", lw=1.0, ls=":", zorder=3)

    ax.set_xlim(*xlim_hours)
    ax.set_ylim(y0, y1)
    ax.grid(alpha=0.3)

    ax.text(
        0.985, 0.99, "\n".join(text_lines),
        transform=ax.transAxes, ha="right", va="top", fontsize=8,
        family="monospace",
        bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.6", alpha=0.95),
        zorder=5,
    )
    return stats_by_window


def make_figure(
    npz_path: str,
    baseline_npz_path: str | None,
    out_path: str,
) -> dict:
    """Build the 5-row trend-fit figure and save it to `out_path`.

    Returns a nested dict {quantity: {"A": stats, "B": stats}} for reporting.
    """
    d = np.load(npz_path, allow_pickle=True)
    times = d["times"]
    z_km = d["z_km"]
    mf_zt = d["mf_zt"]
    amec_w_t = d["amec_w_t"]
    du_max_t = d["du_max_t"]
    u_bar = float(d["u_bar"])

    baseline = None
    if baseline_npz_path is not None and os.path.exists(baseline_npz_path):
        baseline = np.load(baseline_npz_path, allow_pickle=True)

    win_A = (7200.0, 10800.0, WINDOW_A_COLOR)
    win_B = (float(times[-1] - 3600.0), float(times[-1]), WINDOW_B_COLOR)
    windows = [win_A, win_B]

    iz_levels = [int(np.argmin(np.abs(z_km - lv))) for lv in MF_LEVELS_KM]

    fig, axes = plt.subplots(5, 1, figsize=(11, 13), sharex=True)

    all_stats = {}

    # Row 1: AMEC hill mean |w|, with 3 h baseline overlay.
    baseline_t = baseline_w = None
    if baseline is not None:
        baseline_t = baseline["times"]
        baseline_w = baseline["amec_w_t"]
    all_stats["amec_w"] = draw_trendfit_row(
        axes[0], times, amec_w_t, windows, unit_label="m/s",
        baseline_times_s=baseline_t, baseline_series=baseline_w,
    )
    axes[0].set_ylabel("AMEC hill\nmean |w| (m/s)")

    # Rows 2-4: MF at 30, 45, 52 km.
    for row, (iz, lv) in enumerate(zip(iz_levels, MF_LEVELS_KM), start=1):
        mf_mpa = mf_zt[iz, :] * 1e3
        all_stats[f"mf_{int(lv)}km"] = draw_trendfit_row(
            axes[row], times, mf_mpa, windows, unit_label="mPa",
        )
        axes[row].set_ylabel(f"MF at {int(lv)} km\n(mPa)")

    # Row 5: peak mean-flow deviation, with 3 h baseline overlay.
    baseline_t = baseline_du = None
    if baseline is not None:
        baseline_t = baseline["times"]
        baseline_du = baseline["du_max_t"]
    all_stats["du_max"] = draw_trendfit_row(
        axes[4], times, du_max_t, windows, unit_label="m/s",
        baseline_times_s=baseline_t, baseline_series=baseline_du,
    )
    axes[4].set_ylabel("peak <u>_x - U_bar\n(m/s)")

    axes[-1].set_xlabel("time (hours)")

    # Build one combined legend, attached to the top row, covering the
    # elements common to every row (window spans, fits, extensions, the
    # 3 h run-end marker, and the baseline overlay where present).
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D

    legend_handles = [
        Patch(facecolor=WINDOW_A_SPAN_COLOR, alpha=0.6, label="window A (2-3 h)"),
        Patch(facecolor=WINDOW_B_SPAN_COLOR, alpha=0.6, label="window B (last hour)"),
        Line2D([0], [0], color=WINDOW_A_COLOR, lw=2.5, label="fit, window A"),
        Line2D([0], [0], color=WINDOW_B_COLOR, lw=2.5, label="fit, window B"),
        Line2D([0], [0], color="0.4", lw=1.0, ls="--", label="fit extended (0.7 h)"),
        Line2D([0], [0], color="0.5", lw=1.0, label="raw series"),
        Line2D([0], [0], color="0.2", lw=1.0, ls=":", label="3 h run end"),
    ]
    if baseline is not None:
        legend_handles.append(
            Line2D([0], [0], color="black", lw=1.2, ls="--", alpha=0.5,
                   label="3 h baseline run")
        )
    axes[0].legend(handles=legend_handles, loc="upper center",
                    bbox_to_anchor=(0.5, 1.35), ncol=4, fontsize=8, frameon=True)

    fig.suptitle(
        f"U={u_bar:.0f} m/s, h=3 km Stage 1: least-squares trends over the "
        "2-3 h and 5-6 h windows",
        fontsize=13, y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.965))

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)

    return all_stats


def print_stats_table(all_stats: dict) -> None:
    header = f"{'quantity':<12} {'window':<8} {'slope/h':>14} {'frac/h':>10} {'cov':>8}"
    print(header)
    print("-" * len(header))
    for qty, by_window in all_stats.items():
        for win_label, stats in by_window.items():
            slope = stats["slope_per_hour"]
            frac = stats["frac_per_hour"]
            cov = stats["cov"]
            print(f"{qty:<12} {win_label:<8} {slope:>14.4e} {frac:>10.3f} {cov:>8.3f}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--npz", default=DEFAULT_NPZ,
                         help="Path to the extended-run deepdive_series.npz cache")
    parser.add_argument("--baseline-npz", default=DEFAULT_BASELINE_NPZ,
                         help="Path to the 3 h baseline deepdive_series.npz cache")
    parser.add_argument("--out", default=DEFAULT_OUT,
                         help="Output PNG path")
    args = parser.parse_args()

    t_start = time.time()
    all_stats = make_figure(args.npz, args.baseline_npz, args.out)
    duration_s = time.time() - t_start

    print_stats_table(all_stats)
    print(f"\nWrote figure to {args.out}")
    print(f"Wall-clock duration: {duration_s:.2f} s")


if __name__ == "__main__":
    main()
