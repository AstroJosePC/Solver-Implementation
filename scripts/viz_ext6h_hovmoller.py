#!/usr/bin/env python3
"""
viz_ext6h_hovmoller.py

Time-height (Hovmoller) figure for a single Arsia Mons Stage 1 extended
6-hour mountain-wave run (uniform eastward steering wind U_bar = 40 m/s,
hill half-width h = 3 km). The run develops a mean-flow deceleration that
descends in altitude over the 6-hour integration: it is confined above
roughly 52 km near t = 3 h, reaches roughly 42 km by t = 4.5 h, and extends
down to roughly 20 km by t = 6 h, so the descending layer eventually
overlaps the 45-53 km cloud band used elsewhere in this project. Prior
evidence for this descent came from six discrete vertical profiles; this
figure shows it continuously in time and height using the cached
deepdive_series.npz diagnostics (no NetCDF snapshots are read here).

Three stacked panels share the time (hours) and altitude (km) axes,
restricted to 0-60 km because 60-90 km is a sponge layer and is not
physical:

  (a) mean-flow deviation <u>_x(z,t) - U_bar, with the u = 20 m/s contour
      (half the initial wind) and the locus of peak deceleration overlaid
  (b) x-mean momentum flux <rho u'w'>_x(z,t), in mPa
  (c) hill-window RMS vertical velocity w(z,t), in m/s

All three panels mark the 45-53 km cloud band and the t = 2.99 h end time
of the earlier, shorter 3-hour baseline run.

Usage:
  python3 scripts/viz_ext6h_hovmoller.py
  python3 scripts/viz_ext6h_hovmoller.py --npz <path> --out <path>
"""

from __future__ import annotations

import argparse
import os
import time

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm, PowerNorm

DEFAULT_NPZ = (
    "output/arsia_mons_stage1_ext6h_analysis/run_U40.0_h3.0/deepdive_series.npz"
)
DEFAULT_OUT = (
    "output/arsia_mons_stage1_ext6h_analysis/run_U40.0_h3.0/deepdive_hovmoller.png"
)

Z_MAX_KM = 60.0          # sponge base; exclude 60-90 km from the plotted axis
CLOUD_BAND_KM = (45.0, 53.0)
BASELINE_END_H = 2.99    # end time of the earlier 3 h run
HALF_WIND_CONTOUR = 20.0  # m/s, half of U_bar = 40 m/s (secondary contour only)
DU_CONTOUR_LEVELS = [-20.0, -10.0, -5.0, -2.0, -1.0]  # m/s, panel (a) contours
DECEL_THRESHOLD_MS = 5.0  # m/s; defines the "decelerated layer" for the
                           # headline descent curve and the descent table
DECEL_SEARCH_ZMIN_KM = 2.0  # search the descent curve/table only at or above
                             # this altitude, to exclude a near-surface,
                             # single-cell numerical artefact below it
PEAK_DECEL_MASK_MS = 2.0  # hide the peak-deceleration locus while |du_max| is
                           # below this (pre-onset startup noise, not signal)
MF_ARTEFACT_TOP_KM = 2.0  # bottom-boundary levels excluded from panel (b)'s
                           # colour-scale range (see plot_hovmoller docstring)
MF_COLOR_PERCENTILE = 99.5  # percentile of |MF| (above MF_ARTEFACT_TOP_KM)
                             # used to set panel (b)'s symmetric colour limit


def lowest_decel_altitude(
    z_km: np.ndarray, du: np.ndarray, threshold_ms: float, z_min_km: float
) -> np.ndarray:
    """For each time column, return the lowest altitude at or above
    z_min_km where du < -threshold_ms, or np.nan where no level there
    qualifies. z_km and du (z, t) must already share the same z axis.

    This is the direct statement of "how deep has the decelerated layer
    reached" and is used both for the headline curve in panel (a) and for
    the descent table. Restricting the search to z >= z_min_km excludes a
    single-cell near-surface artefact identified in the cached diagnostics
    (see module docstring / plot_hovmoller docstring).
    """
    search_mask = z_km >= z_min_km
    z_search = z_km[search_mask]
    du_search = du[search_mask, :]
    decel = du_search < -threshold_ms
    n_t = du.shape[1]
    z_lowest = np.full(n_t, np.nan)
    for j in range(n_t):
        idxs = np.where(decel[:, j])[0]
        if idxs.size:
            z_lowest[j] = z_search[idxs[0]]  # z_search is ascending
    return z_lowest


def plot_hovmoller(
    times_s: np.ndarray,
    z_km: np.ndarray,
    ubar_zt: np.ndarray,
    mf_zt: np.ndarray,
    wrms_hill_zt: np.ndarray,
    z_du_max_t: np.ndarray,
    du_max_t: np.ndarray,
    u_bar: float,
    out_path: str,
    w_norm_choice: str = "power",
) -> None:
    """Build and save the three-panel time-height figure.

    w_norm_choice selects the normalisation for panel (c):
    "power" uses a PowerNorm(gamma<1) to compress the dynamic range so
    early, weaker structure remains visible alongside late-time growth;
    "log" uses a LogNorm instead; "linear" uses no normalisation.

    Panel (a)'s headline overlay is the lowest altitude (searched at or
    above DECEL_SEARCH_ZMIN_KM) where the mean flow has slowed by more
    than DECEL_THRESHOLD_MS -- the direct visual statement of how deep the
    decelerated layer has descended. A family of du = -1, -2, -5, -10, -20
    m/s contour lines (inline-labelled) gives the shape of the deviation
    field; the u = 20 m/s contour is kept only as a secondary reference
    since it is a much stricter criterion (it only appears above ~55 km).

    Two further data-quality choices are made here rather than plotting
    the raw cached arrays unmodified:

    1) The peak-deceleration locus (z_du_max_t) is only drawn once the
       deviation at that location exceeds PEAK_DECEL_MASK_MS. Before the
       wave has done any real work on the mean flow, du is at the noise
       level everywhere, so the argmax location jumps around and would
       otherwise clutter the first ~1.6 h of panel (a) with a signal that
       is not physically meaningful.
    2) Panel (b)'s colour-scale range is computed from z at or above
       MF_ARTEFACT_TOP_KM only, as the MF_COLOR_PERCENTILE-th percentile
       of |MF| there. The lowest one or two grid levels carry a smoothly
       growing flux 5-10x larger than anywhere else in the free
       atmosphere, consistent with a lower-boundary artefact rather than
       the wave-driven flux this panel is meant to show; without the
       exclusion it swamps the colour scale and hides the structure
       aloft. The data at those levels are still plotted in full (and
       shown as clipped, out-of-range colour via colorbar extend arrows),
       just not used to set vmin/vmax.
    """
    t_h = times_s / 3600.0
    z_mask = z_km <= Z_MAX_KM
    z = z_km[z_mask]

    du = ubar_zt[z_mask, :] - u_bar
    mf_mpa = mf_zt[z_mask, :] * 1e3
    wrms = wrms_hill_zt[z_mask, :]

    fig, axes = plt.subplots(3, 1, figsize=(11, 12), sharex=True)
    ax_a, ax_b, ax_c = axes

    # --- panel (a): mean-flow deviation, diverging, centred on zero ---
    du_vmax = np.nanmax(np.abs(du))
    pa = ax_a.pcolormesh(
        t_h, z, du, cmap="coolwarm", vmin=-du_vmax, vmax=du_vmax, shading="auto"
    )
    cb_a = fig.colorbar(pa, ax=ax_a, pad=0.02)
    cb_a.set_label(r"$\langle u \rangle_x - \bar{U}$ (m/s)")

    # Family of du contour lines with inline labels, giving the shape of
    # the deviation field beyond the single headline threshold below.
    cs = ax_a.contour(
        t_h, z, du, levels=DU_CONTOUR_LEVELS,
        colors="k", linewidths=0.9, linestyles="solid",
    )
    ax_a.clabel(cs, inline=True, fontsize=7, fmt="%d")

    # Secondary reference only: u = 20 m/s is a much stricter criterion
    # than the 5 m/s deceleration threshold below, so it only appears
    # above ~55 km and must not be read as the headline result.
    cs_half = ax_a.contour(
        t_h, z, ubar_zt[z_mask, :], levels=[HALF_WIND_CONTOUR],
        colors="0.4", linewidths=1.2, linestyles="dashed",
    )

    # Headline curve: the lowest altitude (searched at/above
    # DECEL_SEARCH_ZMIN_KM) where the mean flow has slowed by more than
    # DECEL_THRESHOLD_MS. This is the direct visual statement of the
    # descent this figure exists to show.
    z_decel_lowest = lowest_decel_altitude(
        z, du, DECEL_THRESHOLD_MS, DECEL_SEARCH_ZMIN_KM
    )
    ax_a.plot(
        t_h, z_decel_lowest, color="lime", linewidth=2.8, solid_capstyle="round",
        label=f"lowest z with decel > {DECEL_THRESHOLD_MS:.0f} m/s "
              f"(z >= {DECEL_SEARCH_ZMIN_KM:.0f} km)",
    )

    peak_locus = np.where(np.abs(du_max_t) >= PEAK_DECEL_MASK_MS, z_du_max_t, np.nan)
    ax_a.plot(
        t_h, peak_locus, color="darkviolet", marker="o", markersize=3,
        linewidth=1.0, linestyle="-", label="peak deceleration altitude",
    )
    # Proxy handle so the secondary contour (which has no automatic legend
    # entry) appears in the final legend built after all overlays below.
    contour_proxy = plt.Line2D(
        [], [], color="0.4", linewidth=1.2, linestyle="dashed",
        label=f"u = {HALF_WIND_CONTOUR:.0f} m/s contour (secondary, strict)",
    )
    ax_a.set_ylabel("z (km)")
    ax_a.set_title("(a) mean-flow deviation")

    # --- panel (b): momentum flux, diverging, centred on zero ---
    # vmax is the MF_COLOR_PERCENTILE-th percentile of |MF| computed above
    # MF_ARTEFACT_TOP_KM only; see docstring above for why the lowest grid
    # levels are excluded from this range. pcolormesh clips (rather than
    # rescales) values outside [vmin, vmax] to the colormap's end colours.
    mf_range_mask = z >= MF_ARTEFACT_TOP_KM
    mf_vmax = np.nanpercentile(np.abs(mf_mpa[mf_range_mask, :]), MF_COLOR_PERCENTILE)
    pb = ax_b.pcolormesh(
        t_h, z, mf_mpa, cmap="RdBu_r", vmin=-mf_vmax, vmax=mf_vmax, shading="auto"
    )
    cb_b = fig.colorbar(pb, ax=ax_b, pad=0.02, extend="both")
    cb_b.set_label(
        r"MF = $\langle \rho u' w' \rangle_x$ (mPa)"
        f"\nscale set from z >= {MF_ARTEFACT_TOP_KM:.0f} km "
        f"({MF_COLOR_PERCENTILE:g}th pct); near-surface values exceed it"
    )
    ax_b.set_ylabel("z (km)")
    ax_b.set_title("(b) momentum flux (colour scale clipped, see colorbar note)")

    # --- panel (c): hill-window RMS w, sequential ---
    wrms_positive = np.where(wrms > 0, wrms, np.nan)
    if w_norm_choice == "power":
        norm = PowerNorm(gamma=0.5, vmin=0.0, vmax=np.nanmax(wrms))
        pc = ax_c.pcolormesh(t_h, z, wrms, cmap="viridis", norm=norm, shading="auto")
    elif w_norm_choice == "log":
        norm = LogNorm(vmin=np.nanmin(wrms_positive), vmax=np.nanmax(wrms_positive))
        pc = ax_c.pcolormesh(
            t_h, z, wrms_positive, cmap="viridis", norm=norm, shading="auto"
        )
    else:
        pc = ax_c.pcolormesh(t_h, z, wrms, cmap="viridis", shading="auto")
    cb_c = fig.colorbar(pc, ax=ax_c, pad=0.02)
    cb_c.set_label("hill-window RMS w (m/s)")
    ax_c.set_ylabel("z (km)")
    ax_c.set_xlabel("time (h)")
    ax_c.set_title("(c) hill-window RMS vertical velocity")

    # --- shared overlays: cloud band and baseline end time ---
    for i, ax in enumerate(axes):
        band_label = "cloud band (45-53 km)" if i == 0 else None
        ax.axhline(CLOUD_BAND_KM[0], color="0.15", linestyle="--", linewidth=1.0,
                    label=band_label)
        ax.axhline(CLOUD_BAND_KM[1], color="0.15", linestyle="--", linewidth=1.0)
        baseline_label = "end of 3 h baseline run" if i == 0 else None
        ax.axvline(BASELINE_END_H, color="0.15", linestyle=":", linewidth=1.2,
                    label=baseline_label)
        ax.set_xlim(0, t_h.max())
        ax.set_ylim(0, Z_MAX_KM)

    # Build panel (a)'s legend from the contour proxy plus every labelled
    # overlay on that axis (peak-deceleration line, cloud band, baseline time).
    # Placed lower right: the descending-layer signal this figure is meant
    # to show sits in the upper-right corner (late time, near the sponge),
    # so the legend must not be put there.
    handles, labels = ax_a.get_legend_handles_labels()
    ax_a.legend([contour_proxy] + handles, [contour_proxy.get_label()] + labels,
                loc="lower right", fontsize=8.5, framealpha=0.85)

    fig.suptitle(
        f"U=40 m/s, h=3 km Stage 1: time-height evolution over {t_h.max():.0f} h "
        "(sponge above 60 km excluded)",
        fontsize=13,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def print_descent_table(
    times_s: np.ndarray, z_km: np.ndarray, ubar_zt: np.ndarray, u_bar: float
) -> None:
    """Print the lowest altitude at which the mean flow has decelerated by
    more than DECEL_THRESHOLD_MS, at several requested times, quantifying
    the descent shown by the headline curve in panel (a). Uses the same
    lowest_decel_altitude() helper as the figure, searched only at or
    above DECEL_SEARCH_ZMIN_KM to exclude a single-cell near-surface
    artefact at z = 0.75 km (its immediate neighbours at z = 0.25 and
    1.25 km stay well below threshold, so it is not part of the coherent
    descending layer).
    """
    du = ubar_zt - u_bar
    query_hours = [1.0, 2.0, 3.0, 4.0, 4.5, 5.0, 6.0]
    t_h = times_s / 3600.0
    z_lowest_all = lowest_decel_altitude(
        z_km, du, DECEL_THRESHOLD_MS, DECEL_SEARCH_ZMIN_KM
    )
    print(f"t (h)   lowest z with decel > {DECEL_THRESHOLD_MS:.0f} m/s (km)  "
          f"[search z >= {DECEL_SEARCH_ZMIN_KM:.0f} km]")
    for qh in query_hours:
        idx = int(np.argmin(np.abs(t_h - qh)))
        z_lowest = z_lowest_all[idx]
        if np.isnan(z_lowest):
            print(f"{qh:5.2f}   none below threshold   (actual t = {t_h[idx]:.3f} h)")
        else:
            print(f"{qh:5.2f}   {z_lowest:6.2f}                 "
                  f"(actual t = {t_h[idx]:.3f} h)")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Time-height Hovmoller figure for the U=40 m/s, h=3 km "
        "extended 6 h Stage 1 mountain-wave run."
    )
    parser.add_argument("--npz", default=DEFAULT_NPZ, help="path to deepdive_series.npz")
    parser.add_argument("--out", default=DEFAULT_OUT, help="output PNG path")
    parser.add_argument(
        "--w-norm", default="power", choices=["power", "log", "linear"],
        help="normalisation for panel (c) hill-window RMS w",
    )
    args = parser.parse_args()

    t_start = time.time()

    data = np.load(args.npz, allow_pickle=True)
    times_s = data["times"]
    z_km = data["z_km"]
    ubar_zt = data["ubar_zt"]
    mf_zt = data["mf_zt"]
    wrms_hill_zt = data["wrms_hill_zt"]
    z_du_max_t = data["z_du_max_t"]
    du_max_t = data["du_max_t"]
    u_bar = float(data["u_bar"])

    plot_hovmoller(
        times_s, z_km, ubar_zt, mf_zt, wrms_hill_zt, z_du_max_t, du_max_t, u_bar,
        args.out, w_norm_choice=args.w_norm,
    )

    elapsed = time.time() - t_start
    print(f"wrote {args.out}")
    print(f"script wall-clock duration: {elapsed:.2f} s")
    print()
    print_descent_table(times_s, z_km, ubar_zt, u_bar)


if __name__ == "__main__":
    main()
