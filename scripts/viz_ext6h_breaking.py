"""Wave-breaking indicators for the extended 6 h Arsia Mons Stage 1 run
(U = 40 m/s, h = 3 km).

This figure summarizes the full 323-snapshot deepdive record cached in
deepdive_series.npz, replacing an earlier six-point hand-sampled check.
It shows four time series (statically unstable area fraction, minimum
hill-window zonal wind, maximum hill-window vertical velocity, and minimum
zonal-mean zonal wind) plus a time-height map of the unstable fraction at
each level. Static instability means N^2 < 0, where N is the buoyancy
frequency computed from potential temperature; a value below zero at a grid
column indicates a locally overturned (statically unstable) layer.

Usage:
    python viz_ext6h_breaking.py [--npz PATH] [--out PATH]
"""

import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DEFAULT_NPZ = (
    "output/arsia_mons_stage1_ext6h_analysis/run_U40.0_h3.0/deepdive_series.npz"
)
DEFAULT_OUT = (
    "output/arsia_mons_stage1_ext6h_analysis/run_U40.0_h3.0/deepdive_breaking.png"
)

SHORT_RUN_END_H = 2.99
SPONGE_BASE_KM = 60.0
# Colour scale ceiling for the time-height panel. Set explicitly (well above
# the 98th percentile of the field below 57 km) rather than at the raw data
# maximum, because the raw maximum is set by a narrow band just below the
# sponge base (58 to 60 km) that would otherwise wash out the cloud-band
# structure (45 to 53 km) that is the region of scientific interest. Values
# above this ceiling are shown as the top colour and reported separately.
HOV_VMAX = 0.2


def make_breaking_figure(data, out_path):
    """Build and save the four-row-plus-heatmap wave-breaking figure.

    Parameters
    ----------
    data : dict-like
        Must provide times (s), z_km (km), unstable_frac_t, umin_hill_t,
        wmax_hill_t, umean_min_t, unstable_frac_zt, u_bar, sponge_base_km,
        amec_band_km.
    out_path : str
        Destination PNG path.
    """
    t_h = np.asarray(data["times"]) / 3600.0
    z_km = np.asarray(data["z_km"])
    unstable_frac_t = np.asarray(data["unstable_frac_t"])
    umin_hill_t = np.asarray(data["umin_hill_t"])
    wmax_hill_t = np.asarray(data["wmax_hill_t"])
    umean_min_t = np.asarray(data["umean_min_t"])
    unstable_frac_zt = np.asarray(data["unstable_frac_zt"])
    u_bar = float(data["u_bar"])
    amec_band_km = np.asarray(data["amec_band_km"])

    z_mask = z_km <= SPONGE_BASE_KM

    fig = plt.figure(figsize=(11, 14), dpi=150, constrained_layout=True)
    gs = fig.add_gridspec(5, 1, height_ratios=[1, 1, 1, 1, 2.2], hspace=0.08)

    axes = [fig.add_subplot(gs[i]) for i in range(4)]
    ax_hov = fig.add_subplot(gs[4])

    # Row 1: unstable area fraction
    ax = axes[0]
    ax.plot(t_h, unstable_frac_t, color="tab:red")
    ax.set_ylabel("unstable area\nfraction (below 60 km)")
    ax.grid(alpha=0.3)

    # Row 2: minimum hill-window u, mark reversed flow
    ax = axes[1]
    ax.plot(t_h, umin_hill_t, color="tab:blue")
    ax.axhline(0.0, color="k", linewidth=1, linestyle="-")
    ax.fill_between(
        t_h, umin_hill_t, 0.0, where=(umin_hill_t < 0.0),
        color="tab:blue", alpha=0.15, interpolate=True,
    )
    ax.annotate(
        "below 0: reversed flow",
        xy=(0.02, 0.06), xycoords="axes fraction",
        fontsize=8, color="tab:blue",
    )
    ax.set_ylabel("min hill-window u\n(m/s)")
    ax.grid(alpha=0.3)

    # Row 3: maximum hill-window |w|, mark exceedance of background wind
    ax = axes[2]
    ax.plot(t_h, wmax_hill_t, color="tab:green")
    ax.axhline(u_bar, color="k", linewidth=1, linestyle="--")
    ax.annotate(
        f"initial wind U = {u_bar:.0f} m/s",
        xy=(0.02, 0.88), xycoords="axes fraction",
        fontsize=8, color="k",
    )
    ax.set_ylabel("max hill-window |w|\n(m/s)")
    ax.grid(alpha=0.3)

    # Row 4: minimum zonal-mean wind, self-induced critical level at 0
    ax = axes[3]
    ax.plot(t_h, umean_min_t, color="tab:purple")
    ax.axhline(0.0, color="k", linewidth=1, linestyle="--")
    ax.annotate(
        "self-induced critical level",
        xy=(0.02, 0.06), xycoords="axes fraction",
        fontsize=8, color="k",
    )
    ymin = min(-2.0, 1.1 * np.min(umean_min_t))
    ymax = 1.05 * max(u_bar, np.max(umean_min_t))
    ax.set_ylim(ymin, ymax)
    ax.set_ylabel("min zonal-mean u\n(m/s)")
    ax.grid(alpha=0.3)

    for ax in axes:
        ax.axvline(SHORT_RUN_END_H, color="gray", linewidth=1, linestyle=":")
    axes[-1].set_xlabel("time (h)")
    for ax in axes[:-1]:
        ax.tick_params(labelbottom=False)
    for ax in axes:
        ax.set_xlim(t_h.min(), t_h.max())

    # Bottom panel: time-height map of unstable fraction, restricted to
    # 0-60 km (the region above is a numerical sponge, not physical). The
    # colour scale is clipped at HOV_VMAX (see constant definition above);
    # this does not change the plotted altitude range, only the colouring.
    pc = ax_hov.pcolormesh(
        t_h, z_km[z_mask], unstable_frac_zt[z_mask, :],
        shading="auto", cmap="viridis", vmin=0.0, vmax=HOV_VMAX,
    )
    ax_hov.set_ylim(0.0, SPONGE_BASE_KM)
    ax_hov.set_xlim(t_h.min(), t_h.max())
    ax_hov.set_xlabel("time (h)")
    ax_hov.set_ylabel("altitude z (km)")
    ax_hov.set_title(
        f"colour scale clipped at {HOV_VMAX:.2f}; values near the sponge "
        "base (58 to 60 km) exceed this",
        fontsize=9,
    )
    ax_hov.axvline(SHORT_RUN_END_H, color="white", linewidth=1, linestyle=":")
    for zc in amec_band_km:
        ax_hov.axhline(zc, color="white", linewidth=1, linestyle="--")
    cbar = fig.colorbar(pc, ax=ax_hov, pad=0.02, extend="max")
    cbar.set_label(f"unstable column fraction (clipped at {HOV_VMAX:.2f})")

    fig.suptitle(
        f"U=40 m/s, h=3 km Stage 1: wave-breaking indicators over {t_h.max():.0f} h "
        "(sponge above 60 km excluded)",
        fontsize=12,
    )
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def report_indicators(data):
    """Print the summary tables and first-crossing times used for verification."""
    t_h = np.asarray(data["times"]) / 3600.0
    z_km = np.asarray(data["z_km"])
    unstable_frac_t = np.asarray(data["unstable_frac_t"])
    umin_hill_t = np.asarray(data["umin_hill_t"])
    wmax_hill_t = np.asarray(data["wmax_hill_t"])
    umean_min_t = np.asarray(data["umean_min_t"])
    unstable_frac_zt = np.asarray(data["unstable_frac_zt"])

    print("\nTable: indicators at selected times")
    header = (
        f"{'t (h)':>6} {'unstable_frac':>14} {'umin_hill (m/s)':>17} "
        f"{'wmax_hill (m/s)':>17} {'umean_min (m/s)':>17}"
    )
    print(header)
    for target in [1, 2, 3, 4, 5, 6]:
        idx = int(np.argmin(np.abs(t_h - target)))
        print(
            f"{t_h[idx]:6.3f} {unstable_frac_t[idx]:14.5f} "
            f"{umin_hill_t[idx]:17.3f} {wmax_hill_t[idx]:17.3f} "
            f"{umean_min_t[idx]:17.3f}"
        )

    print("\nFirst-crossing times:")

    def first_time(mask):
        idxs = np.nonzero(mask)[0]
        return t_h[idxs[0]] if idxs.size else None

    t_umin_neg = first_time(umin_hill_t < 0.0)
    t_unstable = first_time(unstable_frac_t > 0.001)
    t_wmax = first_time(wmax_hill_t > 40.0)
    print(f"  umin_hill_t first negative:        t = {t_umin_neg:.3f} h")
    print(f"  unstable_frac_t first > 0.001:     t = {t_unstable:.3f} h")
    print(f"  wmax_hill_t first > 40 m/s:        t = {t_wmax:.3f} h")

    z_mask = z_km <= SPONGE_BASE_KM
    final_col = unstable_frac_zt[z_mask, -1]
    z_argmax = z_km[z_mask][int(np.argmax(final_col))]
    print(
        f"  altitude of max unstable fraction at t=6h (below 60 km): "
        f"z = {z_argmax:.2f} km (value {np.max(final_col):.5f})"
    )

    i59 = int(np.argmin(np.abs(z_km - 59.0)))
    i50 = int(np.argmin(np.abs(z_km - 50.0)))
    print(
        f"  unstable fraction at z={z_km[i59]:.2f} km (near sponge base), "
        f"t=6h: {unstable_frac_zt[i59, -1]:.4f}"
    )
    print(
        f"  unstable fraction at z={z_km[i50]:.2f} km (cloud band), "
        f"t=6h: {unstable_frac_zt[i50, -1]:.4f}"
    )

    print("\nSanity check at t = 6 h (last snapshot):")
    print(f"  unstable_frac_t   = {unstable_frac_t[-1]:.5f}  (expected ~0.086)")
    print(f"  umin_hill_t       = {umin_hill_t[-1]:.3f} m/s  (expected ~-64)")
    print(f"  wmax_hill_t       = {wmax_hill_t[-1]:.3f} m/s  (expected ~51)")
    print(f"  umean_min_t       = {umean_min_t[-1]:.3f} m/s  (expected ~13.9)")


def main():
    parser = argparse.ArgumentParser(
        description="Plot wave-breaking indicators for the ext6h Arsia Mons "
        "Stage 1 U=40 h=3 run."
    )
    parser.add_argument("--npz", default=DEFAULT_NPZ, help="path to deepdive_series.npz")
    parser.add_argument("--out", default=DEFAULT_OUT, help="output PNG path")
    args = parser.parse_args()

    data = np.load(args.npz, allow_pickle=True)
    make_breaking_figure(data, args.out)
    report_indicators(data)
    print(f"\nSaved figure to {args.out}")


if __name__ == "__main__":
    main()
