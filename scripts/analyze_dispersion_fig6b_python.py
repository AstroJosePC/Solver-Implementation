#!/usr/bin/env python3
"""
Figure 6b diagnostics in Python (strict run compatible)

Generates three visualizations from the same simulation outputs:
1) Column windowed FFT in z (per time)
2) Sliding-window FFT in z with max pooling across windows
3) Raw-w' wavelet panel (max-over-z pooling), with density weighting applied
   only after extracting dominant lambda_z (to score ridge strength)

Usage:
  python scripts/analyze_dispersion_fig6b_python.py \
    --indir output/dispersion_test_fig6_strict \
    --outdir output/dispersion_analysis
"""

from __future__ import annotations

import argparse
import glob
import os
from dataclasses import dataclass

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm

from scipy.signal import cwt, morlet2, get_window
from netCDF4 import Dataset


@dataclass
class RunData:
    times: np.ndarray        # (nt,)
    z: np.ndarray            # (nz,)
    x: np.ndarray            # (nx,)
    w_raw: np.ndarray        # (nt, nz) at mid-x
    w_scaled: np.ndarray     # (nt, nz) = w * sqrt(rho/rho0)
    rho_weight: np.ndarray   # (nt, nz) = sqrt(rho/rho0)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--indir", default="output/dispersion_test_fig6_strict")
    p.add_argument("--outdir", default="output/dispersion_analysis")
    p.add_argument("--tmin", type=float, default=1800.0)
    p.add_argument("--lambda-min-km", type=float, default=5.0)
    p.add_argument("--lambda-max-km", type=float, default=40.0)
    p.add_argument("--n-lambda", type=int, default=280)
    p.add_argument("--fft-pad", type=int, default=8)
    p.add_argument("--window-km", type=float, default=30.0)
    p.add_argument("--step-km", type=float, default=2.5)
    p.add_argument("--morlet-w", type=float, default=8.0)
    p.add_argument("--lambda-expected-km", type=float, default=16.7)
    p.add_argument("--env-mult", type=float, default=2.5,
                   help="Envelope RMS window = env_mult * lambda_expected")
    p.add_argument("--env-floor-frac", type=float, default=0.10,
                   help="Minimum envelope floor as fraction of per-profile max envelope")
    return p.parse_args()


def load_midx_profiles(indir: str) -> RunData:
    files = sorted(glob.glob(os.path.join(indir, "snapshot_*.nc")))
    if not files:
        raise FileNotFoundError(f"No snapshot_*.nc found in {indir}")

    times = []
    w_raw_list = []
    w_scaled_list = []
    rho_weight_list = []

    x = z = None
    ix_mid = None

    for fp in files:
        with Dataset(fp, "r") as ds:
            if x is None:
                x = np.array(ds.variables["x"][:], dtype=float)
                z = np.array(ds.variables["z"][:], dtype=float)
                ix_mid = int(np.argmin(np.abs(x - 0.5 * (x.min() + x.max()))))

            w = np.array(ds.variables["w"][:, :], dtype=float).T
            rho = np.array(ds.variables["rho"][:, :], dtype=float).T
            t = float(ds.getncattr("time_s"))

            rho0 = float(np.mean(rho[:, 0]))
            weight = np.sqrt(np.clip(rho[ix_mid, :] / rho0, 0.0, None))
            w_mid = w[ix_mid, :]

            times.append(t)
            w_raw_list.append(w_mid)
            w_scaled_list.append(w_mid * weight)
            rho_weight_list.append(weight)

    times = np.array(times, dtype=float)
    order = np.argsort(times)

    return RunData(
        times=times[order],
        z=z,
        x=x,
        w_raw=np.array(w_raw_list, dtype=float)[order, :],
        w_scaled=np.array(w_scaled_list, dtype=float)[order, :],
        rho_weight=np.array(rho_weight_list, dtype=float)[order, :],
    )


def interp_to_lambda_grid(lambda_native_m: np.ndarray, power_native: np.ndarray, lambda_grid_m: np.ndarray) -> np.ndarray:
    order = np.argsort(lambda_native_m)
    lam = lambda_native_m[order]
    pow_ = power_native[order]
    return np.interp(lambda_grid_m, lam, pow_, left=pow_[0], right=pow_[-1])


def fft_spectrum_lambda(profile: np.ndarray, dz: float, pad_factor: int = 8) -> tuple[np.ndarray, np.ndarray]:
    y = profile - np.mean(profile)
    win = get_window("hann", len(y), fftbins=True)
    y = y * win

    nfft = int(pad_factor * len(y))
    Y = np.fft.rfft(y, n=nfft)
    f_cyc_per_m = np.fft.rfftfreq(nfft, d=dz)

    # remove zero-frequency
    f = f_cyc_per_m[1:]
    p = np.abs(Y[1:]) ** 2

    lam = 1.0 / np.clip(f, 1e-15, None)
    return lam, p


def moving_rms(y: np.ndarray, win_n: int) -> np.ndarray:
    n = int(y.size)
    win_n = max(3, int(win_n))
    win_n = min(win_n, n if n % 2 == 1 else max(1, n - 1))
    if win_n % 2 == 0:
        win_n += 1
    win_n = min(win_n, n)
    kern = np.ones(win_n, dtype=float) / float(win_n)
    return np.sqrt(np.convolve(y * y, kern, mode="same"))


def preprocess_profiles(data: RunData, lambda_expected_km: float, env_mult: float,
                        env_floor_frac: float) -> np.ndarray:
    """
    Preprocess w'(z,t) profiles to suppress stationary/background structure:
      1) subtract time-mean profile over all snapshots
      2) divide each time profile by its smoothed amplitude envelope (running RMS)
    """
    dz = float(data.z[1] - data.z[0])
    lam_exp_m = lambda_expected_km * 1e3
    win_n = max(5, int(round((env_mult * lam_exp_m) / dz)))

    w_mean = np.mean(data.w_raw, axis=0)  # stationary vertical structure
    out = np.zeros_like(data.w_raw)

    for it in range(data.w_raw.shape[0]):
        y = data.w_raw[it, :] - w_mean
        env = moving_rms(y, win_n)
        env_floor = max(1e-12, env_floor_frac * float(np.max(env)))
        out[it, :] = y / np.maximum(env, env_floor)

    return out


def panel_fft_column(profiles_proc: np.ndarray, z: np.ndarray, idx_t: np.ndarray,
                     lambda_grid_m: np.ndarray, pad_factor: int) -> np.ndarray:
    dz = float(z[1] - z[0])
    panel = np.zeros((len(lambda_grid_m), len(idx_t)), dtype=float)

    for col, it in enumerate(idx_t):
        lam, p = fft_spectrum_lambda(profiles_proc[it, :], dz, pad_factor)
        pi = interp_to_lambda_grid(lam, p, lambda_grid_m)
        panel[:, col] = pi

    m = float(np.max(panel))
    if m > 0:
        panel = panel / m
    return np.clip(panel, 0.0, 1.0)


def panel_fft_sliding(profiles_proc: np.ndarray, z: np.ndarray, idx_t: np.ndarray,
                      lambda_grid_m: np.ndarray, pad_factor: int,
                      window_km: float, step_km: float) -> np.ndarray:
    dz = float(z[1] - z[0])
    nz = len(z)

    win_n = max(16, int(round((window_km * 1e3) / dz)))
    step_n = max(1, int(round((step_km * 1e3) / dz)))
    starts = list(range(0, max(1, nz - win_n + 1), step_n))

    panel = np.zeros((len(lambda_grid_m), len(idx_t)), dtype=float)

    for col, it in enumerate(idx_t):
        p_stack = []
        prof = profiles_proc[it, :]
        for s in starts:
            seg = prof[s:s + win_n]
            if len(seg) < win_n:
                continue
            lam, p = fft_spectrum_lambda(seg, dz, pad_factor)
            pi = interp_to_lambda_grid(lam, p, lambda_grid_m)
            p_stack.append(pi)

        if p_stack:
            pmax = np.max(np.stack(p_stack, axis=0), axis=0)
        else:
            pmax = np.zeros_like(lambda_grid_m)

        panel[:, col] = pmax

    m = float(np.max(panel))
    if m > 0:
        panel = panel / m
    return np.clip(panel, 0.0, 1.0)


def panel_wavelet_raw_then_weight(data: RunData, profiles_proc: np.ndarray,
                                  idx_t: np.ndarray, lambda_grid_m: np.ndarray,
                                  morlet_w: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Returns:
      panel_raw_norm(lambda, t): raw-w' wavelet spectra with max-over-z pooling
      lambda_dom_m(t): dominant wavelength from raw panel
      strength_weighted(t): density-weighted dominant-mode strength (post-extraction)
    """
    dz = float(data.z[1] - data.z[0])

    # scipy widths are in samples; for morlet2, pseudo-frequency ~ w/(2*pi*width)
    # => lambda = 1/f = (2*pi*width*dz)/w
    widths = np.linspace(2.0, max(3.0, len(data.z) / 3.0), 200)
    lambda_native = (2.0 * np.pi * widths * dz) / float(morlet_w)

    panel = np.zeros((len(lambda_grid_m), len(idx_t)), dtype=float)
    lambda_dom = np.zeros(len(idx_t), dtype=float)
    strength_weighted = np.zeros(len(idx_t), dtype=float)

    for col, it in enumerate(idx_t):
        y = profiles_proc[it, :]

        # CWT matrix shape (nwidth, nz)
        W = cwt(y, morlet2, widths, w=morlet_w)
        P = np.abs(W) ** 2

        # max-over-altitude pooling
        p_lam = np.max(P, axis=1)
        p_interp = interp_to_lambda_grid(lambda_native, p_lam, lambda_grid_m)

        # dominant lambda from RAW panel
        i_dom = int(np.argmax(p_interp))
        lambda_dom[col] = float(lambda_grid_m[i_dom])

        # post-extraction density weighting: score dominant mode only
        # use scale closest to dominant lambda, then pool weighted-by-rho over z
        i_native = int(np.argmin(np.abs(lambda_native - lambda_dom[col])))
        mode_power_z = P[i_native, :]
        w_rho = np.clip(data.rho_weight[it, :], 0.0, None)
        denom = float(np.sum(w_rho))
        strength_weighted[col] = float(np.sum(mode_power_z * w_rho) / denom) if denom > 0 else 0.0

        panel[:, col] = p_interp

    m = float(np.max(panel))
    if m > 0:
        panel = panel / m

    # normalize weighted strength for overlay
    smax = float(np.max(strength_weighted))
    if smax > 0:
        strength_weighted = strength_weighted / smax

    return np.clip(panel, 0.0, 1.0), lambda_dom, strength_weighted


def dominant_lambda(panel: np.ndarray, lambda_grid_m: np.ndarray) -> np.ndarray:
    idx = np.argmax(panel, axis=0)
    return lambda_grid_m[idx]


def save_three_panel_figure(times_min: np.ndarray, lambda_grid_km: np.ndarray,
                            p1: np.ndarray, p2: np.ndarray, p3: np.ndarray,
                            l1_km: np.ndarray, l2_km: np.ndarray, l3_km: np.ndarray,
                            s3: np.ndarray, out_png: str) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(12, 13), constrained_layout=True)

    def draw(ax, panel, lam_dom_km, title):
        im = ax.imshow(
            panel,
            origin="lower",
            aspect="auto",
            extent=[times_min.min(), times_min.max(), lambda_grid_km.min(), lambda_grid_km.max()],
            cmap="turbo",
            norm=PowerNorm(gamma=0.5, vmin=0.0, vmax=1.0),
        )
        ax.plot(times_min, lam_dom_km, color="white", lw=1.8)
        ax.set_xlabel("Time (min)")
        ax.set_ylabel("Vertical wavelength lambda_z (km)")
        ax.set_title(title)
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label("Normalized power")

    draw(axes[0], p1, l1_km, "1) Windowed FFT in z (global norm)")
    draw(axes[1], p2, l2_km, "2) Sliding-window FFT in z (max over windows, global norm)")

    im3 = axes[2].imshow(
        p3,
        origin="lower",
        aspect="auto",
        extent=[times_min.min(), times_min.max(), lambda_grid_km.min(), lambda_grid_km.max()],
        cmap="turbo",
        norm=PowerNorm(gamma=0.5, vmin=0.0, vmax=1.0),
    )
    axes[2].plot(times_min, l3_km, color="white", lw=1.8, label="dominant lambda_z (raw wavelet)")
    axes[2].set_xlabel("Time (min)")
    axes[2].set_ylabel("Vertical wavelength lambda_z (km)")
    axes[2].set_title("3) Preprocessed raw w' wavelet; density weighting post dominant-lambda")

    # overlay weighted dominant-mode strength on secondary axis
    axr = axes[2].twinx()
    axr.plot(times_min, s3, color="black", lw=1.2, ls="--", alpha=0.9)
    axr.set_ylabel("Weighted dominant-mode strength (0-1)")
    axr.set_ylim(0.0, 1.05)

    cbar3 = fig.colorbar(im3, ax=axes[2])
    cbar3.set_label("Normalized power")

    fig.savefig(out_png, dpi=170)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    data = load_midx_profiles(args.indir)
    idx_t = np.where(data.times >= args.tmin)[0]
    if idx_t.size == 0:
        raise RuntimeError(f"No snapshots satisfy t >= {args.tmin} s")

    lambda_grid_m = np.linspace(args.lambda_min_km * 1e3, args.lambda_max_km * 1e3, args.n_lambda)
    times_min = data.times[idx_t] / 60.0

    profiles_proc = preprocess_profiles(
        data,
        lambda_expected_km=args.lambda_expected_km,
        env_mult=args.env_mult,
        env_floor_frac=args.env_floor_frac,
    )

    panel1 = panel_fft_column(profiles_proc, data.z, idx_t, lambda_grid_m, pad_factor=args.fft_pad)
    panel2 = panel_fft_sliding(
        profiles_proc,
        data.z,
        idx_t,
        lambda_grid_m,
        pad_factor=args.fft_pad,
        window_km=args.window_km,
        step_km=args.step_km,
    )
    panel3, lam3, strength3 = panel_wavelet_raw_then_weight(
        data,
        profiles_proc,
        idx_t,
        lambda_grid_m,
        morlet_w=args.morlet_w,
    )

    lam1 = dominant_lambda(panel1, lambda_grid_m)
    lam2 = dominant_lambda(panel2, lambda_grid_m)

    out_png = os.path.join(args.outdir, "dispersion_fig6b_python_three_methods.png")
    save_three_panel_figure(
        times_min,
        lambda_grid_m / 1e3,
        panel1,
        panel2,
        panel3,
        lam1 / 1e3,
        lam2 / 1e3,
        lam3 / 1e3,
        strength3,
        out_png,
    )

    print(f"Saved figure: {out_png}")
    print(f"Snapshots analyzed: {len(data.times)}, with t>={args.tmin}s: {len(idx_t)}")
    print(
        "Preprocess: time-mean subtraction + envelope normalization "
        f"(lambda_expected={args.lambda_expected_km:.2f} km, env_mult={args.env_mult:.2f}, "
        f"env_floor_frac={args.env_floor_frac:.3f})"
    )


if __name__ == "__main__":
    main()
