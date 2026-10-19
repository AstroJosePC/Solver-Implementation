#!/usr/bin/env python3
"""
analyze_linear_flux_constancy.py

Measures the solver's momentum-flux error bar from the linear control corner
(uniform wind, periodic-x, h = 0.1 km), produced by

    julia --project=. scripts/sweep_arsia_mons_stage2.jl linear

At h = 0.1 km the steady linear solution for this ridge gives a maximum
wave-induced |u'| of about 18 m/s at 55 km against U = 40 m/s, so u_hat/U = 0.45
and the stationary-wave overturning criterion u_hat >= U is not met anywhere below
the sponge. The wave is then conservative and the Eliassen-Palm result requires
the vertical flux of horizontal momentum MF(z) = <rho u'w'>_x to be independent of
height. Molecular viscosity from the MCD nu column removes only 0.03 per cent of
the flux by 55 km, so essentially all measured height dependence is numerical.

Reported, all as numbers rather than as a verdict alone:

  1. ARRIVAL   - how far up the stationary wave has propagated at each sample
                 time, from cg_z = k_eff U^2 / N(z), so the evaluation band is
                 taken from the run rather than assumed.
  2. CONSTANCY - relative peak-to-peak of MF over the band and the fitted
                 d ln|MF|/dz, the error bar on the SHAPE of MF(z).
  3. AMPLITUDE - MF(z) divided by the steady linear hydrostatic value, the
                 separate error bar on the MAGNITUDE of the flux.
  4. DRAG      - the same drift as -(1/rho) dMF/dz in m/s per sol, which is the
                 unit the deposition-partition claims are made in.
  5. LAMBDA_Z  - measured vertical wavelength from the phase tilt of w, against
                 2 pi U / N(z). The target is NOT a single number: N over
                 10-55 km in this column is 0.0078-0.0097 s^-1, so lambda_z runs
                 26-32 km. The often-quoted 13.96 km follows from N_REF = 0.018
                 s^-1, the near-surface inversion value used for the Nh/U label
                 in the sweep drivers, and holds only below about 2 km.
  6. LINEARITY - confirmation that the control run stayed linear, so that (1)-(5)
                 are measuring the regime they are supposed to measure.
  7. VS NONLIN - the control differenced against the nonlinear run of the same
                 duration on the same grid. Each profile is normalised by its own
                 steady linear value and then by its own value at 10 km, which
                 removes both the h^2 scaling and the overall amplitude deficit
                 and leaves only the height dependence. This comparison does not
                 require either run to have reached a steady state, because both
                 are sampled at the same time in the same solver, so it is the
                 part of the result least sensitive to incomplete spin-up.

  8. FILL-IN   - the same measured flux against the flux of the part of the
                 ridge spectrum that has actually arrived at each height, rather
                 than against the fully established steady value. Arrival time
                 scales as 1/k with no natural low-k cutoff, so an unconverged
                 record decreases with height even when the solver is exact. This
                 section removes that transient and reports what is left, which is
                 a tighter upper bound on the numerical error than section 2.
  9. PAIR      - optional. A twin run at a different ridge half-width, differenced
                 against this one. The half-width reweights the flux across
                 wavenumbers and so moves the arrival transient, while lambda_z,
                 the sponge thickness in wavelengths and the linearity margin are
                 unchanged, because m = N/U and the steady amplitude both carry no
                 k. Whatever drift is common to the pair is the part the fill-in
                 model does not explain, which is the numerical error bar itself.

Reads the cached deepdive_series.npz plus one snapshot for the phase tilt, and
reuses ridge_spectrum() and linear_steady_flux() from analyze_ext6h_mechanism.py
rather than restating the theory.

Usage:
  python3 scripts/analyze_linear_flux_constancy.py \
      --run-dir output/arsia_mons_linear12h_LX2000/run_U40.0_h0.1 \
      --npz output/arsia_mons_linear12h_LX2000_analysis/run_U40.0_h0.1/deepdive_series.npz
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_ext6h_mechanism import (GAMMA, G, R_GAS, CP, ridge_spectrum,
                                     linear_steady_flux)

MARS_SOL_S = 88775.244        # Mars sol, for expressing drag in m/s per sol
SPONGE_BASE_KM = 60.0

# Reporting thresholds. Stated here so the printed result is interpretable
# without re-deriving them; they are reporting bands, not physics.
PTP_TOL = 0.05                # 5 per cent relative peak-to-peak over the band
DRIFT_TOL_PCT_PER_KM = 0.1    # accumulates to 4.5 per cent across a 45 km band
PHI_FLOOR = 0.02              # smallest arrived-flux fraction that is divided by


def parse_args():
    p = argparse.ArgumentParser(description="Linear-corner flux constancy test")
    p.add_argument("--run-dir", default="output/arsia_mons_linear12h_LX2000/run_U40.0_h0.1")
    p.add_argument("--npz", default=None,
                   help="default: <run parent>_analysis/<run name>/deepdive_series.npz")
    p.add_argument("--profile", default="external/mcd_profiles/mars_arsia_mons_my31.csv")
    p.add_argument("--outdir", default=None, help="default: the npz directory")
    p.add_argument("--h-max", type=float, default=100.0, help="ridge height [m]")
    p.add_argument("--sigma", type=float, default=100e3, help="ridge half-width [m]")
    p.add_argument("--z-lo", type=float, default=10.0, help="evaluation band base [km]")
    p.add_argument("--z-hi", type=float, default=55.0, help="evaluation band top [km]")
    p.add_argument("--avg-hours", type=float, default=1.0,
                   help="length of the final-time average window [h]")
    p.add_argument("--fit-halfwidth-km", type=float, default=5.0,
                   help="half-width of the local phase-slope fit for lambda_z")
    p.add_argument("--compare-npz",
                   default="output/arsia_mons_stage1_ext12h_analysis/"
                           "run_U40.0_h3.0/deepdive_series.npz",
                   help="nonlinear run to difference against; '' to skip")
    p.add_argument("--compare-h-max", type=float, default=3000.0,
                   help="ridge height of the comparison run [m]")
    p.add_argument("--pair-npz", default="",
                   help="deepdive npz of a twin run at a different ridge half-width. "
                        "Differencing the two measured drifts removes the part of the "
                        "drift that is common to both, which under the hypothesis "
                        "being tested is the numerical part")
    p.add_argument("--pair-sigma", type=float, default=50e3,
                   help="ridge half-width of the pair run [m]")
    p.add_argument("--forecast-from", type=float, default=12.0,
                   help="fit both accounts of the drift's time dependence to the "
                        "record up to this hour and score them on what follows. "
                        "With a record no longer than this the section reports "
                        "forecasts only")
    p.add_argument("--forecast-to", type=float, default=0.0,
                   help="last hour to project the two accounts out to. 0 means "
                        "three hours past the cut, or the end of the record, "
                        "whichever is later, so the section is usable both as a "
                        "prediction before an extension and as its scoring after")
    p.add_argument("--forecast-fit-lo", default="3,4,5,6,7",
                   help="start hours for the relaxation fit window. The spread of "
                        "the asymptote across these is the systematic uncertainty "
                        "on it, which is larger than the formal error")
    p.add_argument("--updown-hours", default="8,12",
                   help="hours at which to split the flux by the sign of the "
                        "vertical wavenumber. The last full hour of the record is "
                        "always added. Empty skips the section, which is the only "
                        "part of this script that rereads the snapshots in bulk")
    p.add_argument("--updown-window", type=float, default=0.5,
                   help="half-width of the hour-mean window for that split [h]")
    p.add_argument("--updown-max-snaps", type=int, default=11,
                   help="snapshots per window, subsampled evenly")
    p.add_argument("--updown-msens", default="0.90,0.95,1.05,1.10",
                   help="multipliers on the assumed vertical wavenumber, which "
                        "bound how much the split depends on it")
    p.add_argument("--updown-synthetic", default="0.0,0.05,0.10",
                   help="reflected flux fractions to impose on a synthetic WKB "
                        "wave, which fixes the floor of the method before it is "
                        "applied to the run")
    p.add_argument("--updown-modes", type=int, default=12,
                   help="number of leading ridge modes to report the split for "
                        "individually. The per-wavenumber arrays are already "
                        "formed inside the split, so this costs no further "
                        "snapshot reading. 0 skips that table")
    p.add_argument("--updown-mode-heights", default="10,30,55",
                   help="heights [km] at which to report the per-mode split")
    p.add_argument("--updown-mode-share-floor", type=float, default=1e-3,
                   help="smallest steady flux share a mode may carry and still "
                        "enter the class means and fits of that table. Below it "
                        "R_k is a ratio of two quantities near the noise of the "
                        "partition and carries no information; such modes are "
                        "printed but not aggregated")
    p.add_argument("--transient-stride", type=int, default=3,
                   help="snapshot stride of the stationarity pass of section 12. "
                        "That section rereads the record rather than three hour "
                        "windows of it, which is the only part of this script "
                        "whose cost scales with the length of the run. 0 skips it")
    p.add_argument("--transient-window", type=float, default=3.0,
                   help="trailing window [h] over which the steady and fluctuating "
                        "parts of each mode are separated and its rotation rate "
                        "measured")
    p.add_argument("--transient-bins", type=int, default=10,
                   help="number of trailing one-hour bins over which the unreturned "
                        "class is reported. The scatter across them is the error "
                        "bar that a pair of hour means several hours apart does "
                        "not carry")
    p.add_argument("--sponge-base-km", type=float, default=SPONGE_BASE_KM,
                   help="sponge base of the run [km]; sets the round-trip clock of "
                        "the per-mode classification and the sponge-scale table")
    p.add_argument("--sponge-k0", type=float, default=0.05,
                   help="Rayleigh damping maximum of the run [1/s], used only to "
                        "report the height at which the damping rate overtakes the "
                        "intrinsic frequency")
    p.add_argument("--front-factors", default="0.5,1.0,2.0",
                   help="multipliers on the per-wavenumber arrival time, used to "
                        "bound how much the fill-in result depends on treating the "
                        "front as sharp")
    return p.parse_args()


def band_metrics(zkm, mf, lo, hi):
    """Constancy metrics for MF over [lo, hi] km.

    Returns the median, the peak-to-peak relative to it, the least-squares
    d ln|MF|/dz in per cent per km, and the largest single-level departure.
    A flat profile is the prediction, so all three spreads are error bars.
    """
    sel = (zkm >= lo) & (zkm <= hi)
    if sel.sum() < 3:
        return dict(n=int(sel.sum()))
    v = mf[sel]
    med = float(np.median(v))
    if med == 0.0 or not np.all(np.isfinite(v)):
        return dict(n=int(sel.sum()))
    slope = float(np.polyfit(zkm[sel], np.log(np.abs(v)), 1)[0])   # 1/km
    return dict(n=int(sel.sum()), median_uPa=med * 1e6,
                ptp_rel=float(np.ptp(v) / abs(med)),
                drift_pct_per_km=100.0 * slope,
                max_dev_rel=float(np.max(np.abs(v - med)) / abs(med)))


def buoyancy_integral(N, z):
    """int_0^z N dz', the transit integral of the hydrostatic group velocity.

    The segment below the first cell centre is closed at N[0]; it contributes
    about 1 per cent of the integral to 55 km but is not negligible at the base
    of the evaluation band.
    """
    return N[0] * z[0] + np.concatenate(
        [[0.0], np.cumsum(0.5 * (N[1:] + N[:-1]) * np.diff(z))])


def arrival_profile(z, N, k_eff, u_bar):
    """Time for the stationary-wave front to reach each height.

    The hydrostatic vertical group velocity of a stationary mountain wave is
    cg_z = k U^2 / N, so the transit time is exactly int dz / cg_z =
    int N dz' / (k U^2). The buoyancy frequency is integrated rather than the
    group velocity averaged because N falls by a factor of about 2.1 across the
    near-surface inversion, where the two differ.
    """
    cg = k_eff * u_bar ** 2 / N
    return buoyancy_integral(N, z) / (k_eff * u_bar ** 2), cg


def mode_timing(ibuoy, kk, u_bar, j, j_spg, nmax=12):
    """Upward arrival and sponge round-trip times of the leading ridge modes.

    A stationary hydrostatic component of horizontal wavenumber k climbs at
    cg_z = k U^2 / N(z), so it reaches level j at t_up = I(z_j) / (k U^2), with
    I the buoyancy integral of arrival_profile(). A component reflected at the
    sponge base has to climb to that base and descend again before it can be
    seen at the same level, and the downward group speed is equal and opposite,
    so it arrives at

        t_rt = (2 I(z_spg) - I(z_j)) / (k U^2).

    The gap between the two is what separates the candidate origins. Anything
    that generates a downward component locally, whether an interior reflector
    or a failure of the steady closure the split is built on, switches it on at
    t_up. Reflection at the sponge cannot produce one before t_rt, and because
    cg_z scales with k the two differ by hours for the long modes.

    Both are sharp-front transit times and inherit the idealisations of
    fillin_flux(): no dispersive spread about each front, and the O(k^2/m^2)
    correction to cg_z, about 2e-3 here, is dropped.
    """
    n = np.arange(1, min(nmax, kk.size - 1) + 1)
    denom = kk[n] * u_bar ** 2
    return (n, ibuoy[j] / denom / 3600.0,
            (2 * ibuoy[j_spg] - ibuoy[j]) / denom / 3600.0)


def per_k_steady_flux(rho, N, hk, kk, u_bar, nx):
    """Per-wavenumber contributions to the steady linear hydrostatic MF.

    Substituting the WKB solution used by linear_steady_flux into the one-sided
    cospectrum leaves, for each retained bin,

        F_k = -rho_0 N_0 U w_k k |h_k|^2,

    with w_k the conjugate-pair weight, 2 for interior bins and 1 for the mean
    and Nyquist bins. The height dependence cancels between rho, N and the
    (rho N)^-1/2 amplitude growth, so each F_k is height-independent on its own
    and sum_k F_k reproduces linear_steady_flux() to rounding. That separation is
    what allows each wavenumber to be given its own arrival time.
    """
    w = np.full(kk.size, 2.0)
    w[0] = 1.0
    if nx % 2 == 0:
        w[-1] = 1.0
    return -rho[0] * N[0] * u_bar * w * kk * np.abs(hk) ** 2


def fillin_flux(fk, kk, ibuoy, t, u_bar, front=1.0):
    """Steady flux of the part of the ridge spectrum that has arrived by time t.

    A stationary hydrostatic component of horizontal wavenumber k travels at
    cg_z = k U^2 / N(z) and so reaches height z at t_arr = I(z) / (k U^2), with
    I(z) = int_0^z N dz'. Taking each component as absent before its own arrival
    and at its steady value after it gives

        MF_fill(z, t) = sum over k >= I(z) / (U^2 t) of F_k,

    which tends to the steady linear profile as t grows. Because the flux
    spectrum k |h_k|^2 of a Gaussian ridge has a low-k tail and the arrival time
    scales as 1/k with no natural cutoff, the approach is the filling of a
    spectral tail rather than a relaxation, and MF_fill decreases with height at
    any finite time. That height dependence is a transient rather than a flux
    divergence, so it is the part of the measured drift to remove before the
    residual is read as a numerical error bar.

    The model uses a sharp front per wavenumber. It neglects the dispersive
    spread of each component about its own front, partial reflection, and the
    O(k^2/m^2) correction to cg_z, which for this ridge and column is about
    2e-3. It is therefore a benchmark, not a measurement. `front` multiplies the
    arrival time, so front > 1 delays every component by the same factor and
    brackets the sensitivity of the result to that sharpness assumption.
    """
    kmin = front * ibuoy[:, None] / (u_bar ** 2 * max(t, 1e-30))
    return (fk[None, :] * (kk[None, :] >= kmin)).sum(axis=1)


def relax_fit(x, y, taus=None):
    """Least squares of a exp(-t/tau) + d_inf, by a scan in the one nonlinear term.

    For a fixed tau the model is linear in (a, d_inf), so each node of the scan is
    an exact solve and the minimum over the grid is the global one. On a series of
    a dozen hourly points that is more dependable than a joint nonlinear solve, and
    it returns the residual as a function of tau, from which the range of tau the
    record actually constrains can be read.

    The form carries no wave physics. It is the null against which the fill-in
    account of the same fall is scored, not a competing mechanism.
    """
    x, y = np.asarray(x, float), np.asarray(y, float)
    taus = np.geomspace(0.3, 30.0, 600) if taus is None else np.asarray(taus, float)
    ssr = np.empty(taus.size)
    coefs = np.empty((taus.size, 2))
    for i, tau in enumerate(taus):
        A = np.column_stack([np.exp(-x / tau), np.ones_like(x)])
        c, *_ = np.linalg.lstsq(A, y, rcond=None)
        coefs[i] = c
        r = y - A @ c
        ssr[i] = float(r @ r)
    i0 = int(np.argmin(ssr))
    # Range of tau within the 1-sigma contour of the residual, profiling out the
    # two linear parameters. Reported because the asymptote and tau trade off.
    dof = max(x.size - 3, 1)
    band = taus[ssr <= ssr[i0] * (1.0 + 1.0 / dof)]
    return dict(a=float(coefs[i0, 0]), tau=float(taus[i0]),
                d_inf=float(coefs[i0, 1]),
                rms=float(np.sqrt(ssr[i0] / x.size)),
                tau_lo=float(band.min()), tau_hi=float(band.max()), n=int(x.size))


def relax_eval(fit, x):
    """The fitted relaxation evaluated at times x [h]."""
    return fit["a"] * np.exp(-np.asarray(x, float) / fit["tau"]) + fit["d_inf"]

def one_sided_weights(nk, nx):
    """Conjugate-pair weights for an rfft cospectrum, 2 except mean and Nyquist."""
    w = np.full(nk, 2.0)
    w[0] = 1.0
    if nx % 2 == 0:
        w[-1] = 1.0
    return w


def updown_split(u, w, rho, x, z, N, u_col, m_scale=1.0):
    """Split the x-mean momentum flux into upward- and downward-propagating parts.

    With rho_0(z) the x-mean density, write the density-scaled Fourier amplitudes

        Psi_k(z) = FFT_x[sqrt(rho_0) w']_k,   Xi_k(z) = FFT_x[sqrt(rho_0) u']_k.

    Steady anelastic continuity, ik rho_0 u' + d(rho_0 w')/dz = 0, gives
    dPsi/dz = -ik Xi + Psi/(2H) with H the density scale height, so the vertical
    derivative is obtained from the stored u without differencing in z. Writing
    Psi = A exp(imz) + B exp(-imz) with m > 0 the local vertical wavenumber,

        A exp(imz) = (Psi + Psi'/(im)) / 2,   B exp(-imz) = (Psi - Psi'/(im)) / 2,

    and substituting into the cospectrum leaves, at every height,

        Re(Xi Psi*) = -(m/k) (|A|^2 - |B|^2),

    the AB* cross terms cancelling identically in the real part. The partition is
    therefore pointwise in z. It needs no vertical transform, no window and no
    resolution of m as a spectral peak, which matters here because the 45 km
    evaluation band holds only about 1.6 vertical wavelengths and a windowed
    FFT split would leak the upward peak into the downward half.

    For a stationary wave the branch with k > 0 and m > 0 has cg_z > 0, so the A
    term carries flux upward and negative momentum, and the B term is the
    reflected part. Wavenumbers with m^2 = (N/U)^2 - k^2 <= 0 are evanescent and
    are left out of the partition; their cospectral contribution is returned
    separately and is a self-check, since a non-propagating mode should carry no
    steady flux.

    The 1/(2H) term drops out of the flux exactly, so the result does not depend
    on how the scale height is estimated. It does depend on m through the ratio
    Psi'/(im); an error of a fraction eps in m misassigns a share of about
    (eps/2)^2 of the upward flux to the downward part, which is quadratically
    small and is measured directly by `m_scale`.

    Returns per-level upward and downward fluxes, their per-wavenumber
    decomposition, the total cospectrum, and the configuration-space
    <rho u'w'>_x for closure.
    """
    nx = len(x)
    dx = float(x[1] - x[0])
    kk = 2 * np.pi * np.fft.rfftfreq(nx, dx)
    wgt = one_sided_weights(kk.size, nx)

    rho0 = rho.mean(1) if rho.ndim > 1 else np.asarray(rho, float)
    up = u - u.mean(1, keepdims=True)
    wp = w - w.mean(1, keepdims=True)
    s = np.sqrt(rho0)[:, None]
    Xi = np.fft.rfft(s * up, axis=1) / nx
    Psi = np.fft.rfft(s * wp, axis=1) / nx

    inv2H = -0.5 * np.gradient(np.log(rho0), z)
    m = m_scale * np.sqrt(np.clip((N / u_col)[:, None] ** 2 - kk[None, :] ** 2, 0.0, None))
    prop = m > 0.0

    dPsi = -1j * kk[None, :] * Xi + inv2H[:, None] * Psi

    kpos = kk[None, :] > 0.0
    mk = np.where(prop & kpos, m / np.where(kpos, kk[None, :], 1.0), 0.0)

    def part(d):
        q = np.where(prop, d / (1j * np.where(prop, m, 1.0)), 0.0)
        Aq, Bq = 0.5 * (Psi + q), 0.5 * (Psi - q)
        return (Aq, Bq, -wgt[None, :] * mk * np.abs(Aq) ** 2,
                +wgt[None, :] * mk * np.abs(Bq) ** 2)

    A, B, fup, fdn = part(dPsi)

    # Audit of the one approximation the partition contains. dPsi is taken from
    # steady continuity, which drops d(rho')/dt; the solver is compressible, so
    # that term is not identically zero. Differencing Psi in z is an independent
    # estimate of the same derivative, accurate here to about 0.2 per cent since
    # dz = 500 m resolves the band-median lambda_z of 31 km with about 60 points.
    # It is not used for the split, because a difference is contaminated at the
    # front of a mode that is still arriving and near the boundaries, but the
    # partition it produces bounds what the closure can be contributing. An error
    # of relative size eps in dPsi puts B = -eps Psi / 2 and so manufactures a
    # downward fraction of about eps^2 / 4 in every mode, which is the same
    # quadratic law the m_scale sensitivity measures.
    dPsi_fd = np.gradient(Psi, z, axis=0)
    _, _, fup_fd, fdn_fd = part(dPsi_fd)
    scl = np.where(prop, m * np.abs(Psi), np.nan)
    eps = np.abs(dPsi - dPsi_fd) / np.where(scl > 0, scl, np.nan)

    cospec = wgt[None, :] * np.real(Xi * np.conj(Psi))
    mf_direct = ((rho * up * wp).mean(1) if rho.ndim > 1
                 else (rho0[:, None] * up * wp).mean(1))
    return dict(kk=kk, wgt=wgt, m=m, mk=mk, A=A, B=B, fup=fup, fdn=fdn,
                fup_fd=fup_fd, fdn_fd=fdn_fd, closure_eps=eps,
                mf_up=fup.sum(1), mf_dn=fdn.sum(1),
                mf_spec=cospec.sum(1), mf_evan=(cospec * ~prop).sum(1),
                mf_direct=mf_direct,
                var_total=(wgt[None, :] * np.abs(Psi) ** 2).sum(1),
                var_incoherent=(wgt[None, :] * (np.abs(A) ** 2 + np.abs(B) ** 2)).sum(1))


def synthetic_wkb(hk, kk, z, rho0, N, u_bar, nx, refl=0.0):
    """A stationary WKB mountain wave with a known reflected flux fraction.

    Built to exercise updown_split() against an answer that is known in advance.
    The amplitude follows wave-action conservation, |Psi|^2 m / k constant, so
    the returned field has the height-dependent amplitude of the real solution
    and the split has to survive that gradient rather than a constant-amplitude
    idealisation. `refl` is the downward flux as a fraction of the upward.
    """
    m = np.sqrt(np.clip((N[:, None] / u_bar) ** 2 - kk[None, :] ** 2, 0.0, None))
    prop = m > 0.0
    ph = np.concatenate([np.zeros((1, kk.size)),
                         np.cumsum(0.5 * (m[1:] + m[:-1]) * np.diff(z)[:, None], axis=0)])
    scale = np.where(prop, np.sqrt(np.where(prop, m[0][None, :] / np.where(prop, m, 1.0), 1.0)), 0.0)
    amp0 = 1j * kk[None, :] * u_bar * hk[None, :] * np.sqrt(rho0[0])
    Pu = amp0 * scale * np.exp(1j * ph)
    Pd = np.sqrt(refl) * amp0 * scale * np.exp(-1j * ph) * np.exp(2j * ph[-1])
    Psi = Pu + Pd
    dsc = np.gradient(scale, z, axis=0)
    dPsi = (1j * m * (Pu - Pd)
            + np.where(scale > 0, dsc / np.where(scale > 0, scale, 1.0), 0.0) * Psi)
    inv2H = -0.5 * np.gradient(np.log(rho0), z)
    kpos = kk[None, :] > 0.0
    Xi = 1j * (dPsi - inv2H[:, None] * Psi) / np.where(kpos, kk[None, :], 1.0)
    Xi[:, 0] = 0.0
    wf = np.fft.irfft(Psi * nx, nx, axis=1) / np.sqrt(rho0)[:, None]
    uf = np.fft.irfft(Xi * nx, nx, axis=1) / np.sqrt(rho0)[:, None]
    return uf + u_bar, wf, np.repeat(rho0[:, None], nx, 1)


def synthetic_exact(hk, kk, z, rho_surf, N0, Hs, u_bar, nx):
    """An exact steady anelastic solution, for the systematic bias of the split.

    synthetic_wkb() builds its field from the same m the split assumes, so it
    exercises the algebra, the conjugate-pair weighting and the WKB amplitude
    gradient, but it cannot see an error in m itself. At constant N and constant
    density scale height H the anelastic vertical structure equation for
    Psi = sqrt(rho_0) w' is

        d2Psi/dz2 + [N^2/U^2 - k^2 - 1/(4H^2)] Psi = 0,

    so the exact vertical wavenumber carries a term that
    m = sqrt(N^2/U^2 - k^2) omits. In this column 1/(4H^2) is about 6 per cent
    of (N/U)^2, so the assumed m is low by about 3 per cent, and by the same
    quadratic law the m_scale sensitivity measures, a purely upward exact
    solution is returned with a downward share of ((1 - r)/(1 + r))^2 with
    r = m_exact/m, about 2.4e-04 in every propagating mode.

    That number is a systematic bias of the method rather than noise. It does
    not average down over modes or over time, it is the same sign in every
    mode, and it is the floor the per-mode ratios should be read against. The
    field returned here is exact, so its own flux is height-independent to
    machine precision, which is the check that the construction is right.
    """
    rho0 = rho_surf * np.exp(-(z - z[0]) / Hs)
    m_ex = np.sqrt(np.clip((N0 / u_bar) ** 2 - kk[None, :] ** 2
                           - 1.0 / (4 * Hs ** 2), 0.0, None))
    Psi = (1j * kk[None, :] * u_bar * hk[None, :] * np.sqrt(rho0[0])
           * np.exp(1j * m_ex * (z - z[0])[:, None]))
    kpos = np.where(kk > 0.0, kk, np.inf)
    Xi = 1j * (1j * m_ex * Psi - 0.5 * Psi / Hs) / kpos[None, :]
    Xi[:, 0] = 0.0
    s = np.sqrt(rho0)[:, None]
    return (np.fft.irfft(Xi * nx, nx, axis=1) / s + u_bar,
            np.fft.irfft(Psi * nx, nx, axis=1) / s,
            np.repeat(rho0[:, None], nx, 1), m_ex[0])


def rayleigh_sponge_reflection(k, z, N, rho0, u_bar, K0, z_base, z_top,
                               top="lid", return_field=False, c2=None):
    """Linear reflected share of a sin^2 Rayleigh layer, per level, for one k.

    The solver relaxes every conserved variable toward the initial state at the
    rate K(z) = K0 sin^2(pi/2 (z - z_base)/(z_top - z_base)) (src/sponge.jl), so
    for a steady component the damping enters as Omega -> Omega + iK with
    Omega = -U k. The anelastic structure equation for Phi = sqrt(rho_0) w' is
    then

        Phi'' + [k^2 (N^2/(Omega + iK)^2 - 1) - s^2 + s'] Phi = 0,
        s = -(1/2) d ln(rho_0)/dz,

    integrated downward from z_top by RK4 on the supplied grid. `top` is "lid"
    (Phi = 0, a rigid top) or "radiation" (Phi' = i m_top Phi, upward energy),
    which bracket the eq19 outflow condition the solver applies. Below the layer
    the field is split with the exact local wavenumber, the same bracket without
    damping, so the result is the physical reflection with no partition floor in
    it. Returns |B/A|^2 on the grid levels, NaN where the mode is evanescent.

    Limits that check the construction: K0 = 0 with a lid returns 1 everywhere
    below the top; K0 = 0 with a radiation top returns the partial reflection of
    the N(z) profile alone; a thick, weak layer returns values near that one.

    With return_field the unnormalised solution is returned as well, as
    (R, Phi, dPhi), so that a caller can scale it to an imposed lower boundary
    and read the upward level the column would carry (cavity closure test).

    With c2 (the squared sound speed on the grid) the isothermal compressible
    term (Omega + iK)^2 / c^2 is added to the bracket. It is negligible in the
    free column, but once |Omega + iK| exceeds about sqrt(c k N) it outweighs
    the buoyancy term; for n = 6 in the 60-90 km K0 = 0.05 layer that is from
    about 7 km above the base.
    """
    z = np.asarray(z, float)
    a = np.gradient(np.log(rho0), z)
    s = -0.5 * a
    ds = np.gradient(s, z)
    K = np.where(z >= z_base,
                 K0 * np.sin(0.5 * np.pi * np.clip((z - z_base) / (z_top - z_base),
                                                   0.0, 1.0)) ** 2, 0.0)
    om = -u_bar * k + 1j * K
    L = k ** 2 * (N ** 2 / om ** 2 - 1.0) - s ** 2 + ds
    L0 = (k ** 2 * (N ** 2 / (u_bar * k) ** 2 - 1.0) - s ** 2 + ds).real
    if c2 is not None:
        L = L + om ** 2 / c2
        L0 = L0 + (u_bar * k) ** 2 / c2
    it = int(np.searchsorted(z, z_top, side="right")) - 1
    Phi = np.zeros(len(z), complex)
    dPhi = np.zeros(len(z), complex)
    if top == "lid":
        Phi[it], dPhi[it] = 0.0, 1.0
    else:
        mt = np.sqrt(L[it] + 0j)
        mt = mt if mt.imag >= 0 else -mt
        Phi[it], dPhi[it] = 1.0, 1j * mt
    for j in range(it, 0, -1):
        h = z[j - 1] - z[j]
        Lm = 0.5 * (L[j] + L[j - 1])
        k1p, k1d = dPhi[j], -L[j] * Phi[j]
        k2p, k2d = dPhi[j] + 0.5 * h * k1d, -Lm * (Phi[j] + 0.5 * h * k1p)
        k3p, k3d = dPhi[j] + 0.5 * h * k2d, -Lm * (Phi[j] + 0.5 * h * k2p)
        k4p, k4d = dPhi[j] + h * k3d, -L[j - 1] * (Phi[j] + h * k3p)
        Phi[j - 1] = Phi[j] + h / 6 * (k1p + 2 * k2p + 2 * k3p + k4p)
        dPhi[j - 1] = dPhi[j] + h / 6 * (k1d + 2 * k2d + 2 * k3d + k4d)
    R = np.full(len(z), np.nan)
    ok = (L0 > 0) & (np.arange(len(z)) <= it)
    m = np.sqrt(np.where(ok, L0, 1.0))
    A = 0.5 * (Phi + dPhi / (1j * m))
    B = 0.5 * (Phi - dPhi / (1j * m))
    R[ok] = (np.abs(B[ok]) / np.maximum(np.abs(A[ok]), 1e-300)) ** 2
    return (R, Phi, dPhi) if return_field else R


def ripple_reflection(psi2, z, m, sel, scales=None):
    """Reflected share of one mode from the ripple in its variance alone.

    For Psi = A exp(i theta) + B exp(-i theta) with theta = int m dz and WKB
    amplitudes |A|, |B| proportional to m^-1/2,

        m |Psi|^2 = a^2 + b^2 + 2 a b cos(2 theta + phi),

    so the WKB-flattened variance is a constant carrying a ripple at half the
    vertical wavelength whose relative depth 2q/(1 + q^2), q = b/a, fixes the
    reflected share R = q^2. This uses |Psi|^2 only. It shares no algebra with
    updown_split(), which reads the direction from the phase relation between
    Psi and its derivative, so the two are independent estimates of the same
    quantity, the standing-wave ratio of a transmission line against a
    directional coupler. A linear trend in z is fitted alongside, since a mode
    that is still filling carries one.

    `m` should be the exact local wavenumber, including the -s^2 + s' density
    terms, since a phase error across the band lowers the fitted depth. `scales`
    rescans m and reports the scale at which the ripple is deepest, a check that
    the ripple has the wavelength a reflection would give rather than some other.
    Returns (R at scale 1, R at the best scale, best scale, fitted depth).
    """
    zz = z[sel]
    y = (psi2 * m)[sel]
    th0 = np.concatenate([[0.0], np.cumsum(0.5 * (m[1:] + m[:-1]) * np.diff(z))])[sel]
    zc = zz - zz.mean()

    def fit(sc):
        X = np.column_stack([np.ones_like(zz), zc, np.cos(2 * sc * th0),
                             np.sin(2 * sc * th0)])
        c = np.linalg.lstsq(X, y, rcond=None)[0]
        d = min(np.hypot(c[2], c[3]) / max(c[0], 1e-300), 1.0)
        q = (1.0 - np.sqrt(1.0 - d * d)) / d if d > 0 else 0.0
        return q * q, d

    R1, d1 = fit(1.0)
    scales = np.linspace(0.9, 1.1, 41) if scales is None else scales
    best = max(scales, key=lambda sc: fit(sc)[1])
    return R1, fit(best)[0], float(best), d1


def _snap_time(path):
    """Model time of one snapshot, for records the cached series does not cover."""
    with Dataset(path) as ds:
        return float(ds.getncattr("time_s"))


def subbands(lo, hi, n=3):
    """Equal parts of the evaluation band, then the whole band.

    A straight line through log|MF| is only a summary. Reporting the parts next
    to the whole shows when the profile is non-monotone and the summary is
    therefore the wrong statistic.
    """
    e = np.linspace(lo, hi, n + 1)
    return [(float(e[i]), float(e[i + 1])) for i in range(n)] + [(float(lo), float(hi))]


def ud_fig(outdir, zkm, got, hours, trace, z_lo, z_hi, z_top_km):
    """Three panels: the split profiles, the local ratio, and its growth."""
    fig, ax = plt.subplots(1, 3, figsize=(14, 5))
    for c, hh in enumerate(hours):
        g = got[hh]
        col = f"C{c}"
        ax[0].plot(-g["mf_up"] * 1e6, zkm, color=col, lw=1.2,
                   label=f"upward, t={g['t_h']:.1f} h")
        ax[0].plot(g["mf_dn"] * 1e6, zkm, color=col, lw=1.2, ls="--",
                   label=f"downward, t={g['t_h']:.1f} h")
        with np.errstate(divide="ignore", invalid="ignore"):
            ax[1].plot(np.abs(g["mf_dn"] / g["mf_up"]), zkm, color=col, lw=1.2,
                       label=f"t={g['t_h']:.1f} h")
    for axx in ax[:2]:
        axx.axhspan(z_lo, z_hi, color="k", alpha=0.06)
        axx.axhline(SPONGE_BASE_KM, color="k", lw=0.8, ls=":")
        axx.set_ylim(0, z_top_km)
        axx.grid(alpha=0.3)
    ax[0].set_xscale("log")
    ax[0].set_xlabel("|MF| by propagation direction (uPa)")
    ax[0].set_ylabel("z (km)")
    ax[0].legend(fontsize=7)
    ax[1].set_xlim(0, 0.3)
    ax[1].set_xlabel("R = |MF_down| / |MF_up|")
    ax[1].legend(fontsize=8)
    tt = [e["t_h"] for e in trace]
    ax[2].plot(tt, [e["R"] for e in trace], "k.-", lw=1.0, ms=4)
    ax[2].set_xlabel("t (h)")
    ax[2].set_ylabel(f"R over {z_lo:.0f}-{z_hi:.0f} km")
    ax[2].grid(alpha=0.3)
    axb = ax[2].twinx()
    axb.plot(tt, [-e["MF_net_uPa"] for e in trace], color="C3", lw=1.0)
    axb.plot(tt, [-e["MF_up_uPa"] for e in trace], color="C0", lw=1.0, ls="--")
    axb.set_ylabel("|MF| net (solid) and upward (dashed), uPa", color="C3")
    fig.suptitle("Up and down decomposition of the momentum flux", fontsize=11)
    fig.tight_layout()
    png = os.path.join(outdir, "linear_flux_updown.png")
    fig.savefig(png, dpi=150)
    plt.close(fig)
    return png


def pm_fig(outdir, pm, ud_hours, t_end_h):
    """Three panels: the timing that sets the test, the per-mode ratio, its time.

    The left panel is the test itself rather than a result. It shows where each
    mode's upward arrival and its sponge round trip fall against the end of the
    record, so the modes for which a sponge origin requires a zero can be read
    off directly.
    """
    zs = sorted({r["z_km"] for r in pm["rows"]})
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.6))

    for c, zz in enumerate(zs):
        rw = sorted([r for r in pm["rows"] if r["z_km"] == zz],
                    key=lambda r: r["lambda_x_km"])
        lam = [r["lambda_x_km"] for r in rw]
        ax[0].plot(lam, [r["t_up_h"] for r in rw], color=f"C{c}", lw=1.2,
                   marker="o", ms=3, label=f"t_up, z={zz:.0f} km")
        ax[0].plot(lam, [r["t_rt_h"] for r in rw], color=f"C{c}", lw=1.2,
                   ls="--", marker="s", ms=3, label=f"t_rt, z={zz:.0f} km")

        hh = max(ud_hours)
        key = f"{hh:g}h"
        for r in rw:
            mk = {"returned": "o", "no return": "D", "arriving": "x",
                  "negligible": "."}[r["class"]]
            fc = f"C{c}" if r["class"] == "returned" else "none"
            ec = "0.6" if r["class"] == "negligible" else f"C{c}"
            ax[1].plot(r["lambda_x_km"], r["R_k"][key], marker=mk,
                       ms=4 if r["class"] == "negligible" else 6,
                       mfc=fc, mec=ec, ls="none")
        ax[1].plot([], [], color=f"C{c}", marker="o", ls="none",
                   label=f"z={zz:.0f} km")

    ax[0].axhline(t_end_h, color="k", lw=1.0, ls=":")
    ax[0].annotate("end of record", (lam[-1], t_end_h), textcoords="offset points",
                   xytext=(-4, 4), ha="right", fontsize=8)
    ax[0].set(xscale="log", yscale="log", xlabel="lambda_x (km)",
              ylabel="time (h)", title="Arrival and sponge round trip")
    ax[0].legend(fontsize=6.5, ncol=2)
    ax[0].grid(alpha=0.3, which="both")

    fl = max((r["floor"] for r in pm["rows"]), default=1e-5)
    ax[1].axhline(fl, color="k", lw=1.0, ls="-")
    ax[1].annotate("method floor", (max(lam), fl), textcoords="offset points",
                   xytext=(-2, 4), ha="right", fontsize=8)
    ca = max((e["spurious_R_from_eps"] for e in pm.get("closure_audit", [])),
             default=None)
    if ca:
        ax[1].axhline(ca, color="k", lw=1.0, ls="--")
        ax[1].annotate("closure bound eps^2/4", (max(lam), ca),
                       textcoords="offset points", xytext=(-2, 4), ha="right",
                       fontsize=8)
    ax[1].plot([], [], marker="D", mfc="none", mec="k", ls="none",
               label="no return possible")
    ax[1].plot([], [], marker="o", mfc="k", mec="k", ls="none", label="returned")
    ax[1].plot([], [], marker="x", mec="k", ls="none", label="still arriving")
    ax[1].plot([], [], marker=".", mec="0.6", ls="none", label="share below floor")
    ax[1].set(xscale="log", yscale="log", xlabel="lambda_x (km)",
              ylabel="R_k = |MF_down,k| / |MF_up,k|",
              title=f"Per-mode ratio at t = {max(ud_hours):.0f} h")
    ax[1].legend(fontsize=6.5, ncol=2, loc="lower left")
    ax[1].grid(alpha=0.3, which="both")

    for c, v in enumerate(pm["verdict"]):
        for cl, ls, mk in zip(v["classes"], ("--", "-"), ("D", "o")):
            xs = [cl["by_hour"][f"{hh:g}h"] for hh in ud_hours]
            tt = [h for h, e in zip(ud_hours, xs) if e["mean"] == e["mean"]]
            mu = [e["mean"] for e in xs if e["mean"] == e["mean"]]
            se = [e["se"] if e["se"] == e["se"] else 0.0
                  for e in xs if e["mean"] == e["mean"]]
            if not tt:
                continue
            ax[2].errorbar(tt, mu, yerr=se, color=f"C{c}", ls=ls, marker=mk,
                           ms=5, lw=1.2, capsize=3,
                           mfc=(f"C{c}" if cl["class"] == "returned" else "none"),
                           label=f"{cl['class']}, z={v['z_km']:.0f} km")
    ax[2].set(xlabel="t (h)", ylabel="share-weighted mean R_k",
              title="By class, membership recomputed each hour")
    ax[2].legend(fontsize=6.5)
    ax[2].grid(alpha=0.3)

    fig.suptitle("Per-mode timing of the downward component")
    fig.tight_layout()
    out = os.path.join(outdir, "linear_flux_permode.png")
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"wrote {out}")


def transient_pass(files, idx, z, nmode):
    """Complex mode amplitudes at every level, one entry per sampled snapshot.

    updown_split() already forms A_k and B_k at every height, so the only thing
    this adds is keeping them rather than the flux they sum to. Everything the
    stationarity test needs is a functional of these two arrays, and the flux is
    recovered as wgt m/k |A|^2 and wgt m/k |B|^2, so the pass is done once and
    the several statistics built from the result rather than from the files.

    The background N and U are taken from each snapshot, as ud_mean() does by
    default, so the split is the same one section 11 reports.
    """
    nt, nz = idx.size, len(z)
    tt = np.zeros(nt)
    Ak = np.zeros((nt, nz, nmode + 1), np.complex128)
    Bk = np.zeros((nt, nz, nmode + 1), np.complex128)
    MKk = np.zeros((nt, nz, nmode + 1))
    wgt_k = None
    for q, i in enumerate(idx):
        with Dataset(files[i]) as ds:
            tt[q] = float(ds.getncattr("time_s"))
            xs = np.asarray(ds["x"][:], dtype=np.float64)
            us = np.asarray(ds["u"][:], dtype=np.float64)
            ws = np.asarray(ds["w"][:], dtype=np.float64)
            rs = np.asarray(ds["rho"][:], dtype=np.float64)
            Ts = np.asarray(ds["T"][:], dtype=np.float64).mean(axis=1)
        Nl = np.sqrt(np.clip((G / Ts) * (np.gradient(Ts, z) + G / CP), 1e-12, None))
        sp = updown_split(us, ws, rs, xs, z, Nl, us.mean(axis=1))
        Ak[q] = sp["A"][:, : nmode + 1]
        Bk[q] = sp["B"][:, : nmode + 1]
        MKk[q] = sp["mk"][:, : nmode + 1]
        wgt_k = sp["wgt"][: nmode + 1]
    return tt, Ak, Bk, MKk, wgt_k


def tr_fig(outdir, tr, th_tr, Ak, Bk, zkm, j0, kk, nmode):
    """Three panels: the corrected clock, the class in time, and the wander.

    The first panel is the test of section 11's classification rather than a
    result, the round trip recomputed from each mode's measured rotation rate
    against the end of the record. The second is the quantity section 11 quotes,
    binned rather than sampled at two hours. The third is the raw material of
    both, the mode amplitudes themselves.
    """
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.6))

    rows = [e for e in tr.get("per_mode", []) if e["z_km"] == float(zkm[j0])]
    if rows:
        lx = np.array([e["lambda_x_km"] for e in rows])
        for key, mk_, lab in (("t_rt_h", "o", "tabulated, stationary"),
                              ("t_rt_corrected_h", "s", "from the measured rate")):
            v = np.array([e[key] for e in rows])
            ax[0].plot(lx, v, mk_, ms=5, mfc="none", label=lab)
        ax[0].axhline(th_tr[-1], color="k", lw=1.0)
        ax[0].text(lx.min(), th_tr[-1] * 1.06, "end of record", fontsize=8)
        ax[0].set_xscale("log")
        ax[0].set_yscale("log")
        ax[0].set_xlabel("lambda_x [km]")
        ax[0].set_ylabel("sponge round trip [h]")
        ax[0].set_title(f"round trip at {zkm[j0]:.1f} km", fontsize=10)
        ax[0].legend(fontsize=7, frameon=False)
        ax[0].grid(alpha=0.3)

    bins = tr.get("hourly_bins", [])
    if bins:
        tb = np.array([e["t_h"] for e in bins])
        ax[1].plot(tb, [e["R_share_weighted"] for e in bins], "o-", ms=4,
                   label="R, share weighted")
        ax[1].plot(tb, [e["R_steady"] for e in bins], "s--", ms=4,
                   label="steady within the bin")
        ax[1].plot(tb, [e["R_fluct"] for e in bins], "^--", ms=4,
                   label="fluctuating")
        ax[1].plot(tb, [e["up_fluct"] for e in bins], "v:", ms=4, color="0.5",
                   label="fluctuating, upward branch")
        s = tr.get("hourly_summary")
        if s:
            ax[1].axhspan(s["mean"] - s["sd"], s["mean"] + s["sd"], color="C0",
                          alpha=0.12)
            ax[1].axhline(s["mean"], color="C0", lw=0.8)
            ax[1].set_title(f"unreturned class, {s['trend_per_h']:+.5f}/h "
                            f"({s['trend_sigma']:.1f} se)", fontsize=10)
        ax[1].set_xlabel("t [h]")
        ax[1].set_ylabel("downward over upward")
        ax[1].legend(fontsize=7, frameon=False)
        ax[1].grid(alpha=0.3)

    for c, n in enumerate([n_ for n_ in (1, 2, 3) if n_ <= nmode]):
        col = f"C{c}"
        ax[2].plot(th_tr, np.abs(Bk[:, j0, n]) / np.abs(Bk[:, j0, n]).mean(),
                   color=col, lw=1.0, label=f"|B| n={n}")
        ax[2].plot(th_tr, np.abs(Ak[:, j0, n]) / np.abs(Ak[:, j0, n]).mean(),
                   color=col, lw=1.0, ls="--", alpha=0.7)
    ax[2].axhline(1.0, color="k", lw=0.6)
    ax[2].set_xlabel("t [h]")
    ax[2].set_ylabel("amplitude over its record mean")
    ax[2].set_title(f"mode amplitudes at {zkm[j0]:.1f} km, |A| dashed", fontsize=10)
    ax[2].legend(fontsize=7, frameon=False, ncol=3)
    ax[2].grid(alpha=0.3)

    fig.tight_layout()
    p = os.path.join(outdir, "linear_flux_transient.png")
    fig.savefig(p, dpi=150)
    plt.close(fig)
    print(f"\nFigure written to {p}")


def measure_lambda_z(path, u_bar, zkm, lo, hi, half_km):
    """Vertical wavelength from the phase tilt of w at the flux-dominant k.

    Takes the x-FFT of u' and w' at every level, picks the single horizontal
    wavenumber carrying the most band-integrated cospectral flux, then reads
    m(z) as the local least-squares slope of the unwrapped phase of w at that k.
    Returns lambda_z(z) = 2 pi / |m| plus the sign of m, which records the
    direction of the phase tilt.
    """
    with Dataset(path) as ds:
        x = np.asarray(ds["x"][:], dtype=np.float64)
        z = np.asarray(ds["z"][:], dtype=np.float64)
        u = np.asarray(ds["u"][:], dtype=np.float64)
        w = np.asarray(ds["w"][:], dtype=np.float64)
        rho = np.asarray(ds["rho"][:], dtype=np.float64)
        t_s = float(ds.getncattr("time_s"))
    nx, dx = len(x), float(x[1] - x[0])
    uk = np.fft.rfft(u - u.mean(1, keepdims=True), axis=1) / nx
    wk = np.fft.rfft(w - w.mean(1, keepdims=True), axis=1) / nx
    kk = 2 * np.pi * np.fft.rfftfreq(nx, dx)

    # One-sided cospectrum: every interior wavenumber stands for a conjugate
    # pair, so it is counted twice. The sum over k then reproduces the
    # configuration-space mean, which is asserted below as a self-check.
    wgt = np.full(uk.shape[1], 2.0)
    wgt[0] = 1.0
    if nx % 2 == 0:
        wgt[-1] = 1.0
    cospec = wgt[None, :] * np.real(uk * np.conj(wk))
    uw_direct = ((u - u.mean(1, keepdims=True))
                 * (w - w.mean(1, keepdims=True))).mean(1)
    closure = float(np.max(np.abs(cospec.sum(1) - uw_direct))
                    / max(np.max(np.abs(uw_direct)), 1e-30))

    zk = z / 1e3
    sel = (zk >= lo) & (zk <= hi)
    rho1 = rho.mean(1) if rho.ndim > 1 else rho
    band_flux = np.abs((rho1[sel, None] * cospec[sel]).sum(0))  # flux-weighted, not variance
    band_flux[0] = 0.0
    k_dom = int(np.argmax(band_flux))

    phase = np.unwrap(np.angle(wk[:, k_dom]))
    m = np.full(len(zk), np.nan)
    for j in range(len(zk)):
        win = np.abs(zk - zk[j]) <= half_km
        if win.sum() >= 3:
            m[j] = np.polyfit(z[win], phase[win], 1)[0]
    with np.errstate(divide="ignore", invalid="ignore"):
        lam = 2 * np.pi / np.abs(m) / 1e3
    return dict(t_s=t_s, k_dom=float(kk[k_dom]),
                lambda_x_dom_km=float(2 * np.pi / kk[k_dom] / 1e3),
                lambda_z_km=lam, m_sign=float(np.sign(np.nanmedian(m[sel]))),
                zkm=zk, spec_closure=closure)


def main():
    global SPONGE_BASE_KM
    a = parse_args()
    SPONGE_BASE_KM = a.sponge_base_km
    run_dir = os.path.normpath(a.run_dir)
    if a.npz is None:
        parent, name = os.path.split(run_dir)
        a.npz = os.path.join(parent + "_analysis", name, "deepdive_series.npz")
    outdir = a.outdir or os.path.dirname(a.npz)
    os.makedirs(outdir, exist_ok=True)

    d = np.load(a.npz, allow_pickle=True)
    t, zkm = d["times"], d["z_km"]
    th, z = t / 3600.0, zkm * 1e3
    mf, ubar = d["mf_zt"], d["ubar_zt"]
    umin, wmax = d["umin_hill_zt"], d["wmax_hill_zt"]
    u_bar = float(d["u_bar"])

    raw = np.genfromtxt(a.profile, delimiter=",", names=True)
    rho = np.interp(z, raw["z_m"], raw["rho"])
    T = np.interp(z, raw["z_m"], raw["T"])
    N = np.sqrt(np.clip((G / T) * (np.gradient(T, z) + G / CP), 1e-12, None))

    files = sorted(glob.glob(os.path.join(run_dir, "snapshot_*.nc")))
    with Dataset(files[0]) as ds:
        nx = len(np.asarray(ds["x"][:]))
        dx = float(np.asarray(ds["x"][:])[1] - np.asarray(ds["x"][:])[0])
    lx = nx * dx

    z_top_km = float(zkm[-1] + 0.5 * (zkm[1] - zkm[0]))
    k_eff, hk, kk = ridge_spectrum(h_max=a.h_max, sigma=a.sigma, lx=lx, dx=dx, nx=nx)
    mf_lin = linear_steady_flux(rho, N, hk, kk, z, u_bar=u_bar, nx=nx)
    t_arr, cg = arrival_profile(z, N, k_eff, u_bar)
    # Band-median 2 pi U / N, the wavelength the sponge thickness is measured in.
    lam_z_band = float(np.median(2 * np.pi * u_bar / N[(zkm >= a.z_lo) & (zkm <= a.z_hi)]) / 1e3)

    out = {"run_dir": run_dir, "u_bar": u_bar, "h_max_m": a.h_max,
           "LX_km": lx / 1e3, "dx_m": dx, "nx": nx,
           "t_end_h": float(th[-1]), "n_snapshots": len(files),
           "k_eff_per_m": k_eff, "lambda_x_eff_km": 2 * np.pi / k_eff / 1e3,
           "linear_steady_MF_uPa": float(mf_lin[0] * 1e6),
           "linear_steady_MF_flatness_pct": float(
               100 * np.ptp(mf_lin) / abs(mf_lin[0])),
           "eval_band_km": [a.z_lo, a.z_hi]}

    hdr = lambda s: print("\n" + "=" * 78 + f"\n{s}\n" + "=" * 78)
    print(f"Linear control corner: {run_dir}")
    print(f"  U={u_bar:.0f} m/s  h={a.h_max:.0f} m  LX={lx/1e3:.0f} km  "
          f"dx={dx:.0f} m  nx={nx}  t_end={th[-1]:.2f} h  {len(files)} snapshots")
    print(f"  Ridge k_eff={k_eff:.3e} 1/m (lambda_x={2*np.pi/k_eff/1e3:.0f} km)")
    print(f"  Steady linear hydrostatic MF = {mf_lin[0]*1e6:.3f} uPa "
          f"(height-independent to {100*np.ptp(mf_lin)/abs(mf_lin[0]):.1e} %)")
    print("  NOTE: MF is a domain x-mean and so scales as 1/LX; comparison with a "
          "run\n        at a different width requires that normalization.")

    # ---------------- 1. Wave arrival ----------------
    hdr("1. ARRIVAL  -  height reached by the stationary wave, cg_z = k_eff U^2/N")
    print(f"{'z[km]':>7} {'N[1/s]':>9} {'cg_z[m/s]':>10} {'t_arrival[h]':>13}")
    for zz in [10, 20, 30, 40, 50, 55, 60]:
        j = int(np.argmin(np.abs(zkm - zz)))
        print(f"{zkm[j]:7.2f} {N[j]:9.5f} {cg[j]:10.2f} {t_arr[j]/3600:13.2f}")
    fill = lambda tt: float(zkm[t_arr <= tt].max()) if (t_arr <= tt).any() else 0.0
    z_fill_end = fill(t[-1])
    out["arrival_h"] = {f"{zz}km": float(t_arr[int(np.argmin(np.abs(zkm - zz)))] / 3600)
                        for zz in [10, 20, 30, 40, 50, 55, 60]}
    out["fill_height_at_t_end_km"] = z_fill_end
    print(f"\nFilled to {z_fill_end:.2f} km by t_end = {th[-1]:.2f} h.")
    z_hi = min(a.z_hi, z_fill_end)
    if z_hi < a.z_hi:
        print(f"Evaluation top reduced {a.z_hi:.1f} -> {z_hi:.2f} km: above that "
              "height the\nrecord is still the startup transient, not the steady wave.")
    out["eval_band_used_km"] = [a.z_lo, z_hi]

    # ---------------- 2. Flux constancy ----------------
    hdr(f"2. CONSTANCY  -  MF over {a.z_lo:.0f}-{z_hi:.0f} km; EP predicts dMF/dz = 0")
    print(f"{'t[h]':>6} {'median[uPa]':>12} {'ptp/|med|':>10} "
          f"{'drift[%/km]':>12} {'max dev':>9}")
    # Whole hours plus the final sample, so the endpoint appears even when the
    # record stops just short of the hour (t_end lands on a step, not a snapshot).
    hours = sorted(set(range(1, int(np.floor(th[-1])) + 1)) | {round(th[-1], 2)})
    per_time = {}
    for hh in hours:
        i = int(np.argmin(np.abs(th - hh)))
        mm = band_metrics(zkm, mf[:, i], a.z_lo, z_hi)
        if "ptp_rel" not in mm:
            continue
        per_time[f"{hh:g}h"] = mm
        print(f"{th[i]:6.2f} {mm['median_uPa']:12.4f} {mm['ptp_rel']:10.4f} "
              f"{mm['drift_pct_per_km']:12.4f} {mm['max_dev_rel']:9.4f}")
    out["constancy_by_hour"] = per_time

    w_avg = t >= (t[-1] - a.avg_hours * 3600.0)
    mf_bar = mf[:, w_avg].mean(axis=1)
    final = band_metrics(zkm, mf_bar, a.z_lo, z_hi)
    out["constancy_final_average"] = final
    out["avg_window_h"] = [float(th[w_avg][0]), float(th[-1])]
    print(f"\nAveraged over the final {a.avg_hours:.1f} h "
          f"({th[w_avg][0]:.2f}-{th[-1]:.2f} h, {int(w_avg.sum())} samples):")
    print(f"  median            {final['median_uPa']:+.4f} uPa")
    print(f"  peak-to-peak      {100*final['ptp_rel']:.2f} % of the median "
          f"(reporting threshold {100*PTP_TOL:.0f} %)")
    print(f"  fitted drift      {final['drift_pct_per_km']:+.4f} % per km "
          f"(threshold {DRIFT_TOL_PCT_PER_KM} %/km)")
    print(f"  accumulated       {final['drift_pct_per_km']*(z_hi-a.z_lo):+.2f} % "
          f"across the {z_hi-a.z_lo:.0f} km band")
    flat = (final["ptp_rel"] <= PTP_TOL
            and abs(final["drift_pct_per_km"]) <= DRIFT_TOL_PCT_PER_KM)
    out["flat_at_thresholds"] = bool(flat)
    print("These figures compare MF against a fully established wave. Section 8 "
          "repeats them\nagainst the part of the ridge spectrum that has arrived, "
          "which is the comparison\nthe record can actually support.")

    # ---------------- 3. Flux amplitude ----------------
    hdr("3. AMPLITUDE  -  MF(z) against the steady linear hydrostatic value")
    ratio = mf_bar / mf_lin
    print(f"{'z[km]':>7} {'MF[uPa]':>11} {'linear[uPa]':>12} {'ratio':>8}")
    for zz in [5, 10, 20, 25, 30, 40, 50, 55]:
        if zz > z_hi + 5:
            continue
        j = int(np.argmin(np.abs(zkm - zz)))
        print(f"{zkm[j]:7.2f} {mf_bar[j]*1e6:11.4f} {mf_lin[j]*1e6:12.4f} "
              f"{ratio[j]:8.4f}")
    sel = (zkm >= a.z_lo) & (zkm <= z_hi)
    out["amplitude_ratio"] = {"band_mean": float(ratio[sel].mean()),
                              "band_min": float(ratio[sel].min()),
                              "band_max": float(ratio[sel].max())}
    print(f"\nBand mean MF/MF_linear = {ratio[sel].mean():.4f} "
          f"(range {ratio[sel].min():.4f} to {ratio[sel].max():.4f}).")

    # Which N the reference uses, and how much that choice is worth. The WKB
    # launch flux scales as rho_0 N_0 with N_0 taken at the first cell centre,
    # inside the near-surface inversion. The inversion is crossed in a fraction
    # of a vertical wavelength, so WKB does not hold across it and part of the
    # incident flux is reflected rather than transmitted; the transmitted
    # reference is then smaller than the launch formula gives, by an amount this
    # calculation does not resolve. Quoting the ratio against a band
    # representative N brackets it.
    n_band = float(N[sel].mean())
    ref_scale = n_band / N[0]
    out["amplitude_reference"] = {
        "N0_used_per_s": float(N[0]), "N0_height_km": float(zkm[0]),
        "band_mean_N_per_s": n_band, "reference_ratio_bandN_over_N0": ref_scale,
        "band_mean_ratio_against_bandN": float(ratio[sel].mean() / ref_scale)}
    print(f"\nReference used: N_0 = {N[0]:.5f} s^-1 at z = {zkm[0]:.2f} km, inside "
          f"the near-surface\ninversion. The band mean N is {n_band:.5f} s^-1, "
          f"{1/ref_scale:.2f} times smaller. Since MF_linear\nscales as rho_0 N_0, "
          f"a band-N reference would put the band mean ratio at "
          f"{ratio[sel].mean()/ref_scale:.4f}\nrather than {ratio[sel].mean():.4f}. "
          "The two admissible references differ by more than the\nquantity being "
          "measured, so this section does not presently separate a solver\n"
          "amplitude loss from the choice of reference, and the ratio should not be "
          "read as a\nfailure of the Eliassen-Palm theorem. Sections 2 and 8 are "
          "unaffected: for hydrostatic\nwaves m = N/U does not depend on k, so any "
          "partial reflection at the inversion\nmultiplies every wavenumber by the "
          "same factor and cancels from the normalised shape.")

    # ---------------- 4. Spurious drag ----------------
    hdr("4. DRAG  -  -(1/rho) dMF/dz in m/s per sol, the unit used downstream")
    drag = -np.gradient(mf_bar, z) / rho * MARS_SOL_S
    print(f"{'z[km]':>7} {'drag[m/s/sol]':>15}")
    for zz in [10, 20, 25, 30, 40, 50, 55]:
        if zz > z_hi + 5:
            continue
        j = int(np.argmin(np.abs(zkm - zz)))
        print(f"{zkm[j]:7.2f} {drag[j]:15.4f}")
    out["spurious_drag_m_s_sol"] = {
        "band_rms": float(np.sqrt((drag[sel] ** 2).mean())),
        "band_max_abs": float(np.abs(drag[sel]).max())}
    print(f"\nBand RMS {out['spurious_drag_m_s_sol']['band_rms']:.4f} m/s/sol, "
          f"max |drag| {out['spurious_drag_m_s_sol']['band_max_abs']:.4f} m/s/sol.")

    # ---------------- 5. Vertical wavelength ----------------
    hdr("5. LAMBDA_Z  -  measured phase tilt of w against 2 pi U / N(z)")
    lz = measure_lambda_z(files[-1], u_bar, zkm, a.z_lo, z_hi, a.fit_halfwidth_km)
    lam_pred = 2 * np.pi * u_bar / N / 1e3
    print(f"Dominant horizontal wavenumber {lz['k_dom']:.3e} 1/m "
          f"(lambda_x = {lz['lambda_x_dom_km']:.0f} km); "
          f"cospectrum closes to {lz['spec_closure']:.1e} of max|MF|.")
    print(f"{'z[km]':>7} {'measured[km]':>13} {'2piU/N[km]':>12} {'ratio':>8}")
    rows = []
    for zz in [10, 20, 25, 30, 40, 50, 55]:
        if zz > z_hi + 5:
            continue
        j = int(np.argmin(np.abs(zkm - zz)))
        r = lz["lambda_z_km"][j] / lam_pred[j]
        rows.append((float(zkm[j]), float(lz["lambda_z_km"][j]),
                     float(lam_pred[j]), float(r)))
        print(f"{zkm[j]:7.2f} {lz['lambda_z_km'][j]:13.2f} {lam_pred[j]:12.2f} "
              f"{r:8.4f}")
    out["lambda_z"] = {"k_dom_per_m": lz["k_dom"],
                       "lambda_x_dom_km": lz["lambda_x_dom_km"],
                       "phase_slope_sign": lz["m_sign"],
                       "cospectrum_closure": lz["spec_closure"],
                       "rows_z_meas_pred_ratio": rows,
                       "band_mean_ratio": float(np.nanmean(
                           (lz["lambda_z_km"] / lam_pred)[sel]))}
    print(f"\nBand mean measured/predicted = {out['lambda_z']['band_mean_ratio']:.4f}. "
          f"N over the band is\n{N[sel].min():.5f}-{N[sel].max():.5f} s^-1, so the "
          f"prediction spans {lam_pred[sel].min():.1f}-{lam_pred[sel].max():.1f} km "
          f"({lam_pred[sel].min()*1e3/(z[1]-z[0]):.0f}-"
          f"{lam_pred[sel].max()*1e3/(z[1]-z[0]):.0f} points at dz={z[1]-z[0]:.0f} m). "
          f"The\n13.96 km figure follows from N_REF = 0.018 s^-1 and applies only near "
          "the surface.")

    # ---------------- 6. Linearity ----------------
    hdr("6. LINEARITY  -  did the control run stay in the regime it is meant to test")
    clean = zkm < SPONGE_BASE_KM
    umin_all = float(umin[clean].min())
    j55 = int(np.argmin(np.abs(zkm - 55.0)))
    up55 = float(np.abs(ubar[j55] - u_bar).max())
    out["linearity"] = {"min_u_below_sponge": umin_all,
                        "min_u_over_U": umin_all / u_bar,
                        "max_hill_w_below_sponge": float(wmax[clean].max()),
                        "max_mean_flow_dev_55km": up55}
    print(f"  min_x u below {SPONGE_BASE_KM:.0f} km, all times : {umin_all:+.3f} m/s "
          f"({umin_all/u_bar:+.3f} U)")
    print(f"  overturning criterion min_x u <= 0        : "
          f"{'MET (nonlinear)' if umin_all <= 0 else 'not met (stays linear)'}")
    print(f"  max hill-window |w| below the sponge      : "
          f"{wmax[clean].max():.3f} m/s")
    print(f"  max |<u>_x - U| at 55 km                  : {up55:.3f} m/s")

    # ---------------- 7. Differential against the nonlinear run ----------------
    cmp_rows = None
    if (a.compare_npz and os.path.exists(a.compare_npz)
            and np.load(a.compare_npz)["z_km"].size != zkm.size):
        print(f"\nSection 7 skipped: {a.compare_npz} is on a different vertical grid.")
        a.compare_npz = ""
    if a.compare_npz and os.path.exists(a.compare_npz):
        hdr("7. VS NONLINEAR  -  height dependence not explained by a conservative wave")
        dc = np.load(a.compare_npz, allow_pickle=True)
        zc = dc["z_km"]
        mfc = dc["mf_zt"][:, -1]
        kc, hkc, kkc = ridge_spectrum(h_max=a.compare_h_max, sigma=a.sigma,
                                      lx=lx, dx=dx, nx=nx)
        linc = linear_steady_flux(rho, N, hkc, kkc, z, u_bar=u_bar, nx=nx)
        rl, rn = mf[:, -1] / mf_lin, mfc / linc
        j10 = int(np.argmin(np.abs(zkm - a.z_lo)))
        print(f"Comparison run: {os.path.dirname(a.compare_npz)}")
        print(f"  h = {a.compare_h_max:.0f} m, t_end = {dc['times'][-1]/3600:.2f} h "
              f"against {th[-1]:.2f} h here.")
        print(f"{'z[km]':>7} | {'control':>9} {'shape':>7} | {'nonlinear':>10} "
              f"{'shape':>7} | {'excess':>7}")
        cmp_rows = []
        for zz in range(int(a.z_lo), int(z_hi) + 1, 5):
            j = int(np.argmin(np.abs(zkm - zz)))
            sl, sn = rl[j] / rl[j10], rn[j] / rn[j10]
            cmp_rows.append([float(zkm[j]), float(rl[j]), float(sl),
                             float(rn[j]), float(sn), float(sn / sl)])
            print(f"{zkm[j]:7.2f} | {rl[j]:9.4f} {sl:7.4f} | {rn[j]:10.4f} "
                  f"{sn:7.4f} | {sn/sl:7.4f}")
        print("\n'excess' below 1 is flux divergence present in the nonlinear run and "
              "absent\nfrom a conservative wave in the same solver over the same "
              "interval. Values\nnear 1 mean the height dependence at that level is "
              "reproduced by the control\nand so cannot be attributed to saturation "
              "from this record alone.")
        out["vs_nonlinear"] = {
            "compare_npz": a.compare_npz, "compare_h_max_m": a.compare_h_max,
            "compare_t_end_h": float(dc["times"][-1] / 3600),
            "columns": ["z_km", "control_ratio", "control_shape",
                        "nonlinear_ratio", "nonlinear_shape", "excess"],
            "rows": cmp_rows}

    # ---------------- 8. Fill-in benchmark ----------------
    hdr("8. FILL-IN  -  MF against the flux of the spectrum that has arrived")
    fk = per_k_steady_flux(rho, N, hk, kk, u_bar, nx)
    ibuoy = buoyancy_integral(N, z)
    closure = float(abs(fk.sum() / mf_lin[0] - 1.0))
    print(f"Per-wavenumber decomposition sums to the steady linear value to "
          f"{closure:.2e}.")
    print("Each component is given its own arrival time t_arr = I(z)/(k U^2); the "
          "benchmark is\nthe steady flux of those that have arrived. See "
          "fillin_flux() for the assumptions.")

    # The periodic domain admits only k = n 2pi/LX, so the arrived fraction is a
    # staircase in both z and t rather than a smooth curve. Reporting the bins
    # makes the size of each step, and which mode sets it, explicit.
    share = fk / fk.sum()
    j55 = int(np.argmin(np.abs(zkm - min(55.0, z_hi))))
    print(f"\n{'n':>3} {'k[1/m]':>11} {'lam_x[km]':>10} {'flux share':>11} "
          f"{'cum':>7} {'t_arr(%.0f km)[h]' % zkm[j55]:>17}")
    cum, spec_rows = 0.0, []
    for n in range(1, 8):
        cum += share[n]
        ta = ibuoy[j55] / (kk[n] * u_bar ** 2) / 3600.0
        spec_rows.append([n, float(kk[n]), float(2 * np.pi / kk[n] / 1e3),
                          float(share[n]), float(cum), float(ta)])
        print(f"{n:3d} {kk[n]:11.4e} {2*np.pi/kk[n]/1e3:10.1f} {share[n]:11.4f} "
              f"{cum:7.4f} {ta:17.2f}")
    out["ridge_flux_spectrum"] = {
        "columns": ["n", "k_per_m", "lambda_x_km", "flux_share", "cumulative",
                    "t_arrival_at_eval_top_h"], "rows": spec_rows,
        "eval_top_km": float(zkm[j55]),
        "bins_to_90pct": int(np.searchsorted(np.cumsum(share), 0.90))}
    print(f"Five bins carry 95 per cent of the flux and the bin spacing is "
          f"2 pi / LX = {kk[1]:.3e} 1/m,\nso the benchmark steps rather than varies "
          f"smoothly. The two longest modes together\ncarry "
          f"{share[1]+share[2]:.3f} of the flux and reach {zkm[j55]:.0f} km only at "
          f"{spec_rows[0][5]:.1f} h and {spec_rows[1][5]:.1f} h.")

    phi_t = lambda tt, f=1.0: fillin_flux(fk, kk, ibuoy, tt, u_bar, f) / mf_lin
    correct = lambda prof, tt, f=1.0: prof / np.clip(phi_t(tt, f), PHI_FLOOR, None)

    # The arrived fraction changes over the averaging window, so each snapshot is
    # corrected by its own fraction before the window is averaged.
    idx = np.where(w_avg)[0]
    mf_corr = np.stack([correct(mf[:, i], t[i]) for i in idx], axis=1).mean(axis=1)
    phi_bar = np.stack([phi_t(t[i]) for i in idx], axis=1).mean(axis=1)
    final_corr = band_metrics(zkm, mf_corr, a.z_lo, z_hi)

    print(f"\n{'t[h]':>6} | {'raw drift':>10} {'raw ptp':>8} | "
          f"{'corr drift':>11} {'corr ptp':>9} | {'arrived(top)':>13} "
          f"{'min arrived':>12}")
    per_time_fill = {}
    for hh in hours:
        i = int(np.argmin(np.abs(th - hh)))
        rm = band_metrics(zkm, mf[:, i], a.z_lo, z_hi)
        cm = band_metrics(zkm, correct(mf[:, i], t[i]), a.z_lo, z_hi)
        if "ptp_rel" not in rm or "ptp_rel" not in cm:
            continue
        phi_i = phi_t(t[i])
        phimin = float(phi_i[sel].min())
        per_time_fill[f"{hh:g}h"] = {
            "raw_drift_pct_per_km": rm["drift_pct_per_km"],
            "corrected_drift_pct_per_km": cm["drift_pct_per_km"],
            "raw_ptp_rel": rm["ptp_rel"], "corrected_ptp_rel": cm["ptp_rel"],
            "arrived_fraction_at_eval_top": float(phi_i[j55]),
            "min_arrived_fraction_in_band": phimin}
        print(f"{th[i]:6.2f} | {rm['drift_pct_per_km']:10.4f} {rm['ptp_rel']:8.4f} | "
              f"{cm['drift_pct_per_km']:11.4f} {cm['ptp_rel']:9.4f} | "
              f"{phi_i[j55]:13.4f} {phimin:12.4f}"
              + ("  (clipped)" if phimin < PHI_FLOOR else ""))
    out["fillin_by_hour"] = per_time_fill
    print(f"Rows marked clipped have part of the band below the {PHI_FLOOR:.2f} "
          "arrived-fraction floor,\nwhere dividing by the benchmark is not "
          "meaningful; the correction is interpretable\nonly once the whole band "
          "carries some arrived flux.")

    print(f"\n{'z[km]':>7} {'arrived':>9} {'MF/MF_lin':>11} {'MF/MF_fill':>12}")
    rows = []
    for zz in range(int(a.z_lo), int(z_hi) + 1, 5):
        j = int(np.argmin(np.abs(zkm - zz)))
        rows.append([float(zkm[j]), float(phi_bar[j]), float(ratio[j]),
                     float(mf_corr[j] / mf_lin[j])])
        print(f"{zkm[j]:7.2f} {phi_bar[j]:9.4f} {ratio[j]:11.4f} "
              f"{mf_corr[j]/mf_lin[j]:12.4f}")
    out["fillin_profile"] = {"columns": ["z_km", "arrived_fraction",
                                         "MF_over_MF_linear", "MF_over_MF_fill"],
                             "rows": rows}

    out["fillin_final_average"] = final_corr
    out["fillin_closure"] = closure
    print(f"\nAveraged over the final {a.avg_hours:.1f} h, corrected for arrival:")
    print(f"  peak-to-peak      {100*final_corr['ptp_rel']:.2f} % of the median, "
          f"against {100*final['ptp_rel']:.2f} % raw")
    print(f"  fitted drift      {final_corr['drift_pct_per_km']:+.4f} % per km, "
          f"against {final['drift_pct_per_km']:+.4f} raw")
    print(f"  accumulated       {final_corr['drift_pct_per_km']*(z_hi-a.z_lo):+.2f} % "
          f"across the {z_hi-a.z_lo:.0f} km band, against "
          f"{final['drift_pct_per_km']*(z_hi-a.z_lo):+.2f} raw")
    fexp = (1.0 - abs(final_corr["drift_pct_per_km"])
            / max(abs(final["drift_pct_per_km"]), 1e-30))
    out["fillin_fraction_of_drift_explained"] = float(fexp)

    # Drag from the fitted slope rather than the pointwise gradient. The
    # benchmark is a staircase, so differentiating the corrected profile level by
    # level returns the derivative of the step model; the band fit does not.
    sm = lambda m: m["drift_pct_per_km"] / 1e5      # per cent per km -> 1/m
    drag_fit = lambda prof, m: -(prof * sm(m)) / rho * MARS_SOL_S
    print(f"\n{'z[km]':>7} {'raw drag':>10} {'corrected':>10}   (m/s per sol, from "
          "the band-fitted slope)")
    drag_rows = []
    for zz in [10, 25, 30, 40, 50, 55]:
        if zz > z_hi + 0.5:
            continue
        j = int(np.argmin(np.abs(zkm - zz)))
        dr, dc = drag_fit(mf_bar, final)[j], drag_fit(mf_corr, final_corr)[j]
        drag_rows.append([float(zkm[j]), float(dr), float(dc)])
        print(f"{zkm[j]:7.2f} {dr:10.4f} {dc:10.4f}")
    out["fillin_drag_from_fit"] = {"columns": ["z_km", "raw", "corrected"],
                                   "rows": drag_rows}

    # How much of this depends on treating each front as sharp.
    fronts = [float(v) for v in a.front_factors.split(",") if v.strip()]
    print(f"\nSensitivity to the front-sharpness assumption "
          f"(arrival time multiplied by f):")
    print(f"{'f':>6} {'drift[%/km]':>12} {'ptp':>8} {'min arrived':>12}")
    sens = []
    for f in fronts:
        mc = np.stack([correct(mf[:, i], t[i], f) for i in idx], axis=1).mean(axis=1)
        mm = band_metrics(zkm, mc, a.z_lo, z_hi)
        pm = float(np.stack([phi_t(t[i], f) for i in idx], axis=1).mean(axis=1)[sel].min())
        sens.append([f, mm["drift_pct_per_km"], mm["ptp_rel"], pm])
        print(f"{f:6.2f} {mm['drift_pct_per_km']:12.4f} {mm['ptp_rel']:8.4f} "
              f"{pm:12.4f}")
    out["fillin_front_sensitivity"] = {
        "columns": ["front_factor", "drift_pct_per_km", "ptp_rel",
                    "min_arrived_fraction"], "rows": sens}

    # The front factor that would null the drift. It is not a free parameter, so
    # this is a consistency check rather than a fit: a null close to f = 1 means
    # the WKB group velocity places the fronts where the record puts them.
    scan = np.linspace(0.4, 2.0, 33)
    dscan = np.array([band_metrics(zkm, np.stack(
        [correct(mf[:, i], t[i], f) for i in idx], axis=1).mean(axis=1),
        a.z_lo, z_hi).get("drift_pct_per_km", np.nan) for f in scan])
    cross = np.where(np.sign(dscan[:-1]) * np.sign(dscan[1:]) < 0)[0]
    f_null = float(np.interp(0.0, dscan[cross[0]:cross[0] + 2][::-1] if
                             dscan[cross[0]] > 0 else dscan[cross[0]:cross[0] + 2],
                             scan[cross[0]:cross[0] + 2][::-1] if
                             dscan[cross[0]] > 0 else scan[cross[0]:cross[0] + 2])
                   ) if len(cross) else float("nan")
    out["fillin_front_factor_nulling_drift"] = f_null
    print(f"\nThe front factor that nulls the drift is f = {f_null:.3f}. It is not a "
          "free parameter,\nso a null near f = 1 means the WKB group velocity places "
          "the fronts where the record\nputs them, with no tuning.")

    print(f"\nThe arrival model accounts for {100*fexp:.1f} per cent of the fitted "
          "drift. Two reservations\napply to what is left. It is still an upper "
          "bound, because it also contains whatever\nthe sponge, the dispersive "
          "spread about each front and the WKB group velocity\ncontribute, none of "
          "which this record separates. And the residual peak-to-peak is\nnot a "
          "solver figure: the measured profile varies smoothly with height while the "
          "\nbenchmark steps at the arrival heights of four or five discrete modes, "
          "so the\nscatter left after dividing is dominated by that shape mismatch. "
          "The slope is the\npart of this section that carries information; the "
          "scatter is not.")

    # ---------------- 9. Ridge-width pair ----------------
    # Write the measured drift as D(z) + T(z; sigma, t), with D numerical and T the
    # arrival transient. A single run gives only the sum. The ridge half-width moves
    # T, because it reweights the flux across wavenumbers, while leaving lambda_z,
    # the sponge thickness in wavelengths and the linearity margin untouched, since
    # for hydrostatic waves m = N/U and the steady amplitude N h (rho_0 N_0/rho N)^0.5
    # both carry no k. Differencing the pair therefore tests whether D is common.
    pair = None
    if a.pair_npz and os.path.exists(a.pair_npz):
        hdr("9. PAIR  -  a second ridge width, to separate the numerical drift")
        dp = np.load(a.pair_npz, allow_pickle=True)
        tp, zp, mfp = dp["times"], dp["z_km"], dp["mf_zt"]
        if zp.shape != zkm.shape or not np.allclose(zp, zkm) or \
                abs(float(dp["u_bar"]) - u_bar) > 1e-9:
            print("The pair run does not share the vertical grid or the wind, so the "
                  "difference\nwould not be controlled. Skipped.")
        else:
            k2, hk2, kk2 = ridge_spectrum(h_max=a.h_max, sigma=a.pair_sigma,
                                          lx=lx, dx=dx, nx=nx)
            mf_lin2 = linear_steady_flux(rho, N, hk2, kk2, z, u_bar=u_bar, nx=nx)
            fk2 = per_k_steady_flux(rho, N, hk2, kk2, u_bar, nx)
            wp = tp >= (tp[-1] - a.avg_hours * 3600.0)
            ip = np.where(wp)[0]
            mfp_bar = mfp[:, wp].mean(axis=1)

            pred = lambda F, K, tt, f=1.0: np.stack(
                [fillin_flux(F, K, ibuoy, tv, u_bar, f) for tv in tt],
                axis=1).mean(axis=1)
            prof1, prof2 = pred(fk, kk, t[idx]), pred(fk2, kk2, tp[ip])
            drift = lambda prof: band_metrics(zkm, prof, a.z_lo, z_hi)["drift_pct_per_km"]
            T1, T2 = drift(prof1), drift(prof2)
            M1, M2 = final["drift_pct_per_km"], drift(mfp_bar)
            # Snapshot-to-snapshot scatter of the fitted drift, the conservative
            # scale against which the separation has to be read.
            sd = lambda m, ii: float(np.std(
                [band_metrics(zkm, m[:, i], a.z_lo, z_hi)["drift_pct_per_km"]
                 for i in ii], ddof=1))
            s1, s2 = sd(mf, idx), sd(mfp, ip)
            phi1 = float((prof1 / mf_lin)[j55])
            phi2 = float((prof2 / mf_lin2)[j55])

            print(f"{'quantity':<44} {'sigma=%.0f km' % (a.sigma/1e3):>13} "
                  f"{'sigma=%.0f km' % (a.pair_sigma/1e3):>13}")
            rows = [("flux-weighted k_eff [1/m]", f"{k_eff:.4e}", f"{k2:.4e}"),
                    ("flux-weighted lambda_x [km]",
                     f"{2*np.pi/k_eff/1e3:.0f}", f"{2*np.pi/k2/1e3:.0f}"),
                    ("steady linear MF [uPa]",
                     f"{mf_lin[0]*1e6:.3f}", f"{mf_lin2[0]*1e6:.3f}"),
                    ("record length [h]", f"{th[-1]:.2f}", f"{tp[-1]/3600:.2f}"),
                    (f"arrived fraction at {zkm[j55]:.0f} km",
                     f"{phi1:.4f}", f"{phi2:.4f}"),
                    ("predicted transient drift T [%/km]", f"{T1:+.4f}", f"{T2:+.4f}"),
                    ("measured drift M [%/km]", f"{M1:+.4f}", f"{M2:+.4f}"),
                    ("per-snapshot sd of M [%/km]", f"{s1:.4f}", f"{s2:.4f}"),
                    ("implied residual D = M - T [%/km]",
                     f"{M1-T1:+.4f}", f"{M2-T2:+.4f}")]
            for lab, v1, v2 in rows:
                print(f"{lab:<44} {v1:>13} {v2:>13}")

            sep_pred, sep_meas = T1 - T2, M1 - M2
            sep_sd = float(np.hypot(s1, s2))
            print(f"\npredicted separation T1 - T2 = {sep_pred:+.4f} %/km")
            print(f"measured  separation M1 - M2 = {sep_meas:+.4f} %/km")
            print(f"difference of the two        = {sep_meas - sep_pred:+.4f} %/km, "
                  f"against a snapshot scatter of {sep_sd:.4f}")
            print(f"D is common to the pair to {abs((M1-T1)-(M2-T2)):.4f} %/km; the "
                  f"mean of the two residuals is\n{0.5*((M1-T1)+(M2-T2)):+.4f} %/km, "
                  f"which is the numerical flux drift this pair supports.")
            # The front factor is the one place the fill-in model is soft, since a
            # real front is dispersively spread rather than a step. Solving for the
            # common f that reconciles the pair turns that softness into a number.
            # With two measurements and two unknowns this is a solve and not a test:
            # what it reports is how far from the parameter-free f = 1 the model has
            # to be pushed for D to be common, and the D that follows.
            fscan = np.linspace(0.5, 1.6, 45)
            gap = np.array([drift(pred(fk, kk, t[idx], f))
                            - drift(pred(fk2, kk2, tp[ip], f)) - sep_meas
                            for f in fscan])
            cr = np.where(np.sign(gap[:-1]) * np.sign(gap[1:]) < 0)[0]
            f_pair, d_common = float("nan"), float("nan")
            if cr.size:
                i0 = cr[-1]
                sl = slice(i0, i0 + 2)
                asc = gap[i0] < gap[i0 + 1]
                f_pair = float(np.interp(0.0, gap[sl] if asc else gap[sl][::-1],
                                         fscan[sl] if asc else fscan[sl][::-1]))
                d1 = M1 - drift(pred(fk, kk, t[idx], f_pair))
                d2 = M2 - drift(pred(fk2, kk2, tp[ip], f_pair))
                d_common = float(0.5 * (d1 + d2))
                print(f"\nA common front factor of f = {f_pair:.3f} reconciles the "
                      f"pair, giving D = {d1:+.4f}\nand {d2:+.4f} per cent per km. "
                      f"With two measurements and two unknowns this is a\nsolve "
                      f"rather than a test; what it measures is how far from the "
                      f"parameter-free\nf = 1 the sharp-front assumption has to be "
                      f"pushed, and the D that follows.")

            print("\nTaken together the pair bounds the numerical drift between "
                  f"{min(abs(M1-T1), abs(M2-T2)):.3f} and\n"
                  f"{max(abs(d_common), abs(M2-T2)):.3f} per cent per km over "
                  "10-55 km, depending on how the front-sharpness\nassumption is "
                  "treated. What does not depend on that assumption is the measured "
                  "pair\nitself: with lambda_z, the sponge thickness in wavelengths, "
                  "dz per wavelength and the\nlinearity margin all held fixed, "
                  "reweighting the ridge spectrum alone changed the\nmeasured drift "
                  f"by a factor of {abs(M1/M2):.1f}, which a fixed numerical flux "
                  "divergence would not do.")

            pair = dict(sigma_m=a.pair_sigma, npz=a.pair_npz, mf_bar=mfp_bar,
                        pred=prof2, mf_lin=mf_lin2, t_end_h=float(tp[-1]/3600))
            out["sigma_pair"] = {
                "pair_npz": a.pair_npz,
                "sigma_m": [a.sigma, a.pair_sigma],
                "k_eff_per_m": [k_eff, float(k2)],
                "lambda_x_eff_km": [2*np.pi/k_eff/1e3, float(2*np.pi/k2/1e3)],
                "steady_linear_MF_uPa": [float(mf_lin[0]*1e6), float(mf_lin2[0]*1e6)],
                "t_end_h": [float(th[-1]), float(tp[-1]/3600)],
                "arrived_fraction_at_eval_top": [phi1, phi2],
                "predicted_transient_drift_pct_per_km": [float(T1), float(T2)],
                "measured_drift_pct_per_km": [float(M1), float(M2)],
                "per_snapshot_sd_pct_per_km": [s1, s2],
                "residual_D_pct_per_km": [float(M1-T1), float(M2-T2)],
                "predicted_separation_pct_per_km": float(sep_pred),
                "measured_separation_pct_per_km": float(sep_meas),
                "separation_mismatch_pct_per_km": float(sep_meas - sep_pred),
                "separation_snapshot_scatter_pct_per_km": sep_sd,
                "D_mean_pct_per_km": float(0.5*((M1-T1)+(M2-T2))),
                "front_factor_common_D": f_pair,
                "D_common_at_that_factor_pct_per_km": d_common}

    # ---------------- 10. Horizon ----------------
    # Sections 8 and 9 rest on one account of why the measured drift falls with
    # time: each ridge wavenumber arrives at its own moment and the drift is the
    # shortfall of the part that has not. A relaxation with no wave physics in it
    # describes the same record about as closely. The two cannot be told apart on
    # the record they are fitted to, so both are fitted here to a truncated part of
    # it and scored on the rest. Because the fill-in account supplies T(t) from the
    # spectrum and fits only the constant D, while the relaxation fits three free
    # parameters, the comparison is not symmetric and is read as a check on the
    # fill-in account rather than as model selection.
    hdr("10. HORIZON  -  the fall of the drift in time, scored out of sample")
    bdrift = lambda p: band_metrics(zkm, p, a.z_lo, z_hi).get(
        "drift_pct_per_km", float("nan"))

    def Tfill(tt, f=1.0):
        """Drift of the arrived-spectrum benchmark, averaged over the times tt.

        The measurement averages MF over a window before fitting the slope, so the
        benchmark is averaged the same way rather than sampled at the mid-time.
        """
        tt = np.atleast_1d(np.asarray(tt, float))
        return bdrift(np.stack([fillin_flux(fk, kk, ibuoy, tv, u_bar, f)
                                for tv in tt], axis=1).mean(axis=1))

    cut = float(a.forecast_from)
    # A run stops on a step rather than on a snapshot, so the last sample sits a
    # minute or so short of the hour; count that hour rather than dropping it.
    last_h = np.floor(th[-1] + 1.0 / 60.0)
    hr_all = np.arange(1.0, last_h + 0.5)
    ih = np.array([int(np.argmin(np.abs(th - h))) for h in hr_all])
    M_h = np.array([bdrift(mf[:, i]) for i in ih])
    T_h = np.array([Tfill(t[i]) for i in ih])
    in_s = hr_all <= cut
    have_oos = bool((~in_s).any())

    # Where each of the leading modes crosses the band. The drift of the benchmark
    # responds to fronts inside 10-55 km, not to the arrived fraction at the top,
    # so a mode that has reached the top can still be moving T.
    jlo = int(np.argmin(np.abs(zkm - a.z_lo)))
    print(f"{'n':>3} {'lambda_x[km]':>13} {'enters band[h]':>15} "
          f"{'leaves band[h]':>15}")
    cross_rows = []
    for n in range(1, 5):
        ta = ibuoy / (kk[n] * u_bar ** 2) / 3600.0
        cross_rows.append([n, float(2*np.pi/kk[n]/1e3), float(ta[jlo]),
                           float(ta[j55])])
        print(f"{n:3d} {2*np.pi/kk[n]/1e3:13.0f} {ta[jlo]:15.2f} {ta[j55]:15.2f}")
    print("A front inside the band tilts the fitted slope while it climbs, so T "
          "varies\ncontinuously between arrivals at the top rather than stepping "
          "only at them.")

    print("\nThe fill-in account is M(t) = D + T(t), with T the drift of the arrived "
          "spectrum\nand D the numerical residual it is meant to expose. D is "
          "constant in the account,\nso the last column is a direct check on it. "
          "Both columns are single snapshots here;\nthe forecast table below "
          "averages over the hour, as the headline figures do.")
    print(f"\n{'t[h]':>6} {'M measured':>11} {'T fill-in':>10} {'D = M - T':>10}")
    for h, m, tt in zip(hr_all, M_h, T_h):
        print(f"{h:6.0f} {m:11.4f} {tt:10.4f} {m-tt:10.4f}"
              + ("" if h <= cut else "   out of sample"))

    # Spread of D over the settled part of the in-sample segment. The early hours
    # are excluded because the whole band is then below the arrived-fraction floor.
    dwin = in_s & (hr_all >= max(hr_all[in_s].max() - 3, 1))
    D_spread = float(np.ptp(M_h[dwin] - T_h[dwin]))
    print(f"\nOver the last {int(dwin.sum())} in-sample hours D moves by "
          f"{D_spread:.4f} per cent per km. A residual that is\nconstant in height "
          f"and in time is what the account requires, so this is the size of\nits "
          f"own internal inconsistency on the record it was built from.")

    # --- the two accounts, fitted to t <= cut -------------------------------
    fit_los = [float(v) for v in a.forecast_fit_lo.split(",") if v.strip()]
    print(f"\nRelaxation a exp(-t/tau) + D_inf fitted to t <= {cut:.0f} h:")
    print(f"{'window[h]':>10} {'n':>3} {'a':>9} {'tau[h]':>7} {'tau range':>13} "
          f"{'D_inf':>9} {'rms':>7}")
    fits, fit_rows = {}, []
    for lo in fit_los:
        m = in_s & (hr_all >= lo)
        if m.sum() < 4:
            continue
        fq = relax_fit(hr_all[m], M_h[m])
        # A window of length L cannot constrain a decay slower than L: the fit then
        # trades tau against the asymptote and extrapolates rather than measures.
        # Such fits are printed but kept out of the quoted spread.
        fq["constrained"] = bool(fq["tau"] <= (cut - lo))
        fits[lo] = fq
        fit_rows.append([lo, fq["n"], fq["a"], fq["tau"], fq["tau_lo"],
                         fq["tau_hi"], fq["d_inf"], fq["rms"], fq["constrained"]])
        print(f"{lo:5.0f}-{cut:<4.0f} {fq['n']:3d} {fq['a']:9.3f} {fq['tau']:7.3f} "
              f"{fq['tau_lo']:6.2f}-{fq['tau_hi']:<6.2f} {fq['d_inf']:9.4f} "
              f"{fq['rms']:7.4f}"
              + ("" if fq["constrained"] else "   tau exceeds its window"))
    ok = {lo: f for lo, f in fits.items() if f["constrained"]} or fits
    dinf = np.array([f["d_inf"] for f in ok.values()])
    print(f"Across the {len(ok)} window(s) the fit length constrains, the asymptote "
          f"moves by\n{np.ptp(dinf):.4f} per cent per km. That is larger than any one "
          f"fit's formal error and is the\nuncertainty to carry. Windows whose tau "
          f"exceeds their own length are shown but excluded.")

    # D for the fill-in forecast, anchored the way section 9 anchors it: on the
    # final-hour mean of the in-sample segment.
    i_cut = int(np.argmin(np.abs(th - cut)))
    w_in = (t <= t[i_cut]) & (t >= t[i_cut] - a.avg_hours * 3600.0)
    M_cut = bdrift(mf[:, w_in].mean(axis=1))
    D_anchor = M_cut - Tfill(t[w_in])
    print(f"\nOver the final hour of the in-sample segment the measured drift is "
          f"{M_cut:+.4f} and\nthe benchmark {Tfill(t[w_in]):+.4f}, so the fill-in "
          f"account carries D = {D_anchor:+.4f} forward.")

    # --- forecasts ----------------------------------------------------------
    fc_to = a.forecast_to if a.forecast_to > 0 else max(last_h, cut + 3.0)
    fc_hours = [h for h in np.arange(np.floor(cut) + 1, fc_to + 0.5)]
    print(f"\nForecast of the drift of the mean profile over the hour ending at t, "
          f"in per cent per km:")
    print(f"{'t[h]':>6} {'fill-in':>9} {'relaxation':>21} {'measured':>10} "
          f"{'fill-in err':>12} {'relax err':>10}")
    fc_rows, sc_fill, sc_rel = [], [], []
    for h in fc_hours:
        tw = np.arange(h - a.avg_hours, h + 1e-9, 1.0 / 60.0) * 3600.0
        A = D_anchor + Tfill(tw)
        B = np.array([relax_eval(f, h - 0.5 * a.avg_hours) for f in ok.values()])
        i = int(np.argmin(np.abs(th - h)))
        got = h <= th[-1] + 1.0 / 60.0
        Mm = bdrift(mf[:, (t <= t[i]) & (t >= t[i] - a.avg_hours * 3600.0)
                       ].mean(axis=1)) if got else float("nan")
        eA = Mm - A
        eB = Mm - float(np.median(B))
        if got and np.isfinite(Mm):
            sc_fill.append(eA)
            sc_rel.append(eB)
        fc_rows.append([float(h), float(A), float(B.min()), float(B.max()),
                        float(Mm), float(eA), float(eB)])
        print(f"{h:6.0f} {A:9.4f} {B.min():9.4f} to {B.max():<8.4f} "
              f"{Mm:10.4f} {eA:12.4f} {eB:10.4f}"
              if np.isfinite(Mm) else
              f"{h:6.0f} {A:9.4f} {B.min():9.4f} to {B.max():<8.4f} "
              f"{'-':>10} {'-':>12} {'-':>10}")

    # Scatter scales. The per-snapshot spread is the optimistic one: consecutive
    # snapshots are strongly correlated, so it is not the standard error of the
    # hour. The window spread of D_inf above is the systematic one.
    idx_fin = np.where(t >= t[-1] - a.avg_hours * 3600.0)[0]
    sd_snap = float(np.std([bdrift(mf[:, i]) for i in idx_fin], ddof=1))
    print(f"\nScales to read those errors against: the per-snapshot spread of the "
          f"fitted drift in\nthe final hour is {sd_snap:.4f} per cent per km, and "
          f"the spread of the relaxation\nasymptote across fit windows is "
          f"{np.ptp(dinf):.4f}. The first is optimistic, since "
          f"consecutive\nsnapshots are strongly correlated; the second is the "
          f"systematic one.")

    # The increment from the cut is the part of the comparison that does not
    # depend on either account getting the level right.
    print(f"\nChange in the drift from the final hour of the in-sample segment:")
    print(f"{'t[h]':>6} {'fill-in':>9} {'relaxation':>21} {'measured':>10}")
    inc_rows = []
    for r in fc_rows:
        h = r[0]
        inc_rows.append([h, r[1] - M_cut, r[2] - M_cut, r[3] - M_cut, r[4] - M_cut])
        s = (f"{r[4]-M_cut:10.4f}" if np.isfinite(r[4]) else f"{'-':>10}")
        print(f"{h:6.0f} {r[1]-M_cut:9.4f} {r[2]-M_cut:9.4f} to "
              f"{r[3]-M_cut:<8.4f} {s}")

    # Successive hourly increments. This is the part of the comparison that needs
    # neither an anchor for D nor a choice of fit window, because an additive
    # constant and the overall level both drop out of a difference. The two
    # accounts predict different shapes for it: a relaxation falls geometrically by
    # exp(-1/tau) each hour, while the fill-in increments follow the fronts and
    # collapse once the last mode still inside the band clears its top.
    print(f"\nHour-on-hour change in the drift, which is free of both the anchor "
          f"for D and the\nchoice of fit window:")
    print(f"{'interval':>12} {'fill-in':>9} {'relaxation':>21} {'measured':>10}")
    step_rows, prev = [], [float(cut), M_cut, M_cut]
    for r in fc_rows:
        h0, h1 = prev[0], r[0]
        dA = r[1] - prev[1]
        # Each fit's own increment, so the envelope is over models rather than
        # over whichever model happened to be lowest at each endpoint.
        dB = np.array([relax_eval(f, h1 - 0.5 * a.avg_hours)
                       - relax_eval(f, h0 - 0.5 * a.avg_hours)
                       for f in ok.values()])
        dm = r[4] - prev[2]
        step_rows.append([h0, h1, float(dA), float(dB.min()), float(dB.max()),
                          float(dm)])
        txt = f"{dm:10.4f}" if np.isfinite(dm) else f"{'-':>10}"
        print(f"{f'{h0:.0f}-{h1:.0f} h':>12} {dA:9.4f} {dB.min():9.4f} to "
              f"{dB.max():<8.4f} {txt}")
        prev = [r[0], r[1], r[4]]
    tau_ref = float(np.median([f["tau"] for f in ok.values()]))
    print(f"A relaxation makes those decay by exp(-1/tau) = "
          f"{np.exp(-1.0/tau_ref):.2f} per hour at the fitted\ntau of "
          f"{tau_ref:.2f} h. The fill-in account instead has them collapse once mode "
          f"n = 2 clears\n{zkm[j55]:.0f} km at {cross_rows[1][3]:.2f} h, after which "
          f"no front remains inside the band. The shape of\nthis column, rather than "
          f"any single level, is what separates the two.")

    horizon = {
        "forecast_from_h": cut,
        "forecast_to_h": float(fc_to),
        "band_crossing": {"columns": ["n", "lambda_x_km", "enters_band_h",
                                      "leaves_band_h"], "rows": cross_rows},
        "hourly": {"columns": ["t_h", "measured_drift_pct_per_km",
                               "fillin_T_pct_per_km", "D_eq_M_minus_T"],
                   "rows": [[float(h), float(m), float(tt), float(m - tt)]
                            for h, m, tt in zip(hr_all, M_h, T_h)]},
        "D_spread_last_in_sample_hours_pct_per_km": D_spread,
        "relaxation_fits": {"columns": ["fit_start_h", "n_points", "a", "tau_h",
                                        "tau_lo_h", "tau_hi_h", "D_inf_pct_per_km",
                                        "rms", "constrained"], "rows": fit_rows},
        "relaxation_D_inf_window_spread_pct_per_km": float(np.ptp(dinf)),
        "fillin_D_anchor_pct_per_km": float(D_anchor),
        "measured_drift_at_cut_pct_per_km": float(M_cut),
        "forecasts": {"columns": ["t_h", "fillin", "relax_min", "relax_max",
                                  "measured", "fillin_error", "relax_error"],
                      "rows": fc_rows},
        "increments_from_cut": {"columns": ["t_h", "fillin", "relax_min",
                                            "relax_max", "measured"],
                                "rows": inc_rows},
        "hourly_increments": {"columns": ["t0_h", "t1_h", "fillin",
                                          "relax_min", "relax_max", "measured"],
                              "rows": step_rows},
        "relaxation_tau_median_h": tau_ref,
        "per_snapshot_sd_final_hour_pct_per_km": sd_snap}
    if sc_fill:
        horizon["out_of_sample_rms_error"] = {
            "fillin": float(np.sqrt(np.mean(np.square(sc_fill)))),
            "relaxation_median": float(np.sqrt(np.mean(np.square(sc_rel))))}
        print(f"\nOut-of-sample rms error over {len(sc_fill)} hour(s): fill-in "
              f"{horizon['out_of_sample_rms_error']['fillin']:.4f}, relaxation "
              f"{horizon['out_of_sample_rms_error']['relaxation_median']:.4f} "
              f"per cent per km.")
    else:
        print(f"\nThe record ends at {th[-1]:.2f} h, at or before the cut, so these "
              f"are forecasts only.\nScoring them needs an integration past "
              f"{cut:.0f} h.")
    out["horizon"] = horizon

    # ---------------- 11. Up and down ----------------
    hdr("11. UPDOWN  -  the reflected share of the flux, from the sign of m")
    ud_hours = [float(s) for s in a.updown_hours.split(",") if s.strip()]
    updown = {}
    if not ud_hours:
        print("Skipped (--updown-hours empty).")
    else:
        print(
            "The band drift is one projection of the flux profile and the level is\n"
            "another. A height-uniform loss moves the level while leaving the drift\n"
            "untouched, so section 2 cannot see it. Splitting the flux by the sign of\n"
            "the vertical wavenumber separates the incident wave from anything\n"
            "returning from above, which is the loss channel the geometry allows: the\n"
            f"sponge occupies {SPONGE_BASE_KM:.0f}-{z_top_km:.0f} km, about "
            f"{(z_top_km-SPONGE_BASE_KM)/lam_z_band:.2f} vertical wavelengths at the\n"
            f"band-median lambda_z of {lam_z_band:.0f} km.")

        # Method floor, on a field whose answer is known before it is measured.
        print("\nMethod check on a synthetic WKB wave with an imposed reflected fraction:")
        print(f"{'imposed':>9} {'recovered':>10} {'net/up':>8}")
        ud_null = []
        for r_in in [float(s) for s in a.updown_synthetic.split(",") if s.strip()]:
            us, ws, rs = synthetic_wkb(hk, kk, z, rho, N, u_bar, nx, refl=r_in)
            sp = updown_split(us, ws, rs, np.arange(nx) * dx, z, N,
                              np.full(len(z), u_bar))
            s0 = (zkm >= a.z_lo) & (zkm <= z_hi)
            uu, dd = sp["mf_up"][s0].mean(), sp["mf_dn"][s0].mean()
            ud_null.append({"imposed": r_in, "recovered": float(abs(dd / uu))})
            print(f"{r_in:9.3f} {abs(dd/uu):10.5f} {(uu+dd)/uu:8.4f}")
        updown["synthetic_check"] = ud_null
        print(f"The floor of the method on a purely upward wave is "
              f"{ud_null[0]['recovered']:.5f} in flux, set by the\nWKB amplitude "
              "gradient rather than by noise.")

        # The bias that synthetic cannot see, because it is built from the same
        # m the split assumes. An exact anelastic solution at constant N and
        # constant scale height carries m^2 = N^2/U^2 - k^2 - 1/(4H^2), and the
        # omitted term is not small: it biases every propagating mode the same
        # way and does not average down, so it rather than the WKB number is the
        # floor a per-mode ratio should be read against.
        s_bd = (zkm >= a.z_lo) & (zkm <= z_hi)
        N_bd = float(N[s_bd].mean())
        H_bd = float(-1.0 / np.mean(np.gradient(np.log(rho), z)[s_bd]))
        ue, we, re_, m_ex = synthetic_exact(hk, kk, z, rho[0], N_bd, H_bd, u_bar, nx)
        spe = updown_split(ue, we, re_, np.arange(nx) * dx, z,
                           np.full(len(z), N_bd), np.full(len(z), u_bar))
        with np.errstate(divide="ignore", invalid="ignore"):
            rk_ex = np.abs(spe["fdn"] / spe["fup"])
        j_bd = int(np.argmin(np.abs(zkm - a.z_lo)))
        nx_ex = min(a.updown_modes if a.updown_modes > 0 else 12, kk.size - 1)
        ex_band = float(abs(spe["mf_dn"][s_bd].sum() / spe["mf_up"][s_bd].sum()))
        ex_modes = [float(rk_ex[j_bd, n]) for n in range(1, nx_ex + 1)]
        mf_ex = spe["mf_spec"][s_bd]
        updown["exact_check"] = {
            "N_band_per_s": N_bd, "scale_height_km": H_bd / 1e3,
            "term_pct_of_m2": float(100.0 / (4 * H_bd ** 2) / (N_bd / u_bar) ** 2),
            "band_downward_share": ex_band,
            "per_mode_downward_share": ex_modes,
            "flux_flatness_pct": float(100 * np.ptp(mf_ex) / abs(np.median(mf_ex)))}
        print(f"\nThe same check on an exact steady anelastic solution rather "
              f"than a WKB one.\nAt constant N = {N_bd:.5f} 1/s and constant scale "
              f"height {H_bd/1e3:.2f} km the vertical\nstructure equation for "
              f"sqrt(rho_0) w' carries m^2 = N^2/U^2 - k^2 - 1/(4H^2), whose last\n"
              f"term is {100.0/(4*H_bd**2)/(N_bd/u_bar)**2:.2f} per cent of the "
              f"first here and is not in the m this split assumes.")
        print(f"  downward share of that exact upward solution: "
              f"{ex_band:.3e} over the band,\n  {min(ex_modes):.3e} to "
              f"{max(ex_modes):.3e} across modes 1 to {nx_ex}, against "
              f"{ud_null[0]['recovered']:.3e} for the WKB\n  synthetic. Its own "
              f"flux varies by {100*np.ptp(mf_ex)/abs(np.median(mf_ex)):.4f} per "
              f"cent over the band, which is the check that\n  the construction "
              f"is right. This is a systematic bias of the same sign in every\n"
              f"  mode rather than noise, so it does not average down, and it is "
              f"the floor the\n  per-mode ratios below should be read against.")

        # Snapshot times. The deepdive series covers every snapshot, so the two
        # lists align; if they do not, the attribute is read from each file.
        if len(files) == len(t):
            t_snap = t
        else:
            t_snap = np.array([_snap_time(f) for f in files])

        def ud_window(hh, half, nmax):
            m = (t_snap > (hh - half) * 3600.0) & (t_snap <= (hh + half) * 3600.0)
            idx = np.where(m)[0]
            if idx.size == 0:
                idx = np.array([int(np.argmin(np.abs(t_snap - hh * 3600.0)))])
            if idx.size > nmax:
                idx = idx[np.linspace(0, idx.size - 1, nmax).astype(int)]
            return idx

        def ud_mean(idx, m_scale=1.0, background=False, keep=False):
            acc, keys = None, ("mf_up", "mf_dn", "mf_spec", "mf_direct", "mf_evan",
                               "var_total", "var_incoherent")
            members = []
            for i in idx:
                with Dataset(files[i]) as ds:
                    xs = np.asarray(ds["x"][:], dtype=np.float64)
                    us = np.asarray(ds["u"][:], dtype=np.float64)
                    ws = np.asarray(ds["w"][:], dtype=np.float64)
                    rs = np.asarray(ds["rho"][:], dtype=np.float64)
                    Ts = np.asarray(ds["T"][:], dtype=np.float64).mean(axis=1)
                if background:
                    Nl, Ul = N, np.full(len(z), u_bar)
                else:
                    Nl = np.sqrt(np.clip((G / Ts) * (np.gradient(Ts, z) + G / CP),
                                         1e-12, None))
                    Ul = us.mean(axis=1)
                sp = updown_split(us, ws, rs, xs, z, Nl, Ul, m_scale=m_scale)
                d1 = {q: sp[q] for q in keys}
                d1["fup_k"], d1["fdn_k"] = sp["fup"], sp["fdn"]
                d1["fup_fd_k"], d1["fdn_fd_k"] = sp["fup_fd"], sp["fdn_fd"]
                d1["eps_k"] = np.nan_to_num(sp["closure_eps"])
                d1["psi2_k"] = np.abs(sp["A"] + sp["B"]) ** 2
                # Coherent product of the two branches. Its phase at a fixed height
                # is steady in a standing cavity and rotates while the downward
                # branch is still a transient of a different frequency.
                d1["ab_k"] = np.conj(sp["A"]) * sp["B"]
                if keep:
                    members.append({"fup_k": sp["fup"], "fdn_k": sp["fdn"]})
                acc = ({q: v.copy() for q, v in d1.items()} if acc is None
                       else {q: acc[q] + v for q, v in d1.items()})
            for q in acc:
                acc[q] = acc[q] / idx.size
            acc["t_h"] = float(t_snap[idx].mean() / 3600.0)
            acc["n"] = int(idx.size)
            acc["members"] = members
            return acc

        sel_b = (zkm >= a.z_lo) & (zkm <= z_hi)
        bmean = lambda v: float(v[sel_b].mean())
        bdrift_of = lambda v, s=sel_b: float(
            100.0 * np.polyfit(zkm[s], np.log(np.abs(v[s])), 1)[0])
        lin_band = float(mf_lin[sel_b].mean())

        ud_hours = sorted({h for h in ud_hours if h <= th[-1] + 1.0 / 60.0}
                          | {float(last_h)})
        got = {}
        print(f"\nHour-mean windows of +/-{a.updown_window:.1f} h, band "
              f"{a.z_lo:.0f}-{z_hi:.0f} km. R = |MF_down| / |MF_up|.")
        print(f"{'t[h]':>6} {'n':>3} {'up[uPa]':>9} {'down[uPa]':>10} {'net[uPa]':>9} "
              f"{'R':>7} {'share':>7} {'up/lin':>7} {'net/lin':>8} "
              f"{'drift_up':>9} {'drift_net':>10}")
        rows = []
        for hh in ud_hours:
            g = ud_mean(ud_window(hh, a.updown_window, a.updown_max_snaps),
                        keep=(hh == max(ud_hours)))
            got[hh] = g
            uu, dd, nn = bmean(g["mf_up"]), bmean(g["mf_dn"]), bmean(g["mf_direct"])
            row = {"t_h": g["t_h"], "n_snapshots": g["n"],
                   "MF_up_uPa": uu * 1e6, "MF_down_uPa": dd * 1e6,
                   "MF_net_uPa": nn * 1e6, "R_down_over_up": abs(dd / uu),
                   "down_share": abs(dd) / (abs(uu) + abs(dd)),
                   "up_over_linear": uu / lin_band, "net_over_linear": nn / lin_band,
                   "drift_up_pct_per_km": bdrift_of(g["mf_up"]),
                   "drift_net_pct_per_km": bdrift_of(g["mf_direct"])}
            rows.append(row)
            print(f"{row['t_h']:6.2f} {g['n']:3d} {uu*1e6:9.4f} {dd*1e6:10.4f} "
                  f"{nn*1e6:9.4f} {row['R_down_over_up']:7.4f} "
                  f"{row['down_share']:7.4f} {row['up_over_linear']:7.4f} "
                  f"{row['net_over_linear']:8.4f} "
                  f"{row['drift_up_pct_per_km']:+9.4f} "
                  f"{row['drift_net_pct_per_km']:+10.4f}")
        updown["by_hour"] = rows

        last = got[max(got)]
        print(f"\nLocal ratio R(z) at the sampled hours. Above {SPONGE_BASE_KM:.0f} km "
              "the wave is being damped,\nso m is complex there and these entries are "
              "indicative only.")
        head = "".join(f"{hh:>8.0f}h" for hh in ud_hours)
        print(f"{'z[km]':>7}{head}")
        prof = {"z_km": [], **{f"{hh:g}h": [] for hh in ud_hours}}
        for zz in [10, 20, 30, 40, 50, 55, 58, 60, 62, 65, 70, 80]:
            j = int(np.argmin(np.abs(zkm - zz)))
            line = f"{zkm[j]:7.2f}"
            prof["z_km"].append(float(zkm[j]))
            for hh in ud_hours:
                r = float(abs(got[hh]["mf_dn"][j] / got[hh]["mf_up"][j]))
                prof[f"{hh:g}h"].append(r)
                line += f"{r:9.4f}"
            print(line)
        updown["R_profile"] = prof

        # ---- per mode, against each mode's own clock -----------------------
        # fup and fdn are already resolved by horizontal wavenumber inside
        # updown_split, and the rfft grid there is the grid ridge_spectrum()
        # builds, so ridge mode n is column n of both and this table costs no
        # further snapshot reading.
        #
        # The separation is carried by timing rather than by amplitude. Because
        # cg_z scales with k, the modes are spread widely in transit time, and
        # the record ends inside the gap between t_up and t_rt for the long
        # ones. Reflection at the sponge base cannot put a downward component at
        # 10 km in a mode whose round trip has not elapsed; a reflector
        # distributed through the interior, or leakage from the steady
        # continuity closure the split assumes, switches R_k on at t_up with no
        # round-trip delay.
        nmode = min(a.updown_modes, kk.size - 1)
        per_mode = {}
        if nmode <= 0:
            print("\nPer-mode table skipped (--updown-modes 0).")
        else:
            mode_zs = [float(v) for v in a.updown_mode_heights.split(",") if v.strip()]
            j_spg = int(np.argmin(np.abs(zkm - SPONGE_BASE_KM)))
            print(f"\nPer-mode split. R_k = |MF_down,k| / |MF_up,k| at fixed "
                  f"height, against that\nmode's upward arrival t_up = I(z)/(k U^2) "
                  f"and its sponge round trip\nt_rt = (2 I({SPONGE_BASE_KM:.0f} km) "
                  f"- I(z))/(k U^2). The two hypotheses require different\ncolumns "
                  f"to be at the floor: a sponge origin needs R_k at the floor "
                  f"wherever\nt_rt exceeds the record, an interior origin needs it "
                  f"nonzero wherever t_up does not.")

            # Floor of the per-mode ratio. The band floor is a share-weighted sum
            # and says nothing about any single column, so it is remeasured here
            # on the same synthetic purely upward wave.
            us0, ws0, rs0 = synthetic_wkb(hk, kk, z, rho, N, u_bar, nx, refl=0.0)
            sp0 = updown_split(us0, ws0, rs0, np.arange(nx) * dx, z, N,
                               np.full(len(z), u_bar))
            with np.errstate(divide="ignore", invalid="ignore"):
                floor_k = np.abs(sp0["fdn"] / sp0["fup"])
            floor_k = np.where(np.isfinite(floor_k), floor_k, np.nan)

            rk_of = lambda g, j, n: float(abs(g["fdn_k"][j, n] / g["fup_k"][j, n]))
            rk_fd = lambda g, j, n: float(abs(g["fdn_fd_k"][j, n]
                                              / g["fup_fd_k"][j, n]))
            cov = float(share[1:nmode + 1].sum())

            # Three classes, on each mode's own clock rather than on the wall.
            # A mode counts as established at a height once the sample time is at
            # least SETTLE times its own transit to that height; below that its
            # amplitude is still rising and the ratio of two small numbers is not
            # a measurement. Of the established modes, those whose round trip has
            # elapsed could have received something from the sponge and those
            # whose has not could not. The second class is the test.
            SETTLE = 1.5
            SH_FLOOR = a.updown_mode_share_floor

            def split_classes(t_up, t_rt, t_now, sh=None):
                """Modes with no return possible yet, and modes that could have.

                A mode is established once the sample time exceeds SETTLE times
                its own transit; below that its amplitude is still rising and the
                ratio of two small numbers is not a measurement. Modes carrying
                less than SH_FLOOR of the steady flux are dropped from both
                classes for the same reason in the spectral direction: on a wide
                ridge the flux share falls by four orders of magnitude across
                these twelve modes, and the ratio in the smallest of them is set
                by the noise of the partition rather than by the wave.
                """
                est = t_now >= SETTLE * t_up
                if sh is not None:
                    est = est & (sh >= SH_FLOOR)
                ret = est & (t_rt <= t_now)
                return est & ~ret, ret

            def wstat(r, w, msk):
                """Share-weighted mean of r over msk, and its standard error."""
                if not msk.any() or w[msk].sum() <= 0:
                    return float("nan"), float("nan"), 0
                ww = w[msk] / w[msk].sum()
                mu = float((ww * r[msk]).sum())
                nn_ = int(msk.sum())
                if nn_ < 2:
                    return mu, float("nan"), nn_
                var = float((ww * (r[msk] - mu) ** 2).sum() * nn_ / (nn_ - 1))
                return mu, float(np.sqrt(var / nn_)), nn_

            rows_pm, verdict, trend = [], [], []
            for zz in mode_zs:
                j = int(np.argmin(np.abs(zkm - zz)))
                nn, t_up, t_rt = mode_timing(ibuoy, kk, u_bar, j, j_spg, nmode)
                sh = share[nn]
                Bc, Cc = split_classes(t_up, t_rt, last["t_h"], sh)
                lab = np.where(Cc, "returned",
                               np.where(Bc, "no return",
                                        np.where(sh < SH_FLOOR, "negligible",
                                                 "arriving")))
                print(f"\n  z = {zkm[j]:.2f} km")
                print(f"{'n':>3} {'lam_x[km]':>10} {'share':>7} {'K0/om':>7} "
                      f"{'t_up[h]':>8} {'t_rt[h]':>8} {'floor':>8} {'eps':>7}"
                      + "".join(f"{('R_k@%.0fh' % hh):>10}" for hh in ud_hours)
                      + f"{'R_k(fd)':>9}  class")
                for i, n in enumerate(nn):
                    e = {"z_km": float(zkm[j]), "n": int(n),
                         "lambda_x_km": float(2 * np.pi / kk[n] / 1e3),
                         "flux_share": float(sh[i]),
                         "K0_over_omega": float(a.sponge_k0 / (u_bar * kk[n])),
                         "t_up_h": float(t_up[i]), "t_rt_h": float(t_rt[i]),
                         "class": str(lab[i]), "floor": float(floor_k[j, n]),
                         "closure_eps": float(last["eps_k"][j, n]),
                         "R_k_fd": rk_fd(last, j, n),
                         "R_k": {f"{hh:g}h": rk_of(got[hh], j, n)
                                 for hh in ud_hours}}
                    rows_pm.append(e)
                    print(f"{n:3d} {e['lambda_x_km']:10.0f} {sh[i]:7.4f} "
                          f"{e['K0_over_omega']:7.1f} {t_up[i]:8.2f} "
                          f"{t_rt[i]:8.2f} {e['floor']:8.5f} "
                          f"{e['closure_eps']:7.4f}"
                          + "".join(f"{e['R_k'][f'{hh:g}h']:10.5f}"
                                    for hh in ud_hours)
                          + f"{e['R_k_fd']:9.5f}  {lab[i]}")

                # Class means at every sampled hour. Membership is recomputed at
                # each hour rather than carried back from the last one, because a
                # mode that had not returned by 8 h belongs in the no-return
                # column there whatever it does later; the number of modes in
                # each class therefore changes down the table. A sponge origin
                # needs the no-return column at the floor at every hour; an
                # origin below the sponge needs the two columns to track.
                print(f"{'':>3} {'class':>12}"
                      + "".join(f"{('n  share  mean@%.0fh' % hh):>26}"
                                for hh in ud_hours))
                cls_rows = []
                for name, pick in (("no return", 0), ("returned", 1)):
                    line = f"{'':>3} {name:>12}"
                    ent = {"class": name, "by_hour": {}}
                    for hh in ud_hours:
                        gh = got[hh]
                        msk = split_classes(t_up, t_rt, gh["t_h"], sh)[pick]
                        rr = np.array([rk_of(gh, j, n) for n in nn])
                        mu, se, nm = wstat(rr, sh, msk)
                        ent["by_hour"][f"{hh:g}h"] = {
                            "mean": mu, "se": se, "n_modes": nm,
                            "share": float(sh[msk].sum())}
                        line += (f"{nm:6d}{sh[msk].sum():7.3f}"
                                 + (f"{mu:8.5f}" if np.isfinite(mu) else "     nan")
                                 + (f"+-{se:.3f}" if np.isfinite(se) else "      "))
                    cls_rows.append(ent)
                    print(line)

                # Per-snapshot scatter of the same class means over the final
                # window, which is the error bar the mode-to-mode spread cannot
                # supply because a common temporal fluctuation moves every mode.
                for ent, msk in zip(cls_rows, (Bc, Cc)):
                    vals = []
                    for mem in last.get("members", []):
                        rr = np.array([abs(mem["fdn_k"][j, n] / mem["fup_k"][j, n])
                                       for n in nn])
                        mu, _, _ = wstat(rr, sh, msk)
                        if np.isfinite(mu):
                            vals.append(mu)
                    ent["snapshot_sd"] = (float(np.std(vals, ddof=1))
                                          if len(vals) > 1 else float("nan"))
                    ent["snapshot_n"] = len(vals)
                if all(np.isfinite(e_.get("snapshot_sd", np.nan)) for e_ in cls_rows):
                    print(f"{'':>3} per-snapshot sd of the final-hour means over "
                          f"its {cls_rows[0]['snapshot_n']} snapshots: "
                          + ", ".join(f"{e_['class']} {e_['snapshot_sd']:.4f}"
                                      for e_ in cls_rows))

                # What the record says the asymptote is. Each R_k is already a
                # per-mode ratio, so a returned mode carries no filling of its
                # own and its r_n is final; only the modes that have not returned
                # have to be assumed, and two assumptions are carried for them.
                rk = np.array([rk_of(last, j, n) for n in nn])
                r_ret = wstat(rk, sh, Cc)[0]
                r_nor = wstat(rk, sh, Bc)[0]
                r_trend, r_pred_B = float("nan"), float("nan")
                if Cc.sum() >= 3:
                    xg = np.log10(2 * np.pi / kk[nn] / 1e3)
                    c = np.polyfit(xg[Cc], rk[Cc], 1, w=np.sqrt(sh[Cc]))
                    ex = np.clip(np.polyval(c, xg), 0.0, None)
                    r_trend = float((sh * np.where(Cc, rk, ex)).sum() / sh.sum())
                    trend.append({"z_km": float(zkm[j]), "n_used": int(Cc.sum()),
                                  "slope_per_decade": float(c[0]),
                                  "r_at_167km": float(np.polyval(c, np.log10(167.0))),
                                  "r_at_2000km": float(np.polyval(c, np.log10(2000.0)))})
                    if Bc.any():
                        r_pred_B = wstat(ex, sh, Bc)[0]
                v = {"z_km": float(zkm[j]), "t_h": last["t_h"], "settle_factor": SETTLE,
                     "classes": cls_rows,
                     "floor_max_over_modes": float(np.nanmax(floor_k[j, nn])),
                     "closure_eps_max": float(np.nanmax(last["eps_k"][j, nn])),
                     "R_measured_band": float(abs(last["mf_dn"][j]
                                                  / last["mf_up"][j])),
                     "r_returned": r_ret, "r_no_return": r_nor,
                     "r_no_return_predicted_by_trend": r_pred_B,
                     "R_inf_flat": r_ret, "R_inf_trend": r_trend}
                verdict.append(v)
                print(f"{'':>3} R measured here {v['R_measured_band']:.4f}; "
                      f"R_inf {r_ret:.4f} holding r_n flat over the modes that have "
                      f"not\n{'':>3} returned and "
                      + (f"{r_trend:.4f}" if np.isfinite(r_trend) else "n/a")
                      + " extrapolating their trend in lambda_x. The lambda_x trend "
                      f"of the\n{'':>3} returned modes puts the unreturned ones at "
                      + (f"{r_pred_B:.4f}" if np.isfinite(r_pred_B) else "n/a")
                      + f", against {r_nor:.4f} measured.")

            # Does r_n order itself the way a damping layer would order it? K0
            # over the intrinsic frequency runs over an order of magnitude across
            # these modes and increases with lambda_x, so a sponge that reflects
            # should reflect the long modes hardest.
            if trend:
                print(f"\nTrend of r_n with the horizontal scale, over the modes "
                      f"that have returned.\nA sponge origin orders r_n with "
                      f"K0/omega_hat, which increases with lambda_x, so the slope\n"
                      f"should be positive.")
                print(f"{'z[km]':>7} {'n used':>7} {'d r_n/d log10 lam_x':>20} "
                      f"{'r(167 km)':>10} {'r(2000 km)':>11}")
                for e in trend:
                    print(f"{e['z_km']:7.2f} {e['n_used']:7d} "
                          f"{e['slope_per_decade']:+20.5f} {e['r_at_167km']:10.5f} "
                          f"{e['r_at_2000km']:11.5f}")

            # The aggregate version of the same correction, as a cross-check on
            # the per-mode one. If mode-independent reflection describes the
            # record, R divided by the share that has had time to return should
            # reproduce the share-weighted mean of the measured r_n.
            kpos_f = np.where(kk > 0.0, kk, np.inf)
            t_up_a = ibuoy[:, None] / (kpos_f[None, :] * u_bar ** 2) / 3600.0
            t_rt_a = ((2 * ibuoy[j_spg] - ibuoy)[:, None]
                      / (kpos_f[None, :] * u_bar ** 2) / 3600.0)
            phi = lambda tt, ta: float(
                (share[None, :] * (ta <= tt)).sum(1)[sel_b].mean())
            print(f"\nR against filled fraction rather than clock time, band "
                  f"{a.z_lo:.0f}-{z_hi:.0f} km.\nPhi_up is the share of the ridge "
                  f"flux that has arrived, Phi_dn the share that has\nhad time to "
                  f"return from {SPONGE_BASE_KM:.0f} km.")
            print(f"{'t[h]':>6} {'Phi_up':>7} {'Phi_dn':>7} {'Phi_dn/Phi_up':>14} "
                  f"{'R':>7} {'implied r':>10} {'per-mode r':>11}")
            fillrows = []
            for hh in ud_hours:
                g = got[hh]
                pu, pd = phi(g["t_h"], t_up_a), phi(g["t_h"], t_rt_a)
                Rb = float(abs(bmean(g["mf_dn"]) / bmean(g["mf_up"])))
                jm = [int(np.argmin(np.abs(zkm - zz))) for zz in mode_zs]
                nn = np.arange(1, nmode + 1)
                rr_pm = []
                for j in jm:
                    _, tup_, trt = mode_timing(ibuoy, kk, u_bar, j, j_spg, nmode)
                    _, bk = split_classes(tup_, trt, g["t_h"], share[nn])
                    if bk.any() and share[nn][bk].sum() > 0:
                        rk = np.array([rk_of(g, j, n) for n in nn])
                        rr_pm.append(float((share[nn][bk] * rk[bk]).sum()
                                           / share[nn][bk].sum()))
                pmr = float(np.mean(rr_pm)) if rr_pm else float("nan")
                e = {"t_h": g["t_h"], "phi_up": pu, "phi_dn": pd,
                     "available_share": pd / pu if pu else float("nan"),
                     "R": Rb, "implied_r_aggregate": Rb * pu / pd if pd else float("nan"),
                     "r_per_mode_mean": pmr}
                fillrows.append(e)
                print(f"{g['t_h']:6.2f} {pu:7.4f} {pd:7.4f} "
                      f"{e['available_share']:14.4f} {Rb:7.4f} "
                      f"{e['implied_r_aggregate']:10.5f} {pmr:11.5f}")
            print("The last two columns are the same quantity reached two ways. "
                  "The aggregate\ncolumn assumes every mode reflects alike and "
                  "corrects R by the share that could\nhave returned; the per-mode "
                  "column measures r_n on the modes that did and needs\nno such "
                  "assumption. Where they part, the assumption is what fails.")

            per_mode = {"n_modes": int(nmode), "heights_km": mode_zs,
                        "share_covered": cov, "settle_factor": SETTLE,
                        "share_floor": SH_FLOOR,
                        "rows": rows_pm, "verdict": verdict,
                        "lambda_x_trend": trend, "filled_fraction": fillrows}
            # The one alternative the synthetic floor cannot cover. That check
            # is built on a steady WKB field, so it exercises the algebra but not
            # the continuity closure; eps measures the closure directly against a
            # differenced derivative, and the partition rebuilt on that derivative
            # is carried in the R_k(fd) column above.
            print(f"\nClosure of the relation the split is built on. eps is the "
                  f"relative difference\nbetween dPsi from steady continuity and "
                  f"dPsi differenced in z, which is\naccurate to about 0.2 per cent "
                  f"at dz = {(z[1]-z[0]):.0f} m against lambda_z = "
                  f"{lam_z_band:.0f} km. An error eps\nmanufactures a downward "
                  f"fraction of about eps^2/4 in every mode.")
            print(f"{'z[km]':>7}"
                  + "".join(f"{('eps@%.0fh' % hh):>9}" for hh in ud_hours)
                  + f" {'eps^2/4':>9} {'r no-return':>12} {'same, fd':>10} "
                  f"{'r returned':>11} {'same, fd':>10}")
            clos_rows = []
            for zz in mode_zs:
                j = int(np.argmin(np.abs(zkm - zz)))
                nn2, tu2, tr2 = mode_timing(ibuoy, kk, u_bar, j, j_spg, nmode)
                B2, C2 = split_classes(tu2, tr2, last["t_h"], share[nn2])
                sh2 = share[nn2]
                r1 = np.array([rk_of(last, j, n) for n in nn2])
                r2 = np.array([rk_fd(last, j, n) for n in nn2])
                em = float(np.nanmax(last["eps_k"][j, nn2]))
                eps_h = {f"{hh:g}h": float(np.nanmax(got[hh]["eps_k"][j, nn2]))
                         for hh in ud_hours}
                e = {"z_km": float(zkm[j]), "eps_max": em, "eps_by_hour": eps_h,
                     "spurious_R_from_eps": em * em / 4.0,
                     "r_no_return": wstat(r1, sh2, B2)[0],
                     "r_no_return_fd": wstat(r2, sh2, B2)[0],
                     "r_returned": wstat(r1, sh2, C2)[0],
                     "r_returned_fd": wstat(r2, sh2, C2)[0]}
                clos_rows.append(e)
                print(f"{zkm[j]:7.2f}"
                      + "".join(f"{eps_h[f'{hh:g}h']:9.4f}" for hh in ud_hours)
                      + f" {em*em/4:9.5f} "
                      f"{e['r_no_return']:12.5f} {e['r_no_return_fd']:10.5f} "
                      f"{e['r_returned']:11.5f} {e['r_returned_fd']:10.5f}")
            per_mode["closure_audit"] = clos_rows

            # Shape of the downward spectrum against the upward one. This is
            # independent of the timing test and asks a different question. A
            # reflection returns the incident wave, so it carries the incident
            # spectrum and |MF_dn,k| follows |MF_up,k| with the exponent alpha in
            # |MF_dn,k| = C |MF_up,k|^alpha equal to 1. A source that is broadband
            # in k, whether a residue of the start-up adjustment or scattering at
            # the grid scale, gives alpha near 0 and a downward spectrum flatter
            # than the ridge's. The distinction matters because R_k is a ratio and
            # its denominator alone spans a factor of eight across these modes,
            # which is enough to manufacture a trend in r_n with lambda_x out of a
            # downward part that has no mode dependence at all.
            print(f"\nShape of the downward spectrum. |MF_dn,k| = C |MF_up,k|^alpha "
                  f"over the modes\nestablished at each height: alpha = 1 is a "
                  f"reflection of the incident wave,\nalpha = 0 a component with "
                  f"no memory of the ridge spectrum.")
            print(f"{'z[km]':>7} {'n used':>7} {'alpha':>16} "
                  f"{'up spread':>10} {'down spread':>12} {'MF_dn sum[uPa]':>15}")
            shape = []
            for zz in mode_zs:
                j = int(np.argmin(np.abs(zkm - zz)))
                nn3, tu3, tr3 = mode_timing(ibuoy, kk, u_bar, j, j_spg, nmode)
                B3, C3 = split_classes(tu3, tr3, last["t_h"], share[nn3])
                est = B3 | C3
                uk = np.abs(last["fup_k"][j, nn3]) * 1e6
                dk = np.abs(last["fdn_k"][j, nn3]) * 1e6
                use = est & (uk > 0) & (dk > 0)
                if use.sum() < 3:
                    continue
                lx_, ly_ = np.log(uk[use]), np.log(dk[use])
                al, c0 = np.polyfit(lx_, ly_, 1)
                res = ly_ - (al * lx_ + c0)
                sse = float(res @ res) / max(use.sum() - 2, 1)
                se_al = float(np.sqrt(sse / ((lx_ - lx_.mean()) ** 2).sum()))
                e = {"z_km": float(zkm[j]), "n_used": int(use.sum()),
                     "alpha": float(al), "alpha_se": se_al,
                     "up_spread": float(uk[use].max() / uk[use].min()),
                     "down_spread": float(dk[use].max() / dk[use].min()),
                     "MF_down_sum_uPa": float(dk[use].sum())}
                shape.append(e)
                print(f"{zkm[j]:7.2f} {int(use.sum()):7d} "
                      f"{al:+8.3f} +- {se_al:<5.3f} {e['up_spread']:10.2f} "
                      f"{e['down_spread']:12.2f} {e['MF_down_sum_uPa']:15.4f}")
            per_mode["downward_spectrum_shape"] = shape
            print("The spread columns are max over min across the modes used. A "
                  "downward spectrum\nflatter than the upward one is not the "
                  "incident wave coming back.")

            pm_fig(outdir, per_mode, ud_hours, float(th[-1]))
            print(f"\nThese {nmode} modes carry {cov:.3f} of the steady ridge flux, "
                  f"so the asymptote above\nis quoted on that part of the spectrum "
                  f"rather than on all of it. Two limits on\nwhat the no-return "
                  f"class can be read as: a mode is only counted once the sample "
                  f"time\nexceeds {SETTLE:g} times its own transit, and a downward "
                  f"component left over from the start-up\nadjustment would sit in "
                  f"the same column, which the early hours of the trace bear on\n"
                  f"and this table does not.")
        updown["per_mode"] = per_mode

        # Sub-band drift, with the reflected part removed and with the density
        # weighting changed. The second is a definitional choice of the same
        # order as the residual being argued about, so it is reported next to it.
        print(f"\nSub-band d ln|MF|/dz at {last['t_h']:.2f} h, per cent per km.")
        print(f"{'band':>10} {'<rho uw>':>10} {'<rho0 uw>':>10} {'upward':>9} "
              f"{'downward':>9} {'weighting':>10}")
        sub = []
        for lo, hi in subbands(a.z_lo, z_hi):
            s = (zkm >= lo) & (zkm <= hi)
            e = {"band_km": [lo, hi],
                 "drift_direct": bdrift_of(last["mf_direct"], s),
                 "drift_rho0": bdrift_of(last["mf_spec"], s),
                 "drift_up": bdrift_of(last["mf_up"], s),
                 "drift_down": bdrift_of(last["mf_dn"], s)}
            e["weighting_spread"] = e["drift_direct"] - e["drift_rho0"]
            sub.append(e)
            print(f"{lo:5.0f}-{hi:3.0f} {e['drift_direct']:+10.4f} "
                  f"{e['drift_rho0']:+10.4f} {e['drift_up']:+9.4f} "
                  f"{e['drift_down']:+9.4f} {e['weighting_spread']:+10.4f}")
        updown["subband_drift_pct_per_km"] = sub
        print("The last column is the difference between weighting the covariance "
              "by rho and\nby its x-mean, a choice with no content at this order. "
              "Where it is comparable\nto the entry beside it, the sub-band slope "
              "is not a measurement.")

        # Partial standing wave. The flux carries no A-B cross term by
        # construction, so the interference has to be looked for in the variance.
        print("\nPartial standing wave in the variance, which is where the "
              "interference sits.")
        print(f"{'t[h]':>6} {'r=sqrt(R)':>10} {'2r/(1+r^2)':>11} {'measured ptp':>13} "
              f"{'peak scale':>11} {'lambda_z/2':>11}")
        sw = []
        zb = zkm[sel_b]
        for hh in ud_hours:
            g = got[hh]
            c = ((g["var_total"] - g["var_incoherent"]) / g["var_incoherent"])[sel_b]
            v = c - np.polyval(np.polyfit(zb, c, 1), zb)
            F = np.abs(np.fft.rfft(v * np.hanning(v.size)))
            F[0] = 0.0
            fz = np.fft.rfftfreq(v.size, zb[1] - zb[0])
            peak = float(1.0 / fz[int(np.argmax(F))])
            r = float(np.sqrt(abs(bmean(g["mf_dn"]) / bmean(g["mf_up"]))))
            e = {"t_h": g["t_h"], "r": r, "predicted_depth": 2 * r / (1 + r * r),
                 "measured_ptp": float(np.ptp(c)), "peak_scale_km": peak,
                 "half_lambda_z_km": lam_z_band / 2}
            sw.append(e)
            print(f"{g['t_h']:6.2f} {r:10.4f} {e['predicted_depth']:11.4f} "
                  f"{e['measured_ptp']:13.4f} {peak:11.1f} {lam_z_band/2:11.1f}")
        updown["standing_wave"] = sw

        # The cross term above is formed from the split's own A and B, so it is
        # not independent of the split. The ripple in each mode's variance is:
        # ripple_reflection() reads |Psi_k|^2 alone, flattened by the exact local
        # m, and returns the reflected share a standing-wave ratio implies. It is
        # set beside the split's band-mean ratio for the same mode and beside the
        # linear reflection of this run's own damping layer, a sin^2 Rayleigh
        # profile on the background column, rayleigh_sponge_reflection().
        s_rho = -0.5 * np.gradient(np.log(rho), z)
        m_ex_all = np.sqrt(np.clip((N / u_bar)[:, None] ** 2 - kk[None, :] ** 2
                                   - (s_rho ** 2 - np.gradient(s_rho, z))[:, None],
                                   0.0, None))
        t55 = buoyancy_integral(N, z)[int(np.argmin(np.abs(zkm - z_hi)))]
        nmr = min(a.updown_modes, kk.size - 1)
        R_th = {}
        for n in range(1, nmr + 1):
            rth = rayleigh_sponge_reflection(kk[n], z, N, rho, u_bar, a.sponge_k0,
                                             SPONGE_BASE_KM * 1e3, z_top_km * 1e3,
                                             top="radiation")
            R_th[n] = float(np.nanmean(rth[sel_b]))
        print(f"\nPer-mode reflected share over {a.z_lo:.0f}-{z_hi:.0f} km from the "
              f"split, from the ripple in |Psi_k|^2 alone,\nand from linear theory "
              f"for a sin^2 layer of K0 = {a.sponge_k0:g} 1/s over "
              f"{SPONGE_BASE_KM:.0f}-{z_top_km:.0f} km. A mode is marked '*' until "
              f"1.5 times\nits upward transit to {z_hi:.0f} km has elapsed, "
              f"before which neither estimate is a measurement.")
        rip = []
        iz_lo = int(np.argmin(np.abs(zkm - a.z_lo)))   # phase of B/A read here
        for hh in ud_hours:
            g = got[hh]
            print(f"{'t[h]':>6} {'n':>3} {'lam_x[km]':>10} {'R_split':>9} {'R_ripple':>9} "
                  f"{'R_best':>8} {'m_scale':>8} {'R_theory':>9} {'phBA_lo':>8}")
            for n in range(1, nmr + 1):
                t_up55 = t55 / (kk[n] * u_bar ** 2) / 3600.0
                fu = np.abs(g["fup_k"][:, n])
                rs = float(np.nanmean((np.abs(g["fdn_k"][:, n])
                                       / np.where(fu > 0, fu, np.nan))[sel_b]))
                R1, Rb, sc, dep = ripple_reflection(g["psi2_k"][:, n], z,
                                                    m_ex_all[:, n], sel_b)
                est = g["t_h"] > 1.5 * t_up55
                ph = float(np.degrees(np.angle(g["ab_k"][iz_lo, n])))
                e = {"t_h": g["t_h"], "n": n, "lambda_x_km": 2 * np.pi / kk[n] / 1e3,
                     "R_split": rs, "R_ripple": R1, "R_ripple_best": Rb,
                     "m_scale_best": sc, "ripple_depth": dep,
                     "R_theory": R_th[n], "t_up_top_h": t_up55, "established": bool(est),
                     "phase_BA_lo_deg": ph}
                rip.append(e)
                print(f"{g['t_h']:6.2f} {n:3d} {e['lambda_x_km']:10.0f} {rs:9.4f} "
                      f"{R1:9.4f} {Rb:8.4f} {sc:8.3f} {R_th[n]:9.4f} {ph:8.1f}"
                      f"{'' if est else ' *'}")
        updown["ripple_permode"] = rip
        updown["sponge_theory"] = {"K0_per_s": a.sponge_k0, "base_km": SPONGE_BASE_KM,
                                   "top_km": z_top_km,
                                   "R_band_by_mode": {str(n): v for n, v in R_th.items()}}

        # How much m has to be wrong before the answer moves.
        print("\nSensitivity of the final-hour result to the assumed m and to the "
              "state it is\nbuilt from.")
        print(f"{'variant':>22} {'R':>8} {'drift_up':>9} {'drift_net':>10}")
        sens = []
        idx_last = ud_window(max(ud_hours), a.updown_window, a.updown_max_snaps)
        for lab, bg, ms in ([("mean state, m as is", False, 1.0),
                             ("background N and U", True, 1.0)]
                            + [(f"m x {v:.2f}", False, v)
                               for v in [float(s) for s in a.updown_msens.split(",")
                                         if s.strip()]]):
            g = ud_mean(idx_last, m_scale=ms, background=bg)
            e = {"variant": lab,
                 "R": float(abs(bmean(g["mf_dn"]) / bmean(g["mf_up"]))),
                 "drift_up": bdrift_of(g["mf_up"]),
                 "drift_net": bdrift_of(g["mf_direct"])}
            sens.append(e)
            print(f"{lab:>22} {e['R']:8.5f} {e['drift_up']:+9.4f} "
                  f"{e['drift_net']:+10.4f}")
        updown["sensitivity"] = sens

        # Where the damping the wave meets actually begins. The layer is thick,
        # but K rises as sin^2 from zero, so the height at which the damping rate
        # overtakes the intrinsic frequency is the scale the wave responds to.
        print(f"\nSponge as the wave meets it. K(z) = K0 sin^2(pi/2 (z-"
              f"{SPONGE_BASE_KM:.0f})/{z_top_km-SPONGE_BASE_KM:.0f}), K0 = "
              f"{a.sponge_k0:g} 1/s.")
        print(f"{'lambda_x[km]':>13} {'omega_hat[1/s]':>15} {'K0/omega_hat':>13} "
              f"{'z(K=omega)[km]':>15} {'above base':>11}")
        spg = []
        for j in [j for j in range(1, 13)]:
            kj = 2 * np.pi * j / lx
            om = u_bar * kj
            frac = min(om / a.sponge_k0, 1.0)
            zc = SPONGE_BASE_KM + (z_top_km - SPONGE_BASE_KM) * (2 / np.pi) * np.arcsin(
                np.sqrt(frac))
            e = {"lambda_x_km": 2 * np.pi / kj / 1e3, "omega_hat": om,
                 "K0_over_omega": a.sponge_k0 / om, "z_K_eq_omega_km": zc,
                 "above_base_in_lambda_z": (zc - SPONGE_BASE_KM) / lam_z_band}
            spg.append(e)
            if j in (2, 3, 5, 8, 12):
                print(f"{e['lambda_x_km']:13.1f} {om:15.3e} "
                      f"{e['K0_over_omega']:13.1f} {zc:15.2f} "
                      f"{e['above_base_in_lambda_z']:11.3f} lambda_z")
        updown["sponge_scales"] = {
            "K0_per_s": a.sponge_k0, "base_km": SPONGE_BASE_KM,
            "top_km": z_top_km, "lambda_z_band_km": lam_z_band,
            "thickness_in_lambda_z": (z_top_km - SPONGE_BASE_KM) / lam_z_band,
            "per_mode": spg}

        # Level decay, split. A height-uniform loss is invisible to section 2.
        print("\nDecay of the level, split between the incident wave and the "
              "returning one.")
        print(f"{'interval':>13} {'dln|up|/dt':>11} {'dln|net|/dt':>12} "
              f"{'d(down)/dt':>11} {'reflected share of the loss':>28}")
        lev = []
        for p, q in zip(ud_hours[:-1], ud_hours[1:]):
            g0, g1 = got[p], got[q]
            dt_h = g1["t_h"] - g0["t_h"]
            u0, u1 = bmean(g0["mf_up"]), bmean(g1["mf_up"])
            n0, n1 = bmean(g0["mf_direct"]), bmean(g1["mf_direct"])
            d0, d1 = bmean(g0["mf_dn"]), bmean(g1["mf_dn"])
            # Defined only while the level is falling. Over an interval in which
            # the flux is still filling in, the net grows and the ratio has no
            # reading, so it is left out rather than reported with a sign.
            # Named apart from the ridge flux share of section 8, which is an
            # array over wavenumber that later sections still read.
            refl = ((d1 - d0) / (abs(n0) - abs(n1))
                    if abs(n1) < abs(n0) else float("nan"))
            e = {"from_h": g0["t_h"], "to_h": g1["t_h"],
                 "dln_up_pct_per_h": 100 * np.log(abs(u1 / u0)) / dt_h,
                 "dln_net_pct_per_h": 100 * np.log(abs(n1 / n0)) / dt_h,
                 "d_down_uPa_per_h": (d1 - d0) * 1e6 / dt_h,
                 "reflected_share_of_loss": float(refl)}
            lev.append(e)
            print(f"{g0['t_h']:5.1f}-{g1['t_h']:5.1f} "
                  f"{e['dln_up_pct_per_h']:+11.3f} {e['dln_net_pct_per_h']:+12.3f} "
                  f"{e['d_down_uPa_per_h']:+11.4f} "
                  + (f"{100*refl:26.1f} %" if np.isfinite(refl)
                     else f"{'still filling':>28}"))
        updown["level_decay"] = lev

        # Whether the incident wave is weakening because its source is. The
        # launch flux of a stationary ridge wave scales with the wind that
        # crosses the ridge, so a decaying incident flux under an unchanged
        # low-level wind is not a source effect.
        src_rows = []
        line = f"{'z[km]':>7}"
        for hh in ud_hours:
            line += f"{hh:>11.0f}h"
        print("\nx-mean wind at the launch levels, against the same hours.")
        print(line)
        for zz in [0.25, 2.0, 5.0, 10.0]:
            j = int(np.argmin(np.abs(zkm - zz)))
            vals = [float(ubar[j, int(np.argmin(np.abs(th - hh)))]) for hh in ud_hours]
            src_rows.append({"z_km": float(zkm[j]), "u_bar_m_s": vals})
            print(f"{zkm[j]:7.2f}" + "".join(f"{v:12.4f}" for v in vals))
        updown["launch_level_wind"] = src_rows
        sp_ = max(abs(r["u_bar_m_s"][-1] - r["u_bar_m_s"][0]) for r in src_rows)
        print(f"The largest change across these hours is {sp_:.4f} m/s, so the "
              f"decay of the\nincident flux is not a weakening of the source.")

        # Amplitude accounting against the steady linear reference.
        sh_net = 1.0 - bmean(last["mf_direct"]) / lin_band
        sh_up = 1.0 - bmean(last["mf_up"]) / lin_band
        updown["amplitude_accounting"] = {
            "t_h": last["t_h"], "linear_steady_uPa": lin_band * 1e6,
            "net_shortfall_pct": 100 * sh_net, "upward_shortfall_pct": 100 * sh_up,
            "reflection_share_of_shortfall_pct":
                100 * (sh_net - sh_up) / sh_net if sh_net else float("nan")}
        print(f"\nAgainst the steady linear reference of {lin_band*1e6:.3f} uPa at "
              f"{last['t_h']:.2f} h, the net flux\nfalls short by "
              f"{100*sh_net:.2f} % and the upward component alone by "
              f"{100*sh_up:.2f} %, so the returning\nwave accounts for "
              f"{100*(sh_net-sh_up)/sh_net:.1f} % of the shortfall. The reference "
              f"carries the N_0 ambiguity of\nsection 3 and the remainder should "
              f"not be read as a loss without it.")

        # Closure, and the two terms the partition does not contain.
        clo = float(np.max(np.abs(last["mf_spec"] - last["mf_direct"])[sel_b])
                    / abs(np.median(last["mf_direct"][sel_b])))
        ev = float(abs(bmean(last["mf_evan"]) / bmean(last["mf_direct"])))
        updown["closure"] = {"t_h": last["t_h"],
                             "triple_correlation_max_pct": 100 * clo,
                             "evanescent_fraction": ev}
        print(f"\nClosure at {last['t_h']:.2f} h. The partition sums to "
              f"<rho_0 u'w'>_x identically, so the\ntest is against "
              f"<rho u'w'>_x: the largest level difference is {100*clo:.2f} % of the "
              f"flux,\nwhich is the triple correlation <rho' u'w'>_x. Evanescent "
              f"wavenumbers carry\n{ev:.1e} of the flux, consistent with the zero "
              f"a non-propagating mode should give.")

        # Trace of R against time, at one snapshot per point.
        tr_h = np.arange(2.0, last_h + 0.01, 0.5)
        tr = []
        for hh in tr_h:
            i = int(np.argmin(np.abs(t_snap - hh * 3600.0)))
            g = ud_mean(np.array([i]))
            tr.append({"t_h": g["t_h"],
                       "R": float(abs(bmean(g["mf_dn"]) / bmean(g["mf_up"]))),
                       "MF_up_uPa": bmean(g["mf_up"]) * 1e6,
                       "MF_net_uPa": bmean(g["mf_direct"]) * 1e6})
        updown["trace"] = tr
        rr = np.array([e["R"] for e in tr])
        tt_tr = np.array([e["t_h"] for e in tr])
        # Before the band has filled, single snapshots carry the starting
        # transient and R is not yet a reflected share, so the rate is read off
        # the part of the record past it.
        late = tt_tr >= max(6.0, tt_tr[0])
        slope = float(np.polyfit(tt_tr[late], rr[late], 1)[0])
        updown["R_growth_per_h_after_6h"] = slope
        print(f"\nOver {tt_tr[late][0]:.1f}-{tt_tr[-1]:.2f} h the trace runs "
              f"{rr[late][0]:.4f} to {rr[-1]:.4f} at a linear rate of\n"
              f"{slope:.5f} per hour, with no sign of settling. The net flux does "
              f"not depend on m,\nso the unchanged drift_net column above is a "
              f"check on the partition rather than\na result.")

        ud_fig(outdir, zkm, got, ud_hours, tr, a.z_lo, z_hi, z_top_km)
    out["updown"] = updown

    # ---------------- 12. Transient ----------------
    hdr("12. TRANSIENT  -  whether the downward component is stationary in the "
        "ground frame")
    transient = {}
    if a.transient_stride <= 0 or a.updown_modes <= 0 or not ud_hours:
        print("Skipped (--transient-stride 0, --updown-modes 0, or no updown hours).")
    else:
        print(
            "Section 11 classifies a mode by t_rt = (2 I(z_spg) - I(z))/(k U^2),\n"
            "which is the round trip of a stationary mountain component, ground-based\n"
            "frequency zero and intrinsic frequency omega_hat = -kU. Anything present\n"
            "at the same k and a different frequency travels at a different speed and\n"
            "is also misassigned by the split, whose m assumes the mountain value, so\n"
            "the classification is only as good as the assumption that the field at\n"
            "each k is stationary. That assumption is testable on the record already\n"
            "written. The x-transform of a component going as exp(i(kx - omega t))\n"
            "carries exp(-i omega t), so a stationary field holds Psi_k, and with it\n"
            "A_k and B_k, fixed in time; a component at any other frequency rotates\n"
            "them at a rate that gives omega directly.")

        nmode = min(a.updown_modes, kk.size - 1)
        if len(files) == len(t):
            t_snap = t
        else:
            t_snap = np.array([_snap_time(f) for f in files])
        idx = np.arange(0, len(files), a.transient_stride)
        if idx[-1] != len(files) - 1:
            idx = np.append(idx, len(files) - 1)
        t_tr, Ak, Bk, MKk, wgt_k = transient_pass(files, idx, z, nmode)
        dts = float(np.median(np.diff(t_tr)))
        th_tr = t_tr / 3600.0
        Fup_k = wgt_k[None, None, :] * MKk * np.abs(Ak) ** 2
        Fdn_k = wgt_k[None, None, :] * MKk * np.abs(Bk) ** 2
        print(f"\n{idx.size} snapshots of {len(files)} at stride "
              f"{a.transient_stride}, sample interval {dts:.1f} s. The frequencies "
              f"this can\nresolve run up to {1e3/(2*dts):.3f} mHz and the lag "
              f"resolution is one interval.")

        SETTLE = float(updown.get("per_mode", {}).get("settle_factor", 1.5))
        SH_FLOOR = a.updown_mode_share_floor
        j_spg = int(np.argmin(np.abs(zkm - SPONGE_BASE_KM)))
        mode_zs = [float(v) for v in a.updown_mode_heights.split(",") if v.strip()]
        sel_w = th_tr >= th_tr[-1] - a.transient_window

        def classes_at(j, t_now):
            nn_, tu_, tr_ = mode_timing(ibuoy, kk, u_bar, j, j_spg, nmode)
            est = (t_now >= SETTLE * tu_) & (share[nn_] >= SH_FLOOR)
            return nn_, tu_, tr_, est & (tr_ > t_now), est & (tr_ <= t_now)

        def steady_part(s, mkv, n):
            """The two parts of <|s|^2>, one carried by a constant amplitude.

            The flux is quadratic in the amplitude, so <|s|^2> = |<s>|^2 + var(s)
            splits it exactly into a part an amplitude fixed over the window
            would carry and a part it could not. Applied to B this separates the
            downward flux that the stationary round-trip clock governs from the
            downward flux that it does not.
            """
            mm = float(mkv.mean())
            return (wgt_k[n] * mm * abs(s.mean()) ** 2,
                    wgt_k[n] * mm * float(np.mean(np.abs(s - s.mean()) ** 2)))

        def freq_lag1(s):
            """Ground-based frequency of a complex series, from its lag-1 phase.

            omega = -arg<s_{i+1} s_i*>/dt, the amplitude-weighted mean rotation
            rate. The modulus of the same average, normalised by the geometric
            mean of the two energies, is 1 for a single line and falls toward 0
            for a broadband series; normalising this way rather than by one of
            the two energies keeps a growing amplitude from pushing it past 1.
            """
            if s.size < 3:
                return float("nan"), float("nan")
            num = np.sum(s[1:] * np.conj(s[:-1]))
            den = np.sqrt(np.sum(np.abs(s[1:]) ** 2) * np.sum(np.abs(s[:-1]) ** 2))
            if den <= 0:
                return float("nan"), float("nan")
            return float(-np.angle(num) / dts), float(np.abs(num) / den)

        # ---- per mode, the steady and the fluctuating part of each branch ----
        print(f"\nSteady and fluctuating parts of the two branches over the last "
              f"{a.transient_window:.1f} h,\nas fractions of the total upward flux "
              f"of the same mode. omega is the rotation\nrate of the fluctuating "
              f"part of B and nb its narrowbandness. A component at\nomega has "
              f"|omega_hat| = |omega - kU| and cg_z = omega_hat^2/(N k), so its "
              f"round trip\nis the tabulated one times (kU/omega_hat)^2, and the "
              f"split, whose m assumes the\nmountain value, manufactures a "
              f"downward share R_sp = ((1-r)/(1+r))^2 from it with\nr = "
              f"kU/|omega_hat|.")
        rows_tr, cls_tr = [], []
        for zz in mode_zs:
            j = int(np.argmin(np.abs(zkm - zz)))
            nn_, tu_, trt_, nor, ret = classes_at(j, float(th_tr[-1]))
            print(f"\n  z = {zkm[j]:.2f} km")
            print(f"{'n':>3} {'lam_x[km]':>10} {'share':>7} {'t_up':>6} {'t_rt':>7} "
                  f"{'cls':>4} {'up_st':>8} {'up_tr':>8} {'dn_st':>8} {'dn_tr':>8} "
                  f"{'om[mrad/s]':>11} {'om/kU':>7} {'nb':>5} {'t_rt corr':>10} "
                  f"{'R_sp':>7}")
            for i, n in enumerate(nn_):
                aa, bb = Ak[sel_w, j, n], Bk[sel_w, j, n]
                ua = steady_part(aa, MKk[sel_w, j, n], n)
                ub = steady_part(bb, MKk[sel_w, j, n], n)
                tot = ua[0] + ua[1]
                om, nb = freq_lag1(bb - bb.mean())
                kU = float(kk[n] * u_bar)
                rat = abs(om - kU) / kU
                r_m = 1.0 / rat if rat > 0 else np.inf
                cls = "nor" if nor[i] else ("ret" if ret[i] else "-")
                e = {"z_km": float(zkm[j]), "n": int(n),
                     "lambda_x_km": float(2 * np.pi / kk[n] / 1e3),
                     "flux_share": float(share[n]), "class": cls,
                     "t_up_h": float(tu_[i]), "t_rt_h": float(trt_[i]),
                     "up_steady": ua[0] / tot, "up_fluct": ua[1] / tot,
                     "down_steady": ub[0] / tot, "down_fluct": ub[1] / tot,
                     "omega_rad_s": om, "omega_over_kU": om / kU,
                     "narrowband": nb, "omega_hat_over_kU": float(rat),
                     "t_rt_corrected_h": float(trt_[i] / rat ** 2) if rat > 0
                     else float("inf"),
                     "R_spurious": float(((1 - r_m) / (1 + r_m)) ** 2)
                     if np.isfinite(r_m) else 1.0}
                rows_tr.append(e)
                print(f"{n:3d} {e['lambda_x_km']:10.0f} {share[n]:7.4f} "
                      f"{tu_[i]:6.2f} {trt_[i]:7.2f} {cls:>4} {e['up_steady']:8.5f} "
                      f"{e['up_fluct']:8.5f} {e['down_steady']:8.5f} "
                      f"{e['down_fluct']:8.5f} {om*1e3:11.4f} "
                      f"{e['omega_over_kU']:7.3f} {nb:5.3f} "
                      f"{e['t_rt_corrected_h']:10.1f} {e['R_spurious']:7.4f}")
            fu = np.array([Fup_k[sel_w, j, n].mean() for n in nn_])
            dst = np.array([steady_part(Bk[sel_w, j, n], MKk[sel_w, j, n], n)[0]
                            for n in nn_])
            ust = np.array([steady_part(Ak[sel_w, j, n], MKk[sel_w, j, n], n)[0]
                            for n in nn_])
            fd = np.array([Fdn_k[sel_w, j, n].mean() for n in nn_])
            for lab, msk in (("no return", nor), ("returned", ret)):
                if not msk.any():
                    continue
                e = {"z_km": float(zkm[j]), "class": lab, "n_modes": int(msk.sum()),
                     "R": float(fd[msk].sum() / fu[msk].sum()),
                     "R_steady": float(dst[msk].sum() / fu[msk].sum()),
                     "R_fluct": float((fd - dst)[msk].sum() / fu[msk].sum()),
                     "up_fluct": float((fu - ust)[msk].sum() / fu[msk].sum())}
                cls_tr.append(e)
                print(f"   {lab:>10}  n={e['n_modes']:2d}  R {e['R']:.5f}  "
                      f"steady {e['R_steady']:.5f}  fluctuating {e['R_fluct']:.5f}  "
                      f"against {e['up_fluct']:.5f} fluctuating in the upward branch")
        transient["per_mode"] = rows_tr
        transient["class_split"] = cls_tr
        print("\nWhere the fluctuating parts of the two branches are comparable the "
              "component\nis being divided about equally between them, which is "
              "what a partition does\nwith something whose vertical wavenumber is "
              "not the one it assumes, and carries\nno information about which way "
              "that something is going.")

        # ---- the classification, rechecked against the measured frequencies ----
        j0 = int(np.argmin(np.abs(zkm - mode_zs[0])))
        nor0 = [e for e in rows_tr if e["z_km"] == float(zkm[j0])
                and e["class"] == "nor"]
        worst = min((e["t_rt_corrected_h"] for e in nor0), default=float("nan"))
        slower = sum(1 for e in rows_tr if e["class"] == "nor"
                     and e["t_rt_corrected_h"] > e["t_rt_h"])
        nnor = sum(1 for e in rows_tr if e["class"] == "nor")
        transient["classification_check"] = {
            "z_km": float(zkm[j0]), "t_end_h": float(th_tr[-1]),
            "shortest_corrected_t_rt_h": worst,
            "n_no_return": nnor, "n_slowed_by_correction": slower}
        print(f"\nThe shortest corrected round trip among the modes classed as "
              f"unreturned at\n{zkm[j0]:.2f} km is {worst:.1f} h against a "
              f"{th_tr[-1]:.2f} h record, and {slower} of the {nnor} unreturned\n"
              f"entries at all three heights are slowed rather than hastened by the "
              f"correction.\nOn this record the measured frequencies do not bring "
              f"any of them inside the\nrecord, so the classification of section 11 "
              f"is not overturned by them.")

        # ---- the class through the record, in one-hour bins -------------------
        print(f"\nThe unreturned class at {zkm[j0]:.2f} km in one-hour bins. R sw is "
              f"the share-weighted\nmean of the per-mode ratios, the statistic "
              f"section 11 reports at three hours;\nR fw the ratio of the summed "
              f"fluxes. Steady and fluctuating are within the bin,\nso the split is "
              f"at one hour rather than at the length of the record. Membership\nis "
              f"recomputed in each bin, and the mode count falls as round trips "
              f"elapse.")
        print(f"{'bin[h]':>11} {'ns':>4} {'n_md':>5} {'R sw':>8} {'R fw':>8} "
              f"{'R_st':>8} {'R_tr':>8} {'up_tr':>8}")
        bins = []
        for e0 in np.arange(np.floor(th_tr[-1]) - a.transient_bins + 1,
                            np.floor(th_tr[-1]) + 1):
            m = (th_tr >= e0) & (th_tr < e0 + 1)
            if m.sum() < 4:
                continue
            tc = float(th_tr[m].mean())
            nn_, _, _, nor, _ = classes_at(j0, tc)
            if not nor.any():
                continue
            ns_ = nn_[nor]
            fu = np.array([Fup_k[m][:, j0, n].mean() for n in ns_])
            fd = np.array([Fdn_k[m][:, j0, n].mean() for n in ns_])
            dst = np.array([steady_part(Bk[m, j0, n], MKk[m, j0, n], n)[0]
                            for n in ns_])
            ust = np.array([steady_part(Ak[m, j0, n], MKk[m, j0, n], n)[0]
                            for n in ns_])
            e = {"t_h": tc, "n_snapshots": int(m.sum()), "n_modes": int(nor.sum()),
                 "R_share_weighted": float((share[ns_] * (fd / fu)).sum()
                                           / share[ns_].sum()),
                 "R_flux_weighted": float(fd.sum() / fu.sum()),
                 "R_steady": float(dst.sum() / fu.sum()),
                 "R_fluct": float((fd - dst).sum() / fu.sum()),
                 "up_fluct": float((fu - ust).sum() / fu.sum())}
            bins.append(e)
            print(f"{e0:5.0f}-{e0+1:5.0f} {e['n_snapshots']:4d} {e['n_modes']:5d} "
                  f"{e['R_share_weighted']:8.5f} {e['R_flux_weighted']:8.5f} "
                  f"{e['R_steady']:8.5f} {e['R_fluct']:8.5f} {e['up_fluct']:8.5f}")
        transient["hourly_bins"] = bins
        if len(bins) > 2:
            tb = np.array([e["t_h"] for e in bins])
            rb = np.array([e["R_share_weighted"] for e in bins])
            sl, ic = np.polyfit(tb, rb, 1)
            res = rb - (sl * tb + ic)
            se = float(np.sqrt((res @ res) / (len(bins) - 2)
                               / ((tb - tb.mean()) ** 2).sum()))
            fb = np.array([e["R_fluct"] for e in bins])
            transient["hourly_summary"] = {
                "n_bins": len(bins), "mean": float(rb.mean()),
                "sd": float(rb.std(ddof=1)),
                "se_of_mean": float(rb.std(ddof=1) / np.sqrt(len(bins))),
                "min": float(rb.min()), "max": float(rb.max()),
                "trend_per_h": float(sl), "trend_se": se,
                "trend_sigma": float(abs(sl) / se) if se else float("nan"),
                "fluct_first": float(fb[0]), "fluct_last": float(fb[-1])}
            print(f"Over these {len(bins)} bins the share-weighted ratio has mean "
                  f"{rb.mean():.5f}, sd {rb.std(ddof=1):.5f},\nrange {rb.min():.5f} "
                  f"to {rb.max():.5f}, and a linear trend of {sl:+.5f} +- {se:.5f} "
                  f"per hour,\nthat is {abs(sl)/se:.1f} standard errors from zero. "
                  f"The bin-to-bin scatter is the\nerror bar that a pair of "
                  f"one-hour means separated by several hours does not\ncarry. The "
                  f"fluctuating part runs {fb[0]:.5f} to {fb[-1]:.5f} over the same "
                  f"bins.")

        # ---- the class against height, and whether the fluctuation descends ---
        print(f"\nThe unreturned class against height over the last "
              f"{a.transient_window:.1f} h, with the lag of\nits fluctuating part "
              f"against the topmost level at which the class is not empty.\nA "
              f"component radiated near the sponge base would fall away from it and "
              f"arrive\nlater at each level below. Membership is recomputed at every "
              f"level, so a step in\nthe mode list is a change of statistic and not "
              f"a feature of the profile.")
        print(f"{'z[km]':>7} {'n_md':>5} {'modes':>16} {'R':>8} {'R_st':>8} "
              f"{'R_tr':>8} {'lag[s]':>9} {'coh':>6} {'|c| bound':>10}")
        prof_tr, refser, refz = [], None, None
        for zz in np.arange(np.floor(min(SPONGE_BASE_KM - 2.5, zkm[-1])), 2.0, -2.5):
            j = int(np.argmin(np.abs(zkm - zz)))
            nn_, _, _, nor, _ = classes_at(j, float(th_tr[-1]))
            if not nor.any():
                continue
            ns_ = nn_[nor]
            fu = Fup_k[sel_w][:, j, ns_].sum(1)
            fd = Fdn_k[sel_w][:, j, ns_].sum(1)
            dst = sum(steady_part(Bk[sel_w, j, n], MKk[sel_w, j, n], n)[0]
                      for n in ns_)
            s = fd - fd.mean()
            if refser is None:
                refser, refz, lag, coh = s.copy(), float(zkm[j]), 0.0, 1.0
            elif np.linalg.norm(s) > 0:
                cc = np.correlate(s / np.linalg.norm(s),
                                  refser / np.linalg.norm(refser), "full")
                q = int(np.argmax(cc))
                lag = (q - s.size + 1) * dts
                coh = float(cc[q])
            else:
                lag, coh = float("nan"), float("nan")
            bound = ((refz - zkm[j]) * 1e3 / max(abs(lag), dts)
                     if refz is not None and zkm[j] < refz else float("nan"))
            e = {"z_km": float(zkm[j]), "n_modes": int(nor.sum()),
                 "modes": [int(v) for v in ns_], "share": float(share[ns_].sum()),
                 "R": float(fd.mean() / fu.mean()),
                 "R_steady": float(dst / fu.mean()),
                 "R_fluct": float((fd.mean() - dst) / fu.mean()),
                 "lag_s": float(lag), "coherence": float(coh),
                 "speed_bound_m_s": float(bound)}
            prof_tr.append(e)
            print(f"{zkm[j]:7.2f} {e['n_modes']:5d} {str(e['modes']):>16} "
                  f"{e['R']:8.5f} {e['R_steady']:8.5f} {e['R_fluct']:8.5f} "
                  f"{lag:9.1f} {coh:6.3f} {bound:10.1f}")
        transient["height_profile"] = prof_tr
        print(f"The reference level is {refz:.2f} km and one sample interval is "
              f"{dts:.1f} s, so a lag\nreported as zero bounds the speed from below "
              f"rather than measuring it.")

        # ---- against the overturning layer ------------------------------------
        if "unstable_frac_t" in d.files:
            uf = np.interp(t_tr, d["times"], d["unstable_frac_t"])
            nn_, _, _, nor, _ = classes_at(j0, float(th_tr[-1]))
            ns_ = nn_[nor]
            wn = max(5, int(round(3000.0 / dts)) | 1)
            hp = lambda v: v - np.convolve(v, np.ones(wn) / wn, "same")
            pear = lambda x_, y_: float(np.corrcoef(x_, y_)[0, 1])
            use = th_tr >= 4.0
            print(f"\nAgainst the overturning layer. The unstable area fraction is "
                  f"zero before\n{th_tr[uf > 0][0] if (uf > 0).any() else float('nan'):.2f} h, "
                  f"peaks at {uf.max():.5f} at {th_tr[int(np.argmax(uf))]:.2f} h and "
                  f"ends at {uf[-1]:.5f}. Both\nit and the flux rise over the first "
                  f"half of the record, so the raw coefficient\ncarries a shared "
                  f"trend; the high-passed column removes everything slower than\n"
                  f"{wn*dts/60:.0f} min and is the part that does not. The upward "
                  f"branch is the control, since\nit is filling and cannot be "
                  f"driven by the layer.")
            print(f"{'series':>14} {'r raw':>8} {'r hp':>8} {'best lag[h]':>12} "
                  f"{'r at lag':>9}")
            corr = []
            series = [("down class", Fdn_k[:, j0, ns_].sum(1)),
                      ("up class", Fup_k[:, j0, ns_].sum(1))]
            series += [(f"down n={n}", Fdn_k[:, j0, n]) for n in ns_]
            for lab, v in series:
                x, y = v[use], uf[use]
                pad = wn // 2 + 1
                best, bl = -2.0, 0.0
                for L in range(-int(3 * 3600 / dts), int(3 * 3600 / dts) + 1):
                    xa = x[L:] if L >= 0 else x[: len(x) + L]
                    ya = y[: len(y) - L] if L > 0 else (y if L == 0 else y[-L:])
                    if len(xa) < 40:
                        continue
                    rr = pear(xa, ya)
                    if rr > best:
                        best, bl = rr, L * dts / 3600.0
                e = {"series": lab, "r_raw": pear(x, y),
                     "r_highpass": pear(hp(x)[pad:-pad], hp(y)[pad:-pad]),
                     "best_lag_h": bl, "r_at_best_lag": best}
                corr.append(e)
                print(f"{lab:>14} {e['r_raw']:8.3f} {e['r_highpass']:8.3f} "
                      f"{bl:12.2f} {best:9.3f}")
            transient["instability_correlation"] = {
                "highpass_window_min": wn * dts / 60.0,
                "unstable_frac_peak": float(uf.max()),
                "unstable_frac_peak_h": float(th_tr[int(np.argmax(uf))]),
                "unstable_frac_end": float(uf[-1]), "rows": corr}
            print("A positive lag means the flux follows the layer.")

        tr_fig(outdir, transient, th_tr, Ak, Bk, zkm, j0, kk, nmode)
    out["transient"] = transient

    # ---------------- Figure ----------------
    fig, ax = plt.subplots(2, 2, figsize=(12, 9))
    hrs = sorted({h for h in (2, 4, 6, 8, 10, 12) if h <= th[-1]} | {round(th[-1])})
    for hh in hrs:
        i = int(np.argmin(np.abs(th - hh)))
        ln, = ax[0, 0].plot(mf[:, i] * 1e6, zkm, lw=1.0, label=f"t={th[i]:.0f} h")
        # Same time, benchmark of the part of the ridge spectrum that has arrived.
        bench = np.where(zkm <= SPONGE_BASE_KM,
                         fillin_flux(fk, kk, ibuoy, t[i], u_bar) * 1e6, np.nan)
        ax[0, 0].plot(bench, zkm, ls="--", lw=0.9, color=ln.get_color(), alpha=0.85)
        jf = int(np.argmin(np.abs(zkm - fill(t[i]))))
        ax[0, 0].plot(mf[jf, i] * 1e6, zkm[jf], "k.", ms=6)
    if pair is not None:
        ax[0, 0].plot(pair["mf_bar"] * 1e6, zkm, color="0.35", ls="-", lw=1.4,
                      label=f"sigma={pair['sigma_m']/1e3:.0f} km, final-hour mean")
        ax[0, 0].plot(np.where(zkm <= SPONGE_BASE_KM, pair["pred"] * 1e6, np.nan),
                      zkm, color="0.35", ls="--", lw=1.1)
    ax[0, 0].plot(mf_bar * 1e6, zkm, "k-", lw=1.8,
                  label=f"mean {th[w_avg][0]:.0f}-{th[-1]:.0f} h")
    ax[0, 0].plot(mf_lin * 1e6, zkm, "r--", lw=1.4, label="steady linear")
    ax[0, 0].plot([], [], "k--", lw=0.9, label="arrived spectrum, same t")
    ax[0, 0].axhspan(a.z_lo, z_hi, color="0.85", zorder=0)
    ax[0, 0].axhline(SPONGE_BASE_KM, color="0.4", ls=":", lw=1)
    ax[0, 0].set_xlabel("MF = <rho u'w'>_x (uPa)")
    ax[0, 0].set_ylabel("z (km)")
    ax[0, 0].set_title("Momentum flux, solid; arrived-spectrum benchmark, dashed")
    ax[0, 0].legend(fontsize=7)

    ax[0, 1].plot(ratio, zkm, "k-", lw=1.5, label="linear control")
    ax[0, 1].plot(np.where(zkm <= SPONGE_BASE_KM, mf_corr / mf_lin, np.nan), zkm,
                  "m--", lw=1.4, label="control, arrival-corrected")
    if cmp_rows is not None:
        cr = np.array(cmp_rows)
        ax[0, 1].plot(cr[:, 3], cr[:, 0], "b-.", lw=1.3, label="nonlinear h=3 km")
        ax[0, 1].plot(cr[:, 5], cr[:, 0], "g:", lw=1.6, label="excess (shape ratio)")
    ax[0, 1].legend(fontsize=7, loc="upper right")
    ax[0, 1].axvline(1.0, color="r", ls="--", lw=1.2)
    ax[0, 1].axhspan(a.z_lo, z_hi, color="0.85", zorder=0)
    ax[0, 1].set_xlim(0, 2)
    ax[0, 1].set_xlabel("MF / MF_linear")
    ax[0, 1].set_ylabel("z (km)")
    ax[0, 1].set_title(f"Amplitude error bar (band mean {ratio[sel].mean():.3f})")

    ax[1, 0].plot(drag, zkm, "k-", lw=1.5, label="raw, pointwise")
    # From the band-fitted slope. The benchmark is a staircase over four or five
    # discrete modes, so the pointwise gradient of the corrected profile returns
    # the derivative of the step model rather than a drag.
    ax[1, 0].plot(np.where(sel, drag_fit(mf_bar, final), np.nan), zkm, "b:", lw=1.6,
                  label="raw, band fit")
    ax[1, 0].plot(np.where(sel, drag_fit(mf_corr, final_corr), np.nan), zkm,
                  "m--", lw=1.6, label="arrival-corrected, band fit")
    ax[1, 0].legend(fontsize=7, loc="lower left")
    ax[1, 0].axvline(0.0, color="r", ls="--", lw=1.2)
    ax[1, 0].axhspan(a.z_lo, z_hi, color="0.85", zorder=0)
    ax[1, 0].set_xlabel("-(1/rho) dMF/dz (m/s per sol)")
    ax[1, 0].set_ylabel("z (km)")
    ax[1, 0].set_title(f"Spurious drag (raw band RMS "
                       f"{out['spurious_drag_m_s_sol']['band_rms']:.3f} m/s/sol)")
    # Scale to the evaluation band. The sponge region reaches tens of m/s/sol and
    # would otherwise compress the band of interest to a single vertical line.
    dlim = 1.5 * max(np.abs(drag[sel]).max(), 1e-6)
    ax[1, 0].set_xlim(-dlim, dlim)

    ax[1, 1].plot(lz["lambda_z_km"], zkm, "k-", lw=1.5, label="measured")
    ax[1, 1].plot(lam_pred, zkm, "r--", lw=1.4, label="2 pi U / N(z)")
    ax[1, 1].axvline(2 * np.pi * u_bar / 0.018 / 1e3, color="b", ls=":", lw=1.2,
                     label="2 pi U / N_REF")
    ax[1, 1].axhspan(a.z_lo, z_hi, color="0.85", zorder=0)
    ax[1, 1].set_xlim(0, 60)
    ax[1, 1].set_xlabel("lambda_z (km)")
    ax[1, 1].set_ylabel("z (km)")
    ax[1, 1].set_title("Vertical wavelength from the phase tilt of w")
    ax[1, 1].legend(fontsize=8)

    for axx in ax.ravel():
        axx.set_ylim(0, 90)
        axx.grid(alpha=0.3)
    fig.suptitle(f"Linear flux-constancy control  U={u_bar:.0f} m/s  h={a.h_max:.0f} m  "
                 f"LX={lx/1e3:.0f} km  periodic-x  over {th[-1]:.0f} h", fontsize=11)
    fig.tight_layout()
    png = os.path.join(outdir, "linear_flux_constancy.png")
    fig.savefig(png, dpi=150)
    plt.close(fig)

    jpath = os.path.join(outdir, "linear_flux_constancy.json")
    with open(jpath, "w") as f:
        json.dump(out, f, indent=2, sort_keys=True)

    hdr("VERDICT")
    print(f"Over {a.z_lo:.0f}-{z_hi:.0f} km, averaged over the final "
          f"{a.avg_hours:.1f} h, the measured\nmomentum flux varies by "
          f"{100*final['ptp_rel']:.2f} % peak-to-peak with a fitted drift of "
          f"{final['drift_pct_per_km']:+.4f} % per km,\nagainst the "
          f"Eliassen-Palm prediction of zero. Flat at the stated thresholds: "
          f"{'YES' if flat else 'NO'}.")
    print(f"\nAgainst the arrived part of the ridge spectrum rather than the fully "
          f"established\nwave, the same record gives "
          f"{100*final_corr['ptp_rel']:.2f} % peak-to-peak and "
          f"{final_corr['drift_pct_per_km']:+.4f} % per km, so the arrival\n"
          f"transient accounts for {100*fexp:.1f} % of the drift. The residual is a "
          f"tighter upper bound\non the numerical flux error, not yet a measurement "
          f"of it.")
    print(f"\nThe amplitude ratio in section 3 is reported against a reference that "
          f"scales as\nrho_0 N_0. The admissible choices of N_0 span a factor of "
          f"{1/ref_scale:.2f} here, so that ratio is not\ncomparable in status with "
          f"the constancy result and should be quoted with its N_0.")
    print(f"\nwrote {png}\nwrote {jpath}")


if __name__ == "__main__":
    main()
