# Momentum flux constancy in the linear control corner

*Note dated 2026-10-02.* Every run in this report imposes the ground `w` through the ghost cells (`bottom_forcing = :ghost`, the only path before 2026-09-30). Absolute flux levels and amplitudes quoted here are superseded for amplitude: at CFL 0.4 the ghost path loses 10.7 % of the transmitted flux in uniform wind against the face path, by an amount first order in `dt` (`reflection_convergence_tests.md` Section 14, addendum 2 and addendum 3 step 1). Shapes, timings and ratios between runs at equal `dt` are not covered by this note.

Internal validation report, 2026-09-01. Companion to
`docs/ext6h_instability_mechanism.md` and to `docs/TODO_tracker.md`.

Question addressed: how much of the height dependence of the momentum flux
`MF(z) = <rho u'w'>_x` reported in the Arsia Mons runs is produced by the numerics
rather than by the wave. Every flux figure quoted so far, including the 24 to 29
per cent plateau at 25 km, the 58 to 67 per cent deposition partition over 25 to
40 km and the 0.83 to 0.86 Stage 2 deposition fractions, has been reported without
such a bound.

Reproduce with

```
julia --project=. scripts/sweep_arsia_mons_stage2.jl linear
python3 scripts/analyze_arsia_stage1_deepdive.py \
    --run-dir output/arsia_mons_linear12h_LX2000/run_U40.0_h0.1 \
    --outdir  output/arsia_mons_linear12h_LX2000_analysis/run_U40.0_h0.1 \
    --u-bar 40 --label linear12h
python3 scripts/analyze_linear_flux_constancy.py
```

Metrics are written to `linear_flux_constancy.json` and the figures to
`linear_flux_constancy.png`, `linear_flux_updown.png`, `linear_flux_permode.png`
and `linear_flux_transient.png` in the analysis directory. Section 12 of the
script rereads the record at `--transient-stride`, which is the only part of it
whose cost scales with the length of the run; `--transient-stride 0` skips that
section and leaves everything else unchanged.

The 15 h `sigma = 50` km record of Sections 14 to 16 is analysed with

```
python3 scripts/analyze_linear_flux_constancy.py \
    --run-dir output/arsia_mons_linear15h_LX2000_sig50/run_U40.0_h0.1 \
    --h-max 100 --sigma 50e3 --pair-npz '' --forecast-from 12
```

A companion notebook, `scripts/math_mountain_waves.ipynb` with figures in
`output/math_mountain_waves/`, carries the Taylor-Goldstein derivations behind
Sections 15 to 18 and reproduces Nappo Figures 3.11 and 3.12. Two of its results
are used above and were checked independently before being quoted, the
`N(0)/N_band` identity of Section 19 item 6 and the partial-reflection integral
of Section 18; its other contents have not been verified here.

The `h = 30` m amplitude control of Section 18 was produced and analysed with

```
julia --project=. scripts/sweep_arsia_mons_stage2.jl linear sigma=50 hours=15 h=0.03
python3 scripts/analyze_arsia_stage1_deepdive.py \
    --run-dir output/arsia_mons_linear15h_LX2000_sig50/run_U40.0_h0.03 \
    --outdir  output/arsia_mons_linear15h_LX2000_sig50_analysis/run_U40.0_h0.03 \
    --u-bar 40 --label linear15h_sig50_h30m
python3 scripts/analyze_linear_flux_constancy.py \
    --run-dir output/arsia_mons_linear15h_LX2000_sig50/run_U40.0_h0.03 \
    --h-max 30 --sigma 50e3 --pair-npz '' --forecast-from 12
```

and the `sigma = 100` km twin of Section 16 with the same command on
`output/arsia_mons_linear12h_LX2000/run_U40.0_h0.1` at the default `--sigma`.

## 1. Why a linear corner is needed

The production configuration departs from the assumptions of the analytic results
it is compared against in several ways at once, so an agreement or a disagreement
in that configuration is difficult to attribute.

At the launch level `Nh/U = 1.35` for `h = 3` km using the near-surface
`N = 0.018` s^-1, so the wave is at or past overturning where it is forced. The
`:outflow` x boundary of Stage 1 admits lateral influx after the advective transit
time `LX/(2U) = 6.94` h, which is where the x-mean momentum budget was observed to
stop closing. The flux-carrying horizontal scale falls to about 11 points at the
worst moment of the 12 h record. The sponge base at 60 km sits 3 to 7 km above the
deceleration maximum. Under those conditions the observation that the Lindzen
saturation level reproduces the layer base to within 0.5 to 2.0 km cannot be
separated from a coincidence among several biased quantities, because the solver's
own error on the flux has never been measured.

Reducing the ridge height to `h = 0.1` km removes the nonlinearity, and adopting
the Stage 2 periodic x boundary removes the lateral term. At that height the
steady linear solution for this ridge gives a maximum wave-induced `|u'|` of
18.1 m/s at 54.8 km against `U = 40` m/s, so `u_hat/U = 0.45` and the stationary
wave overturning criterion `u_hat >= U` is not met anywhere below the sponge. The
wave is then conservative, and the Eliassen-Palm result requires

    dMF/dz = 0    exactly.

Any measured departure is therefore an error bar rather than a physical signal.
Integrating `2 nu m^2 / cg_z` through the MCD `nu` column gives a cumulative
molecular-viscosity flux loss of 0.03 per cent by 55 km, rising to 4.6 per cent
only at 90 km inside the sponge, so physical dissipation cannot account for a
drift measured below 60 km.

## 2. Why neither existing run serves as the control

Both sweeps already contain an `h = 0.1` km case at `U = 40` m/s, and neither is
usable, for reasons that are visible in their final snapshots.

The Stage 1 case is uniform wind with the outflow x boundary. Its final snapshot
carries `max|w| = 20.8` m/s and `max|u - 40| = 36.6` m/s, against a linear
expectation of 0.65 and 18.1 m/s at 55 km and considerably less below. Its
measured `MF` reaches -80.4 uPa at 16.8 km, or 3.3 times the -24.70 uPa steady
linear value for that domain width, and the profile does not resemble the flat
one theory predicts. The wave field in that configuration is not measuring the
mountain wave.

The Stage 2 case is periodic but uses the real reversing `u_zonal(z)`. Its
amplitudes are clean, `max|u'| = 3.5` m/s and `max|w| = 1.08` m/s, but `MF`
changes sign between 4.8 and 9.8 km at the 4.5 km critical level and settles near
0.2 uPa, about 1 per cent of the linear value, because that critical level absorbs
the wave before it reaches the altitudes of interest.

The combination not previously run is uniform wind together with periodic x, which
is what the `linear` mode of `scripts/sweep_arsia_mons_stage2.jl` provides.

## 3. Configuration

The control differs from the Stage 1 12 h extension in exactly two respects, the
ridge height and the horizontal boundary condition. Everything else is held fixed.

| setting | value | same as |
|---|---|---|
| ridge height `h_max` | 100 m | changed from 3000 m |
| horizontal BC | `:periodic` | changed from `:outflow` |
| background wind `u_init` | uniform 40 m/s | Stage 1 |
| forcing wind `U_bar` | 40 m/s | Stage 1 and Stage 2 |
| domain width `LX` | 2000 km, `nx = 500` | Stage 1 |
| `dx`, `dz` | 4000 m, 500 m | both stages |
| domain height `LZ` | 90 km, `nz = 180` | both stages |
| sponge | 60 to 90 km, `K0 = 0.05` | both stages |
| ridge half-width `sigma` | 100 km | both stages |
| background atmosphere | MCD MY31 Arsia column | both stages |
| viscosity, cooling | on, off | both stages |
| `CFL` | 0.4 | both stages |
| horizon | 12 h | Stage 1 extension |
| snapshot cadence | every 100 steps | both stages |

`LX = 2000` km rather than the Stage 2 default of 4000 km, so the resulting error
bar transfers directly to the Stage 1 flux quotes and applies at the same grid.
It is also the more demanding of the two, since 2000 km spans 2.87 flux-weighted
horizontal wavelengths against 5.7 at 4000 km, so periodic images are closer.

Note that `MF` as defined is a domain x-mean and therefore scales as `1/LX`. The
same wave gives -24.703 uPa at 2000 km and -12.509 uPa at 4000 km. Comparison of
the absolute flux against a run at a different width requires that normalization;
the constancy test itself is unaffected, since it is a statement about the height
dependence at one width.

## 4. Horizon, and why three hours would not have sufficed

The acceptance criterion concerns the flux over roughly 10 to 55 km, and the
stationary wave reaches those heights only slowly. Its hydrostatic vertical group
velocity is `cg_z = k_eff U^2 / N(z)`, which with the flux-weighted
`k_eff = 9.012e-06` 1/m and `N(z)` from the MCD column is 1.5 to 1.8 m/s over the
region of interest. The transit time is `int dz / cg_z = int N dz' / (k U^2)`,
which the script now evaluates by integrating `N` rather than by averaging `cg_z`
between levels; the two differ because `N` falls by a factor of about 2.1 across
the near-surface inversion, and the resulting arrival times shift by under 2 per
cent. On the exact integral the arrival is 19.75 km at 3.55 h, 29.75 km at 5.06 h,
49.75 km at 8.40 h and 54.75 km at 9.33 h.

A 3 h run therefore fills only to about 17 km, and any flux measured above that
height would be the impulsive startup transient of the lower boundary condition
(`src/forcing.jl:115`, no temporal envelope) rather than the steady wave. The 12 h
horizon leaves about 2.6 h of settling after the front passes 55 km. The analysis
script computes this arrival curve from the run's own background column and
restricts the evaluation band to filled heights, so the band is derived rather
than assumed.

## 5. A correction to the vertical wavelength target

The vertical wavelength of a stationary hydrostatic wave is
`lambda_z = 2 pi U / N`, which is often quoted for this configuration as 13.96 km.
That value follows from `N_REF = 0.018` s^-1, the constant used for the `Nh/U`
label in both sweep drivers and described there as a lower-atmosphere mean.
Computing `N(z)` from the MCD column shows that it is instead the near-surface
value, produced by the inversion in which `T` rises from 187.3 K at the ground to
213.0 K at 5 km.

Over 10 to 55 km, `N` is 0.0078 to 0.0097 s^-1, so `lambda_z` is 26.3 to 32.1 km,
or 53 to 64 points at `dz = 500` m. The 13.96 km figure applies only below about
2 km, where the linear solution gives 16.4 km. The measured phase tilt is compared
against `2 pi U / N(z)` rather than against a single number.

The same height dependence bears on the `Nh/U` label itself. At `h = 3` km the
label reads 1.35 using the surface `N`, which is the appropriate value at the
launch level, but the corresponding ratio using the `N` of the propagation region
is about 0.6. This does not change the conclusion that the `h = 3` km run
overturns, which is established directly from `min_x u <= 0`, but it does mean the
sweep label should not be read as a nonlinearity measure aloft.

## 6. The run

| quantity | value |
|---|---|
| wall-clock | 2500.7 s |
| steps | 58564 |
| `final_t` | 43200.0 s, `converged = true` |
| timestep | 0.7370 to 0.7400 s over 117 diagnostic samples, a spread of 0.41 per cent |
| snapshots | 586 at `OUTPUT_INTERVAL = 100`, 2.0 GB |
| analysis chain | under 72 s for extraction, figures and metrics |

The cost came in below the 4500 to 6000 s estimated, because that estimate assumed
the timestep would drift downward as it did in the nonlinear run. It did not.

The constancy of the timestep is a result in its own right. In the 12 h `h = 3` km
extension the timestep fell from about 0.55 to 0.43 s over the second half, which
was attributed to locally supersonic samples appearing from 9.54 h onward, since
the CFL limit scales with `|u| + c_s`. That attribution rested on a single run. The
control holds the grid, horizon, background column, sponge, viscosity settings and
`CFL` fixed and varies only the ridge height and the horizontal boundary, and its
timestep is flat to 0.41 per cent. Under the assumptions of the model, nothing
about integrating this configuration for 12 h reduces the timestep, which supports
reading the reduction in the nonlinear run as a property of the flow. The
`sigma = 50` km twin of Section 13 repeats that behaviour, holding 0.737 to
0.740 s over its own 58609 steps, so the observation is not specific to one ridge
spectrum.

The analysis chain also ran about five times faster than the 350 s scaled from the
earlier extraction, plausibly because the snapshots were still in the page cache
from the run that had just written them.

## 7. The simulated wave is the intended wave

Three independent checks indicate the control is measuring the stationary mountain
wave rather than something else.

The measured vertical wavelength, read from the local slope of the unwrapped phase
of `w` at the flux-dominant horizontal wavenumber, agrees with `2 pi U / N(z)` to
a band mean ratio of 0.9958 over 10 to 55 km, with individual levels between 0.89
and 1.17. The flux-dominant horizontal wavenumber is 9.425e-06 1/m, or
`lambda_x = 667` km, against the flux-weighted `k_eff` of the ridge spectrum at
9.012e-06 1/m or 697 km. The one-sided cospectrum used to select that wavenumber
sums back to the configuration-space `mean_x(u'w')` to 4.4e-16 of its maximum.

The run stayed linear. The minimum `min_x u` anywhere below 60 km over the whole
record is +29.735 m/s, or 0.743 `U`, so the overturning criterion `min_x u <= 0`
is never approached. The largest departure of the zonal-mean wind from 40 m/s at
55 km is 0.227 m/s, so the wave does not modify the flow it propagates through.

The wave amplitude at the end of the record is close to the linear prediction. The
maximum hill-window `|w|` at 11.99 h is 0.85, 0.83, 0.82 and 0.69 of the steady
linear value at 10, 25, 40 and 55 km. The largest `|w|` in the whole record, 8.36
m/s below the sponge, occurs between 0.02 and 0.10 h and is the impulsive
switch-on transient of the lower boundary condition, not the steady wave; taking a
record maximum without reference to time would misrepresent the amplitude by more
than an order of magnitude.

## 8. The flux is not height-independent against a fully established wave

Averaged over the final hour, `MF` over 10 to 55 km varies by 50.58 per cent
peak-to-peak about its median of -16.377 uPa, with a least-squares fitted drift of
-1.2426 per cent per km, accumulating to -55.9 per cent across the 45 km band. The
Eliassen-Palm prediction is zero. At the reporting thresholds of 5 per cent
peak-to-peak and 0.1 per cent per km, the flux is not flat.

That comparison uses the wrong benchmark, and Section 10 replaces it. The steady
linear value assumes every wavenumber in the ridge spectrum has arrived, which at
55 km after 12 h is not the case for the two longest modes. The figures in this
section and in Section 9 are retained because they are what the acceptance
criterion asked for, and because they bound the sum of the transient and the
numerics; they are not the numerical error on their own.

The figure should also not be read as a fixed error, because every measure of it
is still improving at the end of the record.

| t [h] | 1.0 | 3.0 | 5.0 | 7.0 | 9.0 | 11.0 | 12.0 |
|---|---|---|---|---|---|---|---|
| median MF [uPa] | -0.180 | -1.401 | -4.557 | -9.147 | -13.091 | -15.954 | -16.571 |
| peak-to-peak / median | 6.822 | 5.126 | 3.046 | 1.668 | 1.021 | 0.570 | 0.476 |
| drift [% per km] | -7.09 | -7.05 | -5.72 | -3.92 | -2.62 | -1.40 | -1.12 |

The spin-up is strongly height-dependent. Fitting `d ln|MF| / dt` over the final
3 h gives -1.5 per cent per hour at 4.75 km and +1.2 at 9.75 km, which is
consistent with those levels having settled, against +5.4 at 29.75 km, +11.1 at
39.75 km and +20.0 at 54.75 km, which are not settled.

This qualifies the group-velocity estimate used to choose the horizon. Integrating
`cg_z = k_eff U^2 / N(z)` gives the arrival of the wave front at 54.75 km at
9.33 h, and the front does arrive: the run has filled to 68.75 km by 11.99 h on
that measure. Arrival and settling are different quantities, and the record indicates
that settling of the flux takes considerably longer than arrival of the front. A
12 h horizon is sufficient to establish the flux below roughly 10 to 20 km and is
not sufficient above about 30 km.

## 9. Bounds the record does support

Two statements survive the lack of convergence.

The measured drift is an upper bound on the numerical drift. Since every measure
of non-constancy is decreasing monotonically in time, the -1.2426 per cent per km
at 11.99 h is larger than whatever the converged value would be. Expressed as
spurious drag, `-(1/rho) dMF/dz` over the final hour is -0.0095 m/s per sol at
10 km, -0.0119 at 25 km, -0.0393 at 30 km, -0.128 at 40 km, -0.702 at 50 km and
-0.868 at 55 km, with a band RMS of 0.321 m/s per sol. The bound is therefore
tight at the altitudes that have settled and loose at those that have not, which
is the same ordering as the spin-up.

At the settled levels the measured flux falls short of the steady linear value.
Averaged over the final hour, `MF/MF_linear` is 0.7586 at 4.75 km and 0.7653 at
9.75 km, where the growth rate indicates the flux has stopped rising; in the final
snapshot alone the same ratios are 0.7844 and 0.7704. The corresponding `|w|`
ratio at 9.75 km is 0.853, and 0.853 squared is 0.728, so the flux shortfall is
consistent with an amplitude shortfall of about 15 per cent rather than with a
separate effect.

That shortfall is not currently attributable to the solver, because the reference
is not determined to better than the quantity being measured. The steady linear
flux is `-rho_0 N_0 U sum_k w_k k |h_k|^2`, so it scales as `rho_0 N_0`, and the
value used takes `N_0` at the first cell centre, `z = 0.25` km, where
`N_0 = 0.01763` s^-1 inside the near-surface inversion. The band mean over 10 to
55 km is `N = 0.00849` s^-1, smaller by a factor of 2.08, and a reference built on
it would put the band mean ratio at 1.32 rather than 0.637. The launch formula is
the WKB result, and WKB does not hold across the inversion here: `N` falls by that
factor over about 4 km while the local vertical wavelength is 14 to 31 km, so part
of the incident flux is reflected rather than transmitted and the transmitted
reference lies below the launch value by an amount this calculation does not
resolve. Net Eliassen-Palm flux is conserved through a non-dissipating partially
reflecting layer, so a shortfall against a single-`N` launch formula is not by
itself a departure from `dMF/dz = 0`.

Two smaller normalization terms were identified earlier and are not the leading
one: the analytic comparison places the surface displacement at the first cell
centre rather than at the boundary, and omits the compressible correction to the
vertical wavenumber, which is about 3.8 per cent of `m^2` for `H = 10.7` km.

The constancy result of Sections 8 and 10 does not inherit this ambiguity. For
hydrostatic waves `m = N/U` does not depend on `k`, and `k^2/m^2` is about
`2 x 10^-3` for this ridge and column, so whatever fraction of the flux the
inversion transmits multiplies every wavenumber by the same factor and cancels
from any profile normalized by its own value at one height. The amplitude ratio
should be quoted with its `N_0`; the shape results need not be.

## 10. Most of the height dependence is an incomplete spectral fill-in

The benchmark used in Sections 8 and 9 is the flux of a fully established wave.
The record does not contain one, and the shortfall is not uniform in height,
because the ridge launches a spectrum whose components arrive at different times.

For a stationary hydrostatic component of horizontal wavenumber `k` the vertical
group velocity is `cg_z = k U^2 / N(z)`, linear in `k`, so the transit time to a
height is `t_arr(k, z) = I(z) / (k U^2)` with `I(z) = int_0^z N dz'`. Arrival time
therefore scales as `1/k`, and the ridge spectrum has no low-wavenumber cutoff, so
what is being watched at late times is the filling of a spectral tail rather than
a relaxation towards a steady state. Extrapolating the drift column of Section 8
as if it were a relaxation is not justified.

The per-wavenumber flux follows from the same WKB solution that gives the steady
value. Substituting it into the one-sided cospectrum leaves

```
F_k = -rho_0 N_0 U w_k k |h_k|^2
```

with `w_k` the conjugate-pair weight, 2 for interior bins and 1 for the mean and
Nyquist bins. Each `F_k` is height-independent on its own, since the height
dependence cancels between `rho`, `N` and the `(rho N)^-1/2` amplitude growth, and
`sum_k F_k` reproduces the steady value to rounding, measured at `0.0e+00` in this
run. The benchmark is then the flux of the components that have arrived,

```
MF_fill(z, t) = sum over k >= I(z) / (U^2 t) of F_k.
```

In a periodic domain of 2000 km only `k = n 2 pi / LX` exist, and the flux is
carried by a handful of them.

| n | lambda_x [km] | flux share | cumulative | t_arr at 55 km [h] |
|---|---|---|---|---|
| 1 | 2000.0 | 0.182 | 0.182 | 26.8 |
| 2 | 1000.0 | 0.271 | 0.452 | 13.4 |
| 3 | 666.7 | 0.248 | 0.700 | 8.9 |
| 4 | 500.0 | 0.166 | 0.866 | 6.7 |
| 5 | 400.0 | 0.085 | 0.951 | 5.4 |

Five bins carry 95 per cent of the flux. At 55 km after 12 h the two longest
modes, together 45.2 per cent of the flux, have not arrived, and they reach that
height only at 13.4 h and 26.8 h. The arrived fraction over the evaluation band in
the final hour is 1.000 up to about 20 km, 0.818 from 25 to 45 km and 0.548 from
50 to 55 km.

Dividing the measured flux by that fraction, snapshot by snapshot before
averaging, removes most of the height dependence.

| quantity over 10 to 55 km, final hour | against the steady value | against the arrived spectrum |
|---|---|---|
| fitted drift [% per km] | -1.2426 | +0.0372 |
| accumulated across the band [%] | -55.9 | +1.7 |
| peak-to-peak / median | 0.506 | 0.279 |

The drift falls by 97.0 per cent. Expressed in the unit the deposition claims are
made in, and taken from the band-fitted slope rather than the pointwise gradient,
the apparent drag drops from -0.008, -0.029, -0.114 and -0.466 m/s per sol at
9.75, 24.75, 39.75 and 54.75 km to +0.0002, +0.001, +0.004 and +0.026 over the
same levels.

The model has one assumption that can be tested against the record, namely that
each component switches on sharply at its own front. Multiplying every arrival
time by a factor `f` and refitting gives -0.877 per cent per km at `f = 0.5`,
+0.037 at `f = 1` and +3.303 at `f = 2`, the last with part of the band below the
arrived-fraction floor and so not interpretable. Solving for the factor that nulls
the drift gives `f = 0.989`. Since `f` is not a free parameter, a null within
about 1 per cent of unity indicates that the WKB group velocity places the fronts
where the record puts them, without tuning.

Two reservations apply to the residual. It remains an upper bound rather than a
measurement, because it still contains whatever the sponge, the dispersive spread
of each component about its own front and the WKB group velocity contribute, none
of which this record separates. And the residual peak-to-peak of 0.279 is not a
solver figure: the measured profile varies smoothly with height while the
benchmark steps at the arrival heights of four or five discrete modes, so the
scatter left after dividing is dominated by that shape mismatch. Within this
section the fitted slope carries information and the scatter does not.

Taken together, the results are consistent with the height dependence of `MF` in
the 12 h control being an arrival transient rather than a numerical flux
divergence, and they do not bound the numerical drift below roughly 1 per cent per
km, which is the spread of the correction under the front-sharpness assumption.

## 11. Differencing against the nonlinear run

The comparison least affected by incomplete spin-up is between the control and the
nonlinear 12 h extension, because both are sampled at the same time in the same
solver on the same grid with the same background column. Normalizing each profile
by its own steady linear value removes the `h^2` scaling, and normalizing again by
its own value at 10 km removes the overall amplitude deficit, leaving the height
dependence.

| z [km] | control ratio | control shape | nonlinear ratio | nonlinear shape | excess |
|---|---|---|---|---|---|
| 9.75 | 0.7704 | 1.0000 | 0.2825 | 1.0000 | 1.000 |
| 14.75 | 0.7571 | 0.9828 | 0.1971 | 0.6976 | 0.710 |
| 19.75 | 0.7299 | 0.9474 | 0.3668 | 1.2984 | 1.370 |
| 24.75 | 0.7166 | 0.9302 | 0.2866 | 1.0144 | 1.091 |
| 29.75 | 0.6900 | 0.8957 | 0.2676 | 0.9473 | 1.058 |
| 34.75 | 0.6540 | 0.8489 | 0.1694 | 0.5997 | 0.706 |
| 39.75 | 0.6180 | 0.8023 | 0.1388 | 0.4914 | 0.613 |
| 44.75 | 0.5783 | 0.7506 | 0.1106 | 0.3916 | 0.522 |
| 49.75 | 0.5084 | 0.6599 | 0.0385 | 0.1361 | 0.206 |
| 54.75 | 0.4509 | 0.5853 | 0.0163 | 0.0578 | 0.099 |

Below about 30 km the excess scatters around unity, between 0.71 and 1.37 with no
trend, so the height dependence of the nonlinear flux in that range is reproduced
by a conservative wave in the same solver over the same interval. Above about
35 km the excess falls monotonically, to 0.61 at 39.75 km, 0.52 at 44.75 km, 0.21
at 49.75 km and 0.099 at 54.75 km, so at those heights the nonlinear run loses
flux by roughly a further factor of ten by 55 km beyond anything the control
produces.

The results are consistent with saturation acting above about 35 km in the
nonlinear run, and do not support attributing the flux decrease below about 30 km
to saturation on the basis of this record. The 58 to 67 per cent deposition
partition previously reported over 25 to 40 km spans the boundary between those
two regimes and should be read with that division in mind.

This table is unchanged by the fill-in correction of Section 10, and for a reason
worth stating. The arrived fraction depends on the launched ridge spectrum, the
background column, the wind and the elapsed time, none of which differ between the
two runs; only the ridge height does, and that cancels in the `h^2` normalization.
Applying the same `MF_fill(z, t)` to both profiles therefore multiplies numerator
and denominator of the excess by the same factor at every level, and the excess is
identical. The one asymmetry runs in the direction that strengthens the reading
rather than weakening it: saturation moves the nonlinear run's flux-carrying scale
to higher `k`, and higher `k` arrives sooner, so if the two spectra differ at all
the nonlinear run is the more completely filled of the two, and the shortfall
attributed to it above 35 km is if anything understated.

## 12. What the 12 h control supports

The measurement asked whether `MF(z)` is constant where theory requires it to be.
It is not, by 50.6 per cent peak-to-peak over 10 to 55 km with a fitted drift of
-1.24 per cent per km, and the acceptance criterion as stated is not met.

The record does not support attributing that departure to the solver. Every
measure of it is decreasing in time, the settled levels below about 10 km show a
much smaller departure than the unsettled levels above 40 km, and the spurious
drag ordering follows the same pattern. The results are more consistent with the
stationary wave still filling the column at 12 h than with a fixed numerical
error, which means the flux error bar the run was intended to produce is bounded
from above but not yet measured.

Section 10 makes that quantitative without a further run. Benchmarking the
measured flux against the part of the ridge spectrum that has arrived, rather than
against a fully established wave, reduces the fitted drift from -1.2426 to +0.0372
per cent per km, and the front factor that would null it is 0.989 against the
parameter-free value of 1. The height dependence of `MF` in this record is
therefore consistent with an arrival transient. What that does not do is bound the
numerical drift tightly: the correction spans roughly 1 per cent per km as the
front-sharpness assumption is varied over its plausible range, so the numerical
component is not separated below about that level from this run alone.

Section 13 tightens that with a second run. Narrowing the ridge to `sigma = 50`
km, which reweights the flux across horizontal wavenumbers while leaving
`lambda_z`, the sponge thickness in wavelengths, the vertical resolution and the
linearity margin unchanged, reduces the measured drift from -1.2426 to -0.1138
per cent per km and the band-fit spurious drag at 54.75 km from -0.466 to -0.080
m/s per sol. A fixed numerical flux divergence would have been reproduced by that
change and was not. Differencing the pair against the fill-in prediction places
the numerical component between 0.02 and 0.26 per cent per km, the spread
reflecting how the front sharpness is treated rather than anything measured.

Two ancillary results are firmer. The wave the model produces has the vertical
wavelength linear theory predicts, to a band mean of 0.996, which supports the
vertical resolution being adequate at 53 to 64 points per wavelength. And the
differential against the nonlinear run localizes genuine saturation to above about
35 km, which is a statement that does not require either run to have converged.

The broader implication is that the flux profiles reported from the 12 h nonlinear
run are also unconverged, since a conservative wave in the same configuration has
not settled above 30 km at that horizon. Flux figures quoted from that record at
those heights describe a filling column rather than a steady state.

## 13. Separating the numerical and transient drifts with a ridge-width pair

Section 10 leaves the measured drift as a sum that a single run cannot resolve.
Writing it as `D(z) + T(z; sigma, t)`, with `D` the numerical part and `T` the
arrival transient, the fill-in model predicts `T` with no numerics in it, but
only `D + T` is observed. A second run at a different ridge half-width separates
them, because `sigma` moves `T` and, under the hypothesis being tested, leaves
`D` alone.

### Why the ridge width moves the transient without moving anything else

An earlier version of this document said that halving `sigma` halves every
arrival time. That is not right and the correction matters for how the pair is
read. Arrival time is `t_arr(k, z) = I(z) / (k U^2)` with `I(z) = int_0^z N dz'`,
which contains `k` and not `sigma`. Narrowing the ridge does not make any given
mode arrive sooner. What it does is move flux out of the slow modes and into the
fast ones, so the flux-weighted mean of an unchanged set of arrival times falls.

Three quantities that would otherwise confound the comparison do not move with
`sigma`, and this is what makes the pair a controlled difference rather than two
separate runs.

1. The vertical structure. For hydrostatic waves `m = N/U` carries no `k`, so
   `lambda_z`, the sponge thickness measured in vertical wavelengths, and the
   number of `dz` per wavelength are identical between the two runs.
2. The linearity margin. The steady wave amplitude is
   `u' = N h (rho_0 N_0 / rho N)^(1/2)`, which also carries no `k`. Computed over
   the evaluation band it is 18.6 m/s in both cases, against `U = 40` m/s, so
   both runs sit at `u_hat/U` of 0.46 to 0.47 and neither approaches overturning.
3. The steady flux itself. For a Gaussian ridge at fixed `h` the discrete sum
   `sum_k k |h_k|^2` is close to independent of `sigma`, and the steady linear
   value moves only from -24.703 to -25.017 uPa, or 1.3 per cent.

### The predicted numbers

Computed from `per_k_steady_flux()` and `fillin_flux()` as they stand, on the
same column, wind, grid and snapshot clock as the existing run, so the only input
that differs is `sigma`.

| quantity | sigma = 100 km | sigma = 50 km |
|---|---|---|
| flux-weighted `k_eff` [1/m] | 9.012e-06 | 1.780e-05 |
| flux-weighted `lambda_x` [km] | 697 | 353 |
| points per `lambda_x` at `dx = 4` km | 174 | 88 |
| peak flux bin | n = 2, 1000 km | n = 5, 400 km |
| flux in the two slowest modes | 0.4524 | 0.1381 |
| 95 per cent of the flux by | n = 5, 400 km, 100 points | n = 11, 182 km, 45 points |
| `k_eff^2/m^2` over the band | 1.4e-03 to 2.1e-03 | 5.4e-03 to 8.3e-03 |
| arrived fraction at 55 km, final hour | 0.5476 | 0.8619 |
| steady linear `MF` [uPa] | -24.703 | -25.017 |
| max `\|u'\|` in band [m/s] | 18.6 | 18.6 |
| **predicted transient drift `T`** | **-1.2637 %/km** | **-0.3137 %/km** |

The predicted separation is `T(100) - T(50) = -0.9500` per cent per km. Against
the measured -1.2426 per cent per km at `sigma = 100` km, the implied numerical
residual is `D = +0.0211` per cent per km, and the drift the `sigma = 50` km twin
should return, if `D` is common to both, is **-0.2926 per cent per km**. This
sits inside the -0.3 to -0.4 per cent per km estimated in advance.

`D` from differencing the fitted slopes here is +0.0211 per cent per km, against
the +0.0372 that Section 10 obtains by dividing each snapshot by its own arrived
fraction and then fitting. The two are different ways of combining the same
model and the same record, and the gap of 0.016 per cent per km is a measure of
how little either procedure resolves at this level.

### Whether the separation is measurable

The fitted drift is well enough determined for a 0.95 per cent per km separation
to be read directly. Over the 49 snapshots of the final hour the per-snapshot
fitted drift has a standard deviation of 0.0854 per cent per km and a range of
-1.4021 to -1.1172. The least-squares slope on the hour-mean profile carries a
formal standard error of 0.0421 per cent per km, with a residual scatter of
5.13 per cent in `ln MF` about the fit. The predicted separation is about
11 times the snapshot-to-snapshot standard deviation, so the two runs are
expected to be distinguishable well outside the noise of the estimate.

### What each outcome would mean

If the `sigma = 50` km run returns a drift near -0.29 per cent per km, then `D`
is common to both and the numerical part of the flux drift is bounded at the
few-hundredths of a per cent per km level, which is the number the control was
built to produce. If it returns something else, the difference is `k`-dependent,
and that is itself the answer to the question the control was meant to settle,
since the nonlinear run's flux-carrying scale collapses from about 980 to 45 km
under saturation. Either result is informative, which is the property the single
longer run does not have.

Two limits on the prediction should be carried with it. First, the `sigma = 50`
km spectrum reaches `k^2/m^2` of 3.1e-02 at the mode where the cumulative flux
passes 95 per cent, against 6.4e-03 for `sigma = 100` km. The wave remains
hydrostatic in both, but the non-hydrostatic correction to `cg_z` is of order
`k^2/m^2`, so predicted arrival times for the fastest flux-carrying modes may be
off by up to about 3 per cent. That is small against a 0.95 per cent per km
separation and is in the same direction for the very modes the narrower ridge
adds. Second, 13.8 per cent of the flux is still unarrived at 55 km at 12 h in
the narrower case, so the `sigma = 50` km run is a better differencing partner
rather than a converged record, and it is not offered as one.

### The twin run

Run on 2026-09-02 with `sigma=50` added to the linear mode of
`scripts/sweep_arsia_mons_stage2.jl`, into
`output/arsia_mons_linear12h_LX2000_sig50/`.

| | sigma = 100 km | sigma = 50 km |
|---|---|---|
| wall-clock [s] | 2500.69 | 2587.66 |
| steps | 58564 | 58609 |
| timestep [s] | steady near 0.72 | steady 0.737 to 0.740 |
| `final_t` [s] | 43200.0 | 43200.0 |
| snapshots | 586 | 587 |

Neither run showed a timestep reduction, which is consistent with the reduction
seen in the nonlinear 12 h run having a nonlinear origin rather than a numerical
one. Extraction of the twin took 28.6 s and the constancy analysis 2.2 s.

The twin stayed in the regime it was meant to test. The minimum `min_x u` below
60 km over the whole record is +21.721 m/s, or 0.543 `U`, so the overturning
criterion is not met; the maximum hill-window `\|w\|` below the sponge is 8.529 m/s;
and the mean flow at 55 km departs from `U` by 0.310 m/s. The measured vertical
wavelength gives a band mean `measured/predicted` of 1.0200 against 0.996 for the
wider ridge, which is the direct check that `lambda_z` did not move with `sigma`.

### An independent check in the time domain

The fill-in argument can be tested at a fixed station without computing a flux at
all, and the twin provides the cleanest version of it. A stationary wave has zero
phase speed, so at a fixed point it produces a constant `w` rather than an
oscillation, and the wavelet of the station record should be empty at every
resolvable frequency once the startup transient has decayed. What should remain,
if the column is still filling, is a slow monotone drift as the long modes arrive,
at an apparent frequency set by the record length rather than by any wave.

Measured at S3, `x_rel` = 298 km, z = 49.75 km, over `t > 6` h, with the cone of
influence excluded.

| band | strongest reliable power | implied amplitude |
|---|---|---|
| 1 to 6.8 mHz, where the transient sat | -48.5 dB | 0.020 m/s |
| 0.3 to 6.8 mHz | -37.9 dB | 0.068 m/s |
| below 0.3 mHz | -22.4 dB at 0.065 mHz | 0.407 m/s |

Against a startup transient that peaks at 5.377 m/s, the wave band is empty to
better than one part in a hundred in amplitude, and the only surviving power sits
at a period of 4.3 h in a 12 h record. The time series says the same thing without
the transform: over `t > 6` h it is 97.7 per cent a linear ramp at -0.0664 m/s per
hour, with a residual standard deviation of 0.0179 m/s and a lag-1 autocorrelation
of 0.998, so what is left after the trend is smooth rather than noisy. The hourly
means run -0.008, +0.043, -0.016 and -0.195 m/s over the 2 to 3, 5 to 6, 8 to 9 and
11 to 12 h windows.

Both features are what the fill-in picture requires. At that height 13.8 per cent
of the flux has not arrived by 12 h, with `n = 1` due at 26.8 h and `n = 2` at
13.4 h, and an arriving mode changes the steady value of `w` at a station rather
than making it oscillate. This is an independent line of evidence for the same
conclusion Section 10 reaches from `MF(z)`, and it uses neither the flux nor the
per-wavenumber benchmark.

### What the twin measured

| over 10-55 km, final hour | sigma = 100 km | sigma = 50 km |
|---|---|---|
| median `MF` [uPa] | -16.376 | -21.194 |
| peak-to-peak / median | 0.5058 | 0.0866 |
| fitted drift [%/km] | -1.2426 | **-0.1138** |
| accumulated across the band [%] | -55.9 | -5.12 |
| per-snapshot sd of the drift [%/km] | 0.0854 | 0.0105 |
| band-fit spurious drag at 54.75 km [m/s/sol] | -0.4662 | -0.0804 |
| pointwise spurious drag, band RMS [m/s/sol] | 0.3207 | 0.1427 |
| `MF/MF_linear`, band mean | 0.6371 | 0.8417 |

The measured drift fell by a factor of 10.9, and the drag in the unit the
deposition claims are made in fell by a factor of 5.8 at the top of the band.
This part of the result depends on no model. With `lambda_z`, the sponge
thickness in vertical wavelengths, the vertical resolution of the wave and the
linearity margin all held fixed, and only the flux weighting across horizontal
wavenumbers changed, a drift that was a fixed numerical flux divergence would
have been reproduced. It was not.

### Where the pair does and does not close

The model-dependent step is less clean. The prediction above was made with the
existing run's snapshot times used for both members of the pair, since the twin
had not been run; recomputing each with its own clock, which ends at 12.00 h
rather than 11.99 h, moves `T` for the narrower ridge from -0.3137 to -0.3130 per
cent per km and the predicted separation from -0.9500 to -0.9508. Against that,
the measured separation is -1.1288, a mismatch of -0.1780, which
is 2.1 times the combined snapshot scatter of 0.0861. Under the parameter-free
front assumption the residual is therefore not common to the pair, at
`D = +0.0211` per cent per km for `sigma = 100` km and `+0.1992` for
`sigma = 50` km.

The natural place for that mismatch is the sharp-front assumption itself, and the
twin makes the reason visible. The modes still unarrived at 55 km at 12 h are the
same two in both runs, `n = 1` at 26.8 h and `n = 2` at 13.4 h, with identical
arrival times; only their flux share differs, 0.4524 against 0.1381. In the
narrower case the whole correction over the band reduces to one step of about
14 per cent at the top, since the arrived fraction takes only three values there,
1.0000 below 20 km, 0.9517 over 24.75 to 44.75 km and 0.8619 above 49.75 km.
Replacing a dispersively spread front by a step does most of the work in that
configuration, and the symptom is visible in the twin's own Section 8: the
correction overshoots, moving the drift from -0.1138 to +0.2003 per cent per km
and raising the peak-to-peak from 8.66 to 13.13 per cent, with a front factor of
0.534 nulling the drift against 0.989 for the wider ridge.

Solving for the common front factor that reconciles the pair gives `f = 1.112`,
at which `D` is `+0.2570` and `+0.2573` per cent per km. With two measurements
and two unknowns this is a solve and not a test, so it should be read only as the
distance from the parameter-free `f = 1` that the front assumption has to be
pushed, together with the `D` that follows.

Taken together the pair places the numerical flux drift over 10 to 55 km between
0.02 and 0.26 per cent per km, depending on how the front sharpness is treated.
That is a tightening from the roughly 1 per cent per km the single run supported,
and the bound now sits below the 0.1 per cent per km reporting threshold at one
end of its range and a little above it at the other.

### What this does and does not transfer

The ordering of the two residuals at `f = 1` is the one result that would be
unfavourable if it survived, since `D` is larger for the narrower ridge, and the
nonlinear runs carry their flux at higher wavenumber still, at a scale that
collapses from about 980 to 45 km under saturation. The pair does not settle
whether that ordering is real or an artifact of the front treatment, because a
common front factor within 12 per cent of unity removes it entirely. Separating
the two requires a measurement that does not depend on the fill-in model at all,
which is what a second `dx` would provide.

## 14. Extending the twin to 15 h, and what it excludes

Sections 10 and 13 both rest on one account of why the measured drift falls with
time. Each ridge wavenumber arrives at its own moment, and the drift is the
shortfall of the part that has not. That account was fitted to the record it was
built from, so it was not tested by it. A relaxation `a exp(-t/tau) + D_inf`,
which carries no wave physics, describes the same twelve points about as closely.
The two agree over the record and separate beyond it, so the twin was integrated
from 12 to 15 h and both were scored out of sample. Section 10 of
`scripts/analyze_linear_flux_constancy.py` carries the comparison, with
`--forecast-from` setting the cut.

### The forecast, fixed before the run

Both accounts were fitted to `t <= 12` h and projected forward before the
integration was launched. The fill-in account supplies `T(t)` from the spectrum
and fits only the constant `D`, taken from the 11 to 12 h mean as in Section 13;
the relaxation fits three parameters over each of five start hours. Drifts are
the slope of the mean profile over the hour ending at `t`, over 10 to 55 km,
which is how the headline figures are reported.

| t [h] | fill-in | relaxation family | measured |
|---|---|---|---|
| 13 | -0.0435 | -0.0887 to -0.0302 | -0.0644 |
| 14 | +0.0317 | -0.0579 to +0.0301 | -0.0383 |
| 15 | +0.0367 | -0.0379 to +0.0762 | -0.0332 |

Two corrections to the framing that motivated the extension are worth recording,
because both were established by running the model rather than by argument. The
first is that the arrived fraction being frozen at 0.8619 between 9 and 13.4 h
holds at the evaluation top and not across the band. Mode `n = 2` enters the band
at 2.84 h and reaches 55 km only at 13.39 h, so its front is climbing through the
column throughout, crossing one `dz = 500` m level about every 0.1 h, and `T`
moves from -0.401 to -0.283 per cent per km over 9 to 12 h rather than standing
still. The second follows from it. The step in `T` as `n = 2` clears the band top
is +0.0145 per cent per km, not the +0.2 that a step in the arrived fraction at
one height would suggest, because the band-fitted slope responds to where fronts
sit inside 10 to 55 km rather than to the fraction arrived at its top. The
separation between the two accounts is therefore carried by the 12 to 15 h
increment, +0.150 against +0.067 per cent per km, and not by a step.

### The run

`julia -t auto --project=. scripts/sweep_arsia_mons_stage2.jl linear u=40 h=0.1
hours=15 sigma=50`, into `output/arsia_mons_linear15h_LX2000_sig50/`. Wall-clock
3694.87 s (61.6 min) for 73266 steps to `t_end = 54000` s, converged, with the
timestep flat at 0.737 to 0.740 s and no reduction. 733 snapshots, 2.5 GB.
Extraction 32.5 s, analysis 1 s. Against 2587.66 s and 58609 steps for the 12 h
twin, the step count scales with the horizon as expected and the cost per step is
about 14 per cent higher, which is consistent with the analysis work sharing the
machine.

The integration is an extension rather than a repetition. All 113 snapshots the
two runs share are bit-identical across every stored field, so hours 0 to 12 are
the same integration and the out-of-sample hours are scored against the same
solver state that produced the forecast.

Controls are unchanged. The minimum `min_x u` below the sponge is +14.529 m/s, or
0.363 `U`, so the overturning criterion is not met; the largest `|w|` below the
sponge is 8.529 m/s; the mean flow at 55 km stays within 0.320 m/s of `U`; and
the static instability at `t_end` is confined to 58.75 to 61.75 km, with zero
unstable area fraction everywhere in 10 to 55 km.

### Neither account survives

The relaxation is refused by the hour-on-hour increments, which are free of both
the anchor for `D` and the choice of fit window.

| interval | fill-in | relaxation family | measured |
|---|---|---|---|
| 12-13 h | +0.0703 | +0.0472 to +0.0793 | +0.0494 |
| 13-14 h | +0.0753 | +0.0307 to +0.0604 | +0.0260 |
| 14-15 h | +0.0050 | +0.0200 to +0.0460 | +0.0052 |

At the fitted `tau` of 2.3 to 2.5 h a relaxation requires successive increments to
fall by `exp(-1/tau) = 0.67` per hour, which over 14 to 15 h is +0.020 to +0.046
per cent per km. The measured increment is +0.0052 against a per-snapshot spread
of 0.0030 in the final hour, and the half-hourly trace over that hour reads
-0.0331, -0.0314 and -0.0349, so the change is not distinguishable from zero. The
drift has stopped falling, and it has stopped at about -0.035 rather than at the
`D_inf` of -0.000 to +0.223 the relaxation family projects.

The fill-in account is refused by the level. Its out-of-sample error is -0.070
per cent per km at both 14 and 15 h, about 23 times the per-snapshot spread and
outside the relaxation family's own range at 14 h. Over the three scored hours
its rms error is 0.0584 against 0.0123 for the median of the relaxation family.
Applied to the 15 h record the correction now moves the drift away from zero
rather than toward it, from -0.0332 to +0.1293, so the fraction of the drift it
accounts for is -290 per cent, against -76 per cent at 12 h. The front factor
that would null the corrected drift falls to 0.567.

That last figure retires the argument the front factor was carrying. The value
was quoted as a consistency check rather than a fit because `f` is not free, and
a null at 0.989 for the wider ridge was read as the group velocity placing the
fronts where the record put them. The same statistic is 0.534 for the twin at
12 h and 0.567 for the same twin at 15 h. A quantity that moves over that range
between runs and between horizons of one run is behaving as a fitted parameter,
so the `f = 1.112` reconciliation of Section 13 should be read as what it was
described as, a solve rather than a test, and the `f = 0.989` null in Section 10
should not be read as support for the model. The statistic that 97.0 per cent of
the drift is arrival transient belongs to the same family and is not quotable on
its own; its companion values are -76 per cent at 12 h and -290 per cent at 15 h.

The one thing the fill-in account gets right out of sample is the 14 to 15 h
increment, +0.0050 predicted against +0.0052 measured. Both the model and the
record have stopped moving by then, for reasons that need not be the same, so
this is weak evidence.

### What the extension does establish

The measured flux is flat at the stated reporting thresholds for the first time.

| over 10-55 km, final hour | 12 h | 15 h |
|---|---|---|
| fitted drift [%/km] | -0.1138 | -0.0332 |
| peak-to-peak / median | 0.0866 | 0.0420 |
| accumulated across the band [%] | -5.12 | -1.49 |
| band median MF [uPa] | -21.194 | -20.246 |
| pointwise drag, band RMS [m/s/sol] | 0.1427 | 0.0815 |
| flat at the 5 % and 0.1 %/km thresholds | no | yes |

Three hours of further integration reduced the drift by a factor of 3.4 and the
peak-to-peak by 2.1, on a spectrum, grid, column, wind and sponge that did not
move. Whatever part of the 12 h drift was numerical was therefore at most
-0.0332 per cent per km, since a numerical flux divergence fixed in height would
not have fallen with the horizon. That is a model-free upper bound and it is
inside the 0.1 per cent per km reporting threshold, which is the result the
extension was run for.

### The divergence is matched by mean-flow deposition

A separate check bears on what the residual is. Over 10 to 55 km, which lies
entirely below the 60 km sponge base, the flux divergence and the measured
acceleration of the `x`-mean flow agree throughout.

| window | band RMS from `dMF/dz` | band RMS from `du/dt` | ratio | correlation | residual |
|---|---|---|---|---|---|
| 9-10 h | 0.2231 | 0.2260 | 1.013 | 0.9964 | 7.6 % |
| 11-12 h | 0.1427 | 0.1504 | 1.054 | 0.9988 | 7.2 % |
| 13-14 h | 0.0986 | 0.1107 | 1.123 | 0.9983 | 13.5 % |
| 14-15 h | 0.0815 | 0.0863 | 1.058 | 0.9991 | 7.1 % |

All in m/s per sol. Over the final hour the residual is 0.0058 m/s per sol rms
with a mean of -0.0013, and the sign is that the acceleration slightly exceeds
the divergence rather than falling short of it. The momentum associated with the
measured flux divergence is arriving in the mean flow, which is still
accelerating at 15 h at a band rms of 0.086 m/s per sol.

This bounds loss and not misplacement, and the distinction matters for how the
diagnostic is described. A scheme that removed momentum from the system would
open a gap between the two columns, and no gap larger than 7 per cent is present.
A scheme that transported momentum to the wrong height would conserve it and
would deposit it into `du/dt` exactly as physical deposition does, so the
comparison cannot detect that mode of error at all. The drift in `MF(z)` remains
the only handle on it, and remains an upper bound rather than a measurement.
On the same reasoning the band-integrated divergence should not be described as
spurious drag while the run is still filling; over 10 to 55 km at 15 h it is
deposition into a mean flow that has not finished responding.

### What is now open

The 15 h record leaves the drift settled at about -0.035 per cent per km with
the run still evolving, since the flux median is still falling at about 1.5 per
cent per hour and the mean flow is still accelerating. Whether that -0.035 is a
numerical floor or the current value of a slower transient is not settled by this
pair of horizons, and neither account of the time dependence is available to
extrapolate it, both having been excluded. A resolution comparison remains the
only probe that does not depend on such an account, and it now has a cleaner
target than before, since at 15 h the quantity to be resolved is a drift already
inside the reporting threshold rather than one several times outside it.

## 15. Up and down decomposition of the 15 h record

Section 2 measures one projection of the flux profile, its height gradient, and
Section 3 measures the other, its level. A loss that is uniform in height moves
the level and leaves the gradient untouched, so the drift is blind to it by
construction. Over 12.0 to 14.7 h the band flux falls at 1.42 per cent per hour
while the drift stalls at -0.035 per cent per km, which is that combination. The
statement in Section 14 that the run has settled therefore rests on one
projection of a two-projection object, and the other one was still moving when
the record ended.

Below the sponge the geometry admits one channel for a height-uniform loss, a
wave returning from above. Splitting the flux by the sign of the vertical
wavenumber measures that channel directly on the record already in hand, at no
simulation cost. It is Section 11 of
`scripts/analyze_linear_flux_constancy.py`. It rereads 60 of the 733 snapshots,
which is the only bulk snapshot access in that script, and takes the whole
analysis from 1 to 2 s to 4.7 s of wall clock.

### How the split is made

With `rho_0(z)` the x-mean density, the density-scaled Fourier amplitudes
`Psi_k = FFT_x[sqrt(rho_0) w']_k` and `Xi_k = FFT_x[sqrt(rho_0) u']_k` satisfy the
steady anelastic continuity relation `dPsi/dz = -i k Xi + Psi/(2H)`, so the
vertical derivative is read from the stored `u` rather than differenced in `z`.
Writing `Psi = A exp(imz) + B exp(-imz)` with `m > 0` the local vertical
wavenumber gives `A exp(imz) = (Psi + Psi'/(im))/2` and
`B exp(-imz) = (Psi - Psi'/(im))/2`, and substituting into the cospectrum leaves

```
Re(Xi Psi*) = -(m/k) (|A|^2 - |B|^2)
```

at every height, the `AB*` cross terms cancelling identically in the real part.
Three consequences matter for what follows.

The partition is pointwise in `z`. It needs no vertical transform and no window,
which is what makes it usable here: the evaluation band holds about 1.6 vertical
wavelengths at the band-median `lambda_z` of 30.7 km, and a windowed FFT split
over 45 km would leak the upward peak into the downward half at a level
comparable to the signal.

The `1/(2H)` term drops out of the flux exactly, so the result does not depend on
how the density scale height is estimated.

The flux carries no `AB*` cross term. A partial standing wave therefore produces
no modulation in `MF(z)` at all, and the interference appears only in the
variance. Structure in the flux residual at `lambda_z/2` is consequently not a
reflection signature, and the test has to be run on `|Psi|^2` instead.

### What fixes the floor of the method

The section builds a synthetic WKB mountain wave with the same ridge spectrum,
density and `N(z)` and an imposed reflected fraction, and measures it before it
is applied to the run. On a purely upward wave the method returns
`8.4e-05` in flux, set by the WKB amplitude gradient rather than by noise. Imposed
fractions of 0.05 and 0.10 are returned as 0.05004 and 0.10002. That synthetic is
built from the same `m` the split assumes and so cannot see an error in `m`
itself; Section 18 repeats the check on an exact anelastic solution and finds a
systematic `2.4e-04`, which is the floor the later sections are read against.

The remaining dependence is on `m`. An error of a fraction `eps` misassigns about
`(eps/2)^2` of the upward flux, which is quadratically small, and the run
confirms it: taking `m` from the instantaneous mean state gives `R = 0.0486` at
the final hour, taking it from the background profile with `U = 40` m/s gives
0.0491, and multiplying `m` by 0.90, 0.95, 1.05 and 1.10 gives 0.0527, 0.0499,
0.0486 and 0.0497. The drift of the upward component moves over +0.035 to +0.055
per cent per km across the same variants. The net flux does not depend on `m` at
all, and its drift is -0.0348 per cent per km in every row, which is a check on
the partition rather than a result.

Two terms sit outside the partition and are reported with it. The sum of the two
parts is `<rho_0 u'w'>_x` identically, so the closure test is against
`<rho u'w'>_x`, and the largest level difference at the final hour is 1.14 per
cent of the flux, which is the triple correlation `<rho' u'w'>_x`. Evanescent
wavenumbers, `m^2 = (N/U)^2 - k^2 <= 0`, carry `1.9e-04` of the flux, against the
zero a non-propagating mode should give.

### A downward component is present and grows through the record

| t [h] | MF up [uPa] | MF down [uPa] | net [uPa] | R = down/up | down share |
|---|---|---|---|---|---|
| 8.00 | -21.216 | +0.527 | -20.681 | 0.0248 | 2.42 % |
| 11.99 | -21.852 | +0.896 | -20.947 | 0.0410 | 3.94 % |
| 14.74 | -21.187 | +1.030 | -20.146 | 0.0486 | 4.64 % |

The half-hourly trace over 6.0 to 14.99 h runs from 0.0187 to 0.0494 at a linear
rate of 0.0037 per hour, with no sign of settling within the record.

Its height structure is that of a source above the band. At 14.74 h the local
ratio rises from 0.037 at 9.75 km to 0.079 at 54.75 km, reaches 0.18 at 59.75 km
just below the sponge base, and tends to 1.00 above 70 km, where the net flux is
0.16 uPa and what remains is a standing residual of negligible amplitude. The
entries above 60 km are indicative only, since `m` is complex where the wave is
being damped.

The independent check in the variance is consistent with it. A superposition of
an upward and a downward wave modulates `|Psi|^2` at `2m`, that is at
`lambda_z/2 = 15.3` km, with a relative depth `2r/(1 + r^2)` for `r = sqrt(R)`.
The measured cross term peaks at a vertical scale of 15.0 km at both 8 and 12 h,
against a predicted depth of 0.31 and 0.39 and a measured peak-to-peak of 0.47
and 0.75. At 14.74 h the dominant scale of that residual moves to 45 km, the band
length, so the signature is present at two of the three sampled hours rather than
at all three.

The `sigma = 100` km twin gives a third handle. At 11.74 h it carries `R = 0.0138`
against 0.0410 for the `sigma = 50` km run at 11.99 h, and its flux has reached
only 0.644 of the steady linear value against 0.837. Less of its spectrum has
reached the sponge by that time, and the reflected share is about a third as
large, which is the direction a top-boundary origin predicts.

### The residual drift changes sign when the returning wave is removed

At 14.74 h the drift of the net flux over 10 to 55 km is -0.0348 per cent per km
and that of the upward component alone is +0.0445. Both are inside the 0.1 per
cent per km reporting threshold, but they differ in sign, so the residual band
drift of Section 14 is at least as consistent with the height structure of the
returning wave as with a flux divergence in the incident one. On this record the
two are not separable.

The sub-band picture is less stable than either.

| band | `<rho u'w'>` | `<rho_0 u'w'>` | upward | weighting spread |
|---|---|---|---|---|
| 10-25 km | -0.0799 | -0.1170 | -0.0828 | +0.0371 |
| 25-40 km | -0.0416 | +0.0228 | +0.0597 | -0.0644 |
| 40-55 km | -0.1951 | -0.3367 | -0.2696 | +0.1416 |
| 10-55 km | -0.0348 | -0.0353 | +0.0445 | +0.0005 |

The last column is the difference between weighting the covariance by `rho` and
by its x-mean, a choice with no physical content at this order. Over the whole
band it is 0.0005 per cent per km, so the -0.035 the previous section quotes is
insensitive to it. Over 40 to 55 km it is 0.142 per cent per km, comparable to
the entry beside it and to the entire residual under discussion. The top sub-band
slope is therefore not a measurement at the level at which it has been read, and
the observation that the full-band slope is smaller in magnitude than each of its
parts is a statement about a summary statistic rather than about the flux.
Removing the returning wave does not remove the non-monotonicity: the upward
component alone still gives -0.083, +0.060 and -0.270 across the three sub-bands.

### Reflection accounts for a fifth of the level loss, not for most of it

The decay of the level splits as follows over 12.0 to 14.7 h. The incident
component falls at 1.123 per cent per hour and the net at 1.417, with the
downward component growing at 0.049 uPa per hour, so the growth of the returning
wave accounts for 16.9 per cent of the loss and the remaining 83 per cent is a
decay of the incident wave itself.

That remainder is not a weakening source. The x-mean wind at 0.25, 1.75, 4.75 and
9.75 km changes by at most 0.0010 m/s across 8, 12 and 15 h, so the launch
condition of the stationary wave is unchanged to one part in `4e4` while its flux
falls by about 1.1 per cent per hour.

The same split applies to the amplitude. Against the steady linear reference of
-25.017 uPa the net flux falls short by 19.47 per cent and the upward component
alone by 15.31, so the returning wave accounts for 21.4 per cent of the
shortfall. The reference carries the factor of 2.08 `N_0` ambiguity of Section 9,
so the remaining 15.31 per cent should not be read as a loss without settling
that reference first. What the decomposition does establish is that reflection is
not large enough to absorb the shortfall on its own.

### What the sponge looks like to the wave

The layer is 30 km thick, 0.98 vertical wavelengths at the band-median
`lambda_z`, which is the number Section 19 item 3 has been comparing against the
usual guidance. The damping the wave meets is sharper than that number suggests.
The profile is `K(z) = K0 sin^2(pi/2 (z - 60)/30)` from `src/sponge.jl`, and the
run uses `K0 = 0.05` s^-1, set at `scripts/sweep_arsia_mons_stage2.jl:129`. The
default in `src/parameters.jl:87` is `0.005` s^-1, recorded there as the
Srivastava optimal value, so the run is an order of magnitude above it.

For a stationary wave the relevant comparison is against the intrinsic frequency
`omega_hat = U k`, which for the flux-carrying modes is `2.5e-04` to `1.5e-03`
s^-1.

| lambda_x [km] | omega_hat [1/s] | K0/omega_hat | z where K = omega_hat [km] | above base |
|---|---|---|---|---|
| 1000.0 | 2.513e-04 | 198.9 | 61.36 | 0.044 `lambda_z` |
| 666.7 | 3.770e-04 | 132.6 | 61.66 | 0.054 `lambda_z` |
| 400.0 | 6.283e-04 | 79.6 | 62.15 | 0.070 `lambda_z` |
| 250.0 | 1.005e-03 | 49.7 | 62.72 | 0.089 `lambda_z` |
| 166.7 | 1.508e-03 | 33.2 | 63.33 | 0.109 `lambda_z` |

The wave crosses from undamped to strongly overdamped within 1.4 to 3.3 km of the
sponge base, that is within 0.04 to 0.11 of a vertical wavelength. On the scale
the wave responds to, the layer is closer to a step than to the gradual 30 km
ramp its thickness suggests, which is the configuration in which a Rayleigh layer
is expected to reflect rather than absorb. At `K0 = 0.005` s^-1 the same crossing
would sit at 64.3 to 71.1 km, 0.14 to 0.36 `lambda_z` above the base. These are
scale arguments rather than a reflection calculation, and they identify a
parameter to vary rather than establish that it is responsible.

### What this changes and what it leaves open

The bound of Section 14 is unchanged as a statement about the profile shape.
Whether it should be read as a bound on the solver is now weaker, because the
upward component alone drifts at +0.0445 rather than -0.0348 per cent per km,
so the sign of the residual depends on which of the two waves is being asked
about.

The claim that the run has converged does not survive. One projection settled and
the other did not: the band flux was still falling at 1.42 per cent per hour at
the end of the record, of which about a sixth is the reflected share growing and
the rest a decay of the incident wave under an unchanged source. A height-uniform
multiplicative loss is invisible to `d ln|MF|/dz` by construction, and this is
one.

Neither the reflected share nor the incident decay had settled by 15 h, so the
horizon does not resolve them, and the headroom for extending it is limited: by
15 h `min_u_over_U` below the sponge has fallen to 0.363 from 0.543 at 12 h, and
the base of the unstable layer has descended from 59.75 km at 12 h to 58.75 km at
14 and 15 h.

## 16. Per-mode timing of the downward component

Section 15 measures a downward-propagating component and attributes it, on
geometry and on scale arguments about the damping layer, to reflection at the
sponge base. That attribution was not tested. It can be, at no simulation cost
and without any new quantity, because `updown_split()` already forms `fup` and
`fdn` per horizontal wavenumber and only their sum over `k` was being reported.
The section is `--updown-modes` and `--updown-mode-heights` of
`scripts/analyze_linear_flux_constancy.py`. It rereads no snapshots beyond the 60
Section 15 already reads. The whole analysis takes about 13 s of wall clock
against the 4.7 s recorded for Section 15; the per-mode reporting itself accounts
for none of that, 13.1 s with it against 13.9 s with `--updown-modes 0`, and the
increase is the second partition below, which is formed for every snapshot. The
same table is then run on the `sigma = 100` km twin of Section 13, which shares
the `k` grid but weights it very differently, at the same zero simulation cost.

### Why the timing separates the candidates and the amplitude does not

The rfft grid inside `updown_split()` is the grid `ridge_spectrum()` builds, so
ridge mode `n` is column `n` of both partitions. Each mode has its own clock. A
stationary hydrostatic component of wavenumber `k` climbs at `cg_z = k U^2 / N`,
so it reaches height `z` at `t_up = I(z)/(k U^2)` with `I` the buoyancy integral
already used by the fill-in benchmark, and a component reflected at the sponge
base reaches the same height only after climbing to that base and descending
again,

    t_rt = (2 I(z_spg) - I(z)) / (k U^2).

Because `cg_z` scales with `k`, the two are separated by hours for the long
modes, and the record ends inside the gap. At 9.75 km the round trips of `n = 1`,
2 and 3 are 53.15, 26.57 and 17.72 h against a record of 14.99 h, while their
upward arrivals are 5.67, 2.84 and 1.89 h. Reflection at 60 km cannot place a
downward component in those three modes at that height within this record. Any
process that generates one in the interior, or any leakage from the steady
continuity closure the partition is built on, would appear at `t_up` with no
round-trip delay.

The predictions were fixed before the table was produced. A sponge origin
requires `R_k` at the method floor for `n = 1, 2, 3` at 10 km and `R_k > 0` for
`n >= 4`, and it orders `R_k` with `K0/omega_hat`, which runs from 33.2 at
`n = 12` to 397.9 at `n = 1` and so increases with `lambda_x`. An interior origin
requires `R_k` nonzero wherever `t_up` has elapsed, with no step at `t_rt`.

### What fixes the floor for a single mode

The floor quoted in Section 15 is a share-weighted sum over the spectrum and
does not bound any single column, so it was remeasured per mode on the same
synthetic purely upward WKB wave. It is `2e-05` at 9.75 and 29.75 km and `1e-05`
at 54.75 km. Those numbers understate the floor, because that synthetic shares the
split's `m`; the exact-solution floor of Section 18 is `2.4e-04` in every
propagating mode, about two orders of magnitude below the values reported here
rather than three.

That synthetic is a steady field, so it exercises the algebra but not the one
approximation the partition contains. `dPsi/dz` is taken from steady continuity,
which drops `d(rho')/dt`, and the solver is compressible, so that term is not
identically zero. Differencing `Psi` in `z` gives an independent estimate of the
same derivative, accurate to about 0.2 per cent at `dz = 500` m against the
band-median `lambda_z` of 31 km. Writing the relative difference between the two
estimates as `eps`, an error of that size puts `B = -eps Psi / 2` and so
manufactures a downward fraction of about `eps^2/4` in every mode, the same
quadratic law the `m_scale` sensitivity of Section 15 measures.

| z [km] | eps at 8 h | eps at 12 h | eps at 15 h | eps^2/4 |
|---|---|---|---|---|
| 9.75 | 0.169 | 0.100 | 0.108 | 0.0029 |
| 29.75 | 0.102 | 0.094 | 0.121 | 0.0037 |
| 54.75 | 0.202 | 0.161 | 0.135 | 0.0045 |

The partition rebuilt on the differenced derivative is carried alongside the
first one. At 9.75 km at 14.74 h it gives 0.0160 against 0.0202 for the modes
that cannot have received a reflection, and 0.0301 against 0.0345 for those that
can, so the two derivative estimates agree on the downward amplitude to about 12
per cent. Neither the floor nor the closure accounts for the values below.

### The measurement

Modes are classed on their own clock at each sampled hour rather than on the
wall. A mode counts as established at a height once the sample time exceeds 1.5
times its own transit to that height, below which its amplitude is still rising
and the ratio of two small numbers is not a measurement. Of the established
modes, those whose round trip has elapsed could have received something from the
sponge and those whose has not could not. Membership therefore changes from hour
to hour, and the count is given with each entry.

Share-weighted mean `R_k`, with the standard error across the modes in the class:

| z [km] | class | 8 h | 12 h | 14.74 h |
|---|---|---|---|---|
| 9.75 | no return possible | 0.0073 +- 0.002 (5 modes) | 0.0162 +- 0.006 (4) | 0.0202 +- 0.004 (3) |
| 9.75 | returned | 0.0321 +- 0.010 (6) | 0.0288 +- 0.009 (8) | 0.0345 +- 0.008 (9) |
| 29.75 | no return possible | 0.0044 +- 0.002 (3) | 0.0129 +- 0.010 (2) | 0.0515 +- 0.037 (2) |
| 29.75 | returned | 0.0310 +- 0.010 (7) | 0.0330 +- 0.007 (9) | 0.0385 +- 0.007 (9) |
| 54.75 | returned | 0.0520 +- 0.017 (7) | 0.0523 +- 0.007 (9) | 0.0535 +- 0.013 (10) |

The per-snapshot scatter of these means over the 11 snapshots of the final hour
is 0.0002 and 0.0006 at 9.75 km and 0.0015 and 0.0011 at 29.75 km, so the
mode-to-mode spread quoted above is the larger of the two error bars and is the
one to carry. At 54.75 km the round trip is only about 1.2 times the upward
transit, so no established mode remains unreturned and the class is empty; the
test has power only well below the sponge base, which is where the round trip is
long.

### The prediction that a sponge origin makes is not met

The three modes at 9.75 km whose round trips are 53.15, 26.57 and 17.72 h carry
a share-weighted `R_k` of 0.0202 at 14.74 h, against a per-mode floor of `2.4e-04`
and a closure bound `eps^2/4` of 0.0029. The value is a factor of 7 above the
larger of the two bounds and about 5 times the mode-to-mode standard error,
though that error is formed on three modes and should be read as indicative. The
same three modes carry 0.0073 at 8 h and 0.0162 at 12 h, so the observation is
not confined to one sample. Reading the same three columns
individually at 14.74 h gives 0.0263, 0.0117 and 0.0242, so it is not carried by
a single mode.

These results are not consistent with the sponge being the only source of the
downward component. They are consistent with a source that also acts below the
sponge base, whether distributed through the column, sited near the launch
levels, or in the interior scheme. This record does not separate those.

A step at the round trip is nonetheless present. At 9.75 km the returned class
exceeds the unreturned class by 0.0143 at 14.74 h, 0.0345 against 0.0202, and by
0.0248 and 0.0126 at 8 and 12 h. Taken at face value that would put about 41 per
cent of the asymptotic ratio at the round trip and the remainder below it. Two
things weaken that reading. The two classes differ in `lambda_x` as well as in
timing, since the unreturned modes are the long ones, so the step and any scale
dependence are confounded. And at 29.75 km the ordering reverses, 0.0515 against
0.0385, though on two modes and with an error bar of 0.037. The `sigma = 100` km
twin reverses it as well, which is taken up below.

### The scale ordering and the spectral shape

If the damping layer returns the wave, `r_n` should increase with `lambda_x`,
since `K0/omega_hat` does. Fitting `r_n` against `log10 lambda_x` over the
returned modes, weighted by flux share, gives the slope per decade:

| z [km] | sigma = 50 km, modes | slope | sigma = 100 km, modes | slope |
|---|---|---|---|---|
| 9.75 | 9 | -0.107 | 4 | -0.038 |
| 29.75 | 9 | +0.008 | 5 | -0.035 |
| 54.75 | 10 | +0.065 | 5 | -0.274 |

Five of the six entries are zero or negative, which is not the ordering a
damping layer would produce. The two runs agree at 9.75 km and disagree in sign
at the two upper heights, where the `sigma = 100` km fits rest on four or five
modes and the `sigma = 50` km fits at 54.75 km are drawn from modes whose round
trip is barely distinguishable from their arrival. The scale ordering is
therefore evidence against a damping-layer origin below about 30 km and is not
resolved above it.

A caveat limits how far the slope can be read. `R_k` is a ratio and its
denominator follows the ridge spectrum, which falls by a factor of 10.2 across
the `sigma = 50` km modes and 139 across the `sigma = 100` km modes, so a
downward component with no mode dependence at all would produce a rising `r_n`
at short `lambda_x` on that ground alone. Fitting the numerator against the
denominator directly, as `|MF_dn,k| = C |MF_up,k|^alpha`, separates the two. A
reflection carries the incident spectrum and gives `alpha = 1`; a component with
no memory of the ridge gives `alpha = 0`.

| z [km] | sigma = 50 km | sigma = 100 km |
|---|---|---|
| 9.75 | +0.31 +- 0.21 | +0.85 +- 0.20 |
| 29.75 | +0.64 +- 0.34 | +0.74 +- 0.21 |
| 54.75 | +0.57 +- 0.44 | +0.91 +- 0.43 |

Five of these six estimates are within 1.5 standard errors of 1 and the sixth,
`sigma = 50` km at 9.75 km, is 3.2 below it. The `sigma = 100` km fits have the
longer lever arm, 139-fold against 10-fold in the denominator, so on the
better-conditioned pair the downward component carries the incident spectrum.
The data do not support reading the downward component as broadband in `k`, and
the negative slopes in the table above are not attributable to the ratio artefact
alone.

### The same measurement on the sigma = 100 km twin

The twin of Section 13 runs on the same grid and the same `k` grid but weights
those modes very differently, `n = 1` carrying 0.182 of the steady flux against
0.048 for `sigma = 50` km. It therefore replicates the test with the flux moved
onto the modes whose round trip is longest. Its record ends at 11.74 h, so the
comparison is against the 11.99 h column of the narrower ridge.

| quantity, at 9.75 km | sigma = 50 km, 11.99 h | sigma = 100 km, 11.74 h |
|---|---|---|
| no return possible | 0.0162 +- 0.006, 4 modes, share 0.391 | 0.0193 +- 0.008, 4 modes, share 0.866 |
| same at 8 h | 0.0073 +- 0.002 | 0.0079 +- 0.002 |
| returned | 0.0288 +- 0.009, 8 modes, share 0.588 | 0.0137 +- 0.002, 4 modes, share 0.134 |
| closure bound eps^2/4 | 0.0029 | 0.0046 |
| no-return class, differenced derivative | 0.0160 at 14.74 h | 0.0195 |

The unreturned class replicates. It sits at 0.008 at 8 h and 0.016 to 0.019 at
about 12 h on both runs, above the closure bound in both, and the two derivative
estimates agree on the twin to within 1 per cent. A downward component in modes
that cannot have received a reflection from 60 km is present in two runs whose
flux weighting differs by a factor of 3.8 on the longest mode.

The step at the round trip does not replicate. On `sigma = 50` km the returned
class exceeds the unreturned class, 0.0288 against 0.0162; on `sigma = 100` km
the ordering reverses, 0.0137 against 0.0193. Whatever produces the difference
between the classes on one run does not reproduce on the other, so the 41 per
cent attribution above should be treated as an upper bound on what the round trip
contributes rather than as a measurement.

This also bears on how Section 15 read the twin. There the twin's lower band `R`,
0.0138 against 0.0410, was taken as support for a top-boundary origin on the
grounds that less of its spectrum had reached the sponge. The per-mode table
gives the same observation a different reading. The unreturned class is the same
on both runs, and the band ratios differ because the compositions differ: 86.6
per cent of the twin's flux sits in modes that could not have returned by 11.74 h,
against 39.1 per cent for the narrower ridge at 11.99 h. The twin comparison
measures where each spectrum sits relative to the round-trip horizon, which is
what the reassessment of Section 15 anticipated, and it does not by itself
localise the source.

One limit on the twin. At 54.75 km its `eps` reaches 0.277, so the closure bound
`eps^2/4` is 0.019 against a returned-class mean of 0.027, a factor of 1.4 rather
than the factor of 7 available at 9.75 km. The twin's upper-level entries are not
separated from the closure and are not used above.

### R is filling rather than growing, at a lower level than inferred

Because `R_k` is already a per-mode ratio, a mode that has returned carries no
filling of its own and its `r_n` is final. The returned-class mean is therefore
an estimate of the asymptote that needs no correction, and it is close to
constant in time at every height: 0.0321, 0.0288, 0.0345 at 9.75 km, 0.0310,
0.0330, 0.0385 at 29.75 km, and 0.0520, 0.0523, 0.0535 at 54.75 km across 8, 12
and 14.74 h. The reading that the band `R` grows because the column is filling,
rather than because the reflectivity is rising, is supported by that constancy.

The level is lower than the aggregate correction gives. Correcting the band `R`
by the share of the flux that has had time to return, which assumes every mode
reflects alike, gives an implied `r` of 0.0454, 0.0538 and 0.0585 at the three
hours. Measuring `r_n` directly on the modes that did return, and averaging the
three heights, gives 0.0384, 0.0381 and 0.0422. The two part by about 40 per
cent, and the direct measurement is the one that carries no assumption about the
modes still in transit.

The height dependence is also weaker than the band ratio suggests. Measured `R`
at 14.74 h is 0.0371, 0.0480 and 0.0786 at the three heights, a factor of 2.12
across the band. The returned-class asymptote is 0.0345, 0.0385 and 0.0535, a
factor of 1.55. Part of the difference is that modes which have not yet arrived
at the upper levels contribute to the numerator of the band ratio while
contributing almost nothing to its denominator. At 54.75 km `n = 1` has not
arrived at all, its upward transit being 26.77 h, and its `R_k` of 2.30 is a
ratio of two small numbers; `n = 2` arrived 1.35 h before the sample and reads
0.51. Both enter the band sum. The `R(z)` profile of Section 15 should be read
with that in mind above about 40 km.

A residual height dependence survives the correction. A single reflector at 60
km would give an `R` independent of height in density-scaled WKB variables, and
0.0345 to 0.0535 over 45 km is not that. The observation raised at the end of
Section 15 therefore stands in weakened form. The corrected gradient is about
0.98 per cent per km in `r`, which under reciprocity would require a symmetric
attenuation of about 0.49 per cent per km on each leg, against the 0.67 per cent
per km the uncorrected profile implies.

### What this changes

The claim that the returning wave originates at the damping layer is not
supported by the timing on either run, and the scale ordering below 30 km runs
the other way on both. What Section 15 established is unchanged: a downward
component is present, it is measured well above the floor of the method, and the
band ratio grows because the column is filling rather than because the
reflectivity is rising. What is not established is where it comes from.

The `K0` variation proposed in Section 19 therefore addresses at most part
of the signal, and its prediction should be restated before it is run. If the
sponge is the source, `R_k` in the returned class should fall while `R_k` in the
unreturned class at 9.75 km is unchanged. A run that moves both, or neither,
would separate the sponge from whatever produces the 0.0202. Comparing the two
runs at equal filled fraction rather than equal clock remains necessary, and the
per-mode table supplies that directly, since the class means are already indexed
on each mode's own transit rather than on the wall.

Two limits on the reading. A mode is counted only once the sample time exceeds
1.5 times its own transit, and the classification inherits the sharp-front
idealisation of the fill-in benchmark, which neglects the dispersive spread of
each component about its own front. And a downward component left over from the
start-up adjustment would sit in the unreturned column; the growth of that column
at 9.75 km, from 0.0073 to 0.0202 on the same three modes between 8 and 14.74 h,
is not the decay a settling transient would show, but it is also not the
step-then-plateau an interior source switching on at `t_up` would give, the
latest of those three arrivals being 5.67 h. Neither candidate accounts for the
time dependence of that column on this record. Section 17 revisits that
column with the whole record binned at one hour and class membership recomputed
in each bin, rather than sampled at two hours, and does not resolve a trend in it
on either run, so the time dependence quoted here is better read as the scatter
of a quantity that fluctuates.

The figure is `linear_flux_permode.png`. Its left panel is the test rather than a
result, showing where each mode's arrival and round trip fall against the end of
the record.

## 17. Whether the downward component is stationary in the ground frame

### The assumption the per-mode classification rests on

Section 16 classes a mode by `t_rt = (2 I(z_spg) - I(z))/(k U^2)`, the round trip
of a stationary mountain component. Such a component has ground-based frequency
`omega = 0` and intrinsic frequency `omega_hat = -kU`, and both the transit time
and the vertical wavenumber `m = sqrt(N^2/U^2 - k^2)` that the split of Section 15
applies follow from that. Anything present at the same `k` and a different
frequency departs from both. Its vertical group speed is `omega_hat^2/(N k)`
rather than `k U^2 / N`, so its round trip is the tabulated one multiplied by
`(kU/omega_hat)^2`, and the split, applying the mountain `m` to it, reports a
downward share of

    R_sp = ((1 - r)/(1 + r))^2,    r = kU/|omega_hat|,

from it even when the component is purely upward, since a pure upward wave with
true wavenumber `m_t` gives `A = Psi (1 + m_t/m)/2` and `B = Psi (1 - m_t/m)/2`.
`R_sp` rises toward 1 as `r` moves away from 1 in either direction, so the
misassignment is not a marginal effect: a component far enough from the mountain
frequency is divided about equally between the two branches.

The unreturned class of Section 16 is therefore unreturned only for waves at the
mountain frequency. That is a structural limit of the classification rather than
a defect in its measurement, and it is the one place where the argument of
Section 16 could fail without any of its numbers being wrong. The record already
written is sufficient to test it, at no simulation cost.

### What the transform already carries

A component of the physical field going as `exp(i(kx - omega t))` has an rfft
coefficient in `x` going as `exp(-i omega t)`. A stationary field therefore holds
`Psi_k` fixed in time, and with it `A_k` and `B_k`, while a component at any other
frequency rotates them at `omega = -d arg / dt`. The complex amplitudes are
already formed at every level inside `updown_split()`, so the only addition is to
keep them rather than the flux they sum to. Two statistics follow.

The first is a partition of each branch. The flux is quadratic in the amplitude,
so over any window

    <|B|^2> = |<B>|^2 + var(B),

which splits the downward flux exactly into a part that an amplitude fixed over
the window would carry and a part that it could not. The stationary round-trip
clock governs the first and says nothing about the second. The same partition
applied to `A` is the control, since the incident wave is still filling and the
mean flow still adjusting, so some non-steadiness is common to both branches.

The second is the rotation rate itself, taken as `omega = -arg<B_{i+1} B_i*>/dt`
over the trailing window, with the modulus of the same average normalised by the
geometric mean of the two energies reported alongside as a narrowbandness. That
normalisation rather than a division by one of the energies keeps a growing
amplitude from pushing the diagnostic past 1.

This is Section 12 of `scripts/analyze_linear_flux_constancy.py`, added with
`--transient-stride`, `--transient-window` and `--transient-bins`. It is the only
part of that script whose cost scales with the length of the record rather than
with a fixed number of hour windows, since it reads the whole run at a stride
rather than three windows of it.

### The measured rates do not bring any unreturned mode inside the record

At 9.75 km on the `sigma = 50` km run the three modes of the unreturned class
rotate at `omega/kU` of -0.287, 1.470 and 0.692, that is at ground-based periods
of order tens of hours rather than at the fractions of an hour a fast transient
would give. The correction runs the wrong way for the loophole in two of the
three: `|omega_hat|/kU` is 1.279, 0.470 and 0.301, so `n = 2` and `n = 3` descend
more slowly than the tabulated round trip assumes, not faster. Their corrected
round trips are 32.1, 120.5 and 187.3 h against a 14.99 h record.

| run | shortest corrected round trip | record | unreturned entries slowed |
|---|---|---|---|
| `sigma = 50` km | 32.1 h (`n = 1` at 9.75 km) | 14.99 h | 3 of 4 |
| `sigma = 100` km | 21.4 h (`n = 1` at 9.75 km) | 11.99 h | 4 of 6 |

The closest case is `n = 1` on the twin, whose corrected round trip is 1.8 times
its record rather than the 4.4 the tabulated value gives. On neither run does the
correction bring any unreturned mode inside the record, so the classification of
Section 16 is not overturned by the frequencies actually present.

### About half the unreturned ratio is carried by an amplitude that does not change

Splitting each branch over the last three hours separates the part of the ratio
the stationary clock governs from the part it does not.

| run | z | class | R | steady | fluctuating | fluctuating, upward branch |
|---|---|---|---|---|---|---|
| `sigma = 50` km | 9.75 km | no return | 0.01463 | 0.00969 | 0.00494 | 0.00511 |
| `sigma = 50` km | 9.75 km | returned | 0.03266 | 0.02796 | 0.00471 | 0.00672 |
| `sigma = 50` km | 29.75 km | no return | 0.01023 | 0.00707 | 0.00317 | 0.00873 |
| `sigma = 100` km | 9.75 km | no return | 0.01390 | 0.00949 | 0.00442 | 0.00958 |
| `sigma = 100` km | 29.75 km | no return | 0.00610 | 0.00270 | 0.00340 | 0.01006 |

The steady column is the result. At 9.75 km it is 0.0097 and 0.0095, against a
per-mode method floor and a closure bound `eps^2/4` of 0.0029 and 0.0045
respectively, so it stands at 3.3 and 2.1 times the closure bound and about 40
times the floor, the floor being the `2.4e-04` of Section 18 rather than the
`1.7e-05` this section originally quoted. The conclusion of Section 16 survives the stationarity
test at roughly half the level Section 16 quotes for the same class, the earlier
0.0202 being a share-weighted mean at one hour rather than a flux-weighted mean
over three.

The last two columns are the reason the reduction is not simply a loss of signal.
Across the ten class entries measured on the two runs the fluctuating part of the
downward branch is 0.21 to 0.97 of the fluctuating part of the upward branch,
with a median of 0.42, while the steady parts stand in a ratio of 0.0027 to
0.0504. Entry by entry the fluctuating component is divided between the two
branches 8 to 125 times more evenly than the steady one, with a median of 39. A
reflection returns a small fraction of the incident wave and cannot produce that;
a component whose
vertical wavenumber is not the one the partition assumes produces exactly that,
by the `R_sp` relation above. The fluctuating part therefore carries no
information about which way anything is going, and removing it from the ratio is
a correction rather than a truncation.

### The growth reported in Section 16 is not a trend

Section 16 recorded, as unexplained by either candidate origin, that the
unreturned class at 9.75 km ran 0.0073 at 8 h to 0.0202 at 14.74 h. Those are two
`+/- 0.5` h means of eleven snapshots each. The third panel of the new figure
shows why that is not enough: `|B|` at 9.75 km in `n = 1`, `2` and `3` wanders by
factors of several on hour timescales throughout the record, while `|A|` in the
same modes rises smoothly and settles. Binning the whole record at one hour and
recomputing class membership in each bin gives the following.

| run | bins | mean | sd | range | linear trend |
|---|---|---|---|---|---|
| `sigma = 50` km | 10 | 0.01419 | 0.00545 | 0.00783 to 0.02281 | +0.00103 +/- 0.00052 per h, 2.0 se |
| `sigma = 100` km | 10 | 0.01686 | 0.00614 | 0.00878 to 0.02569 | -0.00061 +/- 0.00068 per h, 0.9 se |

The `sigma = 50` km sequence is not monotone: it runs 0.0123, 0.0095, 0.0078,
0.0097, 0.0121, 0.0220, 0.0228, 0.0114, 0.0143, 0.0200 across 5-6 h to 14-15 h.
The two hours Section 16 quotes are the second-lowest and the highest of the ten.
On the twin the fitted trend is negative and within one standard error of zero.
The data are therefore consistent with a quantity that fluctuates with a standard
deviation of about 0.0055 about a constant mean of about 0.014 to 0.017, and the
apparent growth is better read as two draws from that fluctuation than as a time
dependence. The item left open at the end of Section 16 is closed on that
reading, and the bin-to-bin scatter, `+/- 0.0017` on the mean of ten bins, is the
error bar that a pair of hour means several hours apart does not carry.

The one part of the record that does show a monotone time dependence is the
fluctuating part, which falls from 0.00304 to 0.00014 across the ten bins on
`sigma = 50` km and from 0.01110 to 0.00106 on the twin, tracking the fluctuating
part of the upward branch throughout. That is the settling of the start-up
adjustment, and it is settling rather than growing.

### No descent and no association with the overturning layer

If the fluctuating part were radiated near the sponge base it would fall away
from it with a measurable delay. Cross-correlating its class-summed downward flux
against the topmost level at which the class is not empty gives a lag of zero
within one sample interval over 39.25 to 29.25 km on `sigma = 50` km, at
coherences of 0.96 to 0.98. With a 221 s interval and a 10 km separation that
bounds any propagation from below at 45 m/s, against the 1.16 to 1.27 m/s
stationary group speed of `n = 2`, the only mode in that class over that depth.
The fluctuation is therefore quasi-simultaneous across a 10 km layer rather than
descending through it. Below 29 km the coherence falls to 0.33 to 0.86, the lags
run to thousands of seconds, and on the twin they change sign between adjacent
levels, so no propagation speed is measurable there.

Correlating the same downward flux at 9.75 km against the unstable area fraction
over the whole record gives, with everything slower than 55 min removed, 0.052 on
`sigma = 50` km and 0.112 on the twin for the class sum, and -0.062 to 0.160 and
-0.385 to 0.186 for the individual modes. The raw coefficients are larger, 0.334
and 0.293, but the upward branch is the control here, since it is filling and
cannot be driven by the layer, and it correlates more strongly than the downward
branch does on both runs, 0.705 and 0.889. What raw correlation there is, is
therefore a shared trend. The unstable fraction is also not monotone, peaking at
0.00367 at 13.27 h and ending at 0.00207 on `sigma = 50` km, a 44 per cent
decline over which the downward ratio does not follow it.

### What this changes

The classification of Section 16 stands, and its central result stands at a
revised level. Modes whose sponge round trip has not elapsed, on either run and
after correcting each mode's round trip for its own measured rotation rate, carry
a time-independent downward component of 0.0097 and 0.0095 at 9.75 km, two to three
times the closure bound and about 40 times the per-mode floor of Section 18. Neither a sponge
reflection, on the timing, nor the closure of the split, on that margin, accounts
for it.

Three things previously carried as open or as inference are now measured. The
growth of the unreturned class is not resolved as a trend on either run, so the
one feature Section 16 could not attribute to any candidate is removed from the
list rather than explained. Roughly a third of the ratio Section 16 quotes is a
fluctuating component divided nearly evenly between the two branches, which the
partition cannot classify and which should not be read as returning flux. And
neither the descent a source near the sponge base requires nor an association
with the overturning layer is present in the record at the resolution available.

What remains to be accounted for is narrower than what Section 16 left open: a
downward component that is steady in the ground frame, present at about the same
level on two runs whose spectral weighting of the affected modes differs by a
factor of 3.8, not growing, not descending, and not tracking the instability.
Those are the properties of something in the steady solution the solver settles
into rather than of a propagating return, which is what makes an amplitude
variation the discriminating next run and is why Section 18 reports it. That
section also corrects the floor these margins are quoted against.

The figure is `linear_flux_transient.png`. Its left panel is the test of Section
16's classification, the round trip of each mode recomputed from its measured
rotation rate against the end of the record; its middle panel is the quantity
Section 16 quotes, binned rather than sampled at two hours; its right panel is
the raw material of both.

## 18. The amplitude control, and a correction to the floor

### A correction to the floor of the method

The floor quoted in Sections 15 to 17, `8.4e-05` in the band and `1.7e-05` per
mode, comes from `synthetic_wkb()`, which builds its field from the same
`m = sqrt(N^2/U^2 - k^2)` that the split assumes. It therefore exercises the
algebra, the conjugate-pair weighting and the WKB amplitude gradient, but it
cannot see an error in `m` itself, and there is one.

At constant `N` and constant density scale height `H` the anelastic vertical
structure equation for `Psi = sqrt(rho_0) w'` is

    d2Psi/dz2 + [N^2/U^2 - k^2 - 1/(4H^2)] Psi = 0,

so the exact vertical wavenumber carries a term the split omits. Over 10 to 55 km
in this column the band-mean `N` is 0.00849 s^-1 and the band-mean scale height
9.61 km, which puts `1/(4H^2)` at 6.00 per cent of `(N/U)^2`. The assumed `m` is
then low by 3.0 per cent, and by the same quadratic law the `m_scale` sensitivity
of Section 15 measures, a purely upward exact solution is reported with a
downward share of `((1-r)/(1+r))^2` with `r = m_exact/m`.

Building that exact solution and passing it through the split gives `2.438e-04`
over the band and `2.394e-04` to `2.556e-04` across modes 1 to 12, matching the
analytic value to five digits, while the field's own flux is height-independent
to 0.0000 per cent, which is the check that the construction is right. The check
is now reported alongside the WKB one as `updown/exact_check`.

This is a systematic bias, the same sign in every mode, so it does not average
down over modes or over time, and it is the floor a per-mode ratio should be read
against rather than the `1.7e-05`. Nothing in Sections 15 to 17 reverses, but the
margins quoted there should be read as follows: the 0.0202 of Section 16 stands
at 84 times this floor rather than three orders above it, and the 0.0097 of
Section 17 at 40 times rather than 570. Both remain above the closure bound of
0.0029, which is the binding limit in any case, and the correction leaves the
timing arguments untouched, since they turn on when a mode can have received
anything rather than on how small the floor is.

### The prediction, and the run

`R` is a dimensionless ratio, and every candidate still on the list after Section
17 except one is a linear operator acting on the wave. Sponge reflection, the
lower boundary condition, the residue of the start-up adjustment and truncation
error on either axis all scale the numerator and the denominator together, so
they leave `R` unchanged when the forcing amplitude changes. A source that is
nonlinear in the wave amplitude does not. That makes the amplitude the one
variable that separates the two classes of explanation without a second grid, a
second sponge or a longer record.

The run is `h = 30` m at `sigma = 50` km, everything else identical, launched
2026-09-03 as

```
julia --project=. scripts/sweep_arsia_mons_stage2.jl linear sigma=50 hours=15 h=0.03
```

and analysed with `--h-max 30`. It took 3684.6 s of wall clock over 73055 steps
against 3695.1 s and 73266 steps for the `h = 100` m run, confirming that the
timestep is set by the vertical acoustic limit and not by the wave. The
prediction was fixed before it finished and is recorded in Section 19 of the
version of this document written that morning: that the steady unreturned-class
ratio at 9.75 km either holds at 0.0097 within the `+/- 0.0017` the bin scatter
of Section 17 supports, in which case the origin is linear, or falls toward the
0.0029 closure bound, in which case it is not.

### The controls behave as they should

Before reading the result, the quantities that ought to be amplitude-independent
are. Writing each as `h^p` with the incident flux scaling as `h^2`, so that `p = 2`
means the ratio is unchanged:

| quantity | h = 100 m | h = 30 m | implied p |
|---|---|---|---|
| fluctuating downward part, 9.75 km no return | 0.00494 | 0.00360 | 2.26 |
| fluctuating upward part, 9.75 km no return | 0.00511 | 0.00470 | 2.07 |
| fluctuating downward part, 9.75 km returned | 0.00471 | 0.00492 | 1.96 |
| fluctuating upward part, 9.75 km returned | 0.00672 | 0.00654 | 2.02 |
| fluctuating downward part, 54.75 km returned | 0.00502 | 0.00469 | 2.06 |
| fluctuating upward part, 54.75 km returned | 0.01298 | 0.01291 | 2.00 |

Five of the six sit within 0.07 of `p = 2`. The fluctuating component that
Section 17 identified as something the partition divides evenly between its
branches is therefore amplitude-independent, which is consistent with reading it
as a property of the partition and of the start-up rather than as a source, and
it establishes that the comparison has the sensitivity to detect a departure from
`p = 2` when there is one.

The steady linear reference scales as `h^2` exactly, `-25.017` to `-2.252` uPa, a
ratio of 11.111. The measured incident flux scales as 10.28, `-21.187` to
`-2.061` uPa at 14.74 h, so the shortfall against that reference falls from 15.3
to 8.5 per cent when the amplitude is reduced. About half the amplitude shortfall
of Section 9 is therefore itself amplitude-dependent and is not a numerical loss.

### The steady unreturned component is not amplitude independent

| z, class, quantity | h = 100 m | h = 30 m | ratio | implied p |
|---|---|---|---|---|
| 9.75 km, no return, steady | 0.00969 | 0.00199 | 4.87 | 3.31 |
| 29.75 km, no return, steady | 0.00707 | 0.00100 | 7.04 | 3.62 |
| 9.75 km, returned, steady | 0.02796 | 0.02372 | 1.18 | 2.14 |
| 29.75 km, returned, steady | 0.03306 | 0.01737 | 1.90 | 2.53 |
| 54.75 km, returned, steady | 0.05040 | 0.03134 | 1.61 | 2.39 |

The prediction is answered, and against the branch recorded as the alternative.
The steady part of the unreturned class falls by a factor of 4.87 at 9.75 km and
7.04 at 29.75 km for a 3.33-fold reduction in amplitude, giving exponents of 3.31
and 3.62 against the 2.00 that every linear mechanism requires. At `h = 30` m it
reads 0.00199 and 0.00100, both below the closure bounds of 0.0029 and 0.0038 at
those heights, so it is no longer measurable there.

Two amplitudes cannot resolve the exponent uniquely, but they do fix a
two-parameter reading. Writing the ratio as an amplitude-independent part plus a
part quadratic in `h`, which is what a source quadratic in the wave field would
give since its radiated flux would then go as `h^4` against the incident `h^2`,

| z | amplitude-independent part | quadratic part at h = 100 m |
|---|---|---|
| 9.75 km | 0.00123 +/- 0.00112 | 0.00846 |
| 29.75 km | 0.00040 | 0.00667 |

The offset is `1.0989 R(30) - 0.0989 R(100)`, so it is carried almost entirely by
the smaller-amplitude run and inherits its scatter. Propagating the bin-to-bin
standard errors of the steady part, 0.00101 at `h = 30` m and 0.00156 at
`h = 100` m over ten one-hour bins each, gives the error bar above, which puts
the offset 1.1 standard errors from zero. It is therefore not resolved, and the
29.75 km entry, whose scatter was not separately estimated but is of the same
order, is not resolved either. What is resolved is the quadratic part: on this
reading 87 and 94 per cent of what Section 17 measured at `h = 100` m is
amplitude-dependent, and that share is insensitive to where inside its error bar
the offset actually sits.

The returned class behaves differently, and the difference is informative. At
9.75 km its exponent is 2.14, that is amplitude-independent to 14 per cent, which
is what a linear return channel gives. At 29.75 and 54.75 km it is 2.53 and 2.39,
so the two contributions are mixed there, which is where the wave is larger and
the sponge base nearer. The band ratio `R` falls from 0.0486 to 0.0269 and the
local ratio at 10, 30 and 55 km from 0.0371, 0.0480 and 0.0786 to 0.0257, 0.0219
and 0.0739, so the height dependence within the band largely disappears while the
value just below the sponge is almost unchanged.

### The linearity of the pair, and one prediction that does not hold

The run was also intended to supply the linear control this record has been
treating itself as.

| | h = 100 m | h = 30 m | linear projection |
|---|---|---|---|
| minimum `u/U` below 60 km, all times | 0.363 | 0.762 | 0.809 |
| peak mean-flow deviation [m/s] | -2.150 | -0.206 | -0.194 as `h^2` |
| mean-flow deviation at 55 km [m/s] | 0.320 | 0.033 | 0.029 as `h^2` |
| 40 to 55 km weighting spread [per cent per km] | 0.142 | 0.057 | 0.043 as `h` |

The mean-flow response follows `h^2` to within 6 per cent, and the minimum wind
over the hill recovers to 0.762 of `U` against a linear projection of 0.809, the
6 per cent shortfall being the sense a saturating `h = 100` m run would give. The
weighting spread, which is the direct measure of the triple correlation
`<rho' u'w'>_x` and so of the nonlinearity of the diagnostic itself, falls by 2.5
against the 3.3 a linear scaling would give, so it too is closer to linear at the
smaller amplitude without being exactly so.

One prediction fails, and it matters for the reading of Section 17. The static
instability fraction was expected to go to zero and instead rises, from 0.00207
to 0.00750 at the end of the record. Resolving it in height explains why and
withdraws the interpretation it had been given: in both runs the `N^2 < 0` points
sit in a narrow layer at 58.75 to 61.75 km straddling the sponge base at 60 km,
seven levels at `h = 100` m and five at `h = 30` m, and there are none below it.
That diagnostic is therefore measuring a static instability produced at the top
of the damping layer rather than wave overturning in the free column, and it does
not scale with the wave. The null result of Section 17, that the downward flux
does not track it, is unaffected, but the layer should not be described as the
overturning of the mountain wave, and it is not the nonlinearity this section
detects.

Note added 2026-09-23. `reflection_convergence_tests.md` Sections 9 and 11 reclassify
this layer. In the wave-free Z0 run under `pbar_mode = :total`, the x-mean `N^2` at
60.25 km falls to -152 % of its initial value by 6 h, with no wave present, and in Z1
(`:pressure`) the column holds to round-off. The data are consistent with the
58.75 to 61.75 km instability being a product of the kinetic-energy term in `P_bar`
rather than of the damping layer or the wave.

### What linear theory accounts for in the offset

Two contributions to the unresolved offset can be computed rather than fitted,
and both are linear and so belong there rather than in the quadratic part.

The first is the partition bias of the first subsection, `2.44e-04` in every
propagating mode. The second is partial reflection off the represented `N(z)`
itself, which is a distributed interior reflector and was left on the candidate
list by Section 16 because the timing does not exclude it. Integrating
`Phi'' + [N^2/U^2 - k^2 - s^2 + s'] Phi = 0` for `Phi = sqrt(rho_0) w'`, with
`s = -d ln(rho_0)/dz / 2`, downward from the sponge base with a pure upgoing
radiation condition and decomposing at the base of the band gives a reflected
share of `4.87e-04` over the three modes of the unreturned class at 9.75 km,
with no solver and no sponge in the calculation. The two together come to
`7.31e-04`, which sits inside the one standard error of the fitted offset. The
offset is therefore consistent with being entirely linear, though the fit does
not have the precision to say so.

Two limits on that reading. The reflected share is not a cumulative quantity: at
a base of 29.75 km it is `7.89e-04`, larger than over the longer path from 9.75
km, because the contributions from different parts of the column interfere. The
fitted offsets run the other way, 0.00123 at 9.75 km against 0.00040 at 29.75 km,
so the two do not track each other in height, which is a further reason not to
read the agreement at 9.75 km as an accounting. And the reflected share depends
on how `N(z)` is represented rather than only on the data: the `4.87e-04` uses
the centred-difference `N` on the 500 m grid, which is the `N` the split itself
forms, whereas linearly interpolating the same 500 m MCD column onto a finer grid
makes `dT/dz` a staircase and raises the flux-weighted share to `1.84e-03`,
converged from 250 m to 62.5 m. The quantity is therefore of order `5e-04` to
`2e-03` depending on the representation, and only the first of those is the one
the run can be compared against.

### What this changes

The result removes a class of explanations rather than naming one. A downward
component whose steady part scales as `h^3.3` to `h^3.6` cannot be produced by
sponge reflection, by the lower boundary condition, by the residue of the
start-up adjustment or by truncation error on either the horizontal or the
vertical axis, since each of those is a linear operator on the wave and would
leave the ratio unchanged. That is the whole of the candidate list Section 16
left open, and it is the reason this run was placed ahead of the grid and sponge
variations rather than after them.

What it leaves is a nonlinearity, not yet located. The obvious candidate, wave
overturning in the free column, is not supported by the one diagnostic available,
since the static instability in both runs lives at the sponge base and moves the
wrong way with amplitude. A second reading consistent with the numbers is that
the returned and unreturned classes are not separate channels at all below about
30 km, the returned class being 2.14 at 9.75 km and so nearly linear while the
unreturned class is 3.3, which would put the nonlinear part in the interior and
the linear part above.

The practical consequence for the rest of the record is narrower and firmer. The
`h = 30` m run is the linear control, and on it the band ratio is 0.0269 rather
than 0.0486, the unreturned class is unmeasurable against the closure bound at
both heights, the mean-flow response is 0.21 m/s rather than 2.15, and the
weighting spread over 40 to 55 km is 0.057 per cent per km rather than 0.142. Any
statement about the solver's flux error bar should be made on that run rather
than on the `h = 100` m one, whose residual is now known to be substantially a
property of its own amplitude.

## 19. Proposed next checks

1. Done on 2026-09-02, recorded in Sections 13, 14 and 15, and extended on
   2026-09-03 in Sections 16, 17 and 18. The `sigma = 50` km twin was run first,
   then extended from 12 to 15 h, and the 15 h record was split by the sign of
   the vertical wavenumber at no further simulation cost. The extension left the
   drift at -0.0332 per cent per km and the peak-to-peak at 4.20 per cent, both
   inside the reporting thresholds, and excluded both accounts of the time
   dependence that had been available to extrapolate it. The decomposition then
   showed that the band flux was still losing 1.42 per cent of itself per hour in
   a way the drift cannot see. Section 16 resolved that decomposition by
   horizontal wavenumber, again at no simulation cost, and found that the
   attribution of the downward component to the damping layer is not supported by
   the timing on either run. Section 17 tested the assumption that classification
   rests on and found that the measured rotation rates do not bring any
   unreturned mode inside either record, that about a third of the ratio is a
   fluctuating component the partition divides nearly evenly between its two
   branches, and that the remaining steady part neither grows, descends, nor
   tracks the static instability. Section 18 then varied the forcing amplitude
   and found that steady part is not amplitude-independent, scaling as `h^3.3` to
   `h^3.6` against the `h^2` every linear mechanism requires, which removes the
   whole of the candidate list Section 16 left open and leaves a nonlinearity
   that is not yet located. Section 18 also corrected the floor of the method
   from `1.7e-05` to `2.4e-04`, the earlier synthetic having shared the split's
   own `m`.

2. Locate the nonlinearity, which is now the open question and is the cheapest
   remaining item. Three things belong in it and none needs a new configuration.

   A third amplitude constrains the exponent, which two points cannot, and the
   useful direction is downward. `h = 10` m at `sigma = 50` km and 15 h costs
   about 3690 s, the timestep being set by `dz`. There the two readings of
   Section 18 part: an amplitude-independent 0.00123 plus a part quadratic in `h`
   predicts a steady unreturned-class ratio at 9.75 km of 0.00131, while a single
   exponent, `R` going as `h^1.32`, predicts 0.00047, a factor of 2.8. Both lie
   above the `2.4e-04` exact-solution floor and below the 0.0029 closure bound,
   so the reading turns on how tight that bound is, and the `h = 30` m point
   indicates it is not tight, its measured 0.00199 already sitting below it.
   Going the other way is not available. Extrapolating the measured perturbation,
   which grows slightly more slowly than `h` between 30 and 100 m, the minimum
   wind over the hill reaches about 0.04 `U` near `h = 150` m and passes through
   zero shortly above, so any amplitude much above 100 m overturns and the
   comparison stops being a clean amplitude variation. At `h = 30` and 100 m the
   two accounts agree by construction and at `h = 60` m they differ by 14 per
   cent, so an intermediate point is not worth its hour.

   A wave-overturning diagnostic that excludes the damping layer is needed
   because the one in hand does not. Section 18 finds the `N^2 < 0` points of
   `unstable_frac_t` confined to 58.75 to 61.75 km in both runs, straddling the
   sponge base, with none below, and rising rather than falling when the
   amplitude is reduced. Restricting the same fraction to below 55 km, and
   reporting the minimum `N^2` rather than the fraction below zero so that a
   column that is stabilising is distinguishable from one that never overturned,
   would make it a measurement of the wave rather than of the layer. That is an
   analysis change with no simulation cost.

   And the nonlinear term should be looked for where it would act. The triple
   correlation `<rho' u'w'>_x` is already reported as the difference between
   weighting the covariance by `rho` and by its x-mean, 0.142 per cent per km
   over 40 to 55 km at `h = 100` m against 0.057 at `h = 30` m. Resolving that
   difference by horizontal wavenumber, as Section 16 resolved the split, would
   show whether it sits in the same modes as the unexplained downward component
   and would cost one further pass over the snapshots.


3. Vary the sponge at fixed grid, same `sigma`, same 15 h, and on the `h = 30` m
   configuration rather than on `h = 100` m. Section 18 changes what this run is
   for. The unreturned class is now known to be substantially amplitude-dependent
   and so is not a sponge product, but the returned class at 9.75 km has an
   amplitude exponent of 2.14, that is amplitude-independent to 14 per cent,
   which is what a linear return channel gives and is the signature the sponge
   would have. The returned class, at 0.0237 on the `h = 30` m run, is therefore
   the quantity this variation should be read on, and running it at the smaller
   amplitude removes the nonlinear contribution from the comparison. The cheapest
   variant changes
   `K0_sponge` from the `0.05` s^-1 of `scripts/sweep_arsia_mons_stage2.jl:129`
   to the `0.005` s^-1 default of `src/parameters.jl:87`, which moves the height
   at which the damping rate overtakes the intrinsic frequency from 1.4 to 3.3 km
   above the sponge base to 4.3 to 11.1 km, that is from 0.04 to 0.11 vertical
   wavelengths to 0.14 to 0.36. It costs about the same as the 3695 s of the 15 h
   run, since neither the grid nor the timestep changes.

   Sections 16 to 18 narrow what this run can settle. The prediction is not that
   `R` should fall, but that the returned class should fall, since only it can
   carry anything the damping layer returns and only it has the amplitude
   signature of a linear channel. A run that moves both would indicate
   that `K0` is changing something other than the reflection at the layer, which
   it plausibly does, since it also changes the absorption depth and so the
   mean-flow response and the descending unstable base, which makes it a
   two-variable change rather than the clean probe it looks like. A run that
   moves neither would remove the sponge as a candidate. The comparison should be
   made at equal filled fraction rather than equal clock, because a 15 h run at
   the new `K0` is also mid-fill, and the per-mode class means already supply
   that, being indexed on each mode's own transit. The bin scatter of Section 17
   fixes what counts as a change: `+/- 0.0017` on a ten-bin mean, so a pair of
   hour means is not sufficient to compare the two runs.

4. Refine the vertical grid before the horizontal one, `dz = 250` m at
   `sigma = 50` km and 15 h. Halving `dz` halves the acoustic timestep as well as
   doubling the cell count, so this costs about four times the 3695 s. It is
   placed ahead of the `dx` pair because the vertical axis is the less resolved
   of the two for exactly the modes that carry the unexplained signal: at
   `dx = 4` km, `n = 1`, `2` and `3` sit at 500, 250 and 167 points per horizontal
   wavelength, while `lambda_z = 30.7` km on `dz = 500` m is 61 points.

5. Repeat at `dx = 2` km last. Two arguments in the record now place it there
   rather than second. Horizontal truncation error enters as `(k dx)^2`, so it
   would put a downward spectrum proportional to `(k dx)^2 |MF_up,k|`, whose
   spread across `n = 1` to `12` exceeds the upward spread by a factor of
   `12^2 = 144`. The measured downward spreads are 7.81, 17.75 and 26.75 at 9.75,
   29.75 and 54.75 km against upward spreads of 10.16, 10.51 and 18.38, so the
   largest excess anywhere in the record is a factor of 1.7. Vertical truncation
   error enters as `(m dz)^2`, and since `k^2/m^2` is at most 0.029 across these
   modes, `m` varies by 1.5 per cent among them and that error is
   spectrum-preserving, so the same test does not exclude it. Second, the cost is
   about twice rather than four times the 3695 s, because the timestep is set by
   the vertical acoustic limit: with `c` of 173 to 230 m/s over the column,
   `dz/c` is 2.17 s and `dx/(c + U)` is 14.8 s, against the run's `dt` of 0.737 s,
   so halving `dx` doubles the cell count at unchanged `dt`. That makes it cheap
   but not informative on present evidence, and a null result would be the
   expected one.

   If it is run, the prediction is that an interior horizontal-truncation origin
   drops the steady unreturned-class ratio at 9.75 km by about a factor of four,
   from 0.0097 to the closure bound, while a boundary, start-up or physical
   origin leaves it unchanged. The comparison target is a drift of -0.033 per
   cent per km rather than the -1.24 the single 12 h run at `sigma = 100` km
   supported, and the momentum budget is known to close to 7 per cent, so a
   change in the drift cannot be confused with a change in what the flux is
   depositing.

6. Settle the amplitude reference by solving the steady linear problem with the
   actual `N(z)` rather than the WKB launch formula, integrating the vertical
   structure equation for each `k` with a radiation condition aloft. That would
   remove the factor of 2.08 identified in Section 9 and is the only item here
   that would make the `MF/MF_linear` ratio quotable.

   The source of that factor is now identified and it is a single number.
   `linear_steady_flux()` normalises the WKB amplitude as `(rho_0 N)^-1/2` at the
   launch level, so the reference carries `N` at the first grid point, where the
   MCD column has a temperature inversion of `+11.23` K/km against a dry adiabat
   of `-4.66` K/km. Evaluating that function with the observed `N(z)` and with a
   constant band-mean `N` gives a ratio of 2.0759, and `N(0)/N_band` is 0.017631
   over 0.008493, which is 2.0759 to four decimals. The ambiguity is therefore
   entirely the choice of launch `N` and not the WKB step, which narrows this
   item to fixing that choice rather than replacing the method. Section 15 raises the value
   of doing so, because the shortfall now splits into 4.2 points that the
   returning wave accounts for and 15.3 points that it does not, and the second
   number cannot be read as a loss until the reference is fixed. Evaluating the
   comparison at the boundary rather than at the first cell centre, and including
   the compressible correction to `m`, are smaller terms within the same
   calculation.

   Item 3 was run on 2026-09-22 together with a thickness and a base-height variant,
   and is recorded in `docs/reflection_convergence_tests.md`, which also carries a
   linear prediction of the layer's reflection and an estimate of the reflected share
   from the variance ripple that is independent of the split.

7. What no variant of items 3 to 5 can settle on its own is the level. Even a
   null result on the reflected share leaves the incident flux decaying at 1.12
   per cent per hour under a launch wind that is constant to 0.0010 m/s, and that
   decay is height-uniform, which is the one error mode `d ln|MF|/dz` cannot
   register. Tracking `MF_up` rather than the drift is the projection that sees
   it, and it is reported by hour.
