# Mechanism of the descending instability layer, Stage 1 U=40/h=3 extensions

*Note dated 2026-10-02.* Every run in this report imposes the ground `w` through the ghost cells (`bottom_forcing = :ghost`, the only path before 2026-09-30). Absolute flux levels and amplitudes quoted here are superseded for amplitude: at CFL 0.4 the ghost path loses 10.7 % of the transmitted flux in uniform wind against the face path, by an amount first order in `dt` (`reflection_convergence_tests.md` Section 14, addendum 2 and addendum 3 step 1). Shapes, timings and ratios between runs at equal `dt` are not covered by this note.

Internal analysis report, 2026-08-31, extended 2026-09-01 with the 12 h record.
Sections 1 to 10 describe the 6 h extension and are unchanged. Sections 11 to 16
report the 12 h extension, which tests the Section 7 prediction and revises part
of the Section 9 reading. The `ext6h` stem in the file and script names is
historical and now covers both horizons.

Original report header follows. Companion to `docs/TODO_tracker.md`
(section "Arsia Mons Stage 1 U=40/h=3 time extension") and to
`docs/arsia_stage1_U57h3_deep_dive.md`.

Question addressed: the static-instability fraction in the hill window below
60 km rises from 0.000 at 1 h to 0.086 at 6 h and the region in which it
appears moves downward from about 58 km to about 30 km. This report asks what
produces that pattern, whether it is a property of the resolved dynamics or of
the configuration, and how far it can be reproduced from analytical results.

Reproduce with `python3 scripts/analyze_ext6h_mechanism.py [--npz PATH]`
(about 1.5 s wall-clock; reads `deepdive_series.npz` and the MCD background
column, does not touch the NetCDF snapshots). The sampled hours are derived from
the series, so the same script serves either horizon. Metrics are written to
`mechanism_metrics.json` in the matching `..._analysis/run_U40.0_h3.0/`
directory, `arsia_mons_stage1_ext6h_analysis` or `..._ext12h_analysis`.

## 1. Framework

For a stationary wave the intrinsic phase speed is the mean wind itself, so
three quantities reduce to one criterion. Wave-action conservation gives the
displacement amplitude zeta ~ (rho N)^-1/2, the polarization relation gives
u_hat = N zeta, and convective overturning occurs when |d(zeta)/dz| >= 1.
Combining these,

    overturning  <=>  u_hat >= U  <=>  min_x u <= 0.

The corresponding saturated momentum flux (Lindzen 1981) is

    F_sat(z) = rho k U^3 / (2 N),

so a wave carrying flux F is saturated above the level where |F| = F_sat. Both
forms are evaluated below with the solver's own constants (gamma = 1.31,
g = 3.72 m/s^2, R_gas = 188.9 J/(kg K)), the run's MCD background column, and a
flux-weighted ridge wavenumber k_eff = 9.01e-6 1/m (lambda_x = 697 km) obtained
from the k|h(k)|^2 spectrum of the Gaussian ridge.

## 2. Two wave populations are present

The steady hydrostatic mountain wave has vertical group velocity
cg_z = k U^2 / N, which for this configuration is 1.5 to 1.8 m/s. Integrating
that through the background column, the stationary wave front reaches 20 km at
3.6 h, 30 km at 5.1 h and 50 km only at 8.5 h. The run ends at 6 h, so the
stationary wave has not established above about 37 km at any point in the
record.

The measured flux fill-in times run in the opposite direction with height. The
first time |MF(z)| exceeds 25% of its own 6 h value is 1.55 h at 10 km, 2.01 h
at 20 km, 1.53 h at 50 km and 1.16 h at 57 km. An inverted arrival sequence is
consistent with the upper-column signal at t < 2 h being the startup transient
of the impulsively switched-on lower boundary condition
(`src/forcing.jl:115`, no temporal envelope), which is broadband and carries
higher intrinsic frequencies and therefore larger cg_z, rather than the
stationary wave.

The steady linear hydrostatic solution for this ridge carries
|MF| = 22.23 mPa, height-independent to machine precision. The simulation
reaches 10.87 mPa at 5 km by 6 h, about 49% of that value, so the lower
boundary is still spinning up at t_end.

## 3. Amplitude growth is consistent with conservative WKB scaling

Ratios of hill-window RMS w between band endpoints, against the same ratio for
(rho N)^-1/2:

| t (h) | 10-30 km | 30-45 km | 45-59 km |
|---|---|---|---|
| 1.0 | 0.51 | 0.55 | 1.01 |
| 2.0 | 0.60 | 0.66 | 0.70 |
| 3.0 | 0.66 | 0.55 | 0.72 |

At 1 h the 45-59 km band grows at the conservative rate to within 1%. Below
45 km the growth is 50 to 66% of that rate, which is the expected signature of
a wave field that has not yet filled the column rather than of dissipation.
After about 1.9 h the 45-59 km agreement degrades, which coincides with the
onset of overturning in that band (Section 4) and therefore with the point at
which conservative scaling ceases to apply.

## 4. Onset progresses downward from the top of the sponge-free column

First time each level satisfies the overturning criterion, and first time it
registers N^2 < 0:

| z (km) | 59.75 | 58.75 | 56.75 | 54.75 | 52.75 | 49.75 | 44.75 | 39.75 | 29.75 |
|---|---|---|---|---|---|---|---|---|---|
| t, min u <= 0 (h) | 1.68 | 1.84 | 1.92 | 1.94 | 2.07 | 2.24 | 2.51 | 2.78 | 4.00 |
| t, N^2 < 0 (h) | 2.49 | 2.72 | 3.21 | 3.67 | 3.97 | 4.32 | 4.41 | 3.63 | 4.12 |

The sequence is monotonic in height above 30 km. The N^2 < 0 signal follows the
reversal at each level by 0.6 to 1.9 h, which is consistent with the reversal
criterion being the leading indicator and the finite-amplitude overturning
developing afterwards.

The 0.75 km level does not belong to this sequence. It reverses at 3.44 h and
registers N^2 < 0 at 2.17 h, out of order with respect to the levels above it,
which supports treating it as the separate near-surface feature already
excluded from the descent metric.

## 5. The deceleration equals the resolved wave-momentum flux divergence

The sponge coefficient is identically zero below 60 km (`src/sponge.jl`), so
below that level the only mean-flow forcing available to the solver is the flux
divergence. Comparing d<u>_x/dt against -(1/rho) d<rho u'w'>_x/dz, both in
m/s per hour, at the level of peak deceleration:

| t (h) | z of peak (km) | d<u>/dt | wave drag | ratio |
|---|---|---|---|---|
| 2 | 59.75 | -4.66 | -4.75 | 0.98 |
| 3 | 56.75 | -9.13 | -7.81 | 1.17 |
| 4 | 54.75 | -15.30 | -13.02 | 1.18 |
| 5 | 49.75 | -12.62 | -12.23 | 1.03 |

Closure is within 2 to 18% from 2 h to 5 h. Two systematic departures are
present. At mid-levels the residual is positive (+0.4 to +1.7 m/s per hour),
which is the sign expected from the lateral term: Stage 1 uses
`bc_x_type = :outflow`, so the domain is not periodic and undecelerated 40 m/s
air enters through the inflow edge. The x-mean budget is therefore open by
construction. At 6 h the closure fails above 50 km, where the mean flow is
re-accelerating (+0.8 to +3.2 m/s per hour) while the flux divergence remains
negative; that region is discussed in Section 8.

## 6. The descent matches two independent predictions

The base of the decelerated layer, defined as the lowest level above 2 km where
<u>_x has slowed by more than 5 m/s, is compared against (a) the Lindzen
saturation level, the lowest level at which |MF| >= F_sat, computed from the
background and the instantaneous <u>_x, and (b) the kinematic descent speed of a
mean-wind contour, -(d<u>/dt)/(d<u>/dz), which is an identity rather than a
model.

| t (h) | observed base (km) | z where \|MF\| = F_sat (km) | observed dz/dt (km/h) | -(du/dt)/(du/dz) (km/h) |
|---|---|---|---|---|
| 3.0 | 57.25 | 56.75 | -4.71 | -5.04 |
| 3.5 | 55.25 | 54.25 | -3.50 | -3.24 |
| 4.0 | 53.75 | 52.75 | -3.43 | -3.44 |
| 4.5 | 50.75 | 48.75 | -12.19 | -11.62 |
| 5.0 | 41.25 | 43.25 | -16.35 | -13.90 |
| 5.5 | 37.25 | 40.75 | -7.92 | -9.11 |
| 6.0 | 35.25 | 40.75 | -1.97 | -3.84 |

The saturation level tracks the observed base to within 0.5 to 2.0 km from 3 h
to 5 h, using a formula with no free parameter beyond k_eff. The kinematic
speed reproduces the observed rate to within 10 to 30% throughout, including
the acceleration of the descent near 4.5 to 5 h. The agreement supports reading
the descent as the downward migration of a mean-wind contour driven by drag,
rather than as a structure propagating downward.

The two predictions diverge from the observation after 5.5 h, by 5.5 km at 6 h.
That interval is also where the momentum budget stops closing and where the
horizontal resolution becomes marginal (Section 8), so the late part of the
record is the least constrained.

## 7. The descent has a predicted lower bound

Because F_sat scales as U^3 and rho decreases upward, the level at which a given
flux saturates falls as the flux grows. Evaluating F_sat(z) at U = 40 m/s
against the flux actually supplied at 5 km:

| supplied flux | value (mPa) | saturation level (km) |
|---|---|---|
| measured at t = 1 h | 2.60 | 45.2 |
| measured at t = 2 h | 6.22 | 37.8 |
| measured at t = 3 h | 9.84 | 33.8 |
| measured at t = 6 h | 10.87 | 33.2 |
| full steady linear solution | 22.23 | 25.8 |

The low-level flux is close to its plateau by 3 h, so the saturation level
implied by the flux alone is near 33 km. The observed base reaches 35.25 km at
6 h and its descent rate falls from -16.4 km/h at 5 h to -2.0 km/h at 6 h. Both
are consistent with the layer approaching a level near 33 km rather than
continuing toward the surface. The ascending stationary-wave front is near
37 km at 6 h, so the descending saturation level and the ascending flux front
meet at about that time, which offers an alternative and not separable
explanation for the same slowdown.

This yields a testable prediction. Under the assumptions above, extending the
integration should show the layer base stalling in the 26 to 34 km range rather
than descending further, and the low-level flux approaching the 22.2 mPa linear
value as the stationary wave completes its spin-up.

## 8. Configuration influences that are not part of the mechanism

**Sponge reference state.** `src/run.jl:57` sets `Q_bar = copy(Q)`, the initial
condition, which carries the uniform 40 m/s wind. The sponge therefore relaxes
the zonal wind toward 40 m/s at all times. Measured <u>_x at 6 h is 40.02 m/s at
69.75 km and 40.01 m/s at 79.75 km, confirming the pinning, against 14.4 m/s at
56.75 km. A shear layer of about 2 m/s per km is maintained between the
decelerated interior and the pinned layer. The mean Richardson number there is
8 to 15, above the 1/4 threshold, so the mean shear alone does not satisfy the
Miles-Howard condition for instability.

**Instability excess near the sponge base.** Normalizing the unstable fraction
by a slope proxy S = w_rms N / (k U^2), which is the area-consistent analogue of
|d(zeta)/dz|:

| t (h) | 45-57 km | 57.5-60 km | ratio |
|---|---|---|---|
| 4 | 0.0157 | 0.0484 | 3.1 |
| 5 | 0.0067 | 0.0252 | 3.8 |
| 6 | 0.0054 | 0.0252 | 4.7 |

The 57.5-60 km band carries 3 to 5 times more static instability per unit
resolved wave slope than the rest of the column. The 0.604 peak at 58.25 km at
6 h therefore should not be read as a measurement of wave breaking, and the
previously quoted 3.9x instability contrast between 58.75 and 49.75 km is at
least partly attributable to the band's proximity to the sponge rather than to
the wave field.

**Sponge strength.** The run uses `K0_SPONGE = 0.05`
(`scripts/sweep_arsia_mons_stage1.jl:65`), ten times the value documented in
`src/sponge.jl` as the Srivastava optimum for a reflection coefficient near
0.55%. This bears on the outstanding r = 3.74% item in the tracker.

**Onset altitude is bounded by the domain.** First overturning occurs at
59.75 km, the highest sponge-free level. Since the amplitude increases
monotonically with height, a taller free column would place the first
overturning higher and earlier, and the descent would begin from there. The
asymptotic saturation level of Section 7 is set by the low-level flux and the
background column and is not domain-bounded, but the trajectory toward it is.

**Deposition partition.** The fraction of the flux entering at 25 km that is
absorbed in each band:

| t (h) | 25-40 km | 40-45 km | 45-53 km | 53-56 km | 56-60 km |
|---|---|---|---|---|---|
| 3 | 61.4% | 11.6% | 9.1% | 2.5% | 4.8% |
| 4 | 62.5% | 12.1% | 9.5% | 5.1% | 3.5% |
| 5 | 57.7% | 10.0% | 22.3% | 2.4% | 0.2% |
| 6 | 67.1% | 10.6% | 13.7% | -0.2% | 1.3% |

Between 58 and 67% of the deposition occurs between 25 and 40 km at all times,
and the sponge-adjacent 56-60 km band accounts for 0.2 to 4.8%. The
deceleration profile appears top-heavy because the drag carries a 1/rho factor,
not because the momentum is deposited there. This limits how much of the
column-integrated response the sponge-adjacent band can be responsible for.

**Resolution.** The vertical wavelength 2 pi U / N remains resolved by 19 to 54
points at dz = 500 m throughout. The horizontal scale required to carry the
measured flux at saturation amplitude, lambda_x = pi rho U^3 / (|MF| N), falls
from about 980 km at 55 km at 3 h to 45 km at 55 km at 6 h, which is 11 points
at dx = 4 km. The state after about 5 h is therefore marginally resolved
horizontally, consistent with the existing caveat that the overturning eddies
are not resolved.

**Omitted physics.** Radiative cooling is disabled and the kinematic viscosity
is capped at 1e4 m^2/s, so no thermal damping limits the amplitude growth. The
outstanding Eckermann et al. (2011) CO2 damping item would be expected to
reduce the flux reaching the AMEC band and therefore to lower the saturation
level reached in a given time. The model is also two-dimensional and
non-rotating.

## 9. Interpretation

The observed pattern is reproducible from a single criterion applied to the
run's own background. The wave amplitude grows at the conservative
(rho N)^-1/2 rate, first meets u_hat >= U at the top of the sponge-free column,
and the level at which that criterion is met then descends as the flux supplied
from below grows and as the deposited momentum reduces U, which lowers
F_sat ~ U^3 further. Both the position of that level and the rate at which it
moves are reproduced by independent analytical expressions, and the associated
deceleration matches the resolved flux divergence to within 2 to 18% over
2 to 5 h. Under the assumptions of the model, the results are consistent with
wave saturation and its mean-flow feedback rather than with a numerical
instability.

Two qualifications limit how far that reading extends. The altitude and time at
which the sequence begins are set by the top of the sponge-free column and are
not independent of the domain. The instability concentration in the uppermost
2.5 km of that column exceeds what the resolved wave amplitude accounts for by a
factor of 3 to 5 and is not usable as a breaking measurement. Neither
qualification affects the bulk deposition between 25 and 45 km or the
descent of the layer through the AMEC band, both of which close against the
momentum budget.

## 10. Proposed checks

Ordered by cost.

1. Sponge sensitivity. Re-run 3 h at `K0_SPONGE = 0.005`, the documented
   optimum, and separately with the sponge base at 70 km on a 100 km column.
   Test whether the 57.5-60 km instability excess and the onset height and time
   change. Approximately 1200 s each based on the 2434 s cost of the 6 h run.
2. Extend to 9 or 12 h to test the Section 7 prediction that the layer base
   stalls near 26 to 34 km. Completed 2026-09-01 at 12 h; see Sections 11 to 16.
   The cost estimate given here was low by roughly a factor of three.
3. Horizontal resolution. Re-run at dx = 2 km for 6 h to test whether the
   post-5 h state, which is marginally resolved, is converged. Cost scales
   approximately as dx^-2 in the CFL-limited step count.
4. Domain top. Raise to 120 km using `mars_arsia_mons_my31_cap156km.csv` with a
   raised viscosity cap, to test whether the first overturning altitude is
   domain-bounded as Section 8 argues. This is the existing thermospheric
   extension item and carries the known blowup risk.
5. Repeat the analysis on the Stage 2 configuration, which is periodic in x and
   therefore closes the x-mean momentum budget without a lateral term.

## 11. The 12 h extension and its cost

The configuration of Section 1 was integrated to 12 h with no change other than
`t_end`. The extension mode of `scripts/sweep_arsia_mons_stage1.jl` now derives
its output directory from the requested horizon, so the 6 h record is retained
alongside the new one.

```
julia --project=. scripts/sweep_arsia_mons_stage1.jl ext u=40 h=3 hours=12
```

| quantity | 6 h run | 12 h run |
|---|---|---|
| wall-clock | 2434.1 s | 6253.5 s |
| steps | 32253 | 72922 |
| snapshots | 323 | 730 |
| output size | 1.1 GB | 2.5 GB |
| final t | 21600.0 s | 43200.0 s |
| converged | true | true |

The timestep held between 0.43 and 0.55 s through the final hours and the
`dt < 1e-4` abort guard at `src/run.jl:126` did not trigger, so no critical
level with zero mean wind formed in the resolved flow. Mean throughput was
11.7 steps per second against 21.3 steps per second measured over the first
183 s. Since the per-step cost does not depend on dt, the difference is
attributable to snapshot output rather than to the solver, which suggests that
the output cadence rather than the timestep sets throughput for long runs on
this volume.

The deep-dive extraction over 730 snapshots took 382.9 s, against 193.3 s for
323 snapshots, which is consistent with linear scaling in snapshot count. The
three figures and the mechanism report regenerate in about 6 s combined.

The solver contains no `@threads`, `@spawn` or `Threads.` call in `src/`, so a
single-mode run is serial and deterministic. The 12 h series is bit-identical to
the 6 h series over the whole 5.991 h overlap, on `ubar_zt`, `mf_zt`,
`wrms_hill_zt`, `umin_hill_zt` and `unstable_frac_zt` and on the time axis, with
a maximum absolute difference of 0.000e+00 in every field. The 6 h results
reported in Sections 2 to 9 therefore carry over without restatement.

## 12. The base stalls in the predicted band between 6.5 and 10.5 h

Section 7 predicted that the base of the decelerated layer, measured at a 5 m/s
threshold with z < 2 km excluded, would stop descending in a band near 26 to
34 km rather than continuing downward. Over 6.5 to 10.5 h the observed base
remains within that band and its descent rate falls to near zero, reversing sign
in four of the nine samples.

| t [h] | base [km] | dz/dt [km/h] | z(\|MF\| = F_sat) [km] |
|---|---|---|---|
| 6.0 | 35.25 | -2.38 | 40.25 |
| 6.5 | 33.25 | -4.35 | 39.75 |
| 7.0 | 31.25 | -2.24 | 37.25 |
| 7.5 | 29.25 | -3.15 | 36.25 |
| 8.0 | 28.25 | +0.40 | 38.25 |
| 8.5 | 28.25 | +0.72 | 37.75 |
| 9.0 | 28.25 | +1.51 | 36.25 |
| 9.5 | 28.25 | -1.85 | 37.75 |
| 10.0 | 27.25 | +1.17 | 38.25 |
| 10.5 | 26.25 | -3.02 | 38.75 |
| 11.0 | 24.75 | -14.36 | 44.75 |
| 11.5 | 15.75 | -0.67 | 44.25 |
| 12.0 | 11.25 | -46.25 | 40.25 |

The plateau at 26 to 33 km spanning about four hours is consistent with the
Section 7 prediction. The record does not, however, show that plateau as a
terminal state. After about 10.5 h the base falls to 11.25 km, while the
saturation-level diagnostic moves in the opposite direction, from 36.25 km at
9 h to 44.75 km at 11 h. The two quantities tracked each other to within 0.5 to
6 km from 3 to 7.5 h and decouple thereafter. Section 15 gives grounds for
treating the post-7 h portion of this table as unreliable rather than as a
measurement of the predicted behaviour.

## 13. After about 5 h the layer deepens rather than descends

The altitude of peak deceleration, plotted in panel (a) of
`deepdive_hovmoller.png`, remains between 50 and 55 km from 5 h to 12 h while
the base of the 5 m/s deficit contour moves downward. Over the same interval the
levels that had decelerated most strongly recover.

| z [km] | 3 h | 6 h | 8 h | 10 h | 12 h |
|---|---|---|---|---|---|
| 54.75 | 38.50 | 14.54 | 21.37 | 24.01 | 26.00 |
| 49.75 | 40.97 | 21.18 | 21.44 | 23.34 | 22.79 |
| 39.75 | 40.83 | 30.70 | 28.73 | 32.51 | 36.44 |
| 29.75 | 40.81 | 37.11 | 34.34 | 33.15 | 31.77 |
| 24.75 | 40.72 | 38.53 | 36.53 | 36.23 | 33.42 |
| 9.75 | 40.78 | 41.29 | 41.25 | 39.09 | 35.75 |

Domain-mean zonal wind in m/s. The 39.75 and 54.75 km levels recover by 5.7 and
11.5 m/s respectively between 6 h and 12 h, while 9.75 and 24.75 km lose 5.5 and
5.1 m/s. Both bands lie well below the sponge base at 60 km.

The 6 h record was read in Section 6 as a descending layer. The 12 h record is
more consistent with a layer whose upper boundary is approximately stationary
near 50 to 55 km and which extends downward with time, with partial recovery
above. The overturning criterion is met progressively lower over the same
interval, at 29.75 km at 4.00 h, 19.75 km at 6.23 h and 9.75 km at 9.91 h, and
the hill-window unstable fraction becomes distributed through the column rather
than concentrated in a band, reaching 0.204 at 5 km, 0.132 at 10 km, 0.156 at
15 km and 0.124 at 25 km by 12 h.

## 14. The supplied flux does not approach the steady linear value

The Section 7 lower bound of 25.8 km assumed the low-level flux would grow
toward the steady linear hydrostatic value of 22.23 mPa. That assumption is not
supported by the 12 h record. The flux at 25 km plateaus between 5.3 and
6.4 mPa from 5 h onward, which is 24 to 29% of the steady linear value, and the
flux at 5 km peaks near 11.9 mPa at 4 h and 8 h before falling to 3.66 mPa at
12 h. The corresponding equilibrium crossing height stays between 32.2 and
36.2 km from 4 h to 11 h rather than approaching 25.8 km.

A linear control corner run on 2026-09-01 qualifies the reading of these ratios
(`docs/linear_flux_constancy.md`). At h = 0.1 km with the same grid, background
column and horizon but a conservative wave and a periodic x boundary, the flux
after 12 h reaches only 0.77 of its own steady linear value at 9.75 km and 0.42
at 54.75 km, and is still rising at 20% per hour at 54.75 km. A shortfall against
the steady linear value at 12 h is therefore not by itself evidence of
saturation. A follow-up analysis on the same record identifies the mechanism:
arrival time for a stationary hydrostatic component scales as 1/k, and at 55 km
after 12 h the two longest modes of the ridge spectrum, together 45% of the flux,
have not arrived, reaching that height only at 13.4 h and 26.8 h. Benchmarking the
control against the arrived part of the spectrum rather than against a fully
established wave reduces its fitted drift over 10-55 km from -1.24 to +0.04% per
km. Flux ratios against the steady linear value at any finite time should be read
with the corresponding arrived fraction, which over 10-55 km in the final hour is
1.00 below about 20 km, 0.82 from 25 to 45 km and 0.55 from 50 to 55 km. The part of the shortfall here that the control does not reproduce is
quantified in Section 16.

Two factors plausibly contribute. The wind at the levels that supply the flux
falls over the record, from 40.7 m/s at 24.75 km at 3 h to 33.4 m/s at 12 h, so
F_sat at those levels declines by a factor of about (33.4/40.7)^3 = 0.55 and the
column absorbs a larger share of the flux at low altitude. Separately, the
deposition partition remains bottom-heavy at all times, with 49 to 69% of the
flux entering at 25 km absorbed between 25 and 40 km and only 0.2 to 5.9%
absorbed in the sponge-adjacent 56 to 60 km band, so the sponge remains
unimplicated in the bulk response.

Two of the Section 8 qualifications weaken over the extended record. The
instability excess per unit resolved wave slope in the uppermost 2.5 km of the
sponge-free column falls from a factor of 4.6 at 6 h to 1.06 at 12 h, as the
interior instability grows to match it, which suggests the excess was a feature
of the startup phase rather than a persistent property. The horizontal
resolution of the flux-carrying scale at 55 km improves from 12.8 points per
wavelength at 6 h to 98.6 at 12 h as the wind there recovers, so the late-time
state is better resolved than the state at 6 h.

## 15. The domain-mean budget stops closing after about 7 h

The x-mean momentum budget of Section 5 closes to within 20% from 2 h to 6 h and
then degrades.

| t [h] | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ratio d⟨u⟩/dt to drag | 1.07 | 1.20 | 1.18 | 1.20 | 1.00 | 0.56 | 0.57 | 0.25 | 0.38 | 0.51 | 0.70 |

Evaluated at the level of peak deceleration. From 7 h onward the resolved flux
divergence would produce a deceleration two to four times larger than the one
observed, so a term of the opposite sign is supplying momentum at those levels.

Two properties of the Stage 1 configuration bear on this and were not
quantified in Section 8.

The x boundary condition is `:outflow` and therefore not periodic, so the
x-mean budget is open by construction. The advective transit time from the
inflow edge to the ridge is LX / (2U) = 6.94 h and across the full domain
13.89 h. The budget begins to fail between the 6 h and 7 h samples, which
coincides with the first arrival at the ridge of air that entered through the
inflow boundary after the deceleration was established. The coincidence does not
by itself establish causation, but the sign of the residual is the sign expected
from lateral influx of undecelerated air.

The domain is 2000 km wide against a flux-weighted horizontal wavelength of
697 km, so it spans 2.87 wavelengths. Once the wave amplitude approaches U, an
x-mean over that width does not cleanly separate wave from mean flow. Comparing
the mean over a 700 km window centred on the ridge with the mean over the
remainder gives differences that reach 26.2 m/s at 49.75 km and 18.0 m/s at
39.75 km at 12 h, alternating in sign with height on a scale comparable to the
vertical wavelength of about 15 km. At 3 h the same differences are 6.5 and
9.7 m/s. The quantity called the mean wind throughout this report is therefore
progressively less well defined as the run proceeds.

A third indication concerns the local extrema rather than the mean. The minimum
hill-window zonal velocity reaches -113.5 m/s at 11 h, or Mach 0.58 against a
local sound speed of 195.7 m/s, and from 9.54 h onward isolated samples exceed
the local sound speed, with excursions near -250 m/s visible in panel (b) of
`deepdive_breaking.png`. Local flow approaching and briefly exceeding the sound
speed in a 40 m/s background is not a plausible amplitude for a Mars mountain
wave and is more consistent with grid-scale noise produced where the breaking is
unresolved. This is also the likely origin of the timestep reduction from about
0.55 to 0.43 s over the same interval, since the CFL limit scales with
\|u\| + c_s. The domain-mean quantities remain smooth, so this affects the extrema
rather than the budget, but it independently marks the second half of the record
as unsuitable for quantitative claims.

By contrast, the minimum zonal-mean wind below 60 km, shown in panel (d) of the
same figure, falls from 40 m/s to 13.4 m/s at 5.7 h and then recovers to about
22 m/s and holds that value from 8 h to 12 h. It does not approach zero, so the
run does not develop a self-induced critical level in the mean flow. That the
mean-flow deceleration reaches a floor and stabilizes is consistent with a
saturation equilibrium, although at a mean-wind value rather than at the layer
base altitude that Section 7 predicted.

## 16. What the 12 h record supports

The saturation account of Sections 1 to 9 continues to describe the run over the
interval where the diagnostics are well posed. Through about 7 h the base tracks
the Lindzen saturation level to within 0.5 to 6 km, the budget closes to within
20%, and the deposition stays bottom-heavy and away from the sponge. Two of the
Section 8 qualifications weaken rather than strengthen with time.

The Section 7 prediction is partially supported. The base does decelerate and
occupy the predicted 26 to 34 km band for about four hours, from 6.5 to 10.5 h,
with a descent rate that reverses sign in four samples. It does not remain
there, and the premise behind the 25.8 km asymptotic value, that the supplied
flux grows toward the steady linear flux, is not borne out. The stated
prediction should be regarded as describing a transient plateau rather than an
equilibrium.

The clearest result of the extension may be a bound on the configuration rather
than on the physics. Two independent diagnostics, the budget residual and the
ridge-to-far-field wind contrast, both indicate that domain-mean quantities in
this Stage 1 setup become unreliable after roughly 7 h, on a timescale set by
LX / (2U). Conclusions drawn from the 6 h record fall inside that window.

A linear control corner run on 2026-09-01 (`docs/linear_flux_constancy.md`) adds
a height bound to the saturation attribution itself. Differencing this record
against a conservative wave in the same solver, on the same grid and at the same
time, with each profile normalized by its own steady linear value and then by its
own 10 km value, leaves an excess that scatters around unity below 30 km, between
0.71 and 1.37 with no trend, and falls monotonically above 35 km to 0.61 at
39.75 km, 0.52 at 44.75 km, 0.21 at 49.75 km and 0.099 at 54.75 km. The results
are consistent with saturation acting above about 35 km. This excess is
unaffected by the arrival correction noted in Section 14, because the arrived
fraction depends on the launched ridge spectrum, background column, wind and
elapsed time, all of which the two runs share, so it multiplies numerator and
denominator alike. They do not support
attributing the flux decrease below about 30 km to saturation from this record,
since a conservative wave in the same configuration reproduces that part of the
height dependence. The 58 to 67% deposition partition over 25 to 40 km reported
in Section 9 spans the boundary between the two regimes and should be read with
that division in mind.

The same control indicates that a conservative wave has not settled above 30 km
at 12 h in this configuration, so the flux profiles reported here above that
height describe a filling column rather than a steady state.
Conclusions about the post-10.5 h collapse of the base do not, and are better
treated as a diagnostic artifact of a domain-mean taken over 2.87 wavelengths in
a laterally open domain than as a physical result.

Under the assumptions of the model, the results are consistent with wave
saturation and its mean-flow feedback over 2 to 7 h, with the layer thereafter
deepening downward from a roughly stationary top near 50 to 55 km rather than
descending as a whole. Testing whether that deepening is physical requires a
configuration whose x-mean budget closes. The Stage 2 configuration is periodic
in x and is the natural next test, now item 1 rather than item 5 in the
Section 10 list.
