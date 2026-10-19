# Reflection from the damping layer and convergence in the linear control corner

Internal validation report, started 2026-09-22. Companion to
`docs/linear_flux_constancy.md` (Sections 15 to 19) and to
`scripts/math_mountain_waves.ipynb`.

Question addressed: whether the downward-propagating component measured in the Arsia
Mons linear control corner is returned by the model-top damping layer, how that return
depends on the layer's strength, thickness and base height, and whether the record is
long enough for the downward share to have settled.

The linear theory in the notebook fixes what the downward share should be for a clean
top: a partition floor of `2.4e-04` from the dropped `1/4H_s^2` term, and a partial
reflection off the represented `N(z)` of order `5e-04` to `2e-03`. The start-up
transient is the only candidate that carries a clock. A plateau above those floors is
therefore a property of the solver, and the damping layer is the first candidate for
it. This document records what has been tested on that question, what had not been,
and the runs made to close it.

## 1. Inventory of earlier tests

### In the Arsia linear control corner

All of these use `h = 0.1` or `0.03` km, uniform `U = 40` m/s, periodic x,
`LX = 2000` km, `dx = 4` km, `dz = 500` m, a `sin^2` Rayleigh layer over 60 to 90 km
with `K0 = 0.05` s^-1, and the `:eq19` outflow top.

| test | record | main result | reference |
|---|---|---|---|
| flux constancy, `sigma = 100` km | 12 h, `h = 0.1` km | drift -1.24 %/km, consistent with incomplete spectral fill-in | `linear_flux_constancy.md` §8-12 |
| ridge-width pair and extension | 12 and 15 h, `sigma = 50` km | drift -0.033 %/km at 15 h, inside the 0.1 %/km threshold; both time-dependence models excluded out of sample | §13-14 |
| up and down split | 15 h | band `R = 0.049`, rising; returning wave is 17 % of the level loss; incident flux decays 1.1 % per hour under a constant source | §15 |
| variance cross term | 15 h | peak scale 15.0 km against `lambda_z/2 = 15.3` km; built from the split's own A and B, so not independent of it | §15 |
| per-mode timing | 15 h, both widths | downward share present in modes whose sponge round trip has not elapsed; scale ordering runs against a damping-layer origin below 30 km | §16 |
| ground-frame stationarity | 15 h | measured rotation rates do not bring any unreturned mode inside the record; about a third of `R` is a fluctuating part the split divides evenly | §17 |
| amplitude control | `h = 0.03` km, 15 h | unreturned steady part scales as `h^3.3` to `h^3.6` (nonlinear); returned class at 9.75 km scales as `h^2.14` (linear); static instability confined to 58.75-61.75 km, on the sponge base | §18 |
| linear theory | none | closed forms, the `2.4e-04` floor, `N(z)` reflection by shooting | notebook §1-10 |

### Elsewhere in the repository

None of these ran in the Arsia configuration.

| test | configuration | varied | result | reference |
|---|---|---|---|---|
| Srivastava sponge validation | Earth, isothermal 239 K, wave packet, 20 km sponge | `K0` from 0 to 0.5 | `r` about 20 % for `K0 >= 0.005`, set by thickness | `output/sponge_validation/` |
| thickness and parameter sweep | Earth 239 and 289 K | 20/40/60 km, `K0` 0 and 0.005, viscosity | `r` = 19.7 / 7.4 / 3.7 % at 239 K; the paper's 0.55 % not reached | `TODO_tracker.md:13-28`, `suggestions_compilation.md` |
| well-balanced against standard scheme | Earth | scheme, `K0` | low paper reflections attributed to numerical diffusion of the standard scheme | `suggestions_compilation.md:214-345` |
| R3 thickness sweep | Mars isothermal 210 K, acoustic, `lambda_z = 92` km | 30 to 120 km | `r` levels at about 1.3 %, attributed to the eq19 top | `sponge_sweep_R3.md` |
| characteristic and Sommerfeld tops | Earth | top BC | both fail within 14 to 400 steps | `TODO_tracker.md:531-533` |
| Klemp-Durran radiation top | Tier 2 R2/R3 Mars MCS | top BC, domain height | clean envelope on the 120 km domain with steady `w` 66 % larger; clamp active on 200 km; R3 unstable | `klemp_durran_bc_results.md` |
| cavity interference fits | Tier 2 R2/R3 | time envelope at 100 km | no clean window; decay dominates the interference sum | `cavity_interference_analysis.md` |
| grid, boundary and low-frequency sweeps | Earth and idealized Mars | `dx = dz`, `bc_type`, thickness, `K0` | scripts exist, results not written up | `sweep_dx_reflection_scaling.jl`, `sweep_wb_boundary_sponge.jl`, `test_ultra_low_freq_sponge.jl` |

`ext6h_instability_mechanism.md` proposed re-running at `K0 = 0.005` and with the sponge
base at 70 km. Neither had been run.

### What had not been tested

1. No variation of the damping layer in the Arsia corner. Strength, thickness and base
   height were all fixed, and the driver did not expose them.
2. No record longer than 15 h, while mode `n = 1` reaches 55 km only at 26.8 h.
3. No estimate of the downward share that is independent of `updown_split()`.
4. No theoretical value for how much this layer should reflect.
5. No `dz` or `dx` refinement in the Arsia configuration.

Items 1, 3 and 4 are addressed below and item 2 in part. Item 5 is deferred (Section 6).

## 2. Changes made for these tests

`scripts/sweep_arsia_mons_stage2.jl` accepts `k0=` (1/s), `lz=` (km) and `spbase=`
(km) in `linear` mode, with defaults 0.05, 90 and 60. Non-default values are appended
to the output directory name, so no existing record can be overwritten. The layer
always extends from its base to the model top. `sweep_results.csv` is appended to
rather than rewritten. A new file carries the sponge columns. A 0.2 h run at the
defaults reproduces all ten snapshots it shares with
`arsia_mons_linear15h_LX2000_sig50/run_U40.0_h0.03`, bit for bit in every stored field.

`scripts/analyze_linear_flux_constancy.py` accepts `--sponge-base-km`, which sets the
round-trip clock of the per-mode classification. It also carries the two functions
described next.

## 3. How much the layer should reflect

`rayleigh_sponge_reflection()` integrates the anelastic structure equation for
`Phi = sqrt(rho_0) w'` downward from the model top,

```
Phi'' + [k^2 (N^2 / (Omega + iK)^2 - 1) - s^2 + s'] Phi = 0,
Omega = -U k,   s = -(1/2) d ln(rho_0)/dz,
K(z) = K0 sin^2(pi/2 (z - z_base)/(z_top - z_base)),
```

on the MCD column and the solver's 500 m grid. It then splits the field at each level
with the exact local wavenumber. The damping enters as `Omega -> Omega + iK` because
the solver relaxes all four conserved variables toward the initial state at the same
rate. This happens in both the predictor and the corrector of the z sweep
(`src/euler.jl:149-189`). The top is either a rigid lid or a radiation condition,
which bracket the eq19 outflow. For every configuration below, the two agree to within
5 per cent, because the wave is returned or absorbed well below the top.

The construction reproduces its limits. With `K0 = 0` and a lid it returns 1 at every
level. With `K0 = 1000` s^-1 and a lid it returns 0.99999999. With `K0 = 0` and a
radiation top at 90 km it returns `2.2e-04`, `1.3e-03` and `5.2e-04` at 9.75, 29.75
and 54.75 km for `n = 3`, the `N(z)` reflection alone. The notebook puts that same
quantity at `5e-04` to `2e-03`. With `K0 = 0` and a radiation top at 120 km it
returns `6.9e-03` at 9.75 km for `n = 3`. That is a separate result: the MCD column above 90 km,
where the temperature rises again, reflects on its own. A higher model top therefore
adds a reflector that a thicker layer has to absorb before reaching it.

Predicted band-mean reflected share over 10 to 55 km, by mode. Fixed at 18:35 on
2026-09-22, before any of the runs below had produced a result.

| run | layer | `K0` [1/s] | n=1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| S0 | 60-90 km | 0.05 | 0.605 | 0.495 | 0.424 | 0.373 | 0.333 | 0.301 | 0.275 | 0.253 | 0.234 | 0.217 | 0.203 | 0.191 |
| S2 | 60-90 km | 0.0158 | 0.414 | 0.290 | 0.220 | 0.175 | 0.142 | 0.119 | 0.100 | 0.086 | 0.074 | 0.065 | 0.058 | 0.051 |
| S1 | 60-90 km | 0.005 | 0.211 | 0.110 | 0.066 | 0.043 | 0.030 | 0.021 | 0.015 | 0.012 | 0.009 | 0.007 | 0.006 | 0.005 |
| T1 | 60-120 km | 0.005 | 0.045 | 0.012 | 0.005 | 0.003 | 0.002 | 0.002 | 0.002 | 0.002 | 0.002 | 0.002 | 0.002 | 0.002 |
| T2 | 70-100 km | 0.05 | 0.639 | 0.539 | 0.474 | 0.427 | 0.389 | 0.359 | 0.334 | 0.313 | 0.295 | 0.279 | 0.265 | 0.254 |

On this theory the present layer is close to a lid at its base. `K0/omega_hat` is 33
to 398 across the flux-carrying modes, so the damping rate overtakes the intrinsic
frequency within 1.4 to 3.3 km of the base (`linear_flux_constancy.md` §15). That
is 0.04 to 0.11 of a vertical wavelength. The theory also orders the return with
`lambda_x`: long waves reflect more, since `K0/omega_hat` is larger for them.

## 4. An estimate of the reflected share that does not use the split

`updown_split()` reads the direction of each mode from the phase relation between
`Psi` and its derivative. `ripple_reflection()` uses `|Psi_k|^2` alone. For
`Psi = A e^{i theta} + B e^{-i theta}` with WKB amplitudes proportional to `m^-1/2`,

```
m |Psi|^2 = a^2 + b^2 + 2ab cos(2 theta + phi),   theta = int m dz,
```

so the density-scaled, WKB-flattened variance carries a ripple at `lambda_z/2`. The
ripple's relative depth `2q/(1 + q^2)` fixes `R = q^2`. This is the standing-wave
ratio of a transmission line, set beside a directional coupler. A linear trend in `z`
is fitted alongside, since a mode that is still filling carries one. The
exact local `m`, with the `-s^2 + s'` terms, is used for the phase. The fit is repeated
over `m` scaled from 0.9 to 1.1, and the scale at which the ripple is deepest is
reported, as a check that the ripple has the wavelength a reflection gives.

On the synthetic WKB wave, imposed reflections of 0.00, 0.05 and 0.10 are returned as
0.0000, 0.0500 and 0.1000 in every mode tested, and the exact anelastic solution is
returned as 0.

### Applied to the existing 15 h records

Band-mean over 10 to 55 km at 14.74 h, modes whose round trip from a 60 km base has
elapsed at every level of the band:

| n | `lambda_x` [km] | `h = 30` m split | `h = 30` m ripple | `h = 100` m split | `h = 100` m ripple | theory |
|---|---|---|---|---|---|---|
| 6 | 333 | 0.0304 | 0.0292 | 0.0190 | 0.0178 | 0.301 |
| 7 | 286 | 0.0276 | 0.0254 | 0.0297 | 0.0223 | 0.275 |
| 8 | 250 | 0.0290 | 0.0170 | 0.0302 | 0.0208 | 0.253 |
| 9 | 222 | 0.0365 | 0.0346 | 0.0304 | 0.0276 | 0.234 |
| 10 | 200 | 0.0426 | 0.0388 | 0.0339 | 0.0286 | 0.217 |
| 11 | 182 | 0.0686 | 0.0533 | 0.0675 | 0.0584 | 0.203 |
| 12 | 167 | 0.0559 | 0.0532 | 0.0531 | 0.0604 | 0.191 |

The best-fit `m` scale for these modes lies between 0.955 and 1.02.

Three observations follow.

The two estimators agree. Across the fourteen entries the ripple lies between 41 per
cent below and 14 per cent above the split, with a median ratio of 0.91, and its
wavelength is that of a reflected wave. The downward
component in these modes is therefore present in the variance of the field and is not
an artefact of the partition. This settles the tension noted in the inventory, where
the variance cross term and the split appeared to disagree: that cross term was
built from the split's own `A` and `B`.

The measured share is 3 to 17 times smaller than the linear Rayleigh-layer value. It
also rises with `n` (0.02 to 0.07 from `n = 6` to 12), where the theory falls (0.30 to
0.19). A layer that reflected as the theory describes, scaled by a mode-independent
factor, would give the opposite ordering.

The share is close to amplitude independent between `h = 30` and 100 m in these modes,
which agrees with the `h^2.14` exponent Section 18 of `linear_flux_constancy.md`
reports for the returned class at 9.75 km.

These results are consistent with either of two readings, which the runs of Section 5
separate. Either the solver's layer reflects much less than the continuous linear
theory predicts, and in a different mode ordering. Or the downward component in these
modes is not mainly produced by the layer.

## 5. Runs and the predictions fixed for them

All at `h = 30` m, `sigma = 50` km, `U = 40` m/s, periodic `LX = 2000` km,
`dx = 4` km, `dz = 500` m, 15 h, with only the layer changed. Each was launched as
its own process. The solver is serial, and the runs shared the 24-core machine, so
wall-clock values carry some contention.

| id | layer | model top | command suffix | output |
|---|---|---|---|---|
| S0 | 60-90 km, `K0 = 0.05` | 90 km | existing | `arsia_mons_linear15h_LX2000_sig50/run_U40.0_h0.03` |
| S1 | 60-90 km, `K0 = 0.005` | 90 km | `k0=0.005` | `..._sig50_K0_0.005/` |
| S2 | 60-90 km, `K0 = 0.0158` | 90 km | `k0=0.0158` | `..._sig50_K0_0.0158/` |
| T1 | 60-120 km, `K0 = 0.005` | 120 km | `k0=0.005 lz=120` | `..._sig50_K0_0.005_LZ120/` |
| T2 | 70-100 km, `K0 = 0.05` | 100 km | `lz=100 spbase=70` | `..._sig50_LZ100_spb70/` |

Base command:
`julia --project=. scripts/sweep_arsia_mons_stage2.jl linear sigma=50 hours=15 h=0.03 <suffix>`.

T2 translates the present layer upward by 10 km with its strength and thickness
unchanged, so that the theory predicts almost the same reflected share (Section 3)
but a later return. At 9.75 km the round trip grows from 13.35, 10.68 and 8.90 h to
16.09, 12.87 and 10.73 h for `n = 4`, 5 and 6. The 70 km base was chosen over 75 km on
the linearity margin: at `h = 30` m the steady wave reaches `u'/U` of about 0.38 at
70 km, against 0.50 at 75 km.

Readings fixed before the results:

1. If the layer is the source of the returned-class share, `R` in modes 6 to 12
   falls from S0 to S2 to S1 to T1, by roughly the theoretical ratios. S1/S0 would
   then be 0.07 at `n = 6` and 0.03 at `n = 12`, and S2/S0 0.39 and 0.27. The
   unreturned modes stay at or below the closure bound `eps^2/4`. In T2, `n = 4` and
   5 at 9.75 km move from the returned to the unreturned class.
2. If the layer is not the source, `R` in modes 6 to 12 does not change across S0,
   S1, S2 and T1 beyond the mode-to-mode scatter. Its origin is then below the layer,
   and the linear theory of Section 3 overstates what this layer returns.
3. A partial response, where `R` falls but by less than the theoretical ratios,
   indicates a mixture. The level that `R` approaches as `K0` falls estimates the
   part that does not come from the layer.
4. Split and ripple agree within about 20 per cent in every run. A run in which they
   part would indicate a change in the character of the downward component.

A convergence run, L1, follows these. It is a longer record on the configuration that
returns least, chosen on the measured result rather than on the theory. Its purpose is
to follow the downward share by one-hour bins until the slowest modes have arrived and
settled. Its criteria are: bin slopes within one standard error of zero over the last
6 h; the band drift under 0.1 %/km with peak-to-peak under 5 %; and whether the 1.1 %
per hour decay of the incident flux (`linear_flux_constancy.md` §19 item 7) stops.

## 6. Deferred

- `h = 10` m amplitude point (`linear_flux_constancy.md` §19 item 2), which constrains
  the nonlinear exponent of the unreturned class.
- `dz = 250` m (about 4 times the cost) and `dx = 2` km (about 2 times).
- Klemp-Durran radiation top, which has failed on the acoustic R3 case and was left
  out of scope here.
- Settling the `N_0` reference of the amplitude ratio (`linear_flux_constancy.md` §19
  item 6).

## 7. Results

### The runs

| id | steps | wall clock [s] | `min u/U` below 60 km | unstable points, last hour |
|---|---|---|---|---|
| S0 | 73055 | 3684.6 (alone) | 0.762 | 46 % of the row at 59.75 km, none elsewhere |
| S2 | 73054 | 10446.5 | 0.877 | none from 50 to 90 km |
| S1 | 73050 | 10559.9 | 0.883 | none from 50 to 90 km |
| T1 | 73465 | 14416.8 | 0.883 | none |
| T2 | 73054 | 11928.1 | 0.689 | 25 % of the row at 69.75 km, none elsewhere |

All converged to `t = 54000` s. The four new runs were integrated at the same time as
separate processes, all writing snapshots to the network-mounted drive. From about
20 minutes in, three of the four sat at 0 per cent CPU with the machine 85 per cent
idle, and the snapshot rate fell from about 13 to about 6 per minute per run. Their
wall-clock figures are therefore set mainly by I/O and are not a compute cost. They
should not be compared with the 3685 s of S0, which ran alone. The analysis chain took
33 to 53 s for extraction and about 25 s for the flux analysis per run. Commands and
logs are in `output/reflection_tests/logs/`, and the comparison figures and
`compare_sponge_reflection.json` are in `output/reflection_tests/`
(`python3 scripts/compare_sponge_reflection.py`). T2 was analysed with
`--sponge-base-km 70`, and each run with its own `--sponge-k0`.

### The returned share does not follow the linear prediction

Final-hour band-mean over 10 to 55 km, averaged over modes 6 to 12, whose round trip
from the layer base has elapsed at every level of the band.

| run | layer | split | ripple | theory | split / theory | split / S0 | theory / S0 |
|---|---|---|---|---|---|---|---|
| S0 | 60-90 km, `K0 = 0.05` | 0.0415 | 0.0359 | 0.2390 | 0.17 | 1.00 | 1.00 |
| S2 | 60-90 km, `K0 = 0.0158` | 0.0748 | 0.0682 | 0.0790 | 0.95 | 1.80 | 0.33 |
| S1 | 60-90 km, `K0 = 0.005` | 0.0583 | 0.0504 | 0.0107 | 5.4 | 1.40 | 0.045 |
| T1 | 60-120 km, `K0 = 0.005` | 0.0172 | 0.0114 | 0.0018 | 9.7 | 0.42 | 0.007 |
| T2 | 70-100 km, `K0 = 0.05` | 0.1418 | 0.1339 | 0.2998 | 0.47 | 3.42 | 1.25 |

Against the readings fixed in Section 5:

1. The sponge-origin reading is not met in the form stated. The theory predicts the
   returned share to fall by factors of 3, 22 and 130 from S0 to S2, S1 and T1. The
   measured share rises by 1.8 and 1.4 for the two weaker layers and falls by 2.4 for
   the thick one. Moving the layer up by 10 km, for which the theory predicts 1.25,
   raises it by 3.4.
2. The not-the-sponge reading is not met either. The measured share spans a factor of
   8.2 across the five layers (0.017 to 0.142), well outside the mode-to-mode scatter.
   The layer configuration therefore changes what returns, but not in the way the
   continuous linear theory describes.
3. The data fit a partial response in the sense of reading 3. Because the measured
   share is not monotonic in `K0`, the level it approaches as damping is reduced
   cannot be estimated from these runs.
4. Split and ripple agree in four of five runs, the ripple being 0.86 to 0.94 of the
   split. In T1, where the share is smallest, it is 0.66, which is the regime in
   which the fluctuating component of Section 17 is largest relative to the signal.

The theory overstates the return of the two strong layers (S0 by about 6, T2 by about
2) and understates that of the three weak ones (S2 about equal, S1 by 5, T1 by 10).
The measured values span a much narrower range than the theory does. The per-mode
shape is similar across the layers (`reflection_vs_theory.png`). It peaks near
`lambda_x = 333` km, close to the flux-weighted `lambda_eff = 353` km of the ridge,
whereas the theory falls monotonically with `n`.

### Timing

The returned share appears near the round trip the stationary clock predicts. At
9.75 km in S0, `n = 6` (`t_rt = 8.9` h) reads 0.002, 0.015 and 0.042 at 8, 12 and
15 h. In T2 (`t_rt = 10.7` h) it reads 0.001, 0.016 and 0.035. That supports energy
returning from above, on roughly the group-velocity clock. In T2, `n = 4` cannot
have returned from a 70 km base within the record (`t_rt = 16.0` h), yet it reads
0.032 at 15 h against 0.002 in S0. That is a single point in a quantity whose
single-level values vary by factors of several between samples (Section 17), so the
translation test is not resolved by it.

The per-mode returned share was still rising between 12 and 15 h in most returned
modes of S1, S2 and T2. The five configurations are therefore compared while the
return is still building, whereas the theory gives its steady value. That is one
candidate for part of the difference, and it is what L1 addresses.

### What does not change with the layer

The steady part of the unreturned class at 9.75 km, averaged over ten one-hour bins,
is 0.0103, 0.0092, 0.0102, 0.0092 and 0.0108 in S0, S2, S1, T1 and T2, with bin
standard deviations of 0.002 to 0.005. It does not respond to a 10-fold change in
`K0`, a doubling of the layer thickness or a 10 km change in its base. This is
consistent with Section 18 of `linear_flux_constancy.md`, which places that component
below the layer and finds it nonlinear in amplitude.

### What else changes with the layer

The incident flux level changes. The final-hour `|MF_up|` over 10 to 55 km, relative
to the steady linear value, is 0.915, 0.885, 0.791, 0.805 and 1.166 in S0, S2, S1, T1
and T2. The launch wind is the same in every run. The runs are deterministic (L1 and the
smoke test reproduce S0 bit for bit), so these differences are responses to the layer
and not sampling noise. They are ratios between runs on the same reference, so the
factor 2.08 ambiguity in the launch `N_0` (`linear_flux_constancy.md` §9) cancels
from them. Relative to S0 the level changes by -3 % (S2), -12 % (T1), -14 % (S1)
and +27 % (T2). T2 is the one layer that fails the linearity and flatness checks
(`min u/U = 0.689`, an instability at its base, drift -0.364 %/km, peak-to-peak
16 %), and it changes the base and the model top together. Without it the four
linear layers span 0.791 to 0.915, a max/min ratio of 1.16 and a standard deviation
of about 7 % of the mean.

At 15 h the incident flux level below the layer therefore depends on the layer
configuration by -14 % to 0 % relative to S0 across the four linear layers (+27 % in
T2, which fails the linearity thresholds). Whether this is a steady cavity effect or
a transient difference is not settled by 15 h records, in which the returned share
still rises in S1, S2 and T2. (An earlier version of this paragraph said the S0 level
still decays by 1.1 % per hour. That rate was measured on the `h = 100` m record of
`linear_flux_constancy.md` §15, whose linear reference is -25.0 uPa. In S0, at
`h = 30` m, the incident level rose by 0.64 % per hour over 12 to 15 h; see L1 below.)

A cavity between the ground, whose kinematic condition reflects a downward wave back
upward, and the returning layer is one candidate. With a fully reflecting ground the
upward amplitude scales as `1/|1 + r e^{i phi}|^2`. For the measured `R = r^2` of
0.017 to 0.04 that allows a level range of 0.78 to 1.32 up to 0.69 to 1.56, and the
observed spread lies inside it. The mechanism is thus plausible in magnitude but not
demonstrated. It is also the interference channel at the launch level that the flux
bookkeeping of `linear_flux_constancy.md` §15 cannot see. If it holds, the 15 per
cent incident shortfall reported in §9 and §15 depends on the damping layer and is
not a fixed property of the solver. If the least reflecting layers approximate the
answer without a top, the reference level is the cluster near 0.80 (S1, T1), and the
S0 layer used in every earlier record would bias the incident level high by about
14 %.

The static instability found at the sponge base in Section 18 follows the base of the
strong layer. It is at 59.75 km in S0 and at 69.75 km in T2, and absent in the three
weak layers. That supports reading it as produced by a strongly damped layer rather
than by the wave.

Every run meets the flatness thresholds on the net-flux drift except T2: final-hour
drift -0.021, -0.068, +0.024, +0.067 and -0.364 per cent per km, and peak-to-peak
1.9, 3.9, 3.2, 4.0 and 16.3 per cent, in S0, S2, S1, T1 and T2. T2 carries both the
largest return and the largest incident excess.

### What this changes

The continuous linear theory of a `sin^2` Rayleigh layer does not predict the return
of this solver's layer, in magnitude or in mode ordering, under the evaluated
conditions. It should not be used to size the layer. Of the five layers, T1 (60 to
120 km, `K0 = 0.005` s^-1) returns the least, 0.011 to 0.017 in modes 6 to 12 against
0.036 to 0.042 for the layer used in every earlier record, and it produces no
instability at its base. Its cost per step is higher. Under contention T1 took 1.4
times as long as S1, and a single 0.2 h smoke run took 2.5 times as long as the 90 km
column. The cause was identified on 2026-09-23. With `nz = 240` the triangular solves of
the implicit diffusion step (`src/diffusion.jl`, a sparse LU of 500 x 240 unknowns,
solved three times per step) run on multithreaded OpenBLAS. That process then uses
about 13 cores (1343 % CPU) while advancing at about 5 steps per second. With
`OPENBLAS_NUM_THREADS=1` the same column advances at 14.8 steps per second on one
core. The 90 km column does not show the effect (99 % CPU). Runs with `lz > 90` should
therefore be launched with `OPENBLAS_NUM_THREADS=1`. The T1 wall clock above
includes this overhead.

At 15 h the shape bound survives the layer variations. In the four linear layers the
final-hour drift is -0.068 to +0.067 %/km, at most about 3 % accumulated over 45 km.
A percentage change in level and a gradient in %/km are not commensurable, though.
The drift `d ln|MF|/dz` is blind by construction to a loss that is uniform in height
(§9 and §15 of `linear_flux_constancy.md`), so the drift work never bounded the level.
The layer dependence therefore adds a separate row to the error budget. It does not
replace the -0.033 %/km bound of §14 (about 1.5 % accumulated).

The flux error budget of the linear corner, as separate rows that are not summed
(superseded by the post-fix budget of Section 11, since these rows come from `:total`
records):

| component | size | status |
|---|---|---|
| shape (band drift) | 15 h, four linear layers: within +-0.07 %/km, at most about 3 % over 45 km. 36 h, S0 layer (L1): up to 0.15 %/km, peak-to-peak 9 to 19 % | insensitive to the layer at 15 h; not steady in time (L1) |
| level, from the damping layer and top | about 10 to 15 % across linear layers; 27 % in T2 | provisional, 15 h records |
| level, from the linear reference | factor 2.08 from the choice of `N_0` | open (§19 item 6) |
| level, from the transient | `h = 100` m: -1.1 % per hour at 15 h. `h = 30` m (S0, L1): +0.64 % per hour at 12 to 15 h, then about 0.945 +- 0.01 of the linear value over 26 to 36 h | L1, below |

Scope of the level row. It applies to level metrics in configurations where a
substantial share of the wave reaches the layer. In the Stage 2 production runs the
flux is absorbed at the 4.5 km critical level (`TODO_tracker.md`, 2026-09-01 entry),
so the top may matter much less there; this has to be checked against those records
before being generalised. Where the waves do reach the top, the cavity phase depends
on `lambda_z = 2 pi U / N`. A suite with a fixed top would then not cancel the effect
across a sweep in `U`, and it could bias the trend.

Check against the Stage 2 records (`output/arsia_mons_stage2/`, real reversing wind,
`LX = 4000` km, 3 h). The x-mean `<rho u'w'>` was averaged over the final hour, from
every third snapshot, in uPa:

| run | 1.25 km | 2.75 km | 10.25 km | 30.25 km | 50.25 km | 55.25 km | `|MF(55)/MF(2.75)|` | `w'` rms at 55 km [m/s] |
|---|---|---|---|---|---|---|---|---|
| U20 h0.1 | -0.085 | -0.028 | +0.024 | +0.037 | +0.022 | +0.007 | 0.26 | 0.006 |
| U20 h1.0 | -62.2 | -8.38 | +1.68 | +3.60 | +2.20 | +0.72 | 0.086 | 0.063 |
| U20 h3.0 | -966 | -132 | -3.00 | +26.5 | +19.3 | +5.22 | 0.039 | 0.18 |
| U40 h0.1 | -0.621 | -0.140 | +0.093 | +0.149 | +0.087 | +0.030 | 0.21 | 0.013 |
| U40 h1.0 | -356 | -49.1 | +2.79 | +13.2 | +8.77 | +2.64 | 0.054 | 0.13 |
| U40 h3.0 | -4473 | -169 | -96.9 | +68.1 | +55.7 | +11.3 | 0.067 | 0.34 |
| U57 h0.1 | -1.74 | -0.331 | +0.184 | +0.302 | +0.178 | +0.060 | 0.18 | 0.018 |
| U57 h1.0 | -848 | -115 | -1.30 | +24.5 | +17.6 | +4.83 | 0.042 | 0.17 |
| U57 h3.0 | -8527 | +1200 | -310 | +80.8 | +73.1 | +13.2 | 0.011 | 0.51 |

The flux launched below the 4.5 km critical level does not reach 55 km with its sign.
Above about 8 km the x-mean flux has the opposite sign in every run, and at 55 km its
magnitude is 0.2 to 8 % of the 1.25 km value. Consistent with linear theory, the
stationary launched wave is absorbed near 4.5 km. What reaches 55 km is a component
of opposite sign, plausibly non-stationary in a 3 h record. It is 25 to 35 % of the
flux at 50 km in the AMEC band (18 % for U57 h3.0), and it falls by about 70 % between
50 and 55 km, across the 53.3 km critical level.

For these records, then, the criterion of the steering plan (the fraction of the launch
flux reaching 55 km) does not justify a top swap on the production configuration. Two
qualifications apply. First, the component that does reach the layer is a substantial
share of the flux in the band that Stage 2 analyses, so a layer return of the size
measured here (`R` of 0.02 to 0.04, amplitude ratio 0.14 to 0.2) could alter the band
at the level of its own amplitude ratio. Second, the 4.5 km and 53.3 km critical
levels interrupt any ground-to-layer path for stationary waves, so the cavity of the
linear corner has no direct counterpart there. A top swap on Stage 2 would be
informative only on records longer than 3 h, and only for the non-stationary part.

### L1, the long record

Launched at 21:51 on 2026-09-22 on the S0 configuration for 36 h. The choice was
made on the measured returned share of S0, S1, S2 and T2 before T1 had finished. On
the same rule T1 would now be chosen, and it is listed under next steps. The first
200 snapshots of L1 reproduce S0 bit for bit.

Predictions for the incident level, fixed on 2026-09-23 after the run ended but
before any of its output was analysed:

- Cavity reading. The final-hour-to-hour `|MF_up|` over 10 to 55 km, relative to the
  steady linear value, stops decaying by about 25 to 30 h. Operationally, its slope
  over the last 6 h is below 0.2 % per hour, against 1.1 % per hour at 15 h. The
  phase of the coherent product `conj(A) B` at the band base (10 km) in the returned
  modes 6 to 12 stops rotating: it changes by less than 30 degrees between the 28 h,
  32 h and final-hour windows.
- Transient reading. The level keeps decaying at 0.5 % per hour or more through 36 h,
  and the phase at 10 km keeps rotating.

Either outcome decides whether the 15 h level differences between layers can be
quoted as properties of the layers. The convergence criteria of Section 5 (bin
slopes, band drift, peak-to-peak) are read as well.

The phase column (`phBA_lo`) was added to the per-mode table of
`analyze_linear_flux_constancy.py` for this reading. It is the phase of the window
mean of `conj(A) B`, so a branch that rotates within the one-hour window also lowers
its modulus.

A premise of the predictions above does not hold for this run. The 1.1 % per hour
decay was measured on the `h = 100` m record of §15, not on S0 (`h = 30` m). In S0 the
incident level rose over 12 to 15 h. The predictions are kept as they were written.

Run: 175341 steps to `t = 129600` s, converged, 15515 s of wall clock (4.31 h), about
half of it alone on the machine. One snapshot (`snapshot_0101900.nc`, `t = 20.9` h)
had HDF5 metadata that crashed netCDF4 on close. Its data read correctly and lie within
9 to 14 per cent of the midpoint of their neighbours, relative to the neighbour
difference. It was rewritten into the structure of the neighbouring file, and the
original is kept as `snapshot_0101900.nc.corrupt`. The probable cause is a write
over the SMB mount.

Incident level, `|MF_up|` over 10 to 55 km relative to the steady linear value, from
the half-hourly trace (slopes in per cent of the window mean per hour, one standard
error):

| window [h] | mean `up/lin` | slope [%/h] | band `R` | slope of `R` [1/h] |
|---|---|---|---|---|
| 12-15 | 0.910 | +0.64 +- 0.02 | 0.0269 | +0.00008 +- 0.00015 |
| 15-20 | 0.921 | +0.21 +- 0.07 | 0.0268 | -0.00012 +- 0.00008 |
| 20-26 | 0.930 | +0.12 +- 0.04 | 0.0253 | -0.00002 +- 0.00009 |
| 24-30 | 0.940 | +0.40 +- 0.05 | 0.0251 | -0.00006 +- 0.00008 |
| 26-32 | 0.943 | -0.01 +- 0.09 | 0.0249 | +0.00009 +- 0.00007 |
| 30-36 | 0.947 | +0.46 +- 0.10 | 0.0259 | +0.00042 +- 0.00004 |

The hourly means over 28 to 36 h are 0.947, 0.950, 0.945, 0.939, 0.939, 0.941, 0.949
and 0.961. The level rises by about 3.5 % from 15 to 28 h, which matches in size the
arrival of mode `n = 1` (4.8 % of the flux, top of the band at 26.8 h). After that it
varies by about +-1 % around 0.945, and it rises again over the last 2 h.

Against the readings fixed above:

- The transient reading is not met. The level does not decay at any time in this run.
- The cavity reading is met in part. The level stops trending over 26 to 32 h (slope
  -0.01 +- 0.09 % per hour). Over the last 6 h the slope is +0.46 +- 0.10 % per hour,
  above the 0.2 % per hour criterion, and in the direction opposite to the decay the
  criterion was written against. The phase of `conj(A) B` at 10 km locks in modes 8
  to 12, which change by 4 to 26 degrees across the 28 h, 32 h and final windows and
  lie between 7 and 49 degrees from 20 h on. Mode 7 changes by 31 degrees and mode 6
  by 55 degrees. Modes 2 and 3 lock near +-170 degrees.
- Band `R` is flat at 0.025 to 0.027 from 8 to 36 h. The returned-class mean in modes
  6 to 12 falls from 0.041 at 15 h to 0.034 to 0.037 over 20 to 36 h (split), and
  from 0.035 to 0.030 to 0.033 (ripple). The rise of the returned share seen at 15 h
  in S1, S2 and T2 is not present in S0 over 15 to 36 h.

A common phase across modes 7 to 12 fits a single reflecting height. All the
hydrostatic modes share `m = N/U`, so a reflection at the same level returns to 10 km
with the same phase whatever `k` is. It does not identify the height.

Against the convergence criteria of Section 5:

- Bin slopes. The unreturned class holds only mode 1 after 27.5 h. Its twenty one-hour
  bins give a trend of `-2e-05 +- 3.5e-04` per hour. That is within one standard error,
  but the last six bins carry a single mode that rises from 0.0004 to 0.030, so the
  criterion carries little information in this record.
- Band drift and peak-to-peak. Not met. The final-hour drift is -0.147 %/km and the
  peak-to-peak 11.5 %, against -0.007 %/km and 1.6 % at 15 h. From 17 h onward the
  peak-to-peak stays between 9 and 19 %, and the drift varies between -0.15 and +0.13
  %/km.

The mode decomposition of the final-hour net flux shows where the peak-to-peak comes
from. At 15 h the individual modes carry peak-to-peak variations of 2 to 10 per cent
of the linear value across the band (`n = 4` 9.8 %), but these cancel in the sum to
2.1 %. At 36 h the individual variations are smaller, 2 to 4 %, and they add to
11.2 %. Modes 1 and 2 are not the cause: removing them leaves 9.7 %. For a steady,
conservative, linear wave each mode's net flux is independent of height, whatever the
reflection. A per-mode variation of 2 to 4 % therefore indicates a part of each mode
that is not steady or not conservative. The alignment of those parts after 20 h is in
step with the phase locking above.

What this changes. The 15 h shape bound (-0.033 %/km at `h = 100` m; -0.021 %/km and
1.9 % in S0) was measured at a time when the per-mode variations cancelled. On this
record it is not a steady property of the solver. The shape row of the error budget
above should carry the 36 h figures of the S0 layer (drift up to 0.15 %/km,
peak-to-peak up to 19 %) until a longer record on another layer shows otherwise. The
incident level at 36 h, 0.95 +- 0.01 of the linear value, lies above the 0.79 to 0.92
of the 15 h layer comparison. That comparison was made while the level was still
rising in S0, so the 15 h level differences between layers include a part that comes
from the time of sampling.

### Next steps

1. A 36 h record on the T1 layer, about 6 h of wall clock if run alone, to follow the
   return of the layer that returns least until the slowest modes have settled.
2. A test that separates the layer from the rest of the solver. Force a single
   horizontal wavenumber, with a ridge replaced by one Fourier mode, so the returned
   share is measured for one `k` without the ridge spectrum, the fill-in or the
   mode-to-mode scatter. Then vary `K0` at fixed thickness. This isolates whether the
   solver's layer departs from the continuous theory, or whether the difference comes
   from the cavity and the time dependence.
3. The cavity reading predicts that the incident level oscillates with the layer's
   depth at a period of `lambda_z/2`. Two further base heights, 65 and 75 km, would
   test it.

## 8. Single-mode cavity closure test

Added on 2026-09-23 after the steering review of the Section 7 results. It is next
step 2 above, run on two layers.

### Design

`sweep_arsia_mons_stage2.jl linear mode=N` replaces the Gaussian ridge with the single
Fourier mode `h sin(2 pi N x / LX)`. The kinematic condition is then
`w = U k h cos(k x)`, applied through the existing `:linear` forcing path with
`omega = 0`. The environment variable `MAGNUSP_OUTPUT` sends the output to a local
disk.

| id | mode | `h` [m] | layer | top | horizon | directory (run under `~/MAGNUSP_runs`, moved to `output/` afterwards) |
|---|---|---|---|---|---|---|
| C0 | 6 (`lambda_x = 333` km) | 5 | 60-90 km, `K0 = 0.05` | 90 km | 18 h | `arsia_mons_linear18h_LX2000_mode6/` |
| C1 | 6 | 5 | 60-120 km, `K0 = 0.005` | 120 km | 18 h | `arsia_mons_linear18h_LX2000_mode6_K0_0.005_LZ120/` |

Mode 6 sits near the flux-weighted `lambda_eff = 353` km of the ridge. Its round trip
from a 60 km base to 9.75 km takes 8.9 h, so 18 h covers two. At `h = 5` m the steady
wave has `u'/U` of about 0.02 at 55 km. The bottom `w` of the first snapshots carries
mode 6 only, with amplitude 0.0034 m/s in the first cell against the imposed 0.0038 m/s
at the face.

Analysis: `python3 scripts/compare_sponge_reflection.py --cavity <run dirs>`. For each
hour it takes the hour-mean complex split amplitudes of mode 6 and reports four
quantities:

- The band reflected share, from the split and from the ripple.
- The measured upward flux over 10 to 55 km, divided by the flux the same forcing gives
  in the linear column with no damping layer and a radiating top (`level`).
- Two predictions of that level from the reflection measured in the run. The first
  follows the steering note. It reads `r0 = B/A` at the lowest level and predicts
  `|1 + r_open,0|^2 / |1 + r0|^2`. The second avoids the near-surface structure. At the
  band base (9.75 km) it takes `r_top = B/A` from the run, by the split and
  independently by the ripple fit referred to that level. It takes the reflection of
  the column below, `r_low`, from the undamped structure equation integrated upward
  from a rigid ground. The prediction is then
  `|1 - r_low r_top,open|^2 / |1 - r_low r_top|^2`.

The second form is needed because the MCD column is not uniform near the ground. In
the open column, `|B/A|^2` is `7e-02` to `1e-01` at the lowest level against
`8e-04` in the band. The near-surface `N(z)` structure therefore reflects between the
ground and 10 km about as much as the layer does. `r_low` has modulus 1.000 at
9.75 km (a lossless lower column) and phase 59 degrees. Put the theory's own `r_top`
into the second form and it returns a level of 0.42, against 0.41 from the forced
theory column directly, which checks the construction.

The split-based closure at the ground is close to an identity, since `A + B = Psi` by
construction and `Psi(0)` is imposed. It tests only that the solver keeps the ground
condition and carries the upward branch conservatively. The ripple-based closure at
9.75 km shares no algebra with the split and is the test proper.

### Predictions fixed before the runs had produced a usable hour

Linear Rayleigh-layer theory for mode 6, forced as in the runs (radiating top):

| id | band `R` | `|r_top|^2` at 9.75 km | phase | upward level against the open column |
|---|---|---|---|---|
| C0 | 0.301 | 0.296 | 121 deg | 0.411 |
| C1 | 0.0019 | (open column: `8e-03` band) | | 0.951 |

The measured ridge share for `n = 6` in S0 was 0.02 to 0.03 (Section 7, L1), about 10
times below the C0 theory. Readings:

1. Closure. Over the last 6 h (12 to 18 h) the measured level agrees with the
   ripple-based prediction at 9.75 km within 10 per cent in both runs. If it does, the
   level differences between layers are accounted for by interference between the
   layer's return and the lower column, and the cavity reading of Section 7 stands for
   a single `k`.
2. Magnitude. If the C0 band `R` approaches the theory (about 0.3), the low values in
   the ridge runs come from the ridge spectrum, the fill-in or the time dependence and
   not from the layer. If it stays near 0.02 to 0.04, the solver's layer returns about
   ten times less than the continuous theory for one `k`, and the difference lies in
   the discrete layer or the eq19 top.
3. Layer contrast. The C0/C1 ratio of measured levels matches the ratio of the two
   predicted levels built from each run's own measured `r_top`.

### Results

| id | steps | wall clock [s] | notes |
|---|---|---|---|
| C0 | 87551 | 3871.7 | local disk, beside C1 |
| C1 | 87576 | 5230.0 | local disk, `OPENBLAS_NUM_THREADS=1` (Section 7, cost) |

Both reached `t = 64800` s. Both stay linear throughout. The largest `|u'|` at any
height is 2.1 m/s (`u'/U` 0.05, at 65 km in C1), `min u` stays above 37.8 m/s, and
no grid point is statically unstable at 9 or 18 h. In C1 the wave is absorbed by about
95 km (`|u'|` about 0.1 m/s above). The imposed ground amplitude is kept: `|Psi(0)|`
over the imposed value is 0.94 to 1.03 over 12 to 18 h.

Means over the hour windows centred on 12 to 17 h (six windows):

| id | band `R`, split | band `R`, ripple | theory band `R` | measured level | predicted, ripple at 9.75 km | predicted, split at 9.75 km | predicted, split at ground | theory level | net-flux peak-to-peak |
|---|---|---|---|---|---|---|---|---|---|
| C0 | 0.042 | 0.033 | 0.301 | 0.942 (0.923-0.966) | 1.047 | 0.916 | 0.986 | 0.411 | 31-54 % |
| C1 | 0.047 | 0.036 | 0.0019 | 0.733 (0.652-0.768) | 0.620 | 0.595 | 0.726 | 0.951 | 46-63 % |

The band share rises in both runs from about 0.007 at 5 to 7 h to 0.03 to 0.05 at 12
to 17 h. The rise begins between 8 and 10 h in both, close to the 8.9 h round trip
from a 60 km base to 9.75 km. The two layers share that base. In C0 the share peaks at
0.048 at 13 h and falls to 0.033 at 17 h. In C1 it reaches 0.053 at 16 h. In C0 the
level stays between 0.92 and 0.97 from 10 h on. In C1 it falls from 0.97 at 8 h to
0.65 at 17 h.

Against the readings fixed above:

1. Closure. Not met at the 10 per cent level. With `r_top` from the ripple, the
   hour-by-hour misfit is 1 to 19 % in C0 and 1 to 25 % in C1. The window means
   differ by 10 % (C0, prediction high) and 18 % (C1, prediction low). With `r_top`
   from the split at 9.75 km, the means differ by 3 % and 23 %. The ground form
   differs by 5 % and 1 %, but as noted above it is close to an identity. The
   interference account does reproduce the direction and most of the size of the
   level change: C1 falls to 0.73 while the prediction from its own measured return
   gives 0.60 to 0.73, and C0 stays near 0.94 while its prediction gives 0.92 to
   1.05. The account is therefore consistent in sign and order of magnitude, but the
   fields are not steady (net-flux peak-to-peak of 31 to 63 per cent across the band)
   and it is not closed to 10 per cent.
2. Magnitude. For a single `k` the measured share is 0.03 to 0.05 in both layers. The
   theory gives 0.30 for C0 and 0.0019 for C1, a ratio of about 160. The measured
   shares differ by less than 15 per cent between the two layers. The share is about
   seven times below the theory in C0 and about 25 times above it in C1. The ridge
   spectrum and the fill-in are absent here, so they do not account for the
   difference from the theory. The layer strength and thickness do not set the
   measured share at this `k` either.
3. Layer contrast. Not met. The measured C0/C1 level ratio is 1.29. The predictions
   built from each run's own `r_top` give 1.69 (ripple) and 1.54 (split at 9.75 km).
   The theory gives 0.43, in the opposite direction.

Interpretation, hedged. The single-mode runs are consistent with the ridge results of
Section 7. The downward share in the band does not follow the Rayleigh-layer theory of
the layer, and its size is nearly the same for layers whose theoretical returns differ
by two orders of magnitude. The timing of its growth fits energy returning from about
the height of the common 60 km base. Neither layer produces instability or a nonlinear
amplitude. Candidates that remain are:

- A start-up transient. The forcing starts impulsively (`omega = 0` at full amplitude
  from `t = 0`), so the field carries non-stationary components. The lower column
  reflects them losslessly (`|r_low| = 1`), and the stationary split does not separate
  them. The net-flux peak-to-peak of 31 to 63 per cent at 12 to 17 h shows how much
  non-stationary content remains.
- A departure of the discrete layer from the continuous theory that does not depend
  on `K0` or thickness. Reading `src/sponge.jl` and `src/euler.jl:140-190`, the source
  enters once per step at rate `K` on all four conserved variables. That is the form
  the theory assumes, so no factor in the code as written explains the difference.
- Structure of the MCD column near 60 km that the theory represents but the split
  misreads while the field is still transient.

The level change between layers in C1 (to 0.73 of the open column) shows that a
reflected share of only 0.05 can move the upward level by about 25 per cent when its
phase is unfavourable. That supports the level row of the error budget in Section 7
as a real sensitivity, while the source of the return remains open.

### Next steps from this test

1. The same pair with the forcing ramped on (`:linear_ramp`, `w = A cos(k x) tanh(t/t0)`,
   already in `src/forcing.jl`), with `t0` of a few hours. If the share falls toward
   the theory in C0 and toward the open-column value in C1, the start-up transient is
   the source. About 1 h and 1.5 h of wall clock on local disk.
2. A longer single-mode record (36 h) on one layer, to see whether the net-flux
   peak-to-peak decays and the closure tightens as the transient leaves.
3. Only if both leave the share unchanged: an isolated test of the discrete layer, a
   single `k` in an isothermal column where the theory has no `N(z)` structure.

## 9. Background audit

A review of Sections 7 and 8 (2026-09-23) pointed out that the domain integrals in the
C0 and C1 logs change during the first hours in a way the wave cannot account for.
E falls by 0.57 per cent and KE by 0.73 per cent over 5 to 6 h, and `drho_L2` settles
at 0.019 (C0) and 0.024 (C1). The review proposed that the free column adjusts away
from the initial state while the damping layer holds its region at that state, and
that the mismatch at the 60 km base returns part of the wave. This section records
the audit, one step at a time.

### What drho_L2 measures

`print_diagnostics` (`src/diagnostics.jl:102-112`) calls `check_hydrostatic_drift`
(`src/diagnostics.jl:36-59`), which returns

    drho_L2 = sqrt( mean over physical cells of ((rho - rho_bar) / rho_bar)^2 ).

The mean is an unweighted count over all `nx * nz` cells. `rho_bar` is `Q_bar[1]`, and
`Q_bar = copy(Q)` is taken right after `init_state` (`src/run.jl:53-54`), so the
reference is the initial column with the uniform wind. The damping layer relaxes
toward the same `Q_bar` (`src/euler.jl:153-154, 185-190`). The diagnostic therefore
includes both the wave and any change in the x-mean column. Because it is relative
and unweighted, it weights the thin upper levels as much as the dense lower ones.

The C0 and C1 snapshots were split into an x-mean part and an x-varying (wave) part
(`rel = (rho - rho_bar)/rho_bar`, with the x-mean removed for the second part):

| Run | t [h] | drho_L2 | x-mean part | x-varying part | share of drho_L2^2 above 60 km |
|---|---|---|---|---|---|
| C0 | 0.10 | 0.0217 | 0.0217 | 1.6e-4 | 0.95 |
| C0 | 1.03 | 0.0171 | 0.0171 | 1.2e-4 | 0.06 |
| C0 | 9.99 | 0.0191 | 0.0191 | 1.3e-3 | 0.02 |
| C0 | 17.98 | 0.0189 | 0.0189 | 1.1e-3 | 0.02 |
| C1 | 0.10 | 0.133 | 0.133 | 2.0e-4 | 1.00 |
| C1 | 0.94 | 0.051 | 0.051 | 1.7e-4 | 0.94 |
| C1 | 9.90 | 0.0245 | 0.0244 | 1.8e-3 | 0.55 |
| C1 | 17.90 | 0.0235 | 0.0235 | 1.6e-3 | 0.51 |

The x-mean part accounts for more than 99.7 per cent of `drho_L2^2` at every time
sampled. The departure is present 6 min after the start, when the wave has not left
the lowest levels, so it develops on an acoustic time scale and not with the wave.
`drho_L2` in these runs measures the adjustment of the background column. It does not
measure the wave.

### Profile of the departure at 18 h

x-mean departures from the initial column at the last snapshot:

| z [km] | C0 rho [%] | C0 P [%] | C0 T [%] | C1 rho [%] | C1 P [%] | C1 T [%] |
|---|---|---|---|---|---|---|
| 0.25 | -0.55 | -0.05 | +0.50 | -0.55 | -0.05 | +0.50 |
| 9.75 | -0.58 | -0.46 | +0.11 | -0.58 | -0.47 | +0.11 |
| 19.75 | -1.05 | -0.95 | +0.10 | -1.05 | -0.96 | +0.09 |
| 39.75 | -2.46 | -2.08 | +0.39 | -2.46 | -2.16 | +0.31 |
| 54.75 | -4.01 | -1.06 | +3.07 | -3.98 | -1.60 | +2.48 |
| 57.75 | -4.34 | +0.01 | +4.55 | -4.28 | -0.80 | +3.65 |
| 59.75 | -2.02 | +0.96 | +3.04 | -4.50 | +0.05 | +4.76 |
| 61.75 | +0.68 | +1.07 | +0.39 | -2.16 | +0.92 | +3.15 |
| 64.75 | +0.59 | +0.83 | +0.24 | +0.45 | +1.21 | +0.75 |
| 79.75 | -0.26 | -0.26 | -0.01 | -0.14 | +0.18 | +0.32 |

The mean wind changes by less than 0.13 m/s at every level. Below 50 km the two runs
agree to about 0.03 percentage points in density, although their layers differ in
`K0` by a factor of 10 and in thickness by a factor of 2. Density decreases steadily with
height, and the departure ends within 2 to 5 km of the 60 km base in both runs. There
it gives way to a density step of about 4.5 per cent and a warm layer of +3 to +5 per
cent in T across 55 to 62 km. The column loses 0.74 per cent of its mass. With the
wind unchanged, that matches the 0.73 per cent fall in KE. The path by which the mass
leaves has not been identified.

### A term in the well-balanced flux that is consistent with the drift

The vertical flux subtracts a background pressure, `rho*w^2 + (P - P_bar)`
(`src/euler.jl:27-34`). `P_bar` is formed as `(gamma - 1) * Q_bar[4]` at
`src/euler.jl:143-144` and `171-172`, and in the same way at `src/boundary.jl:138`.
`Q_bar[4]` is the total energy, `P0/(gamma - 1) + rho0 U^2/2`. The subtracted
pressure is therefore `P0 + (gamma - 1) rho0 U^2 / 2`, not `P0`. With `U = 0` the two
agree, and they differ once there is a mean wind. For the initial state
the unbalanced vertical force per unit mass is

    (gamma - 1) U^2 / (2 H_rho),

which amounts to an effective gravity `g (1 + eps)` with
`eps = (gamma - 1) U^2 / (2 g H_rho)`. With `gamma = 1.31`, `g = 3.72 m/s^2`,
`U = 40 m/s` and the MCD `H_rho`, `eps` is 0.006 to 0.011 over 0 to 55 km. At fixed
temperature, a column in balance under `g (1 + eps)` has its `d ln rho / dz` changed
by `-eps / H_p`:

| Quantity, 5 to 50 km | C0 | C1 |
|---|---|---|
| Observed `d ln rho / dz` change [%/km] | -0.069 | -0.069 |
| Observed `d ln P / dz` change [%/km] | -0.048 | -0.052 |
| Predicted `-eps / H_p` [%/km] | -0.072 | -0.072 |

The observed density gradient change is within 5 per cent of the estimate, which
contains no fitted parameter. This agreement is consistent with the offset in
`P_bar` driving the adjustment, but it does not establish it. The estimate ignores
the temperature change of 0.1 to 0.5 per cent, and the pressure gradient
matches less closely. A direct test is a short `h = 0` run, once with the current
`P_bar` and once with `P_bar = (gamma - 1) (Q_bar[4] - (Q_bar[2]^2 + Q_bar[3]^2) / (2 Q_bar[1]))`.
The corrected run should keep `drho_L2` near zero. The same term is present in every
run with a mean wind, including Stage 2 and the earlier linear-corner records.

In the damping layer the relaxation toward `Q_bar` holds the unadjusted state (80 to
90 km in C0 departs by less than 0.3 per cent). That gives the 55 to 62 km interface
the review proposed. At 4.5 per cent it is larger than the review's order-1-per-cent
estimate, and it is the same in both layers. Whether this interface returns the
observed share at `n = 6` is the question for the next audit steps.

## 10. Summary figures

`scripts/meeting_figures.py` writes two figures to `output/meeting_figures/`. They
are drawn from existing output and need no new runs. `--gif` adds an animation of
the second figure's run (144 frames at 15 min spacing, about 13 s to render).

### Figure 1, benchmark of the top boundary and convergence

`fig1_benchmark_sponge_convergence.png`. All runs use `h = 30 m`, `sigma = 50 km`,
`U = 40 m/s`, the MCD MY31 column, and the steady linear flux `MF_lin = -2.25 uPa` as
the normaliser.

- (a) Final-hour flux profiles. Every layer removes the flux within about 10 km of
  its base: `MF/MF_lin` at 85 km is below 0.002 in all runs. At 70 km it is 0.001
  (S0), 0.008 (S2), 0.045 (S1) and 0.31 (T1, a weaker and thicker layer), and 0.86 in
  T2, whose base is at 70 km. The band profile is flat to 1.9 to 4.0 per cent
  peak-to-peak at 15 h in the four linear layers, and to 11.4 per cent in L1 at 36 h.
- (b) Band-mean flux against time. All layers follow the same curve up to about 7 h
  (half of the linear value at 3.8 h), which is before any return from the top reaches
  the band. After that they separate. At 15 h the band-mean levels are 0.89 (S0),
  0.84 (S2), 0.79 (T1), 0.75 (S1) and 1.03 (T2, outside the linear range). L1
  continues S0 and reaches 0.935 at 35 to 36 h.
- (c) Measured returned share (modes 6 to 12, split estimator) against the linear
  Rayleigh-layer prediction. The measured values span 0.017 to 0.14, while the
  predictions span 0.0018 to 0.30. The single-mode runs C0 and C1 (stars) give 0.042
  and 0.047 against predictions of 0.30 and 0.0019.

Points to raise with this figure:

1. The test design. A linear wave with a known steady flux makes both the flatness of
   the flux profile and its level quotable against theory. The damping layer is then
   varied in strength, thickness and base height at fixed everything else.
2. Result. The layers absorb what reaches them, and the band settles to within 25 per
   cent of the linear value. The level depends on the layer by up to about 15 per cent
   among the linear layers, which enters the error budget as a separate row (Section 7).
3. Open issue. The returned share does not scale with the layer theory. Section 9 finds
   a drift in the background column that the layer does not share (a 4.5 per cent
   density step at 55 to 62 km). The data are consistent with this being the common
   source, and the next tests are designed to check it.
4. Caveats. The theory rows for `K0 = 0.05` lie outside the range of the anelastic
   calculation. The flux level also carries the `N0` reference factor (Section 7).

### Figure 2, the showcase run

`fig2_showcase_L1_wave.png` (and `.gif`). L1 is chosen because it is the longest record
(36 h), it is linear and has been checked against theory at each stage, and it is the
run the benchmark refers to. The Stage 2 runs use the real ridge height, but they are
3 h records whose flux in the 45 to 53 km band has the opposite sign to the launched
wave. That result is not yet explained, so those runs are less suited to showing what
the model does.

- Top row: `w'` scaled by `(rho/rho_s)^(1/2)` at 1, 3, 6, 12 and 36 h, with the x-mean
  removed. The scaling offsets the growth of amplitude with height, so the phase lines
  keep a similar colour depth. The x-mean is removed because at 1 h it is dominated by
  the background adjustment of Section 9, a domain-wide downward `w` of the same order
  as the wave. The wave fills the column from the ground upward in about 6 h. The
  phase lines tilt upstream with height, with `lambda_z` near 30 km, and they fade
  above the 60 km layer base.
- Bottom: time-height record of the x-mean flux over 36 h. The front rises to 60 km
  by about 7 h. The curved features after 15 h are arrivals and returns of the longer
  modes, the non-stationary content behind the 9 to 19 per cent peak-to-peak
  variation in Section 7.

### Time evolution and the N^2 interface

`compare_sponge_reflection.py --background` differences the x-mean rho, T, u and N^2
of each run against its initial column. N^2 is formed as `g d ln(theta)/dz` from the
x-mean T and P. The audit was run on C0 and C1 at 0.5, 6 and 18 h and on L1 (S0 layer,
ridge forcing, h = 30 m) at 0.5, 6, 18 and 36 h. Table in
`output/reflection_tests/logs/background_audit.txt`, profiles in
`output/reflection_tests/background_audit.json`, figure
`output/reflection_tests/background_audit.png`.

The review's prediction was a departure of order 1 per cent below 60 km that stops
at the layer base in both runs, leaving an N^2 feature across 58 to 62 km. The 18 h
C0 and C1 density profiles were already known from the previous step, so the test
is not blind for those two. The time evolution, N^2, u and L1 were not known.

x-mean N^2 departure near the base, last snapshot of each run:

| z [km] | C0, 18 h | C1, 18 h | L1, 36 h |
|---|---|---|---|
| 50.25 | +8.5 % | +7.1 % | +8.7 % |
| 55.25 | +13.6 % | +11.1 % | +14.1 % |
| 58.25 | +24.1 % | +19.2 % | +9.4 % |
| 59.25 | -92.0 % | +16.3 % | -92.3 % |
| 59.75 | -91.4 % | +5.4 % | -98.6 % |
| 60.25 | -80.0 % | -15.1 % | -93.3 % |
| 60.75 | -86.0 % | -38.6 % | -50.2 % |
| 61.75 | +6.9 % | -58.5 % | -2.6 % |
| 62.75 | +1.0 % | -41.2 % | +0.7 % |
| 63.75 | +0.5 % | -18.9 % | +0.5 % |
| 64.75 | +0.3 % | -7.3 % | +0.4 % |

Observations:

1. **Below 50 km the adjustment is common to all three runs and is complete by 6 h.**
   Density departures at 6 h are within 0.01 percentage points of those at 18 h
   (C0) and 36 h (L1), and C0, C1 and L1 agree to 0.03 points or better. At 0.5 h
   roughly half of the 18 h departure is in place below 40 km. Ridge forcing at
   h = 30 m (L1) and single-mode forcing at h = 5 m (C0) give the same profile, which
   is consistent with a background adjustment that does not depend on the wave.
2. **The layer base carries a thin near-neutral layer.** In C0 and L1, N^2 falls to
   2 to 9 per cent of its initial value over 59.0 to 61.0 km (four to five grid
   levels), under a region 50 to 58 km where N^2 is raised by 8 to 24 per cent. The
   same interface is already present at 0.5 h (-44 per cent at 61.75 km) and is
   fully formed by 6 h. In C1, where `K` rises more slowly above the base, the
   deficit is shallower (-59 per cent at most) and broader (60 to 65 km), and it
   sits about 2 km higher. It forms between 0.4 and 6 h.
3. **Temperature carries a warm bump of 3 to 5 per cent peaking just below the base**
   (C0 +4.6 per cent at 57.75 km, C1 +4.8 per cent at 59.75 km), which rises
   gradually from 35 km and drops sharply at the base. The N^2 structure in item 2
   is the gradient signature of this bump.
4. **The mean wind is nearly unchanged in the single-mode runs.** |Δu| < 0.01 m/s
   below 50 km and < 0.13 m/s at the base. In L1 the base region decelerates by
   0.2 m/s at 36 h, which is consistent with wave deposition at the layer base and
   is not part of the background adjustment.
5. **In C1 the upper part of the thick layer is not held at the initial state.** At
   0.4 h density at 120 km is +47 per cent; by 6 h it settles to -5 per cent and
   continues the free-column trend from 90 to 120 km. With `K0 = 0.005` the
   relaxation time is 200 s at the top, so this departure is maintained against
   the relaxation, possibly through the top boundary. It was not examined further.

Against the review's prediction, the departure stops at the base in both layers and
leaves an N^2 feature, as predicted. The amplitude is larger than predicted: 4.3 per
cent in density and 4.8 per cent in T just below the base, not about 1 per cent. The
feature is not the same in C0 and C1: it is deeper and thinner in C0 and shallower
and broader in C1. Its N^2 deficit, integrated over the base region, is similar in
the two runs. The C0 deficit is equivalent to a 2.0 km layer with N^2 = 0, and the
C1 deficit to a 1.7 km layer.

A rough scale for what such a layer returns: for a layer thin against half a
vertical wavelength (13 km at `m = N/U = 2.4e-4 1/m`), the Born estimate is
`|r| ≈ m d_eq / 2`, with `d_eq` the equivalent zero-N^2 thickness. This gives
`R ≈ 0.056` for C0 and `R ≈ 0.041` for C1. The measured returned shares are 0.042 and
0.047. The estimate ignores the positive N^2 lobe over 50 to 58 km, which is not
thin against the wavelength and can add to or cancel the deficit's contribution
depending on phase. It is an order-of-magnitude guide and does not replace the
calculation on the adjusted column in the next step. Taken with the uniform
adjustment below 50 km, the interface is a candidate source that would be the
same in form, and similar in size, in every layer with a 60 km base. It is linear
in the wave amplitude and in place within the first hours. That matches the
properties of the return listed in Section 8.

Because the interface is one to ten grid levels thick, any return it produces would
depend on vertical resolution. That is a further reason to remove it at the source
rather than to characterize it.

### Theory on the adjusted column

The damping-layer theory of Section 3 and the cavity test of Section 8 used the
initial MCD column. Two changes are tested here, for `n = 6` and the layers of C0
and C1:

- **Adjusted column.** N^2 and rho are taken from the x-mean of the last snapshot
  of each run, in place of the initial column. N^2 is formed as
  `g ((1/gamma) d ln P/dz - d ln rho/dz)`, which does not assume hydrostatic balance.
  That matters because the adjusted column is in balance under `g (1 + eps)`, not `g`.
  The sponge term `K (Q_bar - Q)` splits into a steady part that holds the adjusted
  background and a part `-K q'` on the wave, so the linear theory on the adjusted
  column with the same `K(z)` is the consistent one.
- **Compressible term.** Uniform Rayleigh damping on every prognostic variable
  replaces every `-i Omega` by `-i (Omega + i K)`. The isothermal compressible
  structure equation then gains `(Omega + i K)^2 / c^2` beside the anelastic terms,
  with `c^2 = gamma P / rho`. For `K = 0.05 1/s` this term exceeds the buoyancy term
  `k^2 N^2 / (Omega + i K)^2` by more than three orders of magnitude, so in the strong
  layers the anelastic form may not describe the wave.

Predictions fixed before the calculation:

1. **Interface.** With `K0 = 0` on the adjusted column (the interface is then the
   only feature that departs from the initial column), R in 10 to 55 km should be
   0.02 to 0.06 for both runs if the interface accounts for the measured return.
   The Born estimate of the previous subsection gives 0.056 (C0) and 0.041 (C1).
2. **Full model.** Adjusted column, own layer, compressible term: the interface
   hypothesis passes if both runs give 0.02 to 0.06, within a factor of 1.5 of the
   measured 0.042 (C0) and 0.047 (C1).
3. **Compressibility.** On the initial column with the C0 layer, the compressible
   term accounts for the overstatement of the strong layer if R falls from 0.30 to
   below about 0.1. For the C1 layer (`K0 = 0.005`) the change should be small.

#### Results

`compare_sponge_reflection.py --theory` (log `output/reflection_tests/logs/adjusted_theory.txt`,
JSON `output/reflection_tests/adjusted_theory.json`). R is the mean over 10 to 55 km, and
the top is radiating in every case. The initial-column values reproduce Section 8
(C0 layer 0.3012 against 0.3011, C1 layer 0.0019, C0 open column 0.00077). The C1 open
column gives 0.0086 against 0.0084 there, a difference that comes from the N^2 formula.

| Column | C0, K0 = 0 | C0 layer | C0 layer, compressible | C1, K0 = 0 | C1 layer | C1 layer, compressible |
|---|---|---|---|---|---|---|
| Initial | 0.0008 | 0.301 | 0.333 | 0.0086 | 0.0019 | 0.0019 |
| Adjusted, 6 h | 0.047 | 0.448 | 0.493 | 0.0015 | 0.0071 | 0.0070 |
| Adjusted, 18 h | 0.065 | 0.443 | 0.491 | 0.034 | 0.031 | 0.031 |
| Measured, 12 to 17 h (split / ripple) | | 0.042 / 0.033 | | | 0.047 / 0.036 | |

Hourly theory on the adjusted column against the measured band share:

| t [h] | C0 theory, K0 = 0 | C0 measured, split / ripple | C1 theory, own layer | C1 measured, split / ripple |
|---|---|---|---|---|
| 3 | 0.050 | 0.020 / 0.000 | 0.003 | 0.020 / 0.000 |
| 6 | 0.047 | 0.008 / 0.001 | 0.007 | 0.007 / 0.000 |
| 9 | 0.048 | 0.016 / 0.008 | 0.014 | 0.010 / 0.006 |
| 12 | 0.045 | 0.045 / 0.032 | 0.020 | 0.034 / 0.027 |
| 14 | 0.045 | 0.047 / 0.041 | 0.024 | 0.047 / 0.036 |
| 16 | 0.065 | 0.039 / 0.025 | 0.028 | 0.053 / 0.043 |

The split values before about 4 h come from the start-up and have no ripple
counterpart. They are not read as reflection.

Reading against the predictions:

1. **Interface: met for C1, marginal for C0.** On the adjusted column with no
   damping, C1 gives 0.034 at 18 h. C0 gives 0.045 to 0.052 from 2 to 14 h and
   0.065 from 15 h, when the deficit deepens by one grid level. The range fixed
   in advance was 0.02 to 0.06.
2. **Full model: met for C1, not met for C0.** C1 gives 0.031, a factor of 1.5
   below the split and 1.2 below the ripple. C0 gives 0.49, twelve times the
   measured share.
3. **Compressibility: not met.** The compressible term raises the C0 layer's value
   from 0.30 to 0.33. It does not lower it. For C1 it changes nothing. In the
   `K0 = 0.05` layer, `K` passes `|Omega| = 7.5e-4 1/s` 2.3 km above the base,
   while the compressible term dominates only from about 7 km above it. A plausible
   reading is that the theory's reflection forms in the lower few km, before the
   compressible term matters. So the anelastic form does not appear to be why the
   theory overstates the strong layer.

Interpretation. Under the assumptions of the linear steady theory, the interface
that the background adjustment leaves at the layer base produces a returned share
of the measured order in both runs (0.03 to 0.065 against 0.03 to 0.05), with no
fitted parameter. In C1 the interface forms slowly: the theory rises from 0.0003 at
5 h to 0.034 at 18 h, and the measured share rises over the same hours. That is
consistent with the interface being the main source of the C1 return and of its
growth after 8 h. The measured share runs ahead of the theory by a factor of about
2 at 12 to 16 h. That is not what a causal lag would give, since the wave takes
about 4 to 5 h to travel from the base to the band. The steady theory may therefore
understate the return of an interface that is still forming.

In C0 the interface-only value matches the measured share, while the continuous
theory of the strong layer predicts 0.44 to 0.49 on the same column. So the
difference between the continuous theory and the solver in the strong layers,
first seen in Section 7, remains and is not due to compressibility. The data are
consistent with the discrete `K0 = 0.05` layer reflecting much less than the
continuous theory, while the interface below it reflects about as the theory says.
This section does not identify why.

The row "theory overstates by 6 to 7 times" in Section 7 therefore stands for the
strong layers. The weak-layer excess (C1, T1, S1) that Section 7 left unexplained
has a candidate source in the adjusted interface that is quantitative to within
a factor of about 1.5 to 2.

### Background-pressure switch test

**Design.** A switch `pbar_mode` in `SimConfig` selects the background pressure
subtracted in the vertical flux:

- `:total` is the historical behaviour and stays the default, `(gamma - 1) Q_bar[4]`.
- `:pressure` subtracts the background pressure alone,
  `(gamma - 1) (Q_bar[4] - (Q_bar[2]^2 + Q_bar[3]^2) / (2 Q_bar[1]))`.

The sponge target is not changed: relaxing E toward the initial E is correct. The
switch covers every site in `src/` that forms a pressure from `Q_bar[4]`:

- the vertical flux in `src/euler.jl`,
- the Klemp-Durran top in `src/boundary.jl`,
- the characteristic top in `src/boundary.jl`.

The positivity floors in `src/run.jl` and `src/diffusion.jl` are not changed; they
use `Q_bar[4]` only as a scale for a floor. With `:total` every earlier record
(S0 to T2, L1, C0, C1, Stage 1 and 2) remains reproducible bit for bit, and this is
checked against C0 before the test runs.

Pair, wave-free (`h = 0`), C0 configuration (60 to 90 km, K0 = 0.05, `mode=6`),
6 h, run in parallel on local disk with `OPENBLAS_NUM_THREADS=1`:

- Z0, `pbar_mode = :total`.
- Z1, `pbar_mode = :pressure`.

A mass audit reports, for both, the domain mass and the time-integrated mass
exchange through the eq19 top, through the bottom and through the damping layer
(`sum K (rho_bar - rho)`). A `P_bar` error changes the momentum balance and not
continuity, so it should redistribute mass but not remove it. If Z1 removes the
drift but the 0.74 per cent mass loss remains, there is a second issue, tied to
a boundary or to the layer.

**Criteria fixed before launch** (5 to 50 km for the gradient; x-mean profiles at 6 h):

| Quantity | Z0, reproduces the drift without the wave | Z1, pass |
|---|---|---|
| change in `d ln rho/dz`, 5 to 50 km | about -0.07 %/km | magnitude below 0.01 %/km |
| minimum of N^2 / N^2_initial over 59 to 61 km | below 0.10 | above 0.80 |
| `drho_L2` at 6 h | about 0.019 | below 0.002 |

A Z1 result between these columns counts as a partial pass. It would point to a
remaining discrete imbalance in the initial column, and the next step would be
to initialize the column by discrete hydrostatic integration.

If Z1 passes, fixed-`P_bar` reruns of C0 and C1 (18 h each) follow in parallel:

- **C1 fixed.** The interface hypothesis passes if the band share over 12 to 17 h is
  0.01 or less. It is rejected if the share stays near 0.04. Values between 0.002
  and 0.008 cannot be told apart, since the initial 120 km column gives 0.0086
  with no damping and the split reads about 0.007 in the transient field.
- **C0 fixed.** Both outcomes are informative. A share near 0.30 would support the
  continuous layer theory. A share near 0.01 would indicate that the solver's
  `K0 = 0.05` layer absorbs almost completely, and that the difference from the
  theory lies in the discrete layer.

**Flag for production.** In a sheared wind the spurious term becomes
`d/dz[(gamma - 1) rho_bar u^2 / 2]`, which includes `(gamma - 1) rho_bar u du/dz`. With
Stage 2 shears of order `4e-3 1/s` that term can exceed the uniform-wind term several
times, and it changes sign across the critical levels. Before any Stage 2 rerun, an
estimate of the effective-gravity error on the Stage 2 wind profile should show
whether the 4.5 km and 53.3 km critical-level structure, or the opposite-sign
AMEC-band flux, could be affected.

**Mass audit of C0 and C1, before the pair finished.** `compare_sponge_reflection.py --mass`
on the existing 18 h records (percentages of the initial domain mass, cumulative):

| t [h] | C0 dM | C0 top out | C0 layer in | C0 bottom out | C1 dM | C1 top out | C1 layer in | C1 bottom out |
|---|---|---|---|---|---|---|---|---|
| 0.5 | -0.189 | -0.0005 | +0.023 | +0.213 | -0.188 | 0.0000 | +0.004 | +0.192 |
| 2 | -0.532 | -0.002 | +0.031 | +0.564 | -0.538 | 0.0000 | +0.011 | +0.549 |
| 6 | -0.721 | -0.005 | +0.032 | +0.758 | -0.726 | 0.0000 | +0.022 | +0.748 |
| 18 | -0.735 | -0.016 | +0.023 | +0.774 | -0.738 | +0.0001 | +0.032 | +0.769 |

x is periodic and the density floor does not act. The bottom-face exchange is
therefore obtained by elimination (`-(dM + top out - layer in)`); cell-centre
snapshots cannot resolve it directly. The eq19 top passes almost nothing: it
returns 0.016 per cent inward in C0. The layer relaxation adds 0.02 to 0.03 per
cent. About 0.77 per cent leaves through the bottom face, with a time profile
matching that of the drift.

Under the mirror condition `rho*w` cancels at the face. The Lax-Wendroff half step
adds `-(dt / 2 dz)` times the difference in `(P - P_bar)` across the face. With
`P_bar` carrying `(gamma - 1) rho_bar U^2 / 2`, that difference is nonzero between the
ghost and the first level, and this term is a candidate path for the loss. This
suggests one further prediction, written before the Z runs were read: **in Z1 the
mass loss at 6 h falls below 0.07 per cent** (a tenth of Z0). If it does not, the
leak is a second issue of the bottom boundary, independent of `P_bar`.

**Bit-for-bit check.** With the switch in place and the default `:total`, a 0.1 h
rerun of the C0 command reproduced C0 snapshots 100 to 400 bit for bit (all of
rho, u, w, P and T, and the times). Wall clock was 22 s.

**Runs.**

| id | command (under `MAGNUSP_OUTPUT=$HOME/MAGNUSP_runs OPENBLAS_NUM_THREADS=1`, `julia -t 1`) | steps | wall clock [s] | directory under `output/` |
|---|---|---|---|---|
| Z0 | `linear mode=6 h=0 hours=6` | 29184 | 1256.5 | `arsia_mons_linear6h_LX2000_mode6/run_U40.0_h0.0` |
| Z1 | `linear mode=6 h=0 hours=6 pbar=1` | 29165 | 1255.6 | `arsia_mons_linear6h_LX2000_mode6_pbarP/run_U40.0_h0.0` |

The two ran in parallel at about 23 steps/s each. Logs:
`output/reflection_tests/logs/Z0_h0_6h.log`, `Z1_h0_6h.log`. Profiles:
`output/reflection_tests/zpair_background_audit.{json,png}`. Mass:
`zpair_mass_audit.json`, with the C0 and C1 audit in `mass_audit_C0C1.json`.

**Results against the criteria** (6 h):

| Quantity | Z0 | Z0 criterion | Z1 | Z1 criterion |
|---|---|---|---|---|
| change in `d ln rho/dz`, 5 to 50 km | -0.069 %/km | about -0.07 | +0.0000 %/km | below 0.01 in magnitude |
| minimum of N^2 / N^2_initial, 59 to 61 km | -0.52 | below 0.10 | 1.000 | above 0.80 |
| `drho_L2` | 0.0190 | about 0.019 | 5.2e-14 | below 0.002 |
| domain mass change | -0.720 % | | -0.0000 % | below 0.07 % (mass prediction) |
| largest x-mean density departure | 4.46 % | | 3.5e-11 % | |

All criteria are met, Z0 at its expected values and Z1 with a large margin. In Z1
the initial column stays a fixed point of the discrete scheme to round-off
(`drho_L2` grows from 5e-15 to 5e-14 over 6 h). So the adjusted profile, the N^2
interface and the mass loss all come from the kinetic-energy part of `P_bar`. No
separate discrete imbalance in the initial column remains, and the fallback of a
discrete hydrostatic initialization is not needed. Without a wave, Z0 reproduces
C0's mass loss (-0.720 against -0.721 per cent at 6 h) and its gradient change
(-0.069 against -0.069 %/km). This supports the earlier finding that the wave plays
no part in the drift.

The mass audit identifies the path of the loss. In Z0 the bottom face carries 0.76
per cent out by 6 h, and in Z1 it carries nothing. Under `:total` the Lax-Wendroff
half step at the bottom face sees a pressure difference across the mirror that is
not balanced, which is consistent with the loss leaving there.

**The interface can be statically unstable.** At 6 h the x-mean N^2 in Z0 is
negative at 60.25 km (-152 per cent of its initial value) and at 60.75 km (-112 per
cent): the background itself is statically unstable in the lowest kilometre of the
layer, with no wave present. The instability that `linear_flux_constancy.md` §18
found at 58.75 to 61.75 km "on the sponge base", and that Section 7 found to follow
the strong layer, sits on this interface. The data are consistent with that
instability being a product of the `P_bar` term rather than of the wave. That
reading needs checking in the fixed C0 rerun, where the same diagnostics can be
repeated.

**Stage 2 exposure estimate.** Under `:total`, the unbalanced vertical force per unit
mass in the initial column is `(1/rho_bar) d/dz[(gamma - 1) rho_bar u^2 / 2]`
`= (gamma - 1) [u du/dz + (u^2 / 2) d ln rho_bar/dz]`. Expressed as a fractional change of
gravity, `eps(z)`, on the MCD `u_zonal` that Stage 2 uses:

| z [km] | u [m/s] | du/dz [m/s per km] | eps, shear part | eps, total | eps, uniform 40 m/s |
|---|---|---|---|---|---|
| 0.25 | +0.2 | +0.5 | -0.001 % | -0.001 % | +1.07 % |
| 7.75 | -9.2 | -3.8 | -0.29 % | -0.26 % | +0.60 % |
| 14.75 | -43.5 | -4.1 | -1.49 % | -0.77 % | +0.61 % |
| 24.75 | -57.4 | +0.0 | +0.02 % | +1.33 % | +0.64 % |
| 44.75 | -43.8 | +1.2 | +0.43 % | +1.36 % | +0.77 % |
| 49.75 | -25.2 | +5.9 | +1.24 % | +1.56 % | +0.80 % |
| 53.25 | -0.4 | +7.2 | +0.02 % | +0.02 % | +0.81 % |
| 57.75 | +29.1 | +5.4 | -1.30 % | -0.86 % | +0.83 % |
| 59.75 | +38.2 | +3.6 | -1.13 % | -0.37 % | +0.84 % |
| 69.75 | +36.2 | -5.0 | +1.52 % | +2.25 % | +0.89 % |

Below 60 km, `|eps|` reaches 2.0 per cent (at 47.8 km, inside the 45 to 53 km AMEC
band). Its rms is 1.6 times that of the uniform 40 m/s case, and it changes sign
near 2, 4, 17 and 53 km. The shear part is of the same order as the uniform-wind
term, not several times larger, and it changes sign across the critical levels as
expected. Near the ground `u` is close to zero, so the bottom-face mass leak seen
here should be small in Stage 2.

Under the assumptions of this estimate, the Stage 2 background is subject to a
spurious effective gravity of the same size as the one that produced the 4.5 per
cent density and 5 per cent temperature departures here. Its sign alternates with
height, and its largest values sit in the AMEC band and across the 53.3 km
critical level. This does not show that the Stage 2 critical-level structure or
the opposite-sign AMEC-band flux are artefacts. It does indicate that neither can
be quoted from `:total` runs without a `:pressure` comparison. Any Stage 2 rerun
should use `pbar=1`.

### Fixed-P_bar reruns of C0 and C1

| id | command (`pbar=1` added to the C0 or C1 command) | steps | wall clock [s] | directory under `output/` |
|---|---|---|---|---|
| C0P | `linear mode=6 h=0.005 hours=18 pbar=1` | 87516 | 3791.3 | `arsia_mons_linear18h_LX2000_mode6_pbarP/run_U40.0_h0.005` |
| C1P | `linear mode=6 h=0.005 hours=18 pbar=1 k0=0.005 lz=120` | 87516 | 5088.5 | `arsia_mons_linear18h_LX2000_mode6_K0_0.005_LZ120_pbarP/run_U40.0_h0.005` |

The two ran in parallel on local disk with `OPENBLAS_NUM_THREADS=1` and `julia -t 1`.
C1P started 2 min after C0P: the first launch passed `k0=0.005 lz=120` as a single
argument, because zsh does not split an unquoted variable, and exited at once. Both
wall clocks are close to those of C0 and C1 (3872 s, 5230 s). Logs:
`output/reflection_tests/logs/C0P_mode6_pbarP.log`, `C1P_mode6_pbarP.log`. Closure:
`logs/fixed_cavity_closure.txt`, `fixed_cavity_closure.{json,png}`. Background:
`fixed_background_audit_C0P.json`, `fixed_background_audit_C1P.json`.

**The background stays at the initial column.** At 18 h the x-mean density departs
by less than 0.07 per cent at every level in C0P and 0.023 per cent in C1P. The
N^2 ratio over 59 to 61 km is 0.99 or above, and nowhere is N^2 negative. The only
mean change is a deceleration of 0.16 m/s (C0P) and 0.15 m/s (C1P) at the layer
base, consistent with the wave depositing momentum there. `drho_L2` is now set by
the wave (1e-3 to 2e-3).

**Band shares and levels, means over 12 to 17 h (± standard deviation over hours):**

| Run | split | ripple | theory (layer) | level measured | level from r0 (split) | level from r_top (ripple) | net peak-to-peak |
|---|---|---|---|---|---|---|---|
| C0, `:total` | 0.042 ± 0.005 | 0.033 ± 0.006 | 0.301 | 0.942 ± 0.016 | 0.986 | 1.047 | 45 ± 7 % |
| C0P, `:pressure` | 0.019 ± 0.001 | 0.016 ± 0.002 | 0.301 | 1.093 ± 0.038 | 1.331 | 1.191 | 35 ± 7 % |
| C1, `:total` | 0.047 ± 0.006 | 0.036 ± 0.007 | 0.0019 | 0.733 ± 0.042 | 0.726 | 0.620 | 54 ± 6 % |
| C1P, `:pressure` | 0.0052 ± 0.0002 | 0.0014 ± 0.0007 | 0.0019 | 0.888 ± 0.016 | 1.039 | 0.923 | 36 ± 6 % |

Reading against the criteria fixed before the reruns:

- **C1P passes.** The share is 0.0052 by the split and 0.0014 by the ripple, below
  the 0.01 threshold and within the band (0.002 to 0.008) that cannot be told apart
  from the damping-layer theory and the open-column floor. The rise after 8 h that
  C1 showed is gone. The split stays at 0.005 to 0.007 from 5 to 17 h, and the
  ripple at 0.0002 to 0.0027. These results support the interpretation that the
  interface left by the `P_bar` term was the source of the weak-layer excess. That
  excess was 25 times the theory in C1, and by extension probably 10 times in T1
  and 5 times in S1.
- **C0P is between the two outcomes.** The share falls by a factor of 2.2 (split)
  and 2.0 (ripple), to 0.019 and 0.016. That is 16 to 19 times below the continuous
  theory of the strong layer, and about twice the 0.01 level the review associated
  with near-complete absorption. The rise at 9 to 11 h remains, from 0.013 to 0.022,
  on the round-trip time from the 60 km base. So a return from the layer region is
  still present, now without the interface. The data are consistent with the
  discrete `K0 = 0.05` layer reflecting a share of about 0.02, well below the 0.30
  the continuous theory gives. Why the discrete layer differs from the theory is
  still open.

**The cavity closure now holds in its ripple form.** The level predicted from the
ripple's `r_top` at 9.75 km with the lossless lower column matches the measured
level within 9 per cent (C0P, 1.19 against 1.09) and 4 per cent (C1P, 0.92 against
0.89). In C0 and C1 the misfit was 11 and 15 per cent. The form that reads `r0`
from the split at the ground misses by 22 and 17 per cent. Section 8's reading 1
(closure within 10 per cent) is therefore met by the ripple form once the interface
is removed. The level contrast between the two layers (1.09 against 0.89, a ratio of
1.23) is now in the direction of the theory's contrast (0.41 against 0.95). It is
not of the theory's size, because the strong layer returns far less than the theory
says.

**The net-flux peak-to-peak does not come from the interface.** It is 35 to 36 per
cent in both fixed runs, against 45 to 54 per cent before. A steady single mode
should have a net flux that does not depend on height, so the steadiness gate
proposed in the review (below 10 per cent) is not met at 12 to 17 h. The remaining
non-stationary content is not explained by the interface. It is the next candidate
for the start-up transient, and the ramped-forcing pair of Section 8 remains the
test for it.

### What the switch test changes

1. The background drift, the 0.74 per cent mass loss through the bottom face and the
   near-neutral or statically unstable interface at the layer base all come from the
   kinetic energy carried in `P_bar`. `pbar_mode = :pressure` removes all three to
   round-off without a wave, and to below 0.1 per cent with one.
2. Every record made with `:total` and a mean wind carries these features. That
   covers S0 to T2, L1, C0, C1, the h = 100 m and h = 30 m linear-corner records of
   `linear_flux_constancy.md`, and Stage 1 and Stage 2. The returned shares of
   Section 7 mix an interface return with the layer's own. The weak-layer rows are
   dominated by the interface (C1: 0.047 with it, 0.005 without). For the
   strong-layer rows the interface accounts for about half the share (C0: 0.042
   against 0.019).
3. The instability reported on the sponge base in `linear_flux_constancy.md` §18,
   and its following the strong layer in Section 7, are consistent with the
   interface, which in Z0 is statically unstable without any wave.
4. The strong-layer difference from the continuous theory remains, at a factor of
   about 16. It is not due to compressibility (Section 9, adjusted-column theory)
   and not to the interface.
5. `:pressure` should be the setting for all further runs. Whether to make it the
   default is a decision for the codebase. Changing the default would break the
   bit-for-bit reproducibility of the earlier records unless those scripts pin
   `:total`.

Section 11 records the decision taken on this point: `:pressure` is now the default,
and the scripts that regenerate mean-wind records pin `:total`.

## 11. Closing the linear corner

*Note dated 2026-10-02.* Every run in Sections 11 to 13 imposes the ground `w` through the
ghost cells (`bottom_forcing = :ghost`). Every absolute flux, level and amplitude quoted in
these sections is superseded for amplitude: the ghost path loses 10.7 % of the flux in
uniform wind and 39 % on the MCD column at CFL 0.4 (Section 14, addendum 2 and addendum 3
step 1). Ratios between ghost-path runs at equal `dt` keep their meaning as comparisons;
on the MCD column the normalized ratios are also shifted (`R10/30` -15 %, `R50/30` +7 %).

Added on 2026-09-23 after a review of Sections 8 and 9. The review recommended
closing the linear corner with one final batch that informs production, and then
moving to production. This section records the re-reading of the fixed runs, the
decisions taken without new compute, the error budget after the fix, and the final
batch with its criteria fixed before launch.

### Re-reading the fixed runs

The hourly series in `output/reflection_tests/fixed_cavity_closure.json` were
re-read directly, rather than through the window means of Section 9.

**The P_bar diagnosis.** The Z pair meets every criterion fixed before it, with no
fitted parameter. The predicted change in the density gradient was -0.072 %/km and
the measured one is -0.069 %/km. In Z1 the initial column is a fixed point of the
discrete scheme to round-off (largest density departure 3.5e-11 %, `drho_L2` 5e-14
at 6 h), and Z0 reproduces the C0 drift without a wave (mass -0.720 against -0.721 %
at 6 h). Within the Arsia column and this discretization, the Z1 test establishes
that the drift, the bottom-face mass loss and the 60 km interface come from the
kinetic-energy part of `P_bar`, and that no second imbalance in the initial column
remains. The mass-loss path at the bottom face is inferred by elimination rather than
measured, which does not affect this conclusion, since Z1 removes the loss.

**C1P and the estimator floor.** The C1P split reads 0.0067, 0.0063, 0.0066, 0.0067,
0.0067, 0.0070 and 0.0069 at 5 to 11 h. The first return from a 60 km base reaches
9.75 km only at about 8.9 h, so the split carries a floor of about 0.006 that is not
reflection. Over the whole record the ripple stays at or below 0.003, and its ground
phase wanders (-97 to -146 deg at 7 and 12 to 17 h, -53 to +39 deg at 8 and 9 h), so
it reads noise. C1P's return cannot be distinguished from the theory's 0.0019. Split
values below about 0.007 are not quoted as reflection from here on.

**C0P's return is coherent.** The C0P ripple is 0.0018, 0.0056 and 0.0153 at 8, 9
and 10 h, rising on the round-trip time. From 10 to 17 h its ground phase stays
within 16 deg of 180 deg (-172, -166, -175, 173, 176, 175, 164, 171 deg), which
indicates a steady, phase-locked return. The discrete `K0 = 0.05` layer therefore
returns a share of about 0.016 (amplitude ratio about 0.13) where the continuous
theory gives 0.30. Compressibility and the interface are excluded (Section 9), and the
cause is open. The sign of the difference means the solver's layer absorbs more than
the theory says.

**The unsteadiness matches between layers.** The hourly net-flux peak-to-peak is
nearly the same in C0P and C1P:

| t [h] | 3 | 5 | 10 | 11 | 12 | 14 | 15 | 16 | 17 |
|---|---|---|---|---|---|---|---|---|---|
| C0P [%] | 433 | 130 | 71 | 53 | 49 | 35 | 31 | 30 | 32 |
| C1P [%] | 470 | 127 | 72 | 53 | 46 | 34 | 33 | 34 | 30 |

The two layers' returns differ by about ten times (ripple 0.016 against 0.0014), so a
cause tied to the layer would not be expected to give matching series. The match is
consistent with the start-up transient of the impulsive forcing. Two further features
fit that reading. The C0P ripple-closure misfit (prediction over measured level, at
9.75 km) shrinks from +26 % at 12 h to +3 to +9 % at 15 to 17 h as the peak-to-peak
falls, and C1P's level swings between 0.84 and 0.99 with an essentially zero return. The
"closure met" reading of Section 9 is therefore conditional: its 9 % window mean
averages a series that is still converging, and every hourly misfit in C0P has the
same sign.

One caveat on the decay. From 10 to 14 h the peak-to-peak roughly halves, but over
14 to 17 h it levels off at 30 to 35 % in both runs. A decay that continues to 10 %
near 30 h is therefore not demonstrated by these records. A floor near 30 % would
point to a non-stationary component that the impulsive start does not explain. The
ramped pair below is designed to separate the two.

For production, the single-mode control has a net-flux peak-to-peak of about 430 % at
3 h. The 3 h Stage 2 records are therefore expected to lie inside the transient, in
addition to their `P_bar` exposure of up to 2 % of gravity (Section 9).

### Decisions

Taken on 2026-09-23, without new compute.

1. `pbar_mode = :pressure` is the default in `SimConfig` (`src/parameters.jl`). With
   zero mean wind the two modes are bit-identical, because the subtracted kinetic
   energy is exactly zero. The scripts that regenerate records made with a mean wind
   pin `pbar_mode = :total`: `sweep_arsia_mons_stage1.jl`, `sweep_mars_mountain.jl`,
   `run_critical_level_fig7.jl` and `run_mountain_benchmark_strict.jl`.
   `sweep_arsia_mons_stage2.jl` passes the mode explicitly and keeps `pbar=0`
   (`:total`) as its default, so the commands recorded for the earlier records
   reproduce them. All new runs pass `pbar=1`.
2. Every snapshot NetCDF carries the global attributes `pbar_mode`, `K0_sponge` and
   `sponge_frac` (`src/io.jl`), so each file records its provenance.
3. Flux levels are quoted against the band-mean `N` (0.00849 s^-1 over 10 to 55 km),
   with the factor of 2.08 to the `N(0)` reference stated alongside
   (`linear_flux_constancy.md` §9 and §19 item 6).
4. The ripple is the estimator of record for the reflected share. The split is kept
   as a secondary estimate with a floor of about 0.006 to 0.007.
5. The production top is the 60 to 120 km layer with `K0 = 0.005`, which is the
   Srivastava default in `parameters.jl`, and where the solver and the theory agree
   (C1P). `K0 = 0.05` is retired.
6. The error budget below replaces the table of Section 7.

### Error budget after the fix

Rows are separate and not summed. Every level and share number from a `:total`
record with a mean wind is superseded.

| component | size | status |
|---|---|---|
| shape (band drift) | within +-0.07 %/km at 15 h across four layers (Section 7) | measured under `:total`; to be re-checked on `:pressure` records |
| level, from the layer, weak thick layer (60 to 120 km, `K0 = 0.005`) | a few per cent (C1P return indistinguishable from theory) | post-fix |
| level, from the layer, S0-type layer (60 to 90 km, `K0 = 0.05`) | up to about +-25 %, depending on the phase of `|r|` about 0.13 | post-fix; layer retired |
| level, from the transient | net peak-to-peak 30 to 35 % at 14 to 17 h with an impulsive start and 13 to 16 % with a 2 h ramp, 430 % at 3 h (single mode); C1R level 6 to 11 % below the undamped-layer value, still rising at 17 h | dominant quantifiable term; ramped pair below |
| level, from the `N_0` reference | factor 2.08 between `N(0)` and band-mean `N` | quoted alongside (decision 3) |
| level, from the ghost-path forcing (Section 14 addendum) | -10.7 % at CFL 0.4 and -4.9 % at CFL 0.2 against the face path (C2R and S14C against S14G1 and S14G2), whose level is `dt`-independent to 0.03 % | added 2026-10-01; new runs use `bcface=1` |
| level, converged (face path) against the linear Rayleigh-layer theory | +11.9 % (1.019 against 0.911, 12-17 h); the `N(z)` treatment of the near-ground inversion moves the theory flux by -0.4 to +0.3 % (Section 14 addendum 2) | added 2026-10-01; mechanism unattributed |
| level, ghost-path forcing on the MCD column | -39 % at `MF(30.25)` at CFL 0.4 (S14M1/S04 = 1.64, 6-9 h); not uniform in height (1.94 at 10.25 km, 1.53 at 50.25 km), so normalized ratios carry it too (`R10/30` -15 %, `R50/30` +7 %) | added 2026-10-02 (Section 14 addendum 3); new runs use `bcface=1` |
| 2 dz boundary mode | ground amplitude about 0.30 of the smooth `w` (0.21 to 0.35 across paths and grids); penetration about `dz^2 / (2 c dt)` in uniform wind (0.36 to 3.7 km); on the MCD column 0.45 to 0.6 of the smooth `w` through the calm layer | added 2026-10-02 (Section 14 addendum 3); descriptive |
| vertical grid spacing on the MCD column (face path, `dz` = 500 m against 250 m, S14M2 against S14S2, 6-9 h) | +8.3 % relative to 250 m at `MF(30.25)` (-7.7 % change on halving `dz`); -12.5 % change at `MF(10.25)`, -6.0 % at `MF(50.25)`; `R10/30` -5.2 %, `R50/30` +1.8 %; the `MF(30.25)` change grows through the window (-5.8, -7.8, -8.8 % by hour) | added 2026-10-19 (Section 14 addendum 3, S2); applies to every absolute Stage 2 flux at `dz` = 500 m until a 250 m record replaces it |

### Ramped-forcing pair

**Design.** C0P and C1P repeated with the forcing switched on as `tanh(t / t0)`
(`:linear_ramp`, `src/forcing.jl`), `t0 = 2` h, 18 h, `pbar=1`, local disk,
`OPENBLAS_NUM_THREADS=1`, `julia -t 1`. The ramp reaches 0.995 of full amplitude by
6 h, so its tail reaches 55 km before about 11 h and lies outside the 12 to 17 h
window. The ramp's frequency, `1/t0 = 1.4e-4` rad/s, is about 20 % of the intrinsic
`kU = 7.5e-4` rad/s, so the forcing stays close to stationary. The script argument is
`ramp=2`, and the directory carries the tag `_ramp2h`.

| id | command (`sweep_arsia_mons_stage2.jl`) |
|---|---|
| C0R | `linear mode=6 h=0.005 hours=18 pbar=1 ramp=2` |
| C1R | `linear mode=6 h=0.005 hours=18 pbar=1 ramp=2 k0=0.005 lz=120` |

**Criteria fixed before launch** (hour windows centred on 12 to 17 h):

1. The start-up transient is confirmed as the source of the unsteadiness if the net
   peak-to-peak is below 10 % in both runs and C1R's level is within 5 % of 0.951.
2. C0R's ripple share should stay at 0.016 +- 0.004 with a ground phase within about
   15 deg of the C0P phase (about 172 deg), since that return belongs to the layer.
3. If the peak-to-peak stays near 30 %, the unsteadiness is logged as an open
   numerical item, not a blocker for production.

The result sets how Stage 2 is forced and how long its records must be.

**Runs.** Launched together with the production A/B at 16:07. Local disk,
`OPENBLAS_NUM_THREADS=1`, `julia -t 1`.

| id | steps | wall clock [s] | end | directory under `output/` |
|---|---|---|---|---|
| C0R | 87513 | 6334.7 | 17:52:58 | `arsia_mons_linear18h_LX2000_mode6_pbarP_ramp2h/run_U40.0_h0.005` |
| C1R | 87514 | 7650.0 | 18:14:56 | `arsia_mons_linear18h_LX2000_mode6_K0_0.005_LZ120_pbarP_ramp2h/run_U40.0_h0.005` |

The wall clocks are 1.67 and 1.50 times those of C0P and C1P (3791 s, 5089 s),
although the step counts match and the ramp adds only one `tanh` per bottom cell per
step. For the first 25 min the machine ran four busy threads (the A/B used two). That
overlap by itself does not account for a slowdown lasting the whole record, so the
cause is not identified. Logs: `output/reflection_tests/logs/C0R_mode6_pbarP_ramp2h.log`,
`C1R_mode6_K0_0.005_LZ120_pbarP_ramp2h.log`. Closure:
`ramped_cavity_closure.{json,png}` and `logs/ramped_cavity_closure.txt`. Background:
`ramped_background_audit.{json,png}`.

**Background.** Both runs hold the initial column. At 18 h the largest x-mean density
departure is 0.023 % (at 59.75 km, C0R), the change in `d ln rho/dz` over 5 to 50 km
is below 1e-4 %/km, and `min N^2/N^2_0` over 59 to 61 km is 0.993 or higher.
`drho_L2` at 18 h is 1.2e-3 (C0R) and 1.7e-3 (C1R), which is set by the wave.

**Hourly series** (net peak-to-peak in %, with the impulsive runs for comparison):

| t [h] | 3 | 5 | 8 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| C0P, impulsive | 433 | 130 | 61 | 71 | 53 | 49 | 48 | 35 | 31 | 30 | 32 |
| C0R, ramped | 733 | 252 | 56 | 40 | 40 | 31 | 25 | 26 | 21 | 13 | 13 |
| C1P, impulsive | 470 | 127 | 60 | 72 | 53 | 46 | 48 | 34 | 33 | 34 | 30 |
| C1R, ramped | 736 | 244 | 56 | 39 | 42 | 32 | 21 | 23 | 18 | 14 | 16 |

Before about 8 h the ramped flux is still small, so its relative peak-to-peak is
large. Those hours are not comparable with the impulsive runs.

**Means over 12 to 17 h:**

| run | split | ripple | ripple ground phase | level measured | level predicted from ripple `r_top` at 9.75 km | net peak-to-peak, mean of hourly values |
|---|---|---|---|---|---|---|
| C0P | 0.019 | 0.016 | -166 to 176 deg | 1.093 | 1.191 | 38 % |
| C0R | 0.017 | 0.0150 (0.0139 to 0.0162) | -176 to 169 deg | 1.029 (0.954 to 1.103) | 1.206 | 22 % |
| C1P | 0.0052 | 0.0014 | wanders | 0.888 | 0.923 | 38 % |
| C1R | 0.0032 (0.0029 to 0.0035) | 0.0012 (0.0011 to 0.0015) | -113 to -132 deg | 0.868 (0.845 to 0.890) | 0.919 | 21 % |

**Against the criteria.**

1. Not met. Over 12 to 17 h the net peak-to-peak averages 22 % (C0R) and 21 % (C1R).
   It falls below 20 % only at 15 to 17 h, and it is not below 10 % in either run.
   C1R's level is 0.868, which is 8.7 % below 0.951, so the 5 % condition is not met
   either.
2. Met. C0R's ripple share is 0.0150, inside 0.016 +- 0.004. Its ground phase stays
   within 11 deg of 180 deg over 12 to 17 h, as in C0P. The strong layer's return does
   not depend on how the forcing starts.
3. Not the outcome either. The peak-to-peak does not stay near 30 %. It keeps falling
   through 17 h and ends at 13 to 16 %.

**Reading.** The result lies between the pre-registered outcomes. The ramp roughly
halves the unsteadiness at equal times: the 12 to 17 h mean drops from 38 to 21 or
22 % in both layers, and from 30 to 34 % to 13 to 16 % at 16 to 17 h. The decline
continues at the end of the record, whereas the impulsive runs levelled off at 30 to
35 % over 14 to 17 h. The two layers again give nearly the same series, with returns
that differ by more than ten times (ripple 0.015 against 0.0012). These data support
the interpretation that the start-up transient is the main part of the unsteadiness
in the impulsive runs. They do not show that the transient is the whole of it. What
remains at 18 h may be the slower tail of the ramp's own transient, or a separate
non-stationary component. An 18 h record does not separate the two.

Two further observations follow from the ramped pair.

- The split floor in C1R is 0.003, about half the 0.006 of C1P, and it holds at
  0.0018 to 0.0035 through the record. So part of the split floor noted above came
  from the impulsive start. The ripple share (0.0012) agrees with the Rayleigh-layer
  theory (0.0019) within the noise of the estimator, as in C1P.
- The C1R level stays 6 to 11 % below the undamped-layer value of 0.951 and 2 to 8 %
  below the prediction built from the measured `r_top`. Because the return in this
  run is at the theory's value, reflection does not account for the shortfall. The
  level rises from 0.845 at 13 h to 0.89 at 16 to 17 h. The shortfall may therefore
  be partly transient, with a remaining part of unidentified origin. The level row of
  the budget keeps a few per cent for the layer and carries this shortfall as part of
  the transient row until a longer record resolves it.

**Consequences for Stage 2.**

1. Force with a ramp. At equal record length the ramp halves the net-flux
   unsteadiness and does not change the layer's return. A `tanh` ramp on the
   orographic forcing (`MountainParams`) is not yet implemented and is the one code
   change Stage 2 needs.
2. Record length. Quantities that rely on flux steadiness to better than about 15 %
   need records of at least 18 h after a 2 h ramp. The 3 h records remain inside the
   transient.
3. The unsteadiness that remains is logged as an open numerical item. It does not
   block production, but every Stage 2 flux quote should state the net peak-to-peak
   of its own final window.

### Production A/B

**Design.** The Stage 2 main-grid cases U40 h1.0 and U40 h3.0 are repeated with
`pbar=1`, keeping the old 60 to 90 km, `K0 = 0.05` layer, `LX = 4000` km and the 3 h
length, so the `P_bar` fix is the only change. Command: `sweep_arsia_mons_stage2.jl
pbar=1 u=40 h=1,3` (`julia -t 2`, local disk). Comparison against
`output/arsia_mons_stage2/` with `analyze_arsia_stage2.py --compare-dir`, final-hour
mean of `<rho u'w'>_x` at the levels of the Section 7 table.

**Criteria fixed before launch.**

1. The 4.5 km absorption stands if `|MF(10.25)/MF(1.25)|` stays below 0.05 with
   `pbar=1` (0.008 and 0.022 under `:total`) and `MF(1.25)` changes by less than 10 %.
2. The opposite-sign flux in the AMEC band is robust to the fix if its sign at
   50.25 km is kept and its magnitude changes by less than 25 % (the size of the
   S0-type level row). Otherwise it is sensitive to `P_bar`, and the Stage 2 flux
   quotes from `:total` records are withdrawn.
3. `drho_L2` at 3 h is reported for both modes as a measure of the background drift.

**Runs.** One process, `julia -t 2`, beside C0R and C1R, local disk,
`OPENBLAS_NUM_THREADS=1`, started 16:07:21 and ended 16:32:25 (1503 s for the process).

| run | steps | wall clock [s] | directory under `output/` |
|---|---|---|---|
| U40 h1.0, `:pressure` | 15068 | 1444.4 | `arsia_mons_stage2_pbarP/run_U40.0_h1.0` |
| U40 h3.0, `:pressure` | 15509 | 1497.2 | `arsia_mons_stage2_pbarP/run_U40.0_h3.0` |

The `:total` references took 2833 s and 2888 s, as two of nine threads in one process.
Log: `output/reflection_tests/logs/AB_stage2_U40_pbarP.log`. It carries harmless
HDF5-DIAG lines from the existence check before each new snapshot, which appear only
in the threaded process; every snapshot opens. Analysis:
`output/arsia_mons_stage2_pbarP_analysis/` (`mf_compare.{png,json}`,
`background_drift.json`).

**Final-hour flux** (`<rho u'w'>_x` in uPa, every third snapshot counted back from
the last). This estimator gives -361.3 at 1.25 km for the `:total` U40 h1.0 record,
against -356 in the Section 7 table, which took a different subset of the same
snapshots. Both modes are read identically here.

| run | mode | 1.25 km | 2.75 km | 10.25 km | 30.25 km | 50.25 km | 55.25 km | `|MF(10.25)/MF(1.25)|` |
|---|---|---|---|---|---|---|---|---|
| U40 h1.0 | `:total` | -361.3 | -48.87 | +3.085 | +13.57 | +9.056 | +2.697 | 0.009 |
| U40 h1.0 | `:pressure` | -361.1 | -48.81 | +2.868 | +13.48 | +9.098 | +3.726 | 0.008 |
| U40 h3.0 | `:total` | -4515 | -142.8 | -96.46 | +69.58 | +57.25 | +11.40 | 0.021 |
| U40 h3.0 | `:pressure` | -4518 | -139.1 | -97.37 | +69.85 | +59.38 | +17.25 | 0.022 |

**Background at 3 h** (x-mean against the initial column, from the snapshots):

| run | mode | `drho_L2` | density, 50.25 km | density, 53.25 km | `u`, 50.25 km [m/s] | `u`, 53.25 km [m/s] | `u = 0` level near 53 km [km] |
|---|---|---|---|---|---|---|---|
| U40 h1.0 | `:total` | 0.0191 | -4.07 % | -4.31 % | +5.76 | +7.67 | 53.30 to 52.32 |
| U40 h1.0 | `:pressure` | 0.0052 | -0.18 % | -0.18 % | +0.26 | +0.42 | 53.30 to 53.24 |
| U40 h3.0 | `:total` | 0.0264 | -4.80 % | -4.84 % | +7.44 | +10.15 | 53.30 to 52.01 |
| U40 h3.0 | `:pressure` | 0.0160 | -0.89 % | -0.75 % | +1.44 | +2.70 | 53.30 to 52.94 |

**Against the criteria.**

1. Met. The 4.5 km absorption is unchanged: `|MF(10.25)/MF(1.25)|` is 0.008 and
   0.022 with the fix, and `MF(1.25)` changes by 0.06 and 0.07 %.
2. Met in the AMEC band. At 50.25 km the flux keeps its sign and changes by +0.5 %
   (h1.0) and +3.7 % (h3.0). Above the band the fix does matter: at 55.25 km, across
   the 53.3 km critical level, the flux is 38 % (h1.0) and 51 % (h3.0) larger with
   `:pressure`. The profiles in `mf_compare.png` agree below about 53 km and depart
   between 53 and 60 km.
3. Under `:total` the background changes more than the wave can account for. Between
   45 and 60 km the x-mean density falls by 3 to 5 %, and the mean wind at 50 to 55 km
   gains 6 to 10 m/s within 3 h. The `u = 0` level near 53 km drops to 52.5 km within
   the first 0.6 h in both `:total` runs, whatever the wave amplitude. With `:pressure`
   it stays at 53.3 km in h1.0 and descends slowly to 52.9 km in h3.0. That slow
   descent is consistent with deposition of wave momentum below the critical level,
   which the larger wave drives.

**Reading.** For these two cases, the 4.5 km absorption and the sign and size of the
opposite-sign flux inside the 45 to 53 km band hold up under the fix, within 4 %. The
Stage 2 statements that rest on those two quantities therefore survive, with the
3 h transient caveat above. The structure at and above the 53.3 km critical level does
not survive. In the `:total` records the critical level sits about 1 km lower, and
the flux that crosses it is 28 to 34 % smaller, than with the background held. Any
Stage 2 result about the upper edge of the AMEC band, the critical-level position or
the flux reaching 55 km should come from `pbar=1` records.

### Open item, the strong-layer discrepancy

The discrete `K0 = 0.05` layer returns about 0.016 against the continuous theory's
0.30 for mode 6. This is parked as open and non-blocking, since production uses the
weak thick layer. The isolated test is specified here and is run only if a reviewer
asks for it or a configuration needs `K0 = 0.05`: a single `k` in an isothermal
column with a uniform wind, `K0 = 0.005`, 0.0158 and 0.05 on a 60 to 90 km layer, with
the Rayleigh-layer theory recomputed for that column. An isothermal column removes the
`N(z)` structure, so any remaining difference would lie in the discrete layer or the
eq19 top.

## 12. Long Stage 2 records

*Note dated 2026-10-02.* Ghost-path records: absolute fluxes superseded for amplitude, see
the note at the head of Section 11 and Section 14, addendum 2.

Started 2026-09-23 as the first production step after Section 11. The nine cases of
the Stage 2 main grid (`U_bar` of 20, 40 and 57 m/s by `h_max` of 0.1, 1 and 3 km, real
reversing MCD wind, periodic `LX = 4000` km) are rerun with the settings Section 11
arrived at:

- `pbar_mode = :pressure` (`pbar=1`);
- the 60 to 120 km damping layer with `K0 = 0.005` (`lz=120 k0=0.005`);
- the orographic forcing switched on as `tanh(t / 2 h)` (`ramp=2`), through the new
  `MountainParams.t_ramp` field (`src/forcing.jl`). With `t_ramp = 0` the forcing path
  is unchanged, and a 0.1 h rerun of the legacy U40 h1.0 command reproduces
  `output/arsia_mons_stage2/run_U40.0_h1.0` bit for bit (6 snapshots). With the ramp,
  the bottom-row `w` scales as `tanh(t / 7200 s)` (0.0099 against 0.0100 at 72 s);
- 24 h records (`hours=24`). The ramped single-mode control was still settling at
  17 h (Section 11), and 24 h allows a final window to be compared with an early one
  from the same record.

Command, one process per case (`julia -t 1`, `OPENBLAS_NUM_THREADS=1`, local disk):
`sweep_arsia_mons_stage2.jl hours=24 pbar=1 ramp=2 k0=0.005 lz=120 u=<U> h=<h>`.
Output: `arsia_mons_stage2_24h_K0_0.005_LZ120_pbarP_ramp2h/run_U<U>_h<h>`.

These records differ from the 3 h records of Section 7 in four respects (`P_bar`,
layer, top, ramp) and in length. They are therefore not a controlled comparison with
them. The A/B of Section 11 isolates the `P_bar` effect, and the ramped pair isolates
the ramp.

**Caveat on the periodic domain.** Stationary waves are unaffected by the periodic
boundary. Non-stationary content that advects with the mean wind crosses 4000 km in
about 19.5 h at 57 m/s (the jet speed at 20 to 30 km) and in 28 h at 40 m/s. So in the
final window of a 24 h record, part of the non-stationary content in the jet may be
wrapped. The ramp keeps that content small, but it is not zero.

**Windows.** Hour means of `<rho u'w'>_x` from every snapshot. The early window is 6 to
9 h (hours centred on 6.5, 7.5 and 8.5 h), after the ramp has reached 0.995. The
final window is 18 to 24 h (six hours centred on 18.5 to 23.5 h).

**Criteria fixed before launch.**

1. Steadiness. For each run, report the peak-to-peak of the hour-mean `MF(z)` over 15
   to 45 km, relative to its mean over that band, in each hour. Between the 4.5 and
   53.3 km critical levels a steady, non-dissipative wave has a height-independent
   flux. The target is below 10 % in the final window. Following Section 11, a run
   that does not reach it is quoted with its own value.
2. Whether the opposite-sign flux above the 4.5 km critical level is stationary. Take
   `MF(30.25) / |MF(1.25)|` in each window. It is read as a stationary feature if the
   final-window value is within a factor of 1.33 of the early-window value. It is read
   as transient if it falls below 0.25 of the early value. Anything in between is
   reported as partly transient.
3. Absorption at the 53.3 km critical level. Report `MF(55.25) / MF(50.25)` in the
   final window. It is descriptive, and it is compared with the 3 h A/B (0.41 and 0.29
   for U40 h1.0 and h3.0 with `:pressure`).
4. Background integrity. In the three `h = 0.1` km runs, the x-mean `u = 0` level near
   53 km stays within 0.1 km of 53.30 km through 24 h. The largest x-mean density
   departure below 45 km stays below 0.1 %. In the `h = 1` and 3 km runs, the descent
   of the critical level is reported as wave-driven. The 3 h A/B gave 0.36 km for U40
   h3.0 and 0.06 km for U40 h1.0.
5. Numerical health. Every run reaches 24 h. `min_P` and `min_rho` stay positive, and
   the nonlinear runs (`N h / U` above 1) are reported with their largest `|u'|` and the
   fraction of points below 60 km with `d theta/dz < 0`. The ratio `|u'|/|u|` is not
   used, since `u` passes through zero at the critical levels. (This wording was
   settled before any output was read.)

**Runs.** Launched 2026-09-23 at 18:39, nine processes in parallel (the machine has
24 cores, 16 of them performance cores), local disk (`~/MAGNUSP_runs`), `julia -t 1`,
`OPENBLAS_NUM_THREADS=1`. Logs: `S2L_U<U>_h<h>.log`, with start and end stamps in
`S2L.start`. The first minute gave a projected wall clock of 8.0 to 8.6 h per run.
That is about twice what the 90 km top would need, because the time step halves.
`compute_adaptive_dt` (`src/euler.jl`) takes `CFL min(dx, dz) / max(|u| + c)`, and the
MCD zonal wind is -81 m/s at 90 km, -191 m/s at 100 km and -327 m/s at 120 km. So
the 90 to 120 km part of the column sets `dt` at about 0.37 s, against about 0.72 s
with the 90 km top. A dimension-split limit, `min(CFL dx / (|u| + c), CFL dz / (|w| + c))`,
would allow about 1 s here. That is a candidate optimization for later records. It
changes the numerics, so it was not applied to this batch.

The batch ended at 03:53:59 on 2026-09-24. Analysis:
`python3 scripts/analyze_arsia_stage2.py --sweep-dir <batch> --timeseries --outdir
output/arsia_mons_stage2_24h_analysis`, which produced `timeseries.{json,png}`.

| run | `N h / U` | steps | reached [h] | wall clock [s] | outcome |
|---|---|---|---|---|---|
| U20 h0.1 | 0.09 | 223412 | 24.0 | 32976.9 | complete |
| U40 h0.1 | 0.04 | 223415 | 24.0 | 31701.1 | complete |
| U57 h0.1 | 0.03 | 223416 | 24.0 | 33124.1 | complete |
| U20 h1.0 | 0.90 | 223450 | 24.0 | 31769.9 | complete |
| U40 h1.0 | 0.45 | 223636 | 24.0 | 31666.1 | complete |
| U57 h1.0 | 0.32 | 122374 | 12.93 | 20310.9 | aborted, timestep collapse near 42 km |
| U20 h3.0 | 2.70 | 125702 | 13.29 | 20660.9 | aborted, timestep collapse at 100 to 120 km |
| U40 h3.0 | 1.35 | 80038 | 8.29 | 14706.4 | aborted, timestep collapse at 100 to 120 km |
| U57 h3.0 | 0.95 | 55726 | 5.49 | 13195.1 | aborted, timestep collapse at 100 to 120 km |

The directory is `output/arsia_mons_stage2_24h_K0_0.005_LZ120_pbarP_ramp2h/`. Logs:
`output/reflection_tests/logs/S2L_U<U>_h<h>.log`.

### Where the four runs failed

The failure sites were located from the last snapshots of each record (largest `|w|`,
`|u'|` and `|T/T_0 - 1|`).

- **At the model top** (U57 h3.0, U40 h3.0, U20 h3.0). `|w|` in the top cells rises to
  30 to 160 m/s over the one to three hours before the collapse, which then starts
  at 100 to 120 km. In this part of the column the MCD zonal wind is supersonic.
  `|u|/c` is 0.47 at 90 km, 1.13 at 100 km and 1.67 to 1.73 at 110 to 120 km, where
  `u = -289` to `-327` m/s. The damping layer (`K0 = 0.005`) therefore overlaps a
  supersonic background. That configuration was not part of the linear-corner
  validation, which used a uniform 40 m/s wind. The disturbance above 100 km scales
  with the wave amplitude. The largest `|w|` there (snapshots every 25 outputs) is
  0.1 m/s or less in the `h = 0.1` km runs, 0.5 to 20.6 m/s in the `h = 1` km runs,
  and 72 m/s or more in the `h = 3` km runs before their collapse. That is consistent
  with the collapse being driven by wave amplitude reaching the supersonic layer.
  U57 h1.0, which collapsed elsewhere, reached 20.6 m/s at the top at 10.7 h.
- **In the free column** (U57 h1.0). Before the collapse, `|w|` of 14 to 18 m/s and
  `|u'|` of 50 to 65 m/s sit at 37 to 48 km, below the 53 km critical level, which by
  then has descended to 52.1 km. The collapse starts near 42 km. U40 h1.0 carries the
  same structure (`|w|` 12 to 17 m/s at 41 to 45 km, `|u'|` up to 78 m/s) and survives
  to 24 h, with up to 2.9 % of points below 60 km statically unstable.

In all four aborted runs, 1.6 to 4.7 % of points below 60 km have `d theta/dz < 0` in
the last hour. These runs are consistent with nonlinear wave breaking that the scheme
does not survive. The present data do not separate a physical breaking event from a
numerical instability that the breaking triggers.

Note added after Section 13. U57 h1.0 completes 24 h with the 95 km top at mean `dt`
0.66 s and collapses at 15.65 h at `dt` 0.29 s. The wave amplitudes and fluxes of this
batch depend on the time step (Section 13, "CFL test results"), so its fluxes are
specific to `dt` of about 0.38 s.

### The completed runs against the criteria

Final window 18 to 24 h. Early window 6 to 9 h.

| run | MF peak-to-peak 15-45 km, early, final | `MF(30.25)` early, final [uPa] | `MF(1.25)` early, final [uPa] | `MF(30.25)/|MF(1.25)|`, final over early | `MF(55.25)/MF(50.25)`, final | `u = 0` level near 53 km at 24 h [km] | largest x-mean density departure below 45 km |
|---|---|---|---|---|---|---|---|
| U20 h0.1 | 31 %, 7.5 % | 0.045, 0.083 | -0.80, -0.56 | 2.68 | 0.33 | 53.29 | 0.042 % (1.25 km) |
| U40 h0.1 | 30 %, 6.6 % | 0.18, 0.33 | -5.7, -12.2 | 0.85 | 0.31 | 53.27 | 0.156 % (1.25 km) |
| U57 h0.1 | 29 %, 6.1 % | 0.37, 0.67 | -15.5, -36.0 | 0.79 | 0.29 | 53.24 | 0.309 % (1.25 km) |
| U20 h1.0 | 16 %, 24 % | 4.3, 16.4 | -604, -180 | 16.8 | -0.03 | 52.41 | 2.4 % |
| U40 h1.0 | 17 %, 33 % | 17.9, 78.5 | -2545, +589 | 38.9 | 0.00 | 51.16 | 4.8 % |

1. **Steadiness.** Met by the three `h = 0.1` km runs as a window mean (6 to 7.5 %).
   The hourly values are 10.7 to 11.5 % at 18.5 to 19.5 h and 2 to 9 % from 20.5 h.
   The drop at about 20 h occurs in all three runs, whatever the forcing wind. It
   coincides with the time the real -57 m/s jet at 20 to 30 km takes to cross the
   4000 km domain (19.5 h). So the low final values may owe something to the periodic
   wrap (see the caveat above). Not met by the two `h = 1` km runs (24 and 33 %),
   whose unsteadiness grows from the early window.
2. **Whether the opposite-sign flux is stationary.** As written, the criterion reads
   U40 h0.1 (0.85) and U57 h0.1 (0.79) as stationary and U20 h0.1 (2.68) as neither.
   The normalization hides a common trend, though. `MF(30.25)` itself grows by a
   factor of 1.8 to 1.85 between the windows in all three linear runs, and
   `MF(1.25)` also changes (by 0.7 to 2.3 times), because the flux below the 4.5 km
   critical level responds to a slow change of the near-ground mean state. In the
   linear runs the opposite-sign flux above 4.5 km therefore persists, and it grows
   rather than decays over 24 h. This does not support the Section 7 reading that it
   is a transient of the 3 h records. Its growth over 6 to 24 h is consistent with
   slow establishment of the long-wavelength part of the ridge spectrum, as in
   `linear_flux_constancy.md` §10. In the `h = 1` km runs, the ratio is dominated by
   the change of `MF(1.25)`, which in U40 h1.0 reverses sign.
3. **Absorption at the 53.3 km level.** In the linear runs 29 to 33 % of the 50.25 km
   flux reaches 55.25 km in the final window, declining from about 0.35 to 0.23
   through the window. In the `h = 1` km runs the flux at 55.25 km is indistinguishable
   from zero (`-0.03 +- 0.04` of the 50.25 km value). The critical level has descended
   to 52.4 km (U20) and 51.2 km (U40), about 0.09 km per hour in U40 h1.0 over the
   second half of the record.
4. **Background integrity.** The `u = 0` level holds within 0.06 km in all three
   `h = 0.1` km runs (met). The density condition (below 0.1 %) is met by U20 h0.1
   only. In U40 h0.1 and U57 h0.1 the departure sits at 0.25 to 1.25 km and grows
   steadily with time. At fixed `h` it scales close to `U_bar^2` (0.042 : 0.156 :
   0.309 against 1 : 4 : 8.1). That is consistent with a second-order mean-flow
   response to the wave absorbed at the 4.5 km critical level, not with a drift of the
   solver. In the `h = 1` km runs this response is large. In U40 h1.0 the near-ground
   x-mean wind develops additional zero crossings at 1.4 to 4 km from about 10 h, and
   the density at 1.25 km falls by 4.8 %.
5. **Numerical health.** Not met. Four of the nine runs aborted: all three at
   `h = 3` km, plus U57 h1.0.

### Reading

For the linear Stage 2 cases (`h = 0.1` km), the configuration holds the background,
and in the final window the flux between the critical levels is steady to within
about 10 %. The periodic wrap noted above is the main qualification. In these cases
the opposite-sign flux above 4.5 km persists over 24 h, and about 30 % of the flux in
the AMEC band crosses the 53.3 km critical level. For `h = 1` km, the wave drives a
large, progressive change of the mean state near the ground and at the 53 km critical
level over 24 h, so these records have no steady state to quote. For `h = 3` km and
for U57 h1.0, the records end in a collapse before the final window. Three of those
collapses begin in the part of the damping layer where the MCD background wind is
supersonic.

Open decisions before the nonlinear Stage 2 cases can be quoted:

- The treatment of the column above about 95 km, where the background wind is
  supersonic. Possibilities are a lower model top (95 km, `|u|/c` 0.78), a background
  wind limited inside the damping layer, or a stronger damping rate near the top.
  None of these has been validated in the linear corner.
- Whether the free-column collapse near 42 km (U57 h1.0) is numerical. A rerun at a
  smaller CFL number would test that.

## 13. A 95 km top for Stage 2

*Note dated 2026-10-02.* Ghost-path records: absolute fluxes and levels superseded for
amplitude, see the note at the head of Section 11 and Section 14, addendum 2.

Decided 2026-09-24 after Section 12. The top moves to 95 km, with the damping layer at
60 to 95 km and `K0 = 0.005` (`lz=95 k0=0.005`). At 95 km the MCD wind is -135 m/s,
with `|u|/c` of 0.78, so the supersonic part of the column (above about 98 km) is left
out. The time step returns to about 0.65 s. The layer is 35 km thick instead of 60 km,
so it is expected to return more than the 60 to 120 km layer. It is checked in the
linear corner first.

**Theory.** Same Rayleigh-layer calculation as the single-mode predictions of Section 8
(mode 6, U = 40 m/s, MCD column, radiating top; script in the session scratchpad,
reproducing the recorded 0.0019 and 0.951 for 60 to 120 km and 0.301 and 0.411 for the
60 to 90 km, `K0 = 0.05` layer):

| layer | band `R` | open-column band `R` | upward level against the open column |
|---|---|---|---|
| 60 to 120 km, `K0 = 0.005` | 0.0019 | 0.0084 | 0.951 |
| 60 to 95 km, `K0 = 0.005` | 0.0112 | 0.0025 | 0.911 |

### Step 1, single-mode control (C2R)

`sweep_arsia_mons_stage2.jl linear mode=6 h=0.005 hours=18 pbar=1 ramp=2 k0=0.005 lz=95`,
which is C1R with the top at 95 km.

**Criteria fixed before launch** (hour windows centred on 12 to 17 h):

1. The ripple share is 0.02 or less. That is well below the retired layer (0.016 to
   0.03 measured, 0.30 theory), and it is the condition for production use. Agreement
   with the theory's 0.011 within a factor of two is reported alongside.
2. The background holds. The largest x-mean density departure is below 0.1 %, and
   `min N^2/N^2_0` over 59 to 61 km is 0.99 or higher at 18 h.
3. The level ratio C2R / C1R is within 10 % of the theory's 0.911 / 0.951 = 0.958.
   Comparing like with like cancels the shared shortfall of C1R (Section 11).

If criterion 1 or 2 fails, the Stage 2 batch is not launched on this top.

### Step 1 results

**Run.** Launched 2026-09-24 09:48:19 (`julia -t 1`, `OPENBLAS_NUM_THREADS=1`, local
disk), next to the CFL test. 87 514 steps to 18 h, `dt = 0.74 s`, `drho_L2` 1.56e-3 at
the end (1.53e-3 at 17 h), exit 0. Wall clock 11 532.2 s (`real` 11 546 s).

- Output: `output/arsia_mons_linear18h_LX2000_mode6_K0_0.005_LZ95_pbarP_ramp2h/run_U40.0_h0.005`
  (written to `~/MAGNUSP_runs`, 3.1 GB, copied afterwards)
- Log: `output/reflection_tests/logs/C2R_mode6_K0_0.005_LZ95_pbarP_ramp2h.log`
- Analysis: `output/reflection_tests/c2r_cavity_closure.{json,png}` and
  `c2r_background_audit.{json,png}`, both with C1R alongside
  (`compare_sponge_reflection.py --cavity C1R C2R`, `--background C1R=… C2R=…`)

The single-mode corner has a uniform wind, so the time step does not depend on the top,
and C1R took the same 87 514 steps. C2R cost 0.132 s per step against 0.087 s for C1R
(7 650.0 s), although it has 190 levels instead of 240. The C2R process ran at normal
priority (`PRI 31`, nice 0). This matches the unexplained 1.5 to 1.7 fold slowdown of
C0R and C1R against C0P and C1P (Section 11), and it remains open. The CFL test ran
at about 0.21 s per step on the Stage 2 grid, compared with about 0.14 s per step in
the nine-process batch of Section 12.

**Window means over 12 to 17 h** (six hourly windows):

| quantity | C1R (60 to 120 km) | C2R (60 to 95 km) | theory C2R |
|---|---|---|---|
| ripple share | 0.0012 (0.0011 to 0.0015) | 0.0019 (0.0016 to 0.0025) | 0.0112 |
| split share, band | 0.0032 | 0.0039 (0.0037 to 0.0042) | 0.0112 |
| `R_top` at 9.75 km | 0.0010 | 0.0035 (0.0024 to 0.0046) | 0.0114 |
| level | 0.868 (0.845 to 0.890) | 0.910 (0.895 to 0.923) | 0.911 |
| net peak-to-peak | 20.7 % | 20.7 % (15.1 to 32.1 %) | |

**Background at 18 h (C2R).** Largest x-mean density departure 0.026 % (61.75 km),
temperature 0.023 %, wind 0.13 m/s, `N^2` change up to +0.64 % (60.75 km),
`min N^2/N^2_0` over 59 to 61 km 1.0005. The C1R values are 0.018 %, 0.017 %,
0.14 m/s, +0.30 % and 1.0000.

**Scoring against the criteria.**

1. Met. The ripple share of 0.0019 is an order of magnitude below the 0.02 bound.
   It is a factor of 5.8 below the theory's 0.011, so it does not agree within a
   factor of two. The split share (factor 2.9 below) and `R_top` at 9.75 km (factor
   3.3 below) are closer. All three rise from C1R to C2R, in the direction the theory
   predicts, but by 1.6, 1.2 and 3.5 times instead of the predicted 5.9 and 5.0 times.
2. Met. Density departure 0.026 % against 0.1 %, and `min N^2/N^2_0` 1.0005 against 0.99.
3. Met at the edge, with the opposite sign. The ratio C2R / C1R is 0.910 / 0.868 = 1.048,
   which is 9.4 % above the 0.958 target and inside the 10 % band. The theory predicts
   a lower level for the thinner layer. The level instead rises by 4.8 %, and C2R
   matches the absolute theory (0.910 against 0.911), whereas C1R sits 8.7 % below its
   own (0.868 against 0.951). So the shortfall is not shared as assumed. One reading
   consistent with this is that the C1R deficit comes from the 95 to 120 km part of
   its column rather than from a transient common to both runs. This is not tested
   here. The criterion is not a gate.

Criteria 1 and 2 are met, so the Stage 2 batch was launched. Within the linear corner,
the 95 km top returns less than the Rayleigh-layer theory predicts. The measured shares
are at or below the split-estimator floor of Section 11 (0.006 to 0.007 for C1P, about
0.003 for C1R), so the difference from theory is not resolved at this amplitude. For
production use, the bound that matters (criterion 1) holds with margin.

### Step 2, Stage 2 batch on the 95 km top

If step 1 passes, all nine main-grid cases are rerun for 24 h:
`hours=24 pbar=1 ramp=2 k0=0.005 lz=95 u=<U> h=<h>`. That includes the five that
completed with the 120 km top, so the grid shares one configuration, and the
`h = 0.1` and 1 km pairs give a sensitivity to the top for free. The criteria of
Section 12 apply unchanged, with one addition. For the five cases that completed in
Section 12, the final-window `MF(50.25)` and `MF(55.25)/MF(50.25)` are compared
between the two tops. A difference below 10 % means the top does not matter for those
quantities.

**Launch.** 2026-09-24 13:03:44, nine processes (`julia -t 1` each, 2 s stagger,
`OPENBLAS_NUM_THREADS=1`), next to the running CFL test. A dry run confirmed `nz = 190`,
layer 60 to 95 km, `K0 = 0.005`, `pbar_mode` pressure, ramp `tanh(t/2 h)`, CFL 0.40.

- Output: `~/MAGNUSP_runs/arsia_mons_stage2_24h_K0_0.005_LZ95_pbarP_ramp2h/run_U<U>_h<h>`
- Logs: `~/MAGNUSP_runs/logs/S2T95_U<U>_h<h>.log`
- Start and end stamps: `~/MAGNUSP_runs/logs/S2T95.start`

Both are moved to `output/` once the batch finishes.

A throughput check 12 minutes after launch gave 0.27 to 0.29 s per step, which
projected to about 10 h per run. That estimate was taken over the first 2 000 steps,
which include start-up. Over whole runs the batch averaged 0.127 to 0.134 s per step,
against 0.142 s per step for the 120 km batch of Section 12 (240 levels against 190),
so the batch ran at the expected rate. The 1.5 fold per-step cost of C2R against C1R
(Step 1) is measured over whole runs and remains unexplained.

### CFL test for the free-column collapse

`hours=16 pbar=1 ramp=2 k0=0.005 lz=95 cfl=0.2 u=57 h=1`, launched together with
step 1. Its twin at CFL 0.4 is U57 h1.0 in the step 2 batch. Reading:

- If the CFL 0.4 twin collapses near 40 km and the CFL 0.2 run passes the same time
  plus 3 h, the collapse depends on the time step and is numerical.
- If both collapse near the same time and height, the collapse is consistent with a
  breaking event that the scheme does not represent at either time step.
- If neither collapses, the collapse of Section 12 depended on the 120 km top.

### Step 2 results

The batch ended at 18:09:13 on 2026-09-24. Analysis:
`python3 scripts/analyze_arsia_stage2.py --sweep-dir <batch> --timeseries --outdir
output/arsia_mons_stage2_24h_LZ95_analysis`. The comparison between tops is in
`top_comparison.txt` in the same directory.

- Output: `output/arsia_mons_stage2_24h_K0_0.005_LZ95_pbarP_ramp2h/` (74 GB)
- Logs: `output/reflection_tests/logs/S2T95_U<U>_h<h>.log` and `S2T95.start`

| run | `N h / U` | steps | reached [h] | mean `dt` [s] | wall clock [s] | outcome |
|---|---|---|---|---|---|---|
| U20 h0.1 | 0.09 | 131069 | 24.0 | 0.659 | 17528.0 | complete |
| U40 h0.1 | 0.04 | 131072 | 24.0 | 0.659 | 17403.8 | complete |
| U57 h0.1 | 0.03 | 131074 | 24.0 | 0.659 | 17466.2 | complete |
| U20 h1.0 | 0.90 | 131107 | 24.0 | 0.659 | 16635.6 | complete |
| U40 h1.0 | 0.45 | 131186 | 24.0 | 0.659 | 16744.3 | complete |
| U57 h1.0 | 0.32 | 131285 | 24.0 | 0.658 | 16725.8 | complete |
| U20 h3.0 | 2.70 | 131332 | 24.0 | 0.658 | 16700.6 | complete |
| U40 h3.0 | 1.35 | 72313 | 10.11 | 0.503 | 18276.0 | aborted, timestep collapse at 92 to 95 km |
| U57 h3.0 | 0.95 | 46118 | 6.82 | 0.532 | 8840.9 | aborted, timestep collapse at 94 to 95 km |

Seven of nine runs reach 24 h, against five with the 120 km top. The two aborted runs
fail in the top cells. `|w|` above 85 km rises to 34 m/s (U40 h3.0, 9.65 h) and to 61
to 146 m/s (U57 h3.0, 6.4 to 6.6 h) before the collapse, which starts at 92 to 95 km.
At the same stage, near the ground `|u'|` reaches 115 to 190 m/s and `|T/T_0 - 1|`
reaches 3.7 to 9. The top is subsonic in these runs, so the supersonic background of
Section 12 is not required for a failure at the top at `h = 3` km.

**Criteria of Section 12** (early window 6 to 9 h, final window 18 to 24 h):

| run | MF peak-to-peak 15-45 km, early, final | `MF(30.25)` early, final [uPa] | `MF(1.25)` early, final [uPa] | `MF(30.25)/|MF(1.25)|`, final over early | `MF(55.25)/MF(50.25)`, final | `u = 0` level near 53 km at 24 h [km] | largest x-mean density departure below 45 km |
|---|---|---|---|---|---|---|---|
| U20 h0.1 | 28 %, 8.1 % | 0.037, 0.056 | -0.28, -0.22 | 1.99 | 0.34 | 53.30 | 0.031 % |
| U40 h0.1 | 26 %, 7.8 % | 0.15, 0.21 | -1.9, -3.9 | 0.72 | 0.33 | 53.28 | 0.116 % |
| U57 h0.1 | 25 %, 7.4 % | 0.30, 0.41 | -5.0, -14.1 | 0.49 | 0.31 | 53.25 | 0.223 % |
| U20 h1.0 | 15 %, 14 % | 3.2, 13.2 | -218, -33 | 101 | -0.01 | 52.65 | 1.8 % |
| U40 h1.0 | 11 %, 23 % | 10.4, 33.6 | -952, +628 | 8.5 | -0.01 | 51.85 | 3.9 % |
| U57 h1.0 | 18 %, 35 % | 23.7, 60.4 | -2161, +666 | 10.2 | -0.02 | 50.93 | 5.3 % |
| U20 h3.0 | 19 %, 31 % | 24.9, 71.1 | -2511, +494 | 21.9 | -0.10 | 50.62 | 5.6 % |

1. Steadiness. Met by the three `h = 0.1` km runs (7.4 to 8.1 %). As with the 120 km
   top, the hourly values fall from 12 to 14 % at 18.5 to 19.5 h to 3 to 7.5 % from
   20.5 h, at the jet's wrap time. Not met by the `h = 1` and 3 km runs (14 to 35 %).
2. Opposite-sign flux. By the ratio as written, U40 h0.1 (0.72) and U57 h0.1 (0.49)
   are partly transient and U20 h0.1 (1.99) is neither. As in Section 12, `MF(30.25)`
   grows between the windows in all three linear runs (by 1.4 to 1.5 times here,
   against 1.8 there) while `MF(1.25)` changes independently, so the ratio mixes the
   two trends.
3. Absorption at 53.3 km. 31 to 34 % of the 50.25 km flux reaches 55.25 km in the
   linear runs (29 to 33 % with the 120 km top). In the `h >= 1` km runs the value is
   indistinguishable from zero or negative (-0.10 to -0.01).
4. Background integrity. The `u = 0` level holds within 0.05 km in the `h = 0.1` km
   runs (met). The density departure is 0.031, 0.116 and 0.223 % (below 0.1 % for U20
   only), close to the `U_bar^2` scaling seen before (1 : 3.7 : 7.2 against 1 : 4 :
   8.1). In the `h >= 1` km runs the critical level descends to 52.65 to 50.62 km.
5. Numerical health. Not met. Two of nine runs abort, both at `h = 3` km.

**Sensitivity to the top.** The pre-specified comparison fails. For the five cases that
completed with both tops, final-window `MF(50.25)` is 16 to 43 % lower with the 95 km
top (U20 h0.1 -30 %, U40 h0.1 -34 %, U57 h0.1 -37 %, U20 h1.0 -16 %, U40 h1.0 -43 %).
`MF(55.25)/MF(50.25)` differs by +4 to +10 % in the linear runs. The difference is
already present in the early window (-16 to -17 % in the linear runs), and in the
linear runs it extends to `MF(30.25)` (-32 to -38 %) and to `MF(1.25)` (factor 2.5 to
3.2). The two batches differ in the top and in the time step (0.380 s against 0.659 s
mean `dt`). A top whose measured return is about 0.002 (Step 1) is not expected to
change the flux below the 4.5 km critical level by a factor of 2.5, so the time step
is examined next.

### CFL test results

The CFL 0.2 run aborted at 15.65 h (step 192 871, mean `dt` 0.292 s) after 39 672.6 s
of wall clock. The CFL 0.4 twin (mean `dt` 0.658 s) completed 24 h. None of the three
pre-specified readings covers this outcome, because it is the run with the smaller time
step that collapses.

- Output: `output/arsia_mons_stage2_16h_K0_0.005_LZ95_pbarP_ramp2h_cfl0.2/run_U57.0_h1.0`
- Log: `output/reflection_tests/logs/S2T95_U57_h1_cfl0.2.log`
- Analysis: `output/arsia_mons_stage2_24h_LZ95_analysis/cfl0.2/`

The two runs separate long before the collapse. At matched times:

| t [h] | free column `|w|` 20-55 km, `dt` 0.66 / 0.38 / 0.29 s [m/s] | `|w|` above 85 km [m/s] | points below 60 km with `d theta/dz < 0` |
|---|---|---|---|
| 6 | 1.5 / 2.0 / 2.2 | 0.1 / 1.4 / 0.7 | 0.1 / 0.2 / 0.2 % |
| 9 | 8.2 / 13.1 / 8.5 | 0.5 / 10.9 / 5.9 | 0.7 / 1.2 / 1.1 % |
| 12 | 9.6 / 17.6 / 23.2 | 0.4 / 8.8 / 10.9 | 1.4 / 2.5 / 3.2 % |
| 15 | 14.9 / (aborted 12.9) / 20.7 | 0.6 / - / 8.8 | 1.9 / - / 6.8 % |

The middle column is U57 h1.0 with the 120 km top at CFL 0.4. The ordering by `dt`
also holds in the momentum flux (hourly means, uPa):

| t [h] | `MF(1.25)`, `dt` 0.38 / 0.29 / 0.66 s | `MF(30.25)` | `MF(50.25)` |
|---|---|---|---|
| 2 | -456 / -516 / -270 | 3.85 / 3.96 / 3.34 | 1.75 / 1.77 / 1.57 |
| 4 | -2435 / -2807 / -1263 | 11.8 / 12.4 / 10.1 | 8.50 / 8.73 / 7.46 |
| 6 | -3207 / -4926 / -1726 | 26.8 / 26.9 / 14.4 | 20.8 / 22.3 / 11.6 |
| 8 | -6227 / -7733 / -2476 | 46.1 / 49.7 / 34.4 | 25.9 / 24.4 / 23.8 |

Up to 8 h, the 95 km top at `dt` 0.29 s follows the 120 km top at `dt` 0.38 s within
about 8 % in `MF(30.25)` and `MF(50.25)`. The 95 km top at `dt` 0.66 s is 8 to 46 %
lower. This pattern is consistent with the difference between the batches coming from
the time step rather than from the top. Across the three runs, a smaller time step goes
with a larger free-column amplitude at 6 and 12 h (not at 9 h, where the 120 km run
leads), a larger unstable fraction at 12 h, and a collapse
(at 12.9 h and 15.65 h for `dt` 0.38 and 0.29 s, none by 24 h for 0.66 s). Of the
Section 12 collapses, the U57 h1.0 free-column collapse is therefore more consistent
with the smaller `dt` of the 120 km batch than with the supersonic top.

The mechanism is not identified. The bottom forcing is imposed through ghost cells
each step without an explicit `dt`, so the dependence enters through the Lax-Wendroff
update, the sponge source term or the split implicit diffusion. The linear corner of
Sections 8 to 11 ran at one time step (0.74 s) throughout, so it does not constrain
this. The near-ground flux `MF(1.25)`, which differs by the largest factor, sits
between the surface, where the MCD wind is zero, and the 4.5 km critical level. That
region is expected to be sensitive to any dissipation, numerical or physical.

### Open decisions after Section 13

1. Time-step convergence. The Stage 2 fluxes differ by 16 to 43 % between mean time
   steps of 0.38 and 0.66 s. Until this is resolved, Stage 2 amplitudes should be quoted
   with their time step. A pre-specified test is proposed: U40 h0.1 on the 95 km top at
   CFL 0.4, 0.2 and 0.1 for 9 h (early window only). Convergence would mean `MF(30.25)`
   and `MF(1.25)` over 6 to 9 h change by less than 5 % between CFL 0.2 and 0.1. The
   same test in the uniform-wind linear corner (C2R at CFL 0.2) would show whether the
   dependence needs the real wind profile.
2. The `h = 3` km cases. They still fail at the top with a subsonic top, after strongly
   nonlinear near-ground behaviour (`|T/T_0 - 1|` above 3). With the time-step
   dependence open, it is not settled whether these failures are physical breaking.

## 14. Time-step convergence and the launch layer (pre-registration)

Recorded 2026-09-30, after Section 13 and before any Section 14 run. Every numeric
reading in this section was fixed before the suite launched. The runs are launched by
`scripts/run_dt_convergence_S2T95.sh` and scored by `scripts/analyze_dt_convergence.py`.
The driver changes they need (`dx=`, `dz=`, `outint=`, `ublend=`) are in
`scripts/sweep_arsia_mons_stage2.jl` behind defaults that reproduce every earlier record,
which the launcher checks bit for bit before anything else starts. The values from the
MCD column quoted for checks Z4 and Z5 below were computed on 2026-09-29 and 2026-09-30
from the background profile alone, before the suite, and they replace the estimates of
the first draft, which used a synthetic Arsia-shaped column.

### What Section 13 established and what it did not

Section 13 separated the top from the time step by accident. The 95 km batch ran at a
mean `dt` of 0.659 s, the 120 km batch of Section 12 at 0.380 s, and the U57 h1.0 CFL
test at 0.292 s. In the three linear (`h = 0.1` km) cases, final-window `MF(30.25)` is
32–38 % lower at 0.659 s than at 0.380 s (18–20 % lower in the 6–9 h window) and
`MF(1.25)` differs by a factor 2.5–3.2, while the measured top return is about 0.002
(C2R). Up to 8 h the U57 h1.0 record at 0.29 s follows the 120 km record at 0.38 s
within about 8 %, and the 0.66 s record is 8–46 % lower. The ordering is monotone:
smaller `dt`, larger flux, larger free-column `|w|`, earlier collapse.

Two facts from the same tables constrain where the dependence lives. First, most of the
effect is common to all levels above 10 km. Across the linear runs the absolute fluxes
change by 30–42 % between the two batches (final window), while `MF(50.25)/MF(30.25)`
changes by +2.1 to +2.3 % (+2.6 to +2.8 % at 6–9 h) and `MF(55.25)/MF(50.25)` by +4 to
+10 %. `MF(10.25)/MF(30.25)` changes by -2.2 to -6.0 % (-4.0 to -5.9 % at 6–9 h), so a
smaller part of the effect, of a few per cent, accumulates between 10 and 30 km. The
difference in absolute flux roughly doubles from the early to the final window, so it
grows with time. Second, the quantity that moves most, `MF(1.25)`, sits between a surface
where the MCD wind is zero and the 4.5 km crossing, in a layer where `|u|` is at most
0.85 m/s, and it is negative (-1.9 and -3.9 uPa for U40 h0.1, early and final) while
`MF(30.25)` is positive (0.15 and 0.21 uPa).

### The launch-layer hypothesis H_L and the Lax-Wendroff dissipation

The forcing imposes `w(x, 0) = U_bar dh/dx` through ghost cells with `U_bar` = 20, 40 or
57 m/s, but the background it is imposed into is the MCD column, whose wind is zero at
the ground and stays below 1 m/s up to the 4.5 km crossing. The ridge is
`h(x) = h exp(-(x - x_0)^2 / (2 sigma^2))` with `sigma = 100` km (the driver default; the
project baseline of 50 km is not what Stage 2 ran), so the largest imposed velocity is
`U_bar h / (sigma sqrt(e))`, 0.024 m/s for U40 h0.1 and 0.73 m/s for U40 h3.0. The
linearized buoyancy equation for a stationary disturbance,

```
d b'/dt + U d b'/dx + N^2 w' = 0 ,
```

has the steady solution `b' = -N^2 zeta` with `U d zeta/dx = w'` only where `U` is not
zero. Where `U` vanishes and `w'` is imposed, the displacement grows secularly,
`zeta = w' t`, about 88 m per hour for U40 h0.1 and 2.6 km per hour for U40 h3.0. In the
inviscid linear problem nothing stops this growth inside the calm layer, so the amplitude
that the layer holds, and the flux that escapes above 4.5 km, are set by whatever
dissipation acts there. In the solver the only dissipation below 60 km besides molecular
diffusion (`nu` below 1 m^2/s, negligible at these scales) is the truncation dissipation
of the Lax-Wendroff sweeps, and that dissipation depends on `dt`.

The two-step Lax-Wendroff sweep for `q_t + A q_z = 0` reduces, for constant `A`, to

```
q_j^{n+1} = q_j - (dt / 2 dz) A (q_{j+1} - q_{j-1}) + (dt^2 / 2 dz^2) A^2 (q_{j+1} - 2 q_j + q_{j-1}) .
```

The second-order term uses the compact operator `D+D-` where a consistent square of the
centred first derivative would be `D0^2`. Their symbols differ by
`-(4/dz^2) sin^4(k_z dz / 2)`, so each step removes `(A^2 dt^2 dz^2 / 8) k_z^4` of each
Fourier component to leading order, and per unit simulated time the scheme adds a
hyperdiffusion

```
kappa_4 = c^2 dt dz^2 (1 - nu_z^2) / 8 ,   nu_z = c dt / dz ,
```

on the eigencomponents that travel at the acoustic speed `c` in the vertical sweep, which
are the pair (`w`, `p'/(rho_0 c)`). The factor `(1 - nu_z^2)` is 0.84–0.99 at the CFL
numbers used here. The dissipation per unit time is proportional to `dt` at fixed `dz`:
halving the time step halves the numerical damping that the wave experiences over the
same simulated hours. For a hydrostatic stationary mode (`m = N/|U|`, `p' = -rho_0 U u'`),
the acoustic pair carries the fractions `k^2/m^2` (in `w`) and `U^2/c^2` (in `p'`) of the
energy, so the amplitude damping rate is

```
gamma = (dt dz^2 m^4 / 16) (c^2 k^2 / m^2 + U^2) ,
```

and dividing by `c_gz = k U^2 / N` gives the flux attenuation per metre

```
2 gamma / c_gz = dt dz^2 N^5 (1 + c^2 k^2 / N^2) / (8 k U^4) .
```

Integrating toward a critical level with `|U| = |U'| s` gives optical depth one at

```
U* = [ dt dz^2 N^5 (1 + c^2 k^2 / N^2) / (24 k |U'|) ]^(1/3) .
```

Where `|U| < U*` the wave is controlled by numerical dissipation, and `U*` scales as
`dt^(1/3) dz^(2/3)`. The expression is a WKB estimate. It fails within about
`|U|/|U'| < 1/N` of a critical level and where `Ri = N^2/U'^2` is order one, so it is an
order-of-magnitude band, not an edge. On the MCD column (check Z4) at `dt` = 0.659 s and
`dz` = 500 m, the median `U*` below 4.5 km is 11–24 m/s over modes n = 10 to 1 of
`LX = 4000` km (9–19 m/s at `dt` = 0.33 s), against `|u|` of at most 0.85 m/s, so the
controlled fraction of the launch layer is 1.00 for every mode. At the 53.3 km crossing
`U*` is 1.8–3.7 m/s (1.4–2.9 m/s at 0.33 s), a band 245–510 m deep, at most about one
`dz`. The launch-layer values are about five times the synthetic-column estimate of the
first draft (2.5–5 m/s), because `U*` scales as `N^(5/3)` and the MCD near-ground
inversion (`N` of about 0.018 s^-1 at 0.25 km) lies inside the launch layer. That
contrast is the core of H_L: at the resolved critical levels the numerics control only
a sliver about one grid cell thick, whose flux is removed by the physical critical-level
absorption anyway, while the entire calm launch layer is numerics-controlled.

A path version of the same mechanism, applying the hyperdiffusion to the wave in the
sheared free column between 4.5 and 53 km, gives about 0.1 % flux attenuation in the C2R
geometry and a few per cent in the Stage 2 column. A naive model that applies `c^2` to
all components without the energy projection gives about 5 %. Neither can produce the
30–42 % change in absolute flux with a +2 % change in `MF(50.25)/MF(30.25)`, which is why
H_L places the effect in the launch layer. The -2 to -6 % change in
`MF(10.25)/MF(30.25)` is of the size the path version allows, so it is not by itself
evidence against H_L, but it means P2 below is not a formality.

H_L accounts for the other Section 13 regularities without extra assumptions. The largest
x-mean density departure below 45 km scales close to `U_bar^2` (1 : 3.7 : 7.2 against
1 : 4 : 8.1), as expected for a layer pumped by a forcing of amplitude proportional to
`U_bar`. Near-ground `|T/T_0 - 1|` of 3.7–9 in the h = 3 km runs is consistent with the
secular displacement reaching kilometres in a layer only 4.5 km deep. The sign of the
`dt` effect (smaller `dt`, less damping, more flux, earlier collapse) is the sign of the
collapse ordering in the CFL test. The legacy width-sweep figures show the same
structure: `MF` near the ground reaches about -3e-3 Pa against about 2e-5 Pa in the free
column, and the horizontally uniform `|w|` of about 1 m/s at 70–90 km in those figures is
the `k = 0` `pbar :total` adjustment of Section 9, not a wave. The legacy "deposited
fraction" of `stage2_response.md` is not usable (its denominator `MF(10 km)` changes sign,
and `MF` is a domain mean that must be multiplied by `LX` before widths are compared).

### Competing hypotheses and their predictions

| id | mechanism | where `dt` enters | predicted size of the `dt` effect | distinct prediction |
|---|---|---|---|---|
| H_L | calm launch layer with imposed `w`, amplitude set by LW hyperdiffusion | `kappa_4 = c^2 dt dz^2 / 8` below 4.5 km | tens of per cent in absolute `MF`, common to all levels above 4.5 km | normalized ratios invariant; removing the calm layer (`ublend`) removes the effect; halving `dz` at fixed `dt` changes absolute `MF` by 10 % or more |
| H_P | LW hyperdiffusion along the path and at the 53.3 km crossing | same operator, free column | about 0.1 % (C2R) to a few % | trend in `MF(z)` concentrated at 45–53 km |
| H_S | commutator of dimensional (x then z) splitting in shear | `dt U'/2` | about 2e-3, below 1 % | persists with `ublend`; present in the sheared column only |
| H_K | discretization of the sponge source inside the LW z-sweep | `K_0 dt / 2` | about 0.16 % | confined to the layer and the level ratio; C2R `dt` pair shows it |
| H_D | split implicit diffusion | `nu dt / dz^2` | negligible below 60 km (`nu` < 1 m^2/s) | none measurable here |

### Zero-cost checks on existing records (run before the suite)

These read records already on disk and are run first. They cannot change the suite,
only the reading of it. Outputs go to `output/reflection_tests/s14/`.

- Z1, profile ratio. `analyze_dt_convergence.py zprofile --a output/arsia_mons_stage2_24h_K0_0.005_LZ95_pbarP_ramp2h/run_U40.0_h0.1 --b output/arsia_mons_stage2_24h_K0_0.005_LZ120_pbarP_ramp2h/run_U40.0_h0.1 --t0 18 --t1 24`
  (and `--t0 6 --t1 9`). H_L predicts that the ratio `MF_0.66 / MF_0.38` is flat above
  10.25 km, spread below 10 % over 10.25–50.25 km, with its departure from one set below
  10.25 km. H_P predicts a trend concentrated at 45–53 km. The pair also differs in the
  top, whose measured return (about 0.002) is too small to matter at this reading.
- Z2, growth law. `zdensity` on the three LZ95 `h = 0.1` records. Secular pumping of the
  calm layer predicts an exponent near 2 in `t` for the largest x-mean density departure
  below 45 km over 2–9 h; a steady wave-driven mean-flow response predicts near 1.
- Z3, `h = 3` km precursors. `zh3` on the LZ95 U40 h3.0 and U57 h3.0 records. H_L
  predicts that ground `|T/T_0 - 1|` exceeds 0.5 before top `|w|` exceeds 5 m/s, and that
  the top variance is mostly `k = 0` or outside the ridge band n ≤ 12. Top reflection of
  the ridge wave predicts the top crossing first with more than half of the top variance
  in n ≤ 12. Free-column breaking predicts the 20–60 km unstable fraction crossing 1 %
  first.
- Z4, `U*` profile. `zustar --dt 0.659 --dz 500 --lx-km 4000` and again with `--dt 0.33`.
  H_L requires `|u| < U*` over most of 0–4.5 km and a controlled depth at 53.3 km no
  larger than one `dz`. Result (quoted above): met, controlled fraction 1.00 for every
  mode, depth at 53.3 km 194–510 m.
- Z5, blend gate. `zblend --zb-km 8` (then 10, 12 if needed), described with the suite
  below. Result: passes at 8 km (minimum `Ri` 1.23 at 3.99 km, no sign change; 2.12 at
  10 km and 2.17 at 12 km).

### The suite

All cases use the 95 km top, `K0 = 0.005`, `pbar = 1` (pressure mode), the 2 h tanh
ramp, periodic `x`, `LX = 4000` km, `sigma = 100` km, and the Section 13 grid unless
stated. `outint` keeps the snapshot cadence at about 66 s in every case
(`outint = 100 x 0.4 / cfl`, or 200 for the `dz = 250` m cases whose `dt` is halved by the
grid). Wall-clock budgets use the measured 0.13 s per step for the 190-level grid
(Section 13, where nine concurrent processes ran at 0.127 to 0.134 s per step) and twice
that for 380 levels, one process per case with `julia -t 1` and
`OPENBLAS_NUM_THREADS=1`. The CFL 0.4 member S04 is the existing 24 h record (window
6–9 h) and is not rerun; the launcher's `check` mode first reproduces 0.1 h of it, and
0.1 h of C2R, bit for bit with the patched driver.

| id | purpose | driver arguments (after `julia -t 1 --project=. scripts/sweep_arsia_mons_stage2.jl`) | mean `dt` [s] | steps | wall clock | output directory |
|---|---|---|---|---|---|---|
| S04 | CFL 0.4 member (existing) | `hours=24 pbar=1 ramp=2 k0=0.005 lz=95 u=40 h=0.1` | 0.659 | 131 072 | done | `arsia_mons_stage2_24h_K0_0.005_LZ95_pbarP_ramp2h` |
| S14A | CFL 0.2 | `hours=9 pbar=1 ramp=2 k0=0.005 lz=95 u=40 h=0.1 cfl=0.2 outint=200` | 0.330 | 98 330 | 3.6 h | `..._9h_K0_0.005_LZ95_pbarP_ramp2h_cfl0.2` |
| S14B | CFL 0.1 | `hours=9 pbar=1 ramp=2 k0=0.005 lz=95 u=40 h=0.1 cfl=0.1 outint=400` | 0.165 | 196 660 | 7.1 h | `..._9h_K0_0.005_LZ95_pbarP_ramp2h_cfl0.1` |
| S14C | C2R at CFL 0.2 (uniform wind) | `linear mode=6 h=0.005 hours=18 pbar=1 ramp=2 k0=0.005 lz=95 cfl=0.2 outint=200` | 0.37 | about 175 000 | 6.4 h | `arsia_mons_linear18h_LX2000_mode6_K0_0.005_LZ95_pbarP_ramp2h_cfl0.2` |
| S14D | `dz = 250` m at the `dt` of S14A | `hours=9 pbar=1 ramp=2 k0=0.005 lz=95 u=40 h=0.1 dz=250 outint=200` | 0.33 | about 98 000 | 7.9 h | `..._9h_K0_0.005_LZ95_pbarP_ramp2h_dz250` |
| S14E | launch-consistent wind, CFL 0.4 (`ublend` 8 km, from Z5) | `hours=9 pbar=1 ramp=2 k0=0.005 lz=95 u=-40 h=0.1 ublend=8` | 0.659 | 49 165 | 1.8 h | `..._9h_K0_0.005_LZ95_pbarP_ramp2h_ub8` |
| S14F | launch-consistent wind, CFL 0.2 | `hours=9 pbar=1 ramp=2 k0=0.005 lz=95 u=-40 h=0.1 ublend=8 cfl=0.2 outint=200` | 0.330 | 98 330 | 3.6 h | `..._9h_K0_0.005_LZ95_pbarP_ramp2h_ub8_cfl0.2` |
| S14Z | `h = 0` equilibrium regression on the new grid | `hours=1 pbar=1 ramp=2 k0=0.005 lz=95 u=40 h=0 dz=250 outint=200` | 0.33 | about 10 900 | 0.9 h | `arsia_mons_stage2_1h_K0_0.005_LZ95_pbarP_ramp2h_dz250` |

The run directory inside each sweep directory is `run_U40.0_h0.1` (S14E and S14F:
`run_U-40.0_h0.1`; S14C: `run_U40.0_h0.005`; S14Z: `run_U40.0_h0.0`). Tag order in the
sweep name is `..._ramp2h_dx<N>_dz<N>_ub<Z>_cfl<X>`, each present only when it differs
from the default. Total cost is about 31 CPU hours; launched together, the suite
finishes in about 8 h (S14D and S14B set the length), and it needs about 25 GB of local
disk before archiving.

The sequence is `scripts/run_dt_convergence_S2T95.sh dry`, then `check`, then `launch`,
`status` while running, and `archive` at the end.

- `dry` must show `nz = 190` (380 for `dz = 250`), the 60–95 km layer, `K0 = 0.005`,
  pressure `pbar` mode and the tanh ramp.
- `check` compares both 0.1 h reruns with the archived records snapshot by snapshot
  (`analyze_dt_convergence.py bitcheck`). Both must print a maximum difference of exactly
  zero, and `launch` refuses to start unless `check` passed with the current driver and
  `src/`.
- `launch` applies the Z5 rule itself and withdraws S14E and S14F if no height passes.
- `archive` copies to `output/`, compares checksums and prints the `rm` command instead of
  running it.

Scoring: `python3 scripts/analyze_dt_convergence.py suite --run S04=<S04 dir>/run_U40.0_h0.1 --run S14A=... --run S14B=... --run S14C=... --run S14D=... --run S14E=... --run S14F=... --run C2R=output/arsia_mons_linear18h_LX2000_mode6_K0_0.005_LZ95_pbarP_ramp2h/run_U40.0_h0.005`,
which writes `output/reflection_tests/s14/s14_suite.json`. The C2R level and ripple share
of S14C are scored with `compare_sponge_reflection.py --cavity C2R S14C` as in Section 13.

Two design choices need stating. The `ublend` wind is `u(z) = u_L (1 - S) + u_MCD(z) S`
with `S = (1 - cos(pi z / Z)) / 2` below the blend height `Z` and
`u_L = sign(u_MCD(Z)) |U_bar|`. The MCD wind above the 4.5 km crossing is westward
(-10.2 m/s at 8 km), so `u_L = -40` m/s and the case is run with `u=-40`, which makes the
imposed `w = U_bar dh/dx` the lower boundary condition of the same ridge in the flow it
is launched into; the driver refuses a blend whose ground wind and forcing wind differ
in sign. The blend height is fixed by a rule that uses only the background column, never
a wave result: `Z` is the smallest of 8, 10 and 12 km for which the blended wind has no
sign change below `Z` and a minimum Richardson number `N^2/(du/dz)^2` of at least 1
below `Z` (check Z5, `analyze_dt_convergence.py zblend`). The launcher applies this rule
before starting anything. If no height passes, S14E and S14F are withdrawn and P5 is not
read. On the MCD column the rule gives `Z = 8` km (minimum `Ri` 1.23). On a synthetic
column with a uniform 1.5 K/km lapse rate the same blend sat at `Ri` 0.99, so the gate
is not a formality; the MCD near-ground inversion raises it. Once chosen, `Z` is not
varied to improve agreement. The launch-consistent case is a different configuration
from the MCD column, not a correction of it: it removes the 4.5 km critical level and
launches into westward flow, so its fluxes are not compared in absolute value with S04,
only with each other (S14E against S14F). S14D pairs with S14A at the same `dt` and
differs in `dz`; `nu_z` then differs (0.4 against 0.2), which changes the `(1 - nu_z^2)`
factor by 12 %, against the factor 4 from `dz^2`.

### Pre-registered readings

All windows are 6–9 h for the Stage 2 cases and 12–17 h for C2R and S14C. Relative
changes are `(b - a)/|a|`.

- P1, the primary criterion as proposed in Section 13. `MF(30.25)` and `MF(1.25)` change
  by less than 5 % between S14A and S14B. H_L predicts FAIL, with `MF(30.25)` increasing
  by +15 % or more from CFL 0.2 to 0.1. The two extrapolations of the existing 0.38/0.66 s
  pair bracket this: a power law `MF ∝ dt^(-0.82)` gives about +77 %, an exponential in
  `dt` with rate 1.54 s^-1 gives about +29 %. Those two forms are fitted to the 18–24 h
  window. Fitted to the 6–9 h window, where P1 is read (`MF(30.25)` 18–20 % lower at
  0.659 s), the same forms give `dt^(-0.37 to -0.40)`, about +29 to +32 %, and a rate of
  0.73–0.78 s^-1, about +13 to +14 %. The +15 % threshold therefore sits at the lower end
  of the bracket for this window, and a result between +5 and +15 % is recorded as
  neither a pass of P1 nor the change H_L predicts. A pass of P1 falsifies H_L as the
  dominant mechanism at these time steps.
- P2, normalized propagation between S14A and S14B. `MF(10.25)/MF(30.25)` and
  `MF(50.25)/MF(30.25)` change by less than 5 %, and `MF(55.25)/MF(50.25)` by less than
  10 %. H_L predicts PASS. A fail places a `dt` dependence in the free column (H_P or H_S)
  at a level that matters.
- P3, order. With `G1 = MF(S14A)/MF(S04)` and `G2 = MF(S14B)/MF(S14A)`,
  `L = ln G2 / ln G1` near 1 means no finite `dt -> 0` limit over this range, near 0.5
  first-order convergence (Richardson extrapolation admissible), near 0.25 second order.
  Under H_L with a secularly growing launch layer, `L` of 0.7 or more is expected, and
  `L` above 0.8 means extrapolation is not admissible.
- P4, intrinsic solver dependence. C2R against S14C, `MF(10.25)`, `MF(30.25)` and
  `MF(50.25)`: H_L predicts changes below 1 % and a level and ripple share unchanged
  within the Section 11 estimator floor. A change above 2 % in any of them means a `dt`
  dependence intrinsic to the solver in uniform wind (H_K or an unidentified term), and
  it takes precedence over every other reading.
- P5, launch-consistent pair. S14E against S14F: `MF(10.25)` and `MF(30.25)` change by
  less than 5 %. H_L predicts PASS; any mechanism acting above 8 km (H_S, H_P) predicts
  the same change as in the S14A/S14B pair scaled by the `dt` ratio, which is a FAIL.
- P6, resolution at fixed `dt`. S14A against S14D: H_L predicts absolute `MF(30.25)` and
  `MF(1.25)` change by 10 % or more (the hyperdiffusion falls by a factor of about 3.5)
  with `MF(50.25)/MF(30.25)` and `MF(10.25)/MF(30.25)` within 5 %. If the normalized
  ratios also change by more than 5 %, the Stage 2 path is not resolved at `dz = 500` m
  and a resolution pair becomes a requirement for every Stage 2 amplitude, independent
  of H_L.
- Regression. S14Z must satisfy the standing equilibrium test on the `dz = 250` m grid
  (`drho_L2` and energy drift at or below the Section 11 values for the 500 m grid) before
  S14D is read.

### The h = 3 km collapses

The two candidate causes posed for these runs make distinct predictions, and the
existing records already bear on both. Physical breaking that the scheme cannot carry
would start where the wave overturns, in the free column below the critical levels, and
move outward. It did not: the collapses begin in the top cells at 92–95 km, far from any
overturning region, and in the `h ≥ 1` km runs `MF(55.25)/MF(50.25)` is
indistinguishable from zero or negative (-0.10 to -0.01), so there is no stationary
ridge-wave flux to overturn above 55 km. Spurious amplification or reflection of the
ridge wave at the top is disfavoured by the same numbers, since no stationary flux
reaches the damping layer, and by C2R, whose measured return (0.0019) is a factor of 10
below the production bound. What reaches the top is therefore something other than the
ridge wave. A third candidate, which H_L supplies, is consistent with the Section 13
observations: the calm launch layer in the `h = 3` km runs develops `|u'|` of 115–190 m/s
and `|T/T_0 - 1|` of 3.7–9, and a disturbance that strong radiates grid-scale and
acoustic energy (including `k = 0` adjustment) that is not blocked by the stationary
critical levels and grows as `rho^(-1/2)`, a factor of order 100 in amplitude from the
ground to 95 km, until it violates the CFL bound in the top cells. Check Z3
discriminates the three by which threshold is crossed first and by the spectral content
at the top.

The strategy follows from that. Explicit subgrid eddy viscosity (Smagorinsky or an added
hyperdiffusion) is out of scope for this project, and it would also not answer the
question: it would replace a `dt`-dependent numerical damping with a tunable one in
exactly the layer where H_L says the amplitude is set by damping, and the resulting flux
would be a function of the closure coefficient. Localized flux caps are a
parameterization of saturation, not a resolved diagnosis, and they are excluded for the
same reason. Eckermann et al. (2011) scale-dependent radiative damping is in-scope
physics (Newtonian or Eckermann-type damping is allowed), and it belongs in every Stage 2
run once the coefficients are available in `src/damping.jl` (currently missing,
`enable_cooling = false`). It is not a rescue for `h = 3` km: its time scales at 0–95 km
on Mars are of order hours to a sol for the relevant vertical scales, which cannot arrest
a displacement growing at 2.6 km per hour. The `h = 3` km cases are therefore classified
as outside the verified envelope, with no fluxes quoted, until two conditions hold: the
launch configuration passes P5 or an equivalent consistent-launch test, and the
`h = 1` km cases pass their own `dt` convergence test (Section 15, to be written after
these readings). They are then rerun on the resolution pair (`dx = 2` km, `dz = 250` m),
and the report gives only the pre-collapse window, the time and altitude of first
overturning (unstable area fraction and `|u'|` against `|c - U|`), and whether the
collapse site moves with resolution. Breaking is diagnosed, not parameterized.

### Decisions that follow the readings

If P1 fails, P2 and P5 pass and P4 is below 1 %, H_L is supported. The launch-consistent
background (`ublend`) or another launch fix is then proposed for production, subject to
approval, since it changes the configuration. The existing Stage 2 fluxes are re-quoted
normalized at 10.25 km (`MF(z)/MF(10.25)`), which P2 would show to be `dt`-robust, with
the absolute amplitude marked as launch-layer dependent in the error budget. If P1
passes, the standard configuration is converged at CFL 0.2, production moves to CFL 0.2,
and the S04/S14A/S14B sequence gives a Richardson extrapolation row in the budget if P3
gives `L` near 0.5. If P4 exceeds 2 % and P5 fails, the dependence is intrinsic to the
solver and the next test is switch-level: sponge moved out of the LW sweep into a split
exponential step, and the splitting order reversed, each as an on/off pair on C2R. If P6
changes the normalized ratios by more than 5 %, a resolution pair is required for all
Stage 2 amplitudes regardless of the other outcomes. Whichever branch holds, the reading
is recorded against these predictions, and a failed prediction is recorded as failed.

A convergence section for the `h = 1` km cases (Section 15) is written only after these
readings, with its own predictions, because its design depends on which branch holds.

### Launch record

- `dry` (2026-09-30): every case resolved as tabled, `nz` 190 or 380, layer 60–95 km,
  `K0 = 0.005`, pressure mode, `tanh(t/2 h)` ramp, and the eight sweep names match the
  table. Two refusal checks behave as intended: `ublend=8` with `u=40` stops with the
  sign message, and `ublend` in `linear` mode stops before any output.
- `check` (68 s wall clock): the 0.1 h rerun of S04 matches the archived record in all
  6 shared snapshots (steps 0 to 500) and the 0.1 h rerun of C2R in all 5 (steps 0 to
  400), `max |diff|` exactly 0.0 over `rho`, `u`, `w`, `P`, `T` and time. Driver and
  `src/` hash `bb935d8f94affe5f`. Logs: `S14check_S04.log`, `S14check_C2R.log`.
- `launch`: 2026-09-30 00:34:14, all seven cases (S14A to F and S14Z), Z5 blend height
  8 km, `julia -t 1` each with `OPENBLAS_NUM_THREADS=1`, output under `~/MAGNUSP_runs`,
  logs `~/MAGNUSP_runs/logs/S14<id>.log`, stamps `S14.start`. Z1 had finished and Z2 and
  Z3 were still reading archived records at launch; their inputs do not depend on the
  suite.
- S14Z finished at 01:09 (exit 0, 10 974 steps, `dt = 0.328 s`, 2073.9 s of solver wall
  clock, 2082.5 s with start-up).

**How the S14Z regression is read.** The criterion refers to "the Section 11 values for
the 500 m grid". Those values come from Z1 of Section 10, which uses the linear
configuration with a uniform wind (`drho_L2` 5e-14 at 6 h). S14Z runs the main grid with
the sheared MCD column and viscosity on, a configuration that has no h = 0 record at
500 m. To keep the comparison matched, the same command without `dz=250` was run after
S14Z had finished and before any S14D output was read (S14Z0, 01:10, 5 462 steps,
`dt = 0.659 s`, 524.0 s wall clock, log `S14Z0.log`, output
`~/MAGNUSP_runs/arsia_mons_stage2_1h_K0_0.005_LZ95_pbarP_ramp2h/run_U40.0_h0.0`). The
criterion is read as S14Z at or below S14Z0.

| run | `dz` [m] | `drho_L2` peak (time) | `drho_L2` mean, 0.5 to 1 h | `drho_L2` at end | E and KE over 1 h |
|---|---|---|---|---|---|
| S14Z0 | 500 | 4.26e-5 (330 s) | 7.8e-6 | 6.1e-6 | constant to the 6 printed digits |
| S14Z | 250 | 1.13e-5 (492 s) | 2.0e-6 | 1.7e-6 | constant to the 6 printed digits |

S14Z meets the criterion: its departure from the initial column is a factor of about 4
below the 500 m grid at every summary point, and neither run shows a drift in energy at
the precision of the log. The factor of about 4 is consistent with a second-order
truncation error in the initial column, reduced by `(250/500)^2`. In both runs
`drho_L2` rises over the first 5 to 8 min and then decreases, so it measures an initial
adjustment of the sheared, viscous column rather than a growing imbalance. At these
levels (1e-6 to 1e-5) it is two to three orders of magnitude below the wave-driven
`drho_L2` of the Stage 2 records (1e-3). S14D can be read.

### Zero-cost check results

Run 2026-09-30 on archived records (`output/reflection_tests/s14/`, `zchecks.log` with
the JSON outputs `z1_profile_*.json`, `z3_U*_h3.json`, `z4_ustar_dt*.json`; Z1 to Z3 took
916 s of wall clock, reading over the network mount). None of the inputs depends on the
suite.

**Z1, profile ratio `MF(0.659 s) / MF(0.380 s)`, U40 h0.1.** Consistent with H_L in
both windows.

| window | ratio at 1.25 km | at 10.25 km | at 30.25 km | at 50.25 km | median above 10.25 km | spread 10.25–50.25 km |
|---|---|---|---|---|---|---|
| 18–24 h | 0.315 | 0.623 | 0.648 | 0.662 | 0.648 | 6.0 % |
| 6–9 h | 0.329 | 0.769 | 0.810 | 0.832 | 0.810 | 7.7 % |

The spread above 10.25 km is below the 10 % bound, and the departure from one is set
below 10.25 km. The small trend above 10 km is gradual and sits mostly at 10–30 km
(+4.1 and +5.3 %), with +0.8 and +1.0 % over 45–50 km, so it is not the concentration
at 45–53 km that H_P predicts. Between about 4 and 8 km the ratio is irregular (0.37 to
1.29), because the flux passes through small values near the 4.5 km crossing; it is not
read there.

**Z2, growth of the largest x-mean density departure below 45 km (LZ95 `h = 0.1` km).**
The pre-registered 2–9 h exponents are 2.22, 2.20 and 2.16 (U20, U40, U57), near 2 as
secular pumping predicts. Over 9–24 h they fall to 1.27–1.31. The 2–9 h window overlaps
the tanh ramp (0.96 of full forcing at 4 h), which by itself steepens growth, so a
supplementary fit over 6–12 h, after the ramp, was added after the registered reading.
It gives 2.31, 2.19 and 2.07, so the early exponent near 2 is not a product of the ramp.
The slower growth after about 9 h is not attributed here; limitation of the layer
amplitude by the damping that H_L invokes and the jet wrap time (19.5 h) are both
candidates.

**Z3, `h = 3` km precursors (LZ95).** Consistent with H_L on both runs.

| run | ground `|T/T_0 - 1|` > 0.5 | top `|w|` > 5 m/s | 20–60 km unstable fraction > 1 % | top variance at the top crossing (within 0.25 h), `k = 0` / ridge band n ≤ 12 |
|---|---|---|---|---|
| U40 h3.0 | 3.84 h | 5.45 h | 6.64 h | 0.09 / 0.39 |
| U57 h3.0 | 2.82 h | 4.56 h | 5.55 h | 0.07 / 0.24 |

The ground crossing comes first in both runs, 1.6 to 1.7 h before the top crossing, and
the free-column instability comes last. At the top crossing most of the top variance is
outside the ridge band (0.52 to 0.69 in n > 12), and `k = 0` carries 7 to 9 %. Top
reflection of the ridge wave would require the top to cross first with more than half
of the variance in n ≤ 12, and free-column breaking would require the unstable fraction
to cross first; neither is observed. Between the top crossing and the collapse the ridge
fraction varies from 0.01 to 0.89 between snapshots, so the spectral reading is made at
the crossing only. The ground `|T/T_0 - 1|` saturates at 9.00 in the last hour of both
runs, which is consistent with the density floor of `src/run.jl` (`rho >= 0.1 rho_bar`)
acting there, so its value after saturation is not physical.

**Z4 and Z5.** As quoted in the pre-registration: the launch layer is fully
numerics-controlled at both time steps (controlled fraction 1.00), the controlled depth
at 53.3 km is 194–510 m, and the blend gate passes at 8 km (minimum `Ri` 1.23).

Taken together, the five checks agree with the predictions of H_L and do not show the
signatures of H_P (Z1), of top reflection or of free-column breaking (Z3). They do not
test H_L against H_S or H_K, which the suite separates (P4, P5).

### Suite results

All seven cases ran to completion (exit 0, no positivity or timestep aborts). The
suite ran from 00:34:14 to 05:37:09 on 2026-09-30 (18 175 s, 5.05 h of wall clock) and
used 19.4 CPU hours (sum of `real` times, S14Z0 not included), against the 8 h and 31
CPU hours budgeted. The 190-level cases ran at 0.092 to 0.096 s per step, below the
0.13 s of the nine-process Section 13 batch; the 380-level case at 0.184 s per step,
twice that. Local disk before archiving was 25 GB.

| id | steps | mean `dt` [s] | solver wall clock [s] | `real` [s] | s per step | size |
|---|---|---|---|---|---|---|
| S14A | 98 305 | 0.3296 | 9 223.8 | 9 233.0 | 0.094 | 3.5 GB |
| S14B | 196 609 | 0.1648 | 18 143.5 | 18 152.5 | 0.092 | 3.5 GB |
| S14C | 175 029 | 0.3702 | 8 305.7 | 8 314.9 | 0.047 (500 columns) | 3.1 GB |
| S14D | 98 767 | 0.3280 | 18 160.3 | 18 169.5 | 0.184 | 7.0 GB |
| S14E | 49 153 | 0.6592 | 4 705.6 | 4 715.3 | 0.096 | 3.5 GB |
| S14F | 98 307 | 0.3296 | 9 195.6 | 9 204.8 | 0.094 | 3.5 GB |
| S14Z | 10 974 | 0.3280 | 2 073.9 | 2 082.5 | 0.189 | 0.8 GB |
| S14Z0 | 5 462 | 0.6592 | 524.0 | 531.8 | 0.096 | 0.4 GB |

Logs: `output/reflection_tests/logs/S14{A..F,Z,Z0}.log` and `S14.start` (stamps). Outputs
are archived under `output/` with the sweep names of the suite table, plus
`output/arsia_mons_stage2_1h_K0_0.005_LZ95_pbarP_ramp2h` for S14Z0. The scoring read the
new runs from the local copies before they were archived; the archive was verified
against them by checksum (`archive` mode plus S14Z0 by hand, 1 227 s over the network
mount, every file count equal and zero checksum differences in all eight directories).

Scoring (241 s of wall clock): `analyze_dt_convergence.py suite` as registered, output
`output/reflection_tests/s14/s14_suite.json` and `s14_suite.log`. Cavity closure of C2R
against S14C (56 s): `output/reflection_tests/s14c_cavity_closure.{json,png}` and
`s14/s14c_cavity.log`. Supplementary profile ratios (not registered, described below):
`s14/s14_profile_{AB,EF,AD}.json` from `analyze_dt_convergence.py zprofile`.

**Window means** (6–9 h for the Stage 2 cases, 12–17 h for C2R and S14C; `MF` in µPa):

| run | `dt` [s] | `MF(1.25)` | `MF(10.25)` | `MF(30.25)` | `MF(50.25)` | `R10/30` | `R50/30` | `R55/50` |
|---|---|---|---|---|---|---|---|---|
| S04 | 0.659 | -1.878 | 0.1227 | 0.1462 | 0.1080 | 0.839 | 0.739 | 0.400 |
| S14A | 0.330 | -7.521 | 0.1692 | 0.1891 | 0.1352 | 0.894 | 0.715 | 0.382 |
| S14B | 0.165 | -18.882 | 0.1487 | 0.2174 | 0.1520 | 0.684 | 0.699 | 0.371 |
| S14D (`dz` 250 m) | 0.328 | -1.731 | 0.1508 | 0.1771 | 0.1285 | 0.851 | 0.725 | 0.406 |
| S14E (`ublend` 8 km) | 0.659 | 14.466 | 1.1818 | 0.7013 | 0.3362 | 1.685 | 0.479 | 0.035 |
| S14F (`ublend` 8 km) | 0.330 | 16.130 | 1.3488 | 0.7868 | 0.3708 | 1.714 | 0.471 | 0.033 |
| C2R | 0.740 | -1.071 | -1.1280 | -1.1481 | -1.1168 | 0.983 | 0.973 | 1.005 |
| S14C | 0.370 | -1.083 | -1.2016 | -1.2240 | -1.1899 | 0.982 | 0.972 | 1.006 |

**Readings against the registered predictions.** Changes are `(b - a)/|a|` as
registered; for the negative fluxes of C2R and S14C a negative change is an increase in
magnitude.

| reading | measured | H_L predicted | outcome |
|---|---|---|---|
| P1, S14A to S14B | `MF(30.25)` +14.9 %, `MF(1.25)` -151 % (magnitude x2.51) | P1 fails with `MF(30.25)` at +15 % or more | P1 does not pass. `MF(30.25)` falls in the +5 to +15 % band, recorded as neither a pass nor the change H_L predicts |
| P2, S14A to S14B | `R10/30` -23.5 %, `R50/30` -2.2 %, `R55/50` -3.1 % | pass | FAIL, on `R10/30`; prediction failed |
| P3, order | `L` = 0.54 (`MF(30.25)`), 0.52 (`MF(50.25)`), 0.66 (`MF(1.25)`) | `L` of 0.7 or more | prediction failed; values near the first-order reading |
| P4, C2R to S14C | `MF(10.25)`, `MF(30.25)`, `MF(50.25)` -6.5 %, -6.6 %, -6.5 %; level 0.910 to 0.970 (+6.6 %); ripple share 0.0019 in both | below 1 % | prediction failed; above 2 %, so a `dt` dependence intrinsic to the solver in uniform wind |
| P5, S14E to S14F | `MF(10.25)` +14.1 %, `MF(30.25)` +12.2 % | pass | FAIL; prediction failed |
| P6, S14A to S14D | `MF(30.25)` -6.4 %, `MF(1.25)` +77 % (magnitude x0.23); `R50/30` +1.5 %, `R10/30` -4.8 % | absolute changes of 10 % or more, normalized within 5 % | absolute prediction failed for `MF(30.25)`; normalized ratios within 5 % |
| Regression | S14Z below S14Z0 at every summary point | | met (see Launch record) |

None of the predictions that distinguish H_L from the other hypotheses is met. Under the
registered rule, H_L is not supported as the dominant mechanism at these time steps.

**Branch.** P4 exceeds 2 % and P5 fails, so the registered branch is the one for a `dt`
dependence intrinsic to the solver: the next test is switch-level, with the sponge moved
out of the LW sweep into a split exponential step and the splitting order reversed, each
as an on/off pair on C2R. P6 leaves the normalized ratios within 5 %, so the rule does
not require a resolution pair for the Stage 2 amplitudes; `R10/30` at -4.8 % is 0.2
points from the limit. P1 did not pass, so the Richardson row of the budget does not
apply, and production does not move to CFL 0.2 on the strength of this suite.

**Descriptive observations (not registered).**

1. The change per halving of `dt` in `MF(30.25)` depends on the background column. In
   uniform wind (C2R to S14C, 0.74 to 0.37 s) it is +6.6 % in magnitude; with the 8 km
   blend (S14E to S14F, 0.66 to 0.33 s) +12.2 %; on the MCD column +29.4 % (S04 to
   S14A, same `dt` pair as the blend) and +14.9 % (S14A to S14B). A dependence is present
   in all three configurations, and its size increases as the near-ground wind weakens.
   This is consistent with a common term whose effect the calm launch layer amplifies,
   although the suite does not isolate the two.
2. In the uniform-wind pair the top-reflection estimators do not change: 12–17 h means
   of `R_band` 0.0039 in both (+0.4 %), ripple share 0.0019 in both, ripple phase at the
   reference level 22.8 degrees in both, reference-level spectra within 0.1 %. The
   ground-referenced quantities change: `Psi0/imp` 0.985 to 1.055 (+7.1 %), `R_gnd`
   0.318 to 0.374, level +6.6 %. The C2R level agreed with the linear theory of the
   layer (0.910 against 0.911), and S14C sits 6 % above it. The data are consistent
   with a `dt` dependence at or near the forced lower boundary rather than in the sponge
   or along the path.
3. The profile ratios locate the P2 failure. In both S14A and S14B `MF` is negative in
   the launch layer and positive at 10.25 km. The S14B/S14A ratio is 30.7 at 6.25 km and
   negative at 8.25 and 8.75 km, so the two runs cross zero at different heights between
   about 6 and 9.25 km, above the 4.5 km critical level, and 10.25 km lies on the flank
   of that crossing. The S14B/S14A ratio is 0.52 at 9.25 km, 0.88 at 10.25 km
   and 0.99 at 10.75 km, and above 16 km it varies smoothly from 1.156 to 1.124 at
   50.25 km (1.10 at 54.25 km). The failure of `R10/30` reflects this crossing rather
   than a change distributed along the free column. `MF(10.25)` is also not monotonic
   in `dt` (0.123, 0.169, 0.149). Between 11 and 16 km the S14B/S14A ratio alternates
   between adjacent levels by 3 to 5 % (for example 1.089, 1.145, 1.118, 1.165), a
   `2 dz` pattern that is absent from the S14E/S14F and S14A/S14D ratios and is
   consistent with the weaker damping of `2 dz` modes by LW at the smaller Courant
   number of S14B.
4. The S14F/S14E ratio is nearly uniform with height above the lowest level: 0.89 at
   0.25 km, 1.17 to 1.25 between 2 and 6 km, 1.155 at 8.25 km, decreasing smoothly to
   1.103 at 50.25 km. The S14A/S14D ratio is 1.12 to 1.22 between 8.25 and 9.25 km,
   1.097 at 10.25 km, and decreases smoothly to 1.059 at 50.25 km. In both pairs the
   change is set below about 10 km and carried up with a trend of about 4 to 5 % over
   10 to 50 km.
5. At fixed `dt`, the finer grid lowers `|MF(1.25)|` from 7.52 to 1.73 µPa, close to
   the S04 value (1.88), and lowers `MF(30.25)` by 6.4 %. The launch-layer hyperdiffusion
   of H_L falls by a factor of about 3.5 from S14A to S14D, which under H_L would raise
   the fluxes; the measured changes have the opposite sign.
6. Over the `dt` sequence 0.66, 0.33, 0.165 s the increments of `MF(30.25)` shrink
   (+0.043, then +0.028 µPa, ratio 0.66), as do those of `MF(50.25)` (ratio 0.62). A
   first-order Richardson estimate gives `MF(30.25)` 0.246 µPa and an Aitken estimate
   from the observed ratio 0.271 µPa; S04 is 54 to 60 % of these and S14B 80 to 89 %.
   For `MF(50.25)` the two estimates are 0.169 and 0.179 µPa (S04 at 60 to 64 %). These
   figures are descriptive only; the registered rule does not admit them to the budget
   because P1 did not pass. At 1.25 km the increments double (-5.6, then -11.4 µPa), so
   the launch-layer flux shows no approach to a limit over this range.

**Tentative interpretation.** Under the assumptions of the suite, the time-step
dependence of the Stage 2 fluxes is not confined to the calm launch layer. A dependence
of about 6.6 % per halving remains in uniform wind, where the reflection estimators are
unchanged and the ground-referenced launch ratio moves by the same amount. The forced
boundary sets the ghost-cell `w` once per step, before the LW predictor
(`w_ghost = 2 w_forced - w_phys` in `apply_bottom_forcing!`, `src/forcing.jl`), so the
face state that forms the bottom flux carries an `O(dt/dz)` predictor correction. This is
one candidate for a first-order launch error, and it would be consistent with `L` near
0.5; it is not among the hypotheses registered in this section and is recorded as a
candidate only. The registered next step is the switch-level pair on C2R (split sponge,
reversed splitting). A third pair on C2R that imposes the forced `w` on the bottom face
flux directly, independent of the predictor, would test the lower-boundary candidate;
it needs its own predictions written before it runs.

**Consequences for quoted numbers.** The absolute Stage 2 fluxes at CFL 0.4 carry a
`dt` dependence of about 25 to 40 % relative to the smallest time step run here (S04
against S14B: `MF(30.25)` -33 %, `MF(50.25)` -29 %), and of order 6 % in uniform wind.
Ratios taken above about 16 km, such as `R50/30`, change by 3.2 % and 2.2 % over the two
halvings. Ratios normalized at 10.25 km on the MCD column are not `dt`-robust, because
of the zero crossing near 9 km. The `h = 1` km convergence section (Section 15) is
deferred until the switch-level tests identify where the intrinsic dependence enters,
since its design depends on that result.

### Addendum 2026-09-30, ground-forcing test

Written after the Section 14 suite was scored and before any run below. P1 to P6, their
criteria and their recorded outcomes are not changed by this addendum.

**Deviation record.** The Section 14 decision rule names two switch-level pairs on C2R
as the next tests: the sponge moved out of the LW sweep into a split exponential step,
and the splitting order reversed. This addendum runs a ground-forcing test ahead of both.
The reason is the P4 signature. Between C2R and S14C the top-reflection estimators do
not move (12–17 h means from `s14c_cavity_closure.json`: `R_band` +0.2 %, ripple share
-0.4 %, ripple phase at the reference level +0.13 degrees with an hourly scatter of 0.31
degrees), while the ground-referenced amplitude `Psi0/imp` changes by +7.1 % and the
level by +6.6 %. A dependence located in the sponge would appear first in the reflection
estimators, so the P4 data give the sponge pair the lowest prior. The lower boundary is
the part of the solver where an identified discrete path (`apply_bottom_forcing!`)
feeds the imposed `w` into the scheme, so it is tested first. The two registered pairs
are not withdrawn; the criteria below say when they run.

**Near-ground structure used in designing this test.** Two descriptive measurements
were made after scoring and before this addendum, and they informed the predictions:

- The S14C/C2R ratio of `MF(z)` over 12–17 h (`s14/s14_profile_C2R_S14C.json`) is flat
  at 1.065 to 1.066 from 8 to 55 km. Below 8 km it alternates between adjacent levels
  around the same mean: 0.990, 1.112, 1.011, 1.079, 1.062, 1.055, 1.077 at 0.25 to
  3.25 km, with the alternation falling to ±0.1 % by 9 km. The `MF(1.25)` change of
  +1.1 % in the P4 table is a low node of this alternation.
- The mode-6 `w` amplitude over 12–17 h is smooth with height in C2R (3.79, 4.63, 5.14,
  5.60, 6.04, 6.30, 6.54, 6.71 mm/s at 0.25 to 3.75 km) and steps in pairs in S14C
  (4.06, 4.84, 5.45, 5.71, 6.37, 6.40, 6.87, 6.83 mm/s). The imposed amplitude is
  `U k h` = 3.77 mm/s. The linear theory of the C2R column gives `Psi0/imp` = 1.089 at
  the first cell (1.062 with no layer); C2R is 9.5 % below it and S14C 3.1 % below.

The change with `dt` in uniform wind is therefore established within the lowest
kilometre and carried up unchanged, and the smaller time step carries a grid-scale
(`2 dz`) component near the ground that the larger one does not.

**The discrete forcing update, current path (`bottom_forcing = :ghost`).** One step from
`t^n` to `t^{n+1} = t^n + dt`, with ghost rows `j = 1, 2` below the first physical row:

1. `apply_all_bcs!` mirrors `Q^n` into the ghost rows (perturbations of `rho`, `rho u`
   and `E` symmetric, `rho w` antisymmetric).
2. `apply_bottom_forcing!` evaluates `w_f(x, t^n)`, with the ramp `tanh(t^n / t0)`, and
   sets `w_g = 2 w_f - w_1^n` in both ghost rows, keeping the ghost density, `u` and
   pressure.
3. The x-sweep (Eq. 9a, 9b) updates the physical rows only, giving `Q*`. The ghost rows
   keep the state of step 2.
4. The z-predictor (Eq. 10a) forms the ground-face half state
   `Q_{1/2} = (Q_g + Q*_1)/2 - (dt/2dz)(G(Q*_1) - G(Q_g)) + (dt/4)(sources)`.
   The imposed `w` enters here, as a ghost value built at `t^n` from the pre-sweep
   state, combined with the x-swept physical state. The face `w` of `Q_{1/2}` equals
   `w_f` only up to the predictor terms, `(w*_1 - w^n_1)/2` and
   `-(dt/2dz) (G_3(Q*_1) - G_3(Q_g)) / rho` plus the gravity source.
5. The z-corrector (Eq. 10b) of the first physical row takes its ground flux from
   `G(Q_{1/2})`. The imposed `w` reaches the corrector only through step 4.

The imposed `w` therefore enters at the predictor substage, at time level `t^n`, and the
ramp is evaluated at `t^n`.

**Modified path (`bottom_forcing = :face`, driver `bcface=1`, tag `_bcface`).**

1. As above (the mirror of `Q^n` is not used by the z-sweep).
2. `w_f(x, t^{n+1/2})` is evaluated with the ramp `tanh((t^n + dt/2) / t0)` and stored
   per column; no ghost cell is written.
3. The x-sweep is unchanged.
4. After the x-sweep the ghost rows are rebuilt from `Q*`: the reflective mirror of the
   x-swept state, then `w_g = 2 w_f(t^{n+1/2}) - w*_1`. The z-predictor is unchanged in
   form and now combines states of the same time level.
5. In the corrector of the first physical row, the ground-face half state keeps the
   predictor's density, `u` and pressure and takes `w = w_f(t^{n+1/2})`. The mass,
   momentum and energy fluxes and the gravity source at the ground face then carry the
   imposed `w` at the half time level exactly.

The switch removes the time-level mismatch of step 4 and the predictor's modification
of the face `w`. It keeps the reflective mirror (zero `dp'/dz` across the face, an
`O(dz)` term) and the predictor's face density and pressure. In the 12–17 h window C2R
forces with `omega = 0` and `1 - tanh(t / 2 h)` is below 1.3e-5, so the change of time
level contributes less than 1e-5 there; the switch acts in that window through the
spatial and predictor coupling only. With the switch off (default) every earlier record
is reproduced bit for bit, which `check` tests before launch.

**Candidate mechanisms, predicted behaviour and signatures.** Estimates use the C2R
values `dt` = 0.740 s, `dz` = 500 m, `U` = 40 m/s, `k` = 1.885e-5 m^-1, `c` = 228 m/s,
`N(0.25 km)` = 0.0176 s^-1.

| candidate | leading term | estimate at CFL 0.4 | predicted order | where it acts and what responds |
|---|---|---|---|---|
| H_B(a), forcing time level and ramp | `(dt/2) d(w_f)/dt` | below 1e-5 in 12–17 h | any; negligible | nothing measurable in the window |
| H_B(b1), ghost `w` from `Q^n` against x-swept row | `(dt/2) U k` in the face `w` | 2.8e-4, in quadrature with `w` | 1 | ground face; amplitude effect about 4e-8 |
| H_B(b2), mirrored `p'` from `Q^n` against x-swept `p'` | `dt^2 c^2 N / (2 dz U)` in the face `w` | 1.3e-2, in quadrature with `w` | 2 | ground face; amplitude effect about 1e-4 |
| H_B(c), predictor face pressure with a grid-scale `w` component `delta` | `rho c^2 (dt/dz) delta` | `0.08 delta / w_f` of `p'` | 1 through `nu_z`; set by `delta` | first cells; needs `delta / w_f` of order 0.5 to reach 4 % |
| H_B(d), near-ground grid-scale mode | LW damping of `2 dz` per unit time scales as `dt` | not estimated | not fixed; non-convergent if the amplitude scales as the inverse damping | alternation between adjacent levels near the ground; `Psi0/imp`; level through the launched amplitude |
| H_K, sponge inside the LW step | `K0 dt / 2` in the layer | 1.9e-3 in the layer | 2 (in-step), 1 if split | 60–95 km; `R_band`, `R_rip`, ripple phase first; level through `r_top` only, about `2 |r_top| |d r_top|`, below 0.5 % |
| H_S, x-then-z Lie splitting | `(dt/2) [L_x, L_z]` | 3e-3 (`dt N / 2`) to 4e-2 (`dt c m`); the shear term of the registered H_S is zero in uniform wind | 1 | wherever the wave is, including the lowest kilometre; reflection estimators unchanged |
| H_P, LW hyperdiffusion along the path | `kappa_4 m^4` | a trend with height | 1 | a ratio that drifts with height; the flat S14C/C2R ratio above 8 km (1.065 to 1.066) does not show it |

For H_B the identifiable smooth-field terms (a, b1, b2) change the launched amplitude by
1e-4 or less, far below the 6.6 % of P4. Only the grid-scale pathway (c, d) could reach
that size, and its magnitude is not estimated here. The a priori estimates therefore
favour a switch-on result above 2 %; a result below 2 % would place the P4 dependence in
the grid-scale pathway at the forced boundary.

Predicted signatures under each candidate:

- H_B (forcing path). Switch-on pair (G1 against G2): level and `Psi0/imp` differ by
  less than 2 %; the G2/G1 profile ratio of `MF(z)` is flat within 1 % from 0.75 km
  upward, without the alternation of the S14C/C2R ratio; `MF(10.25)`, `MF(30.25)` and
  `MF(50.25)` change by the same fraction as the level. On against off at fixed `dt`
  (G1 against C2R, G2 against S14C): a height-independent shift above about 1 km, with
  the largest difference at the lowest level (0.25 km), whose sign is not predicted.
  `MF(1.25)` is read, but in the switch-off runs it lies on the near-ground alternation
  and is not used for the verdict.
- H_K (sponge). Not tested by these runs. Predicted in its own pair: `R_band`, `R_rip`
  and ripple phase change first, the level by less than 0.5 %. The P4 values already
  place its contribution below 0.5 % per halving.
- H_S (splitting order). Not tested by these runs. Predicted in its own pair: reversing
  the order changes the sign of the leading first-order error, so the reversed-order
  level at CFL 0.4 lies on the other side of the `dt -> 0` limit from C2R by about the
  same distance; reflection estimators unchanged. In the switch-on pair H_S predicts a
  remaining dependence of the same sign as in P4 (above 2 %).

**Convergence order, and the P3 inconsistency.** P3 registered `L = ln G2 / ln G1` with
`G1 = MF(dt/2)/MF(dt)` and `G2 = MF(dt/4)/MF(dt/2)`, reading `L` near 0.5 as first order.
`L` is not the order `p`. For a multiplicative error model `MF = M exp(-a dt^p)`,
`L = 2^-p` exactly. For an additive model `MF = M (1 - a dt^p)`, the increment ratio
`r = (MF(dt/4) - MF(dt/2)) / (MF(dt/2) - MF(dt))` equals `2^-p` exactly, and `L` equals
`2^-p` only in the limit of small changes. With `G1` = 1.294, the first `MF(30.25)`
change on the MCD column, an exactly first-order additive model gives `L` = 0.42,
not 0.5.

On the MCD triple (S04, S14A, S14B), whose mean `dt` ratios are 2.00001 and 2.00001, so
that the adaptive step contributes no ratio error:

| quantity | `L` | `p` multiplicative (`-log2 L`) | increment ratio `r` | `p` additive (`-log2 r`) | reading |
|---|---|---|---|---|---|
| `MF(30.25)` | 0.540 | 0.89 | 0.657 | 0.61 | ambiguous between first order and about 0.6 |
| `MF(50.25)` | 0.521 | 0.94 | 0.617 | 0.70 | ambiguous between first order and about 0.7 |
| `MF(1.25)` | 0.663 | 0.59 | 2.01 (increments grow) | none | not consistent with a convergent order under the additive model |

Three further properties of the MCD triple limit what it can show. The 6–9 h window
lies in the secular growth phase (Z2 exponents 2.1 to 2.3), so a change of `dt` can
shift the growth timing rather than add an asymptotic error term. The MCD response
combines the uniform-wind dependence of P4 (6.6 % per halving) with the
background-dependent part (29 % for the first halving), and a sum of terms of different
order shows an intermediate apparent order. Three points fit one order and leave no
check of the model. The measured 0.52 to 0.66 (`L`) is therefore approximately
consistent with first order for `MF(30.25)` and `MF(50.25)` only under the
multiplicative model, ambiguous overall, and not consistent with a convergent order for
`MF(1.25)` under the additive model. It does not by itself support or exclude a
`dt`-proportional boundary term.

The order is re-measured in uniform wind, where the changes are small enough that the
two models agree: for `G1` = 1.066 an exactly first-order additive model gives `L` =
0.477, that is `p` = 1.07 multiplicative against 1.00 additive. From C2R, S14C and G3
(`dt` 0.740, 0.370, 0.185 s), the order is computed for the level, `Psi0/imp` and
`MF(30.25)` separately, with `p = -log2 r` as the primary value and `-log2 L` alongside:

| band for `p` (additive) | reading | candidates predicting it |
|---|---|---|
| 0.75 to 1.25 | first order | H_B(b1, c), H_S, H_P |
| 1.75 to 2.25 | second order | H_B(b2), H_K in its present in-step form |
| 0.25 or below, including growing increments (`r` of 0.84 or more) | not convergent over this range | H_B(d) |
| any other value, or within 0.05 of a band edge | ambiguous | none |

The paired hourly scatter of the C2R/S14C level difference (0.0598 with a standard
deviation of 0.0010 over six hourly windows) corresponds to about ±0.03 in `p`, below
the band widths.

**Correction to the Section 14 observations (2026-09-30).** In "Suite results", the P3
row of the readings table ("values near the first-order reading") and the "Tentative
interpretation" ("consistent with `L` near 0.5") equate `L` near 0.5 with a measured
order of one. As set out above, that holds only under a multiplicative error model or
for small changes; under the additive model the MCD triple gives `p` of 0.61
(`MF(30.25)`) and 0.70 (`MF(50.25)`), and the order is ambiguous. The statement that a
`dt`-proportional boundary correction "would be consistent with `L` near 0.5" is
therefore withdrawn as evidence; it is tested by the uniform-wind triple instead. The
cavity figures quoted there from the rounded log (`R_band` +0.4 %, ripple share 0.0 %)
are +0.2 % and -0.4 % from the JSON means. P1 to P6 and their outcomes are unchanged.

**Runs.** C2R configuration unchanged (`linear mode=6 h=0.005 hours=18 pbar=1 ramp=2
k0=0.005 lz=95`), all read over 12–17 h as P4, launched in parallel with
`scripts/run_dt_convergence_S2T95.sh launch g` after `check` has passed with the patched
driver (switch off, 0.1 h of S04 and C2R bit for bit).

| id | forcing path | CFL | extra arguments | mean `dt` [s] | steps | wall clock (0.047 to 0.13 s per step) | output directory |
|---|---|---|---|---|---|---|---|
| C2R | ghost | 0.4 | (existing) | 0.740 | 87 514 | done | `arsia_mons_linear18h_LX2000_mode6_K0_0.005_LZ95_pbarP_ramp2h` |
| S14C | ghost | 0.2 | (existing) | 0.370 | 175 029 | done | `..._ramp2h_cfl0.2` |
| S14G1 | face | 0.4 | `bcface=1` | 0.74 | about 87 500 | 1.1 to 3.2 h | `..._ramp2h_bcface` |
| S14G2 | face | 0.2 | `bcface=1 cfl=0.2 outint=200` | 0.37 | about 175 000 | 2.3 to 6.4 h | `..._ramp2h_bcface_cfl0.2` |
| S14G3 | ghost | 0.1 | `cfl=0.1 outint=400` | 0.185 | about 350 000 | 4.6 to 12.8 h | `..._ramp2h_cfl0.1` |

Run directory `run_U40.0_h0.005` in each. The per-step range spans the 0.047 s of S14C
in the seven-process suite and the 0.132 s of C2R (the unexplained per-step cost noted
in Section 13). Disk about 3.1 GB per run. Scoring: `compare_sponge_reflection.py
--cavity` on C2R, S14C, G1, G2, G3 (level, `Psi0/imp`, `R_band`, `R_rip`, ripple phase),
`analyze_dt_convergence.py suite` with the labels C2R, S14C and the three new runs for
`MF(10.25)`, `MF(30.25)`, `MF(50.25)`, and `zprofile` for the G2/G1, G1/C2R and G2/S14C
profile ratios.

**Decision criteria, fixed before launch.** The primary reading is the change in
`Psi0/imp` and in the level between G1 and G2, `(G2 - G1)/|G1|`, both taken as 12–17 h
means:

- Both below 2 % in magnitude (the P4 threshold): the forcing path explains the P4
  dependence.
- The larger of the two between 2 and 5 %: partial explanation; the sponge-step and
  reversed-splitting pairs are run as well.
- The larger above 5 %: the forcing path is not the source; the sponge-step pair runs
  next, then the reversed-splitting pair.

Across all five runs, `R_band` and `R_rip` stay within 1 % of the C2R values and the
ripple phase at the reference level within 3.6 degrees (1 % of a cycle). If any of them
leaves that range, the switch affects more than the lower boundary, and this is
recorded with the reading.

The switch-off order from C2R, S14C and G3 is reported against the bands above for the
level, `Psi0/imp` and `MF(30.25)`.

**Expected level of the switch-on runs, stated before launch.** If H_B holds, the two
paths are consistent discretizations of the same boundary condition, so the switch-on
runs should converge on the `dt -> 0` limit of the switch-off sequence, provided that
limit exists. The first-order extrapolation of C2R and S14C gives a level of 1.030
(expected band 0.99 to 1.07) and `Psi0/imp` of 1.125 (band 1.08 to 1.17), against the
linear theory values 0.911 and 1.089. The G3 result refines the extrapolation but does
not change the band fixed here. The possible outcomes and what each implies for the
Section 13 agreement between C2R and the linear theory (level 0.910 against 0.911):

- Switch-on level near 1.03 (0.99 to 1.07) at both `dt`: the C2R agreement results from
  a `dt`-dependent deficit of about 12 % at CFL 0.4 offsetting a converged excess of
  about 13 % over the linear Rayleigh-layer theory. The Section 11 and 13 level
  closures would be recorded as coincidental at CFL 0.4, and the comparison with theory
  reopened (candidate sources of the excess: vertical resolution at `dz` = 500 m, the
  open-column normalization of the level, the treatment of `N(z)` in the near-ground
  inversion).
- Switch-on level near 0.91 (0.89 to 0.93) at both `dt`: the C2R agreement stands. The
  ghost-path error is small at CFL 0.4 and grows as `dt` decreases, which is the
  behaviour of H_B(d); G3 is then expected to show `p` of 0.25 or below.
- Switch-on level near 0.97 (0.95 to 0.99) at both `dt`: the converged level lies about
  6.5 % above theory, C2R carries the whole ghost-path error at CFL 0.4 and S14C almost
  none, so the ghost-path error is not first order, and the C2R agreement is coincidental
  by about 6.5 %.
- Any other switch-on level at both `dt`: the switch changes the boundary
  representation beyond its `dt` terms (for example at `O(dz)`), and no reading of the
  Section 13 agreement is taken from it.
- Switch-on levels that differ by 2 % or more: the criteria above apply, and the
  expected level is not scored.

Whichever case holds, the readings are recorded against these statements, and a failed
prediction is recorded as failed.

**Launch record.** Patch: `SimConfig.bottom_forcing` (`:ghost` default, `:face`),
`forced_w` and `impose_ghost_w!` factored out of `apply_bottom_forcing!` in
`src/forcing.jl`, the `w_face` keyword of `euler_step!` (`src/euler.jl`), the switch in
the time loop of `src/run.jl`, the NetCDF attribute `bottom_forcing` (`src/io.jl`), the
driver option `bcface=1`, and the case set `g` of the launcher. `check` with the switch
off (64 s): 0.1 h of S04 and C2R reproduced in all 6 and 5 shared snapshots with
`max |diff|` exactly 0.0; driver and `src/` hash `49e29559450b495a`. A 0.05 h smoke run
with the switch on (13.0 s, 244 steps) completed, writes `bottom_forcing = "face"`, and
departs from C2R from the first snapshot after the start (`max |diff|` 1.3e-6 at step 100,
3.7e-6 at step 200), as a changed boundary path should. `launch g`: 2026-09-30 12:49:47,
S14G1 to G3 in parallel, `julia -t 1`, `OPENBLAS_NUM_THREADS=1`, output under
`~/MAGNUSP_runs`, logs `~/MAGNUSP_runs/logs/S14G{1,2,3}.log`, stamps `S14G.start`.

**Wall-clock note for S14G1 to G3.** All three runs stepped at 0.045 to 0.047 s per step
until 13:27:43 and at 0.13 to 0.15 s per step afterwards, with the change occurring
within 2 s in all three (snapshot modification times). Two short returns to the fast
rate (13:30:08 to 13:30:54, 13:46:32 to 13:51:57) also occurred simultaneously in all
three. At 14:51 `macmon` reported the P-core cluster at 702 MHz and the E-core cluster at
912 MHz with 2.6 to 2.7 W of CPU power for three busy processes, CPU temperature 46 C,
load average 4.3, `lowpowermode 0`, normal priority (31, nice 0). `taskpolicy -B` on the
three processes at 14:52:39 left the rate unchanged, so the slowdown is a machine-level
clock reduction rather than a per-process background demotion; its trigger was not
identified. The fields are unaffected (the scheme is deterministic), but the wall clocks
of these runs are not comparable with the Section 14 suite. The per-step cost of C2R
(0.132 s, against 0.047 s for S14C) noted in Section 13 as unexplained is of the same
size and may have the same origin; this is not verified.

### Ground-forcing test results

All three runs completed (exit 0). Solver wall clocks, all affected by the clock
reduction noted above for most of their length: S14G1 7 873.6 s (87 515 steps,
`real` 7 881.8 s), S14G2 18 239.8 s (175 030 steps, 18 247.8 s), S14G3 33 486.8 s
(350 059 steps, 33 494.5 s); last end 2026-09-30 22:08. Mean `dt` 0.7404, 0.3702 and
0.1851 s. Outputs are archived under `output/` with the directory names of the run
table; logs `output/reflection_tests/logs/S14G{1,2,3}.log`. Scoring (12–17 h, as P4):
`compare_sponge_reflection.py --cavity` on the five runs (35 s,
`s14/g_cavity_closure.{json,png}`, `s14/g_cavity.log`), `analyze_dt_convergence.py
suite` (69 s, `s14/g_suite.{json,log}`; the `suite` window rule was extended so that
labels beginning with `S14G` are read over 12–17 h like C2R and S14C), and `zprofile`
ratios `s14/g_profile_{G2_G1,G1_C2R,G2_S14C,G3_S14C}.json`.

**12–17 h means** (`MF` magnitudes in µPa; the fluxes are negative in these runs):

| run | path | `dt` [s] | level | `Psi0/imp` | `R_band` | `R_rip` | ripple phase, ref. level [deg] | `MF(1.25)` | `MF(10.25)` | `MF(30.25)` | `MF(50.25)` |
|---|---|---|---|---|---|---|---|---|---|---|---|
| C2R | ghost | 0.740 | 0.9101 | 0.9848 | 0.00394 | 0.00190 | 22.78 | 1.071 | 1.1280 | 1.1481 | 1.1168 |
| S14C | ghost | 0.370 | 0.9700 | 1.0547 | 0.00394 | 0.00189 | 22.91 | 1.083 | 1.2016 | 1.2240 | 1.1899 |
| S14G3 | ghost | 0.185 | 0.9956 | 1.1187 | 0.00396 | 0.00190 | 22.81 | 1.074 | 1.2412 | 1.2565 | 1.2209 |
| S14G1 | face | 0.740 | 1.0192 | 1.0695 | 0.00396 | 0.00188 | 22.96 | 1.195 | 1.2613 | 1.2868 | 1.2494 |
| S14G2 | face | 0.370 | 1.0196 | 1.0983 | 0.00396 | 0.00188 | 22.97 | 1.137 | 1.2620 | 1.2872 | 1.2500 |

**Primary reading.** From S14G1 to S14G2 the level changes by +0.03 % and `Psi0/imp` by
+2.69 %. The larger of the two lies between 2 and 5 %, so under the registered rule the
outcome is a partial explanation, and the sponge-step and reversed-splitting pairs are
to be run as well. The two quantities disagree: the transmitted flux of the switch-on
pair is independent of `dt` to 0.03 to 0.06 % (`MF(10.25)` +0.06 %, `MF(30.25)` +0.03 %,
`MF(50.25)` +0.05 %), against 6.5 to 6.6 % for the ghost path between the same two time
steps, while the first-cell amplitude still changes by 2.7 %.

**Reflection estimators.** `R_band` stays within 1 % of C2R in all four runs (+0.20,
+0.50, +0.52 and +0.66 % for S14C, G1, G2, G3), and the ripple phase within 0.2 degrees.
`R_rip` departs by -0.44, -1.19, -1.00 and +0.03 %. S14G1 leaves the 1 % range and S14G2
sits at its edge, so under the registered rule this is recorded as the switch affecting
more than the lower boundary. The departure is about 0.4 of the hourly scatter of the
paired difference (2.9 to 3.0 %), and `R_band` and the phase do not move with it.

**Predictions stated before launch, and outcomes.**

| statement | measured | outcome |
|---|---|---|
| H_B: switch-on pair, level changes by less than 2 % | +0.03 % | met |
| H_B: switch-on pair, `Psi0/imp` changes by less than 2 % | +2.69 % | failed |
| H_B: G2/G1 profile ratio flat within 1 % from 0.75 km, without the alternation | alternates 0.951 to 1.036 at 0.75 to 4.75 km; 1.0003 to 1.0005 above 5 km | failed below 5 km |
| H_B: `MF(10.25)`, `MF(30.25)`, `MF(50.25)` change by the same fraction as the level | +0.06, +0.03, +0.05 % against +0.03 % | met |
| H_B: on against off at fixed `dt`, height-independent shift above about 1 km | G1/C2R 1.112 to 1.123 at 0.75 to 50 km; G2/S14C 1.050 to 1.054 | met |
| H_B: largest on/off difference at the lowest level | smallest there (1.067 and 1.036 at 0.25 km) | failed |
| A priori estimate: identifiable smooth-field forcing terms too small for 6.6 %; switch-on result above 2 % favoured | flux: on/off differs by 12 % at CFL 0.4 and the switch-on flux is `dt`-independent; `Psi0/imp`: +2.69 % | failed for the flux; the size estimate of the ghost-path error is low by about two orders of magnitude |
| Expected switch-on level 1.03 (band 0.99 to 1.07) | 1.0192 and 1.0196 | met |
| Expected switch-on `Psi0/imp` 1.125 (band 1.08 to 1.17) | 1.0695 (G1), 1.0983 (G2) | failed for G1, met for G2 |
| Reflection estimators within 1 % and 3.6 degrees | `R_rip` -1.19 % (G1), -1.00 % (G2) | failed for `R_rip` in G1 |

**Order of the ghost path (C2R, S14C, S14G3).** `dt` ratios 2.0003 and 2.0000.

| quantity | values | increment ratio `r` | `p` additive | `p` multiplicative | band reading |
|---|---|---|---|---|---|
| level | 0.9101, 0.9700, 0.9956 | 0.428 | 1.22 (hourly 1.21 to 1.24, sd 0.014) | 1.29 | ambiguous by the edge rule; inside the first-order band, 0.03 below its upper edge |
| `MF(30.25)` | 1.1481, 1.2240, 1.2565 | 0.429 | 1.22 | 1.29 | ambiguous by the edge rule |
| `MF(50.25)` | 1.1168, 1.1899, 1.2209 | 0.423 | 1.24 | 1.31 | ambiguous by the edge rule |
| `MF(10.25)` | 1.1280, 1.2016, 1.2412 | 0.538 | 0.89 | 0.96 | first order |
| `Psi0/imp` | 0.9848, 1.0547, 1.1187 | 0.916 | 0.13 | 0.22 | not convergent over this range |
| `MF(1.25)` | 1.071, 1.083, 1.074 | -0.81 (sign change) | none | none | not defined; near-ground alternation |

The transmitted-flux quantities lie in or at the upper edge of the first-order band,
which the registered table assigns to H_B(b1, c), H_S and H_P; the first-cell amplitude
reads as not convergent, which it assigns to H_B(d). The two additive extrapolations of
the ghost-path level, Aitken from the three runs (1.015) and first-order Richardson from
S14C and S14G3 (1.021), bracket the switch-on level (1.019 to 1.020); for `MF(30.25)` the
Aitken limit is 1.281 µPa against 1.287 in G1 and G2. The ghost path therefore approaches
the switch-on result as `dt -> 0`. This is consistent with the switch removing a
`dt`-dependent error from the ghost path, not with a change of the boundary at `O(dz)`.

**Expected-level outcome and the Section 13 agreement.** The switch-on level (1.019 at
both time steps) falls in the case stated as "near 1.03". As stated before launch, the
C2R agreement with the linear theory (0.910 against 0.911) is then recorded as a
compensation at CFL 0.4: C2R lies 10.7 % below the converged level, which lies 11.9 %
above the linear Rayleigh-layer theory of the layer. The level closures of Sections 11
and 13 are recorded as coincidental at CFL 0.4, and the comparison with theory is
reopened. Candidate sources of the excess named before launch: vertical resolution at
`dz` = 500 m, the open-column normalization of the level, the treatment of `N(z)` in the
near-ground inversion.

**Descriptive observations (not registered).**

1. The near-ground alternation between adjacent levels follows the time step, not the
   forcing path. It is present in S14C, S14G2 and S14G3 (G3/S14C alternates between
   0.983 and 1.059 at 0.25 to 5.25 km, then 1.026 to 1.033 above 10 km) and absent from
   C2R and S14G1; the mode-6 `w` amplitude steps in pairs in S14G2 (6.56, 6.54, 7.07,
   7.00 mm/s at 2.25 to 3.75 km) and is smooth in S14G1 (6.40, 6.66, 6.93, 7.10). The
   on/off ratios at fixed `dt` are smooth. `Psi0/imp` is read at the first cell, which
   lies on this structure, so its +2.7 % in the switch-on pair and its non-convergence on
   the ghost path (`p` 0.13) are properties of the near-ground grid-scale structure at
   smaller Courant number, which the switch does not address and which leaves the
   transmitted flux unchanged. The registered criterion used `Psi0/imp` as one of its two
   primary quantities; with hindsight it mixes the launched amplitude with this structure.
   The criterion is not changed, and the outcome stays "partial".
2. The switch bundles three changes (forcing at `t^{n+1/2}`, ghost rows rebuilt after
   the x-sweep, face `w` imposed in the corrector). The smooth-field estimates assigned
   1e-4 or less to each, so the term responsible for the 12 % on/off difference at CFL 0.4
   is not identified by this test. An ablation with each component alone would
   identify it.
3. Under the ghost path, the transmitted flux at CFL 0.4 in uniform wind is about 11 %
   below the value both paths approach as `dt -> 0`, and about 5 % below at CFL 0.2. The
   Stage 2 records use the ghost path; how much of their larger `dt` dependence on the
   MCD column (29 % for the first halving) the switch removes is not tested here.

**Decision.** The registered rule gives "partial": the sponge-step and reversed-splitting
pairs are run as well. They need code (a split exponential sponge step, a z-then-x sweep
order) and their own predictions written before they run; the P4 and present reflection
data put the expected sponge contribution below 0.5 %. Not registered, and proposed for
a decision: an ablation of the three components of the switch, a Stage 2 MCD pair with
the switch on (the S04 configuration at CFL 0.4 and 0.2), and a test of the near-ground
grid-scale structure.

### Addendum 2 (2026-10-01), face path and the remaining mechanism tests

Written after the ground-forcing test was scored and before any run below was launched.
P1 to P6, the first addendum and their recorded outcomes are not changed.

**Decisions.** All new runs pass `bcface=1` explicitly. The `SimConfig` default stays
`:ghost`, so that every earlier command reproduces its record. No absolute Stage 2 flux
and no `F/F_lin` agreement is quoted until steps 1 to 4 below are complete. Step 4
(the face-path re-verification) precedes any production run and Section 15.

#### Step 1, zero-cost checks on existing records

1. **Hourly level of S14G1 and S14G2** (`s14/g_cavity_closure.json`). The face-path level
   is not stationary over 12–17 h. It rings with the same relative pattern in all five
   runs of the first addendum: 1.119 at 9 h, 1.112 at 10 h, then 1.021, 1.002, 1.009,
   1.028, 1.035 and 1.021 at 12 to 17 h (S14G1; S14G2 differs by +0.03 % in every hour).
   A linear fit over 12–17 h gives +0.33 % per hour (all five runs, +0.32 to +0.33 %);
   the 12–14 h and 15–17 h means are 1.011 and 1.028, and the midpoint of the last trough
   and peak is 1.018. The 12–17 h mean of 1.019 is therefore read with an uncertainty of
   about ±1.5 % from the ringing. Ratios between runs are stationary: G2/G1 is
   1.0003 to 1.0004 in every hour from 4 to 17 h. G1/C2R falls from 1.142 at 4–6 h to
   1.119 to 1.121 at 12–17 h, and G2/S14C from 1.064 to 1.051, so the ghost-path deficit
   was larger during the start-up and has been constant to ±0.1 % since 12 h.
2. **The 2 dz alternation against `dt`** (`analyze_dt_convergence.py zalt`, 12–17 h,
   every second snapshot, 36 s; `s14/g_alternation.json`). The grid-scale part of the
   mode-6 complex `w(z)` is separated by a demodulated second difference with a
   [1, 2, 1] filter. RMS relative amplitude over 0.75–4.75 km:

   | run | path | `dt` [s] | RMS 2 dz / smooth `w` | × `dt` | decay length [km] | × `dt` [km s] | relative amplitude extrapolated to 0.25 km |
   |---|---|---|---|---|---|---|---|
   | C2R | ghost | 0.740 | 0.032 | 0.023 | 0.73 | 0.54 | 0.30 |
   | S14C | ghost | 0.370 | 0.082 | 0.030 | 1.56 | 0.58 | 0.30 |
   | S14G3 | ghost | 0.185 | 0.154 | 0.028 | 3.68 | 0.68 | 0.29 |
   | S14G1 | face | 0.740 | 0.037 | 0.027 | 0.73 | 0.54 | 0.35 |
   | S14G2 | face | 0.370 | 0.085 | 0.032 | 1.56 | 0.58 | 0.32 |

   The decay length is a log-linear fit over 1.25–4.75 km. The RMS over a fixed height
   range scales approximately as `1/dt` (×2.6 and ×1.9 per halving). The scaling sits in
   the decay length (`L dt` 0.54 to 0.68 km s), not in the source: the amplitude
   extrapolated to the first cell is 0.29 to 0.35 of the smooth `w` at every `dt` and on
   both paths. This pattern is consistent with a grid-scale component generated at the
   ground with a `dt`-independent amplitude and damped by the Lax-Wendroff dissipation
   at a rate per unit time proportional to `dt`. The extrapolation to 0.25 km assumes
   exponential decay down to the first cell, where the filter has no value.
3. **Linear theory with the `N(z)` of the near-ground inversion**
   (`analyze_dt_convergence.py ztheory`, `s14/g_theory_N.json`). `N` is 0.0176 s^-1 at
   0.25 km and falls to 0.0092 s^-1 at 3.75 km through the inversion (T from 187 K at
   the ground to 214 K at 3.5 km). Layer-theory flux relative to the value of record
   (V0), and measured S14G1 / theory:

   | variant | theory level | `F_layer / F_layer(V0)` | `Psi(0.25 km)/imp` | S14G1 / theory | C2R / theory |
   |---|---|---|---|---|---|
   | V0 as of record (grid levels, `N` from `dT/dz`, `Phi(0)` extrapolated) | 0.9112 | 1.0000 | 1.0885 | 1.119 | 0.999 |
   | V1 `N^2 = g [(dP0/dz)/(gamma P0) - dln rho0/dz]` of the solver's initial column | 0.8930 | 1.0027 | 1.0892 | 1.116 | 0.996 |
   | V2 25 m grid from `z = 0`, `Phi(0)` imposed at the ground | 0.9113 | 0.9959 | 1.0959 | 1.123 | 1.003 |
   | V3 as V2 with the compressible term | 0.9133 | 0.9983 | 1.0961 | 1.120 | 1.001 |
   | V4 as V2 with `N` of the interpolated `rho` and `P` | 0.9112 | 1.0006 | 1.0962 | 1.118 | 0.998 |

   The treatment of `N(z)` through the inversion changes the predicted flux by -0.4 to
   +0.3 %. It does not account for the +11.9 % excess, which stays at +11.6 to +12.3 %
   across the variants. Of the three candidate sources named before the first addendum,
   this one is not supported. The open-column normalization cancels from the
   measured/theory ratio, which compares absolute fluxes, so it is not a candidate
   either. Vertical resolution remains and is tested by S14R1 below.
4. **Error budget.** Two rows were added to the Section 11 table, not summed with the
   others: the ghost-path forcing (-10.7 % at CFL 0.4, -4.9 % at CFL 0.2) and the
   converged face-path level against theory (+11.9 %, mechanism unattributed).

#### Code for the remaining tests

Each switch defaults to the path of every earlier record:

- `SimConfig.sponge_step = :split` (driver `spsplit=1`, tag `_spsplit`). The LW step
  receives a zero damping matrix. After it, `apply_sponge_exp!` applies the exact
  relaxation `Q = Q_bar + (Q - Q_bar) exp(-K dt)`.
- `SimConfig.sweep_order = :zx` (driver `zx=1`, tag `_zx`). The z-sweep runs before
  the x-sweep. On the face path the ghost rows are rebuilt from the state entering the
  z-sweep, which in this order is `Q^n`.
- `SimConfig.face_probe = N` (driver `probe=N`, no tag). Every `N` steps it records
  the state entering the z-sweep to `face_probe.nc`:
  - the ground-face `w` and `p'` of the z-predictor with the ghost rows as they stand
    (`w_ghost`, `p_ghost`);
  - the same with the ghost rows rebuilt from the current state and `w_f`, as in step 4
    of the face path (`w_rebuilt`, `p_rebuilt`);
  - the plain average of the ghost and first cells (`w_avg`);
  - the forced `w` at the path's own time level (`wf`).

  The probe works on a copy and leaves the solution unchanged.

`euler_step!` is factored into `x_sweep!`, `z_sweep!` and `z_predictor` with the
arithmetic unchanged. The time loop flushes the log at every diagnostic line (every 500
steps), so the per-step cost can be followed live. The NetCDF attributes `sponge_step`
and `sweep_order` are written.

`check` with the switches off (74 s, driver and `src/` hash `79369ec45e1eaa3e`)
reproduced, with `max |diff|` exactly 0.0 in every shared snapshot:
- 0.1 h of S04 (6 snapshots) and C2R (5);
- 0.1 h of the face-path S14G1 (5);
- 0.1 h of C2R with `probe=5` (5), which shows that the probe leaves the solution
  unchanged.

Smoke runs of 0.05 h on the face path departed from the 0.1 h S14G1 check record:
- `zx=1` by 5.4e-7 at step 100 and 8.8e-7 at step 200;
- `spsplit=1` by 4.5e-13 and 8.8e-13 (round-off, since the wave has not reached the
  layer).

**Disclosure.** A 0.2 h smoke run of the probe (ghost path, CFL 0.4, no ramp, record
every 25 steps) was inspected before the predictions below were written. In 3-min
windows from 0 to 12 min:
- `|w_ghost/wf|` at mode 6 was 0.986, 0.959, 0.941 and 0.933, still falling;
- `|w_rebuilt/wf|` was within 0.4 % of it;
- `|w_avg/wf|` was 1.003 to 1.005.

At `t = 0`, the plain average gives `w_avg/wf` = 1.020. This is the density weighting
of the averaged momentum, `1 + (rho_g - rho_1)/(rho_g + rho_1)`, about
`1 + dz/(2 H_rho)`, an `O(dz)` term that does not depend on `dt`.

#### Runs

All runs use `julia -t 1`, `OPENBLAS_NUM_THREADS=1` and output under `~/MAGNUSP_runs`,
and are launched in parallel with `scripts/run_dt_convergence_S2T95.sh launch a2` after
`check` has passed. Logs are `~/MAGNUSP_runs/logs/S14<id>.log`, with stamps in
`S14A2.start`. `monitor a2` writes the per-step cost history and the CPU clocks every
60 s to `S14A2.cost.tsv`.

| id | purpose | arguments | `dt` [s] | steps | est. wall at 0.047 to 0.094 s per step (×2.8 if the clock reduction recurs) | output directory |
|---|---|---|---|---|---|---|
| S14M1 | MCD face pair, CFL 0.4 | `hours=9 pbar=1 ramp=2 k0=0.005 lz=95 u=40 h=0.1 bcface=1` | 0.659 | about 49 200 | 1.3 h (3.6 h) | `arsia_mons_stage2_9h_K0_0.005_LZ95_pbarP_ramp2h_bcface` |
| S14M2 | MCD face pair, CFL 0.2 | as S14M1 plus `cfl=0.2 outint=200` | 0.330 | about 98 300 | 2.6 h (7.2 h) | `..._ramp2h_bcface_cfl0.2` |
| S14K1 | split sponge against S14G1 | `linear mode=6 h=0.005 hours=18 pbar=1 ramp=2 k0=0.005 lz=95 bcface=1 spsplit=1` | 0.740 | about 87 500 | 1.1 h (3.2 h) | `arsia_mons_linear18h_LX2000_mode6_K0_0.005_LZ95_pbarP_ramp2h_bcface_spsplit` |
| S14X1 | z-then-x order against S14G1 | as S14K1 with `zx=1` in place of `spsplit=1` | 0.740 | about 87 500 | 1.1 h (3.2 h) | `..._ramp2h_bcface_zx` |
| S14R1 | `dz` = 250 m against S14G2 at matched `dt` | `linear mode=6 h=0.005 hours=18 pbar=1 ramp=2 k0=0.005 lz=95 bcface=1 dz=250 outint=200` | 0.370 | about 175 000 | 4.6 h (12.8 h) | `..._ramp2h_dz250_bcface` |
| S14P1 | ground-face probe, ghost path, no ramp | `linear mode=6 h=0.005 hours=1 pbar=1 k0=0.005 lz=95 probe=25` | 0.740 | about 4 900 | 4 min | `arsia_mons_linear1h_LX2000_mode6_K0_0.005_LZ95_pbarP` |
| S14P2 | as S14P1 at CFL 0.2 | as S14P1 with `cfl=0.2 probe=50` | 0.370 | about 9 700 | 8 min | `..._LZ95_pbarP_cfl0.2` |
| S14P3 | as S14P1 at CFL 0.1 | as S14P1 with `cfl=0.1 probe=100` | 0.185 | about 19 500 | 15 min | `..._LZ95_pbarP_cfl0.1` |

Run directory `run_U40.0_h0.1` (MCD) or `run_U40.0_h0.005` (single mode). The S14R1
time step equals that of S14G2 because the step is limited by `dz`: CFL 0.4 at 250 m
gives the same `dt` as CFL 0.2 at 500 m. Every probe records every 18.5 s of simulated
time (195 records). Disk is about 21 GB in total.

#### Readings fixed before launch

Window means are taken over 6–9 h for S14M1 and S14M2 (as P1) and over 12–17 h for
S14K1, S14X1 and S14R1 (as P4). Scoring uses:
- `analyze_dt_convergence.py suite` (labels `S14M*` read 6–9 h; `S14K*`, `S14X*` and
  `S14R*` read 12–17 h, by the extended window rule);
- `compare_sponge_reflection.py --cavity` for the single-mode runs;
- `zalt` for S14R1;
- `zprobe` (40–60 min, 10-min sub-windows) for the probes.

A value within 0.5 percentage points of a band edge is reported as at the edge, and the
reading of the higher band applies.

**MCD face pair (S14M1 against S14M2), the decisive Stage 2 test.** The primary
quantity is the `MF(30.25)` change `(M2 - M1)/|M1|`:

- **below 5 %:** Stage 2 production runs at CFL 0.4 on the face path;
- **5 to 15 %:** calm-layer amplification remains, and production waits for its test
  (the reversed-splitting pair on the sheared column, step 4);
- **above 15 %:** the forcing path is not the main Stage 2 term.

Reported alongside, without a band:
- `MF(1.25)`, `MF(10.25)` and `MF(50.25)`;
- the normalized ratios `R10/30`, `R50/30` and `R55/50` against the S04 to S14A changes;
- the fixed-`dt` ratios S14M1/S04 and S14M2/S14A.

*Expected outcome.* Suppose the ghost path loses the same fraction on the MCD column as
in uniform wind (-10.7 % at CFL 0.4, -4.9 % at 0.2). From the 6–9 h S04 and S14A values
(0.1462 and 0.1891 µPa) this transfer estimate gives:
- S14M1 at about 0.164 µPa and S14M2 at about 0.199 µPa;
- a change of about +21 %, in the band above 15 %;
- S14M1/S04 near 1.12 and S14M2/S14A near 1.05.

The P2 outcome (`R10/30` -23.5 % between CFL 0.2 and 0.1) places part of the MCD
dependence in the lower column rather than at the ground, which supports this
expectation. If S14M1/S04 exceeds 1.20, the ghost-path error is larger on the MCD
column than in uniform wind (the MCD ground wind is below 1 m/s, against the 40 m/s
forcing wind), and this is recorded with the reading.

**Split sponge (S14K1 against S14G1), registered null.**
- Primary: the level and `MF(10.25)`, `MF(30.25)` and `MF(50.25)` each change by less
  than 0.5 %.
- Basis: the face pair limits the total `dt` dependence of the transmitted flux to 0.03
  to 0.06 % per halving. The sponge terms are smaller still: the in-step `K0 dt / 2` is
  1.9e-3 in the layer, and the level responds to the layer only through `2 |r_top|
  |d r_top|`.
- Secondary: `R_band` and the ripple share within 1 % of S14G1, the ripple phase within
  1 degree, and `Psi0/imp` within 0.5 %.

A change above 0.5 % in the level or the fluxes fails the null. Because the face pair
shows no `dt` dependence at that size, such a failure would indicate terms that
compensate within the face pair. It is recorded with the reflection estimators.

**Reversed splitting (S14X1 against S14G1), registered null.**
- Primary: the level and the three fluxes each change by less than 0.5 %.
- Basis, under H_S: reversing the order reverses the sign of the leading first-order
  splitting term. The face-pair change bounds that term at about 0.12 % at CFL 0.4, so
  the expected difference is at most about 0.25 %. The shear term of the registered
  H_S is zero in uniform wind.
- Secondary: the reflection estimators within 1 % (phase within 1 degree).
- Recorded without a prediction: `Psi0/imp` and the ground amplitude of the 2 dz
  component, since the order changes the state at which the face is imposed.

A failure is read as for S14K1.

**Vertical resolution (S14R1 against S14G2, `dt` 0.370 s in both).** This tests two
statements.

*(a) The +11.9 % excess.* With the S14G2 level at 1.0196 and the theory at 0.9112,
S14R1 reads as follows:

| S14R1 level | reading |
|---|---|
| 0.999 or above (within 2 % of S14G2, or higher) | vertical resolution at 500 m is not the main source of the excess |
| 0.980 to 0.999 | a `dz` dependence smaller than first order in the excess; recorded without attribution |
| 0.950 to 0.980 (excess halved within about ±2.5 points) | consistent with a first-order `dz` term |
| 0.925 to 0.950 | consistent with a second-order `dz` term |
| below 0.925 | recorded without attribution |

The level is normalized by the theory on each run's own grid, which differs by less
than 0.5 % between the grids (V2 above). `Psi0/imp` is read at 0.125 km in S14R1 and
is not compared with S14G2.

*(b) The 2 dz component shrinks with `dz`.* Over 12–17 h, the decay length of the 2 dz
component is fitted as above, over the same level indices (levels 3 to 10, which are
0.625–2.375 km at 250 m).
- *Setting by the vertical Courant number.* S14R1 has the `nu_z` of S14G1 (0.4). The
  measured scaling in level units (decay length 1.46, 3.12 and 7.36 levels at `nu_z`
  0.4, 0.2 and 0.1) then gives about 1.5 levels, that is 0.37 km (band 0.25 to 0.55
  km). The ground amplitude would be 0.2 to 0.4.
- *Setting by `dt` at fixed physical scale.* The decay length would be near the 1.56 km
  of S14G2.

A decay length below 0.8 km supports the statement that the mode shrinks with `dz`.

**Ground-face probes (S14P1 to S14P3, ghost path, impulsive start, 40–60 min).** The
face-path switch has three components:
- component 1, the forcing at `t^{n+1/2}`, is zero by construction here (`omega` = 0
  and no ramp);
- component 2, the ghost rows rebuilt after the x-sweep;
- component 3, the face `w` imposed in the corrector.

If the flux scales as the square of the face amplitude at fixed face `p'`, the ghost-
to-face flux ratios (0.893, 0.951 and 0.977 at CFL 0.4, 0.2 and 0.1) correspond to
amplitude ratios of 0.945, 0.975 and 0.988. The readings fixed are:

| statement | criterion |
|---|---|
| the face `w` of the ghost path carries the amplitude deficit | `|w_ghost/wf|` at mode 6 in 0.93–0.96 (S14P1), 0.965–0.985 (S14P2), 0.98–0.995 (S14P3) |
| the deficit is first order in `dt` | `1 - |w_ghost/wf|` ratios S14P1/S14P2 and S14P2/S14P3 each 1.6 to 2.4 |
| component 2 does not carry it | `|w_rebuilt/w_ghost|` within 1 % of 1 in all three; component 3 then carries the deficit |
| the deficit lies in the predictor's flux and source terms, not in the averaging | `|w_avg/wf|` within 1 % of 1, and the same in all three within 0.5 % |
| the face pressure is not changed by the rebuild | `|p_ghost/p_rebuilt|` within 1 % of 1 |
| stationarity | the 40–50 and 50–60 min values of `|w_ghost/wf|` differ by less than 0.5 %; otherwise the reading is provisional |

If `|w_ghost/wf|` lies outside its band while the components are as stated, the face
amplitude does not map onto the flux deficit one to one. That points to another face
quantity, or to the near-ground 2 dz state, and is recorded. The 40–60 min state of an
impulsive start is not the 12–17 h state of the ramped runs. The probes identify the
component, and the size is read against the bands above as an approximation.

Whichever case holds, the readings are recorded against these statements, and a failed
prediction is recorded as failed.

#### Results, ground-face probes (S14P1 to S14P3)

Launch `launch a2` at 2026-10-01 14:04:10 (driver and `src/` hash `1fd7c912d485b418`,
after a repeated `check` of 72 s, all four comparisons identical). The probes completed
with exit 0:

| id | `dt` [s] | steps | solver wall clock [s] | `real` [s] | s per step |
|---|---|---|---|---|---|
| S14P1 | 0.740 | 4 862 | 245.7 | 254.4 | 0.051 |
| S14P2 | 0.370 | 9 723 | 481.7 | 490.5 | 0.050 |
| S14P3 | 0.185 | 19 445 | 959.7 | 968.5 | 0.049 |

These were run with seven other processes, and CPU power stayed at 44 to 50 W, so the
clock reduction of the first addendum did not recur here. Outputs are in
`~/MAGNUSP_runs/arsia_mons_linear1h_LX2000_mode6_K0_0.005_LZ95_pbarP{,_cfl0.2,_cfl0.1}`
and the logs are `S14P{1,2,3}.log`. Scoring took 2 s: `analyze_dt_convergence.py
zprobe`, written to `s14/a2_probe.json`.

Mode-6 ratios over 40–60 min (65 records each; magnitude and phase in degrees):

| run | `w_ghost/wf` | `w_rebuilt/wf` | `w_avg/wf` | `p_ghost/p_rebuilt` | `|w_ghost/wf|` 40–50 / 50–60 min |
|---|---|---|---|---|---|
| S14P1 | 0.9128 (+7.4) | 0.9090 (+8.4) | 1.0077 (-0.7) | 1.0068 (-0.8) | 0.8969 / 0.9295 |
| S14P2 | 0.9557 (+3.7) | 0.9544 (+3.9) | 1.0080 (-0.8) | 1.0038 (-0.4) | 0.9479 / 0.9639 |
| S14P3 | 0.9806 (+1.4) | 0.9802 (+1.5) | 1.0082 (-0.9) | 1.0020 (-0.2) | 0.9775 / 0.9838 |

| statement | measured | outcome |
|---|---|---|
| `|w_ghost/wf|` in 0.93–0.96, 0.965–0.985 and 0.98–0.995 | 0.913, 0.956, 0.981 | failed for S14P1 and S14P2 (below the bands); met for S14P3, at the edge |
| deficit first order, ratios 1.6 to 2.4 | deficits 0.087, 0.044 and 0.019; ratios 1.97 and 2.28 | met |
| component 2 does not carry it, `|w_rebuilt/w_ghost|` within 1 % | 0.996, 0.999, 1.000 | met |
| `|w_avg/wf|` within 1 % and equal within 0.5 % | 1.0077 to 1.0082 (spread 0.05 %) | met |
| `|p_ghost/p_rebuilt|` within 1 % | 1.0068, 1.0038, 1.0020 | met |
| stationarity, 40–50 and 50–60 min within 0.5 % | differ by 3.6, 1.7 and 0.6 % | failed; the readings are provisional |

**Reading.**
- Rebuilding the ghost rows after the x-sweep (component 2) changes the predictor's
  face `w` by 0.4 % or less and the face `p'` by 0.7 % or less. The plain average
  carries a +0.8 % offset that is the same at all three `dt`.
- The `dt`-dependent part of the face `w` therefore comes from the predictor's flux and
  source terms. On the face path these are removed only by imposing `w` in the corrector
  (component 3). Component 1 is zero by construction in these runs and below 1e-5 in
  the 12–17 h window of the ramped runs.
- The data are consistent with component 3 carrying the ghost-path deficit.
- The deficit is first order in `dt` in magnitude (ratios 1.97 and 2.28). So is its
  phase lag (+7.4, +3.7 and +1.4 degrees).

**Failed predictions.**
- The face-amplitude deficit is 1.6 to 1.8 times the value implied by the flux ratios
  (0.087 against 0.055 at CFL 0.4, 0.044 against 0.025 at CFL 0.2, 0.019 against 0.012
  at CFL 0.1). The predicted one-to-one mapping from face amplitude onto flux does not
  hold under these conditions.
- None of the three records is stationary at 40–60 min. The ratio was still rising by
  1 to 3 % per 10 min at CFL 0.4.

Two readings remain open, and the probes do not separate them:
- the impulsive start, whose 1 h state differs from the 12–17 h state of the ramped
  runs;
- a flux that responds to the face `w` and the face `p'` together.

The component attribution does not rest on the band values.

#### Results, split sponge and reversed splitting (S14K1, S14X1)

Both completed with exit 0:

| id | steps | solver wall clock [s] | `real` [s] | s per step | ended |
|---|---|---|---|---|---|
| S14K1 | 87 515 | 6 397.7 | 6 406.4 | 0.073 | 15:51:01 |
| S14X1 | 87 515 | 6 379.5 | 6 388.1 | 0.073 | 15:50:44 |

Both wall clocks include the idle-sleep period described below. Scoring took 9 s for
`compare_sponge_reflection.py --cavity` on S14G1, S14K1 and S14X1 (`s14/a2kx/cavity_closure.{json,png}`,
`s14/a2kx_cavity.log`) and 128 s for `analyze_dt_convergence.py suite`, which also covered
S04 and S14M1 (`s14/a2kx/s14_suite.json`, `s14/a2kx_suite.log`). All values are 12–17 h
means, given as changes against S14G1:

| quantity | S14K1 (split sponge) | S14X1 (z then x) | registered bound |
|---|---|---|---|
| level | -0.022 % | -0.042 % | 0.5 % |
| `MF(10.25)` | -0.033 % | -0.029 % | 0.5 % |
| `MF(30.25)` | -0.023 % | -0.046 % | 0.5 % |
| `MF(50.25)` | -0.013 % | -0.049 % | 0.5 % |
| `MF(1.25)` (read, no bound) | -0.041 % | +0.038 % | |
| `Psi0/imp` | +0.001 % | -0.007 % | 0.5 % (S14K1); none (S14X1) |
| `R_band` | -0.015 % | -0.21 % | 1 % |
| ripple share `R_rip` | -0.020 % | +0.71 % | 1 % |
| ripple phase, reference level | +0.38 deg | +0.02 deg | 1 deg |

Both registered nulls are met, with every primary quantity at 0.05 % or less. The
secondary statements are met as well. The S14X1 ripple share (+0.71 %) is the largest
change, and it lies within the 1 % range and well inside the hourly scatter of that
estimator (about 3 %). Within these runs neither the in-step treatment of the damping
layer nor the order of the dimensional splitting carries a measurable part of the flux,
level or first-cell amplitude in uniform wind. The data do not indicate terms that
compensate within the face pair. The H_K and H_S candidates of the first addendum are not
supported in uniform wind at CFL 0.4. The registered H_S concerns the sheared column, for
which step 4 keeps a reversed-splitting pair if the MCD face pair leaves a residual.

#### S14M1 at fixed `dt`, ahead of S14M2

S14M1 completed with exit 0 (49 153 steps, solver 6 824.4 s, `real` 6 832.9 s, 0.139 s
per step including the sleep period, ended 15:58:03). The registered reading needs
S14M2. The fixed-`dt` comparison with the ghost-path S04 (6–9 h, µPa) is recorded now:

| quantity | S04 (ghost) | S14M1 (face) | S14M1 / S04 |
|---|---|---|---|
| `MF(1.25)` | -1.878 | -7.163 | 3.81 |
| `MF(10.25)` | 0.1227 | 0.2377 | 1.94 |
| `MF(30.25)` | 0.1462 | 0.2404 | 1.64 |
| `MF(50.25)` | 0.1080 | 0.1656 | 1.53 |
| `R10/30` | 0.839 | 0.989 | 1.18 |
| `R50/30` | 0.739 | 0.689 | 0.93 |
| `R55/50` | 0.400 | 0.362 | 0.90 |

The expectation stated before launch (S14M1 about 0.164 µPa, S14M1/S04 about 1.12)
fails. S14M1/S04 is 1.64 at 30.25 km, above the 1.20 that the addendum fixed as the
signature of a ghost-path error larger on the MCD column than in uniform wind. This is
recorded as such. S14M1 at CFL 0.4 (0.240 µPa) lies above the ghost-path S14B at CFL 0.1
(0.217 µPa) and below the extrapolated ghost-path limits of the Section 14 suite (0.25 to
0.27 µPa). The ratio of 1.64 is not uniform with height (3.81 at 1.25 km, 1.94 at
10.25 km, 1.53 at 50.25 km), so on the MCD column the switch changes the vertical shape
of the flux as well as its level.

#### Cause of the wall-clock slowdowns (correction to the first addendum's wall-clock note)

`monitor a2` recorded CPU power falling from 33 to 34 W to 0.5 to 4.6 W, with the step
rate falling by a factor of 2 to 8. Between 14:47 and 15:01 there were no rows at all,
because the monitor itself was suspended. `pmset -g log` places the start at "Entering
Sleep state due to 'Maintenance Sleep'" at 14:42:18, followed by DarkWake and Sleep
cycles about every 45 s until 15:03 and a reduced clock until 15:28. The machine is set
to `pmset sleep 1`, so it sleeps after one minute without user input or an idle-sleep
assertion. The Claude Code process holds such an assertion only while it is active
(`caffeinate -i -t 300`). The slowdown of 2026-09-30 matches the same mechanism: "Entering
Sleep" at 13:27:27, 16 s before the step rate fell at 13:27:43 in S14G1 to G3, with
DarkWake at 13:30:09 and 13:46:29 matching the two short returns to the fast rate. The
statement in the first addendum that the trigger was not identified is superseded. The
cause is idle system sleep with DarkWake running, not a thermal or per-process limit.
The unexplained per-step cost of C2R in Section 13 may have the same origin; this is not
verified.

From 16:02, `caffeinate -i -w <monitor pid>` holds an idle-sleep assertion for the
remaining runs. The launcher's `monitor` mode now takes this assertion itself, for as
long as any case of its set is running. The fields are unaffected, since the scheme is
deterministic and the runs only slowed.

#### Results, MCD face pair (S14M1, S14M2)

S14M2 completed with exit 0: 98 305 steps, solver 14 892.7 s, `real` 14 901.4 s, 0.151 s
per step, ended 18:12:34. Its wall clock includes two sleep periods (14:42 to 15:28 and
16:06 to 18:04). The idle-sleep assertion of 16:02 (`caffeinate -i`) was not honoured,
because the machine was already in DarkWake. From 18:04, `caffeinate -s -i`
(PreventSystemSleep, honoured on AC) held, and `pmset -g log` records no "Entering
Sleep" after 18:05. The launcher's `monitor` now uses `-s -i`. Scoring took 53 s:
`analyze_dt_convergence.py suite` on S04, S14A, S14B, S14M1 and S14M2
(`s14/a2m/s14_suite.json`, `s14/a2m_suite.log`). All values are 6–9 h means in µPa.

| run | path | `dt` [s] | `MF(1.25)` | `MF(10.25)` | `MF(30.25)` | `MF(50.25)` | `R10/30` | `R50/30` | `R55/50` |
|---|---|---|---|---|---|---|---|---|---|
| S04 | ghost | 0.659 | -1.878 | 0.1227 | 0.1462 | 0.1080 | 0.839 | 0.739 | 0.400 |
| S14A | ghost | 0.330 | -7.521 | 0.1692 | 0.1891 | 0.1352 | 0.894 | 0.715 | 0.382 |
| S14B | ghost | 0.165 | -18.882 | 0.1487 | 0.2174 | 0.1520 | 0.684 | 0.699 | 0.371 |
| S14M1 | face | 0.659 | -7.163 | 0.2377 | 0.2404 | 0.1656 | 0.989 | 0.689 | 0.362 |
| S14M2 | face | 0.330 | -15.789 | 0.2393 | 0.2412 | 0.1660 | 0.992 | 0.688 | 0.362 |

Changes per halving of `dt`:

| quantity | ghost, S04 to S14A | ghost, S14A to S14B | face, S14M1 to S14M2 |
|---|---|---|---|
| `MF(1.25)` | +300.5 % | +151.1 % | +120.4 % |
| `MF(10.25)` | +37.9 % | -12.1 % | +0.64 % |
| `MF(30.25)` | +29.4 % | +14.9 % | +0.35 % |
| `MF(50.25)` | +25.1 % | +12.4 % | +0.27 % |
| `R10/30` | +6.6 % | -23.5 % | +0.30 % |
| `R50/30` | -3.3 % | -2.2 % | -0.07 % |
| `R55/50` | -4.5 % | -3.1 % | -0.03 % |

**Registered reading.** The `MF(30.25)` change from S14M1 to S14M2 is +0.35 %, in the
band below 5 %. Under the rule fixed before launch, Stage 2 production runs at CFL 0.4
on the face path. Step 4 (re-verification of the face path) still precedes any
production run. The hourly `MF(30.25)` values of the two runs differ by 0.2 to 0.4 % at
6.5, 7.5 and 8.5 h (0.1984 against 0.1988, 0.2302 against 0.2310, 0.2925 against
0.2938). The secular growth of the 6–9 h window (Section 14, Z2) is common to both runs
and does not enter their ratio.

**Statements fixed before launch, and outcomes.**

| statement | measured | outcome |
|---|---|---|
| expected `MF(30.25)` change about +21 % (band above 15 %) | +0.35 % | failed |
| expected S14M1 about 0.164 µPa | 0.240 µPa | failed |
| expected S14M1/S04 near 1.12; above 1.20 recorded as a larger ghost-path error on the MCD column | 1.64 | the larger-error case holds |
| expected S14M2/S14A near 1.05 | 1.275 | failed |

The reasoning given for the expectation also fails. It cited the P2 change of `R10/30`
(-23.5 % from CFL 0.2 to 0.1) as placing part of the MCD dependence in the lower column
rather than at the ground. On the face path `R10/30` changes by +0.30 % per halving, so
that shape change is also a property of the ghost path.

**Descriptive observations (not registered).**

1. Above the 4.5 km critical level, the face path removes the `dt` dependence of the
   transmitted flux on the MCD column to within 0.3 to 0.6 % per halving. The ghost path
   shows 25 to 38 % for the first halving at the same levels. At CFL 0.4 the ghost-path
   loss on this column is 39 % of the face-path value at 30.25 km (S04/S14M1 = 0.608),
   against 10.7 % in uniform wind. The data are consistent with the ghost-path error being
   amplified when the background wind at the forced face is near zero (MCD `u` below 1 m/s
   below 4 km) while the forcing wind is 40 m/s. This is not tested separately.
2. `MF(1.25)`, inside the calm layer below the 4.5 km critical level, still changes by
   +120 % for one halving on the face path. No band was registered for it. A
   `dt`-dependent flux below the critical level remains, without a measurable effect at
   10.25 km and above in this window. Its mechanism is not identified. The `MF(1.25)` of
   the ghost path (S14B, CFL 0.1) exceeds that of the face path at CFL 0.2 in magnitude
   (-18.9 against -15.8 µPa).
3. The ghost-path values of the Section 14 suite do not converge monotonically toward
   the face-path values at every level. S14B lies above S14A at 30.25 and 50.25 km and
   below it at 10.25 km. At 30.25 km the S14B value (0.217) is 10 % below the face path
   (0.241). The Section 14 Aitken limit for `MF(30.25)` (0.271 µPa) overestimates the face
   value by 12 %, consistent with the ambiguous order recorded in the first addendum.

**Consequences for the plan.** The step 4 condition "if the MCD face pair leaves a
residual, add a reversed-splitting pair on the sheared column" is not triggered for the
transmitted flux. `MF(1.25)` is the one residual. Whether it calls for the sheared-column
pair is a decision left open here. It was not registered as a primary quantity, and the
reversed-splitting null in uniform wind does not bear on it.

#### Results, vertical resolution (S14R1)

S14R1 completed with exit 0: 175 038 steps, solver 22 201.4 s, `real` 22 210.0 s, 0.127
s per step, ended 20:14:28. Its wall clock includes both sleep periods. Scoring took 19 s
in total: `--cavity` on S14G2 and S14R1 (`s14/a2r/cavity_closure.{json,png}`,
`s14/a2r_cavity.log`), `suite` (`s14/a2r/s14_suite.json`) and `zalt`
(`s14/a2r_alternation.json`). All values are 12–17 h means at `dt` = 0.3702 s in both
runs.

| quantity | S14G2 (`dz` 500 m) | S14R1 (`dz` 250 m) | change |
|---|---|---|---|
| level (own-grid open-column normalization) | 1.0196 | 0.9999 | -1.93 % |
| layer-theory level on the run's grid | 0.9112 | 0.9072 | |
| measured / own-grid theory | 1.119 | 1.102 | |
| `MF(10.25)`, `MF(30.25)`, `MF(50.25)` [µPa] | -1.2620, -1.2872, -1.2500 | -1.2425, -1.2634, -1.2241 | -1.55, -1.85, -2.07 % |
| `MF(1.25)` [µPa] | -1.137 | -1.233 | +8.4 % |
| `R_band` | 0.00396 | 0.00388 | -2.0 % |
| ripple share | 0.00188 | 0.00188 | 0.0 % |
| ripple phase, reference level [deg] | 22.97 | 20.64 | -2.3 deg |
| `Psi0/imp` (0.25 and 0.125 km) | 1.098 | 1.016 | not compared |

The hourly level ratio S14R1/S14G2 is 0.9801 to 0.9813, so the difference is stationary
over the window.

**(a) Excess over theory.** The S14R1 level of 0.9999 lies 0.09 points above the 0.999
edge. Under the edge rule, the reading of the higher band applies: vertical resolution at
500 m is not the main source of the +11.9 % excess. Halving `dz` lowers the excess over
the own-grid theory from +11.9 to +10.2 %, about one seventh of it, and lowers the
transmitted flux by 1.6 to 2.1 %. A `dz` dependence of this size is present. It is
recorded without attribution, since a first-order `dz` term in the whole excess would
have halved it. The excess remains unattributed: the `N(z)` treatment, the open-column
normalization and the vertical resolution together account for at most about 2 of its
12 points.

**(b) The 2 dz component.** The 2 dz component was fitted over levels 3 to 10.

| run | decay length | in levels | relative amplitude extrapolated to the first cell | RMS over 0.75–4.75 km |
|---|---|---|---|---|
| S14G2 | 1.56 km | 3.13 | 0.32 | 0.085 |
| S14R1 | 0.36 km | 1.38 | 0.22 | 0.007 |

Both registered statements are met: the decay length lies in the band of 0.25 to
0.55 km, and the amplitude in the band of 0.2 to 0.4. The decay length below 0.8 km
supports the statement that the mode shrinks with `dz`. In level units the decay length
(1.38) is close to that of S14G1 and C2R at the same vertical Courant number of 0.4
(1.46). The data are consistent with the 2 dz component being set by the vertical
Courant number, with a ground amplitude of 0.2 to 0.35 of the smooth `w`, as a
grid-scale feature of the Lax-Wendroff treatment of the forced face.

#### Summary of addendum 2

| test | registered reading | outcome |
|---|---|---|
| MCD face pair, `MF(30.25)` S14M1 to S14M2 | below 5 %: Stage 2 production at CFL 0.4 on the face path | +0.35 %; production at CFL 0.4 on the face path, after step 4 |
| split sponge (S14K1) | null, below 0.5 % | met (0.05 % or less) |
| reversed splitting (S14X1) | null, below 0.5 % | met (0.05 % or less) |
| `dz` = 250 m (S14R1), excess | 0.999 or above: `dz` not the main source | at the edge (0.9999); `dz` not the main source; -1.9 % recorded |
| `dz` = 250 m (S14R1), 2 dz mode | decay length 0.25 to 0.55 km; below 0.8 km: shrinks with `dz` | 0.36 km; met |
| probes, component attribution | component 2 below 1 %, so component 3 carries the deficit; first order | met; size bands failed for S14P1 and S14P2; stationarity failed |
| expected MCD outcome (+21 %, S14M1 about 0.164 µPa) | | failed |

Open items carried into step 4 and beyond:
- `MF(1.25)` below the MCD critical level remains `dt`-dependent on the face path (+120 %
  per halving).
- The +11.9 % excess of the converged uniform-wind level over the linear theory is
  unattributed.
- The face-path level rings by about ±1.5 % over 12–17 h.
- Step 4 is next: the H = 0 equilibrium with `w_f = 0` at the face, the analytic
  `F(z)/F_lin` check and the Queney case on the face path. The Queney driver
  (`run_mountain_benchmark_strict.jl`) pins `pbar_mode = :total` and has no
  forcing-path option, so it needs an option and its own criteria before it runs.

**Wall clocks.** The launch ran from 14:04:10 to 20:14:28, 6.17 h of wall clock. The CPU
times summed over the `real` values are 16.2 h for the eight runs (58 451 s). Sleep slowed every run
that was active between 14:42 and 15:28 or between 16:06 and 18:04. The `monitor`
process ended with a syntax error after its loop had completed, because the launcher was
edited while it ran (bash reads a script incrementally). The cost history to 20:14 is
complete in `S14A2.cost.tsv`. The launcher is not to be edited while one of its modes is
running.

**Archive.** `archive a2` (999 s over the network mount) copied all eight run directories
to `output/` under the names in the run table. Every file count is equal and there are
zero checksum differences (493, 493, 877, 877, 877, 51, 100 and 197 files). The logs,
`S14A2.start` and `S14A2.cost.tsv` are in `output/reflection_tests/logs/`. The local copies
in `~/MAGNUSP_runs` are kept.

### Addendum 3 (2026-10-02), face-path verification and the remaining Stage 2 readings

Written on 2026-10-02, before any run below was launched. P1 to P6, the first addendum,
addendum 2 and their recorded outcomes are not changed. Policy, unchanged: every new run
passes `bcface=1` explicitly, and the `SimConfig` default stays `:ghost`. No absolute
Stage 2 flux and no `F/F_lin` agreement is quoted until the readings registered here are
scored.

#### Step 1, zero-cost checks on archived records

Scripts and timing: `analyze_dt_convergence.py zalt`, with the decay fit, the model and an
x-rms form added (44 s for the six uniform-wind records, `s14/a3/zalt_uniform.json`; 93 s
for the four MCD records and three `zprofile` ratios, `s14/a3/zalt_mcd.json`,
`s14/a3/profile_*.json`). The hourly values come from `s14/a2m/s14_suite.json` with no new
snapshot reads.

**a. Hourly fluxes of the MCD face pair and the ghost pair (6 to 9 h, µPa).**

| hour | run | `MF(10.25)` | `MF(30.25)` | `MF(50.25)` | `R10/30` | `R50/30` |
|---|---|---|---|---|---|---|
| 6 to 7 | S04 (ghost, 0.659 s) | 0.1142 | 0.1251 | 0.0854 | 0.913 | 0.683 |
| | S14A (ghost, 0.330 s) | 0.1529 | 0.1589 | 0.1046 | 0.962 | 0.658 |
| | S14M1 (face, 0.659 s) | 0.2011 | 0.1984 | 0.1259 | 1.014 | 0.634 |
| | S14M2 (face, 0.330 s) | 0.2034 | 0.1988 | 0.1261 | 1.023 | 0.634 |
| 7 to 8 | S04 | 0.1206 | 0.1442 | 0.1083 | 0.836 | 0.751 |
| | S14A | 0.1666 | 0.1840 | 0.1346 | 0.906 | 0.732 |
| | S14M1 | 0.2365 | 0.2302 | 0.1638 | 1.027 | 0.711 |
| | S14M2 | 0.2377 | 0.2310 | 0.1642 | 1.029 | 0.711 |
| 8 to 9 | S04 | 0.1333 | 0.1692 | 0.1304 | 0.788 | 0.770 |
| | S14A | 0.1880 | 0.2245 | 0.1663 | 0.838 | 0.741 |
| | S14M1 | 0.2756 | 0.2925 | 0.2071 | 0.942 | 0.708 |
| | S14M2 | 0.2767 | 0.2938 | 0.2078 | 0.942 | 0.707 |

Changes per halving of `dt`, by hour (per cent):

| hour | face `MF(10.25)` | face `MF(30.25)` | face `MF(50.25)` | face `R10/30` | face `R50/30` | ghost `MF(30.25)` | ghost `R10/30` | ghost `R50/30` |
|---|---|---|---|---|---|---|---|---|
| 6 to 7 | +1.15 | +0.20 | +0.15 | +0.95 | -0.05 | +27.0 | +5.4 | -3.6 |
| 7 to 8 | +0.50 | +0.34 | +0.25 | +0.16 | -0.10 | +27.6 | +8.3 | -2.5 |
| 8 to 9 | +0.39 | +0.45 | +0.36 | -0.05 | -0.08 | +32.7 | +6.3 | -3.9 |

- The face-pair `MF(30.25)` change is +0.20, +0.34 and +0.45 % in the three hours. It is
  below the 5 % edge of the registered band in every hour. The window value of +0.35 % is
  the mean of a change that grows by about 0.12 points per hour, at the same time as the
  secular growth of the flux itself (Z2).
- The ghost-pair change at the same level is 27 to 33 % in every hour.
- The magnitude of `MF(1.25)` on the face path rises by 108, 119 and 130 % per halving in
  the three hours, so its `dt` dependence is not stationary either (item c).

Window-mean `MF(z)` ratios over 6 to 9 h (`zprofile`):

| z [km] | 1.25 | 2.25 | 4.75 | 10.25 | 20.25 | 30.25 | 40.25 | 50.25 | 55.25 |
|---|---|---|---|---|---|---|---|---|---|
| S14M1 / S04 (face / ghost, 0.659 s) | 3.82 | 2.10 | 0.74 | 1.94 | 1.70 | 1.645 | 1.60 | 1.53 | 1.39 |
| S14M2 / S14A (face / ghost, 0.330 s) | 2.10 | 2.17 | 1.07 | 1.41 | 1.30 | 1.275 | 1.255 | 1.23 | 1.16 |
| S14M2 / S14M1 (face, 0.330 / 0.659 s) | 2.20 | 3.49 | 2.85 | 1.006 | 1.004 | 1.003 | 1.003 | 1.003 | 1.002 |

The flux changes sign between about 5 and 8 km, above the 4.5 km critical level, so the
ratios at those heights carry no information and are left out. Above 10 km:

- The face path raises the transmitted flux at every level, by more at lower levels. The
  face/ghost ratio falls with height by 0.0066 per km at CFL 0.4 (1.94 at 10.25 km to 1.53
  at 50.25 km) and by 0.0029 per km at CFL 0.2.
- On the face path `MF(10.25)` and `MF(30.25)` are nearly equal (`R10/30` 0.989 over
  the window), and the loss from 30 to 50 km is larger than on the ghost path (`R50/30`
  0.689 against 0.739 at CFL 0.4). The ghost-path `R50/30` moves toward the face value as
  `dt` falls (0.739, 0.715, 0.699). The ghost-path `R10/30` does not (0.839, 0.894, 0.684).
- The two face-path records agree to 0.2 to 0.6 % at every level from 10 to 58 km.

The ghost path therefore changes the vertical shape of the MCD flux profile as well as its
level. Normalized ratios from ghost-path MCD records carry this shape error: at CFL 0.4,
`R10/30` is 15 % low and `R50/30` 7 % high against the face path.

**b. The 2 dz component against a descriptive model** (descriptive, not registered).
The decay length of the 2 dz part of the mode-6 `w(z)` is fitted log-linearly over levels
3 to 10 (1.25 to 4.75 km at 500 m, 0.625 to 2.375 km at 250 m), 12 to 17 h, every second
snapshot, and set beside `L = dz^2 / (2 c dt)` with `c = (gamma R T(0))^(1/2)` = 215.3 m/s
at the MCD ground (187.3 K).

| run | path | `dz` [m] | `dt` [s] | `c dt / dz` | fitted decay [km] | `L` [km] | fitted / `L` | relative amplitude extrapolated to the first cell |
|---|---|---|---|---|---|---|---|---|
| C2R | ghost | 500 | 0.7405 | 0.319 | 0.730 | 0.784 | 0.93 | 0.30 |
| S14C | ghost | 500 | 0.3702 | 0.159 | 1.559 | 1.568 | 0.99 | 0.30 |
| S14G3 | ghost | 500 | 0.1851 | 0.080 | 3.679 | 3.137 | 1.17 | 0.29 |
| S14G1 | face | 500 | 0.7404 | 0.319 | 0.731 | 0.784 | 0.93 | 0.35 |
| S14G2 | face | 500 | 0.3702 | 0.159 | 1.565 | 1.568 | 1.00 | 0.32 |
| S14R1 | face | 250 | 0.3702 | 0.319 | 0.358 | 0.392 | 0.91 | 0.21 |
| S14M1 | face, MCD | 500 | 0.6592 | 0.284 | not defined | 0.881 | | |
| S14M2 | face, MCD | 500 | 0.3296 | 0.142 | not defined | 1.762 | | |

- In uniform wind the model gives the fitted decay length to within -9 to +17 % over a
  factor 4 in `dt` and 2 in `dz`. The fitted length scales as `dz^2.13` (S14G2 to S14R1)
  and grows by 2.14 and 2.36 per halving of `dt` (C2R to S14C to S14G3). The model's
  exponents are 2 and 2.
- The ground amplitude is not `dz`-independent at this precision. At matched `dt` it is
  0.21 at 250 m against 0.32 at 500 m. The 250 m extrapolation spans 1.4 decay lengths
  against 0.64 at 500 m, so it depends more on the assumed exponential. At 500 m the
  amplitude does not depend on `dt` (0.29 to 0.30 on the ghost path, 0.32 to 0.35 on the
  face path).
- On the MCD column the x-rms form (`--field xrms`: rms over x of the 2 dz part over rms
  over x of the smooth part of the window-mean `w(z, x)`, 6 to 9 h) gives 0.45 to 0.61 at
  1.25 to 2.25 km in all four records (S04, S14A, S14M1, S14M2). At CFL 0.4 it falls to
  0.25 at 3.25 km and 0.05 at 4.75 km (S14M1; S04 alike). At CFL 0.2 it stays at 0.46 to
  0.50 up to 4.25 km and is 0.23 at 6.75 km (S14M2; S14A alike). A single exponential does
  not describe these profiles (fits over levels 3 to 8 give 2.0 km for S14M1 and above
  50 km for S14M2), so no decay length is quoted for the MCD column. In the calm layer the
  smooth part is not a single stationary wave, and the ratio is poorly conditioned there.

**c. `MF(1.25)` on the MCD column.** By the model, `L` is 0.88 km at the S14M1 step and
1.76 km at the S14M2 step, so 1.25 km lies outside `L` at CFL 0.4 (1.4 `L`) and inside it
at CFL 0.2 (0.71 `L`). By the measured MCD profiles, 1.25 km lies inside the layer in which
the 2 dz part is 0.45 to 0.5 of the smooth `w` at both steps. The x-rms of the window-mean
`w` at 1.25 km rises from 7.1 to 13.6 mm/s from S14M1 to S14M2 (x 1.9), alongside the
+120 % in `MF(1.25)`.

*Note dated 2026-10-02.* In the MCD calm layer `u` is near 0 below the 4.5 km critical
level, while the imposed ground `w` is `U_bar dh/dx` with `U_bar` = 40 m/s. This forcing
has no steady inviscid solution in that layer, so `MF(1.25)` there has no converged value
to approach, and it lies within the near-ground 2 dz layer whose depth follows `dt`.
`MF(1.25)` on the MCD column is therefore a non-convergent diagnostic by construction. From
this date it is not used as a convergence criterion. It is still reported. The earlier
readings that used it (P1, the addendum 2 descriptive observation 2) are unchanged.

**d. Error budget and superseded amplitudes.** Two rows were added to the Section 11
table, not summed: the ghost-path forcing on the MCD column (-39 % at `MF(30.25)` at CFL
0.4, S14M1/S04 = 1.64) and the 2 dz boundary mode (ground amplitude about 0.30 of the
smooth `w`, penetration about `dz^2 / (2 c dt)`). A dated note at the head of Sections 11,
12 and 13, and of `arsia_stage1_U57_h3_deep_dive.md`, `ext6h_instability_mechanism.md`
and `linear_flux_constancy.md`, marks every absolute amplitude from the ghost path, which
includes every Stage 1 level, as superseded for amplitude, with a pointer to addendum 2.

#### Code for addendum 3

Each new path sits behind a switch whose default is the earlier behaviour.

- `run_mountain_benchmark_strict.jl`: `--bcface=1` selects `bottom_forcing = :face` and
  appends `_bcface` to the default output directory; the default 0 keeps the ghost path.
  The driver keeps its pinned `pbar_mode = :total` (Section 11, decision 1), so both Queney
  runs below use `:total`. It now prints the `pbar_mode`, the ground-`w` path and the
  output directory.
- `analyze_dt_convergence.py`, all reading archived output only:
  - `zalt`: `--field xrms` (ridge runs), `--fit-levels`, and the decay fit with the model
    `L` of step 1b;
  - `zflin`: `F_lin` of the uniform-wind Gaussian ridge, factor by factor, and the measured
    `F/F_lin` (V1);
  - `zequil`: `drho_L2` from the solver log and the energy and mass drift from snapshots
    (V2);
  - `zqueney`: face/ghost amplitude in the Queney case and the Srivastava Table 11 metric
    (V3).
- `run_dt_convergence_S2T95.sh`:
  - set `a3` with a driver per case (`case_driver`), so the V3 cases run the Queney driver;
  - `driver_hash` now covers the Queney driver as well;
  - `S14_REPO` points a copy of the launcher at the repository. The launch runs a copy
    kept in `~/MAGNUSP_runs/logs/S14A3.launcher.sh`, and the stamp records its hash;
  - `check` adds three items: 0.1 h of S14M1 against its archived record; the Queney driver
    on its default path against the copy kept before the option was added
    (`~/MAGNUSP_runs/logs/a3_ref/scripts/`, `dx` = 1000 m, 0.7 h, snapshots at steps 0 and
    2000); and a run of the Queney face path, which must complete.

#### Runs

All runs use `julia -t 1`, `OPENBLAS_NUM_THREADS=1` and output under `~/MAGNUSP_runs`.
They are launched in parallel with the launcher copy, `launch a3`, after `check` has passed.
Logs are `~/MAGNUSP_runs/logs/S14<id>.log`, with stamps in `S14A3.start`. `monitor a3`
holds `caffeinate -s -i` for as long as any case runs and writes the per-step cost and the
CPU power every 60 s to `S14A3.cost.tsv`. The common arguments are
`pbar=1 ramp=2 k0=0.005 lz=95` (written `C` below) unless stated.

| id | purpose | arguments | `dt` [s] | steps | s per step (a2 cost log) | est. wall | disk [GB] | output directory under `~/MAGNUSP_runs` |
|---|---|---|---|---|---|---|---|---|
| S14V1a | ridge `F/F_lin`, face, CFL 0.4 | `linear sigma=50 h=0.03 hours=18 C bcface=1` | 0.740 | 87 600 | 0.050 | 1.2 h | 3.1 | `arsia_mons_linear18h_LX2000_sig50_K0_0.005_LZ95_pbarP_ramp2h_bcface` |
| S14V1b | as V1a, CFL 0.2 | V1a plus `cfl=0.2 outint=200` | 0.370 | 175 000 | 0.050 | 2.4 h | 3.1 | `..._bcface_cfl0.2` |
| S14V2u | H = 0, face, uniform wind (Z1 configuration) | `linear mode=6 h=0 hours=1 pbar=1 bcface=1` | 0.741 | 4 860 | 0.050 | 4 min | 0.2 | `arsia_mons_linear1h_LX2000_mode6_pbarP_bcface` |
| S14V2m | H = 0, face, MCD (S14Z0 configuration) | `hours=1 C u=40 h=0 bcface=1` | 0.659 | 5 460 | 0.097 | 9 min | 0.4 | `arsia_mons_stage2_1h_K0_0.005_LZ95_pbarP_ramp2h_bcface` |
| S14V3g | Queney, ghost | Queney driver `--dx=250 --tend_hr=10` | 0.305 | 118 000 | about 0.12 | 3.9 h | 0.6 | `mountain_benchmark_strict_10h_dx250` |
| S14V3f | Queney, face | V3g plus `--bcface=1` | 0.305 | 118 000 | about 0.12 | 3.9 h | 0.6 | `mountain_benchmark_strict_10h_dx250_bcface` |
| S14S1 | stationarity, S14G1 to 30 h | `linear mode=6 h=0.005 hours=30 C bcface=1` | 0.740 | 146 000 | 0.050 | 2.0 h | 5.1 | `arsia_mons_linear30h_LX2000_mode6_K0_0.005_LZ95_pbarP_ramp2h_bcface` |
| S14S2 | MCD face, `dz` = 250 m | `hours=9 C u=40 h=0.1 bcface=1 dz=250 outint=200` | 0.328 | 98 800 | 0.184 (S14D) | 5.0 h | 7.0 | `arsia_mons_stage2_9h_K0_0.005_LZ95_pbarP_ramp2h_dz250_bcface` |
| S14S3E | launch-consistent face, CFL 0.4 | `hours=9 C u=-40 h=0.1 ublend=8 bcface=1` | 0.659 | 49 200 | 0.097 | 1.3 h | 3.5 | `arsia_mons_stage2_9h_K0_0.005_LZ95_pbarP_ramp2h_ub8_bcface` |
| S14S3F | as S3E, CFL 0.2 | S3E plus `cfl=0.2 outint=200` | 0.330 | 98 300 | 0.097 | 2.6 h | 3.5 | `..._ub8_bcface_cfl0.2` |
| S14S4 | MCD face, z then x | `hours=9 C u=40 h=0.1 bcface=1 zx=1` | 0.659 | 49 200 | 0.097 | 1.3 h | 3.5 | `arsia_mons_stage2_9h_K0_0.005_LZ95_pbarP_ramp2h_bcface_zx` |

Run directories are `run_U40.0_h<h>` (Arsia), `run_U-40.0_h0.1` (S3E, S3F) and
`run_U18_b10` (Queney). The s per step are the awake rates of the a2 batch with eight
processes (14:06 to 14:40 on 2026-10-01): 0.050 on the 2000 km grid, 0.097 on the 4000 km
grid and on the 2000 km grid at `dz` = 250 m, and 0.184 for S14D from its record. The
Queney grid (960 x 240, inviscid) has no measured rate; 0.12 scales the 0.097 by the cell
count. The S2 step is that of S14D, 0.5 % below the 0.3296 s of S14M2, which is as close as
a CFL-limited step at `dz` = 250 m comes. Estimates: about 5.0 h of wall clock, set by
S14S2, and 24 CPU hours in total; 11 processes may add 10 to 30 % of contention. Disk is
about 31 GB; the local disk has 3.4 TB free.

#### Readings fixed before launch

Windows: 6 to 9 h for the MCD runs (S2, S3, S4, as P1), 12 to 17 h for V1 (as P4), 24 to
30 h for S1. A value within 0.5 percentage points of a band edge is reported as at the
edge, and the reading of the higher band applies. Each reading is recorded against these
statements, and a failed prediction is recorded as failed.

**V1, ridge `F/F_lin` on the face path (S14V1a, S14V1b).** This is the Section 10 analytic
check (Figure 1: `h` = 30 m, `sigma` = 50 km, `U` = 40 m/s, `LX` = 2000 km, `MF_lin` =
-2.25 µPa), moved to the current settings (60 to 95 km layer with `K0` = 0.005, `pbar=1`,
2 h ramp, 18 h) and the face path.

*`F_lin`.* For `h = H exp(-x^2 / (2 sigma^2))` the hydrostatic Boussinesq drag per unit
length is `D = rho0 U N H^2` (Nappo 2002, eq. 3.121, as cited for this addendum; that page
is not among `docs/Nappo/`, so the form was checked here: `rho0 N U / pi` times the
integral of `k |h(k)|^2` over `k > 0`, with `h(k) = H sigma (2 pi)^(1/2) exp(-k^2 sigma^2 / 2)`,
gives exactly `rho0 U N H^2`). The domain-mean flux is `F_N = D / LX` in magnitude. With
`rho0` = 7.124e-3 kg m^-3 and `N0` = 0.01763 s^-1 at the first cell (0.25 km),
`F_N` = 2.2609 µPa. Each factor below is the ratio to the line above it, so their product is
the theory used (`zflin`, 139 propagating modes, `s14/a3/flin_theory_V1.json`, 4 s):

| factor | definition | value | theory after it [µPa] |
|---|---|---|---|
| discrete-mode | periodic spectrum, `2 rho0 U N0 sum_n k_n |c_n|^2` | 0.9959 | 2.2515 |
| non-hydrostatic | `m = (N0^2/U^2 - k^2)^(1/2)` | 0.99897 | 2.2492 |
| compressibility | `m^2` plus `k^2 U^2 / c0^2 - 1 / (4 H_rho^2)`, `c0` = 216.9 m/s, `H_rho` = 6.23 km | 0.9833 | 2.2116 |
| column | Taylor-Goldstein with the `N(z)`, `rho(z)` and `c(z)` of the MCD column, `Phi(0)` imposed at `z` = 0, radiation at 59.75 km, net flux over 5.25 to 15.25 km | 1.1018 | 2.4368 |
| layer, steady | the 60 to 95 km layer, radiation at 95 km, every return in place | 0.7966 | 1.9412 |
| layer, time-resolved (primary) | the layer factor of a mode applies once its round trip to the layer base and back to 15.25 km has elapsed after a 2 h ramp delay: 0.9361 of the column value over 12 to 15 h, 0.9048 over 15 to 17 h | 0.924 (window mean) | 2.2506 (window mean) |

Notes on the factors:
- The column factor contains the normalization of `Phi(0)` with `rho(0)` rather than
  `rho0` at the first cell (1.041) and the transmission through the near-ground inversion
  (1.058).
- The radiation level of the column line lies midway between the 500 m nodes of the MCD
  table. At a node the `s'` of the piecewise-linear `rho` has a spike that makes the bracket
  negative, which an earlier version of `zflin` placed at 60 km.
- The per-mode layer factor falls with wavelength (0.904 at n = 6, which is the C2R value
  of 0.911 within the band difference; 0.79 at n = 4; 0.41 at n = 1) because `K0 / (k U)`
  rises as `k` falls. 53 % of the flux is in modes n = 1 to 5, whose return reaches the
  band after 12 h. Classing them with a delay of 0 h or 3 h in place of 2 h moves the window
  theory by -1.4 % and +1.2 %.
- Every mode has arrived at 15.25 km by 10.2 h (n = 1 at 8.2 h after a 2 h delay).

*Measurement.* For each hour of 12 to 17 h, the hour-mean `|MF|` over 5.25 to 15.25 km is
divided by that hour's theory. The primary `F/F_lin` is the mean of the five hourly ratios.

*Bands and discrimination* (each of S14V1a and S14V1b):

| primary `F/F_lin` | reading |
|---|---|
| 0.96 to 1.04 | the face-path ridge closes; the +10 to 12 % C2R excess is attributed to the single-mode cavity theory or its normalization |
| 1.08 to 1.16 | an excess similar to that of C2R (+11.9 % at 500 m); attributed to the solver's boundary closure |
| any other value | mechanism unattributed |

If the classification bounds (delays of 0 and 3 h) put the value in a band other than that
of the primary, the reading is "mechanism unattributed". A second statement: the primary of
S14V1b differs from that of S14V1a by less than 1 %, the face-path `dt` independence on a
multi-mode spectrum (S14G2 against S14G1 +0.03 %, the MCD pair +0.35 %). Reported without a
band: `F/Q4` (column value), `F/Q1` (the discrete Nappo value), and `F(z)/Q4` from 1.25 to
50.25 km.

*Expected outcome.* The archived ghost-path `sigma` = 50 km records at CFL 0.4 (`:total`,
impulsive start, 12 to 15 h; `zflin` with their own tops, 44 s,
`s14/a3/flin_ghost_archived_LZ{90,120}.json`) give `F/Q4` of 0.689 (60 to 90 km layer) and
0.724 (60 to 120 km layer), and `F/F_lin` of 0.80 and 0.74 against their own time-resolved
theories. Their layer theories differ by 10.6 % (0.887 and 0.981 of `Q4`), and the
measured fluxes by 5.1 % in the same direction. Applying the uniform-wind face/ghost ratio
at CFL 0.4 (1/0.893) gives 0.77 to 0.81 of `Q4`. With the 95 km theory (0.924 of `Q4`) the
expected `F/F_lin` is about 0.86 (0.80 to 0.92), in the "mechanism unattributed" range. Under
this expectation the ridge lies about 14 % below its linear theory while the single mode lies
12 % above its own. The expectation carries the `:total` background drift and the
impulsive start of those records, so it is a weak prior.

**V2, H = 0 equilibrium on the face path (S14V2u, S14V2m, 1 h).** With `h` = 0 the face
path imposes `w_f` = 0 on the ground face in the corrector. Read with `zequil`: `drho_L2`
from the solver log at every 500-step line, and the total energy and mass from the snapshot
nearest 1 h against snapshot 0. References (`s14/a3/zequil_refs.json`, 26 s):

| quantity | Z1 (uniform wind, ghost) | S14Z0 (MCD, ghost) |
|---|---|---|
| `drho_L2` peak (time) | 1.32e-14 (2592 s) | 4.26e-5 (330 s) |
| `drho_L2` mean, 0.5 to 1 h | 1.28e-14 | 7.75e-6 |
| `drho_L2` last line | 1.31e-14 | 6.09e-6 |
| energy drift to the snapshot nearest 1 h | +8.9e-16 (3629 s) | -6.07e-7 (3560 s) |
| mass drift | -5.6e-16 | -5.5e-9 |

Criterion, for each run against its reference: the `drho_L2` peak, the 0.5 to 1 h mean and
the last value, and the magnitude of the energy drift, are each at or below the reference.
The uniform-wind values are at round-off. If one of them exceeds its reference while staying
below 1e-13 (`drho_L2`) or 1e-14 (energy), the literal criterion is recorded as failed and
the values are described as round-off. Expected: met in both. S14V2u is Z1 with
`bcface=1` and `hours=1` (60 to 90 km layer, `K0` = 0.05, no ramp, the Z1 configuration).

**V3, Queney case on both paths (S14V3g, S14V3f).** Srivastava et al. (2022) Section 3.3,
Tables 9 to 11: 240 x 60 km, `dx = dz` = 250 m, isothermal 239 K, `U` = 18 m/s, bell hill
`b` = 10 m and `a` = 20 km, top 20 km damping layer with `K0` = 0.005, Earth constants,
inviscid, CFL 0.4, 10 h, `pbar_mode = :total` (pinned, shared by both runs). Both runs use
the same CFL, and `zqueney` reports the ratio of their mean `dt`.

*Primary.* `A` = rms of `w` (face) over rms of `w` (ghost) over 0 to 10 km and
`|x - x0|` up to 100 km, each rms averaged over the snapshots in 8 to 10 h.

*Prediction from the measured ghost-path deficit.* In uniform wind the ghost path carries
0.893 of the face-path flux at CFL 0.4 (C2R against S14G1: `dt` 0.7405 s, `dz` 500 m,
`c0` 215.3 m/s), first order in `dt`. Transferring it needs a scaling variable. Two are
written, with the amplitude ratio taken as the inverse square root of the flux ratio:
- the vertical acoustic Courant number `c dt / dz`: 0.319 in C2R and 0.377 here (`c` =
  309.9 m/s, `dt` = 0.3045 s). Flux deficit 12.7 %, `A` = 1.070;
- `dt` itself, or the advective number `U dt / dz` (0.059 in C2R, 0.022 here): flux
  deficit 4.4 % or 4.0 %, `A` = 1.023 or 1.020.

| `A` | reading |
|---|---|
| 1.05 to 1.09 | consistent with the acoustic Courant-number scaling |
| 1.01 to 1.035 | consistent with a `dt` or advective scaling |
| below 1.01 | the ghost-path deficit does not transfer to this configuration |
| any other value | unattributed |

Expected: the acoustic band, with low confidence. The 2 dz decay length follows `c dt / dz`
(step 1b), while the 39 % MCD loss at a near-zero ground wind is not described by any of
the three variables. Also reported: the Table 11 metric (rms of `u' - u'_Queney` over 0 to
40 km as a percentage of `U`, last snapshot) for each path, against the 0.74 % of
Srivastava at 250 m, and `A` by snapshot over 8 to 10 h. A ghost-path match to Srivastava
is not counted as validation of the boundary, since the ghost path carries a
`dt`-dependent amplitude deficit. The face-path value is reported alongside and is not by
itself a validation either.

**S1, stationarity of the face-path level (S14S1).** The S14G1 configuration to 30 h. Read
with `compare_sponge_reflection.py --cavity`, hourly.
- Criterion: the 24 to 30 h mean level lies within 1.019 +- 1.5 % (1.004 to 1.034), and the
  hourly range over 24 to 30 h (largest minus smallest of the six hourly levels) is at most
  3 % of the mean.
- Readings: both met, the 12 to 17 h level with its ringing stands; mean inside and range
  wider, ringing persists to 30 h; mean outside, the 12 to 17 h level was not stationary
  (the drift is then quoted).
- Also registered: S14S1 reproduces S14G1 bit for bit on every snapshot both records hold
  (0 to 18 h), since only `t_end` differs.
- Expected: both met, with a mean of 1.01 to 1.03 (the 12 to 17 h midpoint of trough and
  peak was 1.018).

**S2, MCD face path at `dz` = 250 m (S14S2 against S14M2).**
- Criteria: the `MF(30.25)` change below 5 %, and the `R50/30` and `R10/30` changes below
  5 %. Met: the face-path Stage 2 flux at `dz` = 500 m is resolved to within 5 % in `dz`.
- Expected: met, with `MF(30.25)` about -2 % (S14R1 against S14G2: -1.85 %) and the ratios
  within 1 %. For reference, the ghost path failed the same test (P6, -6.4 %).
- `MF(1.25)` is reported without a criterion (step 1c).

**S3, launch-consistent column on the face path (S14S3E, S14S3F).** The S14E and S14F
configurations (`ublend=8`, `u=-40`) with `bcface=1`. The prediction if the calm layer
explains the `MF(1.25)` residual: `MF(1.25)`, `MF(10.25)` and `MF(30.25)` each change by
less than 5 % per halving.

| outcome | reading |
|---|---|
| all three below 5 % | consistent with the calm layer as the source of the `MF(1.25)` residual on the MCD column |
| `MF(1.25)` at or above 5 %, the other two below 5 % | the residual is not explained by the calm layer alone |
| `MF(10.25)` or `MF(30.25)` at or above 5 % | the face path is not `dt`-independent on this column |

Expected: all three below 5 % (uniform wind on the face path: 0.03 to 0.06 %; the ghost
path on this column, P5: +12 to +14 %).

**S4, MCD face path with the z-sweep first (S14S4 against S14M1), registered null.**
- Null: `MF(30.25)` and `MF(50.25)` each within 1 % of S14M1.
- A failure would show a splitting term on the sheared column at CFL 0.4 (the shear part of
  H_S), which the uniform-wind S14X1 could not see.
- `MF(10.25)` and `MF(1.25)` are reported without a criterion.
- Expected: met (S14X1 in uniform wind: 0.05 % or less).

#### Check before launch

`check` (72 s, driver, Queney driver and `src/` hash `98f7fe52db0bfcf7`) reproduced with
`max |diff|` exactly 0.0 in every shared snapshot: 0.1 h of S04 (6 snapshots), C2R (5),
S14G1 (5), C2R with `probe=5` (5), the new item S14M1 (6), and the Queney driver on its
default path against its pre-addendum copy (2). The Queney face path completed.

Two earlier failures of `check` led to two corrections, recorded here:
- `get_arg` in the Queney driver returned a `SubString`, which `run_simulation` rejects
  for `outdir`. So `--outdir=` had never worked, and every earlier record used the default
  path. Both the driver and its reference copy now return `String(v)`, so the bit check
  still isolates the `bcface` switch.
- The Queney driver does not create its output directory (NetCDF "Permission denied"). The
  launcher now creates every case directory before it starts the case, and the check
  directories before the check.

#### Launch record

`launch a3` at 2026-10-02 15:02:37 from the copy `~/MAGNUSP_runs/logs/S14A3.launcher.sh`
(launcher hash `944bfc1a90f9d45e`, driver hash `98f7fe52db0bfcf7`), all eleven cases, with
`S14_REPO` set to the repository. `monitor a3` was started from the same copy at 15:03 and
holds `caffeinate -s -i` (PreventSystemSleep and PreventUserIdleSystemSleep asserted). The
launcher in `scripts/` is not edited while the batch runs.

#### Results, H = 0 equilibrium on the face path (S14V2u, S14V2m)

Both completed with exit 0: S14V2u 4 859 steps in 252.0 s (`real`, 0.052 s per step), S14V2m
5 462 steps in 602.3 s (0.110 s per step, with eleven processes running). Scoring took 22 s
(`zequil`, `s14/a3/zequil_V2.json`).

| quantity | Z1 (ghost) | S14V2u (face) | outcome | S14Z0 (ghost) | S14V2m (face) | outcome |
|---|---|---|---|---|---|---|
| `drho_L2` peak (time) | 1.32e-14 (2592 s) | 1.52e-14 (2963 s) | above | 4.260e-5 (330 s) | 4.260e-5 (330 s) | at |
| `drho_L2` mean, 0.5 to 1 h | 1.282e-14 | 1.206e-14 | below | 7.746e-6 | 7.734e-6 | below |
| `drho_L2` last line | 1.31e-14 | 1.21e-14 | below | 6.09e-6 | 6.08e-6 | below |
| energy drift to the snapshot nearest 1 h | +8.9e-16 | +1.3e-15 | above | -6.068e-7 | -6.686e-7 | above (by 10 %) |
| mass drift (no criterion) | -5.6e-16 | -4.4e-16 | | -5.5e-9 | -7.2e-8 | |

**Against the criterion** ("each at or below the reference"):
- S14V2u: failed on the literal criterion. The peak and the energy drift exceed the Z1
  values. Both are below the round-off thresholds stated in advance (1e-13 and 1e-14), so,
  as registered, the values are described as round-off. In both paths the uniform-wind H = 0
  column stays a fixed point of the scheme to round-off for 1 h.
- S14V2m: failed on the energy drift. Its magnitude is 10 % above S14Z0. `drho_L2` matches
  S14Z0 to three or four digits at every logged step (4.26e-5 at 330 s, 6.08e-6 at the end),
  so the initial adjustment of the sheared viscous column does not depend on the path.

**Observation, not registered.** The face path changes the 1 h domain mass by -7.2e-8
against -5.5e-9 on the ghost path (13 times larger, still below 1e-7). With `w_f` = 0
imposed in the corrector, the mass flux through the ground face is zero by construction,
so the difference is expected elsewhere: in the face density flux of the corrector, the
ghost rows rebuilt after the x-sweep, or the top boundary. It is not identified here. It
is about 1e-4 of the 0.72 % mass loss that `:total` produced (Section 10) and is not
expected to affect the fluxes. It is carried as an open item.

*Note dated 2026-10-19, mass and energy audit of S14V2m against S14Z0 (descriptive, zero
cost).* This note addresses the open item above. The observation and the registered
outcome are not changed. Tool: `analyze_dt_convergence.py zequil --budget` (2 s, local
copies of both records, `s14/a3/zequil_V2_budget.json`), 55 snapshots per record (every
100 steps) to 3560 s, time integrals by the trapezoid rule. All terms are relative to the
initial mass.

| term | S14Z0 (ghost) | S14V2m (face) |
|---|---|---|
| sponge source of density, integral of `sum K(z) (rho_bar - rho)` | -4.227e-8 | -4.216e-8 |
| outflow through the top cell, integral of `-sum rho w / dz` | -3.064e-8 | -3.067e-8 |
| sum of the two | -7.291e-8 | -7.283e-8 |
| measured mass drift | -5.53e-9 | -7.221e-8 |
| remainder, drift minus the sum | +6.74e-8 | +6.2e-10 |

- The top and sponge terms agree between the paths to 0.3 %. On the face path they
  account for the drift to within 1 % (remainder 6e-10), which is consistent with zero
  mass flux through the ground face. On that face the corrector uses `rho_f w_f` and
  `w_f (E_f + P_f)` with `w_f = 0`, so both the mass flux and the energy flux there are
  zero whatever density and pressure the predictor and the post-x ghost rebuild give it.
- On the ghost path the drift lies 6.7e-8 above the top and sponge terms, a net inflow
  of that size. The per-level change of face minus ghost, cumulated from the ground,
  reaches half of its total (-6.67e-8) by 6.25 km and 99 % by 44.25 km. Above 60 km the
  per-level difference is below 6e-12. The difference grows linearly in time (-9.7e-9 at
  396 s, -3.8e-8 at 1978 s, -6.7e-8 at 3560 s), as a steady boundary flux does, not as
  part of the initial adjustment.
- Mechanism on the ghost path. At `h = 0` the `mountain` forcing still calls
  `impose_ghost_w!` (`src/forcing.jl`), which sets the ghost velocity `w_g = -w_p`
  rather than the ghost momentum `-(rho w)_p`. In the stratified column the ghost
  background density exceeds that of the first cell, so the predictor half state of
  eq. 10a at the ground face carries
  `rho w = w_p (rho_bar_p - rho_bar_g) / 2 - dt g rho'_p / 2`, which is not zero, and the
  corrector passes it as a mass flux. The x-sum of `w_p` near the ground is negative over
  the hour in both records (the column settles), which gives an inflow. A snapshot
  estimate of the two terms is +1.13e-7 for S14Z0: the same sign and order as the
  +6.7e-8 inferred, and 1.7 times larger. The estimate samples the flux every 66 s and
  extrapolates `rho_bar_g` log-linearly, so it is an order-of-magnitude check only.
- Energy. The energy-drift difference, face minus ghost, is -6.19e-8 of the initial
  energy, which is the 10 % of the V2m reading. It has the same signature: it accumulates
  from the ground (half by 7.25 km, 99 % by 42.75 km), stays below 4e-13 per level above
  70 km, and grows linearly in time. This is consistent with the same ground-face flux
  carrying energy into the ghost-path column. An enthalpy accounting of that flux was not
  attempted.

Reading (descriptive). The 13 times larger mass drift and the 10 % larger energy drift of
S14V2m do not come from the face path. Of the candidates named above, the top outflow is
the same on both paths, and the face density flux of the corrector and the post-x ghost
rebuild add no mass, because the face-path budget closes without a ground term. The
difference enters through the ground face of the ghost reference S14Z0, whose inflow
offsets part of the top and sponge loss and so makes its drift smaller. The S14Z0 energy
drift is therefore a reference that contains this inflow. The V2m outcome stands as
recorded (failed on the energy drift). Section 15 proposes S14V2m as the equilibrium
reference for face-path runs on the MCD column. That proposal is not applied
retroactively.

#### Results, Queney case on both paths (S14V3g, S14V3f)

Both completed with exit 0 at 36 000 s: S14V3g 118 430 steps in 1 995.6 s (`real`, 0.017 s
per step), S14V3f 118 374 steps in 2 064.3 s. The inviscid Earth grid has no implicit
diffusion solve and runs about seven times faster per step than the 0.12 s estimated. Mean
`dt` 0.30398 and 0.30412 s (ratio 1.00047). Scoring took 2 s (`zqueney`,
`s14/a3/queney_V3.json`).

| quantity | S14V3g (ghost) | S14V3f (face) | face / ghost |
|---|---|---|---|
| rms `w`, 0 to 10 km, `|x - x0|` to 100 km, mean over 8 to 10 h [mm/s] | 3.7588 | 3.9463 | **1.0499** |
| regression `sum(w_f w_g) / sum(w_g^2)`, mean over 8 to 10 h | | | 1.0234 |
| Table 11 metric, last snapshot (rms of `u' - u'_Queney` over 0 to 40 km, % of `U`) | 2.494 % (9.96 h) | 2.495 % (9.97 h) | |

By snapshot, `A` stays within 1.0495 to 1.0502 over 8.1 to 10.0 h, and the regression
ratio within 1.0228 to 1.0238, so both are stationary.

**Against the registered bands.** `A` = 1.0499 lies 0.01 points below the 1.05 edge of the
acoustic band. Under the edge rule it is reported at the edge, and the reading of the higher
band applies: consistent with the acoustic Courant-number scaling. The value is 2 points
below the 1.070 that scaling predicted. As a flux deficit it is 9.3 %, against the
predicted 12.7 % and the 10.7 % measured in uniform wind on the Arsia column, both also at
CFL 0.4. The `dt` and advective scalings (deficits of 4.4 and 4.0 %, `A` of 1.023 and
1.020) are not supported by the primary.

**Observations, not registered.**
- The regression ratio (1.023) is lower than the amplitude ratio because the two fields
  differ in phase. If the difference were a pure phase shift it would be about 13 deg, the
  same sign as the +7.4 deg face-`w` lag of the S14P1 probe.
- The ghost-path deficit is of the same size on two columns whose acoustic Courant
  numbers differ by 18 % and whose `dt` differ by a factor 2.4 (9.3 % here, 10.7 % on the
  Arsia column, both at CFL 0.4). This is consistent with a deficit set mainly by the CFL
  number of the run rather than by `dt` itself.
- Neither path reproduces the 0.74 % of Srivastava at 250 m: 2.49 % on both, differing by
  0.001 points. The Table 11 metric is therefore dominated by a term that both paths share
  (for example the `:total` background drift, the impulsive start or the reference
  solution's hydrostatic, Boussinesq form), not by the ground boundary. As registered, the
  ghost-path value is not counted as a validation of the boundary, and neither is the
  face-path value.

#### Results, MCD face path with the z-sweep first (S14S4)

Completed with exit 0: 49 153 steps, 4 884.0 s (`real`, 0.099 s per step), ended 16:24:18.
Scoring with `suite` against S14M1 took 70 s together with S14S3E and S14E
(`s14/a3/s34/s14_suite.json`). 6 to 9 h means, change against S14M1:

| quantity | S14M1 | S14S4 | change |
|---|---|---|---|
| `MF(30.25)` [µPa] | 0.240388 | 0.240308 | -0.033 % |
| `MF(50.25)` [µPa] | 0.165587 | 0.165602 | +0.009 % |
| `MF(10.25)` (no criterion) | 0.237743 | 0.237690 | -0.022 % |
| `MF(1.25)` (no criterion) | -7.1634 | -7.1635 | -0.002 % |
| `R10/30`, `R50/30` | 0.9890, 0.6888 | 0.9891, 0.6891 | +0.011 %, +0.043 % |

The registered null is met: both primary fluxes are within 0.04 % of S14M1 (bound 1 %), and
within 0.04 % in each hour. On the sheared MCD column at CFL 0.4, reversing the order of the
dimensional splitting changes the face-path flux by no more than in uniform wind (S14X1).
The shear term of H_S is not measurable at this size.

#### Results, ridge `F/F_lin` on the face path, CFL 0.4 (S14V1a)

Completed with exit 0: 87 523 steps, 4 403.4 s (`real`, 0.050 s per step). Scoring took 6 s
(`zflin`, `s14/a3/flin_V1a.json`).

| hour | measured `|F|` over 5.25 to 15.25 km [µPa] | theory [µPa] | `F/F_lin` |
|---|---|---|---|
| 12 to 13 | 2.4961 | 2.2811 | 1.0942 |
| 13 to 14 | 2.5097 | 2.2811 | 1.1002 |
| 14 to 15 | 2.5234 | 2.2811 | 1.1062 |
| 15 to 16 | 2.5315 | 2.2048 | 1.1482 |
| 16 to 17 | 2.5495 | 2.2048 | 1.1563 |

- **Primary `F/F_lin` = 1.121.** Classification bounds 1.108 (0 h) and 1.136 (3 h), so all
  three lie in the band 1.08 to 1.16.
- Against the factors: `F/Q4` (column) 1.035, `F/Q3` (uniform `N0`, compressible) 1.140,
  `F/Q1` (discrete Nappo) 1.120, `F/Q5` (steady layer) 1.299.
- `F(z)/Q4` is 0.991 at 1.25 km, 1.043 at 4.75 km, 1.032 at 10.25 km, 1.023 at 20.25 km,
  0.989 at 30.25 km, 0.958 at 40.25 km and 0.918 at 50.25 km.

**Registered reading (S14V1a).** The face-path ridge shows an excess over its linear theory
(+12.1 %) similar to that of the single mode (+11.9 %). Under the discrimination fixed before
launch, the excess is attributed to the solver's boundary closure rather than to the
single-mode cavity theory or its normalization. The second statement (S14V1b within 1 %)
is scored when S14V1b ends.

**Expected outcome, failed.** The expected value was about 0.86 (0.80 to 0.92), from the
archived ghost-path `:total` records. The measured 1.121 is 1.45 to 1.50 times their `F/Q4`
(1.035 against 0.689 and 0.724), where the face/ghost ratio of uniform wind (1/0.893) was
assumed. A factor of 1.45 is larger than any ghost-path loss measured at CFL 0.4 on a
uniform column. The remainder lies in the `:total` background, the impulsive start or the
ridge spectrum of those records, and this addendum does not separate them.

**Observations, not registered.**
- The measured band flux rises steadily by about 0.5 % per hour over 12 to 17 h. It does
  not drop when the theory classes modes 4 and 5 as returned at 15 h (the theory falls by
  3.3 %). The part of the excess that grows from 1.094 to 1.156 comes from that drop. The
  measured flux does not follow the layer factor predicted for the long modes, in line with
  the archived records, whose fluxes differed by half as much as their layer theories.
- Measured against the column value without any layer factor (`Q4`), the ridge is within
  4.3 % of theory at 4.75 to 20.25 km. Measured against the steady layer theory, it is
  30 % above. The size of the excess therefore depends on the layer factor, which is the
  least certain line of the theory. For mode 6 the same comparison gives 1.019 against the
  column and 1.119 against the layer.
- Above 20 km `F(z)/Q4` falls (0.918 at 50.25 km), consistent with the slower modes still
  filling the upper column at 12 to 17 h.

#### Results, ridge `F/F_lin` at CFL 0.2 (S14V1b)

Completed with exit 0: 175 044 steps, 8 664.4 s (`real`, 0.049 s per step). Scoring took 9 s
for both runs (`zflin`, `s14/a3/flin_V1ab.json`).

| quantity | S14V1a (CFL 0.4) | S14V1b (CFL 0.2) | change |
|---|---|---|---|
| primary `F/F_lin` | 1.1210 | 1.1220 | +0.09 % |
| classification bounds | 1.108 to 1.136 | 1.109 to 1.137 | |
| `F/Q4` | 1.0350 | 1.0359 | +0.09 % |
| `F(z)/Q4` at 1.25 km | 0.991 | 0.906 | -8.6 % |
| `F(z)/Q4` at 10.25, 30.25, 50.25 km | 1.032, 0.989, 0.918 | 1.033, 0.989, 0.919 | within 0.1 % |

Both registered statements are met for S14V1b: the primary lies in the band 1.08 to 1.16
(the same reading as S14V1a, boundary closure), and it differs from S14V1a by 0.09 %, below
1 %. On a multi-mode spectrum the face path is `dt`-independent above the ground layer, as
for the single mode. The one level that moves is 1.25 km (-8.6 %). It lies outside the
model `L` at CFL 0.4 (0.78 km) and inside it at CFL 0.2 (1.57 km), in uniform wind, which
places part of the near-ground `dt` dependence of `MF(1.25)` in the 2 dz layer also where
there is no calm layer (descriptive).

#### Results, stationarity of the face-path level (S14S1)

Completed with exit 0: 145 967 steps to 30 h, 7 212.1 s (`real`, 0.049 s per step). Scoring
took 11 s: `compare_sponge_reflection.py --cavity` (6 s, `s14/a3/s1/`,
`s14/a3_s1_cavity.log`) and `bitcheck` against S14G1 (5 s, `s14/a3_s1_bitcheck.log`).

| hour | 12 | 13 | 14 | 15 | 16 | 17 | 18 | 19 | 20 | 21 | 22 | 23 | 24 | 25 | 26 | 27 | 28 | 29 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| level | 1.021 | 1.002 | 1.009 | 1.028 | 1.035 | 1.021 | 1.012 | 1.016 | 1.022 | 1.021 | 1.017 | 1.015 | 1.020 | 1.022 | 1.018 | 1.019 | 1.022 | 1.019 |

- 24 to 30 h mean level 1.0200, inside 1.004 to 1.034. Hourly range 0.0041, or 0.40 % of
  the mean (limit 3 %). Both criteria are met: the 12 to 17 h level with its ringing
  stands. The 12 to 17 h mean of this record is 1.0189, and the 24 to 30 h mean is 0.1 %
  above it.
- The ringing decays after 17 h: from 18 to 29 h every hourly level lies within 1.012 to
  1.022. The +0.33 % per hour trend fitted over 12 to 17 h (addendum 2, step 1) does not
  continue.
- `Psi0/imp` stays at 1.065 to 1.072 throughout.
- S14S1 reproduces S14G1 with `max |diff|` 0.0 on all 876 snapshots both records hold
  (0 to 18 h). Met.

The expected outcome (both met, mean 1.01 to 1.03) holds.

#### Results, launch-consistent column on the face path (S14S3E, S14S3F)

Completed with exit 0: S14S3E 49 153 steps, 4 857.0 s (`real`, 0.099 s per step); S14S3F
98 304 steps, 9 452.9 s (0.096 s per step). Scoring took 7 s (`suite`, `s14/a3/s3/`).
6 to 9 h means (µPa):

| quantity | S14S3E (0.659 s) | S14S3F (0.330 s) | change per halving |
|---|---|---|---|
| `MF(1.25)` | 21.189 | 19.762 | **-6.73 %** |
| `MF(10.25)` | 1.5148 | 1.5191 | +0.28 % |
| `MF(30.25)` | 0.8717 | 0.8719 | +0.02 % |
| `MF(50.25)` | 0.40329 | 0.40356 | +0.07 % |
| `R10/30`, `R50/30` | 1.738, 0.4627 | 1.742, 0.4629 | +0.26 %, +0.05 % |

By hour the `MF(1.25)` change grows: -2.6 % (5 to 6 h), -4.8 %, -6.8 % and -8.5 % (8 to 9 h).
`MF(10.25)` changes by +0.56 to +0.23 %, and `MF(30.25)` by 0.00 to +0.03 %.

**Registered reading.** `MF(1.25)` changes by 6.7 %, at or above 5 %, while `MF(10.25)` and
`MF(30.25)` change by less than 0.3 %. The outcome is the second row: the `MF(1.25)`
residual is not explained by the calm layer alone. The expected outcome (all three below
5 %) fails.

**Observations, not registered.**
- Here the ground wind is -40 m/s, there is no calm layer, and the forcing is
  launch-consistent. `MF(1.25)` still depends on `dt`, by -6.7 % against +120 % on the MCD
  column, and with the opposite sign. With S14V1b (-8.6 % at 1.25 km in uniform wind), this
  is consistent with a common near-ground `dt` dependence tied to the 2 dz layer, to which
  the calm layer adds most of the MCD amplitude.
- The face/ghost ratios at CFL 0.4 (S14S3E against S14E) are 1.465, 1.282, 1.243 and 1.199
  at 1.25, 10.25, 30.25 and 50.25 km. On this column the ghost-path loss at 30.25 km is
  20 %, between the 10.7 % of uniform wind and the 39 % of the MCD column.

#### Results, MCD face path at `dz` = 250 m (S14S2)

Completed with exit 0: 98 761 steps, 18 322.0 s (`real`, 0.186 s per step), ended 20:08:00.
Mean `dt` 0.32804 s against 0.32959 s for S14M2 (-0.47 %). Scoring took 11 s (`suite`,
`s14/a3/s2/`). 6 to 9 h means (µPa):

| quantity | S14M2 (`dz` 500 m) | S14S2 (`dz` 250 m) | change |
|---|---|---|---|
| `MF(30.25)` | 0.24122 | 0.22274 | **-7.66 %** |
| `R10/30` | 0.9919 | 0.9404 | **-5.19 %** |
| `R50/30` | 0.6883 | 0.7009 | +1.82 % |
| `MF(10.25)` (no criterion) | 0.23927 | 0.20947 | -12.5 % |
| `MF(50.25)` (no criterion) | 0.16603 | 0.15611 | -6.0 % |
| `R55/50` (no criterion) | 0.3618 | 0.3851 | +6.4 % |
| `MF(1.25)` (no criterion, step 1c) | -15.789 | -2.893 | magnitude -82 % |

By hour the `MF(30.25)` change is -5.8, -7.8 and -8.8 % (6 to 7, 7 to 8 and 8 to 9 h), and
the `R10/30` change -5.9, -4.6 and -5.3 %.

**Registered reading.** `MF(30.25)` changes by 7.7 %, above 5 %: failed. `R10/30` changes by
5.19 %, within 0.5 points of the 5 % edge, so it is reported at the edge, and the reading of
the higher band (at or above 5 %, failed) applies. `R50/30` (+1.8 %) is met. The face-path
Stage 2 flux at `dz` = 500 m is therefore not resolved to within 5 % in `dz` on the MCD
column. The expected outcome (met, `MF(30.25)` about -2 %) fails. The size is close to that
of the ghost-path P6 test (-6.4 %), so the face path removes the `dt` dependence of the MCD
flux but not its `dz` dependence. That dependence is four times the uniform-wind value
(S14R1: -1.85 %), and it grows through the window.

**Observations, not registered.**
- The change is largest at 10.25 km (-12.5 %) and falls with height (-7.7 % at 30.25 km,
  -6.0 % at 50.25 km). The resolution sensitivity therefore enters in the lower column,
  where the calm layer and the 4.5 km critical level lie.
- Halving `dz` at this `dt` shrinks the model `L` from 1.76 to 0.44 km, and `MF(1.25)` falls
  in magnitude by 82 %. With the S14S3 and S14V1b readings, this is consistent with the
  near-ground 2 dz layer feeding flux through the calm layer into the transmitted wave, by
  an amount that depends on the layer depth. The two refinements act in opposite
  directions: a smaller `dt` deepens the layer, a smaller `dz` makes it shallower, and the
  transmitted flux above 10 km follows `dz` but not `dt` (+0.35 % per halving of `dt`). The
  mechanism is not identified.

#### Summary of addendum 3

| test | registered reading | outcome |
|---|---|---|
| V1a, ridge `F/F_lin`, face, CFL 0.4 | 0.96 to 1.04 closes (cavity theory); 1.08 to 1.16 excess (boundary closure); else unattributed | 1.121 (bounds 1.108 to 1.136): excess similar to C2R, attributed to the boundary closure |
| V1b, same at CFL 0.2, and within 1 % of V1a | as V1a; change below 1 % | 1.122, same reading; +0.09 %, met |
| V1 expected outcome | about 0.86, unattributed | failed |
| V2u, H = 0, uniform wind | each at or below Z1 | failed on the literal criterion (peak, energy), values at round-off |
| V2m, H = 0, MCD | each at or below S14Z0 | failed on the energy drift (+10 %); `drho_L2` met |
| V3, Queney face/ghost amplitude | 1.05 to 1.09 acoustic; 1.01 to 1.035 `dt` or advective | 1.0499, at the edge, acoustic Courant-number reading |
| S1, level 24 to 30 h | mean 1.004 to 1.034, range at most 3 % | 1.0200, 0.40 %: met; S14G1 reproduced bit for bit |
| S2, MCD face `dz` = 250 m | `MF(30.25)`, `R50/30`, `R10/30` within 5 % | -7.66 %, +1.82 %, -5.19 % (edge): failed |
| S3, launch-consistent face pair | all below 5 %: calm layer explains `MF(1.25)` | `MF(1.25)` -6.73 %, others below 0.3 %: not explained by the calm layer alone; expected outcome failed |
| S4, MCD z then x | null within 1 % | -0.03 % and +0.01 %: met |

**What this changes.**
- The face path is `dt`-independent above the near-ground layer in every configuration
  tested: uniform wind, ridge and single mode, the MCD column and the launch-consistent
  column (0.02 to 0.35 % per halving at 10 to 50 km). The order of the dimensional
  splitting does not matter on the sheared column either.
- The excess over linear theory (+12 %) appears for both a single mode and a ridge
  spectrum. Under the discrimination fixed before launch it is attributed to the solver's
  boundary closure. This attribution rests on the layer factor of the theory. Measured
  against the column theory without the layer, the ridge is within 4 % from 4.75 to 20 km
  and the single mode within 2 %, so the excess could equally be read as the solver's layer
  returning less than the theory says. That alternative was not separated by this design.
- On the MCD column the Stage 2 flux at `dz` = 500 m carries a resolution term of about
  -8 % at 30 km (-12 % at 10 km) that does not converge within the window. Absolute Stage 2
  fluxes at `dz` = 500 m need this row in their error budget, or a `dz` = 250 m record.
- `MF(1.25)` depends on `dt` on every column, through the near-ground 2 dz layer, and
  strongly where the calm layer adds to it.
- The ghost-path loss at CFL 0.4 is 9 to 11 % in uniform wind on two columns (Arsia,
  Queney), 20 % on the launch-consistent column and 39 % on the MCD column.

**Wall clocks.** The batch ran from 15:02:37 to 20:08:00 on 2026-10-02, 5.09 h of wall
clock, set by S14S2. The CPU time summed over the `real` values is 17.4 h for the eleven
runs (62 709 s), against the 24 h estimated: the Queney runs took 0.017 s per step instead
of 0.12. `pmset -g log` records no sleep between launch and the end of the batch.
`caffeinate -s -i` was held by `monitor a3` and released at the end. Cost history:
`~/MAGNUSP_runs/logs/S14A3.cost.tsv`.

**Archive.** `archive a3` (1 370 s over the network mount) copied all eleven run
directories to `output/` under the names in the run table. Every file count is equal and
there are zero checksum differences (495, 60, 60, 878, 493, 1 460, 878, 493, 493, 56 and 50
files). The logs, `S14A3.start`, `S14A3.cost.tsv`, the launcher copy and the check log are
in `output/reflection_tests/logs/`. The pre-addendum Queney driver and launcher are in
`output/reflection_tests/logs/a3_ref/`, and the scoring outputs are in
`output/reflection_tests/s14/a3/`. The local copies in `~/MAGNUSP_runs` are kept.

#### Proposals awaiting approval (not launched)

**Section 15, `h` = 1 km on the face path.** MCD column, `u=40 h=1 pbar=1 ramp=2 k0=0.005
lz=95 bcface=1`, 24 h, two runs:
- F1 at `dz` = 500 m, CFL 0.4 (about 131 000 steps, about 3.5 h, 9 GB);
- F2 at `dz` = 250 m, CFL 0.4 (`dt` about 0.33 s, about 263 000 steps, about 14 h, 18 GB).

Because the face path is `dt`-independent above the ground layer, the pair isolates `dz`.
It is the resolution pair that S2 shows is needed.

Diagnostics, hourly:
- the first-overturning time and height (first `min_x u <= 0` relative to the local
  background, and first `d theta/dz < 0` below 60 km);
- the Lindzen saturation height for comparison;
- `MF(10.25)`, `MF(30.25)` and `MF(50.25)`, and `R50/30`;
- the linearity ratio `MF(h = 1) / (100 MF(h = 0.1))` against S14M1 over the pre-collapse
  window.

The pre-collapse window is the last 3 h before the first overturning below 55 km, or 18 to
24 h if none occurs. Registered criterion for the pair: the `MF(30.25)` change below 5 % in
that window. Otherwise the h = 1 km fluxes are quoted with the `dz` row.

**Stage 2 rerun on the face path.** The nine main-grid cases (`U` 20, 40, 57 by `h` 0.1, 1,
3) with `hours=24 pbar=1 ramp=2 k0=0.005 lz=95 bcface=1` at CFL 0.4, plus a CFL 0.2 twin of
U40 h0.1 scored over 18 to 24 h. Registered reading for the twin: `MF(30.25)` and
`MF(50.25)` change by less than 1 %. Cost about 3.6 h per case and 7.3 h for the twin in
parallel, about 90 GB. The h = 3 km cases aborted at the top under the ghost path
(Section 13), and the same is expected here. Recommendation, from S2: either add a
`dz` = 250 m twin of U40 h0.1 (about 14 h, 18 GB), or quote every Stage 2 flux with the S2
resolution row (-8 % at 30 km, -12 % at 10 km, growing through 6 to 9 h).
