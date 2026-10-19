# Arsia Mons Stage 1 deep dive, run U=57 m/s, h=3 km

*Note dated 2026-10-02.* Every run in this report imposes the ground `w` through the ghost cells (`bottom_forcing = :ghost`, the only path before 2026-09-30). Absolute flux levels and amplitudes quoted here are superseded for amplitude: at CFL 0.4 the ghost path loses 10.7 % of the transmitted flux in uniform wind against the face path, by an amount first order in `dt` (`reflection_convergence_tests.md` Section 14, addendum 2 and addendum 3 step 1). Shapes, timings and ratios between runs at equal `dt` are not covered by this note.

Internal analysis report, 2026-07-13. Companion to `docs/TODO_tracker.md`
(Stage 1 results section and Decisions Log rows dated 2026-06-30). Scope per
scoping decisions: Stage 1 only, sponge-free region (z < 60 km) only, full
new-diagnostics set, report in `docs/`.

## 1. Purpose and scope

`run_U57.0_h3.0` is the strongest-response corner of the Stage 1 sweep
(steering wind U_bar = 57 m/s, effective hill height h_max = 3 km,
Nh/U = 0.95). The sweep composite (`sweep_composite_w.png`) shows this panel
dominating the 3x3 grid, with small-scale striations through the AMEC band.
This report (a) inventories the data products for this case, (b) examines the
wave train and its timing, answering when the wave was forced and why the
station oscillations end near 4000 s while the production wavelet shows a
feature near minute 165, (c) characterizes the response in the highest
sponge-free altitudes (53-60 km) relevant to the planned rarefied-atmosphere
(DSMC) coupling, and (d) records what Stage 2 is and the next tasks.

The 60 km sponge base is low relative to the intended coupling altitude
(~120 km); domain extension is deferred (Section 7). All statements here are
restricted to z < 60 km; fields above the sponge base are shaped by the
damping coefficient rather than by the resolved dynamics and are excluded
from metrics.

## 2. Run configuration

From `scripts/sweep_arsia_mons_stage1.jl` and the tracker baseline row
(2026-06-30):

| Item | Value |
|---|---|
| Domain | LX = 2000 km, LZ = 90 km (dx = 4 km, dz = 500 m; nx = 500, nz = 180) |
| Background | Real MCD MY31 column (`external/mcd_profiles/mars_arsia_mons_my31.csv`) |
| Wind | Uniform eastward u_init = 57 m/s (Stage 1 idealization) |
| Forcing | Flat-bottom orographic BC, w = U_bar dh/dx, Gaussian hill sigma = 100 km, h_max = 3 km, centered at x = 1000 km |
| Horizontal BC | Outflow (Dirichlet inflow at the left edge) |
| Sponge | 60-90 km, K0 = 0.05 |
| Integration | CFL = 0.4, t_end = 10800 s (3 h), viscosity on, cooling off |
| Outcome | Converged (final_t = 10794.6 s, 17510 steps), wall-clock 1622 s |

Nonlinearity proxy Nh/U = 0.018 x 3000 / 57 = 0.95, the closest case in the
grid to the linear-theory overturning threshold Nh/U ~ 1.

## 3. Data inventory

Raw solver output, `output/arsia_mons_stage1/run_U57.0_h3.0/`:

- 176 NetCDF snapshots (`snapshot_0000000.nc` ... every 100 steps), t = 0 to
  10794.6 s, median cadence 62.65 s. Total 636 MB.
- Each snapshot: interior grid only (z = 180, x = 500; ghost cells not
  written), variables rho, u, w, P, T plus x, z coordinates, global
  attributes time_s, step, gamma = 1.31, g = 3.72.
- The 62.65 s cadence sets a Nyquist frequency of 7.98 mHz. This is why every
  wavelet spectrogram for this run is empty above ~8 mHz; frequencies above
  that are unresolvable from the snapshot record (the solver timestep,
  ~0.62 s, resolves them internally but they are not archived).

Sweep-level analysis (from `scripts/analyze_mars_mountain_sweep.py`), in
`output/arsia_mons_stage1_analysis/run_U57.0_h3.0/`:

- `snapshot_latest.png` (final-time w and density-scaled w), `station_timeseries_w.png`
  (4 stations at z = 49.75 km, x_rel = 50/150/298/498 km), `wavelet_w.png`
  (station S3, first 120 min trimmed), `mean_flow_evolution.png`
  (u_bar(z,t) - 57), `metrics.json`.

New deep-dive products (this report; `scripts/analyze_arsia_stage1_deepdive.py`,
renamed from `analyze_arsia_U57h3_deepdive.py` on 2026-08-31 when it was
generalized to any Stage 1 run; analysis wall-clock 16.5 s):

- `deepdive_station_timeseries_w_z57.png` - stations at z = 57.25 km
  (highest sponge-free band) overlaid on z = 49.75 km.
- `deepdive_wavelet_w_full_z50.png`, `deepdive_wavelet_w_full_z57.png` -
  full-record (untrimmed) wavelets at S3.
- `deepdive_w_profiles.png` - hill-window max |w|(z) and RMS w(z) at
  t = 0.5/1/2/3 h, z < 60 km.
- `deepdive_mf_ubar_profiles.png` - momentum flux MF(z) = <rho u'w'>_x and
  mean-flow change <u>_x - U_bar at the same times.
- `deepdive_convergence.png` - added 2026-08-31 with the generalized script;
  AMEC-band hill mean |w|, MF at 30/45/52 km and the peak mean-flow deviation
  below the sponge base, all against time.
- `deepdive_metrics.json` - inventory and the quantitative values cited below.
  Since 2026-08-31 it also carries a `convergence_trends` block giving, for
  each of those series, the least-squares slope over the 2-3 h window and over
  the record's last hour, normalized by the window mean magnitude. The figures
  and metrics were regenerated with the renamed script on that date; every
  value quoted in this report reproduced unchanged.

## 4. When the wave was forced

The mountain forcing is a steady lower boundary condition. `mountain_forcing`
(`src/forcing.jl:115`) has no time argument; `apply_bottom_forcing!` imposes
w = U_bar dh/dx in the bottom ghost cells at every step. The forcing is
therefore switched on impulsively at t = 0 and remains constant for the full
10800 s. Nothing in the configuration ramps it up or turns it off. Among the
available forcing types only `:linear_ramp` and the Gaussian-envelope options
carry a temporal envelope; `:mountain` does not.

This impulsive switch-on matters for interpreting every time series in the
run: the model starts from a hydrostatic background with a uniform 57 m/s
wind and a suddenly present obstacle, so the early evolution is the
adjustment ("startup") transient of that initial-value problem, not part of
the stationary mountain-wave solution.

## 5. Why the station oscillations end near 4000 s

The oscillations in `station_timeseries_w.png` are the startup transient, not
the forced wave disappearing. Quantitatively, at station S3 (x_rel = 298 km,
z = 49.75 km):

- The transient peaks at |w| ~ 10.9 m/s near t ~ 350 s and rings down with a
  dominant period of 388 s, or 2.574 mHz, from an FFT of the first 4000 s. The
  FFT bin width over that window is 0.257 mHz, so the estimate spans roughly
  2.32 to 2.83 mHz.
- The envelope falls below 10% of its peak by t ~ 3190 s and is visually flat
  by ~4000 s, which is the "forcing ends at 4000 s" impression. The full-record
  wavelet (`deepdive_wavelet_w_full_z50.png`) shows the same thing in one
  view, a single power burst confined to the first ~25 minutes.

### What the 388 s period is not

An earlier version of this section read the 388 s as the buoyancy period,
comparing it against 2 pi / N_ref = 349 s at N_ref = 0.018 1/s. That comparison
is between a measurement made at z = 49.75 km and a frequency that belongs to
the lowest kilometre or two of the column. Taking N and the acoustic cutoff
omega_a = gamma g / (2 c) from this run's own background at the station height
gives, from `background_frequencies()` in
`scripts/analyze_arsia_stage1_deepdive.py`,

| quantity at z = 49.75 km | value |
|---|---|
| local N | 0.00962 1/s |
| local N / 2 pi | 1.532 mHz (2 pi / N = 653 s) |
| local acoustic cutoff / 2 pi | 1.963 mHz |
| measured ring-down | 2.574 mHz (388 s) |
| omega / N at the station | 1.68 |
| omega / omega_a at the station | 1.31, and 1.18 at the low edge of the bin |

A vertically propagating internal gravity wave requires omega < N. That holds
for 2.574 mHz only below about 1.75 km in this column, inside the near-surface
inversion where N reaches 0.0177 1/s. At the station the measured frequency sits
above both local references, and it stays above the local acoustic cutoff even
at the low edge of the FFT bin, so at that height it lies on the freely
propagating acoustic branch rather than the gravity branch. The agreement with
349 s is therefore better read as a coincidence at the 11 per cent level than as
an identification.

Two readings remain open and this record does not separate them.

1. An acoustic adjustment transient. The `:mountain` boundary condition switches
   on at t = 0 in a compressible model, which excites the acoustic branch
   directly. This is consistent with the frequency at the station, though the
   specific mode is not identified: the local cutoff varies only over 1.69 to
   2.00 mHz through the column, and a ground-to-sponge-base vertical resonance
   at c ~ 210 m/s would give about 571 s rather than 388 s.
2. Advection of small-scale structure past a fixed station. At U = 57 m/s a
   stationary structure of horizontal wavelength lambda_x appears at f = U /
   lambda_x, so 2.574 mHz corresponds to lambda_x = 22.1 km, or 20.1 to 24.6 km
   across the bin. This is the same Doppler mechanism invoked in Section 6 for
   the late 2-8 mHz burst. At dx = 4 km that scale spans 5.5 points, so it would
   be near the grid limit, and the decaying envelope is more suggestive of a
   transient than of steady advection.

What is unchanged by this is the timing conclusion, which is what Section 5 is
for. The oscillation decays, settles by ~3190 s and does not return, while the
forcing runs for the whole simulation. The revision affects the mechanism
attached to the frequency, not the statement that the early oscillations are the
startup transient. The wavelet reference lines have been redrawn at the local N
and local acoustic cutoff for the station height in question, with the surface
value kept as a third, clearly labelled line.
- What remains after the transient is the stationary mountain wave. A
  stationary wave has zero phase speed, so at a fixed station it produces an
  approximately constant w offset rather than oscillations (S3 late-time mean
  ~ +0.14 m/s at this z; larger offsets at other stations). A flat,
  non-oscillating trace is the expected signature of the forced wave, and it
  persists to t_end.

In short, the wave is forced during the whole run; only the transient ends
near 4000 s.

## 6. The wavelet feature near minute 165

The production `wavelet_w.png` and the full-record wavelets differ in two
ways that jointly explain why the late feature is prominent in one and barely
visible in the other:

1. `analyze_mars_mountain_sweep.py` trims the first 120 min
   (`--wavelet-trim-min 120`) precisely to remove the startup transient. With
   the transient excluded, the color scale (normalized to the maximum power
   in the window) stretches to fit the weak late-time signal, making it
   visible.
2. In the full-record wavelet the same late signal is present but appears
   faint because the normalization is dominated by the startup transient,
   which is roughly an order of magnitude stronger (the sweep metrics report
   transient_ratio ~ 8.1 for this run).

The late feature itself is a burst of apparent frequency 2-8 mHz at S3
starting near minute 155-160 and strengthening to the end of the record. Two
observations constrain its interpretation:

- Part of the burst lies above the acoustic cutoff (1.8 mHz) and above
  N_ref/2 pi (2.9 mHz), where vertically propagating gravity waves do not
  exist. However, the observed frequency at a fixed station is
  Doppler-shifted by the mean wind: a quasi-stationary structure of
  horizontal wavelength lambda_x advected at U ~ 57 m/s produces an apparent
  frequency f ~ U / lambda_x, so 2-8 mHz corresponds to lambda_x ~ 7-28 km.
  This matches the ~10-30 km spacing of the vertical striations in the final
  snapshot and in the composite panel. The burst is therefore consistent with
  small-scale structures sweeping past the station rather than with acoustic
  waves.
- The burst timing coincides with the growth of a mean-flow deceleration:
  <u>_x - 57 reaches -22.1 m/s at z = 57.75 km by t_end, growing mainly after
  t ~ 2 h (`mean_flow_evolution.png`, `deepdive_mf_ubar_profiles.png`), while
  the x-mean momentum flux becomes increasingly negative through the column.
  The tracker's Stage 1 results section already flags this run as not settled
  (late CoV 0.37, amplitude still rising at t_end, peak |w| ~ U).

### How much of the burst survives the cone of influence

The full-record wavelets were redrawn on 2026-09-02 with two logarithmic
panels and with the cone of influence marked, after the linear normalisation
was found to render everything outside the startup transient at a few per cent
of full scale. In the linear control at the same station the strongest feature
after 6 h sits 13.7 dB below the record peak, which on a 0 to 1 linear map is
0.04 and reads as empty.

That redraw qualifies the burst rather than removing it. Periods longer than
the cone of influence at a given time are contaminated by the zero padding at
the record ends, and the burst begins about 25 min before `t_end`. At
z = 49.75 km a patch near minute 160 at 3 to 5 mHz lies outside the cone and
survives a blanking threshold of 40 dB below the record peak; at
z = 57.25 km almost the whole feature lies inside it. The burst is therefore
still present in the data at the lower station, but its apparent growth to the
end of the record is not separable from the edge of the transform, and the
strength of the feature at the upper station should not be read from these
panels. A longer record is the only way to settle that, since the cone of
influence is set by the record length and not by the analysis.

Taken together, the data are consistent with progressive nonlinear steepening
of the stationary wave, generation of small horizontal scales, and momentum
deposition concentrated near 52-60 km. This pattern may indicate the early
stage of wave breaking, although the run terminates before any timestep
collapse and the interpretation remains tentative at this resolution
(dx = 4 km barely resolves 7-28 km scales, 2-7 points per wavelength).

## 7. Response at the upper sponge-free altitudes

For the rarefied-atmosphere coupling interest, the highest clean band is
53-60 km (between the AMEC band top and the sponge base). At t = 3 h in the
hill window (|x_rel| <= 500 km):

- Hill-window RMS w in the 53-60 km band is ~7-8 m/s, and max |w| reaches
  ~25-50 m/s (column peak 57.3 m/s at z = 29.75 km). At t = 2 h the same band
  showed RMS ~3 m/s, so the amplitude was still increasing at the end of the
  run.
- Station w at z = 57.25 km tracks the z = 49.75 km series closely during the
  transient (peaks ~ +17/-21 m/s, larger than at 49.75 km, consistent with
  rho^-1/2 amplitude growth) and shows the same post-transient flattening.
- The x-mean momentum flux at t = 3 h is negative through the column and
  still evolving; the mean-flow change peaks at -22.1 m/s at z = 57.75 km,
  immediately below the sponge base.

Two cautions for using these numbers as a coupling-level characterization.
First, the run is not statistically steady; every band-integrated quantity is
still trending at t_end. Second, the deceleration maximum sits against the
sponge base, so some of the structure near 57-60 km may be influenced by the
sponge boundary rather than by free dynamics. Both cautions support raising
the domain top (with realistic thermospheric viscosity) before extracting
DSMC forcing data, as anticipated in the tracker rows of 2026-06-30 on the
155 km attempt and the retained `mars_arsia_mons_my31_cap156km.csv` profile.
The present 60 km ceiling is a Stage 1 constraint, not a target.

## 8. What Stage 2 is

Stage 2 is the follow-on configuration defined in the tracker ("Arsia Mons
Stage 2 plan" section and the 2026-06-30 Decisions Log row). It keeps the
Stage 1 orographic forcing (same U_bar and h_max grids) but replaces the
uniform 57 m/s wind with the real reversing MCD zonal wind u_zonal(z)
(installed as a `u_init` closure) and switches the horizontal boundary to
periodic, which removes the Stage 1 inflow-edge artifact. The science target
is whether the stationary wave deposits momentum at the critical levels of
the real wind (zero crossings at ~4.5, 53.3, 77.7 km; the 53.3 km level
coincides with the top of the AMEC band).

Status at the time of writing: the Stage 2 3x3 sweep has been executed
(`output/arsia_mons_stage2/`, LX = 4000 km, all 9 runs converged at
t = 10800 s, wall-clock 2814-2971 s per run), and a smoke run exists
(`output/arsia_mons_stage2_smoke/`). Update, later on 2026-07-13: the Stage 2
analyzer (`scripts/analyze_arsia_stage2.py`, which adds the MF(z)
critical-level diagnostic) has now been run; results are in
`output/arsia_mons_stage2_analysis/` and summarized in the tracker's "Arsia
Mons Stage 2 results" section. For this corner (U=57, h=3 km) the Stage 2
field is confined below the 53.3 km critical level with hill peak |w| =
16.4 m/s (vs 49.6 m/s here) and about 85% of the jet-level momentum flux
deposited across 45-56 km, so the runaway growth documented in this report
does not recur under the real reversing wind.

## 9. Next steps

Ordered by proximity to the current work; tracker references in parentheses.

1. Run the Stage 2 analysis, `python3 scripts/analyze_arsia_stage2.py`, over
   the existing `output/arsia_mons_stage2/` runs. Priority on the U=57/h=3
   corner for a same-corner comparison with this report; the MF(z) overlay
   against the 4.5/53.3/77.7 km critical levels is the primary product
   (TODO_tracker, "Arsia Mons Stage 2 plan").
2. Done, on a different corner (2026-08-31). The 6 h extension was run on
   U=40/h=3 rather than U=57/h=3, at the user's request, because the U=40
   column at h=1 and h=3 km produces the cleaner wave field. The momentum
   budget of U=40/h=3 was first confirmed to be trending at t_end in both
   Stage 1 and Stage 2, so the same motivation applies to it. Outcome: the
   run does not approach a quasi-steady state. It continues toward
   overturning, with a decelerated layer that descends from ~55 km into and
   below the AMEC band between 3 h and 6 h. Details in the tracker section
   "Arsia Mons Stage 1 U=40/h=3 time extension". The U=57/h=3 corner itself
   has not been extended; the U=40/h=3 result suggests it would follow the
   same route sooner, since it is already further into overturning at 3 h.
3. Optionally re-run this corner with the snapshot cadence shortened (e.g.
   every 25 steps, ~16 s) so apparent frequencies up to ~30 mHz and the
   7-28 km striations are better resolved in station records; the present
   8 mHz Nyquist truncates the late-time burst (Section 3).
4. Prepare the thermospheric extension needed before DSMC coupling: raise the
   domain top above 90 km using `mars_arsia_mons_my31_cap156km.csv`, with a
   raised viscosity cap so molecular viscosity damps the wave physically
   above ~100 km (tracker rows of 2026-06-30 on the 155 km blowup; user goal
   of coupling at rarefied altitudes, for which the 60 km sponge base is
   low).
5. DSMC coupling file writer once a trustworthy coupling level exists
   (tracker Outstanding TODO, "DSMC coupling file writer").
6. Unrelated to this case but still open: the residual sponge reflection
   discrepancy (r = 3.74% vs the 0.55% target) and the Eckermann et al.
   (2011) CO2 radiative damping coefficients (tracker Outstanding TODOs).

## 10. Reproducibility

- Simulation: `julia -t auto --project=. scripts/sweep_arsia_mons_stage1.jl`
  (this run: 17510 steps, wall-clock 1622 s as one of 9 concurrent threads).
- Sweep analysis: `python3 scripts/analyze_mars_mountain_sweep.py --sweep-dir
  output/arsia_mons_stage1 --outdir output/arsia_mons_stage1_analysis
  --run-pattern 'run_U*_h*' --station-z-km 50` (see script header for the
  composite flags).
- Deep dive: `python3 scripts/analyze_arsia_stage1_deepdive.py` (defaults
  point at this run; 16.5 s wall-clock). The script takes `--run-dir`,
  `--compare-run`, `--u-bar`, `--label` and `--profile-hours`, so the same
  diagnostics apply to any Stage 1 run.
