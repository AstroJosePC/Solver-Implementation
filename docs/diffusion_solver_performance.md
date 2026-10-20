# Diffusion solver performance

Campaign doc for the cost of the implicit Crank-Nicolson diffusion step in `src/diffusion.jl`
and for the solver switches that reduce it. Sections are numbered within this doc and are
cited as `diffusion_solver_performance.md §<n>`. The physics, the time loop and the
verification list are those of `reflection_convergence_tests.md`; nothing registered there
is changed here.

## 1. Spectral-tridiagonal diffusion solver, on/off pair on the Stage 2 face path (pre-registration)

Written on 2026-10-19 and 2026-10-20, before any run of the set below was launched
(the calibration pilot is descriptive and is not part of the set). Committed on its own
and tagged `prereg/s01-diffusion-solver`; the tag timestamp precedes the launch.

Policy, unchanged: `SimConfig.diffusion_solver` defaults to `:lu` and
`SimConfig.diffusion_dt_tol` to 0.01, which reproduce the earlier code bit for bit. Every
run here passes `bcface=1` and `pbar=1` explicitly. No production flux is quoted from a
non-default solver before the readings below are scored and the switch is merged to `main`.

The change under test is one degree of freedom, the linear solver of the Crank-Nicolson
diffusion step (Srivastava et al. 2022, Eqs. 12 to 16), at fixed equations, fixed `dt` rule
and fixed refactorization rule.

### What is settled before this registration

**a. The operator is separable and the spectral solve is exact.** The left-hand side of
the Crank-Nicolson step is `I - (dt_ref/2) D`, with `D = kappa_z[j] delta_xx / dx^2 +
d/dz(kappa d/dz)` in divergence form (face `kappa` the mean of the two neighbours, Neumann
rows at `j = 1` and `j = nz`), on a Cartesian grid with uniform `dx` and `dz` and no metric
terms. The x part is `kappa_z[j]` times the periodic second difference, which is circulant
at each level with a coefficient that depends on `j` only. A real FFT in x therefore
diagonalizes it exactly: mode `m` (0 to `nx/2`) has the eigenvalue `-lam_m / dx^2` with
`lam_m = 2 - 2 cos(2 pi m / nx)`, and each mode leaves a real tridiagonal system in z with
diagonal `1 + fx_j lam_m + fz_top + fz_bot` (Neumann rows as in `build_cn_matrix`) and
off-diagonals `-fz_bot`, `-fz_top`. `:fft_tridiag` solves these `nx/2 + 1` systems with the
Thomas algorithm. It is the same discrete system as the global matrix, so the results agree
to round-off and are not bit-identical. The global matrix is periodic in x for every
horizontal boundary type, so the equivalence holds for the Stage 1 `:outflow` runs too.
Validity of the spectral path: uniform `dx` and `dz`, diffusivity a function of `z` only,
no metric terms. An x-dependent diffusivity (for example a Smagorinsky closure) breaks the
separability and would need the global solver.

Checks already run on commit `3ef4ac0` (descriptive, not registered):
- Per-mode algebra against a dense solve of `build_cn_matrix` on `nx` = 1 to 33 (odd and
  even, with and without the Nyquist mode), `nz` = 1 to 9, flat and 1e-4 to 1e6
  diffusivities, `dt` 0.7 and 1e4 s: worst relative solution difference 5.7e-13.
- `scripts/check_diffusion_solvers.jl` (rerun on 2026-10-19 at 22:11, 44 s, julia -t 1,
  alone), three steps `dt0`, `1.004 dt0` (lagged), `0.8 dt0` (refactor) on the Stage 2
  grids:

| grid | case | `:lu_reuse` max\|diff\| | `:fft_tridiag` rel_field | `:fft_tridiag` rel_incr |
|---|---|---|---|---|
| 1000 x 190 (dz 500 m) | random and real snapshot (11.78 h, U40 h1, face) | 0.0 | 3.4e-14 | 5.6e-9 |
| 1000 x 380 (dz 250 m) | random and real snapshot (S15F2, 11.73 h) | 0.0 | 4.0e-14 | 2.9e-11 |

  `rel_field = max|f - f_lu| / max|f_lu|` and `rel_incr` divides the same difference by
  the diffusion increment, over `u'`, `w` and `T'`. The 5.6e-9 is a round-off difference of
  about 1e-13 absolute against a small per-step increment. Log and JSON:
  `~/MAGNUSP_runs/scratch_dfft/check_diffusion_solvers_rerun.{log,json}`.
- Whole run, 0.1 h of the S14M1 command with `diffsolver=fft`: max\|diff\| 3.5e-12
  (absolute, over `rho`, `u`, `w`, `P`, `T`) against the default path.
- Diffusion-only replays of the recorded CFL `dt` sequences (S15F1 at 10.0 h, S15F2 at
  10.6 h and 11.7 h; 2400 to 4000 steps): identical refactorization steps for the three
  solvers, `:lu_reuse` 0.0 exactly, `:fft_tridiag` rel_field at most 5.6e-13 with no growth
  over the replay.

**b. The `dt` lag is kept.** The factorization is built with `dt_ref`, the right-hand side
(`build_cn_rhs!`) with the current `dt`, and the factorization is rebuilt only when
`|dt - dt_ref| / dt > diffusion_dt_tol` (default 0.01). All three solvers apply this rule
in the same way, so they solve the same equations step for step, and `:fft_tridiag`
counts its rebuilds of the Thomas coefficients as factorizations. Under this rule the
effective diffusivity of a step is `kappa (dt + dt_ref) / (2 dt)`, a bias of at most 0.5 %
in the diffusion term. Removing the lag (`difftol=0`) is a different discrete scheme and is
not part of this pair.

**c. Measured cost.** Standalone medians (julia -t 1, BLAS 1, one process), and the
modelled cost per diffusion step from the replayed refactorization rates:

| solver | grid | factorization (nu + alpha) [s] | 3 CN solves [s] | full step, no refactor / refactor [s] | refactor rate per 1000 steps | modelled s per diffusion step |
|---|---|---|---|---|---|---|
| `:lu` | dz 500 | 0.701 | 0.0223 | 0.070 / 0.786 | 21.5 (S15F1, 10.0 to 10.7 h) | 0.085 |
| `:lu_reuse` | dz 500 | 0.435 | 0.0223 | 0.071 / 0.521 | 21.5 | 0.081 |
| `:fft_tridiag` | dz 500 | 0.00014 | 0.0024 | 0.049 / 0.049 | 21.5 | 0.049 |
| `:lu` | dz 250 | 2.068 | 0.0469 | 0.145 / 2.262 | 41.2 (10.6 h); 580 (11.7 h) | 0.232; 1.373 |
| `:lu_reuse` | dz 250 | 1.296 | 0.0469 | 0.146 / 1.483 | 41.2; 580 | 0.200; 0.921 |
| `:fft_tridiag` | dz 250 | 0.00027 | 0.0049 | 0.102 / 0.101 | 41.2; 580 | 0.102 |

The non-diffusion part of a time step costs 0.013 s (dz 500) and 0.026 s (dz 250). About
93 % of the `:fft_tridiag` step is the extraction of primitives, six `zeros(nx, nz)`
allocations and the reconstruction, which are shared by all solvers and are not changed
here (a separate, behaviour-preserving optimisation). In the S15F2 window at 11.7 h, where
54 % of steps change `dt` by more than 1 %, `:fft_tridiag` is modelled at 1.3 h of wall
clock per simulated hour against 14.4 h for `:lu` and 9.8 h for `:lu_reuse`.

**d. Default path.** On commit `3ef4ac0` (driver, Queney driver and `src/` hash
`8683084fb5c6932a`), the launcher `check` (see Check before launch) reproduced the
reference records with max\|diff\| exactly 0.0. Snapshots now carry the global attributes
`diffusion_solver` and `diffusion_dt_tol`, so default-path reproduction is judged on the
variables (`rho`, `u`, `w`, `P`, `T`, `time_s`), not on file checksums.

### Runs

Common arguments `CM` = `pbar=1 ramp=2 k0=0.005 lz=95 bcface=1`
(`sweep_arsia_mons_stage2.jl`, dz 500 m, CFL 0.4, outint 100). The solver argument is
`diffsolver=lu` (default, no tag), `diffsolver=fft` (tag `_dfft`) or
`diffsolver=lu_reuse` (tag `_dlureuse`). All six runs go in one launcher set `d1`, one
process each, `julia -t 1`, `OPENBLAS_NUM_THREADS=1`, concurrency 6, under
`caffeinate -s -i` with the cost monitor. The set writes under
`MAGNUSP_OUTPUT=~/MAGNUSP_runs/d1` and archives to `output/d1/`, so that the `:lu` member
of segment A overwrites neither the local nor the archived copy of S14M1, which has the
same directory name. The set is launched from the worktree of `dev/diffusion-fft-tridiag`,
whose git-ignored `output` is a symbolic link to the archive `output/` of the main
checkout, so `check` reads the archived references and `archive` writes to the mount.

| id | segment | solver | arguments | `dt` [s] | steps | s per step (source) | est. wall | disk [GB] |
|---|---|---|---|---|---|---|---|---|
| S01L1 | A, `h` = 0.1 km, 9 h, `dt` nearly fixed | `:lu` | `hours=9 CM u=40 h=0.1` | 0.659 | 49 150 | 0.103 (PR1-U40-h0.1) | 1.4 h | 3.8 |
| S01F1 | A | `:fft_tridiag` | `hours=9 CM u=40 h=0.1 diffsolver=fft` | 0.659 | 49 150 | about 0.082 (0.103 less 0.021, item c) | 1.1 h | 3.8 |
| S01R1 | A | `:lu_reuse` | `hours=9 CM u=40 h=0.1 diffsolver=lu_reuse` | 0.659 | 49 150 | 0.103 | 1.4 h | 3.8 |
| S01L3 | B, `h` = 3 km, U57, 5 h, shrinking `dt` | `:lu` | `hours=5 CM u=57 h=3` | about 0.57 | about 31 500 | 0.139 (PR1-U57-h3.0) | 1.2 h | 2.4 |
| S01F3 | B | `:fft_tridiag` | `hours=5 CM u=57 h=3 diffsolver=fft` | about 0.57 | about 31 500 | about 0.082 | 0.7 h | 2.4 |
| S01R3 | B | `:lu_reuse` | `hours=5 CM u=57 h=3 diffsolver=lu_reuse` | about 0.57 | about 31 500 | about 0.125 | 1.1 h | 2.4 |

Directories under `~/MAGNUSP_runs/d1`:
`arsia_mons_stage2_9h_K0_0.005_LZ95_pbarP_ramp2h{,_dfft,_dlureuse}_bcface/run_U40.0_h0.1`
and
`arsia_mons_stage2_5h_K0_0.005_LZ95_pbarP_ramp2h{,_dfft,_dlureuse}_bcface/run_U57.0_h3.0`.

Basis of the estimates: the s per step are the awake rates of the s15 batch at concurrency
12. PR1-U57-h3.0 ran at 0.139 s per step against 0.103 for PR1-U40-h0.1; the difference,
0.036 s, is taken as the refactorization surcharge of `:lu` on segment B (about 5 % of
steps at 0.70 s each), scaled by 0.62 for `:lu_reuse` and removed for `:fft_tridiag`.
Segment B ends at 5 h, before the `dt` collapse of PR1-U57-h3.0 at 5.47 h. The wall clock of
the set is set by S01L1 and S01R1, about 1.4 h, and the CPU time summed over the six runs is
about 6.9 h. Disk is about 19 GB.

### Hypotheses

- H0, same discrete system: `:fft_tridiag` and `:lu` differ only by floating-point
  round-off in the diffusion solve. Prediction: field differences stay at round-off level
  (`rel` of `rho`, `u`, `P` and `T` of order 1e-13, absolute `w` difference of order 1e-10
  m/s, as in the calibration pilot below) and grow at most linearly with time wherever the
  flow is not sensitive to its initial state, the refactorization steps and the step counts
  are identical, and window-mean fluxes agree to 1e-7 or better.
- H1, different discrete system: the spectral path solves a different system (for example
  a wrong Neumann row, a wrong Nyquist eigenvalue or a missed `dt_ref`). Prediction: a
  difference of the order of the operator error times the accumulated diffusion increment,
  which appears from the first step and does not depend on the flow regime. The unit check
  bounds such an error at `rel_incr` 5.6e-9 per step; an operator error of 1e-6 or more would
  give field differences of 1e-8 or more within 0.1 h on segment A.
- H2, round-off amplification by the flow: in the breaking flow of segment B a round-off
  difference grows with the flow's own sensitivity. Prediction: the difference stays at
  round-off up to the onset of strong nonlinearity and then grows roughly exponentially,
  in the same way that any round-off perturbation of the `:lu` run would grow. H2 is
  compatible with H0 and is not evidence against the solver.

### Calibration pilot (descriptive)

The first draft of D-1 and D-2 set `rel_max` at most 1e-9, from a 0.1 h difference grown
linearly. A 0.1 h `w` difference that grew faster than linearly over the first 500 steps
put that band in doubt, so the first hour of both segments was run on 2026-10-19 to
calibrate the bands before registration. The pilot is descriptive and is not part of the
set; no reading below is scored on it.

- Runs: `hours=1 CM u=40 h=0.1` and `hours=1 CM u=57 h=3`, each with `diffsolver=lu` and
  `diffsolver=fft`, from the worktree at `b33db15` plus the then uncommitted text of this doc,
  4 processes, concurrency 4, `julia -t 1`, `OPENBLAS_NUM_THREADS=1`, `caffeinate -s -i`,
  22:38:51 to 22:47:31, under `MAGNUSP_OUTPUT=~/MAGNUSP_runs/scratch_dfft/pilot`.
  Step counts 5462 (A) and 5467 (B), equal between solvers; `:fft_tridiag` printed
  2 factorizations on each segment.
- Command: `bitcheck --relative`, 55 shared snapshots each. Products and logs:
  `output/analysis/s01/pilot/` (`rel_A.json`, `rel_B.json`, `logs/`).

| quantity over 0 to 1 h | segment A (U40, `h` 0.1 km) | segment B (U57, `h` 3 km) |
|---|---|---|
| max `rel` of `rho`, `u`, `P`, `T` | 1.2e-13, 2.3e-13, 1.2e-14, 2.1e-13 | 1.1e-13, 2.0e-13, 1.5e-14, 2.0e-13 |
| max `rel_w` (time) | 9.96e-9 (1252 s) | 2.98e-10 (923 s) |
| `rel_w` at 1 h | 4.6e-9 | 1.0e-10 |
| absolute `w` difference, 0.5 to 1 h | 4.2e-11 to 5.8e-11 m/s, no trend | 4.3e-11 to 5.7e-11 m/s, no trend |
| max\|`w_ref`\| at the `rel_w` maximum | 5.9e-3 m/s | 1.0e-1 m/s |
| power-law slope of `rel`, 0.5 to 1 h | `rho` 0.23, `u` 0.43, `T` 0.17, `w` −0.57 | `rho` 0.04, `u` 0.18, `T` 0.12, `w` −0.89 |
| solver wall clock, `:lu` / `:fft_tridiag` | 508.7 s / 403.2 s (ratio 0.79) | 511.8 s / 399.0 s (ratio 0.78) |

Reading. The absolute `w` difference reaches a plateau of about 5e-11 m/s within the first
0.3 h on both segments. `rel_w` is the largest of the five only because it is normalized by
max\|`w`\|, which during the ramp on segment A is 6e-3 to 4e-2 m/s, against about 1e2 to 1e3
for `rho`, `u`, `P` and `T` in their units. The growth that made the first band doubtful was
the approach to this plateau. Over 0.5 to 1 h no variable grows faster than about `t`^0.43.
The pilot covers the ramp only. It does not show the growth after 1 h on segment A, or the
growth on segment B as the flow approaches breaking, which H2 expects.

How the bands follow from the pilot. Each band is the 1 h maximum, extrapolated linearly in
time to the end of the scored interval, times a stated margin. Linear growth is an upper
envelope for the measured slopes of at most 0.43.
- Segment A, 9 h. `rho`, `u`, `P`, `T`: 2.3e-13 × 9 = 2.1e-12, with a margin of 5, gives
  1e-11. `w`: 9.96e-9 × 9 = 9.0e-8, with a margin of about 1.1, gives 1e-7. The strict
  choice (an earlier version of this text had 1e-6, margin about 11) rests on the plateau of the absolute `w`
  difference after 0.3 h and on the 1 h value of `rel_w`, 4.6e-9; it is met only if the
  growth after the ramp stays well below linear, since with the factor 2 rule below a
  value above 5e-8 is read in row 2. `MF(30.25)` is quadratic in the perturbation, and its
  relative change is at most about twice `rel_w` in the free column; its band is 1e-7 as
  well.
- Segment B, `t` at most 2 h. `rho`, `u`, `P`, `T`: 2.0e-13 × 2 = 4e-13, with a margin of 25,
  gives 1e-11. `w`: 2.98e-10 × 2 = 6e-10, with a margin of about 17, gives 1e-8.
- The H1 boundary sits three orders above the H0 band for each variable. The H1 reading is
  taken mainly on `rho`, `u`, `P` and `T`, whose normalization does not depend on the wave
  amplitude.

### Readings fixed before launch

Each reading is recorded against these statements, and a failed prediction is recorded as
failed. `rel` of a variable is `max|f_new - f_ref| / max|f_ref|` at a shared snapshot,
over all cells, for `rho`, `u`, `w`, `P` and `T` separately, with `rel_max` the largest of
the five. A value within a factor 2 of a band edge is reported at the edge, and the reading
of the higher band applies.

**D-0, default path over the full segments.**
- S01L1 against the archived S14M1
  (`output/arsia_mons_stage2_9h_K0_0.005_LZ95_pbarP_ramp2h_bcface/run_U40.0_h0.1`), every
  shared snapshot to 9 h: max\|diff\| exactly 0.0.
- S01L3 against the archived PR1-U57-h3.0
  (`output/arsia_mons_stage2_24h_K0_0.005_LZ95_pbarP_ramp2h_bcface/run_U57.0_h3.0`), every
  shared snapshot to 5 h: max\|diff\| exactly 0.0.

Met: the default path reproduces the earlier code over the full segments, which extends
the 0.1 h `check`. Failed: the switch is not merged until the difference is located; a
failure only in the second item, starting from step 0, would point to a dependence on the
`hours` argument rather than on the solver. Expected: met in both.

**D-1, segment A, `:fft_tridiag` against `:lu` (S01F1 against S01L1).**

| outcome | reading |
|---|---|
| at every shared snapshot, `rel` of each of `rho`, `u`, `P`, `T` at most 1e-11 and `rel_w` at most 1e-7, and the relative \|change\| of the 6 to 9 h `MF(30.25)` at most 1e-7 (1e-5 %) | consistent with H0: the two solvers give the same solution to round-off over 9 h; `:fft_tridiag` can replace `:lu` on this column |
| any of `rho`, `u`, `P`, `T` in (1e-11, 1e-8], or `rel_w` in (1e-7, 1e-4], or the `MF(30.25)` change in (1e-7, 1e-4], and none above | a departure above the calibrated round-off and below any quoted error bar; reported, and the switch is not adopted until the departure is located |
| any of `rho`, `u`, `P`, `T` above 1e-8, or `rel_w` above 1e-4, or the `MF(30.25)` change above 1e-4 | consistent with H1: the solvers differ; not adopted |

Second statement: the factorization counts printed by S01F1 and S01R1 are equal, and the
step counts of S01F1 and S01L1 are equal. Reported without a band: the changes of
`MF(1.25)`, `MF(10.25)`, `MF(50.25)`, `R10/30` and `R50/30` over 6 to 9 h, and `rel` per
variable against time, and the first-hour difference of S01L1 and S01F1 from the pilot runs
of the same solver (the pilot used `hours=1`). Expected: row 1 and the second statement met,
with low to moderate confidence for the `rel_w` item. The pilot shows a plateau of the
absolute `w` difference up to 1 h, but it does not cover 1 to 9 h, and the band leaves a
margin of only about 1.1 over linear growth. Basis of the bands: the calibration pilot above.

**D-2, segment B, `:fft_tridiag` against `:lu` (S01F3 against S01L3).** No `h` = 3 km flux is
quoted (Section 14 classification of `reflection_convergence_tests.md`, unchanged).

| outcome | reading |
|---|---|
| at every shared snapshot with `t` at most 2 h (the ramp), `rel` of each of `rho`, `u`, `P`, `T` at most 1e-11 and `rel_w` at most 1e-8, and the `zh3` crossing times of the two runs agree within one scored snapshot (stride 2, about 115 s) | consistent with H0, with H2 for any growth after 2 h; `:fft_tridiag` gives the same breaking-flow diagnostics as `:lu` |
| any of `rho`, `u`, `P`, `T` above 1e-11 or `rel_w` above 1e-8 at some `t` at most 2 h, or a `zh3` crossing time differs by more than one scored snapshot | the solvers differ in this regime (H1), or the flow amplifies round-off before 2 h; reported with the growth curve, and the switch is not adopted for breaking runs until the cause is located |

The `zh3` crossings expected inside 5 h are the ground `|T/T_0 - 1| > 0.5` crossing (2.52 h
in PR1-U57-h3.0) and the top `|w| > 5` m/s crossing (4.36 h); the free-column crossing
(5.24 h) and the `dt` collapse (5.47 h) lie after the segment and are not scored.
Reported without a band (descriptive): `rel` per variable against time, the first times at
which `rel_max` exceeds 1e-9, 1e-6 and 1e-3, an e-folding time fitted to the growth, the
first step at which the `dt` sequences differ, the factorization counts of S01F3 and S01R3
per 1000 steps, and the hourly `MF(30.25)` of both runs. Expected: row 1, with low
confidence for the 4.36 h crossing, since a round-off difference may have grown to a
visible amplitude by then (H2).

**D-3, the interim `:lu_reuse` (S01R1 against S01L1, S01R3 against S01L3).** Criterion:
max\|diff\| exactly 0.0 at every shared snapshot of both segments. Met: `:lu_reuse` is
bit-identical to `:lu` on these segments and can serve as the interim switch for runs that
must keep the global LU. Failed: UMFPACK's reuse of the symbolic analysis changed the
pivoting for some `dt_ref`; `:lu_reuse` is then judged by the D-1 and D-2 bands. Expected:
met, with moderate confidence (bit-identical in the unit check and in all 1643 replayed
refactorizations, but the reuse keeps the pivot order of the first matrix).

**D-4, cost.** Reported from the run logs and `d1.cost.tsv`: solver wall clock, `real`
time and s per step of each run at concurrency 6, and the ratios `:fft_tridiag`/`:lu` and
`:lu_reuse`/`:lu` per segment. Expected: about 0.80 on segment A and about 0.60 on segment B
for `:fft_tridiag`, and about 1.0 and 0.9 for `:lu_reuse`. No band: the cost does not decide
adoption, which rests on D-1 to D-3.

### Scoring commands

All read archived output only (`analyze_dt_convergence.py` unless stated; the run
directories of this set are under `output/d1/`), with the products written to
`output/analysis/s01/`:
- D-0: `bitcheck --new <S01L1 dir> --ref output/arsia_mons_stage2_9h_K0_0.005_LZ95_pbarP_ramp2h_bcface/run_U40.0_h0.1`
  and `bitcheck --new <S01L3 dir> --ref output/arsia_mons_stage2_24h_K0_0.005_LZ95_pbarP_ramp2h_bcface/run_U57.0_h3.0`.
- D-1: `bitcheck --relative --new <S01F1 dir> --ref <S01L1 dir> --out output/analysis/s01/rel_A_fft.json`;
  `suite --run S01L1=<dir> --run S01F1=<dir> --window 6 9 --pair S01L1:S01F1 --outdir output/analysis/s01/mf_A`
  (the change is read from `s14_suite.json` at full precision, not from the printed table);
  `grep -h "factorizations" <S01F1 log> <S01R1 log>` and the step counts from the run logs.
- D-2: `bitcheck --relative --new <S01F3 dir> --ref <S01L3 dir> --out output/analysis/s01/rel_B_fft.json`;
  `zh3 --run <S01L3 dir> --stride 2 --out output/analysis/s01/z3_B_lu.json` and the same for
  S01F3 (`z3_B_fft.json`); `suite --run S01L3=<dir> --run S01F3=<dir> --window 1 5 --pair S01L3:S01F3 --outdir output/analysis/s01/mf_B`
  for the hourly `MF` (descriptive).
- D-3: `bitcheck --new <S01R1 dir> --ref <S01L1 dir>` and `bitcheck --new <S01R3 dir> --ref <S01L3 dir>`.
- D-4: the `d1` rows of `docs/run_matrix.csv` and `~/MAGNUSP_runs/d1/logs/d1.cost.tsv`
  (archived to `output/d1/logs/`).

### Code for Section 1

Already on `dev/diffusion-fft-tridiag` (`18c96ef`, `5a9cf57`, `3ef4ac0`):
- `SimConfig.diffusion_solver` (`:lu` default, `:lu_reuse`, `:fft_tridiag`) and
  `SimConfig.diffusion_dt_tol` (0.01 default), checked at the start of `run_simulation`;
  `:lu_reuse` refactors with `lu!(F, A; reuse_symbolic = true)`; `:fft_tridiag` builds the
  per-mode Thomas coefficients and solves through FFTW `rfft`/`irfft` plans (ESTIMATE);
  non-default solves fail on NaN or Inf.
- Snapshot attributes `diffusion_solver` and `diffusion_dt_tol`; log lines for a
  non-default solver or tolerance, including the factorization count at the end of a run.
- Driver options `diffsolver=` and `difftol=` in `sweep_arsia_mons_stage2.jl` (tags
  `_dfft` or `_dlureuse`, then `_dtol<X>`, after the `_dz` tag) and `--diffsolver`,
  `--difftol` in `run_mountain_benchmark_strict.jl`.
- `scripts/check_diffusion_solvers.jl`, the unit equivalence on the Stage 2 grids.

Added afterwards, each as its own commit, none of them read by the solver except the
line-ending restore:
- `7c7d7d7`: `src/diffusion.jl` back to CRLF line endings (commit `18c96ef` had converted
  it to LF); content unchanged, so the change against `94da3cf` reads as its 281 lines.
- `fb532ed`: `.gitignore` whitelists `docs/diffusion_solver_performance.md`.
- `c349f7f`: `analyze_dt_convergence.py bitcheck --relative [--out <json>]`: per shared
  snapshot, `time_s` of both runs and, per variable, `max|diff|` and `rel`; exit status 0 in
  this mode. Without `--relative` the subcommand is unchanged.
- `b33db15`: `run_dt_convergence_S2T95.sh` set `d1` with the six cases above, run IDs
  `S01L1` to `S01R3`, stamps `d1.start`, `MAGNUSP_OUTPUT` defaulting to `~/MAGNUSP_runs/d1`
  and the archive to `output/d1` for this set only; every earlier set keeps its paths.
- `8abf18e`: merge of `main` (`479e8e1`, the Section 16 splitting switches of
  `reflection_convergence_tests.md`, `sweep_order` `:alternate` and `:strang`,
  `operator_split`, set `s16`) into this branch, so that the set runs on the code that
  lands on `main`. Conflicts were confined to the driver header and the launcher (union of
  the `s16` and `d1` arms); `src/diffusion.jl` and `src/io.jl` are those of `b33db15`, and
  `src/euler.jl` is that of `main`. The combined tag order is
  `..._dz<N>_dfft_dtol<X>_ub<Z>_bcface_spsplit_{zx|alt|strang}_opS_cfl<X>`. Every `d1` case
  uses the default splitting (`sweep_order` `:xz`, `operator_split` `:lie`).

### Check before launch

`dry d1` (2026-10-19, `~/MAGNUSP_runs/d1/logs/d1.dry.log`) resolved the six cases to
`pbar_mode` pressure, the face ground path, CFL 0.40 and the output directories listed
under Runs. `check d1` on commit `b33db15`, 2026-10-19 22:35:12 (67 s wall, nine
processes, under `caffeinate -s -i`), driver, Queney driver and `src/` hash
`5b1718cc227961b4`, all new switches at their defaults, from the worktree with
`MAGNUSP_OUTPUT=~/MAGNUSP_runs/d1` and the archived references of the main checkout:
max\|diff\| exactly 0.0 in every shared snapshot, 0.1 h of S14G1 (5 snapshots), C2R with
`probe=5` (5), S04 (6), C2R (5), S14M1 (6), and the Queney driver on its default path
against its pre-addendum-3 copy (2), 29 snapshots in all. The Queney face path and the
Queney face path with `pbar=1` completed. Summary
`~/MAGNUSP_runs/d1/logs/d1.check.out`, logs `~/MAGNUSP_runs/d1/logs/S14check_*.log`. This
replaces the earlier `check` on `3ef4ac0` (hash `8683084fb5c6932a`, also 29 snapshots at
0.0), whose hash predates the line-ending restore.

Repeated after the merge `8abf18e` (2026-10-19 23:59:40 to 2026-10-20 00:01:43, together
with the two items below, under `caffeinate -s -i`). `dry d1`
(`~/MAGNUSP_runs/d1/logs/d1.dry.merged.log`) resolved the same six directories. `check d1`
on `8abf18e`, driver, Queney driver and `src/` hash `6345878dbd21bf05`: max\|diff\| exactly
0.0 in every shared snapshot of the 29 above and of 0.1 h of S14V1a (5), 34 in all; the
Queney face paths completed; and the Section 16 Strang and alternating paths with Strang
operators reached `t_end`. Summary `~/MAGNUSP_runs/d1/logs/d1.check.merged.out`. This hash
is the one the launch requires. `scripts/check_diffusion_solvers.jl` on the merged code
reproduced its earlier maxima (`:lu_reuse` 0.0; `:fft_tridiag` `rel_field` 3.4e-14 and
4.0e-14, `rel_incr` 5.6e-9 and 2.9e-11 at dz 500 and 250 m). Descriptive, not a reading: a
0.1 h smoke pair of the face-path ridge with `opstrang=1` (the diffusion solve called twice
per step over dt/2), `:fft_tridiag` against `:lu`, gave equal step counts (487),
2 factorizations, `rel` of `rho`, `u`, `P`, `T` at most 1.5e-14 and `rel_w` at most
9.5e-10 (6 snapshots). Copies in `output/analysis/s01/merge_smoke/`. The calibration pilot
ran on `b33db15`, before the merge; the merged `check` reproduces the default paths it
used bit for bit.

### Variant L, a quantized dt ladder (registered separately, not launched with this pair)

A ladder `dt_q = dt0 (1 - q)^k`, with `k` the smallest integer such that `dt_q` is at most
the CFL `dt`, would keep the global LU but reuse a small cache of factorizations. Offline
estimates on the replayed `dt` sequences (descriptive, from `count_rule.py` and
`ladder_cost.py`):

| ladder | S15F2 at 10.6 h, factorizations / extra steps | S15F2 at 11.7 h, factorizations / extra steps |
|---|---|---|
| current 1 % lagged rule | 165 / 0 | 1392 / 0 |
| q 0.01, no cache | 185 / +0.5 % | 1647 / +0.5 % |
| q 0.05, cache 4 | 19 / +2.7 % | 294 / +2.5 % |
| q 0.05, hysteresis 1, cache 4 | 18 / +5.8 % | 167 / +4.8 % |

A 4-rung cache holds about 2.3 GB per run at dz 250 m (about 28 GB for twelve runs). Even
the best ladder leaves `:lu` about 2.4 times slower than `:fft_tridiag` with the unchanged
rule in the 11.7 h window, while `:fft_tridiag` refactors in 0.3 ms and gains nothing from
a ladder.

The ladder differs from the reference in a second degree of freedom: it changes the `dt`
sequence and removes the lag (left- and right-hand sides both with `dt_q`). It therefore
cannot reproduce the reference bit for bit, and it is judged as a time-step change, not as
a solver change. If it is ever launched, its registration carries:
- a switch `dt_ladder_q` with default 0 (off) that reproduces the reference, and the cache
  size as a separate option;
- the time-step criteria of the standing table: free-column `MF` change below 5 % per CFL
  halving (production) and below 2 % (mechanism tests), here applied to the ladder against
  the reference at CFL 0.4, with the stated expectation that a mean `dt` shortened by 2.5 to
  6 % changes `MF(30.25)` by about 0.03 % (scaled from the +0.35 % per halving of S14M1/M2),
  and a band of 0.5 % on the 6 to 9 h `MF(30.25)`;
- the step-count increase and the peak resident memory per run as reported quantities.

Registered decision rule: Variant L is not implemented unless D-1 fails, or unless a later
configuration needs the global LU (an x-dependent diffusivity).

### Naming

Sections of this doc are numbered within it and cited as `diffusion_solver_performance.md
§<n>`. Run IDs use the two-digit section form `S01<letter><n>`, since the project rule
forbids IDs that start with `S` plus a single digit; analysis products go to
`output/analysis/s01/` and the tags are `prereg/s01-diffusion-solver` and
`run/d1-<YYYYMMDD>`.

### Results

Scored on 2026-10-20 with the registered commands on the archived output, after the
registration (`93a716a`, tag `prereg/s01-diffusion-solver`, 00:14:32) and the launch tag
`run/d1-20261020`.

Run record. Set `d1`, launched 2026-10-20 00:14:42 from `93a716a` (clean tree), driver and
`src/` hash `6345878dbd21bf05`, launcher copy `~/MAGNUSP_runs/d1/logs/d1.launcher.sh`
(hash `23449f8fc0f9248c`), concurrency 6, `julia -t 1`, `OPENBLAS_NUM_THREADS=1`, cost
monitor holding `caffeinate -s -i`. Last case ended 01:39:40 (1.42 h); no sleep or DarkWake
entries in `pmset -g log` over the batch. Archived 01:40:00 to 01:56:29 to `output/d1/`, 0
checksum differences and equal file counts in all six sweep directories; logs, stamps,
cost history and launcher copy in `output/d1/logs/`. The local copies are kept; the `rm`
command is printed in `~/MAGNUSP_runs/d1/logs/d1.archive.log`. Products in
`output/analysis/s01/`.

| id | solver | steps | `t_end` [s] | mean `dt` [s] | `real` [s] | s per step | factorizations |
|---|---|---|---|---|---|---|---|
| S01L1 | `:lu` | 49 153 | 32 400 | 0.659166 | 4758.78 | 0.0968 | not printed (default path) |
| S01F1 | `:fft_tridiag` | 49 153 | 32 400 | 0.659166 | 3718.07 | 0.0756 | 2 |
| S01R1 | `:lu_reuse` | 49 153 | 32 400 | 0.659166 | 5093.67 | 0.1036 | 2 |
| S01L3 | `:lu` | 29 545 | 18 000 | 0.609240 | 3386.70 | 0.1146 | not printed (default path) |
| S01F3 | `:fft_tridiag` | 29 545 | 18 000 | 0.609240 | 2257.14 | 0.0763 | 400 |
| S01R3 | `:lu_reuse` | 29 545 | 18 000 | 0.609240 | 3256.94 | 0.1102 | 400 |

Note on the scoring run. The ten scoring commands were first started concurrently; seven
of them stopped at the opening of a snapshot on the mount with `NetCDF: HDF error`, at
files that two processes were opening at the same time. The same files opened without
error one at a time, and their checksums matched the local copies. The seven commands were
rerun unchanged, one after another (01:58:49 to 02:02:59), and every reading below comes
from a command that completed. No command or argument was altered.

**D-0, met.** S01L1 against the archived S14M1: 492 shared snapshots to 9 h, max\|diff\|
exactly 0.0. S01L3 against the archived PR1-U57-h3.0: 296 shared snapshots to 5 h,
max\|diff\| exactly 0.0. The default path of the merged code reproduces the earlier records
over the full segments.

**D-1, row 1 met; second statement met.** S01F1 against S01L1, 492 shared snapshots to
32 365 s.

| quantity | value | band of row 1 |
|---|---|---|
| max `rel_rho` | 1.39e-13 (1.54 h) | at most 1e-11 |
| max `rel_u` | 1.11e-12 (6.04 h) | at most 1e-11 |
| max `rel_P` | 1.72e-14 (3.53 h) | at most 1e-11 |
| max `rel_T` | 2.71e-13 (5.88 h) | at most 1e-11 |
| max `rel_w` | 9.96e-9 (0.35 h) | at most 1e-7 |
| relative change of the 6 to 9 h `MF(30.25)` | −2.55e-12 | at most 1e-7 |

Every value lies more than a factor 2 below its edge, so the factor 2 rule does not apply.
The factorization counts of S01F1 and S01R1 are equal (2 and 2), and the step counts of
S01F1 and S01L1 are equal (49 153). Reported without a band: the changes over 6 to 9 h
of `MF(1.25)` −2.7e-12, `MF(10.25)` 5.3e-12, `MF(50.25)` 4.7e-13, `R10/30` 7.8e-12 and
`R50/30` 3.0e-12 (relative). `rel_w` reaches its maximum during the ramp, at the same
snapshot and value as in the calibration pilot, then falls to 3.9e-9 at 1 h, 1.0e-9 at
4 h and 4.1e-10 at 9 h; the absolute `w` difference stays between 2.7e-11 and 5.3e-11 m/s
after 1 h, with no trend. The snapshot times of the two runs differ from step 300 on, by at
most 1.8e-11 s, a round-off difference of the CFL `dt`. The first hour of S01L1 and S01F1 is
bit for bit that of the pilot runs of the same solver (55 snapshots each, max\|diff\|
0.0). The reading of row 1 is that the data are consistent with H0 over 9 h on this
column: the two solvers give the same solution to round-off.

**D-2, row 1 met.** S01F3 against S01L3, 296 shared snapshots to 5 h. Up to 2 h: max
`rel_rho` 3.96e-13, `rel_u` 5.03e-13, `rel_P` 2.77e-14, `rel_T` 2.72e-13 (all at most 1e-11)
and `rel_w` 2.98e-10 (0.26 h; at most 1e-8). The `zh3` crossing times agree: ground
`|T/T_0 - 1| > 0.5` at 2.5232375 h in both runs, top `|w| > 5` m/s at 4.3646013 h in both
(differences of order 1e-12 h, far below one scored snapshot); no free-column crossing
inside 5 h in either run, as expected. The 4.36 h crossing, registered with low confidence,
is met. Reported without a band (descriptive): `rel_max` first exceeds 1e-9 at 4.65 h and
does not reach 1e-6 inside 5 h; its maximum is 3.25e-9 (`w`, 4.75 h), with `rel_u` 1.1e-9,
`rel_T` 6.8e-10, `rel_rho` 2.2e-10 and `rel_P` 5.2e-11 at their maxima, all after 4.6 h. A
fit of `ln rel_max` over 2 to 5 h gives a growth rate of 1.63 per hour (e-folding time
0.61 h). The absolute `w` difference is 3e-11 m/s up to 3 h, then 2.9e-10 m/s at 4 h and
3.5e-8 m/s at 5 h. The snapshot times first differ at step 300 (equal at step 200), so
the `dt` sequences part between steps 201 and 300. S01F3 and S01R3 each refactorized 400
times, 13.5 per 1000 steps. Hourly `MF(30.25)` [µPa], S01L3 and S01F3: 2.5086 (0.5 h),
12.966 (1.5 h), 28.471 (2.5 h), 52.112 (3.5 h) and 180.47 (4.5 h), with relative differences
of at most 1.5e-12 up to 3.5 h and 3.3e-11 at 4.5 h. The growth after 2 h, in the run that
approaches breaking and nowhere in segment A, is consistent with H2, round-off
amplification by the flow, and is not by itself evidence on the solver.

**D-3, met.** S01R1 against S01L1, 492 shared snapshots, and S01R3 against S01L3, 296
shared snapshots: max\|diff\| exactly 0.0 in both. `:lu_reuse` is bit-identical to `:lu`
over both segments, including the 400 refactorizations of segment B.

**D-4, cost (no band).** From the run logs, with the `real` time of each run, and from
`d1.cost.tsv` over 00:16:06 to 00:52:05 (2159 s), while all six runs were active:

| ratio | expected | whole run (`real`) | all six active |
|---|---|---|---|
| `:fft_tridiag` / `:lu`, segment A | about 0.80 | 0.781 | 0.789 |
| `:fft_tridiag` / `:lu`, segment B | about 0.60 | 0.666 | 0.719 |
| `:lu_reuse` / `:lu`, segment A | about 1.0 | 1.070 | 1.071 |
| `:lu_reuse` / `:lu`, segment B | about 0.9 | 0.962 | 1.000 |

The whole-run ratios mix concurrency levels, since the segment B runs ended between 00:52
and 01:11 and the segment A runs then ran at concurrency 3 to 1; the all-six column holds
concurrency fixed. On segment B the saving of `:fft_tridiag` is smaller than expected
(0.67 to 0.72 against 0.60): up to 5 h the refactorizations are 13.5 per 1000 steps, fewer
than in the S15F2 breaking window that set the estimate. On segment A, where both
`:lu` variants refactorize twice, `:lu_reuse` ran 7 % slower than `:lu` at fixed
concurrency, which the refactorization count does not explain; it may reflect a
difference in the solve path or core placement and is not followed further here. The mean
CPU power over the all-six interval was 39 W.

Summary of outcomes. D-0, D-1 (row 1 and the second statement), D-2 (row 1) and D-3 are
met; no prediction failed. Under the registration, `:fft_tridiag` can replace `:lu` on the
segment A column and gives the same breaking-flow diagnostics on segment B, and `:lu_reuse`
is a bit-identical interim. The default stays `:lu`: changing the default, or quoting a
production flux from `:fft_tridiag`, is a separate decision. Variant L is not implemented,
by its registered rule, since D-1 did not fail.
