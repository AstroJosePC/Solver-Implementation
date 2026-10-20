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

**Note added 2026-10-20, the diffusion pack refactor and the `operator_split` attribute.**
Two of the open items after these results were closed on `dev/fft-production`; neither
changes a reading above, and the readings stay those of the code on `93a716a`.
- Open item b, `0f6ff62` (`src: preallocate the diffusion pack and unpack buffers`). The six
  `[nx, nz]` fields of `diffusion_step!` are preallocated `DiffusionSolver` buffers, the pack
  and unpack loops are typed kernels behind a function barrier (the unpack loop had
  dispatched dynamically on `solve_cn!` results typed `Any`), and `build_cn_rhs!` selects
  the periodic neighbours without `mod1`; every floating-point expression and its order are
  unchanged, and there is no switch. One `diffusion_step!` on the dz 500 m grid (median of
  30 calls, no refactorization, load average about 17): `:lu` 0.0710 to 0.0250 s,
  `:lu_reuse` 0.0765 to 0.0248 s, `:fft_tridiag` 0.0529 to 0.0040 s per call, with
  allocations of 97.3 MB to 9.1 MB (`:lu`, `:lu_reuse`) and 88.2 MB to 0 B
  (`:fft_tridiag`). Bitwise evidence: a scratch harness comparing the old and new
  conservative arrays (36 arrays, all three solvers, steps with and without
  refactorization, max\|diff\| 0.0); and, on hash `1b99274a03cc2a37`, 0.1 h of the S01F1,
  S01R1, S01F3 and S01R3 commands against the archived `output/d1` records and of S14M1 in
  `check`, max\|diff\| exactly 0.0 in every shared snapshot (Section 2, Check before
  launch). Measurements and logs are in Section 2 item b.
- Open item c, `c5e17e2` (`src: write operator_split to the snapshot attributes`).
  Snapshots carry `operator_split` (string) and `enable_viscosity`, `enable_cooling`
  (integers, 1 on, 0 off), and `face_probe.nc` carries `operator_split`. Data variables are
  unchanged: `check` on the same hash reproduced 34 snapshots at 0.0. A new snapshot reads
  `operator_split = "lie"` on the default path and `"strang"` with `opstrang=1`.

The cost ratios of D-4 were measured before `0f6ff62`, which speeds up `:lu` as well, and
are not the ratios of the current code.


## 2. Production use of the spectral-tridiagonal solver, fft twins of the PR1 face-path cases (pre-registration)

Written on 2026-10-20, before any run of the set below was launched. To be committed on its
own and tagged `prereg/s02-fft-production`; the tag timestamp precedes the launch. Nothing
registered in Section 1 is changed.

Policy, unchanged and confirmed for this section: `SimConfig.diffusion_solver` keeps its
default `:lu` and `SimConfig.diffusion_dt_tol` its default 0.01, so the default path keeps
reproducing every earlier record and the launcher `check` keeps judging it. Production use
of `:fft_tridiag` means that production launchers pass `diffsolver=fft` explicitly, in the
same way as `bcface=1` and `pbar=1`. No flux is quoted from a run of this set before its
readings are scored, and the decision rule below says which readings allow production
launchers to pass `diffsolver=fft`.

The change under test is the one degree of freedom of Section 1, the linear solver of the
Crank-Nicolson diffusion step, at fixed equations, `dt` rule and refactorization rule.
Section 1 scored it on two segments of the dz 500 m column, up to 9 h (U40, `h` = 0.1 km)
and 5 h (U57, `h` = 3 km). This section extends the comparison to the configurations that
production uses: the 24 h PR1 face-path cases at the other `(U, h)` pairs, and the 250 m
grid that Section 11 of `reflection_convergence_tests.md` requires as the resolution twin.
Each `:fft_tridiag` run here is scored against an archived `:lu` record of the same command
without `diffsolver`, so the set contains no new `:lu` run.

### What is settled before this registration

**a. Section 1 results.** Scored on 2026-10-20 (`53558a8`). On segment A, 492 shared
snapshots to 9 h, the largest `rel` values were `rho` 1.39e-13, `u` 1.11e-12, `P` 1.72e-14,
`T` 2.71e-13 and `w` 9.96e-9 (0.35 h, during the ramp), and the 6 to 9 h `MF(30.25)`
changed by −2.55e-12; step counts were equal (49 153). The absolute `w` difference stayed
between 2.7e-11 and 5.3e-11 m/s after 1 h with no trend, on both segments, at forcing
amplitudes `U h` of 4 and 171 km m/s. On segment B, every variable stayed at round-off up
to 2 h, the `zh3` crossings at 2.52 h and 4.36 h agreed to about 1e-12 h, and `rel_max`
grew after 2 h with an e-folding time of 0.61 h, first exceeding 1e-9 at 4.65 h
(consistent with H2 of Section 1). Cost, whole run: `:fft_tridiag`/`:lu` 0.78 (A) and 0.67
(B), at concurrency 6.

**b. The shared pack and unpack, refactored without a switch (`0f6ff62`).** The six
`zeros(nx, nz)` per call and the reconstruction loop, which dispatched dynamically on
`solve_cn!` results typed `Any`, were replaced by preallocated `DiffusionSolver` buffers and
two typed kernels, and `build_cn_rhs!` selects the periodic neighbours without `mod1` and
computes the face diffusivities once per level. Every floating-point expression and its
evaluation order are unchanged,
so the change has no switch. Evidence, from a scratch harness that saves the conservative
arrays of the old code and compares the new ones element by element (`:lu`, `:lu_reuse` and
`:fft_tridiag`; steps `dt0`, `1.004 dt0`, `0.8 dt0`; dz 500 m random state with the CFL
`dt` and with `dt0` = 60 s, the U40 `h` 1 km face snapshot at 11.78 h, dz 250 m random
state): 36 arrays bitwise identical, max\|diff\| 0.0. `scripts/check_diffusion_solvers.jl`
on `0f6ff62` reproduces the Section 1 maxima (`:lu_reuse` 0.0; `:fft_tridiag` `rel_field`
3.4e-14 and 4.0e-14, `rel_incr` 5.6e-9 and 2.9e-11 at dz 500 and 250 m).

Measured cost of one `diffusion_step!` on the dz 500 m grid (snapshot state, no
refactorization, median of 30 calls after warm-up, julia -t 1, `OPENBLAS_NUM_THREADS=1`,
with about 16 cores busy with other jobs, load average about 17, old and new measured about
2 min apart):

| solver | old [s] | new [s] | ratio | allocated per call, old / new |
|---|---|---|---|---|
| `:lu` | 0.0710 | 0.0250 | 0.353 | 97.3 MB / 9.1 MB (the three backsolve results) |
| `:lu_reuse` | 0.0765 | 0.0248 | 0.324 | 97.3 MB / 9.1 MB |
| `:fft_tridiag` | 0.0529 | 0.0040 | 0.076 | 88.2 MB / 0 B |

The unpack loop fell from 0.0440 s to 0.00024 s per call (it was 83 % of the old
`:fft_tridiag` step and 62 % of the old `:lu` step). The three solves are unchanged
(0.0229 s `:lu`, 0.0027 s `:fft_tridiag`) and are now 92 % of the `:lu` step and 66 % of the
`:fft_tridiag` step. The pieces sum to the measured step within 0.0006 s. At dz 250 m,
`check_diffusion_solvers.jl` gives 0.008 to 0.010 s of wall clock per `:fft_tridiag` step
including the field extraction, against 0.050 to 0.059 s for `:lu` without refactorization
and 2.4 to 2.6 s with one. The refactor speeds up `:lu` as well, so the cost ratio of the two solvers
on the new code is not the Section 1 ratio. The six buffers add `6 nx nz 8` bytes per
`DiffusionSolver` (9 MB at 1000 x 190, 18 MB at 1000 x 380), allocated also when viscosity
is off. Logs: `~/MAGNUSP_runs/d2/prereg/pack/` (`harness.jl`, `harness_{old,new}.log`,
`check_diffusion_solvers_new.{log,json}`).

**c. Snapshot attributes (`c5e17e2`).** Snapshots now also carry `operator_split`,
`enable_viscosity` and `enable_cooling` (the last two as Int, 1 on, 0 off). Data variables,
dimensions and the solution are unchanged. Reproduction is judged on the variables (`rho`,
`u`, `w`, `P`, `T`, `time_s`), as in Section 1.

**d. The off members.** Every `:lu` member is an archived record. Its archived log header
was read on 2026-10-20 (`output/reflection_tests/logs/<run ID>.log`) and matches the twin
in every option except `diffsolver`: CFL 0.40, K0 0.005, LZ 95 km, `pbar_mode` pressure,
`ramp` tanh(t/2 h), ground `w` on the face, `LX` 4000 km, `sigma` 100 km, `dx` 4000 m, sweep
order `xz`, sponge inline, face probe off, and the `dz`, `outint` and `t_end` of the twin.
The s15 records ran on `671bf32` (clean, driver and `src/` hash `c503a6c106887e07`) at
concurrency 12; S14S2 ran in set `a3`, before the provenance headers, and its header block is
the driver's own configuration printout. The diffusion solver of all of them is the global
LU of the code of that time, which the merged default path reproduces: Section 1 D-0 (S14M1
to 9 h and PR1-U57-h3.0 to 5 h, max\|diff\| 0.0) and, for S14S2, the check below.

S14S2 as the off member of the dz 250 m pair. The registered pair is an on/off pair at
dz 250 m on the face path, `h` = 0.1 km, 9 h. The archived S14S2
(`output/arsia_mons_stage2_9h_K0_0.005_LZ95_pbarP_ramp2h_dz250_bcface/run_U40.0_h0.1`,
command `hours=9 pbar=1 ramp=2 k0=0.005 lz=95 u=40 h=0.1 bcface=1 dz=250 outint=200`,
98 768 steps, 18 322 s `real` at concurrency 11) is exactly that configuration on `:lu`.
0.1 h of its command on `0f6ff62`, run on 2026-10-20 at 05:20 to scratch (1098 steps, 104 s
`real`, julia -t 1, `OPENBLAS_NUM_THREADS=1`, load average about 17), reproduced it with
max\|diff\| exactly 0.0 in all 6 shared snapshots
(`~/MAGNUSP_runs/d2/prereg/S14S2_repro_0.1h.{log,bitcheck.log}`). S14S2 is therefore reused
as the off member, and no new `:lu` run is made. The 0.1 h reproduction is repeated on the
launch hash (Check before launch); if it fails there, the set gains the `:lu` member S02L1
under a dated note before launch.

**e. First overturning in the off members.** `zoverturn --stride 2` (definitions of
Section 15 of `reflection_convergence_tests.md`) on the archived records, from the s15
products (`output/analysis/s15/zoverturn_F.json`, `zoverturn_PR1_h1.json`) and from new runs
of the same command on 2026-10-20 (`~/MAGNUSP_runs/d2/prereg/zoverturn_{S14S2,PR1_U20_h0.1,PR1_U57_h0.1,PR1_h3}.{json,log}`,
one process at a time). `t_c` is the earliest first onset over all bands and both kinds.

| off member | twin | `t_c` [h] | event at `t_c` | later onsets (descriptive) |
|---|---|---|---|---|
| S14S2 | S02F1 | none to 9 h | none | none |
| PR1-U20-h0.1 | S02F2 | none to 24 h | none | none |
| PR1-U57-h0.1 | S02F3 | 10.69 | convective, 0 to 5 km, at 1.0 km | convective 45 to 55 km at 15.75 h (54.5 km); none in 5 to 45 km |
| PR1-U20-h1.0 | S02F4 | 4.58 | convective, 0 to 5 km, at 1.0 km | kinematic 45 to 55 km 6.77 h; free column 15.34 h |
| S15F1 (PR1-U40-h1.0) | S02F5 | 3.18 | convective, 0 to 5 km, at 1.0 km | kinematic 45 to 55 km 4.03 h; free column 13.02 h |
| PR1-U57-h1.0 | S02F6 | 2.67 | convective, 0 to 5 km, at 1.0 km | kinematic 45 to 55 km 3.15 h; free column 10.15 h |
| PR1-U20-h3.0 | S02F7 | 2.64 | convective, 0 to 5 km, at 1.0 km | kinematic 45 to 55 km 3.04 h; free column 9.52 h |
| PR1-U40-h3.0 | S02F8 | 1.90 | convective, 0 to 5 km, at 1.0 km | kinematic 45 to 55 km 2.12 h; free column 5.66 h |
| S15F2 | S02F9 | 2.70 | convective, 0 to 5 km, at 0.5 km | kinematic 45 to 55 km 3.64 h; free column 10.92 h |

The `t_c` event lies at 2 dz in every breaking record, inside the near-ground layer, and
precedes the `zh3` ground crossing of the `h` = 3 km records (5.27 h for U20, 3.40 h for U40).
Two of the three `h` = 0.1 km records show no overturning; PR1-U57-h0.1 shows the near-ground
event at 10.69 h and one above the 53.3 km critical level at 15.75 h, so its 18 to 24 h
window lies after `t_c`.

### Runs

Common arguments `P24` = `hours=24 pbar=1 ramp=2 k0=0.005 lz=95 bcface=1`
(`sweep_arsia_mons_stage2.jl`, CFL 0.4). Every run passes `diffsolver=fft` (tag `_dfft`,
after the `_dz` tag and before `_bcface`) and `bcface=1` explicitly. All nine go in one
launcher set `d2`, one process each, `julia -t 1`, `OPENBLAS_NUM_THREADS=1`, concurrency 9,
under `caffeinate -s -i` with the cost monitor. The set writes under
`MAGNUSP_OUTPUT=~/MAGNUSP_runs/d2` and archives to `output/d2/`. It is launched from the
worktree of `dev/fft-production`, whose git-ignored `output` is to be a symbolic link to the
archive `output/` of the main checkout (as for `d1`), so `check` reads the archived
references and `archive` writes to the mount. Archive paths of the off members below are
relative to the main checkout; `S24` = `output/arsia_mons_stage2_24h_K0_0.005_LZ95_pbarP_ramp2h_bcface`.

| id | regime | arguments | duration | off member | off member archive | off member steps, `real` [s] |
|---|---|---|---|---|---|---|
| S02F1 | non-breaking, dz 250 m | `hours=9 pbar=1 ramp=2 k0=0.005 lz=95 u=40 h=0.1 bcface=1 dz=250 outint=200 diffsolver=fft` | 9 h | S14S2 | `output/arsia_mons_stage2_9h_K0_0.005_LZ95_pbarP_ramp2h_dz250_bcface/run_U40.0_h0.1` | 98 768, 18 322 |
| S02F2 | non-breaking | `P24 u=20 h=0.1 diffsolver=fft` | 24 h | PR1-U20-h0.1 | `S24/run_U20.0_h0.1` | 131 070, 13 462 |
| S02F3 | non-breaking to 10.69 h | `P24 u=57 h=0.1 diffsolver=fft` | 24 h | PR1-U57-h0.1 | `S24/run_U57.0_h0.1` | 131 079, 13 612 |
| S02F4 | breaking, `h` = 1 km | `P24 u=20 h=1 diffsolver=fft` | 24 h | PR1-U20-h1.0 | `S24/run_U20.0_h1.0` | 131 223, 13 555 |
| S02F5 | breaking, `h` = 1 km | `P24 u=40 h=1 diffsolver=fft` | 24 h or `dt` collapse (23.81 h in the off member) | S15F1 | `S24/run_U40.0_h1.0` | 160 044, 23 815 |
| S02F6 | breaking, `h` = 1 km | `P24 u=57 h=1 diffsolver=fft` | to `dt` collapse (11.99 h) | PR1-U57-h1.0 | `S24/run_U57.0_h1.0` | 78 071, 12 850 |
| S02F7 | breaking, `h` = 3 km | `P24 u=20 h=3 diffsolver=fft` | to `dt` collapse (11.45 h) | PR1-U20-h3.0 | `S24/run_U20.0_h3.0` | 72 106, 10 302 |
| S02F8 | breaking, `h` = 3 km | `P24 u=40 h=3 diffsolver=fft` | to `dt` collapse (7.02 h) | PR1-U40-h3.0 | `S24/run_U40.0_h3.0` | 45 200, 6 593 |
| S02F9 | breaking, `h` = 1 km, dz 250 m | `P24 u=40 h=1 dz=250 outint=200 diffsolver=fft` | 24 h or `dt` collapse; off member to 11.73 h only | S15F2 | `output/arsia_mons_stage2_24h_K0_0.005_LZ95_pbarP_ramp2h_dz250_bcface/run_U40.0_h1.0` | 142 400 to 11.73 h (stopped), 35 227 |

Directories under `~/MAGNUSP_runs/d2` (from `dry d2`):
`arsia_mons_stage2_9h_K0_0.005_LZ95_pbarP_ramp2h_dz250_dfft_bcface/run_U40.0_h0.1`,
`arsia_mons_stage2_24h_K0_0.005_LZ95_pbarP_ramp2h_dfft_bcface/run_U<U>.0_h<h>` for S02F2 to
S02F8, and `arsia_mons_stage2_24h_K0_0.005_LZ95_pbarP_ramp2h_dz250_dfft_bcface/run_U40.0_h1.0`.
None of them is the name of an archived record, so no record is overwritten.

Not repeated here: PR1-U40-h0.1 (Section 1 segment A, to 9 h, on the same column and
command apart from `hours`) and PR1-U57-h3.0 (Section 1 segment B, to 5 h of its 5.47 h).
S15T1 (CFL 0.2) and S14V3p (Queney driver) are not production cases and are not twinned.
S02F9 is the one case of this set outside the PR1 list. It is the `h` = 1 km twin at
dz 250 m that Section 15 needs and that `:lu` could not finish (Section 15 note on stopping
S15F2: about 1.5 s per step at 11.5 h, 50 h or more for the remaining 12.3 h).

### Cost table

Basis. The `:fft_tridiag` rows of `d1` ran at 0.0756 s per step (S01F1) and 0.0763 s per
step (S01F3) at concurrency 6. Their standalone model was 0.013 s outside diffusion plus
0.0529 s of `diffusion_step!`, 0.066 s, so the batch ran at 1.15 times the standalone rate.
After `0f6ff62` the standalone model is 0.013 + 0.0040 = 0.017 s at dz 500 m and
0.026 + 0.008 = 0.034 s at dz 250 m (item b; non-diffusion parts from Section 1 item c).
With the same factor 1.15: about 0.020 s per step at dz 500 m and 0.039 s at dz 250 m. The
`:fft_tridiag` refactorization costs 0.14 to 0.27 ms, so breaking runs, which refactorize
often, are costed at the same rate. Step counts are those of the off members.

| id | steps | s per step (est.) | est. wall | disk [GB] |
|---|---|---|---|---|
| S02F1 | 98 768 | 0.039 | 1.07 h | 7.0 |
| S02F2 | 131 070 | 0.020 | 0.73 h | 9.3 |
| S02F3 | 131 079 | 0.020 | 0.73 h | 9.3 |
| S02F4 | 131 223 | 0.020 | 0.73 h | 9.3 |
| S02F5 | 160 044 | 0.020 | 0.89 h | 11.0 |
| S02F6 | 78 071 | 0.020 | 0.43 h | 5.5 |
| S02F7 | 72 106 | 0.020 | 0.40 h | 5.1 |
| S02F8 | 45 200 | 0.020 | 0.25 h | 3.2 |
| S02F9 | 142 400 to 11.73 h; about 585 000 to 24 h if `dt` stays near 0.10 s | 0.039 | 1.54 h to 11.73 h; about 6.3 h to 24 h | 10 to 11.73 h; up to about 44 |

The `dz` 250 m `h` = 1 km step is modelled at about 0.40 h of wall clock per simulated hour
in the 10.6 to 11.7 h window, where `dt` is about 0.1 s. The wall clock of the set is set by
S02F9, between about 1.5 h (a `dt` collapse soon after 11.73 h) and 6.3 h. The CPU time
summed over the nine runs is 6.8 to 11.5 h, against 41 h of `real` for the nine off
members. Disk is 70 to 104 GB, against 3.2 TB free locally and 1.4 TB free on the archive
mount; the archive copy takes 0.9 to 1.3 h at the a3 rate (31 GB in 1370 s).

Concurrency. The machine had a load average of about 17 from other jobs on 2026-10-20.
Nine more processes would then exceed the 16 performance cores, and the estimates above
would not hold. The launch waits until the load average is below about 7, or the set is
launched with the measured load recorded in `d2.start`, in which case the cost readings are
descriptive of that load only. The concurrency and the load average at launch are recorded
either way.

### Hypotheses

- H0, same discrete system, now at dz 250 m and over 24 h: `:fft_tridiag` and `:lu` differ
  only by round-off in the diffusion solve. Prediction: wherever the flow has not overturned,
  `rel` of `rho`, `u`, `P` and `T` stays of order 1e-13 to 1e-12 and the absolute `w`
  difference at about 5e-11 m/s with no trend, so that `rel_w` scales as `1/(U h)` (about
  2e-8 at most for U20 `h` 0.1 km, half the forcing amplitude of S01F1), the step counts are
  equal, and window-mean fluxes agree to 1e-7 or better.
- H1, different discrete system in a configuration Section 1 did not cover (for example a
  dz 250 m row, a state reached only after 9 h, or a refactorization pattern of the 24 h
  runs). Prediction: a difference of the size of the operator error times the accumulated
  diffusion increment, present from the first steps and independent of the flow regime. The
  unit check bounds it at `rel_incr` 2.9e-11 per step at dz 250 m.
- H2, round-off amplification by the flow: in a run that overturns, a round-off difference
  grows with the flow's own sensitivity. Prediction: round-off up to `t_c`, then roughly
  exponential growth with an e-folding time of order the 0.61 h of S01F3, as any round-off
  perturbation of the `:lu` run would grow. Compatible with H0.
- H3, slow secular growth in a non-breaking run: the round-off difference keeps growing over
  9 to 24 h without an overturning onset. Prediction: in S02F2 the power-law slope of `rel`
  of `rho`, `u` and `T` over 9 to 24 h is 1 or more, so that `rel_u` reaches about 4e-12 by
  24 h (1.11e-12 at 6 h grown linearly), against a slope of at most about 0.43 under H0
  (Section 1 pilot) and about 2e-12 at 24 h. H3 is compatible with production use if the
  bands hold; it bounds the length of runs over which the equivalence is shown.

### Readings fixed before launch

Each reading is recorded against these statements, and a failed prediction is recorded as
failed. `rel` of a variable is `max|f_new - f_ref| / max|f_ref|` at a shared snapshot, over
all cells, for `rho`, `u`, `w`, `P` and `T` separately, with `rel_max` the largest of the
five; `t_c` is the off member's value in item e; the `MF(30.25)` change is the relative
change of the window mean read from `s14_suite.json` at full precision; windows are the hour
rows of `suite` inside the window, and every window lies outside the 2 h ramp. No ratio of
fluxes is banded; any ratio reported is normalized above 16 km. A value within a factor 2
of a band edge is reported at the edge, and the reading of the higher band applies. All bands
are those of Section 1 D-1 and D-2, unchanged.

**D-0, prerequisites, judged before launch.** (i) `check` on the launch hash reproduces the
reference records with max\|diff\| exactly 0.0; (ii) 0.1 h of the S14S2 command on the launch
hash reproduces S14S2 with max\|diff\| exactly 0.0; (iii) `dry d2` resolves the nine
directories above. Met: the set may be launched with the off members as registered. Failed:
no launch; a failure of (ii) alone adds S02L1, the `:lu` member at dz 250 m, under a dated
note before launch. Expected: met (item b, item d).

**D-1, non-breaking pairs (S02F1 against S14S2, S02F2 against PR1-U20-h0.1, S02F3 against
PR1-U57-h0.1).**

| outcome | reading |
|---|---|
| at every shared snapshot with `t` at most `t_c` (to the end of S02F1 and S02F2; to 10.69 h in S02F3), `rel` of each of `rho`, `u`, `P`, `T` at most 1e-11 and `rel_w` at most 1e-7, and the relative \|change\| of `MF(30.25)` at most 1e-7 over 6 to 9 h (all three) and over 18 to 24 h (S02F2, S02F3) | consistent with H0: the two solvers give the same solution to round-off on this case, at dz 250 m for S02F1 and over 24 h for S02F2 and S02F3 |
| any of `rho`, `u`, `P`, `T` in (1e-11, 1e-8], or `rel_w` in (1e-7, 1e-4], or an `MF(30.25)` change in (1e-7, 1e-4], and none above | a departure above the calibrated round-off and below any quoted error bar; reported, and the case is read under the decision rule |
| any of `rho`, `u`, `P`, `T` above 1e-8, or `rel_w` above 1e-4, or an `MF(30.25)` change above 1e-4 | consistent with H1: the solvers differ on this case |

The 18 to 24 h `MF(30.25)` item of S02F3 lies after its `t_c`. If it alone falls outside row
1 while the field bands hold to `t_c` and `rel_max` grows after `t_c` (descriptive fit), the
outcome is recorded as failed for that item and read as H2, not H1.

Second statement: the step counts of each twin and its off member are equal (98 768,
131 070, 131 079). Reported without a band: `rel` per variable against time over the whole
run; for S02F3 after 10.69 h, the first times at which `rel_max` exceeds 1e-9 and 1e-6 and
an e-folding time fitted to any growth; the power-law slope of `rel` of `rho`, `u` and `T`
over 9 to 24 h in S02F2 and S02F3 (H3); the changes of `MF(1.25)`, `MF(10.25)`, `MF(50.25)`,
`R10/30` and `R50/30` over both windows; the factorization counts printed by the twins; the
first-overturning times of the twins (`zoverturn`). Expected: row 1 and the second statement
met for S02F1 and S02F2, with moderate confidence for S02F1, the first whole run of
`:fft_tridiag` at dz 250 m. For S02F3, row 1 with moderate confidence for the fields to
10.69 h and low to moderate confidence for the 18 to 24 h `MF(30.25)` item, which follows a
near-ground overturning (10.69 h) and one above the critical level (15.75 h). H3 is expected
not to hold (slope below 1), with low confidence, since Section 1 covered only 9 h.

**D-2, breaking pairs (S02F4 against PR1-U20-h1.0, S02F5 against S15F1, S02F6 against
PR1-U57-h1.0, S02F7 against PR1-U20-h3.0, S02F8 against PR1-U40-h3.0, S02F9 against S15F2).**
No `h` = 3 km flux is quoted (Section 14 classification of `reflection_convergence_tests.md`,
unchanged).

| outcome | reading |
|---|---|
| at every shared snapshot with `t` at most `t_c`, `rel` of each of `rho`, `u`, `P`, `T` at most 1e-11 and `rel_w` at most 1e-8, and the `t_c` event of the twin (same band, kind and height) at a time within one scored snapshot (stride 2, about 132 s on both grids) of the off member's | consistent with H0, with H2 for any growth after `t_c`; `:fft_tridiag` gives the same flow as `:lu` up to the first overturning |
| any of `rho`, `u`, `P`, `T` above 1e-11 or `rel_w` above 1e-8 at some `t` at most `t_c`, or the `t_c` event of the twin absent or more than one scored snapshot away | the solvers differ in this regime (H1), or the flow amplifies round-off before the first overturning; reported with the growth curve |

Reported without a band (descriptive), per pair: `rel` per variable against time over the
shared snapshots; the first times at which `rel_max` exceeds 1e-9, 1e-6 and 1e-3; an
e-folding time fitted to the growth after `t_c`, compared with the 0.61 h of S01F3; the first
step at which the `dt` sequences part; the factorization counts per 1000 steps of the twin;
every `zoverturn` onset of both runs, and for the `h` = 1 km cases `t_ov` and `dz_L` of both
runs and the window means of `MF(30.25)` over the off member's pre-collapse window
(12.34 to 15.34 h for PR1-U20-h1.0, 10.02 to 13.02 h for S15F1, 7.15 to 10.15 h for
PR1-U57-h1.0, 7.92 to 10.92 h for S15F2) with their relative difference and the hourly range
of each; for the `h` = 3 km cases the `zh3` crossings and the end of both runs, against the
H3 table of Section 15; the end time and step count of every twin. S02F9 is compared with
S15F2 only on the shared snapshots to 11.73 h. Its record after 11.73 h has no off member and
is reported (end, onsets, hourly `MF(30.25)`) without a reading; using it for a flux or for
the F-b resolution row of Section 15 needs a separate registration. Expected: row 1 in all
six pairs, with moderate confidence for S02F8, whose `t_c` (1.90 h) lies inside the ramp,
and for S02F9 at dz 250 m. After `t_c` the differences are expected to grow (H2), so the
later onsets, the pre-collapse window means and the collapse times may differ by more than
one snapshot or by more than round-off; no expectation is registered for their size.

**D-3, cost (no band).** From the run logs, `d2.cost.tsv` and the `d2` rows of
`docs/run_matrix.csv`: `real` time and s per step of each twin at concurrency 9, with the
load average at launch; the ratio of each twin's `real` to its off member's `real` (which
mixes `0f6ff62`, which also speeds up `:lu`, concurrency 9 against 12 and the load, and is
therefore descriptive of the production gain only); and for S02F9 and S15F2 the wall clock
per simulated hour over 10.6 to 11.7 h, from the modification times and `time_s` of the
first snapshots at or after 10.6 h and 11.7 h in each record. Expected: about 0.020 s per step
at dz 500 m and 0.039 s at dz 250 m on an otherwise idle machine, twin/off member `real`
ratios of about 0.2 for the 24 h `h` = 0.1 km cases, and about 0.4 h per simulated hour for
S02F9 in the 10.6 to 11.7 h window against about 15 h for S15F2 (1.5 s per step at about
0.1 s of `dt`). No band: cost does not decide production use, which rests on D-0 to D-2.

### Decision rule for production use

The `SimConfig` default stays `:lu` in every outcome. What changes is whether new production
launchers pass `diffsolver=fft`.
- D-0 met, D-1 row 1 and the second statement met in all three pairs, and D-2 row 1 met in all
  six pairs: new production sets pass `diffsolver=fft` together with `bcface=1` and `pbar=1`,
  on both grids and for every `h`. This is implemented as its own `launcher:` commit after a
  dated note in this section; no existing set changes, and `check` keeps judging the `:lu`
  default.
- As above, except that the 18 to 24 h `MF(30.25)` item of S02F3 fails while its fields hold
  to `t_c` (read as H2): the same adoption. Late-window fluxes from runs that have overturned
  near the ground then carry the measured S02F3/PR1-U57-h0.1 change as a solver row of the
  error budget (Section 11 table of `reflection_convergence_tests.md`; rows are not summed).
- D-1 met in all three pairs and D-2 row 2 in any pair: production passes `diffsolver=fft`
  for runs with no overturning (the `h` = 0.1 km cases, with S02F3's rule above), and runs
  that overturn keep `:lu` until the difference is located.
- D-1 row 2 or 3 in any pair: production keeps `:lu`, and the difference is located before
  any further use. Variant L of Section 1 is then reconsidered under its own registration.
- Quoting, in every outcome. Values already quoted from `:lu` records (Section 15 P, F-a to
  F-c, H3) stay as quoted. A value quoted from a `:fft_tridiag` record names the solver along
  with the forcing path, `dt` or CFL, `dz`, window and run ID. For a run that overturns, the
  descriptive D-2 differences after `t_c` (onset times, pre-collapse window means, collapse
  times) are recorded as a solver row of the error budget, kept separate from the sampling
  scatter.

### Scoring commands

All read archived output only (`analyze_dt_convergence.py` unless stated), with the products
written to `output/analysis/s02/`. The run directories of this set are under `output/d2/`;
`F24` = `output/d2/arsia_mons_stage2_24h_K0_0.005_LZ95_pbarP_ramp2h_dfft_bcface` and `S24`
as under Runs. The commands run one after another in one shell, with
`PY=~/.venvs/magnusp/bin/python`, so that no two processes open a file on the mount at the
same time (the d1 scoring lost seven commands to `NetCDF: HDF error` when two processes
opened the same file). A command that still stops with an HDF error is rerun unchanged after
the previous one has ended, and the incident is recorded.
- D-1:
  `bitcheck --relative --new output/d2/arsia_mons_stage2_9h_K0_0.005_LZ95_pbarP_ramp2h_dz250_dfft_bcface/run_U40.0_h0.1 --ref output/arsia_mons_stage2_9h_K0_0.005_LZ95_pbarP_ramp2h_dz250_bcface/run_U40.0_h0.1 --out output/analysis/s02/rel_F1.json`;
  `bitcheck --relative --new F24/run_U20.0_h0.1 --ref S24/run_U20.0_h0.1 --out output/analysis/s02/rel_F2.json`;
  `bitcheck --relative --new F24/run_U57.0_h0.1 --ref S24/run_U57.0_h0.1 --out output/analysis/s02/rel_F3.json`;
  `suite --run S14S2=<off> --run S02F1=<twin> --window 6 9 --pair S14S2:S02F1 --outdir output/analysis/s02/mf_F1_6_9`;
  `suite --run PR1-U20-h0.1=S24/run_U20.0_h0.1 --run S02F2=F24/run_U20.0_h0.1 --window 6 9 --pair PR1-U20-h0.1:S02F2 --outdir output/analysis/s02/mf_F2_6_9`
  and the same with `--window 18 24 --outdir output/analysis/s02/mf_F2_18_24`; the same two
  for S02F3 against PR1-U57-h0.1 (`mf_F3_6_9`, `mf_F3_18_24`); step counts from
  `grep -h "Simulation complete" output/d2/logs/S02F[123].log` against item d.
- D-2: `bitcheck --relative --new <twin> --ref <off member> --out output/analysis/s02/rel_F<n>.json`
  for n = 4 to 9 (pairs as under Runs).
- Onsets, every twin: `zoverturn --run S02F<n>=<twin> --stride 2 --out output/analysis/s02/zoverturn_F<n>.json`
  for n = 1 to 9; the off members' values are those of item e. `zh3 --run F24/run_U20.0_h3.0 --stride 2 --out output/analysis/s02/z3_F7.json`
  and `zh3 --run F24/run_U40.0_h3.0 --stride 2 --out output/analysis/s02/z3_F8.json`,
  against `output/analysis/s15/z3_face_U20_h3.json` and `z3_face_U40_h3.json`.
- D-2 descriptive window means: `suite --run PR1-U20-h1.0=S24/run_U20.0_h1.0 --run S02F4=F24/run_U20.0_h1.0 --window 12.34 15.34 --pair PR1-U20-h1.0:S02F4 --outdir output/analysis/s02/pre_F4`;
  `suite --run S15F1=S24/run_U40.0_h1.0 --run S02F5=F24/run_U40.0_h1.0 --window 10.02 13.02 --pair S15F1:S02F5 --outdir output/analysis/s02/pre_F5`;
  `suite --run PR1-U57-h1.0=S24/run_U57.0_h1.0 --run S02F6=F24/run_U57.0_h1.0 --window 7.15 10.15 --pair PR1-U57-h1.0:S02F6 --outdir output/analysis/s02/pre_F6`;
  `suite --run S15F2=<off> --run S02F9=<twin> --window 7.92 10.92 --pair S15F2:S02F9 --outdir output/analysis/s02/pre_F9`.
- Factorizations: `grep -h "factorizations" output/d2/logs/S02F*.log`.
- D-3: the `d2` rows of `docs/run_matrix.csv`, `output/d2/logs/d2.cost.tsv` and the logs;
  for the S02F9/S15F2 window, `stat -f '%m %N'` on the snapshots and `ncdump -h` for their
  `time_s`.
- The pre-registration products of items d and e are copied from `~/MAGNUSP_runs/d2/prereg/`
  to `output/analysis/s02/prereg/` with the archive.

### Code for Section 2

On `dev/fft-production`, stacked on `main` (`52c035a`):
- `c5e17e2`: `src/io.jl` writes `operator_split`, `enable_viscosity` and `enable_cooling` to
  the snapshot attributes (Section 1 open item c). No change to the solution.
- `0f6ff62`: `src/diffusion.jl` preallocates the pack and unpack buffers and restructures
  the `build_cn_rhs!` loop (item b; Section 1 open item b). Bitwise identical for all three solvers; no switch.
- `dd3877f`: `run_dt_convergence_S2T95.sh` set `d2` with the nine cases above, run IDs
  `S02F1` to `S02F9`, stamps `d2.start`, `MAGNUSP_OUTPUT` defaulting to `~/MAGNUSP_runs/d2`
  and the archive to `output/d2`; `check` is shared and every earlier set keeps its
  behaviour. `dry d2` (2026-10-20, `~/MAGNUSP_runs/d2/logs/d2.dry.log`) resolved the nine
  cases to `pbar_mode` pressure, the face ground path, CFL 0.40, the `dz`, `outint` and
  `t_end` above, and the directories listed under Runs, with the tag order
  `..._ramp2h_dz250_dfft_bcface` of the driver header.

### Naming

Run IDs use the two-digit section form `S02F<n>` (`F` for `:fft_tridiag`, `n` the case), and
`S02L1` is reserved for a `:lu` member that D-0 (ii) might require. No ID in
`docs/run_matrix.csv` or the docs starts with `S02`, and none collides with the legacy `S0`
to `S2`. The set is `d2`, analysis products go to `output/analysis/s02/`, and the tags are
`prereg/s02-fft-production` and `run/d2-<YYYYMMDD>`.

### Check before launch

All on commit `dd3877f` (tree clean apart from this uncommitted Section 2 text), driver,
Queney driver and `src/` hash `1b99274a03cc2a37`, from the worktree with
`MAGNUSP_OUTPUT=~/MAGNUSP_runs/d2`, `output` a symbolic link to the `output/` of the main
checkout, Python from `~/.venvs/magnusp`, and the pre-addendum-3 Queney copy taken from
`~/MAGNUSP_runs/d1/logs/a3_ref` into `~/MAGNUSP_runs/d2/logs/a3_ref`, on 2026-10-20 under
`caffeinate -s -i` with a load average of about 17 to 18 from other jobs.

D-0 (i), met. `check d2`, 05:35:44 to 05:36:32 (48 s wall, twelve processes): max\|diff\|
exactly 0.0 in every shared snapshot, 0.1 h of S14G1 (5 snapshots), C2R with `probe=5` (5),
S04 (6), C2R (5), S14M1 (6), the Queney driver on its default path against its
pre-addendum-3 copy (2) and S14V1a (5), 34 in all. The Queney face path and the Queney face
path with `pbar=1` completed, and the Section 16 Strang and alternating paths with Strang
operators reached `t_end` (487 steps). Summary `~/MAGNUSP_runs/d2/logs/d2.check.out`, logs
`~/MAGNUSP_runs/d2/logs/S14check_*.log` and `S16check_*.log`; the launch marker
`~/MAGNUSP_runs/d2/logs/S14.check_ok` holds the hash above.

D-0 (ii), met. 0.1 h of the S14S2 command (05:36:32 to 05:38:21, 1098 steps, 98.5 s solver
wall clock), against the archived S14S2: max\|diff\| exactly 0.0 in all 6 shared snapshots.
S14S2 stays the off member of S02F1.

D-0 (iii), met. `dry d2` repeated on `dd3877f` (05:48,
`~/MAGNUSP_runs/d2/logs/d2.dry.dd3877f.log`) is identical to the earlier `d2.dry.log` and
resolves the nine directories under Runs.

Further reproductions on the same hash, run together with the S14S2 command (eight
processes) into `~/MAGNUSP_runs/d2/check_d1/` and compared one at a time with the archives
(05:47:17 to 05:47:24), each max\|diff\| exactly 0.0 in 6 shared snapshots: 0.1 h of the
S01F1 (`:fft_tridiag`) and S01R1 (`:lu_reuse`) commands against `output/d1`, the same for
S01F3 and S01R3 (U57 `h` = 3 km), and 0.1 h of the S15F2 command against S15F2, the off
member of S02F9. Each of the four `d1` runs printed 2 factorizations, so the refactorization
path after the first `dt` change is covered. With S14M1 above (`:lu`), the refactored
`diffusion_step!` of item b reproduces records of all three solvers written by the code
before it. Descriptive, not a reading: the 0.1 h face-path ridge with `opstrang=1` (the
diffusion solve called twice per step over dt/2), `:fft_tridiag` against `:lu`, both on the
new code, gives per-variable `max|diff|` and `rel` equal, number for number (70 of 70 over
the 5 shared snapshots), to those of the Section 1 merge smoke pair
(`output/analysis/s01/merge_smoke/rel_opS.json`, `rel_max` 9.46e-10). A new snapshot
(S01F1 command, step 500) carries `operator_split = "lie"`, `enable_viscosity = 1`,
`enable_cooling = 0` (32-bit integers, as `step`) and `diffusion_solver = "fft_tridiag"`;
the Strang ridge snapshots carry `operator_split = "strang"`. Summary
`~/MAGNUSP_runs/d2/check_d1/bitcheck.out`, run logs `~/MAGNUSP_runs/d2/check_d1/*_0.1h.log`,
`opS_fft.log`, `opS_lu.log`, and `rel_opS_new.json`.

Any edit to the driver or `src/` before the launch requires a new `check`.

**Note added 2026-10-20, S02F9 dropped from set `d2` before launch.** By decision of
2026-10-20, written before any run of the set was started, S02F9 (the dz 250 m, U40,
`h` = 1 km twin of S15F2) is removed from set `d2` (`0caebc7`). Its command and directory
stay defined in the launcher, so a later registration can run it. The driver and `src/` are
unchanged (hash `1b99274a03cc2a37`), so the check of D-0 (i) and (ii) stands; `dry d2` on
`0caebc7` (`~/MAGNUSP_runs/d2/logs/d2.dry.drop9.log`) is identical to `d2.dry.dd3877f.log`
apart from the S02F9 block. The launch stamp now also records the 1, 5 and 15 min load
averages (`aef8f86`). This note corrects the following statements of this section and
changes nothing else.
- Runs: the set is S02F1 to S02F8, eight processes, concurrency 8 (not 9). The directory
  `arsia_mons_stage2_24h_K0_0.005_LZ95_pbarP_ramp2h_dz250_dfft_bcface/run_U40.0_h1.0` is not
  written. D-0 (iii) is read against the eight remaining directories.
- Cost table: the S02F9 row does not apply. The wall clock of the set is set by S02F1, about
  1.07 h at the registered rate on an otherwise idle machine. The CPU time summed over the
  eight runs is about 5.2 h, against 31.3 h of `real` for their eight off members. Disk is
  about 60 GB, and the archive copy takes about 0.7 h at the a3 rate.
- D-2: five pairs (S02F4 to S02F8). The S02F9 items are not scored: its `bitcheck`, its
  `zoverturn` (n = 9), the `pre_F9` window means and, under D-3, the wall clock per
  simulated hour of S02F9 and S15F2 over 10.6 to 11.7 h.
- Decision rule: "D-2 row 1 met in all six pairs" and "D-2 row 2 in any pair" read as five
  pairs. The rule is otherwise applied as registered.
- Scope, a limitation of the reduced set. At dz 250 m the set now covers only the
  non-breaking S02F1 (U40, `h` = 0.1 km, 9 h). Adoption "on both grids and for every `h`"
  would therefore rest, for dz 250 m runs that overturn, on no full-run pair: only on the
  unit check (`rel_incr` 2.9e-11 per step at dz 250 m), on the 0.1 h reproduction of the
  S15F2 command on the default path, and on the dz 500 m breaking pairs. This may bias the
  reading for that regime, because the refactorization pattern and the state near `t_c` of
  a dz 250 m overturning run are not compared. The reading of the decision rule records this
  limitation next to its outcome.

### Results

Scored on 2026-10-20 with the registered commands on the archived output, after the
registration (`bd42b47`, tag `prereg/s02-fft-production`), the note of 2026-10-20 above
(`df52f0e`) and the launch tag `run/d2-20261020`.

Run record. Set `d2`, launched 2026-10-20 07:03:44 from `df52f0e` (clean tree) in the
worktree of `dev/fft-production`, driver and `src/` hash `1b99274a03cc2a37`, launcher copy
`~/MAGNUSP_runs/d2/logs/d2.launcher.sh` (hash `e5b991de0590dccd`), concurrency 8,
`julia -t 1`, `OPENBLAS_NUM_THREADS=1`, cost monitor holding `caffeinate -s -i`. Load
average at launch 5.98, 13.73, 15.73 (1, 5, 15 min): a 16-process DSMC spin-up of the
MarsAtmos queue had ended at 07:02:48, and the launch waited for the 1 min value to fall
below 7. No other heavy process ran during the batch; a hold job in the shared Task Spooler
queue kept the next DSMC job waiting until the last case ended. The last case ended
08:01:41 (0.97 h); no sleep or DarkWake entries in `pmset -g log` over the batch. Mean CPU
power while all eight ran was 59 W (`d2.cost.tsv`; its `P_use` column reads 0.00 throughout,
a parsing fault of the monitor, not a reading). Archived 08:02:37 to 08:57:34 to `output/d2/`,
0 checksum differences and equal file counts in both sweep directories (495 and 6835 files);
logs, stamps, cost history and launcher copy in `output/d2/logs/`. The local copies are kept;
the `rm` command is printed in `~/MAGNUSP_runs/d2/logs/d2.archive.log`. The 26 scoring
commands ran one after another, 08:57:47 to 09:38:49, every one with exit 0 and no HDF error.
Products in `output/analysis/s02/`; the extracted readings are in
`output/analysis/s02/readings_d2.json`.

| id | off member | steps (off member) | `t_end` [h] (off member) | mean `dt` [s] | `real` [s] | s per step | factorizations (per 1000 steps) |
|---|---|---|---|---|---|---|---|
| S02F1 | S14S2 | 98 768 (98 768) | 9.00 (9.00) | 0.328041 | 3476.99 | 0.0352 | 2 |
| S02F2 | PR1-U20-h0.1 | 131 070 (131 070) | 24.00 (24.00) | 0.659189 | 2339.65 | 0.0178 | 2 |
| S02F3 | PR1-U57-h0.1 | 131 079 (131 079) | 24.00 (24.00) | 0.659144 | 2344.11 | 0.0179 | 2 |
| S02F4 | PR1-U20-h1.0 | 131 223 (131 223) | 24.00 (24.00) | 0.658421 | 2333.88 | 0.0178 | 3 (0.02) |
| S02F5 | S15F1 | 94 845 (160 044) | 15.35 (23.81) | 0.582695 | 1713.81 | 0.0181 | 5677 (59.9) |
| S02F6 | PR1-U57-h1.0 | 75 700 (78 071) | 11.91 (11.99) | 0.566580 | 1365.00 | 0.0180 | 4477 (59.1) |
| S02F7 | PR1-U20-h3.0 | 74 152 (72 106) | 11.49 (11.45) | 0.557832 | 1339.90 | 0.0181 | 4491 (60.6) |
| S02F8 | PR1-U40-h3.0 | 44 959 (45 200) | 7.01 (7.02) | 0.561354 | 816.09 | 0.0182 | 1921 (42.7) |

**D-0, met.** (i) and (ii) as recorded under Check before launch, on the unchanged hash;
(iii) `dry d2` on `0caebc7` resolves the eight directories (note of 2026-10-20).

**D-1, row 1 met in all three pairs; second statement met.** Maxima over all shared
snapshots up to `t_c` (all of S02F1 and S02F2; to 10.69 h in S02F3), with the time of each
maximum:

| quantity | S02F1 (494 snapshots, dz 250 m) | S02F2 (1311) | S02F3 (585 to `t_c`) | band of row 1 |
|---|---|---|---|---|
| max `rel_rho` | 3.30e-13 (1.59 h) | 1.69e-13 (23.60 h) | 1.71e-13 (10.67 h) | at most 1e-11 |
| max `rel_u` | 1.88e-12 (7.02 h) | 1.20e-12 (22.96 h) | 9.87e-13 (5.42 h) | at most 1e-11 |
| max `rel_P` | 2.67e-14 (7.45 h) | 1.87e-14 (19.10 h) | 1.72e-14 (1.57 h) | at most 1e-11 |
| max `rel_T` | 5.71e-13 (8.33 h) | 3.78e-13 (17.74 h) | 3.25e-13 (5.80 h) | at most 1e-11 |
| max `rel_w` | 2.55e-8 (0.26 h) | 1.38e-8 (0.35 h) | 7.05e-9 (0.35 h) | at most 1e-7 |
| `MF(30.25)` change, 6 to 9 h | 2.32e-12 | 4.16e-12 | 2.28e-12 | at most 1e-7 |
| `MF(30.25)` change, 18 to 24 h | | 1.02e-12 | 3.04e-12 | at most 1e-7 |

Every value lies more than a factor 2 below its edge (the closest, `rel_w` of S02F1, by a
factor 3.9), so the factor 2 rule does not apply. The 18 to 24 h `MF(30.25)` item of S02F3,
which lies after its `t_c`, is met as well. Step counts are equal in each pair (98 768,
131 070, 131 079). The data are consistent with H0 on all three cases: the two solvers give
the same solution to round-off at dz 250 m over 9 h and at dz 500 m over 24 h. As H0
predicted, `rel_w` peaks during the ramp and falls with the forcing amplitude on the
dz 500 m column (1.38e-8 at `U h` = 2 km m/s, 7.05e-9 at 5.7), and the absolute `w`
difference after 1 h stays between 2.4e-11 and 5.7e-11 m/s at dz 500 m and between 4.2e-11
and 9.6e-11 m/s at dz 250 m, with no trend.

Reported without a band. Window means of `MF(30.25)` [µPa], off member and twin equal to
the digits shown: S02F1 0.222736 (6 to 9 h, hourly range 0.187 to 0.268), S02F2 0.0586574
(6 to 9 h) and 0.153513 (18 to 24 h, hourly range 0.133 to 0.173), S02F3 0.499126 (6 to 9 h)
and 1.57874 (18 to 24 h, hourly range 1.341 to 1.822). The changes of `MF(1.25)`,
`MF(10.25)`, `MF(50.25)`, `R10/30` and `R50/30` over both windows lie between 3.4e-14 and
8.7e-12 in magnitude. H3: the 9 to 24 h power-law slopes of `rel` in S02F2 are `rho` 0.52,
`u` 0.087 and `T` 0.090, below 1, and `rel_u` at 24 h is 1.05e-12, against about 4e-12 under
H3; H3 is not supported, as expected. The `rho` slope exceeds the 0.43 envelope of the
Section 1 pilot while `rel_rho` stays between 1.5e-13 and 1.7e-13. In S02F3, after its
`t_c`, the slopes are 1.37 (`rho`), 1.54 (`u`) and 1.43 (`T`), and `rel_u` reaches 5.35e-12
at 23.97 h; `rel_max` stays below 2.9e-10 after `t_c`, so no growth is fitted (the fitted
rate of `ln rel_max` over 10.69 to 24 h is −0.05 per hour, set by the decay of `rel_w`).
The onsets of the twins equal those of the off members: none in S02F1 and S02F2; in S02F3
the 0 to 5 km convective onset at 10.6929 h (1.0 km) and the 45 to 55 km one at 15.7464 h
(54.5 km). Snapshot times differ from step 300 (S02F1 from step 400) by at most 1.8e-10 s,
a round-off difference of the CFL `dt`.

**D-2, row 1 met in all five pairs.** Up to `t_c` the largest `rel` of `rho`, `u`, `P`, `T`
over the five pairs is 8.95e-13 (`u`, S02F4, 3.61 h), against 1e-11; the largest `rel_w`
is 2.33e-9 (S02F4, 0.35 h), against 1e-8. The `t_c` event of every twin is the off member's
event, a convective onset in 0 to 5 km at 1.0 km, at the same snapshot (time differences
1e-15 to 2e-14 h, against one scored snapshot of about 132 s).

| pair | snapshots to `t_c` | `t_c` [h] | max thermo `rel` to `t_c` | max `rel_w` to `t_c` | `rel_max` first above 1e-6, 1e-3 [h] | e-folding time after `t_c` [h] |
|---|---|---|---|---|---|---|
| S02F4 / PR1-U20-h1.0 | 251 | 4.5769 | 8.95e-13 | 2.33e-9 | 19.83, not reached | 1.21 |
| S02F5 / S15F1 | 175 | 3.1850 | 7.47e-13 | 1.16e-9 | 11.13, 12.84 | 0.57 |
| S02F6 / PR1-U57-h1.0 | 147 | 2.6721 | 5.98e-13 | 8.77e-10 | 9.08, 10.74 | 0.45 |
| S02F7 / PR1-U20-h3.0 | 145 | 2.6355 | 6.77e-13 | 8.33e-10 | 8.99, 10.36 | 0.45 |
| S02F8 / PR1-U40-h3.0 | 105 | 1.9026 | 4.71e-13 | 4.36e-10 | 6.40, 6.98 | 0.33 |

The data are consistent with H0 up to the first overturning in every breaking pair: the
two solvers give the same flow to round-off up to `t_c`, for `h` = 1 km and `h` = 3 km.

Reported without a band (descriptive). The e-folding times are fitted to `ln rel_max` from
`t_c` to the first time `rel_max` exceeds 1e-3 (or the last shared snapshot), and bracket
the 0.61 h of S01F3 (0.33 to 1.21 h). `rel_max` first exceeds 1e-9 during the ramp in
S02F4 and S02F5 (0.22 h and 0.31 h, `rel_w` only) and at 5.10 to 6.99 h in the others.
Every later onset of every twin is the off member's, same band, kind and height, within
3.6e-6 h (13 ms), except one: the first wind reversal in 20 to 45 km of S02F5 at 14.1083 h
against 14.1369 h (44.25 km in both), 103 s earlier, inside one scored snapshot and after
`rel_max` passed 1e-3. The first free-column onsets `t_ov` and `dz_L` (from the `zoverturn`
reports) are equal in each `h` = 1 km pair: 15.3367 h and −6.25 km (S02F4), 13.0174 h and
+6.75 km (S02F5), 10.1512 h and +9.75 km (S02F6). The `zh3` crossings of the `h` = 3 km
twins equal those of `output/analysis/s15/z3_face_U*_h3.json` to about 1e-12 h: ground
`|T/T_0 - 1| > 0.5` at 5.2699 h (U20) and 3.4018 h (U40), top `|w| > 5` m/s at 8.7439 h and
5.0785 h, free-column instability above 1 % at 8.5172 h and 6.1274 h. Pre-collapse window
means of `MF(30.25)` [µPa], off member and twin, with the relative change and the hourly
range of the off member: S02F4 35.1961 and 35.1961 (12.34 to 15.34 h, −1.96e-11, 18.2 to
52.4); S02F5 50.9989 and 51.0002 (10.02 to 13.02 h, 2.55e-5, 27.5 to 86.8); S02F6 79.9813
and 79.9813 (7.15 to 10.15 h, 2.13e-8, 50.0 to 127.3). Each change is at least 1e4 times
smaller than the relative hourly range. The ends of the runs: S02F6, S02F7 and S02F8 end at
11.91, 11.49 and 7.01 h against 11.99, 11.45 and 7.02 h, and S02F5 ends at 15.35 h against
23.81 h for S15F1, after `rel_max` passed 1e-3 at 12.84 h. S02F4 ran 24 h with the off
member's step count and stays below 1e-6 to 19.8 h, although it overturns near the ground
at 4.58 h. The growth after `t_c`, and the different end of S02F5, are consistent with H2,
round-off amplification by the overturning flow, which any perturbation of the `:lu` run
would also show; they are not by themselves evidence on the solver. The snapshot times
first differ at step 300 in every pair, as in Section 1.

**D-3, cost (no band).** At concurrency 8 on an otherwise idle machine, the twins ran at
0.0178 to 0.0182 s per step on the dz 500 m column and 0.0352 s on dz 250 m, against the
estimates 0.020 and 0.039. The cost per step is the same in breaking and non-breaking runs:
the breaking twins refactorized 43 to 61 times per 1000 steps at no visible cost, where their
`:lu` off members ran at 0.143 to 0.165 s per step against 0.103 for the non-breaking ones.

| pair | twin / off member, `real` | twin / off member, s per step |
|---|---|---|
| S02F1 / S14S2 (dz 250 m, 9 h) | 0.190 | 0.190 |
| S02F2 / PR1-U20-h0.1 | 0.174 | 0.174 |
| S02F3 / PR1-U57-h0.1 | 0.172 | 0.172 |
| S02F4 / PR1-U20-h1.0 | 0.172 | 0.172 |
| S02F5 / S15F1 | 0.072 | 0.121 |
| S02F6 / PR1-U57-h1.0 | 0.106 | 0.110 |
| S02F7 / PR1-U20-h3.0 | 0.130 | 0.127 |
| S02F8 / PR1-U40-h3.0 | 0.124 | 0.124 |

The ratios mix three changes and are descriptive of the production gain only: the solver,
the pack refactor `0f6ff62` (which also speeds up `:lu`), and concurrency 8 on an idle
machine against 11 or 12 for the off members. The expected `real` ratio of about 0.2 for the
24 h `h` = 0.1 km cases is met at 0.17. The eight twins used 4.4 h of CPU in all, against
31.3 h of `real` for their off members. The S02F9 item is not scored (note of 2026-10-20).

**Decision rule.** D-0 met, D-1 row 1 and the second statement met in all three pairs, and
D-2 row 1 met in all five pairs: the first case of the rule applies. New production sets
pass `diffsolver=fft` together with `bcface=1` and `pbar=1`. The `SimConfig` default stays
`:lu`, no existing set changes, and `check` keeps judging the `:lu` default. Limitation, as
recorded in the note of 2026-10-20: on the dz 250 m grid this rests on S02F1 alone, which
does not overturn. For dz 250 m runs that overturn, no full-run pair compares the two
solvers, and a value quoted from such a run carries this as a solver row of the error budget
until a registered pair covers it. Quoting: a value from a `:fft_tridiag` record names the
solver along with the forcing path, `dt` or CFL, `dz`, window and run ID; the descriptive
D-2 differences after `t_c` (the S02F5 end time and its 2.55e-5 window-mean change) are the
solver row of the error budget for runs that overturn, kept separate from the sampling
scatter.

**Note added 2026-10-20, after scoring: the solver's share of the D-3 ratios (descriptive).**
The D-3 ratios of 0.11 to 0.19 (5 to 9 times) are not the gain of the solver alone. They
split into three factors with the standalone medians already registered (Section 1 item c
for the non-diffusion part and the factorization, Section 2 item b for `diffusion_step!`
after `0f6ff62`); no new measurement is made. On the dz 500 m column without breaking:
(i) the pack refactor takes the `:lu` step from 0.013 + 0.0710 = 0.084 s to
0.013 + 0.0250 = 0.038 s, a factor 2.2; (ii) the solver takes it from 0.038 s to
0.013 + 0.0040 = 0.017 s, a factor 2.2; (iii) the batch conditions add a factor 1.17, the
off members having run at 0.103 / 0.084 = 1.23 times the old standalone rate at
concurrency 11 or 12 and the twins at 0.0178 / 0.017 = 1.05 times the new one at
concurrency 8. The product, 2.2 x 2.2 x 1.17 = 5.7, matches the measured 1 / 0.174 = 5.75.
At dz 250 m without breaking the solver alone gives (0.026 + 0.050 to 0.059) / 0.034, 2.2
to 2.5 times. In runs that overturn, `:lu` adds one factorization of 0.701 s (dz 500 m)
at the rate the twins counted, 42.7 to 60.6 per 1000 steps, so the solver alone gives
(0.038 + 0.030 to 0.042) / 0.017, 4.0 to 4.7 times; the same model on the old code with
the factor 1.23 gives 0.140 to 0.156 s per step, against the 0.143 to 0.165 s measured for
the off members, whose own rates are not counted here. The production gain of
`diffsolver=fft` on the current code is therefore about 2 times before breaking and about
4 to 5 times after it at dz 500 m. The other factor of about 2.6 belongs to the pack
refactor, which `:lu` runs keep, and to the quieter machine. A dz 250 m run that overturns
would gain more, with a 2.068 s factorization, but no pair measures it (open item (e)).
