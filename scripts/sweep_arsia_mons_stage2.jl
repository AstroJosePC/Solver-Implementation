# =============================================================================
# Arsia Mons Stage 2 mountain-wave sweep (real reversing shear, critical levels)
#
# Stage 2 keeps the Stage 1 orographic forcing but replaces the uniform steering
# wind with the REAL reversing zonal wind u(z) from the MY31 MCD column, and
# switches the horizontal boundary to PERIODIC. The aim is to test whether a
# stationary mountain wave launched into this shear deposits momentum near the
# 45-53 km AMEC cloud band.
#
# What changed from Stage 1 (and why):
#   1. Background wind: u_init is now a Function closure interpolating the CSV
#      u_zonal column onto grid.z_full (the real reversing jet), instead of a
#      uniform scalar. The sheared background is a discrete steady state of the
#      solver (x-uniform => zero x-sweep tendency; w=0, P=P_bar => zero z-sweep
#      flux/source), and the implicit diffusion splits u = u_bar(z) + u' and
#      diffuses only u', so the jet is preserved.
#   2. Horizontal BC: :periodic instead of :outflow. Stage 1's Dirichlet inflow
#      edge produced a stationary perturbation that amplified as rho^-1/2 with
#      height (a spurious feature ~900 km upwind). Periodic-x removes that
#      special edge entirely.
#   3. Forcing wind: U_bar in the orographic BC w = U_bar * dh/dx stays the
#      Stage 1 sweep value (eastward +20/40/57). The real surface wind is
#      U(0)=0, so the forcing amplitude is a deliberate EFFECTIVE launching wind
#      (decoupled from the background u(z) the wave then propagates through).
#      Keeping it identical to Stage 1 isolates the effect of the shear alone.
#   4. Domain width LX is exposed as a sensitivity axis. With periodic-x the
#      wakes eventually wrap; a wider domain delays that. The main U x h grid
#      runs at a single baseline width; a separate `width` mode sweeps LX at one
#      representative (U, h) to check convergence.
#
# Critical levels (zero-crossings of u_zonal; critical for a stationary c=0 wave):
#   ~4.5 km, ~53.3 km (top of the AMEC band, below the 60 km sponge), ~77.7 km
#   (inside the sponge, damped). The 53.3 km level is the science target.
#
# Usage (run on the Mac Studio M2; threads via -t auto):
#   julia -t auto --project=. scripts/sweep_arsia_mons_stage2.jl smoke
#       -> single short run (U=40, h=1 km, t_end=1 h, baseline width) to confirm
#          the sheared + periodic 90 km config is stable; reports wall-clock.
#   julia -t auto --project=. scripts/sweep_arsia_mons_stage2.jl
#       -> main 3x3 sweep (U_bar x h_max) at the baseline width.
#   julia -t auto --project=. scripts/sweep_arsia_mons_stage2.jl width
#       -> width sensitivity {2000,3000,4000} km at the representative (U=40,h=1).
#   julia --project=. scripts/sweep_arsia_mons_stage2.jl linear [u=40 h=0.1 hours=12
#                                                                 sigma=100]
#       -> LINEAR CONTROL corner: uniform wind (u_init = U_bar) with periodic-x at
#          h = 0.1 km, into output/arsia_mons_linear<hours>h/. The wave stays below
#          the overturning criterion everywhere under the sponge, so Eliassen-Palm
#          requires the momentum flux to be height-independent and any measured
#          drift is the solver's flux error bar. LX defaults to 2000 km (the Stage 1
#          width) rather than 4000 km, so the result transfers to the Stage 1 runs.
#          `sigma=NN` (km) reweights the ridge spectrum without moving lambda_z or
#          any per-mode arrival time, which is what makes a sigma pair difference
#          out the numerical part of the measured flux drift.
#          `k0=X` (1/s), `lz=NNN` (km) and `spbase=NN` (km) vary the sponge
#          strength, the model top and the sponge base (defaults 0.05, 90, 60);
#          non-default values are appended to the output directory name.
#          `mode=N` replaces the ridge by the single Fourier mode h sin(2 pi N x/LX)
#          of amplitude `h` (the cavity closure test); the directory is tagged
#          _mode<N> in place of the sigma tag.
#   `ramp=T` (hours, any mode) switches the forcing on as tanh(t/T) instead of
#   impulsively: `:linear_ramp` for the single mode, MountainParams.t_ramp for the
#   ridge. The directory is tagged _ramp<T>h.
#   `pbar=1` selects SimConfig.pbar_mode = :pressure (the fixed background
#   pressure, docs/reflection_convergence_tests.md Section 9) and tags the output
#   directory _pbarP. The script default stays pbar=0 (:total) so the commands
#   recorded for the earlier records reproduce them; all new work passes pbar=1.
#   In the main grid, `u=` and `h=` (km) accept comma lists to run a subset,
#   e.g. `pbar=1 u=40 h=1,3`, and `hours=H` replaces the 3 h record (directory
#   tagged _<H>h when H != 3). `cfl=X` sets the CFL number (default 0.4, tag
#   _cfl<X> otherwise). The sponge, top, pbar and ramp tags apply to every
#   mode, and every mode writes under MAGNUSP_OUTPUT when it is set.
#   Section 14 options (docs/reflection_convergence_tests.md), each a no-op at its
#   default: `dx=N` and `dz=N` (m, defaults 4000 and 500, tags _dx<N>, _dz<N>) set
#   the grid; `outint=N` (steps, default 100, no tag) sets the snapshot cadence, so
#   a smaller CFL can keep about 66 s between snapshots; `ublend=Z` (km, main grid
#   only, tag _ub<Z>) replaces the MCD wind below Z by the launch-consistent blend
#   u = u_L (1 - S) + u_MCD S, S = (1 - cos(pi z/Z))/2, u_L = sign(u_MCD(Z)) |U_bar|,
#   and refuses a forcing wind `u=` whose sign differs from u_L. `bcface=1` selects
#   SimConfig.bottom_forcing = :face (forced w at t^{n+1/2} imposed on the ground face
#   of the z-sweep; Section 14 addendum, tag _bcface). Addendum 2 options:
#   `spsplit=1` (SimConfig.sponge_step = :split, exact relaxation after the LW step,
#   tag _spsplit), `zx=1` (sweep_order = :zx, z-sweep before x-sweep, tag _zx) and
#   `probe=N` (face_probe = N: ground-face w and p' of the z-predictor every N steps
#   to face_probe.nc; leaves the solution unchanged, no tag). Tag order is
#   ..._ramp<T>h_dx<N>_dz<N>_ub<Z>_bcface_spsplit_zx_cfl<X>.
#   Section 16 options (docs/reflection_convergence_tests.md): `alt=1` (sweep_order
#   = :alternate, tag _alt) or `strang=1` (:strang, x(dt/2) z(dt) x(dt/2), tag
#   _strang) in place of `zx=1`, and `opstrang=1` (operator_split = :strang: split sponge,
#   diffusion and cooling over dt/2 on both sides of the step, tag _opS, after the
#   sweep tag). `final=1` also writes the snapshot at t_end (no tag).
#   Optional `lx=NNNN` (km) overrides the baseline width of any mode, e.g.
#   `julia ... sweep_arsia_mons_stage2.jl lx=3000`. The solver itself is serial (no
#   @threads in src/), so -t auto only helps the multi-combination sweep modes.
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using Base.Threads
using MAGNUSP
using Printf
using DelimitedFiles

# --- Run mode (ARGS) ---------------------------------------------------------
const SMOKE = "smoke" in ARGS
const WIDTH_SWEEP = "width" in ARGS
const LINEAR = "linear" in ARGS   # uniform-wind linear control (see mode block)

"""Parse an optional `key=value` numeric argument; fall back to `default`."""
function getarg(key::String, default::Float64)
    i = findfirst(a -> startswith(a, key * "="), ARGS)
    return i === nothing ? default : parse(Float64, split(ARGS[i], "=")[2])
end

"""Parse an optional comma-separated `key=a,b,...` list; fall back to `default`."""
function getlist(key::String, default::Vector{Float64})
    i = findfirst(a -> startswith(a, key * "="), ARGS)
    return i === nothing ? default :
           parse.(Float64, split(split(ARGS[i], "=")[2], ","))
end

"""Domain width [m] from an optional `lx=NNNN` (km) argument."""
parse_lx_arg(default_km::Float64) = getarg("lx", default_km) * 1e3

# --- Background atmosphere (real MY31 MCD, 0-120 km) -------------------------
const ATM_FILE = abspath(joinpath(@__DIR__, "..", "external", "mcd_profiles",
                                   "mars_arsia_mons_my31.csv"))

isfile(ATM_FILE) || error("Background atmosphere not found: $ATM_FILE\n" *
    "Generate it with: python scripts/prep_mars_arsia_mons_atmosphere.py " *
    "--output external/mcd_profiles/mars_arsia_mons_my31.csv")

# --- Real reversing zonal wind u(z) -----------------------------------------
# Closure interpolating the CSV u_zonal column (col 6). Clamps to the endpoints
# (below ground -> u(0); above the table top -> u(120 km)), so ghost cells are
# well defined. This is the BACKGROUND the wave propagates through, NOT the
# orographic forcing wind.
function make_u_zonal_closure(atm_file::String)
    data, _ = readdlm(atm_file, ',', Float64, header=true)
    z_data = data[:, 1]
    u_data = data[:, 6]
    return function (z::Float64)
        if z <= z_data[1]
            return u_data[1]
        elseif z >= z_data[end]
            return u_data[end]
        else
            idx = searchsortedlast(z_data, z)
            t = (z - z_data[idx]) / (z_data[idx + 1] - z_data[idx])
            return u_data[idx] * (1.0 - t) + u_data[idx + 1] * t
        end
    end
end

const U_ZONAL = make_u_zonal_closure(ATM_FILE)

# Launch-consistent background (Section 14, `ublend=Z`): below the blend height the
# MCD wind is replaced by u = u_L (1 - S) + u_MCD S with S = (1 - cos(pi z/Z))/2, so
# the ground wind is u_L = sign(u_MCD(Z)) |U_bar| and the imposed w = U_bar dh/dx is
# the kinematic condition of the flow the wave is launched into. Above Z the MCD
# wind is unchanged. The forcing wind must have the sign of u_L, otherwise the ridge
# would be forced against the wind it sits in.
function make_blended_wind(u_mcd, zb_m::Float64, U_bar::Float64)
    uL = sign(u_mcd(zb_m)) * abs(U_bar)
    sign(U_bar) == sign(uL) || error(@sprintf(
        "ublend: forcing wind u=%.1f and blended ground wind u_L=%.1f differ in sign; pass u=%.0f",
        U_bar, uL, uL))
    return function (z::Float64)
        z >= zb_m && return u_mcd(z)
        S = 0.5 * (1.0 - cos(pi * max(z, 0.0) / zb_m))
        return uL * (1.0 - S) + u_mcd(z) * S
    end
end

# --- Domain / numerics (vertical column identical to Stage 1) ----------------
# Vertical extent [m]. Settable with `lz=NNN` (km) for the sponge thickness and
# model-top tests of docs/reflection_convergence_tests.md; the MCD column reaches
# 120 km, so lz <= 120.
const LZ = getarg("lz", 90.0) * 1e3
LZ <= 120e3 || error("lz = $(LZ / 1e3) km exceeds the 120 km MCD column")
const DX = getarg("dx", 4000.0)   # [m]
const DZ = getarg("dz", 500.0)    # [m] -> nz = LZ / DZ (190 for the 95 km top)
DX > 0 && DZ > 0 || error("dx and dz must be positive")
# Gaussian hill half-width [m] (Arsia-scale). Settable with `sigma=NN` (km) so a
# twin of the linear control can be run at a narrower ridge. The half-width sets
# only the flux weighting across horizontal wavenumbers: the per-mode arrival time
# I(z)/(k U^2) depends on k and not on sigma, and the steady wave amplitude
# u' = N h (rho_0 N_0 / rho N)^(1/2) carries no k either, so the linearity margin
# and the vertical structure are unchanged. See docs/linear_flux_constancy.md.
const SIGMA = getarg("sigma", 100.0) * 1e3
nz = Int(round(LZ / DZ))

# Sponge base 60 km by default (60-90 km layer), 7 km above the 53 km AMEC
# level. `spbase=NN` (km) moves the base and `k0=X` (1/s) sets the peak damping
# rate of the sin^2 profile; the layer always extends from the base to LZ.
const SPONGE_BASE = getarg("spbase", 60.0) * 1e3
0 < SPONGE_BASE < LZ || error("spbase must lie inside (0, lz)")
const SPONGE_FRAC = (LZ - SPONGE_BASE) / LZ
const K0_SPONGE = getarg("k0", 0.05)
# Single Fourier mode in place of the Gaussian ridge (LINEAR mode only). With
# `mode=N` the terrain is h(x) = h sin(2 pi N x / LX), so the kinematic condition
# is w = U k h cos(k x) with k = 2 pi N / LX, one stationary horizontal wavenumber
# with no ridge spectrum and no fill-in. 0 keeps the ridge.
const SINGLE_MODE = Int(getarg("mode", 0.0))
# Background pressure in the well-balanced vertical flux (SimConfig.pbar_mode).
# `pbar=1` subtracts the background pressure alone; the default 0 keeps the
# historical (gamma-1)*Q_bar[4], which includes the mean-wind kinetic energy.
const PBAR_MODE = getarg("pbar", 0.0) > 0 ? :pressure : :total
# Forcing ramp time [h]; 0 keeps the impulsive start.
const RAMP_H = getarg("ramp", 0.0)
RAMP_H >= 0 || error("ramp must be >= 0")
# CFL number of the adaptive time step (`cfl=X`); 0.4 is the default of every record.
const CFL_NUM = getarg("cfl", 0.4)
# Launch-consistent blend height [km] (`ublend=Z`, main grid only); 0 keeps the MCD wind.
const UBLEND_KM = getarg("ublend", 0.0)
# Path of the imposed ground w (`bcface=1`: :face); 0 keeps the ghost-cell path.
const BOTTOM_FORCING = getarg("bcface", 0.0) > 0 ? :face : :ghost
# Mechanism switches of the Section 14 addendum 2; 0 keeps the default path.
const SPONGE_STEP = getarg("spsplit", 0.0) > 0 ? :split : :inline
const SWEEP_ORDER = getarg("strang", 0.0) > 0 ? :strang :
                    getarg("alt", 0.0) > 0 ? :alternate :
                    getarg("zx", 0.0) > 0 ? :zx : :xz
# Symmetric split operators (Section 16); 0 keeps the Lie sequence.
const OPERATOR_SPLIT = getarg("opstrang", 0.0) > 0 ? :strang : :lie
# `final=1` also writes the state at t_end (SimConfig.final_snapshot); no tag.
const FINAL_SNAPSHOT = getarg("final", 0.0) > 0
const FACE_PROBE = Int(getarg("probe", 0.0))
UBLEND_KM >= 0 || error("ublend must be >= 0")
UBLEND_KM > 0 && (LINEAR || SMOKE || WIDTH_SWEEP) &&
    error("ublend applies to the main grid only (the linear control has a uniform wind)")
# Output root. MAGNUSP_OUTPUT points it at a local disk, which keeps the snapshot
# writes off a network mount; the default is the repository's output/.
const OUT_ROOT = get(ENV, "MAGNUSP_OUTPUT", joinpath(@__DIR__, "..", "output"))
# Directory tags for non-default layer, top, background pressure and ramp. Each is
# appended only when it departs from the default, so the directories of the earlier
# records keep their names and a variant never overwrites them.
numlabel(v) = v == round(v) ? string(Int(round(v))) : string(v)
const SPGLABEL = (K0_SPONGE == 0.05 ? "" : "_K0_$(K0_SPONGE)") *
                 (LZ == 90e3 ? "" : "_LZ$(Int(round(LZ / 1e3)))") *
                 (SPONGE_BASE == 60e3 ? "" : "_spb$(Int(round(SPONGE_BASE / 1e3)))") *
                 (PBAR_MODE === :total ? "" : "_pbarP") *
                 (RAMP_H == 0 ? "" : "_ramp$(numlabel(RAMP_H))h") *
                 (DX == 4000.0 ? "" : "_dx$(numlabel(DX))") *
                 (DZ == 500.0 ? "" : "_dz$(numlabel(DZ))") *
                 (UBLEND_KM == 0 ? "" : "_ub$(numlabel(UBLEND_KM))") *
                 (BOTTOM_FORCING === :ghost ? "" : "_bcface") *
                 (SPONGE_STEP === :inline ? "" : "_spsplit") *
                 (SWEEP_ORDER === :xz ? "" : "_" * Dict(:zx => "zx", :alternate => "alt",
                                                        :strang => "strang")[SWEEP_ORDER]) *
                 (OPERATOR_SPLIT === :lie ? "" : "_opS") *
                 (CFL_NUM == 0.4 ? "" : "_cfl$(CFL_NUM)")

# Representative N for the nonlinearity diagnostic Nh/U (lower-atmosphere mean).
const N_REF = 0.018

const OUTPUT_INTERVAL = Int(getarg("outint", 100.0))   # snapshot cadence [steps]
OUTPUT_INTERVAL > 0 || error("outint must be a positive number of steps")
const DIAG_INTERVAL = 500     # diagnostics / blowup-check cadence

# Zero-crossings of u_zonal: critical levels for a stationary mountain wave.
# (Recorded as metadata; the analyzer re-detects these from the CSV.)
const CRITICAL_LEVELS_KM = [4.5, 53.3, 77.7]

# --- Sweep axes --------------------------------------------------------------
if SMOKE
    U_bar_vals = [40.0]
    h_max_vals = [1000.0]
    LX_vals = [parse_lx_arg(4000.0)]
    t_end = 3600.0                 # 1 hour
    sweep_base = abspath(joinpath(OUT_ROOT, "arsia_mons_stage2_smoke$(SPGLABEL)"))
elseif LINEAR
    # Linear control corner: UNIFORM wind (u_init = U_bar, as in Stage 1) with the
    # Stage 2 periodic-x boundary, at h = 0.1 km. At that height the steady linear
    # solution for this ridge gives u_hat/U = 0.45 at 55 km, so the stationary-wave
    # overturning criterion u_hat >= U is not met anywhere below the sponge, and
    # molecular viscosity from the MCD nu column removes only 0.03% of the flux by
    # 55 km. Eliassen-Palm then requires d<rho u'w'>_x/dz = 0, so any measured drift
    # is the solver's own flux error bar. LX defaults to the Stage 1 width so the
    # drift transfers directly to the Stage 1 flux quotes.
    U_bar_vals = [getarg("u", 40.0)]
    h_max_vals = [getarg("h", 0.1) * 1e3]
    lin_hours = getarg("hours", 12.0)
    t_end = lin_hours * 3600.0
    LX_vals = [parse_lx_arg(2000.0)]
    # Horizon- and width-named output directory, so successive horizons do not
    # overwrite each other or their recorded wall-clock. The width is in the name
    # because the diagnostic <rho u'w'>_x is a domain x-mean and so scales as 1/LX:
    # the same wave gives -24.70 uPa at 2000 km and -12.51 uPa at 4000 km.
    hlabel = lin_hours == round(lin_hours) ? string(Int(round(lin_hours))) :
             string(lin_hours)
    # The half-width is appended only when it departs from the 100 km default, so
    # the directory names of the runs already on disk do not move.
    siglabel = SIGMA == 100e3 ? "" : "_sig$(Int(round(SIGMA / 1e3)))"
    modelabel = SINGLE_MODE > 0 ? "_mode$(SINGLE_MODE)" : siglabel
    sweep_base = abspath(joinpath(OUT_ROOT,
        "arsia_mons_linear$(hlabel)h_LX$(Int(round(LX_vals[1] / 1e3)))$(modelabel)$(SPGLABEL)"))
elseif WIDTH_SWEEP
    # Width sensitivity at one representative (U, h); LX is the swept axis.
    U_bar_vals = [40.0]
    h_max_vals = [1000.0]
    LX_vals = [2000e3, 3000e3, 4000e3]
    t_end = 10800.0                # 3 hours
    sweep_base = abspath(joinpath(OUT_ROOT, "arsia_mons_stage2_width$(SPGLABEL)"))
else
    # Main grid: mirror Stage 1 (U_bar x h_max) at a single baseline width.
    # `u=` and `h=` (km) lists select a subset of the grid.
    U_bar_vals = getlist("u", [20.0, 40.0, 57.0])       # m/s eastward effective launching wind
    h_max_vals = getlist("h", [0.1, 1.0, 3.0]) .* 1e3   # m (linear -> nonlinear onset)
    LX_vals = [parse_lx_arg(4000.0)]       # baseline width (override with lx=NNNN)
    main_hours = getarg("hours", 3.0)   # record length [h]
    t_end = main_hours * 3600.0
    hourlabel = main_hours == 3.0 ? "" : "_$(numlabel(main_hours))h"
    sweep_base = abspath(joinpath(OUT_ROOT, "arsia_mons_stage2$(hourlabel)$(SPGLABEL)"))
end

# Tag run dirs with the width only when LX is the swept axis (keeps the main
# grid's names matching the Stage 1 convention run_U<u>_h<h>).
const TAG_LX = length(LX_vals) > 1

run_dir_name(U, h, lx) = TAG_LX ?
    "run_LX$(Int(round(lx / 1e3)))_U$(U)_h$(h / 1e3)" :
    "run_U$(U)_h$(h / 1e3)"

struct RunResult
    U_bar::Float64
    h_max::Float64
    LX::Float64
    Nh_over_U::Float64
    duration_s::Float64
    final_t::Float64
    steps::Int
    converged::Bool
    outdir::String
end

results = Vector{RunResult}()
results_lock = ReentrantLock()

mode_label = SMOKE ? "SMOKE" : (LINEAR ? "LINEAR CONTROL" :
              (WIDTH_SWEEP ? "WIDTH-SWEEP" : "MAIN U x h"))
@printf("Arsia Mons Stage 2 sweep (%s)\n", mode_label)
@printf("  atm_file : %s\n", ATM_FILE)
if LINEAR
    @printf("  wind     : UNIFORM u_init = U_bar = %.1f m/s (linear control; no critical levels)\n",
            U_bar_vals[1])
else
    @printf("  wind     : REAL reversing u_zonal(z) (closure); U(0)=%.2f m/s; crit. levels ~%s km\n",
            U_ZONAL(0.0), join(string.(CRITICAL_LEVELS_KM), "/"))
    if UBLEND_KM > 0
        # Resolve every blend now, so a sign mismatch stops the run before any output.
        for U in U_bar_vals
            ub = make_blended_wind(U_ZONAL, UBLEND_KM * 1e3, U)
            @printf("  ublend   : %.0f km | u(0)=%.1f m/s -> u_MCD(%.0f km)=%.2f m/s (forcing u=%.0f)\n",
                    UBLEND_KM, ub(0.0), UBLEND_KM, U_ZONAL(UBLEND_KM * 1e3), U)
        end
    end
end
@printf("  bc_x     : periodic | bottom: reflective | top: eq19 outflow + sponge\n")
@printf("  vertical : nz=%d (dz=%.0f) LZ=%.0f km | sponge %.0f-%.0f km (K0=%.4f)\n",
        nz, DZ, LZ / 1e3, SPONGE_BASE / 1e3, LZ / 1e3, K0_SPONGE)
@printf("  grid     : dx=%.0f m | snapshot every %d steps\n", DX, OUTPUT_INTERVAL)
@printf("  widths   : %s km | sigma=%.0f km | t_end=%.0f s | threads=%d\n",
        join(string.(Int.(round.(LX_vals ./ 1e3))), ","), SIGMA / 1e3, t_end,
        Threads.nthreads())
@printf("  pbar_mode: %s | forcing ramp: %s | ground w: %s | CFL=%.2f\n",
        PBAR_MODE, RAMP_H > 0 ? "tanh(t/$(numlabel(RAMP_H)) h)" : "none (impulsive)",
        BOTTOM_FORCING, CFL_NUM)
@printf("  sponge   : %s | sweep order: %s | split ops: %s | face probe: %s\n  output   : %s\n",
        SPONGE_STEP, SWEEP_ORDER, OPERATOR_SPLIT, FACE_PROBE > 0 ? "every $(FACE_PROBE) steps" : "off",
        sweep_base)
"dryrun" in ARGS && exit(0)   # print the resolved configuration only

mkpath(sweep_base)

combinations_arr = collect(Iterators.product(U_bar_vals, h_max_vals, LX_vals))

# Pre-create run directories sequentially (safe with external-drive metadata).
for (U_bar, h_max, lx) in combinations_arr
    mkpath(normpath(joinpath(sweep_base, run_dir_name(U_bar, h_max, lx))))
end

Threads.@threads for params in combinations_arr
    U_bar, h_max, lx = params
    outdir = normpath(joinpath(sweep_base, run_dir_name(U_bar, h_max, lx)))
    nx_local = Int(round(lx / DX))
    Nh_over_U = N_REF * h_max / U_bar

    config = mars_config(
        nx = nx_local,
        nz = nz,
        dx = DX,
        dz = DZ,
        CFL = CFL_NUM,
        t_end = t_end,
        output_interval = OUTPUT_INTERVAL,
        enable_viscosity = true,
        enable_cooling = false,
        sponge_frac = SPONGE_FRAC,
        K0_sponge = K0_SPONGE,
        pbar_mode = PBAR_MODE,
        bottom_forcing = BOTTOM_FORCING,
        sponge_step = SPONGE_STEP,
        sweep_order = SWEEP_ORDER,
        operator_split = OPERATOR_SPLIT,
        final_snapshot = FINAL_SNAPSHOT,
        face_probe = FACE_PROBE,
    )

    mountain_params = MAGNUSP.MountainParams(
        h_max = h_max,
        sigma = SIGMA,
        x0 = 0.5 * lx,
        U_bar = U_bar,                 # effective launching wind (NOT u_init)
        t_ramp = RAMP_H * 3600.0,      # 0 = impulsive (bit-identical to earlier runs)
    )
    single = LINEAR && SINGLE_MODE > 0
    k_mode = 2pi * SINGLE_MODE / lx
    # With a ramp, t0 is the tanh time scale of :linear_ramp; unused otherwise.
    forcing_params = single ?
        MAGNUSP.WaveParams(amplitude = U_bar * k_mode * h_max, omega = 0.0,
                           kx = k_mode, t0 = RAMP_H > 0 ? RAMP_H * 3600.0 : 1200.0) :
        mountain_params
    forcing_type = !single ? :mountain : (RAMP_H > 0 ? :linear_ramp : :linear)

    local state
    duration = @elapsed begin
        state = run_simulation(config;
            forcing_type = forcing_type,
            wave_params = forcing_params,
            # uniform control wind, real shear, or the launch-consistent blend
            u_init = LINEAR ? U_bar :
                     (UBLEND_KM > 0 ? make_blended_wind(U_ZONAL, UBLEND_KM * 1e3, U_bar) : U_ZONAL),
            bc_x_type = :periodic,         # removes the Stage 1 inflow-edge artifact
            atm_file = ATM_FILE,
            outdir = outdir,
            diag_interval = DIAG_INTERVAL,
        )
    end

    converged = state.t >= t_end * 0.999
    res = RunResult(U_bar, h_max, lx, Nh_over_U, duration, state.t, state.step,
                    converged, outdir)

    lock(results_lock) do
        push!(results, res)
        @printf("Done U=%.0f h=%.2fkm LX=%dkm | Nh/U=%.2f | t=%.0f/%.0f s steps=%d %s | wall=%.1fs\n",
                U_bar, h_max / 1e3, Int(round(lx / 1e3)), Nh_over_U, state.t, t_end,
                state.step, converged ? "OK" : "BLEW-UP", duration)
    end
end

# --- Save summary (with per-run wall-clock duration) -------------------------
# Appended rather than overwritten, so successive invocations into one sweep
# directory (e.g. two ridge heights) keep every run's wall-clock. The sponge
# columns are written only to a new file, so older CSVs keep their layout.
out_csv = joinpath(sweep_base, "sweep_results.csv")
new_csv = !isfile(out_csv)
open(out_csv, "a") do io
    new_csv && write(io, "U_bar,h_max,LX_km,Nh_over_U,duration_s,final_t,steps," *
                         "converged,outdir,K0_sponge,sponge_base_km,LZ_km,pbar_mode\n")
    for r in results
        write(io, "$(r.U_bar),$(r.h_max),$(r.LX / 1e3),$(r.Nh_over_U)," *
                  "$(r.duration_s),$(r.final_t),$(r.steps),$(r.converged),$(r.outdir)" *
                  (new_csv ? ",$(K0_SPONGE),$(SPONGE_BASE / 1e3),$(LZ / 1e3),$(PBAR_MODE)" : "") *
                  "\n")
    end
end

total_wall = sum(r.duration_s for r in results)
@printf("\nSweep complete: %d run(s), total wall-clock %.1f s. Results: %s\n",
        length(results), total_wall, out_csv)
