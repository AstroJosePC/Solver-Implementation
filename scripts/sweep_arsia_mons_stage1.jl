# =============================================================================
# Arsia Mons Stage 1 mountain-wave sweep (AMEC mid-atmosphere focus)
#
# Launches orographic gravity waves over an Arsia-scale Gaussian hill using the
# real MY31 MCD background atmosphere and resolves the wave field through the
# 45-53 km Arsia Mons Elongated Cloud (AMEC) formation level.
#
# Domain top is 90 km with a sponge over 60-90 km. An earlier 155 km attempt
# blew up: a stationary mountain wave amplifies as rho^-1/2 into the thermosphere
# (supersonic ~800 m/s by 120 km) and breaks, because the viscosity cap (nu<=1e4)
# throttles the molecular viscosity that physically damps these waves above
# ~100 km. The mid-atmosphere (<=80 km) is clean, so we cut the domain below the
# breaking zone; the 120 km DSMC coupling product is deferred to a later stage
# (would need realistic thermospheric viscosity / a raised nu cap).
#
# Stage 1 baseline (per plan): uniform EASTWARD steering wind (u_init = U_bar)
# with the outflow x-BC (the validated config). Sweeps the steering wind U_bar
# and the effective hill height h_max. The real reversing u_zonal shear and its
# critical levels are deferred to Stage 2 (needs periodic-x + wider domain).
#
# Usage:
#   julia -t auto --project=. scripts/sweep_arsia_mons_stage1.jl smoke
#       -> single short config (U_bar=40, h_max=1 km, t_end=1 h) to confirm
#          stability of the 90 km domain and measure wall-clock time.
#   julia -t auto --project=. scripts/sweep_arsia_mons_stage1.jl
#       -> full 3x3 sweep (U_bar x h_max) at the production t_end.
#   julia --project=. scripts/sweep_arsia_mons_stage1.jl ext [u=40 h=3 hours=6]
#       -> single-run time extension (default U=40, h=3 km, 6 h) into
#          output/arsia_mons_stage1_ext<hours>h/, for modes whose momentum budget
#          was still trending at the 3 h sweep t_end. The output directory is
#          named after the horizon, so `hours=12` writes to ..._ext12h/ and leaves
#          the 6 h tree intact. The solver itself is serial (no @threads in src/),
#          so -t auto only helps the multi-combination sweep modes.
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using Base.Threads
using MAGNUSP
using Printf
using DelimitedFiles

const SMOKE = "smoke" in ARGS
const EXT = "ext" in ARGS   # single-run time extension (see mode block below)

# --- Background atmosphere (real MCD MY31, 0-120 km) ---
# The 90 km domain (+ ghosts) lies entirely within the real MCD column, so the
# isothermal-capped 156 km file is not needed here (that file is reserved for a
# later thermospheric stage).
const ATM_FILE = abspath(joinpath(@__DIR__, "..", "external", "mcd_profiles",
                                   "mars_arsia_mons_my31.csv"))

# --- Domain / numerics ---
const LX = 2000e3       # horizontal extent [m]
const LZ = 90e3         # vertical extent [m] (90 km; sponge 60-90 km)
const DX = 4000.0       # [m] -> nx = 500
const DZ = 500.0        # [m] -> nz = 180
const SIGMA = 100e3     # Gaussian hill half-width [m] (Arsia-scale), fixed

nx = Int(round(LX / DX))
nz = Int(round(LZ / DZ))

# Sponge occupies the top 30 km (60-90 km), ~1.4 vertical wavelengths; its base
# (K=0) sits 7 km above the 53 km AMEC level, leaving the cloud region clean.
const SPONGE_FRAC = 30e3 / LZ
const K0_SPONGE = 0.05

# Representative N for the nonlinearity diagnostic Nh/U (lower-atmosphere mean).
const N_REF = 0.018

# --- Sweep grid ---
if SMOKE
    U_bar_vals = [40.0]
    h_max_vals = [1000.0]
    t_end = 3600.0                 # 1 hour
    sweep_base = abspath(joinpath(@__DIR__, "..", "output", "arsia_mons_stage1_smoke"))
elseif EXT
    # Time extension of a single mode: the 3 h sweep left the momentum budget
    # (MF growth, mean-flow drift) still trending even where the band-amplitude
    # CoV looked settled. Defaults U=40, h=3 km, 6 h; override with u=, h= (km),
    # hours= args, e.g. `... sweep_arsia_mons_stage1.jl ext u=40 h=3 hours=6`.
    getarg(k, d) = (i = findfirst(a -> startswith(a, k * "="), ARGS);
                    i === nothing ? d : parse(Float64, split(ARGS[i], "=")[2]))
    U_bar_vals = [getarg("u", 40.0)]
    h_max_vals = [getarg("h", 3.0) * 1e3]
    ext_hours = getarg("hours", 6.0)
    t_end = ext_hours * 3600.0
    # Output directory is named after the horizon so successive extensions do not
    # overwrite each other (hours=6 still resolves to .../arsia_mons_stage1_ext6h).
    hlabel = ext_hours == round(ext_hours) ? string(Int(round(ext_hours))) : string(ext_hours)
    sweep_base = abspath(joinpath(@__DIR__, "..", "output",
                                  "arsia_mons_stage1_ext$(hlabel)h"))
else
    U_bar_vals = [20.0, 40.0, 57.0]        # m/s (eastward magnitude; westward jet mirrored)
    h_max_vals = [100.0, 1000.0, 3000.0]   # m (linear -> nonlinear onset)
    t_end = 10800.0                # 3 hours (90 km column reaches quasi-steady fast)
    sweep_base = abspath(joinpath(@__DIR__, "..", "output", "arsia_mons_stage1"))
end

const OUTPUT_INTERVAL = 100   # snapshot cadence (~73 s at dt~0.74 s)
const DIAG_INTERVAL = 500     # diagnostics / blowup-check cadence

struct RunResult
    U_bar::Float64
    h_max::Float64
    Nh_over_U::Float64
    duration_s::Float64
    final_t::Float64
    steps::Int
    converged::Bool
    outdir::String
end

results = Vector{RunResult}()
results_lock = ReentrantLock()

@printf("Arsia Mons Stage 1 sweep (%s)\n", SMOKE ? "SMOKE" : "FULL")
@printf("  atm_file : %s\n", ATM_FILE)
@printf("  grid     : nx=%d nz=%d (dx=%.0f dz=%.0f) | LX=%.0f km LZ=%.0f km\n",
        nx, nz, DX, DZ, LX/1e3, LZ/1e3)
@printf("  sponge   : top %.0f km (frac=%.3f, K0=%.3f) | AMEC level ~45-53 km\n",
        SPONGE_FRAC*LZ/1e3, SPONGE_FRAC, K0_SPONGE)
@printf("  t_end    : %.0f s | threads=%d\n", t_end, Threads.nthreads())

isfile(ATM_FILE) || error("Background atmosphere not found: $ATM_FILE\n" *
    "Generate it with: python scripts/prep_mars_arsia_mons_atmosphere.py " *
    "--output external/mcd_profiles/mars_arsia_mons_my31.csv")

mkpath(sweep_base)

combinations_arr = collect(Iterators.product(U_bar_vals, h_max_vals))

# Pre-create run directories sequentially (safe with external-drive metadata).
for (U_bar, h_max) in combinations_arr
    mkpath(normpath(joinpath(sweep_base, "run_U$(U_bar)_h$(h_max/1e3)")))
end

Threads.@threads for params in combinations_arr
    U_bar, h_max = params
    outdir = normpath(joinpath(sweep_base, "run_U$(U_bar)_h$(h_max/1e3)"))
    Nh_over_U = N_REF * h_max / U_bar

    config = mars_config(
        nx = nx,
        nz = nz,
        dx = DX,
        dz = DZ,
        CFL = 0.4,
        t_end = t_end,
        output_interval = OUTPUT_INTERVAL,
        enable_viscosity = true,
        enable_cooling = false,
        sponge_frac = SPONGE_FRAC,
        K0_sponge = K0_SPONGE,
        pbar_mode = :total,   # legacy records; see reflection_convergence_tests.md §9
    )

    mountain_params = MAGNUSP.MountainParams(
        h_max = h_max,
        sigma = SIGMA,
        x0 = 0.5 * LX,
        U_bar = U_bar,
    )

    local state
    duration = @elapsed begin
        state = run_simulation(config;
            forcing_type = :mountain,
            wave_params = mountain_params,
            u_init = U_bar,                # uniform eastward steering wind
            bc_x_type = :outflow,          # prevents wrap-around (unidirectional wind)
            atm_file = ATM_FILE,
            outdir = outdir,
            diag_interval = DIAG_INTERVAL,
        )
    end

    converged = state.t >= t_end * 0.999
    res = RunResult(U_bar, h_max, Nh_over_U, duration, state.t, state.step,
                    converged, outdir)

    lock(results_lock) do
        push!(results, res)
        @printf("Done U=%.0f h=%.2fkm | Nh/U=%.2f | t=%.0f/%.0f s steps=%d %s | wall=%.1fs\n",
                U_bar, h_max/1e3, Nh_over_U, state.t, t_end, state.step,
                converged ? "OK" : "BLEW-UP", duration)
    end
end

# --- Save summary (with per-run wall-clock duration) ---
out_csv = joinpath(sweep_base, "sweep_results.csv")
open(out_csv, "w") do io
    write(io, "U_bar,h_max,Nh_over_U,duration_s,final_t,steps,converged,outdir\n")
    for r in results
        write(io, "$(r.U_bar),$(r.h_max),$(r.Nh_over_U),$(r.duration_s)," *
                  "$(r.final_t),$(r.steps),$(r.converged),$(r.outdir)\n")
    end
end

total_wall = sum(r.duration_s for r in results)
@printf("\nSweep complete: %d run(s), total wall-clock %.1f s. Results: %s\n",
        length(results), total_wall, out_csv)
