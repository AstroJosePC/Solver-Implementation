# =============================================================================
# Strict Mountain-Wave Benchmark Driver (Srivastava et al. 2022, Sec 3.3)
# Matches Tables 9, 10, 11 exactly.
#
# Options (--key=value): --dx, --dz (m), --tend_hr, --bcface=1, --pbar=1,
# --diffsolver=lu|lu_reuse|fft, --difftol=X and --outdir. The default outdir is
# output/mountain_benchmark_strict_dx<dx> with tags appended only when an option
# departs from its default, in the order _bcface_pbarP_dfft_dtol<X>
# (--diffsolver=lu_reuse tags _dlureuse in place of _dfft). --diffsolver and
# --difftol set SimConfig.diffusion_solver and diffusion_dt_tol (src/diffusion.jl);
# they act only when the implicit viscosity is on, and this benchmark runs inviscid
# (Table 10), so they leave its solution unchanged.
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using MAGNUSP
using Printf

function get_arg(name::String, default)
    key = "--" * name * "="
    for a in ARGS
        if startswith(a, key)
            v = split(a, "=", limit=2)[2]
            return typeof(default) <: AbstractFloat ? parse(Float64, v) :
                   typeof(default) <: Integer ? parse(Int, v) :
                   String(v)  # a SubString would not match run_simulation(outdir::String)
        end
    end
    return default
end

# Table 9: 240 x 60 km
const LX = 240e3
const LZ = 60e3
dx_in = get_arg("dx", 1000.0) # Start with 1000m for fast verification, then 250m
dz_in = get_arg("dz", dx_in)
tend_hr = get_arg("tend_hr", 10.0)
# Path of the imposed ground w (reflection_convergence_tests.md Section 14, addendum 3 V3):
# --bcface=1 selects SimConfig.bottom_forcing = :face and tags the default outdir _bcface;
# the default 0 keeps the ghost-cell path of every earlier record.
bcface = get_arg("bcface", 0) > 0
# Background pressure (reflection_convergence_tests.md Section 15, S14V3p): --pbar=1 selects
# SimConfig.pbar_mode = :pressure and tags the default outdir _pbarP; the default 0 keeps
# the pinned :total of every earlier record (Section 11, decision 1).
pbar = get_arg("pbar", 0) > 0
# Implicit diffusion method and refactorization tolerance (header); the defaults lu and
# 0.01 are the path of every earlier record.
diffsolver_in = get_arg("diffsolver", "lu")
diffsolver = diffsolver_in == "lu" ? :lu : diffsolver_in == "lu_reuse" ? :lu_reuse :
             diffsolver_in == "fft" ? :fft_tridiag :
             error("--diffsolver must be lu, lu_reuse or fft, got $(diffsolver_in)")
difftol = get_arg("difftol", 0.01)
outdir = get_arg("outdir", "output/mountain_benchmark_strict_dx$(Int(dx_in))" * (bcface ? "_bcface" : "") *
                           (pbar ? "_pbarP" : "") *
                           (diffsolver === :lu ? "" : diffsolver === :fft_tridiag ? "_dfft" : "_dlureuse") *
                           (difftol == 0.01 ? "" : "_dtol$(difftol)"))

nx = Int(round(LX / dx_in))
nz = Int(round(LZ / dz_in))

# Table 10/11 parameters
config = earth_config(
    nx = nx,
    nz = nz,
    dx = dx_in,
    dz = dz_in,
    CFL = 0.4,
    T0 = 239.0,
    t_end = tend_hr * 3600.0,
    output_interval = 2000, # save some intermediate steps
    enable_viscosity = false,
    enable_cooling = false,
    sponge_frac = 20.0 / 60.0, # 20 km sponge in 60 km domain
    K0_sponge = 0.005,
    pbar_mode = pbar ? :pressure : :total, # :total for the legacy records; see reflection_convergence_tests.md §9
    bottom_forcing = bcface ? :face : :ghost,
    diffusion_solver = diffsolver,
    diffusion_dt_tol = difftol,
)

# Bell-shaped hill (Srivastava Sec 3.3)
# h(x) = b / (1 + (x/a)^2) where b=10m, a=20km.
mountain_params = MAGNUSP.MountainBellParams(
    h_max = 10.0,            # b = 10 m
    a = 20e3,                # a = 20 km
    x0 = 120e3,              # centered in 240km domain
    U_bar = 18.0,            # 18 m/s
)

println("=" ^ 88)
println("STRICT MOUNTAIN-WAVE BENCHMARK (Table 9/10/11)")
println("=" ^ 88)
@printf("Domain: %.1f x %.1f km, dx=%.2f km, dz=%.2f km\n", LX/1e3, LZ/1e3, dx_in/1e3, dz_in/1e3)
@printf("Run: t_end=%.1f hr, CFL=%.2f, sponge=top 20.0 km (frac=%.3f)\n",
        config.t_end/3600, config.CFL, config.sponge_frac)
@printf("pbar_mode=%s, ground w: %s, output: %s\n", config.pbar_mode, config.bottom_forcing, outdir)
@printf("Hill: b=%.1f m, a=%.1f km, center=%.1f km, U_bar=%.1f m/s\n\n",
        mountain_params.h_max, mountain_params.a/1e3, mountain_params.x0/1e3, mountain_params.U_bar)

state = run_simulation(config;
    forcing_type = :mountain_bell,
    wave_params = mountain_params,
    u_init = mountain_params.U_bar,
    outdir = outdir,
    diag_interval = 500,
)

println()
println("Simulation complete.")
@printf("Final time: %.1f s, steps: %d\n", state.t, state.step)
println("Snapshots: $(outdir)/")
