# =============================================================================
# Orographic Mountain-Wave Benchmark Driver
#
# Srivastava-style benchmark over an isolated hill using steady bottom forcing:
#   w(x, z=0) = U_bar * dh/dx
#   h(x) = h_max * exp(-(x-x0)^2 / (2*sigma^2))
#
# This driver uses run_simulation with forcing_type=:mountain.
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
                   v
        end
    end
    return default
end

# Domain inspired by Srivastava mountain test scale
const LX = 240e3
const LZ = 80e3      # includes 20 km sponge at top
dx_in = get_arg("dx", 250.0)
dz_in = get_arg("dz", dx_in)
tend_hr = get_arg("tend_hr", 10.0)
outdir = get_arg("outdir", "output/mountain_benchmark")

nx = Int(round(LX / dx_in))
nz = Int(round(LZ / dz_in))

config = earth_config(
    nx = nx,
    nz = nz,
    dx = dx_in,
    dz = dz_in,
    CFL = 0.4,
    T0 = 239.0,
    t_end = tend_hr * 3600.0,
    output_interval = 3000,
    enable_viscosity = false,
    enable_cooling = false,
    sponge_frac = 20e3 / LZ, # top 20 km buffer
    K0_sponge = 0.005,
)

mountain_params = MAGNUSP.MountainParams(
    h_max = 10.0,            # m
    sigma = 20e3,            # m (half-width scale)
    x0 = 0.5 * LX,           # centered hill
    U_bar = 18.0,            # m/s
)

println("=" ^ 88)
println("MOUNTAIN-WAVE BENCHMARK")
println("=" ^ 88)
@printf("Domain: %.1f x %.1f km, dx=%.2f km, dz=%.2f km\n", LX/1e3, LZ/1e3, dx_in/1e3, dz_in/1e3)
@printf("Run: t_end=%.1f hr, CFL=%.2f, sponge=top %.1f km\n", config.t_end/3600, config.CFL, config.sponge_frac*LZ/1e3)
@printf("Hill: h_max=%.1f m, sigma=%.1f km, x0=%.1f km, U_bar=%.1f m/s\n\n",
        mountain_params.h_max, mountain_params.sigma/1e3, mountain_params.x0/1e3, mountain_params.U_bar)

state = run_simulation(config;
    forcing_type = :mountain,
    wave_params = mountain_params,
    outdir = outdir,
    diag_interval = 500,
)

println()
println("Simulation complete.")
@printf("Final time: %.1f s, steps: %d\n", state.t, state.step)
println("Snapshots: $(outdir)/")
println()
println("Suggested next checks:")
println("  1) Compare density-scaled w contours against Queney-style analytical pattern")
println("  2) Compute correlation / normalized RMSE of simulated vs analytical field")
println("  3) Verify phase-line tilt and vertical wavelength consistency")
