# =============================================================================
# Critical Layer Benchmark Driver (Srivastava et al. 2022, Sec 3.2)
# Matches Tables 7 and 8 exactly to reproduce Figure 7.
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using MAGNUSP
using Printf

# Table 7: Domain Parameters
# L_x = 400 km, X in [-200, 200] km -> dx = 1 km -> nx = 400
# L_z = 160 km, Z in [0, 160] km -> dz = 0.5 km -> nz = 320
const LX = 400e3
const LZ = 160e3
const DX = 1000.0
const DZ = 500.0
const T_MAX = 10000.0
const CFL = 0.8

nx = Int(round(LX / DX))
nz = Int(round(LZ / DZ))

# Background atmosphere: Earth isothermal T=239K
# Viscosity and thermal conduction are ON
config = earth_config(
    nx = nx,
    nz = nz,
    dx = DX,
    dz = DZ,
    CFL = CFL,
    T0 = 239.0,
    t_end = T_MAX,
    output_interval = 240,   # Every 240 steps (~300s, ~5min) to capture 50, 80, 120, 150 min exactly
    enable_viscosity = true,
    enable_cooling = false,  # Only molecular viscosity and thermal conduction according to text
    sponge_frac = 60.0 / 160.0, # Implement a 60km sponge layer at the top to absorb any escaped waves and suppress vacuum instability
    K0_sponge = 0.05,           # Stronger sponge damping
    pbar_mode = :total,         # legacy records; see reflection_convergence_tests.md §9
)

# Table 8: Time and Forcing Parameters
# Note: Forcing is w(x, z=0, t) = w_tilde * cos(omega*(t-t0) - kx*x) * exp(-(t-t0)^2/(2*sigma_t^2)) * exp(-(x-x0)^2/(2*sigma_x^2))
# Notice the sign on kx is negative in Eq 33 compared to Eq 31, which makes it travel rightward.
lambda_x = 20e3
kx = -2.0 * pi / lambda_x  # Negative so it travels positive-x
omega = 0.01

wave_params = MAGNUSP.WaveParams(
    amplitude = 0.001,
    omega = omega,
    kx = kx,
    sigma_t = 600.0,
    t0 = 1200.0,
    x0 = 80e3,      # Shifted from -120km to 80km because domain is [0, 400km]
    sigma_x = 10e3
)

# Gaussian Mean Wind Profile (Eq 34)
# u_bar(z) = u_0 * exp(-(z-z_0)^2 / (2*sigma_u^2))
u_0 = 100.0
z_0 = 80e3
sigma_u = 10e3

function background_wind(z::Float64)
    return u_0 * exp(-(z - z_0)^2 / (2.0 * sigma_u^2))
end

outdir = "output/critical_layer_test"

println("=" ^ 88)
println("CRITICAL LAYER BENCHMARK (Table 7/8)")
println("=" ^ 88)
@printf("Domain: %.1f x %.1f km, dx=%.2f km, dz=%.2f km\n", LX/1e3, LZ/1e3, DX/1e3, DZ/1e3)
@printf("Run: t_end=%.1f s, CFL=%.2f\n", config.t_end, config.CFL)
@printf("Wind: u_0=%.1f m/s, z_0=%.1f km, sigma_u=%.1f km\n", u_0, z_0/1e3, sigma_u/1e3)
println("Viscosity: ON, Conduction: ON")

state = run_simulation(config;
    forcing_type = :critical,
    wave_params = wave_params,
    u_init = background_wind,
    bc_x_type = :periodic,   # testing if outflow caused instability
    outdir = outdir,
    diag_interval = 50,
)

println()
println("Simulation complete.")
@printf("Final time: %.1f s, steps: %d\n", state.t, state.step)
println("Snapshots: $(outdir)/")
