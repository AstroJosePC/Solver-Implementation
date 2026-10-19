# =============================================================================
# Strict Figure-6 parameter-matching run
# Targets Srivastava Fig. 6 dispersion setup with N = 0.0182 rad/s.
# Uses same domain/forcing as Tables 5/6 and adjusts T0 so N matches exactly.
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using MAGNUSP

const GAMMA = 1.4
const G = 9.81
const R_GAS = 287.0
const N_TARGET = 0.0182

# From N = sqrt(gamma-1) * g / cs  and  cs^2 = gamma*R*T0
cs_target = sqrt(GAMMA - 1.0) * G / N_TARGET
T0_strict = cs_target^2 / (GAMMA * R_GAS)

config = earth_config(
    nx = 80,
    nz = 320,
    dx = 500.0,
    dz = 500.0,
    CFL = 0.8,
    T0 = T0_strict,
    t_end = 8000.0,
    output_interval = 50,   # dense snapshots for clean travel diagram
    enable_viscosity = false,
    enable_cooling = false,
    sponge_frac = 20.0 / 160.0,
    K0_sponge = 0.005,
)

lambda_x = 40e3
kx = 2.0 * pi / lambda_x
omega = 0.007

wave_params = MAGNUSP.WaveParams(
    amplitude = 0.001,
    omega = omega,
    kx = kx,
    sigma_t = 600.0,
    t0 = 1200.0,
)

cs = MAGNUSP.sound_speed(config.gamma, config.R_gas, config.T0)
N = MAGNUSP.brunt_vaisala(config.gamma, config.g, cs)
m2 = kx^2 * (N^2 / omega^2 - 1.0)
lz = m2 > 0 ? 2pi / sqrt(m2) : NaN

println("Strict Figure-6 run setup:")
println("  T0 = $(round(config.T0, digits=3)) K (chosen from N target)")
println("  cs = $(round(cs, digits=3)) m/s")
println("  N  = $(round(N, digits=6)) rad/s (target $(N_TARGET))")
println("  Predicted lambda_z = $(round(lz/1e3, digits=3)) km")
println("  Output interval = $(config.output_interval)")

state = run_simulation(config;
    forcing_type = :gaussian,
    wave_params = wave_params,
    outdir = "output/dispersion_test_fig6_strict",
    diag_interval = 200,
)

println("Done. Strict snapshots saved to output/dispersion_test_fig6_strict/")
println("Total steps: $(state.step)")
