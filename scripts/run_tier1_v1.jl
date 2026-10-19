# =============================================================================
# Tier 1, Case V1: Linear GW Propagation without Dissipation (Mars Isothermal)
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using MAGNUSP

config = mars_config(
    nx = 300,
    nz = 400,
    dx = 500.0,
    dz = 500.0,
    CFL = 0.8,
    t_end = 10000.0,
    output_interval = 50,
    enable_viscosity = false,
    enable_cooling = false,
    sponge_frac = 30.0 / 200.0, # 30 km sponge out of 200 km domain
    K0_sponge = 0.005
)

# Wave forcing parameters
lambda_x = 150e3
kx = 2.0 * pi / lambda_x
omega = 0.005

wave_params = MAGNUSP.WaveParams(
    amplitude = 0.001,
    omega = omega,
    kx = kx,
    sigma_t = 600.0,
    t0 = 1200.0
)

println("Running Tier 1 Case V1: Isothermal Mars (No Dissipation)")
println("  Lx = $(config.nx * config.dx / 1e3) km")
println("  Lz = $(config.nz * config.dz / 1e3) km")
println("  omega = $(omega) rad/s")
println("  lambda_x = $(lambda_x / 1e3) km")

state = run_simulation(config;
    forcing_type = :gaussian,
    wave_params = wave_params,
    outdir = "output/tier1_v1",
    diag_interval = 100
)

println("Tier 1 Case V1 Complete. Outputs saved to output/tier1_v1/")
