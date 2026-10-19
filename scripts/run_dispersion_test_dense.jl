# =============================================================================
# Earth Isothermal Linear Dispersion Test (Dense Output)
# Same physics as run_dispersion_test.jl, but denser snapshot cadence for
# wavelet/FFT panel convergence diagnostics.
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using MAGNUSP

config = earth_config(
    nx = 80,
    nz = 320,
    dx = 500.0,
    dz = 500.0,
    CFL = 0.8,
    T0 = 239.0,
    t_end = 8000.0,
    output_interval = 50,
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

println("Starting dense-output dispersion simulation...")
println("  output_interval = $(config.output_interval) (denser snapshots)")

state = run_simulation(config;
    forcing_type = :gaussian,
    wave_params = wave_params,
    outdir = "output/dispersion_test_dense",
    diag_interval = 200,
)

println("Done. Dense snapshots saved to output/dispersion_test_dense/")
println("Total steps: $(state.step)")
