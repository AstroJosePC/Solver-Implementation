# =============================================================================
# Eckermann Figure 14 Replication (Fast Wave)
#
# Validates the CO2 radiative damping implementation in MAGNUS-P by structurally
# replicating the vertical displacement amplitude profiles of a faster gravity
# wave (c = 19.0 m/s, lambda_z ~ 15 km), allowing observation over the full
# 200 km atmospheric column response in a computationally feasible timeframe.
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using MAGNUSP
using Printf

function run_damping_case(enable_cooling::Bool, outdir::String)
    mkpath(outdir)

    # Domain parameters
    LX = 100e3
    LZ = 200e3
    dx_in = 1000.0
    dz_in = 500.0

    nx = Int(round(LX / dx_in))
    nz = Int(round(LZ / dz_in))

    config = mars_config(
        nx = nx,
        nz = nz,
        dx = dx_in,
        dz = dz_in,
        CFL = 0.8,
        t_end = 65000.0, # Run long enough to reach steady state
        output_interval = 1000,
        enable_viscosity = true,
        enable_cooling = enable_cooling,
        sponge_frac = 40e3 / LZ, # Top 40 km buffer
        K0_sponge = 0.005,
        T0 = 210.0 # Isothermal Mars
    )

    # Wave parameters
    c = 19.0 # m/s
    lambda_x = 100e3
    omega = 2 * pi * c / lambda_x # ~ 1.1938e-3 rad/s

    # Using 1.0m displacement source -> w' = omega * zeta
    # Keeping it at 1.0m ensures the wave stays in the strictly linear regime
    # preventing non-linear shocks from corrupting the WKB profile.
    w_tilde = omega * 1.0 # ~ 0.0011938 m/s

    wave_params = MAGNUSP.WaveParams(
        amplitude = w_tilde,
        omega = omega,
        kx = 2 * pi / lambda_x,
        t0 = 5000.0 # Used as the tanh ramp-up scaling factor in linear_ramp forcing
    )

    println("=" ^ 88)
    @printf("ECKERMANN FIG 14 GW - DAMPING %s\n", enable_cooling ? "ON" : "OFF")
    println("=" ^ 88)

    state = run_simulation(config;
        forcing_type = :linear_ramp,
        wave_params = wave_params,
        outdir = outdir,
        diag_interval = 1000,
    )
    return state
end

println("Starting Eckermann Fig 14 Damping Validation...")

out_off = "output/eckermann_fig14/cooling_off"
out_on  = "output/eckermann_fig14/cooling_on"

# Run with viscosity only (no cooling)
time_off = @elapsed run_damping_case(false, out_off)

# Run with viscosity AND cooling
time_on = @elapsed run_damping_case(true, out_on)

println()
println("Validation runs complete.")
@printf("Time without cooling: %.2f s\n", time_off)
@printf("Time with cooling: %.2f s\n", time_on)
println("Run analyze_eckermann_fig14_fft.py to compare.")
