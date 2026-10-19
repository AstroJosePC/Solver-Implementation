# =============================================================================
# Eckermann Figure 12 Replication (Mars Radiative Damping Validation)
#
# Simulates a slow, short-vertical-wavelength gravity wave on Mars:
# lambda_x = 100 km, c = 6.6 m/s, lambda_z ~ 5 km
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
    dx_in = 500.0
    dz_in = 250.0

    nx = Int(round(LX / dx_in))
    nz = Int(round(LZ / dz_in))

    config = mars_config(
        nx = nx,
        nz = nz,
        dx = dx_in,
        dz = dz_in,
        CFL = 0.8,
        t_end = 50000.0, # Run to steady state
        output_interval = 2000, 
        enable_viscosity = true,
        enable_cooling = enable_cooling,
        sponge_frac = 40e3 / LZ, # Top 40 km buffer
        K0_sponge = 0.005,
        T0 = 210.0 # Isothermal Mars
    )

    # Wave parameters
    c = 6.6 # m/s
    lambda_x = 100e3
    omega = 2 * pi * c / lambda_x # ~ 4.1469e-4
    
    # Using 100m displacement source -> w' = omega * zeta
    w_tilde = omega * 100.0 # ~ 0.04147 m/s

    wave_params = MAGNUSP.WaveParams(
        amplitude = w_tilde,
        omega = omega,
        kx = 2 * pi / lambda_x,
        # Other params unused for :linear forcing
    )

    println("=" ^ 88)
    @printf("ECKERMANN FIG 12 GW - DAMPING %s\n", enable_cooling ? "ON" : "OFF")
    println("=" ^ 88)
    
    state = run_simulation(config;
        forcing_type = :linear,
        wave_params = wave_params,
        outdir = outdir,
        diag_interval = 1000,
    )
    return state
end

println("Starting Eckermann Fig 12 Damping Validation...")

out_off = "output/eckermann_fig12/cooling_off"
out_on  = "output/eckermann_fig12/cooling_on"

# Run with viscosity only (no cooling)
time_off = @elapsed run_damping_case(false, out_off)

# Run with viscosity AND cooling
time_on = @elapsed run_damping_case(true, out_on)

println()
println("Validation runs complete.")
@printf("Time without cooling: %.2f s\n", time_off)
@printf("Time with cooling: %.2f s\n", time_on)
println("Run analyze_eckermann_fig12.py to compare.")
