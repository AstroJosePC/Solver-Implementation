# =============================================================================
# Mars Orographic GW Damping Validation
#
# Runs a mountain wave benchmark using a Martian atmosphere, comparing the
# result with and without the altitude-dependent CO2 radiative damping 
# parameterization.
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using MAGNUSP
using Printf
using DelimitedFiles

# Use standard Mars atmosphere if available or isothermal
function run_damping_case(enable_cooling::Bool, outdir::String)
    mkpath(outdir)
    # Define domain for Mars GW propagation (large enough to reach thermosphere)
    LX = 600e3
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
        CFL = 0.6,
        t_end = 12 * 3600.0, # 12 hours
        output_interval = 1800, # every 30 mins
        enable_viscosity = true,
        enable_cooling = enable_cooling,
        sponge_frac = 40e3 / LZ, # 40 km sponge at top
        K0_sponge = 0.005,
    )

    mountain_params = MAGNUSP.MountainParams(
        h_max = 500.0,           # m (larger mountain for Mars)
        sigma = 30e3,            # m 
        x0 = 0.5 * LX,           # centered
        U_bar = 20.0,            # m/s
    )

    println("=" ^ 88)
    @printf("MARS OROGRAPHIC GW - DAMPING %s\n", enable_cooling ? "ON" : "OFF")
    println("=" ^ 88)
    
    state = run_simulation(config;
        forcing_type = :mountain,
        wave_params = mountain_params,
        outdir = outdir,
        diag_interval = 500,
    )
    return state
end

println("Starting Mars Orographic Damping Validation...")

out_off = "output/mars_mountain_damping/cooling_off"
out_on  = "output/mars_mountain_damping/cooling_on"

# Run without cooling
time_off = @elapsed run_damping_case(false, out_off)

# Run with cooling
time_on = @elapsed run_damping_case(true, out_on)

println()
println("Validation runs complete.")
@printf("Time without cooling: %.2f s\n", time_off)
@printf("Time with cooling: %.2f s\n", time_on)
println("Check the output directories to compare the density-scaled vertical velocity profiles.")
