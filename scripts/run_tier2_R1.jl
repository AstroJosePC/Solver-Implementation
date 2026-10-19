# =============================================================================
# Tier 2 DSMC Handoff - Case R1 (Baseline Gravity Wave)
#
# Generates realistic fluid boundary conditions for the HARRAH DSMC model.
# This simulation propagates a mid-range Martian gravity wave through a realistic
# background atmosphere up to the exobase, utilizing our validated compressible
# Navier-Stokes solver with molecular viscosity and thermal conduction.
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using MAGNUSP
using Printf

function run_tier2_R1()
    outdir = "output/tier2_R1"
    mkpath(outdir)

    # Domain & Grid parameters
    LX = 200e3
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
        t_end = 50000.0,
        output_interval = 2000, # save full state snapshot every ~2000 steps
        enable_viscosity = true,
        enable_cooling = false,
        sponge_frac = 30e3 / LZ, 
        K0_sponge = 0.005,
        T0 = 210.0 # Placeholder, will be overridden by tabulated atmosphere
    )

    # Wave parameters
    omega = 0.005 # rad/s (~21 min period)
    lambda_x = 200e3 # m
    
    # 0.01 m/s amplitude for linear regime validation
    w_tilde = 0.01

    wave_params = MAGNUSP.WaveParams(
        amplitude = w_tilde,
        omega = omega,
        kx = 2 * pi / lambda_x,
        t0 = 5000.0 # Smooth start
    )

    println("=" ^ 88)
    @printf("TIER 2 - DSMC HANDOFF: R1 (Baseline Gravity Wave)\n")
    println("=" ^ 88)

    atm_file = joinpath(@__DIR__, "..", "external/digitized_figure_11/mars_mcs_uniform.csv")

    state = run_simulation(config;
        forcing_type = :linear_ramp,
        wave_params = wave_params,
        outdir = outdir,
        diag_interval = 20, # Extracts coupling data every 20 steps (~34s)
        coupling_altitude = 100e3,
        atm_file = atm_file
    )
    
    # Export DSMC Handoff Data
    println("\nProcessing DSMC coupling buffer...")
    nt = length(state.coupling_buffer)
    nx_cells = state.grid.nx
    
    times = zeros(nt)
    n_mat = zeros(nt, nx_cells)
    T_mat = zeros(nt, nx_cells)
    u_mat = zeros(nt, nx_cells)
    w_mat = zeros(nt, nx_cells)
    
    for (i, record) in enumerate(state.coupling_buffer)
        times[i] = record[:time][1]
        n_mat[i, :] .= record[:n]
        T_mat[i, :] .= record[:T]
        u_mat[i, :] .= record[:u]
        w_mat[i, :] .= record[:w]
    end
    
    coupling_data = Dict(
        :time => times,
        :x => state.grid.x,
        :n => n_mat,
        :T => T_mat,
        :u => u_mat,
        :w => w_mat
    )
    
    outpath = joinpath(outdir, "dsmc_handoff_100km.nc")
    MAGNUSP.write_coupling_file_nc(coupling_data, config, outpath)
    println("Saved DSMC coupling data to $outpath")

    return state
end

# Execute the benchmark
run_tier2_R1()