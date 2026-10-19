# =============================================================================
# Tier 2 DSMC Handoff - Case R3 with Klemp-Durran Radiation BC
#
# Same wave parameters as run_tier2_R3.jl, but with the Klemp-Durran (1983)
# radiation boundary condition at the top of the domain instead of Eq. 19.
# The sponge is retained at reduced strength to absorb acoustic residuals.
#
# NOTE: The KD BC is derived from the gravity wave dispersion relation.
# R3 is an acoustic wave (omega=0.02 > omega_ac~0.011), so the KD radiation
# condition is not physically exact for this mode. This run tests whether
# KD's transmissive behavior still outperforms Eq.19's partial-wall behavior
# for acoustic waves.
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using MAGNUSP
using Printf

function run_tier2_R3_kd()
    outdir = "output/tier2_R3_kd"
    mkpath(outdir)

    # Domain & Grid parameters (identical to R3)
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
        t_end = 25000.0,
        output_interval = 30,
        enable_viscosity = true,
        enable_cooling = false,
        sponge_frac = 30e3 / LZ,
        K0_sponge = 0.002,  # reduced sponge: acoustic residual insurance only
        T0 = 210.0
    )

    # Wave parameters (identical to R3)
    omega = 0.02     # rad/s (~5.2 min period, acoustic regime)
    lambda_x = 200e3 # m
    w_tilde = 0.01   # linear acoustic regime

    wave_params = MAGNUSP.WaveParams(
        amplitude = w_tilde,
        omega = omega,
        kx = 2 * pi / lambda_x,
        t0 = 1200.0  # shorter ramp — acoustic waves arrive fast
    )

    println("=" ^ 88)
    @printf("TIER 2 - DSMC HANDOFF: R3 with Klemp-Durran Radiation BC\n")
    println("=" ^ 88)

    atm_file = joinpath(@__DIR__, "..", "external/digitized_figure_11/mars_mcs_uniform.csv")

    state = run_simulation(config;
        forcing_type = :linear_ramp,
        wave_params = wave_params,
        bc_type = :klemp_durran,
        outdir = outdir,
        diag_interval = 20,
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

# Execute
run_tier2_R3_kd()
