# =============================================================================
# Tier 2 DSMC Handoff - Case R2 with Klemp-Durran BC on a 120 km domain
#
# Same wave parameters as run_tier2_R2.jl, but with:
# - Domain height reduced from 200 km to 120 km
# - Klemp-Durran radiation BC at the top (120 km)
# - Thin sponge from 110-120 km (acoustic residual insurance)
#
# Rationale: At 200 km, rho ~ 8e-13 kg/m^3 and the KD pressure perturbation
# exceeds P_bar, requiring aggressive clamping that neutralizes the BC.
# At 120 km, rho ~ 1e-7 kg/m^3 (6 orders denser), so KD can operate with
# much less clamping. The handoff altitude (100 km) is 10 km below the
# sponge onset, preserving clean data extraction.
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using MAGNUSP
using Printf

function run_tier2_R2_kd_120km()
    outdir = "output/tier2_R2_kd_120km"
    mkpath(outdir)

    # Domain & Grid parameters
    LX = 200e3      # horizontal: same as original
    LZ = 120e3      # vertical: reduced from 200 km to 120 km
    dx_in = 1000.0
    dz_in = 500.0

    nx = Int(round(LX / dx_in))   # 200
    nz = Int(round(LZ / dz_in))   # 240

    config = mars_config(
        nx = nx,
        nz = nz,
        dx = dx_in,
        dz = dz_in,
        CFL = 0.8,
        t_end = 50000.0,
        output_interval = 100,
        enable_viscosity = true,
        enable_cooling = false,
        sponge_frac = 10e3 / LZ,   # sponge from 110-120 km (10 km thick)
        K0_sponge = 0.002,          # reduced: acoustic residual insurance
        T0 = 210.0
    )

    # Wave parameters (identical to R2)
    omega = 0.005     # rad/s (gravity wave regime)
    lambda_x = 200e3  # m
    w_tilde = 0.05    # nonlinear regime

    wave_params = MAGNUSP.WaveParams(
        amplitude = w_tilde,
        omega = omega,
        kx = 2 * pi / lambda_x,
        t0 = 5000.0
    )

    println("=" ^ 88)
    @printf("TIER 2 - DSMC HANDOFF: R2 with Klemp-Durran BC (120 km domain)\n")
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
run_tier2_R2_kd_120km()
