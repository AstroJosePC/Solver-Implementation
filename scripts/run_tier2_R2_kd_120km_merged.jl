# =============================================================================
# Tier 2 DSMC Handoff - Case R2 with Klemp-Durran BC on a 120 km domain
# using merged MCS-MCD background atmospheres for dayside and nightside.
#
# Tracks execution wall-clock time for performance monitoring.
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using MAGNUSP
using Printf
using Dates

function run_case(case_name::String, atm_file::String, outdir::String, w_tilde::Float64)
    mkpath(outdir)
    
    # Domain & Grid parameters
    LX = 200e3      # horizontal: same as original
    LZ = 120e3      # vertical: 120 km
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

    # Wave parameters
    # omega lowered from 0.005 to 0.003 rad/s so the forced mode propagates
    # vertically through the whole column. The turning-point frequency (m^2=0) for
    # the merged atmosphere with lambda_x=200 km is ~0.00363 rad/s at z~10 km (the
    # minimum-N^2 level); at 0.005 the mode was evanescent over ~95% of the column.
    # omega=0.003 gives m^2>0 everywhere with ~17% margin, lambda_z ~ 60-150 km,
    # and stays well below the acoustic cutoff (~0.0108 rad/s).
    omega = 0.003     # rad/s (internal gravity wave, propagating regime)
    lambda_x = 200e3  # m

    wave_params = MAGNUSP.WaveParams(
        amplitude = w_tilde,
        omega = omega,
        kx = 2 * pi / lambda_x,
        t0 = 5000.0
    )

    println("=" ^ 80)
    @printf("RUNNING CASE: %s\n", case_name)
    @printf("Atmosphere Profile: %s\n", atm_file)
    println("=" ^ 80)

    start_time = now()
    
    state = run_simulation(config;
        forcing_type = :linear_ramp,
        wave_params = wave_params,
        bc_type = :klemp_durran,
        outdir = outdir,
        diag_interval = 20,
        coupling_altitude = 100e3,
        atm_file = atm_file
    )

    end_time = now()
    elapsed_seconds = (end_time - start_time).value / 1000.0
    @printf("Case %s completed in %.2f seconds (wall-clock).\n", case_name, elapsed_seconds)

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
    println("Saved DSMC coupling data to $outpath\n")
    
    return elapsed_seconds
end

function run_all_cases()
    atm_dir = joinpath(@__DIR__, "..", "external/digitized_figure_11")
    
    # Dayside run (calibrated w_tilde = 0.15 for 20-40% target)
    day_atm = joinpath(atm_dir, "mars_mcs_mcd_merged_day.csv")
    day_out = "output/tier2_R2_kd_120km_merged_day"
    day_elapsed = run_case("Dayside (DD2-like)", day_atm, day_out, 0.15)
    
    # Nightside run (calibrated w_tilde = 0.25 for up to 70% target)
    night_atm = joinpath(atm_dir, "mars_mcs_mcd_merged_night.csv")
    night_out = "output/tier2_R2_kd_120km_merged_night"
    night_elapsed = run_case("Nightside (DD6-like)", night_atm, night_out, 0.25)
    
    println("=" ^ 80)
    @printf("SUMMARY OF SIMULATION RUNS\n")
    @printf("Dayside Case duration:   %.2f s\n", day_elapsed)
    @printf("Nightside Case duration: %.2f s\n", night_elapsed)
    println("=" ^ 80)
end

run_all_cases()
