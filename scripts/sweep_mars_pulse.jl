# =============================================================================
# Parameter Sweep: Pressure Pulse on Mars
#
# Explores pulse spatial widths (sigma_x, sigma_z) to test if a broader source
# can excite gravity modes (frequencies < 1.8 mHz) in the Martian atmosphere.
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using Base.Threads
using MAGNUSP
using Printf
using DelimitedFiles

# Sweep grid
sigma_x_vals = [500.0, 2000.0, 5000.0]  # m
sigma_z_vals = [1000.0, 2000.0, 5000.0] # m
p_tilde_vals = [44.83, 44.83/2, 2*44.83] # Pa
z0_vals = [5e3, 10e3, 20e3, 30e3] # m

const LX = 200e3
const LZ = 40e3
const DX = 1000.0
const DZ = 500.0

nx = Int(round(LX / DX))
nz = Int(round(LZ / DZ))

# Acoustic cutoff for Mars
omega_ac_mHz = 1.8

struct PulseRunResult
    sigma_x::Float64
    sigma_z::Float64
    p_tilde::Float64
    z0::Float64
    duration_s::Float64
    outdir::String
end

results = Vector{PulseRunResult}()
results_lock = ReentrantLock()

println("Starting parallel sweep for Mars Pressure Pulse...")
println("Number of threads: ", Threads.nthreads())

# Base output directory
sweep_base = abspath(joinpath(@__DIR__, "..", "output", "sweep_mars_pulse"))
mkpath(sweep_base)

combinations = Iterators.product(sigma_x_vals, sigma_z_vals, p_tilde_vals, z0_vals)
combinations_arr = collect(combinations)

# Pre-create run directories sequentially to be safe with external drive metadata
for params in combinations_arr
    sigma_x, sigma_z, p_tilde, z0 = params
    outdir = normpath(joinpath(sweep_base, "run_sx$(sigma_x)_sz$(sigma_z)_pt$(p_tilde)_z0$(z0)"))
    mkpath(outdir)
end

Threads.@threads for params in combinations_arr
    sigma_x, sigma_z, p_tilde, z0 = params

    outdir = normpath(joinpath(sweep_base, "run_sx$(sigma_x)_sz$(sigma_z)_pt$(p_tilde)_z0$(z0)"))

    # Keep outputs clean when re-running the sweep
    for f in filter(f -> endswith(f, ".nc"), readdir(outdir; join=true))
        rm(f; force=true)
    end

    atm_file = normpath(joinpath(@__DIR__, "..", "external", "digitized_figure_11", "mars_mcs_uniform.csv"))

    config = mars_config(
        nx = nx,
        nz = nz,
        dx = DX,
        dz = DZ,
        CFL = 0.4,
        t_end = 4000.0, # Run for ~1 hour simulation time to capture lower frequencies
        output_interval = 5, # in model steps; dense output needed for Fig 12-15 analysis
        enable_viscosity = true,
        enable_cooling = false,
        sponge_frac = 0.0,
        K0_sponge = 0.05
    )

    # Initial Pressure Pulse Parameters (Eq 43)
    pulse_params = MAGNUSP.PulseParams(
        P_tilde = p_tilde,
        sigma_x = sigma_x,
        sigma_z = sigma_z,
        x0 = 100e3,
        z0 = z0
    )

    duration = @elapsed begin
        state = run_simulation(config;
            forcing_type = :none,
            pulse_params = pulse_params,
            bc_x_type = :outflow,
            bc_bottom_type = :outflow,
            atm_file = atm_file,
            outdir = outdir,
            diag_interval = 2000
        )
    end

    res = PulseRunResult(sigma_x, sigma_z, p_tilde, z0, duration, outdir)

    lock(results_lock) do
        push!(results, res)
        @printf("Completed sigma_x=%.1f, sigma_z=%.1f, p_tilde=%.1f, z0=%.1f | Time=%.1fs\n",
                sigma_x, sigma_z, p_tilde, z0, duration)
    end
end

# Save summary
out_csv = joinpath(@__DIR__, "..", "output", "mars_pulse_sweep_results.csv")
mkpath(dirname(out_csv))
open(out_csv, "w") do io
    write(io, "sigma_x,sigma_z,p_tilde,z0,duration_s,outdir\n")
    for r in results
        write(io, "$(r.sigma_x),$(r.sigma_z),$(r.p_tilde),$(r.z0),$(r.duration_s),$(r.outdir)\n")
    end
end
println("Sweep complete! Results saved to ", out_csv)
