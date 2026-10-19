# =============================================================================
# Parameter Sweep: Mountain Wave on Mars
#
# Explores U_bar and sigma to find proper AGWs (gravity modes) with frequencies 
# lower than the acoustic cut off in a Martian atmosphere.
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using Base.Threads
using MAGNUSP
using Printf
using DelimitedFiles

# Sweep grid
U_bar_vals = [5.0, 10.0, 20.0]        # m/s
sigma_vals = [5e3, 20e3, 50e3]        # m
h_max_vals = [10.0]                   # m

# Domain config (expanded to prevent boundary cliffs for 50km mountain)
const LX = 600e3
const LZ = 60e3
const DX = 1000.0
const DZ = 500.0

nx = Int(round(LX / DX))
nz = Int(round(LZ / DZ))

# Basic properties for Mars
# omega_ac approx 1.8 mHz
omega_ac_mHz = 1.8

struct RunResult
    U_bar::Float64
    sigma::Float64
    h_max::Float64
    duration_s::Float64
    kx::Float64
    intrinsic_omega_mHz::Float64
    is_gravity_mode::Bool
    outdir::String
end

results = Vector{RunResult}()
results_lock = ReentrantLock()

println("Starting parallel sweep for Mars Mountain Wave...")
println("Number of threads: ", Threads.nthreads())

# Base output directory
sweep_base = abspath(joinpath(@__DIR__, "..", "output", "sweep_mars_mountain"))
mkpath(sweep_base)

# Iterate over parameter grid
combinations = Iterators.product(U_bar_vals, sigma_vals, h_max_vals)
combinations_arr = collect(combinations)

# Pre-create run directories sequentially to be safe with external drive metadata
for params in combinations_arr
    U_bar, sigma, h_max = params
    outdir = normpath(joinpath(sweep_base, "run_U$(U_bar)_sig$(sigma/1e3)"))
    mkpath(outdir)
end

Threads.@threads for params in combinations_arr
    U_bar, sigma, h_max = params
    
    # Calculate expected parameters
    kx = 2.0 * pi / (4.0 * sigma)  # rough dominant kx for a gaussian hill of half-width sigma
    omega_intrinsic = kx * U_bar
    omega_mHz = (omega_intrinsic / (2 * pi)) * 1000.0
    
    is_gravity_mode = omega_mHz < omega_ac_mHz
    
    outdir = normpath(joinpath(sweep_base, "run_U$(U_bar)_sig$(sigma/1e3)"))
    
    config = mars_config(
        nx = nx,
        nz = nz,
        dx = DX,
        dz = DZ,
        CFL = 0.4,
        t_end = 10 * 3600.0, # 10 hours for gravity wave to fully establish
        output_interval = 20, # Output every ~35s for proper wavelet sampling
        enable_viscosity = true,
        enable_cooling = false,
        sponge_frac = 20e3 / LZ,
        K0_sponge = 0.05,
        pbar_mode = :total,   # legacy records; see reflection_convergence_tests.md §9
    )
    
    mountain_params = MAGNUSP.MountainParams(
        h_max = h_max,
        sigma = sigma,
        x0 = 0.5 * LX,
        U_bar = U_bar
    )
    
    # Measure execution time
    duration = @elapsed begin
        state = run_simulation(config;
            forcing_type = :mountain,
            wave_params = mountain_params,
            u_init = U_bar,           # CRITICAL: atmosphere must flow over mountain
            bc_x_type = :outflow,     # CRITICAL: prevent wrap-around interference
            outdir = outdir,
            diag_interval = 2000
        )
    end
    
    res = RunResult(U_bar, sigma, h_max, duration, kx, omega_mHz, is_gravity_mode, outdir)
    
    lock(results_lock) do
        push!(results, res)
        @printf("Completed U_bar=%.1f, sigma=%.1fkm | Omega=%.3f mHz (<%.1f? %s) | Time=%.1fs\n", 
                U_bar, sigma/1e3, omega_mHz, omega_ac_mHz, is_gravity_mode, duration)
    end
end

# Save summary
out_csv = joinpath(@__DIR__, "..", "output", "mars_mountain_sweep_results.csv")
mkpath(dirname(out_csv))
open(out_csv, "w") do io
    write(io, "U_bar,sigma,h_max,duration_s,kx,intrinsic_omega_mHz,is_gravity_mode,outdir\n")
    for r in results
        write(io, "$(r.U_bar),$(r.sigma),$(r.h_max),$(r.duration_s),$(r.kx),$(r.intrinsic_omega_mHz),$(r.is_gravity_mode),$(r.outdir)\n")
    end
end
println("Sweep complete! Results saved to ", out_csv)
