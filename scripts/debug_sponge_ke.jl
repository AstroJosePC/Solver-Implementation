# =============================================================================
# Debug script: run K0=0 sponge case with detailed KE diagnostics
# Prints KE, max|w|, location of max|w|, and density extremes every 100s
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using MAGNUSP
using MAGNUSP: SVector, primitive_from_conservative, phys_x, phys_z
using Printf
using NCDatasets

const NX = 160
const NZ = 240
const DX = 500.0
const DZ = 500.0
const T_ISOTHERMAL = 239.0
const T_END = 10000.0
const DIAG_DT = 100.0   # every 100 s for diagnostics

const W_TILDE = 0.001
const T0_FORCING = 1200.0
const SIGMA_T = 600.0
const OMEGA = 0.01
const LAMBDA_X = 80e3
const KX = 2.0 * pi / LAMBDA_X

K0 = 0.005  # enable sponge to verify stabilization
sponge_frac = 20.0 / 120.0

config = earth_config(
    nx = NX, nz = NZ, dx = DX, dz = DZ,
    CFL = 0.8, T0 = T_ISOTHERMAL, t_end = T_END,
    output_interval = 999999,
    enable_viscosity = true,
    enable_cooling = false,
    sponge_frac = sponge_frac,
    K0_sponge = K0,
)

wave_params = MAGNUSP.WaveParams(
    amplitude = W_TILDE,
    omega = OMEGA,
    kx = KX,
    sigma_t = SIGMA_T,
    t0 = T0_FORCING,
)

grid = Grid2D(config)
atm = isothermal_atmosphere(config, Float64.(grid.z_full))
Q = MAGNUSP.init_state(grid, atm, config)
Q_bar = copy(Q)

sponge_layer = MAGNUSP.init_sponge(grid, config)
diff_solver = MAGNUSP.DiffusionSolver(grid.nx, grid.nz)
Q_half_x = similar(Q)
Q_half_z = similar(Q)

rho_surface = atm.rho0[grid.ng + 1]

# --- Time loop with diagnostics ---
t = 0.0
step_count = 0
next_diag = 0.0

# Save KE time series
times = Float64[]
ke_values = Float64[]

# Open live log file for incremental writes
outdir = joinpath(@__DIR__, "..", "output", "sponge_validation")
mkpath(outdir)
logfile = joinpath(outdir, "debug_ke_live.txt")
logf = open(logfile, "w")

header = @sprintf("%8s  %12s  %12s  %8s  %12s  %12s\n",
        "t [s]", "KE [J/m]", "max|w| [m/s]", "z_maxw[km]", "min_rho", "max_rho_pert")
println("Debug run: K0 = $(K0)")
println("=" ^ 70)
print(header)
print(logf, "# K0=$(K0)  NZ=$(NZ)  T_END=$(T_END)\n")
print(logf, header)
flush(logf)

while t < config.t_end
    global t, step_count, next_diag
    dt = MAGNUSP.compute_adaptive_dt(Q, grid, config.CFL, config.gamma, config.R_gas)
    dt = min(dt, config.t_end - t)

    MAGNUSP.apply_all_bcs!(Q, Q_bar, grid, config.gamma)
    MAGNUSP.apply_bottom_forcing!(Q, grid, config, t, wave_params, :gaussian)
    MAGNUSP.euler_step!(Q, Q_bar, Q_half_x, Q_half_z, grid, dt, config.gamma, config.g,
                        sponge_layer.K)
    if config.enable_viscosity
        MAGNUSP.diffusion_step!(Q, Q_bar, grid, dt, config, diff_solver)
    end

    t += dt
    step_count += 1

    if t >= next_diag
        # Compute KE
        ke = 0.0
        for j in phys_z(grid), i in phys_x(grid)
            grid.z_full[j] >= 100_000.0 && continue
            rho_bg = Q_bar[1, i, j]
            u_pert = Q[2, i, j] / Q[1, i, j]
            w_pert = Q[3, i, j] / Q[1, i, j]
            ke += 0.5 * rho_bg * (u_pert^2 + w_pert^2) * DX * DZ
        end

        # Find max |w| and its location
        max_w = 0.0
        z_max_w = 0.0
        min_rho = Inf
        max_rho_pert = 0.0
        for j in phys_z(grid)
            for i in phys_x(grid)
                rho = Q[1, i, j]
                w = Q[3, i, j] / rho
                rho_bg = Q_bar[1, i, j]

                aw = abs(w)
                if aw > max_w
                    max_w = aw
                    z_max_w = grid.z_full[j]
                end

                if rho < min_rho
                    min_rho = rho
                end

                rho_rel = abs(rho - rho_bg) / rho_bg
                if rho_rel > max_rho_pert
                    max_rho_pert = rho_rel
                end
            end
        end

        push!(times, t)
        push!(ke_values, ke)

        row = @sprintf("%8.1f  %12.4e  %12.6e  %8.1f  %12.4e  %12.6e\n",
                t, ke, max_w, z_max_w/1e3, min_rho, max_rho_pert)
        print(row)
        print(logf, row)
        flush(logf)

        next_diag += DIAG_DT
    end
end
close(logf)

println("\n--- Summary ---")
@printf("Peak KE = %.6e at t = %.1f s\n", maximum(ke_values), times[argmax(ke_values)])
@printf("Final KE = %.6e\n", ke_values[end])
@printf("Steps: %d, avg dt = %.3f s\n", step_count, T_END / step_count)

println("Live log written to: $logfile")
