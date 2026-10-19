using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using MAGNUSP
using MAGNUSP: phys_x, phys_z, primitive_from_conservative, SVector
using Printf

const NX = 160; const NZ = 240; const DX = 500.0; const DZ = 500.0

function check_state(Q, grid, gamma, label)
    min_P = Inf
    worst_i = 0; worst_j = 0
    for j in phys_z(grid), i in phys_x(grid)
        Qi = SVector(Q[1,i,j], Q[2,i,j], Q[3,i,j], Q[4,i,j])
        rho, u, w, P = primitive_from_conservative(Qi, gamma)
        if P < min_P
            min_P = P; worst_i = i; worst_j = j
        end
    end
    if min_P < 0
        @printf("  %s: min_P = %.4e at (i=%d, j=%d)\n", label, min_P, worst_i, worst_j)
    end
    return min_P
end

function run_bc_test(; bc_type::Symbol)
    config = earth_config(
        nx = NX, nz = NZ, dx = DX, dz = DZ, CFL = 0.8, T0 = 239.0,
        t_end = 8000.0, output_interval = 999999,
        enable_viscosity = true, enable_cooling = false,
        sponge_frac = 20.0 / 120.0, K0_sponge = 0.0,
    )
    wave_params = MAGNUSP.WaveParams(
        amplitude = 0.001, omega = 0.01, kx = 2*pi/80e3,
        sigma_t = 600.0, t0 = 1200.0,
    )
    grid = Grid2D(config)
    atm  = isothermal_atmosphere(config, Float64.(grid.z_full))
    Q    = MAGNUSP.init_state(grid, atm, config)
    Q_bar = copy(Q)
    sponge = MAGNUSP.init_sponge(grid, config)
    diff_solver = MAGNUSP.DiffusionSolver(grid.nx, grid.nz)
    Q_hx = similar(Q); Q_hz = similar(Q)

    t = 0.0; step_count = 0
    ke_peak = 0.0
    next_print = 500.0

    println("K0=0 (no sponge), bc_type=$bc_type")
    @printf("  R_gas=%.1f, T0=%.1f, gamma=%.2f\n", config.R_gas, config.T0, config.gamma)

    function compute_ke()
        ke = 0.0
        for j in phys_z(grid), i in phys_x(grid)
            u = Q[2, i, j] / Q[1, i, j]
            w = Q[3, i, j] / Q[1, i, j]
            ke += 0.5 * Q_bar[1, i, j] * (u^2 + w^2) * DX * DZ
        end
        return ke
    end

    while t < config.t_end
        min_P = check_state(Q, grid, config.gamma, "step$(step_count)")
        if min_P < 0
            ng = grid.ng
            for j in phys_z(grid), i in phys_x(grid)
                Qi = SVector(Q[1,i,j], Q[2,i,j], Q[3,i,j], Q[4,i,j])
                rho, u, w, P = primitive_from_conservative(Qi, config.gamma)
                if P < 0
                    z_km = (j - ng - 1) * DZ / 1e3
                    @printf("  NEG P at j=%d (z=%.1f km): rho=%.4e, u=%.4e, w=%.4e, P=%.4e\n",
                            j, z_km, rho, u, w, P)
                    break
                end
            end
            println("CRASH at step $step_count, t=$t")
            return
        end

        dt = MAGNUSP.compute_adaptive_dt(Q, grid, config.CFL, config.gamma, config.R_gas)
        dt = min(dt, config.t_end - t)

        MAGNUSP.apply_all_bcs!(Q, Q_bar, grid, config.gamma;
                               bc_type = bc_type,
                               R_gas = config.R_gas,
                               T0 = config.T0)
        MAGNUSP.apply_bottom_forcing!(Q, grid, config, t, wave_params, :gaussian)
        MAGNUSP.euler_step!(Q, Q_bar, Q_hx, Q_hz, grid, dt, config.gamma, config.g, sponge.K)
        if config.enable_viscosity
            MAGNUSP.diffusion_step!(Q, Q_bar, grid, dt, config, diff_solver)
        end
        t += dt; step_count += 1

        if step_count % 50 == 0
            ke = compute_ke()
            @printf("  step %d: t=%.1f, KE=%.4e\n", step_count, t, ke)
        end
    end

    println("Completed successfully.")
end

println("=== Characteristic BC (debug, 10 steps) ===")
run_bc_test(bc_type = :characteristic)
