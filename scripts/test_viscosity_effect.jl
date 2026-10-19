using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using MAGNUSP
using MAGNUSP: phys_x, phys_z
using Printf

const NX = 160; const NZ = 240; const DX = 500.0; const DZ = 500.0

function run_k0_test(; enable_viscosity::Bool)
    config = earth_config(
        nx = NX, nz = NZ, dx = DX, dz = DZ, CFL = 0.8, T0 = 239.0,
        t_end = 8000.0, output_interval = 999999,
        enable_viscosity = enable_viscosity, enable_cooling = false,
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
    next_print = 1000.0

    visc_label = enable_viscosity ? "ON (altitude-dependent)" : "OFF"
    println("K0=0, viscosity=$visc_label")
    println("t [s]      KE           KE_peak")
    println("-" ^ 45)

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
        dt = MAGNUSP.compute_adaptive_dt(Q, grid, config.CFL, config.gamma, config.R_gas)
        dt = min(dt, config.t_end - t)
        MAGNUSP.apply_all_bcs!(Q, Q_bar, grid, config.gamma)
        MAGNUSP.apply_bottom_forcing!(Q, grid, config, t, wave_params, :gaussian)
        MAGNUSP.euler_step!(Q, Q_bar, Q_hx, Q_hz, grid, dt, config.gamma, config.g, sponge.K)
        if enable_viscosity
            MAGNUSP.diffusion_step!(Q, Q_bar, grid, dt, config, diff_solver)
        end
        t += dt; step_count += 1

        if t >= next_print
            ke = compute_ke()
            ke_peak = max(ke_peak, ke)
            @printf("t = %6.0f s:  KE = %.4e   peak = %.4e\n", t, ke, ke_peak)
            next_print += 1000.0
        end
    end

    println()
    @printf("Final KE_peak = %.4e\n", ke_peak)
end

println("=== Test 1: Viscosity ON (altitude-dependent nu(z)) ===")
run_k0_test(enable_viscosity=true)
println()
println("=== Test 2: Viscosity OFF (reference) ===")
run_k0_test(enable_viscosity=false)
println()
println("If viscosity-ON KE matches viscosity-OFF (no growth after t~3000s), the fix works.")
