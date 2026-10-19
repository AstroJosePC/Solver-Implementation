# Ultra-low frequency sponge validation
# Focus on omega = 0.1 mHz (6.28e-4 rad/s)
# lambda_x = 200 km (sigma = 50 km)

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using MAGNUSP
using MAGNUSP: SVector, phys_x, phys_z
using Printf
using NCDatasets

const NX = 200
const NZ = 480      # 240 km / 0.5 km (deeper domain)
const DX = 1000.0   # 1 km
const DZ = 500.0    # 0.5 km
const T_ISOTHERMAL = 210.0  # Mars-like
const T_END = 50000.0       # 14 hours
const DIAG_DT = 100.0

const W_TILDE = 0.001
const T0_FORCING = 10000.0
const SIGMA_T = 5000.0
const OMEGA = 0.000628      # 0.1 mHz
const LAMBDA_X = 200e3      # 200 km
const KX = 2.0 * pi / LAMBDA_X

function run_case(K0)
    config = mars_config(
        nx = NX,
        nz = NZ,
        dx = DX,
        dz = DZ,
        CFL = 0.8,
        T0 = T_ISOTHERMAL,
        t_end = T_END,
        output_interval = 999999,
        enable_viscosity = true,
        enable_cooling = false,
        sponge_frac = 40.0 / (NZ * DZ / 1000.0), # 40 km sponge
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

    times = Float64[]
    ke_values = Float64[]
    t = 0.0
    next_diag = 0.0

    while t < config.t_end
        dt = MAGNUSP.compute_adaptive_dt(Q, grid, config.CFL, config.gamma, config.R_gas)
        dt = min(dt, config.t_end - t)
        MAGNUSP.apply_all_bcs!(Q, Q_bar, grid, config.gamma)
        MAGNUSP.apply_bottom_forcing!(Q, grid, config, t, wave_params, :gaussian)
        MAGNUSP.euler_step!(Q, Q_bar, Q_half_x, Q_half_z, grid, dt, config.gamma, config.g, sponge_layer.K)
        if config.enable_viscosity
            MAGNUSP.diffusion_step!(Q, Q_bar, grid, dt, config, diff_solver)
        end
        t += dt
        if t >= next_diag
            ke = 0.0
            for j in phys_z(grid), i in phys_x(grid)
                rho_bg = Q_bar[1, i, j]
                u_pert = Q[2, i, j] / Q[1, i, j]
                w_pert = Q[3, i, j] / Q[1, i, j]
                ke += 0.5 * rho_bg * (u_pert^2 + w_pert^2) * DX * DZ
            end
            push!(times, t)
            push!(ke_values, ke)
            next_diag += DIAG_DT
        end
    end
    return times, ke_values
end

K0_test = 0.05
times, ke = run_case(K0_test)
println("Peak KE: ", maximum(ke))
