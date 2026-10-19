# =============================================================================
# Sanity Check: Standard LW at dx=dz=125 m
#
# Focused follow-up run for the stable STD trend from the grid sweep.
# Runs only dx=125 m cases:
#   1) STD, K0=0.0
#   2) STD, K0=0.005 (same 20 km sponge as baseline setup)
#
# Reflection metric:
#   r(%) = sqrt(KE(t=5000) / KE_peak) * 100
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using MAGNUSP
using MAGNUSP: phys_x, phys_z, SVector, primitive_from_conservative
using Printf

const LX = 80_000.0
const LZ = 120_000.0
const DX = 125.0
const DZ = 125.0
const NX = Int(round(LX / DX))   # 640
const NZ = Int(round(LZ / DZ))   # 960

const T_ISO = 239.0
const CFL = 0.8
const T_END = 8000.0
const DIAG_DT = 10.0
const SPONGE_KM = 20.0

const W_TILDE = 0.001
const T0_FORCING = 1200.0
const SIGMA_T = 600.0
const OMEGA = 0.01
const LAMBDA_X = 80e3
const KX = 2.0 * pi / LAMBDA_X

function compute_flux_z_standard(Q::SVector{4,Float64}, gamma::Float64)
    rho, u, w, P = primitive_from_conservative(Q, gamma)
    E = Q[4]
    return SVector(rho * w, rho * u * w, rho * w^2 + P, w * (E + P))
end

function compute_source_standard(Q::SVector{4,Float64}, g::Float64)
    rho = Q[1]
    w = Q[3] / rho
    return SVector(0.0, 0.0, -rho * g, -rho * g * w)
end

function euler_step_standard!(Q::Array{Float64,3}, Q_bar::Array{Float64,3},
                              Q_half_x::Array{Float64,3}, Q_half_z::Array{Float64,3},
                              grid::Grid2D, dt::Float64, gamma::Float64,
                              g::Float64, sponge_K::Matrix{Float64})
    dtdx = dt / grid.dx
    dtdz = dt / grid.dz
    ix = phys_x(grid)
    iz = phys_z(grid)

    for j in iz
        for i in (first(ix)-1):(last(ix))
            Qi   = SVector(Q[1,i,j], Q[2,i,j], Q[3,i,j], Q[4,i,j])
            Qip1 = SVector(Q[1,i+1,j], Q[2,i+1,j], Q[3,i+1,j], Q[4,i+1,j])
            Fi   = MAGNUSP.compute_flux_x(Qi, gamma)
            Fip1 = MAGNUSP.compute_flux_x(Qip1, gamma)
            @inbounds for k in 1:4
                Q_half_x[k, i, j] = 0.5 * (Qi[k] + Qip1[k]) - 0.5 * dtdx * (Fip1[k] - Fi[k])
            end
        end
    end

    for j in iz
        for i in ix
            Qh_right = SVector(Q_half_x[1,i,j], Q_half_x[2,i,j], Q_half_x[3,i,j], Q_half_x[4,i,j])
            Qh_left  = SVector(Q_half_x[1,i-1,j], Q_half_x[2,i-1,j], Q_half_x[3,i-1,j], Q_half_x[4,i-1,j])
            Fh_right = MAGNUSP.compute_flux_x(Qh_right, gamma)
            Fh_left  = MAGNUSP.compute_flux_x(Qh_left, gamma)
            @inbounds for k in 1:4
                Q[k, i, j] -= dtdx * (Fh_right[k] - Fh_left[k])
            end
        end
    end

    for j in (first(iz)-1):(last(iz))
        for i in ix
            Qj   = SVector(Q[1,i,j], Q[2,i,j], Q[3,i,j], Q[4,i,j])
            Qjp1 = SVector(Q[1,i,j+1], Q[2,i,j+1], Q[3,i,j+1], Q[4,i,j+1])
            Gj   = compute_flux_z_standard(Qj, gamma)
            Gjp1 = compute_flux_z_standard(Qjp1, gamma)
            Sj   = compute_source_standard(Qj, g)
            Sjp1 = compute_source_standard(Qjp1, g)
            Kj   = sponge_K[i, j]
            Kjp1 = sponge_K[i, j+1]
            @inbounds for k in 1:4
                sp_j   = Kj   * (Q_bar[k, i, j]   - Qj[k])
                sp_jp1 = Kjp1 * (Q_bar[k, i, j+1] - Qjp1[k])
                Q_half_z[k, i, j] = 0.5 * (Qj[k] + Qjp1[k]) -
                                    0.5 * dtdz * (Gjp1[k] - Gj[k]) +
                                    0.25 * dt * (Sj[k] + Sjp1[k] + sp_j + sp_jp1)
            end
        end
    end

    for j in iz
        for i in ix
            Qh_top = SVector(Q_half_z[1,i,j], Q_half_z[2,i,j], Q_half_z[3,i,j], Q_half_z[4,i,j])
            Qh_bot = SVector(Q_half_z[1,i,j-1], Q_half_z[2,i,j-1], Q_half_z[3,i,j-1], Q_half_z[4,i,j-1])
            Gh_top = compute_flux_z_standard(Qh_top, gamma)
            Gh_bot = compute_flux_z_standard(Qh_bot, gamma)
            Sh_top = compute_source_standard(Qh_top, g)
            Sh_bot = compute_source_standard(Qh_bot, g)
            K_top = 0.5 * (sponge_K[i, j]   + sponge_K[i, j+1])
            K_bot = 0.5 * (sponge_K[i, j-1] + sponge_K[i, j])
            @inbounds for k in 1:4
                Qbar_half_top = 0.5 * (Q_bar[k, i, j]   + Q_bar[k, i, j+1])
                Qbar_half_bot = 0.5 * (Q_bar[k, i, j-1] + Q_bar[k, i, j])
                sp_top = K_top * (Qbar_half_top - Qh_top[k])
                sp_bot = K_bot * (Qbar_half_bot - Qh_bot[k])
                Q[k, i, j] -= dtdz * (Gh_top[k] - Gh_bot[k])
                Q[k, i, j] += dt * 0.5 * (Sh_top[k] + Sh_bot[k] + sp_top + sp_bot)
            end
        end
    end

    return nothing
end

function run_std_case(K0::Float64)
    sponge_frac = SPONGE_KM / (LZ / 1e3)
    config = earth_config(
        nx = NX, nz = NZ, dx = DX, dz = DZ,
        CFL = CFL, T0 = T_ISO, t_end = T_END,
        output_interval = 999999,
        enable_viscosity = false,
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
    Q_hx = similar(Q)
    Q_hz = similar(Q)

    times = Float64[]
    ke_arr = Float64[]

    t = 0.0
    steps = 0
    next_diag = 0.0
    crashed = false

    @printf("\n[STD dx=125m] K0=%.4f\n", K0)

    while t < config.t_end
        local dt::Float64
        try
            dt = MAGNUSP.compute_adaptive_dt(Q, grid, config.CFL, config.gamma, config.R_gas)
            if !isfinite(dt) || dt <= 0.0
                crashed = true
                @printf("  crashed at t=%.1f s (invalid dt=%.3e)\n", t, dt)
                break
            end
            dt = min(dt, config.t_end - t)

            MAGNUSP.apply_all_bcs!(Q, Q_bar, grid, config.gamma)
            MAGNUSP.apply_bottom_forcing!(Q, grid, config, t, wave_params, :gaussian)
            euler_step_standard!(Q, Q_bar, Q_hx, Q_hz, grid, dt, config.gamma, config.g,
                                 sponge_layer.K)
        catch err
            crashed = true
            @printf("  crashed at t=%.1f s (%s)\n", t, sprint(showerror, err))
            break
        end

        t += dt
        steps += 1

        if t >= next_diag
            ke = 0.0
            for j in phys_z(grid), i in phys_x(grid)
                rho = Q[1, i, j]
                u = Q[2, i, j] / rho
                w = Q[3, i, j] / rho
                ke += 0.5 * Q_bar[1, i, j] * (u^2 + w^2) * DX * DZ
            end
            push!(times, t)
            push!(ke_arr, ke)
            next_diag += DIAG_DT

            if steps % 2500 == 0
                @printf("  t=%6.1f s, KE=%.3e, steps=%d\n", t, ke, steps)
            end
        end
    end

    if crashed || isempty(ke_arr)
        println("  Result: CRASH")
        return
    end

    ke_peak = maximum(ke_arr)
    idx5000 = argmin(abs.(times .- 5000.0))
    ke_5000 = ke_arr[idx5000]
    r = sqrt(ke_5000 / ke_peak) * 100.0

    max_drift = 0.0
    for j in phys_z(grid), i in phys_x(grid)
        rho_bar = Q_bar[1, i, j]
        if rho_bar > 1e-20
            drift = abs(Q[1, i, j] - rho_bar) / rho_bar
            max_drift = max(max_drift, drift)
        end
    end

    @printf("  Done: steps=%d, r=%.3f%%, KE_peak=%.3e, KE@5000=%.3e, drift=%.2e\n",
            steps, r, ke_peak, ke_5000, max_drift)
end

println("=" ^ 76)
println("STD dx=125m SANITY CHECK")
println("=" ^ 76)
@printf("Domain: %.0f x %.0f km, grid: nx=%d nz=%d\n", LX/1e3, LZ/1e3, NX, NZ)

run_std_case(0.0)
run_std_case(0.005)

println("\nDone.")
