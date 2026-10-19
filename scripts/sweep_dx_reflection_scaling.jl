# =============================================================================
# Grid Refinement Sweep: Reflection Scaling with Δx=Δz
#
# Purpose:
#   Test whether reflection coefficient scales with grid spacing, which would
#   indicate strong numerical diffusion in the advection scheme.
#
# Main comparison:
#   1) Well-balanced solver   (current MAGNUSP.euler_step!)
#   2) Standard LW solver     (total fields in z-flux and source)
#
# Key metric (same as prior scripts):
#   r(%) = sqrt( KE(t=5000) / KE_peak ) * 100
#
# Expected interpretation:
#   - Standard LW: strong Δx sensitivity (numerical diffusion dominated)
#   - Well-balanced: weak Δx sensitivity (less artificial dissipation)
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using MAGNUSP
using MAGNUSP: phys_x, phys_z, SVector, primitive_from_conservative
using Printf
using Statistics

# --- Physical/domain constants ---
const LX = 80_000.0           # m
const LZ = 120_000.0          # m (includes top sponge region)
const T_ISO = 239.0           # K
const CFL = 0.8
const T_END = 8000.0          # s
const DIAG_DT = 10.0          # s

# Wave forcing (Srivastava Table 3 / Eq. 29 style)
const W_TILDE = 0.001
const T0_FORCING = 1200.0
const SIGMA_T = 600.0
const OMEGA = 0.01
const LAMBDA_X = 80e3
const KX = 2.0 * pi / LAMBDA_X

# Sponge settings used for the core scaling test
const SPONGE_KM = 20.0
const K0_MAIN = 0.0

# Optional third case to compare with sponge-enabled standard run
const RUN_STANDARD_WITH_SPONGE = true
const K0_STANDARD_SPONGE = 0.005

# Grid sweep (must exactly divide both LX and LZ)
# 1000m => nx=80,  nz=120
#  800m => nx=100, nz=150
#  500m => nx=160, nz=240
#  400m => nx=200, nz=300
#  250m => nx=320, nz=480 (more expensive)
const DX_VALUES = [1000.0, 800.0, 500.0, 400.0, 250.0]

# =============================================================================
# Standard (non-well-balanced) z-flux/source and Euler step
# =============================================================================

function compute_flux_z_standard(Q::SVector{4,Float64}, gamma::Float64)
    rho, u, w, P = primitive_from_conservative(Q, gamma)
    E = Q[4]
    return SVector(
        rho * w,
        rho * u * w,
        rho * w^2 + P,
        w * (E + P),
    )
end

function compute_source_standard(Q::SVector{4,Float64}, g::Float64)
    rho = Q[1]
    w = Q[3] / rho
    return SVector(0.0, 0.0, -rho * g, -rho * g * w)
end

function euler_step_standard!(Q::Array{Float64,3},
                              Q_bar::Array{Float64,3},
                              Q_half_x::Array{Float64,3},
                              Q_half_z::Array{Float64,3},
                              grid::Grid2D,
                              dt::Float64,
                              gamma::Float64,
                              g::Float64,
                              sponge_K::Matrix{Float64})
    dtdx = dt / grid.dx
    dtdz = dt / grid.dz
    ix = phys_x(grid)
    iz = phys_z(grid)

    # X-sweep (same as well-balanced)
    for j in iz
        for i in (first(ix)-1):(last(ix))
            Qi   = SVector(Q[1,i,j], Q[2,i,j], Q[3,i,j], Q[4,i,j])
            Qip1 = SVector(Q[1,i+1,j], Q[2,i+1,j], Q[3,i+1,j], Q[4,i+1,j])
            Fi   = MAGNUSP.compute_flux_x(Qi, gamma)
            Fip1 = MAGNUSP.compute_flux_x(Qip1, gamma)

            @inbounds for k in 1:4
                Q_half_x[k, i, j] = 0.5 * (Qi[k] + Qip1[k]) -
                                    0.5 * dtdx * (Fip1[k] - Fi[k])
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

    # Z-sweep (standard total-field)
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

# =============================================================================
# Experiment drivers
# =============================================================================

struct CaseResult
    dx::Float64
    nx::Int
    nz::Int
    label::String
    K0::Float64
    times::Vector{Float64}
    ke::Vector{Float64}
    ke_peak::Float64
    ke_5000::Float64
    r::Float64
    hydro_drift::Float64
    steps::Int
    crashed::Bool
end

function run_case(dx::Float64; use_standard::Bool, K0::Float64, label::String)
    nx = Int(round(LX / dx))
    nz = Int(round(LZ / dx))

    if abs(nx * dx - LX) > 1e-6 || abs(nz * dx - LZ) > 1e-6
        error("dx=$(dx) does not exactly divide domain lengths LX/LZ")
    end

    sponge_frac = SPONGE_KM / (LZ / 1e3)

    config = earth_config(
        nx = nx,
        nz = nz,
        dx = dx,
        dz = dx,
        CFL = CFL,
        T0 = T_ISO,
        t_end = T_END,
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

    @printf("  [%s] dx=%6.1f m (nx=%3d, nz=%3d), K0=%.4f\n", label, dx, nx, nz, K0)

    while t < config.t_end
        local dt::Float64
        try
            dt = MAGNUSP.compute_adaptive_dt(Q, grid, config.CFL, config.gamma, config.R_gas)
            if !isfinite(dt) || dt <= 0.0
                crashed = true
                @printf("    crashed at t=%.1f s (invalid dt=%.3e)\n", t, dt)
                break
            end
            dt = min(dt, config.t_end - t)

            MAGNUSP.apply_all_bcs!(Q, Q_bar, grid, config.gamma)
            MAGNUSP.apply_bottom_forcing!(Q, grid, config, t, wave_params, :gaussian)

            if use_standard
                euler_step_standard!(Q, Q_bar, Q_hx, Q_hz, grid, dt, config.gamma, config.g,
                                     sponge_layer.K)
            else
                MAGNUSP.euler_step!(Q, Q_bar, Q_hx, Q_hz, grid, dt, config.gamma, config.g,
                                    sponge_layer.K)
            end
        catch err
            crashed = true
            @printf("    crashed at t=%.1f s (%s)\n", t, sprint(showerror, err))
            break
        end

        t += dt
        steps += 1

        if steps % 200 == 0
            probe = Q[1, grid.ng + 2, grid.ng + 2]
            if isnan(probe) || isinf(probe)
                crashed = true
                @printf("    crashed at t=%.1f s (NaN/Inf detected)\n", t)
                break
            end
        end

        if t >= next_diag
            ke = 0.0
            for j in phys_z(grid), i in phys_x(grid)
                rho = Q[1, i, j]
                u = Q[2, i, j] / rho
                w = Q[3, i, j] / rho
                ke += 0.5 * Q_bar[1, i, j] * (u^2 + w^2) * dx * dx
            end

            push!(times, t)
            push!(ke_arr, ke)
            next_diag += DIAG_DT
        end
    end

    if crashed || isempty(ke_arr)
        return CaseResult(dx, nx, nz, label, K0, Float64[], Float64[], NaN, NaN, NaN, NaN, steps, true)
    end

    ke_peak = maximum(ke_arr)
    idx5000 = argmin(abs.(times .- 5000.0))
    ke_5000 = ke_arr[idx5000]
    r = ke_peak > 0 ? sqrt(ke_5000 / ke_peak) * 100.0 : NaN

    max_drift = 0.0
    for j in phys_z(grid), i in phys_x(grid)
        rho_bar = Q_bar[1, i, j]
        if rho_bar > 1e-20
            drift = abs(Q[1, i, j] - rho_bar) / rho_bar
            max_drift = max(max_drift, drift)
        end
    end

    @printf("    done: steps=%6d, r=%.3f%%, KE_peak=%.3e, drift=%.2e\n",
            steps, r, ke_peak, max_drift)

    return CaseResult(dx, nx, nz, label, K0, times, ke_arr, ke_peak, ke_5000, r, max_drift, steps, false)
end

function scaling_exponent(results::Vector{CaseResult})
    good = filter(r -> !r.crashed && isfinite(r.r) && r.r > 0.0, results)
    if length(good) < 2
        return NaN
    end
    xs = log10.([r.dx for r in good])
    ys = log10.([r.r for r in good])
    mx = mean(xs)
    my = mean(ys)
    denom = sum((x - mx)^2 for x in xs)
    if denom <= 0
        return NaN
    end
    slope = sum((x - mx) * (y - my) for (x, y) in zip(xs, ys)) / denom
    return slope
end

# =============================================================================
# Main sweep
# =============================================================================

println("=" ^ 88)
println("GRID REFINEMENT SWEEP: Reflection scaling with Δx=Δz")
println("=" ^ 88)

cs = MAGNUSP.sound_speed(1.4, 287.0, T_ISO)
N_bv = MAGNUSP.brunt_vaisala(1.4, 9.81, cs)
m_sq = KX^2 * (N_bv^2 / OMEGA^2 - 1.0)
lambda_z = m_sq > 0 ? (2.0 * pi / sqrt(m_sq)) : NaN
@printf("T=%.0f K, cs=%.1f m/s, N=%.5f rad/s, λz≈%.1f km\n",
        T_ISO, cs, N_bv, lambda_z / 1e3)
@printf("Sweep dx values [m]: %s\n", string(DX_VALUES))
@printf("Core test uses K0=%.4f, sponge=%.0f km\n\n", K0_MAIN, SPONGE_KM)

wb_results = CaseResult[]
std_results = CaseResult[]
std_sponge_results = CaseResult[]

for dx in DX_VALUES
    push!(wb_results, run_case(dx; use_standard = false, K0 = K0_MAIN, label = "WB"))
    push!(std_results, run_case(dx; use_standard = true,  K0 = K0_MAIN, label = "STD"))
    if RUN_STANDARD_WITH_SPONGE
        push!(std_sponge_results,
              run_case(dx; use_standard = true, K0 = K0_STANDARD_SPONGE, label = "STD+SP"))
    end
    println()
end

# Summary table
println("=" ^ 110)
println("SUMMARY TABLE")
println("=" ^ 110)
@printf("%-8s %-6s %-6s %-10s %-11s %-13s %-11s %-9s\n",
        "Case", "dx(m)", "nx", "nz", "r(%)", "KE_peak", "KE@5000", "drift")
println("-" ^ 110)

function print_rows(rows::Vector{CaseResult})
    for r in rows
        if r.crashed
            @printf("%-8s %-6.0f %-6d %-10d %-11s %-13s %-11s %-9s\n",
                    r.label, r.dx, r.nx, r.nz, "CRASH", "-", "-", "-")
        else
            @printf("%-8s %-6.0f %-6d %-10d %-11.3f %-13.3e %-11.3e %-9.2e\n",
                    r.label, r.dx, r.nx, r.nz, r.r, r.ke_peak, r.ke_5000, r.hydro_drift)
        end
    end
end

print_rows(wb_results)
print_rows(std_results)
if RUN_STANDARD_WITH_SPONGE
    print_rows(std_sponge_results)
end

# Scaling fits: r ~ dx^p
p_wb = scaling_exponent(wb_results)
p_std = scaling_exponent(std_results)
p_std_sp = RUN_STANDARD_WITH_SPONGE ? scaling_exponent(std_sponge_results) : NaN

println()
println("=" ^ 110)
println("SCALING FIT (log10 r vs log10 dx)")
println("=" ^ 110)
@printf("WB      : p = %.3f\n", p_wb)
@printf("STD     : p = %.3f\n", p_std)
if RUN_STANDARD_WITH_SPONGE
    @printf("STD+SP  : p = %.3f\n", p_std_sp)
end

# Numerical diffusion estimate for standard LW at each dx
println()
println("=" ^ 110)
println("STANDARD LW NUMERICAL DIFFUSION ESTIMATE")
println("=" ^ 110)
@printf("%-6s %-14s %-14s %-12s\n", "dx(m)", "nu_num(m^2/s)", "gamma_num(1/s)", "r_pred(%)")
println("-" ^ 110)

m = sqrt(m_sq)
for dx in DX_VALUES
    nu_num = (cs * dx / 2.0) * (1.0 - CFL^2)
    gamma_num = nu_num * m^2
    ke_ratio_5000 = exp(-2.0 * gamma_num * 5000.0)
    r_pred = sqrt(ke_ratio_5000) * 100.0
    @printf("%-6.0f %-14.3e %-14.3e %-12.3f\n", dx, nu_num, gamma_num, r_pred)
end

println()
println("Done.")
println("Interpretation guide:")
println("  - If STD r decreases strongly as dx increases (coarser grid), numerical diffusion dominates.")
println("  - If WB r stays relatively high and weakly dx-dependent, well-balanced formulation suppresses")
println("    artificial dissipation and preserves wave energy.")
