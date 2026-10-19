# =============================================================================
# Well-Balanced vs Standard Lax-Wendroff: Numerical Diffusion Test
#
# Hypothesis: our well-balanced Euler formulation (P-P_bar in z-flux,
# (rho-rho_bar)*g in source) dramatically reduces numerical diffusion
# compared to the standard total-field formulation (P in z-flux, rho*g
# in source). The paper likely uses standard LW, which has inherent
# numerical viscosity ~ (c*dx/2)*(1-CFL^2) that damps waves with
# r ~ 7% for K0=0, matching their Table 4 result.
#
# This script runs K0=0 (no sponge damping) with both formulations and
# compares reflection coefficients. If the standard scheme gives r ~ 7%
# while well-balanced gives r ~ 91%, the hypothesis is confirmed.
#
# Runs:
#   1. Well-balanced (current solver)   K0=0, visc=OFF
#   2. Standard (total fields)          K0=0, visc=OFF
#   3. Standard (total fields)          K0=0.005, visc=OFF, 60km sponge
#
# Run 3 checks whether the standard scheme + sponge approaches r=0.55%.
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using MAGNUSP
using MAGNUSP: phys_x, phys_z, SVector, primitive_from_conservative
using Printf

# --- Fixed parameters (paper Tables 2 & 3) ---
const NX = 160
const NZ = 240
const DX = 500.0
const DZ = 500.0
const DOMAIN_KM = NZ * DZ / 1e3  # 120 km
const T_END = 8000.0
const DIAG_DT = 10.0
const T_ISO = 239.0

# Wave forcing
const W_TILDE = 0.001
const T0_FORCING = 1200.0
const SIGMA_T = 600.0
const OMEGA = 0.01
const LAMBDA_X = 80e3
const KX = 2.0 * pi / LAMBDA_X

# =============================================================================
# Standard (non-well-balanced) flux and source functions
# =============================================================================

"""
Standard z-flux: uses total pressure P (not P - P_bar).
G = [rho*w, rho*u*w, rho*w^2 + P, w*(E+P)]
"""
function compute_flux_z_standard(Q::SVector{4,Float64}, gamma::Float64)
    rho, u, w, P = primitive_from_conservative(Q, gamma)
    E = Q[4]
    return SVector(rho * w,
                   rho * u * w,
                   rho * w^2 + P,
                   w * (E + P))
end

"""
Standard source: uses total density rho (not rho - rho_bar).
S = [0, 0, -rho*g, -rho*g*w]
"""
function compute_source_standard(Q::SVector{4,Float64}, g::Float64)
    rho = Q[1]
    w   = Q[3] / rho
    return SVector(0.0, 0.0, -rho * g, -rho * g * w)
end

"""
    euler_step_standard!(Q, Q_bar, Q_half_x, Q_half_z, grid, dt, gamma, g, sponge_K)

Non-well-balanced Euler step. Identical structure to the well-balanced version
but uses total P in z-flux and total rho*g in source. The x-sweep is unchanged
(no gravity in x-direction). Sponge damping is still included.
"""
function euler_step_standard!(Q::Array{Float64,3},
                              Q_bar::Array{Float64,3},
                              Q_half_x::Array{Float64,3},
                              Q_half_z::Array{Float64,3},
                              grid::Grid2D, dt::Float64,
                              gamma::Float64, g::Float64,
                              sponge_K::Matrix{Float64})
    dx = grid.dx
    dz = grid.dz
    dtdx = dt / dx
    dtdz = dt / dz
    ix = phys_x(grid)
    iz = phys_z(grid)

    # =========================================================================
    # X-SWEEP: identical to well-balanced (no gravity in x)
    # =========================================================================
    for j in iz
        for i in (first(ix)-1):(last(ix))
            Qi   = SVector(Q[1,i,j], Q[2,i,j], Q[3,i,j], Q[4,i,j])
            Qip1 = SVector(Q[1,i+1,j], Q[2,i+1,j], Q[3,i+1,j], Q[4,i+1,j])
            Fi   = MAGNUSP.compute_flux_x(Qi, gamma)
            Fip1 = MAGNUSP.compute_flux_x(Qip1, gamma)

            for k in 1:4
                Q_half_x[k, i, j] = 0.5 * (Qi[k] + Qip1[k]) - 0.5 * dtdx * (Fip1[k] - Fi[k])
            end
        end
    end

    for j in iz
        for i in ix
            Qh_right = SVector(Q_half_x[1,i,j], Q_half_x[2,i,j],
                               Q_half_x[3,i,j], Q_half_x[4,i,j])
            Qh_left  = SVector(Q_half_x[1,i-1,j], Q_half_x[2,i-1,j],
                               Q_half_x[3,i-1,j], Q_half_x[4,i-1,j])
            Fh_right = MAGNUSP.compute_flux_x(Qh_right, gamma)
            Fh_left  = MAGNUSP.compute_flux_x(Qh_left, gamma)

            for k in 1:4
                Q[k, i, j] -= dtdx * (Fh_right[k] - Fh_left[k])
            end
        end
    end

    # =========================================================================
    # Z-SWEEP: STANDARD (non-well-balanced) — total P in flux, total rho*g
    # =========================================================================

    # Predictor
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

            for k in 1:4
                sp_j   = Kj   * (Q_bar[k, i, j]   - Qj[k])
                sp_jp1 = Kjp1 * (Q_bar[k, i, j+1] - Qjp1[k])
                Q_half_z[k, i, j] = 0.5 * (Qj[k] + Qjp1[k]) -
                                    0.5 * dtdz * (Gjp1[k] - Gj[k]) +
                                    0.25 * dt * (Sj[k] + Sjp1[k] + sp_j + sp_jp1)
            end
        end
    end

    # Corrector
    for j in iz
        for i in ix
            Qh_top = SVector(Q_half_z[1,i,j], Q_half_z[2,i,j],
                             Q_half_z[3,i,j], Q_half_z[4,i,j])
            Qh_bot = SVector(Q_half_z[1,i,j-1], Q_half_z[2,i,j-1],
                             Q_half_z[3,i,j-1], Q_half_z[4,i,j-1])
            Gh_top = compute_flux_z_standard(Qh_top, gamma)
            Gh_bot = compute_flux_z_standard(Qh_bot, gamma)
            Sh_top = compute_source_standard(Qh_top, g)
            Sh_bot = compute_source_standard(Qh_bot, g)

            K_top = 0.5 * (sponge_K[i, j]   + sponge_K[i, j+1])
            K_bot = 0.5 * (sponge_K[i, j-1] + sponge_K[i, j])

            for k in 1:4
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
# Run a single case
# =============================================================================

struct TestResult
    label::String
    times::Vector{Float64}
    ke::Vector{Float64}
    ke_peak::Float64
    r::Float64
    hydro_drift::Float64   # max |rho - rho_bar| / rho_bar at t_end
end

function run_case(; label::String, K0::Float64, sponge_km::Float64,
                    use_standard::Bool)
    sponge_frac = sponge_km / DOMAIN_KM

    config = earth_config(
        nx = NX, nz = NZ, dx = DX, dz = DZ, CFL = 0.8, T0 = T_ISO,
        t_end = T_END, output_interval = 999999,
        enable_viscosity = false, enable_cooling = false,
        sponge_frac = sponge_frac, K0_sponge = K0,
    )

    wave_params = MAGNUSP.WaveParams(
        amplitude = W_TILDE, omega = OMEGA, kx = KX,
        sigma_t = SIGMA_T, t0 = T0_FORCING,
    )

    grid = Grid2D(config)
    atm = isothermal_atmosphere(config, Float64.(grid.z_full))
    Q = MAGNUSP.init_state(grid, atm, config)
    Q_bar = copy(Q)
    sponge_layer = MAGNUSP.init_sponge(grid, config)
    Q_hx = similar(Q); Q_hz = similar(Q)

    times = Float64[]
    ke_arr = Float64[]

    t = 0.0
    step_count = 0
    next_diag = 0.0
    crashed = false

    scheme_str = use_standard ? "STANDARD" : "WELL-BALANCED"
    @printf("  [%s] K0=%.4f, sponge=%.0fkm ...\n", scheme_str, K0, sponge_km)

    while t < config.t_end
        dt = MAGNUSP.compute_adaptive_dt(Q, grid, config.CFL, config.gamma, config.R_gas)
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

        t += dt
        step_count += 1

        # Check for NaN / crash
        if step_count % 100 == 0
            test_val = Q[1, div(NX,2)+2, div(NZ,2)+2]
            if isnan(test_val) || isinf(test_val)
                @printf("  CRASH at step %d, t=%.1f\n", step_count, t)
                crashed = true
                break
            end
        end

        if t >= next_diag
            ke = 0.0
            for j in phys_z(grid), i in phys_x(grid)
                u = Q[2, i, j] / Q[1, i, j]
                w = Q[3, i, j] / Q[1, i, j]
                ke += 0.5 * Q_bar[1, i, j] * (u^2 + w^2) * DX * DZ
            end
            push!(times, t)
            push!(ke_arr, ke)
            next_diag += DIAG_DT

            # Print at 1000s intervals
            if abs(t - round(t / 1000.0) * 1000.0) < DIAG_DT
                @printf("    t=%5.0fs  KE=%.4e  steps=%d\n", t, ke, step_count)
            end
        end
    end

    if crashed || isempty(ke_arr)
        return TestResult(label, Float64[], Float64[], NaN, NaN, NaN)
    end

    # Reflection coefficient
    ke_peak = maximum(ke_arr)
    t5000_idx = argmin(abs.(times .- 5000.0))
    r = ke_peak > 0 ? sqrt(ke_arr[t5000_idx] / ke_peak) * 100.0 : NaN

    # Hydrostatic drift: max relative density error vs background
    max_drift = 0.0
    ng = grid.ng
    for j in phys_z(grid), i in phys_x(grid)
        rho_bar = Q_bar[1, i, j]
        if rho_bar > 1e-20
            drift = abs(Q[1, i, j] - rho_bar) / rho_bar
            max_drift = max(max_drift, drift)
        end
    end

    @printf("  Done: %d steps, KE_peak=%.4e, r=%.2f%%, hydro_drift=%.2e\n\n",
            step_count, ke_peak, r, max_drift)

    return TestResult(label, times, ke_arr, ke_peak, r, max_drift)
end

# =============================================================================
# Main
# =============================================================================
println("=" ^ 80)
println("WELL-BALANCED vs STANDARD LAX-WENDROFF: Numerical Diffusion Test")
println("=" ^ 80)

cs = MAGNUSP.sound_speed(1.4, 287.0, T_ISO)
N_bv = MAGNUSP.brunt_vaisala(1.4, 9.81, cs)
m = sqrt(KX^2 * (N_bv^2 / OMEGA^2 - 1.0))
lz = 2 * pi / m / 1e3
@printf("T=%.0fK, cs=%.1f m/s, N=%.5f rad/s, lambda_z=%.1f km\n", T_ISO, cs, N_bv, lz)

# Back-of-envelope numerical diffusion estimate
CFL_est = 0.8
nu_num = (cs * DX / 2.0) * (1.0 - CFL_est^2)
rate_num = nu_num * m^2
ke_ratio_5000 = exp(-2.0 * rate_num * 5000.0)
r_predicted = sqrt(ke_ratio_5000) * 100.0
@printf("LW numerical viscosity estimate: nu_num=%.0f m^2/s, damping rate=%.4e /s\n",
        nu_num, rate_num)
@printf("Predicted r(t=5000) from numerical diffusion alone: %.1f%%\n\n", r_predicted)

# --- Run cases ---
results = TestResult[]

push!(results, run_case(label = "Well-balanced, K0=0",
                        K0 = 0.0, sponge_km = 20.0, use_standard = false))

push!(results, run_case(label = "Standard, K0=0",
                        K0 = 0.0, sponge_km = 20.0, use_standard = true))

push!(results, run_case(label = "Standard, K0=0.005, 60km sponge",
                        K0 = 0.005, sponge_km = 60.0, use_standard = true))

# --- Summary ---
println("=" ^ 80)
println("SUMMARY")
println("=" ^ 80)
@printf("%-40s  %-10s  %-10s  %-12s\n", "Case", "r (%)", "KE_peak", "Hydro drift")
println("-" ^ 80)
for r in results
    if isnan(r.r)
        @printf("%-40s  CRASHED\n", r.label)
    else
        @printf("%-40s  %-10.2f  %-10.3e  %-12.2e\n",
                r.label, r.r, r.ke_peak, r.hydro_drift)
    end
end

println()
println("Paper targets: K0=0 -> r=7.6%, K0=0.005 -> r=0.55%")
println()
println("INTERPRETATION:")
println("  If standard K0=0 gives r ~ 7% while well-balanced gives r ~ 91%,")
println("  the paper's K0=0 result is dominated by LW numerical diffusion.")
println("  Our well-balanced scheme is more accurate (less artificial dissipation).")
println()
println("  If standard K0=0 CRASHES (hydrostatic drift -> negative P), that")
println("  confirms well-balancing is necessary for stability in deep atmospheres,")
println("  and the paper may use a different stabilization mechanism.")
