# =============================================================================
# WB-only Boundary + Sponge Sweep
#
# Purpose:
#   Diagnose near-total reflection in the well-balanced (WB) solver by jointly
#   sweeping top boundary condition and sponge design.
#
# Sweep dimensions:
#   - Top BC type:      [:eq19, :characteristic]
#   - Sponge thickness: [1.0, 1.5, 2.0] * lambda_z
#   - Sponge K0:        [0.0, 0.003, 0.005, 0.01]
#
# Diagnostics:
#   - r(%) = sqrt(KE_total(t=5000) / KE_total_peak) * 100
#   - t_top_peak: time of peak KE in sponge region
#   - refl_idx: post/pre ratio of low-domain KE around top interaction
#               (>1 indicates stronger downward returned energy)
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using MAGNUSP
using MAGNUSP: phys_x, phys_z
using Printf

# --- Base setup (aligned with prior validation scripts) ---
const LX = 80_000.0
const LZ = 120_000.0
const DX = 500.0
const DZ = 500.0
const NX = Int(round(LX / DX))
const NZ = Int(round(LZ / DZ))

const T_ISO = 239.0
const CFL = 0.8
const T_END = 8000.0
const DIAG_DT = 10.0

# Forcing
const W_TILDE = 0.001
const T0_FORCING = 1200.0
const SIGMA_T = 600.0
const OMEGA = 0.01
const LAMBDA_X = 80e3
const KX = 2.0 * pi / LAMBDA_X

# Sweep
const BC_TYPES = [:eq19, :characteristic]
const THICKNESS_LZ = [1.0, 1.5, 2.0]
const K0_VALUES = [0.0, 0.003, 0.005, 0.01]

struct SweepResult
    bc_type::Symbol
    sponge_km::Float64
    sponge_lz::Float64
    K0::Float64
    steps::Int
    crashed::Bool
    r::Float64
    t_top_peak::Float64
    refl_idx::Float64
    ke_peak::Float64
    ke_5000::Float64
end

function mean_window(t::Vector{Float64}, y::Vector{Float64}, t1::Float64, t2::Float64)
    idx = findall((t .>= t1) .& (t .<= t2))
    isempty(idx) && return NaN
    return sum(y[idx]) / length(idx)
end

function run_case(; bc_type::Symbol, sponge_km::Float64, sponge_lz::Float64, K0::Float64)
    sponge_frac = sponge_km / (LZ / 1e3)

    config = earth_config(
        nx = NX,
        nz = NZ,
        dx = DX,
        dz = DZ,
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

    z_phys = [grid.z_full[j] for j in phys_z(grid)]
    z_sponge_start = (LZ / 1e3 - sponge_km) * 1e3
    top_idx_local = findall(z_phys .>= z_sponge_start)
    low_idx_local = findall(z_phys .<= 80_000.0)

    times = Float64[]
    ke_total = Float64[]
    ke_top = Float64[]
    ke_low = Float64[]

    t = 0.0
    steps = 0
    next_diag = 0.0
    crashed = false

    @printf("[bc=%s, sponge=%.1f km (%.1f lz), K0=%.4f] ...\n",
            String(bc_type), sponge_km, sponge_lz, K0)

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

            MAGNUSP.apply_all_bcs!(Q, Q_bar, grid, config.gamma;
                                   bc_type = bc_type,
                                   R_gas = config.R_gas,
                                   T0 = config.T0)
            MAGNUSP.apply_bottom_forcing!(Q, grid, config, t, wave_params, :gaussian)
            MAGNUSP.euler_step!(Q, Q_bar, Q_hx, Q_hz, grid, dt, config.gamma, config.g,
                                sponge_layer.K)
        catch err
            crashed = true
            @printf("  crashed at t=%.1f s (%s)\n", t, sprint(showerror, err))
            break
        end

        t += dt
        steps += 1

        if t >= next_diag
            ke_t = 0.0
            ke_t_top = 0.0
            ke_t_low = 0.0

            for (j_local, j) in enumerate(phys_z(grid))
                row_ke = 0.0
                for i in phys_x(grid)
                    rho = Q[1, i, j]
                    u = Q[2, i, j] / rho
                    w = Q[3, i, j] / rho
                    row_ke += 0.5 * Q_bar[1, i, j] * (u^2 + w^2) * DX * DZ
                end
                ke_t += row_ke
                if j_local in top_idx_local
                    ke_t_top += row_ke
                end
                if j_local in low_idx_local
                    ke_t_low += row_ke
                end
            end

            push!(times, t)
            push!(ke_total, ke_t)
            push!(ke_top, ke_t_top)
            push!(ke_low, ke_t_low)
            next_diag += DIAG_DT
        end
    end

    if crashed || isempty(ke_total)
        return SweepResult(bc_type, sponge_km, sponge_lz, K0, steps, true,
                           NaN, NaN, NaN, NaN, NaN)
    end

    i_top = argmax(ke_top)
    t_top_peak = times[i_top]

    ke_peak = maximum(ke_total)
    i5000 = argmin(abs.(times .- 5000.0))
    ke_5000 = ke_total[i5000]
    r = ke_peak > 0 ? sqrt(ke_5000 / ke_peak) * 100.0 : NaN

    pre = mean_window(times, ke_low, max(0.0, t_top_peak - 1500.0), max(0.0, t_top_peak - 500.0))
    post = mean_window(times, ke_low, t_top_peak + 500.0, min(T_END, t_top_peak + 1500.0))
    refl_idx = (isfinite(pre) && pre > 0 && isfinite(post)) ? (post / pre) : NaN

    @printf("  done: steps=%d, r=%.2f%%, t_top=%.1fs, refl_idx=%.3f\n",
            steps, r, t_top_peak, refl_idx)

    return SweepResult(bc_type, sponge_km, sponge_lz, K0, steps, false,
                       r, t_top_peak, refl_idx, ke_peak, ke_5000)
end

# --- Derived wave scale for thickness conversion ---
cs = MAGNUSP.sound_speed(1.4, 287.0, T_ISO)
N_bv = MAGNUSP.brunt_vaisala(1.4, 9.81, cs)
m_sq = KX^2 * (N_bv^2 / OMEGA^2 - 1.0)
lambda_z = m_sq > 0 ? 2.0 * pi / sqrt(m_sq) : NaN
lambda_z_km = lambda_z / 1e3

println("=" ^ 96)
println("WB BOUNDARY + SPONGE SWEEP")
println("=" ^ 96)
@printf("T=%.0fK, dx=dz=%.0f m, lambda_z=%.1f km, domain=%.0f x %.0f km\n\n",
        T_ISO, DX, lambda_z_km, LX/1e3, LZ/1e3)

results = SweepResult[]
for bc in BC_TYPES
    for th in THICKNESS_LZ
        sponge_km = min(th * lambda_z_km, LZ / 1e3 - 5.0)  # keep at least 5 km non-sponge region
        for K0 in K0_VALUES
            push!(results, run_case(bc_type = bc, sponge_km = sponge_km, sponge_lz = th, K0 = K0))
        end
    end
end

println()
println("=" ^ 120)
println("SUMMARY")
println("=" ^ 120)
@printf("%-14s %-14s %-8s %-8s %-10s %-10s %-12s %-12s\n",
        "bc", "sponge(km/lz)", "K0", "crash", "r(%)", "t_top(s)", "refl_idx", "KE_peak")
println("-" ^ 120)

for r in results
    sk = @sprintf("%.1f/%.1f", r.sponge_km, r.sponge_lz)
    if r.crashed
        @printf("%-14s %-14s %-8.4f %-8s %-10s %-10s %-12s %-12s\n",
                String(r.bc_type), sk, r.K0, "YES", "-", "-", "-", "-")
    else
        @printf("%-14s %-14s %-8.4f %-8s %-10.2f %-10.1f %-12.3f %-12.3e\n",
                String(r.bc_type), sk, r.K0, "NO", r.r, r.t_top_peak, r.refl_idx, r.ke_peak)
    end
end

println()
println("Best (lowest r) per boundary type:")
for bc in BC_TYPES
    sub = filter(x -> (!x.crashed) && x.bc_type == bc && isfinite(x.r), results)
    if isempty(sub)
        println("  $(bc): no successful runs")
        continue
    end
    best = sub[argmin([x.r for x in sub])]
    @printf("  %-14s r=%.2f%%, sponge=%.1f km (%.1f lz), K0=%.4f, refl_idx=%.3f\n",
            String(bc), best.r, best.sponge_km, best.sponge_lz, best.K0, best.refl_idx)
end

println("\nDone.")
