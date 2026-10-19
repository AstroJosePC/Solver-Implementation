# =============================================================================
# Long WB Run (eq19 only): Late Reflection + Below-80 km Fidelity Check
#
# Purpose:
#   1) Extend simulation beyond 8000 s to detect late re-entry reflections.
#   2) Verify wave amplitude/phase fidelity below 80 km is not over-damped.
#
# Configuration (fixed):
#   - Well-balanced solver only
#   - Top BC: eq19 (characteristic excluded)
#   - Sponge: 2.0 * lambda_z, K0=0.0100 (best-r case from sweep)
#
# Diagnostics:
#   - r_5000, r_8000, r_10000, r_12000
#   - KE below 80 km (pre/post/top-late windows)
#   - Mid-x w and density-scaled w at z=[20,40,60,80] km
#   - Period consistency via lag-correlation estimate around forcing period
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using MAGNUSP
using MAGNUSP: phys_x, phys_z
using Printf
using Statistics

# --- Base setup ---
const LX = 80_000.0
const LZ = 120_000.0
const DX = 500.0
const DZ = 500.0
const NX = Int(round(LX / DX))
const NZ = Int(round(LZ / DZ))

const T_ISO = 239.0
const CFL = 0.8
const T_END = 12_000.0
const DIAG_DT = 10.0

# Forcing
const W_TILDE = 0.001
const T0_FORCING = 1200.0
const SIGMA_T = 600.0
const OMEGA = 0.01
const LAMBDA_X = 80e3
const KX = 2.0 * pi / LAMBDA_X
const T_FORCE = 2.0 * pi / OMEGA

# Sponge config from WB sweep best-r case
const SPONGE_LZ = 2.0
const K0_SPONGE = 0.0100

# Altitudes for fidelity probes
const Z_PROBE_KM = [20.0, 40.0, 60.0, 80.0]

function nearest_index(vals::Vector{Float64}, x::Float64)
    return argmin(abs.(vals .- x))
end

function mean_window(t::Vector{Float64}, y::Vector{Float64}, t1::Float64, t2::Float64)
    idx = findall((t .>= t1) .& (t .<= t2))
    isempty(idx) && return NaN
    return mean(y[idx])
end

function rms_window(t::Vector{Float64}, y::Vector{Float64}, t1::Float64, t2::Float64)
    idx = findall((t .>= t1) .& (t .<= t2))
    isempty(idx) && return NaN
    return sqrt(mean(y[idx] .^ 2))
end

function est_period_lagcorr(t::Vector{Float64}, y::Vector{Float64};
                            t1::Float64, t2::Float64,
                            pmin::Float64 = 300.0,
                            pmax::Float64 = 900.0)
    idx = findall((t .>= t1) .& (t .<= t2))
    if length(idx) < 20
        return NaN
    end

    yw = y[idx]
    yw = yw .- mean(yw)
    dt = median(diff(t[idx]))
    lmin = max(1, Int(round(pmin / dt)))
    lmax = min(length(yw) - 2, Int(round(pmax / dt)))
    if lmax <= lmin
        return NaN
    end

    best_lag = lmin
    best_corr = -Inf
    for lag in lmin:lmax
        a = @view yw[1:end-lag]
        b = @view yw[1+lag:end]
        denom = sqrt(sum(a .^ 2) * sum(b .^ 2))
        if denom <= 0
            continue
        end
        c = sum(a .* b) / denom
        if c > best_corr
            best_corr = c
            best_lag = lag
        end
    end

    return best_lag * dt
end

let
    # --- Derived wave scale ---
    cs = MAGNUSP.sound_speed(1.4, 287.0, T_ISO)
    N_bv = MAGNUSP.brunt_vaisala(1.4, 9.81, cs)
    m_sq = KX^2 * (N_bv^2 / OMEGA^2 - 1.0)
    lambda_z = m_sq > 0 ? 2.0 * pi / sqrt(m_sq) : NaN
    lambda_z_km = lambda_z / 1e3
    sponge_km = SPONGE_LZ * lambda_z_km
    sponge_km = min(sponge_km, LZ / 1e3 - 5.0)
    sponge_frac = sponge_km / (LZ / 1e3)

    println("=" ^ 96)
    println("LONG WB EQ19 FIDELITY RUN")
    println("=" ^ 96)
    @printf("T=%.0fK, dx=dz=%.0fm, t_end=%.0fs, forcing period=%.1fs\n", T_ISO, DX, T_END, T_FORCE)
    @printf("Top BC=eq19, sponge=%.1f km (%.1f lz), K0=%.4f\n\n", sponge_km, sponge_km / lambda_z_km, K0_SPONGE)

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
        K0_sponge = K0_SPONGE,
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

    # Probe mapping
    z_phys = [grid.z_full[j] for j in phys_z(grid)]
    j_phys = collect(phys_z(grid))
    probe_j_local = Int[]
    probe_j = Int[]
    for zk in Z_PROBE_KM
        jl = nearest_index(z_phys, zk * 1e3)
        push!(probe_j_local, jl)
        push!(probe_j, j_phys[jl])
    end

    # Region indices
    top_idx_local = findall(z_phys .>= (LZ / 1e3 - sponge_km) * 1e3)
    low80_idx_local = findall(z_phys .<= 80_000.0)

    # Diagnostics arrays
    times = Float64[]
    ke_total = Float64[]
    ke_top = Float64[]
    ke_low80 = Float64[]

    nprobe = length(probe_j)
    w_probe = [Float64[] for _ in 1:nprobe]
    wscaled_probe = [Float64[] for _ in 1:nprobe]

    i_mid = grid.ng + div(grid.nx, 2)
    rho_surface = Q_bar[1, i_mid, grid.ng + 1]

    t = 0.0
    steps = 0
    next_diag = 0.0
    crashed = false

    while t < config.t_end
    local dt::Float64
    try
        dt = MAGNUSP.compute_adaptive_dt(Q, grid, config.CFL, config.gamma, config.R_gas)
        if !isfinite(dt) || dt <= 0.0
            crashed = true
            @printf("CRASH at t=%.1f s (invalid dt=%.3e)\n", t, dt)
            break
        end
        dt = min(dt, config.t_end - t)

        MAGNUSP.apply_all_bcs!(Q, Q_bar, grid, config.gamma;
                               bc_type = :eq19,
                               R_gas = config.R_gas,
                               T0 = config.T0)
        MAGNUSP.apply_bottom_forcing!(Q, grid, config, t, wave_params, :gaussian)
        MAGNUSP.euler_step!(Q, Q_bar, Q_hx, Q_hz, grid, dt, config.gamma, config.g,
                            sponge_layer.K)
    catch err
        crashed = true
        @printf("CRASH at t=%.1f s (%s)\n", t, sprint(showerror, err))
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
            if j_local in low80_idx_local
                ke_t_low += row_ke
            end
        end

        push!(times, t)
        push!(ke_total, ke_t)
        push!(ke_top, ke_t_top)
        push!(ke_low80, ke_t_low)

        # Mid-x probes
        for p in 1:nprobe
            j = probe_j[p]
            rho = Q[1, i_mid, j]
            w = Q[3, i_mid, j] / rho
            rho_bg = Q_bar[1, i_mid, j]
            w_sc = w * sqrt(rho_bg / rho_surface)
            push!(w_probe[p], w)
            push!(wscaled_probe[p], w_sc)
        end

        next_diag += DIAG_DT

        if steps % 2500 == 0
            @printf("  t=%7.1f s, KE_total=%.3e\n", t, ke_t)
        end
    end
    end

    if crashed || isempty(ke_total)
        error("Run failed; no diagnostics produced.")
    end

    # --- Summary metrics ---
    ke_peak = maximum(ke_total)
    function r_at(tq)
        iq = argmin(abs.(times .- tq))
        return sqrt(ke_total[iq] / ke_peak) * 100.0
    end

    r5000 = r_at(5000.0)
    r8000 = r_at(8000.0)
    r10000 = r_at(10000.0)
    r12000 = r_at(12000.0)

    i_top = argmax(ke_top)
    t_top = times[i_top]

    pre_low = mean_window(times, ke_low80, max(0.0, t_top - 1500.0), max(0.0, t_top - 500.0))
    post_low = mean_window(times, ke_low80, t_top + 500.0, min(T_END, t_top + 1500.0))
    late_low = mean_window(times, ke_low80, 10000.0, 12000.0)

    refl_idx = (isfinite(pre_low) && pre_low > 0 && isfinite(post_low)) ? (post_low / pre_low) : NaN
    late_idx = (isfinite(pre_low) && pre_low > 0 && isfinite(late_low)) ? (late_low / pre_low) : NaN

    println()
    println("--- Reflection / late re-entry diagnostics ---")
    @printf("steps=%d, KE_peak=%.3e, t_top_peak=%.1f s\n", steps, ke_peak, t_top)
    @printf("r@5000=%.2f%%, r@8000=%.2f%%, r@10000=%.2f%%, r@12000=%.2f%%\n", r5000, r8000, r10000, r12000)
    @printf("KE_low80 pre=%.3e, post=%.3e, late=%.3e\n", pre_low, post_low, late_low)
    @printf("refl_idx(post/pre)=%.3f, late_idx(10000-12000/pre)=%.3f\n", refl_idx, late_idx)

    println()
    println("--- Below-80 km amplitude/phase fidelity (mid-x probes) ---")
    @printf("Reference forcing period: %.1f s\n", T_FORCE)
    @printf("%-8s %-10s %-10s %-10s %-10s %-10s\n",
            "z(km)", "rms_pre", "rms_late", "late/pre", "T_pre(s)", "T_late(s)")

    for p in 1:nprobe
        zkm = z_phys[probe_j_local[p]] / 1e3
        sig = wscaled_probe[p]

        rms_pre = rms_window(times, sig, 2500.0, 4500.0)
        rms_late = rms_window(times, sig, 9000.0, 11500.0)
        amp_ratio = (isfinite(rms_pre) && rms_pre > 0 && isfinite(rms_late)) ? rms_late / rms_pre : NaN

        T_pre = est_period_lagcorr(times, sig; t1 = 2500.0, t2 = 4500.0)
        T_late = est_period_lagcorr(times, sig; t1 = 9000.0, t2 = 11500.0)

        @printf("%-8.1f %-10.3e %-10.3e %-10.3f %-10.1f %-10.1f\n",
                zkm, rms_pre, rms_late, amp_ratio, T_pre, T_late)
    end

    println()
    println("Done.")
end
