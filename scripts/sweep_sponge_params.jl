# =============================================================================
# Parameter Sweep: Sponge Validation
#
# Systematically explores parameter combinations to identify what reproduces
# the paper's Table 4 reflection coefficients (K0=0 -> r=7.6%, K0=0.005 -> 0.55%).
#
# Sweep dimensions:
#   T:               [239, 289] K         (tests N/T inconsistency)
#   Sponge thickness: [20, 40, 60] km     (tests Klemp minimum thickness rule)
#   K0:              [0, 0.005]           (isolates BC vs sponge effect)
#   Viscosity:       [ON, OFF]           (tests physical damping contribution)
#   KE region:       [full, below-sponge] (post-processing, no extra runs)
#
# K0=0 is independent of sponge thickness -> 2T x 2visc = 4 runs
# K0=0.005: 2T x 3sponge x 2visc = 12 runs
# Total: 16 runs, each t_end=8000s
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using MAGNUSP
using MAGNUSP: phys_x, phys_z
using Printf

# --- Fixed parameters (paper Tables 2 & 3) ---
const NX = 160       # 80 km / 0.5 km
const NZ = 240       # 120 km / 0.5 km
const DX = 500.0     # 0.5 km
const DZ = 500.0     # 0.5 km
const DOMAIN_KM = NZ * DZ / 1e3  # 120 km
const T_END = 8000.0
const DIAG_DT = 10.0

# Wave forcing (Table 3, Eq. 29)
const W_TILDE = 0.001
const T0_FORCING = 1200.0
const SIGMA_T = 600.0
const OMEGA = 0.01
const LAMBDA_X = 80e3
const KX = 2.0 * pi / LAMBDA_X

# --- Sweep dimensions ---
T_values = [239.0, 289.0]
sponge_km_values = [20.0, 40.0, 60.0]
K0_values = [0.0, 0.005]
visc_values = [false, true]

# --- Result storage ---
struct SweepResult
    T::Float64
    sponge_km::Float64
    K0::Float64
    viscosity::Bool
    times::Vector{Float64}
    ke_full::Vector{Float64}       # KE over full domain (0 to z_max)
    ke_below::Vector{Float64}      # KE below sponge only (0 to z_sponge)
    ke_peak_full::Float64
    ke_peak_below::Float64
    r_full::Float64                # reflection coeff using full-domain KE
    r_below::Float64               # reflection coeff using below-sponge KE
    lambda_z_km::Float64
    N_bv::Float64
    cs::Float64
end

"""
    run_sweep_case(T, sponge_km, K0, enable_viscosity)

Run a single sweep case. Returns SweepResult with both KE diagnostics.
"""
function run_sweep_case(T::Float64, sponge_km::Float64, K0::Float64,
                        enable_viscosity::Bool)
    sponge_frac = sponge_km / DOMAIN_KM

    config = earth_config(
        nx = NX, nz = NZ, dx = DX, dz = DZ, CFL = 0.8, T0 = T,
        t_end = T_END, output_interval = 999999,
        enable_viscosity = enable_viscosity, enable_cooling = false,
        sponge_frac = sponge_frac, K0_sponge = K0,
    )

    wave_params = MAGNUSP.WaveParams(
        amplitude = W_TILDE, omega = OMEGA, kx = KX,
        sigma_t = SIGMA_T, t0 = T0_FORCING,
    )

    # Atmosphere info
    cs = MAGNUSP.sound_speed(config.gamma, config.R_gas, T)
    N_bv = MAGNUSP.brunt_vaisala(config.gamma, config.g, cs)
    m_sq = KX^2 * (N_bv^2 / OMEGA^2 - 1.0)
    lambda_z_km = m_sq > 0 ? 2.0 * pi / sqrt(m_sq) / 1e3 : NaN

    # Setup
    grid = Grid2D(config)
    atm = isothermal_atmosphere(config, Float64.(grid.z_full))
    Q = MAGNUSP.init_state(grid, atm, config)
    Q_bar = copy(Q)
    sponge_layer = MAGNUSP.init_sponge(grid, config)
    diff_solver = MAGNUSP.DiffusionSolver(grid.nx, grid.nz)
    Q_hx = similar(Q); Q_hz = similar(Q)

    # Sponge boundary z-index: z_sponge = (1 - sponge_frac) * domain_height
    z_sponge_m = (1.0 - sponge_frac) * DOMAIN_KM * 1e3
    ng = grid.ng
    # j index where z >= z_sponge (first sponge cell)
    j_sponge_start = ng + 1 + Int(round(z_sponge_m / DZ))

    # Diagnostic arrays
    times = Float64[]
    ke_full_arr = Float64[]
    ke_below_arr = Float64[]

    t = 0.0
    step_count = 0
    next_diag = 0.0

    while t < config.t_end
        dt = MAGNUSP.compute_adaptive_dt(Q, grid, config.CFL, config.gamma, config.R_gas)
        dt = min(dt, config.t_end - t)

        MAGNUSP.apply_all_bcs!(Q, Q_bar, grid, config.gamma)
        MAGNUSP.apply_bottom_forcing!(Q, grid, config, t, wave_params, :gaussian)
        MAGNUSP.euler_step!(Q, Q_bar, Q_hx, Q_hz, grid, dt, config.gamma, config.g,
                            sponge_layer.K)
        if enable_viscosity
            MAGNUSP.diffusion_step!(Q, Q_bar, grid, dt, config, diff_solver)
        end

        t += dt
        step_count += 1

        if t >= next_diag
            ke_full = 0.0
            ke_below = 0.0
            for j in phys_z(grid), i in phys_x(grid)
                u = Q[2, i, j] / Q[1, i, j]
                w = Q[3, i, j] / Q[1, i, j]
                ke_cell = 0.5 * Q_bar[1, i, j] * (u^2 + w^2) * DX * DZ
                ke_full += ke_cell
                if j < j_sponge_start
                    ke_below += ke_cell
                end
            end
            push!(times, t)
            push!(ke_full_arr, ke_full)
            push!(ke_below_arr, ke_below)
            next_diag += DIAG_DT
        end
    end

    # Compute reflection coefficients
    ke_peak_full = maximum(ke_full_arr)
    ke_peak_below = maximum(ke_below_arr)
    t5000_idx = argmin(abs.(times .- 5000.0))

    r_full = ke_peak_full > 0 ? sqrt(ke_full_arr[t5000_idx] / ke_peak_full) * 100.0 : NaN
    r_below = ke_peak_below > 0 ? sqrt(ke_below_arr[t5000_idx] / ke_peak_below) * 100.0 : NaN

    return SweepResult(T, sponge_km, K0, enable_viscosity,
                       times, ke_full_arr, ke_below_arr,
                       ke_peak_full, ke_peak_below,
                       r_full, r_below, lambda_z_km, N_bv, cs)
end

# =============================================================================
# Main sweep
# =============================================================================
println("=" ^ 80)
println("PARAMETER SWEEP: Sponge Validation")
println("=" ^ 80)

# Print atmosphere summary for each T
for T in T_values
    cs = MAGNUSP.sound_speed(1.4, 287.0, T)
    H = MAGNUSP.scale_height(T, 9.81, 287.0)
    N = MAGNUSP.brunt_vaisala(1.4, 9.81, cs)
    m_sq = KX^2 * (N^2 / OMEGA^2 - 1.0)
    lz = 2.0 * pi / sqrt(m_sq) / 1e3
    @printf("  T=%3.0fK: cs=%.1f m/s, H=%.2f km, N=%.5f rad/s, lambda_z=%.1f km\n",
            T, cs, H/1e3, N, lz)
end
println()

# Build run list — K0=0 doesn't depend on sponge thickness, so run only once per (T, visc)
run_list = Tuple{Float64, Float64, Float64, Bool}[]

for T in T_values, visc in visc_values
    # K0=0: sponge thickness irrelevant, run with 20 km (arbitrary)
    push!(run_list, (T, 20.0, 0.0, visc))
    # K0=0.005: sweep sponge thickness
    for sp in sponge_km_values
        push!(run_list, (T, sp, 0.005, visc))
    end
end

results = SweepResult[]
total = length(run_list)

for (idx, (T, sp, K0, visc)) in enumerate(run_list)
    visc_str = visc ? "ON" : "OFF"
    @printf("[%2d/%2d] T=%3.0fK, sponge=%2.0fkm, K0=%.4f, visc=%s\n",
            idx, total, T, sp, K0, visc_str)
    result = run_sweep_case(T, sp, K0, visc)
    push!(results, result)
    @printf("        r_full=%.2f%%, r_below=%.2f%%, peak_full=%.2e, peak_below=%.2e\n",
            result.r_full, result.r_below, result.ke_peak_full, result.ke_peak_below)
end

# =============================================================================
# Summary tables
# =============================================================================
println()
println("=" ^ 100)
println("SUMMARY: Reflection Coefficients r (%)")
println("Paper targets: K0=0 -> 7.6%, K0=0.005 -> 0.55%")
println("=" ^ 100)

# Table header
@printf("%-5s  %-10s  %-6s  %-5s  %-9s  %-10s  %-10s  %-10s  %-10s  %-10s\n",
        "T(K)", "Sponge", "K0", "Visc", "lam_z(km)", "r_full(%)", "r_below(%)",
        "KE_pk_full", "KE_pk_bel", "KE@5000_f")
println("-" ^ 100)

for r in results
    visc_str = r.viscosity ? "ON" : "OFF"
    sp_str = @sprintf("%.0fkm (%.1flz)", r.sponge_km, r.sponge_km / r.lambda_z_km)
    t5000_idx = argmin(abs.(r.times .- 5000.0))
    @printf("%-5.0f  %-10s  %-6.3f  %-5s  %-9.1f  %-10.2f  %-10.2f  %-10.3e  %-10.3e  %-10.3e\n",
            r.T, sp_str, r.K0, visc_str,
            r.lambda_z_km, r.r_full, r.r_below,
            r.ke_peak_full, r.ke_peak_below, r.ke_full[t5000_idx])
end

# =============================================================================
# Grouped analysis: what combination gets closest to paper?
# =============================================================================
println()
println("=" ^ 80)
println("CLOSEST TO PAPER TARGETS")
println("=" ^ 80)

# K0=0 target: r=7.6%
k0_zero = filter(r -> r.K0 == 0.0, results)
if !isempty(k0_zero)
    println("\nK0=0 (paper target: r=7.6%):")
    for r in sort(k0_zero, by=x -> abs(x.r_full - 7.6))
        visc_str = r.viscosity ? "ON" : "OFF"
        @printf("  T=%3.0fK, visc=%-3s: r_full=%.2f%%, r_below=%.2f%%  [delta=%.2f%%]\n",
                r.T, visc_str, r.r_full, r.r_below, abs(r.r_full - 7.6))
    end
end

# K0=0.005 target: r=0.55%
k0_005 = filter(r -> r.K0 == 0.005, results)
if !isempty(k0_005)
    println("\nK0=0.005 (paper target: r=0.55%):")
    for r in sort(k0_005, by=x -> min(abs(x.r_full - 0.55), abs(x.r_below - 0.55)))
        visc_str = r.viscosity ? "ON" : "OFF"
        sp_lz = r.sponge_km / r.lambda_z_km
        @printf("  T=%3.0fK, sponge=%2.0fkm(%.1flz), visc=%-3s: r_full=%.2f%%, r_below=%.2f%%  [delta_f=%.2f%%, delta_b=%.2f%%]\n",
                r.T, r.sponge_km, sp_lz, visc_str,
                r.r_full, r.r_below,
                abs(r.r_full - 0.55), abs(r.r_below - 0.55))
    end
end

println()
println("Done. Total runs: $total")
