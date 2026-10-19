# =============================================================================
# Sponge Thickness Sweep for R3 Acoustic Wave Parameters
#
# Purpose:
#   Measure reflection coefficient r(%) as a function of sponge thickness
#   at the R3 wave parameters (omega=0.02, acoustic regime) on Mars.
#   Uses Gaussian wave packet forcing to isolate incident vs reflected energy.
#
# Sweep dimensions:
#   Sponge thickness: [30, 60, 90, 120, 150] km
#   K0:              [0.0, 0.005]
#
# Domain: 200 km vertical (dz=500m, nz=400), 200 km horizontal (dx=1km, nx=200)
# Atmosphere: Mars isothermal (gamma=1.31, g=3.72, R=188.9, T=210K)
# Note: 200 km domain is proven stable (used for actual R3 run). Sponge onset
# altitude varies with thickness, but reflection diagnostic remains valid.
#
# Diagnostics:
#   r(%) = sqrt(KE(t_late) / KE_peak) * 100
#   KE measured in full domain and below-sponge region separately
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using MAGNUSP
using MAGNUSP: phys_x, phys_z
using Printf

# --- Domain ---
const LX = 200_000.0    # 200 km (one horizontal wavelength)
const LZ = 200_000.0    # 200 km (proven stable for Mars isothermal)
const DX = 1000.0       # 1 km horizontal
const DZ = 500.0        # 0.5 km vertical
const NX = Int(round(LX / DX))   # 200
const NZ = Int(round(LZ / DZ))   # 400

# --- Simulation control ---
const T_ISO = 210.0     # Mars isothermal temperature [K]
const CFL = 0.8
const T_END = 10000.0   # enough for packet to travel up and reflected energy to settle
const DIAG_DT = 10.0    # diagnostic sampling interval [s]

# --- R3 wave parameters ---
const W_TILDE = 0.001   # small amplitude for clean linear test
const OMEGA = 0.02      # rad/s (acoustic regime, > omega_ac ~ 0.011)
const LAMBDA_X = 200e3  # 200 km
const KX = 2.0 * pi / LAMBDA_X
const SIGMA_T = 600.0   # Gaussian temporal width [s]
const T0_FORCING = 1200.0  # Gaussian center time [s]

# --- Sweep dimensions ---
# Max 120 km sponge in 200 km domain (leaves 80 km physical domain)
const SPONGE_KM_VALUES = [30.0, 60.0, 90.0, 120.0]
const K0_VALUES = [0.0, 0.005]
const DOMAIN_KM = LZ / 1e3  # 300 km

# --- Result storage ---
struct SweepResult
    sponge_km::Float64
    K0::Float64
    times::Vector{Float64}
    ke_full::Vector{Float64}
    ke_below::Vector{Float64}
    ke_peak_full::Float64
    ke_peak_below::Float64
    r_full::Float64
    r_below::Float64
    lambda_z_km::Float64
    sponge_over_lz::Float64
end

"""
Compute vertical wavelength for acoustic branch dispersion relation.
kz² = (ω² - ωac²)/cs² - kx²(1 - N²/ω²)
"""
function acoustic_lambda_z(gamma, g, cs, omega, kx)
    omega_ac = gamma * g / (2.0 * cs)
    N_bv = sqrt((gamma - 1.0) / gamma) * g / cs
    kz_sq = (omega^2 - omega_ac^2) / cs^2 - kx^2 * (1.0 - N_bv^2 / omega^2)
    if kz_sq > 0
        return 2.0 * pi / sqrt(kz_sq)
    else
        return NaN
    end
end

function run_sweep_case(sponge_km::Float64, K0::Float64, lambda_z_km::Float64)
    sponge_frac = sponge_km / DOMAIN_KM

    config = mars_config(
        nx = NX, nz = NZ, dx = DX, dz = DZ, CFL = CFL, T0 = T_ISO,
        t_end = T_END, output_interval = 999999,
        enable_viscosity = false, enable_cooling = false,
        sponge_frac = sponge_frac, K0_sponge = K0,
    )

    wave_params = MAGNUSP.WaveParams(
        amplitude = W_TILDE, omega = OMEGA, kx = KX,
        sigma_t = SIGMA_T, t0 = T0_FORCING,
    )

    # Setup
    grid = Grid2D(config)
    atm = isothermal_atmosphere(config, Float64.(grid.z_full))
    Q = MAGNUSP.init_state(grid, atm, config)
    Q_bar = copy(Q)
    sponge_layer = MAGNUSP.init_sponge(grid, config)
    Q_hx = similar(Q); Q_hz = similar(Q)

    # Sponge boundary index
    z_sponge_m = (1.0 - sponge_frac) * LZ
    ng = grid.ng
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
        if !isfinite(dt) || dt <= 0.0
            @printf("  CRASH at t=%.1f s (invalid dt)\n", t)
            break
        end
        dt = min(dt, config.t_end - t)

        MAGNUSP.apply_all_bcs!(Q, Q_bar, grid, config.gamma)
        MAGNUSP.apply_bottom_forcing!(Q, grid, config, t, wave_params, :gaussian)
        MAGNUSP.euler_step!(Q, Q_bar, Q_hx, Q_hz, grid, dt, config.gamma, config.g,
                            sponge_layer.K)

        t += dt
        step_count += 1

        if t >= next_diag
            ke_full = 0.0
            ke_below = 0.0
            for j in phys_z(grid), i in phys_x(grid)
                rho = Q[1, i, j]
                u = Q[2, i, j] / rho
                w = Q[3, i, j] / rho
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

    if isempty(ke_full_arr)
        return SweepResult(sponge_km, K0, Float64[], Float64[], Float64[],
                           NaN, NaN, NaN, NaN, lambda_z_km, sponge_km / lambda_z_km)
    end

    # Compute reflection coefficients
    ke_peak_full = maximum(ke_full_arr)
    ke_peak_below = maximum(ke_below_arr)

    # Measure residual KE at t=8000s (well after packet has interacted with sponge)
    t_late = 8000.0
    i_late = argmin(abs.(times .- t_late))
    r_full = ke_peak_full > 0 ? sqrt(ke_full_arr[i_late] / ke_peak_full) * 100.0 : NaN
    r_below = ke_peak_below > 0 ? sqrt(ke_below_arr[i_late] / ke_peak_below) * 100.0 : NaN

    @printf("  done: steps=%d, r_full=%.2f%%, r_below=%.2f%%\n", step_count, r_full, r_below)

    return SweepResult(sponge_km, K0,
                       times, ke_full_arr, ke_below_arr,
                       ke_peak_full, ke_peak_below,
                       r_full, r_below, lambda_z_km, sponge_km / lambda_z_km)
end

# =============================================================================
# Main
# =============================================================================
println("=" ^ 90)
println("SPONGE THICKNESS SWEEP — R3 ACOUSTIC WAVE (Mars)")
println("=" ^ 90)

# Compute wave properties
cs = MAGNUSP.sound_speed(1.31, 188.9, T_ISO)
lz_m = acoustic_lambda_z(1.31, 3.72, cs, OMEGA, KX)
lz_km = lz_m / 1e3
omega_ac = 1.31 * 3.72 / (2.0 * cs)
N_bv = sqrt((1.31 - 1.0) / 1.31) * 3.72 / cs
c_gz = cs^2 * (2*pi/lz_m) / OMEGA  # vertical group velocity

@printf("Mars isothermal: T=%.0f K, cs=%.1f m/s, H=%.1f km\n", T_ISO, cs, 188.9*T_ISO/3.72/1e3)
@printf("omega_ac=%.4f rad/s, N=%.4f rad/s, omega=%.3f rad/s (acoustic branch)\n", omega_ac, N_bv, OMEGA)
@printf("lambda_z=%.1f km, c_gz=%.1f m/s\n", lz_km, c_gz)
@printf("Domain: %.0f x %.0f km, dx=%.0f m, dz=%.0f m, nx=%d, nz=%d\n",
        LX/1e3, LZ/1e3, DX, DZ, NX, NZ)
println()

# Build run list
run_list = Tuple{Float64, Float64}[]
# K0=0 baseline (sponge thickness irrelevant, run once with 30 km)
push!(run_list, (30.0, 0.0))
# K0=0.005 at each sponge thickness
for sp in SPONGE_KM_VALUES
    push!(run_list, (sp, 0.005))
end

results = SweepResult[]
total = length(run_list)
wall_start = time()

for (idx, (sp, K0)) in enumerate(run_list)
    @printf("[%d/%d] sponge=%.0f km (%.2f λz), K0=%.4f\n",
            idx, total, sp, sp / lz_km, K0)
    t0 = time()
    result = run_sweep_case(sp, K0, lz_km)
    elapsed = time() - t0
    push!(results, result)
    @printf("  wall-clock: %.1f s\n", elapsed)
end

wall_total = time() - wall_start

# =============================================================================
# Summary
# =============================================================================
println()
println("=" ^ 90)
println("SUMMARY")
println("=" ^ 90)
@printf("lambda_z = %.1f km\n\n", lz_km)

@printf("%-12s  %-12s  %-6s  %-10s  %-10s  %-12s  %-12s\n",
        "Sponge(km)", "Sponge/λz", "K0", "r_full(%)", "r_below(%)",
        "KE_pk_full", "KE_pk_below")
println("-" ^ 90)

for r in results
    @printf("%-12.0f  %-12.2f  %-6.3f  %-10.2f  %-10.2f  %-12.3e  %-12.3e\n",
            r.sponge_km, r.sponge_over_lz, r.K0,
            r.r_full, r.r_below,
            r.ke_peak_full, r.ke_peak_below)
end

# Best configuration
k0_005 = filter(x -> x.K0 == 0.005 && isfinite(x.r_below), results)
if !isempty(k0_005)
    best = k0_005[argmin([x.r_below for x in k0_005])]
    println()
    @printf("Best: sponge=%.0f km (%.2f λz), r_below=%.2f%%\n",
            best.sponge_km, best.sponge_over_lz, best.r_below)
end

@printf("\nTotal wall-clock: %.1f s (%.1f min)\n", wall_total, wall_total / 60.0)
println("Done.")
