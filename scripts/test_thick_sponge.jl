# =============================================================================
# Thick Sponge Layer Test
# Tests whether a thicker sponge layer eliminates wave reflections.
#
# The wave has lambda_z ~ 46 km, so:
#   - 20 km sponge (current) = 0.43 lambda_z  (less than half a wavelength)
#   - 40 km sponge           = 0.87 lambda_z  (nearly one wavelength)
#   - 60 km sponge           = 1.30 lambda_z  (more than one wavelength)
#
# Uses K0=0.005 (paper reference case) with three sponge thicknesses.
# Domain is always 120 km total; only sponge_frac changes.
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using MAGNUSP
using MAGNUSP: SVector, primitive_from_conservative, phys_x, phys_z
using Printf

# --- Configuration matching paper Tables 2 & 3 ---
const NX = 160
const NZ = 240
const DX = 500.0
const DZ = 500.0
const T_ISOTHERMAL = 239.0
const T_END = 8000.0
const DIAG_DT = 10.0

# Wave forcing parameters (Table 3)
const W_TILDE = 0.001
const T0_FORCING = 1200.0
const SIGMA_T = 600.0
const OMEGA = 0.01
const LAMBDA_X = 80e3
const KX = 2.0 * pi / LAMBDA_X

# K0 for sponge damping (paper reference case)
const K0 = 0.005

# Print atmosphere info
cs = MAGNUSP.sound_speed(1.4, 287.0, T_ISOTHERMAL)
H = MAGNUSP.scale_height(T_ISOTHERMAL, 9.81, 287.0)
N_bv = MAGNUSP.brunt_vaisala(1.4, 9.81, cs)
m_squared = KX^2 * (N_bv^2 / OMEGA^2 - 1.0)
m = sqrt(m_squared)
lambda_z = 2.0 * pi / m

println("Thick Sponge Layer Test")
println("=" ^ 60)
println("Atmosphere: Earth isothermal T = $(T_ISOTHERMAL) K")
@printf("  cs = %.1f m/s, H = %.2f km\n", cs, H/1e3)
@printf("  N = %.5f rad/s\n", N_bv)
@printf("  lambda_z = %.1f km\n", lambda_z/1e3)
println("\nWave: K0 = $(K0), omega = $(OMEGA) rad/s, lambda_x = $(LAMBDA_X/1e3) km")
println("Domain: $(NX*DX/1e3) km x $(NZ*DZ/1e3) km, t_end = $(T_END) s")
println()

# --- Sponge thickness cases ---
sponge_thicknesses_km = [20.0, 40.0, 60.0]
domain_height_km = NZ * DZ / 1e3  # 120 km

"""
    run_thick_sponge_case(sponge_km)

Run a single case with the given sponge thickness (in km).
Returns (times, ke_values).
"""
function run_thick_sponge_case(sponge_km::Float64)
    sponge_frac = sponge_km / domain_height_km

    @printf("Running sponge = %.0f km (%.1f%% of domain, %.2f lambda_z)...\n",
            sponge_km, sponge_frac * 100.0, sponge_km / (lambda_z / 1e3))

    config = earth_config(
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

    # --- Setup ---
    grid = Grid2D(config)
    atm = isothermal_atmosphere(config, Float64.(grid.z_full))
    Q = MAGNUSP.init_state(grid, atm, config)
    Q_bar = copy(Q)

    sponge_layer = MAGNUSP.init_sponge(grid, config)
    diff_solver = MAGNUSP.DiffusionSolver(grid.nx, grid.nz)
    Q_half_x = similar(Q)
    Q_half_z = similar(Q)

    # --- Diagnostic arrays ---
    times = Float64[]
    ke_values = Float64[]

    # --- Time loop ---
    t = 0.0
    step_count = 0
    next_diag = 0.0

    while t < config.t_end
        dt = MAGNUSP.compute_adaptive_dt(Q, grid, config.CFL, config.gamma, config.R_gas)
        dt = min(dt, config.t_end - t)

        # BCs
        MAGNUSP.apply_all_bcs!(Q, Q_bar, grid, config.gamma)

        # Bottom forcing (Gaussian envelope)
        MAGNUSP.apply_bottom_forcing!(Q, grid, config, t, wave_params, :gaussian)

        # Euler step
        MAGNUSP.euler_step!(Q, Q_bar, Q_half_x, Q_half_z, grid, dt, config.gamma, config.g,
                            sponge_layer.K)

        # Diffusion
        if config.enable_viscosity
            MAGNUSP.diffusion_step!(Q, Q_bar, grid, dt, config, diff_solver)
        end

        t += dt
        step_count += 1

        # Diagnostics at regular time intervals
        if t >= next_diag
            # Compute perturbation KE over all gridpoints (including sponge zone)
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

            # Print KE at every 1000s
            if t >= 1000.0 && abs(t - round(t / 1000.0) * 1000.0) < DIAG_DT
                @printf("  t = %6.0f s, KE = %12.6e\n", t, ke)
            end
        end
    end

    @printf("  Done: %d steps, peak KE = %.6e\n\n", step_count, maximum(ke_values))
    return (times, ke_values)
end

# --- Run all three cases ---
results = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}()

for sponge_km in sponge_thicknesses_km
    times, ke = run_thick_sponge_case(sponge_km)
    results[sponge_km] = (times, ke)
end

# --- Summary ---
println("=" ^ 60)
println("RESULTS: Thick Sponge Test (K0 = $(K0))")
println("=" ^ 60)
@printf("  lambda_z = %.1f km\n\n", lambda_z/1e3)

@printf("  %-12s  %-12s  %-14s  %-10s  %s\n",
        "Sponge (km)", "Frac (%)", "Sponge/lam_z", "KE_peak", "r (%)")
println("  " * "-" ^ 68)

for sponge_km in sponge_thicknesses_km
    times_k, ke_k = results[sponge_km]
    ke_peak = maximum(ke_k)

    # Reflection coefficient: r = sqrt(KE_t5000 / KE_peak) * 100%
    t5000_idx = argmin(abs.(times_k .- 5000.0))
    ke_t5000 = ke_k[t5000_idx]
    r = sqrt(ke_t5000 / ke_peak) * 100.0

    sponge_frac_pct = sponge_km / domain_height_km * 100.0
    sponge_over_lambda = sponge_km / (lambda_z / 1e3)

    @printf("  %-12.0f  %-12.1f  %-14.2f  %-10.4e  %.2f\n",
            sponge_km, sponge_frac_pct, sponge_over_lambda, ke_peak, r)
end

println()
println("Lower r (%) = less reflection = better sponge absorption.")
println("If r decreases significantly with thicker sponge, the 20 km")
println("sponge is simply too thin relative to lambda_z ~ $(round(lambda_z/1e3, digits=1)) km.")
