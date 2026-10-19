# =============================================================================
# Sponge Layer Validation Test
# Reproduces Srivastava et al. (2022) Figure 4 (KE evolution)
# and Figure 3 (wave travel diagram)
#
# Parameters from Tables 2 and 3:
#   Domain: 80 km horizontal x 120 km vertical (20 km sponge at top)
#   dx = dz = 0.5 km, T = 239 K, Earth isothermal
#   Gaussian wave packet: w_tilde=0.001 m/s, t0=1200s, sigma_t=600s
#   omega=0.01 rad/s, lambda_x=80 km
#   t_end = 10000 s, CFL = 0.8
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using MAGNUSP
using MAGNUSP: SVector, primitive_from_conservative, phys_x, phys_z
using Printf
using NCDatasets

# --- Configuration matching paper Tables 2 & 3 ---
# Domain: Lx=80km (x from -40 to +40 km), Lz=120km (0 to 120 km)
# 20 km sponge at top: sponge starts at z=100 km, domain top at 120 km
const NX = 160      # 80 km / 0.5 km
const NZ = 240      # 120 km / 0.5 km
const DX = 500.0    # 0.5 km [m]
const DZ = 500.0    # 0.5 km [m]
const T_ISOTHERMAL = 239.0  # K (Srivastava Section 3)
const T_END = 15000.0       # s
const DIAG_DT = 10.0        # save diagnostics every ~10 s

# Wave forcing parameters (Table 3, Eq. 29)
const W_TILDE = 0.001       # m/s
const T0_FORCING = 1200.0   # s
const SIGMA_T = 600.0       # s
const OMEGA = 0.01           # rad/s
const LAMBDA_X = 80e3       # m
const KX = 2.0 * pi / LAMBDA_X

# Print atmosphere info
cs = MAGNUSP.sound_speed(1.4, 287.0, T_ISOTHERMAL)
H = MAGNUSP.scale_height(T_ISOTHERMAL, 9.81, 287.0)
N_bv = MAGNUSP.brunt_vaisala(1.4, 9.81, cs)
m_squared = KX^2 * (N_bv^2 / OMEGA^2 - 1.0)
m = sqrt(m_squared)
lambda_z = 2.0 * pi / m

println("Sponge Layer Validation (Srivastava Fig. 3-4)")
println("=" ^ 60)
println("Atmosphere: Earth isothermal T = $(T_ISOTHERMAL) K")
@printf("  cs = %.1f m/s, H = %.2f km\n", cs, H/1e3)
@printf("  N = %.5f rad/s, omega_ac = %.5f rad/s\n", N_bv,
        MAGNUSP.acoustic_cutoff(1.4, 9.81, cs))
println("\nWave parameters:")
@printf("  omega = %.3f rad/s (gravity wave: omega < N = %.3f)\n", OMEGA, N_bv)
@printf("  lambda_x = %.0f km, lambda_z = %.1f km\n", LAMBDA_X/1e3, lambda_z/1e3)
@printf("  w_tilde = %.4f m/s, t0 = %.0f s, sigma_t = %.0f s\n", W_TILDE, T0_FORCING, SIGMA_T)
println("\nDomain: $(NX*DX/1e3) km x $(NZ*DZ/1e3) km, dx=dz=$(DX) m")
println("Sponge: top 20 km (z_sponge=100 km to z_max=120 km)")
println()

"""
    run_sponge_case(K0; label="")

Run a single sponge test case with given K0 value.
Returns (time_array, ke_array, travel_diagram_data).
"""
function run_sponge_case(K0::Float64; label::String="")
    sponge_frac = 20.0 / 120.0  # 20 km sponge out of 120 km domain

    config = earth_config(
        nx = NX,
        nz = NZ,
        dx = DX,
        dz = DZ,
        CFL = 0.8,
        T0 = T_ISOTHERMAL,
        t_end = T_END,
        output_interval = 999999,  # don't write NetCDF snapshots
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

    # Surface density for density-scaling
    rho_surface = atm.rho0[grid.ng + 1]

    # --- Diagnostic arrays ---
    times = Float64[]
    ke_values = Float64[]
    # Travel diagram: density-scaled w at mid-x vs (t, z)
    i_mid = grid.ng + div(grid.nx, 2)
    nz_phys = grid.nz
    travel_w = Float64[]  # will reshape later

    # --- Time loop ---
    t = 0.0
    step_count = 0
    next_diag = 0.0

    @printf("Running K0 = %.4f %s...\n", K0, label)

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
            # Compute perturbation KE over ALL gridpoints including sponge zone
            # Matches paper Eq. 28: E_k = sum 0.5*rho_0*(u'^2+w'^2) at all gridpoints (Section 2.2.6)
            ke = 0.0
            for j in phys_z(grid), i in phys_x(grid)
                rho_bg = Q_bar[1, i, j]
                u_pert = Q[2, i, j] / Q[1, i, j]
                w_pert = Q[3, i, j] / Q[1, i, j]
                ke += 0.5 * rho_bg * (u_pert^2 + w_pert^2) * DX * DZ
            end

            push!(times, t)
            push!(ke_values, ke)

            # Travel diagram: density-scaled w at mid-x
            for j in phys_z(grid)
                rho = Q[1, i_mid, j]
                w = Q[3, i_mid, j] / rho
                rho_bg = Q_bar[1, i_mid, j]
                # Density-scaled: w * sqrt(rho / rho_surface)
                w_scaled = w * sqrt(rho_bg / rho_surface)
                push!(travel_w, w_scaled)
            end

            next_diag += DIAG_DT

            if step_count % 1000 == 0
                @printf("  step %6d, t = %8.1f s, KE = %12.6e\n", step_count, t, ke)
            end
        end
    end

    nt = length(times)
    travel_diagram = reshape(travel_w, nz_phys, nt)

    @printf("  Done: %d steps, peak KE = %.6e\n\n", step_count, maximum(ke_values))

    return (times, ke_values, travel_diagram)
end

# --- Run cases ---
K0_values = [0.0, 0.0005, 0.005, 0.05, 0.5]
results = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}, Matrix{Float64}}}()

for K0 in K0_values
    times, ke, travel = run_sponge_case(K0; label="(case $(findfirst(==(K0), K0_values))/$(length(K0_values)))")
    results[K0] = (times, ke, travel)
end

# --- Save results to NetCDF ---
outdir = joinpath(@__DIR__, "..", "output", "sponge_validation")
mkpath(outdir)
outfile = joinpath(outdir, "sponge_ke_evolution.nc")

NCDataset(outfile, "c") do ds
    # Use times from K0=0.005 case as reference
    ref_times = results[0.005][1]
    nt = length(ref_times)

    defDim(ds, "time", nt)
    defDim(ds, "K0", length(K0_values))

    ds.attrib["description"] = "Sponge layer validation - KE evolution (Srivastava Fig. 4)"
    ds.attrib["T_isothermal"] = T_ISOTHERMAL
    ds.attrib["omega"] = OMEGA
    ds.attrib["lambda_x"] = LAMBDA_X
    ds.attrib["w_tilde"] = W_TILDE

    t_var = defVar(ds, "time", Float64, ("time",), attrib=Dict("units" => "s"))
    t_var[:] = ref_times

    k0_var = defVar(ds, "K0_values", Float64, ("K0",))
    k0_var[:] = K0_values

    for (idx, K0) in enumerate(K0_values)
        times_k, ke_k, _ = results[K0]
        # Interpolate to reference time grid if needed
        ke_interp = if length(times_k) == nt
            ke_k
        else
            # Simple nearest-neighbor interpolation
            [ke_k[argmin(abs.(times_k .- t_ref))] for t_ref in ref_times]
        end
        varname = @sprintf("KE_K0_%.4f", K0)
        v = defVar(ds, varname, Float64, ("time",),
                   attrib=Dict("units" => "J/m", "K0" => K0))
        v[:] = ke_interp
    end
end
println("KE evolution saved to: $outfile")

# --- Save travel diagram for K0=0.005 ---
travel_file = joinpath(outdir, "travel_diagram_K0_0.005.nc")
times_ref, _, travel_ref = results[0.005]

# Physical z grid
config_ref = earth_config(nx=NX, nz=NZ, dx=DX, dz=DZ, T0=T_ISOTHERMAL)
grid_ref = Grid2D(config_ref)

NCDataset(travel_file, "c") do ds
    nt = length(times_ref)
    nz = grid_ref.nz

    defDim(ds, "time", nt)
    defDim(ds, "z", nz)

    ds.attrib["description"] = "Wave travel diagram - density-scaled w (Srivastava Fig. 3)"
    ds.attrib["K0"] = 0.005

    t_var = defVar(ds, "time", Float64, ("time",), attrib=Dict("units" => "s"))
    z_var = defVar(ds, "z", Float64, ("z",), attrib=Dict("units" => "m"))
    w_var = defVar(ds, "w_scaled", Float64, ("z", "time"),
                   attrib=Dict("units" => "m/s", "long_name" => "w * sqrt(rho/rho_surface)"))

    t_var[:] = times_ref
    z_var[:] = grid_ref.z
    w_var[:, :] = travel_ref
end
println("Travel diagram saved to: $travel_file")

# Also save K0=0 travel diagram
travel_file_nosponge = joinpath(outdir, "travel_diagram_K0_0.nc")
times_ns, _, travel_ns = results[0.0]

NCDataset(travel_file_nosponge, "c") do ds
    nt = length(times_ns)
    nz = grid_ref.nz

    defDim(ds, "time", nt)
    defDim(ds, "z", nz)

    ds.attrib["description"] = "Wave travel diagram - no sponge (Srivastava Fig. 3 left)"
    ds.attrib["K0"] = 0.0

    t_var = defVar(ds, "time", Float64, ("time",), attrib=Dict("units" => "s"))
    z_var = defVar(ds, "z", Float64, ("z",), attrib=Dict("units" => "m"))
    w_var = defVar(ds, "w_scaled", Float64, ("z", "time"),
                   attrib=Dict("units" => "m/s", "long_name" => "w * sqrt(rho/rho_surface)"))

    t_var[:] = times_ns
    z_var[:] = grid_ref.z
    w_var[:, :] = travel_ns
end
println("No-sponge travel diagram saved to: $travel_file_nosponge")

# --- Print summary for quick comparison with Fig. 4 ---
println("\n" * "=" ^ 60)
println("KE Peak Values (compare with Srivastava Fig. 4 peak ~ 0.02)")
println("=" ^ 60)
for K0 in K0_values
    _, ke, _ = results[K0]
    peak_ke = maximum(ke)
    peak_t = results[K0][1][argmax(ke)]
    @printf("  K0 = %6.4f: peak KE = %12.6e at t = %.0f s\n", K0, peak_ke, peak_t)
end

# Reflection coefficient — matches paper Table 4 definition exactly (Section 2.2.6):
# r = sqrt(KE_t5000 / KE_peak) * 100%
# KE is summed over ALL gridpoints (including sponge zone), at exactly t=5000s.
# Paper peak ~0.02 kg/ms² (our units J/m = kg/ms² * dx*dz).
println("\nReflection coefficients at t=5000s (paper Table 4: K0=0→7.6%, K0=0.005→0.55%):")
for K0 in K0_values
    times_k, ke, _ = results[K0]
    ke_peak = maximum(ke)
    t5000_idx = argmin(abs.(times_k .- 5000.0))
    ke_t5000 = ke[t5000_idx]
    r = sqrt(ke_t5000 / ke_peak) * 100.0
    @printf("  K0 = %6.4f: r = %.2f%%  (KE_peak=%.1f, KE_t5000=%.4f)\n",
            K0, r, ke_peak, ke_t5000)
end
