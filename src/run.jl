# =============================================================================
# Main simulation time loop orchestrator
# =============================================================================

"""
    SimulationState

Holds all mutable state for a running simulation.
"""
mutable struct SimulationState
    Q::Array{Float64,3}          # current state
    Q_bar::Array{Float64,3}      # reference/background state
    grid::Grid2D
    config::SimConfig
    atm::AtmosphereProfile
    sponge::SpongeLayer
    diff_solver::DiffusionSolver
    t::Float64
    step::Int
    coupling_buffer::Vector{Dict{Symbol, Vector{Float64}}}  # time series at z_couple
end

"""
    run_simulation(config; forcing_type=:linear, wave_params=WaveParams(),
                   outdir="output", diag_interval=100, coupling_altitude=nothing)

Main entry point: run the MAGNUS-P simulation from t=0 to t_end.

Returns the final SimulationState for post-processing.
"""
function run_simulation(config::SimConfig;
                        forcing_type::Symbol = :linear,
                        wave_params = WaveParams(),
                        pulse_params::Union{Nothing, PulseParams} = nothing,
                        u_init::Union{Float64, Vector{Float64}, Function} = 0.0,
                        w_init::Union{Float64, Vector{Float64}, Function} = 0.0,
                        bc_type::Symbol = :eq19,
                        bc_x_type::Symbol = :periodic,
                        bc_bottom_type::Symbol = :reflective,
                        outdir::String = "output",
                        diag_interval::Int = 100,
                        coupling_altitude::Union{Nothing, Float64} = nothing,
                        atm_file::Union{Nothing, String} = nothing)
    # --- Setup ---
    grid = Grid2D(config)

    atm = if atm_file !== nothing
        tabulated_atmosphere(atm_file, Float64.(grid.z_full))
    else
        isothermal_atmosphere(config, Float64.(grid.z_full))
    end

    Q = init_state(grid, atm, config, u0=u_init, w0=w_init)
    Q_bar = copy(Q)

    if pulse_params !== nothing
        pressure_pulse_ic!(Q, grid, config, pulse_params)
    end

    sponge_layer = init_sponge(grid, config)
    diff_solver = DiffusionSolver(grid.nx, grid.nz)

    # Workspace arrays for Euler solver
    Q_half_x = similar(Q)
    Q_half_z = similar(Q)

    # Coupling buffer
    coupling_buffer = Dict{Symbol, Vector{Float64}}[]

    # Find coupling altitude index (if specified)
    j_couple = nothing
    if coupling_altitude !== nothing
        j_couple = argmin(abs.(grid.z .- coupling_altitude)) + grid.ng
    end

    # --- Klemp-Durran radiation BC initialization ---
    kd_state = nothing
    if bc_type == :klemp_durran
        j_top = grid.nz + grid.ng
        T_top = atm.T0[j_top]
        T_below = atm.T0[j_top - 1]
        rho_top = atm.rho0[j_top]
        N_local = brunt_vaisala_local(T_below, T_top, grid.dz,
                                      config.g, config.gamma, config.R_gas)
        kd_state = init_klemp_durran(grid, N_local, rho_top)
        @printf("  Klemp-Durran BC: N_local=%.4f rad/s, rho_top=%.3e kg/m^3\n",
                N_local, rho_top)
    end

    # --- Initial diagnostics ---
    @printf("MAGNUS-P Solver initialized\n")
    @printf("  Grid: %d x %d (dx=%.0f m, dz=%.0f m)\n", grid.nx, grid.nz, grid.dx, grid.dz)
    @printf("  Domain: %.1f km x %.1f km\n", domain_width(grid)/1e3, domain_height(grid)/1e3)
    @printf("  Planet: %s (gamma=%.2f, g=%.2f, R=%.1f)\n",
            config.planet, config.gamma, config.g, config.R_gas)

    cs = sound_speed(config.gamma, config.R_gas, config.T0)
    H = scale_height(config.T0, config.g, config.R_gas)
    N_bv = brunt_vaisala(config.gamma, config.g, cs)
    omega_ac = acoustic_cutoff(config.gamma, config.g, cs)
    @printf("  c_s=%.1f m/s, H=%.1f km, N=%.4f rad/s, omega_ac=%.4f rad/s\n",
            cs, H/1e3, N_bv, omega_ac)
    @printf("  Forcing: %s, t_end=%.1f s\n", forcing_type, config.t_end)
    @printf("  Viscosity: %s, Cooling: %s\n",
            config.enable_viscosity ? "ON" : "OFF",
            config.enable_cooling ? "ON" : "OFF")
    println("="^80)

    # --- Time loop ---
    t = 0.0
    step_count = 0
    wall_start = time()

    # Initial snapshot
    write_snapshot_nc(Q, grid, config, t, step_count, outdir)

    while t < config.t_end
        # Adaptive timestep
        dt = compute_adaptive_dt(Q, grid, config.CFL, config.gamma, config.R_gas)
        dt = min(dt, config.t_end - t)  # don't overshoot t_end

        # 1. Apply boundary conditions
        apply_all_bcs!(Q, Q_bar, grid, config.gamma;
                       bc_type=bc_type, bc_x_type=bc_x_type,
                       bc_bottom_type=bc_bottom_type, kd_state=kd_state)

        # 2. Apply bottom forcing
        apply_bottom_forcing!(Q, grid, config, t, wave_params, forcing_type)

        # 3. Euler step (hyperbolic) — sponge damping is integrated inside as source term
        euler_step!(Q, Q_bar, Q_half_x, Q_half_z, grid, dt, config.gamma, config.g,
                    sponge_layer.K)

        # Enforce positivity to prevent floating point crashes at 160km
        for j in 1:grid.nz_total
            for i in 1:grid.nx_total
                rho_old = Q[1, i, j]
                rho_min = Q_bar[1, i, j] * 0.1

                if rho_old < rho_min
                    # Floor density and damp momentum to prevent velocity spikes
                    Q[1, i, j] = rho_min
                    Q[2, i, j] = Q_bar[2, i, j]
                    Q[3, i, j] = Q_bar[3, i, j]
                    Q[4, i, j] = Q_bar[4, i, j]
                end

                rho = Q[1, i, j]
                u = Q[2, i, j] / rho
                w = Q[3, i, j] / rho
                E = Q[4, i, j]
                P = (config.gamma - 1.0) * (E - 0.5 * rho * (u^2 + w^2))

                P_min = Q_bar[4, i, j] * (config.gamma - 1.0) * 0.1
                if P < P_min
                    Q[4, i, j] = P_min / (config.gamma - 1.0) + 0.5 * rho * (u^2 + w^2)
                end
            end
        end

        # 4. Diffusion step (implicit, if enabled)
        if config.enable_viscosity
            diffusion_step!(Q, Q_bar, grid, dt, config, diff_solver)
        end

        # 6. Newtonian cooling (if enabled)
        if config.enable_cooling
            newtonian_cooling_step!(Q, grid, config, atm, dt)
        end

        # Advance time
        t += dt
        step_count += 1

        # Diagnostics
        if step_count % diag_interval == 0
            print_diagnostics(step_count, t, dt, Q, Q_bar, grid, config)

            # Check for blowup
            min_rho, min_P = check_positivity(Q, grid, config.gamma)
            if min_rho < 0 || min_P < 0 || isnan(min_rho) || isnan(min_P)
                @printf("ERROR: Simulation blowup at step %d, t=%.3f s\n", step_count, t)
                @printf("  min_rho = %e, min_P = %e\n", min_rho, min_P)
                break
            end
        end

        # Output snapshots
        if step_count % config.output_interval == 0
            write_snapshot_nc(Q, grid, config, t, step_count, outdir)
        end

        # Extract coupling data at z_couple
        if j_couple !== nothing && step_count % diag_interval == 0
            record = Dict{Symbol, Vector{Float64}}()
            record[:time] = [t]
            ng = grid.ng
            nx = grid.nx
            n_vals = zeros(nx)
            T_vals = zeros(nx)
            u_vals = zeros(nx)
            w_vals = zeros(nx)
            for i in 1:nx
                ig = i + ng
                Qi = SVector(Q[1,ig,j_couple], Q[2,ig,j_couple],
                             Q[3,ig,j_couple], Q[4,ig,j_couple])
                rho, uv, wv, P = primitive_from_conservative(Qi, config.gamma)
                T_val = P / (rho * config.R_gas)
                n_vals[i] = rho / (M_CO2 / N_A)  # number density
                T_vals[i] = T_val
                u_vals[i] = uv
                w_vals[i] = wv
            end
            record[:n] = n_vals
            record[:T] = T_vals
            record[:u] = u_vals
            record[:w] = w_vals
            push!(coupling_buffer, record)
        end
    end

    wall_end = time()
    duration = wall_end - wall_start
    @printf("\nSimulation complete: %d steps, t=%.3f s\n", step_count, t)
    @printf("Wall-clock duration: %.2f seconds\n", duration)

    return SimulationState(Q, Q_bar, grid, config, atm, sponge_layer, diff_solver,
                           t, step_count, coupling_buffer)
end
