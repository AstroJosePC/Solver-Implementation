# =============================================================================
# Diagnostic monitors for simulation health and validation
# =============================================================================

"""
    compute_total_energy(Q, grid, config)

Compute volume-integrated total energy over the physical domain.
"""
function compute_total_energy(Q::Array{Float64,3}, grid::Grid2D)
    E_total = 0.0
    dV = grid.dx * grid.dz
    for j in phys_z(grid), i in phys_x(grid)
        E_total += Q[4, i, j] * dV
    end
    return E_total
end

"""
    compute_kinetic_energy(Q, grid)

Compute volume-integrated kinetic energy: sum of 0.5*rho*(u^2+w^2)*dV.
"""
function compute_kinetic_energy(Q::Array{Float64,3}, grid::Grid2D)
    KE = 0.0
    dV = grid.dx * grid.dz
    for j in phys_z(grid), i in phys_x(grid)
        rho = Q[1, i, j]
        u = Q[2, i, j] / rho
        w = Q[3, i, j] / rho
        KE += 0.5 * rho * (u^2 + w^2) * dV
    end
    return KE
end

"""
    check_hydrostatic_drift(Q, Q_bar, grid)

Compute L2 and Linf norms of perturbation (Q - Q_bar) relative to Q_bar.
Returns (L2_norm, Linf_norm) for density field.
"""
function check_hydrostatic_drift(Q::Array{Float64,3}, Q_bar::Array{Float64,3},
                                  grid::Grid2D)
    L2 = 0.0
    Linf = 0.0
    count = 0
    for j in phys_z(grid), i in phys_x(grid)
        rho = Q[1, i, j]
        rho_bar = Q_bar[1, i, j]
        if rho_bar > 1e-30
            rel_err = abs(rho - rho_bar) / rho_bar
            L2 += rel_err^2
            Linf = max(Linf, rel_err)
            count += 1
        end
    end
    L2 = sqrt(L2 / max(count, 1))
    return (L2, Linf)
end

"""
    cfl_monitor(Q, grid, gamma)

Compute maximum CFL number across the domain.
"""
function cfl_monitor(Q::Array{Float64,3}, grid::Grid2D, gamma::Float64, R_gas::Float64)
    cfl_max = 0.0
    for j in phys_z(grid), i in phys_x(grid)
        Qi = SVector(Q[1,i,j], Q[2,i,j], Q[3,i,j], Q[4,i,j])
        rho, u, w, P = primitive_from_conservative(Qi, gamma)
        cs = sqrt(gamma * P / rho)
        cfl_x = (abs(u) + cs) / grid.dx
        cfl_z = (abs(w) + cs) / grid.dz
        cfl_max = max(cfl_max, cfl_x, cfl_z)
    end
    return cfl_max  # multiply by dt to get actual CFL number
end

"""
    check_positivity(Q, grid)

Check that density and pressure remain positive. Returns (min_rho, min_P).
"""
function check_positivity(Q::Array{Float64,3}, grid::Grid2D, gamma::Float64)
    min_rho = Inf
    min_P = Inf
    for j in phys_z(grid), i in phys_x(grid)
        rho = Q[1, i, j]
        KE = 0.5 * (Q[2,i,j]^2 + Q[3,i,j]^2) / rho
        P = (gamma - 1.0) * (Q[4,i,j] - KE)
        min_rho = min(min_rho, rho)
        min_P = min(min_P, P)
    end
    return (min_rho, min_P)
end

"""
    print_diagnostics(step, t, dt, Q, Q_bar, grid, config)

Print a one-line diagnostic summary.
"""
function print_diagnostics(step::Int, t::Float64, dt::Float64,
                           Q::Array{Float64,3}, Q_bar::Array{Float64,3},
                           grid::Grid2D, config::SimConfig)
    E_total = compute_total_energy(Q, grid)
    KE = compute_kinetic_energy(Q, grid)
    L2, Linf = check_hydrostatic_drift(Q, Q_bar, grid)
    min_rho, min_P = check_positivity(Q, grid, config.gamma)

    @printf("Step %7d | t=%10.3f s | dt=%8.2e s | E=%12.5e | KE=%12.5e | drho_L2=%8.2e | min_rho=%8.2e | min_P=%8.2e\n",
            step, t, dt, E_total, KE, L2, min_rho, min_P)
end
