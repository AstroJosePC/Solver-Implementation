# =============================================================================
# Hyperbolic solver: Dimensionally-split 2-step Richtmyer Lax-Wendroff
# Srivastava et al. (2022), Equations 9-10
# =============================================================================

"""
    compute_flux_x(Q, gamma)

Compute x-direction flux vector F from conservative state Q.
F = [rho*u, rho*u^2 + P, rho*u*w, u*(E+P)]
"""
function compute_flux_x(Q::AbstractVector, gamma::Float64)
    rho, u, w, P = primitive_from_conservative(Q, gamma)
    E = Q[4]
    return SVector(rho * u,
                   rho * u^2 + P,
                   rho * u * w,
                   u * (E + P))
end

"""
    compute_flux_z(Q, P_bar, gamma)

Compute z-direction flux vector G from conservative state Q.
G = [rho*w, rho*u*w, rho*w^2 + P, w*(E+P)]
"""
function compute_flux_z(Q::AbstractVector, P_bar::Float64, gamma::Float64)
    rho, u, w, P = primitive_from_conservative(Q, gamma)
    E = Q[4]
    return SVector(rho * w,
                   rho * u * w,
                   rho * w^2 + (P - P_bar),
                   w * (E + P))
end

"""
    compute_source(Q, rho_bar, g)

Compute source vector S for gravitational terms.
S = [0, 0, -rho*g, -rho*g*w]
"""
function compute_source(Q::AbstractVector, rho_bar::Float64, g::Float64)
    rho = Q[1]
    w   = Q[3] / rho
    return SVector(0.0, 0.0, -(rho - rho_bar) * g, -rho * g * w)
end

"""
    compute_adaptive_dt(Q, grid, CFL, gamma)

Compute adaptive timestep from CFL condition.
dt = CFL * min(dx, dz) / max(|u| + cs, |w| + cs)
"""
function compute_adaptive_dt(Q::Array{Float64,3}, grid::Grid2D,
                              CFL::Float64, gamma::Float64, R_gas::Float64)
    c_max = 0.0
    for j in phys_z(grid), i in phys_x(grid)
        Qi = SVector(Q[1,i,j], Q[2,i,j], Q[3,i,j], Q[4,i,j])
        rho, u, w, P = primitive_from_conservative(Qi, gamma)
        # Floor P to prevent sqrt(negative) from floating point errors at high altitude
        cs = sqrt(gamma * max(P, 0.0) / rho)
        c_local = max(abs(u) + cs, abs(w) + cs)
        c_max = max(c_max, c_local)
    end
    if c_max < 1e-15
        return 1e-6  # fallback for zero-velocity state
    end
    return CFL * min(grid.dx, grid.dz) / c_max
end

"""
    euler_step!(Q, Q_bar, Q_half_x, Q_half_z, grid, dt, gamma, g, sponge_K)

Perform one full Euler timestep using dimensionally-split Richtmyer Lax-Wendroff.

Step 1: X-sweep (Eqs. 9a, 9b) - updates Q with x-fluxes
Step 2: Z-sweep (Eqs. 10a, 10b) - updates Q with z-fluxes, source terms, and sponge

Sponge damping K*(Q_bar - Q) is included as a source term inside the LW integration
(both predictor and corrector) rather than applied as a separate operator-split step.
Q_half_x and Q_half_z are pre-allocated workspace arrays of shape (4, nx_total, nz_total).
"""
function euler_step!(Q::Array{Float64,3},
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
    # X-SWEEP: Richtmyer Lax-Wendroff (Eq. 9a, 9b)
    # =========================================================================

    # Predictor (Eq. 9a): compute half-step values at i+1/2 interfaces
    # Q_{i+1/2}^{n+1/2} = 0.5*(Q_i + Q_{i+1}) - (dt/2dx)*(F_{i+1} - F_i)
    for j in iz
        for i in (first(ix)-1):(last(ix))  # need i and i+1 in physical range
            Qi   = SVector(Q[1,i,j], Q[2,i,j], Q[3,i,j], Q[4,i,j])
            Qip1 = SVector(Q[1,i+1,j], Q[2,i+1,j], Q[3,i+1,j], Q[4,i+1,j])
            Fi   = compute_flux_x(Qi, gamma)
            Fip1 = compute_flux_x(Qip1, gamma)

            for k in 1:4
                Q_half_x[k, i, j] = 0.5 * (Qi[k] + Qip1[k]) - 0.5 * dtdx * (Fip1[k] - Fi[k])
            end
        end
    end

    # Corrector (Eq. 9b): update Q using half-step fluxes
    # Q_i^{n+1} = Q_i^n - (dt/dx)*(F_{i+1/2} - F_{i-1/2})
    for j in iz
        for i in ix
            Qh_right = SVector(Q_half_x[1,i,j], Q_half_x[2,i,j],
                               Q_half_x[3,i,j], Q_half_x[4,i,j])
            Qh_left  = SVector(Q_half_x[1,i-1,j], Q_half_x[2,i-1,j],
                               Q_half_x[3,i-1,j], Q_half_x[4,i-1,j])
            Fh_right = compute_flux_x(Qh_right, gamma)
            Fh_left  = compute_flux_x(Qh_left, gamma)

            for k in 1:4
                Q[k, i, j] -= dtdx * (Fh_right[k] - Fh_left[k])
            end
        end
    end

    # =========================================================================
    # Z-SWEEP: Richtmyer Lax-Wendroff with source terms (Eq. 10a, 10b)
    # =========================================================================

    # Predictor (Eq. 10a):
    # Q_{j+1/2}^{n+1/2} = 0.5*(Q_j + Q_{j+1}) - (dt/2dz)*(G_{j+1} - G_j) + (dt/2)*S_avg
    for j in (first(iz)-1):(last(iz))
        for i in ix
            Qj   = SVector(Q[1,i,j], Q[2,i,j], Q[3,i,j], Q[4,i,j])
            Qjp1 = SVector(Q[1,i,j+1], Q[2,i,j+1], Q[3,i,j+1], Q[4,i,j+1])
            P_bar_j   = (gamma - 1.0) * Q_bar[4, i, j]
            P_bar_jp1 = (gamma - 1.0) * Q_bar[4, i, j+1]
            Gj   = compute_flux_z(Qj, P_bar_j, gamma)
            Gjp1 = compute_flux_z(Qjp1, P_bar_jp1, gamma)
            Sj   = compute_source(Qj, Q_bar[1, i, j], g)
            Sjp1 = compute_source(Qjp1, Q_bar[1, i, j+1], g)
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

    # Corrector (Eq. 10b):
    # Q_j^{n+1} = Q_j^n - (dt/dz)*(G_{j+1/2} - G_{j-1/2}) + dt*S^{n+1/2}
    # S^{n+1/2} is evaluated from predictor half-step values
    for j in iz
        for i in ix
            Qh_top = SVector(Q_half_z[1,i,j], Q_half_z[2,i,j],
                             Q_half_z[3,i,j], Q_half_z[4,i,j])
            Qh_bot = SVector(Q_half_z[1,i,j-1], Q_half_z[2,i,j-1],
                             Q_half_z[3,i,j-1], Q_half_z[4,i,j-1])
            P_bar_top = 0.5 * (gamma - 1.0) * (Q_bar[4,i,j] + Q_bar[4,i,j+1])
            P_bar_bot = 0.5 * (gamma - 1.0) * (Q_bar[4,i,j-1] + Q_bar[4,i,j])
            rho_bar_top = 0.5 * (Q_bar[1,i,j] + Q_bar[1,i,j+1])
            rho_bar_bot = 0.5 * (Q_bar[1,i,j-1] + Q_bar[1,i,j])
            Gh_top = compute_flux_z(Qh_top, P_bar_top, gamma)
            Gh_bot = compute_flux_z(Qh_bot, P_bar_bot, gamma)

            # Source at half-step: average of predictor boundary values
            Sh_top = compute_source(Qh_top, rho_bar_top, g)
            Sh_bot = compute_source(Qh_bot, rho_bar_bot, g)

            # Sponge source at half-nodes (j+1/2 and j-1/2)
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
