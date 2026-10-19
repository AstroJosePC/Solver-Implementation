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
    background_pressure_energy(Q_bar, i, j, pbar_mode)

The energy whose (gamma - 1) multiple is taken as the background pressure P_bar:
Q_bar[4] itself for :total, Q_bar[4] minus the background kinetic energy for
:pressure. Returning Q_bar[4] unchanged for :total keeps earlier runs bit-identical.
"""
@inline function background_pressure_energy(Q_bar::Array{Float64,3}, i::Int, j::Int,
                                            pbar_mode::Symbol)
    pbar_mode === :total && return Q_bar[4, i, j]
    return Q_bar[4, i, j] - 0.5 * (Q_bar[2, i, j]^2 + Q_bar[3, i, j]^2) / Q_bar[1, i, j]
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
    z_predictor(Q, Q_bar, sponge_K, i, j, dt, dtdz, gamma, g, pbar_mode)

Half state of Eq. 10a on the face j+1/2 of column i:
Q_{j+1/2}^{n+1/2} = 0.5*(Q_j + Q_{j+1}) - (dt/2dz)*(G_{j+1} - G_j) + (dt/2)*S_avg,
with S the gravity source plus the in-step sponge K*(Q_bar - Q).
"""
@inline function z_predictor(Q::Array{Float64,3}, Q_bar::Array{Float64,3},
                             sponge_K::Matrix{Float64}, i::Int, j::Int, dt::Float64,
                             dtdz::Float64, gamma::Float64, g::Float64, pbar_mode::Symbol)
    Qj   = SVector(Q[1,i,j], Q[2,i,j], Q[3,i,j], Q[4,i,j])
    Qjp1 = SVector(Q[1,i,j+1], Q[2,i,j+1], Q[3,i,j+1], Q[4,i,j+1])
    P_bar_j   = (gamma - 1.0) * background_pressure_energy(Q_bar, i, j, pbar_mode)
    P_bar_jp1 = (gamma - 1.0) * background_pressure_energy(Q_bar, i, j+1, pbar_mode)
    Gj   = compute_flux_z(Qj, P_bar_j, gamma)
    Gjp1 = compute_flux_z(Qjp1, P_bar_jp1, gamma)
    Sj   = compute_source(Qj, Q_bar[1, i, j], g)
    Sjp1 = compute_source(Qjp1, Q_bar[1, i, j+1], g)
    Kj   = sponge_K[i, j]
    Kjp1 = sponge_K[i, j+1]
    return SVector(ntuple(Val(4)) do k
        sp_j   = Kj   * (Q_bar[k, i, j]   - Qj[k])
        sp_jp1 = Kjp1 * (Q_bar[k, i, j+1] - Qjp1[k])
        0.5 * (Qj[k] + Qjp1[k]) -
        0.5 * dtdz * (Gjp1[k] - Gj[k]) +
        0.25 * dt * (Sj[k] + Sjp1[k] + sp_j + sp_jp1)
    end)
end

"""
    x_sweep!(Q, Q_half_x, grid, dt, gamma)

Richtmyer Lax-Wendroff x-sweep (Eqs. 9a, 9b) over the physical rows.
"""
function x_sweep!(Q::Array{Float64,3}, Q_half_x::Array{Float64,3}, grid::Grid2D,
                  dt::Float64, gamma::Float64)
    dtdx = dt / grid.dx
    ix = phys_x(grid)
    iz = phys_z(grid)

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
    return nothing
end

"""
    z_sweep!(Q, Q_bar, Q_half_z, grid, dt, gamma, g, sponge_K; pbar_mode, w_face)

Richtmyer Lax-Wendroff z-sweep with gravity and sponge sources (Eqs. 10a, 10b).
With `w_face` (the :face path) the bottom ghost rows are first rebuilt from the
current state (reflective mirror, then the forced w), and the ground-face half state
of the first physical row's corrector takes w = w_face with its density, u and
pressure kept.
"""
function z_sweep!(Q::Array{Float64,3}, Q_bar::Array{Float64,3},
                  Q_half_z::Array{Float64,3}, grid::Grid2D, dt::Float64,
                  gamma::Float64, g::Float64, sponge_K::Matrix{Float64};
                  pbar_mode::Symbol = :total,
                  w_face::Union{Nothing, Vector{Float64}} = nothing)
    dtdz = dt / grid.dz
    ix = phys_x(grid)
    iz = phys_z(grid)

    # Face forcing: ghost rows consistent with the state entering the z-sweep
    if w_face !== nothing
        apply_bottom_bc!(Q, Q_bar, grid)
        for i in ix
            impose_ghost_w!(Q, i, grid.ng, gamma, w_face[i - grid.ng])
        end
    end

    # Predictor (Eq. 10a)
    for j in (first(iz)-1):(last(iz))
        for i in ix
            q = z_predictor(Q, Q_bar, sponge_K, i, j, dt, dtdz, gamma, g, pbar_mode)
            for k in 1:4
                Q_half_z[k, i, j] = q[k]
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
            if w_face !== nothing && j == first(iz)
                # Ground face: impose w, keep the predictor's density, u and pressure
                rho_f, u_f, _, P_f = primitive_from_conservative(Qh_bot, gamma)
                wf = w_face[i - grid.ng]
                Qh_bot = SVector(rho_f, rho_f * u_f, rho_f * wf,
                                 P_f / (gamma - 1.0) + 0.5 * rho_f * (u_f^2 + wf^2))
            end
            Eb_jm1 = background_pressure_energy(Q_bar, i, j-1, pbar_mode)
            Eb_j   = background_pressure_energy(Q_bar, i, j,   pbar_mode)
            Eb_jp1 = background_pressure_energy(Q_bar, i, j+1, pbar_mode)
            P_bar_top = 0.5 * (gamma - 1.0) * (Eb_j + Eb_jp1)
            P_bar_bot = 0.5 * (gamma - 1.0) * (Eb_jm1 + Eb_j)
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

"""
    FaceProbe(grid)

Workspace and per-column record of the ground face (z = 0) as the z-predictor sees it
(SimConfig.face_probe; Section 14 addendum 2). `wf` is set by the caller to the forced
w of the step. probe_ground_face! fills, for the state entering the z-sweep:
`w_ghost`, `p_ghost`: w and p' (minus the face P_bar) of the predictor half state with
the ghost rows as they stand (the :ghost path); `w_avg`: w of the plain average of the
ghost and first physical cells; `w_rebuilt`, `p_rebuilt`: the same predictor with the
ghost rows rebuilt from the current state and wf (step 4 of the :face path).
"""
struct FaceProbe
    wf::Vector{Float64}
    w_ghost::Vector{Float64}
    w_avg::Vector{Float64}
    w_rebuilt::Vector{Float64}
    p_ghost::Vector{Float64}
    p_rebuilt::Vector{Float64}
    Qc::Array{Float64,3}
end
FaceProbe(grid::Grid2D) = FaceProbe([zeros(grid.nx) for _ in 1:6]...,
                                    zeros(4, grid.nx_total, grid.nz_total))

function probe_ground_face!(pr::FaceProbe, Q::Array{Float64,3}, Q_bar::Array{Float64,3},
                            grid::Grid2D, dt::Float64, gamma::Float64, g::Float64,
                            sponge_K::Matrix{Float64}, pbar_mode::Symbol)
    ng = grid.ng
    dtdz = dt / grid.dz
    copyto!(pr.Qc, Q)
    apply_bottom_bc!(pr.Qc, Q_bar, grid)
    for i in phys_x(grid)
        impose_ghost_w!(pr.Qc, i, ng, gamma, pr.wf[i - ng])
    end
    for i in phys_x(grid)
        c = i - ng
        P_bar_face = 0.5 * (gamma - 1.0) * (background_pressure_energy(Q_bar, i, ng, pbar_mode) +
                                            background_pressure_energy(Q_bar, i, ng + 1, pbar_mode))
        qa = z_predictor(Q, Q_bar, sponge_K, i, ng, dt, dtdz, gamma, g, pbar_mode)
        qb = z_predictor(pr.Qc, Q_bar, sponge_K, i, ng, dt, dtdz, gamma, g, pbar_mode)
        _, _, wa, Pa = primitive_from_conservative(qa, gamma)
        _, _, wb, Pb = primitive_from_conservative(qb, gamma)
        pr.w_ghost[c], pr.p_ghost[c] = wa, Pa - P_bar_face
        pr.w_rebuilt[c], pr.p_rebuilt[c] = wb, Pb - P_bar_face
        pr.w_avg[c] = (Q[3, i, ng] + Q[3, i, ng + 1]) / (Q[1, i, ng] + Q[1, i, ng + 1])
    end
    return nothing
end

"""
    euler_step!(Q, Q_bar, Q_half_x, Q_half_z, grid, dt, gamma, g, sponge_K;
                pbar_mode, w_face, order, probe)

Perform one full Euler timestep using dimensionally-split Richtmyer Lax-Wendroff:
the x-sweep (Eqs. 9a, 9b) then the z-sweep (Eqs. 10a, 10b), or the reverse with
`order = :zx` (SimConfig.sweep_order).

Sponge damping K*(Q_bar - Q) is included as a source term inside the LW integration
(both predictor and corrector); the :split sponge of SimConfig.sponge_step is applied
by the caller after this step, with a zero `sponge_K` passed here.
Q_half_x and Q_half_z are pre-allocated workspace arrays of shape (4, nx_total, nz_total).

`w_face` (length nx, SimConfig.bottom_forcing = :face) is the imposed ground w at
t^{n+1/2}; see z_sweep!. `nothing` (the default) leaves the ghost rows as set before
the step. `probe` (a FaceProbe) records the ground face of the state entering the
z-sweep without changing the solution.

`order = :strang` runs x(dt/2) z(dt) x(dt/2). `refresh_x` (a function of Q, e.g. the
periodic x boundary) is called before every x-sweep that follows another sweep; it is
required for :strang and, as `nothing`, keeps the :zx record, whose x-sweep reads the
x ghosts set before the z-sweep.
"""
function euler_step!(Q::Array{Float64,3},
                     Q_bar::Array{Float64,3},
                     Q_half_x::Array{Float64,3},
                     Q_half_z::Array{Float64,3},
                     grid::Grid2D, dt::Float64,
                     gamma::Float64, g::Float64,
                     sponge_K::Matrix{Float64};
                     pbar_mode::Symbol = :total,
                     w_face::Union{Nothing, Vector{Float64}} = nothing,
                     order::Symbol = :xz,
                     probe::Union{Nothing, FaceProbe} = nothing,
                     refresh_x::Union{Nothing, Function} = nothing)
    if order === :strang
        refresh_x === nothing && error("order = :strang needs refresh_x")
        x_sweep!(Q, Q_half_x, grid, 0.5 * dt, gamma)
        refresh_x(Q)
    end
    order === :xz && x_sweep!(Q, Q_half_x, grid, dt, gamma)
    probe === nothing ||
        probe_ground_face!(probe, Q, Q_bar, grid, dt, gamma, g, sponge_K, pbar_mode)
    z_sweep!(Q, Q_bar, Q_half_z, grid, dt, gamma, g, sponge_K;
             pbar_mode = pbar_mode, w_face = w_face)
    order === :xz || refresh_x === nothing || refresh_x(Q)
    order === :zx && x_sweep!(Q, Q_half_x, grid, dt, gamma)
    order === :strang && x_sweep!(Q, Q_half_x, grid, 0.5 * dt, gamma)
    return nothing
end
