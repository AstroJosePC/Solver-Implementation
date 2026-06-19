# =============================================================================
# Crank-Nicholson implicit diffusion solver
# Srivastava et al. (2022), Equations 12-16
#
# Solves: du/dt = div(nu(z) * grad(u)) implicitly with altitude-dependent
# kinematic viscosity nu(z) = mu(T) / rho_bar(z).
#
# Using a domain-average nu is WRONG for stratified atmospheres:
# velocity amplitude grows as exp(z/2H), so uniform-nu diffusion transfers
# energy from near-vacuum top cells to dense bottom cells, injecting KE.
# =============================================================================

"""
    DiffusionSolver

Pre-assembled sparse matrices and workspace for the implicit diffusion solve.
Stores altitude-dependent viscosity profile and cached LU factorizations.
"""
mutable struct DiffusionSolver
    nx::Int
    nz::Int
    N::Int                      # total unknowns = nx * nz
    A_lu_nu::Any                # cached LU for momentum diffusion
    A_lu_alpha::Any             # cached LU for thermal diffusion
    rhs::Vector{Float64}
    nu_z::Vector{Float64}       # kinematic viscosity at each z-level [nz]
    alpha_z::Vector{Float64}    # thermal diffusivity at each z-level [nz]
    dt_ref::Float64             # dt used to build current LU
    initialized::Bool           # whether nu_z/alpha_z have been computed
end

function DiffusionSolver(nx::Int, nz::Int)
    N = nx * nz
    DiffusionSolver(nx, nz, N, nothing, nothing, zeros(N),
                    zeros(nz), zeros(nz), 0.0, false)
end

"""
    linear_index(i, j, nx)

Map 2D physical index (i,j) to 1D index for the linear system.
i in 1:nx, j in 1:nz.
"""
linear_index(i::Int, j::Int, nx::Int) = (j - 1) * nx + i

"""
    build_cn_matrix(nx, nz, dx, dz, kappa_z, dt)

Build the Crank-Nicholson left-hand side matrix:
  (I - dt/2 * D) where D is the 2D diffusion operator with altitude-dependent
  coefficient kappa_z[j].

x-direction: kappa_z[j] * d²u/dx² (periodic)
z-direction: d/dz(kappa_z * du/dz) in divergence form (Neumann at boundaries)

Returns sparse matrix A of size (nx*nz) x (nx*nz).
"""
function build_cn_matrix(nx::Int, nz::Int, dx::Float64, dz::Float64,
                         kappa_z::Vector{Float64}, dt::Float64)
    N = nx * nz
    c = dt / 2.0

    I_idx = Int[]
    J_idx = Int[]
    V_val = Float64[]
    sizehint!(I_idx, 5 * N)
    sizehint!(J_idx, 5 * N)
    sizehint!(V_val, 5 * N)

    for j in 1:nz, i in 1:nx
        row = linear_index(i, j, nx)

        # x-direction coefficient (uniform per z-level)
        fx = c * kappa_z[j] / dx^2

        # z-direction face-centered coefficients (divergence form)
        if j < nz
            kappa_top_face = 0.5 * (kappa_z[j] + kappa_z[j+1])
        else
            kappa_top_face = kappa_z[j]  # Neumann: ghost = interior
        end
        if j > 1
            kappa_bot_face = 0.5 * (kappa_z[j] + kappa_z[j-1])
        else
            kappa_bot_face = kappa_z[j]  # Neumann: ghost = interior
        end
        fz_top = c * kappa_top_face / dz^2
        fz_bot = c * kappa_bot_face / dz^2

        # Diagonal
        diag_val = 1.0 + 2.0 * fx + fz_top + fz_bot
        # Neumann adjustments: remove ghost-cell contributions
        if j == 1
            diag_val -= fz_bot  # ghost = interior cancels the bot term
        end
        if j == nz
            diag_val -= fz_top  # ghost = interior cancels the top term
        end

        push!(I_idx, row); push!(J_idx, row); push!(V_val, diag_val)

        # x-neighbors (periodic)
        i_left  = mod1(i - 1, nx)
        i_right = mod1(i + 1, nx)
        push!(I_idx, row); push!(J_idx, linear_index(i_left, j, nx)); push!(V_val, -fx)
        push!(I_idx, row); push!(J_idx, linear_index(i_right, j, nx)); push!(V_val, -fx)

        # z-neighbors
        if j > 1
            push!(I_idx, row); push!(J_idx, linear_index(i, j-1, nx)); push!(V_val, -fz_bot)
        end
        if j < nz
            push!(I_idx, row); push!(J_idx, linear_index(i, j+1, nx)); push!(V_val, -fz_top)
        end
    end

    return sparse(I_idx, J_idx, V_val, N, N)
end

"""
    build_cn_rhs!(rhs, field, nx, nz, dx, dz, kappa_z, dt)

Build the Crank-Nicholson right-hand side vector:
  (I + dt/2 * D) * u^n
with altitude-dependent diffusion coefficient kappa_z[j].
"""
function build_cn_rhs!(rhs::Vector{Float64}, field::AbstractMatrix{Float64},
                       nx::Int, nz::Int, dx::Float64, dz::Float64,
                       kappa_z::Vector{Float64}, dt::Float64)
    c = dt / 2.0

    for j in 1:nz, i in 1:nx
        row = linear_index(i, j, nx)

        u_c = field[i, j]

        # x-neighbors (periodic)
        u_left  = field[mod1(i - 1, nx), j]
        u_right = field[mod1(i + 1, nx), j]

        # z-neighbors (Neumann: ghost = interior value)
        u_below = j > 1  ? field[i, j-1] : field[i, j]
        u_above = j < nz ? field[i, j+1] : field[i, j]

        # x-Laplacian with local kappa
        lap_x = kappa_z[j] * (u_left - 2.0 * u_c + u_right) / dx^2

        # z-divergence form: d/dz(kappa * du/dz)
        if j < nz
            kappa_top_face = 0.5 * (kappa_z[j] + kappa_z[j+1])
        else
            kappa_top_face = kappa_z[j]
        end
        if j > 1
            kappa_bot_face = 0.5 * (kappa_z[j] + kappa_z[j-1])
        else
            kappa_bot_face = kappa_z[j]
        end
        lap_z = (kappa_top_face * (u_above - u_c) - kappa_bot_face * (u_c - u_below)) / dz^2

        rhs[row] = u_c + c * (lap_x + lap_z)
    end
    return nothing
end

"""
    diffusion_step!(Q, Q_bar, grid, dt, config, solver)

Apply one implicit diffusion step to velocity (u, w) and temperature fields.
Uses altitude-dependent viscosity nu(z) = mu(T) / rho_bar(z).

On first call, computes and caches nu_z from the background state Q_bar.
LU factorizations are rebuilt when dt changes by more than 1%.
"""
function diffusion_step!(Q::Array{Float64,3}, Q_bar::Array{Float64,3},
                         grid::Grid2D, dt::Float64,
                         config::SimConfig, solver::DiffusionSolver)
    nx = grid.nx
    nz = grid.nz
    ng = grid.ng

    # Compute altitude-dependent viscosity on first call
    if !solver.initialized
        T0 = config.T0
        for j in 1:nz
            jg = j + ng
            # Use x-averaged background density (should be uniform in x for isothermal)
            rho_bar_j = Q_bar[1, ng + 1, jg]
            solver.nu_z[j] = kinematic_viscosity(T0, rho_bar_j, config.planet)
            solver.alpha_z[j] = thermal_diffusivity(solver.nu_z[j], config.Pr)
        end
        solver.initialized = true
    end

    # Rebuild LU factorizations if dt changed significantly
    if solver.A_lu_nu === nothing || abs(dt - solver.dt_ref) / max(dt, 1e-30) > 0.01
        A_nu = build_cn_matrix(nx, nz, grid.dx, grid.dz, solver.nu_z, dt)
        solver.A_lu_nu = lu(A_nu)

        A_alpha = build_cn_matrix(nx, nz, grid.dx, grid.dz, solver.alpha_z, dt)
        solver.A_lu_alpha = lu(A_alpha)

        solver.dt_ref = dt
    end

    # Extract physical domain fields (perturbations only)
    u_pert_field = zeros(nx, nz)
    w_field = zeros(nx, nz)
    T_pert_field = zeros(nx, nz)
    rho_field = zeros(nx, nz)
    
    u_bar_field = zeros(nx, nz)
    T_bar_field = zeros(nx, nz)

    for j in 1:nz, i in 1:nx
        ig = i + ng
        jg = j + ng
        Qi = SVector(Q[1,ig,jg], Q[2,ig,jg], Q[3,ig,jg], Q[4,ig,jg])
        rho, u, w, P = primitive_from_conservative(Qi, config.gamma)
        
        Qbar_i = SVector(Q_bar[1,ig,jg], Q_bar[2,ig,jg], Q_bar[3,ig,jg], Q_bar[4,ig,jg])
        rho_b, u_b, w_b, P_b = primitive_from_conservative(Qbar_i, config.gamma)
        T_b = P_b / (rho_b * config.R_gas)
        
        u_bar_field[i, j] = u_b
        T_bar_field[i, j] = T_b
        
        u_pert_field[i, j] = u - u_b
        w_field[i, j] = w   # w_bar is assumed 0
        T_pert_field[i, j] = P / (rho * config.R_gas) - T_b
        rho_field[i, j] = rho
    end

    # Solve for u' (momentum diffusion)
    build_cn_rhs!(solver.rhs, u_pert_field, nx, nz, grid.dx, grid.dz, solver.nu_z, dt)
    u_pert_new = solver.A_lu_nu \ solver.rhs

    # Solve for w (momentum diffusion)
    build_cn_rhs!(solver.rhs, w_field, nx, nz, grid.dx, grid.dz, solver.nu_z, dt)
    w_new = solver.A_lu_nu \ solver.rhs

    # Solve for T' (thermal diffusion)
    build_cn_rhs!(solver.rhs, T_pert_field, nx, nz, grid.dx, grid.dz, solver.alpha_z, dt)
    T_pert_new = solver.A_lu_alpha \ solver.rhs

    # Reconstruct conservative variables from (background + diffused perturbation)
    for j in 1:nz, i in 1:nx
        ig = i + ng
        jg = j + ng
        idx = linear_index(i, j, nx)

        rho = rho_field[i, j]
        u = u_bar_field[i, j] + u_pert_new[idx]
        w = w_new[idx]
        T = T_bar_field[i, j] + T_pert_new[idx]
        P = rho * config.R_gas * T
        
        # Prevent numerical negative pressures
        P = max(P, Q_bar[4, ig, jg] * (config.gamma - 1.0) * 1e-4)
        
        E = P / (config.gamma - 1.0) + 0.5 * rho * (u^2 + w^2)

        Q[1, ig, jg] = rho
        Q[2, ig, jg] = rho * u
        Q[3, ig, jg] = rho * w
        Q[4, ig, jg] = E
    end

    return nothing
end
