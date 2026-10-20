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

Pre-assembled operators and workspace for the implicit Crank-Nicolson diffusion solve
of Srivastava et al. (2022), Eqs. 12-16. Stores the altitude-dependent diffusivity
profiles and the factorized left-hand sides of the method chosen by
SimConfig.diffusion_solver:

- `:lu`: global sparse LU (UMFPACK) of the (nx*nz) x (nx*nz) operator, rebuilt from
  scratch at every refactorization (every record so far).
- `:lu_reuse`: the same global matrices, refactorized in place with `lu!`, which keeps
  the symbolic analysis of the first factorization (the sparsity pattern is checked to
  be identical at every rebuild).
- `:fft_tridiag`: the same discrete operator diagonalized by a real FFT in x, one real
  tridiagonal system in z per horizontal mode (FFTTridiag).

Fields: `nu_z`, `alpha_z` [nz] in m^2/s; `dt_ref` [s] is the step size the current
factorization was built with; `rhs` [nx*nz] with linear index (j - 1) nx + i;
`n_factor` counts the (re)factorizations; `kind` is the method fixed at the first call.
Field buffers reused by every diffusion_step! call, each [nx, nz] (x first, physical
cells only): `u_pert` u' = u - u_bar and `w_pert` w in m/s, `T_pert` T' = T - T_bar
in K, `rho_f` rho in kg/m^3, `u_bar` in m/s and `T_bar` in K of the background Q_bar.
They are separate arrays from `rhs` and the FFTTridiag workspace (`out`, `rhs_mat`).
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
    kind::Symbol                # SimConfig.diffusion_solver, fixed at the first call
    n_factor::Int               # number of (re)factorizations of the CN left-hand sides
    pattern::Any                # :lu_reuse: (colptr, rowval) of the first global matrix
    fft::Any                    # :fft_tridiag: FFTTridiag workspace
    u_pert::Matrix{Float64}     # u' [nx, nz]
    w_pert::Matrix{Float64}     # w [nx, nz]
    T_pert::Matrix{Float64}     # T' [nx, nz]
    rho_f::Matrix{Float64}      # rho [nx, nz]
    u_bar::Matrix{Float64}      # background u [nx, nz]
    T_bar::Matrix{Float64}      # background T [nx, nz]
end

function DiffusionSolver(nx::Int, nz::Int)
    N = nx * nz
    DiffusionSolver(nx, nz, N, nothing, nothing, zeros(N),
                    zeros(nz), zeros(nz), 0.0, false, :lu, 0, nothing, nothing,
                    zeros(nx, nz), zeros(nx, nz), zeros(nx, nz),
                    zeros(nx, nz), zeros(nx, nz), zeros(nx, nz))
end

"""Implicit diffusion methods accepted by SimConfig.diffusion_solver."""
const DIFFUSION_SOLVERS = (:lu, :lu_reuse, :fft_tridiag)

"""
    check_diffusion_config(config)

Fail loudly on an unknown SimConfig.diffusion_solver or a negative or non-finite
SimConfig.diffusion_dt_tol (relative dt change that triggers a refactorization).
"""
function check_diffusion_config(config::SimConfig)
    config.diffusion_solver in DIFFUSION_SOLVERS ||
        error("diffusion_solver must be one of $(DIFFUSION_SOLVERS), got :$(config.diffusion_solver)")
    isfinite(config.diffusion_dt_tol) && config.diffusion_dt_tol >= 0 ||
        error("diffusion_dt_tol must be finite and >= 0, got $(config.diffusion_dt_tol)")
    return nothing
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
Shapes: field [nx, nz] (x first), rhs [nx*nz] with linear index (j - 1) nx + i,
kappa_z [nz] in m^2/s; dx, dz in m, dt in s.
"""
function build_cn_rhs!(rhs::Vector{Float64}, field::AbstractMatrix{Float64},
                       nx::Int, nz::Int, dx::Float64, dz::Float64,
                       kappa_z::Vector{Float64}, dt::Float64)
    c = dt / 2.0
    size(field) == (nx, nz) && length(rhs) == nx * nz && length(kappa_z) == nz ||
        error("build_cn_rhs!: field $(size(field)), rhs $(length(rhs)), kappa_z $(length(kappa_z)) for nx = $nx, nz = $nz")

    @inbounds for j in 1:nz
        # z-divergence form face coefficients (Neumann: ghost = interior value)
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
        j_below = j > 1  ? j - 1 : j
        j_above = j < nz ? j + 1 : j

        for i in 1:nx
            row = linear_index(i, j, nx)

            u_c = field[i, j]

            # x-neighbors (periodic), mod1(i - 1, nx) and mod1(i + 1, nx)
            u_left  = field[i == 1 ? nx : i - 1, j]
            u_right = field[i == nx ? 1 : i + 1, j]

            # z-neighbors (Neumann: ghost = interior value)
            u_below = field[i, j_below]
            u_above = field[i, j_above]

            # x-Laplacian with local kappa
            lap_x = kappa_z[j] * (u_left - 2.0 * u_c + u_right) / dx^2

            # z-divergence form: d/dz(kappa * du/dz)
            lap_z = (kappa_top_face * (u_above - u_c) - kappa_bot_face * (u_c - u_below)) / dz^2

            rhs[row] = u_c + c * (lap_x + lap_z)
        end
    end
    return nothing
end

"""
    FFTTridiag

Workspace of the `:fft_tridiag` diffusion solve (SimConfig.diffusion_solver).

The CN left-hand side (I - dt/2 D) of build_cn_matrix is separable on the Cartesian
grid: its x part is kappa_z[j] times the periodic second difference, circulant at
each level j with a coefficient that depends on j only. A real FFT in x
diagonalizes it exactly. Mode m = 0, ..., nx/2 (array index k = m + 1) has the
eigenvalue -lam[k]/dx^2 of the periodic second difference, lam[k] = 2 - 2 cos(2 pi m/nx),
which leaves one real tridiagonal system in z per mode,

    -fz_bot[j] X[k, j-1] + (1 + fx[j] lam[k] + fz_top[j] + fz_bot[j]) X[k, j] - fz_top[j] X[k, j+1]

with fx = (dt/2) kappa_z[j]/dx^2, fz = (dt/2) kappa_face/dz^2, face kappa the mean of
the two neighbours, and the Neumann rows of build_cn_matrix at j = 1 and j = nz.
The systems are symmetric and strictly diagonally dominant, so the Thomas algorithm
needs no pivoting. The solution is the same discrete system as the global LU, equal
to it to round-off but not bit for bit. m = 0 (lam = 0) and, for even nx, the
Nyquist mode m = nx/2 (lam = 4) are real for real data and stay real, since every
coefficient is real.

Validity: uniform dx and dz with no metric terms, a diffusivity that depends on z only
(nu_z, alpha_z), and the periodic x coupling of build_cn_matrix, which the global LU
also applies whatever the hyperbolic x boundary is.

Fields: `lam` [nk]; Thomas coefficients `cp_*` (c'_j) and `inv_*` (1/denominator)
[nk, nz] and lower diagonals `lo_*` [nz] for kappa = nu (momentum) and alpha
(heat); `xhat` [nk, nz] complex spectral buffer; `out` three [nx*nz] solution
vectors (u', w, T') with the linear index of build_cn_rhs!, and `out_mat` their
[nx, nz] views; `pf`/`pb` the FFTW rfft/irfft plans along dimension 1 (x) of an
[nx, nz] array, nk = nx/2 + 1 (integer division).
"""
struct FFTTridiag{PF, PB}
    nk::Int
    lam::Vector{Float64}
    cp_nu::Matrix{Float64}
    inv_nu::Matrix{Float64}
    lo_nu::Vector{Float64}
    cp_alpha::Matrix{Float64}
    inv_alpha::Matrix{Float64}
    lo_alpha::Vector{Float64}
    xhat::Matrix{ComplexF64}
    out::NTuple{3, Vector{Float64}}
    out_mat::NTuple{3, Matrix{Float64}}
    rhs_mat::Matrix{Float64}
    pf::PF
    pb::PB
end

"""
    FFTTridiag(nx, nz, rhs)

Allocate the `:fft_tridiag` workspace for an nx x nz grid. `rhs` [nx*nz] is the
DiffusionSolver right-hand side, viewed as [nx, nz] (x first) by the forward plan.
Plans use FFTW.ESTIMATE, so the transform algorithm and its round-off are the same
in every run.
"""
function FFTTridiag(nx::Int, nz::Int, rhs::Vector{Float64})
    nk = nx ÷ 2 + 1
    lam = [2.0 - 2.0 * cos(2pi * (k - 1) / nx) for k in 1:nk]
    rhs_mat = reshape(rhs, nx, nz)
    xhat = zeros(ComplexF64, nk, nz)
    out = (zeros(nx * nz), zeros(nx * nz), zeros(nx * nz))
    out_mat = map(v -> reshape(v, nx, nz), out)
    pf = plan_rfft(rhs_mat, 1; flags = FFTW.ESTIMATE)
    pb = plan_irfft(xhat, nx, 1; flags = FFTW.ESTIMATE)
    FFTTridiag(nk, lam, zeros(nk, nz), zeros(nk, nz), zeros(nz),
               zeros(nk, nz), zeros(nk, nz), zeros(nz), xhat, out, out_mat, rhs_mat, pf, pb)
end

"""
    build_mode_tridiag!(cp, inv_den, lo, lam, nx, nz, dx, dz, kappa_z, dt)

Thomas-algorithm factorization of the per-mode tridiagonal CN systems (FFTTridiag)
for diffusivity kappa_z [nz] (m^2/s), spacings dx, dz (m) and step dt (s), with the
coefficients of build_cn_matrix (Srivastava et al. 2022, Eqs. 12-16). Writes, for
mode index k and level j, cp[k, j] = c_j / den and inv_den[k, j] = 1/den with
den = b_j - a_j cp[k, j-1], and lo[j] = a_j = -fz_bot[j] (lo[1] = 0). Arrays are
[nk, nz] with the mode index first.
"""
function build_mode_tridiag!(cp::Matrix{Float64}, inv_den::Matrix{Float64},
                             lo::Vector{Float64}, lam::Vector{Float64},
                             nx::Int, nz::Int, dx::Float64, dz::Float64,
                             kappa_z::Vector{Float64}, dt::Float64)
    c = dt / 2.0
    nk = length(lam)
    @inbounds for j in 1:nz
        fx = c * kappa_z[j] / dx^2
        kappa_top_face = j < nz ? 0.5 * (kappa_z[j] + kappa_z[j+1]) : kappa_z[j]
        kappa_bot_face = j > 1  ? 0.5 * (kappa_z[j] + kappa_z[j-1]) : kappa_z[j]
        fz_top = c * kappa_top_face / dz^2
        fz_bot = c * kappa_bot_face / dz^2
        bz = fz_top + fz_bot
        j == 1  && (bz -= fz_bot)   # Neumann: ghost = interior
        j == nz && (bz -= fz_top)
        a = j > 1  ? -fz_bot : 0.0
        up = j < nz ? -fz_top : 0.0
        lo[j] = a
        for k in 1:nk
            b = 1.0 + fx * lam[k] + bz
            den = j == 1 ? b : b - a * cp[k, j-1]
            inv_den[k, j] = 1.0 / den
            cp[k, j] = up / den
        end
    end
    return nothing
end

"""
    fft_tridiag_solve!(ws, which, slot)

Solve (I - dt_ref/2 D) x = rhs for the right-hand side held in ws.rhs_mat [nx, nz]
(the DiffusionSolver rhs): rfft in x, Thomas solve in z per mode for kappa = nu
(`which = :nu`) or alpha (`:alpha`), irfft. Returns the solution vector ws.out[slot]
[nx*nz], linear index (j - 1) nx + i.
"""
function fft_tridiag_solve!(ws::FFTTridiag, which::Symbol, slot::Int)
    cp, inv_den, lo = which === :nu ? (ws.cp_nu, ws.inv_nu, ws.lo_nu) :
                                      (ws.cp_alpha, ws.inv_alpha, ws.lo_alpha)
    X = ws.xhat
    nk, nz = size(X)
    mul!(X, ws.pf, ws.rhs_mat)
    @inbounds begin
        for k in 1:nk
            X[k, 1] *= inv_den[k, 1]
        end
        for j in 2:nz
            a = lo[j]
            for k in 1:nk
                X[k, j] = (X[k, j] - a * X[k, j-1]) * inv_den[k, j]
            end
        end
        for j in nz-1:-1:1
            for k in 1:nk
                X[k, j] -= cp[k, j] * X[k, j+1]
            end
        end
    end
    mul!(ws.out_mat[slot], ws.pb, X)   # c2r overwrites X, which is rebuilt next call
    return ws.out[slot]
end

"""
    factorize_cn!(solver, grid, dt)

(Re)build the CN left-hand sides (I - dt/2 D) for nu_z and alpha_z at step dt [s]
with the method in solver.kind, and count it in solver.n_factor.
"""
function factorize_cn!(solver::DiffusionSolver, grid::Grid2D, dt::Float64)
    nx, nz = solver.nx, solver.nz
    if solver.kind === :fft_tridiag
        ws = solver.fft
        build_mode_tridiag!(ws.cp_nu, ws.inv_nu, ws.lo_nu, ws.lam, nx, nz,
                            grid.dx, grid.dz, solver.nu_z, dt)
        build_mode_tridiag!(ws.cp_alpha, ws.inv_alpha, ws.lo_alpha, ws.lam, nx, nz,
                            grid.dx, grid.dz, solver.alpha_z, dt)
    elseif solver.kind === :lu_reuse && solver.A_lu_nu !== nothing
        A_nu = build_cn_matrix(nx, nz, grid.dx, grid.dz, solver.nu_z, dt)
        A_alpha = build_cn_matrix(nx, nz, grid.dx, grid.dz, solver.alpha_z, dt)
        colptr, rowval = solver.pattern
        for A in (A_nu, A_alpha)
            A.colptr == colptr && A.rowval == rowval ||
                error("diffusion_solver = :lu_reuse: CN matrix sparsity pattern changed")
        end
        lu!(solver.A_lu_nu, A_nu; reuse_symbolic = true)
        lu!(solver.A_lu_alpha, A_alpha; reuse_symbolic = true)
    else
        A_nu = build_cn_matrix(nx, nz, grid.dx, grid.dz, solver.nu_z, dt)
        solver.A_lu_nu = lu(A_nu)

        A_alpha = build_cn_matrix(nx, nz, grid.dx, grid.dz, solver.alpha_z, dt)
        solver.A_lu_alpha = lu(A_alpha)

        if solver.kind === :lu_reuse
            A_nu.colptr == A_alpha.colptr && A_nu.rowval == A_alpha.rowval ||
                error("diffusion_solver = :lu_reuse: nu and alpha CN patterns differ")
            solver.pattern = (copy(A_nu.colptr), copy(A_nu.rowval))
        end
    end
    solver.n_factor += 1
    return nothing
end

"""
    solve_cn!(solver, which, slot)

Solve the CN system for solver.rhs with the left-hand side of `which` (:nu or
:alpha). :lu and :lu_reuse return a new vector (LU backsolve), :fft_tridiag the
workspace vector of `slot`. The non-default methods fail on NaN or Inf.
"""
function solve_cn!(solver::DiffusionSolver, which::Symbol, slot::Int)
    if solver.kind === :lu
        return (which === :nu ? solver.A_lu_nu : solver.A_lu_alpha) \ solver.rhs
    end
    x = solver.kind === :fft_tridiag ? fft_tridiag_solve!(solver.fft, which, slot) :
        (which === :nu ? solver.A_lu_nu : solver.A_lu_alpha) \ solver.rhs
    all(isfinite, x) ||
        error("diffusion_solver = :$(solver.kind): non-finite value in the CN solution (slot $slot)")
    return x
end

"""
    diffusion_step!(Q, Q_bar, grid, dt, config, solver)

Apply one implicit Crank-Nicolson diffusion step (Srivastava et al. 2022, Eqs. 12-16)
to the perturbations u' = u - u_bar, w and T' = T - T_bar, with
nu(z) = mu(T0)/rho_bar(z) (m^2/s) for u', w and alpha = nu/Pr for T'; rho is
unchanged and P = rho R T is floored at 1e-4 of the background pressure.
Q, Q_bar [4, nx + 2ng, nz + 2ng] are conservative states (rho, rho u, rho w, E) in SI
units, dt in s.

On first call, computes and caches nu_z from the background state Q_bar and fixes
the method of SimConfig.diffusion_solver. The left-hand sides are built with the step
dt_ref of the last factorization and rebuilt when |dt - dt_ref|/dt exceeds
SimConfig.diffusion_dt_tol (default 0.01); the right-hand side always uses dt.
"""
function diffusion_step!(Q::Array{Float64,3}, Q_bar::Array{Float64,3},
                         grid::Grid2D, dt::Float64,
                         config::SimConfig, solver::DiffusionSolver)
    nx = grid.nx
    nz = grid.nz
    ng = grid.ng

    # Compute altitude-dependent viscosity on first call
    if !solver.initialized
        check_diffusion_config(config)
        T0 = config.T0
        for j in 1:nz
            jg = j + ng
            # Use x-averaged background density (should be uniform in x for isothermal)
            rho_bar_j = Q_bar[1, ng + 1, jg]
            solver.nu_z[j] = kinematic_viscosity(T0, rho_bar_j, config.planet)
            solver.alpha_z[j] = thermal_diffusivity(solver.nu_z[j], config.Pr)
        end
        solver.kind = config.diffusion_solver
        solver.kind === :fft_tridiag && (solver.fft = FFTTridiag(nx, nz, solver.rhs))
        solver.initialized = true
    end
    config.diffusion_solver === solver.kind ||
        error("diffusion_solver changed from :$(solver.kind) to :$(config.diffusion_solver)")

    # Rebuild the factorizations if dt changed by more than diffusion_dt_tol
    if solver.n_factor == 0 || abs(dt - solver.dt_ref) / max(dt, 1e-30) > config.diffusion_dt_tol
        factorize_cn!(solver, grid, dt)
        solver.dt_ref = dt
    end

    # Extract physical domain fields (perturbations only) into the solver buffers
    pack_perturbations!(solver, Q, Q_bar, grid, config)

    # Solve for u' (momentum diffusion)
    build_cn_rhs!(solver.rhs, solver.u_pert, nx, nz, grid.dx, grid.dz, solver.nu_z, dt)
    u_pert_new = solve_cn!(solver, :nu, 1)

    # Solve for w (momentum diffusion)
    build_cn_rhs!(solver.rhs, solver.w_pert, nx, nz, grid.dx, grid.dz, solver.nu_z, dt)
    w_new = solve_cn!(solver, :nu, 2)

    # Solve for T' (thermal diffusion)
    build_cn_rhs!(solver.rhs, solver.T_pert, nx, nz, grid.dx, grid.dz, solver.alpha_z, dt)
    T_pert_new = solve_cn!(solver, :alpha, 3)

    # Reconstruct conservative variables from (background + diffused perturbation)
    unpack_conservative!(Q, Q_bar, solver, u_pert_new, w_new, T_pert_new, grid, config)

    return nothing
end

"""
    pack_perturbations!(solver, Q, Q_bar, grid, config)

Fill the DiffusionSolver buffers [nx, nz] (x first, physical cells) from the
conservative states Q, Q_bar [4, nx + 2ng, nz + 2ng] (rho, rho u, rho w, E; SI):
u_pert = u - u_bar, w_pert = w (w_bar assumed 0), T_pert = P/(rho R) - T_bar with
T_bar = P_bar/(rho_bar R), rho_f = rho, and u_bar, T_bar of the background, with
P = (gamma - 1)(E - rho (u^2 + w^2)/2) (primitive_from_conservative).
"""
function pack_perturbations!(solver::DiffusionSolver, Q::Array{Float64,3},
                             Q_bar::Array{Float64,3}, grid::Grid2D, config::SimConfig)
    nx, nz, ng = grid.nx, grid.nz, grid.ng
    u_pert_field, w_field, T_pert_field = solver.u_pert, solver.w_pert, solver.T_pert
    rho_field, u_bar_field, T_bar_field = solver.rho_f, solver.u_bar, solver.T_bar
    for A in (u_pert_field, w_field, T_pert_field, rho_field, u_bar_field, T_bar_field)
        size(A) == (nx, nz) || error("DiffusionSolver buffer $(size(A)), grid is $((nx, nz))")
    end
    size(Q) == size(Q_bar) == (4, nx + 2ng, nz + 2ng) ||
        error("diffusion_step!: Q $(size(Q)), Q_bar $(size(Q_bar)), grid $((nx, nz, ng))")

    @inbounds for j in 1:nz, i in 1:nx
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
    return nothing
end

"""
    unpack_conservative!(Q, Q_bar, solver, u_pert_new, w_new, T_pert_new, grid, config)

Write the physical cells of Q [4, nx + 2ng, nz + 2ng] from the diffused perturbations
u_pert_new, w_new (m/s) and T_pert_new (K), each [nx*nz] with linear index
(j - 1) nx + i, and the solver buffers rho_f, u_bar, T_bar [nx, nz]: rho unchanged,
u = u_bar + u', T = T_bar + T', P = rho R T floored at (gamma - 1) Q_bar[4] 1e-4,
E = P/(gamma - 1) + rho (u^2 + w^2)/2. A separate function so that the three
solution vectors, which solve_cn! returns through untyped solver fields, have a
concrete type inside the loop.
"""
function unpack_conservative!(Q::Array{Float64,3}, Q_bar::Array{Float64,3},
                              solver::DiffusionSolver, u_pert_new::Vector{Float64},
                              w_new::Vector{Float64}, T_pert_new::Vector{Float64},
                              grid::Grid2D, config::SimConfig)
    nx, nz, ng = grid.nx, grid.nz, grid.ng
    rho_field, u_bar_field, T_bar_field = solver.rho_f, solver.u_bar, solver.T_bar
    for A in (rho_field, u_bar_field, T_bar_field)
        size(A) == (nx, nz) || error("DiffusionSolver buffer $(size(A)), grid is $((nx, nz))")
    end
    length(u_pert_new) == length(w_new) == length(T_pert_new) == nx * nz ||
        error("unpack_conservative!: solution lengths differ from nx*nz = $(nx * nz)")
    size(Q) == size(Q_bar) == (4, nx + 2ng, nz + 2ng) ||
        error("diffusion_step!: Q $(size(Q)), Q_bar $(size(Q_bar)), grid $((nx, nz, ng))")

    @inbounds for j in 1:nz, i in 1:nx
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
