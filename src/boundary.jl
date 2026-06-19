# =============================================================================
# Boundary conditions
# Srivastava et al. (2022), Equations 18-19
# Klemp & Durran (1983) radiation BC for gravity waves
#
# For a stratified atmosphere, BCs must preserve hydrostatic balance.
# We reflect/extrapolate *perturbations* from the background state Q_bar,
# not raw state values, to maintain discrete hydrostatic equilibrium.
# =============================================================================

# ---- Klemp-Durran radiation BC state and functions -------------------------

"""
    KlempDurranState

Persistent state for the Klemp-Durran (1983) radiation boundary condition.

The KD BC relates pressure perturbation to vertical velocity perturbation at
the top boundary via the gravity wave dispersion relation in Fourier space:
    p̂'(kx) = ρ₀ (N / |k_w|) ŵ'(kx)
where k_w is the FD-corrected effective wavenumber (KD Appendix, Eq. A2).

Fields:
- `N_bv`: Brunt-Väisälä frequency at the top boundary [rad/s]
- `multiplier`: precomputed ρ₀ N / |k_w| for each Fourier mode
- `phi_prev`: previous-timestep pressure perturbation (for time smoothing)
- `fft_plan`: precomputed real FFT plan (FFTW)
- `initialized`: false until first call (disables time smoothing on first step)
"""
mutable struct KlempDurranState
    N_bv::Float64
    multiplier::Vector{Float64}
    phi_prev::Vector{Float64}
    fft_plan::FFTW.rFFTWPlan
    initialized::Bool
end

"""
    init_klemp_durran(grid, N_bv, rho_top)

Allocate and precompute the KD radiation BC state.

The FD-corrected wavenumber for a cell-centered periodic grid is
    k_w(m) = (2/Δx) sin(π m / nx)    for mode m = 0, …, nx/2
which is the standard effective wavenumber for centered finite differences
(matches the solver's spatial discretization accuracy).

The multiplier for each mode is ρ₀ N / k_w, where ρ₀ is the background
density at the top boundary. The k=0 mode (domain-mean) is zeroed out
since we only radiate perturbations.
"""
function init_klemp_durran(grid::Grid2D, N_bv::Float64, rho_top::Float64)
    nx = grid.nx
    dx = grid.dx

    n_freq = nx ÷ 2 + 1  # number of positive-frequency rfft bins
    multiplier = zeros(n_freq)
    for m in 1:n_freq
        if m == 1
            # k=0 mode: no mean pressure perturbation
            multiplier[m] = 0.0
        else
            # FD-corrected effective wavenumber for centered differences
            kw = (2.0 / dx) * sin(π * (m - 1) / nx)
            multiplier[m] = rho_top * N_bv / kw
        end
    end

    phi_prev = zeros(nx)
    w_tmp = zeros(nx)
    fft_plan = plan_rfft(w_tmp)

    return KlempDurranState(N_bv, multiplier, phi_prev, fft_plan, false)
end

"""
    apply_radiation_top_klemp_durran!(Q, Q_bar, grid, gamma, kd)

Apply the Klemp-Durran (1983) radiation BC at the top boundary.

Algorithm:
1. Extract w' = w - w̄ at the topmost physical row
2. FFT → multiply by ρ₀N/|k_w| → IFFT to get radiation-condition p'
3. Time-smooth p' with previous step (KD Eq. 40, r=0.2 for stability)
4. Set ghost cell pressure from p_bar + p'; density and velocity extrapolated

The time smoothing prevents high-frequency noise feedback. KD recommend
r = 0.2 (80% weight on previous timestep), which we adopt here.
"""
function apply_radiation_top_klemp_durran!(
    Q::Array{Float64,3}, Q_bar::Array{Float64,3},
    grid::Grid2D, gamma::Float64,
    kd::KlempDurranState)

    ng = grid.ng
    nx = grid.nx
    j_top = grid.nz + ng  # last physical cell in z

    # 1. Extract w perturbation at top physical row
    w_pert = Vector{Float64}(undef, nx)
    for i in 1:nx
        ii = i + ng
        rho = Q[1, ii, j_top]
        w = Q[3, ii, j_top] / rho
        w_bar = Q_bar[3, ii, j_top] / max(Q_bar[1, ii, j_top], 1e-20)
        w_pert[i] = w - w_bar
    end

    # 2. FFT -> multiply by rho*N/|k_w| -> IFFT  (KD Eq. A3)
    w_hat = kd.fft_plan * w_pert
    phi_hat = w_hat .* kd.multiplier
    phi_pert = irfft(phi_hat, nx)

    # 3. Time smoothing (KD Eq. 40, r = 0.2)
    r_smooth = 0.2
    if kd.initialized
        phi_pert .= r_smooth .* phi_pert .+ (1.0 - r_smooth) .* kd.phi_prev
    end
    kd.phi_prev .= phi_pert
    kd.initialized = true

    # 4. Set ghost cells using radiation-condition pressure
    #
    # Deep-atmosphere safety: In a 200 km Mars domain, background pressure at
    # the top is O(1e-8) Pa while wave amplitude grows as sqrt(rho_0/rho_top).
    # The KD perturbation can exceed P_bar where the linear assumption breaks
    # down. We clamp phi_pert to a fraction of P_bar to prevent catastrophic
    # instability while preserving KD physics in the well-resolved region.
    alpha_clamp = 0.3  # max |phi_pert| / P_bar_g

    for g_idx in 1:ng
        j_ghost = j_top + g_idx
        for i in 1:nx
            ii = i + ng

            # Background state at ghost cell
            rho_bar_g = Q_bar[1, ii, j_ghost]
            P_bar_g = (gamma - 1.0) * Q_bar[4, ii, j_ghost]

            # Clamp perturbation to maintain linear regime validity
            phi_max = alpha_clamp * P_bar_g
            phi_clamped = clamp(phi_pert[i], -phi_max, phi_max)

            # Radiation-condition pressure: background + clamped KD perturbation
            P_ghost = P_bar_g + phi_clamped

            # Density: scale perturbation by sqrt(rho_bar ratio) (Eq. 19 convention)
            rho_bar_t = Q_bar[1, ii, j_top]
            rho_t = Q[1, ii, j_top]
            rho_pert_scaled = (rho_t - rho_bar_t) * sqrt(rho_bar_g / max(rho_bar_t, 1e-20))
            rho_g = max(rho_bar_g + rho_pert_scaled, rho_bar_g * 0.1)

            # Velocity: zero-order extrapolation from top physical cell, clamped
            rho_phys = max(Q[1, ii, j_top], rho_bar_t * 0.1)
            cs_top = sqrt(gamma * P_bar_g / max(rho_bar_g, 1e-30))
            u_g = clamp(Q[2, ii, j_top] / rho_phys, -cs_top, cs_top)
            w_g = clamp(Q[3, ii, j_top] / rho_phys, -cs_top, cs_top)

            # Reconstruct conservative variables
            Q[1, ii, j_ghost] = rho_g
            Q[2, ii, j_ghost] = rho_g * u_g
            Q[3, ii, j_ghost] = rho_g * w_g
            Q[4, ii, j_ghost] = P_ghost / (gamma - 1.0) + 0.5 * rho_g * (u_g^2 + w_g^2)
        end
    end
    return nothing
end

"""
    apply_periodic_x!(Q, grid)

Apply periodic boundary conditions in x-direction.
Eq. 18: Copy physical edge cells into ghost cells on opposite side.
Uses 2 ghost cells (ng=2).
"""
function apply_periodic_x!(Q::Array{Float64,3}, grid::Grid2D)
    ng = grid.ng
    nx = grid.nx

    for j in 1:grid.nz_total
        for k in 1:4
            # Left ghosts <- right physical edge
            Q[k, 1, j]  = Q[k, nx + 1, j]
            Q[k, 2, j]  = Q[k, nx + 2, j]
            # Right ghosts <- left physical edge
            Q[k, nx + ng + 1, j] = Q[k, ng + 1, j]
            Q[k, nx + ng + 2, j] = Q[k, ng + 2, j]
        end
    end
    return nothing
end

"""
    apply_outflow_x!(Q, Q_bar, grid)

Apply mixed outflow boundary conditions in x-direction.
Left boundary is set to background (incoming wind).
Right boundary uses zero-order extrapolation (outgoing wind).
"""
function apply_outflow_x!(Q::Array{Float64,3}, Q_bar::Array{Float64,3}, grid::Grid2D)
    ng = grid.ng
    nx = grid.nx

    for j in 1:grid.nz_total
        for k in 1:4
            # Left boundary: strict Dirichlet (background state) because u > 0
            Q[k, 1, j] = Q_bar[k, 1, j]
            Q[k, 2, j] = Q_bar[k, 2, j]
            # Right boundary: zero-order extrapolation
            Q[k, nx + ng + 1, j] = Q[k, nx + ng, j]
            Q[k, nx + ng + 2, j] = Q[k, nx + ng, j]
        end
    end
    return nothing
end

"""
    apply_bottom_bc!(Q, Q_bar, grid)

Apply reflective (closed wall) boundary condition at the bottom.

Perturbations from background are reflected: density and pressure perturbations
are symmetric, vertical momentum perturbation is antisymmetric (w -> -w).
Background values are preserved at ghost cell positions to maintain hydrostatic balance.
"""
function apply_bottom_bc!(Q::Array{Float64,3}, Q_bar::Array{Float64,3}, grid::Grid2D)
    ng = grid.ng
    for i in 1:grid.nx_total
        for g_idx in 1:ng
            j_ghost  = ng + 1 - g_idx   # ghost cell index
            j_mirror = ng + g_idx        # corresponding interior cell

            # Reflect perturbations from background state
            for k in 1:4
                perturbation = Q[k, i, j_mirror] - Q_bar[k, i, j_mirror]
                Q[k, i, j_ghost] = Q_bar[k, i, j_ghost] + perturbation
            end

            # Antisymmetric reflection for vertical momentum (w -> -w)
            w_pert = Q[3, i, j_mirror] - Q_bar[3, i, j_mirror]
            Q[3, i, j_ghost] = Q_bar[3, i, j_ghost] - w_pert
        end
    end
    return nothing
end

"""
    apply_outflow_top!(Q, Q_bar, grid)

Apply outflow boundary condition at the top of the domain.
Eq. 19: perturbation extrapolated with density scaling.

Q_ghost = Q_bar_ghost + (Q_phys - Q_bar_phys) * sqrt(rho_bar_ghost / rho_bar_phys)

This allows waves to exit without significant reflection while preserving
the hydrostatic background in the ghost cells.

Note: a Sommerfeld-like linear gradient extrapolation was tested but caused
acoustic instability at T=239K because the near-vacuum top cells (rho~1e-10)
amplify floating-point noise into acoustic waves that traverse the domain in
~400 steps. Eq. 19 zero-order is numerically stable.
"""
function apply_outflow_top!(Q::Array{Float64,3}, Q_bar::Array{Float64,3}, grid::Grid2D)
    ng = grid.ng
    nz = grid.nz
    j_top = nz + ng  # last physical cell

    for i in 1:grid.nx_total
        for g_idx in 1:ng
            j_ghost = j_top + g_idx

            # Density ratio for amplitude scaling
            rho_bar_ghost = Q_bar[1, i, j_ghost]
            rho_bar_phys  = Q_bar[1, i, j_top]

            if rho_bar_phys > 1e-30
                scale = sqrt(rho_bar_ghost / rho_bar_phys)
            else
                scale = 1.0
            end

            for k in 1:4
                perturbation = Q[k, i, j_top] - Q_bar[k, i, j_top]
                Q[k, i, j_ghost] = Q_bar[k, i, j_ghost] + perturbation * scale
            end
        end
    end
    return nothing
end

"""
    apply_outflow_top_characteristic!(Q, Q_bar, grid, gamma, R_gas, T0)

Apply characteristic-based non-reflecting outflow BC at the top boundary.
Based on Thompson (1987) / Poinsot & Lele (1992) NSCBC formulation.

For the 2D Euler equations in the z-direction, the characteristic wave speeds are:
- L1: w - cs  (incoming acoustic wave for w < cs)
- L2: w       (entropy wave, outgoing for w >= 0)
- L3: w       (vorticity/shear wave, outgoing for w >= 0)
- L4: w + cs  (outgoing acoustic wave)

For gravity waves where w << cs, the L1 wave is incoming. Setting its amplitude
to zero prevents spurious reflections from the top boundary.

Uses BACKGROUND state for sound speed and characteristic decomposition to ensure
numerical stability at near-vacuum densities (rho ~ 1e-10 to 1e-8 kg/m^3).
"""
function apply_outflow_top_characteristic!(Q::Array{Float64,3}, Q_bar::Array{Float64,3},
                                           grid::Grid2D, gamma::Float64,
                                           R_gas::Float64, T0::Float64)
    ng = grid.ng
    nz = grid.nz
    dz = grid.dz
    j_top  = nz + ng      # last physical cell
    j_top1 = j_top - 1    # one cell below top

    # Background sound speed (constant for isothermal atmosphere)
    cs = sqrt(gamma * R_gas * T0)

    for i in 1:grid.nx_total
        # --- Extract PERTURBATION primitives at two topmost physical cells ---
        # The background hydrostatic gradient must be removed before characteristic
        # decomposition, otherwise the huge exp(-z/H) pressure gradient gets
        # amplified by cs and creates unphysical ghost cell values.
        rho_bar_t  = Q_bar[1, i, j_top]
        rho_bar_1  = Q_bar[1, i, j_top1]

        rho_t = Q[1, i, j_top]
        u_t   = Q[2, i, j_top] / rho_t
        w_t   = Q[3, i, j_top] / rho_t
        E_t   = Q[4, i, j_top]
        P_t   = (gamma - 1.0) * (E_t - 0.5 * rho_t * (u_t^2 + w_t^2))

        rho_1 = Q[1, i, j_top1]
        u_1   = Q[2, i, j_top1] / rho_1
        w_1   = Q[3, i, j_top1] / rho_1
        E_1   = Q[4, i, j_top1]
        P_1   = (gamma - 1.0) * (E_1 - 0.5 * rho_1 * (u_1^2 + w_1^2))

        # Background primitives (u_bar = w_bar = 0)
        P_bar_t = (gamma - 1.0) * (Q_bar[4, i, j_top] - 0.0)   # no KE in background
        P_bar_1 = (gamma - 1.0) * (Q_bar[4, i, j_top1] - 0.0)

        # Perturbation gradients (background gradient subtracted)
        drho_dz = ((rho_t - rho_bar_t) - (rho_1 - rho_bar_1)) / dz
        du_dz   = (u_t - u_1) / dz       # u_bar = 0 everywhere
        dw_dz   = (w_t - w_1) / dz       # w_bar = 0 everywhere
        dP_dz   = ((P_t - P_bar_t) - (P_1 - P_bar_1)) / dz

        # Use background density for stable characteristic scaling
        rho_ref = rho_bar_t

        # --- Characteristic wave amplitudes (LODI) on perturbation gradients ---
        # L1 = (w - cs)(dP'/dz - rho_ref*cs*dw'/dz)  incoming acoustic
        # L2 = w*(cs^2*drho'/dz - dP'/dz)             entropy
        # L3 = w*du'/dz                                vorticity
        # L4 = (w + cs)(dP'/dz + rho_ref*cs*dw'/dz)   outgoing acoustic
        L1 = (w_t - cs) * (dP_dz - rho_ref * cs * dw_dz)
        L2 = w_t * (cs^2 * drho_dz - dP_dz)
        L3 = w_t * du_dz
        L4 = (w_t + cs) * (dP_dz + rho_ref * cs * dw_dz)

        # --- Zero out incoming wave amplitudes ---
        if (w_t - cs) < 0.0
            L1 = 0.0
        end
        if w_t < 0.0
            L2 = 0.0
            L3 = 0.0
        end
        if (w_t + cs) < 0.0
            L4 = 0.0
        end

        # --- Filtered perturbation spatial gradients from outgoing-only L's ---
        # For w ~ 0 (gravity waves), L2=L3=0 always. Only L4 survives.
        # Perturbation gradients from acoustic-only:
        dP_dz_f  = 0.5 * (L4 + L1)
        dw_dz_f  = (1.0 / (2.0 * rho_ref * cs)) * (L4 - L1)
        drho_dz_f = dP_dz_f / cs^2    # acoustic relation: drho' = dP'/cs^2
        du_dz_f  = 0.0                # no shear extrapolation for w ~ 0

        # --- Fill ghost cells: background + filtered perturbation extrapolation ---
        # Perturbation at top physical cell
        rho_pert_t = rho_t - rho_bar_t
        P_pert_t   = P_t - P_bar_t

        for g_idx in 1:ng
            j_ghost = j_top + g_idx
            dist = Float64(g_idx) * dz

            # Background at ghost cell (from Q_bar, preserves hydrostatic balance)
            rho_bar_g = Q_bar[1, i, j_ghost]
            P_bar_g   = (gamma - 1.0) * Q_bar[4, i, j_ghost]

            # Extrapolate perturbations using filtered gradients
            rho_pert_g = rho_pert_t + drho_dz_f * dist
            u_g        = u_t + du_dz_f * dist
            w_g        = w_t + dw_dz_f * dist
            P_pert_g   = P_pert_t + dP_dz_f * dist

            # Reconstruct total = background + perturbation
            rho_g = rho_bar_g + rho_pert_g
            P_g   = P_bar_g + P_pert_g

            # Floor for near-vacuum stability
            rho_g = max(rho_g, rho_bar_g * 0.01)
            P_g   = max(P_g, P_bar_g * 0.01)

            # Convert to conservative
            E_g = P_g / (gamma - 1.0) + 0.5 * rho_g * (u_g^2 + w_g^2)

            Q[1, i, j_ghost] = rho_g
            Q[2, i, j_ghost] = rho_g * u_g
            Q[3, i, j_ghost] = rho_g * w_g
            Q[4, i, j_ghost] = E_g
        end
    end
    return nothing
end

"""
    apply_outflow_bottom!(Q, Q_bar, grid)

Apply outflow boundary condition at the bottom of the domain.
Eq. 19 style: perturbation extrapolated with density scaling.
"""
function apply_outflow_bottom!(Q::Array{Float64,3}, Q_bar::Array{Float64,3}, grid::Grid2D)
    ng = grid.ng
    j_bot = ng + 1  # first physical cell

    for i in 1:grid.nx_total
        for g_idx in 1:ng
            j_ghost = ng + 1 - g_idx

            rho_bar_ghost = Q_bar[1, i, j_ghost]
            rho_bar_phys  = Q_bar[1, i, j_bot]

            if rho_bar_phys > 1e-30
                scale = sqrt(rho_bar_ghost / rho_bar_phys)
            else
                scale = 1.0
            end

            for k in 1:4
                perturbation = Q[k, i, j_bot] - Q_bar[k, i, j_bot]
                Q[k, i, j_ghost] = Q_bar[k, i, j_ghost] + perturbation * scale
            end
        end
    end
    return nothing
end

function apply_all_bcs!(Q::Array{Float64,3}, Q_bar::Array{Float64,3},
                        grid::Grid2D, gamma::Float64;
                        bc_type::Symbol = :eq19,
                        bc_x_type::Symbol = :periodic,
                        bc_bottom_type::Symbol = :reflective,
                        R_gas::Float64 = 287.0,
                        T0::Float64 = 239.0,
                        kd_state::Union{Nothing, KlempDurranState} = nothing)
    if bc_x_type == :outflow
        apply_outflow_x!(Q, Q_bar, grid)
    else
        apply_periodic_x!(Q, grid)
    end

    if bc_bottom_type == :outflow
        apply_outflow_bottom!(Q, Q_bar, grid)
    else
        apply_bottom_bc!(Q, Q_bar, grid)
    end

    if bc_type == :klemp_durran && kd_state !== nothing
        apply_radiation_top_klemp_durran!(Q, Q_bar, grid, gamma, kd_state)
    elseif bc_type == :characteristic
        apply_outflow_top_characteristic!(Q, Q_bar, grid, gamma, R_gas, T0)
    else
        apply_outflow_top!(Q, Q_bar, grid)
    end
    return nothing
end
