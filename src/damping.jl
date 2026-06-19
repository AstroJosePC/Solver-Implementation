# =============================================================================
# Thermal damping: Newtonian cooling
# Srivastava et al. (2022), Equations 25-26
#
# Placeholder for full CO2 radiative damping (Eckermann et al. 2011).
# See docs/TODO_tracker.md for status on obtaining coefficients.
# =============================================================================

"""
    compute_tau_r(z, config)

Compute an altitude-dependent Newtonian cooling timescale (tau_r) [s]
as a proxy for CO2 15 µm radiative damping (e.g., Eckermann et al. 2011; Wang et al. 2024).
"""
function compute_tau_r(z::Float64, config::SimConfig)
    # Background/large value far from the cooling peak
    tau_r_max = 1e6
    
    # Calculate a sharper super-Gaussian profile matching the parameters
    # The timescale dips to tau_r_min around z_cooling_peak.
    # We use a power of 8 to create a sharper wall, making the lower atmosphere
    # completely transparent to the waves.
    decay_arg = -0.5 * ((z - config.z_cooling_peak) / config.cooling_width)^8
    
    # We want tau_r(z) to be tau_r_min at z = z_cooling_peak, and larger elsewhere.
    # Inverse Gaussian mapping on the rate (1/tau_r):
    # 1/tau_r(z) = (1/tau_r_min - 1/tau_r_max) * exp(decay_arg) + 1/tau_r_max
    
    inv_tau_r_min = 1.0 / config.tau_r_min
    inv_tau_r_max = 1.0 / tau_r_max
    
    inv_tau = (inv_tau_r_min - inv_tau_r_max) * exp(decay_arg) + inv_tau_r_max
    return 1.0 / inv_tau
end

"""
    newtonian_cooling_step!(Q, grid, config, atm, dt)

Apply Newtonian cooling via Crank-Nicholson time discretization.

Eq. 25-26:
  T^{n+1} = [(1 - dt/(2*tau_r)) * T^n + T_a * dt/tau_r] / (1 + dt/(2*tau_r))

where T_a is the ambient (background) temperature and tau_r is the local cooling timescale.

This modifies the energy equation while keeping rho and momentum unchanged,
effectively relaxing temperature perturbations toward the background state.
"""
function newtonian_cooling_step!(Q::Array{Float64,3}, grid::Grid2D,
                                 config::SimConfig, atm::AtmosphereProfile,
                                 dt::Float64)
    ng = grid.ng

    for j in phys_z(grid)
        # Get ambient temperature at this altitude
        # atm is evaluated on z_full, so index j maps directly
        T_ambient = atm.T0[j]
        z = grid.z_full[j]
        
        # Compute level-specific tau_r and Crank-Nicholson coefficients
        tau_r = compute_tau_r(z, config)
        coeff_a = 1.0 - dt / (2.0 * tau_r)
        coeff_b = dt / tau_r
        coeff_c = 1.0 / (1.0 + dt / (2.0 * tau_r))

        for i in phys_x(grid)
            rho = Q[1, i, j]
            u   = Q[2, i, j] / rho
            w   = Q[3, i, j] / rho
            KE  = 0.5 * rho * (u^2 + w^2)
            P   = (config.gamma - 1.0) * (Q[4, i, j] - KE)
            T   = P / (rho * config.R_gas)

            # Crank-Nicholson update
            T_new = (coeff_a * T + coeff_b * T_ambient) * coeff_c

            # Reconstruct energy with updated temperature
            P_new = rho * config.R_gas * T_new
            Q[4, i, j] = P_new / (config.gamma - 1.0) + KE
        end
    end

    return nothing
end
