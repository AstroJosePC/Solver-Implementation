# =============================================================================
# Wave forcing functions for the bottom boundary
# Srivastava et al. (2022), Equations 29, 31, 33, 43, 49
# =============================================================================

"""
    WaveParams

Parameters for wave forcing at the bottom boundary.
"""
@kwdef struct WaveParams
    amplitude::Float64 = 1.0    # w amplitude [m/s]
    omega::Float64 = 0.01       # angular frequency [rad/s]
    kx::Float64 = 2e-4          # horizontal wavenumber [1/m]
    sigma_t::Float64 = 600.0    # temporal Gaussian width [s]
    t0::Float64 = 1200.0        # Gaussian center time [s]
    x0::Float64 = 0.0           # spatial Gaussian center [m]
    sigma_x::Float64 = 10000.0  # spatial Gaussian width [m]
end

"""
    MountainParams

Parameters for orographic (mountain) forcing (Gaussian hill).
"""
@kwdef struct MountainParams
    h_max::Float64 = 1000.0   # mountain height [m]
    sigma::Float64 = 10000.0  # mountain half-width [m]
    x0::Float64 = 100000.0    # mountain center [m]
    U_bar::Float64 = 20.0     # background wind speed [m/s]
end

"""
    MountainBellParams

Parameters for orographic (mountain) forcing (Bell-shaped hill, Queney 1948).
h(x) = h_max / (1 + ((x - x0) / a)^2)
Used in Srivastava et al. (2022) Section 3.3.
"""
@kwdef struct MountainBellParams
    h_max::Float64 = 10.0     # peak height [m]
    a::Float64 = 20000.0      # half-width scale [m]
    x0::Float64 = 120000.0    # center position [m]
    U_bar::Float64 = 18.0     # background wind speed [m/s]
end

"""
    PulseParams

Parameters for pressure pulse initial condition.
"""
@kwdef struct PulseParams
    P_tilde::Float64 = 100.0   # pressure perturbation amplitude [Pa]
    sigma_x::Float64 = 10000.0 # horizontal width [m]
    sigma_z::Float64 = 5000.0  # vertical width [m]
    x0::Float64 = 100000.0     # center x [m]
    z0::Float64 = 20000.0      # center z [m]
end

"""
    wave_forcing_linear_ramp(x, t, params)

Continuous harmonic plane wave forcing with a smooth temporal ramp-up.
w(x, t) = A * cos(omega*t + kx*x) * tanh(t / t0)

Used to prevent shock transients in steady-state envelope validations.
"""
function wave_forcing_linear_ramp(x::Float64, t::Float64, p::WaveParams)
    return p.amplitude * cos(p.omega * t + p.kx * x) * tanh(t / p.t0)
end

"""
    wave_forcing_linear(x, t, params)

Monochromatic plane wave forcing (Eq. 31).
w(x, t) = A * cos(omega*t + kx*x)

Used for linear dispersion validation.
"""
function wave_forcing_linear(x::Float64, t::Float64, p::WaveParams)
    return p.amplitude * cos(p.omega * t + p.kx * x)
end

"""
    wave_forcing_gaussian(x, t, params)

Gaussian-enveloped wave packet (Eq. 29).
w(x, t) = A * cos(omega*(t-t0) + kx*x) * exp(-(t-t0)^2 / (2*sigma_t^2))

Used for sponge layer testing.
"""
function wave_forcing_gaussian(x::Float64, t::Float64, p::WaveParams)
    envelope = exp(-(t - p.t0)^2 / (2.0 * p.sigma_t^2))
    return p.amplitude * cos(p.omega * (t - p.t0) + p.kx * x) * envelope
end

"""
    wave_forcing_critical(x, t, params)

Gaussian wave packet with spatial localization (Eq. 33).
Adds a Gaussian envelope in x for critical layer testing.
"""
function wave_forcing_critical(x::Float64, t::Float64, p::WaveParams)
    t_env = exp(-(t - p.t0)^2 / (2.0 * p.sigma_t^2))
    x_env = exp(-(x - p.x0)^2 / (2.0 * p.sigma_x^2))
    return p.amplitude * cos(p.omega * (t - p.t0) + p.kx * x) * t_env * x_env
end

"""
    mountain_forcing(x, params)

Orographic forcing (Gaussian hill).
w(x) = U_bar * dh/dx where h(x) = h_max * exp(-(x-x0)^2 / (2*sigma^2))
"""
function mountain_forcing(x::Float64, p::MountainParams)
    dhdx = -p.h_max * (x - p.x0) / p.sigma^2 * exp(-(x - p.x0)^2 / (2.0 * p.sigma^2))
    return p.U_bar * dhdx
end

"""
    mountain_bell_forcing(x, params)

Orographic forcing (Bell-shaped hill, Eq. 41/49 style with Queney profile).
h(x) = b / (1 + (x/a)^2)
dh/dx = -2*b*(x/a^2) / (1 + (x/a)^2)^2
w(x) = U_bar * dh/dx
"""
function mountain_bell_forcing(x::Float64, p::MountainBellParams)
    X = (x - p.x0)
    denom = 1.0 + (X / p.a)^2
    dhdx = -2.0 * p.h_max * (X / p.a^2) / (denom^2)
    return p.U_bar * dhdx
end

"""
    pressure_pulse_ic!(Q, grid, config, params)

Apply Gaussian pressure pulse as initial condition (Eq. 43).
P(x,z,0) = P_bg + P_tilde * exp(-(z-z0)^2/(2*sigma_z^2)) * exp(-(x-x0)^2/(2*sigma_x^2))

Modifies only the energy component of Q.
"""
function pressure_pulse_ic!(Q::Array{Float64,3}, grid::Grid2D,
                            config::SimConfig, p::PulseParams)
    ng = grid.ng
    for j in phys_z(grid)
        z = grid.z_full[j]
        for i in phys_x(grid)
            x = grid.x_full[i]

            dP = p.P_tilde * exp(-(z - p.z0)^2 / (2.0 * p.sigma_z^2)) *
                              exp(-(x - p.x0)^2 / (2.0 * p.sigma_x^2))

            # Add pressure perturbation to energy
            Q[4, i, j] += dP / (config.gamma - 1.0)
        end
    end
    return nothing
end

"""
    apply_bottom_forcing!(Q, grid, config, t, wave_params, forcing_type)

Apply wave forcing at the bottom boundary by setting w in the bottom ghost cells.

forcing_type: :linear, :gaussian, :critical, :mountain, :mountain_bell, :none
"""
function apply_bottom_forcing!(Q::Array{Float64,3}, grid::Grid2D,
                               config::SimConfig, t::Float64,
                               forcing_params,
                               forcing_type::Symbol)
    if forcing_type == :none
        return nothing
    end

    ng = grid.ng

    for i in phys_x(grid)
        x = grid.x_full[i]

        # Compute forced vertical velocity
        w_forced = if forcing_type == :linear
            p = forcing_params isa WaveParams ? forcing_params : WaveParams()
            wave_forcing_linear(x, t, p)
        elseif forcing_type == :linear_ramp
            p = forcing_params isa WaveParams ? forcing_params : WaveParams()
            wave_forcing_linear_ramp(x, t, p)
        elseif forcing_type == :gaussian
            p = forcing_params isa WaveParams ? forcing_params : WaveParams()
            wave_forcing_gaussian(x, t, p)
        elseif forcing_type == :critical
            p = forcing_params isa WaveParams ? forcing_params : WaveParams()
            wave_forcing_critical(x, t, p)
        elseif forcing_type == :mountain
            p = forcing_params isa MountainParams ? forcing_params : MountainParams()
            mountain_forcing(x, p)
        elseif forcing_type == :mountain_bell
            p = forcing_params isa MountainBellParams ? forcing_params : MountainBellParams()
            mountain_bell_forcing(x, p)
        else
            0.0
        end

        # Set w in ghost cells so the face at z=0 carries w_forced.
        # We use 1st-order enforcement: (w_ghost + w_phys)/2 = w_forced
        # => w_ghost = 2*w_forced - w_phys
        j_phys = ng + 1
        rho_p = Q[1, i, j_phys]
        w_p   = Q[3, i, j_phys] / rho_p

        w_g = 2.0 * w_forced - w_p

        for g in 1:ng
            j_ghost = ng + 1 - g   # j=ng (g=1), j=ng-1 (g=2)
            rho   = Q[1, i, j_ghost]
            u_g   = Q[2, i, j_ghost] / rho
            w_old = Q[3, i, j_ghost] / rho
            P     = (config.gamma - 1.0) * (Q[4, i, j_ghost] - 0.5 * rho * (u_g^2 + w_old^2))
            P     = max(P, 0.0)
            Q[3, i, j_ghost] = rho * w_g
            Q[4, i, j_ghost] = P / (config.gamma - 1.0) + 0.5 * rho * (u_g^2 + w_g^2)
        end
    end

    return nothing
end
