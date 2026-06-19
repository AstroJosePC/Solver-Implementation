# =============================================================================
# Background atmosphere profiles
# Hydrostatic equilibrium initialization for isothermal and tabulated cases
# =============================================================================

"""
    AtmosphereProfile

Background state arrays evaluated on the z-grid.
"""
struct AtmosphereProfile
    z::Vector{Float64}      # altitude [m]
    rho0::Vector{Float64}   # density [kg/m^3]
    P0::Vector{Float64}     # pressure [Pa]
    T0::Vector{Float64}     # temperature [K]
end

"""
    scale_height(T, g, R_gas)

Pressure scale height H = R_gas * T / g [m].
"""
scale_height(T::Float64, g::Float64, R_gas::Float64) = R_gas * T / g

"""
    sound_speed(gamma, R_gas, T)

Adiabatic sound speed c_s = sqrt(gamma * R_gas * T) [m/s].
"""
sound_speed(gamma::Float64, R_gas::Float64, T::Float64) = sqrt(gamma * R_gas * T)

"""
    brunt_vaisala(gamma, g, cs)

Brunt-Vaisala (buoyancy) frequency N = sqrt(gamma-1) * g / cs [rad/s].
For isothermal atmosphere.
"""
brunt_vaisala(gamma::Float64, g::Float64, cs::Float64) = sqrt(gamma - 1.0) * g / cs

"""
    acoustic_cutoff(gamma, g, cs)

Acoustic cutoff frequency omega_ac = gamma * g / (2 * cs) [rad/s].
"""
acoustic_cutoff(gamma::Float64, g::Float64, cs::Float64) = gamma * g / (2.0 * cs)

"""
    isothermal_atmosphere(config, z_grid)

Initialize a hydrostatic isothermal atmosphere profile.

For isothermal gas:
  P(z) = P_surface * exp(-z / H)
  rho(z) = P(z) / (R_gas * T0)
"""
function isothermal_atmosphere(config::SimConfig, z_grid::AbstractVector{Float64})
    H = scale_height(config.T0, config.g, config.R_gas)
    nz = length(z_grid)

    P0   = Vector{Float64}(undef, nz)
    rho0 = Vector{Float64}(undef, nz)
    T0   = fill(config.T0, nz)

    for j in 1:nz
        P0[j]   = config.P_surface * exp(-z_grid[j] / H)
        rho0[j] = P0[j] / (config.R_gas * config.T0)
    end

    return AtmosphereProfile(collect(z_grid), rho0, P0, T0)
end

"""
    brunt_vaisala_local(T_below, T_above, dz, g, gamma, R_gas)

Brunt-Väisälä frequency from a finite-difference temperature gradient.
N² = (g/T)(dT/dz + g/cp) where cp = gamma*R_gas/(gamma-1).

For non-isothermal atmospheres (e.g., MCS profile), this gives the local N
at the interface between two cells, rather than the isothermal approximation.
Returns zero if the layer is convectively unstable (N² < 0).
"""
function brunt_vaisala_local(T_below::Float64, T_above::Float64,
                             dz::Float64, g::Float64,
                             gamma::Float64, R_gas::Float64)
    T_avg = 0.5 * (T_below + T_above)
    dTdz = (T_above - T_below) / dz
    cp = gamma * R_gas / (gamma - 1.0)
    N2 = (g / T_avg) * (dTdz + g / cp)
    N2 = max(N2, 0.0)  # convective stability guard
    return sqrt(N2)
end

using DelimitedFiles

"""
    tabulated_atmosphere(filepath::String, z_target::AbstractVector{Float64})

Loads a background atmosphere from a CSV file and interpolates it onto the
solver's `z_target` grid (which includes ghost cells).

The CSV must have a header and the following columns:
z_m, rho, P, T, nu
"""
function tabulated_atmosphere(filepath::String, z_target::AbstractVector{Float64})
    data, header = readdlm(filepath, ',', Float64, header=true)
    
    # Extract columns
    z_data = data[:, 1]
    rho_data = data[:, 2]
    P_data = data[:, 3]
    T_data = data[:, 4]

    # Linear interpolation function
    function interpolate(z_val, z_arr, val_arr)
        if z_val <= z_arr[1]
            return val_arr[1]
        elseif z_val >= z_arr[end]
            return val_arr[end]
        else
            idx = searchsortedlast(z_arr, z_val)
            t = (z_val - z_arr[idx]) / (z_arr[idx+1] - z_arr[idx])
            return val_arr[idx] * (1.0 - t) + val_arr[idx+1] * t
        end
    end

    nz = length(z_target)
    rho0 = zeros(nz)
    P0 = zeros(nz)
    T0 = zeros(nz)

    for j in 1:nz
        z = z_target[j]
        rho0[j] = interpolate(z, z_data, rho_data)
        P0[j]   = interpolate(z, z_data, P_data)
        T0[j]   = interpolate(z, z_data, T_data)
    end

    return AtmosphereProfile(collect(z_target), rho0, P0, T0)
end
