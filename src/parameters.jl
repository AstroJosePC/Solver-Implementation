# =============================================================================
# Physical constants and simulation configuration
# Based on Srivastava et al. (2022), Appendix A
# =============================================================================

# --- Physical Constants ---
const k_B = 1.380649e-23      # Boltzmann constant [J/K]
const N_A = 6.02214076e23     # Avogadro number [1/mol]
const R_universal = 8.314462  # Universal gas constant [J/(mol*K)]

# --- Molecular masses [kg/mol] ---
const M_CO2 = 44.01e-3
const M_N2  = 28.01e-3
const M_O2  = 32.00e-3
const M_Ar  = 39.95e-3
const M_O   = 16.00e-3

# --- Planet enumeration ---
@enum Planet Earth Mars Venus

"""
    SimConfig

Complete simulation configuration struct.

Fields:
- `nx, nz`: number of physical grid cells in x, z
- `dx, dz`: grid spacing [m]
- `ng`: number of ghost cells per side (default 2)
- `CFL`: CFL number for adaptive timestep (default 0.8)
- `gamma`: ratio of specific heats
- `g`: gravitational acceleration [m/s^2]
- `R_gas`: specific gas constant [J/(kg*K)]
- `Pr`: Prandtl number
- `T0`: reference/surface temperature [K]
- `P_surface`: surface pressure [Pa]
- `t_end`: simulation end time [s]
- `output_interval`: steps between NetCDF snapshots
- `planet`: Planet enum (Earth, Mars, Venus)
- `enable_viscosity`: enable viscous diffusion
- `enable_cooling`: enable Newtonian cooling
- `tau_r_min`: peak Newtonian cooling timescale [s] (default 2500.0)
- `z_cooling_peak`: altitude of maximum cooling [m] (default 70000.0)
- `cooling_width`: vertical spread of cooling region [m] (default 25000.0)
- `sponge_frac`: fraction of domain height for sponge layer
- `K0_sponge`: maximum sponge damping coefficient [1/s]
"""
@kwdef struct SimConfig
    # Grid
    nx::Int = 200
    nz::Int = 80
    dx::Float64 = 1000.0     # [m]
    dz::Float64 = 500.0      # [m]
    ng::Int = 2

    # Numerics
    CFL::Float64 = 0.8

    # Gas properties
    gamma::Float64 = 1.4
    g::Float64 = 9.81
    R_gas::Float64 = 287.0   # J/(kg*K), specific gas constant
    Pr::Float64 = 0.7

    # Thermodynamic state
    T0::Float64 = 300.0      # Surface/reference temperature [K]
    P_surface::Float64 = 1.01325e5  # [Pa]

    # Simulation control
    t_end::Float64 = 3600.0  # [s]
    output_interval::Int = 100

    # Planet
    planet::Planet = Earth

    # Physics toggles
    enable_viscosity::Bool = true
    enable_cooling::Bool = false
    
    # Altitude-dependent CO2 radiative damping parameters
    tau_r_min::Float64 = 2500.0       # [s] Peak damping timescale
    z_cooling_peak::Float64 = 70000.0 # [m] Altitude of maximum cooling
    cooling_width::Float64 = 25000.0  # [m] Spread of the cooling region

    # Sponge layer
    sponge_frac::Float64 = 0.15   # top 15% of domain
    K0_sponge::Float64 = 0.005    # [1/s], Srivastava optimal
end

"""
    kinematic_viscosity(T, rho, planet)

Compute kinematic viscosity nu [m^2/s] using planet-specific power law.
Srivastava Eq. A3: mu = C * T^0.69, then nu = mu / rho.
"""
function kinematic_viscosity(T::Float64, rho::Float64, planet::Planet)
    if planet == Mars
        mu = 4.2e-7 * T^0.69   # CO2-dominated
    elseif planet == Venus
        mu = 3.38e-7 * T^0.69
    else  # Earth
        mu = 3.56e-7 * T^0.69
    end
    nu = mu / rho
    # Cap kinematic viscosity to prevent implicit matrix ill-conditioning
    # at extreme altitudes (e.g., above 130 km where continuum breaks down anyway)
    return min(nu, 1e4)
end

"""
    thermal_diffusivity(nu, Pr)

Thermal diffusivity alpha = nu / Pr.
"""
thermal_diffusivity(nu::Float64, Pr::Float64) = nu / Pr

"""
    earth_config(; kwargs...)

Convenience constructor for Earth isothermal atmosphere.
"""
function earth_config(; kwargs...)
    SimConfig(;
        gamma = 1.4,
        g = 9.81,
        R_gas = 287.0,
        T0 = 300.0,
        P_surface = 1.01325e5,
        planet = Earth,
        kwargs...
    )
end

"""
    mars_config(; kwargs...)

Convenience constructor for Mars CO2 atmosphere.
"""
function mars_config(; kwargs...)
    SimConfig(;
        gamma = 1.31,
        g = 3.72,
        R_gas = R_universal / M_CO2,  # ~188.9 J/(kg*K)
        T0 = 210.0,
        P_surface = 636.0,  # ~636 Pa mean surface pressure
        planet = Mars,
        Pr = 0.73,
        kwargs...
    )
end
