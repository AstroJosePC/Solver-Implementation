# =============================================================================
# Earth Isothermal Linear Dispersion Test
# Validates against Srivastava et al. (2022) Fig. 6-7
#
# Parameters from Tables 5 and 6:
#   Domain: 40 km horizontal x 160 km vertical (20 km sponge at top)
#   dx = dz = 0.5 km, T0 = 300 K (Earth isothermal)
#   Gaussian wave packet: w_tilde=0.001 m/s, t0=1200s, sigma_t=600s
#   omega=0.007 rad/s, lambda_x=40 km
#   t_end = 8000 s, CFL = 0.8
#   Expected: lambda_z = 16.7 km (paper), N = 0.0182 rad/s
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using MAGNUSP

# --- Earth isothermal configuration (Table 5) ---
# Domain: Lx=40km (x from -20 to +20 km), Lz=160km (0 to 160 km)
# 20 km sponge at top: z_sponge=140 km, z_max=160 km
config = earth_config(
    nx = 80,            # 40 km / 0.5 km
    nz = 320,           # 160 km / 0.5 km
    dx = 500.0,         # 0.5 km [m]
    dz = 500.0,         # 0.5 km [m]
    CFL = 0.8,
    T0 = 239.0,         # K — paper Section 3 opening: "Earth isothermal T=239K" for ALL Earth cases
    t_end = 8000.0,     # Table 6
    output_interval = 500,
    enable_viscosity = false,  # inviscid for dispersion test (Section 3.1)
    enable_cooling = false,
    sponge_frac = 20.0 / 160.0,   # 20 km sponge out of 160 km domain
    K0_sponge = 0.005,
)

# --- Wave parameters (Table 6) ---
cs = MAGNUSP.sound_speed(config.gamma, config.R_gas, config.T0)
H  = MAGNUSP.scale_height(config.T0, config.g, config.R_gas)
N  = MAGNUSP.brunt_vaisala(config.gamma, config.g, cs)
omega_ac = MAGNUSP.acoustic_cutoff(config.gamma, config.g, cs)

println("Earth isothermal atmosphere:")
println("  T0 = $(config.T0) K")
println("  cs = $(round(cs, digits=1)) m/s")
println("  H  = $(round(H/1e3, digits=2)) km")
println("  N  = $(round(N, digits=5)) rad/s (period = $(round(2*pi/N, digits=0)) s)")
println("  omega_ac = $(round(omega_ac, digits=5)) rad/s")
println()

lambda_x = 40e3          # horizontal wavelength [m] (Table 6)
kx = 2.0 * pi / lambda_x
omega = 0.007            # rad/s (Table 6, gravity wave branch)

# Compute expected vertical wavelength from dispersion relation (Eq. 32)
m_squared = kx^2 * (N^2 / omega^2 - 1.0)
if m_squared > 0
    m = sqrt(m_squared)
    lambda_z = 2.0 * pi / m
    println("Gravity wave parameters (Table 6):")
    println("  lambda_x = $(lambda_x/1e3) km, kx = $(round(kx*1e3, digits=4)) 1/km")
    println("  omega = $omega rad/s (period = $(round(2*pi/omega, digits=0)) s)")
    println("  m = $(round(m*1e3, digits=4)) 1/km, lambda_z = $(round(lambda_z/1e3, digits=2)) km")
    println("  Paper reports lambda_z = 16.7 km (Section 3.1)")
    println("  Vertical phase speed: $(round(omega/m, digits=2)) m/s")
    println("  Horizontal phase speed: $(round(omega/kx, digits=2)) m/s")
else
    println("WARNING: omega = $omega is not on the GW branch (m^2 < 0).")
end

wave_params = MAGNUSP.WaveParams(
    amplitude = 0.001,  # w_tilde = 0.001 m/s (Table 6, linear regime)
    omega = omega,
    kx = kx,
    sigma_t = 600.0,    # Table 6
    t0 = 1200.0,        # Table 6
)

# --- Run simulation ---
println("\nStarting simulation...")
println("="^80)

state = run_simulation(config;
    forcing_type = :gaussian,   # Eq. 31 (Gaussian envelope, same as sponge test)
    wave_params = wave_params,
    outdir = "output/dispersion_test",
    diag_interval = 200,
)

println("\nDone! Output files in output/dispersion_test/")
println("Total steps: $(state.step)")
