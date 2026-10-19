# =============================================================================
# Acoustic Wave Propagation on Mars Benchmark (Srivastava et al. 2022, Sec 3.4)
# Matches Tables 12 and 13 exactly to reproduce Figures 12, 13, 14, 15.
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using MAGNUSP
using Printf

# Table 12: Domain Parameters
# L_x = 200 km, L_z = 40 km
# To match the "pulse centered at x=0" and monitoring stations at x=20, 40, 60, 80 km
# We define physical domain [0, 200] km and put pulse at x0 = 100 km.
# The monitoring stations will then be at 120, 140, 160, 180 km.
const LX = 200e3
const LZ = 40e3
const DX = 1000.0
const DZ = 500.0

nx = Int(round(LX / DX))
nz = Int(round(LZ / DZ))

# Load interpolated MCS profile
atm_file = "external/digitized_figure_11/mars_mcs_uniform.csv"

# Table 13: Time Parameters
# t_max = 8000 s, outputs every 5 s
config = mars_config(
    nx = nx,
    nz = nz,
    dx = DX,
    dz = DZ,
    CFL = 0.4,
    t_end = 8000.0,
    output_interval = 5,      # Outputs every 5s for Nyquist 100 mHz
    enable_viscosity = true,  # Text says classical absorption is modeled
    enable_cooling = false,
    sponge_frac = 0.0,         # No sponge layer for this simulation (paper spec)
    K0_sponge = 0.05
)

# Initial Pressure Pulse Parameters (Eq 43)
# Centered at z0=20km. x0=0km in paper, but we use x0=100km to stay in [0, 200] domain.
# P_tilde is half the ambient pressure at 20km; they report 44.83 Pa
pulse_params = MAGNUSP.PulseParams(
    P_tilde = 44.83,          # ~Half of ambient pressure at 20km
    sigma_x = 500.0,          # 500 m
    sigma_z = 1000.0,         # 1 km
    x0 = 100e3,               # Center of domain
    z0 = 20e3                 # 20 km altitude
)

outdir = "output/mars_acoustic_test"

println("=" ^ 88)
println("MARS ACOUSTIC WAVE BENCHMARK (Table 12/13)")
println("=" ^ 88)
@printf("Domain: %.1f x %.1f km, dx=%.2f km, dz=%.2f km\n", LX/1e3, LZ/1e3, DX/1e3, DZ/1e3)
@printf("Run: t_end=%.1f s, CFL=%.2f, dt_out=%.1f s\n", config.t_end, config.CFL, config.output_interval)
println("Viscosity: ON, Sponge: OFF (as per paper)")

state = run_simulation(config;
    forcing_type = :none,            # No steady bottom forcing, just initial pulse
    pulse_params = pulse_params,
    bc_x_type = :outflow,            # Side boundaries set to outflow (paper specified)
    bc_bottom_type = :outflow,       # Bottom boundary set to outflow
    atm_file = atm_file,             # Digitized MCS profile
    outdir = outdir,
    diag_interval = 200,
)

println()
println("Simulation complete.")
@printf("Final time: %.1f s, steps: %d\n", state.t, state.step)
println("Snapshots: $(outdir)/")
