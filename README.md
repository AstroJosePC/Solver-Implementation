# Model for Acoustic Gravity wave Numerical Simulation in Planetary atmospheres (MAGNUS-P)

A nonlinear numerical model for simulating the propagation of Acoustic-Gravity Waves (AGWs) in planetary atmospheres (e.g., Earth and Mars). 

This codebase is heavily inspired by and validated against the numerical approaches detailed in:
> **Srivastava et al. 2022**: *"A nonlinear numerical model for comparative study of acoustic-gravity wave propagation in planetary atmospheres - Application to Earth and Mars"* (JGR Planets).

## Repository Scope
This repository is configured strictly for core solver development and collaboration. **It tracks only the Julia source code (`src/`) and environment configuration files.** 

Documentation, validation reports, output data, and plotting scripts are maintained separately to keep the primary codebase clean and lightweight.

## Getting Started

This project is written in [Julia](https://julialang.org/). The dependencies and exact package versions are tracked via `Project.toml` and `Manifest.toml` to ensure complete reproducibility.

### Installation
1. Clone the repository:
   ```bash
   git clone git@github.com:AstroJosePC/Solver-Implementation.git
   cd Solver-Implementation
   ```
2. Launch Julia in the project directory:
   ```bash
   julia --project
   ```
3. Instantiate the environment to download exactly matching dependencies:
   ```julia
   julia> using Pkg
   julia> Pkg.instantiate()
   ```

### Running the Solver

The solver is designed to be highly modular. You define the background atmospheric configuration (`earth_config` or `mars_config`), specify the wave forcing parameters (`WaveParams` or `MountainParams`), and pass them to `run_simulation`.

**Example 1: Gaussian Wave Pulse on Mars**
```julia
using MAGNUSP

# 1. Configure the background atmosphere and grid
config = mars_config(
    nx = 300, nz = 400, 
    dx = 500.0, dz = 500.0, 
    CFL = 0.8, t_end = 10000.0,
    output_interval = 50
)

# 2. Define the localized wave forcing (e.g., a pulse)
wave_params = MAGNUSP.WaveParams(
    amplitude = 0.001,
    omega = 0.005,
    kx = 2.0 * pi / 150e3,
    sigma_t = 600.0,
    t0 = 1200.0
)

# 3. Run the simulation
state = run_simulation(config;
    forcing_type = :gaussian,
    wave_params = wave_params,
    outdir = "output/tier1_v1",
    diag_interval = 100
)
```

**Example 2: Orographic Mountain-Wave on Earth**
```julia
using MAGNUSP

# 1. Configure the background atmosphere and grid
config = earth_config(
    nx = 960, nz = 320, 
    dx = 250.0, dz = 250.0,
    CFL = 0.4, t_end = 36000.0
)

# 2. Define the mountain forcing
mountain_params = MAGNUSP.MountainParams(
    h_max = 10.0,            # Hill height (m)
    sigma = 20e3,            # Hill half-width (m)
    x0 = 120e3,              # Center position (m)
    U_bar = 18.0             # Background wind (m/s)
)

# 3. Run the simulation
state = run_simulation(config;
    forcing_type = :mountain,
    wave_params = mountain_params,
    outdir = "output/mountain_benchmark"
)
```

## Codebase Structure
* `src/`: Core Julia implementation.
  * `MAGNUSP.jl`: Main module definition.
  * `euler.jl` / `diffusion.jl` / `forcing.jl`: Core physics and fluid dynamics equations.
  * `boundary.jl` / `sponge.jl` / `damping.jl`: Boundary conditions and numerical damping layers.
  * `grid.jl` / `atmosphere.jl`: Spatial discretization and background atmospheric states.
* `Project.toml` & `Manifest.toml`: Julia package environment configuration.

## License
This project is licensed under the [MIT License](LICENSE).
