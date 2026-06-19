"""
    MAGNUSP

2D nonlinear compressible fluid solver for acoustic-gravity wave propagation
in planetary atmospheres. Based on Srivastava et al. (2022).

Uses dimensionally-split Richtmyer Lax-Wendroff for hyperbolic terms and
Crank-Nicholson implicit scheme for viscous/thermal diffusion.
"""
module MAGNUSP

using LinearAlgebra
using SparseArrays
using StaticArrays
using Printf
using FFTW
using Base.Threads: ReentrantLock

# Global lock for thread-safe NetCDF I/O (libnetcdf/libhdf5 are generally not thread-safe)
const IO_LOCK = ReentrantLock()

include("parameters.jl")
include("atmosphere.jl")
include("grid.jl")
include("boundary.jl")
include("euler.jl")
include("diffusion.jl")
include("sponge.jl")
include("damping.jl")
include("forcing.jl")
include("diagnostics.jl")
include("io.jl")
include("run.jl")

export SimConfig, earth_config, mars_config
export Grid2D, init_state
export isothermal_atmosphere, AtmosphereProfile
export brunt_vaisala_local, KlempDurranState, init_klemp_durran
export run_simulation

end # module
