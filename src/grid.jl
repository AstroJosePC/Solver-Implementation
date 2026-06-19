# =============================================================================
# 2D Cartesian grid with ghost cells
# =============================================================================

"""
    Grid2D

Uniform 2D Cartesian grid with ghost cells for boundary conditions.

Physical cells: indices [ng+1 : nx+ng] in x, [ng+1 : nz+ng] in z
Total cells including ghosts: (nx + 2*ng) x (nz + 2*ng)

The state array Q has shape (4, nx_total, nz_total) where:
  Q[1,:,:] = rho       (density)
  Q[2,:,:] = rho * u   (x-momentum)
  Q[3,:,:] = rho * w   (z-momentum)
  Q[4,:,:] = E         (total energy per unit volume)
"""
struct Grid2D
    nx::Int              # physical cells in x
    nz::Int              # physical cells in z
    ng::Int              # ghost cells per side
    dx::Float64          # grid spacing x [m]
    dz::Float64          # grid spacing z [m]
    nx_total::Int        # nx + 2*ng
    nz_total::Int        # nz + 2*ng
    x::Vector{Float64}  # x coordinates of cell centers (physical only)
    z::Vector{Float64}  # z coordinates of cell centers (physical only)
    x_full::Vector{Float64}  # x coordinates including ghosts
    z_full::Vector{Float64}  # z coordinates including ghosts
end

function Grid2D(config::SimConfig)
    ng = config.ng
    nx_total = config.nx + 2 * ng
    nz_total = config.nz + 2 * ng

    # Physical cell centers
    x = [(i - 0.5) * config.dx for i in 1:config.nx]
    z = [(j - 0.5) * config.dz for j in 1:config.nz]

    # Full grid including ghosts (ghost cells extend beyond physical domain)
    x_full = [(i - ng - 0.5) * config.dx for i in 1:nx_total]
    z_full = [(j - ng - 0.5) * config.dz for j in 1:nz_total]

    Grid2D(config.nx, config.nz, ng, config.dx, config.dz,
           nx_total, nz_total, x, z, x_full, z_full)
end

# Index ranges for physical domain
phys_x(g::Grid2D) = (g.ng + 1):(g.nx + g.ng)
phys_z(g::Grid2D) = (g.ng + 1):(g.nz + g.ng)

# Domain extents
domain_width(g::Grid2D) = g.nx * g.dx
domain_height(g::Grid2D) = g.nz * g.dz

"""
    conservative_from_primitive(rho, u, w, P, gamma)

Convert primitive variables to conservative state vector Q = [rho, rho*u, rho*w, E].
E = P/(gamma-1) + 0.5*rho*(u^2 + w^2)
"""
function conservative_from_primitive(rho::Float64, u::Float64, w::Float64,
                                     P::Float64, gamma::Float64)
    E = P / (gamma - 1.0) + 0.5 * rho * (u^2 + w^2)
    return SVector(rho, rho * u, rho * w, E)
end

"""
    primitive_from_conservative(Q, gamma)

Convert conservative variables to primitives (rho, u, w, P).
"""
function primitive_from_conservative(Q::AbstractVector, gamma::Float64)
    rho = Q[1]
    u   = Q[2] / rho
    w   = Q[3] / rho
    KE  = 0.5 * rho * (u^2 + w^2)
    P   = (gamma - 1.0) * (Q[4] - KE)
    return (rho, u, w, P)
end

"""
    temperature_from_conservative(Q, gamma, R_gas)

Extract temperature from conservative variables: T = P / (rho * R_gas).
"""
function temperature_from_conservative(Q::AbstractVector, gamma::Float64, R_gas::Float64)
    rho, u, w, P = primitive_from_conservative(Q, gamma)
    return P / (rho * R_gas)
end

"""
    init_state(grid, atm, config; u0=0.0, w0=0.0)

Initialize the conservative state array Q with background density/pressure
and optional background velocities.
"""
function init_state(grid::Grid2D, atm::AtmosphereProfile, config::SimConfig;
                    u0::Union{Float64, Vector{Float64}, Function} = 0.0, 
                    w0::Union{Float64, Vector{Float64}, Function} = 0.0)
    Q = zeros(Float64, 4, grid.nx_total, grid.nz_total)

    @assert length(atm.z) == grid.nz_total "Atmosphere profile must cover full grid including ghost cells"

    for j in 1:grid.nz_total
        rho = atm.rho0[j]
        P   = atm.P0[j]
        
        uj = if u0 isa Float64
            u0
        elseif u0 isa Vector{Float64}
            u0[j]
        else
            u0(grid.z_full[j])
        end
        
        wj = if w0 isa Float64
            w0
        elseif w0 isa Vector{Float64}
            w0[j]
        else
            w0(grid.z_full[j])
        end

        q = conservative_from_primitive(rho, uj, wj, P, config.gamma)
        for i in 1:grid.nx_total
            Q[1, i, j] = q[1]
            Q[2, i, j] = q[2]
            Q[3, i, j] = q[3]
            Q[4, i, j] = q[4]
        end
    end

    return Q
end
