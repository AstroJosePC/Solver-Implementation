# =============================================================================
# Rayleigh damping sponge layer
# Srivastava et al. (2022), Equation 27 (Klemp & Lilly 1978)
# =============================================================================

"""
    SpongeLayer

Pre-computed sponge damping coefficients for each grid cell.
"""
struct SpongeLayer
    K::Matrix{Float64}  # damping coefficient at each (i,j), size (nx_total, nz_total)
end

"""
    init_sponge(grid, config)

Initialize sponge layer coefficients.

K(z) = K0 * sin^2( pi/2 * (z - z_sponge) / (z_max - z_sponge) )  for z >= z_sponge
K(z) = 0  for z < z_sponge

The sponge occupies the top `sponge_frac` fraction of the domain.
Optimal K0 = 0.005 gives reflection coefficient r ~ 0.55% (Srivastava).
"""
function init_sponge(grid::Grid2D, config::SimConfig)
    z_max = domain_height(grid)
    z_sponge = z_max * (1.0 - config.sponge_frac)

    K = zeros(grid.nx_total, grid.nz_total)

    for j in 1:grid.nz_total
        z_j = grid.z_full[j]
        if z_j >= z_sponge && z_max > z_sponge
            eta = (z_j - z_sponge) / (z_max - z_sponge)
            K_val = config.K0_sponge * sin(0.5 * pi * eta)^2
            for i in 1:grid.nx_total
                K[i, j] = K_val
            end
        end
    end

    return SpongeLayer(K)
end

"""
    apply_sponge!(Q, Q_bar, sponge, dt)

Apply Rayleigh damping: relax Q toward reference state Q_bar.
Q -= K * dt * (Q - Q_bar)
"""
function apply_sponge!(Q::Array{Float64,3}, Q_bar::Array{Float64,3},
                       sponge::SpongeLayer, dt::Float64)
    K = sponge.K
    for j in axes(Q, 3), i in axes(Q, 2)
        if K[i, j] > 0.0
            damping = K[i, j] * dt
            for k in 1:4
                Q[k, i, j] -= damping * (Q[k, i, j] - Q_bar[k, i, j])
            end
        end
    end
    return nothing
end
