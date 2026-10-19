# =============================================================================
# KE(z,t) Hovmoller Diagnostics for Reflection Analysis
#
# Runs selected cases and builds Hovmoller fields:
#   KEz(z,t) = sum_x [0.5 * rho_bar * (u^2 + w^2)] * dx
#
# Default cases (dx=500 m):
#   1) WB,  K0=0.0
#   2) STD, K0=0.0
#   3) STD, K0=0.005
#
# Outputs:
#   - NetCDF files with KE(z,t) per case
#   - PNG heatmaps if Plots.jl is available
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using MAGNUSP
using MAGNUSP: phys_x, phys_z, SVector, primitive_from_conservative
using Printf
using NCDatasets

const LX = 80_000.0
const LZ = 120_000.0
const DX = 500.0
const DZ = 500.0
const NX = Int(round(LX / DX))
const NZ = Int(round(LZ / DZ))

const T_ISO = 239.0
const CFL = 0.8
const T_END = 8000.0
const DIAG_DT = 20.0
const SPONGE_KM = 20.0

const W_TILDE = 0.001
const T0_FORCING = 1200.0
const SIGMA_T = 600.0
const OMEGA = 0.01
const LAMBDA_X = 80e3
const KX = 2.0 * pi / LAMBDA_X

struct CaseCfg
    label::String
    use_standard::Bool
    K0::Float64
end

const CASES = [
    CaseCfg("WB_K0_0", false, 0.0),
    CaseCfg("STD_K0_0", true, 0.0),
    CaseCfg("STD_K0_0005", true, 0.005),
]

function compute_flux_z_standard(Q::SVector{4,Float64}, gamma::Float64)
    rho, u, w, P = primitive_from_conservative(Q, gamma)
    E = Q[4]
    return SVector(rho * w, rho * u * w, rho * w^2 + P, w * (E + P))
end

function compute_source_standard(Q::SVector{4,Float64}, g::Float64)
    rho = Q[1]
    w = Q[3] / rho
    return SVector(0.0, 0.0, -rho * g, -rho * g * w)
end

function euler_step_standard!(Q::Array{Float64,3}, Q_bar::Array{Float64,3},
                              Q_half_x::Array{Float64,3}, Q_half_z::Array{Float64,3},
                              grid::Grid2D, dt::Float64, gamma::Float64,
                              g::Float64, sponge_K::Matrix{Float64})
    dtdx = dt / grid.dx
    dtdz = dt / grid.dz
    ix = phys_x(grid)
    iz = phys_z(grid)

    for j in iz
        for i in (first(ix)-1):(last(ix))
            Qi   = SVector(Q[1,i,j], Q[2,i,j], Q[3,i,j], Q[4,i,j])
            Qip1 = SVector(Q[1,i+1,j], Q[2,i+1,j], Q[3,i+1,j], Q[4,i+1,j])
            Fi   = MAGNUSP.compute_flux_x(Qi, gamma)
            Fip1 = MAGNUSP.compute_flux_x(Qip1, gamma)
            @inbounds for k in 1:4
                Q_half_x[k, i, j] = 0.5 * (Qi[k] + Qip1[k]) - 0.5 * dtdx * (Fip1[k] - Fi[k])
            end
        end
    end

    for j in iz
        for i in ix
            Qh_right = SVector(Q_half_x[1,i,j], Q_half_x[2,i,j], Q_half_x[3,i,j], Q_half_x[4,i,j])
            Qh_left  = SVector(Q_half_x[1,i-1,j], Q_half_x[2,i-1,j], Q_half_x[3,i-1,j], Q_half_x[4,i-1,j])
            Fh_right = MAGNUSP.compute_flux_x(Qh_right, gamma)
            Fh_left  = MAGNUSP.compute_flux_x(Qh_left, gamma)
            @inbounds for k in 1:4
                Q[k, i, j] -= dtdx * (Fh_right[k] - Fh_left[k])
            end
        end
    end

    for j in (first(iz)-1):(last(iz))
        for i in ix
            Qj   = SVector(Q[1,i,j], Q[2,i,j], Q[3,i,j], Q[4,i,j])
            Qjp1 = SVector(Q[1,i,j+1], Q[2,i,j+1], Q[3,i,j+1], Q[4,i,j+1])
            Gj   = compute_flux_z_standard(Qj, gamma)
            Gjp1 = compute_flux_z_standard(Qjp1, gamma)
            Sj   = compute_source_standard(Qj, g)
            Sjp1 = compute_source_standard(Qjp1, g)
            Kj   = sponge_K[i, j]
            Kjp1 = sponge_K[i, j+1]
            @inbounds for k in 1:4
                sp_j   = Kj   * (Q_bar[k, i, j]   - Qj[k])
                sp_jp1 = Kjp1 * (Q_bar[k, i, j+1] - Qjp1[k])
                Q_half_z[k, i, j] = 0.5 * (Qj[k] + Qjp1[k]) -
                                    0.5 * dtdz * (Gjp1[k] - Gj[k]) +
                                    0.25 * dt * (Sj[k] + Sjp1[k] + sp_j + sp_jp1)
            end
        end
    end

    for j in iz
        for i in ix
            Qh_top = SVector(Q_half_z[1,i,j], Q_half_z[2,i,j], Q_half_z[3,i,j], Q_half_z[4,i,j])
            Qh_bot = SVector(Q_half_z[1,i,j-1], Q_half_z[2,i,j-1], Q_half_z[3,i,j-1], Q_half_z[4,i,j-1])
            Gh_top = compute_flux_z_standard(Qh_top, gamma)
            Gh_bot = compute_flux_z_standard(Qh_bot, gamma)
            Sh_top = compute_source_standard(Qh_top, g)
            Sh_bot = compute_source_standard(Qh_bot, g)
            K_top = 0.5 * (sponge_K[i, j]   + sponge_K[i, j+1])
            K_bot = 0.5 * (sponge_K[i, j-1] + sponge_K[i, j])
            @inbounds for k in 1:4
                Qbar_half_top = 0.5 * (Q_bar[k, i, j]   + Q_bar[k, i, j+1])
                Qbar_half_bot = 0.5 * (Q_bar[k, i, j-1] + Q_bar[k, i, j])
                sp_top = K_top * (Qbar_half_top - Qh_top[k])
                sp_bot = K_bot * (Qbar_half_bot - Qh_bot[k])
                Q[k, i, j] -= dtdz * (Gh_top[k] - Gh_bot[k])
                Q[k, i, j] += dt * 0.5 * (Sh_top[k] + Sh_bot[k] + sp_top + sp_bot)
            end
        end
    end

    return nothing
end

function run_case(cfg::CaseCfg)
    sponge_frac = SPONGE_KM / (LZ / 1e3)

    config = earth_config(
        nx = NX, nz = NZ, dx = DX, dz = DZ,
        CFL = CFL, T0 = T_ISO, t_end = T_END,
        output_interval = 999999,
        enable_viscosity = false,
        enable_cooling = false,
        sponge_frac = sponge_frac,
        K0_sponge = cfg.K0,
    )

    wave_params = MAGNUSP.WaveParams(
        amplitude = W_TILDE,
        omega = OMEGA,
        kx = KX,
        sigma_t = SIGMA_T,
        t0 = T0_FORCING,
    )

    grid = Grid2D(config)
    atm = isothermal_atmosphere(config, Float64.(grid.z_full))
    Q = MAGNUSP.init_state(grid, atm, config)
    Q_bar = copy(Q)
    sponge_layer = MAGNUSP.init_sponge(grid, config)
    Q_hx = similar(Q)
    Q_hz = similar(Q)

    z_vec = [grid.z_full[j] for j in phys_z(grid)]
    times = Float64[]
    kez_snapshots = Vector{Vector{Float64}}()

    t = 0.0
    steps = 0
    next_diag = 0.0
    crashed = false

    @printf("\n[%s] dx=%.0f m, K0=%.4f\n", cfg.label, DX, cfg.K0)

    while t < config.t_end
        local dt::Float64
        try
            dt = MAGNUSP.compute_adaptive_dt(Q, grid, config.CFL, config.gamma, config.R_gas)
            if !isfinite(dt) || dt <= 0.0
                crashed = true
                @printf("  crashed at t=%.1f s (invalid dt=%.3e)\n", t, dt)
                break
            end
            dt = min(dt, config.t_end - t)

            MAGNUSP.apply_all_bcs!(Q, Q_bar, grid, config.gamma)
            MAGNUSP.apply_bottom_forcing!(Q, grid, config, t, wave_params, :gaussian)

            if cfg.use_standard
                euler_step_standard!(Q, Q_bar, Q_hx, Q_hz, grid, dt, config.gamma, config.g,
                                     sponge_layer.K)
            else
                MAGNUSP.euler_step!(Q, Q_bar, Q_hx, Q_hz, grid, dt, config.gamma, config.g,
                                    sponge_layer.K)
            end
        catch err
            crashed = true
            @printf("  crashed at t=%.1f s (%s)\n", t, sprint(showerror, err))
            break
        end

        t += dt
        steps += 1

        if t >= next_diag
            ke_z = zeros(Float64, grid.nz)
            for (j_local, j) in enumerate(phys_z(grid))
                acc = 0.0
                for i in phys_x(grid)
                    rho = Q[1, i, j]
                    u = Q[2, i, j] / rho
                    w = Q[3, i, j] / rho
                    acc += 0.5 * Q_bar[1, i, j] * (u^2 + w^2) * DX
                end
                ke_z[j_local] = acc
            end

            push!(times, t)
            push!(kez_snapshots, ke_z)
            next_diag += DIAG_DT
        end
    end

    if crashed || isempty(kez_snapshots)
        return (false, z_vec, Float64[], zeros(Float64, length(z_vec), 0), NaN)
    end

    nt = length(kez_snapshots)
    nz = length(z_vec)
    kez = zeros(Float64, nz, nt)
    for n in 1:nt
        kez[:, n] = kez_snapshots[n]
    end

    ke_total = vec(sum(kez, dims=1))
    ke_peak = maximum(ke_total)
    idx5000 = argmin(abs.(times .- 5000.0))
    r = sqrt(ke_total[idx5000] / ke_peak) * 100.0

    @printf("  done: steps=%d, nt=%d, r=%.3f%%\n", steps, nt, r)
    return (true, z_vec, times, kez, r)
end

function save_case_netcdf(path::String, label::String, z::Vector{Float64}, t::Vector{Float64},
                          kez::Matrix{Float64}, r::Float64)
    NCDataset(path, "c") do ds
        defDim(ds, "z", length(z))
        defDim(ds, "time", length(t))

        ds.attrib["description"] = "KE(z,t) Hovmoller"
        ds.attrib["case_label"] = label
        ds.attrib["reflection_r_percent"] = r

        zvar = defVar(ds, "z", Float64, ("z",), attrib = Dict("units" => "m"))
        tvar = defVar(ds, "time", Float64, ("time",), attrib = Dict("units" => "s"))
        kvar = defVar(ds, "KE_z_t", Float64, ("z", "time"),
                      attrib = Dict("units" => "J/m^2", "long_name" => "column-integrated KE at each z"))

        zvar[:] = z
        tvar[:] = t
        kvar[:, :] = kez
    end
end

function maybe_plot_hovmoller_png(outpath::String, label::String, z::Vector{Float64}, t::Vector{Float64},
                                  kez::Matrix{Float64}, r::Float64)
    try
        @eval import Plots
    catch
        println("Plots.jl not available in this environment. NetCDF outputs were saved.")
        return
    end

    z_km = z ./ 1e3
    kez_plot = log10.(kez .+ maximum(kez) * 1e-12)

    p = Base.invokelatest(
        Plots.heatmap,
        t,
        z_km,
        kez_plot;
        xlabel = "time (s)",
        ylabel = "z (km)",
        title = "$(label)  |  log10 KE(z,t), r=$(round(r, digits=2))%",
        colorbar_title = "log10 KE",
        ylims = (minimum(z_km), maximum(z_km)),
    )

    Base.invokelatest(Plots.savefig, p, outpath)
end

println("=" ^ 84)
println("KE(z,t) HOVMOLLER DIAGNOSTICS")
println("=" ^ 84)

outdir = joinpath(@__DIR__, "..", "output", "hovmoller_ke")
mkpath(outdir)

for cfg in CASES
    ok, z, t, kez, r = run_case(cfg)
    if !ok
        println("Case $(cfg.label) failed; skipping outputs.")
        continue
    end

    nc_path = joinpath(outdir, "KE_zt_$(cfg.label).nc")
    png_path = joinpath(outdir, "KE_zt_$(cfg.label).png")

    save_case_netcdf(nc_path, cfg.label, z, t, kez, r)
    maybe_plot_hovmoller_png(png_path, cfg.label, z, t, kez, r)

    @printf("Saved: %s\n", nc_path)
    if isfile(png_path)
        @printf("Saved: %s\n", png_path)
    end
end

println("\nDone.")
