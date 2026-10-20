# =============================================================================
# Equivalence of the implicit diffusion solvers (SimConfig.diffusion_solver) against
# the global sparse LU on the Stage 2 grids.
#
# Builds the 1000 x 190 (dz 500 m) and 1000 x 380 (dz 250 m) grids of
# scripts/sweep_arsia_mons_stage2.jl (dx 4000 m, LZ 95 km, MY31 MCD background with
# the real u_zonal, Mars constants of mars_config), and applies diffusion_step! with
# :lu, :lu_reuse and :fft_tridiag to (a) seeded random perturbations u', w, T' and
# (b) the fields of a real snapshot of each grid. Each solver runs the same three
# steps from the same state: dt0 (first factorization), 1.004 dt0 (inside the 1 %
# tolerance, so the left-hand side keeps dt_ref = dt0 while the right-hand side uses
# the new dt), and 0.8 dt0 (refactorization, the lu! path of :lu_reuse). dt0 is the
# CFL 0.4 step of the state; a second random case uses dt0 = 60 s to load the
# high-wavenumber modes. Reported per step and field (u', w, T'):
#   rel_field = max|f_s - f_lu| / max|f_lu|
#   rel_incr  = max|f_s - f_lu| / max|f_lu - f_before|   (the diffusion increment)
# plus the factorization counts and the wall time of each step.
#
# Usage (from the repository root):
#   OPENBLAS_NUM_THREADS=1 julia -t 1 --project=. scripts/check_diffusion_solvers.jl \
#       [snap500=<snapshot.nc>] [snap250=<snapshot.nc>] [json=<out.json>]
# Defaults: the U40 h1 face-path 24 h records under $MAGNUSP_OUTPUT (default
# ~/MAGNUSP_runs) at about 11.7 h, snapshot_0065800.nc (dz 500 m) and
# snapshot_0142400.nc (dz 250 m, S15F2).
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using MAGNUSP
using NCDatasets
using Printf
using DelimitedFiles
using SparseArrays
using LinearAlgebra
BLAS.set_num_threads(1)

const ATM_FILE = abspath(joinpath(@__DIR__, "..", "external", "mcd_profiles",
                                   "mars_arsia_mons_my31.csv"))
const OUT_ROOT = get(ENV, "MAGNUSP_OUTPUT", joinpath(homedir(), "MAGNUSP_runs"))
const SOLVERS = (:lu, :lu_reuse, :fft_tridiag)

function strarg(key::String, default::String)
    i = findfirst(a -> startswith(a, key * "="), ARGS)
    return i === nothing ? default : String(split(ARGS[i], "=", limit = 2)[2])
end

const SNAP = Dict(
    500.0 => strarg("snap500", joinpath(OUT_ROOT,
        "arsia_mons_stage2_24h_K0_0.005_LZ95_pbarP_ramp2h_bcface/run_U40.0_h1.0/snapshot_0065800.nc")),
    250.0 => strarg("snap250", joinpath(OUT_ROOT,
        "arsia_mons_stage2_24h_K0_0.005_LZ95_pbarP_ramp2h_dz250_bcface/run_U40.0_h1.0/snapshot_0142400.nc")))
const JSON_OUT = strarg("json", "")

"""Linear interpolation of the MCD u_zonal column (as in sweep_arsia_mons_stage2.jl)."""
function u_zonal_closure(atm_file::String)
    data, _ = readdlm(atm_file, ',', Float64, header = true)
    zd, ud = data[:, 1], data[:, 6]
    return z -> z <= zd[1] ? ud[1] : z >= zd[end] ? ud[end] : begin
        k = searchsortedlast(zd, z)
        t = (z - zd[k]) / (zd[k + 1] - zd[k])
        ud[k] * (1 - t) + ud[k + 1] * t
    end
end

"""Seeded uniform numbers in [-1, 1) (SplitMix64), so the random case is reproducible."""
function splitmix_uniform(n::Int, seed::UInt64)
    out = Vector{Float64}(undef, n)
    s = seed
    for k in 1:n
        s += 0x9e3779b97f4a7c15
        z = s
        z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
        z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
        z ⊻= z >> 31
        out[k] = 2.0 * ((z >> 11) * 2.0^-53) - 1.0
    end
    return out
end

"""Stage 2 configuration, grid and background Q_bar [4, nx + 2ng, nz + 2ng] for dz [m]."""
function stage2_setup(dz::Float64; diffusion_solver::Symbol = :lu)
    LZ = 95e3
    config = mars_config(nx = 1000, nz = Int(round(LZ / dz)), dx = 4000.0, dz = dz,
                         CFL = 0.4, enable_viscosity = true, enable_cooling = false,
                         sponge_frac = (LZ - 60e3) / LZ, K0_sponge = 0.005,
                         pbar_mode = :pressure, bottom_forcing = :face,
                         diffusion_solver = diffusion_solver)
    grid = Grid2D(config)
    atm = MAGNUSP.tabulated_atmosphere(ATM_FILE, Float64.(grid.z_full))
    Q_bar = init_state(grid, atm, config, u0 = u_zonal_closure(ATM_FILE))
    return config, grid, Q_bar
end

"""Perturbations (u', w, T') [nx, nz] of Q against Q_bar on the physical cells."""
function perturbations(Q, Q_bar, grid, config)
    nx, nz, ng = grid.nx, grid.nz, grid.ng
    up, wp, Tp = zeros(nx, nz), zeros(nx, nz), zeros(nx, nz)
    g1 = config.gamma - 1.0
    for j in 1:nz, i in 1:nx
        ig, jg = i + ng, j + ng
        rho = Q[1, ig, jg]; u = Q[2, ig, jg] / rho; w = Q[3, ig, jg] / rho
        P = g1 * (Q[4, ig, jg] - 0.5 * rho * (u^2 + w^2))
        rb = Q_bar[1, ig, jg]; ub = Q_bar[2, ig, jg] / rb; wb = Q_bar[3, ig, jg] / rb
        Pb = g1 * (Q_bar[4, ig, jg] - 0.5 * rb * (ub^2 + wb^2))
        up[i, j] = u - ub
        wp[i, j] = w
        Tp[i, j] = P / (rho * config.R_gas) - Pb / (rb * config.R_gas)
    end
    return (up, wp, Tp)
end

"""State from background plus perturbations (u', w, T') [nx, nz] at the background density."""
function state_from_perturbations(Q_bar, grid, config, up, wp, Tp)
    Q = copy(Q_bar)
    g1 = config.gamma - 1.0
    for j in 1:grid.nz, i in 1:grid.nx
        ig, jg = i + grid.ng, j + grid.ng
        rb = Q_bar[1, ig, jg]; ub = Q_bar[2, ig, jg] / rb
        Tb = g1 * (Q_bar[4, ig, jg] - 0.5 * rb * ub^2) / (rb * config.R_gas)
        u = ub + up[i, j]; w = wp[i, j]; P = rb * config.R_gas * (Tb + Tp[i, j])
        Q[1, ig, jg] = rb; Q[2, ig, jg] = rb * u; Q[3, ig, jg] = rb * w
        Q[4, ig, jg] = P / g1 + 0.5 * rb * (u^2 + w^2)
    end
    return Q
end

"""State with the physical cells of a snapshot (rho, u, w, P [nx, nz]); ghosts from Q_bar."""
function state_from_snapshot(path, Q_bar, grid, config)
    Q = copy(Q_bar)
    NCDataset(path) do ds
        rho, u, w, P = (Array{Float64}(ds[v][:, :]) for v in ("rho", "u", "w", "P"))
        size(rho) == (grid.nx, grid.nz) ||
            error("snapshot $(path) is $(size(rho)), grid is $((grid.nx, grid.nz))")
        for j in 1:grid.nz, i in 1:grid.nx
            ig, jg = i + grid.ng, j + grid.ng
            Q[1, ig, jg] = rho[i, j]; Q[2, ig, jg] = rho[i, j] * u[i, j]
            Q[3, ig, jg] = rho[i, j] * w[i, j]
            Q[4, ig, jg] = P[i, j] / (config.gamma - 1.0) + 0.5 * rho[i, j] * (u[i, j]^2 + w[i, j]^2)
        end
    end
    return Q
end

"""Run the three-step sequence with every solver from Q0; return the comparison rows."""
function compare_case(label, dz, Q0, dt0)
    dts = (dt0, 1.004 * dt0, 0.8 * dt0)
    res = Dict{Symbol, Any}()
    for s in SOLVERS
        config, grid, Q_bar = stage2_setup(dz; diffusion_solver = s)
        solver = MAGNUSP.DiffusionSolver(grid.nx, grid.nz)
        Q = copy(Q0)
        fields, walls = Any[perturbations(Q, Q_bar, grid, config)], Float64[]
        for dt in dts
            push!(walls, @elapsed MAGNUSP.diffusion_step!(Q, Q_bar, grid, dt, config, solver))
            all(isfinite, Q) || error("$(label) $(s): non-finite state")
            push!(fields, perturbations(Q, Q_bar, grid, config))
        end
        res[s] = (fields = fields, walls = walls, n_factor = solver.n_factor, dt_ref = solver.dt_ref)
    end
    rows = []
    ref = res[:lu]
    for s in SOLVERS
        r = res[s]
        r.n_factor == ref.n_factor == 2 ||
            error("$(label) $(s): $(r.n_factor) factorizations, expected 2")
        r.dt_ref == ref.dt_ref || error("$(label) $(s): dt_ref differs")
        for n in 1:3, (f, name) in enumerate(("u'", "w", "T'"))
            a, b, before = r.fields[n + 1][f], ref.fields[n + 1][f], ref.fields[n][f]
            d = maximum(abs.(a .- b))
            push!(rows, (case = label, dz = dz, solver = s, step = n, dt = dts[n], field = name,
                         maxabs = d, rel_field = d / maximum(abs.(b)),
                         rel_incr = d / maximum(abs.(b .- before)),
                         wall = r.walls[n], n_factor = r.n_factor))
        end
    end
    return rows
end

# Warm-up on a small grid so the timings below exclude compilation
let config = mars_config(nx = 16, nz = 12, dx = 4000.0, dz = 500.0, enable_viscosity = true)
    for s in SOLVERS
        c = mars_config(nx = 16, nz = 12, dx = 4000.0, dz = 500.0, enable_viscosity = true,
                        diffusion_solver = s)
        grid = Grid2D(c)
        atm = MAGNUSP.isothermal_atmosphere(c, Float64.(grid.z_full))
        Qb = init_state(grid, atm, c, u0 = 10.0)
        sol = MAGNUSP.DiffusionSolver(grid.nx, grid.nz)
        for dt in (1.0, 1.004, 0.8)
            MAGNUSP.diffusion_step!(copy(Qb), Qb, grid, dt, c, sol)
        end
    end
end

all_rows = []
for dz in (500.0, 250.0)
    config, grid, Q_bar = stage2_setup(dz)
    @printf("\n=== grid %d x %d, dx 4000 m, dz %.0f m ===\n", grid.nx, grid.nz, dz)

    # Sparsity pattern of the CN matrix does not depend on the values: sparse() keeps
    # numerical zeros, so a zero diffusivity gives the same structure (lu_reuse premise)
    A1 = MAGNUSP.build_cn_matrix(grid.nx, grid.nz, grid.dx, grid.dz, ones(grid.nz), 1.0)
    A0 = MAGNUSP.build_cn_matrix(grid.nx, grid.nz, grid.dx, grid.dz, zeros(grid.nz), 1.0)
    (A0.colptr == A1.colptr && A0.rowval == A1.rowval) ||
        error("CN sparsity pattern depends on the diffusivity values")
    @printf("CN pattern: nnz = %d for kappa = 1 and kappa = 0 (identical)\n", nnz(A1))

    n = grid.nx * grid.nz
    rnd = reshape.((10.0 .* splitmix_uniform(n, UInt64(1)), 1.0 .* splitmix_uniform(n, UInt64(2)),
                    5.0 .* splitmix_uniform(n, UInt64(3))), grid.nx, grid.nz)
    Q_rnd = state_from_perturbations(Q_bar, grid, config, rnd...)
    dt_cfl = MAGNUSP.compute_adaptive_dt(Q_rnd, grid, 0.4, config.gamma, config.R_gas)

    snap = SNAP[dz]
    isfile(snap) || error("snapshot not found: $(snap)")
    Q_snap = state_from_snapshot(snap, Q_bar, grid, config)
    dt_snap = MAGNUSP.compute_adaptive_dt(Q_snap, grid, 0.4, config.gamma, config.R_gas)
    t_snap = NCDataset(ds -> ds.attrib["time_s"], snap)
    @printf("snapshot %s (t = %.3f h)\n", snap, t_snap / 3600)

    for (label, Q0, dt0) in (("random_cfl", Q_rnd, dt_cfl), ("random_dt60", Q_rnd, 60.0),
                             ("snapshot_cfl", Q_snap, dt_snap))
        rows = compare_case(label, dz, Q0, dt0)
        append!(all_rows, rows)
        @printf("\n%-13s dt0 = %.4g s\n", label, dt0)
        @printf("  %-11s %4s %3s %10s %10s %10s %8s\n", "solver", "step", "f",
                "max|diff|", "rel_field", "rel_incr", "wall[s]")
        for r in rows
            r.solver === :lu && continue
            @printf("  %-11s %4d %3s %10.3e %10.3e %10.3e %8.3f\n", r.solver, r.step,
                    r.field, r.maxabs, r.rel_field, r.rel_incr, r.wall)
        end
        for s in SOLVERS
            w = [r.wall for r in rows if r.solver === s && r.field == "u'"]
            @printf("  wall per step %-11s %s s (factorizations %d)\n", s,
                    join((@sprintf("%.3f", x) for x in w), " / "),
                    first(r.n_factor for r in rows if r.solver === s))
        end
    end
end

println("\nMaximum over steps and fields:")
for dz in (500.0, 250.0), s in (:lu_reuse, :fft_tridiag)
    sel = [r for r in all_rows if r.dz == dz && r.solver === s]
    @printf("  dz %3.0f m %-11s rel_field %.3e  rel_incr %.3e\n", dz, s,
            maximum(r.rel_field for r in sel), maximum(r.rel_incr for r in sel))
end

if !isempty(JSON_OUT)
    open(JSON_OUT, "w") do io
        println(io, "{\"snapshots\": {\"dz500\": \"$(SNAP[500.0])\", \"dz250\": \"$(SNAP[250.0])\"},")
        println(io, " \"rows\": [")
        for (k, r) in enumerate(all_rows)
            @printf(io, "  {\"case\": \"%s\", \"dz_m\": %.1f, \"solver\": \"%s\", \"step\": %d, \"dt_s\": %.17g, \"field\": \"%s\", \"maxabs\": %.6e, \"rel_field\": %.6e, \"rel_incr\": %.6e, \"wall_s\": %.6f, \"n_factor\": %d}%s\n",
                    r.case, r.dz, r.solver, r.step, r.dt, r.field, r.maxabs, r.rel_field,
                    r.rel_incr, r.wall, r.n_factor, k == length(all_rows) ? "" : ",")
        end
        println(io, " ]}")
    end
    println("JSON: ", JSON_OUT)
end
