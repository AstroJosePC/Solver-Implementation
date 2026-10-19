# =============================================================================
# Mountain Benchmark Post-Processing (Srivastava Sec. 3.3 style)
#
# Inputs: snapshots in output directory (default output/mountain_benchmark)
# Outputs:
#   - Figure-10-style comparison plot (density-scaled w): sim vs theory
#   - Table-11-style scalar metric: RMS(u_sim - u_theory) / Ubar [%]
#   - CSV metrics for report ingestion
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using NCDatasets
using Printf
using Statistics
using Plots

function get_arg(name::String, default)
    key = "--" * name * "="
    for a in ARGS
        if startswith(a, key)
            v = split(a, "=", limit=2)[2]
            return typeof(default) <: AbstractFloat ? parse(Float64, v) :
                   typeof(default) <: Integer ? parse(Int, v) :
                   v
        end
    end
    return default
end

const GAMMA = 1.4
const G = 9.81
const R_GAS = 287.0
const T0 = 239.0

 indir = get_arg("indir", joinpath("output", "mountain_benchmark"))
outdir = get_arg("outdir", joinpath("output", "mountain_benchmark_analysis"))
mkpath(outdir)

files = filter(f -> endswith(f, ".nc"), readdir(indir, join = true))
sort!(files)
isempty(files) && error("No snapshot files found in $(indir)")

f_last = files[end]
ds = NCDataset(f_last, "r")
x = ds["x"][:]
z = ds["z"][:]
rho = ds["rho"][:, :]
u = ds["u"][:, :]
w = ds["w"][:, :]
t_sim = ds.attrib["time_s"]
close(ds)

nx = length(x)
nz = length(z)

# Center x coordinate around hill center (assume center in mid-domain)
x_center = 0.5 * (minimum(x) + maximum(x))
xx = x .- x_center

# Isothermal density ratio for theory scaling
H = R_GAS * T0 / G
rho_ratio = exp.(-z ./ H)      # rho/rho0
scale = sqrt.(rho_ratio)       # sqrt(rho/rho0)

# Queney-like analytic fields (as listed in paper equations)
ks = sqrt(GAMMA - 1.0) * G / sqrt(GAMMA * R_GAS * T0) / Ubar  # N/U

zeta = zeros(nx, nz)
for j in 1:nz, i in 1:nx
    X = xx[i]
    Z = z[j]
    zeta[i, j] = b * sqrt(1.0 / rho_ratio[j]) * (cos(ks * Z) - (X / a) * sin(ks * Z)) / (1.0 + (X / a)^2)
end

# Numerical derivatives
function ddx(field::Matrix{Float64}, x::Vector{Float64})
    nx, nz = size(field)
    out = zeros(nx, nz)
    dx = x[2] - x[1]
    for j in 1:nz
        out[1, j] = (field[2, j] - field[1, j]) / dx
        for i in 2:nx-1
            out[i, j] = (field[i+1, j] - field[i-1, j]) / (2dx)
        end
        out[nx, j] = (field[nx, j] - field[nx-1, j]) / dx
    end
    return out
end

function ddz(field::Matrix{Float64}, z::Vector{Float64})
    nx, nz = size(field)
    out = zeros(nx, nz)
    dz = z[2] - z[1]
    for i in 1:nx
        out[i, 1] = (field[i, 2] - field[i, 1]) / dz
        for j in 2:nz-1
            out[i, j] = (field[i, j+1] - field[i, j-1]) / (2dz)
        end
        out[i, nz] = (field[i, nz] - field[i, nz-1]) / dz
    end
    return out
end

zeta_x = ddx(zeta, xx)
zeta_z = ddz(zeta, z)

w_theory = Ubar .* zeta_x
u_theory = Ubar .* (1 .- zeta_z)

# Density-scaled w for Figure-10-style comparison
rho0_sim = mean(rho[:, 1])
w_sim_scaled = similar(w)
for j in 1:nz, i in 1:nx
    w_sim_scaled[i, j] = w[i, j] * sqrt(rho[i, j] / rho0_sim)
end

w_theory_scaled = similar(w_theory)
for j in 1:nz, i in 1:nx
    w_theory_scaled[i, j] = w_theory[i, j] * scale[j]
end

# RMS metric (Table-11 style): % of Ubar
u_rms = sqrt(mean((u .- u_theory).^2))
u_rms_pct = 100.0 * u_rms / Ubar

# Plot side-by-side contours
x_km = xx ./ 1e3
z_km = z ./ 1e3

clim_abs = max(maximum(abs.(w_sim_scaled)), maximum(abs.(w_theory_scaled)))

p1 = heatmap(x_km, z_km, w_sim_scaled',
    xlabel="x (km)", ylabel="z (km)", title="Simulation (density-scaled w)",
    color=:balance, clims=(-clim_abs, clim_abs))

p2 = heatmap(x_km, z_km, w_theory_scaled',
    xlabel="x (km)", ylabel="z (km)", title="Queney theory (density-scaled w)",
    color=:balance, clims=(-clim_abs, clim_abs))

p = plot(p1, p2, layout=(1,2), size=(1400, 500),
    plot_title=@sprintf("Mountain benchmark at t=%.1f s | u_RMS/Ubar = %.3f%%", t_sim, u_rms_pct))

fig_path = joinpath(outdir, "mountain_fig10_style_comparison.png")
savefig(p, fig_path)

# Write CSV metrics
csv_path = joinpath(outdir, "mountain_metrics.csv")
open(csv_path, "w") do io
    println(io, "snapshot_file,time_s,dx_m,dz_m,u_rms_pct")
    dx = x[2] - x[1]
    dz = z[2] - z[1]
    println(io, "$(basename(f_last)),$(t_sim),$(dx),$(dz),$(u_rms_pct)")
end

println("=" ^ 88)
println("MOUNTAIN BENCHMARK METRICS")
println("=" ^ 88)
println("Snapshot: $(f_last)")
@printf("t = %.1f s, dx = %.1f m, dz = %.1f m\n", t_sim, x[2]-x[1], z[2]-z[1])
@printf("u_RMS as %% of Ubar = %.3f%%\n", u_rms_pct)
println("Saved figure: $(fig_path)")
println("Saved metrics: $(csv_path)")
