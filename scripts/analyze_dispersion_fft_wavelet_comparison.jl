# =============================================================================
# Dispersion panel comparison: FFT proxy vs true wavelet spectrum
# Uses dense-output linear dispersion run and produces side-by-side figure.
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using NCDatasets
using FFTW
using Plots
using ContinuousWavelets
using Statistics
using Printf

const INDIR = joinpath(@__DIR__, "..", "output", "dispersion_test_dense")
const OUTDIR = joinpath(@__DIR__, "..", "output", "dispersion_analysis")
mkpath(OUTDIR)

function list_snapshots(indir::String)
    files = filter(f -> endswith(f, ".nc"), readdir(indir, join = true))
    sort!(files)
    return files
end

function fftfreq(n::Int, d::Float64)
    f = zeros(Float64, n)
    n2 = fld(n, 2)
    for i in 0:n-1
        k = i <= n2 ? i : i - n
        f[i+1] = k / (n * d)
    end
    return f
end

function interp_sorted(x::Vector{Float64}, y::Vector{Float64}, xq::Vector{Float64})
    n = length(x)
    n == 0 && return fill(NaN, length(xq))
    n == 1 && return fill(y[1], length(xq))
    yq = similar(xq)
    i = 1
    for k in eachindex(xq)
        xv = xq[k]
        while i < n - 1 && x[i+1] < xv
            i += 1
        end
        if xv <= x[1]
            yq[k] = y[1]
        elseif xv >= x[end]
            yq[k] = y[end]
        else
            x1, x2 = x[i], x[i+1]
            y1, y2 = y[i], y[i+1]
            α = (xv - x1) / (x2 - x1)
            yq[k] = (1 - α) * y1 + α * y2
        end
    end
    return yq
end

function hanning(n::Int)
    n <= 1 && return ones(Float64, n)
    return [0.5 - 0.5*cos(2pi*(i-1)/(n-1)) for i in 1:n]
end

function enhance_contrast(A::AbstractMatrix{<:Real}; p_lo::Float64 = 0.05,
                          p_hi::Float64 = 0.995, gamma::Float64 = 0.8)
    vals = vec(Float64.(A))
    lo = quantile(vals, p_lo)
    hi = quantile(vals, p_hi)
    scale = hi > lo ? (hi - lo) : 1.0
    B = clamp.((Float64.(A) .- lo) ./ scale, 0.0, 1.0)
    return B .^ gamma
end

files = list_snapshots(INDIR)
isempty(files) && error("No NetCDF snapshots in $(INDIR). Run scripts/run_dispersion_test_dense.jl first.")

# Grid from first snapshot
let ds0 = NCDataset(files[1], "r")
    global x = ds0["x"][:]
    global z = ds0["z"][:]
    close(ds0)
end

ix_mid = argmin(abs.(x .- 0.5 * (minimum(x) + maximum(x))))

times = Float64[]
profiles = Vector{Vector{Float64}}()

for f in files
    ds = NCDataset(f, "r")
    w = ds["w"][:, :]
    rho = ds["rho"][:, :]
    t = ds.attrib["time_s"]
    close(ds)

    rho_surface = mean(rho[:, 1])
    wscaled_mid = [w[ix_mid, j] * sqrt(rho[ix_mid, j] / rho_surface) for j in eachindex(z)]

    push!(times, t)
    push!(profiles, wscaled_mid)
end

perm = sortperm(times)
times = times[perm]
profiles = profiles[perm]

tsel = findall(t -> t >= 1800.0, times)
isempty(tsel) && error("No snapshots with t >= 1800 s. Increase output density or run length.")

lambda_grid = collect(range(5e3, 40e3, length = 260))
fft_panel = zeros(length(lambda_grid), length(tsel))
wav_panel = zeros(length(lambda_grid), length(tsel))

# FFT settings
dz = z[2] - z[1]
kz = 2pi .* fftfreq(length(z), dz)
kz_pos_idx = findall(kz .> 1e-8)
wz = hanning(length(z))

# Wavelet settings (true CWT)
c = CWT(morse)
f_wave = getMeanFreq(length(z), c, 1.0 / dz)
lambda_wave = 1.0 ./ f_wave

for (col, it) in enumerate(tsel)
    p = profiles[it] .- mean(profiles[it])

    # FFT proxy spectrum -> lambda
    Y = fft(p .* wz)
    pz = abs.(Y).^2
    kz_pos = kz[kz_pos_idx]
    pz_pos = pz[kz_pos_idx]
    lambda_fft = 2pi ./ kz_pos
    ordf = sortperm(lambda_fft)
    lf = lambda_fft[ordf]
    pf = pz_pos[ordf]
    pfi = interp_sorted(lf, pf, lambda_grid)
    maxf = maximum(pfi)
    fft_panel[:, col] .= maxf > 0 ? pfi ./ maxf : pfi

    # True wavelet spectrum -> lambda
    W = cwt(p, c)
    pw = vec(mean(abs2.(W), dims = 1))
    ordw = sortperm(lambda_wave)
    lw = lambda_wave[ordw]
    pw_s = pw[ordw]
    pwi = interp_sorted(lw, pw_s, lambda_grid)
    maxw = maximum(pwi)
    wav_panel[:, col] .= maxw > 0 ? pwi ./ maxw : pwi
end

# Dominant mode tracks
lambda_fft_dom = [lambda_grid[argmax(fft_panel[:, i])] for i in 1:size(fft_panel, 2)]
lambda_wav_dom = [lambda_grid[argmax(wav_panel[:, i])] for i in 1:size(wav_panel, 2)]

tmin = times[tsel] ./ 60.0

# Contrast-enhanced versions for visualization (preserve raw normalized arrays above)
# fft_panel_vis = enhance_contrast(fft_panel; p_lo = 0.06, p_hi = 0.995, gamma = 0.75)
# wav_panel_vis = enhance_contrast(wav_panel; p_lo = 0.9, p_hi = 0.99999, gamma = 0.5)
fft_panel_vis = fft_panel
wav_panel_vis = wav_panel

p_fft = heatmap(tmin, lambda_grid ./ 1e3, fft_panel_vis,
    xlabel = "Time (min)", ylabel = "Vertical wavelength lambda_z (km)",
    title = "FFT proxy (normalized)",
    color = :turbo, clims = (0.0, 1.0),
    colorbar_title = "Normalized power (enhanced)",
    left_margin = 8Plots.mm, right_margin = 6Plots.mm, bottom_margin = 8Plots.mm)
plot!(p_fft, tmin, lambda_fft_dom ./ 1e3, color = :white, lw = 2, label = false)

p_wav = heatmap(tmin, lambda_grid ./ 1e3, wav_panel_vis,
    xlabel = "Time (min)", ylabel = "Vertical wavelength lambda_z (km)",
    title = "True wavelet (CWT, Morse)",
    color = :turbo, clims = (0.0, 1.0),
    colorbar_title = "Normalized power (enhanced)",
    left_margin = 6Plots.mm, right_margin = 8Plots.mm, bottom_margin = 8Plots.mm,
    legend = false)
plot!(p_wav, tmin, lambda_wav_dom ./ 1e3, color = :white, lw = 2, label = false)

p_cmp = plot(p_fft, p_wav, layout = (1, 2), size = (1520, 600),
    plot_title = "Dispersion right-panel comparison on same dense run")

outfile = joinpath(OUTDIR, "dispersion_fft_vs_wavelet_comparison.png")
savefig(p_cmp, outfile)

@printf("Dense snapshots analyzed: %d (t >= 1800 s: %d)\n", length(times), length(tsel))
println("Saved figure: " * outfile)
