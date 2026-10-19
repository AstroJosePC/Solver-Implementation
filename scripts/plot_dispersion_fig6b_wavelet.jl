# =============================================================================
# Plot Figure 6b-style panel using true wavelet transform (CWT, Morse)
# from a given dispersion output directory.
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using NCDatasets
using Statistics
using Plots
using ContinuousWavelets

function parse_args(args)
    indir = "output/dispersion_test_fig6_strict"
    outfile = "output/dispersion_analysis/dispersion_fig6b_wavelet_strict_2x2.png"
    wavelet = "morlet"
    morlet_omega0 = 6.0
    morse_beta = 20.0
    morse_gamma = 3.0
    for arg in args
        if startswith(arg, "--indir=")
            indir = split(arg, "=", limit=2)[2]
        elseif startswith(arg, "--outfile=")
            outfile = split(arg, "=", limit=2)[2]
        elseif startswith(arg, "--wavelet=")
            wavelet = lowercase(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--morlet-omega0=")
            morlet_omega0 = parse(Float64, split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--morse-beta=")
            morse_beta = parse(Float64, split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--morse-gamma=")
            morse_gamma = parse(Float64, split(arg, "=", limit=2)[2])
        end
    end
    return indir, outfile, wavelet, morlet_omega0, morse_beta, morse_gamma
end

function hanning(n::Int)
    n <= 1 && return ones(Float64, n)
    return [0.5 - 0.5*cos(2pi*(i-1)/(n-1)) for i in 1:n]
end

function detrend_linear(y::Vector{Float64})
    n = length(y)
    n <= 1 && return copy(y)
    x = collect(1.0:n)
    x̄ = mean(x)
    ȳ = mean(y)
    num = sum((x .- x̄) .* (y .- ȳ))
    den = sum((x .- x̄).^2)
    b = den > 0 ? num / den : 0.0
    a = ȳ - b * x̄
    return y .- (a .+ b .* x)
end

function process_profile(raw::Vector{Float64}; detrend_mode::Symbol = :mean,
                         apply_window::Bool = false)
    y = copy(raw)
    if detrend_mode == :mean
        y .-= mean(y)
    elseif detrend_mode == :linear
        y = detrend_linear(y)
    else
        error("Unsupported detrend_mode=$(detrend_mode)")
    end
    if apply_window
        y .*= hanning(length(y))
    end
    return y
end

function compute_wavelet_panel(profiles::Vector{Vector{Float64}}, times::Vector{Float64},
                               z::Vector{Float64}, tsel::Vector{Int},
                               lambda_grid::Vector{Float64};
                               zmin::Float64 = minimum(z), zmax::Float64 = maximum(z),
                               wavelet_name::Symbol = :morlet,
                               morlet_omega0::Float64 = 6.0,
                               morse_beta::Float64 = 20.0,
                               morse_gamma::Float64 = 3.0,
                               detrend_mode::Symbol = :mean,
                               apply_window::Bool = false,
                               pool_mode::Symbol = :max,
                               norm_mode::Symbol = :per_time,
                               floor_frac::Float64 = 0.0)

    idx = findall(zz -> zz >= zmin && zz <= zmax, z)
    isempty(idx) && error("No z indices in requested range $(zmin)-$(zmax)")

    zsub = z[idx]
    dz = zsub[2] - zsub[1]
    c = if wavelet_name == :morlet
        CWT(Morlet(morlet_omega0))
    elseif wavelet_name == :morse
        CWT(Morse(morse_gamma, morse_beta, 1.0))
    else
        error("Unsupported wavelet_name=$(wavelet_name). Use :morlet or :morse")
    end
    f_wave = getMeanFreq(length(zsub), c, 1.0 / dz)
    lambda_wave = 1.0 ./ f_wave
    ord = sortperm(lambda_wave)
    lw = lambda_wave[ord]

    panel = zeros(length(lambda_grid), length(tsel))
    rawmax = 0.0

    for (col, it) in enumerate(tsel)
        y = process_profile(profiles[it][idx]; detrend_mode = detrend_mode,
                            apply_window = apply_window)
        W = cwt(y, c)
        Wp = abs2.(W)
        for j in axes(Wp, 2)
            Wp[:, j] .*= f_wave[j]
        end
        pw = if pool_mode == :max
            vec(maximum(Wp, dims = 1))[ord]
        elseif pool_mode == :mean
            vec(mean(Wp, dims = 1))[ord]
        else
            error("Unsupported pool_mode=$(pool_mode)")
        end
        pwi = interp_sorted(lw, pw, lambda_grid)

        pmax = maximum(pwi)
        if floor_frac > 0 && pmax > 0
            pwi[pwi .< floor_frac * pmax] .= 0.0
            pmax = maximum(pwi)
        end

        panel[:, col] .= pwi
        rawmax = max(rawmax, pmax)
    end

    if norm_mode == :per_time
        for col in 1:size(panel, 2)
            pmax = maximum(panel[:, col])
            if pmax > 0
                panel[:, col] ./= pmax
            end
        end
    elseif norm_mode == :global
        if rawmax > 0
            panel ./= rawmax
        end
    else
        error("Unsupported norm_mode=$(norm_mode)")
    end

    panel = clamp.(panel, 0.0, 1.0)
    lambda_dom = [lambda_grid[argmax(panel[:, i])] for i in 1:size(panel, 2)]
    return panel, lambda_dom
end

function list_snapshots(indir::String)
    files = filter(f -> endswith(f, ".nc"), readdir(indir, join=true))
    sort!(files)
    return files
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

let
    indir_rel, outfile_rel, wavelet_str, morlet_omega0, morse_beta, morse_gamma = parse_args(ARGS)
    wavelet_name = Symbol(wavelet_str)
    indir = isabspath(indir_rel) ? indir_rel : joinpath(@__DIR__, "..", indir_rel)
    outfile = isabspath(outfile_rel) ? outfile_rel : joinpath(@__DIR__, "..", outfile_rel)
    mkpath(dirname(outfile))

    files = list_snapshots(indir)
    isempty(files) && error("No snapshots found in $(indir)")

    ds0 = NCDataset(files[1], "r")
    x = ds0["x"][:]
    z = ds0["z"][:]
    close(ds0)

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
    isempty(tsel) && error("No snapshots with t >= 1800 s in $(indir)")

    lambda_grid = collect(range(5e3, 40e3, length = 260))

    # 2x2 denoising/tuning variants (no display gamma/contrast transforms)
    panelA, domA = compute_wavelet_panel(profiles, times, z, tsel, lambda_grid;
        zmin = minimum(z), zmax = maximum(z),
        wavelet_name = wavelet_name,
        morlet_omega0 = morlet_omega0,
        morse_beta = morse_beta,
        morse_gamma = morse_gamma,
        detrend_mode = :mean, apply_window = false,
        pool_mode = :max,
        norm_mode = :per_time, floor_frac = 0.00)

    panelB, domB = compute_wavelet_panel(profiles, times, z, tsel, lambda_grid;
        zmin = 10e3, zmax = 140e3,
        wavelet_name = wavelet_name,
        morlet_omega0 = morlet_omega0,
        morse_beta = morse_beta,
        morse_gamma = morse_gamma,
        detrend_mode = :mean, apply_window = true,
        pool_mode = :max,
        norm_mode = :per_time, floor_frac = 0.00)

    panelC, domC = compute_wavelet_panel(profiles, times, z, tsel, lambda_grid;
        zmin = 10e3, zmax = 140e3,
        wavelet_name = wavelet_name,
        morlet_omega0 = morlet_omega0,
        morse_beta = morse_beta,
        morse_gamma = morse_gamma,
        detrend_mode = :linear, apply_window = true,
        pool_mode = :max,
        norm_mode = :per_time, floor_frac = 0.10)

    panelD, domD = compute_wavelet_panel(profiles, times, z, tsel, lambda_grid;
        zmin = 10e3, zmax = 140e3,
        wavelet_name = wavelet_name,
        morlet_omega0 = morlet_omega0,
        morse_beta = morse_beta,
        morse_gamma = morse_gamma,
        detrend_mode = :linear, apply_window = true,
        pool_mode = :max,
        norm_mode = :global, floor_frac = 0.10)

    tplot = times[tsel] ./ 60.0
    yplot = lambda_grid ./ 1e3

    common_kwargs = (
        xlabel = "Time (min)",
        ylabel = "Vertical wavelength lambda_z (km)",
        color = :turbo,
        clims = (0.0, 1.0),
        colorbar_title = "Normalized power",
        legend = false,
        left_margin = 6Plots.mm,
        right_margin = 6Plots.mm,
        bottom_margin = 7Plots.mm,
        top_margin = 5Plots.mm,
    )

    pA = heatmap(tplot, yplot, panelA;
        title = "A) Baseline: full-z, MAX pooling, per-time norm",
        common_kwargs...)
    plot!(pA, tplot, domA ./ 1e3, color = :white, lw = 2, label = false)

    pB = heatmap(tplot, yplot, panelB;
        title = "B) 10-140 km + Hann window + MAX pooling",
        common_kwargs...)
    plot!(pB, tplot, domB ./ 1e3, color = :white, lw = 2, label = false)

    pC = heatmap(tplot, yplot, panelC;
        title = "C) B + linear detrend + 10% floor + MAX pooling",
        common_kwargs...)
    plot!(pC, tplot, domC ./ 1e3, color = :white, lw = 2, label = false)

    pD = heatmap(tplot, yplot, panelD;
        title = "D) C + global normalization + MAX pooling",
        common_kwargs...)
    plot!(pD, tplot, domD ./ 1e3, color = :white, lw = 2, label = false)

    pfig = plot(pA, pB, pC, pD,
        layout = (2, 2),
        size = (1700, 1200),
        plot_title = "Figure 6b wavelet tuning grid (strict run) | wavelet=$(wavelet_name)")

    savefig(pfig, outfile)
    println("Saved figure: " * outfile)
    println("Snapshots analyzed: $(length(times)), with t>=1800s: $(length(tsel))")
    if wavelet_name == :morlet
        println("Wavelet config: Morlet(omega0=$(morlet_omega0))")
    else
        println("Wavelet config: Morse(gamma=$(morse_gamma), beta=$(morse_beta), center=1.0)")
    end
end
