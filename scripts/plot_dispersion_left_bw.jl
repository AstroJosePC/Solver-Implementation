# =============================================================================
# Plot grayscale left panel (wave travel diagram) for a given dispersion output
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using NCDatasets
using Statistics
using Plots

function parse_args(args)
    indir = "output/dispersion_test_fig6_strict"
    outfile = "output/dispersion_analysis/dispersion_fig6_left_bw_strict.png"
    for arg in args
        if startswith(arg, "--indir=")
            indir = split(arg, "=", limit=2)[2]
        elseif startswith(arg, "--outfile=")
            outfile = split(arg, "=", limit=2)[2]
        end
    end
    return indir, outfile
end

function list_snapshots(indir::String)
    files = filter(f -> endswith(f, ".nc"), readdir(indir, join=true))
    sort!(files)
    return files
end

function interp_linear(t::Vector{Float64}, y::Vector{Float64}, tu::Vector{Float64})
    yu = similar(tu)
    i = 1
    nt = length(t)
    for k in eachindex(tu)
        tk = tu[k]
        while i < nt-1 && t[i+1] < tk
            i += 1
        end
        if tk <= t[1]
            yu[k] = y[1]
        elseif tk >= t[end]
            yu[k] = y[end]
        else
            t1, t2 = t[i], t[i+1]
            y1, y2 = y[i], y[i+1]
            α = (tk - t1) / (t2 - t1)
            yu[k] = (1 - α) * y1 + α * y2
        end
    end
    return yu
end

let
    indir_rel, outfile_rel = parse_args(ARGS)
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

    w_mid = hcat(profiles...)

    # Upsample for smoother display
    tfine = collect(range(times[1], times[end], length = 4 * length(times)))
    zfine = collect(range(z[1], z[end], length = 4 * length(z)))

    w_time_interp = zeros(length(z), length(tfine))
    for j in eachindex(z)
        w_time_interp[j, :] .= interp_linear(times, vec(w_mid[j, :]), tfine)
    end

    w_fine = zeros(length(zfine), length(tfine))
    for it in eachindex(tfine)
        w_fine[:, it] .= interp_linear(z, vec(w_time_interp[:, it]), zfine)
    end

    wplot = 1e4 .* w_fine
    cabs = maximum(abs, wplot)
    cmax = cabs > 0 ? ceil(cabs) : 1.0
    cb_ticks = collect(-cmax:1.0:cmax)

    p = heatmap(tfine ./ 60.0, zfine ./ 1e3, wplot,
        xlabel = "Time (min)",
        ylabel = "Altitude z (km)",
        title = "Wave travel diagram at mid-x (strict Fig. 6 params, grayscale)",
        color = :grays,
        clims = (-cmax, cmax),
        colorbar_ticks = cb_ticks,
        colorbar_title = "m/s\n×10^-4",
        right_margin = 6Plots.mm,
        left_margin = 8Plots.mm,
        bottom_margin = 8Plots.mm,
        top_margin = 6Plots.mm,
        size = (820, 620)
    )

    savefig(p, outfile)
    println("Saved figure: " * outfile)
end
