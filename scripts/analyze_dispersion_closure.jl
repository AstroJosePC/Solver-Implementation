# =============================================================================
# Linear Dispersion Closure Post-Processing (Spectral)
#
# Robust spectral method inspired by Srivastava's spectral analysis approach:
#   - Temporal FFT at probe altitude -> ω_sim, T_sim
#   - Phase at same ω for two x-probes -> kx_sim, cx_sim
#   - 2D FFT (x-z) at energetic snapshot, isolate forced |kx| -> kz_sim, λz_sim
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using NCDatasets
using Printf
using Statistics
using FFTW
using Plots

# Match run_dispersion_test.jl settings (adjust if changed there)
const GAMMA = 1.4
const G = 9.81
const R_GAS = 287.0
const T0 = 239.0
const OMEGA = 0.007
const LAMBDA_X = 40e3
const KX = 2pi / LAMBDA_X

const INDIR = joinpath(@__DIR__, "..", "output", "dispersion_test")
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

function hanning(n::Int)
    n <= 1 && return ones(Float64, n)
    return [0.5 - 0.5*cos(2pi*(i-1)/(n-1)) for i in 1:n]
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
    files = list_snapshots(INDIR)
    if isempty(files)
        error("No NetCDF snapshots found in $(INDIR). Run scripts/run_dispersion_test.jl first.")
    end

    # Grid from first snapshot
    ds0 = NCDataset(files[1], "r")
    x = ds0["x"][:]
    z = ds0["z"][:]
    close(ds0)

    # Probes
    ix1 = argmin(abs.(x .- (minimum(x) + 0.25 * (maximum(x)-minimum(x)))))
    ix2 = argmin(abs.(x .- (minimum(x) + 0.50 * (maximum(x)-minimum(x)))))
    jp  = argmin(abs.(z .- 60e3))

    times = Float64[]
    p1 = Float64[]
    p2 = Float64[]
    energ = Float64[]
    fields = Matrix{Float64}[]

    # Read snapshots
    for f in files
        ds = NCDataset(f, "r")
        w = ds["w"][:, :]
        rho = ds["rho"][:, :]
        t = ds.attrib["time_s"]
        close(ds)

        push!(times, t)
        push!(p1, w[ix1, jp])
        push!(p2, w[ix2, jp])

        # Build density-scaled w field for spectral use
        rho_surface = mean(rho[:, 1])
        wscaled = similar(w)
        for j in axes(w, 2), i in axes(w, 1)
            wscaled[i, j] = w[i, j] * sqrt(rho[i, j] / rho_surface)
        end
        push!(fields, wscaled)
        push!(energ, sqrt(mean(abs2, wscaled)))
    end

    # Sort by time
    perm = sortperm(times)
    times = times[perm]
    p1 = p1[perm]
    p2 = p2[perm]
    fields = fields[perm]
    energ = energ[perm]

    # Uniform time grid for FFT
    dtu = median(diff(times))
    tu = collect(range(times[1], times[end], step = dtu))
    y1 = interp_linear(times, p1, tu)
    y2 = interp_linear(times, p2, tu)

    # Window + FFT
    wt = hanning(length(tu))
    y1w = (y1 .- mean(y1)) .* wt
    y2w = (y2 .- mean(y2)) .* wt
    Y1 = fft(y1w)
    Y2 = fft(y2w)
    f = fftfreq(length(tu), dtu)
    ω = 2pi .* f

    # positive frequencies only
    pos = findall(ω .> 0)
    ωpos = ω[pos]
    Pω = abs.(Y1[pos]).^2

    # pick dominant ω near target forcing (±40%)
    cand = findall((ωpos .>= 0.6*OMEGA) .& (ωpos .<= 1.4*OMEGA))
    if isempty(cand)
        iω = pos[argmax(Pω)]
    else
        ic = cand[argmax(Pω[cand])]
        iω = pos[ic]
    end
    ω_sim = ω[iω]
    T_sim = 2pi / ω_sim

    # Phase-based kx and cx from two probes at same ω
    φ1 = angle(Y1[iω])
    φ2 = angle(Y2[iω])
    dφ = φ2 - φ1
    while dφ > pi
        dφ -= 2pi
    end
    while dφ < -pi
        dφ += 2pi
    end
    dx12 = abs(x[ix2] - x[ix1])
    kx_sim = abs(dφ) > 0 ? abs(dφ) / dx12 : NaN
    cx_sim = (isfinite(kx_sim) && kx_sim > 0) ? (ω_sim / kx_sim) : NaN

    # 2D FFT in (x,z): choose energetic snapshot after transients
    i_start = findfirst(t -> t >= 1800.0, times)
    i_start = isnothing(i_start) ? 1 : i_start
    i_best = i_start - 1 + argmax(energ[i_start:end])
    W = fields[i_best]

    # Restrict z analysis band to 20-120 km
    jlo = argmin(abs.(z .- 20e3))
    jhi = argmin(abs.(z .- 120e3))
    Wb = W[:, jlo:jhi]

    # Detrend/center + windows
    Wb = Wb .- mean(Wb)
    wx = hanning(size(Wb, 1))
    wz = hanning(size(Wb, 2))
    for j in 1:size(Wb,2), i in 1:size(Wb,1)
        Wb[i,j] *= wx[i] * wz[j]
    end

    F2 = fft(Wb)
    kx_vec = 2pi .* fftfreq(size(Wb,1), x[2]-x[1])
    kz_vec = 2pi .* fftfreq(size(Wb,2), z[2]-z[1])

    # Isolate |kx| near forcing kx and find dominant kz by power
    ikx = argmin(abs.(abs.(kx_vec) .- KX))
    pkz = abs.(F2[ikx, :]).^2
    # exclude kz ~ 0
    nz_idx = findall(abs.(kz_vec) .> 1e-8)
    ikz = nz_idx[argmax(pkz[nz_idx])]
    kz_sim = kz_vec[ikz]
    lambda_z_sim = 2pi / abs(kz_sim)

    # Track dominant lambda_z(t) from 2D FFT at forced |kx|
    lambda_z_t = fill(NaN, length(times))
    for it in eachindex(times)
        Wti = fields[it][:, jlo:jhi]
        Wti = Wti .- mean(Wti)
        for j in 1:size(Wti,2), i in 1:size(Wti,1)
            Wti[i,j] *= wx[i] * wz[j]
        end
        Fti = fft(Wti)
        pkzi = abs.(Fti[ikx, :]).^2
        ikzi = nz_idx[argmax(pkzi[nz_idx])]
        kz_i = kz_vec[ikzi]
        lambda_z_t[it] = abs(kz_i) > 0 ? 2pi / abs(kz_i) : NaN
    end

    # Theory
    cs = sqrt(GAMMA * R_GAS * T0)
    N = sqrt((GAMMA - 1.0) * G^2 / cs^2)
    m2 = KX^2 * (N^2 / OMEGA^2 - 1.0)
    lambda_z_theory = m2 > 0 ? 2pi / sqrt(m2) : NaN
    T_theory = 2pi / OMEGA
    cx_theory = OMEGA / KX

    err_lz = isfinite(lambda_z_sim) ? 100 * abs(lambda_z_sim - lambda_z_theory) / lambda_z_theory : NaN
    err_T = isfinite(T_sim) ? 100 * abs(T_sim - T_theory) / T_theory : NaN
    err_cx = isfinite(cx_sim) ? 100 * abs(cx_sim - cx_theory) / cx_theory : NaN

    println("=" ^ 92)
    println("LINEAR DISPERSION CLOSURE SUMMARY (SPECTRAL)")
    println("=" ^ 92)
    @printf("Files analyzed: %d, energetic snapshot t=%.1f s\n", length(files), times[i_best])
    @printf("Probe z = %.1f km, x1 = %.1f km, x2 = %.1f km\n", z[jp]/1e3, x[ix1]/1e3, x[ix2]/1e3)
    println()
    @printf("%-24s %-14s %-14s %-12s\n", "Metric", "Theory", "Simulated", "|err| %")
    println("-" ^ 92)
    @printf("%-24s %-14.3f %-14.3f %-12.2f\n", "lambda_z (km)", lambda_z_theory/1e3, lambda_z_sim/1e3, err_lz)
    @printf("%-24s %-14.3f %-14.3f %-12.2f\n", "Period T (s)", T_theory, T_sim, err_T)
    @printf("%-24s %-14.3f %-14.3f %-12.2f\n", "c_x (m/s)", cx_theory, cx_sim, err_cx)
    println()
    @printf("Recovered ω_sim = %.6f rad/s, kx_sim = %.6e 1/m, kz_sim = %.6e 1/m\n", ω_sim, kx_sim, kz_sim)
    println("Heuristic closure pass: all three errors <= 10%")
    pass = all(isfinite.([err_lz, err_T, err_cx])) && err_lz <= 10 && err_T <= 10 && err_cx <= 10
    println(pass ? "PASS" : "CHECK")

    # Save CSV metrics for report tables
    csv_path = joinpath(OUTDIR, "dispersion_closure_metrics.csv")
    open(csv_path, "w") do io
        println(io, "metric,theory,simulated,error_percent")
        println(io, "lambda_z_km,$(lambda_z_theory/1e3),$(lambda_z_sim/1e3),$(err_lz)")
        println(io, "period_s,$(T_theory),$(T_sim),$(err_T)")
        println(io, "cx_mps,$(cx_theory),$(cx_sim),$(err_cx)")
        println(io, "omega_rads,$(OMEGA),$(ω_sim),$(100*abs(ω_sim-OMEGA)/OMEGA)")
    end

    # Plot 1: λz evolution (FFT-based) vs theory
    p1 = plot(times ./ 60.0, lambda_z_t ./ 1e3,
              xlabel = "time (min)", ylabel = "lambda_z (km)",
              label = "FFT dominant lambda_z", lw = 2,
              title = "Dispersion closure: dominant vertical wavelength vs time")
    hline!(p1, [lambda_z_theory/1e3], label = "Theory lambda_z", ls = :dash)
    savefig(p1, joinpath(OUTDIR, "dispersion_lambda_z_vs_time.png"))

    # Plot 2: k_z power slice at forced |k_x| for energetic snapshot
    p2 = plot(kz_vec .* 1e3, pkz ./ maximum(pkz),
              xlabel = "kz (1/km)", ylabel = "normalized power",
              label = "|kx| ~ forcing", lw = 2,
              title = @sprintf("kz spectrum at t=%.1f s", times[i_best]))
    vline!(p2, [kz_sim*1e3], label = "kz_sim", ls = :dash)
    savefig(p2, joinpath(OUTDIR, "dispersion_kz_spectrum.png"))

    # Figure-6-style panel (a): wave travel diagram at domain mid-x using density-scaled w
    ix_mid = argmin(abs.(x .- 0.5 * (minimum(x) + maximum(x))))
    w_mid = hcat([fields[it][ix_mid, :] for it in eachindex(times)]...)

    # Upsample for smoother display (paper-like rendering while preserving raw data values)
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

    # Scale for cleaner colorbar tick labels (~integer) with explicit 10^-4 factor
    w_mid_scaled = 1e4 .* w_fine
    cabs = maximum(abs, w_mid_scaled)
    cmax = cabs > 0 ? ceil(cabs) : 1.0
    cb_ticks = collect(-cmax:1.0:cmax)

    pa = heatmap(tfine ./ 60.0, zfine ./ 1e3, w_mid_scaled,
                 xlabel = "Time (min)", ylabel = "Altitude z (km)",
                 title = "Wave travel diagram at mid-x",
                 color = :balance,
                 clims = (-cmax, cmax),
                 colorbar_ticks = cb_ticks,
                 colorbar_title = "m/s\n×10^-4",
                 right_margin = 6Plots.mm,
                 left_margin = 8Plots.mm,
                 bottom_margin = 8Plots.mm,
                 top_margin = 6Plots.mm)

    # Grayscale variant for apples-to-apples visual comparison with paper styling
    pa_bw = heatmap(tfine ./ 60.0, zfine ./ 1e3, w_mid_scaled,
                    xlabel = "Time (min)", ylabel = "Altitude z (km)",
                    title = "Wave travel diagram at mid-x (grayscale)",
                    color = :grays,
                    clims = (-cmax, cmax),
                    colorbar_ticks = cb_ticks,
                    colorbar_title = "m/s\n×10^-4",
                    right_margin = 6Plots.mm,
                    left_margin = 8Plots.mm,
                    bottom_margin = 8Plots.mm,
                    top_margin = 6Plots.mm)

    # Figure-6-style panel (b): normalized FFT power vs vertical wavelength and time (proxy to wavelet)
    tsel = findall(t -> t >= 1800.0, times)
    dz = z[2] - z[1]
    kz_full = 2pi .* fftfreq(length(z), dz)
    nz_idx_full = findall(abs.(kz_full) .> 1e-8)
    kz_pos_idx = [k for k in nz_idx_full if kz_full[k] > 0]

    lambda_grid = collect(range(5e3, 40e3, length = 220))
    spec_lambda_t = fill(0.0, length(lambda_grid), length(tsel))
    wz_full = hanning(length(z))

    for (col, it) in enumerate(tsel)
        profile = w_mid[:, it] .- mean(w_mid[:, it])
        Yz = fft(profile .* wz_full)
        pz = abs.(Yz).^2

        kz_pos = kz_full[kz_pos_idx]
        pz_pos = pz[kz_pos_idx]
        lambda_pos = 2pi ./ kz_pos

        ord = sortperm(lambda_pos)
        lambda_sorted = lambda_pos[ord]
        pz_sorted = pz_pos[ord]

        pz_interp = interp_sorted(lambda_sorted, pz_sorted, lambda_grid)
        maxp = maximum(pz_interp)
        spec_lambda_t[:, col] .= maxp > 0 ? pz_interp ./ maxp : pz_interp
    end

    # Dominant λz(t) curve from panel (b)
    lambda_dom = similar(tsel, Float64)
    for col in eachindex(tsel)
        imax = argmax(spec_lambda_t[:, col])
        lambda_dom[col] = lambda_grid[imax]
    end

    pb = heatmap(times[tsel] ./ 60.0, lambda_grid ./ 1e3, spec_lambda_t,
                 xlabel = "Time (min)", ylabel = "Vertical wavelength lambda_z (km)",
               title = "Normalized vertical-wavelength power (FFT proxy)",
                 color = :viridis, clims = (0.0, 1.0), colorbar_title = "Normalized power",
               left_margin = 6Plots.mm,
                 right_margin = 8Plots.mm,
                 bottom_margin = 8Plots.mm,
               top_margin = 6Plots.mm,
               legend = false)
    plot!(pb, times[tsel] ./ 60.0, lambda_dom ./ 1e3, color = :white, lw = 2,
           label = false)
    hline!(pb, [lambda_z_theory / 1e3], color = :red, ls = :dash, lw = 2,
            label = false)

    pfig6 = plot(pa, pb, layout = (1, 2), size = (1500, 600),
                 plot_title = "Linear dispersion replication (Figure 6 style)")
    savefig(pfig6, joinpath(OUTDIR, "dispersion_fig6_proxy.png"))

    # Save grayscale comparison outputs
    savefig(pa_bw, joinpath(OUTDIR, "dispersion_fig6_left_bw.png"))
    pfig6_bwleft = plot(pa_bw, pb, layout = (1, 2), size = (1500, 600),
                        plot_title = "Linear dispersion replication (left panel grayscale)")
    savefig(pfig6_bwleft, joinpath(OUTDIR, "dispersion_fig6_proxy_bwleft.png"))

    println("Saved metrics: " * csv_path)
    println("Saved figure: " * joinpath(OUTDIR, "dispersion_lambda_z_vs_time.png"))
    println("Saved figure: " * joinpath(OUTDIR, "dispersion_kz_spectrum.png"))
    println("Saved figure: " * joinpath(OUTDIR, "dispersion_fig6_proxy.png"))
    println("Saved figure: " * joinpath(OUTDIR, "dispersion_fig6_left_bw.png"))
    println("Saved figure: " * joinpath(OUTDIR, "dispersion_fig6_proxy_bwleft.png"))
end
