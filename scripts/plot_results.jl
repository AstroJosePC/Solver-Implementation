# =============================================================================
# Visualization helpers for MAGNUS-P output
# Requires: Plots.jl (install separately: ] add Plots)
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using NCDatasets
using Printf

"""
    plot_snapshot(filepath; variable=:w)

Load and display a NetCDF snapshot. Requires Plots.jl to be loaded separately.
This script prints the data summary; actual plotting requires:
  using Plots; heatmap(x/1e3, z/1e3, data')
"""
function load_snapshot(filepath::String)
    ds = NCDataset(filepath, "r")
    x = ds["x"][:]
    z = ds["z"][:]
    rho = ds["rho"][:, :]
    u = ds["u"][:, :]
    w = ds["w"][:, :]
    P = ds["P"][:, :]
    T = ds["T"][:, :]
    t = ds.attrib["time_s"]
    close(ds)
    return Dict(:x => x, :z => z, :rho => rho, :u => u, :w => w, :P => P, :T => T, :time => t)
end

"""
    summarize_snapshot(filepath)

Print a text summary of a snapshot file.
"""
function summarize_snapshot(filepath::String)
    data = load_snapshot(filepath)
    t = data[:time]
    @printf("Snapshot at t = %.3f s\n", t)
    @printf("  Grid: %d x %d\n", length(data[:x]), length(data[:z]))
    @printf("  x range: [%.1f, %.1f] km\n", data[:x][1]/1e3, data[:x][end]/1e3)
    @printf("  z range: [%.1f, %.1f] km\n", data[:z][1]/1e3, data[:z][end]/1e3)
    for var in [:rho, :u, :w, :P, :T]
        vals = data[var]
        @printf("  %3s: min=%12.5e, max=%12.5e, mean=%12.5e\n",
                var, minimum(vals), maximum(vals), sum(vals)/length(vals))
    end
end

# If run as script, summarize all snapshots in a directory
if abspath(PROGRAM_FILE) == @__FILE__
    outdir = length(ARGS) > 0 ? ARGS[1] : "output/dispersion_test"
    files = filter(f -> endswith(f, ".nc"), readdir(outdir, join=true))
    sort!(files)
    for f in files
        summarize_snapshot(f)
        println()
    end
end
