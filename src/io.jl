# =============================================================================
# NetCDF I/O for simulation snapshots and DSMC coupling output
# =============================================================================

using NCDatasets

"""
    write_snapshot_nc(Q, grid, config, t, step, outdir)

Write a 2D snapshot of primitive variables to a NetCDF file.

File: outdir/snapshot_XXXXXXX.nc
Variables: rho, u, w, P, T on the physical grid.
"""
function write_snapshot_nc(Q::Array{Float64,3}, grid::Grid2D, config::SimConfig,
                           t::Float64, step::Int, outdir::String)
    # Ensure filename is clean and absolute
    fname = abspath(joinpath(outdir, @sprintf("snapshot_%07d.nc", step)))

    # NetCDF/HDF5 libraries are generally not thread-safe. Use global lock for I/O.
    # Serializing snapshots avoids metadata corruption and "No such file or directory" race conditions.
    lock(IO_LOCK) do
        nx = grid.nx
        nz = grid.nz
        ng = grid.ng

        # Extract primitive variables on physical grid
        rho = zeros(nx, nz)
        u   = zeros(nx, nz)
        w   = zeros(nx, nz)
        P   = zeros(nx, nz)
        T   = zeros(nx, nz)

        for j in 1:nz, i in 1:nx
            ig = i + ng
            jg = j + ng
            Qi = SVector(Q[1,ig,jg], Q[2,ig,jg], Q[3,ig,jg], Q[4,ig,jg])
            rho_v, u_v, w_v, P_v = primitive_from_conservative(Qi, config.gamma)
            rho[i, j] = rho_v
            u[i, j]   = u_v
            w[i, j]   = w_v
            P[i, j]   = P_v
            T[i, j]   = P_v / (rho_v * config.R_gas)
        end

        # Simple retry logic to handle transient flakiness on external drives/concurrent VFS metadata issues
        for attempt in 1:3
            try
                NCDataset(fname, "c") do ds
                    # Dimensions
                    defDim(ds, "x", nx)
                    defDim(ds, "z", nz)

                    # Coordinate variables
                    x_var = defVar(ds, "x", Float64, ("x",),
                                   attrib = Dict("units" => "m", "long_name" => "horizontal coordinate"))
                    z_var = defVar(ds, "z", Float64, ("z",),
                                   attrib = Dict("units" => "m", "long_name" => "altitude"))
                    x_var[:] = grid.x
                    z_var[:] = grid.z

                    # Time
                    ds.attrib["time_s"] = t
                    ds.attrib["step"] = step
                    ds.attrib["gamma"] = config.gamma
                    ds.attrib["g"] = config.g

                    # Data variables
                    rho_var = defVar(ds, "rho", Float64, ("x", "z"),
                                     attrib = Dict("units" => "kg/m3", "long_name" => "density"))
                    u_var = defVar(ds, "u", Float64, ("x", "z"),
                                   attrib = Dict("units" => "m/s", "long_name" => "horizontal velocity"))
                    w_var = defVar(ds, "w", Float64, ("x", "z"),
                                   attrib = Dict("units" => "m/s", "long_name" => "vertical velocity"))
                    P_var = defVar(ds, "P", Float64, ("x", "z"),
                                   attrib = Dict("units" => "Pa", "long_name" => "pressure"))
                    T_var = defVar(ds, "T", Float64, ("x", "z"),
                                   attrib = Dict("units" => "K", "long_name" => "temperature"))

                    rho_var[:, :] = rho
                    u_var[:, :]   = u
                    w_var[:, :]   = w
                    P_var[:, :]   = P
                    T_var[:, :]   = T
                end
                break # Success, exit retry loop
            catch e
                if attempt < 3
                    @warn "Snapshot retry $(attempt)/3: Failed to write $(fname). Retrying..."
                    sleep(0.5) # Give the drive/VFS a breather
                else
                    rethrow(e) # Give up after 3 attempts
                end
            end
        end
    end

    return fname
end

"""
    write_coupling_file_nc(coupling_data, config, outpath)

Write the 1D time series at the coupling altitude for DSMC ingestion.
Per planA Phase 2 format.

coupling_data: Dict with keys:
  :time => Vector{Float64}
  :x    => Vector{Float64}
  :n    => Matrix{Float64} (num_steps x num_x_cells) number density
  :T    => Matrix{Float64} temperature
  :u    => Matrix{Float64} horizontal velocity
  :w    => Matrix{Float64} vertical velocity
"""
function write_coupling_file_nc(coupling_data::Dict, config::SimConfig,
                                outpath::String)
    times = coupling_data[:time]
    x = coupling_data[:x]
    nt = length(times)
    nx = length(x)

    lock(IO_LOCK) do
        NCDataset(outpath, "c") do ds
            defDim(ds, "time", nt)
            defDim(ds, "x", nx)

            ds.attrib["description"] = "MAGNUS-P coupling output for DSMC boundary"
            ds.attrib["dt_fluid"] = length(times) > 1 ? times[2] - times[1] : 0.0
            ds.attrib["dx"] = length(x) > 1 ? x[2] - x[1] : 0.0

            t_var = defVar(ds, "time", Float64, ("time",),
                           attrib = Dict("units" => "s"))
            x_var = defVar(ds, "x", Float64, ("x",),
                           attrib = Dict("units" => "m"))
            n_var = defVar(ds, "n", Float64, ("time", "x"),
                           attrib = Dict("units" => "1/m3", "long_name" => "number density"))
            T_var = defVar(ds, "T", Float64, ("time", "x"),
                           attrib = Dict("units" => "K", "long_name" => "temperature"))
            u_var = defVar(ds, "u", Float64, ("time", "x"),
                           attrib = Dict("units" => "m/s", "long_name" => "horizontal velocity"))
            w_var = defVar(ds, "w", Float64, ("time", "x"),
                           attrib = Dict("units" => "m/s", "long_name" => "vertical velocity"))

            t_var[:] = times
            x_var[:] = x
            n_var[:, :] = coupling_data[:n]
            T_var[:, :] = coupling_data[:T]
            u_var[:, :] = coupling_data[:u]
            w_var[:, :] = coupling_data[:w]
        end
    end

    return outpath
end
