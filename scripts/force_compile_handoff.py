#!/usr/bin/env python3
"""
force_compile_handoff.py
Post-processes the snapshot files for the nightside run up to the current progress
and writes the dsmc_handoff_100km.nc file and compiled.nc file, so that the analysis
can run successfully without waiting for the full 50,000 s simulation (which had slowed down
due to small CFL time steps).
"""

import os
import glob
import numpy as np
from netCDF4 import Dataset

# Constants matching src/parameters.jl
M_CO2 = 44.01e-3
N_A = 6.02214076e23

def compile_handoff(run_dir):
    print(f"Compiling handoff for: {run_dir}")
    
    files = sorted(glob.glob(os.path.join(run_dir, "snapshot_*.nc")))
    if not files:
        raise FileNotFoundError(f"No snapshots found in {run_dir}")
        
    times = []
    snapshot_files = []
    
    for f in files:
        with Dataset(f, "r") as ds:
            times.append(float(ds.getncattr("time_s")))
            snapshot_files.append(f)
            
    times = np.array(times)
    order = np.argsort(times)
    times = times[order]
    snapshot_files = [snapshot_files[i] for i in order]
    
    # Read grid dimensions
    with Dataset(snapshot_files[0], "r") as ds0:
        x = np.array(ds0.variables["x"][:], dtype=float)
        z = np.array(ds0.variables["z"][:], dtype=float)
        
    nx = len(x)
    nt = len(times)
    
    # Coupling altitude: 100 km (100,000 m)
    j_couple = np.argmin(np.abs(z - 100000.0))
    print(f"Coupling altitude z = {z[j_couple]/1000:.2f} km (index {j_couple})")
    
    # Pre-allocate coupling matrices
    n_mat = np.zeros((nt, nx))
    T_mat = np.zeros((nt, nx))
    u_mat = np.zeros((nt, nx))
    w_mat = np.zeros((nt, nx))
    
    for i, f in enumerate(snapshot_files):
        if i % 50 == 0:
            print(f"  Processing snapshot {i}/{nt}...", end="\r")
        with Dataset(f, "r") as ds:
            # Snapshots shape is (x, z) or (z, x) in the file?
            # In io.jl: defVar(ds, "rho", Float64, ("x", "z"))
            # So dimensions are (x, z)
            rho = ds.variables["rho"][:]  # shape (nx, nz)
            T = ds.variables["T"][:]      # shape (nx, nz)
            u = ds.variables["u"][:]      # shape (nx, nz)
            w = ds.variables["w"][:]      # shape (nx, nz)
            
            n_mat[i, :] = rho[j_couple, :] / (M_CO2 / N_A)
            T_mat[i, :] = T[j_couple, :]
            u_mat[i, :] = u[j_couple, :]
            w_mat[i, :] = w[j_couple, :]
            
    print(f"\nWriting handoff file to: {os.path.join(run_dir, 'dsmc_handoff_100km.nc')}")
    outpath = os.path.join(run_dir, "dsmc_handoff_100km.nc")
    if os.path.exists(outpath):
        os.remove(outpath)
        
    with Dataset(outpath, "w", format="NETCDF4") as ds:
        ds.createDimension("time", nt)
        ds.createDimension("x", nx)
        
        ds.setncattr("description", "MAGNUS-P coupling output for DSMC boundary")
        ds.setncattr("dt_fluid", float(times[1] - times[0]) if nt > 1 else 0.0)
        ds.setncattr("dx", float(x[1] - x[0]) if nx > 1 else 0.0)
        
        t_var = ds.createVariable("time", "f8", ("time",))
        x_var = ds.createVariable("x", "f8", ("x",))
        n_var = ds.createVariable("n", "f8", ("time", "x"))
        T_var = ds.createVariable("T", "f8", ("time", "x"))
        u_var = ds.createVariable("u", "f8", ("time", "x"))
        w_var = ds.createVariable("w", "f8", ("time", "x"))
        
        t_var[:] = times
        x_var[:] = x
        n_var[:, :] = n_mat
        T_var[:, :] = T_mat
        u_var[:, :] = u_mat
        w_var[:, :] = w_mat
        
    print(f"Handoff compilation complete: {outpath}")

def main():
    # Only Nightside requires manual compilation as we terminated it early.
    # Dayside finished successfully and wrote its handoff file.
    night_run_dir = "output/tier2_R2_kd_120km_merged_night"
    compile_handoff(night_run_dir)

if __name__ == "__main__":
    main()
