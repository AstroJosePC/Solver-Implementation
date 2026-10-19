import argparse
import glob
import os
import numpy as np
from netCDF4 import Dataset

def compile_run(run_dir: str, outdir: str):
    run_name = os.path.basename(os.path.normpath(run_dir))
    files = sorted(glob.glob(os.path.join(run_dir, "snapshot_*.nc")))
    if not files:
        print("No snapshots found")
        return

    times = []
    for f in files:
        with Dataset(f, "r") as ds:
            times.append(float(ds.getncattr("time_s")))
    times = np.array(times)
    order = np.argsort(times)
    files = [files[i] for i in order]
    times = times[order]

    nt = len(files)
    with Dataset(files[0], "r") as ds0:
        x_m = np.asarray(ds0.variables["x"][:], dtype=np.float64)
        z_m = np.asarray(ds0.variables["z"][:], dtype=np.float64)
        field_names = [name for name, var in ds0.variables.items() if var.ndim == 2]

    coord_names = {"x", "z"}
    field_names = [n for n in field_names if n not in coord_names]

    nz, nx = len(z_m), len(x_m)
    x_km = x_m / 1000.0
    z_km = z_m / 1000.0

    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, f"{run_name}.nc")

    with Dataset(outpath, "w", format="NETCDF4") as out:
        out.createDimension("time", nt)
        out.createDimension("z", nz)
        out.createDimension("x", nx)

        vt = out.createVariable("time", "f8", ("time",))
        vt[:] = times
        vz = out.createVariable("z", "f4", ("z",))
        vz[:] = z_km.astype(np.float32)
        vx = out.createVariable("x", "f4", ("x",))
        vx[:] = x_km.astype(np.float32)

        field_handles = {}
        for name in field_names:
            field_handles[name] = out.createVariable(name, "f4", ("time", "z", "x"))

        for i, f in enumerate(files):
            print(f"Processing {i}/{nt}", end="\r")
            with Dataset(f, "r") as ds:
                for name in field_names:
                    field_handles[name][i] = ds.variables[name][:]

    print(f"\nDone: {outpath}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir")
    parser.add_argument("--outdir", default="output/compiled")
    args = parser.parse_args()
    compile_run(args.run_dir, args.outdir)
