#!/usr/bin/env python3
"""
compile_mountain_run.py

Consolidate all per-snapshot NetCDF files from one Mars mountain wave run
into a single structured NetCDF4 file.

Output layout
-------------
Dimensions  : time, z, x
Coordinates : time(time)  [s]
              z(z)        [km]
              x(x)        [km, absolute]
              x_rel(x)    [km, relative to mountain center]
Fields      : w(time, z, x)    [m/s]   vertical velocity
              u(time, z, x)    [m/s]   horizontal velocity
              rho(time, z, x)  [kg/m3] density
              (any additional 2-D fields present in the snapshots)
Attributes  : run_name, u_bar_mps, sigma_km, n_snapshots,
              source_dir, x_center_m

Chunking    : (chunk_t, nz, nx)
              One full z-x plane per time-chunk -> O(1) snapshot reads,
              O(nt / chunk_t) reads for a full point time series.
Compression : zlib level 4 (float32 fields, float64 time coord)

Usage
-----
    python scripts/compile_mountain_run.py output/sweep_mars_mountain/run_U5.0_sig5.0

    # compile all runs in a sweep directory (PowerShell)
    Get-ChildItem output/sweep_mars_mountain -Directory | ForEach-Object {
        python scripts/compile_mountain_run.py $_.FullName
    }

    # lower memory use on a small machine
    python scripts/compile_mountain_run.py ... --chunk-t 50 --batch-t 50
"""

from __future__ import annotations

import argparse
import glob
import os
import re

import numpy as np
from netCDF4 import Dataset


RUN_DIR_PATTERN = re.compile(r"run_U(?P<u>[-+0-9.]+)_sig(?P<sig>[-+0-9.]+)$")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Compile Mars mountain wave run snapshots into one NetCDF4 file"
    )
    p.add_argument(
        "run_dir",
        help="Path to run directory, e.g. output/sweep_mars_mountain/run_U5.0_sig5.0",
    )
    p.add_argument(
        "--outdir",
        default="output/sweep_mars_mountain_compiled",
        help="Directory for the compiled file (default: output/sweep_mars_mountain_compiled)",
    )
    p.add_argument(
        "--chunk-t",
        type=int,
        default=100,
        help="Time steps per NetCDF chunk (default: 100)",
    )
    p.add_argument(
        "--batch-t",
        type=int,
        default=200,
        help="Time steps read into RAM at once; reduce to lower peak memory (default: 200)",
    )
    p.add_argument(
        "--complevel",
        type=int,
        default=4,
        help="zlib compression level 1-9 (default: 4)",
    )
    return p.parse_args()


def discover_snapshots(run_dir: str) -> tuple[list[str], np.ndarray]:
    files = sorted(glob.glob(os.path.join(run_dir, "*.nc")))
    if not files:
        raise FileNotFoundError(f"No .nc snapshot files found in {run_dir!r}")
    times = np.empty(len(files), dtype=np.float64)
    for i, fpath in enumerate(files):
        with Dataset(fpath, "r") as ds:
            times[i] = float(ds.getncattr("time_s"))
    order = np.argsort(times)
    return [files[i] for i in order], times[order]


def get_2d_field_names(sample_file: str) -> list[str]:
    """Return names of all 2-D (nz x nx) variables found in a snapshot file."""
    with Dataset(sample_file, "r") as ds:
        return [name for name, var in ds.variables.items() if var.ndim == 2]


def compile_run(run_dir: str, outdir: str, chunk_t: int, batch_t: int, complevel: int) -> str:
    run_name = os.path.basename(os.path.normpath(run_dir))
    m = RUN_DIR_PATTERN.match(run_name)
    if m is None:
        raise ValueError(
            f"Cannot parse U/sigma from run directory name: {run_name!r}\n"
            "Expected pattern: run_U<float>_sig<float>"
        )
    u_bar    = float(m.group("u"))
    sigma_km = float(m.group("sig"))

    print(f"Run     : {run_name}  (U={u_bar} m/s, sigma={sigma_km} km)")

    files, times = discover_snapshots(run_dir)
    nt = len(files)
    print(f"Snapshots: {nt}  |  t = {times[0]:.0f} s ... {times[-1]:.0f} s")

    # Read grid and field inventory from first snapshot
    with Dataset(files[0], "r") as ds0:
        x_m         = np.asarray(ds0.variables["x"][:], dtype=np.float64)
        z_m         = np.asarray(ds0.variables["z"][:], dtype=np.float64)
        field_names = get_2d_field_names(files[0])

    # Remove coordinate variables from the field list
    coord_names = {"x", "z"}
    field_names = [n for n in field_names if n not in coord_names]

    nz, nx = len(z_m), len(x_m)
    x_center_m = 0.5 * (x_m.min() + x_m.max())
    x_rel_km   = (x_m - x_center_m) / 1000.0
    x_km       = x_m / 1000.0
    z_km       = z_m / 1000.0

    ct  = min(chunk_t, nt)
    bt  = min(batch_t, nt)
    enc = dict(zlib=True, complevel=complevel)

    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, f"{run_name}.nc")

    print(f"Output  : {outpath}")
    print(f"Grid    : nt={nt}, nz={nz}, nx={nx}")
    print(f"Fields  : {field_names}")
    print(f"Chunks  : ({ct}, {nz}, {nx})   batch_t={bt}   complevel={complevel}")

    with Dataset(outpath, "w", format="NETCDF4") as out:
        # Global attributes
        out.run_name    = run_name
        out.u_bar_mps   = u_bar
        out.sigma_km    = sigma_km
        out.n_snapshots = nt
        out.source_dir  = os.path.abspath(run_dir)
        out.x_center_m  = float(x_center_m)

        # Dimensions
        out.createDimension("time", nt)
        out.createDimension("z",    nz)
        out.createDimension("x",    nx)

        # Coordinate variables
        vt            = out.createVariable("time",  "f8", ("time",), **enc)
        vt.units      = "s"
        vt.long_name  = "simulation time"
        vt[:]         = times

        vz            = out.createVariable("z",     "f4", ("z",), **enc)
        vz.units      = "km"
        vz.long_name  = "altitude"
        vz[:]         = z_km.astype(np.float32)

        vx            = out.createVariable("x",     "f4", ("x",), **enc)
        vx.units      = "km"
        vx.long_name  = "x position (absolute)"
        vx[:]         = x_km.astype(np.float32)

        vxr           = out.createVariable("x_rel", "f4", ("x",), **enc)
        vxr.units     = "km"
        vxr.long_name = "x relative to mountain center"
        vxr[:]        = x_rel_km.astype(np.float32)

        # Field variables — pre-create with chunking
        field_handles: dict[str, object] = {}
        for name in field_names:
            v = out.createVariable(
                name, "f4", ("time", "z", "x"),
                chunksizes=(ct, nz, nx),
                **enc,
            )
            field_handles[name] = v

        # Fill fields in time batches to limit peak RAM
        n_batches = (nt + bt - 1) // bt
        for b in range(n_batches):
            i0   = b * bt
            i1   = min(i0 + bt, nt)
            blen = i1 - i0
            pct  = 100.0 * i1 / nt
            print(f"  {pct:5.1f}%  [{i1}/{nt}]  t-indices {i0}-{i1 - 1} ...", end="\r", flush=True)

            buf: dict[str, np.ndarray] = {
                name: np.empty((blen, nz, nx), dtype=np.float32)
                for name in field_names
            }
            for j, fpath in enumerate(files[i0:i1]):
                with Dataset(fpath, "r") as ds:
                    for name in field_names:
                        buf[name][j] = np.asarray(ds.variables[name][:], dtype=np.float32)

            for name in field_names:
                field_handles[name][i0:i1] = buf[name]

    size_mb = os.path.getsize(outpath) / 1024**2
    print(f"\nDone. {outpath}  ({size_mb:.1f} MB)")
    return outpath


def main() -> None:
    args = parse_args()
    compile_run(
        run_dir=args.run_dir,
        outdir=args.outdir,
        chunk_t=args.chunk_t,
        batch_t=args.batch_t,
        complevel=args.complevel,
    )


if __name__ == "__main__":
    main()
