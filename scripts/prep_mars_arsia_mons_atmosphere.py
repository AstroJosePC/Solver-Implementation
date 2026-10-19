#!/usr/bin/env python3
"""
prep_mars_arsia_mons_atmosphere.py

Builds a solver-ready background atmosphere for the Arsia Mons Elongated Cloud
(AMEC) formation scenario from a single Mars Climate Database (MCD) vertical
column.

Scenario (MY31 baseline, AMEC formation window):
    Latitude       -9.5  deg N
    Longitude     237.0  deg E
    Solar long.   226.0  deg (Ls)
    Local time      6.0  Martian hours
    Dust / EUV    Martian Year 31

The MCD column (Temperature, Pressure, Density, W-E wind vs height above the
local surface, 0-120 km) is fetched by mars_mcd_helper.fetch_web_form_data(),
which presets exactly these form fields. This script consumes the resulting
ASCII file (default: the already-cached fetch) and produces a profile on the
solver's uniform z-grid.

Physics handling (matches prep_mars_mcs_mcd_atmosphere.py conventions):
  * T(z) and U(z) are interpolated onto the solver grid with a monotone cubic
    (PCHIP) interpolant to avoid overshoot near the wind reversals.
  * Above the MCD top (120 km), T and U are held flat (isothermal + constant
    wind cap). This lets the solver grid extend past the MCD column (e.g. to
    ~156 km, keeping the 120 km coupling level below the sponge) without
    relying on unreliable PCHIP cubic extrapolation. The 120-156 km values are
    a modeling assumption, not MCD data, and lie inside/above the sponge layer.
  * A defensive static-stability guard clips any superadiabatic layer so
    dT/dz >= -g/cp (the solver reads T/rho/P directly and does not re-stabilize).
  * Pressure is re-integrated in DISCRETE hydrostatic balance on the solver grid
    from the MCD surface pressure, so the background is well-balanced. Density is
    then rho = P/(R*T). Raw MCD rho/P are used only as a cross-check.
  * nu = mu/rho with mu = 4.2e-7 * T^0.69, capped at 1e4 (src/parameters.jl).

Output CSV columns: z_m, rho, P, T, nu, u_zonal
  - The solver's tabulated_atmosphere() reads columns 1-4 (z_m, rho, P, T).
  - u_zonal is the MCD W-E wind on the grid, intended to be passed to
    run_simulation as `u_init` (a Vector{Float64} over grid.z_full); its
    near-surface value is the natural choice for the mountain steering wind
    U_bar in the bottom orographic forcing.
"""

import os
import argparse
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

# --- Mars constants (matching src/parameters.jl mars_config) ---
G_MARS = 3.72
R_GAS_MARS = 188.92210861167916   # J/(kg*K) = R_universal / M_CO2
GAMMA = 1.31                      # CO2-dominated ratio of specific heats
CP = GAMMA * R_GAS_MARS / (GAMMA - 1.0)

# Default cached MCD fetch for the Arsia Mons MY31 scenario (produced by
# mars_mcd_helper.fetch_web_form_data()).
DEFAULT_MCD_FILE = os.path.join(
    "external", "mars-mcd-helper", "tmp_mcd_browser_helper",
    "marsdata_web_3639f2e1c98b.txt",
)
OUTPUT_PATH = os.path.join("external", "mcd_profiles", "mars_arsia_mons_my31.csv")

# Map an MCD "Column 2 is ..." label substring to our internal variable key.
_LABEL_MAP = [
    ("Temperature", "T"),
    ("Pressure", "P"),
    ("Density", "rho"),
    ("W-E wind", "u"),
]


def parse_mcd_sections(path):
    """Parse a multi-section MCD ASCII file into {key: (z[m], value)} arrays.

    Sections are delimited by '#' comment headers; each section's variable is
    identified by its '### Column 2 is <name> (<unit>)' line.
    """
    sections = {}
    cur_key = None
    cur_z, cur_v = [], []

    def _flush():
        if cur_key is not None and cur_z:
            sections[cur_key] = (np.array(cur_z), np.array(cur_v))

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                if "Column 2 is" in s:
                    # Starting a new section: flush the previous one.
                    _flush()
                    cur_z, cur_v = [], []
                    cur_key = None
                    for sub, key in _LABEL_MAP:
                        if sub in s:
                            cur_key = key
                            break
                continue
            parts = s.split()
            if len(parts) >= 2 and cur_key is not None:
                cur_z.append(float(parts[0]))
                cur_v.append(float(parts[1]))
    _flush()
    return sections


def correct_surface_skin_temperature(z, T, P, rho, R=R_GAS_MARS, tol=0.03):
    """Replace an anomalous ground-level temperature with the air temperature.

    MCD reports the surface (skin) temperature at z=0, which at early local
    times is colder than the near-surface air and is inconsistent with the
    surface density/pressure. We detect this via the per-level implied gas
    constant R_eff = P/(rho*T): if the surface value deviates from the column
    median by more than `tol`, T(z=0) is reset to the air-consistent value
    P(0)/(rho(0)*R). All other levels are left untouched.
    """
    assert len(z) == len(P) == len(rho), "T/P/rho must share the MCD grid"
    R_eff = P / (rho * T)
    R_med = float(np.median(R_eff[1:]))
    T_out = T.copy()
    if abs(R_eff[0] - R_med) / R_med > tol:
        T_air = P[0] / (rho[0] * R)
        print(f"  surface skin-temp correction: T(z=0) {T[0]:.2f} K -> {T_air:.2f} K "
              f"(R_eff[0]={R_eff[0]:.1f} vs column median {R_med:.1f} J/kg/K)")
        T_out[0] = T_air
    return T_out


def enforce_static_stability(z_grid, T_grid):
    """Clip superadiabatic layers so dT/dz >= -g/cp, integrating upward."""
    min_dTdz = -G_MARS / CP
    T_adj = T_grid.copy()
    n_clipped = 0
    for i in range(len(z_grid) - 1):
        dz = z_grid[i + 1] - z_grid[i]
        dTdz = (T_adj[i + 1] - T_adj[i]) / dz
        if dTdz < min_dTdz:
            T_adj[i + 1] = T_adj[i] + min_dTdz * dz
            n_clipped += 1
    return T_adj, n_clipped


def integrate_hydrostatic(z_grid, T_grid, P_surf):
    """Discrete hydrostatic pressure on the solver grid."""
    P = np.zeros_like(z_grid)
    P[0] = P_surf
    for i in range(len(z_grid) - 1):
        dz = z_grid[i + 1] - z_grid[i]
        T_avg = 0.5 * (T_grid[i] + T_grid[i + 1])
        P[i + 1] = P[i] * np.exp(-G_MARS * dz / (R_GAS_MARS * T_avg))
    return P


def brunt_vaisala_sq(z_grid, T_grid):
    dTdz = np.gradient(T_grid, z_grid)
    return (G_MARS / T_grid) * (dTdz + G_MARS / CP)


def find_zero_crossings(z, u):
    """Altitudes where u changes sign (linearly interpolated)."""
    crossings = []
    for i in range(len(u) - 1):
        if u[i] == 0.0:
            crossings.append(z[i])
        elif u[i] * u[i + 1] < 0.0:
            t = -u[i] / (u[i + 1] - u[i])
            crossings.append(z[i] + t * (z[i + 1] - z[i]))
    return crossings


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mcd-file", default=DEFAULT_MCD_FILE,
                    help="MCD ASCII column file (default: cached Arsia Mons fetch)")
    ap.add_argument("--output", default=OUTPUT_PATH)
    ap.add_argument("--ztop", type=float, default=120000.0, help="top altitude [m]")
    ap.add_argument("--dz", type=float, default=500.0, help="grid spacing [m]")
    args = ap.parse_args()

    if not os.path.exists(args.mcd_file):
        raise FileNotFoundError(
            f"MCD file not found: {args.mcd_file}\n"
            "Re-fetch with mars_mcd_helper.fetch_web_form_data() "
            "(presets the Arsia Mons MY31 form fields)."
        )

    print(f"Reading MCD column: {args.mcd_file}")
    sec = parse_mcd_sections(args.mcd_file)
    for key in ("T", "P", "rho", "u"):
        if key not in sec:
            raise ValueError(f"MCD file missing required section '{key}'. Found: {list(sec)}")

    z_T, T_raw = sec["T"]
    z_P, P_raw = sec["P"]
    z_rho, rho_raw = sec["rho"]
    z_u, u_raw = sec["u"]
    print(f"  MCD levels: {len(z_T)}  span: {z_T.min():.0f}-{z_T.max():.0f} m")

    # Correct the anomalous ground-level (skin) temperature to air-consistent T.
    T_raw = correct_surface_skin_temperature(z_T, T_raw, P_raw, rho_raw)

    # Solver grid (matches src/grid.jl uniform spacing; 0 -> ztop inclusive).
    z_grid = np.arange(0.0, args.ztop + 0.5 * args.dz, args.dz)

    # Monotone-cubic interpolation onto the solver grid (no overshoot).
    T_grid = PchipInterpolator(z_T, T_raw, extrapolate=True)(z_grid)
    u_grid = PchipInterpolator(z_u, u_raw, extrapolate=True)(z_grid)

    # Above the MCD top, hold T and U flat (isothermal + constant-wind cap).
    # PCHIP cubic extrapolation above the data is unreliable; the flat cap keeps
    # the hydrostatic reconstruction well-balanced when the solver grid extends
    # past the MCD column (e.g. ztop=156 km).
    z_cap = float(z_T.max())
    above = z_grid > z_cap
    if np.any(above):
        T_grid[above] = T_raw[-1]
        u_grid[above] = u_raw[-1]
        print(f"  cap above {z_cap/1000:.0f} km: isothermal T={T_raw[-1]:.2f} K, "
              f"constant U={u_raw[-1]:.2f} m/s on {int(above.sum())} level(s)")

    # Defensive static-stability guard.
    T_grid, n_clipped = enforce_static_stability(z_grid, T_grid)
    if n_clipped:
        print(f"  stability guard clipped {n_clipped} superadiabatic layer(s)")

    # Well-balanced hydrostatic reconstruction from MCD surface pressure.
    P_surf = float(P_raw[0])
    P_grid = integrate_hydrostatic(z_grid, T_grid, P_surf)
    rho_grid = P_grid / (R_GAS_MARS * T_grid)
    mu = 4.2e-7 * T_grid ** 0.69
    nu_grid = np.minimum(mu / rho_grid, 1e4)

    # --- Diagnostics ---
    N2 = brunt_vaisala_sq(z_grid, T_grid)
    i0 = int(np.argmin(N2))
    z_top = float(z_grid[-1])
    # Cross-check hydrostatic rho against raw MCD rho only within the MCD span;
    # above z_cap the MCD column is extrapolated (cubic) and not a meaningful
    # reference. Interpolate MCD rho linearly in log to avoid PCHIP overshoot.
    in_mcd = (z_grid <= z_cap) & (z_grid >= z_rho.min())
    rho_mcd_on_grid = np.exp(np.interp(z_grid, z_rho, np.log(rho_raw)))
    rel_rho = np.abs(rho_grid[in_mcd] - rho_mcd_on_grid[in_mcd]) / rho_mcd_on_grid[in_mcd]
    crossings = find_zero_crossings(z_grid, u_grid)

    def at(z_q):
        return (np.interp(z_q, z_grid, T_grid),
                np.exp(np.interp(z_q, z_grid, np.log(rho_grid))),
                np.interp(z_q, z_grid, u_grid))

    print("\n--- Derived profile (solver grid) ---")
    print(f"  surface pressure          : {P_surf:.3f} Pa")
    print(f"  grid: {len(z_grid)} levels, 0 - {z_top/1000:.1f} km "
          f"(MCD data 0 - {z_cap/1000:.0f} km, capped above)")
    for zq in (0.0, 60000.0, 120000.0, z_top):
        Tq, rhoq, uq = at(zq)
        print(f"  z={zq/1000:6.1f} km : T={Tq:7.2f} K  rho={rhoq:.3e} kg/m3  U={uq:8.2f} m/s")
    print(f"  min N^2                    : {N2[i0]:.3e} s^-2 at z={z_grid[i0]/1000:.1f} km "
          f"({'STABLE' if N2[i0] >= 0 else 'UNSTABLE'})")
    print(f"  hydrostatic rho vs MCD rho : max rel diff = {rel_rho.max()*100:.2f}% "
          f"(within MCD span 0 - {z_cap/1000:.0f} km)")
    if crossings:
        print("  zonal-wind critical levels : "
              + ", ".join(f"{z/1000:.1f} km" for z in crossings))

    # --- Write ---
    df = pd.DataFrame({
        "z_m": z_grid,
        "rho": rho_grid,
        "P": P_grid,
        "T": T_grid,
        "nu": nu_grid,
        "u_zonal": u_grid,
    })
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    df.to_csv(args.output, index=False)
    print(f"\nWrote {len(z_grid)} levels to {args.output}")


if __name__ == "__main__":
    main()
