#!/usr/bin/env python3
"""
prep_mars_mcs_mcd_atmosphere.py
Constructs realistic Martian background atmosphere profiles (dayside and nightside)
by blending low-altitude Mars Climate Sounder (MCS) temperature profiles with
high-altitude Mars Climate Database (MCD) profiles. It ensures strict discrete
numerical hydrostatic balance by integrating the temperature profile on the solver's grid.
"""

import os
import numpy as np
import pandas as pd
import scipy.interpolate as interp
# Note: `requests` and `bs4` are imported lazily inside fetch_mcd_profile() so the
# script can regenerate profiles from cached MCD files without the network libraries.

# Mars Constants (matching src/parameters.jl mars_config)
G_MARS = 3.72
R_GAS_MARS = 188.92210861167916  # J/(kg*K)
GAMMA = 1.31                      # ratio of specific heats for CO2-dominated atmosphere
CP = GAMMA * R_GAS_MARS / (GAMMA - 1.0)  # specific heat at constant pressure [J/(kg*K)]

def fetch_mcd_profile(var, ls, localtime, cache_path):
    """
    Fetches the profile for the given variable from MCD, or loads it from cache if network fails.
    """
    url = "http://www-mars.lmd.jussieu.fr/mcd_python/cgi-bin/mcdcgi.py"
    p = {
        "datekeyhtml": 1,
        "ls": ls,
        "localtime": localtime,
        "latitude": 0.0,
        "longitude": 0.0,
        "altitude": "0 120000",
        "zkey": 2,
        "isfixedlt": "off",
        "dust": 1,
        "hrkey": 1,
        "zonmean": "off",
        "var1": var,
        "var2": "none",
        "var3": "none",
        "var4": "none",
        "dpi": 80,
        "islog": "off",
        "colorm": "Blues",
        "proj": "cyl",
        "iswind": "off",
        "betatest": "on",  # Selects MCD v6.1 database
        "averaging": "off",
        "diumean": "off",
        "animation": "off"
    }
    
    raw_text = None
    try:
        import requests
        from bs4 import BeautifulSoup
        print(f"Querying MCD CGI for variable={var}, Ls={ls}, LT={localtime}...")
        r = requests.get(url, params=p, timeout=20)
        r.raise_for_status()
        
        soup = BeautifulSoup(r.text, "html.parser")
        a_tag = soup.find("a")
        if not a_tag:
            raise ValueError(f"No download link found in CGI response. Response:\n{r.text[:500]}")
            
        href = a_tag["href"]
        download_url = "http://www-mars.lmd.jussieu.fr/mcd_python/" + href.replace("../", "")
        print(f"Downloading data from {download_url}...")
        r_dl = requests.get(download_url, timeout=20)
        r_dl.raise_for_status()
        raw_text = r_dl.text
        
        # Save to cache
        os.makedirs(os.path.dirname(cache_path), exist_ok=True)
        with open(cache_path, "w", encoding="utf-8") as f:
            f.write(raw_text)
        print(f"Saved profile to cache: {cache_path}")
        
    except Exception as e:
        print(f"Network request or parsing failed: {e}")
        if os.path.exists(cache_path):
            print(f"Loading from cache fallback: {cache_path}")
            with open(cache_path, "r", encoding="utf-8") as f:
                raw_text = f.read()
        else:
            raise RuntimeError(f"No cached file found at {cache_path} and network request failed.")
            
    # Parse the raw text
    altitudes = []
    values = []
    for line in raw_text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 2:
            altitudes.append(float(parts[0]))
            values.append(float(parts[1]))
            
    return np.array(altitudes), np.array(values)

def read_digitized_mcs(filepath):
    # Column 0: value, Column 1: altitude (km)
    df = pd.read_csv(filepath, header=None)
    val = df.iloc[:, 0].values
    z_km = df.iloc[:, 1].values
    
    # Sort by altitude ascending
    idx = np.argsort(z_km)
    return z_km[idx] * 1000.0, val[idx]

def blend_profiles(z_grid, z_mcs, T_mcs, z_mcd, T_mcd):
    import scipy.ndimage as ndimage
    
    # 1. Interpolate MCS profile onto our uniform grid
    f_mcs_interp = interp.interp1d(
        z_mcs, T_mcs, kind='linear', bounds_error=False,
        fill_value=(T_mcs[0], T_mcs[-1])
    )
    T_mcs_grid = f_mcs_interp(z_grid)
    
    # Below the lowest digitized MCS altitude, extend the MCS profile downward
    # using its own near-surface lapse rate (from the two lowest MCS samples).
    #
    # The previous behaviour grafted the MCD surface temperature (~232 K) onto z=0
    # and ramped down to the MCS value (~205 K) over <1 km. That produced a lapse
    # rate of ~37-41 K/km (the dry adiabatic value is ~4.66 K/km), i.e. a
    # superadiabatic, convectively unstable surface layer (N^2 < 0). Because the
    # solver reads T/rho/P directly without re-stabilizing, this seeded spurious
    # convection. Extending the MCS lapse keeps the surface near ~205 K (matching
    # Srivastava et al. 2022, Fig. 11) and the layer statically stable.
    z_min_mcs = z_mcs[0]
    idx_below = z_grid < z_min_mcs
    if np.any(idx_below):
        dTdz_surf = (T_mcs[1] - T_mcs[0]) / (z_mcs[1] - z_mcs[0])  # MCS near-surface lapse
        T_mcs_grid[idx_below] = T_mcs[0] + dTdz_surf * (z_grid[idx_below] - z_min_mcs)
        
    # 2. Smooth the MCD profile to represent a background in equilibrium without thermal tides / wavy artifacts
    T_mcd_smoothed = ndimage.gaussian_filter1d(T_mcd, sigma=3.0)
        
    # 3. Interpolate smoothed MCD profile onto our uniform grid
    f_mcd_interp = interp.interp1d(
        z_mcd, T_mcd_smoothed, kind='cubic', bounds_error=False,
        fill_value="extrapolate"
    )
    T_mcd_grid = f_mcd_interp(z_grid)
    
    # 4. Blend them over 25-35 km using a cosine transition
    # (Since digitized temps.csv only goes up to ~39.4 km)
    T_blended = np.zeros_like(z_grid)
    z_blend_start = 25000.0
    z_blend_end = 35000.0
    
    for i, z in enumerate(z_grid):
        if z < z_blend_start:
            T_blended[i] = T_mcs_grid[i]
        elif z > z_blend_end:
            T_blended[i] = T_mcd_grid[i]
        else:
            w = 0.5 - 0.5 * np.cos(np.pi * (z - z_blend_start) / (z_blend_end - z_blend_start))
            T_blended[i] = (1.0 - w) * T_mcs_grid[i] + w * T_mcd_grid[i]
            
    return T_blended

def enforce_static_stability(z_grid, T_grid):
    """
    Defensive guard: clip any superadiabatic layer so the temperature lapse rate
    never exceeds the dry adiabatic value (dT/dz >= -g/cp), integrating upward
    from the surface. Returns the adjusted profile and the count of clipped layers.
    With a correct blend this should clip zero layers.
    """
    min_dTdz = -G_MARS / CP  # most negative (steepest) statically stable gradient
    T_adj = T_grid.copy()
    n_clipped = 0
    for i in range(len(z_grid) - 1):
        dz = z_grid[i + 1] - z_grid[i]
        dTdz = (T_adj[i + 1] - T_adj[i]) / dz
        if dTdz < min_dTdz:
            T_adj[i + 1] = T_adj[i] + min_dTdz * dz
            n_clipped += 1
    return T_adj, n_clipped

def min_brunt_vaisala(z_grid, T_grid):
    """Return (min N^2 [s^-2], altitude of minimum [m]) for the temperature profile."""
    dTdz = np.gradient(T_grid, z_grid)
    N2 = (G_MARS / T_grid) * (dTdz + G_MARS / CP)
    i = int(np.argmin(N2))
    return N2[i], z_grid[i]

def integrate_hydrostatic(z_grid, T_grid, P_surf):
    P_grid = np.zeros_like(z_grid)
    P_grid[0] = P_surf
    for i in range(len(z_grid) - 1):
        dz = z_grid[i+1] - z_grid[i]
        T_avg = 0.5 * (T_grid[i] + T_grid[i+1])
        P_grid[i+1] = P_grid[i] * np.exp(-G_MARS * dz / (R_GAS_MARS * T_avg))
    return P_grid

def process_case(case_name, ls, localtime, z_grid, z_mcs, T_mcs):
    cache_dir = "scratch"
    t_cache = os.path.join(cache_dir, f"mcd_cache_{case_name}_t.txt")
    p_cache = os.path.join(cache_dir, f"mcd_cache_{case_name}_p.txt")
    
    # Fetch MCD profiles
    z_mcd_t, T_mcd = fetch_mcd_profile("t", ls, localtime, t_cache)
    z_mcd_p, P_mcd = fetch_mcd_profile("p", ls, localtime, p_cache)
    
    # Surface pressure from MCD at z=0 (first element of profile)
    P_surface_mcd = P_mcd[0]
    print(f"[{case_name}] Surface Pressure from MCD: {P_surface_mcd:.3f} Pa")
    
    # Blend MCS and MCD temperatures
    T_blended = blend_profiles(z_grid, z_mcs, T_mcs, z_mcd_t, T_mcd)

    # Defensive static-stability guard against superadiabatic layers
    T_blended, n_clipped = enforce_static_stability(z_grid, T_blended)
    if n_clipped > 0:
        print(f"[{case_name}] WARNING: stability guard clipped {n_clipped} superadiabatic layer(s).")

    # Report minimum static stability of the final temperature profile
    N2_min, z_N2min = min_brunt_vaisala(z_grid, T_blended)
    status = "STABLE" if N2_min >= 0.0 else "UNSTABLE"
    print(f"[{case_name}] Min N^2 = {N2_min:.3e} s^-2 at z = {z_N2min/1000:.1f} km ({status}); "
          f"T_surface = {T_blended[0]:.2f} K")

    # Integrate hydrostatic pressure
    P_blended = integrate_hydrostatic(z_grid, T_blended, P_surface_mcd)
    
    # Compute density rho = P / (R * T)
    rho_blended = P_blended / (R_GAS_MARS * T_blended)
    
    # Compute kinematic viscosity nu = mu / rho where mu = 4.2e-7 * T^0.69
    mu = 4.2e-7 * T_blended**0.69
    nu_blended = mu / rho_blended
    # Cap kinematic viscosity to 1e4 to match src/parameters.jl capping
    nu_blended = np.minimum(nu_blended, 1e4)
    
    # Create DataFrame
    df = pd.DataFrame({
        'z_m': z_grid,
        'rho': rho_blended,
        'P': P_blended,
        'T': T_blended,
        'nu': nu_blended
    })
    
    # Save to CSV
    output_path = f"external/digitized_figure_11/mars_mcs_mcd_merged_{case_name}.csv"
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    df.to_csv(output_path, index=False)
    print(f"[{case_name}] Successfully generated and saved to {output_path}")

def main():
    # Grid spacing: 0 to 120 km with 500m resolution (matching solver grid)
    z_grid = np.arange(0, 120001.0, 500.0)
    
    # Read MCS digitized profile
    mcs_filepath = "external/digitized_figure_11/temps.csv"
    print(f"Loading low-altitude MCS digitized profile from {mcs_filepath}...")
    z_mcs, T_mcs = read_digitized_mcs(mcs_filepath)
    
    # Case 1: Dayside (DD2-like)
    # Ls = 327, LT = 12.0
    process_case("day", 327.0, 12.0, z_grid, z_mcs, T_mcs)
    
    # Case 2: Nightside (DD6-like)
    # Ls = 195, LT = 0.0
    process_case("night", 195.0, 0.0, z_grid, z_mcs, T_mcs)

if __name__ == "__main__":
    main()
