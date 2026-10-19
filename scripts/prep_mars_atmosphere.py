#!/usr/bin/env python3
"""
prep_mars_atmosphere.py
Enforces strict numerical hydrostatic balance by integrating the temperature profile.
Output is a CSV file with columns: z_m, rho, P, T, nu
"""

import numpy as np
import pandas as pd
import scipy.interpolate as interp

# Mars Constants (matching src/parameters.jl)
G_MARS = 3.72
R_GAS_MARS = 188.92  # J/(kg*K)
P_SURFACE_MARS = 636.0

def read_digitized(filepath):
    # Column 0: value, Column 1: altitude (km)
    df = pd.read_csv(filepath, header=None)
    val = df.iloc[:, 0].values
    z_km = df.iloc[:, 1].values
    
    # Sort by altitude ascending
    idx = np.argsort(z_km)
    return z_km[idx] * 1000.0, val[idx]

def main():
    # Desired uniform grid: 0 to 220 km with 500m spacing
    z_uniform = np.arange(0, 220001.0, 500.0)
    
    # Read ONLY Temperature
    z_T, T_raw = read_digitized("external/digitized_figure_11/temps.csv")
    
    # Interpolate Temperature onto uniform grid
    # Clamp the values at boundaries (isothermal thermosphere approximation above max z)
    f_T = interp.interp1d(z_T, T_raw, kind='linear', bounds_error=False, fill_value=(T_raw[0], T_raw[-1]))
    T_uni = f_T(z_uniform)
    T_uni = np.maximum(T_uni, 10.0) # Safety floor
    
    # Integrate Hydrostatic Equation: dP/dz = -rho * g = -(P/(RT)) * g
    # dP/P = - (g / (R*T)) dz
    # Discrete integration: P(i+1) = P(i) * exp(-g * dz / (R * T_avg))
    P_uni = np.zeros_like(z_uniform)
    P_uni[0] = P_SURFACE_MARS
    
    for i in range(len(z_uniform) - 1):
        dz = z_uniform[i+1] - z_uniform[i]
        T_avg = 0.5 * (T_uni[i] + T_uni[i+1])
        P_uni[i+1] = P_uni[i] * np.exp(-G_MARS * dz / (R_GAS_MARS * T_avg))
        
    # Compute Density via Ideal Gas Law: rho = P / (R*T)
    rho_uni = P_uni / (R_GAS_MARS * T_uni)
    
    # Calculate kinematic viscosity nu = mu / rho
    # Using formula: mu = 4.2e-7 * T^0.69 (from Appendix A)
    nu_uni = (4.2e-7 * T_uni**0.69) / rho_uni
    
    out_df = pd.DataFrame({
        'z_m': z_uniform,
        'rho': rho_uni,
        'P': P_uni,
        'T': T_uni,
        'nu': nu_uni
    })
    
    out_df.to_csv("external/digitized_figure_11/mars_mcs_uniform.csv", index=False)
    print("Successfully generated external/digitized_figure_11/mars_mcs_uniform.csv with numerical hydrostatic balance.")

if __name__ == "__main__":
    main()