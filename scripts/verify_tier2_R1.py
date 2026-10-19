import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

def verify_handoff():
    ds = nc.Dataset('output/tier2_R1/dsmc_handoff_100km.nc', 'r')
    
    times = ds.variables['time'][:]
    w = ds.variables['w'][:]
    T = ds.variables['T'][:]
    n = ds.variables['n'][:]
    x = ds.variables['x'][:]
    
    # Extract the time series at the center of the domain (x = 100 km)
    mid_idx = len(x) // 2
    
    # NetCDF dims from Julia are ("time", "x") which load as ("x", "time") in Python netCDF4
    w_mid = w[mid_idx, :]
    T_mid = T[mid_idx, :]
    n_mid = n[mid_idx, :]
    
    fig, axs = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    
    axs[0].plot(times, w_mid, color='blue')
    axs[0].set_ylabel('Vertical Velocity (m/s)')
    axs[0].set_title('DSMC Handoff Data at z = 100 km (Tier 2 - R1)')
    axs[0].grid(True)
    
    axs[1].plot(times, T_mid, color='red')
    axs[1].set_ylabel('Temperature (K)')
    axs[1].grid(True)
    
    axs[2].plot(times, n_mid, color='green')
    axs[2].set_ylabel('Number Density (m^-3)')
    axs[2].set_xlabel('Time (s)')
    axs[2].grid(True)
    
    plt.tight_layout()
    plt.savefig('output/tier2_R1/dsmc_handoff_verification.png')
    print("Verification plot saved to output/tier2_R1/dsmc_handoff_verification.png")
    
    # Print some stats
    print(f"Max w perturbation: {np.max(np.abs(w_mid)):.4f} m/s")
    print(f"Mean T: {np.mean(T_mid):.2f} K, Max T perturbation: {np.max(T_mid) - np.mean(T_mid):.2f} K")
    print(f"Mean n: {np.mean(n_mid):.2e} m^-3, Fractional n pert: {(np.max(n_mid) - np.mean(n_mid))/np.mean(n_mid)*100:.2f}%")

if __name__ == '__main__':
    verify_handoff()
