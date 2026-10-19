import netCDF4 as nc
import numpy as np

def check():
    outdir_off = "output/eckermann_fig14/cooling_off"
    outdir_on = "output/eckermann_fig14/cooling_on"
    
    import glob, os
    files_off = sorted(glob.glob(os.path.join(outdir_off, "snapshot_*.nc")))
    files_on = sorted(glob.glob(os.path.join(outdir_on, "snapshot_*.nc")))
    
    print(f"Num files OFF: {len(files_off)}, ON: {len(files_on)}")
    
    ds = nc.Dataset(files_off[-1])
    w = ds.variables['w'][:]
    z = ds.variables['z'][:] / 1e3
    w_mid = w[w.shape[0]//2, :]
    print("Cooling OFF max w:", np.max(np.abs(w_mid)))
    
    ds_on = nc.Dataset(files_on[-1])
    w_on = ds_on.variables['w'][:]
    w_mid_on = w_on[w_on.shape[0]//2, :]
    print("Cooling ON max w:", np.max(np.abs(w_mid_on)))

if __name__ == "__main__":
    check()
