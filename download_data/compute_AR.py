from netCDF4 import Dataset
import numpy as np
import os
import re
import sys

g0 = 9.81

def pleaseRun(cmd):
    print(">> %s" % cmd)
    os.system(cmd)


def computeAR(in_filename, out_filename, varnames=dict(u="UGRD", v="VGRD", q="q", lat="lat", lon="lon", lev="lev")):


    with Dataset(in_filename, "r") as ds:

        lat = ds.variables[varnames["lat"]][:]
        lon = ds.variables[varnames["lon"]][:]
        lev = ds.variables[varnames["lev"]][:]

        lev_crop_idx = (lev >= 200) & (lev <= 1000)
        lev = lev[lev_crop_idx]

        q = ds.variables[varnames["q"]][:, lev_crop_idx, :, :]
        U = ds.variables[varnames["u"]][:, lev_crop_idx, :, :]
        V = ds.variables[varnames["v"]][:, lev_crop_idx, :, :]

        lev_weight = np.zeros((len(lev),))
        dlev = lev[1:] - lev[:-1]
        lev_weight[1:-1] = (dlev[1:] + dlev[:-1]) / 2
        lev_weight[1]    = dlev[1]  / 2
        lev_weight[-1]   = dlev[-1] / 2
        lev_weight *= 100.0 / g0  # convert into density
        
        lev_weight = lev_weight[None, :, None, None]
      
        AR_vars = {
            'IWV' : np.sum(q * lev_weight, axis=1),
            'IVT' : np.sqrt(np.sum(q * U * lev_weight, axis=1)**2 + np.sum(q * V * lev_weight, axis=1)**2),
            'IWVKE' : np.sum(q * (U**2 + V**2) * lev_weight, axis=1),
        }

        with Dataset(out_filename, mode='w', format='NETCDF4_CLASSIC') as ds_out: 

            lat_dim = ds_out.createDimension('lat', len(lat))
            lon_dim = ds_out.createDimension('lon', len(lon)) 
            time_dim = ds_out.createDimension('time', None)

            var_lat = ds_out.createVariable('lat', np.float32, ('lat',))
            var_lon = ds_out.createVariable('lon', np.float32, ('lon',))
            
            var_lat[:] = lat
            var_lon[:] = lon

            for k, d in AR_vars.items():
                
                _var = ds_out.createVariable(k, np.float32, ('time', 'lat', 'lon'))
                _var[0:d.shape[0], :, :] = d

if __name__ == "__main__":

    if len(sys.argv) < 3:
        print("Needs two arguments.")
   
    computeAR(sys.argv[1], sys.argv[2]) 
