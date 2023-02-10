import netCDF4
import numpy as np
import xarray as xr

ref_filename = "./data/ECCO_LLC/ECCO_L4_POSTPROC_MXLANA_LLC0090GRID_DAILY_V4R4/MXLANA_day_mean_2001-01-01_ECCO_V4r4_native_llc0090.nc"
out_filename = "./mask_ECCO.nc"


ds = xr.open_dataset(ref_filename)
    

mask = xr.ones_like(ds.coords["XC"])
mask = xr.where(ds.MLGs_ttl[0, :, :, :].isnull() & ds.MLGs_hadv[0, :, :, :].isnull(), 0.0, 1.0).rename("mask")

print("Output: ", out_filename)
mask.to_netcdf(out_filename)
