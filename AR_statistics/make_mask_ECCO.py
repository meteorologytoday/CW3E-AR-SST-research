import netCDF4
import numpy as np
import xarray as xr

ref_filename = "./data/ECCO_LLC/ECCO_L4_POSTPROC_MXLANA_LLC0090GRID_DAILY_V4R4/MXLANA_day_mean_2001-01-01_ECCO_V4r4_native_llc0090.nc"
out_filename = "./mask_ECCO.nc"


ds = xr.open_dataset(ref_filename)
    

mask = xr.ones_like(ds.coords["XC"])

varnames = [
    'dMLTdt',
    'MLG_ttl',
    'MLG_frc_sw',
    'MLG_frc_lw',
    'MLG_frc_sh',
    'MLG_frc_lh',
    'MLG_frc_fwf',
    'MLG_rescale',
    'MLG_hadv',
    'MLG_vadv',
    'MLG_hdiff',
    'MLG_vdiff',
    'MLG_ent',
]


null_idx = None
for varname in varnames:
    print("Overlapping mask with variable: ", varname)
    if null_idx is None:
        null_idx = ds[varname][0, :, :, :].isnull()

    else:
        null_idx = null_idx & ds[varname][0, :, :, :].isnull()


mask = xr.where(null_idx, 0.0, 1.0).rename("mask")

print("Output: ", out_filename)
mask.to_netcdf(out_filename)
