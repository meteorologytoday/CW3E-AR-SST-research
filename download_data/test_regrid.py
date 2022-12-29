import netCDF4
import numpy as np
import grid_tools

print("start 2")
ds = netCDF4.Dataset("data/ERA5/sfc/ERA5_sfc_2012-10-18.nc", "r")
lat = ds.variables['latitude'][:]
lon = ds.variables['longitude'][:]
data = ds.variables['sst'][:]
imask = np.logical_not(data.mask).astype(np.int32)

grid_tools.genSCRIPFileWithCornerMissing_regular(lat=lat, lon=lon, imask=imask, output_filename="grid_ECMWF.nc")

print("done")

print("start 1")

ds = netCDF4.Dataset("data/ORA5/processed/ORA5_NillerKrausMixedLayerDynamics_somxl010_1990-10.nc", "r")
lat = ds.variables['lat'][:]
lon = ds.variables['lon'][:]
data = ds.variables['MLD'][:]
imask = np.logical_not(data.mask).astype(np.int32)

grid_tools.genSCRIPFileWithCornerMissing_curvlinear(lat=lat, lon=lon, imask=imask, output_filename="grid_NEMO.nc")

print("done")


