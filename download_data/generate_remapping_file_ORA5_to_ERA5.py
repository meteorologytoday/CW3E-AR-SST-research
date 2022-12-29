import netCDF4
import numpy as np
import grid_tools
from cmd_tools import pleaseRun

print(" *** This program generates SCRIP Grid formatted files and regrid weighting file of ORA5 and ERA5 data. *** ")

ORA5_sample_file = "data/ORA5/processed/ORA5_NillerKrausMixedLayerDynamics_somxl010_1990-10.nc"
ERA5_sample_file = "data/ERA5/sfc/ERA5_sfc_2001-01-01.nc"

ERA5_grid_file = "grid_ERA5.nc"
ORA5_grid_file = "grid_ORA5.nc"
algo = "bilinear"

remapping_file = "remap_ORA5_to_ERA5_%s.nc" % (algo,)

print("ERA5 sample file: %s => %s" % (ERA5_sample_file, ERA5_grid_file,))
print("ORA5 sample file: %s => %s" % (ORA5_sample_file, ORA5_grid_file,))

print("Remapping algo: %s" % (algo,))
print("Remapping file: %s" % (remapping_file,))

print("Making grid file of ERA5 in SCRIP format... ")
ds = netCDF4.Dataset(ERA5_sample_file, "r")
lat = ds.variables['latitude'][:]
lon = ds.variables['longitude'][:]
data = ds.variables['sst'][:]
imask = np.logical_not(data.mask).astype(np.int32)
grid_tools.genSCRIPFileWithCornerMissing_regular(lat=lat, lon=lon, imask=imask, output_filename=ERA5_grid_file, rev_lat=True)

print("Making grid file of ORA5 in SCRIP format... ")
ds = netCDF4.Dataset("data/ORA5/processed/ORA5_NillerKrausMixedLayerDynamics_somxl010_1990-10.nc", "r")
lat = ds.variables['lat'][:]
lon = ds.variables['lon'][:]
data = ds.variables['MLD'][:]
imask = np.logical_not(data.mask).astype(np.int32)
grid_tools.genSCRIPFileWithCornerMissing_curvlinear(lat=lat, lon=lon, imask=imask, output_filename=ORA5_grid_file)

print("Making regrid weighting file")

pleaseRun("time ncremap -a %s --src_grd=%s --dst_grd=%s --msk_src=imask --msk_dst=imask -m %s" % (algo, ORA5_grid_file, ERA5_grid_file, remapping_file))

