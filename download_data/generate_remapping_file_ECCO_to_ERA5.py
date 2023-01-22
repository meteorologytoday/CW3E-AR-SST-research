import netCDF4
import numpy as np
import grid_tools
from cmd_tools import pleaseRun

print(" *** This program generates SCRIP Grid formatted files and regrid weighting file of ECCO and ERA5 data. *** ")

ECCO_sample_file = "data/ECCO/ECCO_L4_TEMP_SALINITY_05DEG_DAILY_V4R4/OCEAN_TEMPERATURE_SALINITY_day_mean_2004-12-26_ECCO_V4r4_latlon_0p50deg.nc"
ERA5_sample_file = "data/ERA5/sfc/ERA5_sfc_2001-01-01.nc"

ERA5_grid_file = "grid_ERA5.nc"
ECCO_grid_file = "grid_ECCO.nc"
algo = "bilinear"

remapping_file = "remap_ECCO_to_ERA5_%s.nc" % (algo,)

print("ERA5 sample file: %s => %s" % (ERA5_sample_file, ERA5_grid_file,))
print("ECCO sample file: %s => %s" % (ECCO_sample_file, ECCO_grid_file,))

print("Remapping algo: %s" % (algo,))
print("Remapping file: %s" % (remapping_file,))

print("Making grid file of ERA5 in SCRIP format... ")
ds = netCDF4.Dataset(ERA5_sample_file, "r")
lat = ds.variables['latitude'][:]
lon = ds.variables['longitude'][:]
data = ds.variables['sst'][:]
imask = np.logical_not(data.mask).astype(np.int32)
grid_tools.genSCRIPFileWithCornerMissing_regular(lat=lat, lon=lon, imask=imask, output_filename=ERA5_grid_file, rev_lat=True)

print("Making grid file of ECCO in SCRIP format... ")
ds = netCDF4.Dataset(ECCO_sample_file, "r")
lat = ds.variables['latitude'][:]
lon = ds.variables['longitude'][:]
data = ds.variables['THETA'][0, 0, :, :]
imask = np.logical_not(data.mask).astype(np.int32)
grid_tools.genSCRIPFileWithCornerMissing_regular(lat=lat, lon=lon, imask=imask, output_filename=ECCO_grid_file)

print("Making regrid weighting file")

pleaseRun("time ncremap -a %s --src_grd=%s --dst_grd=%s --msk_src=imask --msk_dst=imask -m %s" % (algo, ECCO_grid_file, ERA5_grid_file, remapping_file))

