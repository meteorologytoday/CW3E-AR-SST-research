from multiprocessing import Pool
import numpy as np
import os.path as path
#import fmon_tools, watertime_tools
import anomalies
import ARstat_tool
import xarray as xr
import watertime_tools
import traceback
from pathlib import Path
import argparse
from datetime import (timedelta, datetime, timezone)

parser = argparse.ArgumentParser(
                    prog = 'plot_skill',
                    description = 'Plot prediction skill of GFS on AR.',
)

parser.add_argument('--input-dir', type=str, help='Input file', required=True)
parser.add_argument('--beg-year', type=int, help='Input file', required=True)
parser.add_argument('--end-year', type=int, help='Input file', required=True)
parser.add_argument('--watermonths', type=int, nargs="+", help='Input file', default=[1, 2, 3, 4, 5, 6])

parser.add_argument('--output-dir', type=str, help='Output dir', default="")
parser.add_argument('--title', type=str, help='Output title', default="")
parser.add_argument('--no-display', action="store_true")
args = parser.parse_args()
print(args)

yrs = list(range(args.beg_year, args.end_year+1))
ds = ARstat_tool.loadDatasets(args.input_dir, yrs)


if args.output_dir == "":
    args.output_dir = args.input_dir

output_filename = "%s/AR_simple_statistics_%04d-%04d.nc" % (args.output_dir, args.beg_year, args.end_year, )

print("Planned output file: ", output_filename)

output_varnames = ["IVT", "IWV", "MEAN_VEL", "mtpr", "absU10", "u10", "v10"]

MEAN_VEL = (ds.IVT / ds.IWV).rename("MEAN_VEL")
absU10 = ((ds["u10"]**2 + ds["v10"]**2)**0.5).rename("absU10")

ds = xr.merge([ds.IWV, ds.IVT, ds.u10, ds.v10, ds.mtpr, MEAN_VEL, absU10])


wms = args.watermonths
stat = ["clim", "std", "AR"]

empty_arr = np.zeros((len(wms), len(stat), len(ds.coords["lat"]), len(ds.coords["lon"])))

dims = ["watermonth", "stat", "lat", "lon"]


data_vars = {
    varname : (dims, empty_arr.copy())
    for varname in output_varnames
}

output_ds = xr.Dataset(

    data_vars = data_vars,

    coords = dict(
        watermonth = wms,
        stat       = stat,
        lat        = ds.coords["lat"],
        lon        = ds.coords["lon"],
    ),

)


print(list(ds.keys()))

for wm in wms:


    print("Doing statistics of watermonth: ", wm)


    cond_clim = ds.time.dt.month.isin(watertime_tools.wm2m(wm))
    cond_AR = (ds.IVT >= 250) & (ds.IWV >= 20.0) & ds.time.dt.month.isin(watertime_tools.wm2m(wm))
    ds_AR = ds.where(cond_AR)
    ds_clim = ds.where(cond_clim)
    # Compute the mean IVT, IWV and v = ds.IVT / ds.IWV
 

    for varname in output_varnames:

        #_mean = ds_subset[varname].mean(dim=("time", "lon"), skipna=True)
        output_ds[varname][wm-1, 0, :, :] = ds_clim[varname].mean(dim=("time", ), skipna=True)
        output_ds[varname][wm-1, 1, :, :] = ds_clim[varname].std(dim=("time", ), skipna=True)
        output_ds[varname][wm-1, 2, :, :] = ds_AR[varname].mean(dim=("time", ), skipna=True)



print("Generating output file: ", output_filename)
output_ds.to_netcdf(output_filename)
