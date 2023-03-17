import numpy as np
import os.path as path
#import fmon_tools, watertime_tools
import anomalies
import ARstat_tool

print("Loading xarray")
import xarray as xr
print("done")

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
parser.add_argument('--beg-year', type=int, help='Input file', default=0)
parser.add_argument('--end-year', type=int, help='Input file', default=0)
parser.add_argument('--selected-years', type=int, nargs="*", help='separate years. If nonempty then will use this instead of --beg-year and --end-year', default=[])

parser.add_argument('--output-dir', type=str, help='Output dir', default="")
parser.add_argument('--suffix', type=str, help='Suffix for output netcdf', default="")
parser.add_argument('--title', type=str, help='Output title', default="")
parser.add_argument('--no-display', action="store_true")
args = parser.parse_args()
print(args)


if len(args.selected_years) != 0:
    yrs = args.selected_years
else:
    yrs = list(range(args.beg_year, args.end_year+1))

ds = ARstat_tool.loadDatasets(args.input_dir, yrs)


if args.output_dir == "":
    args.output_dir = args.input_dir

if len(args.selected_years) != 0:
    output_filename = "%s/AR_interannual_statistics_selected_years%s.nc" % (args.output_dir, args.suffix)
else:
    output_filename = "%s/AR_interannual_statistics_%04d-%04d%s.nc" % (args.output_dir, args.beg_year, args.end_year, args.suffix)

print("Planned output file: ", output_filename)

ds = xr.merge([ds.IWV, ds.IVT, ])

valid_wms = [1, 2, 3, 4, 5, 6]
        
empty_arr = np.zeros((len(ds.coords["lat"]), len(ds.coords["lon"])))

dims = ["lat", "lon"]

output_ds = xr.Dataset(

    data_vars = dict(
        IVT = (dims, empty_arr.copy()),
        IWV = (dims, empty_arr.copy()),
        MEAN_VEL = (dims, empty_arr.copy()),
        count = (dims, empty_arr.copy()),
    ),

    coords = dict(
        lat        = ds.coords["lat"],
        lon        = ds.coords["lon"],
    ),

)



for i, yr in enumerate(yrs):


    print("Doing statistics of year: ", yr)

    cond = (ds.IVT >= 250) & (ds.IWV >= 20.0) & ds.time.dt.month.isin(watertime_tools.wm2m(wms)) & (ds.time.dt.year == yr )
    ds_subset = ds.where(cond)

    # Compute the mean IVT, IWV and v = ds.IVT / ds.IWV
    
    MEAN_VEL = (ds_subset.IVT / ds_subset.IWV).rename("MEAN_VEL")
    
    count = cond.sum(dim=("time",), skipna=True).rename("count")

    ds_subset = xr.merge([ds_subset, MEAN_VEL, count])

    for varname in ["IVT", "IWV", "MEAN_VEL"]:
        #_mean = ds_subset[varname].mean(dim=("time", "lon"), skipna=True)
        output_ds[varname][i, :, :] = ds_subset[varname].mean(dim=("time", ), skipna=True)

    output_ds["count"][i, :, :] = cond.sum(dim=("time",), skipna=True)


print("Generating output file: ", output_filename)
output_ds.to_netcdf(output_filename)
