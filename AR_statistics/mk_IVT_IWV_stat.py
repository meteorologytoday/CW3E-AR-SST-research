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

parser.add_argument('--output-dir', type=str, help='Output dir', default="")
parser.add_argument('--title', type=str, help='Output title', default="")
parser.add_argument('--no-display', action="store_true")
args = parser.parse_args()
print(args)

yrs = list(range(args.beg_year, args.end_year+1))
ds = ARstat_tool.loadDatasets(args.input_dir, yrs)


if args.output_dir == "":
    args.output_dir = args.input_dir

output_filename = "%s/lat_statistics_%04d-%04d.nc" % (args.output_dir, args.beg_year, args.end_year, )

print("Planned output file: ", output_filename)

ds = xr.merge([ds.IWV, ds.IVT, ])
wms = [1, 2, 3, 4, 5, 6]
        
empty_arr = np.zeros((len(wms), len(ds.coords["lat"]),))

dims = ["watermonth", "lat"]
output_ds = xr.Dataset(

    data_vars = dict(
        IVT = (dims, empty_arr.copy()),
        IWV = (dims, empty_arr.copy()),
        MEAN_VEL = (dims, empty_arr.copy()),
    ),

    coords = dict(
        watermonth = wms,
        lat        = ds.coords["lat"],
    ),

)



for wm in wms:


    print("Doing statistics of watermonth: ", wm)

    cond = (ds.IVT >= 250) & (ds.IWV >= 20.0) & ds.time.dt.month.isin(watertime_tools.wm2m(wm))
    ds_subset = ds.where(cond)

    # Compute the mean IVT, IWV and v = ds.IVT / ds.IWV
    
    MEAN_VEL = (ds_subset.IVT / ds_subset.IWV).rename("MEAN_VEL")
    
    ds_subset = xr.merge([ds_subset, MEAN_VEL])

    latmeans = {}
    for varname in ["IVT", "IWV", "MEAN_VEL"]:
        _mean = ds_subset[varname].mean(dim=("time", "lon"), skipna=True)
        output_ds[varname][wm-1, :] = ds_subset[varname].mean(dim=("time", "lon"), skipna=True)



print("Generating output file: ", output_filename)
output_ds.to_netcdf(output_filename)
