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
import pandas as pd
import re


parser = argparse.ArgumentParser(
                    prog = 'plot_skill',
                    description = 'Plot prediction skill of GFS on AR.',
)

parser.add_argument('--input-dir', type=str, help='Input file', required=True)
parser.add_argument('--beg-year', type=int, help='Input file', required=True)
parser.add_argument('--end-year', type=int, help='Input file', required=True)
parser.add_argument('--ncpu', type=int, help='Number of CPUs.', default=4)
parser.add_argument('--overwrite', help='If we overwrite the output', action="store_true")
parser.add_argument('--AR-algo', type=str, required=True, choices=["Rutz2017", "ANOM_LEN"])

parser.add_argument('--output-dir', type=str, help='Output dir', default="")
parser.add_argument('--title', type=str, help='Output title', default="")
parser.add_argument('--no-display', action="store_true")
args = parser.parse_args()
print(args)

if args.output_dir == "":
    args.output_dir = "%s/climanom_%04d-%04d" % (args.input_dir, args.beg_year, args.end_year, )
    
output_dir_stat = "%s/%s" % (args.output_dir, args.AR_algo)

print("Planned output dir 1: %s" % (args.output_dir,))
print("Planned output dir 2: %s" % (output_dir_stat,))

print("Create dir: %s" % (args.output_dir,))
Path(args.output_dir).mkdir(parents=True, exist_ok=True)
Path(output_dir_stat).mkdir(parents=True, exist_ok=True)

filename_clim = "%s/clim.nc" % (args.output_dir,)
filename_anom = "%s/anom.nc" % (args.output_dir,)
filename_ttl  = "%s/ttl.nc" % (args.output_dir,)

yrs = list(range(args.beg_year, args.end_year+1))

ds = ARstat_tool.loadDatasets(args.input_dir, yrs)

MLG_frc = (ds['MLG_frc_sw'] + ds['MLG_frc_lw'] + ds['MLG_frc_sh']  + ds['MLG_frc_lh'] + ds['MLG_frc_fwf'] + ds['MLG_rescale']).rename('MLG_frc')
MLG_adv = (ds['MLG_hadv'] + ds['MLG_vadv']).rename('MLG_adv')
MLG_diff = (ds['MLG_vdiff'] + ds['MLG_hdiff']).rename('MLG_diff')
MLG_nonfrc = (MLG_adv + MLG_diff + ds['MLG_ent']).rename('MLG_nonfrc')

#dMLDdt = ( ds['MLG_ent'] / ( ds[''] ) ).rename('dMLDdt')

ds = xr.merge(
    [
        ds,
        MLG_frc,
        MLG_nonfrc,
        MLG_adv,
        MLG_diff,
    ]
)

# Compute mean and anomaly

#def convertdt64todt(dt64):
#    return datetime.fromtimestamp( (dt64 - np.datetime64("1970-01-01")) / np.timedelta64(1, 's' ))

 
#ts_np64 = ds.time.to_numpy().astype("datetime64[ns]")
#ts_dt   = ts_np64.astype("datetime64[s]").astype(object) # need to be datetime64[s] first so that the object casting will result in a datetime object.


#tm_np64 = np.array([ np.datetime64('2022-01-01') + np.timedelta64(1,'D') * d for d in range(365) ]).astype("datetime64[ns]")

tm = pd.date_range("2021-01-01", "2021-12-31", freq="D", inclusive="both")
ts = pd.DatetimeIndex(ds.time.to_numpy())

ds_clim = []
ds_anom = []
ds_ttl  = []

target_varnames = list(ds.keys())

print("Target varnames: ", target_varnames)

# For testing:
# target_varnames = ["IWV", "IVT", "dMLTdt", "MLG_frc", "MLG_nonfrc", "MXLDEPTH"]


def doStat(varname):

    print("Doing stat of variable: %s" % (varname,))

    global tm
    _da_mean = xr.DataArray(
        name = varname,
        data = np.zeros((len(tm), len(ds.coords["lat"]), len(ds.coords["lon"]))),
        dims = ["time", "lat", "lon"],
        coords = {
            "time" : tm,
            "lat"  : ds.coords["lat"],
            "lon"  : ds.coords["lon"],
        }
    )
   

 
    _da_anom = xr.DataArray(
        name = varname,
        data = np.zeros((len(ts), len(ds.coords["lat"]), len(ds.coords["lon"]))),
        dims = ["time", "lat", "lon"],
        coords = {
            "time" : ts,
            "lat"  : ds.coords["lat"],
            "lon"  : ds.coords["lon"],
        }
    )
 
    _da_ttl = xr.DataArray(
        name = varname,
        data = np.zeros((len(ts), len(ds.coords["lat"]), len(ds.coords["lon"]))),
        dims = ["time", "lat", "lon"],
        coords = {
            "time" : ts,
            "lat"  : ds.coords["lat"],
            "lon"  : ds.coords["lon"],
        }
    )
   
    for i in range(len(ds.coords["lon"])):
        for j in range(len(ds.coords["lat"])):
            
            _var = ds[varname][:, j, i]

            xs = _var.to_numpy()
            
            tm, xm, xa, cnt, _ = anomalies.decomposeClimAnom(ts, xs)

            _da_mean[:, j, i] = xm
            _da_anom[:, j, i] = xa[:]
            _da_ttl[:, j, i] = xs[:]

            #print(ts[0], ";" , ts[-1])
    return varname, _da_mean, _da_anom, _da_ttl


if args.overwrite or (not path.exists(filename_clim)) or (not path.exists(filename_anom) or (not path.exists(filename_ttl))):

    print("Ready to multiprocess the statistical job.")
    with Pool(processes=args.ncpu) as pool:

        it = pool.imap(doStat, target_varnames)
        for (varname, _da_mean, _da_anom, _da_ttl) in it:

            ds_clim.append(_da_mean)
            ds_anom.append(_da_anom)

            if re.match('^map_', varname):
                print("Need total field of AR object mask variable: %s" % (varname,))
                ds_ttl.append(_da_ttl)



    print("Stat all done. Merge the outcome")

    ds_clim = xr.merge(ds_clim)
    ds_anom = xr.merge(ds_anom)
    ds_ttl  = xr.merge(ds_ttl)

    ds_clim.to_netcdf(filename_clim)
    ds_anom.to_netcdf(filename_anom)
    ds_ttl.to_netcdf(filename_ttl)

else:
    print("Files %s and %s already exists. Skip the computation." % (filename_clim, filename_anom, filename_ttl))

    ds_clim = xr.open_dataset(filename_clim)
    ds_anom = xr.open_dataset(filename_anom)
    ds_ttl  = xr.open_dataset(filename_ttl)
    

# Construct
time_constrains = [
    [1,],
    [2,],
    [3,],
    [4,],
    [5,],
    [6,],
    [1, 2, 3, 4, 5, 6],
]

time_labels = [
    "Oct",
    "Nov",
    "Dec",
    "Jan",
    "Feb",
    "Mar",
    "Oct-Mar",
]

    
if args.AR_algo == "Rutz2017":

    AR_cond  = (ds.IVT >= 250) & (ds.IWV >= 20.0)
    ARf_cond = (ds.IVT < 250) | ( (ds.IVT >= 250) & (ds.IWV < 20.0) )
  
elif args.AR_algo == "ANOM_LEN":

    AR_cond  = np.isfinite(ds.map_ANOM_LEN) & (ds.map_ANOM_LEN > 0)
    ARf_cond  = np.isfinite(ds.map_ANOM_LEN) & (ds.map_ANOM_LEN == 0)
 
elif args.AR_algo == "TOTIVT250":

    AR_cond  = np.isfinite(ds.map_TOTIVT250) & (ds.map_TOTIVT250 > 0)
    ARf_cond  = np.isfinite(ds.map_TOTIVT250) & (ds.map_TOTIVT250 == 0)


ds_stats = {}

for condition_name in ["clim", "ARf", "AR", "AR+ARf"]:

    print("Process condition: ", condition_name)

    _tmp = {}
    for varname, _ in ds_anom.items():
        _tmp[varname] = (["time", "lat", "lon", "stat"], np.zeros((len(time_constrains), len(ds.coords["lat"]), len(ds.coords["lon"]), 4)) )

    ds_stat = xr.Dataset(
        _tmp,

        coords = {
            "time" : time_labels,
            "lat"  : ds.coords["lat"],
            "lon"  : ds.coords["lon"],
            "stat" : ["mean", "std", "var", "cnt"],
        }
    )

    ds_stats[condition_name] = ds_stat


    for m, wm in enumerate(time_constrains): 

        time_cond = ds.time.dt.month.isin(watertime_tools.wm2m(wm))
        time_clim_cond = ds_clim.time.dt.month.isin(watertime_tools.wm2m(wm))

        
        if condition_name == "clim":
            
            _ds_ref = ds
            total_cond = time_cond

        elif condition_name == "AR+ARf":
        
            _ds_ref = ds_anom
            total_cond = time_cond
        
        elif condition_name == "AR":

            _ds_ref = ds_anom
            total_cond = time_cond & AR_cond
 
        elif condition_name == "ARf":

            _ds_ref = ds_anom
            total_cond = time_cond & ARf_cond
 
        else:
            raise Exception("Unknown condition_name: ", condition_name) 

        
        # Construct n-days in a row selection
        #ds.time.dt.month.isin(watertime_tools.wm2m(wm))

        _ds = _ds_ref.where(total_cond)

        for varname, _ in ds_stat.items():

            _data = _ds[varname].to_numpy()

            ds_stat[varname][m, :, :, 0] = np.nanmean(_data, axis=0) #_data.mean( dim="time", skipna=True)
            ds_stat[varname][m, :, :, 1] = np.nanstd(_data,  axis=0) #_data.std(  dim="time", skipna=True)
            ds_stat[varname][m, :, :, 2] = np.nanvar(_data,  axis=0) #_data.std(  dim="time", skipna=True)
            ds_stat[varname][m, :, :, 3] = np.nansum(np.isfinite(_data),  axis=0)#_data.std(  dim="time", skipna=True)
            

ds_stats["AR-ARf"] = ds_stats["AR"] - ds_stats["ARf"]

for k in ["clim", "AR", "ARf", "AR+ARf", "AR-ARf"]:
    output_filename = "%s/stat_%s.nc" % (output_dir_stat, k,)
    print("Writing output file: %s" % (output_filename,))
    
    ds_stats[k].to_netcdf(output_filename)

