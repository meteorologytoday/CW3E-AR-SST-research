import numpy as np
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

if args.output_dir == "":
    args.output_dir = args.input_dir

print("Planned output dir: %s" % (args.output_dir,))

yrs = list(range(args.beg_year, args.end_year+1))

ds = ARstat_tool.loadDatasets(args.input_dir, yrs)

MLG_frc = (ds['MLG_frc_sw'] + ds['MLG_frc_lw'] + ds['MLG_frc_sh']  + ds['MLG_frc_lh'] + ds['MLG_frc_fwf'] + ds['MLG_rescale']).rename('MLG_frc')
MLG_adv = (ds['MLG_hadv'] + ds['MLG_vadv']).rename('MLG_adv')
MLG_diff = (ds['MLG_vdiff'] + ds['MLG_hdiff']).rename('MLG_diff')
MLG_nonfrc = (MLG_adv + MLG_diff + ds['MLG_ent']).rename('MLG_nonfrc')

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
def convertdt64todt(dt64):
    return datetime.fromtimestamp( (dt64 - np.datetime64("1970-01-01")) / np.timedelta64(1, 's' ))
 
ts_np64 = ds.time.to_numpy().astype("datetime64[ns]")
ts_dt   = ts_np64.astype("datetime64[s]").astype(object) # need to be datetime64[s] first so that the object casting will result in a datetime object.


#for t in ts_dt:
#    print(t)

tm_np64 = np.array([ np.datetime64('2022-01-01') + np.timedelta64(1,'D') * d for d in range(365) ]).astype("datetime64[ns]")


ds_mean = []
ds_anom = []

#target_varnames = ["IVT", "dMLTdt", "MLG_frc", "MLG_nonfrc"]
target_varnames = ds.keys()

assist = None # For climatology decomposition speed up





for _, varname in enumerate(target_varnames):

    print("Doing statistics of ", varname)
    
    _da_mean = xr.DataArray(
        name = varname,
        data = np.zeros((len(tm_np64), len(ds.coords["lat"]), len(ds.coords["lon"]))),
        dims = ["time", "lat", "lon"],
        coords = {
            "time" : tm_np64,
        }
    )
   

 
    _da_anom = xr.DataArray(
        name = varname,
        data = np.zeros((len(ts_np64), len(ds.coords["lat"]), len(ds.coords["lon"]))),
        dims = ["time", "lat", "lon"],
        coords = {
            "time" : ts_np64,
        }
    )
    
    for i in range(len(ds.coords["lon"])):
        for j in range(len(ds.coords["lat"])):
            
            _var = ds[varname][:, j, i]

            xs = _var.to_numpy()
            
            tm, xm, xa, cnt, assist = anomalies.decomposeClimAnom(ts_dt, xs, assist=assist)

            _da_mean[:, j, i] = xm
            _da_anom[:, j, i] = xa[:]

    ds_anom.append(_da_anom)
    ds_mean.append(_da_mean)

ds_anom = xr.merge(ds_anom)
ds_mean = xr.merge(ds_mean)


ds_anom.to_netcdf("ds_anom.nc")
ds_mean.to_netcdf("ds_mean.nc")

# Construct
t_months = np.arange(1, 7)

ds_stats = {}

for condition_name, (IVT_min, IVT_max), in [
    ("clim",  (0, np.inf)   ),
    ("ARf",   (0, 250)      ),
    ("AR",    (250, np.inf )),
    ("AR+ARf", (0, np.inf  )),
]:

    print("Process condition: ", condition_name)

    _tmp = {}
    for varname, _ in ds_anom.items():
        _tmp[varname] = (["time", "lat", "lon", "stat"], np.zeros((len(t_months), len(ds.coords["lat"]), len(ds.coords["lon"]), 4)) )

    ds_stat = xr.Dataset(
        _tmp,

        coords = {
            "time" : t_months,
            "lat"  : ds.coords["lat"],
            "lon"  : ds.coords["lon"],
            "stat" : ["mean", "std", "var", "cnt"],
        }
    )

    ds_stats[condition_name] = ds_stat

    if condition_name == "clim":
        IVT_cond = np.isfinite(ds.IVT)

    elif condition_name == "AR":
        IVT_cond = (ds.IVT >= 250) & (ds.IWV >= 20.0)
    
    elif condition_name == "ARf":
        IVT_cond = (ds.IVT < IVT_min) | ( (ds.IVT >= 250) & (ds.IWV < 20.0) )
  
    elif condition_name == "AR+ARf":
        IVT_cond = np.isfinite(ds.IVT)
          
   
    #IVT_cond_ndays = ifNdaysInARow(IVT_cond)
 
    #print("SHAPE OF IVT_COND: ", IVT_cond.shape)
 
    for m, wm in enumerate(t_months): 

        time_cond = ds.time.dt.month.isin(watertime_tools.wm2m(wm))
        time_clim_cond = ds_mean.time.dt.month.isin(watertime_tools.wm2m(wm))

        
        if condition_name == "clim":
            #_ds_ref = ds_mean
            #total_cond = time_clim_cond
            
            _ds_ref = ds
            total_cond = time_cond
            
            #print( "sum time_clim_cond: ",  np.nansum( time_clim_cond.astype(int).to_numpy() ) )
            #print( "sum time_cond:      ",  np.nansum( time_cond.astype(int).to_numpy() ) )


        elif condition_name == "AR+ARf":
            _ds_ref = ds_anom
            total_cond = time_cond
        
        elif condition_name in ["AR", "ARf"]:
            _ds_ref = ds_anom
            total_cond = time_cond & IVT_cond
 
        else:
            raise Exception("Unknown condition_name: ", condition_name) 

        
        # Construct n-days in a row selection
        #ds.time.dt.month.isin(watertime_tools.wm2m(wm))

        _ds = _ds_ref.where(total_cond)


        #if condition_name == "clim":
        #    print("MEAN = ")
        #    print( np.nanmean( _ds["dMLTdt"].to_numpy(), axis=0) ) 
        #    print("Valid points: ", np.sum(np.isfinite(np.nanmean( _ds["dMLTdt"].to_numpy(), axis=0) ))) 

        for varname, _ in ds_stat.items():

            _data = _ds[varname].to_numpy()
            ds_stat[varname][m, :, :, 0] = np.nanmean(_data, axis=0) #_data.mean( dim="time", skipna=True)
            ds_stat[varname][m, :, :, 1] = np.nanstd(_data,  axis=0) #_data.std(  dim="time", skipna=True)
            ds_stat[varname][m, :, :, 2] = np.nanvar(_data,  axis=0) #_data.std(  dim="time", skipna=True)
            ds_stat[varname][m, :, :, 3] = np.nansum(np.isfinite(_data),  axis=0)#_data.std(  dim="time", skipna=True)
            

ds_stats["AR-ARf"] = ds_stats["AR"] - ds_stats["ARf"]

for k in ["clim", "AR", "ARf", "AR+ARf", "AR-ARf"]:
    output_filename = "%s/stat_%s.nc" % (args.output_dir, k,)
    print("Writing output file: %s" % (output_filename,))
    
    ds_stats[k].to_netcdf(output_filename)

