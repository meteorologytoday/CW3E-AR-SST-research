import numpy as np
import load_data
import netCDF4
from datetime import (datetime, timedelta, timezone)
import traceback
from pathlib import Path
import argparse

parser = argparse.ArgumentParser(
                    prog = 'plot_skill',
                    description = 'Plot prediction skill of GFS on AR.',
)

parser.add_argument('--beg-date', type=str, help='Date string: yyyy-mm-dd', required=True)
parser.add_argument('--end-date', type=str, help='Date string: yyyy-mm-dd', required=True)
parser.add_argument('--output', type=str, help='Output directory', default="")
parser.add_argument('--products', type=str, nargs="+", help='Forcast products.', default=["GFS",])
parser.add_argument('--lat-rng', type=float, nargs=2, help='Latitude  range', required=True)
parser.add_argument('--lon-rng', type=float, nargs=2, help='Longitude range. 0-360', required=True)
parser.add_argument('--save', action="store_true", help='Save the timeseries into an npz file')


args = parser.parse_args()

print(args)

# Configuration

beg_date = datetime.strptime(args.beg_date, "%Y-%m-%d").replace(tzinfo=timezone.utc)
end_date = datetime.strptime(args.end_date, "%Y-%m-%d").replace(tzinfo=timezone.utc)

total_days = (end_date - beg_date).days


print("Beg: ", beg_date)
print("End: ", end_date)
print("Total days: ", total_days)

if total_days <= 0:
    raise Exception("No days are avaiable.")



AR_vars = ["IWV", "IVT", "IWVKE", "sst", "mslhf", "msshf", "msnlwrf", "msnswrf", "mtpr", "mvimd", "t2m", "u10", "v10"]

aug_AR_vars = AR_vars + ['U',]

data_good = np.ones((total_days, len(aug_AR_vars)))

ts = {}

for AR_var in aug_AR_vars:
    ts[AR_var] = np.zeros((total_days,)) 

# Prep

lat = None
lon = None
lat_idx = None
lon_idx = None
wgt = None

# Load data
for d in range(total_days):
        
    _t = beg_date + timedelta(days=d)
        
            
    _data = {}

    for i, varname in enumerate(AR_vars):

        try:

            load_varname = varname

            # Load observation (the 'truth')
            info = load_data.getFileAndIndex("ERA5", _t, root_dir="data", varname=varname)

            print("Load file: ", info['filename'])


            with netCDF4.Dataset(info['filename'], "r") as ds:
                # decide range
                if lat is None:
                    
                    lat = ds.variables[info['varnames']['lat']][:]
                    lon = ds.variables[info['varnames']['lon']][:] % 360.0


                    lat_rng = np.array(args.lat_rng)
                    lon_rng = np.array(args.lon_rng) % 360
                    
                    if lat_rng[1] < lat_rng[0]:  # across lon=0
                        raise Exception("Latitude range should be lat_min, lat_max")

                    lat_idx = (lat_rng[0] < lat) & (lat < lat_rng[1])

                    if lon_rng[1] >= lon_rng[0]:
                        lon_idx = (lon_rng[0] < lon) & (lon < lon_rng[1])
                    
                    else:  # across lon=0
                        print("Across lon=0")
                        lon_idx = (lon_rng[0]) < lon | (lon < lon_rng[1])

                    lat = lat[lat_idx]
                    lon = lon[lon_idx]
                    wgt = np.cos(lat * np.pi / 180)


                _data[load_varname] = ds.variables[load_varname][info['idx'], lat_idx, lon_idx]
                #ts[varname][d] = np.average(np.average( _data[load_varname], axis=1), weights=wgt)

            keep_going = True

        except Exception as e:

            print(traceback.format_exc()) 
            print("Someting wrong happened when loading date: %s" % (_t.strftime("%Y-%m-%d"),))

            data_good[d, i] = 0.0

            keep_going = False


    if keep_going:
        _data['U'] = np.sqrt(_data['u10']**2 + _data['v10']**2)
        for varname in aug_AR_vars:
            ts[varname][d] = np.average(np.average( _data[load_varname], axis=1), weights=wgt)



for i, AR_var in enumerate(aug_AR_vars):
    data_good_idx = data_good[:, i] == 0.0
    ts[AR_var][data_good_idx] = np.nan


t_vec = np.array([ beg_date + timedelta(days=d) for d in range(total_days)], dtype="datetime64[s]")


if args.output != "":
            
    print("Output: %s" % (args.output,))
    #Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    with netCDF4.Dataset(args.output, 'w', format='NETCDF4') as ds:

        dim_time = ds.createDimension('time', len(t_vec))
        
        var_time = ds.createVariable('time', 'f4', ('time',))

        var_time[:] = t_vec
        
        for varname in aug_AR_vars:

            _var = ds.createVariable(varname, 'f4', ('time',))
            #value.units = 'Unknown'

            _var[:] = ts[varname]
            


