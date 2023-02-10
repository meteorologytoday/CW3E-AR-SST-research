import numpy as np
import load_data
import netCDF4
from datetime import (datetime, timedelta)
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
parser.add_argument('--save', action="store_true", help='Save the timeseries into an npz file')


args = parser.parse_args()

print(args)

# Configuration

beg_date = datetime.strptime(args.beg_date, "%Y-%m-%d")
end_date = datetime.strptime(args.end_date, "%Y-%m-%d")

total_days = (end_date - beg_date).days


print("Beg: ", beg_date)
print("End: ", end_date)
print("Total days: ", total_days)

if total_days <= 0:
    raise Exception("No days are avaiable.")



AR_vars = ["IVT",]
data_good = np.ones((total_days, len(AR_vars)))

# Prep

lat = None
lon = None
lat_idx = None
lon_idx = None
wgt = None

AR_days = None

# Load data
for d in range(total_days):
        
    _t = beg_date + timedelta(days=d)
        
    for i, varname in enumerate(AR_vars):

        try:

            load_varname = varname

            # Load observation (the 'truth')
            info = load_data.getFileAndIndex("ERA5", _t, root_dir="data", varname=varname)

            print("Load file: ", info['filename'])

            _data = {}
            with netCDF4.Dataset(info['filename'], "r") as ds:
                # decide range
                if lat is None:
                    
                    lat = ds.variables[info['varnames']['lat']][:]
                    lon = ds.variables[info['varnames']['lon']][:] % 360.0

                    AR_days = np.zeros((len(lat), len(lon),), dtype=int)


                _data[load_varname] = ds.variables[load_varname][info['idx'], :, :]


                if load_varname == "IVT":
                    AR_days[_data["IVT"] > 250] += 1

        except Exception as e:

            print(traceback.format_exc()) 
            print("Someting wrong happened when loading date: %s" % (_t.strftime("%Y-%m-%d"),))

t_vec = np.array([ beg_date + timedelta(days=d) for d in range(total_days)], dtype="datetime64[s]")


if args.output != "":
            
    print("Output: %s" % (args.output,))
    #Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    with netCDF4.Dataset(args.output, 'w', format='NETCDF4') as ds:

        dim_time = ds.createDimension('time', len(t_vec))
        dim_lat  = ds.createDimension('lat',  len(lat))
        dim_lon  = ds.createDimension('lon',  len(lon))
        
        var_time = ds.createVariable('time', 'f4', ('time',))
        var_lat  = ds.createVariable('lat', 'f4', ('lat',))
        var_lon  = ds.createVariable('lon', 'f4', ('lon',))
        
        var_time[:] = t_vec
        var_lat[:]  = lat
        var_lon[:]  = lon
        
        var_AR_days = ds.createVariable('AR_days', 'i4', ('lat', 'lon', ))
        var_AR_days[:, :] = AR_days

