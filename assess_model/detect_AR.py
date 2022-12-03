import numpy as np
import load_data
import netCDF4
from datetime import (datetime, timedelta)

from pathlib import Path
import argparse

parser = argparse.ArgumentParser(
                    prog = 'plot_skill',
                    description = 'Plot prediction skill of GFS on AR.',
)

parser.add_argument('--beg-date', type=str, help='Date string: yyyy-mm-dd', required=True)
parser.add_argument('--end-date', type=str, help='Date string: yyyy-mm-dd', required=True)
parser.add_argument('--lat-rng', type=float, nargs=2, help='Latitude  range', required=True)
parser.add_argument('--lon-rng', type=float, nargs=2, help='Longitude range. 0-360', required=True)
parser.add_argument('--save-npz', type=str, help='Save the timeseries into an npz file', default="")


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

data_good = np.ones((total_days,))

AR_vars = ["IWV", "IVT", "IWVKE"]

regavg = {}

for AR_var in AR_vars:
    regavg[AR_var] = np.zeros((total_days,))


# prep

lat = None
lon = None
lat_idx = None
lon_idx = None
wgt = None

# Load data
for d in range(total_days):

    try:

        _t = beg_date + timedelta(days=d)

        info = load_data.getFileAndIndex("ERA5", _t, root_dir="data")

        print("Load file: ", info['filename'])

        _data = {}
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

            for varname in AR_vars:
                _data[varname] = ds.variables[varname][info['idx'], lat_idx, lon_idx]
    

        # compute skill
        for AR_var in AR_vars:
            regavg[AR_var][d] = np.average(np.average( _data[AR_var], axis=1), weights=wgt)

    except Exception as e:

        print(e)
        print("Someting wrong happened when loading date: %s" % (_t.strftime("%Y-%m-%d"),))

        data_good[d] = 0.0

for AR_var in AR_vars:
    regavg[AR_var][data_good == 0.0] = np.nan


t_vec = np.array([ beg_date + timedelta(days=d) for d in range(total_days)], dtype="datetime64[s]")


if args.save_npz != "":
    np.savez(args.save_npz, time=t_vec, **regavg)


