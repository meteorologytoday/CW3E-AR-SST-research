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
parser.add_argument('--output-dir', type=str, help='Output directory', default="output_ts")
parser.add_argument('--fcst',     type=int, help='Forecast hours', default=240)
parser.add_argument('--products', type=str, nargs="+", help='Forcast products.', default=["GFS",])
parser.add_argument('--lat-rng', type=float, nargs=2, help='Latitude  range', required=True)
parser.add_argument('--lon-rng', type=float, nargs=2, help='Longitude range. 0-360', required=True)
parser.add_argument('--save', action="store_true", help='Save the timeseries into an npz file')


args = parser.parse_args()

print(args)

# Configuration

beg_date = datetime.strptime(args.beg_date, "%Y-%m-%d")
end_date = datetime.strptime(args.end_date, "%Y-%m-%d")

fcst_time = timedelta(hours=args.fcst)

fcst_beg_date = beg_date - fcst_time
fcst_end_date = end_date - fcst_time

total_days = (end_date - beg_date).days


print("Beg: ", beg_date)
print("End: ", end_date)
print("Total days: ", total_days)

if total_days <= 0:
    raise Exception("No days are avaiable.")



AR_vars = ["IWV", "IVT", "IWVKE", "HGT_850mb", "HGT_500mb", "LHTFL", "SHTFL", "PRATE"]

data_good = np.ones((total_days, len(AR_vars)))


ts = {}

for AR_var in AR_vars:
    ts[AR_var] = np.zeros((2, total_days,)) # first index: [0=obs, 1=fcst], second index: days

# Prep
datas = {}


lat = None
lon = None
lat_idx = None
lon_idx = None
wgt = None





# Load data
for d in range(total_days):
        
    _t = beg_date + timedelta(days=d)
        
    for i, varname in enumerate(AR_vars):

        try:

            if varname in [ "HGT_500mb" , "HGT_850mb" ]:
                load_varname = "HGT"
            else:
                load_varname = varname

            # Load observation (the 'truth')
            infos = {
                'obs'  : load_data.getFileAndIndex("GFS", _t,             root_dir="data", varname=varname, fcst=0        ),
                'fcst' : load_data.getFileAndIndex("GFS", _t - fcst_time, root_dir="data", varname=varname, fcst=args.fcst),
            }

            _data = {}
            
            for k, info in infos.items():
                
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

                    if load_varname == "HGT":
                        _data[k] = ds.variables[load_varname][info['idx'], 0, lat_idx, lon_idx]
                    else:    
                        _data[k] = ds.variables[load_varname][info['idx'], lat_idx, lon_idx]
            

            # compute skill
            #_diff = _data['fcst'] - _data['obs']

            ts[varname][0, d] = np.average(np.average( _data['obs'], axis=1), weights=wgt)
            ts[varname][1, d] = np.average(np.average( _data['fcst'], axis=1), weights=wgt)

        except Exception as e:

            print(traceback.format_exc()) 
            print("Someting wrong happened when loading date: %s" % (_t.strftime("%Y-%m-%d"),))

            data_good[d, i] = 0.0

for i, AR_var in enumerate(AR_vars):
    data_good_idx = data_good[:, i] == 0.0
    ts[AR_var][:, data_good_idx] = np.nan


t_vec = np.array([ beg_date + timedelta(days=d) for d in range(total_days)], dtype="datetime64[s]")


if args.save:

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    filename = "output_ts/%s_%s_f%03d.npz" % (beg_date.strftime("%Y-%m-%d"), end_date.strftime("%Y-%m-%d"), args.fcst)
    print("Output: %s" % (filename,))
    np.savez(filename, time=t_vec, **ts)

"""
#print(lat)
#print(lon)

# Plot data
print("Loading Matplotlib...")
import matplotlib
if args.no_display is False:
    matplotlib.use('TkAgg')
else:
    matplotlib.use('Agg')
    
import matplotlib.pyplot as plt
from matplotlib import cm
import cartopy.crs as ccrs
print("done")



fig, ax = plt.subplots(len(AR_vars), 1, figsize=(6, 4))

for i, AR_var in enumerate(AR_vars):
    _ax = ax[i]

    _ax.plot(t_vec, skills[AR_var])
    _ax.set_title("Forecast error of %s" % (AR_var, ))


if not args.no_display:
    plt.show()

if args.output_dir != "":
    
    print("Create dir: %s" % (args.output_dir,))
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    fig.savefig("%s/skills.png" % (args.output_dir,), dpi=200)

"""
