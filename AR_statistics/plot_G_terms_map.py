import numpy as np
#import fmon_tools, watertime_tools
import anomalies
import ARstat_tool
import xarray as xr

import traceback
from pathlib import Path
import argparse
from datetime import (timedelta, datetime, timezone)

parser = argparse.ArgumentParser(
                    prog = 'plot_skill',
                    description = 'Plot prediction skill of GFS on AR.',
)

parser.add_argument('--input-dir', type=str, help='Input file', required=True)
parser.add_argument('--output', type=str, help='Output file', default="")
parser.add_argument('--title', type=str, help='Output title', default="")
parser.add_argument('--no-display', action="store_true")
args = parser.parse_args()
print(args)

t_months = np.array([1, 2, 3, 4, 5, 6])

ds_stat = {}
for k in ["clim", "AR", "ARf", "AR-ARf", "AR+ARf"]:
    ds_stat[k] = xr.open_dataset("%s/stat_%s.nc" % (args.input_dir, k))

plot_infos_scnario = {

    "clim" : {
        "title" : "All",
    },

    "AR" : {
        "title" : "AR",
    },

    "ARf" : {
        "title" : "AR free",
    },

    "AR-ARf" : {
        "title" : "AR minus AR free",
    }

}

plot_infos = {
    "dMLTdt" : {
        "levels": np.linspace(-1, 1, 21)
    },

    "MLG_frc" : {
        "levels": np.linspace(-1, 1, 21)
    }, 
    "MLG_nonfrc" : {
        "levels": np.linspace(-1, 1, 21)
    }, 

}


plot_ylim = {

    "atmocn" : {
        "mean" : [-1.5, 1.5],
        "anom" : [-0.7, 0.7],
    },

    "atm" : {
        "mean" : [-1.5, 1.5],
        "anom" : [-0.3, 0.7],
    },

    "ocn" : {
        "mean" : [-0.5, 0.5],
        #"anom" : [-0.06, 0.01],
        "anom" : [-0.5, 0.5],
    },

}




# Plot data
print("Loading Matplotlib...")
import matplotlib as mpl
if args.no_display is False:
    mpl.use('TkAgg')
else:
    mpl.use('Agg')
    #mpl.rc('font', size=20)
    #mpl.rc('axes', labelsize=15)
     
 
  
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Rectangle
import matplotlib.transforms as transforms
from matplotlib.dates import DateFormatter
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


print("done")

cent_lon = 180.0

plot_lon_l = 110.0
plot_lon_r = 250.0
plot_lat_b = 5.0
plot_lat_t = 65.0

proj = ccrs.PlateCarree(central_longitude=cent_lon)
proj_norm = ccrs.PlateCarree()

varnames = ["dMLTdt", "MLG_frc", "MLG_nonfrc"]
fig, ax = plt.subplots(3, len(t_months), figsize=(4*len(t_months), 2*3), subplot_kw=dict(projection=proj))

coords = ds_stat["clim"].coords
cmap = "bwr"
for m, mon in enumerate(t_months):

    _ax = ax[:, m]
  
    s = "AR-ARf"

    for i, varname in enumerate(varnames):
        mappable = _ax[i].contourf(coords["lon"], coords["lat"], ds_stat["AR-ARf"][varname][m, :, :, 0] * 1e6, levels=plot_infos["dMLTdt"]["levels"], cmap=cmap, extend="both", transform=proj_norm)

        _diff = np.abs(ds_stat["AR-ARf"][varname][m, :, :, 0].to_numpy())
        _std1 = ds_stat["AR"][varname][m, :, :, 1].to_numpy()
        _std2 = ds_stat["ARf"][varname][m, :, :, 1].to_numpy()
        
        _std3 = ds_stat["AR+ARf"][varname][m, :, :, 1].to_numpy()

        _dot = _diff * 0 
        _significant_idx =  (_diff > _std1) | (_diff > _std2)
        #_significant_idx =  (_diff > _std3)

        print("_diff: ", _diff[0, 0:20])
        print("_std3: ", _std1[0, 0:20])

        _dot[ _significant_idx                 ] = 0.75
        _dot[ np.logical_not(_significant_idx) ] = 0.25

        _ax[i].contourf(coords["lon"], coords["lat"], _dot, colors='none', levels=[0, 0.5, 1], hatches=[None, "..."], transform=proj_norm)
       
    for __ax in _ax: 
        __ax.set_global()
        __ax.gridlines()
        __ax.coastlines()

        gl = __ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=1, color='gray', alpha=0.5, linestyle='--')
        gl.xlabels_top = False
        gl.ylabels_left = False
        #gl.xlocator = mticker.FixedLocator(np.arange(-180, 181, 60))
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 10, 'color': 'black'}
        gl.ylabel_style = {'size': 10, 'color': 'black'}

        __ax.set_extent([plot_lon_l, plot_lon_r, plot_lat_b, plot_lat_t], crs=proj_norm)

if not args.no_display:
    plt.show()

if args.output != "":
    
    fig.savefig(args.output, dpi=200)

