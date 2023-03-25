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

parser.add_argument('--input', type=str, help='Input file', required=True)
parser.add_argument('--output', type=str, help='Input file', default="")
parser.add_argument('--title', type=str, help='Output title', default="")
parser.add_argument('--no-display', action="store_true")
parser.add_argument('--markers', action="store_true")
args = parser.parse_args()
print(args)


ds = xr.open_dataset(args.input)

AR_mean= ds["count"].mean(dim=("time",))
AR_std = ds["count"].std(dim=("time",))


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

plot_lon_l = 120.0
plot_lon_r = 240.0
plot_lat_b = 10.0
plot_lat_t = 60.0

proj = ccrs.PlateCarree(central_longitude=cent_lon)
proj_norm = ccrs.PlateCarree()

fig, ax = plt.subplots(
    1, 1,
    figsize=(6, 4),
    subplot_kw=dict(projection=proj),
    gridspec_kw=dict(hspace=0, wspace=0.2),
    constrained_layout=False,
)

ax.set_title("AR days per year (shading) and variability (contour)")

coords = ds.coords
cmap = cm.get_cmap("GnBu")

#cmap.set_over("yellow")
#cmap.set_under("red")

mean_levels = np.linspace(0, 70, 15)
std_levels = np.linspace(0, 20, 11)

mappable = ax.contourf(coords["lon"], coords["lat"], AR_mean, levels=mean_levels, cmap=cmap, extend="max", transform=proj_norm)
cs = ax.contour(coords["lon"], coords["lat"], AR_std, levels=std_levels, colors="k", transform=proj_norm, linewidths=1)

ax.clabel(cs, fmt="%d")


ax.set_global()
ax.coastlines()
ax.set_extent([plot_lon_l, plot_lon_r, plot_lat_b, plot_lat_t], crs=proj_norm)

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.5, linestyle='--')

gl.xlabels_top   = False
gl.ylabels_right = False

#gl.xlocator = mticker.FixedLocator(np.arange(-180, 181, 30))
gl.xlocator = mticker.FixedLocator([120, 150, 180, -150, -120])#np.arange(-180, 181, 30))
gl.ylocator = mticker.FixedLocator([10, 20, 30, 40, 50])

gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 10, 'color': 'black'}
gl.ylabel_style = {'size': 10, 'color': 'black'}

    
cb = plt.colorbar(mappable, ax=ax, orientation="vertical", pad=0.01, shrink=0.5)
cb.ax.set_ylabel("AR days per wateryear")

if not args.no_display:
    plt.show()

if args.output != "":
    
    fig.savefig(args.output, dpi=200)

