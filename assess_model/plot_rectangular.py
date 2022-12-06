import argparse
import numpy as np

parser = argparse.ArgumentParser(
                    prog = 'plot_rectangular',
                    description = 'Plot a box on world map',
)

parser.add_argument('--no-display', action="store_true")
parser.add_argument('--lat-rng', type=float, nargs=2, help='Latitude  range', required=True)
parser.add_argument('--lon-rng', type=float, nargs=2, help='Longitude range. 0-360', required=True)
parser.add_argument('--plot-lon-rng', type=float, nargs=2, help='Longitude range. 0-360', default=[0, 360])
parser.add_argument('--plot-lat-rng', type=float, nargs=2, help='Longitude range. 0-360', default=[-90, 90])
parser.add_argument('--output', type=str, help='Output file.', default="")
args = parser.parse_args()

print(args)

import matplotlib
if args.no_display is False:
    matplotlib.use('TkAgg')
else:
    matplotlib.use('Agg')
    
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.patches as mpatches

import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


print("done")


lat_t = args.lat_rng[0]
lat_b = args.lat_rng[1]

lon_l = args.lon_rng[0]
lon_r = args.lon_rng[1]

dlon = lon_r - lon_l
dlat = lat_t - lat_b

cent_lon = 180.0

proj = ccrs.PlateCarree(central_longitude=cent_lon)
proj_norm = ccrs.PlateCarree()


    
fig, ax = plt.subplots(1, 1, subplot_kw=dict(projection=proj))
ax.set_global()

ax.add_patch(mpatches.Rectangle(xy=[lon_l, lat_b], width=dlon, height=dlat,
                                facecolor='blue',
                                alpha=0.2,
                                transform=proj_norm)
)

ax.gridlines()
ax.coastlines()

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--')
gl.xlabels_top = False
gl.ylabels_left = False
#gl.xlocator = mticker.FixedLocator(np.arange(-180, 181, 60))
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 10, 'color': 'black'}
gl.ylabel_style = {'size': 10, 'color': 'black'}

plot_lat_b = args.plot_lat_rng[0]
plot_lat_t = args.plot_lat_rng[1]

plot_lon_l = args.plot_lon_rng[0]
plot_lon_r = args.plot_lon_rng[1]

ax.set_extent([plot_lon_l, plot_lon_r, plot_lat_b, plot_lat_t], crs=proj_norm)

if not args.no_display:
    plt.show()

if args.output != "":
    
    fig.savefig(args.output, dpi=200)

