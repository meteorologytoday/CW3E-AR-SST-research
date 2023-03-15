import argparse
import numpy as np

parser = argparse.ArgumentParser(
                    prog = 'plot_skill',
                    description = 'Plot prediction skill of GFS on AR.',
)

parser.add_argument('--output', type=str, help='Output file', default="")
parser.add_argument('--no-display', action="store_true")
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


lat_rng = [5, 65]
lon_rng = [115, 245]



lat_t = lat_rng[0]
lat_b = lat_rng[1]

lon_l = lon_rng[0]
lon_r = lon_rng[1]

dlon = lon_r - lon_l
dlat = lat_t - lat_b

cent_lon = 180.0

proj = ccrs.PlateCarree(central_longitude=cent_lon)
proj_norm = ccrs.PlateCarree()

fig, ax = plt.subplots(1, 1, subplot_kw=dict(projection=proj))
ax.set_global()

"""
for b, box in enumerate(boxes):
    
    p = box['polygon']

    lon_l = p['lon_bnds'][0]
    lat_b = p['lat_bnds'][0]

    dlon = p['lon_bnds'][1] - p['lon_bnds'][0]
    dlat = p['lat_bnds'][1] - p['lat_bnds'][0]

    ax.add_patch(mpatches.Rectangle(xy=[lon_l, lat_b], width=dlon, height=dlat,
                                    facecolor='blue',
                                    edgecolor='black',
                                    alpha=0.2,
                                    transform=proj_norm)

    )
    
    mid_lon = np.mean(p['lon_bnds'])
    mid_lat = np.mean(p['lat_bnds'])

    ax.text(mid_lon, mid_lat, "%d" % b, transform=proj_norm, size=12, ha="center", va="center")
    
"""

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
