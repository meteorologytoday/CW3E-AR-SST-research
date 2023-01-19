import numpy as np
import netCDF4
import argparse

parser = argparse.ArgumentParser(description='Description')
parser.add_argument('--input-file', type=str, )
parser.add_argument('--no-display', action="store_true")
parser.add_argument('--scale-days', type=int, default=1)
parser.add_argument('--output', type=str, default="")
parser.add_argument('--title', type=str, default="")

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

if args.no_display is True:  # Using TkAgg
    matplotlib.rc('axes', labelsize=25)
    matplotlib.rc('font', size=25)



print("done")


    
scaled_varnames = [
    "TOTTTEND_mean", "TEND_SFCFLX", "TEND_VDIFF", "TEND_HDIFF",
    "TEND_ENT_dhdt", "TEND_RES", "TEND_SUM",
    "TEND_TRANSPORT", "TEND_SW", "TEND_CORRFLX",
]

nonscaled_varnames = [
    "EXFuwind", "EXFvwind", "oceFWflx",
]


data = {}
lat = None
lon = None
# load file
with netCDF4.Dataset(args.input_file, "r") as ds:

    lat = ds.variables["lat"][:].flatten() 
    lon = ds.variables["lon"][:].flatten()
    
    for varname in scaled_varnames:
        print("Loading %s" % (varname,))
        data[varname] = ds.variables[varname][0, :, :] * 86400.0 * args.scale_days

    for varname in nonscaled_varnames:
        print("Loading %s" % (varname,))
        data[varname] = ds.variables[varname][0, :, :]


    data["TEND_SFCFLX_ALL"]    = data["TEND_SW"] + data["TEND_SFCFLX"]
    data["TEND_TRANSPORT_ALL"] = data["TEND_TRANSPORT"] + data["TEND_CORRFLX"]
    data["TEND_BOTTOM_ALL"]    = data["TEND_ENT_dhdt"]  + data["TEND_VDIFF"]
    data["wind"] = np.sqrt(data["EXFuwind"]**2 + data["EXFvwind"]**2)

threshold=0.000 
PmE = data["oceFWflx"] * 86400.0 / 1000.0 * 1000 # mm / day
PmE[PmE <= threshold]  = -1
PmE[PmE > threshold] = 1

print("Sum of PmE: ", np.sum(PmE == 1))

cent_lon = 180.0

proj = ccrs.PlateCarree(central_longitude=cent_lon)
proj_norm = ccrs.PlateCarree()

plot_data = [
    "weather",  "TOTTTEND_mean",  "TEND_SUM",       "TEND_RES",
    "TEND_SW",        "TEND_SFCFLX",    "TEND_SFCFLX_ALL", None, 
    "TEND_TRANSPORT", "TEND_CORRFLX",   "TEND_TRANSPORT_ALL", None,
    "TEND_ENT_dhdt",  "TEND_VDIFF",     "TEND_BOTTOM_ALL",  None,
]
rows = 4

if len(plot_data) % rows != 0:
    cols = len(plot_data) // rows + 1
else:    
    cols = len(plot_data) // rows



fig, ax = plt.subplots(rows, cols, figsize=(6*cols, 5*rows),
    squeeze=False,
    gridspec_kw = dict(hspace=0.3, wspace=0.4),
    subplot_kw=dict(projection=proj)
)

if args.title == "":
    fig.suptitle(args.input_file)
else:
    fig.suptitle(args.title)


for _ax in ax.flatten():
    _ax.set_global()
    _ax.coastlines()
    gl = _ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    #gl.xlocator = mticker.FixedLocator(np.arange(-180, 181, 60))
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}

    _ax.set_extent([-132, -118, 30, 44], crs=proj_norm)
    #_ax.gridlines()
    _ax.coastlines()



levs = np.linspace(-1, 1, 21) * 1.0

cmap = plt.cm.get_cmap("bwr")
cmap.set_under("yellow")
cmap.set_over("green")

for i, varname in enumerate(plot_data): 
    
    _ax = ax.flatten()[i]
    
    if varname is None:
        fig.delaxes(_ax)
        continue


    
    if varname == "weather":
        
        #mappable_wind = _ax.contourf(lon, lat, data["wind"], np.linspace(0, 20, 11), cmap="hot_r", transform=proj_norm, extend="max")
        #mappable_wind = _ax.contourf(lon, lat, data["oceFWflx"] * 86400, np.linspace(-5, 5.1, 0.5), cmap="hot_r", transform=proj_norm, extend="max")
        #_ax.contour(lon, lat, data["oceFWflx"]*86400, np.linspace(-5, 5, 11), colors="black", transform=proj_norm)
       

        mappable_wind = _ax.contourf(lon, lat, PmE, [0.5, 1.1], colors=["red",], alpha=0.25, transform=proj_norm)
        #_ax.contour(lon, lat, data["wind"], np.linspace(0, 20, 11), colors="black", transform=proj_norm)

    else:

        mappable = _ax.contourf(lon, lat, data[varname], levs, cmap=cmap, transform=proj_norm, extend="both")


    _ax.set_title(varname)

    
cbar = plt.colorbar(mappable, ax=ax.ravel().tolist(), orientation="vertical")
cbar.ax.set_ylabel("[ degC / %d days ]" % (args.scale_days,))

cbar_wind = plt.colorbar(mappable_wind, ax=ax[0, 0], orientation="vertical")
#cbar_wind.ax.set_ylabel("[ m / s ]")

if not args.no_display:
    plt.show()

if args.output != "":
    
    fig.savefig(args.output, dpi=200)

