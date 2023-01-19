import numpy as np
import netCDF4
import argparse

parser = argparse.ArgumentParser(description='Description')
parser.add_argument('--input-file', type=str, )

args = parser.parse_args()
print(args)

import matplotlib
matplotlib.use('TkAgg')
    
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.patches as mpatches

import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

#matplotlib.rc('axes', labelsize=25)
#matplotlib.rc('font', size=25)



print("done")

varnames = [
    "TEND_SFCFLX", "TEND_VDIFF", "TEND_HDIFF",
    "TEND_RES", "TEND_SUM",
    "TEND_TRANSPORT", "TEND_SW", "TEND_CORRFLX",
]

data = {}
lat = None
lon = None
# load file
with netCDF4.Dataset(args.input_file, "r") as ds:

    lat = ds.variables["lat"][:].flatten() 
    lon = ds.variables["lon"][:].flatten()
    
    for varname in varnames:
        print("Loading %s" % (varname,))
        data[varname] = ds.variables["%s_3D" % (varname,)][0, 0, :, :] 


ratio = data["TEND_SFCFLX"] / data["TEND_RES"]

ratio[np.isnan(ratio)] = np.nan

print("Max  : %f" % (np.amax(ratio),))
print("Min  : %f" % (np.amin(ratio),))
print("Std  : %f" % (np.mean(ratio),))
print("mean : %f" % (np.std(ratio),))


#for r in ratio[:]:
#    print(r)

cent_lon = 180.0

proj = ccrs.PlateCarree(central_longitude=cent_lon)
proj_norm = ccrs.PlateCarree()

fig, ax = plt.subplots(2, 3, figsize=(5*3, 5*2),
    squeeze=False,
    gridspec_kw = dict(hspace=0.3, wspace=0.4),
    subplot_kw=dict(projection=proj)
)


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
    _ax.coastlines()


levs=np.linspace(-1, 1, 21)
cmap = plt.cm.get_cmap("bwr")
cmap.set_under("yellow")
cmap.set_over("green")


ax[0, 0].set_title("TEND_RES")
mappable = ax[0, 0].contourf(lon, lat, data["TEND_RES"], levs * 2e-6, cmap=cmap, transform=proj_norm, extend="both")
cbar = plt.colorbar(mappable, ax=ax[0, 0], orientation="vertical")

ax[0, 1].set_title("TEND_SFCFLX")
mappable = ax[0, 1].contourf(lon, lat, data["TEND_SFCFLX"], levs * 1e-4, cmap=cmap, transform=proj_norm, extend="both")
cbar = plt.colorbar(mappable, ax=ax[0, 1], orientation="vertical")

ax[0, 2].set_title("TEND_VDIFF")
mappable = ax[0, 2].contourf(lon, lat, data["TEND_VDIFF"], levs * 1e-4, cmap=cmap, transform=proj_norm, extend="both")
cbar = plt.colorbar(mappable, ax=ax[0, 2], orientation="vertical")

ax[1, 0].set_title("TEND_CORRFLX")
mappable = ax[1, 0].contourf(lon, lat, data["TEND_CORRFLX"], levs * 1e-5, cmap=cmap, transform=proj_norm, extend="both")
cbar = plt.colorbar(mappable, ax=ax[1, 0], orientation="vertical")


ax[1, 1].set_title("TEND_SFCFLX / TEND_RES")
mappable = ax[1, 1].contourf(lon, lat, ratio, np.linspace(30, 80, 101), cmap=cmap, transform=proj_norm, extend="both")
cbar = plt.colorbar(mappable, ax=ax[1, 1], orientation="vertical")


ax[1, 2].set_title("- TEND_VDIFF / TEND_RES")
mappable = ax[1, 2].contourf(lon, lat, - data["TEND_VDIFF"] / data["TEND_RES"], np.linspace(10, 50, 101), cmap=cmap, transform=proj_norm, extend="both")
cbar = plt.colorbar(mappable, ax=ax[1, 2], orientation="vertical")



#cbar.ax.set_ylabel("[ degC / %d days ]" % (args.scale_days,))

plt.show()
