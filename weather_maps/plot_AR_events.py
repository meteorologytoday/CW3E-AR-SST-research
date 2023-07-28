import numpy as np
import xarray as xr
import argparse
import pandas as pd
import load_data
from pathlib import Path
import tool_fig_config

parser = argparse.ArgumentParser(
                    prog = 'plot_skill',
                    description = 'Plot prediction skill of GFS on AR.',
)

parser.add_argument('--date-rng', type=str, nargs=2, help='Date range.', required=True)
parser.add_argument('--output-dir', type=str, help='Output dir', default="output_figure")
parser.add_argument('--lat', type=float, nargs=2, help='Latitudes in degree', default=[30, 45])
parser.add_argument('--lon', type=float, nargs=2, help='Longitudes in degree', default=[360-130, 360-120])
parser.add_argument('--no-display', action="store_true")
args = parser.parse_args()
print(args)


dts = pd.date_range(args.date_rng[0], args.date_rng[1], freq="D", inclusive="both")

args.lon = np.array(args.lon) % 360.0

lat_n, lat_s = np.amax(args.lat), np.amin(args.lat)
lon_w, lon_e = np.amin(args.lon), np.amax(args.lon)


print("==================================")
print("Date range: ", dts[0], " to ", dts[-1])
print("Latitude  box: %.2f %.2f" % (lat_s, lat_n))
print("Longitude box: %.2f %.2f" % (lon_w, lon_e))
print("==================================")





import matplotlib as mpl
if args.no_display is False:
    mpl.use('TkAgg')
else:
    mpl.use('Agg')
    #mpl.rc('font', size=20)
 
# This program plots the AR event
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Rectangle
import matplotlib.transforms as transforms
from matplotlib.dates import DateFormatter
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

print("done")

print("Create dir: %s" % (args.output_dir,))
Path(args.output_dir).mkdir(parents=True, exist_ok=True)

for dt in dts:

    print("Doing date: ", dt)

    output_filename = "%s/weather_%s.png" % (args.output_dir, dt.strftime("%Y-%m-%d"))
    
    _data = []
    
    for varname in ["sst", "u10", "v10", "t2m", "IWV", "IVT", "mtpr"]:
        info = load_data.getFileAndIndex(product="ERA5", date=dt, varname=varname)
        _data.append(xr.load_dataset(info['filename'])[varname])
        _time = _data[-1].coords["time"]
 
    for varname in ["map",]:
        info = load_data.getFileAndIndex(product="ERA5_ARobj", date=dt, varname=varname, method="ANOM_LEN")

        # these two datasets do not have a uniform time for the same day
        _tmp = xr.load_dataset(info['filename'])[varname]
        _tmp = _tmp.assign_coords(time=_time)
        
        _data.append(_tmp)
 
    _data = xr.merge(_data).isel(time=0)

    cent_lon = 180.0

    plot_lon_l = 360 - 130.0
    plot_lon_r = 360 - 120.0
    plot_lat_b = 30.0
    plot_lat_t = 45.0

    proj = ccrs.PlateCarree(central_longitude=cent_lon)
    proj_norm = ccrs.PlateCarree()

    figsize, gridspec_kw = tool_fig_config.calFigParams(
        w = 5,
        h = 5,
        wspace = 1.5,
        hspace = 0.5,
        w_left = 1.0,
        w_right = 1.5,
        h_bottom = 1.0,
        h_top = 1.0,
        ncol = 2,
        nrow = 1,
    )

    fig, ax = plt.subplots(
        1, 2,
        figsize=figsize,
        subplot_kw=dict(projection=proj, aspect="auto"),
        gridspec_kw=gridspec_kw,
        constrained_layout=False,
        squeeze=False,
    )

    coords = _data.coords
    #cmap = cm.get_cmap("bwr")
    #cmap.set_over("green")
    #cmap.set_under("yellow")

    fig.suptitle(dt.strftime("%Y-%m-%d"))

    _ax = ax[0, 0]
    mappable = _ax.contourf(coords["lon"], coords["lat"], _data.IVT, levels=np.linspace(0, 600, 13), transform=proj_norm, extend="max", cmap="GnBu")
    
    _ax.quiver(coords["lon"], coords["lat"], _data.u10.to_numpy(), _data.v10.to_numpy(), scale=200, transform=proj_norm)
 
    cs = _ax.contourf(coords["lon"], coords["lat"], _data['map'], colors='none', levels=[0, 0.5, np.inf], hatches=[None, "."], transform=proj_norm)

    # Remove the contour lines for hatches 
    for _, collection in enumerate(cs.collections):
        collection.set_edgecolor("red")
     
    cax = tool_fig_config.addAxesNextToAxes(fig, _ax, "right", thickness=0.03, spacing=0.05)
    cb = plt.colorbar(mappable, cax=cax, orientation="vertical", pad=0.0)
    cb.ax.set_ylabel(" IVT [ $ \\mathrm{kg} \\, / \\, \\mathrm{m} \\, / \\, \\mathrm{s} $ ]")
 
    _ax = ax[0, 1]
    mappable = _ax.contourf(coords["lon"], coords["lat"], _data.sst - 273.15, levels=np.linspace(10, 20, 21), transform=proj_norm, extend="both", cmap="rainbow")
   
    print(np.nanmean(_data.mtpr * 86400.0)) 
    cs = _ax.contourf(coords["lon"], coords["lat"], _data.mtpr * 86400.0, colors='none', levels=[0, 10, 2e5], hatches=[None, "//"], transform=proj_norm)

    # Remove the contour lines for hatches 
    for _, collection in enumerate(cs.collections):
        collection.set_edgecolor("yellow")
        #collection.set_linewidth(0.)


    cax = tool_fig_config.addAxesNextToAxes(fig, _ax, "right", thickness=0.03, spacing=0.05)
    cb = plt.colorbar(mappable, cax=cax, orientation="vertical", pad=0.0)
    cb.ax.set_ylabel(" SST [ $ \\mathrm{C} $ ]")
    
           
    for __ax in ax.flatten(): 

        __ax.set_global()
        #__ax.gridlines()
        __ax.coastlines(color='gray')
        __ax.set_extent([plot_lon_l, plot_lon_r, plot_lat_b, plot_lat_t], crs=proj_norm)

        gl = __ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=1, color='gray', alpha=0.5, linestyle='--')

        gl.xlabels_top   = False
        gl.ylabels_right = False

        #gl.xlocator = mticker.FixedLocator(np.arange(-180, 181, 30))
        #gl.xlocator = mticker.FixedLocator([120, 150, 180, -150, -120])#np.arange(-180, 181, 30))
        #gl.ylocator = mticker.FixedLocator([10, 20, 30, 40, 50])
        
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 10, 'color': 'black'}
        gl.ylabel_style = {'size': 10, 'color': 'black'}


    print("Output file: ", output_filename)
    
    fig.savefig(output_filename, dpi=200)
    if not args.no_display:
        plt.show()
