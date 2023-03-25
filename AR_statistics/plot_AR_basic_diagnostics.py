import numpy as np
import netCDF4
import fmon_tools, watertime_tools
import anomalies
import pandas as pd

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
parser.add_argument('--output', type=str, help='Output file', default="")
parser.add_argument('--no-display', action="store_true")
args = parser.parse_args()
print(args)


ds = xr.open_dataset(args.input)

ds = ds.mean(dim="lon", skipna=True)

ds['mtpr'] *= 3600 * 24


plot_infos = {

    'IWV' : {
        'var'  : "$ \\mathrm{IWV}$",
        'unit' : "$ \\mathrm{kg} / \\mathrm{m}^2 $",
        'levs' : np.linspace(20, 50, 7),
        'cmap' : "GnBu",
    },

    'IVT' : {
        'var'  : "$ \\mathrm{IVT}$",
        'unit' : "$ \\mathrm{kg} / \\mathrm{m} / \\mathrm{s} $",
        'levs' : np.linspace(250, 500, 6),
        'cmap' : "hot_r",
        'fmt'  : "%d",
    },

    'MEAN_VEL' : {
        'var'  : "$ \\overline{V}$",
        'unit' : "$ \\mathrm{m} / \\mathrm{s} $",
        'levs' : np.linspace(0, 30, 7),
        'anomratio_levs' : np.linspace(-2, 2, 11),
        'cmap' : "hot_r",
        'fmt'  : "%.1f",
        'extend' : "both",
    },

    'absU10' : {
        'var'  : "U10",
        'unit' : "$ \\mathrm{m} / \\mathrm{s} $",
        'levs' : np.linspace(0, 30, 7),
        'anomratio_levs' : np.linspace(-1, 1, 11),
        'cmap' : "bwr",
        'fmt'  : "%.1f",
        'extend' : "both",
    },


    'u10' : {
        'var'  : "Zonal 10m wind",
        'unit' : "$ \\mathrm{m} / \\mathrm{s} $",
        'levs' : np.linspace(0, 30, 7),
        'anom_levs' : np.linspace(-5, 5, 11),
        'anomratio_levs' : np.linspace(-1, 1, 11),
        'cmap' : "bwr",
        'fmt'  : "%.1f",
        'extend' : "both",
    },

    'v10' : {
        'var'  : "Meridional 10m wind",
        'unit' : "$ \\mathrm{m} / \\mathrm{s} $",
        'levs' : np.linspace(0, 30, 7),
        'anom_levs' : np.linspace(-5, 5, 11),
        'anomratio_levs' : np.linspace(-1, 1, 11),
        'cmap' : "bwr",
        'fmt'  : "%.1f",
        'extend' : "both",
    },




    'mtpr' : {
        'var'  : "Total precipitation",
        'unit' : "$ \\mathrm{mm} / \\mathrm{day} $",
        'levs' : np.linspace(0, 30, 7),
        'anom_levs' : np.linspace(0, 20, 11),
        'extend' : "max",
        'cmap' : "hot_r",
        'fmt'  : "%.1f",
    },




}

# Plot data
print("Loading Matplotlib...")
import matplotlib as mpl
if args.no_display is False:
    mpl.use('TkAgg')
else:
    mpl.use('Agg')
    mpl.rc('font', size=20)
    mpl.rc('axes', labelsize=15)
     
 
  
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Rectangle
import matplotlib.transforms as transforms
from matplotlib.dates import DateFormatter

print("done")


plot_vars = [
    ("clim", "IWV", "IVT"),
    ("anom",  "u10", "v10"),
    ("anom",  "mtpr", None),
    ("anomratio",  "absU10", "MEAN_VEL"),
]

fig, ax = plt.subplots(len(plot_vars), 1, sharex=True, figsize=(5, len(plot_vars) * 3), squeeze=False, gridspec_kw = dict(hspace=0.3, wspace=0.4))


for (i, (measure, ctrf_varname, ctr_varname)) in enumerate(plot_vars):
    
    _ax = ax[i, 0]
    
    title = "(%s)" % ("abcdefghijklmn"[i],)

    if measure == "clim":

        if ctrf_varname is not None:

            pinfo = plot_infos[ctrf_varname]
            
            title += " %s (shading)" % pinfo["var"]
            
            print(ds[ctrf_varname])

            mappable = _ax.contourf(ds.coords['lat'], ds.coords['watermonth'], ds[ctrf_varname].sel(stat="AR"), levels=pinfo["levs"], cmap=pinfo["cmap"], extend="max")
            cb = plt.colorbar(mappable, ax=_ax, ticks=pinfo["levs"], orientation="vertical")
            cb.set_label("%s [ %s ]" % (pinfo["var"], pinfo["unit"]))

        if ctr_varname is not None:

            print("ctr_varname: ", ctr_varname)

            pinfo = plot_infos[ctr_varname]
            
            title += ", %s (contour)" % pinfo["var"]
            cs = _ax.contour(ds.coords['lat'], ds.coords['watermonth'], ds[ctr_varname].sel(stat="AR"), levels=pinfo["levs"], colors="k")
            _ax.clabel(cs, fmt=pinfo["fmt"])

    elif measure == "anom":

        if ctrf_varname is not None:

            pinfo = plot_infos[ctrf_varname]
            
            title += " %s (shading)" % pinfo["var"]

            data_anom = (ds[ctrf_varname].sel(stat="AR") - ds[ctrf_varname].sel(stat="clim"))

            mappable = _ax.contourf(ds.coords['lat'], ds.coords['watermonth'], data_anom, levels=pinfo["anom_levs"], cmap=pinfo["cmap"], extend=pinfo["extend"])
            cb = plt.colorbar(mappable, ax=_ax, ticks=pinfo["anom_levs"], orientation="vertical")
            cb.set_label("%s anomaly [ %s ]" % (pinfo["var"], pinfo["unit"]))

        if ctr_varname is not None:

            print("ctr_varname: ", ctr_varname)

            pinfo = plot_infos[ctr_varname]
            
            title += ", %s (contour)" % pinfo["var"]



            data_anom = (ds[ctr_varname].sel(stat="AR") - ds[ctr_varname].sel(stat="clim"))

            cs = _ax.contour(ds.coords['lat'], ds.coords['watermonth'], data_anom, levels=pinfo["anom_levs"], colors="k")
            _ax.clabel(cs, fmt=pinfo["fmt"])


    elif measure == "anomratio":

        if ctrf_varname is not None:

            pinfo = plot_infos[ctrf_varname]
            
            title += " %s (shading)" % pinfo["var"]

            data_anomratio = (ds[ctrf_varname].sel(stat="AR") - ds[ctrf_varname].sel(stat="clim")) / ds[ctrf_varname].sel(stat="std")

            mappable = _ax.contourf(ds.coords['lat'], ds.coords['watermonth'], data_anomratio, levels=pinfo["anomratio_levs"], cmap=pinfo["cmap"], extend=pinfo["extend"])
            cb = plt.colorbar(mappable, ax=_ax, ticks=pinfo["anomratio_levs"], orientation="vertical")
            cb.set_label("%s anomaly normalized by std" % (pinfo["var"],))

        if ctr_varname is not None:

            print("ctr_varname: ", ctr_varname)

            pinfo = plot_infos[ctr_varname]
            
            title += ", %s (contour)" % pinfo["var"]



            data_anomratio = (ds[ctr_varname].sel(stat="AR") - ds[ctr_varname].sel(stat="clim")) / ds[ctrf_varname].sel(stat="std")

            cs = _ax.contour(ds.coords['lat'], ds.coords['watermonth'], data_anomratio, levels=pinfo["anomratio_levs"], colors="k")
            _ax.clabel(cs, fmt=pinfo["fmt"])


    _ax.set_title(title)
    _ax.set_xticks(np.linspace(10, 60, 6)) 
    _ax.set_yticks(np.arange(1, 1+len(ds.coords['watermonth'])))    
    
    _ax.set_xlabel("Latitude [ deg ]")
 
if args.output != "":
   
    print("Output filename: %s" % (args.output,))
    fig.savefig(args.output, dpi=200)


if not args.no_display:
    print("Show figure")
    plt.show()


