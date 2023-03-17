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



plot_infos = {

    'IWV' : {
        'var'  : "$ \\mathrm{IWV}$",
        'unit' : "$ \\mathrm{kg} / \\mathrm{m}^2 $",
        'levs' : np.linspace(20, 50, 6),
        'cmap' : "bone_r",
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
        'cmap' : "hot_r",
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
    ("IWV", "IVT"),
    ("MEAN_VEL", None),
]

fig, ax = plt.subplots(len(plot_vars), 1, sharex=True, figsize=(5, len(plot_vars) * 3), squeeze=False, gridspec_kw = dict(hspace=0.3, wspace=0.4))


for (i, (ctrf_varname, ctr_varname)) in enumerate(plot_vars):
    
    _ax = ax[i, 0]
    
    title = "(%s)" % ("abcdefghijklmn"[i],)

    if ctrf_varname is not None:

        pinfo = plot_infos[ctrf_varname]
        
        title += " %s (shading)" % pinfo["var"]

        mappable = _ax.contourf(ds.coords['lat'], ds.coords['watermonth'], ds[ctrf_varname], levels=pinfo["levs"], cmap=pinfo["cmap"], extend="max")
        cb = plt.colorbar(mappable, ax=_ax, ticks=pinfo["levs"], orientation="vertical")
        cb.set_label("%s [ %s ]" % (pinfo["var"], pinfo["unit"]))

    if ctr_varname is not None:

        print("ctr_varname: ", ctr_varname)

        pinfo = plot_infos[ctr_varname]
        
        title += ", %s (contour)" % pinfo["var"]
        cs = _ax.contour(ds.coords['lat'], ds.coords['watermonth'], ds[ctr_varname], levels=pinfo["levs"], colors="k")
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


