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
        'levs' : np.linspace(0, 50, 6),
        'cmap' : "rainbow",
    },

    'IVT' : {
        'var'  : "$ \\mathrm{IVT}$",
        'unit' : "$ \\mathrm{kg} / \\mathrm{m} / \\mathrm{s} $",
        'levs' : np.linspace(0, 500, 6),
        'cmap' : "hot_r",
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


plot_vars = ["IWV", "IVT", "MEAN_VEL"]

fig, ax = plt.subplots(len(plot_vars), 1, sharex=True, figsize=(5, len(plot_vars) * 3), squeeze=False, gridspec_kw = dict(hspace=0.3, wspace=0.4))



for (i, varname) in enumerate(plot_vars):
    
    _ax = ax[i]

    pinfo = plot_infos[varname]

    _ax.set_title("(%s) %s" % ("abcdefghijklmn"[i], pinfo["var"])

    mappable = _ax.contourf(ds.coord['lat'], ds.coord['watermonth'], ds[varname], levels=pinfo["levs"], cmap=pinfo["cmap"])
    
    cb = plt.colorbar(mappable, ax=_ax, levels=pinfo["levs"], orientation="vertical")
    cb.set_label("%s [ %s ]" % (pinfo["var"], pinfo["unit"])
    
    
if args.output != "":
   
    print("Output filename: %s" % (args.output,))
    fig.savefig(args.output, dpi=200)


if not args.no_display:
    print("Show figure")
    plt.show()


