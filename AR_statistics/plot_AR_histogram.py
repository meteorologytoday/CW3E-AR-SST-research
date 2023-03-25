from multiprocessing import Pool
import numpy as np
import os.path as path
#import fmon_tools, watertime_tools
import anomalies
import ARstat_tool
import xarray as xr
import watertime_tools
import traceback
from pathlib import Path
import argparse
from datetime import (timedelta, datetime, timezone)

parser = argparse.ArgumentParser(
                    prog = 'plot_skill',
                    description = 'Plot prediction skill of GFS on AR.',
)

parser.add_argument('--climanom-input-dir', type=str, help='Input file', required=True)
parser.add_argument('--input-dir', type=str, help='Input file', required=True)
parser.add_argument('--output', type=str, help='Input file', default="")
parser.add_argument('--beg-year', type=int, help='Input file', required=True)
parser.add_argument('--end-year', type=int, help='Input file', required=True)
parser.add_argument('--ncpu', type=int, help='Number of CPUs.', default=4)
parser.add_argument('--overwrite', help='If we overwrite the output', action="store_true")
parser.add_argument('--IWV-cond-off', action="store_true")

parser.add_argument('--title', type=str, help='Output title', default="")
parser.add_argument('--no-display', action="store_true")
args = parser.parse_args()
print(args)

yrs = list(range(args.beg_year, args.end_year+1))
ds = ARstat_tool.loadDatasets(args.input_dir, yrs)

filename_clim = "%s/clim.nc" % (args.climanom_input_dir,)
filename_anom = "%s/anom.nc" % (args.climanom_input_dir,)

ds_clim = xr.open_dataset(filename_clim)
ds_anom = xr.open_dataset(filename_anom)


# Construct
t_months = np.arange(1, 7)

ds_stats = {}

print("Constructing AR condition")
AR_cond = (ds.IVT >= 250) & (ds.IWV >= 20.0)
time_cond = ds.time.dt.month.isin(watertime_tools.wm2m(t_months))

spatial_cond = ds.coords["lat"] > 20.0

total_cond = AR_cond & time_cond & spatial_cond

ds_anom = ds_anom.where(total_cond)


nethf_anom = (ds_anom.mslhf + ds_anom.msshf + ds_anom.msnswrf + ds_anom.msnlwrf).rename("nethf")

ds_anom = xr.merge([ds_anom, nethf_anom])


def doHist(X, Y, X_edges, Y_edges):
    
    finite_idx = np.isfinite(Y) & np.isfinite(X)
    
    Y = Y[finite_idx]
    X = X[finite_idx]

    hist, X_edges, Y_edges = np.histogram2d(X, Y, bins=[X_edges, Y_edges], density=True)

    X_mids = (X_edges[:-1] + X_edges[1:]) / 2
    Y_mids = (Y_edges[:-1] + Y_edges[1:]) / 2

    return dict(hist=hist, X_edges=X_edges, Y_edges=Y_edges, X_mids=X_mids, Y_mids=Y_mids)



hist_pairs = [
    ("MLG_nonfrc", "u10"),
    ("MLG_nonfrc", "mtpr"),
    ("MLG_nonfrc", "nethf"),
]

edges = {
    "MLG_nonfrc" : np.linspace(-1, 1, 201) * 1e-6,
    "u10"        : np.linspace(-1, 1, 201) * 20,
    "mtpr"       : np.linspace(-1, 1, 201) * 1e-4,
    "nethf"      : np.linspace(-1, 1, 201) * 400,
}


data = []

print("Doing histogram2d...")
for i, (Y_varname, X_varname) in enumerate(hist_pairs):
    
    Y = ds_anom[Y_varname].to_numpy().flatten()
    X = ds_anom[X_varname].to_numpy().flatten()

    X_edges = edges[X_varname]
    Y_edges = edges[Y_varname]
    
    result = doHist(X, Y, X_edges, Y_edges)

    result["varnames"] = dict(X=X_varname, Y=Y_varname)

    data.append(result)


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

print("done")



fig, ax = plt.subplots(1, len(data), figsize=(4*len(data), 4))

for i, _data in enumerate(data):
    _ax = ax[i] 
    _hist = _data["hist"].transpose()
    _hist[_hist == 0] = np.nan
    mappable = _ax.contourf(_data["X_mids"], _data["Y_mids"], _hist, cmap="hot_r")
    #cb = plt.colorbar(mappable, ax=_ax)
    _ax.set_title("(%s) %s v.s. %s" % ("abcdefg"[i], _data["varnames"]["X"], _data["varnames"]["Y"]))

if not args.no_display:
    plt.show()

if args.output != "":
    
    fig.savefig(args.output, dpi=200)

