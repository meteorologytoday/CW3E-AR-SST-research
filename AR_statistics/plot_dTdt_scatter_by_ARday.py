import numpy as np
import netCDF4
import AR_tools, NK_tools, fmon_tools, watertime_tools
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
parser.add_argument('--IVT-threshold', type=float, help='Threshold of IVT to determin AR condition.', default=250.0)
#parser.add_argument('--output-database', type=str, help='CSV file.', default="")
parser.add_argument('--no-display', action="store_true")
args = parser.parse_args()
print(args)

ds = xr.open_dataset(args.input)
ds = ds.where(ds.IVT >= args.IVT_threshold)

MLG_frc = (ds['MLG_frc_sw'] + ds['MLG_frc_lw'] + ds['MLG_frc_sh']  + ds['MLG_frc_lh'] + ds['MLG_frc_fwf']).rename('MLG_frc')

ds = xr.merge([ds, MLG_frc])
    
watertime = np.vectorize( lambda d, t: watertime_tools.getWatertime(datetime.combine(d, t)))(ds.time.dt.date, ds.time.dt.time)




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
from scipy.stats import linregress

print("done")


var_infos = {

    'MLG_ttl' : {
        'var'  : "$\\dot{T}_\\mathrm{ttl}$",
        'unit' : "$ \\mathrm{T} / \\mathrm{s} $",
    },

    'MLG_frc' : {
        'var'  : "$\\dot{T}_\\mathrm{shf}$",
        'unit' : "$ \\mathrm{T} / \\mathrm{s} $",
    },

}

def plot_linregress(ax, X, Y, eq_x=0.1, eq_y=0.9, transform=None):


    idx = np.isfinite(X) & np.isfinite(Y)
    
    X = X[idx]
    Y = Y[idx]

    if transform is None:
        transform = ax.transAxes

    res = linregress(X, Y)

    X_min, X_max = np.amin(X), np.amax(X)

    regressline_X = np.array([X_min, X_max])
    regressline_Y = res.slope * regressline_X + res.intercept

    ax.plot(regressline_X, regressline_Y, 'k--')
    
    ax.text(eq_x, eq_y, "$ y = %.2ex + %.2e $\n$R=%.2f$" % (res.slope, res.intercept, res.rvalue), transform=transform, ha="left", va="top", size=8)


    print("Number of data points: %d" % (len(X),))



plot_data = [
    ('MLG_ttl', 'MLG_frc'), 
]

rows = 1


if len(plot_data) % rows != 0:
    cols = len(plot_data) // rows + 1
else:    
    cols = len(plot_data) // rows



color_info = {

    'watertime': {
        'cmap' : 'rainbow',
        'varname'  : 'watertime', 
        'factor'   : 6.0,  # 6 months
        'label' : "Watertime [ mon ]",
        'bnd'   : [0, 6],
    },

    'waterdate': {
        'cmap' : 'rainbow',
        'varname'  : 'waterdate', 
        'factor'   : 12.0,  # 6 months
        'label' : "Watertime [ mon ]",
        'bnd'   : [0, 6],
    },


}['watertime']

fig, ax = plt.subplots(rows, cols, figsize=(6*cols, 5*rows), squeeze=False, gridspec_kw = dict(hspace=0.3, wspace=0.4))

ax_flat = ax.flatten()
    
#fig.suptitle("AR duration time range to do linear regression: %d - %d days\n# of cases = %d" % (args.AR_dt_rng[0], args.AR_dt_rng[1], np.sum(test_data['picked']==True)))

for i, _plot_data in enumerate(plot_data):

    if _plot_data is None:
        continue

    print("Plotting data: ", _plot_data)

    var_info_x = var_infos[_plot_data[0]]
    var_info_y = var_infos[_plot_data[1]]

    _ax = ax_flat[i]
    
    data_needed = {
        'X': _plot_data[0],
        'Y': _plot_data[1],
        'Z': color_info['varname'],
        'picked': 'do_linregress',
    }


    X = ds[_plot_data[0]]
    Y = ds[_plot_data[1]]
    Z = np.mod(watertime, 1)
    
    mappable =  _ax.scatter(X, Y, c = Z * color_info['factor'], s=8, cmap=color_info['cmap'], vmin=color_info['bnd'][0], vmax=color_info['bnd'][1])

    plot_linregress(_ax, X, Y)

    _ax.set_xlabel("%s [ %s ]" % (var_info_x['var'], var_info_x['unit']))
    _ax.set_ylabel("%s [ %s ]" % (var_info_y['var'], var_info_y['unit']))
    
    cbar = plt.colorbar(mappable, ax=_ax, orientation='vertical')
    cbar.ax.set_ylabel(color_info['label'])

for i, _ax in enumerate(ax_flat):
    
    if i > len(plot_data)-1 or plot_data[i] is None: 
        fig.delaxes(_ax)

if not args.no_display:
    print("Show figure")
    plt.show()

if args.output != "":
    print("Output figure: %s" % args.output)
    fig.savefig(args.output, dpi=200)

