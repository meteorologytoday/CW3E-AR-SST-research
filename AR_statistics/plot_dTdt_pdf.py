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
parser.add_argument('--IVT-rng', type=float, nargs=2, help='Threshold of IVT to determin AR condition.', default=[250.0, np.inf])
parser.add_argument('--no-display', action="store_true")
args = parser.parse_args()
print(args)


ds = xr.open_dataset(args.input)

MLG_frc = (ds['MLG_frc_sw'] + ds['MLG_frc_lw'] + ds['MLG_frc_sh']  + ds['MLG_frc_lh'] + ds['MLG_frc_fwf']).rename('MLG_frc')
MLG_nonfrc = (ds['MLG_ttl'] - MLG_frc).rename('MLG_nonfrc')
MLG_adv = (ds['MLG_hadv'] + ds['MLG_vadv']).rename('MLG_adv')
MLG_nonadv = (MLG_nonfrc - MLG_adv).rename('MLG_nonadv')
MLG_diff = (ds['MLG_vdiff'] + ds['MLG_hdiff']).rename('MLG_diff')
MLG_nondiff = (MLG_nonfrc - MLG_diff).rename('MLG_nondiff')
MLG_ent = (MLG_nonfrc - MLG_adv - MLG_diff).rename('MLG_ent')
MLG_MLB = (ds['MLG_vdiff'] + MLG_ent).rename('MLG_MLB')

MLG_res = (ds['MLG_ttl'] - MLG_frc - MLG_ent - MLG_diff - MLG_adv).rename('MLG_res')

ds = xr.merge(
    [
        ds,
        MLG_frc,
        MLG_nonfrc,
        MLG_adv,
        MLG_nonadv,
        MLG_diff,
        MLG_ent,
        MLG_nondiff,
        MLG_MLB,
        MLG_res,
    ]
)

var_infos = {

    'MLG_MLB' : {
        'var'  : "$\\dot{T}_\\mathrm{ent} + \\dot{T}_\\mathrm{vdiff}$",
        'unit' : "$ \\mathrm{T} / \\mathrm{s} $",
        'lc' : "green",
    },


    'MLG_ent' : {
        'var'  : "$\\dot{T}_\\mathrm{ent}$",
        'unit' : "$ \\mathrm{T} / \\mathrm{s} $",
    },


    'MLG_nondiff' : {
        'var'  : "$\\dot{T}_\\mathrm{nondiff}$",
        'unit' : "$ \\mathrm{T} / \\mathrm{s} $",
    },

    'MLG_vdiff' : {
        'var'  : "$\\dot{T}_\\mathrm{vdiff}$",
        'unit' : "$ \\mathrm{T} / \\mathrm{s} $",
        'lc' : 'pink'
    },


    'MLG_diff' : {
        'var'  : "$\\dot{T}_\\mathrm{diff}$",
        'unit' : "$ \\mathrm{T} / \\mathrm{s} $",
    },

    'MLG_hdiff' : {
        'var'  : "$\\dot{T}_\\mathrm{hdiff}$",
        'unit' : "$ \\mathrm{T} / \\mathrm{s} $",
    },



    'MLG_adv' : {
        'var'  : "$\\dot{T}_\\mathrm{adv}$",
        'unit' : "$ \\mathrm{T} / \\mathrm{s} $",
        'lc' : "darkorange",
    },

    'MLG_nonadv' : {
        'var'  : "$\\dot{T}_\\mathrm{nonadv}$",
        'unit' : "$ \\mathrm{T} / \\mathrm{s} $",
    },

    'MLG_nonfrc' : {
        'var'  : "$\\dot{T}_\\mathrm{nonfrc}$",
        'unit' : "$ \\mathrm{T} / \\mathrm{s} $",
        'lc' : "blue",
    },



    'MLG_ttl' : {
        'var'  : "$\\dot{T}_\\mathrm{ttl}$",
        'unit' : "$ \\mathrm{T} / \\mathrm{s} $",
        'lc' : "black",
    },

    'MLG_frc' : {
        'var'  : "$\\dot{T}_\\mathrm{shf}$",
        'unit' : "$ \\mathrm{T} / \\mathrm{s} $",
        'lc' : "red",
    },

    'MLG_res' : {
        'var'  : "$\\dot{T}_\\mathrm{res}$",
        'unit' : "$ \\mathrm{T} / \\mathrm{s} $",
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
from scipy.stats import linregress

print("done")


mag = 10.0
bin_edges = np.linspace(-1, 1, 501) * mag
bin_mids = (bin_edges[:-1] + bin_edges[1:])/2
months = [[10, 11, 12], [1, 2, 3]]
titles = ["Oct-Dec", "Jan-Mar"]


rows = 2
cols = len(months)

fig, ax = plt.subplots(rows, cols, sharey='col', figsize=(6*cols, 5*rows), squeeze=False, gridspec_kw = dict(hspace=0.3, wspace=0.4))


fig.suptitle("PDF of T tendency")

ax_flat = ax.flatten()
    

for i, m in enumerate(months):

    if m is None:
        continue

    _ax = ax[:, i]

    _ds = ds.where(ds.time.dt.month.isin(m)).where((ds.IVT >= args.IVT_rng[0]) & (ds.IVT < args.IVT_rng[1]))


    for j, (_bin_edges, varnames) in enumerate([
        (bin_edges, ["MLG_ttl", "MLG_frc", "MLG_nonfrc"]),
        (bin_edges, ["MLG_nonfrc", "MLG_adv", "MLG_vdiff", "MLG_MLB"]),
    ]):
    

        for varname in varnames:
            if np.any((_ds[varname] > _bin_edges[-1]) | (_ds[varname] < _bin_edges[0]) ):
                raise Exception("[Figure %d] Variable %s exceeds bin edges" % (i, varname)) 

            var_info = var_infos[varname]
            _d = _ds[varname] * 1e6
            hist_cnts, _ = np.histogram(_d, bin_edges, density=True)
            #_ax[0].hist(hist_cnt, hist_edges, alpha=0.5, label=varname, histtype='stepfilled', density=True)
            kwargs = {}
            if 'lc' in var_info:
                kwargs['color'] = var_info['lc']

            _ax[j].plot(bin_mids, hist_cnts, label=var_info['var'], **kwargs)
            #_ax[j].hist(_ds[varname], bin_edges, label=varname, histtype='step', fill=False, density=True)
    
        _ax[j].set_title(titles[i])


    unit = "[ $ \\times 10^{-6} \\, \\mathrm{K} \\, / \\, \\mathrm{s} $ ]"
    for __ax in _ax.flatten():
        __ax.legend(fontsize=10)
        __ax.set_xlim(np.array([-3, 1]))
        __ax.set_xlabel(unit)

if args.output != "":
   
    print("Output filename: %s" % (args.output,))
    fig.savefig(args.output, dpi=200)


if not args.no_display:
    print("Show figure")
    plt.show()


