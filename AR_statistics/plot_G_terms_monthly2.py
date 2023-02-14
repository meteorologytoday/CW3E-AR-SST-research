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

#MLG_frc = (ds['MLG_frc_sw'] + ds['MLG_frc_lw'] + ds['MLG_frc_sh']  + ds['MLG_frc_lh'] + ds['MLG_frc_fwf']).rename('MLG_frc')
MLG_frc = (ds['MLG_frc_sw'] + ds['MLG_frc_lw'] + ds['MLG_frc_sh']  + ds['MLG_frc_lh']).rename('MLG_frc')
MLG_nonfrc = (ds['MLG_ttl'] - MLG_frc).rename('MLG_nonfrc')
MLG_adv = (ds['MLG_hadv'] + ds['MLG_vadv']).rename('MLG_adv')
MLG_nonadv = (MLG_nonfrc - MLG_adv).rename('MLG_nonadv')
MLG_diff = (ds['MLG_vdiff'] + ds['MLG_hdiff']).rename('MLG_diff')
MLG_nondiff = (MLG_nonfrc - MLG_diff).rename('MLG_nondiff')
#MLG_MLB = (ds['MLG_vdiff'] + ds["MLG_ent"]).rename('MLG_MLB')

#MLG_res = (ds['MLG_ttl'] - MLG_frc - ds["MLG_ent"] - MLG_diff - MLG_adv).rename('MLG_res')
MLG_res2 = (ds['MLG_ttl'] - MLG_frc - ds["MLG_vdiff"] - MLG_adv).rename('MLG_res2')

ds = xr.merge(
    [
        ds,
        MLG_frc,
        MLG_nonfrc,
        MLG_adv,
        MLG_nonadv,
        MLG_diff,
        MLG_nondiff,
#        MLG_MLB,
#        MLG_res,
        MLG_res2,
    ]
)

# Construct
t_months = np.arange(1, 7)
ds_stat = xr.Dataset(

    {
        "MLG_ttl"   : (["time", "stat"], np.zeros((len(t_months), 5)) ),
        "MLG_frc"   : (["time", "stat"], np.zeros((len(t_months), 5)) ),
        "MLG_nonfrc"   : (["time", "stat"], np.zeros((len(t_months), 5)) ),
        "MLG_adv"   : (["time", "stat"], np.zeros((len(t_months), 5)) ),
        "MLG_vdiff" : (["time", "stat"], np.zeros((len(t_months), 5)) ),
        "MLG_res2"   : (["time", "stat"], np.zeros((len(t_months), 5)) ),
    },

    coords = {
        "time" : t_months,
        "stat" : ["mean", "q25", "q75", "q10", "q90"]
    }
)

for t, m in enumerate(t_months): 
    
    print("m = ", m, "; wm = ", watertime_tools.wm2m(m))
    
    _ds = ds.where(ds.time.dt.month.isin(watertime_tools.wm2m(m))).where((ds.IVT >= args.IVT_rng[0]) & (ds.IVT < args.IVT_rng[1]) & (ds.MLG_ttl < 5e-6))
        
    for varname, _ in ds_stat.items():

        #print(var[np.isfinite(var)])
        #print(np.nanquantile(var, [0.0, .25, .5, .75, 1.0]))
        _data = _ds[varname] * 1e6
        ds_stat[varname][t, 0] = np.nanmean(_data)
        ds_stat[varname][t, 1] = np.nanquantile(_data, .25)
        ds_stat[varname][t, 2] = np.nanquantile(_data, .75)
        ds_stat[varname][t, 3] = np.nanquantile(_data, .10)
        ds_stat[varname][t, 4] = np.nanquantile(_data, .90)
        
print(ds_stat)


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


fig, ax = plt.subplots(1, 1, figsize=(6, 4), squeeze=False, gridspec_kw = dict(hspace=0.3, wspace=0.4))


fig.suptitle(args.input)

bar_width = 0.15

#for i, varname in enumerate(["MLG_ttl", "MLG_frc", "MLG_adv", "MLG_vdiff", "MLG_res2"]):
for i, varname in enumerate(["MLG_ttl", "MLG_frc", "MLG_nonfrc"]):

    _data = ds_stat[varname]
    ax[0, 0].bar(t_months + i*bar_width, _data[:, 0],   bar_width, label=varname)

    for m, t_month in enumerate(t_months): 
        ax[0, 0].plot([t_month + i*bar_width]*2, [_data[m, 1], _data[m, 2]], "k-") 
        ax[0, 0].scatter([t_month + i*bar_width]*2, [_data[m, 3], _data[m, 4]], s=10, c="red") 

ax[0, 0].set_xticks(t_months)
ax[0, 0].set_xticklabels(np.vectorize(watertime_tools.wm2m)(t_months))

ax[0, 0].set_xlabel("Month")
ax[0, 0].set_ylabel("[ $ 1 \\times 10^{-6} \\mathrm{K} \\, / \\, \\mathrm{s} $ ]")

if args.output != "":
   
    print("Output filename: %s" % (args.output,))
    fig.savefig(args.output, dpi=200)


if not args.no_display:
    print("Show figure")
    plt.show()


