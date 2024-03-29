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
parser.add_argument('--IVT-rng', type=float, nargs=2, help='Threshold of IVT to determin AR condition.', default=[250.0, np.inf])
parser.add_argument('--watermonths', type=float, nargs='+', help='The months you need.', default=[1, 2, 3, 4, 5, 6])
#parser.add_argument('--output-database', type=str, help='CSV file.', default="")
parser.add_argument('--no-display', action="store_true")
args = parser.parse_args()
print(args)

ds = xr.open_dataset(args.input)
ds = ds.where((ds.IVT >= args.IVT_rng[0]) & (ds.IVT < args.IVT_rng[1]))

MLG_frc = (ds['MLG_frc_sw'] + ds['MLG_frc_lw'] + ds['MLG_frc_sh']  + ds['MLG_frc_lh'] + ds['MLG_frc_fwf'] + ds['MLG_rescale']).rename('MLG_frc')

MLG_adv = (ds['MLG_hadv'] + ds['MLG_vadv']).rename('MLG_adv')
MLG_diff = (ds['MLG_vdiff'] + ds['MLG_hdiff']).rename('MLG_diff')
MLG_nonfrc = (MLG_adv + MLG_diff + ds['MLG_ent']).rename('MLG_nonfrc')

MLG_nonfrc = (ds['dMLTdt'] - MLG_frc).rename('MLG_nonfrc')

dTdz_b_over_MLD = (ds['dTdz_b'] / ds['MLD']).rename("dTdz_b_over_MLD")

MLG_res2 = (ds['dMLTdt'] - (
      ds['MLG_frc_sw']
    + ds['MLG_frc_lw']
    + ds['MLG_frc_sh']
    + ds['MLG_frc_lh']
    + ds['MLG_frc_fwf']
    + ds['MLG_rescale']
    + ds["MLG_vdiff"]
    + ds["MLG_hdiff"]
    + ds["MLG_vadv"]
    + ds["MLG_hadv"]
    + ds["MLG_ent"]
)).rename('MLG_res2')

print(MLG_res2)

print("RESIDUE: ", np.amax(np.abs(ds['MLG_residue'])))

ds = xr.merge(
    [
        ds,
        MLG_frc,
        MLG_nonfrc,
        MLG_adv,
        MLG_diff,
        MLG_res2,
        dTdz_b_over_MLD,
    ]
)

ds = xr.merge([ds, MLG_frc])

m2wm = np.vectorize( lambda m: ((m - 10) % 12 + 1) )
wm2m = np.vectorize( lambda wm:  ((wm + 9) - 1) % 12 + 1)

ds = ds.where(ds.time.dt.month.isin(wm2m(args.watermonths)))
    
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

    'MLD' : {
        'var'  : "h",
        'unit' : "$ \\mathrm{m} $",
    },

    'MXLDEPTH' : {
        'var'  : "MXLDEPTH",
        'unit' : "$ \\mathrm{m} $",
    },


    'dTdz_b' : {
        'var'  : "$\\partial T_\\mathrm{b} / \\partial z$",
        'unit' : "$ \\mathrm{K} \\,/\\, \\mathrm{m} $",
    },

    'dTdz_b_over_MLD' : {
        'var'  : "$\\partial T_\\mathrm{b} / \\partial z \\, / \\, h$",
        'unit' : "$ \\mathrm{K} \\,/\\, \\mathrm{m} $",
    },




    'MLG_MLB' : {
        'var'  : "$\\dot{T}_\\mathrm{ent} + \\dot{T}_\\mathrm{vdiff}$",
        'unit' : "$ \\mathrm{K} / \\mathrm{s} $",
    },


    'MLG_ent' : {
        'var'  : "$\\dot{T}_\\mathrm{ent}$",
        'unit' : "$ \\mathrm{K} / \\mathrm{s} $",
    },


    'MLG_nondiff' : {
        'var'  : "$\\dot{T}_\\mathrm{nondiff}$",
        'unit' : "$ \\mathrm{K} / \\mathrm{s} $",
    },

    'ttl_frc' : {
        'var'  : "$F_\\mathrm{sfc}$",
        'unit' : "$ \\mathrm{W} / \\mathrm{m}^2 $",
    },


    'MLG_vdiff' : {
        'var'  : "$\\dot{T}_\\mathrm{vdiff}$",
        'unit' : "$ \\mathrm{K} / \\mathrm{s} $",
    },


    'MLG_diff' : {
        'var'  : "$\\dot{T}_\\mathrm{diff}$",
        'unit' : "$ \\mathrm{K} / \\mathrm{s} $",
    },

    'MLG_hdiff' : {
        'var'  : "$\\dot{T}_\\mathrm{hdiff}$",
        'unit' : "$ \\mathrm{K} / \\mathrm{s} $",
    },


    'MLG_frc' : {
        'var'  : "$\\dot{T}_\\mathrm{frc}$",
        'unit' : "$ \\mathrm{K} / \\mathrm{s} $",
    },

    'MLG_res' : {
        'var'  : "$\\dot{T}_\\mathrm{res}$",
        'unit' : "$ \\mathrm{K} / \\mathrm{s} $",
    },

    'MLG_adv' : {
        'var'  : "$\\dot{T}_\\mathrm{adv}$",
        'unit' : "$ \\mathrm{K} / \\mathrm{s} $",
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
    ('MLG_vdiff', 'dTdz_b'),            ('MLG_frc', 'dTdz_b'),          ('MLG_frc', 'MLG_adv'),
    ('MLG_vdiff', 'MLD'),               ('MLG_frc', 'MLD'),             None,
    ('MLG_vdiff', 'dTdz_b_over_MLD'),   ('MLG_frc', 'dTdz_b_over_MLD'), None,
]

rows = 3


if len(plot_data) % rows != 0:
    cols = len(plot_data) // rows + 1
else:    
    cols = len(plot_data) // rows



color_info = {

    'watertime': {
        'cmap' : 'rainbow',
        'varname'  : 'watertime', 
        'factor'   : 12.0,  # 6 months
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

def m2str(m):

    return {
        '10' : 'Oct',
        '11' : 'Nov',
        '12' : 'Dec',
        '1'  : 'Jan',
        '2'  : 'Feb',
        '3'  : 'Mar',
    }["%d" % m]


def wm2str(wm):
    return m2str(wm2m(wm))


aspect_ratio = 0.8

fig, ax = plt.subplots(rows, cols, figsize=(4*cols, 6*aspect_ratio*rows), squeeze=False, constrained_layout=True)#, subplot_kw={'aspect': aspect_ratio})#, constrained_layout=True)

ax_flat = ax.flatten()
    
#fig.suptitle("AR duration time range to do linear regression: %d - %d days\n# of cases = %d" % (args.AR_dt_rng[0], args.AR_dt_rng[1], np.sum(test_data['picked']==True)))

#fig.suptitle("IVT range: [%d, %d]" % (args.IVT_rng[0], args.IVT_rng[1]))

#ax[0, 0].set_title("%s-%s" % (wm2str(np.amin(args.watermonths[0])), wm2str(np.amax(args.watermonths[-1]))))



for i, _plot_data in enumerate(plot_data):

    if _plot_data is None:
        continue

    print("Plotting data: ", _plot_data)

    var_X = _plot_data[0]
    var_Y = _plot_data[1]

    var_info_x = var_infos[var_X]
    var_info_y = var_infos[var_Y]

    _ax = ax_flat[i]
    

    X = ds[var_X]
    Y = ds[var_Y]
    Z = np.mod(watertime, 1)
    
    mappable =  _ax.scatter(X, Y, c = Z * color_info['factor'], s=1, cmap=color_info['cmap'], vmin=color_info['bnd'][0], vmax=color_info['bnd'][1], alpha=0.5)
    #_ax.scatter(X, Y, c = "black", s=4)

    #plot_linregress(_ax, X, Y)

    _ax.set_xlabel("%s [ %s ]" % (var_info_x['var'], var_info_x['unit']))
    _ax.set_ylabel("%s [ %s ]" % (var_info_y['var'], var_info_y['unit']))
    
    cbar = plt.colorbar(mappable, ax=_ax, orientation='vertical')
    cbar.ax.set_ylabel(color_info['label'])

    #if i==0:
    #    shared_limx = np.array([-2, 0.5])
    #else:
    #    shared_limx = np.array([-3.2, 0.5])
        
    #shared_limy = np.array([-1, 0]) * (shared_limx[1] - shared_limx[0]) * aspect_ratio + 0.5
    #shared_ticks = np.arange(-3, 1)
    #_ax.set_xlim(shared_limx)
    #_ax.set_ylim(shared_limy)
    #_ax.set_xticks(shared_ticks)
    #_ax.set_yticks(shared_ticks)

    if _plot_data[1] in ["MXLDEPTH", "MLD"]:
        _ax.invert_yaxis()

for i, _ax in enumerate(ax_flat):
    
    if i > len(plot_data)-1 or plot_data[i] is None: 
        fig.delaxes(_ax)

if not args.no_display:
    print("Show figure")
    plt.show()

if args.output != "":
    print("Output figure: %s" % args.output)
    fig.savefig(args.output, dpi=200)

