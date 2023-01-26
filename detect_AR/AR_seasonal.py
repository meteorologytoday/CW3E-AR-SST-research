import numpy as np
import netCDF4
import AR_tools, NK_tools, fmon_tools, watertime_tools
import anomalies
import pandas as pd

import traceback
from pathlib import Path
import argparse
from datetime import (timedelta, datetime, timezone)

def findfirst(a):
    return np.where(a)[0][0]

def findlast(a):
    return np.where(a)[0][-1]

def within(a, m, M):

    return m <= a and a < M

parser = argparse.ArgumentParser(
                    prog = 'plot_skill',
                    description = 'Plot prediction skill of GFS on AR.',
)

parser.add_argument('--input', type=str, help='Input file', required=True)
parser.add_argument('--output', type=str, help='Output file', default="")
parser.add_argument('--AR-dt-rng', type=float, nargs=2, help='Days of AR duration to do linear regression.', default=[0.0, 30.0])
parser.add_argument('--IVT-threshold', type=float, help='Threshold of IVT to determin AR condition.', default=250.0)
parser.add_argument('--output-database', type=str, help='CSV file.', default="")
parser.add_argument('--no-display', action="store_true")
args = parser.parse_args()
print(args)

data = {
    "ttl" : {},
    "clim" : {},
    "anom" : {},
}

data_dim = {}

AR_varnames = ["IWV", "IVT", "IWVKE", "sst", "mslhf", "msshf", "msnlwrf", "msnswrf", "mtpr", "mer", "mvimd", "t2m", "u10", "v10", "U", "MLD", "dT", "db", "dTdt", "dTdt_sfchf", "dTdt_no_sfchf", "w_deepen", "dTdt_deepen", "net_sfc_hf", "net_conv_wv", "vort10", "curltau", "dTdt_Ekman"]



#AR_varnames = ["IVT", "sst", "mslhf", "msshf"]

with netCDF4.Dataset(args.input, "r") as ds:
 
    for varname in ["time", "time_clim",]:
        t = ds.variables[varname][:]
        data_dim[varname] = [ datetime(1970, 1, 1) + _t * timedelta(days=1) for _t in t]

    for k, subdata in data.items():

        for AR_varname in AR_varnames:
            subdata[AR_varname] = ds.variables["%s_%s" % (AR_varname, k)][:]



    #data['time_clim'] = data['time_clim'][:]
    #data['time'] = data['time'][:]

data['ttl']['IVT'][np.isnan(data['ttl']['IVT'])] = 0.0



# chop timeseries into year
watertimes = [watertime_tools.getWatertime(t) for t in data_dim['time']]
wateryears = [np.floor(wt) for wt in watertimes]
waterdates = [watertimes[t] - wateryears[t] for t in range(len(watertimes))]
unique_wateryears = np.unique(wateryears)

wateryear_idxes = []

for i, wateryear in enumerate(unique_wateryears):
   
    flags = wateryears == wateryear
    wateryear_beg_idx = findfirst(flags)
    wateryear_end_idx = findlast(flags)
        
    idx = slice(wateryear_beg_idx, wateryear_end_idx+1)
    
    wateryear_idxes.append(
        (wateryear, idx)
    )            


AR_evts_by_wateryear = {}

for wateryear, idx in wateryear_idxes:

    AR_evts = []
    AR_evts_by_wateryear["%04d" % wateryear] = AR_evts


    AR_t_segs, AR_t_inds = AR_tools.detectAbove(data_dim['time'][idx], data['ttl']['IVT'][idx], args.IVT_threshold, glue_threshold=timedelta(hours=24))


    for k, AR_t_seg in enumerate(AR_t_segs):

        AR_evts.append({
            'watertime_beg' : watertime_tools.getWatertime(AR_t_seg[0]),
            'watertime_end' : watertime_tools.getWatertime(AR_t_seg[1]),
        })


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
from scipy.stats import linregress

print("done")


var_infos = {

    'dTdt' : {
        'var'  : "$\\dot{T}_\\mathrm{ttl}$",
        'unit' : "$ \\mathrm{T} / \\mathrm{s} $",
    },

    'dTdt_sfchf' : {
        'var'  : "$\\dot{T}_\\mathrm{shf}$",
        'unit' : "$ \\mathrm{T} / \\mathrm{s} $",
    },

    'dTdt_no_sfchf' : {
        'var'  : "$\\dot{T}_\\mathrm{ttl} - \\dot{T}_\\mathrm{shf}$",
        'unit' : "$ \\mathrm{T} / \\mathrm{s} $",
    },

    'dTdt_deepen' : {
        'var'  : "$\\dot{T}_\\mathrm{e-deep}$",
        'unit' : "$ \\mathrm{T} / \\mathrm{s} $",
    },

    'w_deepen' : {
        'var'  : "$w_\\mathrm{deep}$",
        'unit' : "$ \\mathrm{m} / \\mathrm{s} $",
    },

    'vort10' : {
        'var'  : "$\\cdot\\zeta$",
        'unit' : "$ 1 / \\mathrm{s} $",
    },

    'curltau' : {
        'var'  : "$\\hat{k}\\cdot\\nabla \\times \\vec{\\tau}$",
        'unit' : "$ 1 / \\mathrm{s} $",
    },


    'dTdt_Ekman' : {
        'var'  : "$\\dot{T}_{\\mathrm{e-Ekman-pump}}$",
        'unit' : "$ \\mathrm{K} / \\mathrm{s} $",
    },

    'MLD' : {
        'var'  : "Mixed-layer depth",
        'unit' : "$ \\mathrm{m}$",
    },



    'U' : {
        'var'  : "$\\left|\\vec{U}_{10}\\right|$",
        'unit' : "$ \\mathrm{m} / \\mathrm{s} $",
    },

    'Delta' : {
        'var'  : "$ \\Delta b w_e$",
        'unit' : "$ \\mathrm{m}^2 / \\mathrm{s}^3 $",
    },

    'DeltaOnlyU' : {
        'var'  : "$ \\Delta b w_e$ -- U only",
        'unit' : "$ \\mathrm{m}^2 / \\mathrm{s}^3 $",
    },


    'dt' : {
        'var'  : "$\\Delta t_{\\mathrm{AR}}$",
        'unit' : "$ \\mathrm{s} $",
    },

    'mean_IVT' : {
        'var'  : "$\\mathrm{IVT}_{\\mathrm{mean}}$",
        'unit' : "$ \\mathrm{kg} \\, \\mathrm{m} / \\mathrm{s} $",
    },

    'max_IVT' : {
        'var'  : "$\\mathrm{IVT}_{\\mathrm{max}}$",
        'unit' : "$ \\mathrm{kg} \\, \\mathrm{m} / \\mathrm{s} $",
    },



}

rows = len(unique_wateryears)


color_info = {
    'AR_duration': {
        'cmap' : 'bone_r',
        'varname'  : 'dt',
        'factor'   : 1 / 86400.0,
        'label' : "AR duration [ day ]", 
        'bnd'   : [0, 10],
    },

    'watertime': {
        'cmap' : 'rainbow',
        'varname'  : 'watertime', 
        'factor'   : 6.0,  # 6 months
        'label' : "Watertime [ mon ]",
        'bnd'   : [0, 6],
    },

}['watertime']

fig, ax = plt.subplots(rows, 1, figsize=(8, 2*rows), squeeze=False, gridspec_kw = dict(hspace=0.3, wspace=0.4))

ax_flat = ax.flatten()
    

fig.suptitle("")

for i, wateryear in enumerate(unique_wateryears):

    _ax = ax_flat[i]
    wateryear_idx = wateryear_idxes[i][1]

    t = waterdates[wateryear_idx]
    sst = data['ttl']['sst'][wateryear_idx]

    _ax.plot(t, sst, color="black")
    print(sst)


    AR_evts = AR_evts_by_wateryear["%04d" % wateryear]

    trans = transforms.blended_transform_factory(_ax.transData, _ax.transAxes)

    for j, AR_evt in enumerate(AR_evts):
   
        print("Evt %d: %.3f to %.3f" % (j, AR_evt['watertime_beg'], AR_evt['watertime_end']))
 
        dwatertime = AR_evt['watertime_end'] - AR_evt['watertime_beg']
        rect = Rectangle((AR_evt['watertime_beg'] % 1, 0), dwatertime, 1, transform=trans, color='gray', alpha=0.5)
        _ax.add_patch(rect)
                



    
ax_flat[-1].set_xlabel("waterdate")
#_ax.set_title("Wateryear: %d" % (wateryear,))
    

if not args.no_display:
    plt.show()

if args.output != "":
    print("Output figure: %s" % args.output)
    fig.savefig(args.output, dpi=200)

