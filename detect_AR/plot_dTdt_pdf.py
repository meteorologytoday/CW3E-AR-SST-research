import numpy as np
import netCDF4
import AR_analysis_tools, NK_tools, fmon_tools, watertime_tools
import anomalies
import pandas as pd

import traceback
from pathlib import Path
import argparse
from datetime import (timedelta, datetime, timezone)

parser = argparse.ArgumentParser(
                    prog = 'plot_skill',
                    description = 'Plot prediction skill of GFS on AR.',
)

parser.add_argument('--input', type=str, help='Input file', required=True)
parser.add_argument('--output-dir', type=str, help='Output file', default="")
parser.add_argument('--IVT-threshold', type=float, help='Threshold of IVT to determin AR condition.', default=250.0)
parser.add_argument('--no-display', action="store_true")
args = parser.parse_args()
print(args)

data = {
    "ttl" : {},
    "clim" : {},
    "anom" : {},
}

data_dim = {}

AR_varnames = ["IWV", "IVT", "IWVKE", "sst", "mslhf", "msshf", "msnlwrf", "msnswrf", "mtpr", "mer", "mvimd", "t2m", "u10", "v10", "U", "MLD", "dT", "db", "dTdt", "dTdt_sfchf", "dTdt_no_sfchf", "w_deepen", "dTdt_deepen", "net_sfc_hf", "net_conv_wv", "vort10", "curltau", "ent_Ekman", "EkmanAdv", "EkmanAdv_adj", "ent_ADV", "ttl_ADV", "geo_ADV", "ageo_ADV"]


AR_evts = AR_analysis_tools.constructCases(args.input, args.IVT_threshold, method="isolated")

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


    'ent_Ekman' : {
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


    'U_mean' : {
        'var'  : "$\\left|\\vec{U}_{10}\\right|_{\\mathrm{mean}}$",
        'unit' : "$ \\mathrm{m} / \\mathrm{s} $",
    },


    'u10' : {
        'var'  : "$u_{10}$",
        'unit' : "$ \\mathrm{m} / \\mathrm{s} $",
    },



    'v10' : {
        'var'  : "$v_{10}$",
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


mag = 10e-6
bin_edges = np.linspace(-1, 1, 501) * mag
months = [[10, 11], [12, 1], [2, 3], [10, 11, 12, 1, 2, 3]]
titles = ["Oct-Nov", "Dec-Jan", "Feb-Mar", "Oct-Mar"]


rows = 4
cols = len(months)



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

    'waterdate': {
        'cmap' : 'rainbow',
        'varname'  : 'waterdate', 
        'factor'   : 12.0,  # 6 months
        'label' : "Watertime [ mon ]",
        'bnd'   : [0, 6],
    },


}['AR_duration']#waterdate']

fig, ax = plt.subplots(rows, cols, figsize=(6*cols, 5*rows), squeeze=False, gridspec_kw = dict(hspace=0.3, wspace=0.4))


fig.suptitle("PDF of T tendency")

ax_flat = ax.flatten()
    

for i, m in enumerate(months):

    if m is None:
        continue

    _ax = ax[:, i]
    
    _data = AR_evts.loc[AR_evts["month"].isin(m)] 

    for varname in ["dTdt", "dTdt_sfchf", "dTdt_no_sfchf"]:
    #for varname in ["dTdt_deepen"]:

        if np.any(np.abs(_data[varname]) > mag):
            raise Exception("Variable %s exceeds the magnitude" % varname) 

        _ax[0].hist(_data[varname], bin_edges, alpha=0.5, label=varname, histtype='stepfilled', density=True)

    for varname in ["dTdt_no_sfchf", "ttl_ADV", "geo_ADV", "ageo_ADV"]:
        if np.any(np.abs(_data[varname]) > mag):
            raise Exception("Variable %s exceeds the magnitude" % varname) 

        _ax[1].hist(_data[varname], bin_edges, alpha=0.5, label=varname, histtype='stepfilled', density=True)

    for varname in ["dTdt_no_sfchf", "ent_Ekman", "ent_ADV", "ent_slope",]:
        if np.any(np.abs(_data[varname]) > mag):
            raise Exception("Variable %s exceeds the magnitude" % varname) 

        _ax[2].hist(_data[varname], bin_edges, alpha=0.5, label=varname, histtype='stepfilled', density=True)


    for varname in ["dTdt_no_sfchf", "OCN_nonADV",]:

        if np.any(np.abs(_data[varname]) > mag):
            raise Exception("Variable %s exceeds the magnitude" % varname) 

        _ax[3].hist(_data[varname], bin_edges, alpha=0.5, label=varname, histtype='stepfilled', density=True)


    _ax[0].set_title(titles[i])
    #_ax.set_xlabel("%s [ %s ]" % (var_info_x['var'], var_info_x['unit']))
    #_ax.set_ylabel("%s [ %s ]" % (var_info_y['var'], var_info_y['unit']))


    for __ax in _ax.flatten():
        __ax.legend(fontsize=12)
        __ax.set_xlim(np.array([-1, 1]) * 2e-6)

if args.output_dir != "":
   
    output_filename = "%s/AR_dTdt_pdf.png" % (args.output_dir,)
    print("Output filename: %s" % (output_filename,))
    fig.savefig(output_filename, dpi=200)


if not args.no_display:
    print("Show figure")
    plt.show()


