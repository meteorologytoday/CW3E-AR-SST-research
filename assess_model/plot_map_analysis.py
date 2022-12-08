import numpy as np
import load_data
import netCDF4
from datetime import (datetime, timedelta)
import AR_tools

from pathlib import Path
import argparse

parser = argparse.ArgumentParser(
                    prog = 'plot_skill',
                    description = 'Plot prediction skill of GFS on AR.',
)
parser.add_argument('--AR-glue-threshold', type=int, help='AR glue threshold in hours.', default=0)

parser.add_argument('--AR-npz', type=str, help='Save the timeseries into an npz file', required=True)
parser.add_argument('--ts-npz', type=str, nargs="+", help='Have to be a pair: label, filename.', required=True)
parser.add_argument('--lines', type=str, nargs="+", help='Have to be a pair: linecolor, linestyle.', default=["red", "solid"])
parser.add_argument('--output', type=str, help='Output directory.', default="")
parser.add_argument('--no-display', action="store_true")

args = parser.parse_args()

print(args)

if len(args.ts_npz) % 2 != 0:
    print("Error: --ts-npz has to be pairs")

if len(args.ts_npz) != len(args.lines):
    print("Error: --lines has to be consistent with --ts-npz")



AR_data    = np.load(args.AR_npz)


ts_data = {}

for i in range(len(args.ts_npz) // 2):
    _ts_data = np.load(args.ts_npz[i*2+1])
    
    # Reconstruct a normal dict so that I can calculate water budget
    _ts_data = dict( zip( (k for k in _ts_data.keys()) , (_ts_data[k] for k in _ts_data.keys()) ) )
    _ts_data['time'] = _ts_data['time'].tolist()

    _ts_data['WBGT'] = ( _ts_data['LHTFL'] / 2.5e6 - _ts_data['PRATE'] ) * 86400.0

    _ts_data['line'] = {
        'lc' : args.lines[i*2],
        'ls' : args.lines[i*2+1],
    }

    ts_data[args.ts_npz[i*2]] = _ts_data

# test if times are the same
t_AR = AR_data['time'].tolist()


#if len(t_AR) != len(t_sk):
#    raise Exception("Time vectors have different lengths.")




"""
for i in range(len(t_AR)):
    if t_AR[i] != t_sk[i]:
        raise Exception("Time vectors have different values.")
"""
# Plot data
print("Loading Matplotlib...")
import matplotlib
if args.no_display is False:
    matplotlib.use('TkAgg')
else:
    matplotlib.use('Agg')
    
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Rectangle
import matplotlib.transforms as transforms
from matplotlib.dates import DateFormatter

import cartopy.crs as ccrs
print("done")

#AR_vars = ["IWV", "IVT", "IWVKE", "HGT_850mb"]
AR_vars = ["IWV", "IVT", "HGT_500mb", "HGT_850mb", "PRATE", "LHTFL", "WBGT"]
#AR_vars = ["IWV", "IVT", ]

var_infos = {

    'IWV' : {
        'label' : "[ $ \\mathrm{kg} / \\mathrm{m}^2 $ ]",
    },

    'IVT' : {
        'label' : "[ $ \\mathrm{kg} \\, \\mathrm{m} / \\mathrm{s} / \\mathrm{m}^2 $ ]",
    },


    'HGT_850mb' : {
        'label' : "[ $ \\mathrm{m} $ ]",
    },
    
    'HGT_500mb' : {
        'label' : "[ $ \\mathrm{m} $ ]",
    },

    'PRATE' : {
        'label' : "[ $ \\mathrm{kg} / \\mathrm{s} / \\mathrm{m}^2 $ ]",
    },


    'LHTFL' : {
        'label' : "[ $ \\mathrm{W} / \\mathrm{m}^2 $ ]",
    },

    'SHTFL' : {
        'label' : "[ $ \\mathrm{W} / \\mathrm{m}^2 $ ]",
    },


    'WBGT' : {
        'label' : "[ $ \\mathrm{kg} / \\mathrm{day} / \\mathrm{m}^2 $ ]",
    },


}



fig, ax = plt.subplots(len(AR_vars), 1, figsize=(12, 4*len(AR_vars)), sharex=True)

t_segs = AR_tools.detectAbove(t_AR, AR_data["IVT"], 250, glue_threshold=timedelta(hours=args.AR_glue_threshold))

analysis_plotted = np.zeros((len(ts_data.keys()),))
for i, AR_var in enumerate(AR_vars):

    print("Plot %s" % (AR_var,))
    _ax = ax[i]

    trans = transforms.blended_transform_factory(_ax.transData, _ax.transAxes)
        
    for k, t_seg in enumerate(t_segs):
        _ax.add_patch( Rectangle((t_seg[0], 0), t_seg[1] - t_seg[0], 1, fc ='#cccccc', ec='none', transform=trans, zorder=1)) 

    if AR_var == "IVT":
        _ax.plot(t_AR, AR_data[AR_var], color="gray", zorder=3, linewidth=1)


    for k, data in ts_data.items():

        t = data["time"]

        if analysis_plotted[i] == 0:
            analysis_plotted[i] = 1
            _ax.plot(t, data[AR_var][0, :], "k-", linewidth=2, zorder=4, label="analysis")

        _ax.plot(t, data[AR_var][1, :], color=data['line']['lc'], linestyle=data['line']['ls'], linewidth=2, zorder=5, label=k)

  
    #_ax.set_title("Forecast error of %s" % (AR_var, ))
    _ax.set_ylabel("%s\n%s" % (AR_var, var_infos[AR_var]['label']))

    _ax.xaxis.set_major_formatter(DateFormatter('%Y-%m-%d'))

    _ax.legend()

if not args.no_display:
    plt.show()

if args.output != "":
    
    fig.savefig(args.output, dpi=200)


