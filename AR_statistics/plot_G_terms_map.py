import numpy as np
import fmon_tools, watertime_tools
import anomalies
import ARstat_tool
import xarray as xr

import traceback
from pathlib import Path
import argparse
from datetime import (timedelta, datetime, timezone)

parser = argparse.ArgumentParser(
                    prog = 'plot_skill',
                    description = 'Plot prediction skill of GFS on AR.',
)

parser.add_argument('--input-dir', type=str, help='Input file', required=True)
parser.add_argument('--beg-year', type=int, help='Input file', required=True)
parser.add_argument('--output', type=str, help='Output file', default="")
parser.add_argument('--title', type=str, help='Output title', default="")
parser.add_argument('--breakdown', type=str, help='Output title', default="atmocn", choices=["atmocn", "atm", "ocn"])
parser.add_argument('--no-display', action="store_true")
args = parser.parse_args()
print(args)

plotted_varnames = {
    "atmocn" : ["dMLTdt", "MLG_frc", "MLG_nonfrc"],
    "atm" : ["MLG_frc", "MLG_frc_sw", "MLG_frc_lw", "MLG_frc_sh", "MLG_frc_lh", "MLG_frc_fwf"],
    "ocn" : ["MLG_nonfrc", "MLG_adv", "MLG_vdiff", "MLG_ent", "MLG_hdiff", "MLG_res2"],
}[args.breakdown]

print(plotted_varnames)


ds = xr.open_dataset(args.input).astype(np.float64)

MLG_frc = (ds['MLG_frc_sw'] + ds['MLG_frc_lw'] + ds['MLG_frc_sh']  + ds['MLG_frc_lh'] + ds['MLG_frc_fwf'] + ds['MLG_rescale']).rename('MLG_frc')

MLG_adv = (ds['MLG_hadv'] + ds['MLG_vadv']).rename('MLG_adv')
MLG_diff = (ds['MLG_vdiff'] + ds['MLG_hdiff']).rename('MLG_diff')
MLG_nonfrc = (MLG_adv + MLG_diff + ds['MLG_ent']).rename('MLG_nonfrc')

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
    ]
)

# Compute mean and anomaly
def convertdt64todt(dt64):
    return datetime.fromtimestamp( (dt64 - np.datetime64("1970-01-01")) / np.timedelta64(1, 's' ))
 
ts_np64 = ds.time.to_numpy().astype("datetime64[s]")
ts_dt   = ts_np64.astype(object)

ds_mean = []
ds_anom = []
for _, varname in enumerate(ds.keys()):
    
    _var = ds[varname]
    xs = _var.to_numpy()
    
    tm, xm, xa, cnt = anomalies.decomposeClimAnom(ts_dt, xs)

    ds_mean.append(
        xr.DataArray(
            name = varname,
            data = xm,
            dims = ["time"],
            coords = {
                "time" : np.array(tm).astype("datetime64[ns]"),
            }
        )
    )

    ds_anom.append(
        xr.DataArray(
            name = varname,
            data = xa,
            dims = ["time"],
            coords = {
                "time" : ts_np64,
            }
        )
    )


ds_anom = xr.merge(ds_anom)
ds_mean = xr.merge(ds_mean)

# Construct
t_months = np.arange(1, 7)

ds_stats = {}

for condition_name, (IVT_min, IVT_max) in [
    ("clim",  (0, np.inf)),
    ("ARf",   (0, 250)),
    ("AR",    (250, np.inf)),
    ("AR+ARf",    (None, None)),
]:

    _tmp = {}
    for varname, _ in ds.items():
        _tmp[varname] = (["time", "stat"], np.zeros((len(t_months), 5)) )

    ds_stat = xr.Dataset(
        _tmp,

        coords = {
            "time" : t_months,
            "stat" : ["mean", "q25", "q75", "q10", "q90"]
        }
    )

    ds_stats[condition_name] = ds_stat

    for t, m in enumerate(t_months): 
        
        #print("m = ", m, "; wm = ", watertime_tools.wm2m(m))
        
        if condition_name == "clim":
            _ds_ref = ds
            _cond = _ds_ref.time.dt.month.isin(watertime_tools.wm2m(m))
 
        elif condition_name == "AR+ARf":
            _ds_ref = ds_anom
            _cond = _ds_ref.time.dt.month.isin(watertime_tools.wm2m(m))
        
        elif condition_name in ["AR", "ARf"]:
            _ds_ref = ds_anom
            _cond = (
                  ( ds.IVT >= IVT_min ) 
                & ( ds.IVT <  IVT_max ) 
                & _ds_ref.time.dt.month.isin(watertime_tools.wm2m(m))
            )
   
        else:
            raise Exception("Unknown condition_name: ", condition_name) 

        print("Conditon %s : cnt = %d, wm = %d" % (condition_name, np.sum(_cond), m))

        _ds = _ds_ref.where(_cond)
        #_ds = _ds_ref
            
        for varname, _ in ds_stat.items():

            _data = _ds[varname] * 1e6
            ds_stat[varname][t, 0] = np.nanmean(_data)
            ds_stat[varname][t, 1] = np.nanquantile(_data, .25)
            ds_stat[varname][t, 2] = np.nanquantile(_data, .75)
            ds_stat[varname][t, 3] = np.nanquantile(_data, .10)
            ds_stat[varname][t, 4] = np.nanquantile(_data, .90)
            

    print(IVT_min, " - ", IVT_max)

ds_stats["AR-ARf"] = ds_stats["AR"] - ds_stats["ARf"]



plot_infos = {

    "dMLTdt" : {
        "label" : "$\\dot{T}_\\mathrm{ttl}$",
    },


    "MLG_ttl" : {
        "label" : "$\\dot{T}_\\mathrm{ttl}$",
    },

    "MLG_frc" : {
        "label" : "$\\dot{T}_\\mathrm{frc}$",
    },

    "MLG_nonfrc" : {
        "label" : "$\\dot{T}_\\mathrm{nfrc}$",
    },


    "MLG_frc_lw" : {
        "label" : "$\\dot{T}_\\mathrm{lw}$",
    },

    "MLG_frc_sw" : {
        "label" : "$\\dot{T}_\\mathrm{sw}$",
    },

    "MLG_frc_lh" : {
        "label" : "$\\dot{T}_\\mathrm{lh}$",
    },

    "MLG_frc_sh" : {
        "label" : "$\\dot{T}_\\mathrm{sh}$",
    },

    "MLG_frc_fwf" : {
        "label" : "$\\dot{T}_\\mathrm{fwf}$",
    },


    "MLG_adv" : {
        "label" : "$\\dot{T}_\\mathrm{adv}$",
    },

    "MLG_vdiff" : {
        "label" : "$\\dot{T}_\\mathrm{vdiff}$",
    },

    "MLG_ent" : {
        "label" : "$\\dot{T}_\\mathrm{ent}$",
    },

    "MLG_hdiff" : {
        "label" : "$\\dot{T}_\\mathrm{hdiff}$",
    },


    "MLG_res2" : {
        "label" : "$\\dot{T}_\\mathrm{res}$",
    },



}

plot_infos_scnario = {

    "clim" : {
        "title" : "All",
    },

    "AR" : {
        "title" : "AR",
    },

    "ARf" : {
        "title" : "AR free",
    },

    "AR-ARf" : {
        "title" : "AR minus AR free",
    }

}

plot_ylim = {

    "atmocn" : {
        "mean" : [-1.5, 1.5],
        "anom" : [-0.7, 0.7],
    },

    "atm" : {
        "mean" : [-1.5, 1.5],
        "anom" : [-0.3, 0.7],
    },

    "ocn" : {
        "mean" : [-0.5, 0.5],
        #"anom" : [-0.06, 0.01],
        "anom" : [-0.5, 0.5],
    },

}

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

fig, ax = plt.subplots(2, 1, figsize=(6, 8), squeeze=False, constrained_layout=True)# gridspec_kw = dict(hspace=0.3, wspace=0.4))

if args.title == "":
    fig.suptitle(args.input)
else:
    fig.suptitle(args.title)

bar_width = 0.8 / len(plotted_varnames) #0.15


for s, sname in enumerate(["clim", "AR-ARf"]):

    ds_stat = ds_stats[sname]
    
    _ax = ax[s, 0]

    for i, varname in enumerate(plotted_varnames):

        _data = ds_stat[varname].to_numpy()
        _ax.bar(t_months + i*bar_width, _data[:, 0],   bar_width, label=plot_infos[varname]['label'])


        _anom_ARpARf_data = ds_stats["AR+ARf"][varname].to_numpy()
        _anom_AR_data = ds_stats["AR"][varname].to_numpy()
        _anom_ARf_data = ds_stats["ARf"][varname].to_numpy()
        _mean_data = ds_stats["clim"][varname].to_numpy()

        # error bar
        _error_bar_lower = _anom_ARpARf_data[:, 1]
        _error_bar_upper = _anom_ARpARf_data[:, 2]
        _offset = _data[:, 0]

        for m, t_month in enumerate(t_months): 
            _ax.plot(
                [t_month + i*bar_width]*2,
                np.array([_error_bar_lower[m], _error_bar_upper]) + _offset[m],
                "k-",
            )

    _ax.set_xticks(t_months)
    _ax.set_xticklabels(np.vectorize(watertime_tools.wm2m)(t_months))

    _ax.set_xlabel("Month")
    _ax.set_ylabel("[ $ 1 \\times 10^{-6} \\mathrm{K} \\, / \\, \\mathrm{s} $ ]")

    _ax.set_title("(%s) %s " % (
        "abcdefghijklmn"[s],
        plot_infos_scnario[sname]['title'],
    ))

    
    _ax.set_xlim([0.5, 8.5])

    _ax.legend(loc="center right", borderpad=0.1, labelspacing=0.1)
    
    if s in [0, ]:
        _ax.set_ylim(plot_ylim[args.breakdown]['mean'])
    else:
        _ax.set_ylim(plot_ylim[args.breakdown]['anom'])

    _ax.grid(True)

if args.output != "":
   
    print("Output filename: %s" % (args.output,))
    fig.savefig(args.output, dpi=200)


if not args.no_display:
    print("Show figure")
    plt.show()


