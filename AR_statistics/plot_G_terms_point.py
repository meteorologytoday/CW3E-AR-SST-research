import numpy as np
#import fmon_tools, watertime_tools
import anomalies
import ARstat_tool
import xarray as xr
import watertime_tools
import traceback
from pathlib import Path
import argparse
from datetime import (timedelta, datetime, timezone)
from scipy.stats import ttest_ind_from_stats

parser = argparse.ArgumentParser(
                    prog = 'plot_skill',
                    description = 'Plot prediction skill of GFS on AR.',
)

parser.add_argument('--input-dir', type=str, help='Input file', required=True)
parser.add_argument('--lat', type=float, help='Latitude', required=True)
parser.add_argument('--lon', type=float, help='Longitude', required=True)
parser.add_argument('--output', type=str, help='Output file', default="")
parser.add_argument('--title', type=str, help='Output title', default="")
parser.add_argument('--title-style', type=str, help='Output title', default="folder", choices=["folder", "latlon"])
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


t_months = np.array([1, 2, 3, 4, 5, 6])

ds_stat = {}
for k in ["clim", "AR", "ARf", "AR-ARf", "AR+ARf"]:
    ds_stat[k] = xr.open_dataset("%s/stat_%s.nc" % (args.input_dir, k))
   
    ds = ds_stat[k] 
    
    MLG_res2 = (ds['dMLTdt'] - (
          ds['MLG_frc']
        + ds['MLG_nonfrc']
        + ds['MLG_rescale']
    )).rename('MLG_res2')

    print("RESIDUE: ", np.amax(np.abs(ds['MLG_residue'])))

    ds = xr.merge(
        [
            ds,
            MLG_res2,
        ]
    )

    ds_stat[k] = ds

    ds = None 

print(ds_stat["AR"].coords["lat"])
print(ds_stat["AR"].coords["lon"])
    
print("Selecting data at (lat, lon) = (%.2f, %.2f)" % (args.lat, args.lon))

for k, _ in ds_stat.items():
    ds_stat[k] = ds_stat[k].sel(lat=args.lat, lon=args.lon, method="nearest")


final_lat = ds_stat["AR"].coords["lat"].data
final_lon = ds_stat["AR"].coords["lon"].data

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


shared_levels = np.linspace(-1, 1, 11) * 0.5
plot_infos = {
    "dMLTdt" : {
        "levels": shared_levels,
        "label" : "$ G_{\mathrm{ttl}} $",
        "color" : "gray",
        "hatch" : '///',
    },

    "MLG_frc" : {
        "levels": shared_levels,
        "label" : "$ G_{\mathrm{frc}} $",
        "color" : "orangered",
        "hatch" : '///',
    }, 

    "MLG_nonfrc" : {
        "levels": shared_levels,
        "label" : "$ G_{\mathrm{nfrc}} $",
        "color" : "dodgerblue",
        "hatch" : '///',
    }, 

    "MLG_adv" : {
        "levels": shared_levels,
        "label" : "$ G_{\mathrm{adv}} $",
    }, 

    "MLG_vdiff" : {
        "levels": shared_levels,
        "label" : "$ G_{\mathrm{vdiff}} $",
    }, 

    "MLG_ent" : {
        "levels": shared_levels,
        "label" : "$ G_{\mathrm{ent}} $",
    }, 

    "MLG_hdiff" : {
        "levels": shared_levels,
        "label" : "$ G_{\mathrm{hdiff}} $",
    }, 

    "MLG_res2" : {
        "levels": shared_levels,
        "label" : "$ G_{\mathrm{res}} $",
    }, 

    "MLG_frc_sw" : {
        "levels": shared_levels,
        "label" : "$ G_{\mathrm{sw}} $",
    }, 

    "MLG_frc_lw" : {
        "levels": shared_levels,
        "label" : "$ G_{\mathrm{lw}} $",
    }, 

    "MLG_frc_lh" : {
        "levels": shared_levels,
        "label" : "$ G_{\mathrm{lh}} $",
    }, 

    "MLG_frc_sh" : {
        "levels": shared_levels,
        "label" : "$ G_{\mathrm{sh}} $",
    }, 

    "MLG_frc_fwf" : {
        "levels": shared_levels,
        "label" : "$ G_{\mathrm{fwf}} $",
    }, 


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


def pretty_latlon(lat, lon):

    lon %= 360
    
    if lon > 180:
        lon = 360 - lon
        lon_EW = "W"
    else:
        lon_EW = "E"

    if lat > 0:
        lat_NS = "N"
    elif lat == 0:
        lat_NS = "E"
    else:
        lat_NS = "S"

    lat = abs(lat)

    if lat % 1 == 0 and lon % 1 == 0:
        return "%d%s" % (lat, lat_NS), "%d%s" % (lon, lon_EW)

    else:
        return "%.2f%s" % (lat, lat_NS), "%.2f%s" % (lon, lon_EW)



if args.title == "":

    if args.title_style == "folder":
        fig.suptitle(args.input_dir)
    elif args.title_style == "latlon":
        title = ("%s %s : " % pretty_latlon(final_lat, final_lon)) + args.breakdown
        fig.suptitle(title)


else:
    fig.suptitle(args.title)

bar_width = 0.8 / len(plotted_varnames) #0.15


#for s, sname in enumerate(["clim", "AR-ARf"]):
for s, sname in enumerate(["clim", "AR"]):

    ds = ds_stat[sname]
   
    _ax = ax[s, 0]

    for i, varname in enumerate(plotted_varnames):
        _data = ds[varname].to_numpy() * 1e6

        #print("data shape: ", _data.shape)

        kwargs = {}
        if 'color' in plot_infos[varname]:
            kwargs['color'] = plot_infos[varname]['color']

        if 'hatch' in plot_infos[varname]:
            kwargs['hatch'] = plot_infos[varname]['hatch']


        if sname == "AR":
            _data -= ds_stat["clim"][varname].to_numpy() * 1e6

        _ax.bar(t_months + i*bar_width, _data[:, 0],   bar_width, label=plot_infos[varname]['label'], **kwargs)

        #print(ds)
        _anom_ARpARf_data = ds_stat["AR+ARf"][varname].to_numpy() * 1e6
        _anom_AR_data = ds_stat["AR"][varname].to_numpy() * 1e6
        _anom_ARf_data = ds_stat["ARf"][varname].to_numpy() * 1e6
        _mean_data = ds_stat["clim"][varname].to_numpy() * 1e6

        # error bar
        _error_bar_lower = - _anom_AR_data[:, 1]
        _error_bar_upper =   _anom_AR_data[:, 1]
        _offset = _data[:, 0]


        for m, t_month in enumerate(t_months): 
            _ax.plot(
                [t_month + i*bar_width] * 2,
                np.array([_error_bar_lower[m], _error_bar_upper[m]]) + _offset[m],
                "k-",
                zorder=99,
            )

            result = ttest_ind_from_stats(
                _anom_AR_data[m, 0], _anom_AR_data[m, 1], _anom_AR_data[m, 3],
                _anom_ARf_data[m, 0], _anom_ARf_data[m, 1], _anom_ARf_data[m, 3],
                equal_var=False,
                alternative='two-sided',
            )

            print("[m=%d] Result of T-test: " % (m,), result)


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

    """    
    if s in [0, ]:
        _ax.set_ylim(plot_ylim[args.breakdown]['mean'])
    else:
        _ax.set_ylim(plot_ylim[args.breakdown]['anom'])
    """

    _ax.grid(True)
    _ax.set_axisbelow(True)

        
ax[0, 0].set_ylim([-1.2, 1.2])
ax[1, 0].set_ylim([-.8, .8])


if args.output != "":
   
    print("Output filename: %s" % (args.output,))
    fig.savefig(args.output, dpi=200)


if not args.no_display:
    print("Show figure")
    plt.show()

