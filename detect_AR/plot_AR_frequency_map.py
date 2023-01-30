import argparse
import numpy as np
import netCDF4

parser = argparse.ArgumentParser(
                    prog = 'plot_rectangular',
                    description = 'Plot a box on world map',
)

parser.add_argument('--no-display', action="store_true")
parser.add_argument('--lat-rng', type=float, nargs=2, help='Latitude  range', required=True)
parser.add_argument('--lon-rng', type=float, nargs=2, help='Longitude range. 0-360', required=True)
parser.add_argument('--input',  nargs=2,    type=str, help='Input file.', required=True)
parser.add_argument('--output-dir', type=str, help='Output direrctory.', default="")

args = parser.parse_args()

print(args)

data = []

varnames  = [ "IVT", "IWV", "net_sfc_hf", "mslhf", "msshf", "msnswrf", "msnlwrf", "dTdt", "U"]
subgroups = [ "avg", "std", "cnt", "cnt_ttl"]

for i, filename in enumerate(args.input):

    print("Loading file : %s" % (filename,))

    _data = { subgroup : {} for subgroup in subgroups }
    _data["frq"] = {}

    with netCDF4.Dataset(filename, "r") as f:

        nbox = f.dimensions["box"].size
        nlat = f.dimensions["lat"].size
        nlon = f.dimensions["lon"].size
        
        lat = f.variables["lat"][:]
        lon = f.variables["lon"][:]

         
        for varname in varnames:

            for subgroup in subgroups:
                _data[subgroup][varname] = f.groups[subgroup].variables[varname][:]
  
            if "frq" not in _data["avg"]:
                _data["avg"]["frq"] = _data["cnt"][varname] / _data["cnt_ttl"][varname]
   
        boxnames = netCDF4.chartostring(f.variables["boxname"][:])


    print(boxnames) 
    data.append(_data)

    

import matplotlib
if args.no_display is False:
    matplotlib.use('TkAgg')
else:
    matplotlib.use('Agg')
    
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.patches as mpatches

import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


print("done")

plot_info = {

    "IVT" : {
        "name" : "IVT",
        "unit" : "[ $ \\mathrm{kg} \\, \\mathrm{m} \\, / \\, \\mathrm{s} \\, / \\, \\mathrm{m}^2$ ]",
        "diff" : {
            "levs" : np.linspace(-500, 500, 21),
            "cmap" : cm.get_cmap("bwr"),
        },
    },

    "IWV" : {
        "name" : "IWV",
        "unit" : "[ $ \\mathrm{kg} \\, / \\, \\mathrm{m}^2$ ]",
        "diff" : {
            "levs" : np.linspace(-20, 20, 11),
            "cmap" : cm.get_cmap("bwr"),
        },
    },



    "msnswrf" : {
        "name" : "SW Rad HF",
        "unit" : "[ $ \\mathrm{W} \\, / \\, \\mathrm{s} $ ]",
        "diff" : {
            "levs" : np.linspace(-150, 150, 31),
            "cmap" : cm.get_cmap("bwr"),
        },
    },



    "msnlwrf" : {
        "name" : "LW Rad HF",
        "unit" : "[ $ \\mathrm{W} \\, / \\, \\mathrm{s} $ ]",
        "diff" : {
            "levs" : np.linspace(-150, 150, 31),
            "cmap" : cm.get_cmap("bwr"),
        },
    },


    "msshf" : {
        "name" : "Sensible HF",
        "unit" : "[ $ \\mathrm{W} \\, / \\, \\mathrm{s} $ ]",
        "diff" : {
            "levs" : np.linspace(-150, 150, 31),
            "cmap" : cm.get_cmap("bwr"),
        },
    },


    "mslhf" : {
        "name" : "Latent HF",
        "unit" : "[ $ \\mathrm{W} \\, / \\, \\mathrm{s} $ ]",
        "diff" : {
            "levs" : np.linspace(-150, 150, 31),
            "cmap" : cm.get_cmap("bwr"),
        },
    },


    "net_sfc_hf" : {
        "name" : "Net SHF",
        "unit" : "[ $ \\mathrm{W} \\, / \\, \\mathrm{s} $ ]",
        "diff" : {
            "levs" : np.linspace(-150, 150, 31),
            "cmap" : cm.get_cmap("bwr"),
        },
    },

    "U" : {
        "name" : "Wind speed",
        "unit" : "[ $ \\mathrm{m} \\, / \\, \\mathrm{s} $ ]",
        "diff" : {
            "levs" : np.linspace(-10, 10, 11),
            "cmap" : cm.get_cmap("bwr"),
        },
    },

    "frq" : {
        "name" : "Frequency",
        "unit" : "",
        "avg" : {
            "levs" : np.linspace(0, 1, 11),
            "cmap" : cm.get_cmap("hot_r"),
        },
    },
}

def plot_diff(ax, box, varname):

    pinfo = plot_info[varname]["diff"]
    factor = plot_info[varname]["factor"] if "factor" in plot_info[varname] else 1.0

    diff = (data[0]["avg"][varname][box, :, :] - data[1]["avg"][varname][box, :, :]) * factor
    std  = (data[0]["std"][varname][box, :, :] + data[1]["std"][varname][box, :, :]) / 2 * factor
    sigf = (np.abs(diff) >= std).astype(np.float32) + 0.25  # 0.25 = false, 0.75 = true
    
    mappable = ax.contourf(lon, lat, diff, pinfo["levs"], cmap=pinfo["cmap"], transform=proj_norm, extend="both")
    cb = plt.colorbar(mappable, ax=ax, orientation="vertical")

    cs = ax.contourf(lon, lat, sigf, [0, 0.5, np.inf], colors="none", hatches=[None, "..", ], transform=proj_norm)
    plt.setp(cs.collections , linewidth=.5, edgecolor="black")

    ax.set_title(plot_info[varname]["name"])
    cb.ax.set_ylabel(plot_info[varname]["unit"])



def plot_avg(ax, data_i, box, varname):

    _data = data[data_i]["avg"]

    pinfo = plot_info[varname]["avg"]
    factor = plot_info[varname]["factor"] if "factor" in plot_info[varname] else 1.0

    mappable = ax.contourf(lon, lat, _data[varname][box, :, :] * factor, pinfo["levs"], cmap=pinfo["cmap"], transform=proj_norm)
    cb = plt.colorbar(mappable, ax=ax, orientation="vertical")


    ax.set_title(plot_info[varname]["name"])
    cb.ax.set_ylabel(plot_info[varname]["unit"])


lat_t = args.lat_rng[0]
lat_b = args.lat_rng[1]

lon_l = args.lon_rng[0]
lon_r = args.lon_rng[1]

dlon = lon_r - lon_l
dlat = lat_t - lat_b

plot_lat_b = args.lat_rng[0]
plot_lat_t = args.lat_rng[1]

plot_lon_l = args.lon_rng[0]
plot_lon_r = args.lon_rng[1]

cent_lon = 180.0

proj = ccrs.PlateCarree(central_longitude=cent_lon)
proj_norm = ccrs.PlateCarree()

for box in range(0, nbox, 1):

    boxname = boxnames[box]
    print("Box: %d => %s" % (box, boxname))
    
    fig, ax = plt.subplots(3, 3, figsize=(15, 10), subplot_kw=dict(projection=proj))

    fig.suptitle(boxname)
    

    plot_avg(ax[0, 0], 0, box, "frq")
    plot_diff(ax[0, 1], box, "net_sfc_hf")
    plot_diff(ax[1, 0], box, "msshf")
    plot_diff(ax[1, 1], box, "mslhf")
    plot_diff(ax[2, 0], box, "msnlwrf")
    plot_diff(ax[2, 1], box, "msnswrf")

    plot_diff(ax[0, 2], box, "IVT")
    plot_diff(ax[1, 2], box, "IWV")
    plot_diff(ax[2, 2], box, "U")



    for _ax in ax.flatten():
        _ax.set_global()



        _ax.gridlines()
        _ax.coastlines()

        gl = _ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=2, color='gray', alpha=0.5, linestyle='--')
        gl.xlabels_top = False
        gl.ylabels_right = False
        #gl.xlocator = mticker.FixedLocator(np.arange(-180, 181, 60))
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 10, 'color': 'black'}
        gl.ylabel_style = {'size': 10, 'color': 'black'}


        _ax.set_extent([plot_lon_l, plot_lon_r, plot_lat_b, plot_lat_t], crs=proj_norm)


    
    if args.output_dir != "":
       
        output_filename = "%s/AR_freq_map_%s.png" % (args.output_dir, boxname,)
        print("Output filename: %s" % (output_filename,))
        fig.savefig(output_filename, dpi=200)


    


if not args.no_display:
    plt.show()


