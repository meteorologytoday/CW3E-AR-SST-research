import numpy as np
import load_data
import netCDF4
from datetime import datetime

from pathlib import Path
import argparse

parser = argparse.ArgumentParser(
                    prog = 'plot_SST',
                    description = 'Plot SST comparison',
)

parser.add_argument('--date', type=str, help='Date string: yyyy-mm-dd', required=True)
parser.add_argument('--output-dir', type=str, help='Output directory.', default="")
parser.add_argument('--no-display', action="store_true")
parser.add_argument('--rows', type=int, default=1)

args = parser.parse_args()


print(args)

# Configuration
ref_product = "ERA5"
compare_products = [ "OSTIA", "CFSR", "AVHRR", "OISST",]


data_date = datetime.strptime(args.date, "%Y-%m-%d")#, args.date)

# Prep
all_products = [ref_product,] + compare_products
datas = {}

# Load data
for product in all_products:

    print("Loading data of %s" % (product,))

    info = load_data.getFileAndIndex(product, data_date, root_dir="data")
    data = {}

    offset = 0.0

    if info['unit'] == 'K':

        offset = 273.15

    with netCDF4.Dataset(info['filename'], "r") as ds:

        data['sst'] = ds.variables[info['varnames']['sst']][info['idx'], :, :] - offset
        data['lat'] = ds.variables[info['varnames']['lat']][:]
        data['lon'] = ds.variables[info['varnames']['lon']][:]


        if product == "ERA5":
            data['sst'] = np.flip( data['sst'], axis=0 )
            data['lat'] = np.flip( data['lat'], axis=0 )


        print("[%s] Rearranging data" % (product,))
        new_lon, new_arrs = load_data.rearrangeLonAndField(data['lon'], [data['sst'],], axis=1)
        data['lon'] = new_lon
        data['sst'] = new_arrs[0]

    datas[product] = data


extra_data = {}
with netCDF4.Dataset("data/ERA5/AR_processed/ERA5_AR_%s.nc" % (data_date.strftime("%Y-%m-%d"),), "r") as ds:

    extra_data['IWV'] = np.flip(ds.variables['IWV'][0, :, :], axis=0)
    extra_data['lat'] = np.flip(ds.variables['lat'][:], axis=0)
    extra_data['lon'] = ds.variables['lon'][:]



#lnd_mask_idx = np.isnan(datas[ref_product]['sst'])
#lnd_mask_idx = datas[ref_product]['sst'] == -32767
#lnd_mask_idx = np.isfinite(datas[ref_product]['sst'])

lnd_mask_idx = datas[ref_product]['sst'].mask

#print(lnd_mask_idx.size)
#print(np.sum(lnd_mask_idx))

for product in compare_products:

    print("Interpolation compare product %s onto reference product %s" % (product, ref_product,) )

    data = datas[product]
    ref_data = datas[ref_product]
    #print(ref_data['lat'])
    
    if np.any((data['lon'][1:] - data['lon'][:-1]) <= 0):
        raise Exception("Not right")

    if np.any((ref_data['lon'][1:] - ref_data['lon'][:-1]) <= 0):

        for i, lon in enumerate(ref_data['lon']):
            print("lon[%d] = %.2f" % (i, lon))

        raise Exception("Ref Lon Not right")

    if np.any((ref_data['lat'][1:] - ref_data['lat'][:-1]) <= 0):
        raise Exception("Ref Lat Not right")



    data['interp_sst'] = load_data.interpolate(data['lat'], data['lon'], data['sst'], ref_data['lat'], ref_data['lon'], mask_idx=lnd_mask_idx)
    


# Plot data
print("Loading Matplotlib...")
import matplotlib
if args.no_display is False:
    matplotlib.use('TkAgg')
else:
    matplotlib.use('Agg')
    
import matplotlib.pyplot as plt
from matplotlib import cm
import cartopy.crs as ccrs
print("done")

cent_lon = 180.0
proj = ccrs.PlateCarree(central_longitude=cent_lon, globe=None)


rows = args.rows
cols = int(np.ceil(len(all_products) / rows)) 

fig, gridded_ax = plt.subplots(rows, cols, subplot_kw = dict(projection=proj), figsize=(6*cols, 4*rows))

ax = gridded_ax.flatten()

fig.suptitle("Date: %s" % (data_date.strftime("%Y-%m-%d"),))

for i, product in enumerate(all_products):
 
    print("Plotting product: %s" % (product,) )
    
    _ax = ax[i]

    data = datas[product]
    ref_data = datas[ref_product]

    if product == ref_product:

        sst_levs  = np.linspace(-2, 35, 38) 
        sst_ticks = [-2, 0, 5, 10, 15, 20, 25, 30, 35] 
        cmap = cm.get_cmap("nipy_spectral")
        _plot_data = data['sst']
        cb_label = "SST [${}^\\circ\\mathrm{C}$]"
        title_text = product

    else:

        sst_levs = np.linspace(-0.5, 0.5, 11) 
        sst_ticks = np.arange(-0.5, 0.6, 0.1)
        cmap = cm.get_cmap("bwr")
        cmap.set_over("brown")
        cmap.set_under("navy")
        _plot_data = data['interp_sst'] - ref_data['sst']
        cb_label = "SST bias [${}^\\circ\\mathrm{C}$]"
        title_text = "%s minus %s" % (product, ref_product)
        

    mappable = _ax.contourf(ref_data['lon'], ref_data['lat'], _plot_data, levels=sst_levs,
                transform=ccrs.PlateCarree(),
                cmap=cmap,
                extend="both")

    #cs = _ax.contour(extra_data['lon'], extra_data['lat'], extra_data['IWV'], levels=[20, 40, 60, 80],
    #            transform=ccrs.PlateCarree(),
    #            colors = "white",
    #            linewidths=2,
    #)

    _ax.coastlines()
    _ax.set_global()

    _ax.set_title(title_text)
    _ax.set_extent((cent_lon - 60, cent_lon + 60, 0, 65), ccrs.PlateCarree())

    cb = plt.colorbar(mappable, orientation="horizontal", ticks=sst_ticks) 
    cb.ax.set_xlabel(cb_label)

for _ax in ax[len(all_products):]:

    fig.delaxes(_ax)


if not args.no_display:
    plt.show()

if args.output_dir != "":
    
    print("Create dir: %s" % (args.output_dir,))
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    fig.savefig("%s/SST_compare_%s.png" % (args.output_dir, data_date.strftime("%Y-%m-%d")), dpi=200)


