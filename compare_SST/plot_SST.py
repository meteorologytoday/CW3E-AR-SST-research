import numpy as np
import load_data
import netCDF4
from datetime import datetime
from scipy.interpolate import RectBivariateSpline

# Configuration
ref_product = "ERA5"
compare_products = [ "CFSR", "OSTIA", "OISST",]
data_date = datetime(2017, 1, 7)

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
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import cm
import cartopy.crs as ccrs
print("done")

cent_lon = 180.0
proj = ccrs.PlateCarree(central_longitude=cent_lon, globe=None)

fig, ax = plt.subplots(1, len(all_products), subplot_kw = dict(projection=proj), gridspec_kw=dict(width_ratios=[1]*len(all_products)), figsize=(16, 4))

fig.suptitle("Date: %s" % (data_date.strftime("%Y-%m-%d"),))

for i, product in enumerate(all_products):
 
    print("Plotting product: %s" % (product,) )

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
        
    _ax = ax[i]
    mappable = _ax.contourf(ref_data['lon'], ref_data['lat'], _plot_data, levels=sst_levs,
                transform=ccrs.PlateCarree(),
                cmap=cmap,
                extend="both")

    _ax.coastlines()
    _ax.set_global()

    _ax.set_title(title_text)
    _ax.set_extent((cent_lon - 60, cent_lon + 60, 0, 65), ccrs.PlateCarree())

    cb = plt.colorbar(mappable, orientation="horizontal", ticks=sst_ticks) 
    cb.ax.set_xlabel("SST [${}^\\circ\\mathrm{C}$]")

plt.show()


