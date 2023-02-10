import numpy as np
import netCDF4
import AR_tools
import anomalies

import traceback
from pathlib import Path
import argparse
from datetime import (timedelta, datetime, timezone)

parser = argparse.ArgumentParser(
                    prog = 'plot_skill',
                    description = 'Plot prediction skill of GFS on AR.',
)

parser.add_argument('--input', type=str, help='Input file', required=True)
parser.add_argument('--no-display', action="store_true")
args = parser.parse_args()
print(args)

data = {
    "ttl" : {},
    "clim" : {},
    "anom" : {},
}
AR_varnames = ["IWV", "IVT", "IWVKE", "sst", "mslhf", "msshf", "msnlwrf", "msnswrf", "mtpr", "mvimd", "t2m", "u10", "v10", "U"]

AR_varnames = ["IVT", "sst", "U", "mtpr", ]#"mer"]

with netCDF4.Dataset(args.input, "r") as ds:
    
    for AR_varname in AR_varnames:
        for k, subdata in data.items():
            subdata[AR_varname] = ds.variables["%s_%s" % (AR_varname, k)][:]

    for varname in ["time", "time_clim",]:
        t = ds.variables[varname][:]
        data[varname] = [ datetime.fromtimestamp(_t) for _t in t]

    #data['time_clim'] = data['time_clim'][:]
    #data['time'] = data['time'][:]

data['ttl']['IVT'][np.isnan(data['ttl']['IVT'])] = 0.0


AR_t_segs, AR_t_inds = AR_tools.detectAbove(data['time'], data['ttl']['IVT'], 250.0, glue_threshold=timedelta(hours=24))


# Find correct indices
for i, AR_t_seg in enumerate(AR_t_segs):
    print("[%d] : %s to %s" % (i, AR_t_seg[0].strftime("%Y-%m-%d"), AR_t_seg[1].strftime("%Y-%m-%d")))
    
    
print("List detected AR events:")
for i, AR_t_seg in enumerate(AR_t_segs):
    print("[%d] : %s to %s" % (i, AR_t_seg[0].strftime("%Y-%m-%d"), AR_t_seg[1].strftime("%Y-%m-%d")))


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

print("done")


var_infos = {

    'IWV' : {
        'label' : "[ $ \\mathrm{kg} / \\mathrm{m}^2 $ ]",
    },

    'IVT' : {
        'label' : "[ $ \\mathrm{kg} \\, \\mathrm{m} / \\mathrm{s} / \\mathrm{m}^2 $ ]",
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



fig, ax = plt.subplots(len(AR_varnames), 1, figsize=(12, 4*len(AR_varnames)), sharex=True)

for i, AR_varname in enumerate(AR_varnames):

    print("Plot %s" % (AR_varname,))
    _ax = ax[i]

    trans = transforms.blended_transform_factory(_ax.transData, _ax.transAxes)
        
    for k, t_seg in enumerate(AR_t_segs):

        ind = AR_t_inds[k, :] == True
       
        dt = t_seg[1] - t_seg[0] 
        if dt.total_seconds() < 86400*3:
            continue

        _IVT = data['ttl']['IVT'][ind]

        _time = np.array(data["time"])[ind]
        _ind0 = np.argmax(_IVT)
        _time0 = _time[_ind0]

        _time_aligned_day = np.array([ _t.total_seconds() for _t in (_time - _time0) ])/86400.0

        if AR_varname in ["IVT",]:

            _ax.plot(_time_aligned_day, data['ttl'][AR_varname][ind], color='gray', linestyle='solid', linewidth=1)

        else: 
            #_ax.plot(_time_aligned_day, data['ttl'][AR_varname][ind], color='gray', linestyle='solid', linewidth=1)
            
            _ax.plot(_time_aligned_day, data['ttl'][AR_varname][ind], color='gray', linestyle='solid', linewidth=1)
            

  
    #_ax.set_title("Forecast error of %s" % (AR_var, ))
    #_ax.set_ylabel("%s\n%s" % (AR_var, var_infos[AR_var]['label']))
    _ax.set_ylabel("%s" % (AR_varname,))
    #_ax.xaxis.set_major_formatter(DateFormatter('%Y-%m-%d'))

    _ax.legend()

if not args.no_display:
    plt.show()

if args.output != "":
    
    fig.savefig(args.output, dpi=200)
























