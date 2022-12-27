import numpy as np
import netCDF4
import AR_tools, NK_tools
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
parser.add_argument('--output', type=str, help='Output file', default="")
parser.add_argument('--AR-dt-rng', type=float, nargs=2, help='Days of AR duration to do linear regression.', default=[0.0, 30.0])
parser.add_argument('--no-display', action="store_true")
args = parser.parse_args()
print(args)

data = {
    "ttl" : {},
    "clim" : {},
    "anom" : {},
}

data_dim = {}

AR_varnames = ["IWV", "IVT", "IWVKE", "sst", "mslhf", "msshf", "msnlwrf", "msnswrf", "mtpr", "mer", "mvimd", "t2m", "u10", "v10", "U", "MLD"]

zeta = 23.0

#AR_varnames = ["IVT", "sst", "mslhf", "msshf"]

with netCDF4.Dataset(args.input, "r") as ds:
 
    for varname in ["time", "time_clim",]:
        t = ds.variables[varname][:]
        data_dim[varname] = [ datetime.fromtimestamp(_t) for _t in t]

    for k, subdata in data.items():

        for AR_varname in AR_varnames:
            subdata[AR_varname] = ds.variables["%s_%s" % (AR_varname, k)][:]

        print(subdata.keys())    
        subdata['net_sfc_hf']  = subdata['msnswrf'] + subdata['msnlwrf'] + subdata['msshf'] + subdata['mslhf']
        subdata['net_conv_wv'] = subdata['mtpr'] + subdata['mer'] + subdata['mvimd']
   
        wb_prime = - ( subdata['msnlwrf'] + subdata['msshf'] + subdata['mslhf'] )
        subdata['Delta'] = NK_tools.calDelta(subdata["MLD"], subdata["U"], wb_prime, zeta, subdata['msnswrf']) 



def findfirst(a):
    return np.where(a)[0][0]

def findlast(a):
    return np.where(a)[0][-1]

def within(a, m, M):

    return m <= a and a < M

    #data['time_clim'] = data['time_clim'][:]
    #data['time'] = data['time'][:]

data['ttl']['IVT'][np.isnan(data['ttl']['IVT'])] = 0.0


AR_t_segs, AR_t_inds = AR_tools.detectAbove(data_dim['time'], data['ttl']['IVT'], 250.0, glue_threshold=timedelta(hours=24))


for i, AR_t_seg in enumerate(AR_t_segs):
    print("[%d] : %s to %s" % (i, AR_t_seg[0].strftime("%Y-%m-%d"), AR_t_seg[1].strftime("%Y-%m-%d")))
   

AR_evts = []
for k, t_seg in enumerate(AR_t_segs):
    print("Processing the %d-th AR time segment." % (k, ))
    AR_evt = {}

    ind = AR_t_inds[k, :] == True
    ind_first = findfirst(ind)
    ind_last  = findlast(ind)

    #print("(%d, %d) " % (ind_first, ind_last,))

    sst  = data['ttl']['sst']
    time = data_dim['time']
    AR_evt['dT']   = sst[ind_last] - sst[ind_first] 
    AR_evt['dt']   = (time[ind_last] - time[ind_first]).total_seconds()
    AR_evt['dTdt'] = AR_evt['dT'] / AR_evt['dt']
    
    AR_evt['U']   = np.mean(data['ttl']['U'][ind])
    AR_evt['MLD']   = np.mean(data['ttl']['MLD'][ind])
    
    AR_evt['net_sfc_hf']   = np.mean(data['ttl']['net_sfc_hf'][ind])
    AR_evt['net_conv_wv']   = np.mean(data['ttl']['net_conv_wv'][ind])

    
    AR_evt['dTdt_sfchf'] = AR_evt['net_sfc_hf'] / (3996*1026 * AR_evt['MLD'])
    AR_evt['dTdt_no_sfchf'] = AR_evt['dTdt'] - AR_evt['dTdt_sfchf']
    AR_evt['dTdt_ratio_from_sfchf']    = AR_evt['dTdt_sfchf'] / AR_evt['dTdt']
    AR_evt['dTdt_ratio_from_no_sfchf'] = 1.0 - AR_evt['dTdt_ratio_from_sfchf'] 

    AR_evt['t2m']   = np.mean(data['ttl']['t2m'][ind])
    AR_evt['ao_Tdiff']  = np.mean(data['ttl']['t2m'][ind] - data['ttl']['sst'][ind])
    
    AR_evt['U*ao_Tdiff']  = np.mean( (data['ttl']['t2m'][ind] - data['ttl']['sst'][ind]) * data['ttl']['U'][ind])
    
    AR_evt['Delta']   = np.mean(data['ttl']['Delta'][ind])
    
    AR_evt['mean_IVT'] = np.mean(data['ttl']['IVT'][ind])
    AR_evt['max_IVT']  = np.amax(data['ttl']['IVT'][ind])
    

    
    if AR_evt['dt'] == 0:
        AR_evt = None

    else:

        if within(AR_evt['dt'], args.AR_dt_rng[0]*86400, args.AR_dt_rng[1]*86400):
            AR_evt['do_linregress'] = True
        else:
            AR_evt['do_linregress'] = False


    AR_evts.append(AR_evt)


def collectData(AR_evts, varnames, ignorenan=True):

    data = {}

    if ignorenan:
        _AR_evts = []
        for i, AR_evt in enumerate(AR_evts):

            if AR_evt is not None:
                _AR_evts.append(AR_evt)


    AR_evts = _AR_evts

    if len(AR_evts) == 0:
        raise Exception("No available data")
    
    for k, varname in varnames.items():
        data[k] = np.zeros((len(AR_evts),))
        data[k][:] = np.nan


    for i, AR_evt in enumerate(AR_evts):

        for k, varname in varnames.items():
            
            if AR_evt is None:
                continue
            
            data[k][i] = AR_evt[varname]


    return data

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
from scipy.stats import linregress

print("done")


var_infos = {

    'dTdt' : {
        'var'  : "$\dot{T}_\\mathrm{ttl}$",
        'unit' : "$ \\mathrm{T} / \\mathrm{s} $",
    },

    'dTdt_sfchf' : {
        'var'  : "$\dot{T}_\\mathrm{hf}$",
        'unit' : "$ \\mathrm{T} / \\mathrm{s} $",
    },

    'dTdt_no_sfchf' : {
        'var'  : "$\dot{T}_\\mathrm{res}$",
        'unit' : "$ \\mathrm{T} / \\mathrm{s} $",
    },

    'U' : {
        'var'  : "$\\left|\\vec{U}_{10}\\right|$",
        'unit' : "$ \\mathrm{m} / \\mathrm{s} $",
    },

    'Delta' : {
        'var'  : "$ h\\Delta b w_e$",
        'unit' : "$ \\mathrm{m}^3 / \\mathrm{s}^3 $",
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

def plot_linregress(ax, X, Y, eq_x=0.1, eq_y=0.9, transform=None):

    if transform is None:
        transform = ax.transAxes

    res = linregress(X, Y)

    X_min, X_max = np.amin(X), np.amax(X)

    regressline_X = np.array([X_min, X_max])
    regressline_Y = res.slope * regressline_X + res.intercept

    ax.plot(regressline_X, regressline_Y, 'k--')
    
    ax.text(eq_x, eq_y, "$ y = %.2ex + %.2e $\n$R=%.2f$" % (res.slope, res.intercept, res.rvalue), transform=transform, ha="left", va="top", size=8)


    print("Number of data points: %d" % (len(X),))

#plot_data = [
#    ('dTdt', 'dTdt_sfchf'),    ('dTdt_sfchf', 'dTdt_no_sfchf'), ('U',     'dTdt_no_sfchf'), ('dt',      'mean_IVT'),
#    ('dTdt', 'dTdt_no_sfchf'), None,                            ('Delta', 'dTdt_no_sfchf'), ('max_IVT', 'mean_IVT')
#]

plot_data = [
    ('dTdt', 'dTdt_sfchf'),    ('dTdt_sfchf', 'dTdt_no_sfchf'), ('U',     'dTdt_no_sfchf'), 
    ('dTdt', 'dTdt_no_sfchf'), None,                            ('Delta', 'dTdt_no_sfchf'), 
]


rows = 2


if len(plot_data) % rows != 0:
    cols = len(plot_data) // rows + 1
else:    
    cols = len(plot_data) // rows


fig, ax = plt.subplots(rows, cols, figsize=(6*cols, 5*rows), squeeze=False, gridspec_kw = dict(hspace=0.3, wspace=0.4))

ax_flat = ax.flatten()
    

test_data = collectData(AR_evts, dict(picked='do_linregress'))

fig.suptitle("AR duration time range to do linear regression: %d - %d days\n# of cases = %d" % (args.AR_dt_rng[0], args.AR_dt_rng[1], np.sum(test_data['picked']==True)))

for i, _plot_data in enumerate(plot_data):

    if _plot_data is None:
        continue

    var_info_x = var_infos[_plot_data[0]]
    var_info_y = var_infos[_plot_data[1]]

    _ax = ax_flat[i]
    _data = collectData(AR_evts, dict(X=_plot_data[0], Y=_plot_data[1], Z='dt', picked='do_linregress'))
    mappable =  _ax.scatter(_data['X'], _data['Y'], c=_data['Z']/86400, s=10, cmap='bone_r', vmin=0, vmax=5)
    plot_linregress(_ax, _data['X'][_data['picked']==True], _data['Y'][_data['picked']==True])


    _ax.set_xlabel("%s [ %s ]" % (var_info_x['var'], var_info_x['unit']))
    _ax.set_ylabel("%s [ %s ]" % (var_info_y['var'], var_info_y['unit']))
    
    cbar = plt.colorbar(mappable, ax=_ax, orientation='vertical')
    cbar.ax.set_ylabel("AR duration [ day ]")

for i, _ax in enumerate(ax_flat):
    
    if i > len(plot_data)-1 or plot_data[i] is None: 
        fig.delaxes(_ax)


if not args.no_display:
    plt.show()

if args.output != "":
    print("Output figure: %s" % args.output)
    fig.savefig(args.output, dpi=200)

