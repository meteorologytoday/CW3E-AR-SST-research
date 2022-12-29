import numpy as np
import netCDF4
import AR_tools, NK_tools
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
parser.add_argument('--output', type=str, help='Output file', required=True)
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
    
    mid_time = time[ind_first] + (time[ind_last] - time[ind_first]) / 2
    AR_evt['mid_time'] = mid_time.timestamp()
    AR_evt['month'] = mid_time.month
    AR_evt['year'] = mid_time.year + 1 if within(mid_time.month, 10, 12) else mid_time.year
    
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

    if AR_evt is not None:
        AR_evts.append(AR_evt)

    
AR_database = pd.DataFrame.from_records(AR_evts)

print("Output file: %s" % (args.output_file,))
AR_database.to_csv(args.output_file)
