import numpy as np
import netCDF4
import AR_tools, NK_tools, fmon_tools, watertime_tools
import anomalies
import pandas as pd
import traceback
from datetime import (timedelta, datetime, timezone)

def findfirst(a):
    return np.where(a)[0][0]

def findlast(a):
    return np.where(a)[0][-1]

def within(a, m, M):

    return m <= a and a < M

def constructCases(filename, IVT_threshold):

    data = {
        "ttl" : {},
        "clim" : {},
        "anom" : {},
    }

    data_dim = {}

    AR_varnames = ["IWV", "IVT", "IWVKE", "sst", "mslhf", "msshf", "msnlwrf", "msnswrf", "mtpr", "mer", "mvimd", "t2m", "u10", "v10", "U", "MLD", "dT", "db", "dTdt", "dTdt_sfchf", "dTdt_no_sfchf", "w_deepen", "dTdt_deepen", "net_sfc_hf", "net_conv_wv", "vort10", "curltau", "dTdt_Ekman"]


    with netCDF4.Dataset(filename, "r") as ds:
     
        for varname in ["time", "time_clim",]:
            t = ds.variables[varname][:]
            data_dim[varname] = [ datetime(1970, 1, 1) + _t * timedelta(days=1) for _t in t]

        for k, subdata in data.items():

            for AR_varname in AR_varnames:
                subdata[AR_varname] = ds.variables["%s_%s" % (AR_varname, k)][:]


    data['ttl']['IVT'][np.isnan(data['ttl']['IVT'])] = 0.0

    AR_t_segs, AR_t_inds = AR_tools.detectAbove(data_dim['time'], data['ttl']['IVT'], IVT_threshold, glue_threshold=timedelta(hours=24))


    for i, AR_t_seg in enumerate(AR_t_segs):
        print("[%d] : %s to %s" % (i, AR_t_seg[0].strftime("%Y-%m-%d"), AR_t_seg[1].strftime("%Y-%m-%d")))
       

    AR_evts = []
    for k, t_seg in enumerate(AR_t_segs):

        print("Processing the %d-th AR time segment." % (k, ))

        AR_evt = {}

        ind = AR_t_inds[k, :] == True
        ind_first = findfirst(ind)
        ind_last  = findlast(ind)

        #sst  = data['ttl']['sst']
        time = data_dim['time']
        AR_evt['dt']   = (time[ind_last] - time[ind_first]).total_seconds()
        AR_evt['dTdt'] = np.mean(data['ttl']['dTdt'][ind])
        
        mid_time = time[ind_first] + (time[ind_last] - time[ind_first]) / 2
        AR_evt['datetime_beg'] = time[ind_first]
        AR_evt['datetime_end'] = time[ind_last]
        AR_evt['mid_time'] = mid_time.timestamp()
        AR_evt['month'] = mid_time.month
        AR_evt['year'] = mid_time.year + 1 if within(mid_time.month, 10, 12) else mid_time.year
        AR_evt['watertime'] = watertime_tools.getWatertime(datetime.fromtimestamp(AR_evt['mid_time']))
        AR_evt['wateryear'] = np.floor(AR_evt['watertime'])
        AR_evt['waterdate'] = AR_evt['watertime'] - AR_evt['wateryear']
        
        AR_evt['vort10']   = np.mean(data['ttl']['vort10'][ind])
        AR_evt['curltau']   = np.mean(data['ttl']['curltau'][ind])
        AR_evt['U']   = np.mean(data['ttl']['U'][ind])
        AR_evt['MLD']   = np.mean(data['ttl']['MLD'][ind])
        
        AR_evt['u10']   = np.mean(data['ttl']['u10'][ind])
        AR_evt['v10']   = np.mean(data['ttl']['v10'][ind])
        
        AR_evt['U_mean']   = (AR_evt['u10']**2 + AR_evt['v10']**2)**0.5
        
        AR_evt['net_sfc_hf']   = np.mean(data['ttl']['net_sfc_hf'][ind])
        AR_evt['net_conv_wv']   = np.mean(data['ttl']['net_conv_wv'][ind])

        AR_evt['dTdt_sfchf']    = np.mean(data['ttl']['dTdt_sfchf'][ind])
        AR_evt['dTdt_no_sfchf'] = np.mean(data['ttl']['dTdt_no_sfchf'][ind])
        
        AR_evt['dTdt_Ekman']    = np.mean(data['ttl']['dTdt_Ekman'][ind])

        AR_evt['t2m']   = np.mean(data['ttl']['t2m'][ind])
        AR_evt['ao_Tdiff']  = np.mean(data['ttl']['t2m'][ind] - data['ttl']['sst'][ind])
        
        AR_evt['U*ao_Tdiff']  = np.mean( (data['ttl']['t2m'][ind] - data['ttl']['sst'][ind]) * data['ttl']['U'][ind])
        
        AR_evt['dTdt_deepen'] = np.mean(data['ttl']['dTdt_deepen'][ind])
        AR_evt['w_deepen']    = np.mean(data['ttl']['w_deepen'][ind])
        
        AR_evt['mean_IVT'] = np.mean(data['ttl']['IVT'][ind])
        AR_evt['max_IVT']  = np.amax(data['ttl']['IVT'][ind])
        
        if AR_evt['dt'] == 0:
            AR_evt = None


        """
        else:

            if within(AR_evt['dt'], args.AR_dt_rng[0]*86400, args.AR_dt_rng[1]*86400):
                AR_evt['do_linregress'] = True
            else:
                AR_evt['do_linregress'] = False
        """

        AR_evts.append(AR_evt)


    _AR_evts = []
    for AR_evt in AR_evts:
        if AR_evt is not None:
            _AR_evts.append(AR_evt)

    AR_evts = _AR_evts
    _AR_evts = None

    # Convert to a dataframe
    colnames = AR_evts[0].keys()
    database = { colname : [] for colname in colnames }

    for i, AR_evt in enumerate(AR_evts):
        
        for colname in colnames:
            database[colname].append(AR_evt[colname]) 

    
    df = pd.DataFrame.from_dict(database, orient='columns')
    #df.to_csv(output_database)

    df.attrs["wateryear_beg"] = watertime_tools.getWateryear(data_dim["time"][0])
    df.attrs["wateryear_end"] = watertime_tools.getWateryear(data_dim["time"][-1])

    return df
        
