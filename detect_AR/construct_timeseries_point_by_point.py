import numpy as np
import load_data
import netCDF4
from datetime import (datetime, timedelta, timezone)
import traceback
import anomalies
import date_tools, fmon_tools, domain_tools, NK_tools
import earth_constants as ec
from pathlib import Path
import argparse
import buoyancy_linear

parser = argparse.ArgumentParser(
                    prog = 'plot_skill',
                    description = 'Plot prediction skill of GFS on AR.',
)

parser.add_argument('--beg-year', type=int, help='Date string: yyyy-mm-dd', required=True)
parser.add_argument('--end-year', type=int, help='Date string: yyyy-mm-dd', required=True)
parser.add_argument('--output-dir', type=str, help='Output directory', default="")
parser.add_argument('--lat-rng', type=float, nargs=2, help='Latitude  range', required=True)
parser.add_argument('--lon-rng', type=float, nargs=2, help='Longitude range. 0-360', required=True)
parser.add_argument('--mask', type=str, help='Mask file. Land=0, Ocean=1', required=True)
parser.add_argument('--mld', type=str, help='The mixed layer specifier.', choices=['somxl010', 'somxl030'], required=True)

args = parser.parse_args()

print(args)

# Configuration

beg_date = datetime(args.beg_year-1, 10,  1 )
end_date = datetime(args.end_year,    4,  1 )

total_days = (end_date - beg_date).days
t_vec = [ beg_date + timedelta(days=d) for d in range(total_days) ]
t_vec_npdatetime = np.array(t_vec, dtype="datetime64[s]")



print("Beg: ", beg_date)
print("End: ", end_date)
print("Total days: ", total_days)

if total_days <= 0:
    raise Exception("No days are avaiable.")



ERA5_varnames = ["IWV", "IVT", "IWVKE", "sst", "mslhf", "msshf", "msnlwrf", "msnswrf", "mtpr", "mer", "mvimd", "t2m", "u10", "v10", "vort10", "curltau", ]
ORA5_varnames = ["MLD", "db", "dT"]
            
ignored_months = [4, 5, 6, 7, 8, 9]


domain_check_tolerance = 1e-10
ERA5_lat_raw = None
ERA5_lon_raw = None
ORA5_lat_raw = None
ORA5_lon_raw = None
lat_idx = None
lon_idx = None

lat = None
lon = None
wgt = None
f_co = None

computed_vars = ['U', 'db', 'dT', 'hdb', 'dTdt', 'w_deepen', 'dTdt_deepen', 'net_sfc_hf', 'net_conv_wv', 'sfhf_wosw', 'pme', 'dTdt_sfchf', 'dTdt_no_sfchf', 'dT_contribution_to_db', 'dTdt_Ekman', 'w_Ekman']

ts = { varname : np.zeros((total_days,), dtype=np.float32) 
    for varname in (ERA5_varnames + ORA5_varnames + computed_vars) 
}

data_good = np.zeros((total_days,), dtype=bool)


def magicalExtension(_data):

    #_data["MLD"] *= 2

    _data['U']   = np.sqrt(_data['u10']**2 + _data['v10']**2)
    #_data['dT']  = _data['T_upper'] - _data['T_lower']
    _data['hdb'] = _data['MLD'] * _data['db']

    #                           shortwave            longwave          sensible            latent
    _data['net_sfc_hf']  = _data['msnswrf'] + _data['msnlwrf'] + _data['msshf'] + _data['mslhf']
    _data['net_conv_wv'] = _data['mtpr'] + _data['mer'] + _data['mvimd']

    _data['sfhf_wosw'] = - ( _data['msnlwrf'] + _data['msshf'] + _data['mslhf'] ) # positive upwards
    _data['pme']       = _data['mtpr'] + _data['mer']                             # positive means raining
 
    _data['w_deepen'] = NK_tools.cal_we(_data["MLD"], _data["U"], _data['sfhf_wosw'], _data['pme'], _data['msnswrf'], _data['db']) 
    _data['dTdt_deepen'] = NK_tools.cal_dSSTdt_e(_data["MLD"], _data["U"], _data['sfhf_wosw'], _data['pme'], _data['msnswrf'], _data['db'], _data['dT'])
 
    #_data['DeltaOnlyU'] = NK_tools.calDeltaOnlyU(_data["MLD"], _data["U"])

    
    _data['dTdt_sfchf'] = _data['net_sfc_hf'] / (3996*1026 * _data['MLD'])
    _data['dTdt_no_sfchf'] = _data['dTdt'] - _data['dTdt_sfchf']


    _data['dT_contribution_to_db'] = buoyancy_linear.g0 * buoyancy_linear.alpha_T * _data['dT'] / _data['db']

    _data['w_Ekman']  = _data['curltau'] / f_co[:, np.newaxis] / ec.rho_sw
    _data['dTdt_Ekman']  = _data['w_Ekman'] * _data['dT'] / _data['MLD']

    _data['dTdt_Ekman'][_data['dTdt_Ekman'] < 0] = 0
 

with netCDF4.Dataset(args.mask, "r") as ds:
       
    print("Loading mask file: %s" % (args.mask,))
    landsea_mask = np.ma.array(ds.variables["mask"][:], keep_mask=False)  # 1=ocean, 0=land
    mask = (1 - landsea_mask).astype(bool) # in masked_array, unused grid are denoted with "true"
    #print("Sum of mask (): ", np.sum(mask.mask))

for d, _t in enumerate(t_vec):
        
    _data = {}
            
    I_have_all_data_for_today = True
    
    if _t.month in ignored_months:

        print("We do not need this time of data.")
        I_have_all_data_for_today = False
        
    else:

        # Load ERA5 data
        for i, varname in enumerate(ERA5_varnames):

            try:

                load_varname = varname

                # Load observation (the 'truth')
                info = load_data.getFileAndIndex("ERA5", _t, root_dir="data", varname=varname)

                print("Load `%s` from file: %s" % ( varname, info['filename'] ))


                with netCDF4.Dataset(info['filename'], "r") as ds:
                    
                    if ERA5_lat_raw is None:
                       
                        print("Coordinate loading...") 
                        ERA5_lat_raw = ds.variables[info['varnames']['lat']][:]
                        ERA5_lon_raw = ds.variables[info['varnames']['lon']][:] % 360.0

                        lat_rng = np.array(args.lat_rng)
                        lon_rng = np.array(args.lon_rng) % 360
                        
                        lat_idx, lon_idx, wgt = domain_tools.detectIndexRange(ERA5_lat_raw, ERA5_lon_raw, lat_rng, lon_rng)

                        lat = ERA5_lat_raw[lat_idx]
                        lon = ERA5_lon_raw[lon_idx]
                        wgt = np.cos(lat * np.pi / 180)

                        mask = mask[lat_idx, :][:, lon_idx]

                        print("Shape of mask: ", mask.shape)

                        f_co = 2 * ec.Omega * np.sin(np.pi / 180 * lat)


                    _data[load_varname] = ds.variables[load_varname][info['idx'], lat_idx, lon_idx]

            except Exception as e:

                print(traceback.format_exc()) 
                print("Someting wrong happened when loading date: %s" % (_t.strftime("%Y-%m-%d"),))

                I_have_all_data_for_today = False


        # Special: compute dSST/dt as variable "dTdt"
        try:

            varname = "sst"
            load_varname = varname

            # Load observation (the 'truth')
            info_l = load_data.getFileAndIndex("ERA5", _t + timedelta(days=-1), root_dir="data", varname=varname)
            info_r = load_data.getFileAndIndex("ERA5", _t + timedelta(days=1), root_dir="data", varname=varname)

            print("Load `%s` from file: %s" % (load_varname, info_l['filename']))
            with netCDF4.Dataset(info_l['filename'], "r") as ds:
                _data_l = ds.variables[load_varname][info_l['idx'], lat_idx, lon_idx]

            print("Load `%s` from file: %s" % (load_varname, info_r['filename']))
            with netCDF4.Dataset(info_r['filename'], "r") as ds:
                _data_r = ds.variables[load_varname][info_r['idx'], lat_idx, lon_idx]

            _data['dTdt'] = (_data_r - _data_l) / (2 * 86400.0)

        except Exception as e:

            print(traceback.format_exc()) 
            print("Someting wrong happened when loading date: %s" % (_t.strftime("%Y-%m-%d"),))

            I_have_all_data_for_today = False

        #del info

        ############ Loading ORA5 data ############ 
        # Determine boundaries
        if _t.day < 15:
            _t_l = fmon_tools.fmon2datetime(fmon_tools.datetime2fmon(_t)-1)
            _t_r = fmon_tools.fmon2datetime(fmon_tools.datetime2fmon(_t))
        else:
            _t_l = fmon_tools.fmon2datetime(fmon_tools.datetime2fmon(_t))
            _t_r = fmon_tools.fmon2datetime(fmon_tools.datetime2fmon(_t)+1)
        
        _t_l = _t_l.replace(day=15)
        _t_r = _t_r.replace(day=15)

        Delta_t_all = _t_r - _t_l
        Delta_t     = _t - _t_l


        for i, varname in enumerate(ORA5_varnames):

            try:

                load_varname = varname
                
                info_l = load_data.getFileAndIndex("ORA5", _t_l, root_dir="data", varname=varname, mxl_algo=args.mld)
                info_r = load_data.getFileAndIndex("ORA5", _t_r, root_dir="data", varname=varname, mxl_algo=args.mld)
                tmp_l = None
                tmp_r = None

                print("Load file [left]: ", info_l['filename'])
                with netCDF4.Dataset(info_l['filename'], "r") as ds:
                    tmp_l = ds.variables[load_varname][info_l['idx'], lat_idx, lon_idx]

                    if ORA5_lat_raw is None:
                        
                        ORA5_lat_raw = ds.variables[info_l['varnames']['lat']][:]
                        ORA5_lon_raw = ds.variables[info_l['varnames']['lon']][:] % 360.0

                        # Compare coordinate
                        if np.any(np.abs(ORA5_lat_raw - ERA5_lat_raw) > domain_check_tolerance) or np.any(np.abs(ORA5_lon_raw - ERA5_lon_raw) > domain_check_tolerance):
                            raise Error("Fatal error: ERA5 and ORA5 has different domain")


                print("Load file [right]: ", info_r['filename'])
                with netCDF4.Dataset(info_r['filename'], "r") as ds:
                    tmp_r = ds.variables[load_varname][info_r['idx'], lat_idx, lon_idx]
                    

                _data[varname] = tmp_l + (tmp_r - tmp_l) * (Delta_t / Delta_t_all)

            except Exception as e:

                print(traceback.format_exc()) 
                print("Someting wrong happened when loading date: %s" % (_t.strftime("%Y-%m-%d"),))

                raise e


        data_good[d] = I_have_all_data_for_today

    if not I_have_all_data_for_today:

        print("Missing data for date: ", _t)
        continue

    # Add other vairables inside
    magicalExtension(_data)


    for varname, var_data in _data.items():
        _tmp_data = np.ma.array(
            _data[varname],
            keep_mask=False,
            mask=mask,
        )

        ts[varname][d] = np.average(np.average( _tmp_data, axis=1), weights=wgt)


print("Exclude non-consecutive years")
data_good_t = t_vec_npdatetime[data_good == True]
missing_dates = date_tools.findMissingDatetime(data_good_t, beg_date, end_date, timedelta(days=1))


# I let the months [ (y-1).11, (y-1).12, y.1, y.2, y.3, y.4 ] be the winter of the year `y`. 
rm_years = []
needed_missing_dates = np.zeros((len(missing_dates),), dtype=bool)
for i, missing_date in enumerate(missing_dates):

    if missing_date.month in [11, 12]:#[10, 11, 12]:
        
        rm_years.append(missing_date.year+1)
        needed_missing_dates[i] = True


    elif missing_date.month in [1, 2]:#, 3]:
        
        rm_years.append(missing_date.year)
        needed_missing_dates[i] = True
        
rm_years = np.unique(rm_years)

ii=0
for i, missing_date in enumerate(missing_dates):
    if needed_missing_dates[i]:
        ii+=1
        print("[%d] Missing date needed: %s" % (ii, missing_date.strftime("%Y-%m-%d"),))


print("Original year range: %d-%d (%d years)." % (args.beg_year, args.end_year, args.end_year - args.beg_year + 1, ))

if len(rm_years) == 0:
    print("Congratulations! No missing data.")
else:
    print("I am going to remove %d years: %s" % ( len(rm_years), ",".join(["%d" % year for year in rm_years]), ))


marked_to_remove = np.ones((len(t_vec),), dtype=bool)
for i, t in enumerate(t_vec):
    marked_to_remove[i] = t.year in rm_years
    
data_good[marked_to_remove == True] = False
 


print("Clear the data of the date if we do not have all the data of that date, or if any date in that winter is missing.")
data_not_good_idx = (data_good == False)
for varname, var_data in ts.items():
    var_data[data_not_good_idx] = np.nan

# produce decomposed data
data_decompose = {
    'clim' : {},
    'anom' : {},
}

for varname in ts.keys():
    time_clim, data_decompose['clim'][varname], data_decompose['anom'][varname], cnt = anomalies.decomposeClimAnom(t_vec, ts[varname])


if args.output_dir != "":
            
    print("Output directory: %s" % (args.output_dir,))
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    with netCDF4.Dataset("%s/AR_timeseries_climanom.nc" % (args.output_dir,), 'w', format='NETCDF4') as ds:

        dim_time_clim = ds.createDimension('time', len(t_vec))
        var_time = ds.createVariable('time', 'f4', ('time',))
        var_time.units    = 'days since 1970-01-01 00:00:00'
        var_time.calendar = 'standard'
        var_time[:] = (t_vec_npdatetime - np.datetime64('1970-01-01')) / np.timedelta64(1, 'D')
        
        dim_time_clim = ds.createDimension('time_clim', len(time_clim))
        var_time_clim = ds.createVariable('time_clim', 'f4', ('time_clim',))
        var_time_clim.units    = 'days since 0001-01-01 00:00:00'
        var_time_clim.calendar = 'noleap'
        var_time_clim[:] = time_clim
        
        for varname in ts.keys():

            _var_ttl     = ds.createVariable("%s_ttl" % (varname, ), 'f4', ('time',))
            _var_ttl[:]  = ts[varname]
 
            _var_clim = ds.createVariable("%s_clim" % (varname, ), 'f4', ('time_clim',))
            _var_clim[:] = data_decompose['clim'][varname]
 
            _var_anom = ds.createVariable("%s_anom" % (varname, ), 'f4', ('time',))
            _var_anom[:] = data_decompose['anom'][varname]
            

    with netCDF4.Dataset("%s/AR_timeseries.nc" % (args.output_dir,), 'w', format='NETCDF4') as ds:

        dim_time = ds.createDimension('time', len(t_vec_npdatetime))
        var_time = ds.createVariable('time', 'f4', ('time',))
        var_time[:] = t_vec_npdatetime
        
        for varname in ts.keys():

            _var = ds.createVariable(varname, 'f4', ('time',))
            #value.units = 'Unknown'

            _var[:] = ts[varname]


    with open("%s/rm_years.txt" % (args.output_dir,), "w") as f:
        for y in rm_years:
            f.write("%d\n" % y)           



