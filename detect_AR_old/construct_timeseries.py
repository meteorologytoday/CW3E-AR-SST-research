import numpy as np
import load_data
import netCDF4
from datetime import (datetime, timedelta, timezone)
import traceback
import anomalies
import date_tools, fmon_tools
from pathlib import Path
import argparse

parser = argparse.ArgumentParser(
                    prog = 'plot_skill',
                    description = 'Plot prediction skill of GFS on AR.',
)

parser.add_argument('--beg-year', type=int, help='Date string: yyyy-mm-dd', required=True)
parser.add_argument('--end-year', type=int, help='Date string: yyyy-mm-dd', required=True)
parser.add_argument('--output-dir', type=str, help='Output directory', default="")
parser.add_argument('--products', type=str, nargs="+", help='Forcast products.', default=["GFS",])
parser.add_argument('--lat-rng', type=float, nargs=2, help='Latitude  range', required=True)
parser.add_argument('--lon-rng', type=float, nargs=2, help='Longitude range. 0-360', required=True)
parser.add_argument('--mld', type=str, help='The mixed layer specifier.', choices=['somxl010', 'somxl030'], required=True)
parser.add_argument('--save', action="store_true", help='Save the timeseries into an npz file')


args = parser.parse_args()

print(args)

# Configuration

beg_date = datetime(args.beg_year-1, 10,  1 )
end_date = datetime(args.end_year,    4,  1 )

total_days = (end_date - beg_date).days


print("Beg: ", beg_date)
print("End: ", end_date)
print("Total days: ", total_days)

if total_days <= 0:
    raise Exception("No days are avaiable.")



AR_varnames = ["IWV", "IVT", "IWVKE", "sst", "mslhf", "msshf", "msnlwrf", "msnswrf", "mtpr", "mer", "mvimd", "t2m", "u10", "v10"]

aug_AR_varnames = AR_varnames + ['U', 'MLD', 'hdb', 'dT', 'db']

data_good = np.ones((total_days, len(aug_AR_varnames)))

ts = {}

for AR_var in aug_AR_varnames:
    ts[AR_var] = np.zeros((total_days,)) 

# npdatetime is more like an array of numbers

t_vec = [ beg_date + timedelta(days=d) for d in range(total_days) ]
t_vec_npdatetime = np.array(t_vec, dtype="datetime64[s]")


# Prep

lat = None
lon = None
lat_idx = None
lon_idx = None
wgt = None

# Load MLD data
beg_fmon = fmon_tools.datetime2fmon(beg_date) - 1
end_fmon = fmon_tools.datetime2fmon(end_date) + 1
total_fmon = end_fmon - beg_fmon + 1

fmon_data = {}
fmon_varnames = ["MLD", "db", "T_upper", "T_lower"]  # loaded from file
aug_fmon_varnames = fmon_varnames + ["hdb", "dT"]         # generate from above
for varname in aug_fmon_varnames:
    fmon_data[varname] = np.zeros((total_fmon,)) 



for k, _fmon in enumerate(range(beg_fmon, end_fmon+1)):
        
    y, m = fmon_tools.getYearMonthFromfmon(_fmon)
    
    _t = datetime(y, m, 1)
    
    _data = {}
            
    for i, varname in enumerate(fmon_varnames):

        try:

            load_varname = varname

            info = load_data.getFileAndIndex("ORA5", _t, root_dir="data", varname=varname, mxl_algo=args.mld)

            print("Load file: ", info['filename'])
            
            with netCDF4.Dataset(info['filename'], "r") as ds:

                # decide range
                if lat is None:
                    
                    lat = ds.variables[info['varnames']['lat']][:]
                    lon = ds.variables[info['varnames']['lon']][:] % 360.0

                    lat_rng = np.array(args.lat_rng)
                    lon_rng = np.array(args.lon_rng) % 360
                    
                    if lat_rng[1] < lat_rng[0]:  # across lon=0
                        raise Exception("Latitude range should be lat_min, lat_max")

                    lat_idx = (lat_rng[0] < lat) & (lat < lat_rng[1])

                    if lon_rng[1] >= lon_rng[0]:
                        lon_idx = (lon_rng[0] < lon) & (lon < lon_rng[1])
                    
                    else:  # across lon=0
                        print("Across lon=0")
                        lon_idx = (lon_rng[0]) < lon | (lon < lon_rng[1])

                    lat = lat[lat_idx]
                    lon = lon[lon_idx]
                    wgt = np.cos(lat * np.pi / 180)


                _data[load_varname] = ds.variables[load_varname][info['idx'], lat_idx, lon_idx]

        except Exception as e:

            print(traceback.format_exc()) 
            print("Someting wrong happened when loading date: %s" % (_t.strftime("%Y-%m-%d"),))

            raise e


    _data['dT'] = _data['T_upper'] - _data['T_lower']
    _data['hdb'] = _data['MLD'] * _data['db']
    for varname in aug_fmon_varnames:
        fmon_data[varname][k] = np.average(np.average( _data[varname], axis=1), weights=wgt)

# Interpolate ORA5 data
ts_fmon = np.zeros((total_fmon,), dtype='datetime64[s]')
for i, _fmon in enumerate(range(beg_fmon, end_fmon+1)):
    y, m = fmon_tools.getYearMonthFromfmon(_fmon)
    ts_fmon[i] = np.datetime64(datetime(y, m, 15, 12, 0, 0))

for varname in aug_fmon_varnames:
    ts[varname] = np.interp(t_vec_npdatetime.astype(np.float64), ts_fmon.astype(np.float64), fmon_data[varname]) 


# Load ERA5 data
for d in range(total_days):
        
    _t = beg_date + timedelta(days=d)
        
            
    _data = {}
            
    keep_going = True

    for i, varname in enumerate(AR_varnames):

        try:

            if _t.month in [4, 5, 6, 7, 8, 9]:
                raise Exception("We do not this time of data.")

            load_varname = varname

            # Load observation (the 'truth')
            info = load_data.getFileAndIndex("ERA5", _t, root_dir="data", varname=varname)

            print("Load file: ", info['filename'])


            with netCDF4.Dataset(info['filename'], "r") as ds:
                # decide range
                if lat is None:
                    
                    lat = ds.variables[info['varnames']['lat']][:]
                    lon = ds.variables[info['varnames']['lon']][:] % 360.0


                    lat_rng = np.array(args.lat_rng)
                    lon_rng = np.array(args.lon_rng) % 360
                    
                    if lat_rng[1] < lat_rng[0]:  # across lon=0
                        raise Exception("Latitude range should be lat_min, lat_max")

                    lat_idx = (lat_rng[0] < lat) & (lat < lat_rng[1])

                    if lon_rng[1] >= lon_rng[0]:
                        lon_idx = (lon_rng[0] < lon) & (lon < lon_rng[1])
                    
                    else:  # across lon=0
                        print("Across lon=0")
                        lon_idx = (lon_rng[0]) < lon | (lon < lon_rng[1])

                    lat = lat[lat_idx]
                    lon = lon[lon_idx]
                    wgt = np.cos(lat * np.pi / 180)


                _data[load_varname] = ds.variables[load_varname][info['idx'], lat_idx, lon_idx]
                #ts[varname][d] = np.average(np.average( _data[load_varname], axis=1), weights=wgt)



        except Exception as e:

            print(traceback.format_exc()) 
            print("Someting wrong happened when loading date: %s" % (_t.strftime("%Y-%m-%d"),))

            data_good[d, i] = 0.0

            keep_going = False




    if keep_going:
        _data['U'] = np.sqrt(_data['u10']**2 + _data['v10']**2)
        for varname in aug_AR_varnames:
            if varname in aug_fmon_varnames:
                print("Already have %s. Skip" % (varname,))
                continue

            ts[varname][d] = np.average(np.average( _data[varname], axis=1), weights=wgt)




   

print("Exclude non-consecutive years")
data_all_good = np.prod(data_good, axis=1) != 0
all_good_t = t_vec_npdatetime[data_all_good]
missing_dates = date_tools.findMissingDatetime(all_good_t, beg_date, end_date, timedelta(days=1))


# I let the months [ (y-1).11, (y-1).12, y.1, y.2, y.3, y.4 ] be the winter of the year `y`. 
rm_years = []
needed_missing_dates = np.zeros((len(missing_dates),), dtype=bool)
for i, missing_date in enumerate(missing_dates):

    if missing_date.month in [10, 11, 12]:
        
        rm_years.append(missing_date.year+1)
        needed_missing_dates[i] = True


    elif missing_date.month in [1, 2, 3]:
        
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
    
data_good[marked_to_remove == True, :] = 0.0
    
print("Filtering data")
for i, AR_var in enumerate(aug_AR_varnames):
    print("Doing: ", AR_var)
    data_good_idx = data_good[:, i] == 0.0
    ts[AR_var][data_good_idx] = np.nan
 




# produce decomposed data
data_decompose = {
    'clim' : {},
    'anom' : {},
}

for AR_varname in aug_AR_varnames:
    time_clim, data_decompose['clim'][AR_varname], data_decompose['anom'][AR_varname], cnt = anomalies.decomposeClimAnom(t_vec, ts[AR_varname])


if args.output_dir != "":
            
    print("Output directory: %s" % (args.output_dir,))
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    with netCDF4.Dataset("%s/AR_timeseries_climanom.nc" % (args.output_dir,), 'w', format='NETCDF4') as ds:

        dim_time_clim = ds.createDimension('time', len(t_vec))
        var_time = ds.createVariable('time', 'f4', ('time',))
        var_time[:] = t_vec_npdatetime
        
        dim_time_clim = ds.createDimension('time_clim', len(time_clim))
        var_time_clim = ds.createVariable('time_clim', 'f4', ('time_clim',))
        var_time_clim.units    = 'days since 0001-01-01 00:00:00'
        var_time_clim.calendar = 'noleap'
        var_time_clim[:] = time_clim
        
        for AR_varname in aug_AR_varnames:

            _var_ttl     = ds.createVariable("%s_ttl" % (AR_varname, ), 'f4', ('time',))
            _var_ttl[:]  = ts[AR_varname]
 
            _var_clim = ds.createVariable("%s_clim" % (AR_varname, ), 'f4', ('time_clim',))
            _var_clim[:] = data_decompose['clim'][AR_varname]
 
            _var_anom = ds.createVariable("%s_anom" % (AR_varname, ), 'f4', ('time',))
            _var_anom[:] = data_decompose['anom'][AR_varname]
            

    with netCDF4.Dataset("%s/AR_timeseries.nc" % (args.output_dir,), 'w', format='NETCDF4') as ds:

        dim_time = ds.createDimension('time', len(t_vec_npdatetime))
        var_time = ds.createVariable('time', 'f4', ('time',))
        var_time[:] = t_vec_npdatetime
        
        for varname in aug_AR_varnames:

            _var = ds.createVariable(varname, 'f4', ('time',))
            #value.units = 'Unknown'

            _var[:] = ts[varname]


    with open("%s/rm_years.txt" % (args.output_dir,), "w") as f:
        for y in rm_years:
            f.write("%d\n" % y)           



