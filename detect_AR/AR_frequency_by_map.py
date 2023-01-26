import numpy as np
import load_data
import netCDF4
from datetime import (datetime, timedelta, timezone)
import traceback
import anomalies
import date_tools, fmon_tools, domain_tools, NK_tools
from pathlib import Path
import argparse

parser = argparse.ArgumentParser(
                    prog = 'plot_skill',
                    description = 'Plot prediction skill of GFS on AR.',
)

parser.add_argument('--beg-year', type=int, help='Date string: yyyy-mm-dd', required=True)
parser.add_argument('--end-year', type=int, help='Date string: yyyy-mm-dd', required=True)
parser.add_argument('--output',   type=str, help='Output file', default="")
parser.add_argument('--lat-rng', type=float, nargs=2, help='Latitude  range', required=True)
parser.add_argument('--lon-rng', type=float, nargs=2, help='Longitude range. 0-360', required=True)
parser.add_argument('--mask', type=str, help='Mask file. Land=0, Ocean=1', required=True)

args = parser.parse_args()

print(args)


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


ERA5_varnames = ["IWV", "IVT", "IWVKE", "sst", "dTdt"]
ignored_months = [4, 5, 6, 7, 8, 9]

lat_idx = None
lon_idx = None

lat = None
lon = None
wgt = None

cnt_boxes_settings = [
    ('Oct', (10, 1), (10, 31)),
    ('Nov', (11, 1), (11, 30)),
    ('Dec', (12, 1), (12, 31)),
    ('Jan', ( 1, 1), (1,  31)),
    ('Feb', ( 2, 1), (2,  28)),
    ('Mar', ( 3, 1), (3,  31)),
]


cnt_boxes = []

for cnt_boxes_setting in enumerate(cnt_boxes_settings):

    cnt_boxes.append({
        'name' : cnt_boxes_setting[0],
        'wd_beg' : watertime_tools.getWaterday(datetime(1, cnt_boxes_setting[1][0], cnt_boxes_setting[1][1]), no_leap=True),
        'wd_end' : watertime_tools.getWaterday(datetime(1, cnt_boxes_setting[2][0], cnt_boxes_setting[2][1]), no_leap=True),
        'total_cnt': 0,
    })


data = {
}


with netCDF4.Dataset(args.mask, "r") as ds:
       
    print("Loading mask file: %s" % (args.mask,))
    landsea_mask = np.ma.array(ds.variables["mask"][:], keep_mask=False)  # 1=ocean, 0=land
    mask = (1 - landsea_mask).astype(bool) # in masked_array, unused grid are denoted with "true"

for d, _t in enumerate(t_vec):
        
    _data = {}
            
    if _t.month in ignored_months or (_t.month == 2 and _t.day == 29):
        print("We do not need this time of data.")
        
    else:


        # Load ERA5 IVT data

        try:
            info = load_data.getFileAndIndex("ERA5", _t, root_dir="data", varname="IVT")

            print("Load `IVT` from file: %s" % ( info['filename'], ))

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


                
                IVT = ds.variables["IVT"][info['idx'], lat_idx, lon_idx]
                IVT_cond_met_idx = IVT >= args.IVT_threshold

        except Exception as e:

            print(traceback.format_exc()) 
            raise Error("Fail to load the day: %s" % (str(_t),))

        


        # Load ERA5 data
        for i, varname in enumerate(ERA5_varnames):

            if not (varname in data):

                data[varname] = {
                    'cnt' : [ np.zeros((len(lat), len(lon))) for _ in range(len(cnt_boxes))],
                    'avg' : [ np.zeros((len(lat), len(lon))) for _ in range(len(cnt_boxes))],
                    'std' : [ np.zeros((len(lat), len(lon))) for _ in range(len(cnt_boxes))],
                }

            try:

                load_varname = varname

                if varname == "dTdt":

                    info_l = load_data.getFileAndIndex("ERA5", _t + timedelta(days=-1), root_dir="data", varname=varname)
                    info_r = load_data.getFileAndIndex("ERA5", _t + timedelta(days=1), root_dir="data", varname=varname)

                    print("Load `%s` from file: %s" % (load_varname, info_l['filename']))
                    with netCDF4.Dataset(info_l['filename'], "r") as ds:
                        _data_l = ds.variables[load_varname][info_l['idx'], lat_idx, lon_idx]

                    print("Load `%s` from file: %s" % (load_varname, info_r['filename']))
                    with netCDF4.Dataset(info_r['filename'], "r") as ds:
                        _data_r = ds.variables[load_varname][info_r['idx'], lat_idx, lon_idx]

                    tmp = (_data_r - _data_l) / (2 * 86400.0)

                else:
                    # Load observation (the 'truth')
                    info = load_data.getFileAndIndex("ERA5", _t, root_dir="data", varname=varname)

                    print("Load `%s` from file: %s" % ( varname, info['filename'] ))

                    with netCDF4.Dataset(info['filename'], "r") as ds:
                        
                        tmp = ds.variables[load_varname][info['idx'], lat_idx, lon_idx]
                    
                for k, cnt_box in enumerate(cnt_boxes):

                    if _t >= cnt_box['wd_beg'] and _t <= cnt_box['wd_end']:

                        data[varname]['cnt'][k][IVT_cond_met_idx] += 1
                        data[varname]['avg'][k][IVT_cond_met_idx] += tmp[IVT_cond_met_idx]
                        data[varname]['std'][k][IVT_cond_met_idx] += tmp[IVT_cond_met_idx]**2

            except Exception as e:

                print(traceback.format_exc()) 
                raise Error("Fail to load the day: %s" % (str(_t),))




for varname in ERA5_varnames:

    _data = data[varname]

    for k, cnt_box in enumeratt(cnt_boxes):

        cnt = _data['cnt'][k]
        _data['avg'][k] /= cnt
        _data['std'][k] = np.sqrt((_data['std'][k][IVT_cond_met_idx] - cnt * tmp[IVT_cond_met_idx]**2) / (cnt - 1))




print("Output file: %s" % (args.output,))
with netCDF4.Dataset(args.output, 'w', format='NETCDF4') as ds:

    dim_lat = ds.createDimension('lat', len(lat))
    dim_lon = ds.createDimension('lon', len(lon))
    dim_box = ds.createDimension('box', len(cnt_boxes))

    var_lat = ds.createVariable('lat', 'f4', ('lat',))
    var_lon = ds.createVariable('lon', 'f4', ('lon',))
  
    var_lat[:] = lat 
    var_lon[:] = lon

    for varname in ERA5_varnames:
       
        for subvarname in ['avg', 'std', 'cnt']: 
            _var = ds.createVariable("%s_%s" % (varname, subvarname), 'f4', ('box', 'lat', 'lon'))
            for k in range(len(cnt_boxes)):
                _var[k, :, :] = data[varname][subvarname][k]
        
 
