import numpy as np
import load_data
import netCDF4
from datetime import (datetime, timedelta, timezone)
import traceback
import anomalies
import date_tools, fmon_tools, domain_tools, NK_tools, KPP_tools, watertime_tools
import earth_constants as ec
from pathlib import Path
import argparse
import buoyancy_linear

import xarray as xr
import ECCO_helper

def weightedAvg(var_data, wgts):

    d = var_data.to_numpy()

    idx = np.isfinite(d)
    d = d[idx]
    w = wgts.to_numpy()[idx]

    print("Weight: ", np.mean(d))

    return np.mean(d)
    #return np.sum(d * w) / np.sum(w)


print("Loading libraries completed.")

parser = argparse.ArgumentParser(
                    prog = 'plot_skill',
                    description = 'Plot prediction skill of GFS on AR.',
)

parser.add_argument('--beg-year', type=int, help='Date string: yyyy-mm-dd', required=True)
parser.add_argument('--end-year', type=int, help='Date string: yyyy-mm-dd', required=True)
parser.add_argument('--output-dir', type=str, help='Output directory', default="")
parser.add_argument('--lat-rng', type=float, nargs=2, help='Latitude  range', required=True)
parser.add_argument('--lon-rng', type=float, nargs=2, help='Longitude range. 0-360', required=True)
parser.add_argument('--mask-ERA5', type=str, help='mask file of ERA5', required=True)
parser.add_argument('--mask-ECCO', type=str, help='mask file of ECCO', required=True)
#parser.add_argument('--mask', type=str, help='Mask file. Land=0, Ocean=1', required=True)

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


ERA5_varnames = ["IWV", "IVT", "IWVKE", "sst", "mslhf", "msshf", "msnlwrf", "msnswrf", ]
ECCO_varnames = [
    "dMLTdt",
    "MLT",
    "MXLDEPTH",
    "MLG_ttl",
    "MLG_frc_sw",
    "MLG_frc_lw",
    "MLG_frc_sh",
    "MLG_frc_lh",
    "MLG_frc_fwf",
    "MLG_hadv",
    "MLG_vadv",
    "MLG_hdiff",
    "MLG_vdiff",
    "MLG_ent",
    "MLG_rescale",
    "MLD",
    "dTdz_b",
    "MLU",
    "MLV",
    "U_g",
    "V_g",
    "dMLTdx",
    "dMLTdy",
]
            
ignored_months = [4, 5, 6, 7, 8, 9]


tendency_residue_tolerance = 1e-10

domain_check_tolerance = 1e-10
ERA5_lat_raw = None
ERA5_lon_raw = None

ecco_grid = None

lat = None
lon = None
f_co = None

computed_LLC_vars  = ["MLG_geo", "MLG_ageo", "dTdz_b_over_h", "MLG_residue"]
computed_ERA5_vars = ["ERA5_MLG_ttl",]
#        'ERA5_MLG_ttl', 'ERA5_MLG_frc', 'ERA5_sfc_hf', 'ERA5_MLG_ttl_exp', 'ERA5_MLG_ttl_uexp'
#]


ts_ds = xr.Dataset(
    { 
        varname : (['time',], np.zeros((total_days,), dtype=np.float64)) 
        for varname in (ERA5_varnames + ECCO_varnames + computed_LLC_vars + computed_ERA5_vars + ['data_good',] ) 
    },

    coords = {
        'time' : t_vec_npdatetime, 
    },

)

lat_rng = np.array(args.lat_rng)
lon_rng = np.array(args.lon_rng) % 360

def magicalExtension(_data):
    
    #_data['ERA5_sfc_hf']  = _data['msnswrf'] + _data['msnlwrf'] + _data['msshf'] + _data['mslhf']
    #_data['ERA5_MLG_ttl_exp']  = _data['ERA5_sfc_hf'] / (3996*1026 * _data['MLD'])
    #_data['ERA5_MLG_ttl_uexp'] = _data['ERA5_MLG_ttl'] - _data['ERA5_MLG_frc']
    
    _data["MLG_geo"]  = - ( _data["MLU"] * _data["dMLTdx"] + _data["MLV"] * _data["dMLTdy"] )
    _data["MLG_ageo"] = - ( (_data["MLU"] - _data["U_g"]) * _data["dMLTdx"] + (_data["MLV"] - _data["V_g"]) * _data["dMLTdy"] )
    _data["dTdz_b_over_h"] = _data["dTdz_b"] / _data["MLD"]
    
    _data['MLG_residue'] = _data['dMLTdt'] - (
          _data['MLG_frc_sw']
        + _data['MLG_frc_lw']
        + _data['MLG_frc_sh']
        + _data['MLG_frc_lh']
        + _data['MLG_frc_fwf']
        + _data['MLG_rescale']
        + _data['MLG_hadv']
        + _data['MLG_vadv']
        + _data['MLG_hdiff']
        + _data['MLG_vdiff']
        + _data['MLG_ent']
    )
    
    res = _data["MLG_residue"].to_numpy()
    res_max = np.amax(np.abs(res[np.isfinite(res)]))
    print("Max of abs(MLG_residue): ", res_max)
    

#mask_ERA5 = ds


ditch_this_wateryear = np.nan
current_wateryear = np.nan

for d, _t in enumerate(t_vec):

    print("# Processing date: ", _t)
        
    _data = {}
            
    I_have_all_data_for_today = True
    
    if _t.month in ignored_months:
        ts_ds.data_good[d] = 0
        print("We do not need this time of data: ", _t)
        continue
        

    current_wateryear = watertime_tools.getWateryear(_t)

    if np.isfinite(ditch_this_wateryear):

        if ditch_this_wateryear == current_wateryear:
            print("Ditch this wateryear: ", _t)
            continue

        elif ditch_this_wateryear == current_wateryear + 1: # Just moved into the next water year. Reset the flag.
            ditch_this_wateryear = np.nan
        
        else:
            raise Exception("Wrong counting of the year. Please check")

        

    # Load ERA5 data
    for i, varname in enumerate(ERA5_varnames):

        try:

            load_varname = varname

            # Load observation (the 'truth')
            info = load_data.getFileAndIndex("ERA5", _t, root_dir="data", varname=varname)

            print("Load `%s` from file: %s" % ( varname, info['filename'] ))


            ds_ERA5 = xr.open_dataset(info["filename"])
            _var = ds_ERA5[varname].isel(time=0)

            if ERA5_lat_raw is None:
              
                print("Coordinate loading...")

                mask_ERA5 = xr.open_dataset(args.mask_ERA5).mask.to_numpy()


                ERA5_lat_raw = ds_ERA5.coords["lat"]
                ERA5_lon_raw = ds_ERA5.coords["lon"] % 360

                ERA5_lat, ERA5_lon = np.meshgrid(ERA5_lat_raw.to_numpy(), ERA5_lon_raw.to_numpy(), indexing='ij')

                ERA5_wgts = np.cos(ERA5_lat * np.pi / 180)

                ERA5_subset_idx = (
                      ( ERA5_lat >= lat_rng[0] )
                    & ( ERA5_lat <= lat_rng[1] )
                    & ( ERA5_lon >= lon_rng[0] )
                    & ( ERA5_lon <= lon_rng[1] )
                    & ( mask_ERA5 == 1)
                ) 

                ERA5_grid = xr.Dataset(
                    { 
                        "llat" : (['lat', 'lon'], ERA5_lat), 
                        "llon" : (['lat', 'lon'], ERA5_lon), 
                        "wgts" : (['lat', 'lon'], ERA5_wgts),
                    },

                    coords = {
                        'lat' : ERA5_lat_raw,
                        'lon' : ERA5_lon_raw,
                    },
                )

                ERA5_wgts = ERA5_grid.wgts.where(ERA5_subset_idx, other=0.0)



            _var = _var.where(ERA5_subset_idx)

            """
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

            """
            _data[load_varname] = _var.load()

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
        _var_l = (xr.open_dataset(info_l["filename"])[varname]).isel(time=0).where(ERA5_subset_idx)
        
        print("Load `%s` from file: %s" % (load_varname, info_r['filename']))
        _var_r = (xr.open_dataset(info_r["filename"])[varname]).isel(time=0).where(ERA5_subset_idx)

        dvardt = (_var_r - _var_l) / (2 * 86400.0)
        _data['ERA5_MLG_ttl'] = (_var_r - _var_l) / (2 * 86400.0)

    except Exception as e:

        print(traceback.format_exc()) 
        print("Someting wrong happened when loading date: %s" % (_t.strftime("%Y-%m-%d"),))

        I_have_all_data_for_today = False

    #del info

    ############ Loading ECCOv4 data ############

    try:

        for varname in ECCO_varnames:

            ecco_filename = ECCO_helper.getECCOFilename(varname, "DAILY", _t)
            ecco_filename = "data/ECCO_LLC/%s/%s" % ecco_filename

            print("Load `%s` from file: %s" % ( varname, ecco_filename, ))

            if varname == "MLD":
                ds_ECCO = xr.open_dataset(ecco_filename).isel(time_snp=0)
                
            else:
                ds_ECCO = xr.open_dataset(ecco_filename).isel(time=0)

            if ecco_grid is None:
                
                mask_ECCO = xr.open_dataset(args.mask_ECCO).mask.to_numpy()

                ecco_grid = ECCO_helper.getECCOGrid()

                ecco_lat = ecco_grid["YC"]
                ecco_lon = ecco_grid["XC"] % 360
                
                ecco_subset_idx = (
                    (ecco_lat >= lat_rng[0])
                    & (ecco_lat <= lat_rng[1])
                    & (ecco_lon <= lon_rng[0])
                    & (ecco_lon <= lon_rng[1])
                    & (mask_ECCO == 1)
                )

                ecco_wgts = ecco_grid.rA.where(ecco_subset_idx, other=0.0)
                
            ds_ECCO = ds_ECCO.where(ecco_subset_idx)
            _data[varname] = ds_ECCO[varname].load()

    except Exception as e:

        print(traceback.format_exc()) 
        print("ECCO: Someting wrong happened when loading date: %s" % (_t.strftime("%Y-%m-%d"),))

        I_have_all_data_for_today = False


    if I_have_all_data_for_today:
        ts_ds['data_good'][d] = 1

    else:
        ts_ds['data_good'][d] = 0
        print("Missing data for date: ", _t)
        ditch_this_wateryear = current_wateryear
        continue

    # Add other vairables inside
    magicalExtension(_data)


    for varname, var_data in _data.items():

        if (varname in ERA5_varnames) or (varname in computed_ERA5_vars): 
            ts_ds[varname][d] = weightedAvg(var_data, ERA5_wgts)
            
        elif (varname in ECCO_varnames) or (varname in computed_LLC_vars):
            ts_ds[varname][d] = weightedAvg(var_data, ecco_wgts)

        else:
            raise Exception("Unknown variable : %s" % (varname,) )

    print("Averaged Residue: ", ts_ds["MLG_residue"][d])

    MLG_res2 = (ts_ds['dMLTdt'][d] - (
          ts_ds['MLG_frc_sw'][d]
        + ts_ds['MLG_frc_lw'][d]
        + ts_ds['MLG_frc_sh'][d]
        + ts_ds['MLG_frc_lh'][d]
        + ts_ds['MLG_frc_fwf'][d]
        + ts_ds['MLG_rescale'][d]
        + ts_ds["MLG_vdiff"][d]
        + ts_ds["MLG_hdiff"][d]
        + ts_ds["MLG_vadv"][d]
        + ts_ds["MLG_hadv"][d]
        + ts_ds["MLG_ent"][d]
    )).rename('MLG_res2')
    
    print("Averaged Residue - forced check: ", MLG_res2)



print("Exclude non-consecutive years")
data_good_t = t_vec_npdatetime[ ts_ds.data_good == 1 ]
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
    
ts_ds.data_good[marked_to_remove == True] = 0
 


print("Clear the data of the date if we do not have all the data of that date, or if any date in that winter is missing.")

ts_ds = ts_ds.where(ts_ds.data_good == 1)
#data_not_good_idx = (ts_ds.data_good ==)
#for varname, var_data in ts.items():
#    var_data[data_not_good_idx] = np.nan

# produce decomposed data
#data_decompose = {
#    'clim' : {},
#    'anom' : {},
#}

#for varname in ts_ds.keys():
#    time_clim, data_decompose['clim'][varname], data_decompose['anom'][varname], cnt = anomalies.decomposeClimAnom(t_vec, ts[varname])


if args.output_dir != "":
            
    print("Output directory: %s" % (args.output_dir,))
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    
    output_filename = "%s/AR_timeseries.nc" % (args.output_dir,)

    print("Output filename: %s" % ( output_filename, ))
    ts_ds.to_netcdf(output_filename)

    with open("%s/rm_years.txt" % (args.output_dir,), "w") as f:
        for y in rm_years:
            f.write("%d\n" % y)           



