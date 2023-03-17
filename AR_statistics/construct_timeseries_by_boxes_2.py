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
import map_divide_tools

import xarray as xr
import ECCO_helper

def weightedAvg(var_data, wgts):

    d = var_data.to_numpy()

    idx = np.isfinite(d)
    d = d[idx]
    w = wgts.to_numpy()[idx]

    return np.sum(d * w) / np.sum(w)


print("Loading libraries completed.")

parser = argparse.ArgumentParser(
                    prog = 'plot_skill',
                    description = 'Plot prediction skill of GFS on AR.',
)

parser.add_argument('--beg-year',   type=int, help='Date string: yyyy-mm-dd', required=True)
parser.add_argument('--end-year',   type=int, help='Date string: yyyy-mm-dd', required=True)
parser.add_argument('--output-dir', type=str, help='Output directory', default="")
parser.add_argument('--output-filename', type=str, help='Output filename', default="full_dataset.nc")
parser.add_argument('--lat-rng',    type=float, nargs=2, help='Latitude  range', required=True)
parser.add_argument('--lon-rng',    type=float, nargs=2, help='Longitude range. 0-360', required=True)
parser.add_argument('--lat-nbox',   type=int, help='Latitude  range', required=True)
parser.add_argument('--lon-nbox',   type=int, help='Longitude range. 0-360', required=True)
parser.add_argument('--mask-ERA5',  type=str, help='mask file of ERA5', required=True)
parser.add_argument('--mask-ECCO',  type=str, help='mask file of ECCO', required=True)
parser.add_argument('--ignore-empty-box',  action="store_true")

args = parser.parse_args()

print(args)

# Configuration

beg_date = datetime(args.beg_year-1, 10,  1 )
end_date = datetime(args.end_year,    4,  1 )


total_days = (end_date - beg_date).days

t_vec = [ beg_date + timedelta(days=d) for d in range(total_days) ]
t_vec_npdatetime = np.array(t_vec, dtype="datetime64[s]")

lat_rng = np.array(args.lat_rng)
lon_rng = np.array(args.lon_rng) % 360

print("Beg: ", beg_date)
print("End: ", end_date)
print("Total days: ", total_days)

if total_days <= 0:
    raise Exception("No days are avaiable.")


lon_bnds = np.array([ lon_rng[0] + (lon_rng[1] - lon_rng[0]) / args.lon_nbox * i for i in range(args.lon_nbox+1)])
lat_bnds = np.array([ lat_rng[0] + (lat_rng[1] - lat_rng[0]) / args.lat_nbox * i for i in range(args.lat_nbox+1)])

boxes = map_divide_tools.makeDividedBoxes(lon_bnds, lat_bnds)

#print("### List of divided boxes: ")
#for box in boxes:
#    print(box)

ERA5_varnames = ["IWV", "IVT", "IWVKE", "sst", "mslhf", "msshf", "msnlwrf", "msnswrf", "u10", "v10", "t2m", "mtpr", "sst"]
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
computed_ERA5_vars = []

all_varnames = ERA5_varnames + ECCO_varnames + computed_LLC_vars + computed_ERA5_vars

print("Construct output datasets: ts_datasets")
ts_datasets = [

    xr.Dataset(
        { 
            varname : (['time',], np.zeros((total_days,), dtype=np.float64)) 
            for varname in all_varnames 
        },

        coords = {
            'time' : t_vec_npdatetime,
        },
    ) 
    
    for b in range(len(boxes))

]

print("Construct output datasets: full_dataset")

full_dataset = xr.Dataset(
    { 
        varname : (['time', 'lat', 'lon', ], np.zeros((total_days, len(lat_bnds)-1, len(lon_bnds)-1), dtype=np.float64)) 
        for varname in (ERA5_varnames + ECCO_varnames + computed_LLC_vars + computed_ERA5_vars) 
    },

    coords = {
        'time' : t_vec_npdatetime,
        'lat'  : (lat_bnds[:-1] + lat_bnds[1:]) / 2,
        'lon'  : (lon_bnds[:-1] + lon_bnds[1:]) / 2,
    },
) 

box_number = np.zeros((len(lat_bnds)-1, len(lon_bnds)-1), dtype=np.int32)

for b, box in enumerate(boxes):
    j = box["j"]
    i = box["i"]
    box_number[j, i] = box["n"]

full_dataset = xr.merge([
    
    full_dataset,
    
    xr.Dataset(
        { 
            "box_number" : (['lat', 'lon'], box_number),
        },
        coords = {
            'lat'  : (lat_bnds[:-1] + lat_bnds[1:]) / 2,
            'lon'  : (lon_bnds[:-1] + lon_bnds[1:]) / 2,
            'lat_bnds' : lat_bnds,
            'lon_bnds' : lon_bnds,
        },
    ),
])


full_dataset.coords["time"].encoding["units"] = "days since 1990-01-01"

# Eventually, data_good will be merge with each dataset of each box
data_good = xr.DataArray(
    name = "data_good",
    data =  np.zeros((total_days,), dtype=np.int32),
    dims=["time",],
    coords=dict(time=t_vec_npdatetime),
)

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
    

ditch_this_wateryear = np.nan
current_wateryear = np.nan

print("Ready to process data.")

for d, _t in enumerate(t_vec):

    print("# Processing date: ", _t)
        
    _data = {}
            
    I_have_all_data_for_today = True
    
    if _t.month in ignored_months:
        data_good[d] = 0
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

                # Create label 
                for b, box in enumerate(boxes):

                    box['ERA5_subset_idx'] = (
                          ( ERA5_lat >= box['polygon']['lat_bnds'][0] )
                        & ( ERA5_lat <  box['polygon']['lat_bnds'][1] )
                        & ( ERA5_lon >= box['polygon']['lon_bnds'][0] )
                        & ( ERA5_lon <  box['polygon']['lon_bnds'][1] )
                        & ( mask_ERA5 == 1)
                    )

                    box['empty_ERA5'] = np.sum(box['ERA5_subset_idx']) == 0

                    if box['empty_ERA5']:
                        if args.ignore_empty_box:
                            print("[ERA5] Ignore empty box: %d" % (b,))
                        else:
                            raise Exception("ERROR: No point is selected in ERA5 latlon grid.")

                    #box['ERA5_wgts'] = ERA5_grid.wgts.where(box['ERA5_subset_idx'], other=0.0).rename("ERA5_wgts")
                    box['ERA5_wgts'] = ERA5_grid.wgts.to_numpy()[box['ERA5_subset_idx']]


            # Subset after
            #_var = _var.where(ERA5_subset_idx)

            _data[load_varname] = _var.load()

        except Exception as e:

            print(traceback.format_exc()) 
            print("Someting wrong happened when loading date: %s" % (_t.strftime("%Y-%m-%d"),))

            I_have_all_data_for_today = False


    """
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
    """
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

            ds_ECCO = ds_ECCO.astype(np.float64)

            if ecco_grid is None:
                
                mask_ECCO = xr.open_dataset(args.mask_ECCO).mask.to_numpy()

                ecco_grid = ECCO_helper.getECCOGrid()

                ecco_lat = ecco_grid["YC"]
                ecco_lon = ecco_grid["XC"] % 360


                for b, box in enumerate(boxes):
                
                    poly = box['polygon']

                    box['ecco_subset_idx'] = (
                        (ecco_lat   >= poly['lat_bnds'][0])
                        & (ecco_lat <  poly['lat_bnds'][1])
                        & (ecco_lon >= poly['lon_bnds'][0])
                        & (ecco_lon <  poly['lon_bnds'][1])
                        & (mask_ECCO == 1)
                    )
                    
                    box['empty_ecco'] = np.sum(box['ecco_subset_idx']) == 0

                    if box['empty_ecco']:
                        if args.ignore_empty_box:
                            print("[ECCO] Ignore empty box: %d" % (b,))
                        else:
                            raise Exception("ERROR: No point is selected in ECCO LLC grid.")

                    
                    #box['ecco_wgts'] = ecco_grid.rA.where(box['ecco_subset_idx'], other=0.0)
                    box['ecco_wgts'] = ecco_grid.rA.to_numpy()[box['ecco_subset_idx']]


            # subset afterwards   
            #ds_ECCO = ds_ECCO.where(ecco_subset_idx)
            _data[varname] = ds_ECCO[varname].load()

    except Exception as e:

        print(traceback.format_exc()) 
        print("ECCO: Someting wrong happened when loading date: %s" % (_t.strftime("%Y-%m-%d"),))

        I_have_all_data_for_today = False


    if I_have_all_data_for_today:
        data_good[d] = 1

    else:
        data_good[d] = 0
        print("Missing data for date: ", _t)
        ditch_this_wateryear = current_wateryear
        continue

    # Add other vairables inside
    print("Do magical extension")
    magicalExtension(_data)


    print("Do average of each box")
    # Make average of each box
    for varname, var_data in _data.items():
        #print("varname: ", varname)
        for b, box in enumerate(boxes):

            if box['empty_ecco'] or box['empty_ERA5']:
                continue

            ts_ds = ts_datasets[b]
            if (varname in ERA5_varnames) or (varname in computed_ERA5_vars):
                subset_idx_varname  = 'ERA5_subset_idx'
                subset_wgts_varname = 'ERA5_wgts'
                
            elif (varname in ECCO_varnames) or (varname in computed_LLC_vars):
                subset_idx_varname  = 'ecco_subset_idx'
                subset_wgts_varname = 'ecco_wgts'

            else:
                raise Exception("Unknown variable : %s" % (varname,) )

            #_masked_data = var_data.where(box[subset_idx_varname])
            #_wgts        = box[subset_wgts_varname]  # already subsetted
            #ts_ds[varname][d] = weightedAvg(_masked_data, _wgts)
 
            idx = box[subset_idx_varname]
            _masked_data = var_data.to_numpy()[idx]
            _wgts        = box[subset_wgts_varname]  # already subsetted
            ts_ds[varname][d] = np.sum(_masked_data * _wgts) / np.sum(_wgts)
            


            
    print("Do recheck of MLG budget")
    for b in range(len(boxes)):

        if box['empty_ecco'] or box['empty_ERA5']:
            continue

        ts_ds = ts_datasets[b]

        MLG_recheck = (ts_ds['dMLTdt'][d] - (
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
        )).rename('MLG_recheck')

        print("[box=%d][v=%s] Double check residue using averaged values: " % (b, varname), ts_ds["MLG_residue"].data[d])



print("Exclude non-consecutive years")
data_good_t = t_vec_npdatetime[ data_good == 1 ]
missing_dates = date_tools.findMissingDatetime(data_good_t, beg_date, end_date, timedelta(days=1))


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
    
data_good[marked_to_remove == True] = 0
 

print("Clear the data of the date if we do not have all the data of that date, or if any date in that winter is missing.")




print("Merge data_good into each dataset")
data_good_idx = data_good == 1
none_is_selected_idx = np.isnan(data_good) # This is an all-false boolean array. It is used to assign an entire dataset as NaN for empty boxes
for i, ts_ds in enumerate(ts_datasets):

    print(boxes[i].keys())

    if boxes[i]['empty_ERA5'] or boxes[i]['empty_ecco']:
        ts_datasets[i] = ts_ds.where(none_is_selected_idx).merge(data_good)
    else: 
        ts_datasets[i] = ts_ds.where(data_good_idx).merge(data_good)

print("Merge each timeseries into a complete file (experimental).")
for b, box in enumerate(boxes):
    
    i = box["i"]
    j = box["j"]
    for varname in all_varnames:
        full_dataset[varname][:, j, i] = ts_datasets[b][varname]





if args.output_dir != "":
            
    print("Output directory: %s" % (args.output_dir,))
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
   
    """ 
    for b in range(len(boxes)):

        box = boxes[b]
        ts_ds = ts_datasets[b]

        output_filename = "%s/AR_timeseries_b%d.nc" % (args.output_dir, b)

        print("[b=%d] Output filename: %s" % ( b, output_filename, ))
        ts_ds.to_netcdf(output_filename)
    """

    output_filename_full_dataset = "%s/%s" % (args.output_dir, args.output_filename)
        
    print("Output filename: %s" % ( output_filename_full_dataset, ))
    full_dataset.to_netcdf(
        output_filename_full_dataset,
        unlimited_dims=["time",]
    )
 
    

    """
    with open("%s/box_info.txt" % (args.output_dir,), "w") as f:

        f.write("Latitude  range: [%f, %f) \n" % (args.lat_rng[0],     args.lat_rng[1])) 
        f.write("Longitude range: [%f, %f) \n" % (args.lon_rng[0]%360, args.lon_rng[1]%360)) 
        
        for b, box in enumerate(boxes):
            f.write("box=%d, label=%s, lat_min=%.2f, lat_max=%.2f, lon_min=%.2f, lon_max=%.2f, nsamples_atm=%d, nsamples_ocn=%d \n" % (
                b,
                box['label'],
                box['polygon']['lat_bnds'][0],
                box['polygon']['lat_bnds'][1],
                box['polygon']['lon_bnds'][0],
                box['polygon']['lon_bnds'][1],
                np.sum(box['ERA5_subset_idx'].astype(int)),
                np.sum(box['ecco_subset_idx'].astype(int)),
            ))


    """

    """
    with open("%s/rm_years.txt" % (args.output_dir,), "w") as f:
        for y in rm_years:
            f.write("%d\n" % y)           
    """


