with open("shared_header.py", "rb") as source_file:
    code = compile(source_file.read(), "shared_header.py", "exec")
exec(code)

import ARdetection
import numpy as np
import xarray as xr
import pandas as pd
import argparse
from earth_constants import r_E

parser = argparse.ArgumentParser(
                    prog = 'make_ERA5_AR_objects.py',
                    description = 'Postprocess ECCO data (Mixed-Layer integrated).',
)

parser.add_argument('--method', required=True, type=str, choices=["ANOM_LEN", "TOTIVT250"])
parser.add_argument('--AR-clim-dir', required=True, type=str)
parser.add_argument('--leftmost-lon', type=float, help='The leftmost longitude on the map. It matters when doing object detection finding connectedness.', default=90.0)
parser.add_argument('--nproc', type=int, default=1)
args = parser.parse_args()
print(args)


input_dir = "data/ERA5/AR_processed"
output_dir = "data/ERA5/ARobjs_%s" % (args.method,) 
input_file_prefix = "ERA5_AR"
input_clim_file_prefix = "ERA5_AR"
output_file_prefix = "ERA5_ARobjs"


beg_time_str = beg_time.strftime("%Y-%m-%d")
end_time_str = end_time.strftime("%Y-%m-%d")
dts = pd.date_range(beg_time_str, end_time_str, freq="D", inclusive="both")


def myARFilter(AR_obj):

    result = True
    
    length_min = 1000e3
    area_min = 1000e3 * 100e3

    if (AR_obj['length'] < length_min) or AR_obj['area'] < area_min:
        result = False
    
    return result



def doJob(dt):

    jobname = dt.strftime("%Y-%m-%d")

    try: 
        
        # Load file
        AR_clim_file = "%s/%s_%s.nc" % (args.AR_clim_dir, input_clim_file_prefix, dt.strftime("%m-%d"))
        AR_full_file = "%s/%s_%s.nc" % (input_dir,        input_file_prefix,      dt.strftime("%Y-%m-%d"))
        output_file = "%s/%s_%s.nc"  % (output_dir, output_file_prefix, dt.strftime("%Y-%m-%d"))

        # Test if file exists

        if os.path.isfile(output_file):
            print("[%s] Warning: File %s already exists. Skip this." % (jobname,))
            return None

        
            
        print("[%s] Making %s with input files: %s, %s." % (jobname, output_file, AR_full_file, AR_clim_file))

        # Make output dir
        Path(os.path.dirname(output_file)).mkdir(parents=True, exist_ok=True)

        # Load anom
        ds_full = xr.open_dataset(AR_full_file)
        ds_clim = xr.open_dataset(AR_clim_file)

        # find the lon=0
        lon_first_zero = np.argmax(ds_full.coords["lon"].to_numpy() >= args.leftmost_lon)

        ds_full = ds_full.roll(lon=-lon_first_zero, roll_coords=True)
        ds_clim = ds_clim.roll(lon=-lon_first_zero, roll_coords=True)
        
        lat = ds_full.coords["lat"].to_numpy() 
        lon = ds_full.coords["lon"].to_numpy()  % 360
      
        # For some reason we need to reassign it otherwise the contourf will be broken... why??? 
        ds_full = ds_full.assign_coords(lon=lon) 
        ds_clim = ds_clim.assign_coords(lon=lon) 
        
        IVT_full = ds_full.IVT[0, :, :].to_numpy()
        IVT_clim = ds_clim.IVT[0, :, :].to_numpy()
        IVT_anom = IVT_full - IVT_clim

        llat, llon = np.meshgrid(lat, lon, indexing='ij')

        dlat = np.deg2rad((lat[0] - lat[1]))  # ERA5 data is werid. Latitude start at north and move southward
        dlon = np.deg2rad((lon[1] - lon[0]))

        area = r_E**2 * np.cos(np.deg2rad(llat)) * dlon * dlat

        print("[%s] Compute AR_objets using method: %s" % (jobname, args.method))

        if args.method == "ANOM_LEN": 

            labeled_array, AR_objs = ARdetection.detectARObjects(
                IVT_anom, llat, llon, area,
                IVT_threshold=250.0,
                weight=IVT_full,
                filter_func = myARFilter,
            )

        elif args.method == "TOTIVT250": 
 
            labeled_array, AR_objs = ARdetection.detectARObjects(
                IVT_full, llat, llon, area,
                IVT_threshold=250.0,
                weight=IVT_full,
                filter_func = myARFilter,
            )
          
        # Convert AR object array into Dataset format
        Nobj = len(AR_objs)
        data_output = dict(
            feature_n    = np.zeros((Nobj,), dtype=int),
            area         = np.zeros((Nobj,),),
            length       = np.zeros((Nobj,),),
            centroid_lat = np.zeros((Nobj,),),
            centroid_lon = np.zeros((Nobj,),),
        )

        for i, AR_obj in enumerate(AR_objs):
            data_output['feature_n'][i] = AR_obj['feature_n']
            data_output['area'][i]      = AR_obj['area']
            data_output['length'][i]    = AR_obj['length']
            data_output['centroid_lat'][i] = AR_obj['centroid'][0]
            data_output['centroid_lon'][i] = AR_obj['centroid'][1]
         
        # Make Dataset
        ds_out = xr.Dataset(

            data_vars=dict(
                map          = (["time", "lat", "lon"], np.reshape(labeled_array, (1, *labeled_array.shape))),
                feature_n    = (["ARobj", ], data_output["feature_n"]),
                length       = (["ARobj", ], data_output["length"]),
                area         = (["ARobj", ], data_output["area"]),
                centroid_lat = (["ARobj", ], data_output["centroid_lat"]),
                centroid_lon = (["ARobj", ], data_output["centroid_lon"]),
            ),

            coords=dict(
                lon=(["lon"], lon),
                lat=(["lat"], lat),
                time=(["time"], [dt,]),
            ),

            attrs=dict(description="AR objects file."),
        )

        ds_out.to_netcdf(output_file, encoding={'time': {'dtype': 'i4'}})

    except Exception as e:

        print("[%s] Error. Now print stacktrace..." % (jobname,))
        import traceback
        traceback.print_exc()


    return output_file




with Pool(processes=args.nproc) as pool:

    it = pool.imap(doJob, dts)


    for result in it:
        print("Task for file %s is compelete." % (result,))


print("Done.")


