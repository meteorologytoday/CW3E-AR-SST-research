from datetime import datetime
import os.path
import numpy as np
import scipy.interpolate
from scipy.interpolate import RectBivariateSpline


def interpolate_deprecated(lat, lon, arr, new_lat, new_lon, mask_idx=None):
    """This is an older version of interpolate 2D data.
        
        It is faster but not good because RectBivariateSpline does not
        provide fill_value option. Data outside of the available range
        are not assigned with a particular value.

    """
    interp_spline = RectBivariateSpline(lat, lon, arr)
    interp_data = interp_spline(new_lat, new_lon)

    if mask_idx is not None:
        interp_data[mask_idx] = np.nan
 
    return interp_data


def interpolate(lat, lon, data, new_lat, new_lon, mask_idx=None):

    ext_lon  = np.zeros((len(lon)+2,), dtype=lon.dtype)
    ext_data = np.zeros((len(lat), len(lon)+2,), dtype=data.dtype)
    
    ext_lon[0]    = lon[-1]   - 360.0
    ext_lon[1:-1] = lon
    ext_lon[-1]   = lon[0]    + 360

    ext_data[:, 0]    = data[:, -1]
    ext_data[:, 1:-1] = data
    ext_data[:, -1]   = data[:,  0]

    interp_func = scipy.interpolate.interp2d(ext_lon, lat, ext_data, kind='linear', fill_value=np.nan)
    interp_data = interp_func(new_lon, new_lat)

    if mask_idx is not None:
        interp_data[mask_idx] = np.nan
 
    return interp_data

def rearrangeLonAndField(lon, arrs, axis=1):
    
    if len(lon.shape) != 1:
        raise Exception("Longitude should be a 1-dim array.")

    new_lon = lon % 360.0
    dlon = new_lon[1:] - new_lon[:-1]
    jmps = np.sum(dlon < 0)

    if jmps == 0:
        
        new_arrs = arrs
 
    elif jmps == 1: 

        new_arrs = []
        
        jmp_idx = np.argwhere(dlon < 0)[0] # the skip happens between arr[skip_idx] and arr[skip_idx+1]
        shift = - ( jmp_idx + 1 )
        new_lon = np.roll(new_lon, shift)
        
        for i, arr in enumerate(arrs):
            new_arrs.append(np.roll(arr, shift, axis=axis))

    else:
        raise Exception("There is more than one jump in longitude")
            



    return new_lon, new_arrs


def getFileAndIndex(product, date, root_dir="data"):
    
    if product == "ERA5":
        
        filename = "ERA5_SST_%s.nc" % (date.strftime("%Y-%m-%d"))
        filename = os.path.join(root_dir, "ERA5", "SST", filename)
        idx = 0
        lat = "latitude"
        lon = "longitude"
        sst = 'sst'   
        unit = 'K'
 
    elif product == "OSTIA":
 
        filename = "OSTIA_SST_%s.nc" % (date.strftime("%Y-%m-%d"),)
        filename = os.path.join(root_dir, "OSTIA", "SST", filename)
        idx = 0
        lat = "lat"
        lon = "lon" 
        sst = 'analysed_sst'
        unit = 'K'


    elif product == "OISST":
        
        filename = "sst.day.mean.%s.nc" % (date.strftime("%Y"))
        filename = os.path.join(root_dir, "OISST", "SST", filename)

        idx = int((datetime(date.year, date.month, date.day) - datetime(date.year, 1, 1)).total_seconds() / 86400)

        lat = "lat"
        lon = "lon"
        sst = 'sst'
        unit = 'C'

    elif product == "CFSR":
        
        filename = "CFSR_SST_%s.nc" % (date.strftime("%Y-%m"))
        filename = os.path.join(root_dir, "CFSR", "SST", filename)

        idx = int((datetime(date.year, date.month, date.day) - datetime(date.year, date.month, 1)).total_seconds() / 86400)

        lat = "lat"
        lon = "lon" 
        sst = 'pt'
        unit = 'K'


    else:
        raise Exception("Unrecognized product: %s" % (product,))

    info = {
        'filename' : filename,
        'idx' : idx,
        'varnames' : {
            'sst': sst,
            'lat': lat,
            'lon': lon,
        },
        'unit' : unit,
    }
    return info
        
    

