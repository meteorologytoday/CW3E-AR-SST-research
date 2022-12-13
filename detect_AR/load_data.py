from datetime import datetime
import os.path
import numpy as np
import scipy.interpolate
from scipy.interpolate import RectBivariateSpline


def getFileAndIndex(product, date, root_dir="data", fcst=-1, varname=""):
    
    if product == "ERA5":

        if varname in ["sst", "mslhf", "msshf", "msnlwrf", "msnswrf", "mtpr", "mvimd", "t2m", "u10", "v10",] :
            subfolder = "sfc"
            filename = "ERA5_sfc_%s.nc" % (date.strftime("%Y-%m-%d"))
        elif varname in ["IWV", "IVT", "IWVKE"]:
            subfolder = "AR_processed"
            filename = "ERA5_AR_%s.nc" % (date.strftime("%Y-%m-%d"))
        else:
            raise Exception("Unrecognized varname: %s " % (varname,) )


        filename = os.path.join(root_dir, "ERA5", subfolder, filename)

        idx = 0
        lat = "lat"
        lon = "lon"
 
    elif product == "GFS":

        if fcst == -1:
            raise Exception("Need parameter `fcst` -- forecast hour.")

        timestr = date.strftime("%Y%m%d")
        if varname in ["IWV", "IVT", "IWVKE"]:
            filename = "GFS_0p25_%s_f%03d.AR.nc" % (timestr, fcst)
            subfolder = "fcst"
        elif varname == "HGT_500mb":
            filename = "GFS_0p25_%s_f%03d.HGT_500mb.nc" % (timestr, fcst)
            subfolder = "fcst"
        elif varname == "HGT_850mb":
            filename = "GFS_0p25_%s_f%03d.HGT_850mb.nc" % (timestr, fcst)
            subfolder = "fcst"
        elif varname in ["LHTFL", "SHTFL", "PRATE", "UFLX", "VFLX"]:
            filename = "GFS_0p25_%s_f%03d.sfcflx.nc" % (timestr, fcst)
            subfolder = "fcst_sfc"
        else:
            raise Exception("Unrecognized varname: %s " % (varname,) )


        filename = os.path.join(root_dir, "GFS", subfolder, filename)
        idx = 0
        lat = "lat"
        lon = "lon"

    else:
        raise Exception("Unrecognized product: %s" % (product,))

    info = {
        'filename' : filename,
        'idx' : idx,
        'varnames' : {
            'lat': lat,
            'lon': lon,
        },
    }
    return info
        
    

