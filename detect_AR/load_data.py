from datetime import datetime
import os.path
import numpy as np
import scipy.interpolate
from scipy.interpolate import RectBivariateSpline


def getFileAndIndex(product, date, root_dir="data", varname="", **kwargs):
    
    if product == "ERA5":

        if varname in ["sst", "mslhf", "msshf", "msnlwrf", "msnswrf", "mtpr", "mer", "mvimd", "t2m", "u10", "v10",] :
            subfolder = "sfc"
            filename = "ERA5_sfc_%s.nc" % (date.strftime("%Y-%m-%d"),)
        elif varname in ["IWV", "IVT", "IWVKE"]:
            subfolder = "AR_processed"
            filename = "ERA5_AR_%s.nc" % (date.strftime("%Y-%m-%d"),)
        elif varname in ["vort10", "curltau"]:
            subfolder = "sfc_processed"
            filename = "ERA5_sfc_processed_%s.nc" % (date.strftime("%Y-%m-%d"),)
        else:
            raise Exception("Unrecognized varname: %s " % (varname,) )


        filename = os.path.join(root_dir, "ERA5", subfolder, filename)

        idx = 0
        lat = "lat"
        lon = "lon"
 
    elif product == "GFS":

        fcst = kwargs['fcst']

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

    elif product == "ECCO":

        if varname in ["SST", "MLT", "SSS", "MLS", "MLD",] :
            filename = "ECCO_mixedlayer_0p50deg_%s.nc" % (date.strftime("%Y-%m-%d"),)
        else:
            raise Exception("Unrecognized varname: %s " % (varname,) )

        filename = os.path.join(root_dir, "ECCO", "processed_0p25deg", filename)

        idx = 0
        lat = "lat"
        lon = "lon"
 
    elif product == "ORA5":

        timestr = date.strftime("%Y-%m")

        if varname in ["MLD", "T_upper", "T_lower", "S_upper", "S_lower", "db", "dT"]:

            mxl_algo = kwargs['mxl_algo']

            if mxl_algo in ['somxl010', 'somxl030']:

                filename = "ORA5_NillerKrausMixedLayerDynamics_%s_%s.nc" % (mxl_algo, timestr)
                subfolder = "processed_remapped" 

            else:
                raise Exception("Unknown mixed-layer algorithm %s" % (mxl_algo,))

        else:
            raise Exception("Unrecognized varname: %s " % (varname,) )


        filename = os.path.join(root_dir, "ORA5", subfolder, filename)
        idx = 0
        lat = "lat"
        lon = "lon"

    elif product == "ORA5-clim":

        timestr = date.strftime("%m")

        if varname in ["MLD", "T_upper", "T_lower", "S_upper", "S_lower", "db", "dT"]:

            mxl_algo = kwargs['mxl_algo']

            if mxl_algo in ['somxl010', 'somxl030']:

                filename = "ORA5_NillerKrausMixedLayerDynamics_%s_%s.nc" % (mxl_algo, timestr)
                kwargs['beg_year'] = 2001
                kwargs['end_year'] = 2014
                subfolder = "processed_remapped_clim_%04d-%04d" % (kwargs['beg_year'], kwargs['end_year']) 

            else:
                raise Exception("Unknown mixed-layer algorithm %s" % (mxl_algo,))

        else:
            raise Exception("Unrecognized varname: %s " % (varname,) )


        filename = os.path.join(root_dir, "ORA5", subfolder, filename)
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
        
    

