import xarray as xr


def loadDatasets(input_dir, yrs, file_fmt="AR_statistics_yr%04d.nc"):


    filenames = [ "%s/%s" % (input_dir, file_fmt % yr) for yr in yrs  ]


    data = xr.load_mfdataset(filenames)

    return data
    
