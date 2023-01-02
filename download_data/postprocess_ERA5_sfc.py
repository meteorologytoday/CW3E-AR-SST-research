import numpy as np
import netCDF4
from buoyancy_linear import TS2b
import vorticity_tools, Ekman_tools
import NK_tools
import earth_constants as ec

fill_value=1e20


def processERA5Sfc(
    input_filename,
    output_filename,
    lat_varname = "latitude",
    lon_varname = "longitude",
):


    #print("INPUT FILE: ", input_filename)
    #print("INPUT FILE: ", input_filename)

    with netCDF4.Dataset(input_filename, "r") as ds:
        u10 = ds.variables["u10"][:, :, :] # (time, y, x)
        v10 = ds.variables["v10"][:, :, :] # (time, y, x)
        sst = ds.variables["sst"][:, :, :] # (time, y, x)
        lat = ds.variables[lat_varname][:]
        lon = ds.variables[lon_varname][:]

    Nlon = len(lon)
    Nlat = len(lat)

    u10 = np.reshape(u10, (Nlat, Nlon))
    v10 = np.reshape(v10, (Nlat, Nlon))
    sst = np.reshape(sst, (Nlat, Nlon))

    _u10_abs = np.abs(u10)
    _v10_abs = np.abs(v10)
    taux = ec.rho_a * NK_tools.calu_star(_u10_abs)**2 * np.sign(u10)
    tauy = ec.rho_a * NK_tools.calu_star(_v10_abs)**2 * np.sign(v10)


    vort10 = vorticity_tools.uv2vort(u10, v10, lat, lon, periodoc_lon=True)
    curltau = vorticity_tools.uv2vort(taux, tauy, lat, lon, periodoc_lon=True)
    EkmanAdv = Ekman_tools.calEkmanAdv(sst, taux, tauy, 50.0, lat, lon, periodoc_lon=True) 

    output_vars = {
        'vort10' : vort10,
        'curltau': curltau,
        'EkmanAdv': EkmanAdv,
    }

    with netCDF4.Dataset(output_filename, mode='w', format='NETCDF4_CLASSIC') as ds: 

        lat_dim    = ds.createDimension('lat', Nlat)
        lon_dim    = ds.createDimension('lon', Nlon)
        time_dim = ds.createDimension('time', None)

        var_lat = ds.createVariable('lat', np.float32, ('lat',))
        var_lon = ds.createVariable('lon', np.float32, ('lon',))
        
        var_lat[:] = lat
        var_lon[:] = lon

        for k, d in output_vars.items():
            _var = ds.createVariable(k, np.float32, ('time', 'lat', 'lon'), fill_value=fill_value)
            _var[0, :, :] = d




   

if __name__ == "__main__" : 

    print("*** This is for testing ***") 
    input_filename = "data/ERA5/sfc/ERA5_sfc_2018-01-07.nc"
    output_filename = "test_convert_ERA5_sfc_to_vort.nc"

    print("Input  file: %s" % (input_filename,)) 
    print("Output file: %s" % (output_filename,)) 
    processERA5Sfc(input_filename, output_filename)
 
