import numpy as np
import netCDF4
from buoyancy_nonlinear import TS2rho
#from buoyancy_linear import TS2rho

default_fill_value = 1e20

def compute_jump(
    T,
    S,
    jump_depth,
    depth,
    sample_above_dist = 10.0,
    sample_below_dist = 20.0,
    fill_value = default_fill_value,
):


    Nt, Nz, Ny, Nx = T.shape

    db = np.zeros((Nt, Ny, Nx), dtype=np.float32)
    db[:] = fill_value 

    T_upper = db.copy()
    T_lower = db.copy()
    S_upper = db.copy()
    S_lower = db.copy()

    for t in range(Nt):
        for j in range(Ny): 
            for i in range(Nx): 
                
                h = jump_depth[t, j, i]

                if h is np.ma.masked:
                    #print("Masked. continue")
                    continue
    
                
                Nz_max = Nz
                for k in range(Nz):
                    if T[t, k, j, i] is np.ma.masked:
                        Nz_max = k
                        break
                
                _depth = depth[:Nz_max]
                depth_min = _depth[0]
                depth_max = _depth[-1]

                ##### Get Themperatrue
                depth_upper = max(depth_min, h - sample_above_dist)
                depth_lower = min(depth_max, h + sample_below_dist)
                sampling_depth = [depth_upper, depth_lower]
    
                #print(sampling_depth)
                #print(depth)

                T_sampled = np.interp(sampling_depth, _depth, T[t, 0:Nz_max, j, i], left=None, right=None, period=None)
                S_sampled = np.interp(sampling_depth, _depth, S[t, 0:Nz_max, j, i], left=None, right=None, period=None)


                T_upper[t, j, i] = T_sampled[0]
                T_lower[t, j, i] = T_sampled[1]
               
                S_upper[t, j, i] = S_sampled[0]
                S_lower[t, j, i] = S_sampled[1]

                db[t, j, i] = TS2b(T_sampled[0], S_sampled[0]) - TS2b(T_sampled[1], S_sampled[1])


    dT = T_upper - T_lower
    dS = S_upper - S_lower

    return T_upper, T_lower, S_upper, S_lower, dT, dS, db

def findMLD_rho(rho, z_T, dev=0.03, mask=None, Nz_bot=None, fill_value=default_fill_value):


    Nz, Ny, Nx = rho.shape

    if mask is None:
        mask = np.ones((Ny, Nx), dtype=np.int32)
 
    if Nz_bot is None:
        Nz_bot = np.zeros((Ny, Nx), dtype=np.int32)
        Nz_bot += Nz
        
        
    MLD = np.zeros((Ny, Nx))
    for j in range(Ny):
        for i in range(Nx):

#            print("(j, i) = (%d, %d)" % (j, i))

            if mask[j, i] == 0:
                MLD[j, i] = fill_value
                continue 

            SSrho = rho[0, j, i]
            rho_threshold = SSrho + dev

            #print("SSrho = ", SSrho)
            #print("dev = ", dev)
            #print("rho_threshold = ", rho_threshold)

            _Nz_bot = Nz_bot[j, i]
            for k in range(_Nz_bot):
                
                if rho[k, j, i] >= rho_threshold:
                
                    MLD[j, i] = - z_T[k-1] + (z_T[k-1] - z_T[k]) * (rho_threshold - rho[k-1, j, i])/(rho[k, j, i] - rho[k-1, j, i]) 
                    break

                if k == _Nz_bot-1:  # cannot find the mixed layer depth
                        MLD[j, i] = - z_T[k]

    # Sanity check
    if np.any(MLD[np.isfinite(MLD)] <= 0):
        throw(ErrorException("Some MLD is negative."))

    return MLD


def processECCO(
    input_filename,
    output_filename,
    varname_TEMP = "THETA",
    varname_SALT = "SALT",
    varname_z_T  = "Z",
    varname_z_W  = "Z_bnds",
    varname_lat  = "latitude",
    varname_lon  = "longitude",
    lon_l = 0.0,
    lon_r = 360.0,
    lat_b = -90.0,
    lat_t = 90.0,
    fill_value   = default_fill_value,
):


    print("INPUT FILE: ", input_filename)

    # 1. convert T and S to density
    # 2. Find the mixed-layer depth
    # 3. Compute mixed-layer temperature
    # 4. Compute mixed-layer salinity
    # 5. Output file.

    with netCDF4.Dataset(input_filename, "r") as ds:
        z_T = ds.variables[varname_z_T][:]
        z_W = ds.variables[varname_z_W][:, :]  # (z_T, nv=2)

        lat = ds.variables[varname_lat][:]
        lon = ds.variables[varname_lon][:] % 360.0

        lat_idx = (lat_b <= lat) & (lat <= lat_t)
        lon_idx = (lon_l <= lon) & (lon <= lon_r)

        lat = lat[lat_idx]
        lon = lon[lon_idx]



    with netCDF4.Dataset(input_filename, "r") as ds:
        TEMP = ds.variables[varname_TEMP][:, :, lat_idx, lon_idx]  # (time, z, y, x)
        SALT = ds.variables[varname_SALT][:, :, lat_idx, lon_idx]  # (time, z, y, x)


    #print(lat.shape)
    #print(lon.shape)
    #print(TEMP.shape)
    #print(SALT.shape)

    Nt, Nz, Ny, Nx = TEMP.shape
    print("Shape = (%d, %d, %d, %d)" % (Nt, Nz, Ny, Nx))
                
    mask = 1 - TEMP[0, 0, :, :].mask.astype(np.int32)
    mask3D = 1 - TEMP[0, :, :, :].mask.astype(np.int32)

    Nz_bot = np.sum(mask3D, axis=0)

    print("Compute density")    
    rho = TS2rho(TEMP, SALT)

    # compute density
    """
    for t in range(Nt):
        for k in range(Nz):
            print("(t, k) = (%d, %d)" % (t, k))
            for j in range(Ny):
                for i in range(Nx):
                    rho[t, k, j, i] = TS2rho(TEMP[t, k, j, i], SALT[t, k, j, i])
    """


    print("Compute mixed-layer depth")    
    # compute mixed-layer depth

    MLD = np.zeros((Nt, Ny, Nx))
    for t in range(Nt):
        MLD[t, :, :] = findMLD_rho(rho[t, :, :, :], z_T, mask=mask, Nz_bot=Nz_bot, dev=0.03)

    output_vars = {
        'MLD'     : MLD,
    }

    print("Writing output file: ", output_filename)
    with netCDF4.Dataset(output_filename, mode='w', format='NETCDF4_CLASSIC') as ds: 

        x_dim    = ds.createDimension('lon', Nx)
        y_dim    = ds.createDimension('lat', Ny)
        time_dim = ds.createDimension('time', None)


        var_lon = ds.createVariable('lon', np.float32, ('lon',))
        var_lat = ds.createVariable('lat', np.float32, ('lat',))
        
        var_lon[:] = lon
        var_lat[:] = lat

        for k, d in output_vars.items():
            _var = ds.createVariable(k, np.float32, ('time', 'lat', 'lon'), fill_value=fill_value)
            _var[:, :, :] = d



if __name__ == "__main__" : 

    print("*** This is for testing ***") 
    input_filename ="data/ECCO/ECCO_L4_TEMP_SALINITY_05DEG_DAILY_V4R4/OCEAN_TEMPERATURE_SALINITY_day_mean_1992-01-01_ECCO_V4r4_latlon_0p50deg.nc"

    output_filename = "test_convert_ECCO5.nc"

    print("Input  file : %s" % (input_filename,))
    print("Output file: %s" % (output_filename,)) 

    processECCO(
        input_filename,
        output_filename,
    ) 
    
