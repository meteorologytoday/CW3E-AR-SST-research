import numpy as np
import netCDF4
from buoyancy_nonlinear import TS2rho, TS2b
#from buoyancy_linear import TS2rho

default_fill_value = 1e20
default_fill_value_int = -1

def detectMLNz(h, z_W, mask=None, fill_value=default_fill_value_int):

    if mask is None:
        mask = np.ones((Ny, Nx), dtype=np.int32)

    Ny, Nx = h.shape
    Nz = len(z_W) - 1
    MLNz = np.zeros((Ny, Nx), dtype=np.int32)
   
    for j in range(Ny):
        for i in range(Nx):
        
            if mask[j, i] == 0:
                MLNz[j, i] = fill_value

            else:

                z = - h[j, i]
                for k in range(Nz):

                    if z_W[k] >= z >= z_W[k+1]:   # Using 2 equalities so that the last grid-box will include z = z_bottom
                        MLNz[j, i] = k
                        break
 
                    elif k == Nz-1:
                        MLNz[j, i] = k
    
    return MLNz



def evalAtMLD_T(fi, h, z_W, mask=None, fill_value=default_fill_value):
  
    Nz, Ny, Nx = fi.shape

    if mask is None:
        mask = np.ones((Ny, Nx), dtype=np.int32)

    fo = np.zeros((Ny, Nx))
    Nz_h = detectMLNz(h, z_W, mask=mask)
    
    for j in range(Ny):
        for i in range(Nx):
            
            if mask[j, i] == 0:
                fo[j, i] = fill_value
            
            else:
                
                _Nz = Nz_h[j, i]
                fo[j, i] = fi[_Nz, j, i]
   
    return fo



def computeMLMean(fi, h, z_W, mask=None, fill_value=default_fill_value):
  
    Nz, Ny, Nx = fi.shape 

    if mask is None:
        mask = np.ones((Ny, Nx), dtype=np.int32)

    dz_T = z_W[:-1] - z_W[1:]
    fo = np.zeros((Ny, Nx))
    Nz_h = detectMLNz(h, z_W, mask=mask)
    
    for j in range(Ny):
        for i in range(Nx):
           

            _Nz = Nz_h[j, i]
            
#            print("(j, i) = (%d, %d); _Nz = %d" % (j, i, _Nz)) 

            if mask[j, i] == 0:
                fo[j, i] = fill_value

            else:
                _tmp = 0.0
                if _Nz > 0:
                    _tmp += np.sum(dz_T[:_Nz] * fi[:_Nz, j, i])

                _tmp += (z_W[_Nz] + h[j, i]) * fi[_Nz, j, i]
                
                fo[j, i] = _tmp / h[j, i]
    
    return fo


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
        z_W_tmp = ds.variables[varname_z_W][:, :]  # (z_T, nv=2)

        z_W = np.zeros((len(z_T)+1,))
        z_W[:-1] = z_W_tmp[:, 0]
        z_W[-1] = z_W_tmp[-1, 1]


        lat = ds.variables[varname_lat][:]
        lon = ds.variables[varname_lon][:] % 360.0

        lat_idx = (lat_b <= lat) & (lat <= lat_t)
        lon_idx = (lon_l <= lon) & (lon <= lon_r)

        lat = lat[lat_idx]
        lon = lon[lon_idx]



    with netCDF4.Dataset(input_filename, "r") as ds:
        TEMP = ds.variables[varname_TEMP][:, :, lat_idx, lon_idx]  # (time, z, y, x)
        SALT = ds.variables[varname_SALT][:, :, lat_idx, lon_idx]  # (time, z, y, x)


    Nt, Nz, Ny, Nx = TEMP.shape
    print("Shape = (%d, %d, %d, %d)" % (Nt, Nz, Ny, Nx))
                
    mask = 1 - TEMP[0, 0, :, :].mask.astype(np.int32)
    mask3D = 1 - TEMP[0, :, :, :].mask.astype(np.int32)

    Nz_bot = np.sum(mask3D, axis=0)

    print("Compute density")    
    rho = TS2rho(TEMP, SALT)

    print("Compute mixed-layer depth")    
    # compute mixed-layer depth

    MLD = np.zeros((Nt, Ny, Nx))
    MLT = np.zeros((Nt, Ny, Nx))
    MLS = np.zeros((Nt, Ny, Nx))
    
    dT = np.zeros((Nt, Ny, Nx))
    dS = np.zeros((Nt, Ny, Nx))
    db = np.zeros((Nt, Ny, Nx))

    for t in range(Nt):
        MLD[t, :, :] = findMLD_rho(rho[t, :, :, :], z_T, mask=mask, Nz_bot=Nz_bot, dev=0.03)
        MLT[t, :, :] = computeMLMean(TEMP[t, :, :, :], MLD[t, :, :], z_W, mask=mask)
        MLS[t, :, :] = computeMLMean(SALT[t, :, :, :], MLD[t, :, :], z_W, mask=mask)

        T_b = evalAtMLD_T(TEMP[t, :, :, :], MLD[t, :, :], z_W, mask=mask)
        S_b = evalAtMLD_T(SALT[t, :, :, :], MLD[t, :, :], z_W, mask=mask)

        MLb = TS2b(MLT[t, :, :], MLS[t, :, :])
        b_b = TS2b(T_b, S_b)
        
        dT[t, :, :] = MLT[t, :, :] - T_b
        dS[t, :, :] = MLS[t, :, :] - S_b
        db[t, :, :] = MLb - b_b
        
               
 



    output_vars = {
        'MLD'     : MLD,
        'MLT'     : MLT,
        'MLS'     : MLS,
        'SST'     : TEMP[:, 0, :, :],
        'SSS'     : SALT[:, 0, :, :],
        'dT'      : dT,
        'dS'      : dS,
        'db'      : db,
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
    
