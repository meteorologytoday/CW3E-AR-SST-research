import numpy as np
import netCDF4
import vorticity_tools, advection_tools
from buoyancy_nonlinear import TS2rho, TS2b
import earth_constants as ec
import calculus_tools
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


def evalAtMLD_W(fi, h, z_W, mask=None, fill_value=default_fill_value):
  
    Nzp1, Ny, Nx = fi.shape

    if mask is None:
        mask = np.ones((Ny, Nx), dtype=np.int32)

    fo = np.zeros((Ny, Nx))
    Nz_h = detectMLNz(h, z_W, mask=mask)
    
    dz_T = z_W[:-1] - z_W[1:]

    for j in range(Ny):
        for i in range(Nx):
            
            if mask[j, i] == 0:
                fo[j, i] = fill_value
            
            else:
                
                _Nz = Nz_h[j, i]
                _h = h[j, i]

                fo[j, i] = fi[_Nz+1, j, i] + (fi[_Nz, j, i] - fi[_Nz+1, j, i]) / dz_T[_Nz] * (- _h - z_W[_Nz+1])
   
    return fo


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
    input_filename_TS,
    input_filename_VEL,
    input_filename_SSH,
    output_filename,
    varname_TEMP = "THETA",
    varname_SALT = "SALT",
    varname_UVEL = "EVEL",
    varname_VVEL = "NVEL",
    varname_WVEL = "WVEL",
    varname_SSH  = "SSH",
    varname_z_T  = "Z",
    varname_z_W  = "Z_bnds",
    varname_lat  = "latitude",
    varname_lon  = "longitude",
    input_filename_MLD = "",
    varname_MLD  = "MXLDEPTH",
    lon_l = 0.0,
    lon_r = 360.0,
    lat_b = -90.0,
    lat_t = 90.0,
    MLD_dev = 0.03,
    fill_value   = default_fill_value,
):


    print("INPUT FILE TS : ", input_filename_TS)
    print("INPUT FILE VEL: ", input_filename_VEL)
    print("INPUT FILE SSH: ", input_filename_SSH)
    print("INPUT FILE MLD: ", input_filename_MLD)

    # 1. convert T and S to density
    # 2. Find the mixed-layer depth
    # 3. Compute mixed-layer temperature
    # 4. Compute mixed-layer salinity
    # 5. Output file.

    with netCDF4.Dataset(input_filename_TS, "r") as ds:
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



    with netCDF4.Dataset(input_filename_TS, "r") as ds:
        TEMP = ds.variables[varname_TEMP][:, :, lat_idx, lon_idx]  # (time, z, y, x)
        SALT = ds.variables[varname_SALT][:, :, lat_idx, lon_idx]  # (time, z, y, x)

    with netCDF4.Dataset(input_filename_VEL, "r") as ds:
        UVEL = ds.variables[varname_UVEL][:, :, lat_idx, lon_idx]  # (time, z, y, x)
        VVEL = ds.variables[varname_VVEL][:, :, lat_idx, lon_idx]  # (time, z, y, x)
        WVEL = ds.variables[varname_WVEL][:, :, lat_idx, lon_idx]  # (time, z, y, x)

    with netCDF4.Dataset(input_filename_SSH, "r") as ds:
        SSH  = ds.variables[varname_SSH][:, lat_idx, lon_idx]  # (time, y, x)

    if input_filename_MLD != "":
        with netCDF4.Dataset(input_filename_MLD, "r") as ds:
            MLD  = ds.variables[varname_MLD][:, lat_idx, lon_idx]  # (time, y, x)



    f_co = 2 * ec.Omega * np.sin(lat * np.pi / 180)
    f_co = f_co[:, None]

    Nt, Nz, Ny, Nx = TEMP.shape
    print("Shape = (%d, %d, %d, %d)" % (Nt, Nz, Ny, Nx))
                
    mask = 1 - TEMP[0, 0, :, :].mask.astype(np.int32)
    mask3D = 1 - TEMP[0, :, :, :].mask.astype(np.int32)

    Nz_bot = np.sum(mask3D, axis=0)

    print("Compute density and buoyancy") 
    rho = TS2rho(TEMP, SALT)
    b   = TS2b(TEMP, SALT)

    print("Compute mixed-layer depth")    
    # compute mixed-layer depth

    if input_filename_MLD == "":
        MLD = np.zeros((Nt, Ny, Nx))

    MLT = np.zeros((Nt, Ny, Nx))
    MLS = np.zeros((Nt, Ny, Nx))
    
    MLU = np.zeros((Nt, Ny, Nx))
    MLV = np.zeros((Nt, Ny, Nx))


    dMLTdx = np.zeros((Nt, Ny, Nx))
    dMLTdy = np.zeros((Nt, Ny, Nx))
    dMLSdx = np.zeros((Nt, Ny, Nx))
    dMLSdy = np.zeros((Nt, Ny, Nx))
    dMLDdx = np.zeros((Nt, Ny, Nx))
    dMLDdy = np.zeros((Nt, Ny, Nx))
 
    dT = np.zeros((Nt, Ny, Nx))
    dS = np.zeros((Nt, Ny, Nx))
    db = np.zeros((Nt, Ny, Nx))


    w_b = np.zeros((Nt, Ny, Nx))
    
    V_g = np.zeros((Nt, Ny, Nx))
    U_g = np.zeros((Nt, Ny, Nx))
    
    dTdz_b = np.zeros((Nt, Ny, Nx))
    dSdz_b = np.zeros((Nt, Ny, Nx))
    dUdz_b = np.zeros((Nt, Ny, Nx))
    dVdz_b = np.zeros((Nt, Ny, Nx))
    N2_b   = np.zeros((Nt, Ny, Nx))
    
    lapT   = np.zeros((Nt, Ny, Nx))
    lapS   = np.zeros((Nt, Ny, Nx))


    print("Compute needed variables...")    
    for t in range(Nt):

        if input_filename_MLD == "":
            MLD[t, :, :] = findMLD_rho(rho[t, :, :, :], z_T, mask=mask, Nz_bot=Nz_bot, dev=MLD_dev)
    
        MLT[t, :, :] = computeMLMean(TEMP[t, :, :, :], MLD[t, :, :], z_W, mask=mask)
        MLS[t, :, :] = computeMLMean(SALT[t, :, :, :], MLD[t, :, :], z_W, mask=mask)
        MLU[t, :, :] = computeMLMean(UVEL[t, :, :, :], MLD[t, :, :], z_W, mask=mask)
        MLV[t, :, :] = computeMLMean(VVEL[t, :, :, :], MLD[t, :, :], z_W, mask=mask)

        T_b = evalAtMLD_T(TEMP[t, :, :, :], MLD[t, :, :], z_W, mask=mask)
        S_b = evalAtMLD_T(SALT[t, :, :, :], MLD[t, :, :], z_W, mask=mask)
        w_b[t, :, :] = evalAtMLD_T(WVEL[t, :, :, :], MLD[t, :, :], z_W, mask=mask)

        MLb = TS2b(MLT[t, :, :], MLS[t, :, :])
        b_b = TS2b(T_b, S_b)
        
        dT[t, :, :] = MLT[t, :, :] - T_b
        dS[t, :, :] = MLS[t, :, :] - S_b
        db[t, :, :] = MLb - b_b
       
        # Compute geostrophic velocity           
        dMLTdx[t, :, :], dMLTdy[t, :, :] = advection_tools.calGrad(MLT[t, :, :], lat, lon, periodoc_lon=False, mask=mask)
        dMLSdx[t, :, :], dMLSdy[t, :, :] = advection_tools.calGrad(MLS[t, :, :], lat, lon, periodoc_lon=False, mask=mask)
        dMLDdx[t, :, :], dMLDdy[t, :, :] = advection_tools.calGrad(MLD[t, :, :], lat, lon, periodoc_lon=False, mask=mask)
        dSSHdx, dSSHdy = advection_tools.calGrad(SSH[t, :, :], lat, lon, periodoc_lon=False, mask=mask)

        V_g[t, :, :] =   ec.g0 / f_co * dSSHdx
        U_g[t, :, :] = - ec.g0 / f_co * dSSHdy


        dTdz = calculus_tools.W_ddz_T(TEMP[t, :, :, :], z_T=z_T)
        dSdz = calculus_tools.W_ddz_T(SALT[t, :, :, :], z_T=z_T)
        dUdz = calculus_tools.W_ddz_T(UVEL[t, :, :, :], z_T=z_T)
        dVdz = calculus_tools.W_ddz_T(VVEL[t, :, :, :], z_T=z_T)
        N2   = calculus_tools.W_ddz_T(b[t, :, :, :], z_T=z_T)
        
        dTdz_b[t, :, :] = evalAtMLD_W(dTdz, MLD[t, :, :], z_W, mask=mask)
        dSdz_b[t, :, :] = evalAtMLD_W(dSdz, MLD[t, :, :], z_W, mask=mask)
        dUdz_b[t, :, :] = evalAtMLD_W(dUdz, MLD[t, :, :], z_W, mask=mask)
        dVdz_b[t, :, :] = evalAtMLD_W(dVdz, MLD[t, :, :], z_W, mask=mask)
        N2_b[t, :, :]   = evalAtMLD_W(N2,   MLD[t, :, :], z_W, mask=mask)
        
        _lapT   = np.zeros((Nz, Ny, Nx))
        _lapS   = np.zeros((Nz, Ny, Nx))

        for k in range(Nz):
            _lapT[k, :, :] = advection_tools.calLaplacian(TEMP[t, k, :, :], lat, lon, periodoc_lon=False, mask=mask)
            _lapS[k, :, :] = advection_tools.calLaplacian(SALT[t, k, :, :], lat, lon, periodoc_lon=False, mask=mask)

        lapT[t, :, :] = computeMLMean(_lapT, MLD[t, :, :], z_W, mask=mask)
        lapS[t, :, :] = computeMLMean(_lapS, MLD[t, :, :], z_W, mask=mask)

    output_vars = {
        'MLD'     : MLD,
        'MLT'     : MLT,
        'MLS'     : MLS,
        'SST'     : TEMP[:, 0, :, :],
        'SSS'     : SALT[:, 0, :, :],
        'dT'      : dT,
        'dS'      : dS,
        'db'      : db,
        'MLU'     : MLU,
        'MLV'     : MLV,
        'w_b'     : w_b,
        'V_g'     : V_g,
        'U_g'     : U_g,
        'dMLTdx'  : dMLTdx,
        'dMLTdy'  : dMLTdy,
        'dMLSdx'  : dMLSdx,
        'dMLSdy'  : dMLSdy,
        'dMLDdx'  : dMLDdx,
        'dMLDdy'  : dMLDdy,
        'dTdz_b'  : dTdz_b,
        'dSdz_b'  : dSdz_b,
        'dUdz_b'  : dUdz_b,
        'dVdz_b'  : dVdz_b,
        'N2_b'    : N2_b,
        'lapT'    : lapT,
        'lapS'    : lapS,
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

            print("Output %s: " % (k, ), "; dim = ", d.shape)

            _var = ds.createVariable(k, np.float32, ('time', 'lat', 'lon'), fill_value=fill_value)
            _var[:, :, :] = d



if __name__ == "__main__" : 

    print("*** This is for testing ***") 
    input_filename_TS ="data/ECCO/ECCO_L4_TEMP_SALINITY_05DEG_DAILY_V4R4/OCEAN_TEMPERATURE_SALINITY_day_mean_1997-12-10_ECCO_V4r4_latlon_0p50deg.nc"
    input_filename_SSH ="data/ECCO/ECCO_L4_SSH_05DEG_DAILY_V4R4B/SEA_SURFACE_HEIGHT_day_mean_1997-12-10_ECCO_V4r4b_latlon_0p50deg.nc"
    input_filename_VEL ="data/ECCO/ECCO_L4_OCEAN_VEL_05DEG_DAILY_V4R4/OCEAN_VELOCITY_day_mean_1997-12-10_ECCO_V4r4_latlon_0p50deg.nc"

    output_filename = "test_convert_ECCO5.nc"

    print("Input  file TS  : %s" % (input_filename_TS,))
    print("Input  file VEL : %s" % (input_filename_VEL,))
    print("Input  file SSH : %s" % (input_filename_SSH,))
    print("Output file: %s" % (output_filename,)) 

    processECCO(
        input_filename_TS,
        input_filename_VEL,
        input_filename_SSH,
        output_filename,
    ) 
    
