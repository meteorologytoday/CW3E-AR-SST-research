
from pathlib import Path
import os.path
import xarray as xr
import ECCO_helper
import ecco_v4_py as ecco

import numpy as np
import netCDF4
import vorticity_tools, advection_tools
from buoyancy_nonlinear import TS2rho, TS2b
import earth_constants as ec
import calculus_tools

from datetime import (datetime, timedelta)

#from buoyancy_linear import TS2rho

default_fill_value = np.nan
default_fill_value_int = -1

# RHO_CONST number is read from
# https://ecco-v4-python-tutorial.readthedocs.io/Thermal_wind.html#Viewing-and-Plotting-Density
RHO_CONST = 1029.0

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

                    if z_W[k] >= z and z >= z_W[k+1]:   # Using 2 equalities so that the last grid-box will include z = z_bottom
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

            if mask[j, i] == 0:
                MLD[j, i] = fill_value
                continue 

            SSrho = rho[0, j, i]
            rho_threshold = SSrho + dev

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
    target_datetime,
    output_filename,
    MLD_dev = 0.03,
):

    print("Target datetime: ", target_datetime)
    print("Output filename: ", output_filename)

    beg_datetime = target_datetime # - timedelta(days=1)

    snp_varnames = ["THETA", "SALT", "ETAN"]
    ave_varnames = ["MXLDEPTH", "G_ttl", "G_hadv", "G_vadv", "G_hdiff", "G_vdiff"]
    #ave_varnames = ["MXLDEPTH", ]

    ds = ECCO_helper.loadECCOData_continuous(
        beg_datetime = beg_datetime,
        ndays = 1,
        snp_varnames = snp_varnames,
        ave_varnames = ave_varnames,
    )

    xgcm_grid = ecco.get_llc_grid(ds)
    ecco_grid = ECCO_helper.getECCOGrid()
    
    sTHETA_snp = ds.THETA_snp * (1 + ds.ETAN_snp/ecco_grid.Depth)
    sSALT_snp  = ds.SALT_snp  * (1 + ds.ETAN_snp/ecco_grid.Depth)

    sTHETA_snp = sTHETA_snp.rename("sTHETA_snp")
    
    Nt, Nz, Nl, Nx, Ny = (ds.dims['time'], ds.dims['k'], ds.dims['tile'], ds.dims['j'], ds.dims['i'])
  
    if Nt != 1:
        raise Exception("Too many records. I only need one.")
 
    Zl = ecco_grid.Zl.load()
    Zu = ecco_grid.Zu.load()

    z_W = np.zeros((len(Zl)+1,),)


    z_W[:-1] = Zl
    z_W[-1] = Zu[-1]
    
    z_T = (z_W[:-1] + z_W[1:]) / 2.0

 
    sample2D_snp = ds.THETA_snp[:, 0, :, :, :]
    sample2D_ave = ds.MXLDEPTH[:, :, :, :]

    ML_snp = {}
    ML_ave = {}

    ML_snp_varnames = ["MLD_snp", "MLT_snp", "MLS_snp"]
    ML_ave_varnames = ["dMLTdt", "MLG_ttl", "MLG_hadv", "MLG_vadv", "MLG_adv", "MLG_vdiff", "MLG_hdiff", "MLG_forcing"]

    for varname in ML_snp_varnames:
        ML_snp[varname] = xr.zeros_like(sample2D_snp).rename(varname)

    for varname in ML_ave_varnames:
        ML_ave[varname] = xr.zeros_like(sample2D_ave).rename(varname)

 
    rho_snp = TS2rho(ds.THETA_snp, ds.SALT_snp).rename('rho_snp')
    
    mask   = [ ds.THETA_snp[0, :, l, :, :].notnull().rename('mask').astype('i4').to_numpy() for l in range(Nl) ]
    mask2D = [ mask[l][0, :, :] for l in range(Nl) ]

    Nz_bot = [ np.sum(mask[l], axis=0) for l in range(Nl) ]

    # Compute variable at the time_bnds
    for s in range(2):
        for l in range(Nl):

            ML_snp["MLD_snp"][s, l, :, :] = findMLD_rho(
                rho_snp[s, :, l, :, :].to_numpy(),
                z_T,
                mask=mask2D[l],
                Nz_bot=Nz_bot[l],
                dev=MLD_dev,
            )

            ML_snp["MLT_snp"][s, l, :, :] = computeMLMean(
                sTHETA_snp[s, :, l, :, :].to_numpy(),
                ML_snp["MLD_snp"][s, l, :, :].to_numpy(),
                z_W,
                mask=mask2D[l]
            )

            ML_snp["MLS_snp"][s, l, :, :] = computeMLMean(
                sSALT_snp[s, :, l, :, :].to_numpy(),
                ML_snp["MLD_snp"][s, l, :, :].to_numpy(),
                z_W,
                mask=mask2D[l]
            )

    dt = 86400.0
    #xgcm_grid.diff(ds.time_snp, 'T', boundary='fill', fill_value=np.nan).astype('f4') / 1e9 # nanosec to sec 
   
    ML_ave["dMLTdt"][0, :, :, :] = (ML_snp["MLT_snp"][1, :, :, :] - ML_snp["MLT_snp"][0, :, :, :]) / dt


    # compute variable in the middle of time_bnds
    MLD_0 = ML_snp["MLD_snp"][0, :, :, :].to_numpy()
    MLD_1 = ML_snp["MLD_snp"][1, :, :, :].to_numpy()
    for varname in ["G_ttl", "G_hadv", "G_vadv", "G_vdiff", "G_hdiff", "G_forcing"]:
        ML_varname = "ML%s" % varname
        for l in range(Nl):
                
            ML_ave[ML_varname][0, l, :, :] = computeMLMean(
                ds[varname][0, :, l, :, :].to_numpy(),
                MLD_1[l, :, :],
                z_W,
                mask=mask2D[l]
            )

    # Compute entrainment term explicitly
    ML_ave["MLG_ent"] = xr.zeros_like(sample2D_ave).rename("MLG_ent") 
    for l in range(Nl):
            
        sTHETA = sTHETA_snp[0, :, l, :, :].to_numpy()
        ML_ave["MLG_ent"][0, l, :, :] = ( computeMLMean(
            sTHETA,
            MLD_1[l, :, :],
            z_W,
            mask=mask2D[l]
        ) -  computeMLMean(
            sTHETA,
            MLD_0[l, :, :],
            z_W,
            mask=mask2D[l]
        )) / dt


    # Additional diagnostic variables 
    ML_ave["MLG_adv"][:, :, :, :] = ML_ave["MLG_hadv"] + ML_ave["MLG_vadv"]
    ML_ave["dMLTdt_res"] = ML_ave["dMLTdt"] - ( ML_ave["MLG_ent"] + ML_ave["MLG_hadv"] + ML_ave["MLG_vadv"] + ML_ave["MLG_hdiff"] + ML_ave["MLG_vdiff"] + ML_ave["MLG_forcing"] )

    ML_ave["dMLTdt_res"] = ML_ave["dMLTdt_res"].rename("dMLTdt_res")

    MLU = np.zeros((Nt, Ny, Nx))
    MLV = np.zeros((Nt, Ny, Nx))


    dMLTdx = np.zeros((Nt, Ny, Nx))
    dMLTdy = np.zeros((Nt, Ny, Nx))
    dMLSdx = np.zeros((Nt, Ny, Nx))
    dMLSdy = np.zeros((Nt, Ny, Nx))
    dMLDdx = np.zeros((Nt, Ny, Nx))
    dMLDdy = np.zeros((Nt, Ny, Nx))
 

    output_data = []

    for k, var in ML_ave.items():
        output_data.append(var)

    for var in output_data:
        for attr in ["valid_min", "valid_max"]:
            if attr in var.attrs:
                del(var.attrs[attr])

    ds_out = xr.merge(output_data, compat='override')


    print("Output: ", output_filename)

    dir_name = os.path.dirname(output_filename)
    if not os.path.isdir(dir_name):
        print("Create dir: %s" % (dir_name,))
        Path(dir_name).mkdir(parents=True, exist_ok=True)


    ds_out.to_netcdf(output_filename)
    # 1. convert T and S to density
    # 2. Find the mixed-layer depth
    # 3. Compute mixed-layer temperature
    # 4. Compute mixed-layer salinity
    # 5. Output file.


    #f_co = 2 * ec.Omega * np.sin(lat * np.pi / 180)
    #f_co = f_co[:, None]

    #Nt, Nz, Ny, Nx = TEMP.shape
    #print("Shape = (%d, %d, %d, %d)" % (Nt, Nz, Ny, Nx))
                
    #mask = 1 - TEMP[0, 0, :, :].mask.astype(np.int32)
    #mask3D = 1 - TEMP[0, :, :, :].mask.astype(np.int32)

    #Nz_bot = np.sum(mask3D, axis=0)

    """
    print("Compute density and buoyancy") 
    rho = TS2rho(TEMP, SALT)
    b   = TS2b(TEMP, SALT)
    """
    """
    N2 = - ec.g0 / RHO_CONST * dRHOdz
     

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
            MLD[t, :, :] = findMLD_rho(RHO[t, :, :, :], z_T, mask=mask, Nz_bot=Nz_bot, dev=MLD_dev)
    
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
        
        dTdz_b[t, :, :] = evalAtMLD_W(dTdz, MLD[t, :, :], z_W, mask=mask)
        dSdz_b[t, :, :] = evalAtMLD_W(dSdz, MLD[t, :, :], z_W, mask=mask)
        dUdz_b[t, :, :] = evalAtMLD_W(dUdz, MLD[t, :, :], z_W, mask=mask)
        dVdz_b[t, :, :] = evalAtMLD_W(dVdz, MLD[t, :, :], z_W, mask=mask)
        N2_b[t, :, :]   = evalAtMLD_T(N2[t, :, :, :], MLD[t, :, :], z_W, mask=mask)
        
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
    """


if __name__ == "__main__" : 

    print("*** This is for testing ***") 

    target_datetime = datetime(2017, 12, 26)
    output_filename = "data/ECCO_LLC/%s/%s" % ECCO_helper.getECCOFilename("MLT", "DAILY", target_datetime)

    processECCO(
        target_datetime,
        output_filename,
    ) 
    
