
from pathlib import Path
import os.path
import xarray as xr
import ECCO_helper
import ecco_v4_py as ecco

import numpy as np
import netCDF4
import earth_constants as ec
import calculus_tools

from datetime import (datetime, timedelta)

from buoyancy_linear import TS2rho

default_fill_value = np.nan
default_fill_value_int = -1

# RHO_CONST number is read from
# https://ecco-v4-python-tutorial.readthedocs.io/Thermal_wind.html#Viewing-and-Plotting-Density
RHO_CONST = 1029.0
OMEGA     = (2*np.pi)/86164

def calddz(Q, xgcm_grid, ecco_grid, interp=False):
    dQdz = xgcm_grid.diff(Q, axis='Z', boundary='fill') #/ ecco_grid.drC

    if interp:
        dQdz = xgcm_grid.interp(dQdz, 'Z')
    
    return dQdz

def calGrad(Q, xgcm_grid, ecco_grid):

    dQdx = (xgcm_grid.diff(Q, axis="X", boundary='extend')) / ecco_grid.dxC
    dQdy = (xgcm_grid.diff(Q, axis="Y", boundary='extend')) / ecco_grid.dyC

    #dQdx.data = dQdx.values
    #dQdy.data = dQdy.values

    #_interp = xgcm_grid.interp_2d_vector({"X": dQdx, "Y": dQdy}, boundary='extend')

    dQdx = xgcm_grid.interp(dQdx, "X")#_interp["X"]
    dQdy = xgcm_grid.interp(dQdy, "Y")
    #dQdy = _interp["Y"]

    return dQdx, dQdy




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
    fixed_MLD = -1.0, # if this option is positive, then MLD will be assignd with this value
):

    print("[processECCO] Target datetime: ", target_datetime)
    print("[processECCO] Output filename: ", output_filename)


    if fixed_MLD <= 0:
        
        print("Notice: the parameter `fixed_MLD` is negative. Later on the mixed-layer depth will be determined dynamically.")

    else:
        
        print("Notice: the parameter `fixed_MLD` = %f is positive. Later on the mixed-layer depth will be assigned to this value." % (fixed_MLD, ))

    beg_datetime = target_datetime # - timedelta(days=1)

    snp_varnames = ["THETA", "SALT", "ETAN"]
    ave_varnames = ["MXLDEPTH", "Gs_ttl", "Gs_hadv", "Gs_vadv", "Gs_hdiff", "Gs_vdiff", "Gs_frc_sw", "Gs_frc_lw", "Gs_frc_sh", "Gs_frc_lh", "Gs_frc_fwf", "THETA", "SALT", "UVEL", "VVEL", "RHOAnoma", "PHIHYDcR"]

    ds = ECCO_helper.loadECCOData_continuous(
        beg_datetime = beg_datetime,
        ndays = 1,
        snp_varnames = snp_varnames,
        ave_varnames = ave_varnames,
    )

    xgcm_grid = ecco.get_llc_grid(ds)
    ecco_grid = ECCO_helper.getECCOGrid()
    
    s_snp = 1 + ds.ETAN_snp / ecco_grid.Depth
    sTHETA_snp = ds.THETA_snp * s_snp
    sSALT_snp  = ds.SALT_snp  * s_snp

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

    ML_snp_varnames = ["MLDs_snp", "MLTs_snp", "MLSs_snp"]
    ML_ave_varnames = [
        "MLTs", "dMLTsdt", "MLGs_ttl",
        "MLGs_hadv", "MLGs_vadv", "MLGs_adv",
        "MLGs_vdiff", "MLGs_hdiff",
        "MLGs_frc_sw", "MLGs_frc_lw", "MLGs_frc_sh", "MLGs_frc_lh", "MLGs_frc_fwf",
    ]

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

            if fixed_MLD <= 0:

                ML_snp["MLDs_snp"][s, l, :, :] = findMLD_rho(
                    rho_snp[s, :, l, :, :].to_numpy(),
                    z_T,
                    mask=mask2D[l],
                    Nz_bot=Nz_bot[l],
                    dev=MLD_dev,
                )
                
            else:

                ML_snp["MLDs_snp"][s, l, :, :] = fixed_MLD

            ML_snp["MLTs_snp"][s, l, :, :] = computeMLMean(
                sTHETA_snp[s, :, l, :, :].to_numpy(),
                ML_snp["MLDs_snp"][s, l, :, :].to_numpy(),
                z_W,
                mask=mask2D[l]
            )

            ML_snp["MLSs_snp"][s, l, :, :] = computeMLMean(
                sSALT_snp[s, :, l, :, :].to_numpy(),
                ML_snp["MLDs_snp"][s, l, :, :].to_numpy(),
                z_W,
                mask=mask2D[l]
            )

    dt = 86400.0
    #xgcm_grid.diff(ds.time_snp, 'T', boundary='fill', fill_value=np.nan).astype('f4') / 1e9 # nanosec to sec 
   
    ML_ave["dMLTsdt"][0, :, :, :] = (ML_snp["MLTs_snp"][1, :, :, :] - ML_snp["MLTs_snp"][0, :, :, :]) / dt
    ML_ave["MLTs"][0, :, :, :] = (ML_snp["MLTs_snp"][0, :, :, :] + ML_snp["MLTs_snp"][1, :, :, :]) / 2.0


    # compute variable in the middle of time_bnds
    MLDs_snp = [ ML_snp["MLDs_snp"][s, :, :, :].to_numpy() for s in range(2) ]
    for phy_proc in ["ttl", "hadv", "vadv", "vdiff", "hdiff", "frc_sw", "frc_lw", "frc_sh", "frc_lh", "frc_fwf"]:
        
        varname = "Gs_%s" % phy_proc
        ML_varname = "ML%s" % varname
        for l in range(Nl):
                
            ML_ave[ML_varname][0, l, :, :] = computeMLMean(
                ds[varname][0, :, l, :, :].to_numpy(),
                MLDs_snp[1][l, :, :],
                z_W,
                mask=mask2D[l]
            )

    # Compute entrainment term explicitly
    ML_ave["MLGs_ent"] = xr.zeros_like(sample2D_ave).rename("MLGs_ent") 
    for l in range(Nl):
            
        sTHETA = sTHETA_snp[0, :, l, :, :].to_numpy()
        ML_ave["MLGs_ent"][0, l, :, :] = ( computeMLMean(
            sTHETA,
            MLDs_snp[1][l, :, :],
            z_W,
            mask=mask2D[l]
        ) -  computeMLMean(
            sTHETA,
            MLDs_snp[0][l, :, :],
            z_W,
            mask=mask2D[l]
        )) / dt


    # Additional diagnostic variables 
    ML_ave["MLGs_adv"][:, :, :, :] = ML_ave["MLGs_hadv"] + ML_ave["MLGs_vadv"]
    ML_ave["dMLTsdt_res"] = ML_ave["dMLTsdt"] - (
              ML_ave["MLGs_ent"] 
            + ML_ave["MLGs_hadv"]
            + ML_ave["MLGs_vadv"]
            + ML_ave["MLGs_hdiff"]
            + ML_ave["MLGs_vdiff"]
            + ML_ave["MLGs_frc_sw"]
            + ML_ave["MLGs_frc_lw"]
            + ML_ave["MLGs_frc_sh"]
            + ML_ave["MLGs_frc_lh"]
            + ML_ave["MLGs_frc_fwf"]
    )

    ML_ave["dMLTsdt_res"] = ML_ave["dMLTsdt_res"].rename("dMLTsdt_res")

    print("Compute physical Z terms.")

    s_0 = s_snp[0, :, : ,:]
    s_1 = s_snp[1, :, : ,:]
    
    for phy_proc in [
        "ttl",
        "ent", "hadv", "vadv", "vdiff", "hdiff",
        "frc_sw", "frc_lw", "frc_sh", "frc_lh", "frc_fwf",
    ]:
        zstar_name = "MLGs_%s" % (phy_proc,)
        z_name     = "MLG_%s" % (phy_proc,)
        ML_ave[z_name] = ML_ave[zstar_name] /  s_0



    MLT_snp = xr.zeros_like(s_snp).rename('MLT_snp')
    MLS_snp = xr.zeros_like(s_snp).rename('MLS_snp')

    for s in range(2):
        for l in range(Nl):

            MLT_snp[s, l, :, :] = computeMLMean(
                ds.THETA_snp[s, :, l, :, :].to_numpy(),
                MLDs_snp[s][l, :, :],
                z_W,
                mask=mask2D[l]
            )

            MLS_snp[s, l, :, :] = computeMLMean(
                ds.SALT_snp[s, :, l, :, :].to_numpy(),
                MLDs_snp[s][l, :, :],
                z_W,
                mask=mask2D[l]
            )

    ML_ave["dMLTdt"] = (MLT_snp[1, :, :, :] - MLT_snp[0, :, :, :]) / dt
    ML_ave["dMLSdt"] = (MLS_snp[1, :, :, :] - MLS_snp[0, :, :, :]) / dt
    
    #ML_ave["MLT"] = (MLT_snp[0, :, :, :] + MLT_snp[1, :, :, :]) / 2.0
    #ML_ave["MLT"] = ML_ave["MLT"].rename()

    ML_ave["MLG_rescale"] = - MLT_snp[1, :, :, :] / s_0 * (s_1 - s_0) / dt
    ML_ave["MLG_rescale"] = ML_ave["MLG_rescale"].rename('MLG_rescale')
    
    ML_ave["dMLTdt_res"] = ML_ave["dMLTdt"] - (
              ML_ave["MLG_ent"] 
            + ML_ave["MLG_hadv"]
            + ML_ave["MLG_vadv"]
            + ML_ave["MLG_hdiff"]
            + ML_ave["MLG_vdiff"]
            + ML_ave["MLG_frc_sw"]
            + ML_ave["MLG_frc_lw"]
            + ML_ave["MLG_frc_sh"]
            + ML_ave["MLG_frc_lh"]
            + ML_ave["MLG_frc_fwf"]
            + ML_ave["MLG_rescale"]
    )

    ML_ave["dMLTdt_res"] = ML_ave["dMLTdt_res"].rename("dMLTdt_res")

    ####################################################
    # compute MLU, MLV, dMLTdx, dMLTdy, dMLSdx, dMLSdy
    #         U_g, V_g
    # with MLD in the end of the time interval
    ####################################################

    def calGrad_wrap(Q):
        return calGrad(Q, xgcm_grid, ecco_grid)
    
    ML_ave2 = {}


    for varname in ["MLT", "MLS", "MLU", "MLV", "U_g", "V_g", "dMLTdx", "dMLTdy", "dMLSdx", "dMLSdy", "dTdz_b", "dSdz_b"]:
         ML_ave2[varname] = xr.zeros_like(ML_ave["MLG_ttl"]).rename(varname)

    
    # Compute geostrophic balance. The code is from 
    # https://ecco-v4-python-tutorial.readthedocs.io/Geostrophic_balance.html#Right-hand-side

    dens  = ds.RHOAnoma + RHO_CONST
    pressanom = ds.PHIHYDcR

    dpdx, dpdy = calGrad_wrap(RHO_CONST * pressanom)

    f_co = 2 * OMEGA * np.sin(ecco_grid.YC * np.pi / 180)
    U_g = - dpdy / RHO_CONST / f_co
    V_g =   dpdx / RHO_CONST / f_co


    # Compute velocity at T grid
    vel_interp = xgcm_grid.interp_2d_vector({'X':ds.UVEL,'Y':ds.VVEL},boundary='extend')
    UVEL = vel_interp['X']
    VVEL = vel_interp['Y']

    for MLvarname, var in {
        "MLT" : ds["THETA"],
        "MLS" : ds["SALT"],
        "MLU" : UVEL,
        "MLV" : VVEL,
        "U_g" : U_g,
        "V_g" : V_g,
    }.items():

        for l in range(Nl):

            MLD = MLDs_snp[1][l, :, :]

            ML_ave2[MLvarname][0, l, :, :] = computeMLMean(
                var[0, :, l, :, :].to_numpy(),
                MLD,
                z_W,
                mask=mask2D[l]
            )

    
    ML_ave2["dMLTdx"], ML_ave2["dMLTdy"] = calGrad_wrap(ML_ave2["MLT"])
    ML_ave2["dMLSdx"], ML_ave2["dMLSdy"] = calGrad_wrap(ML_ave2["MLS"])



    for l in range(Nl):
        MLD = MLDs_snp[1][l, :, :]

        dTdz = calculus_tools.W_ddz_T(ds["THETA"][0, :, l, :, :].to_numpy(), z_W=z_W)
        dSdz = calculus_tools.W_ddz_T(ds["SALT"][0, :, l, :, :].to_numpy(), z_W=z_W)

        ML_ave2["dTdz_b"][0, l, :, :] = evalAtMLD_W(dTdz, MLD, z_W, mask=mask2D[l])
        ML_ave2["dSdz_b"][0, l, :, :] = evalAtMLD_W(dSdz, MLD, z_W, mask=mask2D[l])
    

    
    output_data = []

    for k, var in ML_ave.items():
        var = var.rename(k)
        output_data.append(var)

    for k, var in ML_ave2.items():
        var = var.rename(k)
        output_data.append(var)

    output_data.append(ML_snp["MLDs_snp"][1:2, :, :, :].rename("MLD"))

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

if __name__ == "__main__" : 

    print("*** This is for testing ***") 

    target_datetime = datetime(2006, 12, 29)
    output_filename = "data/ECCO_LLC/%s/%s" % ECCO_helper.getECCOFilename("MLT", "DAILY", target_datetime)
    
    processECCO(
        target_datetime,
        output_filename,
    ) 
    
