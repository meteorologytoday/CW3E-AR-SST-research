from datetime import (datetime, timedelta)
import os.path
import xarray as xr
import xgcm
import ecco_v4_py as ecco
import numpy as np

rhoconst = 1029.0
c_p = 3994.0
R = 0.62
zeta1 = 0.06
zeta2 = 20.0

ECCO_root_dir = "data/ECCO_LLC"
ECCO_grid_dir = "ECCO_L4_GEOMETRY_LLC0090GRID_V4R4"
ECCO_grid_filename = "GRID_GEOMETRY_ECCO_V4r4_native_llc0090.nc"

ECCO_mapping = {

    "TEMP_SALINITY" : {
        "fileprefix": "OCEAN_TEMPERATURE_SALINITY",
        "varnames": ["THETA", "SALT"],
    },

    "OCEAN_3D_TEMPERATURE_FLUX" : {
        "fileprefix": "OCEAN_3D_TEMPERATURE_FLUX",
        "varnames": ["ADVx_TH", "ADVy_TH", "ADVr_TH", "DFxE_TH", "DFyE_TH", "DFrE_TH", "DFrI_TH"],
    },

    "MIXED_LAYER_DEPTH" : {
        "fileprefix": "OCEAN_MIXED_LAYER_DEPTH",
        "varnames": ["MXLDEPTH",],
    },

    
    # oceQnet = EXFls + EXFlh - (EXFlwnet + EXFswnet)
    # oceQsw = - EXFswnet

    # EXFqnet = - oceQnet = EXFlwnet + EXFswnet - EXFlh - EXFhs
    # TFLUX = surForcT + oceQsw + oceFreez + [PmEpR*SST]*Cp
    # oceFWflx = [PmEpR]

    # In our case where sea ice does not involve
    # TFLUX = oceQsw + EXFhl + EXFhs - EXFlwnet + [PmEpR*SST]*Cp
    #                                                  |
    #                                                  +--> the loss/gain of ocean mass * c_p

    "HEAT_FLUX" : {
        "fileprefix" : "OCEAN_AND_ICE_SURFACE_HEAT_FLUX",
        "varnames" : ["oceQsw", "TFLUX", "EXFhs", "EXFhl", "EXFlwnet"],
    },

    "OCEAN_VEL" : {
        "fileprefix" : "OCEAN_VELOCITY",
        "varnames" : ["UVEL", "VVEL", "WVEL"],
    },

    "SSH" : {
        "fileprefix" : "SEA_SURFACE_HEIGHT",
        "varnames" : ["SSH", "ETAN"],
    },

    "POSTPROC_GS_TERMS" : {
        "fileprefix": "GS_TERMS",
        "varnames": ["Gs_ttl", "Gs_hadv", "Gs_vadv", "Gs_hdiff", "Gs_vdiff",
                     "Gs_frc_sw", "Gs_frc_lw", "Gs_frc_sh", "Gs_frc_lh", "Gs_frc_fwf",
                     "Gs_sum", "Gs_res"],
    },


    "POSTPROC_MXLANA" : {
        "fileprefix": "MXLANA",
        "varnames": [

            "MLT", "dMLTdT",

            "MLG_ttl", "MLG_hadv", "MLG_vadv", "MLG_hdiff", "MLG_vdiff",
            "MLG_frc_sw", "MLG_frc_lw", "MLG_frc_sh", "MLG_frc_lh", "MLG_frc_fwf",
            "MLG_sum", "MLG_res",

            "MLGs_ttl", "MLGs_hadv", "MLGs_vadv", "MLGs_hdiff", "MLGs_vdiff",
            "MLGs_frc_sw", "MLGs_frc_lw", "MLGs_frc_sh", "MLGs_frc_lh", "MLGs_frc_fwf",
            "MLGs_sum", "MLGs_res",

        ],
    },


}


map_varname_pathinfo = {}

for dirmidfix, info in ECCO_mapping.items():
    for varname in info["varnames"]:
        map_varname_pathinfo[varname] = {
            "dirmidfix" : dirmidfix,
            "fileprefix" : info["fileprefix"],
        }



grid_mapping = {

    "LATLON" : {
        "dir" : "05DEG",
        "file" : "latlon_0p50deg",
    },

    "LLC" : {
        "dir" : "LLC0090GRID",
        "file" : "native_llc0090",
    },

}

time_character_mapping = {

    "DAILY" : {
        "dir" : "DAILY",
        "file" : "day_mean",
    },

    "SNAPSHOT" : {
        "dir" : "SNAPSHOT",
        "file" : "snap",
    },

}


def getECCOGrid():
        
    print("%s/%s" % (ECCO_root_dir, ECCO_grid_dir))
    print(ECCO_grid_filename)

    ecco_grid = ecco.load_ecco_grid_nc("%s/%s" % (ECCO_root_dir, ECCO_grid_dir), ECCO_grid_filename)
    

    return ecco_grid 


def getECCOFilename(varname, time_character, target_datetime, grid="LLC", version=4, release=4):


    if time_character not in ["DAILY", "SNAPSHOT"]:
        raise Exception("Unknown `time_character`: %s" % (str(time_character),))

    if grid not in ["LLC", "LATLON"]:
        raise Exception("Unknown `grid`: %s" % (str(grid),))


    dirrelease = "V%dR%d" % (version, release)
    dirsuffix = "%s_%s_%s" % (grid_mapping[grid]["dir"], time_character_mapping[time_character]["dir"], dirrelease)
    dirmidfix = map_varname_pathinfo[varname]["dirmidfix"]

    dirname = "ECCO_L4_%s_%s" % (dirmidfix, dirsuffix)

   

    if time_character == "DAILY" :
        time_str_suffix = ""

    elif time_character == "SNAPSHOT":
        time_str_suffix = "T000000"
        
    time_str = "%s%s" % ( target_datetime.strftime("%Y-%m-%d"), time_str_suffix)
    filerelase = "V%dr%d" % (version, release)
    fileprefix = map_varname_pathinfo[varname]["fileprefix"]
    filesuffix = "%s_%s_ECCO_%s_%s" % (time_character_mapping[time_character]["file"], time_str, filerelase, grid_mapping[grid]["file"])
    filename = "%s_%s.nc" % (fileprefix, filesuffix)


    return dirname, filename




def loadECCOData_continuous(
    beg_datetime,
    ndays = 1,
    snp_varnames = [],
    ave_varnames = [],
    return_xgcm_grid = False,
):

    full_list = []
    
    # ndays + 1 SNAPSHOTS are needed
    for _t in range(ndays+1):

        _now_datetime = beg_datetime + timedelta(days=1) * _t

        for varname in snp_varnames:

            dirname, filename = getECCOFilename(varname, "SNAPSHOT", _now_datetime)
            fullpath = "data/ECCO_LLC/%s/%s" % (dirname, filename)

            if not os.path.isfile(fullpath):
                raise Exception("File %s does not exist." % (fullpath,))

            new_varname = "%s_snp" % (varname,)
            _tmp = xr.open_dataset(fullpath)[[varname,]].rename({'time':'time_snp', varname : new_varname})

            full_list.append(_tmp) 
 
    for _t in range(ndays):
        
        _now_datetime = beg_datetime + timedelta(days=1) * _t

        for varname in ave_varnames:

            dirname, filename = getECCOFilename(varname, "DAILY", _now_datetime)
            fullpath = "data/ECCO_LLC/%s/%s" % (dirname, filename)
            if not os.path.isfile(fullpath):
                raise Exception("File %s does not exist." % (fullpath,))


            _tmp = xr.open_dataset(fullpath)[[varname,]]
            
            full_list.append(_tmp)


    ds = xr.merge(full_list)
    ds.time_snp.attrs['c_grid_axis_shift'] = - 0.5

    return ds




# Reference: https://ecco-v4-python-tutorial.readthedocs.io/ECCO_v4_Heat_budget_closure.html
def computeTendency(target_datetime, grid=None):
    snp_varnames = ["THETA", "ETAN"]
    ave_varnames = [ 
        "TFLUX", "oceQsw", "EXFlwnet", "EXFhl", "EXFhs",
        "ADVx_TH", "ADVy_TH", "ADVr_TH", "DFxE_TH", "DFyE_TH", "DFrE_TH", "DFrI_TH", ] 

    ds = loadECCOData_continuous(
        beg_datetime = target_datetime,
        ndays = 1,
        snp_varnames = snp_varnames,
        ave_varnames = ave_varnames,
    )
    xgcm_grid = ecco.get_llc_grid(ds)

 
    delta_t = xgcm_grid.diff(ds.time_snp, 'T', boundary='fill', fill_value=np.nan).astype('f4') / 1e9 # nanosec to sec 

    ecco_grid = getECCOGrid()
    vol = (ecco_grid.rA*ecco_grid.drF*ecco_grid.hFacC).transpose('tile','k','j','i')

    s_star_snap = 1 + ds.ETAN_snp / ecco_grid.Depth
    sTHETA = ds.THETA_snp * s_star_snap
    G_ttl = xgcm_grid.diff(sTHETA, 'T', boundary='fill', fill_value=0.0)/delta_t

    ADVxy_diff = xgcm_grid.diff_2d_vector({'X' : ds.ADVx_TH, 'Y' : ds.ADVy_TH}, boundary = 'fill')

    adv_hConvH = (-(ADVxy_diff['X'] + ADVxy_diff['Y']))

    ADVr_TH = ds.ADVr_TH.transpose('time','tile','k_l','j','i')
    adv_vConvH = xgcm_grid.diff(ADVr_TH, 'Z', boundary='fill')

    G_hadv = adv_hConvH / vol
    G_vadv = adv_vConvH / vol


    DFxyE_diff = xgcm_grid.diff_2d_vector({'X' : ds.DFxE_TH, 'Y' : ds.DFyE_TH}, boundary = 'fill')

    # Convergence of horizontal diffusion (degC m^3/s)
    dif_hConvH = (-(DFxyE_diff['X'] + DFxyE_diff['Y']))

    # Load monthly averages of vertical diffusive fluxes
    DFrE_TH = ds.DFrE_TH.transpose('time','tile','k_l','j','i')
    DFrI_TH = ds.DFrI_TH.transpose('time','tile','k_l','j','i')

    # Convergence of vertical diffusion (degC m^3/s)
    dif_vConvH = xgcm_grid.diff(DFrE_TH, 'Z', boundary='fill') + xgcm_grid.diff(DFrI_TH, 'Z', boundary='fill')

    G_hdiff = dif_hConvH / vol
    G_vdiff = dif_vConvH / vol

    Z = ecco_grid.Z.load()
    RF = np.concatenate([ecco_grid.Zp1.values[:-1],[np.nan]])

    q1 = R*np.exp(1.0/zeta1*RF[:-1]) + (1.0-R)*np.exp(1.0/zeta2*RF[:-1])
    q2 = R*np.exp(1.0/zeta1*RF[1:]) + (1.0-R)*np.exp(1.0/zeta2*RF[1:])


    zCut = np.where(Z < -200)[0][0]
    q1[zCut:] = 0
    q2[zCut-1:] = 0


    q1 = xr.DataArray(q1,coords=[Z.k],dims=['k'])
    q2 = xr.DataArray(q2,coords=[Z.k],dims=['k'])

    mskC = ecco_grid.hFacC.copy(deep=True).load()

    # Change all fractions (ocean) to 1. land = 0
    mskC.values[mskC.values>0] = 1
    forcH_subsurf_sw = ((q1*(mskC==1)-q2*(mskC.shift(k=-1)==1))*ds.oceQsw).transpose('time','tile','k','j','i')


    forcH_surf_sw = ( (q1[0]-q2[0]) * ds.oceQsw
              *mskC[0]).transpose('time','tile','j','i').assign_coords(k=0).expand_dims('k')

    forcH_sw = xr.concat([forcH_surf_sw,forcH_subsurf_sw[:,:,1:]], dim='k').transpose('time','tile','k','j','i')

    forcH_surf_nonsw_shape = (( ds.TFLUX * 0 + 1 )
              *mskC[0]).transpose('time','tile','j','i').assign_coords(k=0).expand_dims('k')
    forcH_subsurf_nonsw_shape = forcH_subsurf_sw * 0

    forcH_nonsw_shape = xr.concat([forcH_surf_nonsw_shape,forcH_subsurf_nonsw_shape[:,:,1:]], dim='k').transpose('time','tile','k','j','i')

    
    EXFfwf = ds.TFLUX - ds.oceQsw - ds.EXFhl - ds.EXFhs + ds.EXFlwnet

    G_frc_sw  = forcH_sw                            / (rhoconst*c_p) / (ecco_grid.hFacC*ecco_grid.drF)
    G_frc_lw  = forcH_nonsw_shape * (- ds.EXFlwnet) / (rhoconst*c_p) / (ecco_grid.hFacC*ecco_grid.drF)
    G_frc_sh  = forcH_nonsw_shape * ds.EXFhs        / (rhoconst*c_p) / (ecco_grid.hFacC*ecco_grid.drF)
    G_frc_lh  = forcH_nonsw_shape * ds.EXFhl        / (rhoconst*c_p) / (ecco_grid.hFacC*ecco_grid.drF)
    G_frc_fwf = forcH_nonsw_shape * EXFfwf          / (rhoconst*c_p) / (ecco_grid.hFacC*ecco_grid.drF)


    G_sum = G_hadv + G_vadv + G_hdiff + G_vdiff + G_frc_sw + G_frc_lw + G_frc_sh + G_frc_lh + G_frc_fwf
    G_res = G_sum - G_ttl

    result = {
        "Gs_ttl"     : G_ttl,
        "Gs_hadv"    : G_hadv,
        "Gs_vadv"    : G_vadv,
        "Gs_frc_sw"  : G_frc_sw,
        "Gs_frc_lw"  : G_frc_lw,
        "Gs_frc_sh"  : G_frc_sh,
        "Gs_frc_lh"  : G_frc_lh,
        "Gs_frc_fwf" : G_frc_fwf,
        "Gs_hdiff"   : G_hdiff,
        "Gs_vdiff"   : G_vdiff,
        "Gs_sum"     : G_sum,
        "Gs_res"     : G_res,
    }


    for k, v in result.items():
        
        result[k] = v.transpose('time', 'k', 'tile', 'j', 'i')


    return result


