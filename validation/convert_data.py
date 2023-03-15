import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import buoyancy_nonlinear

# This program convert a mooring station data from irregular grid
# into regular grid. It also computes all the necessary fluxes which
# will be used to validate ECCO data

def loadData(
    input_dir,
    fn_temp,
    fn_salt,
    fn_sw,
    fn_lw,
    fn_shf,
    fn_lhf,
    fn_qnet,
    vn_temp,
    vn_salt,
    vn_sw,
    vn_lw,
    vn_shf,
    vn_lhf,
    vn_qnet,
    missing_value = 1e35,
):

    fn_z = [
        fn_temp,
        fn_salt,
    ]

    fn_sfc = [
        fn_sw,
        fn_lw,
        fn_shf,
        fn_lhf,
        fn_qnet,
    ]

    vn_sfc = [
        vn_sw,
        vn_lw,
        vn_shf,
        vn_lhf,
        vn_qnet,
    ]


 
    fn_z = ["%s/%s" % (input_dir, fn,) for fn in fn_z]
    fn_sfc = ["%s/%s" % (input_dir, fn,) for fn in fn_sfc]

    
    ds_z   = xr.open_mfdataset(fn_z)
    ds_sfc = xr.open_mfdataset(fn_sfc)

    ds_z   = ds_z.where( (ds_z[vn_temp] != missing_value) & (ds_z[vn_salt] != missing_value) )[dict(lat=0, lon=0)]
    ds_sfc = ds_sfc.where(ds_sfc[vn_qnet] != missing_value)[dict(lat=0, lon=0)]
    
    #chk = ( ds_sfc[vn_sw] - ds_sfc[vn_lw] - ds_sfc[vn_lhf] - ds_sfc[vn_shf]) - ds_sfc[vn_qnet]
    chk = ( ds_sfc[vn_sw] - ds_sfc[vn_lw] - ds_sfc[vn_lhf] - ds_sfc[vn_shf])
  

    rho = buoyancy_nonlinear.TS2b(ds_z[vn_temp], ds_z[vn_salt])
 
    #ds_z = ds_z.where(np.abs(chk) < 1)

    #i = 10
    #print("qnet : ", ds_sfc[vn_qnet][i].to_numpy()[0])
    #print("sw  : ", ds_sfc[vn_sw][i].to_numpy()[0])
    #print("lw  : ", ds_sfc[vn_lw][i].to_numpy()[0])
    #print("shf : ", ds_sfc[vn_shf][i].to_numpy()[0])
    #print("lhf : ", ds_sfc[vn_lhf][i].to_numpy()[0])

    #print(ds_z)
    #print(ds_sfc)
   
    fig, ax = plt.subplots(4, 1)
 
    mappable_T   = ax[0].contourf(ds_z.time, ds_z.depth, ds_z[vn_temp].transpose())
    mappable_S   = ax[1].contourf(ds_z.time, ds_z.depth, ds_z[vn_salt].transpose())
    mappable_rho = ax[2].contourf(ds_z.time, ds_z.depth, rho.transpose())
    #ax[1].plot(ds_sfc.time, ds_sfc[vn_sw])
    #ax[1].plot(ds_sfc.time, ds_sfc[vn_lw])
    #ax[1].plot(ds_sfc.time, ds_sfc[vn_lhf])
    #ax[1].plot(ds_sfc.time, ds_sfc[vn_shf])
    
    ax[3].twinx().plot(ds_sfc.time, chk, "k--", label="chk")
    ax[3].plot(ds_sfc.time, ds_sfc[vn_qnet], "r-", label="qnet")



    cb_T = plt.colorbar(mappable_T, ax = ax[0])
    cb_S = plt.colorbar(mappable_S, ax = ax[1])
    cb_rho = plt.colorbar(mappable_rho, ax = ax[2])

    
    ax[0].invert_yaxis()
    ax[1].invert_yaxis()
    ax[2].invert_yaxis()

    ax[0].set_title("Temperature")
    ax[1].set_title("Salinity")
    ax[2].set_title("Density")
    ax[3].set_title("Surface fluxes")

    plt.show()


def interpProfile(deps, Q):

    
    np.interp(x, xp, fp, left=None, right=None, period=None)


if __name__ == "__main__":

    print("Hi")

    suffix = "50n145w_dy"

    loadData(
        input_dir = "data/papa",
        fn_temp   = "t%s.cdf" % (suffix,),
        fn_salt   = "s%s.cdf" % (suffix,),
        fn_sw     = "swnet%s.cdf" % (suffix,),
        fn_lw     = "lwnet%s.cdf" % (suffix,),
        fn_shf    = "qsen%s.cdf" % (suffix,),
        fn_lhf    = "qlat%s.cdf" % (suffix,),
        fn_qnet   = "qnet%s.cdf" % (suffix,),
        vn_temp = "T_20",
        #vn_temp = "QT_5020",
        vn_salt = "S_41",
        vn_sw   = "SWN_1495",
        vn_lw   = "LWN_1136",
        vn_shf  = "QS_138",
        vn_lhf  = "QL_137",
        vn_qnet = "QT_210",
    )








