import numpy as np
import netCDF4
from buoyancy_linear import TS2b

fill_value=1e20


def processORA5ForNillerKrausMixedLayerDynamics(
    input_MLD_filename,
    input_T_filename,
    input_S_filename,
    output_filename,
    varname_MLD = "somxl",
    varname_T = "votemper",
    varname_S = "vosaline",
    varname_depth='deptht',
    varname_lat='nav_lat',
    varname_lon='nav_lon',
    sample_above_dist = 10.0,
    sample_below_dist = 20.0,
):


    print("INPUT FILE: ", input_MLD_filename)

    # 1. convert T and S to density
    # 2. Find the mixed layer depth
    # 3. Find the depth of ocean used to sample density to represent the jump of density beneath the mixed layer
    # 4. Compute drho, db and (h * db)
    # 5. Output file.

    with netCDF4.Dataset(input_MLD_filename, "r") as ds:
        MLD = ds.variables[varname_MLD][:, :, :] # (time, y, x)

    with netCDF4.Dataset(input_T_filename, "r") as ds:
        T = ds.variables[varname_T][:, :, :, :]  # (time, z, y, x)
        depth = ds.variables[varname_depth][:]
        lat = ds.variables[varname_lat][:]
        lon = ds.variables[varname_lon][:]

    with netCDF4.Dataset(input_S_filename, "r") as ds:
        S = ds.variables[varname_S][:, :, :, :]  # (time, z, y, x)


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
                
                h = MLD[t, j, i]

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




    output_vars = {
        'T_upper' : T_upper,
        'T_lower' : T_lower,
        'S_upper' : S_upper,
        'S_lower' : S_lower,
        'db'      : db,
        'MLD'     : MLD,
    }

    with netCDF4.Dataset(output_filename, mode='w', format='NETCDF4_CLASSIC') as ds: 

        x_dim    = ds.createDimension('y', Ny)
        y_dim    = ds.createDimension('x', Nx)
        time_dim = ds.createDimension('time', None)

        var_lat = ds.createVariable('lat', np.float32, ('y', 'x'))
        var_lon = ds.createVariable('lon', np.float32, ('y', 'x'))
        
        var_lat[:, :] = lat
        var_lon[:, :] = lon

        for k, d in output_vars.items():
            _var = ds.createVariable(k, np.float32, ('time', 'y', 'x'), fill_value=fill_value)
            _var[0, :, :] = d




    
