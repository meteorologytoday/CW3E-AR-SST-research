
import netCDF4
import numpy as np
import os
import re
import sys


# using August–Roche–Magnus formula as described in
# https://en.wikipedia.org/wiki/Clausius%E2%80%93Clapeyron_relation#August-Roche-Magnus_approximation
def getSaturatedVaporPressure(TMP):

    # in this formulation the temperature is in Celcius
    TMP_degC = TMP - 273.15
    p_es = 100 * 6.1094 * np.exp( 17.625 * TMP_degC / (TMP_degC + 243.04) ) # in Pa 
   
    return p_es


def convertRHtoSpecificHumidity(RH, TMP, pres):
    p_w = getSaturatedVaporPressure(TMP) * RH
    w = (p_w * 18) / ( ( pres - p_w ) * 28.9)
    q =  w / (w + 1)
    return q
    

def pleaseRun(cmd):
    print(">> %s" % cmd)
    os.system(cmd)


def convertGFSOutput(in_filename, out_filename, varnames=["UGRD", "VGRD", "HGT", "RH", "TMP"], compute_specific_humidity=True):


    tmp_filename = "%s.tmp" % (out_filename,)
    pleaseRun("wgrib2 -match ':(%s):[0-9]+ mb:' %s -netcdf %s" % ("|".join(varnames), in_filename, tmp_filename,))

    pattern = re.compile("(?P<VARNAME>[A-Za-z]+)_(?P<LEV>[0-9]+)mb")
    
    with netCDF4.Dataset(tmp_filename, "r") as ds_tmp:

        lat = ds_tmp.variables["latitude"][:]
        lon = ds_tmp.variables["longitude"][:]

        all_varnames = ds_tmp.variables.keys()

        lev = []
        print(all_varnames) 
        test_varname = varnames[0]

        # finding levels
        for varname in all_varnames:
            m = pattern.match(varname)
            if m is None:
                print("Varname %s does not match." % (varname,))
                continue

            if m.group("VARNAME") != test_varname:
                continue

            read_lev = int(m.group("LEV"))
            if read_lev < 200:
                continue
            
            lev.append(int(m.group("LEV")))

        lev.sort()
        lev = np.array(lev, dtype=np.float32) 

        if len(lev) == 0:
            raise Exception("Length of lev is 0.")
 
        with netCDF4.Dataset(out_filename, mode='w', format='NETCDF4_CLASSIC') as ds_out:

            time_dim = ds_out.createDimension('time', None)
            lat_dim  = ds_out.createDimension('lat', len(lat))
            lon_dim  = ds_out.createDimension('lon', len(lon)) 
            lev_dim  = ds_out.createDimension('lev', len(lev)) 


            var_lat = ds_out.createVariable('lat', np.float32, ('lat',))
            var_lon = ds_out.createVariable('lon', np.float32, ('lon',))
            var_lev = ds_out.createVariable('lev', np.float32, ('lev',))
            
            var_lat[:] = lat
            var_lon[:] = lon
            var_lev[:] = lev

            var_data = {}
            for varname in varnames:
                
                # Collect layers and put them into a 3D field
                x = np.zeros( (1, len(lev), len(lat), len(lon)), dtype=np.float32)
                for i, _lev in enumerate(lev):
                    x[0, i, :, :] = ds_tmp.variables["%s_%dmb" % (varname, _lev)][0, :, :]
                    
                if compute_specific_humidity and varname in ["TMP", "RH"]:
                    var_data[varname] = x

                _var = ds_out.createVariable(varname, np.float32, ('time', 'lev', 'lat', 'lon'))
                _var[0, :, :, :] = x

            if compute_specific_humidity:

                _q = ds_out.createVariable("q", np.float32, ('time', 'lev', 'lat', 'lon'))

                # The saturated vapor pressure formula is copied from wikipedia.
                # https://en.wikipedia.org/wiki/Clausius%E2%80%93Clapeyron_relation
                _TMP = var_data["TMP"]
                _RH  = var_data["RH"] / 100
                _lev = lev[:, None, None] * 100
                _q[0, :, :, :] = convertRHtoSpecificHumidity(_RH, _TMP, _lev)

    os.remove(tmp_filename)

if __name__ == "__main__":

    if len(sys.argv) < 3:
        print("Needs two arguments.")
   
    convertGFSOutput(sys.argv[1], sys.argv[2]) 
