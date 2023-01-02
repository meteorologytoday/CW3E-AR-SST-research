import numpy as np
import earth_constants as ec 


def calAdv(Q, u, v, lat, lon, periodoc_lon=True, gap=30.0):

    r_E = ec.r_E
    
    Nlat, Nlon = u.shape

    deg2rad = np.pi / 180.0
    vort = np.zeros((Nlat, Nlon))
    _cos  = np.cos(deg2rad * lat)

    dy_2 = (( r_E * deg2rad ) * ( np.roll(lat, -1) - np.roll(lat, 1) ) )[:, np.newaxis]

    dlon_2 =  np.roll(lon, -1) - np.roll(lon, 1)

    # `gap` is just a number that is decidedly large to
    # determin that some dlon_2s are crossing 360
    dlon_2[dlon_2 < - gap] += 360.0
    
    
    if np.any(dlon_2 < - gap):
        raise Exception("Some dlon_2 is more negative than `-gap` = %f" % (gap,))
    
    if np.any(dlon_2 <= 0):
        raise Exception("Some dlon_2 is not positive. Please check.")
    
    if np.any(dlon_2 > gap):
        raise Exception("Some dlon_2 is larger than the gap %f . This should not happen please check." % gap)
    
    
    dlon_2 = dlon_2[np.newaxis, :]
    
    dx_2 = ( r_E * _cos[:, np.newaxis] * deg2rad ) * dlon_2
    
    dQdx = (np.roll(Q, -1, axis=1) - np.roll(Q, 1, axis=1)) / dx_2
    dQdy = (np.roll(Q, -1, axis=0) - np.roll(Q, 1, axis=0)) / dy_2
   
    #print("dQdx:")
    #print(dQdx)
 
    for _var in [dQdx, dQdy]:
        _var[0,  :] = np.nan
        _var[-1, :] = np.nan
        
        if periodoc_lon == False:
            _var[:,  0] = np.nan
            _var[:, -1] = np.nan
    
    adv = - ( u * dQdx + v * dQdy )
    
    return adv

    
