import numpy as np

def makeDividedBoxes(lon_bnds, lat_bnds):
    
    Nx = len(lon_bnds) - 1
    Ny = len(lat_bnds) - 1
    

    boxes = []

    n = 0

    for j in range(Ny):
        for i in range(Nx):
            
            boxes.append({
                "numbering" : n,
                "polygon"   : {
                    "lon_bnds" : (lon_bnds[i]%360, lon_bnds[i+1]%360),
                    "lat_bnds" : (lat_bnds[j], lat_bnds[j+1]),
                },
                "label"     : "(%d, %d)" % (j, i),
            })
            
            n += 1
            
    return boxes
