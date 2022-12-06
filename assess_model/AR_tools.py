import numpy as np
from datetime import datetime

def detectAbove(_t, x, x_threshold):

    dtype_t = type(_t[0])

    if dtype_t is datetime:
        t = [ _tt.timestamp() for _tt in _t ]
    else:
        t = _t

    above = x >= x_threshold

    t_segs  = []
    seg_i = 0

    i     = 0
    flag = False

    t_seg = [np.nan, np.nan]
    for i in range(len(x)):

        if flag is False:
                    
            if x[i] >= x_threshold:
                flag = True
                if i == 0:
                    t_seg[0] = t[0]
                else:
                    t_seg[0] = t[i-1] + (t[i] - t[i-1]) * (x_threshold - x[i-1]) / (x[i] - x[i-1])
                
            #else: # flag is True, meaning it keeps to be above threshold
            #    if i == len(x) - 1: # last element
            #        t_seg[1] = t[-1]
            #        break

            #i#else: # x[i] < x_threshold
            
        else: # Flag is True
            
            if x[i] >= x_threshold:

                if i == len(x) - 1:
            
                    flag = False
                    t_seg[1] = t[-1]

            else:

                flag = False
                t_seg[1] = t[i] - (t[i] - t[i-1]) * (x_threshold - x[i]) / (x[i-1] - x[i])

        if not np.isnan(t_seg[0]) and not np.isnan(t_seg[1]):
            
            if dtype_t is datetime:
                t_seg[0] = datetime.fromtimestamp(t_seg[0])
                t_seg[1] = datetime.fromtimestamp(t_seg[1])

            t_segs.append(t_seg)

            t_seg = [np.nan, np.nan]

    return t_segs 
