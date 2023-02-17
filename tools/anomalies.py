import numpy as np
import datetime
from datetime import timedelta, datetime

def doy_noleap(t):

    ref_t = datetime(2022, t.month, t.day)

    return int(ref_t.strftime('%j'))

def doy_leap(t):
    ref_t = datetime(2020, t.month, t.day)
    return int(ref_t.strftime('%j'))



def decomposeClimAnom(ts: np.ndarray, xs: np.ndarray):

    tm       = np.zeros((365,))
    xm       = np.zeros((len(tm),))
    cnt      = np.zeros((len(tm),), dtype=int)

    if ts.dtype != object:    
        ts_datetime = ts.astype(object)
    else:
        ts_datetime = ts

    #print(ts_datetime)

    for i, t in enumerate(ts_datetime):

        m = t.month
        d = t.day

        if m == 2 and d == 29:
            continue

        doy = doy_noleap(t)
        xm[doy-1]  += xs[i]
        cnt[doy-1] += 1

       
    xm /= cnt
        
    xa = np.zeros((len(xs),))
    for i, t in enumerate(ts_datetime):

        m = t.month
        d = t.day

        if m == 2 and d == 29:

            # Interpolate Feb 29 if Feb 28 and Mar 1 exist.
            if i > 0 and i < (len(ts) - 1):
                xa[i] = (xa[i-1] + xa[i+1]) / 2.0
            else:
                xa[i] = np.nan
        
        else:
            doy = doy_noleap(t)
            xa[i] = xs[i] - xm[doy-1]


    

    tm = np.array([
        datetime(2021, 1, 1) + _t * timedelta(days=1)  for _t in range(len(tm))
    ]).astype("datetime64[s]")

    return tm, xm, xa, cnt
