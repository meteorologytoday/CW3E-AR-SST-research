import numpy as np
import datetime
from datetime import timedelta, datetime

def doy_noleap(t):

    ref_t = datetime(2022, t.month, t.day)

    return int(ref_t.strftime('%j'))

def doy_leap(t):
    ref_t = datetime(2020, t.month, t.day)
    return int(ref_t.strftime('%j'))



def decomposeClimAnom(ts: np.ndarray, xs: np.ndarray, assist = None):

    if ts.dtype != object:    
        ts_datetime = ts.astype(object)
    else:
        ts_datetime = ts


    if assist is None:

        tm = np.array([
            datetime(2021, 1, 1) + _t * timedelta(days=1)  for _t in range(365)
        ]).astype("datetime64[s]")


        doy  = np.zeros(len(ts_datetime), dtype=np.int32)
        skip = np.zeros(len(ts_datetime), dtype=np.bool_)
        cnt  = np.zeros(len(tm), dtype=np.int32)

        doy[:] = -1
        for i, t in enumerate(ts_datetime):

            m = t.month
            d = t.day


            
            if m == 2 and d == 29:
                skip[i] = True
            else:
                doy[i] = doy_noleap(t)
                cnt[doy-1] += 1

        assist = {
            'doy'  : doy,
            'skip' : skip,
            'cnt'  : cnt,
            'tm'   : tm,
        }

    tm = assist['tm']
    doy  = assist['doy']       
    cnt  = assist['cnt']       
    skip = assist['skip']       
    
    xm       = np.zeros((len(tm),))

    for i, t in enumerate(ts_datetime):

        if skip[i] :
            continue
        
        xm[doy[i]-1]  += xs[i]

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
            xa[i] = xs[i] - xm[doy[i]-1]

    tm = np.array([
        datetime(2021, 1, 1) + _t * timedelta(days=1)  for _t in range(len(tm))
    ]).astype("datetime64[s]")

    return tm, xm, xa, cnt, assist
