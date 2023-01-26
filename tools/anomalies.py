import numpy as np
import datetime

def doy_noleap(t):

    ref_t = datetime.datetime(2022, t.month, t.day)

    return int(ref_t.strftime('%j'))

def doy_leap(t):
    ref_t = datetime.datetime(2020, t.month, t.day)
    return int(ref_t.strftime('%j'))



def decomposeClimAnom(ts, xs):

    tm    = np.zeros((365,))
    xm       = np.zeros((len(tm),))
    cnt      = np.zeros((len(tm),), dtype=int)

    for i, t in enumerate(ts):

        m = t.month
        d = t.day

        if m == 2 and d == 29:
            continue

        doy = doy_noleap(t)
        xm[doy-1]  += xs[i]
        cnt[doy-1] += 1

       
    xm /= cnt
        
    xa = np.zeros((len(xs),))
    for i, t in enumerate(ts):

        m = t.month
        d = t.day

        if m == 2 and d == 29:

            if i > 0 and i < (len(ts) - 1):
                xa[i] = (xa[i-1] + xa[i+1]) / 2.0
            else:
                xa[i] = np.nan
        
        else:
            doy = doy_noleap(t)
            xa[i] = xs[i] - xm[doy-1]


   
    return tm, xm, xa, cnt
