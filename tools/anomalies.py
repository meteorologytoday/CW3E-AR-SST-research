import numpy as np
import datetime

def doy_noleap(t):

    ref_t = datetime.datetime(2022, t.month, t.day)

    return int(ref_t.strftime('%j'))


def decomposeClimAnom(ts, xs):

    new_t    = np.zeros((365,))
    xm       = np.zeros((len(new_t),))
    cnt      = np.zeros((len(new_t),), dtype=int)

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
            xa[i] = np.nan
        
        else:
            doy = doy_noleap(t)
            xa[i] = xs[i] - xm[doy-1]


   
    return xm, xa, cnt
