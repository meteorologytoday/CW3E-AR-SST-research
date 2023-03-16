import numpy as np
import datetime
from datetime import timedelta, datetime

def doy_noleap(t):

    ref_t = datetime(2022, t.month, t.day)

    return int(ref_t.strftime('%j'))

def doy_leap(t):
    ref_t = datetime(2020, t.month, t.day)
    return int(ref_t.strftime('%j'))

def total_doy(t):
    return int((datetime(t.year+1, 1, 1) - datetime(t.year, 1, 1) ).total_seconds()/86400)

def fraction_of_year(t):

    return ( doy_leap(dt) - 1 ) / total_doy(dt)


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
                skip[i] = False
                doy[i] = doy_noleap(t)

                cnt[doy[i]-1] += 1

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

    cnt_nz = cnt != 0
    xm[cnt_nz] /= cnt[cnt_nz]
    xm[cnt == 0] = np.nan

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
           
            #print("xs[i]: ", xs[i], "; xm[doy[i]-1]: ", xm[doy[i]-1]) 
            xa[i] = xs[i] - xm[doy[i]-1]

    tm = np.array([
        datetime(2021, 1, 1) + _t * timedelta(days=1)  for _t in range(len(tm))
    ]).astype("datetime64[s]")

    return tm, xm, xa, cnt, assist



if __name__ == "__main__":

    # Creating fake data
    import pandas as pd
   
    beg_year = 2021
    end_year = 2030

    expected_cnt = end_year - beg_year + 1
 
    dates    = pd.date_range(start="%04d-01-01" % (beg_year,), end="%04d-12-31" % (end_year,), freq="D")
    noise    = np.random.randn(len(dates))
    raw_data = np.zeros((len(dates),))
    
    amp = 10.0
    wnm = 1.0

    for (i, _dt) in enumerate(dates):

        dt = pd.to_datetime(_dt)
        frac = fraction_of_year(dt)

        raw_data[i] = np.sin(2*np.pi * wnm * frac ) * amp + noise[i]

        
        #print("%d : %d-%d-%d, fraction of that year: %.4f" % (i, dt.year, dt.month, dt.day, ( doy_leap(dt) - 1 ) / total_doy(dt)))
        


    t_clim, clim, anom, cnt, _ = decomposeClimAnom(pd.to_datetime(dates), raw_data)

    if np.any(cnt != expected_cnt):
        print("[Debug] Expected cnt: ", expected_cnt)
        print("[Debug] Computed cnt = ", cnt)
        raise Exception("Count should be %d but we don't get it correct.")


    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(3, 1)
        
    ax[0].plot(dates, raw_data, label="Raw data")

    ax[1].plot(dates, noise, "k-", label="actual noise")
    ax[1].plot(dates, anom, "r--", label="computed anomalies")
    
    ax[2].plot(t_clim, clim, "k-", label="computed climatology")

    for _ax in ax:
        _ax.legend()

    plt.show()
     


