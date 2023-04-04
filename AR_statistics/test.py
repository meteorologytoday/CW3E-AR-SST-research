import ARstat_tool
import anomalies
import pandas as pd
import numpy as np

ds = ARstat_tool.loadDatasets("output_ECCO/1993-2017_10N-60N-n25_100E-100W-n80", list(range(1994, 2016+1)))

xs = ds.IVT[:, 20, 40].to_numpy()
ts = pd.DatetimeIndex(ds.time.to_numpy())

tm, xm, xa, cnt, _ = anomalies.decomposeClimAnom(ts, xs)



for t, _tm in enumerate(tm):
    print(_tm, " = ", xm[t])

#for t, _ts in enumerate(ts):
#    print(_ts, " = ", xa[t])
