from multiprocessing import Pool
import datetime
from pathlib import Path
import os.path
import os
import netCDF4

def pleaseRun(cmd):
    print(">> %s" % cmd)
    os.system(cmd)

beg_time = datetime.datetime(1992,     9, 30)
end_time = datetime.datetime(2020,    12, 31)

#beg_time = datetime.datetime(1992,   10,  1)
#end_time = datetime.datetime(1993,    1,  1)




g0 = 9.81
