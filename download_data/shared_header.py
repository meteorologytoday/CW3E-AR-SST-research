from multiprocessing import Pool
import datetime
from pathlib import Path
import os.path
import os
import netCDF4

def pleaseRun(cmd):
    print(">> %s" % cmd)
    os.system(cmd)

beg_time = datetime.datetime(1997,    9, 30)
end_time = datetime.datetime(2017,    4,  2)

#beg_time = datetime.datetime(1996,    9, 30)
#end_time = datetime.datetime(1997,   12, 31)




g0 = 9.81
