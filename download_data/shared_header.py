from multiprocessing import Pool
import datetime
from pathlib import Path
import os.path
import os
import netCDF4

def pleaseRun(cmd):
    print(">> %s" % cmd)
    os.system(cmd)

beg_time = datetime.datetime(2016,   11,  1)
end_time = datetime.datetime(2022,   12, 31)

beg_time = datetime.datetime(2017,   11,  1)
end_time = datetime.datetime(2018,   12, 31)



g0 = 9.81
