from multiprocessing import Pool
import datetime
from pathlib import Path
import os.path
import os
import netCDF4

def pleaseRun(cmd):
    print(">> %s" % cmd)
    os.system(cmd)

year_rng = [2017, 2019]
beg_time = datetime.datetime(year_rng[0],   1, 1)
end_time = datetime.datetime(year_rng[1],   1, 1)


g0 = 9.81
