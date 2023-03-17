from multiprocessing import Pool
import numpy as np
from datetime import (datetime, timedelta, timezone)
from pathlib import Path
import os.path
import os
import netCDF4
import argparse
print("Loading libraries completed.")

def pleaseRun(cmd):
    print(">> %s" % cmd)
    os.system(cmd)


parser = argparse.ArgumentParser(
                    prog = 'plot_skill',
                    description = 'Plot prediction skill of GFS on AR.',
)

parser.add_argument('--ncpu',   type=int, help='Date string: yyyy-mm-dd', default=4)
parser.add_argument('--beg-year',   type=int, help='Date string: yyyy-mm-dd', required=True)
parser.add_argument('--end-year',   type=int, help='Date string: yyyy-mm-dd', required=True)
parser.add_argument('--output-dir', type=str, help='Output directory', default="")
parser.add_argument('--lat-rng',    type=float, nargs=2, help='Latitude  range', required=True)
parser.add_argument('--lon-rng',    type=float, nargs=2, help='Longitude range. 0-360', required=True)
parser.add_argument('--lat-nbox',   type=int, help='Latitude  range', required=True)
parser.add_argument('--lon-nbox',   type=int, help='Longitude range. 0-360', required=True)
parser.add_argument('--mask-ERA5',  type=str, help='mask file of ERA5', required=True)
parser.add_argument('--mask-ECCO',  type=str, help='mask file of ECCO', required=True)
parser.add_argument('--ignore-empty-box',  action="store_true")
args = parser.parse_args()

print(args)

class JOB:

    def __init__(self, year):
        self.year = year


    def work(self):
     
        print("##### Doing work of year %04d #####" % (self.year,))
 
        target_filename = "AR_statistics_yr%04d.nc" % (self.year,)
        target_fullname = "%s/%s" % (args.output_dir, target_filename)

        already_exists = os.path.isfile(target_fullname)

        if already_exists:

            print("[%04d] File '%s' already exists. Skip." % (self.year, target_fullname, ))

        else:

            print("[%04d] Now generate file: %s" % (self.year, target_fullname,))


            cmd = [
                "python3", "construct_timeseries_by_boxes_2.py",
                "--beg-year", "%d" % (self.year,),
                "--end-year", "%d" % (self.year,),
                "--output-dir", args.output_dir,
                "--output-filename", target_filename,
                "--lat-rng", "%f %f" % tuple(args.lat_rng),
                "--lon-rng", "%f %f" % tuple(args.lon_rng),
                "--lat-nbox", "%d" % args.lat_nbox,
                "--lon-nbox", "%d" % args.lon_nbox,
                "--mask-ERA5", args.mask_ERA5,
                "--mask-ECCO", args.mask_ECCO,
            ]

            if args.ignore_empty_box:
                cmd.append("--ignore-empty-box")

            cmd = " ".join(cmd)

            pleaseRun(cmd)


def wrap_retrieve(job):

    job.work()


jobs = []
for y in range(args.beg_year, args.end_year+1):

    jobs.append(JOB(y))


with Pool(processes=args.ncpu) as pool:

    result = pool.map(wrap_retrieve, jobs)

