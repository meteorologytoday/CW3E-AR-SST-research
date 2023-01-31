with open("shared_header.py", "rb") as source_file:
    code = compile(source_file.read(), "shared_header.py", "exec")
exec(code)

import numpy as np
import postprocess_ECCO_tools
import argparse

parser = argparse.ArgumentParser(
                    prog = 'postprocess_ECCO.py',
                    description = 'Postprocess ECCO data (Mixed-Layer integrated).',
)

parser.add_argument('--MLD-method', required=True, help="If set then use ECCO MLD instead.", type=str, choices=["RHO", "ECCO"])
args = parser.parse_args()
print(args)

MLD_dev = 0.03


nprocs = 8

download_dir = "data/ECCO"



input_dir_TS = "data/ECCO/ECCO_L4_TEMP_SALINITY_05DEG_DAILY_V4R4"
input_dir_VEL = "data/ECCO/ECCO_L4_OCEAN_VEL_05DEG_DAILY_V4R4"
input_dir_SSH = "data/ECCO/ECCO_L4_SSH_05DEG_DAILY_V4R4B"
input_dir_MLD = "data/ECCO/ECCO_L4_MIXED_LAYER_DEPTH_05DEG_DAILY_V4R4"


if args.MLD_method == "ECCO":
    print("Going to use ECCO provided mixed-layer depth.")
    output_dir_0p50deg = "data/ECCO/processed_0p50deg_ECCO-MLD"
elif args.MLD_method == "RHO":
    print("Going to use my algorithm to detect mixed-layer depth (dRHO = %f)." % (MLD_dev,))
    output_dir_0p50deg = "data/ECCO/processed_0p50deg_RHO-%02d-MLD" % (int(np.floor(MLD_dev * 100)))


input_filename_format_TS = "OCEAN_TEMPERATURE_SALINITY_day_mean_%s_ECCO_V4r4_latlon_0p50deg.nc"
input_filename_format_VEL = "OCEAN_VELOCITY_day_mean_%s_ECCO_V4r4_latlon_0p50deg.nc"
input_filename_format_SSH = "SEA_SURFACE_HEIGHT_day_mean_%s_ECCO_V4r4b_latlon_0p50deg.nc"

# This is only used if ECCO mixed layer depth is used
input_filename_format_MLD = "OCEAN_MIXED_LAYER_DEPTH_day_mean_%s_ECCO_V4r4_latlon_0p50deg.nc"

total_days = (end_time - beg_time).days

print("Going to download %d days of data." % (total_days,))

class JOB:

    def __init__(self, t):
        self.t = t


    def work(self):
        
        y = self.t.year
        m = self.t.month
        d = self.t.day

        time_now_str = "%04d-%02d-%02d" % (y, m, d)

        input_filename_TS = input_filename_format_TS % ( time_now_str, )
        input_filename_TS = "%s/%s" % (input_dir_TS, input_filename_TS)

        input_filename_SSH = input_filename_format_SSH % ( time_now_str, )
        input_filename_SSH = "%s/%s" % (input_dir_SSH, input_filename_SSH)

        input_filename_VEL = input_filename_format_VEL % ( time_now_str, )
        input_filename_VEL = "%s/%s" % (input_dir_VEL, input_filename_VEL)

        input_filename_MLD = input_filename_format_MLD % ( time_now_str, )
        input_filename_MLD = "%s/%s" % (input_dir_MLD, input_filename_MLD)




        output_filename_0p50deg = "%s/ECCO_mixedlayer_0p50deg_%s.nc" % (output_dir_0p50deg, time_now_str, )
        already_exists = os.path.isfile(output_filename_0p50deg)


        for input_filename in [input_filename_TS, input_filename_VEL, input_filename_SSH]:
            if not os.path.isfile(input_filename):
                print("Error: Input file %s does not exist." % (input_filename,))
                return

        if already_exists:
            print("[%s] Data already exists. Skip." % (time_now_str, ))

        else:

            print("[%s] Now postprocessing..." % (time_now_str,))

            postprocess_ECCO_tools.processECCO(
                input_filename_TS,
                input_filename_VEL,
                input_filename_SSH,
                output_filename_0p50deg,
                input_filename_MLD = input_filename_MLD,
                MLD_dev=MLD_dev,
            ) 


            # regrid
            

 
def wrap_retrieve(job):

    job.work()

jobs = []
for d in range(total_days):
    new_d =  beg_time + datetime.timedelta(days=d)

    if 5 <= new_d.month and new_d.month <= 8 :
        continue
 
    # We need extra days to compute dSST/dt
    if new_d.month == 4 and new_d.day != 1:
        continue
 
    if new_d.month == 9 and new_d.day != 30:
        continue
    
    jobs.append(JOB(new_d))


print("Total jobs: %d" % (len(jobs),))

for dir_name in [output_dir_0p50deg, ]:
    print("Create dir: %s" % (dir_name,))
    Path(dir_name).mkdir(parents=True, exist_ok=True)


with Pool(processes=nprocs) as pool:

    result = pool.map(wrap_retrieve, jobs)



