with open("shared_header.py", "rb") as source_file:
    code = compile(source_file.read(), "shared_header.py", "exec")
exec(code)

import postprocess_ECCO_tools


nprocs = 8

download_dir = "data/ECCO"

input_dir = "data/ECCO/ECCO_L4_TEMP_SALINITY_05DEG_DAILY_V4R4"
output_dir_0p50deg = "data/ECCO/processed_0p50deg"

input_filename_format = "OCEAN_TEMPERATURE_SALINITY_day_mean_%s_ECCO_V4r4_latlon_0p50deg"

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

        input_filename = input_filename_format % ( time_now_str, )
        input_filename = "%s/%s.nc" % (input_dir, input_filename)

        output_filename_0p50deg = "%s/ECCO_mixedlayer_0p50deg_%s.nc" % (output_dir_0p50deg, time_now_str, )
        already_exists = os.path.isfile(output_filename_0p50deg)

        if already_exists:
            print("[%s] Data already exists. Skip." % (time_now_str, ))

        else:

            print("[%s] Now process file: %s" % (time_now_str, input_filename,))

            postprocess_ECCO_tools.processECCO(
                input_filename,
                output_filename_0p50deg,
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



