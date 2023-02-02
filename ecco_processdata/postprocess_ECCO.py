with open("shared_header.py", "rb") as source_file:
    code = compile(source_file.read(), "shared_header.py", "exec")
exec(code)

import numpy as np
import argparse
import ECCO_helper
import xarray as xr

parser = argparse.ArgumentParser(
                    prog = 'postprocess_ECCO.py',
                    description = 'Postprocess ECCO data (Mixed-Layer integrated).',
)

parser.add_argument('--MLD-method', required=True, help="If set then use ECCO MLD instead.", type=str, choices=["RHO", "ECCO"])
args = parser.parse_args()
print(args)

output_root_dir = "data/ECCO_LLC"

MLD_dev = 0.03

nprocs = 1

beg_time = datetime.datetime(2017, 12, 25)
end_time = datetime.datetime(2018,  1,  1)

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

        _tmp = ECCO_helper.getECCOFilename("G_total", "DAILY", self.t)
        output_filename = "%s/%s/%s" % (output_root_dir, _tmp[0], _tmp[1])

        dir_name = os.path.dirname(output_filename)
        if not os.path.isdir(dir_name):
            print("Create dir: %s" % (dir_name,))
            Path(dir_name).mkdir(parents=True, exist_ok=True)


        all_exist = True
        for filename in [output_filename,]:
            all_exist = all_exist and os.path.isfile(filename)

        if all_exist:
            print("[%s] Data already exists. Skip." % (time_now_str, ))

        else:

            print("[%s] Now postprocessing..." % (time_now_str,))

            tends = ECCO_helper.computeTendency(self.t)

            ds = xr.Dataset(data_vars={})
            for varname, G in tends.items():
                ds[varname] = G

            ds.time.encoding = {}
            ds.reset_coords(drop=True)

            print("Output: ", output_filename)



            ds.to_netcdf(output_filename, format='NETCDF4')
            

 
def wrap_retrieve(job):

    try:
        job.work()
    
    except Exception as e:
        print("Exception!")
        print(e)

jobs = []
for d in range(total_days):
    new_d =  beg_time + datetime.timedelta(days=d)

    if 5 <= new_d.month and new_d.month <= 8 :
        continue
 
    # We need extra days to compute dSST/dt
    #if new_d.month == 4 and new_d.day != 1:
    #    continue
 
    #if new_d.month == 9 and new_d.day != 30:
    #    continue
    
    jobs.append(JOB(new_d))


print("Total jobs: %d" % (len(jobs),))

with Pool(processes=nprocs) as pool:

    result = pool.map(wrap_retrieve, jobs)



