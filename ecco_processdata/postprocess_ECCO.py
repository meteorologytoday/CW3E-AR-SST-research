with open("shared_header.py", "rb") as source_file:
    code = compile(source_file.read(), "shared_header.py", "exec")
exec(code)

import numpy as np
import argparse
import ECCO_helper, ECCO_computeTendency
import xarray as xr
import postprocess_ECCO_tools

parser = argparse.ArgumentParser(
                    prog = 'postprocess_ECCO.py',
                    description = 'Postprocess ECCO data (Mixed-Layer integrated).',
)

parser.add_argument('--MLD-method', required=True, help="If set then use ECCO MLD instead.", type=str, choices=["RHO", "FIXED500m"])
parser.add_argument('--nproc', type=int, default=2)
args = parser.parse_args()
print(args)

output_root_dir = "data/ECCO_LLC"

MLD_dev = 0.03



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

        global args



        if args.MLD_method == "FIXED500m":
            
            extra_dirsuffix = "_500m"

        else:
            
            extra_dirsuffix = ""
        
        # Phase 1    
        # This one is computing G terms for each grid cell. Does not depend on MLD_method.
        _tmp = ECCO_helper.getECCOFilename("Gs_ttl", "DAILY", self.t)
        output_filename_G_terms = "%s/%s/%s" % (output_root_dir, _tmp[0], _tmp[1])

        if os.path.isfile(output_filename_G_terms):
            print("[%s] File %s already exists. Skip." % (time_now_str, output_filename_G_terms))

        else:
            dir_name = os.path.dirname(output_filename_G_terms)
            if not os.path.isdir(dir_name):
                print("Create dir: %s" % (dir_name,))
                Path(dir_name).mkdir(parents=True, exist_ok=True)
            
            
            print("[%s] Now computing G terms..." % (time_now_str,))
            
            tends = ECCO_computeTendency.computeTendency(self.t)

            ds = xr.Dataset(data_vars={})
            for varname, G in tends.items():
                ds[varname] = G

            ds.time.encoding = {}
            ds.reset_coords(drop=True)

            print("Output: ", output_filename_G_terms)
            ds.to_netcdf(output_filename_G_terms, format='NETCDF4')
           

        # Phase 2
        # This one computes the advection. Does not depend on MLD_method.
        _tmp = ECCO_helper.getECCOFilename("HADV_g", "DAILY", self.t)
        output_filename_ADV = "%s/%s/%s" % (output_root_dir, _tmp[0], _tmp[1])
        
        dir_name = os.path.dirname(output_filename_ADV)
        if not os.path.isdir(dir_name):
            print("Create dir: %s" % (dir_name,))
            Path(dir_name).mkdir(parents=True, exist_ok=True)
        
        if os.path.isfile(output_filename_ADV):
            print("[%s] File %s already exists. Skip." % (time_now_str, output_filename_ADV))

        else:
            print("[%s] Now compute the advection." % (time_now_str, ))
            
            ds = ECCO_computeTendency.computeTendencyAdv(self.t)
            print("Output: ", output_filename_ADV)
            ds.to_netcdf(output_filename_ADV, format='NETCDF4')
 
        # Phase 3
        # This one computes the mixed-layer integrated quantities
        _tmp = ECCO_helper.getECCOFilename("MLT", "DAILY", self.t, extra_dirsuffix=extra_dirsuffix)
        output_filename_MXLANA = "%s/%s/%s" % (output_root_dir, _tmp[0], _tmp[1])
        
        dir_name = os.path.dirname(output_filename_MXLANA)
        if not os.path.isdir(dir_name):
            print("Create dir: %s" % (dir_name,))
            Path(dir_name).mkdir(parents=True, exist_ok=True)
        
        if os.path.isfile(output_filename_MXLANA):
            print("[%s] File %s already exists. Skip." % (time_now_str, output_filename_MXLANA))

        else:
             
            print("[%s] Now compute the mixed-layer integrated quantities. Method = %s" % (time_now_str, args.MLD_method))
            
            if args.MLD_method == "RHO":
                fixed_MLD = -1.0

            elif args.MLD_method == "FIXED500m":
                fixed_MLD = 500.0

            postprocess_ECCO_tools.processECCO(
                self.t,
                output_filename_MXLANA,
                fixed_MLD=fixed_MLD, 
            ) 
 
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
    if new_d.month == 4 and new_d.day != 1:
        continue
 
    if new_d.month == 9 and new_d.day != 30:
        continue

    jobs.append(JOB(new_d))


print("Total jobs: %d" % (len(jobs),))

with Pool(processes=args.nproc) as pool:

    result = pool.map(wrap_retrieve, jobs)



