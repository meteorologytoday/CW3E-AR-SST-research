with open("shared_header.py", "rb") as source_file:
    code = compile(source_file.read(), "shared_header.py", "exec")
exec(code)

import postprocess_ERA5_sfc
import cdsapi

c = cdsapi.Client()

download_dir = "data/ERA5/cloud"
file_prefix = "ERA5_cloud"

processed_dir = "data/ERA5/sfc_processed"

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
        filename = "%s/%s_%s.nc" % (download_dir, file_prefix, time_now_str)
        tmp_filename = "%s/%s_%s.nc.tmp" % (download_dir, file_prefix, time_now_str)

        already_exists = os.path.isfile(filename)

        if already_exists:

            print("[%s] Data already exists. Skip." % (time_now_str, ))


        else:

            print("[%s] Now download file: %s" % (time_now_str, filename,))

            try:
                c.retrieve(
                    'reanalysis-era5-single-levels',
                    {
                        'product_type': 'reanalysis',
                        'format': 'netcdf',
                        'variable': [
                            'high_cloud_cover', 'medium_cloud_cover', 'low_cloud_cover', 'total_cloud_cover',
                        ],
                        'day': [
                                "%02d" % d,
                            ],
                        'month': [
                                "%02d" % m,
                            ],
                        'year': [
                                "%04d" % y,
                            ],
                        'time': [
                            '00:00', '01:00', '02:00',
                            '03:00', '04:00', '05:00',
                            '06:00', '07:00', '08:00',
                            '09:00', '10:00', '11:00',
                            '12:00', '13:00', '14:00',
                            '15:00', '16:00', '17:00',
                            '18:00', '19:00', '20:00',
                            '21:00', '22:00', '23:00',
                            ],
                        'area': [
                            65, -180, 0,
                            180,
                        ],
                    },
                tmp_filename)
                pleaseRun("ncks -O --mk_rec_dmn time %s %s" % (tmp_filename, tmp_filename,))
                pleaseRun("ncra -O %s %s" % (tmp_filename, filename,))
                os.remove(tmp_filename)

            except Exception as e:

                print("Something goes wrong when generating %s" % (filename,))
                print(str(e))
                

 
def wrap_retrieve(job):

    job.work()

#for y in range(year_rng[0], year_rng[1]):
#    for m in range(1, 13):
#        jobs.append((y, m))

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

print("Create dir: %s" % (download_dir,))
Path(download_dir).mkdir(parents=True, exist_ok=True)

print("Create dir: %s" % (processed_dir,))
Path(processed_dir).mkdir(parents=True, exist_ok=True)



with Pool(processes=16) as pool:

    result = pool.map(wrap_retrieve, jobs)


