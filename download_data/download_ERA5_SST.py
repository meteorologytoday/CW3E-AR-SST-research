with open("shared_header.py", "rb") as source_file:
    code = compile(source_file.read(), "shared_header.py", "exec")
exec(code)

import cdsapi

c = cdsapi.Client()

download_dir = "data/ERA5/SST"
file_prefix = "ERA5_SST"

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

            c.retrieve(
                'reanalysis-era5-single-levels',
                {
                    'product_type': 'reanalysis',
                    'format': 'netcdf',
                    'variable': 'sea_surface_temperature',
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

def wrap_retrieve(job):

    job.work()

#for y in range(year_rng[0], year_rng[1]):
#    for m in range(1, 13):
#        jobs.append((y, m))

jobs = []
for d in range(total_days):
    new_d =  beg_time + datetime.timedelta(days=d)

    if 4 <= new_d.month and new_d.month <= 10 :
        continue
    
    jobs.append(JOB(new_d))


print("Total jobs: %d" % (len(jobs),))

print("Create dir: %s" % (download_dir,))
Path(download_dir).mkdir(parents=True, exist_ok=True)


with Pool(processes=2) as pool:

    result = pool.map(wrap_retrieve, jobs)


