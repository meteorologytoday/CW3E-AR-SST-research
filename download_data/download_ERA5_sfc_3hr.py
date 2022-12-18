with open("shared_header.py", "rb") as source_file:
    code = compile(source_file.read(), "shared_header.py", "exec")
exec(code)

import numpy as np
import cdsapi

c = cdsapi.Client()



dhr = 3
if 24 % dhr != 0:
    raise Exception("Not cool. 24 / dhr is not an integer.")

subcycles = int(24 / dhr)


download_dir = "data/ERA5/sfc_%dhr" % (dhr,)
file_prefix = "ERA5_sfc"

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
        tmp_filename = "%s/%s_%s.nc.tmp" % (download_dir, file_prefix, time_now_str)

        filenames = [ 
            "%s/%s_%s_%02d.nc" % (download_dir, file_prefix, time_now_str, dhr * i) for i in range(subcycles)
        ]
 
        existence = [os.path.isfile(filename) for filename in filenames]

        if np.all(existence):

            print("[%s] Data already exists. Skip." % (time_now_str, ))


        else:

            print("[%s] Now producing files: %s" % (time_now_str, filenames,))

            try:
                c.retrieve(
                    'reanalysis-era5-single-levels',
                    {
                        'product_type': 'reanalysis',
                        'format': 'netcdf',
                        'variable': [
                            '10m_u_component_of_wind', '10m_v_component_of_wind', '2m_temperature',
                            'mean_evaporation_rate', 'mean_surface_latent_heat_flux', 'mean_surface_net_long_wave_radiation_flux',
                            'mean_surface_net_short_wave_radiation_flux', 'mean_surface_sensible_heat_flux', 'mean_total_precipitation_rate',
                            'mean_vertically_integrated_moisture_divergence', 'sea_surface_temperature',
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

                for i, filename in enumerate(filenames):
                    pleaseRun("ncra -O -d time,%d,%d %s %s" % (dhr*i, dhr*(i+1)-1, tmp_filename, filename,))
                
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

    if 5 <= new_d.month and new_d.month <= 10 :
        continue
    
    jobs.append(JOB(new_d))


print("Total jobs: %d" % (len(jobs),))

print("Create dir: %s" % (download_dir,))
Path(download_dir).mkdir(parents=True, exist_ok=True)


with Pool(processes=1) as pool:

    result = pool.map(wrap_retrieve, jobs)


