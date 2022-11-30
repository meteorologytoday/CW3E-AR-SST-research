with open("shared_header.py", "rb") as source_file:
    code = compile(source_file.read(), "shared_header.py", "exec")
exec(code)

with open("shared_header_login_OSTIA.py", "rb") as source_file:
    code = compile(source_file.read(), "shared_header_login_OSTIA.py", "exec")
exec(code)


download_dir = "data/OSTIA/SST"

total_days = (end_time - beg_time).days

file_prefix = "OSTIA_SST"

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

            command = "python -m motuclient --motu https://nrt.cmems-du.eu/motu-web/Motu --service-id SST_GLO_SST_L4_NRT_OBSERVATIONS_010_001-TDS --product-id METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2 --longitude-min 100 --longitude-max 260 --latitude-min 0 --latitude-max 65 --date-min '%s 00:00:00' --date-max '%s 23:59:59' --variable analysed_sst --variable analysis_error --variable mask --out-dir '%s' --out-name '%s' --user '%s' --pwd '%s'" % (time_now_str,  time_now_str, download_dir, os.path.basename(tmp_filename), username, password)

            
            pleaseRun(command)
            pleaseRun("ncks -O --mk_rec_dmn time %s %s" % (tmp_filename, filename,))
            #os.system("ncra -O %s %s" % (tmp_filename, filename,))
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


with Pool(processes=1) as pool:

    result = pool.map(wrap_retrieve, jobs)


