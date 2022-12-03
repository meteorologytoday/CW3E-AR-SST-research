with open("shared_header.py", "rb") as source_file:
    code = compile(source_file.read(), "shared_header.py", "exec")
exec(code)

with open("shared_header_login_CFSR.py", "rb") as source_file:
    code = compile(source_file.read(), "shared_header_login_CFSR.py", "exec")
exec(code)


import requests
import shutil
from convert_GFS_output import convertGFSOutput
from compute_AR import computeAR

download_dir = "data/GFS/fcst"
file_prefix = "GFS_0p25"
fcst_hrs = [0, 120, 240]
login_post = {
    'email'  : CFSR_email,
    'passwd' : CFSR_passwd,
    'action' : 'login',
}

shared_session = requests.Session()
#shared_session.auth = (usr, pwd)

def download(url, output, max_attempt=5):


    needs_login = True
    okay = False

    for attempt in range(max_attempt):

        if needs_login:

            print("(Re)Login required. Login now.")
            login_resp = shared_session.post('https://rda.ucar.edu/cgi-bin/login', data=login_post)

        try:

            print("Downloading url: %s to %s" % (url, output,))
            with shared_session.get(url, stream=True) as download_resp:
                
                if download_resp.status_code != 200:
                    print("Request data failed. Login required.")
                    needs_login = True
                    continue
                else:
                    needs_login = False

                with open(output, 'wb') as f:
                    shutil.copyfileobj(download_resp.raw, f)

        except Exception as e:
            print("Someting went wrong: ", e)
            print("Try to redownload again. Current attempt: %d" % (attempt,)) 
            continue

        print("Download success")
        break
    

total_days = (end_time - beg_time).days
print("Going to download %d days of data." % (total_days,))

class JOB:

    def __init__(self, t, fcst_hr, download_dir, file_prefix):
        self.t = t
        self.download_dir = download_dir
        self.file_prefix = file_prefix
        self.fcst_hr = fcst_hr


    def work(self):
        
        y = self.t.year
        m = self.t.month
        d = self.t.day

        time_now_str = "%04d-%02d-%02d_f%03d" % (y, m, d, self.fcst_hr)
        tmp_grb2_filename = "%s/%s_%s.grb2.tmp" % (self.download_dir, self.file_prefix, time_now_str)
        tmp_nc_filename = "%s/%s_%s.nc.tmp" % (self.download_dir, self.file_prefix, time_now_str)
        final_filename_AR  = "%s/%s_%s.AR.nc" % (self.download_dir, self.file_prefix, time_now_str)
        final_filename_HGT500mb = "%s/%s_%s.HGT_500mb.nc" % (self.download_dir, self.file_prefix, time_now_str)
        final_filename_HGT850mb = "%s/%s_%s.HGT_850mb.nc" % (self.download_dir, self.file_prefix, time_now_str)

        already_exists = os.path.isfile(final_filename_AR) and os.path.isfile(final_filename_HGT500mb) and os.path.isfile(final_filename_HGT850mb)

        if already_exists:

            print("[%s] Data already exists. Skip." % (time_now_str, ))


        else:

            print("[%s] Now generate file: %s" % (time_now_str, ", ".join([final_filename_AR, final_filename_HGT500mb, final_filename_HGT850mb],)))

            datetime_str = self.t.strftime("%Y%m%d")
            url = 'https://rda.ucar.edu/data/ds084.1/%04d/%s/gfs.0p25.%s00.f%03d.grib2' % (y, datetime_str, datetime_str, self.fcst_hr)
            
            download(url, tmp_grb2_filename)
            convertGFSOutput(tmp_grb2_filename, tmp_nc_filename)

            pleaseRun("ncks -O -v HGT -d lev,500.0,500.0 %s %s" % (tmp_nc_filename, final_filename_HGT500mb))
            pleaseRun("ncks -O -v HGT -d lev,850.0,850.0 %s %s" % (tmp_nc_filename, final_filename_HGT850mb))

            computeAR(tmp_nc_filename, final_filename_AR)

            os.remove(tmp_grb2_filename)
            os.remove(tmp_nc_filename)

def wrap_retrieve(job):
    job.work()


jobs = []

    
for d in range(total_days):

    new_d =  beg_time + datetime.timedelta(days=d)
    timestr = new_d.strftime("%Y%m%d")

    print("%d : " % (d,), timestr)

    if new_d.month <=10 and new_d.month >=5:
        continue

    for fcst_hr in fcst_hrs:

        jobs.append(JOB(new_d, fcst_hr, download_dir, file_prefix))



print("Total jobs: %d" % (len(jobs),))

print("Create dir: %s" % (download_dir,))
Path(download_dir).mkdir(parents=True, exist_ok=True)


with Pool(processes=2) as pool:
    result = pool.map(wrap_retrieve, jobs)







