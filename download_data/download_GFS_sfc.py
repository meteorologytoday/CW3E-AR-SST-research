with open("shared_header.py", "rb") as source_file:
    code = compile(source_file.read(), "shared_header.py", "exec")
exec(code)

with open("shared_header_login_CFSR.py", "rb") as source_file:
    code = compile(source_file.read(), "shared_header_login_CFSR.py", "exec")
exec(code)


import requests
import shutil

download_dir = "data/GFS/fcst_sfc"
file_prefix = "GFS_0p25"
fcst_hrs = [0, 240]

fcst_hrs_avg = {
    '0'   : [6, 12, 18, 24],   # the first 24 hrs (0-24hr) needs four files fcst_[006, 012, 018, 024] to compute the average fluxes
    '240' : [252, 264],        # the 10-th day (240-264hr) needs two  files fcst_[252, 264] to compute the average fluxes
}


login_post = {
    'email'  : CFSR_email,
    'passwd' : CFSR_passwd,
    'action' : 'login',
}

shared_session = requests.Session()

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

            okay = True

        except Exception as e:
            print("Someting went wrong: ", e)
            print("Try to redownload again. Current attempt: %d" % (attempt,)) 
            continue

        if okay:
            print("Download success")
            break
    

    return okay

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

        time_now_str = "%04d%02d%02d_f%03d" % (y, m, d, self.fcst_hr)

        tmp_nc_filename = "%s/%s_%s.nc.tmp" % (self.download_dir, self.file_prefix, time_now_str)
        tmp_grb2_filenames = []
        tmp_grb2_filename_ave = "%s/%s_%s_ave.grb2.tmp" % (self.download_dir, self.file_prefix, time_now_str)

        final_filename  = "%s/%s_%s.sfcflx.nc" % (self.download_dir, self.file_prefix, time_now_str)

        already_exists = os.path.isfile(final_filename)

        if already_exists:

            print("[%s] Data already exists. Skip." % (time_now_str, ))


        else:

            print("[%s] Now generate file: %s" % (time_now_str, final_filename_,))

            datetime_str = self.t.strftime("%Y%m%d")
            needed_fcst_hrs = fcst_hrs_avg['%d' % self.fcst_hr]
            urls = [
                'https://rda.ucar.edu/data/ds084.1/%04d/%s/gfs.0p25.%s00.f%03d.grib2' % (y, datetime_str, datetime_str, needed_fcst_hr) for _fcst_hr in needed_fcst_hrs
            ]

            okay = True
            needed_fcst_hrs = fcst_hrs_avg['%d' % self.fcst_hr]
            for hr in needed_fcst_hrs:

                tmp_grb2_filename = "%s/%s_%s_%02d.grb2.tmp" % (self.download_dir, self.file_prefix, time_now_str, hr)
                url = 'https://rda.ucar.edu/data/ds084.1/%04d/%s/gfs.0p25.%s%02d.f%03d.grib2' % (y, datetime_str, datetime_str, hr, self.fcst_hr)
                okay = okay and download(url, tmp_grb2_filename)
                
                tmp_grb2_filenames.append(tmp_grb2_filename)

            if okay:

                interval_hrs = needed_fcst_hrs[1] - needed_fcst_hrs[0]

                pleaseRun("gmerge - %s | wgrib2 -  -match ':(LHTFL|SHTFL|PRATE|UFLX|VFLX):surface:' -ave %dhr %s" % (interval_hrs, " ".join(tmp_grb2_filenames), tmp_grb2_filename_ave))
                pleaseRun("wgrib2 %s -netcdf %s" % (tmp_grb2_filename_ave, tmp_nc_filename,))
                pleaseRun("ncks -O -d lat,0.0,65.0 %s %s" % (tmp_nc_filename, final_filename))


            else:

                print("Not okay. Skip.")

            for fn in [tmp_grb2_filename_ave, tmp_nc_filename] + tmp_grb2_filenames:
                try:
                    os.remove(fn)
                except OSError:
                    pass

def wrap_retrieve(job):
    job.work()


jobs = []

    
for cnt in range(total_days):

    new_t =  beg_time + datetime.timedelta(days=cnt)
    timestr = new_t.strftime("%Y%m%d")

    print("cnt=%d, timestr=%s" % (cnt, timestr,))

    if new_t.month <=10 and new_t.month >=5:
        continue

    for fcst_hr in fcst_hrs:

        jobs.append(JOB(new_t, fcst_hr, download_dir, file_prefix))



print("Total jobs: %d" % (len(jobs),))

print("Create dir: %s" % (download_dir,))
Path(download_dir).mkdir(parents=True, exist_ok=True)


with Pool(processes=1) as pool:
    result = pool.map(wrap_retrieve, jobs)







