with open("shared_header.py", "rb") as source_file:
    code = compile(source_file.read(), "shared_header.py", "exec")
exec(code)

with open("shared_header_login_CFSR.py", "rb") as source_file:
    code = compile(source_file.read(), "shared_header_login_CFSR.py", "exec")
exec(code)


import requests
import shutil

download_dir = "data/CFSR/SST"
file_prefix = "CFSR_SST"

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
    
    

total_months = (end_time.year*12 + end_time.month-1) - (beg_time.year*12 + end_time.month-1)

print("Going to download %d months of data." % (total_months,))

class JOB:

    def __init__(self, t, download_dir, file_prefix, url_fmt):
        self.t = t
        self.download_dir = download_dir
        self.file_prefix = file_prefix
        self.url_fmt = url_fmt


    def work(self):
        
        y = self.t.year
        m = self.t.month

        time_now_str = "%04d-%02d" % (y, m)
        tmp_grb2_filename = "%s/%s_%s.grb2.tmp" % (self.download_dir, self.file_prefix, time_now_str)
        tmp_nc_filename = "%s/%s_%s.nc.tmp" % (self.download_dir, self.file_prefix, time_now_str)
        tmp2_nc_filename = "%s/%s_%s.nc.tmp2" % (self.download_dir, self.file_prefix, time_now_str)
        final_filename = "%s/%s_%s.nc" % (self.download_dir, self.file_prefix, time_now_str)

        already_exists = os.path.isfile(final_filename)

        if already_exists:

            print("[%s] Data already exists. Skip." % (time_now_str, ))


        else:

            print("[%s] Now generate file: %s" % (time_now_str, final_filename,))

            
        url = self.url_fmt % ( y, y, m)
        
        download(url, tmp_grb2_filename)
        pleaseRun("cdo -f nc copy %s %s" % (tmp_grb2_filename, tmp_nc_filename,))
        pleaseRun("ncra --mro -d time,,,24,24 -O %s %s" % (tmp_nc_filename, tmp2_nc_filename,))
        pleaseRun("ncwa -a depth -O %s %s" % (tmp2_nc_filename, final_filename,))

        os.remove(tmp_grb2_filename)
        os.remove(tmp_nc_filename)
        os.remove(tmp2_nc_filename)

def wrap_retrieve(job):
    job.work()


jobs = []
for m in range(total_months):

    _m = beg_time.year * 12 + (beg_time.month-1) + m 
    _t = datetime.datetime(_m // 12, _m % 12 + 1, 1)

    print("%d : " % (m,), _t)

    if _t.year >= 1979 and _t.year <= 2010:
        url_fmt = 'https://rda.ucar.edu/data/ds093.1/%04d/ocnsst.gdas.%04d%02d.grb2'
    elif _t.year >= 2011:
        url_fmt = 'https://rda.ucar.edu/data/ds094.1/%04d/ocnsst.cdas1.%04d%02d.grb2'
    else:
        raise Exception("Time frame is out of range. CFSR only provides year 1979 afterwards.")

    jobs.append(JOB(_t, download_dir, file_prefix, url_fmt))


print("Total jobs: %d" % (len(jobs),))

print("Create dir: %s" % (download_dir,))
Path(download_dir).mkdir(parents=True, exist_ok=True)


with Pool(processes=2) as pool:
    result = pool.map(wrap_retrieve, jobs)







