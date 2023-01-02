with open("shared_header.py", "rb") as source_file:
    code = compile(source_file.read(), "shared_header.py", "exec")
exec(code)

# Reset beg_time and end_time for testing
#beg_time = datetime.datetime(1990,   10,  1)
#end_time = datetime.datetime(1990,   12, 31)

import cdsapi
import numpy as np
import postprocess_ORA5

c = cdsapi.Client()

download_dir = "data/ORA5/raw"
processed_dir = "data/ORA5/processed"


def ORA5_filename(varname, dim, y, m):

    dimension_str = "%dD" % (dim,)
    return "%s_control_monthly_highres_%s_%04d%02d_CONS_v0.1.nc" % (varname, dimension_str, y, m)

def datetime2months(t):
    return t.year*12 + (t.month-1)

def months2datetime(m):

    year = m // 12
    mon  = m % 12 + 1

    return datetime.datetime(year, mon, 1)


total_years = end_time.year - beg_time.year + 1

total_months = datetime2months(end_time) - datetime2months(beg_time) + 1

print("Going to download %d months of data." % (total_months,))


class JOB:

    def __init__(self, t):
        self.t = t


    def work(self):
        
        y = self.t.year
        m = self.t.month

        time_now_str = "%04d-%02d" % (y, m)

        filenames = [ 
            "%s/ORA5_votemper_%s.nc" % (download_dir, time_now_str),
            "%s/ORA5_vosaline_%s.nc" % (download_dir, time_now_str),
            "%s/ORA5_somxl010_%s.nc" % (download_dir, time_now_str),
            "%s/ORA5_somxl030_%s.nc" % (download_dir, time_now_str),
        ]
        
        tmp_zip_filename_3D  = "%s/ORA5_3D_%s.zip" % (download_dir, time_now_str)
        tmp_zip_filename_2D = "%s/ORA5_2D_%s.zip" % (download_dir, time_now_str)
        tmp_unzip_folder = "%s/tmp_unzip_%s" % (download_dir, time_now_str)


 
        
        all_files_exist = np.all([ os.path.isfile(filename) for filename in filenames])
        if all_files_exist:

            print("[%s] Data already exists. Skip." % (time_now_str, ))

        else:

            print("[%s] Now download data to produces files:" % (time_now_str,))
            for k, filename in enumerate(filenames):
                print("   (%02d) %s" % (k+1, filename,))
            
            c.retrieve(
                'reanalysis-oras5',
                {
                    'format': 'zip',
                    'product_type': 'consolidated',
                    'variable': [
                        'potential_temperature', 'salinity', 
                    ],
                    'vertical_resolution': 'all_levels',
                    'month': [
                        "%02d" % (m,),
                    ],
                    'year': "%04d" % (y,)
                },
                tmp_zip_filename_3D,
            )

            c.retrieve(
                'reanalysis-oras5',
                {
                    'format': 'zip',
                    'product_type': 'consolidated',
                    'variable': [
                        'mixed_layer_depth_0_01', 'mixed_layer_depth_0_03',
                    ],
                    'vertical_resolution': 'single_level',
                    'month': [
                        "%02d" % (m,),
                    ],
                    'year': "%04d" % (y,)
                },
                tmp_zip_filename_2D,
            )
            
            pleaseRun("unzip -d %s -o %s" % (tmp_unzip_folder, tmp_zip_filename_3D,))
            pleaseRun("unzip -d %s -o %s" % (tmp_unzip_folder, tmp_zip_filename_2D,))
            
            # rename files
            varinfos = [
                ["somxl010", 2],
                ["somxl030", 2],
                ["votemper", 3],
                ["vosaline", 3],
            ]

            for varname, dim in varinfos:
                old_filename = ORA5_filename(varname, dim, y, m)
                new_filename = "ORA5_%s_%s.nc" % (varname, time_now_str)

                pleaseRun("mv %s/%s %s/%s" % (tmp_unzip_folder, old_filename, download_dir, new_filename))


            pleaseRun("rm -rf %s" % (tmp_unzip_folder,))
            pleaseRun("rm -f %s %s" % (tmp_zip_filename_2D, tmp_zip_filename_3D, ))    
 

        MLD_filenames = {
            "somxl010" : "%s/ORA5_%s_%s.nc" % (download_dir, "somxl010", time_now_str),
            "somxl030" : "%s/ORA5_%s_%s.nc" % (download_dir, "somxl030", time_now_str),
        }
            
        for varname_MLD, MLD_filename in MLD_filenames.items():

                      
            processed_filename = "%s/ORA5_NillerKrausMixedLayerDynamics_%s_%s.nc" % (processed_dir, varname_MLD, time_now_str)
            if os.path.isfile(processed_filename):

                print("File %s exists. " % (processed_filename,))

            else:
            
                print("File %s does not exist. Now generating it." % (processed_filename,))
                
                T_filename   = "%s/ORA5_%s_%s.nc" % (download_dir, "votemper", time_now_str)
                S_filename   = "%s/ORA5_%s_%s.nc" % (download_dir, "vosaline", time_now_str)
                output_filename = processed_filename
                postprocess_ORA5.processORA5ForNillerKrausMixedLayerDynamics(
                    MLD_filename,
                    T_filename,
                    S_filename,
                    output_filename,
                    varname_MLD = varname_MLD,
                )
                
        
def wrap_retrieve(job):

    try:
        job.work()

    except Exception as e:
        print(e)
        print("Something wrong with job : ", job.t)


jobs = []
for m in range(total_months):

    new_datetime =  months2datetime(datetime2months(beg_time) + m)

    jobs.append(JOB(new_datetime))


print("Total jobs: %d" % (len(jobs),))

print("Create dir: %s" % (download_dir,))
Path(download_dir).mkdir(parents=True, exist_ok=True)

print("Create dir: %s" % (processed_dir,))
Path(processed_dir).mkdir(parents=True, exist_ok=True)




with Pool(processes=4) as pool:

    result = pool.map(wrap_retrieve, jobs)

print("Done.")
