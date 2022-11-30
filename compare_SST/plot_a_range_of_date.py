from datetime import (datetime, timedelta)
import os

def pleaseRun(cmd):
    print(">> %s" % cmd)
    os.system(cmd)



beg_date = datetime(2018, 1, 1)
end_date = datetime(2018, 3, 5)

total_days = (end_date - beg_date).days

print("Begin date: ", beg_date)
print("End   date: ", end_date)
print("Total days: ", total_days)

for d in range(total_days):
    new_d =  beg_date + timedelta(days=d)
    pleaseRun("python3 plot_SST.py --no-display --output-dir output --date %s" % (new_d.strftime("%Y-%m-%d"), ) )

