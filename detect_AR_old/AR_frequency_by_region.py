import numpy as np
import AR_analysis_tools
import watertime_tools

import argparse
from datetime import (timedelta, datetime, timezone)

parser = argparse.ArgumentParser(
                    prog = 'plot_skill',
                    description = 'Plot prediction skill of GFS on AR.',
)

parser.add_argument('--input', type=str, help='Input file', required=True)
parser.add_argument('--output', type=str, help='Output file', default="")
parser.add_argument('--title', type=str, help='Title of figure', default="")
parser.add_argument('--IVT-threshold', type=float, help='Threshold of IVT to determin AR condition.', default=250.0)
parser.add_argument('--output-database', type=str, help='CSV file.', default="")
parser.add_argument('--no-display', action="store_true")
args = parser.parse_args()
print(args)



# Construct histogram, we select no-leap all the time
waterday_cnt = np.zeros((365,))

df = AR_analysis_tools.constructCases(args.input, args.IVT_threshold)

for index, row in df.iterrows():

    datetime_beg = row['datetime_beg']
    datetime_end = row['datetime_end']

    duration_days = int((datetime_end - datetime_beg).total_seconds() / 86400)
    for i in range(duration_days):
        
        datetime_now = datetime_beg + timedelta(days=i)

        wd = watertime_tools.getWaterday(datetime_now, no_leap=True)

        if np.isnan(wd):
            print("Feb 29 encountered. Skip")
            continue


        waterday_cnt[int(np.floor(wd))] += 1


total_wateryears = df.attrs["wateryear_end"] - df.attrs["wateryear_beg"] + 1
waterday_ratio = waterday_cnt / total_wateryears


print("Begin wateryear: %d" % (df.attrs["wateryear_beg"], ))
print("End   wateryear: %d" % (df.attrs["wateryear_end"], ))


# Plot data
print("Loading Matplotlib...")
import matplotlib as mpl
if args.no_display is False:
    mpl.use('TkAgg')
else:
    mpl.use('Agg')
 
#mpl.rc('font', size=20)
#mpl.rc('axes', labelsize=15)
   
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Rectangle
import matplotlib.transforms as transforms
from matplotlib.dates import DateFormatter
from scipy.stats import linregress

print("done")


fig, ax = plt.subplots(1, 1)

if args.title != "":
    fig.suptitle(args.title)
else:
    fig.suptitle(args.input)

rects1 = ax.bar(range(len(waterday_cnt)), waterday_ratio, color='k', width=1)



plt.show()


