import numpy as np
import load_data
import netCDF4
from datetime import (datetime, timedelta)

from pathlib import Path
import argparse

parser = argparse.ArgumentParser(
                    prog = 'plot_skill',
                    description = 'Plot prediction skill of GFS on AR.',
)

parser.add_argument('--AR-npz', type=str, help='Save the timeseries into an npz file', required=True)
parser.add_argument('--skill-npz', type=str, help='Save the timeseries into an npz file', required=True)
parser.add_argument('--output', type=str, help='Output directory.', default="")
parser.add_argument('--no-display', action="store_true")

args = parser.parse_args()

print(args)

AR_data    = np.load(args.AR_npz)
sk_data = np.load(args.skill_npz)

# test if times are the same
t_AR = AR_data['time']
t_sk = sk_data['time']
if len(t_AR) != len(t_sk):
    raise Exception("Time vectors have different lengths.")

for i in range(len(t_AR)):
    if t_AR[i] != t_sk[i]:
        raise Exception("Time vectors have different values.")

# Plot data
print("Loading Matplotlib...")
import matplotlib
if args.no_display is False:
    matplotlib.use('TkAgg')
else:
    matplotlib.use('Agg')
    
import matplotlib.pyplot as plt
from matplotlib import cm
import cartopy.crs as ccrs
print("done")

#AR_vars = ["IWV", "IVT", "IWVKE", "HGT_850mb"]
AR_vars = ["IWV", "IVT", "HGT_850mb", "HGT_500mb"]

fig, ax = plt.subplots(len(AR_vars), 1, figsize=(6, 4*len(AR_vars)))

for i, AR_var in enumerate(AR_vars):

    print("Plot %s" % (AR_var,))
    _ax = ax[i]
    _axt = _ax.twinx()

    if AR_var not in ["IWV", "IVT", "IWVKE", ]:
        twin_AR_var = "IVT"
    else:
        twin_AR_var = AR_var



    _axt.plot(t_AR, AR_data[twin_AR_var], "k--")
    _ax.plot(t_AR, np.sqrt(sk_data[AR_var]), "r-")

  
    _ax.set_title("Forecast error of %s" % (AR_var, ))


if not args.no_display:
    plt.show()

if args.output != "":
    
    fig.savefig(args.output, dpi=200)


