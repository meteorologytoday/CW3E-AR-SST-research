import xarray as xr
import numpy as np
import pandas as pd

input_file = "nino34.long.data"
output_file = "NINO34.nc"


raw_data = np.loadtxt(input_file, comments='#')


data = raw_data[:, 1:13].reshape(-1, order='C')
year = raw_data[:, 0]

t = pd.date_range('%04d-01-01' % (year[0],), '%04d-01-01' % (year[-1]+1,), freq='MS', inclusive='left')
data[data == -99.99] = np.nan


data_clim_man = np.nanmean(data.reshape((-1, 12), order='C'), axis=0, keepdims=True)
data_anom_man = data.reshape((-1, 12), order='C') - data_clim_man


data = xr.DataArray(
    data=data,
    dims=["time"],
    coords=dict(
        time=t,
    ),
)


# Making climatology
data_clim = data.groupby("time.month").mean(dim="time").rename("clim")
data_anom = (data.groupby("time.month") - data_clim).rename("anom")
data = data.rename("total")


output_ds = xr.merge([data_clim, data_anom, data])

print("Making output: ", output_file)
output_ds.to_netcdf(output_file)


print("Loading matplotlib...")
import matplotlib.pyplot as plt
print("Done")

fig, ax = plt.subplots(3, 1)

ax[0].plot(data.coords['time'], data)

ax[1].plot(data_clim.coords['month'], data_clim, "k-", label="xarray method")
ax[1].plot(data_clim.coords['month'], data_clim_man.flatten(), "r--", label="traditional")

ax[2].plot(data_anom.coords['time'], data_anom, "k-", label="xarray method")
ax[2].plot(data_anom.coords['time'], data_anom_man.flatten(), "r--", label="traditional")

ax[1].legend()
ax[2].legend()

ax[0].set_title("Raw Data")
ax[1].set_title("Climatology")
ax[2].set_title("Anomalies")

plt.show()
