import numpy as np
import earth_constants as ec
import vorticity_tools

deg2rad = np.pi / 180.0


lat = np.linspace(-90, 90, 181)
lon = np.linspace(-180, 180, 361)

llat = lat[:, np.newaxis]
llon = lon[np.newaxis, :]

llat_r = llat * deg2rad
llon_r = llon * deg2rad

samples = []

empty = np.zeros((len(lat), len(lon)))

samples.append({
    'u': np.cos( llat_r ) * np.sin( llon_r ),
    'v': empty,
    'vort_anly': np.sin(llat_r) * np.sin(llon_r) / ec.r_E,
})
samples.append({
    'u': empty,
    'v': np.cos( llat_r ) * np.sin( llon_r ),
    'vort_anly': np.cos(llat_r) * np.cos(llon_r) / ( ec.r_E * np.cos(llat_r) ),
})


for sample in samples:
    sample['vort_comp'] = vorticity_tools.uv2vort(sample['u'], sample['v'], lat, lon) 




print("Importing Matplotlib...")
import matplotlib.pyplot as plt
print("Done.")
levs = np.linspace(-1, 1, 21) * 2e-7
cmap = "bwr"

fig, ax = plt.subplots(len(samples), 2)

for i, sample in enumerate(samples):
    mappable_1 = ax[i, 0].contourf(lon, lat, sample['vort_anly'], levs, cmap=cmap)
    mappable_2 = ax[i, 1].contourf(lon, lat, sample['vort_comp'], levs, cmap=cmap)

    ax[i, 0].set_title("Sample %d : analytical solution" % (i,))
    ax[i, 1].set_title("Sample %d : computed solution" % (i,))


    plt.colorbar(mappable_1, ax=ax[i, 0])
    plt.colorbar(mappable_2, ax=ax[i, 1])

    for _ax in ax[i, :]:
        _ax.set_xlabel("Longitude")
        _ax.set_ylabel("Latitude")


plt.show()

