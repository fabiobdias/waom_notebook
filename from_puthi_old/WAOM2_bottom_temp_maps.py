import xarray as xr
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Agg')

from netCDF4 import Dataset
from netCDF4 import num2date, date2num

output_path = '/scratch/project_2000339/boeiradi/'

waom2 = xr.open_dataset(output_path + 'waom2_mahti/output_13-15mm_fabio_inifile/ocean_avg_0001.nc')
waom2_grd = xr.open_dataset(output_path + 'waom2_frc/waom2_grd.nc')

temp_bot_2km = np.mean(waom2["temp"].isel(s_rho=0), axis=0)

# read grid file for lon/lat coordinates
lat_rho_2km = waom2_grd.variables["lat_rho"]
lon_rho_2km = waom2_grd.variables["lon_rho"]
waom2.coords['lat_rho']=lat_rho_2km
waom2.coords['lon_rho']=lon_rho_2km

print(temp_bot_2km.shape)

def lonlat_labels(ax):
    # latitude labels
    ax.text(120,-80,'80$^{\circ}$S',transform=ccrs.PlateCarree(),color='gray')
    ax.text(120,-70,'70$^{\circ}$S',transform=ccrs.PlateCarree(),color='gray')
    # longitude labels
    ax.text(0,-66,'0$^{\circ}$',transform=ccrs.PlateCarree(),color='gray')
    ax.text(-60,-48,'60$^{\circ}$W',transform=ccrs.PlateCarree(),color='gray')
    ax.text(-120,-48,'120$^{\circ}$W',transform=ccrs.PlateCarree(),color='gray')
    ax.text(180,-60,'180$^{\circ}$',transform=ccrs.PlateCarree(),color='gray')
    return

tmin=-2.5
tmax=2.6
levelsT = np.arange(-2.5,2.5,.05)

# call cartopy projection
proj = ccrs.SouthPolarStereo()
fig = plt.figure(figsize=(10,10))

ax2 = fig.add_subplot(111, projection=proj)
ct=plt.pcolor(lon_rho_2km,lat_rho_2km,temp_bot_2km, transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=tmin, vmax=tmax)
cbart =fig.colorbar(ct, extend='both')
plt.title('Bottom temperature - WAOM 2km')
ax2.gridlines() # draw_labels=True,linewidth=2, color='gray', alpha=0.5, linestyle='--')
ax2.coastlines(resolution='110m')
lonlat_labels(ax2)
ax2.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
ax2.add_feature(cfeature.LAND, zorder=4, edgecolor='black', facecolor='white')

plt.savefig('/users/boeiradi/COLD_project/postprocessing/figs/WAOM_comparison/waom2_bottom_temp_jan.png')
plt.close()
