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

waom4 = xr.open_dataset(output_path + 'waom4_mahti/output_6yr_fabio_inifile/ocean_avg_0001.nc')
waom4_grd = xr.open_dataset(output_path + 'waom4_frc/waom4_grd.nc')

waom2 = xr.open_dataset(output_path + 'waom2_mahti/output_13-15mm_fabio_inifile/ocean_avg_0001.nc')
waom2_grd = xr.open_dataset(output_path + 'waom2_frc/waom2_grd.nc')

temp_bot_4km = np.mean(waom4["temp"].isel(s_rho=0), axis=0)
temp_bot_2km = np.mean(waom2["temp"].isel(s_rho=0), axis=0)

print(temp_bot_4km.shape, temp_bot_2km.shape)

# read grid file for lon/lat coordinates
lat_rho_4km = waom4_grd.variables["lat_rho"]
lon_rho_4km = waom4_grd.variables["lon_rho"]
waom4.coords['lat_rho']=lat_rho_4km#.transpose() # put lat_rho into ds dataset
waom4.coords['lon_rho']=lon_rho_4km#.transpose() # put lon_rho into ds dataset

lat_rho_2km = waom2_grd.variables["lat_rho"]
lon_rho_2km = waom2_grd.variables["lon_rho"]
waom2.coords['lat_rho']=lat_rho_2km
waom2.coords['lon_rho']=lon_rho_2km

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
fig = plt.figure(figsize=(20,10))

ax1 = fig.add_subplot(121, projection=proj)
ct=plt.pcolor(lon_rho_4km,lat_rho_4km,temp_bot_4km, transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=tmin, vmax=tmax)
cbart =fig.colorbar(ct, extend='both')
plt.title('Bottom temperature - WAOM 4km')
ax1.gridlines() # draw_labels=True,linewidth=2, color='gray', alpha=0.5, linestyle='--')
ax1.coastlines(resolution='110m')
lonlat_labels(ax1)
ax1.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
ax1.add_feature(cfeature.LAND, zorder=4, edgecolor='black', facecolor='white')

ax2 = fig.add_subplot(122, projection=proj)
ct=plt.pcolor(lon_rho_2km,lat_rho_2km,temp_bot_2km, transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=tmin, vmax=tmax)
cbart =fig.colorbar(ct, extend='both')
plt.title('Bottom temperature - WAOM 4km')
ax2.gridlines() # draw_labels=True,linewidth=2, color='gray', alpha=0.5, linestyle='--')
ax2.coastlines(resolution='110m')
lonlat_labels(ax1)
ax2.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
ax2.add_feature(cfeature.LAND, zorder=4, edgecolor='black', facecolor='white')

plt.savefig('/users/boeiradi/COLD_project/postprocessing/figs/WAOM_comparison/waom4xwaom2_bottom_temp_jan.png')
plt.close()
