# read nc output from WAOM 4km run

import xarray as xr
import pandas as p
import numpy as np
import numpy.ma as ma
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib.dates as dates
from matplotlib.ticker import NullFormatter
from matplotlib.colors import LinearSegmentedColormap   # for custom colormaps
import matplotlib.path as mpath
import matplotlib as mpl
mpl.use('Agg')

from datetime import datetime, timedelta

from netCDF4 import Dataset
from netCDF4 import num2date, date2num

import iris
import iris.iterate
import iris.coords
import iris.plot as iplt
import gsw

## --------------------- load ROMS avg output
roms_temp = np.empty((0,31,1325,1575))

ds = xr.open_dataset("/scratch/project_2000339/boeiradi/waom4_mahti/output_6yr_fabio_inifile/ocean_avg_00" + mm + ".nc")
temp = ds.variables["temp"]
print('temp size =', temp.shape)
ds = ds.set_coords(['Cs_r', 'Cs_w', 'hc', 'h', 'Vtransform'])

if ds.Vtransform == 1:
    Zo_rho = ds.hc * (ds.s_rho - ds.Cs_r) + ds.Cs_r * ds.h
    z_rho = Zo_rho + ds.zeta * (1 + Zo_rho/ds.h)
    print("Vtransform=1")
elif ds.Vtransform == 2:
    Zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * ds.h) / (ds.hc + ds.h)
    z_rho = ds.zeta + (ds.zeta + ds.h) * Zo_rho
    print("Vtransform=2")
ds.coords['z_rho'] = z_rho.transpose() # put z_rho into ds dataset

# read grid file for lon/lat coordinates
dg = xr.open_dataset("/scratch/project_2000339/boeiradi/waom4_frc/waom4_grd.nc")

lat_rho = dg.variables["lat_rho"]
lon_rho = dg.variables["lon_rho"]
ds.coords['lat_rho']=lat_rho#.transpose() # put lat_rho into ds dataset
ds.coords['lon_rho']=lon_rho#.transpose() # put lon_rho into ds dataset

def lonlat_labels(ax):
    # latitude labels
    ax.text(120,-80,'80$^{\circ}$S',transform=ccrs.PlateCarree(),color='gray')
    ax.text(120,-70,'70$^{\circ}$S',transform=ccrs.PlateCarree(),color='gray')
    # longitude labels
    ax.text(0,-66,'0$^{\circ}$',transform=ccrs.PlateCarree(),color='gray')
    #ax.text(60,-53,'60$^{\circ}$E',transform=ccrs.PlateCarree(),color='gray')
    #ax.text(120,-53,'120$^{\circ}$E',transform=ccrs.PlateCarree(),color='gray')
    ax.text(-60,-48,'60$^{\circ}$W',transform=ccrs.PlateCarree(),color='gray')
    ax.text(-120,-48,'120$^{\circ}$W',transform=ccrs.PlateCarree(),color='gray')
    ax.text(180,-60,'180$^{\circ}$',transform=ccrs.PlateCarree(),color='gray')
    return


# Compute a circle in axes coordinates, which we can use as a boundary
# for the map. We can pan/zoom as much as we like - the boundary will be
# permanently circular.
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)

print('Shape model Temp: ', ds.temp.isel(s_rho=-1).shape)
# time average WAOM4 
#temp_sfc = ds.temp.isel(s_rho=-1).mean(axis=0)
#salt_sfc = ds.salt.isel(s_rho=-1).mean(axis=0)#, ocean_time=1)
#Hsbl = ds.Hsbl.isel(ocean_time=9) #mean(axis=0)
temp_sfc = np.mean(roms_temp[:,-1,:,:], axis=0)
salt_sfc = np.mean(roms_salt[:,-1,:,:], axis=0)
print('temp size=', temp_sfc.shape)
print('lon/lat_rho size =', lon_rho.shape,lat_rho.shape)

print('Shape EN4 Temp: ', EN_temp.shape)
# time average EN4 
EN_temp_sfc = np.mean(EN_temp[:,0,:,:],axis=0)
EN_salt_sfc = np.mean(EN_salt[:,0,:,:],axis=0)

tmin=-2.5
tmax=2.2
smin=34.
smax=35.
levelsT = np.arange(-2.5,2.5,.05)
levelsS = np.arange(34.,35.,.05)

# call cartopy projection
proj = ccrs.SouthPolarStereo()
fig = plt.figure(figsize=(10,10))

ax1 = fig.add_subplot(221, projection=proj)
# sea ice contour
ci=plt.contour(ice_xgrid, ice_ygrid,np.squeeze(ice_mm[8,:,:]), levels=[.75],colors='dimgray', transform=ccrs.Stereographic(**kw))
ci=plt.contour(ice_xgrid, ice_ygrid,np.squeeze(ice_mm[2,:,:]), levels=[.75],colors='white', transform=ccrs.Stereographic(**kw))
ct=plt.pcolor(lon_rho,lat_rho,temp_sfc, transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=tmin, vmax=tmax)
cbart =fig.colorbar(ct, extend='both')
plt.title('Surface temperature - WAOM 4km')
ax1.gridlines() # draw_labels=True,linewidth=2, color='gray', alpha=0.5, linestyle='--')
ax1.coastlines(resolution='110m')
lonlat_labels(ax1)
ax1.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
ax1.add_feature(cfeature.LAND, zorder=4, edgecolor='black', facecolor='white')
    
ax2 = fig.add_subplot(222, projection=proj)
ci=plt.contour(ice_xgrid, ice_ygrid,np.squeeze(ice_mm[8,:,:]), levels=[.75],colors='dimgray', transform=ccrs.Stereographic(**kw)) #
ci=plt.contour(ice_xgrid, ice_ygrid,np.squeeze(ice_mm[2,:,:]), levels=[.75],colors='white', transform=ccrs.Stereographic(**kw)) #
cs=plt.pcolor(lon_rho,lat_rho,salt_sfc, transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=smin, vmax=smax)
cbars =fig.colorbar(cs, extend='both')
plt.title('Surface salinity - WAOM 4km')
ax2.coastlines(resolution='110m')
ax2.gridlines()
lonlat_labels(ax2)
ax2.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
ax2.add_feature(cfeature.LAND, zorder=4, edgecolor='black', facecolor='white')

ax3 = fig.add_subplot(223, projection=proj)
ct2=plt.pcolor(EN_lon[:,:22],EN_lat[:,:22],np.transpose(EN_temp_sfc[:22,:]), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=tmin, vmax=tmax)
cbart =fig.colorbar(ct2, extend='both')
plt.title('Surface temperature - EN4')
ax3.gridlines() # draw_labels=True,linewidth=2, color='gray', alpha=0.5, linestyle='--')
ax3.coastlines(resolution='110m')
lonlat_labels(ax3)
ax3.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
ax3.add_feature(cfeature.LAND, zorder=4, edgecolor='black', facecolor='white')

ax4 = fig.add_subplot(224, projection=proj)
cs2=plt.pcolor(EN_lon[:,:22],EN_lat[:,:22],np.transpose(EN_salt_sfc[:22,:]),transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=smin, vmax=smax)
cbars =fig.colorbar(cs2, extend='both')
plt.title('Surface salinity - EN4')
ax4.gridlines()
ax4.coastlines(resolution='110m')
lonlat_labels(ax4)
ax4.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
ax4.add_feature(cfeature.LAND, zorder=4, edgecolor='black', facecolor='white')

plt.tight_layout()

plt.savefig('figs/WAOM4/waom4xEN4_surface_TS__annual.png')

# 2. Plot bottom temp/salt ROMSxEN4

### plot some maps
    
# time average WAOM4 
#temp_bot = ds.temp.isel(s_rho=0).mean(axis=0)
#salt_bot = ds.salt.isel(s_rho=0).mean(axis=0)#, ocean_time=1)
#Hsbl = ds.Hsbl.isel(ocean_time=9) #mean(axis=0)
temp_bot = np.mean(roms_temp[:,0,:,:], axis=0)
salt_bot = np.mean(roms_salt[:,0,:,:], axis=0)
print('temp size=', temp_bot.shape)

# time average EN4: bottom values are obtained above;
EN_temp_mbot = np.mean(EN_temp_bot,axis=0)
EN_salt_mbot = np.mean(EN_salt_bot,axis=0)

# call cartopy projection
proj = ccrs.SouthPolarStereo()
fig = plt.figure(figsize=(10,10))

ax1 = fig.add_subplot(221, projection=proj)
ct=plt.pcolor(lon_rho,lat_rho,temp_bot, transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=tmin, vmax=tmax)
cbart =fig.colorbar(ct, extend='both')
plt.title('Bottom temperature - WAOM 4km')
ax1.gridlines() # draw_labels=True,linewidth=2, color='gray', alpha=0.5, linestyle='--')
ax1.coastlines(resolution='110m')
lonlat_labels(ax1)
ax1.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
ax1.add_feature(cfeature.LAND, zorder=4, edgecolor='black', facecolor='white')

ax2 = fig.add_subplot(222, projection=proj)
cs=plt.pcolor(lon_rho,lat_rho,salt_bot, transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=smin, vmax=smax)
cbars =fig.colorbar(cs, extend='both')
plt.title('Bottom salinity - WAOM 4km')
ax2.coastlines(resolution='110m')
ax2.gridlines()
lonlat_labels(ax2)
ax2.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
ax2.add_feature(cfeature.LAND, zorder=4, edgecolor='black', facecolor='white')

ax3 = fig.add_subplot(223, projection=proj)
ct2=plt.pcolor(EN_lon[:,:22],EN_lat[:,:22],np.transpose(EN_temp_mbot[:22,:]), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=tmin, vmax=tmax)
cbart =fig.colorbar(ct2, extend='both')
plt.title('Bottom temperature - EN4')
ax3.gridlines() # draw_labels=True,linewidth=2, color='gray', alpha=0.5, linestyle='--')
ax3.coastlines(resolution='110m')
lonlat_labels(ax3)
ax3.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
ax3.add_feature(cfeature.LAND, zorder=4, edgecolor='black', facecolor='white')

ax4 = fig.add_subplot(224, projection=proj)
cs2=plt.pcolor(EN_lon[:,:22],EN_lat[:,:22],np.transpose(EN_salt_mbot[:22,:]),transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=smin, vmax=smax)
cbars =fig.colorbar(cs2, extend='both')
plt.title('Bottom salinity - EN4')
ax4.gridlines()
ax4.coastlines(resolution='110m')
lonlat_labels(ax4)
ax4.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
ax4.add_feature(cfeature.LAND, zorder=4, edgecolor='black', facecolor='white')

plt.tight_layout()

plt.savefig('figs/WAOM4/waom4xEN4_bottom_TS_annual.png')
plt.close()
 
# call cartopy projection
proj = ccrs.SouthPolarStereo()
fig = plt.figure(figsize=(10,10))

ax1 = fig.add_subplot(111, projection=proj)
ct=plt.pcolor(lon_rho,lat_rho,temp_bot, transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=tmin, vmax=tmax)
cbart =fig.colorbar(ct, extend='both')
plt.title('Bottom temperature - WAOM 4km')
ax1.gridlines() # draw_labels=True,linewidth=2, color='gray', alpha=0.5, linestyle='--')
ax1.coastlines(resolution='110m')
lonlat_labels(ax1)
ax1.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
ax1.add_feature(cfeature.LAND, zorder=4, edgecolor='black', facecolor='white')

plt.savefig('figs/WAOM4/waom4_bottom_TS_annual.png')
plt.close()
