# read nc output from WAOM 4km run

import xarray as xr
import pandas as p
import numpy as np
import numpy.ma as ma
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as dates
from matplotlib.ticker import NullFormatter
from matplotlib.colors import LinearSegmentedColormap   # for custom colormaps
from matplotlib import colors
import matplotlib.path as mpath

from datetime import datetime, timedelta

from netCDF4 import Dataset
from netCDF4 import num2date, date2num

import iris
import iris.iterate
import iris.coords
import iris.plot as iplt
import gsw

# load ROMS avg output
for mm  in ['01','02','03','04','05','06','07','08','09','10','11','12']:
    ds = xr.open_dataset('//scratch/project_2000789/boeiradi/waom4_shflim_S_ORAS5em_0.25Q/output_10yr/ocean_avg_00' + mm + '.nc')
    print('MM: size temp and time length: ', mm, ds.temp.shape, len(ds.salt.isel(xi_rho=20, eta_rho=100, s_rho=0)))
    temp_tmp = ds.variables["temp"]
    salt_tmp = ds.variables["salt"]
    temp_sfc_mavg = np.nanmean(temp_tmp[:,-1,:,:], axis=0) # monthly average
    salt_sfc_mavg = np.nanmean(salt_tmp[:,-1,:,:], axis=0)
    temp_bot_mavg = np.nanmean(temp_tmp[:,0,:,:], axis=0) # monthly average
    salt_bot_mavg = np.nanmean(salt_tmp[:,0,:,:], axis=0)
    del temp_tmp, salt_tmp

    # concatenate monthly averages into annual vars:
    if mm == '01':
        temp_sfc_ann = temp_sfc_mavg
        salt_sfc_ann = salt_sfc_mavg
        temp_bot_ann = temp_bot_mavg
        salt_bot_ann = salt_bot_mavg
    elif mm == '02':
        temp_sfc_ann = np.stack((temp_sfc_ann,temp_sfc_mavg), axis=0)
        salt_sfc_ann = np.stack((salt_sfc_ann,salt_sfc_mavg), axis=0)
        temp_bot_ann = np.stack((temp_bot_ann,temp_bot_mavg), axis=0)
        salt_bot_ann = np.stack((salt_bot_ann,salt_bot_mavg), axis=0)
    else:
        temp_sfc_mavg_4thdim = np.expand_dims(temp_sfc_mavg, axis=0)
        temp_sfc_ann = np.concatenate((temp_sfc_ann,temp_sfc_mavg_4thdim), axis=0)
        salt_sfc_mavg_4thdim = np.expand_dims(salt_sfc_mavg, axis=0)
        salt_sfc_ann = np.concatenate((salt_sfc_ann,salt_sfc_mavg_4thdim), axis=0)
        temp_bot_mavg_4thdim = np.expand_dims(temp_bot_mavg, axis=0)
        temp_bot_ann = np.concatenate((temp_bot_ann,temp_bot_mavg_4thdim), axis=0)
        salt_bot_mavg_4thdim = np.expand_dims(salt_bot_mavg, axis=0)
        salt_bot_ann = np.concatenate((salt_bot_ann,salt_bot_mavg_4thdim), axis=0)


# read grid file for lon/lat coordinates
dg = xr.open_dataset("/scratch/project_2000339/boeiradi/waom4_frc/waom4_grd.nc")

lat_rho = dg.variables["lat_rho"]
lon_rho = dg.variables["lon_rho"]
ds.coords['lat_rho']=lat_rho#.transpose() # put lat_rho into ds dataset
ds.coords['lon_rho']=lon_rho#.transpose() # put lon_rho into ds dataset

EN_temp = np.empty((12,42,173,360))
EN_salt = np.empty((12,42,173,360))

count = 0
for mm in ['01','02','03','04','05','06','07','08','09','10','11','12']:
    de_tmp = xr.open_dataset("/scratch/project_2000339/boeiradi/EN/EN.4.2.1.f.analysis.l09.2007" + mm + ".nc")
    EN_temp[count,:] = de_tmp.temperature - 273.15
    EN_salt[count,:] = de_tmp.salinity
    count = count +1

#EN_temp.shape
#EN_salt.shape

de = xr.open_dataset("/scratch/project_2000339/boeiradi/EN/EN.4.2.1.f.analysis.l09.2007" + mm + ".nc")
print(de)

EN_lon, EN_lat = np.meshgrid(de.lon, de.lat, sparse=False, indexing='ij')

ice_mm = np.empty((12,332,316))
count=0
for mm in ['01','02','03','04','05','06','07','08','09','10','11','12']:
    di_tmp = xr.open_dataset("/scratch/project_2000339/boeiradi/NSIDC/seaice_conc_monthly_sh_f13_2007" + mm + "_v03r01.nc")

    ice_mm[count,:,:] = di_tmp.variables["seaice_conc_monthly_cdr"]
    if mm == '01':
        ice_lon = di_tmp.variables["longitude"]
        ice_lat = di_tmp.variables["latitude"]
        ice_xgrid = di_tmp.variables["xgrid"]
        ice_ygrid = di_tmp.variables["ygrid"]
    count = count + 1

### plot some maps
fig_path='/users/boeiradi/COLD_project/postprocessing/figs/Maps_validations/'

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

kw = dict(central_latitude=-90, central_longitude=0, true_scale_latitude=-70)

# time average WAOM10
temp_sfc = np.nanmean(temp_sfc_ann,axis=0)
salt_sfc = np.nanmean(salt_sfc_ann,axis=0)

# time average EN4
EN_temp_sfc = np.mean(EN_temp[:,0,:,:],axis=0)
EN_salt_sfc = np.mean(EN_salt[:,0,:,:],axis=0)

tmin=-2.5
tmax=1.5
smin=34.
smax=35.
levelsT = np.arange(-2.5,1.5,.1)
levelsS = np.arange(34.,35.,.1)

temp_sfc = np.nan_to_num(temp_sfc)
salt_sfc = np.nan_to_num(salt_sfc)

# call cartopy projection
proj = ccrs.SouthPolarStereo()
fig = plt.figure(figsize=(14,10))

ax1 = fig.add_subplot(221, projection=proj)
# sea ice contour
ci=plt.contour(ice_xgrid, ice_ygrid,np.squeeze(ice_mm[8,:,:]), levels=[.75],colors='gray', transform=ccrs.Stereographic(**kw)) #
ci=plt.contour(ice_xgrid, ice_ygrid,np.squeeze(ice_mm[2,:,:]), levels=[.75],colors='white', transform=ccrs.Stereographic(**kw)) #
ct=plt.pcolor(lon_rho,lat_rho,temp_sfc, transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=tmin, vmax=tmax)
#ct=plt.imshow(temp_sfc, extent=[lon_rho.min(), lon_rho.max(), lat_rho.min(), lat_rho.max()], transform=ccrs.PlateCarree(),  vmin=tmin, vmax=tmax, cmap=plt.cm.coolwarm)#, norm=normT_sfc)
cbart =fig.colorbar(ct)# , extend='both')
plt.title('Surface temperature - WAOM 4km')
ax1.gridlines() # draw_labels=True,linewidth=2, color='white', alpha=0.5, linestyle='--')
ax1.coastlines(resolution='110m')
lonlat_labels(ax1)
ax1.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
ax1.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white')

ax2 = fig.add_subplot(222, projection=proj)
ci=plt.contour(ice_xgrid, ice_ygrid,np.squeeze(ice_mm[8,:,:]), levels=[.75],colors='gray', transform=ccrs.Stereographic(**kw)) #
ci=plt.contour(ice_xgrid, ice_ygrid,np.squeeze(ice_mm[2,:,:]), levels=[.75],colors='white', transform=ccrs.Stereographic(**kw)) #
cs=plt.pcolor(lon_rho,lat_rho,salt_sfc, transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=smin, vmax=smax)
cbars =fig.colorbar(cs)# , extend='both')
plt.title('Surface salinity - WAOM 4km')
ax2.coastlines(resolution='110m')
ax2.gridlines()
lonlat_labels(ax2)
ax2.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
ax2.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white')

ax3 = fig.add_subplot(223, projection=proj)
ct2=plt.pcolor(EN_lon[:,:22],EN_lat[:,:22],np.transpose(EN_temp_sfc[:22,:]), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=tmin, vmax=tmax)
cbart =fig.colorbar(ct2)# , extend='both')
plt.title('Surface temperature - EN4')
ax3.gridlines() # draw_labels=True,linewidth=2, color='gray', alpha=0.5, linestyle='--')
ax3.coastlines(resolution='110m')
lonlat_labels(ax3)
ax3.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
ax3.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white')

ax4 = fig.add_subplot(224, projection=proj)
cs2=plt.pcolor(EN_lon[:,:22],EN_lat[:,:22],np.transpose(EN_salt_sfc[:22,:]), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=smin, vmax=smax)
cbars =fig.colorbar(cs2)# , extend='both')
plt.title('Surface salinity - EN4')
ax4.gridlines()
ax4.coastlines(resolution='110m')
lonlat_labels(ax4)
ax4.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
ax4.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white')

#plt.tight_layout()
name_fig = 'waom4_shflim_S_ORAS5em_0.25QxEN4_surface_TS_annual_yr20.png'
plt.savefig(fig_path + name_fig, dpi=300)
plt.close()

# first check the lower value of EN4 depth coordinate which is not-nan ("bottom value"):
EN_nan = np.isnan(EN_temp)
print(EN_nan.shape, EN_temp.shape)

ind_bot = np.empty([12,173,360])
EN_temp_bot = np.empty([12,173,360])
EN_salt_bot = np.empty([12,173,360])

for tt in np.arange(0,12):
    for jj in np.arange(0,173):
        for ii in np.arange(0,360):
            tmp = np.where(np.isfinite(EN_temp[tt,:,jj,ii]))[-1]
            #print('tmp indexes =', EN_temp[0,tmp,jj,ii], tmp)
            if tmp.size == 0:
                ind_bot[tt,jj,ii] = np.nan
                EN_temp_bot[tt,jj,ii] = np.nan
                EN_salt_bot[tt,jj,ii] = np.nan
            else:
                ind_bot[tt,jj,ii] = tmp[-1]
                EN_temp_bot[tt,jj,ii] = np.squeeze(EN_temp[tt,tmp[-1],jj,ii])
                EN_salt_bot[tt,jj,ii] = np.squeeze(EN_salt[tt,tmp[-1],jj,ii])

print(ind_bot.shape, EN_temp_bot.shape, EN_lon.shape)

### plot some maps

# time average WAOM10
temp_bot = np.nanmean(temp_bot_ann,axis=0)
salt_bot = np.nanmean(salt_bot_ann,axis=0)

# time average EN4: bottom values are obtained above;
EN_temp_mbot = np.mean(EN_temp_bot,axis=0)
EN_salt_mbot = np.mean(EN_salt_bot,axis=0)

# call cartopy projection
proj = ccrs.SouthPolarStereo()
fig = plt.figure(figsize=(14,10))

ax1 = fig.add_subplot(221, projection=proj)
ct=plt.pcolor(lon_rho,lat_rho,temp_bot, transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=tmin, vmax=tmax)
cbart =fig.colorbar(ct)# , extend='both')
plt.title('Bottom temperature - WAOM 4km')
ax1.gridlines() # draw_labels=True,linewidth=2, color='gray', alpha=0.5, linestyle='--')
ax1.coastlines(resolution='110m')
lonlat_labels(ax1)
ax1.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
ax1.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white')

ax2 = fig.add_subplot(222, projection=proj)
cs=plt.pcolor(lon_rho,lat_rho,salt_bot, transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=smin, vmax=smax)
cbars =fig.colorbar(cs)# , extend='both')
plt.title('Bottom salinity - WAOM 4km')
ax2.coastlines(resolution='110m')
ax2.gridlines()
lonlat_labels(ax2)
ax2.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
ax2.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white')

ax3 = fig.add_subplot(223, projection=proj)
ct2=plt.pcolor(EN_lon[:,:22],EN_lat[:,:22],np.transpose(EN_temp_mbot[:22,:]), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=tmin, vmax=tmax)
cbart =fig.colorbar(ct2)# , extend='both')
plt.title('Bottom temperature - EN4')
ax3.gridlines() # draw_labels=True,linewidth=2, color='gray', alpha=0.5, linestyle='--')
ax3.coastlines(resolution='110m')
lonlat_labels(ax3)
ax3.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
ax3.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white')

ax4 = fig.add_subplot(224, projection=proj)
cs2=plt.pcolor(EN_lon[:,:22],EN_lat[:,:22],np.transpose(EN_salt_mbot[:22,:]), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=smin, vmax=smax)
cbars =fig.colorbar(cs2)# , extend='both')
plt.title('Bottom salinity - EN4')
ax4.gridlines()
ax4.coastlines(resolution='110m')
lonlat_labels(ax4)
ax4.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
ax4.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white')

#plt.tight_layout()
name_fig ='waom4_shflim_S_ORAS5em_0.25QxEN4_bottom_TS_annual_yr20.png'
plt.savefig(fig_path + name_fig, dpi=300)
plt.close()

# plot some maps

#months_start=[1,32,61,93,124,156,187,219,251,282,314,345]
#months_end=[31,60,92,123,155,186,218,250,281,313,344,365]
#print(months_start,months_end)
#
#
#for mm in np.arange(0,12):
#    swf_range=np.arange(months_start[mm]-1,months_end[mm]-1,1)
#    print(swf_range)
#
#    temp_sfc_mm = ds.temp.isel(s_rho=-1, ocean_time=mm)
#    salt_sfc_mm = ds.salt.isel(s_rho=-1, ocean_time=mm)
#    print('temp,Hsbl size=', temp_sfc_mm.shape)
#    print('EN_temp size=', EN_temp.shape)
#
#    EN_temp_sfc_mm = np.squeeze(EN_temp[mm-1,0,:,:])
#    EN_salt_sfc_mm = np.squeeze(EN_salt[mm-1,0,:,:])
#    print('EN_temp size=', EN_temp_sfc_mm.shape)
#
#    # call cartopy projectione
#    proj = ccrs.SouthPolarStereo()
#    fig = plt.figure(figsize=(14,10))
#
#    ax1 = fig.add_subplot(221, projection=proj)
#    ci=plt.contour(ice_xgrid, ice_ygrid,np.squeeze(ice_mm[mm,:,:]), levels=[.75],colors='dimgray', transform=ccrs.Stereographic(**kw))
#    ct=plt.pcolor(lon_rho,lat_rho,temp_sfc_mm, transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=tmin, vmax=tmax)
#    cbart =fig.colorbar(ct)# , extend='both')
#    plt.title('SST WAOM 4km, m=' + str(mm+1))
#    lonlat_labels(ax1)
#    ax1.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
#
#    ax2 = fig.add_subplot(222, projection=proj)
#    ci=plt.contour(ice_xgrid, ice_ygrid,np.squeeze(ice_mm[mm,:,:]), levels=[.75],colors='dimgray', transform=ccrs.Stereographic(**kw))
#    cs=plt.pcolor(lon_rho,lat_rho,salt_sfc_mm, transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=smin, vmax=smax)
#    cbars =fig.colorbar(cs)# , extend='both')
#    plt.title('SSS WAOM 4km, m=' + str(mm+1))
#    ax1.gridlines()
#    ax2.gridlines()
#    lonlat_labels(ax2)
#    ax2.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
#
#    ax3 = fig.add_subplot(223, projection=proj)
#    ct2=plt.pcolor(EN_lon[:,:22],EN_lat[:,:22],np.transpose(EN_temp_sfc_mm[:22,:]), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=tmin, vmax=tmax)
#    cbart =fig.colorbar(ct2)# , extend='both')
#    plt.title('SST EN4, m=' + str(mm+1))
#    ax3.gridlines()
#    ax3.coastlines(resolution='110m')
#    lonlat_labels(ax3)
#    ax3.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
#
#    ax4 = fig.add_subplot(224, projection=proj)
#    cs2=plt.pcolor(EN_lon[:,:22],EN_lat[:,:22],np.transpose(EN_salt_sfc_mm[:22,:]), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=smin, vmax=smax)#, vmin=32, vmax=35)
#    cbars =fig.colorbar(cs2)# , extend='both')
#    plt.title('SST EN4, m=' + str(mm+1))
#    ax4.gridlines()
#    ax4.coastlines(resolution='110m')
#    lonlat_labels(ax4)
#    ax4.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
#
#    plt.tight_layout()
#
#    month=str(mm+1)
#
#    name_fig='waom4_shflim_S_ORAS5em_0.25QxEN4_surface_TS_m=' + month + '_yr15.png'
#    plt.savefig(fig_path + name_fig)
#
# # plot some maps
#
#for mm in np.arange(0,12):
#    swf_range=np.arange(months_start[mm]-1,months_end[mm]-1,1)
#    print(swf_range)
#
#    temp_bot_mm = ds.temp.isel(s_rho=0, ocean_time=mm)
#    salt_bot_mm = ds.salt.isel(s_rho=0, ocean_time=mm)
#    print('temp,Hsbl size=', temp_bot_mm.shape)
#    print('EN_temp size=', EN_temp.shape)
#
#    EN_temp_bot_mm = np.squeeze(EN_temp_bot[mm-1,:,:])
#    EN_salt_bot_mm = np.squeeze(EN_salt_bot[mm-1,:,:])
#    print('EN_temp size=', EN_temp_bot_mm.shape)
#
#    # call cartopy projectione
#    proj = ccrs.SouthPolarStereo()
#    fig = plt.figure(figsize=(14,10))
#
#    ax1 = fig.add_subplot(221, projection=proj)
#    ci=plt.contour(ice_xgrid, ice_ygrid,np.squeeze(ice_mm[mm,:,:]), levels=[.75],colors='dimgray', transform=ccrs.Stereographic(**kw))
#    ct=plt.pcolor(lon_rho,lat_rho,temp_bot_mm, transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=tmin, vmax=tmax)
#    cbart =fig.colorbar(ct)# , extend='both')
#    plt.title('Temperature WAOM 4km, m=' + str(mm+1))
#    lonlat_labels(ax1)
#    ax1.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
#
#    ax2 = fig.add_subplot(222, projection=proj)
#    ci=plt.contour(ice_xgrid, ice_ygrid,np.squeeze(ice_mm[mm,:,:]), levels=[.75],colors='dimgray', transform=ccrs.Stereographic(**kw))
#    cs=plt.pcolor(lon_rho,lat_rho,salt_bot_mm, transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=smin, vmax=smax)
#    cbars =fig.colorbar(cs)# , extend='both')
#    plt.title('Salinity WAOM 4km, m=' + str(mm+1))
#    ax1.gridlines()
#    ax2.gridlines()
#    lonlat_labels(ax2)
#    ax2.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
#
#    ax3 = fig.add_subplot(223, projection=proj)
#    ct2=plt.pcolor(EN_lon[:,:22],EN_lat[:,:22],np.transpose(EN_temp_bot_mm[:22,:]), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=tmin, vmax=tmax)
#    cbart =fig.colorbar(ct2)# , extend='both')
#    plt.title('Temperature EN4, m=' + str(mm+1))
#    ax3.gridlines()
#    ax3.coastlines(resolution='110m')
#    lonlat_labels(ax3)
#    ax3.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
#
#    ax4 = fig.add_subplot(224, projection=proj)
#    cs2=plt.pcolor(EN_lon[:,:22],EN_lat[:,:22],np.transpose(EN_salt_bot_mm[:22,:]), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=smin, vmax=smax)#, vmin=32, vmax=35)
#    cbars =fig.colorbar(cs2)# , extend='both')
#    plt.title('Temperature EN4, m=' + str(mm+1))
#    ax4.gridlines()
#    ax4.coastlines(resolution='110m')
#    lonlat_labels(ax4)
#    ax4.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
#
#    plt.tight_layout()
#
#    month=str(mm+1)
#
#    
#    name_fig='waom4_shflim_ORAS5em_S_0.25QxEN4_bottom_TS_m=' + month + '_yr15.png'
#    plt.savefig(fig_path + name_fig, dpi=300)


