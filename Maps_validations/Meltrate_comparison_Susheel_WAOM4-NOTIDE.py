# read nc output from WAOM 10km run

import xarray as xr
import pandas as p
import numpy as np
import numpy.ma as ma
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.dates as dates
from matplotlib.ticker import NullFormatter
from matplotlib.colors import LinearSegmentedColormap   # for custom colormaps

from datetime import datetime, timedelta

from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker

from netCDF4 import Dataset
from netCDF4 import num2date, date2num

import gsw

import h5py
from scipy.interpolate import griddata

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# read grid file for lon/lat coordinates
dg = xr.open_dataset("/g/data3/hh5/tmp/access-om/fbd581/ROMS/waom10_frc/waom10extend_grd.nc")
lat_rho_10km= dg.variables["lat_rho"]
lon_rho_10km = dg.variables["lon_rho"]
lat_u_10km= dg.variables["lat_u"]
lon_u_10km = dg.variables["lon_u"]
lat_v_10km= dg.variables["lat_v"]
lon_v_10km = dg.variables["lon_v"]
cor_10km = dg.variables["f"]
pm_10km = dg.variables["pm"]
pn_10km = dg.variables["pn"]
zice_10km = dg.variables["zice"]
h_10km = dg.variables["h"]
dg.close()
print('Print lon/lat_rho shapes',lon_rho_10km.shape, lat_rho_10km.shape)
print('Print lon/lat_rho shapes',lon_rho_10km[0:-1,0:-1].shape, lat_rho_10km[0:-1,0:-1].shape)

dg4 = xr.open_dataset("/g/data3/hh5/tmp/access-om/fbd581/ROMS/waom4_frc/waom4extend_grd.nc")
lat_rho_4km = dg4.variables["lat_rho"]
lon_rho_4km = dg4.variables["lon_rho"]
lat_u_4km= dg4.variables["lat_u"]
lon_u_4km = dg4.variables["lon_u"]
lat_v_4km= dg4.variables["lat_v"]
lon_v_4km = dg4.variables["lon_v"]
cor_4km = dg4.variables["f"]
pm_4km = dg4.variables["pm"]
pn_4km = dg4.variables["pn"]
zice_4km = dg4.variables["zice"]
h_4km = dg4.variables["h"]
dg4.close()

dg2 = xr.open_dataset("/g/data3/hh5/tmp/access-om/fbd581/ROMS/waom2_frc/waom2extend_grd.nc")
lat_rho_2km = dg2.variables["lat_rho"]
lon_rho_2km = dg2.variables["lon_rho"]
lat_u_2km= dg2.variables["lat_u"]
lon_u_2km = dg2.variables["lon_u"]
lat_v_2km= dg2.variables["lat_v"]
lon_v_2km = dg2.variables["lon_v"]
cor_2km = dg2.variables["f"]
pm_2km = dg2.variables["pm"]
pn_2km = dg2.variables["pn"]
zice_2km = dg2.variables["zice"]
h_2km = dg2.variables["h"]
dg2.close()

# load ROMS avg output: 10km
for mm  in ['01','02','03','04','05','06','07','08','09','10','11','12']:
    ds = xr.open_dataset('/g/data3/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom10extend_shflim_S_0.25Q/output_20yr_diag/ocean_avg_00' + mm + '.nc')
    print(ds.variables["temp"].shape)
    m_tmp = np.nanmean(ds.variables["m"], axis=0)
    
    if mm == '01':
        m_10km = m_tmp

    elif mm == '02':
        m_10km = np.stack((m_10km,m_tmp), axis=0)

    else:
        m_tmp_4thdim = np.expand_dims(m_tmp, axis=0)
        m_10km = np.concatenate((m_10km,m_tmp_4thdim), axis=0) 

    ds.close()
    
ds.coords['lat_rho']=lat_rho_10km.transpose() # put lat_rho into ds dataset
ds.coords['lon_rho']=lon_rho_10km.transpose() # put lon_rho into ds dataset

#/scratch/project_2000789/boeiradi/waom4extend_shflim_S_0.25Q/output_yr10_diag/

# load ROMS avg outputWAOM4
for mm  in ['01','02','03','04','05','06','07','08','09','10','11','12']:
    ds = xr.open_dataset('/g/data3/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom4extend_shflim_S_0.25Q/output_yr10_diag/ocean_avg_00' + mm + '.nc')
    print(ds.variables["temp"].shape)
    m_tmp = np.nanmean(ds.variables["m"], axis=0)
    
    if mm == '01':
        m_4km = m_tmp

    elif mm == '02':
        m_4km = np.stack((m_4km,m_tmp), axis=0)

    else:
        m_tmp_4thdim = np.expand_dims(m_tmp, axis=0)
        m_4km = np.concatenate((m_4km,m_tmp_4thdim), axis=0) 

    ds.close()
        
ds.coords['lat_rho']=lat_rho_4km.transpose() # put lat_rho into ds dataset
ds.coords['lon_rho']=lon_rho_4km.transpose() # put lon_rho into ds dataset

#/g/data3/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom2extend_shflim_S_0.25Q/output_yr5

# load ROMS avg output WAOM4-NOTIDE
for mm  in ['01','02','03','04','05','06','07','08','09','10','11','12']:
    ds = xr.open_dataset('/g/data3/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom4extend_shflim_S_0.25Q/output_yr10_notides_diag/ocean_avg_00' + mm + '.nc')
    print(ds.variables["temp"].shape)
    m_tmp = np.nanmean(ds.variables["m"], axis=0)
    
    if mm == '01':
        m_4kmNT = m_tmp

    elif mm == '02':
        m_4kmNT = np.stack((m_4kmNT,m_tmp), axis=0)

    else:
        m_tmp_4thdim = np.expand_dims(m_tmp, axis=0)
        m_4kmNT = np.concatenate((m_4kmNT,m_tmp_4thdim), axis=0) 

    ds.close()
        
# define masks:
mask_zice_10km = ma.masked_where(zice_10km < 0, np.ones(zice_10km.shape))
mask_outice_10km = ma.masked_where(zice_10km >= 0, np.ones(zice_10km.shape))
mask_shelf_10km = ma.masked_where(h_10km > 2000, np.ones(zice_10km.shape))

mask_zice_4km = ma.masked_where(zice_4km < 0, np.ones(zice_4km.shape))
mask_outice_4km = ma.masked_where(zice_4km >= 0, np.ones(zice_4km.shape))
mask_shelf_4km = ma.masked_where(h_4km > 2000, np.ones(zice_4km.shape))


# Load Susheel melta rate data:
# following Susheel github: https://github.com/sioglaciology/ice_shelf_change/blob/master/read_melt_rate_file.ipynb

# filename ='ANT/ANT_iceshelf_melt_rates_CS2_2010-2018_v0.h5'
filename = '/g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/susheel/bb0448974g.h5'
is_wb = h5py.File(filename,'r')

x_wb = np.array(is_wb['/x'])
y_wb = np.array(is_wb['/y'])
wb = np.array(is_wb['/w_b'])

# interpolate Susheel data to WAOM grid:
mr=2
x,y = np.meshgrid(np.arange(-3000,3300+mr/2,mr/2),np.arange(-3000,2600+mr/2,mr/2))
x_rho = x[1::2,1::2]
y_rho = y[1::2,1::2]

x_sus,y_sus = np.meshgrid(x_wb[::2,0],y_wb[::2,0])
x_sus,y_sus = x_sus.flatten()/1000,y_sus.flatten()/1000
wb_sus = wb[::2,::2].flatten()

points = (x_sus,y_sus)
#print(x_rho.shape, y_rho.shape, x_sus.shape, wb_sus.shape)

wb_rho = xr.DataArray(griddata((x_sus,y_sus),wb_sus,(x_rho,y_rho)),dims=('eta_rho','xi_rho'))

### plot some maps
import matplotlib.path as mpath
import cartopy.feature as cfeature
import matplotlib.colors as colors

def lonlat_labels(ax):
    # latitude labels
    ax.text(120,-80,'80$^{\circ}$S',transform=ccrs.PlateCarree(),color='gray')
    ax.text(120,-70,'70$^{\circ}$S',transform=ccrs.PlateCarree(),color='gray')
    # longitude labels
    ax.text(0,-67.5,'0$^{\circ}$',transform=ccrs.PlateCarree(),color='gray')
    #ax.text(60,-53,'60$^{\circ}$E',transform=ccrs.PlateCarree(),color='gray')
    #ax.text(120,-53,'120$^{\circ}$E',transform=ccrs.PlateCarree(),color='gray')
    ax.text(-60,-55,'60$^{\circ}$W',transform=ccrs.PlateCarree(),color='gray')
    ax.text(-120,-55,'120$^{\circ}$W',transform=ccrs.PlateCarree(),color='gray')
    ax.text(180,-65.5,'180$^{\circ}$',transform=ccrs.PlateCarree(),color='gray')
    return

proj = ccrs.SouthPolarStereo()

# plotting surface fluxes:
fig_path = '/g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/figs/Melt_rates/'

bathym = cfeature.NaturalEarthFeature(name='bathymetry_J_1000', scale='10m', category='physical')
# limits for contour of ice front (Ronne-Filchner IS):
xlimit = np.arange(300,500,1)
ylimit = np.arange(100,300,1)

# define your scale, with white at zero
Vmin = -4 
Vmax = 4
#norm = colors.TwoSlopeNorm(vmin=Vmin, vcenter=0, vmax=Vmax)

proj = ccrs.SouthPolarStereo()
fig = plt.figure(figsize=(8,10))
 
ax1 = fig.add_subplot(222, projection=proj)
plt.title('B) WAOM10')
cy=plt.pcolormesh(lon_rho_10km,lat_rho_10km,np.nanmean(m_10km, axis=0)*86400*365.25, transform=ccrs.PlateCarree(), cmap=plt.cm.seismic, vmin=Vmin, vmax=Vmax)#, norm=norm) 
gl = ax1.gridlines(zorder=4,draw_labels=True, dms=False, x_inline=False, y_inline=False)
gl.rotate_labels = False
gl.ylabels_right = False
gl.xlines = True
gl.xlocator = mticker.FixedLocator(np.arange(-180,180,30))
gl.ylocator = mticker.FixedLocator(np.arange(-50,-90,-10))
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
ax1.set_extent([-180, 180, -90, -65], crs=ccrs.PlateCarree())
ax1.add_feature(cfeature.LAND, zorder=3, facecolor='lightgrey')
#plt.contour(lon_rho_10km, lat_rho_10km,h_10km,levels=(2000,2001), transform=ccrs.PlateCarree(), colors='k')
#plt.contour(lon_rho_10km, lat_rho_10km,zice_10km,levels=[-.1], transform=ccrs.PlateCarree(), cmap=plt.cm.binary, linewidths=1)

ax2 = fig.add_subplot(223, projection=proj)
plt.title('C) WAOM4')
cy=plt.pcolormesh(lon_rho_4km,lat_rho_4km,np.nanmean(m_4km, axis=0)*86400*365.25, transform=ccrs.PlateCarree(), cmap=plt.cm.seismic, vmin=Vmin, vmax=Vmax)#, norm=norm) 
gl = ax2.gridlines(zorder=4,draw_labels=True, dms=False, x_inline=False, y_inline=False)
gl.rotate_labels = False
gl.ylabels_right = False
gl.xlines = True
gl.xlocator = mticker.FixedLocator(np.arange(-180,180,30))
gl.ylocator = mticker.FixedLocator(np.arange(-50,-90,-10))
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
ax2.set_extent([-180, 180, -90, -65], crs=ccrs.PlateCarree())
ax2.add_feature(cfeature.LAND, zorder=3, facecolor='lightgrey')
#plt.contour(lon_rho_10km, lat_rho_10km,h_10km,levels=(2000,2001), transform=ccrs.PlateCarree(), colors='k')
#plt.contour(lon_rho_10km, lat_rho_10km,zice_10km,levels=[-.1], transform=ccrs.PlateCarree(), cmap=plt.cm.binary, linewidths=1)

ax3 = fig.add_subplot(224, projection=proj)
plt.title('D) WAOM4-NOTIDE minus WAOM4')
cy=plt.pcolormesh(lon_rho_4km,lat_rho_4km,(np.nanmean(m_4kmNT, axis=0)-np.nanmean(m_4km, axis=0))*86400*365.25, transform=ccrs.PlateCarree(), cmap=plt.cm.seismic, vmin=Vmin, vmax=Vmax)#, norm=norm)
#plt.title('D) WAOM4-NOTIDE')
#cy=plt.pcolormesh(lon_rho_4km,lat_rho_4km,np.nanmean(m_4kmNT, axis=0)*86400*365.25, transform=ccrs.PlateCarree(), cmap=plt.cm.seismic, vmin=Vmin, vmax=Vmax)#, norm=norm) 
gl = ax3.gridlines(zorder=4,draw_labels=True, dms=False, x_inline=False, y_inline=False)
gl.rotate_labels = False
gl.ylabels_right = False
gl.xlines = True
gl.xlocator = mticker.FixedLocator(np.arange(-180,180,30))
gl.ylocator = mticker.FixedLocator(np.arange(-50,-90,-10))
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
ax3.set_extent([-180, 180, -90, -65], crs=ccrs.PlateCarree())
ax3.add_feature(cfeature.LAND, zorder=3, facecolor='lightgrey')
#plt.contour(lon_rho_10km, lat_rho_10km,h_10km,levels=(2000,2001), transform=ccrs.PlateCarree(), colors='k')
#plt.contour(lon_rho_10km, lat_rho_10km,zice_10km,levels=[-.1], transform=ccrs.PlateCarree(), cmap=plt.cm.binary, linewidths=1)

ax = fig.add_subplot(221, projection=proj)
plt.pcolormesh(lon_rho_2km,lat_rho_2km, wb_rho, transform=ccrs.PlateCarree(),  cmap=plt.cm.seismic, vmin=Vmin, vmax=Vmax)#, norm=norm)
plt.title('A) Adusumilli et al. (2020)')
gl = ax.gridlines(zorder=4,draw_labels=True, dms=False, x_inline=False, y_inline=False)
gl.rotate_labels = False
gl.ylabels_right = False
gl.xlines = True
gl.xlocator = mticker.FixedLocator(np.arange(-180,180,30))
gl.ylocator = mticker.FixedLocator(np.arange(-50,-90,-10))
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
ax.set_extent([-180, 180, -90, -65], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.LAND, zorder=3, facecolor='lightgrey')
#plt.contour(lon_rho_10km, lat_rho_10km,h_10km,levels=(2000,2001), transform=ccrs.PlateCarree(), colors='k')
#plt.contour(lon_rho_10km, lat_rho_10km,zice_10km,levels=[-.1], transform=ccrs.PlateCarree(), cmap=plt.cm.binary, linewidths=1)

cbar_ax1 = fig.add_axes([0.13, 0.06, 0.76, 0.01])
fig.colorbar(cy, cax=cbar_ax1, orientation='horizontal')
cbar_ax1.set_xlabel('Basal melt rate (m/yr)',fontsize=12)#, labelpad=-35)
                                                
name_fig="Susheel_waom10x4x4-notides_melt_rates_maps_annual_yr20.png"
plt.savefig(fig_path + name_fig, dpi=300)
plt.close()

# for FRIS region:
fig = plt.figure(figsize=(8,8))

ax1 = fig.add_subplot(222, projection=proj)
plt.title('B) WAOM10')
cy=plt.pcolormesh(lon_rho_10km,lat_rho_10km,np.nanmean(m_10km, axis=0)*86400*365.25, transform=ccrs.PlateCarree(), cmap=plt.cm.seismic, vmin=Vmin, vmax=Vmax)#, norm=norm)
gl = ax1.gridlines(zorder=4,draw_labels=True, dms=False, x_inline=False, y_inline=False)
gl.rotate_labels = False
gl.ylabels_right = False
gl.xlines = True
gl.xlocator = mticker.FixedLocator(np.arange(-180,180,10))
gl.ylocator = mticker.FixedLocator(np.arange(-50,-90,-10))
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
ax1.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
ax1.add_feature(cfeature.LAND, zorder=3, facecolor='lightgrey')
ax1.add_feature(bathym, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='grey')

ax2 = fig.add_subplot(223, projection=proj)
plt.title('C) WAOM4')
cy=plt.pcolormesh(lon_rho_4km,lat_rho_4km,np.nanmean(m_4km, axis=0)*86400*365.25, transform=ccrs.PlateCarree(), cmap=plt.cm.seismic, vmin=Vmin, vmax=Vmax)#, norm=norm)
gl = ax2.gridlines(zorder=4,draw_labels=True, dms=False, x_inline=False, y_inline=False)
gl.rotate_labels = False
gl.ylabels_right = False
gl.xlines = True
gl.xlocator = mticker.FixedLocator(np.arange(-180,180,10))
gl.ylocator = mticker.FixedLocator(np.arange(-50,-90,-10))
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
ax2.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
ax2.add_feature(cfeature.LAND, zorder=3, facecolor='lightgrey')
ax2.add_feature(bathym, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='grey')

ax3 = fig.add_subplot(224, projection=proj)
plt.title('D) WAOM4-NOTIDE minus WAOM4')
cy=plt.pcolormesh(lon_rho_4km,lat_rho_4km,(np.nanmean(m_4kmNT, axis=0)-np.nanmean(m_4km, axis=0))*86400*365.25, transform=ccrs.PlateCarree(), cmap=plt.cm.seismic, vmin=Vmin, vmax=Vmax)#, norm=norm)
gl = ax3.gridlines(zorder=4,draw_labels=True, dms=False, x_inline=False, y_inline=False)
gl.rotate_labels = False
gl.ylabels_right = False
gl.xlines = True
gl.xlocator = mticker.FixedLocator(np.arange(-180,180,10))
gl.ylocator = mticker.FixedLocator(np.arange(-50,-90,-10))
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
ax3.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
ax3.add_feature(cfeature.LAND, zorder=3, facecolor='lightgrey')
ax3.add_feature(bathym, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='grey')

ax = fig.add_subplot(221, projection=proj)
plt.pcolormesh(lon_rho_2km,lat_rho_2km, wb_rho, transform=ccrs.PlateCarree(),  cmap=plt.cm.seismic, vmin=Vmin, vmax=Vmax)#, norm=norm)
plt.title('A) Adusumilli et al. (2020)')
ax.gridlines() # draw_labels=True,linewidth=2, color='white', alpha=0.5, linestyle='--')
lonlat_labels(ax)
ax.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.LAND, zorder=3, facecolor='lightgrey')
ax.add_feature(bathym, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='grey')

cbar_ax1 = fig.add_axes([0.13, 0.08, 0.76, 0.01])
fig.colorbar(cy, cax=cbar_ax1, orientation='horizontal')
cbar_ax1.set_xlabel('Melt rate (m.yr$^{-1}$)',fontsize=12)#, labelpad=-35)

name_fig="Susheel_waom10x4x4-notides_melt_rates_maps_annual_yr20_FRIS.png"

plt.savefig(fig_path + name_fig, dpi=300)
plt.close()

