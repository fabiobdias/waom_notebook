# read nc output from WAOM 4km run and plot TS diagrams for East Antarctica sections

import xarray as xr
import pandas as p
import numpy as np
import numpy.ma as ma
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.dates as dates
from matplotlib.ticker import NullFormatter
from matplotlib.colors import LinearSegmentedColormap   # for custom colormaps
import matplotlib as mpl
#mpl.use('Agg')

from datetime import datetime, timedelta

from netCDF4 import Dataset
from netCDF4 import num2date, date2num
from matplotlib.ticker import NullFormatter
from matplotlib.colors import LinearSegmentedColormap   # for custom colormaps

import iris
import iris.iterate
import iris.coords
import iris.plot as iplt
import gsw

# make grid for density contours
smin = 30 - (0.01 * 30)    #salt_ctrl_subregR.min - (0.01 * salt_ctrl_subregR.min)
smax = 36. + (0.01 * 36.)    #salt_ctrl_subregR.max + (0.01 * salt_ctrl_subregR.max)
tmin = -4. + (0.1 * -4.)       #temp_ctrl_subregR.min - (0.1 * temp_ctrl_subregR.max)
tmax = 5 + (0.1 * 5.)       #temp_ctrl_subregR.max + (0.1 * temp_ctrl_subregR.max)
print('tmin, tmax, smin, smax sizes=,', tmin, tmax, smin, smax)
# Calculate how many gridcells we need in the x and y dimensions
xdim = 30
ydim = 20
# Create empty grid of zeros
dens = np.zeros((ydim,xdim))
# Create temp and salt vectors of appropiate dimensions
ti = np.linspace(-4,5,ydim)
si = np.linspace(31,36,xdim)

Si, Ti = np.meshgrid(si, ti, sparse=False, indexing='ij')
# Loop to fill in grid with densities
for j in range(0,int(ydim)):
    for i in range(0, int(xdim)):
        dens[j,i]=gsw.rho(si[i],ti[j],0)
        # Substract 1000 to convert to sigma-0
dens = dens - 1000
print(np.max(dens),np.min(dens))

# load vars for TS plot
count = 0
monthly_ind_ini = [0, 7, 13, 19, 25, 31, 37, 43, 49, 55, 61, 67]
monthly_ind_end = [6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 73]

# load grid to get lat_rho for sections:
dx = xr.open_dataset("/scratch/project_2000339/boeiradi/waom2_frc/waom2extend_grd.nc")
lat_rho = dx.variables["lat_rho"]
lon_rho = dx.variables["lon_rho"]

# index for eta var sections (WWed, EWed, ERos, WRos)
#xi_pt = [875, 1000, 1275, 1550]
#eta_sec_ini = [1600, 1650,  750,  500]
#eta_sec_end = [2200, 2250, 1075, 1196]
#sec_len = [600, 600, 325, 696]

WWed = '/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM2extend_shflim_S_sec_WWed_5days_annual.nc'
print(WWed)
z_rho_WWed_tmp = iris.load_cube(WWed, 'z_rho')
temp_WWed_tmp = iris.load_cube(WWed, 'time-averaged potential temperature') 
salt_WWed_tmp = iris.load_cube(WWed, 'time-averaged salinity') 
temp_WWed = temp_WWed_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
salt_WWed = salt_WWed_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
z_rho_WWed = z_rho_WWed_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
depth_WWed = z_rho_WWed.coord('bathymetry at RHO-points').points
rho_WWed = gsw.rho(xr.DataArray.from_iris(salt_WWed),xr.DataArray.from_iris(temp_WWed),0) - 1000
lat_WWed_tmp = dx.lat_rho.isel(xi_rho=875, eta_rho=slice(1600,2200))
z_rho_WWed2 = xr.DataArray.from_iris(z_rho_WWed)
z_WWed_mask = ma.array(z_rho_WWed2,mask=np.isnan(z_rho_WWed2))
lat_WWed = np.ones((31,600))
for ii in np.arange(0,31):
    lat_WWed[ii,:] = lat_WWed_tmp
lat_WWed_mask =  ma.array(lat_WWed,mask=np.isnan(lat_WWed))
print(temp_WWed_tmp.shape,salt_WWed_tmp.shape,depth_WWed.shape, temp_WWed.shape, z_rho_WWed.shape, depth_WWed.shape)

WRos = '/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM2extend_shflim_S_sec_WRos_5days_annual.nc'
print(WRos)
z_rho_WRos_tmp = iris.load_cube(WRos, 'z_rho')
temp_WRos_tmp = iris.load_cube(WRos, 'time-averaged potential temperature')
salt_WRos_tmp = iris.load_cube(WRos, 'time-averaged salinity')
temp_WRos = temp_WRos_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
salt_WRos = salt_WRos_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
z_rho_WRos = z_rho_WRos_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
depth_WRos = z_rho_WRos.coord('bathymetry at RHO-points').points
rho_WRos = gsw.rho(xr.DataArray.from_iris(salt_WRos),xr.DataArray.from_iris(temp_WRos),0) - 1000
lat_WRos_tmp = dx.lat_rho.isel(xi_rho=1550, eta_rho=slice(500,1196))
z_rho_WRos2 = xr.DataArray.from_iris(z_rho_WRos)
z_WRos_mask = ma.array(z_rho_WRos2,mask=np.isnan(z_rho_WRos2))
lat_WRos = np.ones((31,696))
for ii in np.arange(0,31):
    lat_WRos[ii,:] = lat_WRos_tmp
lat_WRos_mask =  ma.array(lat_WRos,mask=np.isnan(lat_WRos))
print(temp_WRos_tmp.shape,salt_WRos_tmp.shape,depth_WRos.shape, temp_WRos.shape, z_rho_WRos.shape, depth_WRos.shape)

EWed = '/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM2extend_shflim_S_sec_EWed_5days_annual.nc'
print(EWed)
z_rho_EWed_tmp = iris.load_cube(EWed, 'z_rho')
temp_EWed_tmp = iris.load_cube(EWed, 'time-averaged potential temperature')
salt_EWed_tmp = iris.load_cube(EWed, 'time-averaged salinity')
temp_EWed = temp_EWed_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
salt_EWed = salt_EWed_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
z_rho_EWed = z_rho_EWed_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
depth_EWed = z_rho_EWed.coord('bathymetry at RHO-points').points
rho_EWed = gsw.rho(xr.DataArray.from_iris(salt_EWed),xr.DataArray.from_iris(temp_EWed),0) - 1000
lat_EWed_tmp = dx.lat_rho.isel(xi_rho=1000, eta_rho=slice(1650,2250))
z_rho_EWed2 = xr.DataArray.from_iris(z_rho_EWed)
z_EWed_mask = ma.array(z_rho_EWed2,mask=np.isnan(z_rho_EWed2))
lat_EWed = np.ones((31,600))
for ii in np.arange(0,31):
    lat_EWed[ii,:] = lat_EWed_tmp
lat_EWed_mask =  ma.array(lat_EWed,mask=np.isnan(lat_EWed))
print(temp_EWed_tmp.shape,salt_EWed_tmp.shape,depth_EWed.shape, temp_EWed.shape, z_rho_EWed.shape, depth_EWed.shape)

ERos = '/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM2extend_shflim_S_sec_ERos_5days_annual.nc'
print(ERos)
z_rho_ERos_tmp = iris.load_cube(ERos, 'z_rho')
temp_ERos_tmp = iris.load_cube(ERos, 'time-averaged potential temperature')
salt_ERos_tmp = iris.load_cube(ERos, 'time-averaged salinity')
temp_ERos = temp_ERos_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
salt_ERos = salt_ERos_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
z_rho_ERos = z_rho_ERos_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
depth_ERos = z_rho_ERos.coord('bathymetry at RHO-points').points
rho_ERos = gsw.rho(xr.DataArray.from_iris(salt_ERos),xr.DataArray.from_iris(temp_ERos),0) - 1000
lat_ERos_tmp = dx.lat_rho.isel(xi_rho=1275, eta_rho=slice(750,1075))
z_rho_ERos2 = xr.DataArray.from_iris(z_rho_ERos)
z_ERos_mask = ma.array(z_rho_ERos2,mask=np.isnan(z_rho_ERos2))
lat_ERos = np.ones((31,325))
for ii in np.arange(0,31):
    lat_ERos[ii,:] = lat_ERos_tmp
lat_ERos_mask =  ma.array(lat_ERos,mask=np.isnan(lat_ERos))
print(temp_ERos_tmp.shape,salt_ERos_tmp.shape,depth_ERos.shape, temp_ERos.shape, z_rho_ERos.shape, depth_ERos.shape)

s_rho_WWed = z_rho_WWed.coord('S-coordinate at RHO-points').points
s_rho_WRos = z_rho_WRos.coord('S-coordinate at RHO-points').points
s_rho_EWed = z_rho_EWed.coord('S-coordinate at RHO-points').points
s_rho_ERos = z_rho_ERos.coord('S-coordinate at RHO-points').points

## plots:

fig_path='/users/boeiradi/COLD_project/postprocessing/figs/TS_diagrams/'


# plot ts diagram
fig = plt.figure(figsize=(8, 16))

plt.subplot(4,1,1)
plt.gcf().subplots_adjust(bottom=0.15)
print(salt_WWed.shape, temp_WWed.shape,z_rho_WWed.shape)
for s, t in iris.iterate.izip(salt_WWed, temp_WWed, coords='bathymetry at RHO-points'):
    scat = iplt.scatter(s,t, c=depth_WWed, marker='o', edgecolor='none' ,cmap=plt.cm.gist_rainbow, vmin=0, vmax=2000)
ax = plt.gca()
ax.set_ylabel('Temperature ($^{\circ}$C)', fontsize=16)
#ax.set_xlabel('Salinity', fontsize=16)
ax.set_xlim([33.,35.])
ax.set_xticks([33, 33.5, 34, 34.5, 35])
ax.set_ylim([-2.7,1.5])
plt.title('West Ronne Ice Shelf', fontsize=14)
ax.tick_params(labelsize=16)
CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(26.,28.2,.2),linestyles='dashed', colors='k')
cbw = plt.contour(Si,Ti,dens.transpose(), levels=[36.9],linestyles='dashed', colors='m',linewidth=2)
plt.clabel(CS, fontsize=12, inline=1, fmt='%1.1f')

plt.subplot(4,1,2)
plt.gcf().subplots_adjust(bottom=0.15)
print(salt_EWed.shape, temp_EWed.shape,z_rho_EWed.shape)
for s, t in iris.iterate.izip(salt_EWed, temp_EWed, coords='bathymetry at RHO-points'):
    scat = iplt.scatter(s,t, c=depth_EWed, marker='o', edgecolor='none' ,cmap=plt.cm.gist_rainbow, vmin=0, vmax=2000)
ax = plt.gca()
ax.set_ylabel('Temperature ($^{\circ}$C)', fontsize=16)
#ax.set_xlabel('Salinity', fontsize=16)
ax.set_xlim([33.,35.])
ax.set_xticks([33, 33.5, 34, 34.5, 35])
ax.set_ylim([-2.7,1.5])
plt.title('East Ronne Ice Shelf', fontsize=14)
ax.tick_params(labelsize=16)
CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(26.,28.2,.2),linestyles='dashed', colors='k')
cbw = plt.contour(Si,Ti,dens.transpose(), levels=[36.9],linestyles='dashed', colors='m',linewidth=2)
plt.clabel(CS, fontsize=12, inline=1, fmt='%1.1f')

plt.subplot(4,1,3)
plt.gcf().subplots_adjust(bottom=0.15)
print(salt_ERos.shape, temp_ERos.shape,z_rho_ERos.shape)
for s, t in iris.iterate.izip(salt_ERos, temp_ERos, coords='bathymetry at RHO-points'):
    scat = iplt.scatter(s,t, c=depth_ERos, marker='o', edgecolor='none' ,cmap=plt.cm.gist_rainbow, vmin=0, vmax=2000)
ax = plt.gca()
ax.set_xlim([33.,35.])
ax.set_xticks([33, 33.5, 34, 34.5, 35])
ax.set_ylim([-2.7,1.5])
plt.title('East Ross Ice Shelf', fontsize=14)
ax.tick_params(labelsize=16)
CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(26.,28.2,.2),linestyles='dashed', colors='k')
cbw = plt.contour(Si,Ti,dens.transpose(), levels=[36.9],linestyles='dashed', colors='m',linewidth=2)
plt.clabel(CS, fontsize=12, inline=1, fmt='%1.1f')

plt.subplot(4,1,4)
plt.gcf().subplots_adjust(bottom=0.15)
print(salt_WRos.shape, temp_WRos.shape)
for s, t in iris.iterate.izip(salt_WRos, temp_WRos, coords='bathymetry at RHO-points'):
    scat = iplt.scatter(s,t, c=depth_WRos, marker='o', edgecolor='none' ,cmap=plt.cm.gist_rainbow, vmin=0, vmax=2000)
ax = plt.gca()
#cb.set_label('Depth (m)')
ax.set_xlabel('Salinity', fontsize=16)
ax.set_xlim([33.,35.])
ax.set_xticks([33, 33.5, 34, 34.5, 35])
ax.set_ylim([-2.7,1.5])
plt.title('West Ross Ice Shelf', fontsize=14)
ax.tick_params(labelsize=16)
CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(26.,28.2,.2),linestyles='dashed', colors='k')
cbw = plt.contour(Si,Ti,dens.transpose(), levels=[36.9],linestyles='dashed', colors='m',linewidth=2)
plt.clabel(CS, fontsize=12, inline=1, fmt='%1.1f')

cbar_axim = fig.add_axes([0.92, 0.15, 0.015, 0.7])
cb = fig.colorbar(scat, cax=cbar_axim)
cb.ax.tick_params(labelsize=12)
cb.set_label('Depth (m)',fontsize=14)

name_fig="waom2extend_shflim_S_0.25Q_RF-RISsections_TSdiag_depth_yr20.png"
plt.savefig(fig_path + name_fig, dpi=300)

# fig cross-sections + ts digrams RF-RIS regions

# Plot transects
levelsT = np.arange(-2.5,1.,.1)
levelsTf = np.arange(-2.5,1.,.5)
levelsS = np.arange(33.8,34.8,.05)
levelsSf = np.arange(33.8,34.8,.05)
levelsR = np.arange(27.,28.2,.05)
levelsRf = np.arange(27.,28.2,.05)

#f, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 7))
fig = plt.figure(figsize=(20,15))
ax1 = fig.add_subplot(4,4,1)
ct = plt.contourf(lat_WWed_mask, z_WWed_mask, xr.DataArray.from_iris(temp_WWed), levels=levelsT, cmap=plt.cm.coolwarm)
# plt.colorbar(ct, extend='both')
plt.contour(lat_WWed_mask, z_WWed_mask, xr.DataArray.from_iris(temp_WWed), levels=levelsTf, colors='k', linewidths=.5)
plt.xlim([-78.4,-73])
plt.ylim([-2000,0])
plt.title('West Ronne Ice Shelf')
plt.ylabel('Depth (m)')
# plt.xlabel('Latitude')

ax2 = fig.add_subplot(4,4,2)
#xr.DataArray.from_iris(salt_WWed.plot(levels=levelsS)
cs = plt.contourf(lat_WWed_mask, z_WWed_mask, xr.DataArray.from_iris(salt_WWed), levels=levelsS, cmap=plt.cm.coolwarm)
# plt.colorbar(cs, extend='both')
plt.contour(lat_WWed_mask, z_WWed_mask, xr.DataArray.from_iris(salt_WWed), levels=levelsSf, colors='k', linewidths=.5)
plt.xlim([-78.4,-73])
plt.ylim([-2000,0])

ax3 = fig.add_subplot(4,4,3)
ct = plt.contourf(lat_WWed_mask, z_WWed_mask, rho_WWed, levels=levelsR, cmap=plt.cm.coolwarm)
# plt.colorbar(ct, extend='both')
plt.contour(lat_WWed_mask, z_WWed_mask, rho_WWed, levels=levelsRf, colors='k', linewidths=.5)
plt.xlim([-78.4,-73])
plt.ylim([-2000,0])

ax5 = fig.add_subplot(4,4,9)
ct = plt.contourf(lat_WRos_mask, z_WRos_mask, xr.DataArray.from_iris(temp_WRos), levels=levelsT, cmap=plt.cm.coolwarm)
# plt.colorbar(ct, extend='both')
plt.contour(lat_WRos_mask, z_WRos_mask, xr.DataArray.from_iris(temp_WRos), levels=levelsTf, colors='k', linewidths=.5)
plt.xlim([-84,-72])
plt.ylim([-2000,0])
plt.title('West Ross Ice Shelf')
plt.ylabel('Depth (m)')

ax6 = fig.add_subplot(4,4,10)
cs = plt.contourf(lat_WRos_mask, z_WRos_mask, xr.DataArray.from_iris(salt_WRos), levels=levelsS, cmap=plt.cm.coolwarm)
# plt.colorbar(cs, extend='both')
plt.contour(lat_WRos_mask, z_WRos_mask, xr.DataArray.from_iris(salt_WRos), levels=levelsSf, colors='k', linewidths=.5)
plt.xlim([-84,-72])
plt.ylim([-2000,0])

ax7 = fig.add_subplot(4,4,11)
ct = plt.contourf(lat_WRos_mask, z_WRos_mask, rho_WRos, levels=levelsR, cmap=plt.cm.coolwarm)
# plt.colorbar(ct, extend='both')
plt.contour(lat_WRos_mask, z_WRos_mask, rho_WRos, levels=levelsRf, colors='k', linewidths=.5)
plt.xlim([-84,-72])
plt.ylim([-2000,0])

ax9 = fig.add_subplot(4,4,5)
ct = plt.contourf(lat_EWed_mask, z_EWed_mask, xr.DataArray.from_iris(temp_EWed), levels=levelsT, cmap=plt.cm.coolwarm)
# plt.colorbar(ct, extend='both')
plt.contour(lat_EWed_mask, z_EWed_mask, xr.DataArray.from_iris(temp_EWed), levels=levelsTf, colors='k', linewidths=.5)
plt.xlim([-81,-73.7])
plt.ylim([-2000,0])
plt.title('East Ronne Ice Shelf')
plt.ylabel('Depth (m)') # minus is height not depth!

ax10 = fig.add_subplot(4,4,6)
cs = plt.contourf(lat_EWed_mask, z_EWed_mask, xr.DataArray.from_iris(salt_EWed), levels=levelsS, cmap=plt.cm.coolwarm)
# plt.colorbar(cs, extend='both')
plt.contour(lat_EWed_mask, z_EWed_mask, xr.DataArray.from_iris(salt_EWed), levels=levelsSf, colors='k', linewidths=.5)
plt.xlim([-81,-73.7])
plt.ylim([-2000,0])

ax11 = fig.add_subplot(4,4,7)
ct = plt.contourf(lat_EWed_mask, z_EWed_mask, rho_EWed, levels=levelsR, cmap=plt.cm.coolwarm)
# plt.colorbar(ct, extend='both')
plt.contour(lat_EWed_mask, z_EWed_mask, rho_EWed, levels=levelsRf, colors='k', linewidths=.5)
plt.xlim([-81,-73.7])
plt.ylim([-2000,0])

ax13 = fig.add_subplot(4,4,13)
ct = plt.contourf(lat_ERos_mask, z_ERos_mask, xr.DataArray.from_iris(temp_ERos), levels=levelsT, cmap=plt.cm.coolwarm)
plt.contour(lat_ERos_mask, z_ERos_mask, xr.DataArray.from_iris(temp_ERos), levels=levelsTf, colors='k', linewidths=.5)
plt.xlim([-81.3,-75.7])
plt.ylim([-2000,0])
plt.title('East Ross Ice Shelf')
plt.ylabel('Depth (m)') # minus is height not depth!
plt.xlabel('Latitude')
cbar_ax1 = fig.add_axes([0.125, 0.06, 0.17, 0.01])
cbt = fig.colorbar(ct, cax=cbar_ax1, orientation='horizontal')#plt.colorbar(ct, extend='both')
cbt.ax.set_xlabel('Temperature (degC)')

ax14 = fig.add_subplot(4,4,14)
cs = plt.contourf(lat_ERos_mask, z_ERos_mask, xr.DataArray.from_iris(salt_ERos), levels=levelsS, cmap=plt.cm.coolwarm)
plt.contour(lat_ERos_mask, z_ERos_mask, xr.DataArray.from_iris(salt_ERos), levels=levelsSf, colors='k', linewidths=.5)
plt.xlim([-81.3,-75.7])
plt.ylim([-2000,0])
plt.xlabel('Latitude')
cbar_ax2 = fig.add_axes([0.33, 0.06, 0.17, 0.015])
cbs = fig.colorbar(cs, cax=cbar_ax2, orientation='horizontal')
cbs.ax.set_xlabel('Salinity')

ax15 = fig.add_subplot(4,4,15)
cr = plt.contourf(lat_ERos_mask, z_ERos_mask, rho_ERos, levels=levelsR, cmap=plt.cm.coolwarm)
plt.contour(lat_ERos_mask, z_ERos_mask, rho_ERos, levels=levelsRf, colors='k', linewidths=.5)
plt.xlim([-81.3,-75.7])
plt.ylim([-2000,0])
plt.xlabel('Latitude')
cbar_ax3 = fig.add_axes([0.53, 0.06, 0.17, 0.015])
cbr = fig.colorbar(cr, cax=cbar_ax3, orientation='horizontal')
cbr.ax.set_xlabel('Potential density ($\sigma_2$, kg m$^{-3}$)')

# ts-diagram

plt.subplot(4,4,4)
plt.gcf().subplots_adjust(bottom=0.15)

for s, t in iris.iterate.izip(salt_WWed, temp_WWed, coords='bathymetry at RHO-points'):
    scat = iplt.scatter(s,t, c=depth_WWed, marker='o', edgecolor='none' ,cmap=plt.cm.gist_rainbow, vmin=0, vmax=2000)
ax = plt.gca()
cb.set_label('Depth (m)')
ax.set_ylabel('Temperature ($^{\circ}$C)', fontsize=11)
ax.set_xlim([33,35.])
ax.set_ylim([-2.7,1.5])
plt.title('West Ronne Ice Shelf', fontsize=11)
ax.tick_params(labelsize=11)
CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(26.,28.2,.2),linestyles='dashed', colors='k')
cbw = plt.contour(Si,Ti,dens.transpose(), levels=[36.9],linestyles='dashed', colors='m',linewidth=2)
plt.clabel(CS, fontsize=11, inline=1, fmt='%1.1f')

plt.subplot(4,4,12)
plt.gcf().subplots_adjust(bottom=0.15)

for s, t in iris.iterate.izip(salt_WRos, temp_WRos, coords='bathymetry at RHO-points'):
    scat = iplt.scatter(s,t, c=depth_WRos, marker='o', edgecolor='none' ,cmap=plt.cm.gist_rainbow, vmin=0, vmax=2000)
ax = plt.gca()
# cb = plt.colorbar(scat)
cb.set_label('Depth (m)')
ax.set_ylabel('Temperature ($^{\circ}$C)', fontsize=11)
ax.set_xlim([33,35.])
ax.set_ylim([-2.7,1.5])
plt.title('West Ross Ice Shelf', fontsize=11)
ax.tick_params(labelsize=11)
CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(26.,28.2,.2),linestyles='dashed', colors='k')
cbw = plt.contour(Si,Ti,dens.transpose(), levels=[36.9],linestyles='dashed', colors='m',linewidth=2)
plt.clabel(CS, fontsize=11, inline=1, fmt='%1.1f')

plt.subplot(4,4,8)
plt.gcf().subplots_adjust(bottom=0.15)

for s, t in iris.iterate.izip(salt_EWed, temp_EWed, coords='bathymetry at RHO-points'):
    scat = iplt.scatter(s,t, c=depth_EWed, marker='o', edgecolor='none' ,cmap=plt.cm.gist_rainbow, vmin=0, vmax=2000)
ax = plt.gca()
# cb = plt.colorbar(scat)
cb.set_label('Depth (m)')
ax.set_ylabel('Temperature ($^{\circ}$C)', fontsize=11)
ax.set_xlim([33,35.])
ax.set_ylim([-2.7,1.5])
plt.title('East Ronne Ice Shelf', fontsize=11)
ax.tick_params(labelsize=11)
CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(26.,28.2,.2),linestyles='dashed', colors='k')
cbw = plt.contour(Si,Ti,dens.transpose(), levels=[36.9],linestyles='dashed', colors='m',linewidth=2)
plt.clabel(CS, fontsize=11, inline=1, fmt='%1.1f')

plt.subplot(4,4,16)
plt.gcf().subplots_adjust(bottom=0.15)

for s, t in iris.iterate.izip(salt_ERos, temp_ERos, coords='bathymetry at RHO-points'):
    scat = iplt.scatter(s,t, c=depth_ERos, marker='o', edgecolor='none' ,cmap=plt.cm.gist_rainbow, vmin=0, vmax=2000)
ax = plt.gca()
# cb = plt.colorbar(scat)
cb.set_label('Depth (m)')
ax.set_ylabel('Temperature ($^{\circ}$C)', fontsize=11)
ax.set_xlabel('Salinity', fontsize=11)
ax.set_xlim([33,35.])
ax.set_ylim([-2.7,1.5])
plt.title('East Ross Ice Shelf', fontsize=11)
ax.tick_params(labelsize=11)
CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(26.,28.2,.2),linestyles='dashed', colors='k')
cbw = plt.contour(Si,Ti,dens.transpose(), levels=[36.9],linestyles='dashed', colors='m',linewidth=2)
plt.clabel(CS, fontsize=11, inline=1, fmt='%1.1f')
cbar_ax4 = fig.add_axes([0.735, 0.06, 0.17, 0.015])
cbts = fig.colorbar(scat, cax=cbar_ax4, orientation='horizontal')
cbts.ax.set_xlabel('Depth (m)')
# plt.tight_layout()
#plt.show()

name_fig="waom2extend_shflim_S_0.25Q_RF-RISsections_TSdiag+sections_yr20.png"
plt.savefig(fig_path + name_fig, dpi=300)
