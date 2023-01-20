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
        dens[j,i]=gsw.rho(si[i],ti[j],2000)
        # Substract 1000 to convert to sigma-2
dens = dens - 1000
print(np.max(dens),np.min(dens))

# load vars for TS plot
count = 0
monthly_ind_ini = [0, 7, 13, 19, 25, 31, 37, 43, 49, 55, 61, 67]
monthly_ind_end = [6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 73]

wed = '/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM4_shflim_S_sec_wed_10days_annual.nc'
#sec = iris.load('/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM4_shflim_S_ORAS5em_0.25Q_sec_wed_5days_annual.nc')
print(wed)
z_rho_wed_tmp = iris.load_cube(wed, 'z_rho')
temp_wed_tmp = iris.load_cube(wed, 'time-averaged potential temperature') 
salt_wed_tmp = iris.load_cube(wed, 'time-averaged salinity') 
temp_wed = temp_wed_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
salt_wed = salt_wed_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
z_rho_wed = z_rho_wed_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
depth_wed = z_rho_wed.coord('bathymetry at RHO-points').points

print(temp_wed_tmp.shape,salt_wed_tmp.shape,depth_wed.shape, temp_wed.shape, z_rho_wed.shape, depth_wed.shape)

ros = '/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM4_shflim_S_sec_ros_10days_annual.nc'
#sec = iris.load('/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM4_shflim_S_ORAS5em_0.25Q_sec_ros_5days_annual.nc')
print(ros)
z_rho_ros_tmp = iris.load_cube(ros, 'z_rho')
temp_ros_tmp = iris.load_cube(ros, 'time-averaged potential temperature')
salt_ros_tmp = iris.load_cube(ros, 'time-averaged salinity')
temp_ros = temp_ros_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
salt_ros = salt_ros_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
z_rho_ros = z_rho_ros_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
depth_ros = z_rho_ros.coord('bathymetry at RHO-points').points

print(temp_ros_tmp.shape,salt_ros_tmp.shape,depth_ros.shape, temp_ros.shape, z_rho_ros.shape, depth_ros.shape)

pry = '/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM4_shflim_S_sec_pry_10days_annual.nc'
#sec = iris.load('/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM4_shflim_S_ORAS5em_0.25Q_sec_pry_5days_annual.nc')
print(pry)
z_rho_pry_tmp = iris.load_cube(pry, 'z_rho')
temp_pry_tmp = iris.load_cube(pry, 'time-averaged potential temperature')
salt_pry_tmp = iris.load_cube(pry, 'time-averaged salinity')
temp_pry = temp_pry_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
salt_pry = salt_pry_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
z_rho_pry = z_rho_pry_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
depth_pry = z_rho_pry.coord('bathymetry at RHO-points').points

print(temp_pry_tmp.shape,salt_pry_tmp.shape,depth_pry.shape, temp_pry.shape, z_rho_pry.shape, depth_pry.shape)

neu = '/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM4_shflim_S_sec_neu_10days_annual.nc'
#sec = iris.load('/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM4_shflim_S_ORAS5em_0.25Q_sec_neu_5days_annual.nc')
print(neu)
z_rho_neu_tmp = iris.load_cube(neu, 'z_rho')
temp_neu_tmp = iris.load_cube(neu, 'time-averaged potential temperature')
salt_neu_tmp = iris.load_cube(neu, 'time-averaged salinity')
temp_neu = temp_neu_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
salt_neu = salt_neu_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
z_rho_neu = z_rho_neu_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
depth_neu = z_rho_neu.coord('bathymetry at RHO-points').points

print(temp_neu_tmp.shape,salt_neu_tmp.shape,depth_neu.shape, temp_neu.shape, z_rho_neu.shape, depth_neu.shape)

mai = '/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM4_shflim_S_sec_mai_10days_annual.nc'
#sec = iris.load('/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM4_shflim_S_ORAS5em_0.25Q_sec_mai_5days_annual.nc')
print(mai)
z_rho_mai_tmp = iris.load_cube(mai, 'z_rho')
temp_mai_tmp = iris.load_cube(mai, 'time-averaged potential temperature')
salt_mai_tmp = iris.load_cube(mai, 'time-averaged salinity')
temp_mai = temp_mai_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
salt_mai = salt_mai_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
z_rho_mai = z_rho_mai_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
depth_mai = z_rho_mai.coord('bathymetry at RHO-points').points

print(temp_mai_tmp.shape,salt_mai_tmp.shape,depth_mai.shape, temp_mai.shape, z_rho_mai.shape, depth_mai.shape)

dar = '/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM4_shflim_S_sec_dar_10days_annual.nc'
#sec = iris.load('/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM4_shflim_S_ORAS5em_0.25Q_sec_dar_5days_annual.nc')
print(dar)
z_rho_dar_tmp = iris.load_cube(dar, 'z_rho')
temp_dar_tmp = iris.load_cube(dar, 'time-averaged potential temperature')
salt_dar_tmp = iris.load_cube(dar, 'time-averaged salinity')
temp_dar = temp_dar_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
salt_dar = salt_dar_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
z_rho_dar = z_rho_dar_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
depth_dar = z_rho_dar.coord('bathymetry at RHO-points').points

print(temp_dar_tmp.shape,salt_dar_tmp.shape,depth_dar.shape, temp_dar.shape, z_rho_dar.shape, depth_dar.shape)

mer = '/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM4_shflim_S_sec_mer_10days_annual.nc'
#sec = iris.load('/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM4_shflim_S_ORAS5em_0.25Q_sec_mer_5days_annual.nc')
print(mer)
z_rho_mer_tmp = iris.load_cube(mer, 'z_rho')
temp_mer_tmp = iris.load_cube(mer, 'time-averaged potential temperature')
salt_mer_tmp = iris.load_cube(mer, 'time-averaged salinity')
temp_mer = temp_mer_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
salt_mer = salt_mer_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
z_rho_mer = z_rho_mer_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
depth_mer = z_rho_mer.coord('bathymetry at RHO-points').points

print(temp_mer_tmp.shape,salt_mer_tmp.shape,depth_mer.shape, temp_mer.shape, z_rho_mer.shape, depth_mer.shape)

amu = '/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM4_shflim_S_sec_amu_10days_annual.nc'
#sec = iris.load('/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM4_shflim_S_ORAS5em_0.25Q_sec_amu_5days_annual.nc')
print(amu)
z_rho_amu_tmp = iris.load_cube(amu, 'z_rho')
temp_amu_tmp = iris.load_cube(amu, 'time-averaged potential temperature')
salt_amu_tmp = iris.load_cube(amu, 'time-averaged salinity')
temp_amu = temp_amu_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
salt_amu = salt_amu_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
z_rho_amu = z_rho_amu_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
depth_amu = z_rho_amu.coord('bathymetry at RHO-points').points

print(temp_amu_tmp.shape,salt_amu_tmp.shape,depth_amu.shape, temp_amu.shape, z_rho_amu.shape, depth_amu.shape)


s_rho_wed = z_rho_wed.coord('S-coordinate at RHO-points').points
s_rho_ros = z_rho_ros.coord('S-coordinate at RHO-points').points
s_rho_pry = z_rho_pry.coord('S-coordinate at RHO-points').points
s_rho_neu = z_rho_neu.coord('S-coordinate at RHO-points').points
s_rho_mai = z_rho_mai.coord('S-coordinate at RHO-points').points
s_rho_mer = z_rho_mer.coord('S-coordinate at RHO-points').points
s_rho_dar = z_rho_dar.coord('S-coordinate at RHO-points').points
s_rho_amu = z_rho_amu.coord('S-coordinate at RHO-points').points

## plots:

fig_path='/users/boeiradi/COLD_project/postprocessing/figs/TS_diagrams/'

# plot ts diagram
plt.figure(figsize=(8, 16))
plt.subplot(4,1,1)
plt.gcf().subplots_adjust(bottom=0.15)

print(salt_wed)

for s, t in iris.iterate.izip(salt_wed, temp_wed, coords='bathymetry at RHO-points'):
    scat = iplt.scatter(s,t, c=depth_wed, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=200, vmax=4000)
ax = plt.gca()
cb = plt.colorbar(scat)
cb.set_label('Depth (m)')
ax.set_ylabel('Temperature ($^{\circ}$C)', fontsize=16)
#ax.set_xlabel('Salinity', fontsize=16)
ax.set_xlim([33,35.])
ax.set_ylim([-3,5])
plt.title('Weddell Sea', fontsize=14)
ax.tick_params(labelsize=16)
CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(33.5,38,.2),linestyles='dashed', colors='k')
cbw = plt.contour(Si,Ti,dens.transpose(), levels=[36.9],linestyles='dashed', colors='m',linewidth=2)
plt.clabel(CS, fontsize=12, inline=1, fmt='%1.1f')

plt.subplot(4,1,2)
plt.gcf().subplots_adjust(bottom=0.15)

print(salt_ros.shape, temp_ros.shape)

for s, t in iris.iterate.izip(salt_ros, temp_ros, coords='bathymetry at RHO-points'):
    scat = iplt.scatter(s,t, c=depth_ros, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=200, vmax=4000)
ax = plt.gca()
cb = plt.colorbar(scat)
cb.set_label('Depth (m)')
ax.set_ylabel('Temperature ($^{\circ}$C)', fontsize=16)
#ax.set_xlabel('Salinity', fontsize=16)
ax.set_xlim([33,35.])
ax.set_ylim([-3,5])
plt.title('Ross Sea', fontsize=14)
ax.tick_params(labelsize=16)
CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(33.5,38,.2),linestyles='dashed', colors='k')
cbw = plt.contour(Si,Ti,dens.transpose(), levels=[36.9],linestyles='dashed', colors='m',linewidth=2)
plt.clabel(CS, fontsize=12, inline=1, fmt='%1.1f')

#plt.figure(figsize=(8, 5))
plt.subplot(4,1,3)
plt.gcf().subplots_adjust(bottom=0.15)

print(salt_pry.shape, temp_pry.shape)

for s, t in iris.iterate.izip(salt_dar, temp_dar, coords='bathymetry at RHO-points'):
    scat = iplt.scatter(s,t, c=depth_dar, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=200, vmax=4000)
ax = plt.gca()
cb = plt.colorbar(scat)
cb.set_label('Depth (m)')
ax.set_ylabel('Temperature ($^{\circ}$C)', fontsize=16)
ax.set_xlim([33,35.])
ax.set_ylim([-3,5])
plt.title('Darnley Polynya', fontsize=14)
ax.tick_params(labelsize=16)
CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(33.5,38,.2),linestyles='dashed', colors='k')
cbw = plt.contour(Si,Ti,dens.transpose(), levels=[36.9],linestyles='dashed', colors='m',linewidth=2)
plt.clabel(CS, fontsize=12, inline=1, fmt='%1.1f')

plt.subplot(4,1,4)
plt.gcf().subplots_adjust(bottom=0.15)

print(salt_pry.shape, temp_pry.shape)

for s, t in iris.iterate.izip(salt_mer, temp_mer, coords='bathymetry at RHO-points'):
    scat = iplt.scatter(s,t, c=depth_mer, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=200, vmax=4000)
ax = plt.gca()
cb = plt.colorbar(scat)
cb.set_label('Depth (m)')
ax.set_ylabel('Temperature ($^{\circ}$C)', fontsize=16)
ax.set_xlabel('Salinity', fontsize=16)
ax.set_xlim([33,35.])
ax.set_ylim([-3,5])
plt.title('Merz Glacier', fontsize=14)
ax.tick_params(labelsize=16)
CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(33.5,38,.2),linestyles='dashed', colors='k')
cbw = plt.contour(Si,Ti,dens.transpose(), levels=[36.9],linestyles='dashed', colors='m',linewidth=2)
plt.clabel(CS, fontsize=12, inline=1, fmt='%1.1f')

name_fig="waom4_shflim_S_ORAS5em_0.25Q_DSWsections_TSdiag_depth_yr20.png"
plt.savefig(fig_path + name_fig, dpi=300)


# plot ts diagram
fig = plt.figure(figsize=(20, 10))

plt.subplot(2,4,1)
plt.gcf().subplots_adjust(bottom=0.15)
print(salt_wed.shape, temp_wed.shape,z_rho_wed.shape)
for s, t in iris.iterate.izip(salt_wed, temp_wed, coords='bathymetry at RHO-points'):
    scat = iplt.scatter(s,t, c=depth_wed, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=200, vmax=4000)
ax = plt.gca()
cb.set_label('Depth (m)')
ax.set_ylabel('Temperature ($^{\circ}$C)', fontsize=16)
#ax.set_xlabel('Salinity', fontsize=16)
ax.set_xlim([33.,35.])
ax.set_xticks([33, 33.5, 34, 34.5, 35])
ax.set_ylim([-3,5])
plt.title('Weddell Sea', fontsize=14)
ax.tick_params(labelsize=16)
CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(33.5,38,.2),linestyles='dashed', colors='k')
cbw = plt.contour(Si,Ti,dens.transpose(), levels=[36.9],linestyles='dashed', colors='m',linewidth=2)
plt.clabel(CS, fontsize=12, inline=1, fmt='%1.1f')

plt.subplot(2,4,2)
plt.gcf().subplots_adjust(bottom=0.15)
print(salt_neu.shape, temp_neu.shape,z_rho_neu.shape)
for s, t in iris.iterate.izip(salt_neu, temp_neu, coords='bathymetry at RHO-points'):
    scat = iplt.scatter(s,t, c=depth_neu, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=200, vmax=4000)
ax = plt.gca()
cb.set_label('Depth (m)')
ax.set_xlim([33.,35.])
ax.set_xticks([33, 33.5, 34, 34.5, 35])
ax.set_ylim([-3,5])
plt.title('Maud Land 1', fontsize=14)
ax.tick_params(labelsize=16)
CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(33.5,38,.2),linestyles='dashed', colors='k')
cbw = plt.contour(Si,Ti,dens.transpose(), levels=[36.9],linestyles='dashed', colors='m',linewidth=2)
plt.clabel(CS, fontsize=12, inline=1, fmt='%1.1f')

plt.subplot(2,4,3)
plt.gcf().subplots_adjust(bottom=0.15)
for s, t in iris.iterate.izip(salt_mai, temp_mai, coords='bathymetry at RHO-points'):
    scat = iplt.scatter(s,t, c=depth_mai, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=200, vmax=4000)
ax = plt.gca()
cb.set_label('Depth (m)')
ax.set_xlim([33.,35.])
ax.set_xticks([33, 33.5, 34, 34.5, 35])
ax.set_ylim([-3,5])
plt.title('Maud Land 2', fontsize=14)
ax.tick_params(labelsize=16)
CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(33.5,38,.2),linestyles='dashed', colors='k')
cbw = plt.contour(Si,Ti,dens.transpose(), levels=[36.9],linestyles='dashed', colors='m',linewidth=2)
plt.clabel(CS, fontsize=12, inline=1, fmt='%1.1f')

plt.subplot(2,4,4)
plt.gcf().subplots_adjust(bottom=0.15)
for s, t in iris.iterate.izip(salt_dar, temp_dar, coords='bathymetry at RHO-points'):
    scat = iplt.scatter(s,t, c=depth_dar, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=200, vmax=4000)
ax = plt.gca()
cb.set_label('Depth (m)')
ax.set_xlabel('Salinity', fontsize=16)
ax.set_xlim([33.,35.])
ax.set_xticks([33, 33.5, 34, 34.5, 35])
ax.set_ylim([-3,5])
plt.title('Darnley Polynya', fontsize=14)
ax.tick_params(labelsize=16)
CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(33.5,38,.2),linestyles='dashed', colors='k')
cbw = plt.contour(Si,Ti,dens.transpose(), levels=[36.9],linestyles='dashed', colors='m',linewidth=2)
plt.clabel(CS, fontsize=12, inline=1, fmt='%1.1f')

plt.subplot(2,4,5)
plt.gcf().subplots_adjust(bottom=0.15)
for s, t in iris.iterate.izip(salt_mer, temp_mer, coords='bathymetry at RHO-points'):
    scat = iplt.scatter(s,t, c=depth_mer, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=200, vmax=4000)
ax = plt.gca()
cb.set_label('Depth (m)')
ax.set_ylabel('Temperature ($^{\circ}$C)', fontsize=16)
ax.set_xlabel('Salinity', fontsize=16)
ax.set_xlim([33.,35.])
ax.set_xticks([33, 33.5, 34, 34.5, 35])
ax.set_ylim([-3,5])
plt.title('Merz Glacier', fontsize=14)
ax.tick_params(labelsize=16)
CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(33.5,38,.2),linestyles='dashed', colors='k')
cbw = plt.contour(Si,Ti,dens.transpose(), levels=[36.9],linestyles='dashed', colors='m',linewidth=2)
plt.clabel(CS, fontsize=12, inline=1, fmt='%1.1f')

plt.subplot(2,4,6)
plt.gcf().subplots_adjust(bottom=0.15)
print(salt_ros.shape, temp_ros.shape)
for s, t in iris.iterate.izip(salt_ros, temp_ros, coords='bathymetry at RHO-points'):
    scat = iplt.scatter(s,t, c=depth_ros, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=200, vmax=4000)
ax = plt.gca()
cb.set_label('Depth (m)')
ax.set_xlabel('Salinity', fontsize=16)
ax.set_xlim([33.,35.])
ax.set_xticks([33, 33.5, 34, 34.5, 35])
ax.set_ylim([-3,5])
plt.title('Ross Sea', fontsize=14)
ax.tick_params(labelsize=16)
CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(33.5,38,.2),linestyles='dashed', colors='k')
cbw = plt.contour(Si,Ti,dens.transpose(), levels=[36.9],linestyles='dashed', colors='m',linewidth=2)
plt.clabel(CS, fontsize=12, inline=1, fmt='%1.1f')

plt.subplot(2,4,7)
plt.gcf().subplots_adjust(bottom=0.15)
for s, t in iris.iterate.izip(salt_amu, temp_amu, coords='bathymetry at RHO-points'):
    scat = iplt.scatter(s,t, c=depth_amu, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=200, vmax=4000)
ax = plt.gca()
#cb = plt.colorbar(scat)
cb.set_label('Depth (m)')
ax.set_xlabel('Salinity', fontsize=16)
ax.set_xlim([33.,35.])
ax.set_xticks([33, 33.5, 34, 34.5, 35])
ax.set_ylim([-3,5])
plt.title('Amundsen Sea', fontsize=14)
ax.tick_params(labelsize=16)
CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(33.5,38,.2),linestyles='dashed', colors='k')
cbw = plt.contour(Si,Ti,dens.transpose(), levels=[36.9],linestyles='dashed', colors='m',linewidth=2)
plt.clabel(CS, fontsize=12, inline=1, fmt='%1.1f')

cbar_axim = fig.add_axes([0.92, 0.15, 0.015, 0.7])
cb = fig.colorbar(scat, cax=cbar_axim)
cb.ax.tick_params(labelsize=12)
cb.set_label('Depth (m)',fontsize=14)

name_fig="waom4_shflim_S_ORAS5em_0.25Q_AllDSWsections_TSdiag_depth_yr20.png"
plt.savefig(fig_path + name_fig, dpi=300)

# fig cross-sections + ts digrams DSW regions

# Plot transects
levelsT = np.arange(-2.5,2.5,.1)
levelsTf = np.arange(-2.5,2.5,.5)
levelsS = np.arange(33.5,34.9,.025)
levelsSf = np.arange(33.5,34.9,.05)
levelsR = np.arange(36.,37.3,.05)
levelsRf = np.arange(36.,37.3,.05)

#f, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 7))
fig = plt.figure(figsize=(20,15))
ax1 = fig.add_subplot(4,4,1)
ct = plt.contourf(section_lat_wed_mask, section_z_wed_mask, section_temp_wed, levels=levelsT, cmap=plt.cm.coolwarm)
# plt.colorbar(ct, extend='both')
plt.contour(section_lat_wed_mask, section_z_wed_mask, section_temp_wed, levels=levelsTf, colors='k', linewidths=.5)
plt.xlim([-78,-71])
plt.ylim([-1500,0])
plt.title('Weddell Sea')
plt.ylabel('Depth (m)')
# plt.xlabel('Latitude')

ax2 = fig.add_subplot(4,4,2)
#section_salt_wed.plot(levels=levelsS)
cs = plt.contourf(section_lat_wed_mask, section_z_wed_mask, section_salt_wed, levels=levelsS, cmap=plt.cm.coolwarm)
# plt.colorbar(cs, extend='both')
plt.contour(section_lat_wed_mask, section_z_wed_mask, section_salt_wed, levels=levelsSf, colors='k', linewidths=.5)
plt.xlim([-78,-71])
plt.ylim([-1500,0])

ax3 = fig.add_subplot(4,4,3)
ct = plt.contourf(section_lat_wed_mask, section_z_wed_mask, section_rho_wed, levels=levelsR, cmap=plt.cm.coolwarm)
# plt.colorbar(ct, extend='both')
plt.contour(section_lat_wed_mask, section_z_wed_mask, section_rho_wed, levels=levelsRf, colors='k', linewidths=.5)
plt.xlim([-78,-71])
plt.ylim([-1500,0])

ax5 = fig.add_subplot(4,4,5)
ct = plt.contourf(section_lat_ros_mask, section_z_ros_mask, section_temp_ros, levels=levelsT, cmap=plt.cm.coolwarm)
# plt.colorbar(ct, extend='both')
plt.contour(section_lat_ros_mask, section_z_ros_mask, section_temp_ros, levels=levelsTf, colors='k', linewidths=.5)
plt.xlim([-84,-70])
plt.ylim([-1500,0])
plt.title('Ross Sea')
plt.ylabel('Depth (m)')

ax6 = fig.add_subplot(4,4,6)
cs = plt.contourf(section_lat_ros_mask, section_z_ros_mask, section_salt_ros, levels=levelsS, cmap=plt.cm.coolwarm)
# plt.colorbar(cs, extend='both')
plt.contour(section_lat_ros_mask, section_z_ros_mask, section_salt_ros, levels=levelsSf, colors='k', linewidths=.5)
plt.xlim([-84,-70])
plt.ylim([-1500,0])

ax7 = fig.add_subplot(4,4,7)
ct = plt.contourf(section_lat_ros_mask, section_z_ros_mask, section_rho_ros, levels=levelsR, cmap=plt.cm.coolwarm)
# plt.colorbar(ct, extend='both')
plt.contour(section_lat_ros_mask, section_z_ros_mask, section_rho_ros, levels=levelsRf, colors='k', linewidths=.5)
plt.xlim([-84,-70])
plt.ylim([-1500,0])

ax9 = fig.add_subplot(4,4,9)
ct = plt.contourf(section_lat_dar_mask, section_z_dar_mask, section_temp_dar, levels=levelsT, cmap=plt.cm.coolwarm)
# plt.colorbar(ct, extend='both')
plt.contour(section_lat_dar_mask, section_z_dar_mask, section_temp_dar, levels=levelsTf, colors='k', linewidths=.5)
plt.xlim([-70.5,-65])
plt.ylim([-1500,0])
plt.title('Darnley Polynya')
plt.ylabel('Depth (m)') # minus is height not depth!

ax10 = fig.add_subplot(4,4,10)
cs = plt.contourf(section_lat_dar_mask, section_z_dar_mask, section_salt_dar, levels=levelsS, cmap=plt.cm.coolwarm)
# plt.colorbar(cs, extend='both')
plt.contour(section_lat_dar_mask, section_z_dar_mask, section_salt_dar, levels=levelsSf, colors='k', linewidths=.5)
plt.xlim([-70.5,-65])
plt.ylim([-1500,0])

ax11 = fig.add_subplot(4,4,11)
ct = plt.contourf(section_lat_dar_mask, section_z_dar_mask, section_rho_dar, levels=levelsR, cmap=plt.cm.coolwarm)
# plt.colorbar(ct, extend='both')
plt.contour(section_lat_dar_mask, section_z_dar_mask, section_rho_dar, levels=levelsRf, colors='k', linewidths=.5)
plt.xlim([-70.5,-65])
plt.ylim([-1500,0])

ax13 = fig.add_subplot(4,4,13)
ct = plt.contourf(section_lat_mer_mask, section_z_mer_mask, section_temp_mer, levels=levelsT, cmap=plt.cm.coolwarm)
plt.contour(section_lat_mer_mask, section_z_mer_mask, section_temp_mer, levels=levelsTf, colors='k', linewidths=.5)
plt.xlim([-68.,-65])
plt.ylim([-1500,0])
plt.title('Merz Glacier')
plt.ylabel('Depth (m)') # minus is height not depth!
plt.xlabel('Latitude')
cbar_ax1 = fig.add_axes([0.125, 0.06, 0.17, 0.01])
cbt = fig.colorbar(ct, cax=cbar_ax1, orientation='horizontal')#plt.colorbar(ct, extend='both')
cbt.ax.set_xlabel('Temperature (degC)')

ax14 = fig.add_subplot(4,4,14)
cs = plt.contourf(section_lat_mer_mask, section_z_mer_mask, section_salt_mer, levels=levelsS, cmap=plt.cm.coolwarm)
plt.contour(section_lat_mer_mask, section_z_mer_mask, section_salt_mer, levels=levelsSf, colors='k', linewidths=.5)
plt.xlim([-68.,-65])
plt.ylim([-1500,0])
plt.xlabel('Latitude')
cbar_ax2 = fig.add_axes([0.33, 0.06, 0.17, 0.015])
cbs = fig.colorbar(cs, cax=cbar_ax2, orientation='horizontal')
cbs.ax.set_xlabel('Salinity')

ax15 = fig.add_subplot(4,4,15)
cr = plt.contourf(section_lat_mer_mask, section_z_mer_mask, section_rho_mer, levels=levelsR, cmap=plt.cm.coolwarm)
plt.contour(section_lat_mer_mask, section_z_mer_mask, section_rho_mer, levels=levelsRf, colors='k', linewidths=.5)
plt.xlim([-68.,-65])
plt.ylim([-1500,0])
plt.xlabel('Latitude')
cbar_ax3 = fig.add_axes([0.53, 0.06, 0.17, 0.015])
cbr = fig.colorbar(cr, cax=cbar_ax3, orientation='horizontal')
cbr.ax.set_xlabel('Potential density ($\sigma_2$, kg m$^{-3}$)')

# ts-diagram

plt.subplot(4,4,4)
plt.gcf().subplots_adjust(bottom=0.15)

for s, t in iris.iterate.izip(salt_wed, temp_wed, coords='bathymetry at RHO-points'):
    scat = iplt.scatter(s,t, c=depth_wed, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=200, vmax=4000)
ax = plt.gca()
cb.set_label('Depth (m)')
ax.set_ylabel('Temperature ($^{\circ}$C)', fontsize=11)
ax.set_xlim([33,35.])
ax.set_ylim([-3,5])
plt.title('Weddell Sea', fontsize=11)
ax.tick_params(labelsize=11)
CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(33.5,38,.2),linestyles='dashed', colors='k')
cbw = plt.contour(Si,Ti,dens.transpose(), levels=[36.9],linestyles='dashed', colors='m',linewidth=2)
plt.clabel(CS, fontsize=11, inline=1, fmt='%1.1f')

plt.subplot(4,4,8)
plt.gcf().subplots_adjust(bottom=0.15)

print(salt_ros.shape, temp_ros.shape)

for s, t in iris.iterate.izip(salt_ros, temp_ros, coords='bathymetry at RHO-points'):
    scat = iplt.scatter(s,t, c=depth_ros, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=200, vmax=4000)
ax = plt.gca()
# cb = plt.colorbar(scat)
cb.set_label('Depth (m)')
ax.set_ylabel('Temperature ($^{\circ}$C)', fontsize=11)
ax.set_xlim([33,35.])
ax.set_ylim([-3,5])
plt.title('Ross Sea', fontsize=11)
ax.tick_params(labelsize=11)
CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(33.5,38,.2),linestyles='dashed', colors='k')
cbw = plt.contour(Si,Ti,dens.transpose(), levels=[36.9],linestyles='dashed', colors='m',linewidth=2)
plt.clabel(CS, fontsize=11, inline=1, fmt='%1.1f')

plt.subplot(4,4,12)
plt.gcf().subplots_adjust(bottom=0.15)

print(salt_pry.shape, temp_pry.shape)

for s, t in iris.iterate.izip(salt_dar, temp_dar, coords='bathymetry at RHO-points'):
    scat = iplt.scatter(s,t, c=depth_dar, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=200, vmax=4000)
ax = plt.gca()
# cb = plt.colorbar(scat)
cb.set_label('Depth (m)')
ax.set_ylabel('Temperature ($^{\circ}$C)', fontsize=11)
ax.set_xlim([33,35.])
ax.set_ylim([-3,5])
plt.title('Darnley Polynya', fontsize=11)
ax.tick_params(labelsize=11)
CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(33.5,38,.2),linestyles='dashed', colors='k')
cbw = plt.contour(Si,Ti,dens.transpose(), levels=[36.9],linestyles='dashed', colors='m',linewidth=2)
plt.clabel(CS, fontsize=11, inline=1, fmt='%1.1f')

plt.subplot(4,4,16)
plt.gcf().subplots_adjust(bottom=0.15)

print(salt_pry.shape, temp_pry.shape)

for s, t in iris.iterate.izip(salt_dar, temp_dar, coords='bathymetry at RHO-points'):
    scat = iplt.scatter(s,t, c=depth_dar, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=200, vmax=4000)
ax = plt.gca()
# cb = plt.colorbar(scat)
cb.set_label('Depth (m)')
ax.set_ylabel('Temperature ($^{\circ}$C)', fontsize=11)
ax.set_xlim([33,35.])
ax.set_ylim([-3,5])
plt.title('Darnley Polynya', fontsize=11)
ax.tick_params(labelsize=11)
CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(33.5,38,.2),linestyles='dashed', colors='k')
cbw = plt.contour(Si,Ti,dens.transpose(), levels=[36.9],linestyles='dashed', colors='m',linewidth=2)
plt.clabel(CS, fontsize=11, inline=1, fmt='%1.1f')

plt.subplot(4,4,16)
plt.gcf().subplots_adjust(bottom=0.15)

print(salt_pry.shape, temp_pry.shape)

for s, t in iris.iterate.izip(salt_mer, temp_mer, coords='bathymetry at RHO-points'):
    scat = iplt.scatter(s,t, c=depth_mer, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=200, vmax=4000)
ax = plt.gca()
# cb = plt.colorbar(scat)
cb.set_label('Depth (m)')
ax.set_ylabel('Temperature ($^{\circ}$C)', fontsize=11)
ax.set_xlabel('Salinity', fontsize=11)
ax.set_xlim([33,35.])
ax.set_ylim([-3,5])
plt.title('Merz Glacier', fontsize=11)
ax.tick_params(labelsize=11)
CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(33.5,38,.2),linestyles='dashed', colors='k')
cbw = plt.contour(Si,Ti,dens.transpose(), levels=[36.9],linestyles='dashed', colors='m',linewidth=2)
plt.clabel(CS, fontsize=11, inline=1, fmt='%1.1f')
cbar_ax4 = fig.add_axes([0.735, 0.06, 0.17, 0.015])
cbts = fig.colorbar(scat, cax=cbar_ax4, orientation='horizontal')
cbts.ax.set_xlabel('Depth (m)')
# plt.tight_layout()
plt.show()

name_fig="waom4_shflim_S_ORAS5em_0.25Q_DSW+sections_TSdiag_depth_yr20.png"
plt.savefig(fig_path + name_fig, dpi=300)
