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

# load grid to get lat_rho for sections:
dx = xr.open_dataset("/scratch/project_2000339/boeiradi/waom4_frc/waom4extend_grd.nc")
lat_rho = dx.variables["lat_rho"]
lon_rho = dx.variables["lon_rho"]
#ds.coords['lat_rho']=lat_rho.transpose() # put lat_rho into ds dataset
#ds.coords['lon_rho']=lon_rho.transpose() # put lon_rho into ds dataset

# load section vars:
wed = '/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM4extend_shflim_S_sec_wed_5days_annual.nc'
#sec = iris.load('/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM4extend_shflim_S_ORAS5em_0.25Q_sec_wed_5days_annual.nc')
print(wed)
z_rho_wed_tmp = iris.load_cube(wed, 'z_rho')
temp_wed_tmp = iris.load_cube(wed, 'time-averaged potential temperature') 
salt_wed_tmp = iris.load_cube(wed, 'time-averaged salinity') 
temp_wed = temp_wed_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
salt_wed = salt_wed_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
z_rho_wed = z_rho_wed_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
depth_wed = z_rho_wed.coord('bathymetry at RHO-points').points
print(temp_wed_tmp.shape,salt_wed_tmp.shape,depth_wed.shape, temp_wed.shape, z_rho_wed.shape, depth_wed.shape)
rho_wed = gsw.rho(xr.DataArray.from_iris(salt_wed),xr.DataArray.from_iris(temp_wed),2000) - 1000
lat_wed_tmp = dx.lat_rho.isel(xi_rho=425, eta_rho=slice(775,1200))
z_rho_wed2 = xr.DataArray.from_iris(z_rho_wed)
z_wed_mask = ma.array(z_rho_wed2,mask=np.isnan(z_rho_wed2))
lat_wed = np.ones((31,425))
for ii in np.arange(0,31):
    lat_wed[ii,:] = lat_wed_tmp
lat_wed_mask =  ma.array(lat_wed,mask=np.isnan(lat_wed))

ros = '/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM4extend_shflim_S_sec_ros_5days_annual.nc'
#sec = iris.load('/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM4extend_shflim_S_ORAS5em_0.25Q_sec_ros_5days_annual.nc')
print(ros)
z_rho_ros_tmp = iris.load_cube(ros, 'z_rho')
temp_ros_tmp = iris.load_cube(ros, 'time-averaged potential temperature')
salt_ros_tmp = iris.load_cube(ros, 'time-averaged salinity')
temp_ros = temp_ros_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
salt_ros = salt_ros_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
z_rho_ros = z_rho_ros_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
depth_ros = z_rho_ros.coord('bathymetry at RHO-points').points
print(temp_ros_tmp.shape,salt_ros_tmp.shape,depth_ros.shape, temp_ros.shape, z_rho_ros.shape, depth_ros.shape)
rho_ros = gsw.rho(xr.DataArray.from_iris(salt_ros),xr.DataArray.from_iris(temp_ros),2000) - 1000
lat_ros_tmp = dx.lat_rho.isel(xi_rho=775, eta_rho=slice(125,700))
z_ros_mask = ma.array(xr.DataArray.from_iris(z_rho_ros),mask=np.isnan(xr.DataArray.from_iris(z_rho_ros)))
lat_ros = np.ones((31,575))
for ii in np.arange(0,31):
    lat_ros[ii,:] = lat_ros_tmp
lat_ros_mask =  ma.array(lat_ros,mask=np.isnan(lat_ros))

pry = '/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM4extend_shflim_S_sec_pry_5days_annual.nc'
#sec = iris.load('/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM4extend_shflim_S_ORAS5em_0.25Q_sec_pry_5days_annual.nc')
print(pry)
z_rho_pry_tmp = iris.load_cube(pry, 'z_rho')
temp_pry_tmp = iris.load_cube(pry, 'time-averaged potential temperature')
salt_pry_tmp = iris.load_cube(pry, 'time-averaged salinity')
temp_pry = temp_pry_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
salt_pry = salt_pry_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
z_rho_pry = z_rho_pry_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
depth_pry = z_rho_pry.coord('bathymetry at RHO-points').points
print(temp_pry_tmp.shape,salt_pry_tmp.shape,depth_pry.shape, temp_pry.shape, z_rho_pry.shape, depth_pry.shape)
rho_pry = gsw.rho(xr.DataArray.from_iris(salt_pry),xr.DataArray.from_iris(temp_pry),2000) - 1000
lat_pry_tmp = dx.lat_rho.isel(eta_rho=920, xi_rho=slice(1050,1500))
z_pry_mask = ma.array(xr.DataArray.from_iris(z_rho_pry),mask=np.isnan(xr.DataArray.from_iris(z_rho_pry)))
lat_pry = np.ones((31,450))
for ii in np.arange(0,31):
    lat_pry[ii,:] = lat_pry_tmp
lat_pry_mask =  ma.array(lat_pry,mask=np.isnan(lat_pry))

neu = '/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM4extend_shflim_S_sec_neu_5days_annual.nc'
#sec = iris.load('/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM4extend_shflim_S_ORAS5em_0.25Q_sec_neu_5days_annual.nc')
print(neu)
z_rho_neu_tmp = iris.load_cube(neu, 'z_rho')
temp_neu_tmp = iris.load_cube(neu, 'time-averaged potential temperature')
salt_neu_tmp = iris.load_cube(neu, 'time-averaged salinity')
temp_neu = temp_neu_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
salt_neu = salt_neu_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
z_rho_neu = z_rho_neu_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
depth_neu = z_rho_neu.coord('bathymetry at RHO-points').points
print(temp_neu_tmp.shape,salt_neu_tmp.shape,depth_neu.shape, temp_neu.shape, z_rho_neu.shape, depth_neu.shape)
rho_neu = gsw.rho(xr.DataArray.from_iris(salt_neu),xr.DataArray.from_iris(temp_neu),2000) - 1000
lat_neu_tmp = dx.lat_rho.isel(xi_rho=700, eta_rho=slice(1125,1375))
z_neu_mask = ma.array(xr.DataArray.from_iris(z_rho_neu),mask=np.isnan(xr.DataArray.from_iris(z_rho_neu)))
lat_neu = np.ones((31,250))
for ii in np.arange(0,31):
    lat_neu[ii,:] = lat_neu_tmp
lat_neu_mask =  ma.array(lat_neu,mask=np.isnan(lat_neu))

mai = '/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM4extend_shflim_S_sec_mai_5days_annual.nc'
#sec = iris.load('/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM4extend_shflim_S_ORAS5em_0.25Q_sec_mai_5days_annual.nc')
print(mai)
z_rho_mai_tmp = iris.load_cube(mai, 'z_rho')
temp_mai_tmp = iris.load_cube(mai, 'time-averaged potential temperature')
salt_mai_tmp = iris.load_cube(mai, 'time-averaged salinity')
temp_mai = temp_mai_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
salt_mai = salt_mai_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
z_rho_mai = z_rho_mai_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
depth_mai = z_rho_mai.coord('bathymetry at RHO-points').points
print(temp_mai_tmp.shape,salt_mai_tmp.shape,depth_mai.shape, temp_mai.shape, z_rho_mai.shape, depth_mai.shape)
rho_mai = gsw.rho(xr.DataArray.from_iris(salt_mai),xr.DataArray.from_iris(temp_mai),2000) - 1000
lat_mai_tmp = dx.lat_rho.isel(xi_rho=800, eta_rho=slice(1125,1375))
z_mai_mask = ma.array(xr.DataArray.from_iris(z_rho_mai),mask=np.isnan(xr.DataArray.from_iris(z_rho_mai)))
lat_mai = np.ones((31,250))
for ii in np.arange(0,31):
    lat_mai[ii,:] = lat_mai_tmp
lat_mai_mask =  ma.array(lat_mai,mask=np.isnan(lat_mai))

dar = '/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM4extend_shflim_S_sec_dar_5days_annual.nc'
#sec = iris.load('/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM4extend_shflim_S_ORAS5em_0.25Q_sec_dar_5days_annual.nc')
print(dar)
z_rho_dar_tmp = iris.load_cube(dar, 'z_rho')
temp_dar_tmp = iris.load_cube(dar, 'time-averaged potential temperature')
salt_dar_tmp = iris.load_cube(dar, 'time-averaged salinity')
temp_dar = temp_dar_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
salt_dar = salt_dar_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
z_rho_dar = z_rho_dar_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
depth_dar = z_rho_dar.coord('bathymetry at RHO-points').points
print(temp_dar_tmp.shape,salt_dar_tmp.shape,depth_dar.shape, temp_dar.shape, z_rho_dar.shape, depth_dar.shape)
rho_dar = gsw.rho(xr.DataArray.from_iris(salt_dar),xr.DataArray.from_iris(temp_dar),2000) - 1000
lat_dar_tmp = dx.lat_rho.isel(eta_rho=950, xi_rho=slice(1150,1538))
z_dar_mask = ma.array(xr.DataArray.from_iris(z_rho_dar),mask=np.isnan(xr.DataArray.from_iris(z_rho_dar)))
lat_dar = np.ones((31,388))
for ii in np.arange(0,31):
    lat_dar[ii,:] = lat_dar_tmp
lat_dar_mask =  ma.array(lat_dar,mask=np.isnan(lat_dar))

mer = '/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM4extend_shflim_S_sec_mer_5days_annual.nc'
#sec = iris.load('/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM4extend_shflim_S_ORAS5em_0.25Q_sec_mer_5days_annual.nc')
print(mer)
z_rho_mer_tmp = iris.load_cube(mer, 'z_rho')
temp_mer_tmp = iris.load_cube(mer, 'time-averaged potential temperature')
salt_mer_tmp = iris.load_cube(mer, 'time-averaged salinity')
temp_mer = temp_mer_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
salt_mer = salt_mer_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
z_rho_mer = z_rho_mer_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
depth_mer = z_rho_mer.coord('bathymetry at RHO-points').points
print(temp_mer_tmp.shape,salt_mer_tmp.shape,depth_mer.shape, temp_mer.shape, z_rho_mer.shape, depth_mer.shape)
rho_mer = gsw.rho(xr.DataArray.from_iris(salt_mer),xr.DataArray.from_iris(temp_mer),2000) - 1000
lat_mer_tmp = dx.lat_rho.isel(xi_rho=1088, eta_rho=slice(100,275))
z_mer_mask = ma.array(xr.DataArray.from_iris(z_rho_mer),mask=np.isnan(xr.DataArray.from_iris(z_rho_mer)))
lat_mer = np.ones((31,175))
for ii in np.arange(0,31):
    lat_mer[ii,:] = lat_mer_tmp
lat_mer_mask =  ma.array(lat_mer,mask=np.isnan(lat_mer))

amu = '/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM4extend_shflim_S_sec_amu_5days_annual.nc'
#sec = iris.load('/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM4extend_shflim_S_ORAS5em_0.25Q_sec_amu_5days_annual.nc')
print(amu)
z_rho_amu_tmp = iris.load_cube(amu, 'z_rho')
temp_amu_tmp = iris.load_cube(amu, 'time-averaged potential temperature')
salt_amu_tmp = iris.load_cube(amu, 'time-averaged salinity')
temp_amu = temp_amu_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
salt_amu = salt_amu_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
z_rho_amu = z_rho_amu_tmp.collapsed('averaged time since initialization', iris.analysis.MEAN)
depth_amu = z_rho_amu.coord('bathymetry at RHO-points').points
print(temp_amu_tmp.shape,salt_amu_tmp.shape,depth_amu.shape, temp_amu.shape, z_rho_amu.shape, depth_amu.shape)
rho_amu = gsw.rho(xr.DataArray.from_iris(salt_amu),xr.DataArray.from_iris(temp_amu),2000) - 1000
lat_amu_tmp = dx.lat_rho.isel(eta_rho=563, xi_rho=slice(125,500))
z_amu_mask = ma.array(xr.DataArray.from_iris(z_rho_amu),mask=np.isnan(xr.DataArray.from_iris(z_rho_amu)))
lat_amu = np.ones((31,375))
for ii in np.arange(0,31):
    lat_amu[ii,:] = lat_amu_tmp
lat_amu_mask =  ma.array(lat_amu,mask=np.isnan(lat_amu))

s_rho_wed = z_rho_wed.coord('sea_surface_height_above_reference_ellipsoid').points
s_rho_ros = z_rho_ros.coord('sea_surface_height_above_reference_ellipsoid').points
s_rho_pry = z_rho_pry.coord('sea_surface_height_above_reference_ellipsoid').points
s_rho_neu = z_rho_neu.coord('sea_surface_height_above_reference_ellipsoid').points
s_rho_mai = z_rho_mai.coord('sea_surface_height_above_reference_ellipsoid').points
s_rho_mer = z_rho_mer.coord('sea_surface_height_above_reference_ellipsoid').points
s_rho_dar = z_rho_dar.coord('sea_surface_height_above_reference_ellipsoid').points
s_rho_amu = z_rho_amu.coord('sea_surface_height_above_reference_ellipsoid').points

## plots:

fig_path='/users/boeiradi/COLD_project/postprocessing/figs/TS_diagrams/'

# plot ts diagram
plt.figure(figsize=(8, 16))
plt.subplot(4,1,1)
plt.gcf().subplots_adjust(bottom=0.15)

print(salt_wed)

for s, t, z in iris.iterate.izip(salt_wed, temp_wed, z_rho_wed, coords='S-coordinate at RHO-points'):
    print('z inside irits iterate: ',z)
    sr_tmp = z.coord('S-coordinate at RHO-points').points
    dep_tmp = z.coord('bathymetry at RHO-points').points
    print('s_rho and dep points: ', sr_tmp, dep_tmp, sr_tmp*dep_tmp)
    scat = iplt.scatter(s,t,c=-sr_tmp*dep_tmp, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=0, vmax=4000)
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

for s, t, z in iris.iterate.izip(salt_ros, temp_ros, z_rho_ros, coords='S-coordinate at RHO-points'):
    sr_tmp = z.coord('S-coordinate at RHO-points').points
    dep_tmp = z.coord('bathymetry at RHO-points').points
    scat = iplt.scatter(s,t, c=-sr_tmp*dep_tmp, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=0, vmax=4000)
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

for s, t, z in iris.iterate.izip(salt_dar, temp_dar, z_rho_dar, coords='S-coordinate at RHO-points'):
    sr_tmp = z.coord('S-coordinate at RHO-points').points
    dep_tmp = z.coord('bathymetry at RHO-points').points
    scat = iplt.scatter(s,t, c=-sr_tmp*dep_tmp, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=0, vmax=4000)
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

for s, t, z in iris.iterate.izip(salt_mer, temp_mer, z_rho_mer, coords='S-coordinate at RHO-points'):
    sr_tmp = z.coord('S-coordinate at RHO-points').points
    dep_tmp = z.coord('bathymetry at RHO-points').points
    scat = iplt.scatter(s,t, c=-sr_tmp*dep_tmp, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=0, vmax=4000)
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

name_fig="waom4extend_shflim_S_0.25Q_DSWsections_TSdiag_depth_yr20.png"
plt.savefig(fig_path + name_fig, dpi=300)


# plot ts diagram
fig = plt.figure(figsize=(20, 10))

plt.subplot(2,4,1)
plt.gcf().subplots_adjust(bottom=0.15)
print(salt_wed.shape, temp_wed.shape,z_rho_wed.shape)
for s, t, z in iris.iterate.izip(salt_wed, temp_wed, z_rho_wed, coords='S-coordinate at RHO-points'):
    sr_tmp = z.coord('S-coordinate at RHO-points').points
    dep_tmp = z.coord('bathymetry at RHO-points').points
    scat = iplt.scatter(s,t, c=-sr_tmp*dep_tmp, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=0, vmax=4000)
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
for s, t, z in iris.iterate.izip(salt_neu, temp_neu, z_rho_neu, coords='S-coordinate at RHO-points'):
    sr_tmp = z.coord('S-coordinate at RHO-points').points
    dep_tmp = z.coord('bathymetry at RHO-points').points
    scat = iplt.scatter(s,t, c=-sr_tmp*dep_tmp, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=0, vmax=4000)
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
for s, t, z in iris.iterate.izip(salt_mai, temp_mai, z_rho_mai, coords='S-coordinate at RHO-points'):
    sr_tmp = z.coord('S-coordinate at RHO-points').points
    dep_tmp = z.coord('bathymetry at RHO-points').points
    scat = iplt.scatter(s,t, c=-sr_tmp*dep_tmp, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=0, vmax=4000)
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
for s, t, z in iris.iterate.izip(salt_dar, temp_dar, z_rho_dar, coords='S-coordinate at RHO-points'):
    sr_tmp = z.coord('S-coordinate at RHO-points').points
    dep_tmp = z.coord('bathymetry at RHO-points').points
    scat = iplt.scatter(s,t, c=-sr_tmp*dep_tmp, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=0, vmax=4000)
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
for s, t, z in iris.iterate.izip(salt_mer, temp_mer, z_rho_mer, coords='S-coordinate at RHO-points'):
    sr_tmp = z.coord('S-coordinate at RHO-points').points
    dep_tmp = z.coord('bathymetry at RHO-points').points
    scat = iplt.scatter(s,t, c=-sr_tmp*dep_tmp, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=0, vmax=4000)
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
for s, t, z in iris.iterate.izip(salt_ros, temp_ros, z_rho_ros, coords='S-coordinate at RHO-points'):
    sr_tmp = z.coord('S-coordinate at RHO-points').points
    dep_tmp = z.coord('bathymetry at RHO-points').points
    scat = iplt.scatter(s,t, c=-sr_tmp*dep_tmp, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=0, vmax=4000)
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
for s, t, z in iris.iterate.izip(salt_amu, temp_amu, z_rho_amu, coords='S-coordinate at RHO-points'):
    sr_tmp = z.coord('S-coordinate at RHO-points').points
    dep_tmp = z.coord('bathymetry at RHO-points').points
    scat = iplt.scatter(s,t, c=-sr_tmp*dep_tmp, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=0, vmax=4000)
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

name_fig="waom4extend_shflim_S_0.25Q_AllDSWsections_TSdiag_depth_yr20.png"
plt.savefig(fig_path + name_fig, dpi=300)

# fig cross-sections + ts digrams DSW regions

# Plot transects
levelsT = np.arange(-2.5,2.5,.1)
levelsTf = np.arange(-2.5,2.5,.5)
levelsS = np.arange(33.5,34.9,.025)
levelsSf = np.arange(33.5,34.9,.05)
levelsR = np.arange(36.,37.4,.05)
levelsRf = np.arange(36.,37.4,.05)

#f, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 7))
fig = plt.figure(figsize=(20,15))
ax1 = fig.add_subplot(4,4,1)
ct = plt.contourf(lat_wed_mask, z_wed_mask, xr.DataArray.from_iris(temp_wed), levels=levelsT, cmap=plt.cm.coolwarm)
# plt.colorbar(ct, extend='both')
plt.contour(lat_wed_mask, z_wed_mask, xr.DataArray.from_iris(temp_wed), levels=levelsTf, colors='k', linewidths=.5)
plt.xlim([-78,-71])
plt.ylim([-1500,0])
plt.title('Weddell Sea')
plt.ylabel('Depth (m)')
# plt.xlabel('Latitude')

ax2 = fig.add_subplot(4,4,2)
#salt_wed.plot(levels=levelsS)
cs = plt.contourf(lat_wed_mask, z_wed_mask, xr.DataArray.from_iris(salt_wed), levels=levelsS, cmap=plt.cm.coolwarm)
# plt.colorbar(cs, extend='both')
plt.contour(lat_wed_mask, z_wed_mask, xr.DataArray.from_iris(salt_wed), levels=levelsSf, colors='k', linewidths=.5)
plt.xlim([-78,-71])
plt.ylim([-1500,0])

ax3 = fig.add_subplot(4,4,3)
ct = plt.contourf(lat_wed_mask, z_wed_mask, rho_wed, levels=levelsR, cmap=plt.cm.coolwarm)
# plt.colorbar(ct, extend='both')
plt.contour(lat_wed_mask, z_wed_mask, rho_wed, levels=levelsRf, colors='k', linewidths=.5)
plt.xlim([-78,-71])
plt.ylim([-1500,0])

ax5 = fig.add_subplot(4,4,5)
ct = plt.contourf(lat_ros_mask, z_ros_mask, xr.DataArray.from_iris(temp_ros), levels=levelsT, cmap=plt.cm.coolwarm)
# plt.colorbar(ct, extend='both')
plt.contour(lat_ros_mask, z_ros_mask, xr.DataArray.from_iris(temp_ros), levels=levelsTf, colors='k', linewidths=.5)
plt.xlim([-84,-70])
plt.ylim([-1500,0])
plt.title('Ross Sea')
plt.ylabel('Depth (m)')

ax6 = fig.add_subplot(4,4,6)
cs = plt.contourf(lat_ros_mask, z_ros_mask, xr.DataArray.from_iris(salt_ros), levels=levelsS, cmap=plt.cm.coolwarm)
# plt.colorbar(cs, extend='both')
plt.contour(lat_ros_mask, z_ros_mask, xr.DataArray.from_iris(salt_ros), levels=levelsSf, colors='k', linewidths=.5)
plt.xlim([-84,-70])
plt.ylim([-1500,0])

ax7 = fig.add_subplot(4,4,7)
ct = plt.contourf(lat_ros_mask, z_ros_mask, rho_ros, levels=levelsR, cmap=plt.cm.coolwarm)
# plt.colorbar(ct, extend='both')
plt.contour(lat_ros_mask, z_ros_mask, rho_ros, levels=levelsRf, colors='k', linewidths=.5)
plt.xlim([-84,-70])
plt.ylim([-1500,0])

ax9 = fig.add_subplot(4,4,9)
ct = plt.contourf(lat_dar_mask, z_dar_mask, xr.DataArray.from_iris(temp_dar), levels=levelsT, cmap=plt.cm.coolwarm)
# plt.colorbar(ct, extend='both')
plt.contour(lat_dar_mask, z_dar_mask, xr.DataArray.from_iris(temp_dar), levels=levelsTf, colors='k', linewidths=.5)
plt.xlim([-70.5,-65])
plt.ylim([-1500,0])
plt.title('Darnley Polynya')
plt.ylabel('Depth (m)') # minus is height not depth!

ax10 = fig.add_subplot(4,4,10)
cs = plt.contourf(lat_dar_mask, z_dar_mask, xr.DataArray.from_iris(salt_dar), levels=levelsS, cmap=plt.cm.coolwarm)
# plt.colorbar(cs, extend='both')
plt.contour(lat_dar_mask, z_dar_mask, xr.DataArray.from_iris(salt_dar), levels=levelsSf, colors='k', linewidths=.5)
plt.xlim([-70.5,-65])
plt.ylim([-1500,0])

ax11 = fig.add_subplot(4,4,11)
ct = plt.contourf(lat_dar_mask, z_dar_mask, rho_dar, levels=levelsR, cmap=plt.cm.coolwarm)
# plt.colorbar(ct, extend='both')
plt.contour(lat_dar_mask, z_dar_mask, rho_dar, levels=levelsRf, colors='k', linewidths=.5)
plt.xlim([-70.5,-65])
plt.ylim([-1500,0])

ax13 = fig.add_subplot(4,4,13)
ct = plt.contourf(lat_mer_mask, z_mer_mask, xr.DataArray.from_iris(temp_mer), levels=levelsT, cmap=plt.cm.coolwarm)
plt.contour(lat_mer_mask, z_mer_mask, xr.DataArray.from_iris(temp_mer), levels=levelsTf, colors='k', linewidths=.5)
plt.xlim([-68.,-65])
plt.ylim([-1500,0])
plt.title('Merz Glacier')
plt.ylabel('Depth (m)') # minus is height not depth!
plt.xlabel('Latitude')
cbar_ax1 = fig.add_axes([0.125, 0.06, 0.17, 0.01])
cbt = fig.colorbar(ct, cax=cbar_ax1, orientation='horizontal')#plt.colorbar(ct, extend='both')
cbt.ax.set_xlabel('Temperature (degC)')

ax14 = fig.add_subplot(4,4,14)
cs = plt.contourf(lat_mer_mask, z_mer_mask, xr.DataArray.from_iris(salt_mer), levels=levelsS, cmap=plt.cm.coolwarm)
plt.contour(lat_mer_mask, z_mer_mask, xr.DataArray.from_iris(salt_mer), levels=levelsSf, colors='k', linewidths=.5)
plt.xlim([-68.,-65])
plt.ylim([-1500,0])
plt.xlabel('Latitude')
cbar_ax2 = fig.add_axes([0.33, 0.06, 0.17, 0.015])
cbs = fig.colorbar(cs, cax=cbar_ax2, orientation='horizontal')
cbs.ax.set_xlabel('Salinity')

ax15 = fig.add_subplot(4,4,15)
cr = plt.contourf(lat_mer_mask, z_mer_mask, rho_mer, levels=levelsR, cmap=plt.cm.coolwarm)
plt.contour(lat_mer_mask, z_mer_mask, rho_mer, levels=levelsRf, colors='k', linewidths=.5)
plt.xlim([-68.,-65])
plt.ylim([-1500,0])
plt.xlabel('Latitude')

cbar_ax3 = fig.add_axes([0.53, 0.06, 0.17, 0.015])
cbr = fig.colorbar(cr, cax=cbar_ax3, orientation='horizontal')
cbr.ax.set_xlabel('Potential density ($\sigma_2$, kg m$^{-3}$)')

# ts-diagram

plt.subplot(4,4,4)
plt.gcf().subplots_adjust(bottom=0.15)

for s, t, z in iris.iterate.izip(salt_wed, temp_wed, z_rho_wed, coords='S-coordinate at RHO-points'):
    sr_tmp = z.coord('S-coordinate at RHO-points').points
    dep_tmp = z.coord('bathymetry at RHO-points').points
    scat = iplt.scatter(s,t, c=-sr_tmp*dep_tmp, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=0, vmax=4000)
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

for s, t, z in iris.iterate.izip(salt_ros, temp_ros, z_rho_ros, coords='S-coordinate at RHO-points'):
    sr_tmp = z.coord('S-coordinate at RHO-points').points
    dep_tmp = z.coord('bathymetry at RHO-points').points
    scat = iplt.scatter(s,t, c=-sr_tmp*dep_tmp, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=0, vmax=4000)
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

for s, t, z in iris.iterate.izip(salt_dar, temp_dar, z_rho_dar, coords='S-coordinate at RHO-points'):
    sr_tmp = z.coord('S-coordinate at RHO-points').points
    dep_tmp = z.coord('bathymetry at RHO-points').points
    scat = iplt.scatter(s,t, c=-sr_tmp*dep_tmp, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=0, vmax=4000)
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

for s, t, z in iris.iterate.izip(salt_mer, temp_mer, z_rho_mer, coords='S-coordinate at RHO-points'):
    sr_tmp = z.coord('S-coordinate at RHO-points').points
    dep_tmp = z.coord('bathymetry at RHO-points').points
    scat = iplt.scatter(s,t, c=-sr_tmp*dep_tmp, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=0, vmax=4000)
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
#plt.show()

name_fig="waom4extend_shflim_S_0.25Q_DSW+sections_TSdiag_depth_yr20.png"
plt.savefig(fig_path + name_fig, dpi=300)
