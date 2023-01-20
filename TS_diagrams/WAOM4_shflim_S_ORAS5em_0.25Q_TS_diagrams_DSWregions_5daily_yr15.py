# read nc output from WAOM 10km run

import xarray as xr
import pandas as p
import numpy as np
import numpy.ma as ma
import cartopy.crs as ccrs
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as dates
from matplotlib.ticker import NullFormatter
from matplotlib.colors import LinearSegmentedColormap   # for custom colormaps

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

# load ROMS avg output
for mm  in ['01','02','03','04','05','06','07','08','09','10','11','12']:
    ds = xr.open_dataset('//scratch/project_2000339/boeiradi/waom4_shflim_S_ORAS5em_0.25Q/output_5yr/ocean_avg_00' + mm + '.nc')
    print('MM: size temp and time length: ', mm, ds.temp.shape, len(ds.salt.isel(xi_rho=20, eta_rho=100, s_rho=0)))
    temp_tmp = ds.variables["temp"]
    salt_tmp = ds.variables["salt"]
    temp_mavg = np.nanmean(temp_tmp, axis=0) # monthly average
    salt_mavg = np.nanmean(salt_tmp, axis=0)
    del temp_tmp, salt_tmp

    ds = ds.set_coords(['Cs_r', 'Cs_w', 'hc', 'h', 'Vtransform'])

    if ds.Vtransform == 1:
        Zo_rho = ds.hc * (ds.s_rho - ds.Cs_r) + ds.Cs_r * ds.h
        z_rho_tmp = Zo_rho + ds.zeta * (1 + Zo_rho/ds.h)
        print("Vtransform=1")
    elif ds.Vtransform == 2:
        Zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * ds.h) / (ds.hc + ds.h)
        z_rho_tmp = ds.zeta + (ds.zeta + ds.h) * Zo_rho
        print("Vtransform=2")
    z_rho_mavg = np.nanmean(z_rho_tmp, axis=0)
    del z_rho_tmp

    # concatenate monthly averages into annual vars:
    if mm == '01':
        temp_ann = temp_mavg
        salt_ann = salt_mavg
        z_rho_ann = z_rho_mavg
    elif mm == '02':
        temp_ann = np.stack((temp_ann,temp_mavg), axis=0)
        salt_ann = np.stack((salt_ann,salt_mavg), axis=0)
        z_rho_ann = np.stack((z_rho_ann,z_rho_mavg), axis=0)

    else:
        temp_mavg_4thdim = np.expand_dims(temp_mavg, axis=0)
        temp_ann = np.concatenate((temp_ann,temp_mavg_4thdim), axis=0)
        salt_mavg_4thdim = np.expand_dims(salt_mavg, axis=0)
        salt_ann = np.concatenate((salt_ann,salt_mavg_4thdim), axis=0)
        z_rho_mavg_4thdim = np.expand_dims(z_rho_mavg, axis=0)
        z_rho_ann = np.concatenate((z_rho_ann,z_rho_mavg_4thdim), axis=0)


dg = xr.open_dataset("/scratch/project_2000339/boeiradi/waom4_frc/waom4_grd.nc")
lat_rho = dg.variables["lat_rho"]
lon_rho = dg.variables["lon_rho"]
ds.coords['lat_rho']=lat_rho.transpose() # put lat_rho into ds dataset
ds.coords['lon_rho']=lon_rho.transpose() # put lon_rho into ds dataset

sigma_t_ann = gsw.rho(salt_ann,temp_ann,2000) - 1000
# mask out values ~1e37

# save temp,salt sections
ncdftmp_path = '/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/'

variables = ['temp','salt','zeta','z_rho']
ds[variables].isel(ocean_time=slice(0, 12), #wed
                  xi_rho=425, eta_rho=slice(700,1125)).to_netcdf(ncdftmp_path + 'WAOM4_shflim_S_ORAS5_0.25Q_sec_wed_ann_yr15.nc', mode='w')

ds[variables].isel(ocean_time=slice(0, 12), # ros
                  xi_rho=775, eta_rho=slice(50,625)).to_netcdf(ncdftmp_path + 'WAOM4_shflim_S_ORAS5_0.25Q_sec_ros_ann_yr15.nc', mode='w')

ds[variables].isel(ocean_time=slice(0, 12), # pry
                  xi_rho=slice(1050,1500), eta_rho=845).to_netcdf(ncdftmp_path + 'WAOM4_shflim_S_ORAS5_0.25Q_sec_pry_ann_yr15.nc', mode='w')

ds[variables].isel(ocean_time=slice(0, 12), # neu
                  xi_rho=700, eta_rho=slice(1050,1300)).to_netcdf(ncdftmp_path + 'WAOM4_shflim_S_ORAS5_0.25Q_sec_neu_ann_yr15.nc', mode='w')

ds[variables].isel(ocean_time=slice(0, 12), # mai
                  xi_rho=800, eta_rho=slice(1050,1300)).to_netcdf(ncdftmp_path + 'WAOM4_shflim_S_ORAS5_0.25Q_sec_mai_ann_yr15.nc', mode='w')

ds[variables].isel(ocean_time=slice(0, 12), # dar
                  eta_rho=875, xi_rho=slice(1150,1538)).to_netcdf(ncdftmp_path + 'WAOM4_shflim_S_ORAS5_0.25Q_sec_dar_ann_yr15.nc', mode='w')

ds[variables].isel(ocean_time=slice(0, 12), # mer
                  xi_rho=1088, eta_rho=slice(25,200)).to_netcdf(ncdftmp_path + 'WAOM4_shflim_S_ORAS5_0.25Q_sec_mer_ann_yr15.nc', mode='w')

ds[variables].isel(ocean_time=slice(0, 12), # amu
                  eta_rho=488, xi_rho=slice(125,500)).to_netcdf(ncdftmp_path + 'WAOM4_shflim_S_ORAS5_0.25Q_sec_amu_ann_yr15.nc', mode='w')

print('saved')

# load vars for TS plot
wed=ncdftmp_path + 'WAOM4_shflim_S_ORAS5_0.25Q_sec_wed_ann_yr15.nc'
temp_wed = iris.load_cube(wed, 'time-averaged potential temperature')
salt_wed = iris.load_cube(wed, 'time-averaged salinity')
z_rho_wed = iris.load_cube(wed, 'z_rho')

ros=ncdftmp_path + 'WAOM4_shflim_S_ORAS5_0.25Q_sec_ros_ann_yr15.nc'
temp_ros = iris.load_cube(ros, 'time-averaged potential temperature')
salt_ros = iris.load_cube(ros, 'time-averaged salinity')
z_rho_ros = iris.load_cube(ros, 'z_rho')

pry=ncdftmp_path + 'WAOM4_shflim_S_ORAS5_0.25Q_sec_pry_ann_yr15.nc'
temp_pry = iris.load_cube(pry, 'time-averaged potential temperature')
salt_pry = iris.load_cube(pry, 'time-averaged salinity')
z_rho_pry = iris.load_cube(pry, 'z_rho')

neu=ncdftmp_path + 'WAOM4_shflim_S_ORAS5_0.25Q_sec_neu_ann_yr15.nc'
temp_neu = iris.load_cube(neu, 'time-averaged potential temperature')
salt_neu = iris.load_cube(neu, 'time-averaged salinity')
z_rho_neu = iris.load_cube(neu, 'z_rho')

mai=ncdftmp_path + 'WAOM4_shflim_S_ORAS5_0.25Q_sec_mai_ann_yr15.nc'
temp_mai = iris.load_cube(mai, 'time-averaged potential temperature')
salt_mai = iris.load_cube(mai, 'time-averaged salinity')
z_rho_mai = iris.load_cube(mai, 'z_rho')

mer=ncdftmp_path + 'WAOM4_shflim_S_ORAS5_0.25Q_sec_mer_ann_yr15.nc'
temp_mer = iris.load_cube(mer, 'time-averaged potential temperature')
salt_mer = iris.load_cube(mer, 'time-averaged salinity')
z_rho_mer = iris.load_cube(mer, 'z_rho')

dar=ncdftmp_path + 'WAOM4_shflim_S_ORAS5_0.25Q_sec_dar_ann_yr15.nc'
temp_dar = iris.load_cube(dar, 'time-averaged potential temperature')
salt_dar = iris.load_cube(dar, 'time-averaged salinity')
z_rho_dar = iris.load_cube(dar, 'z_rho')

amu=ncdftmp_path + 'WAOM4_shflim_S_ORAS5_0.25Q_sec_amu_ann_yr15.nc'
temp_amu = iris.load_cube(amu, 'time-averaged potential temperature')
salt_amu = iris.load_cube(amu, 'time-averaged salinity')
z_rho_amu = iris.load_cube(amu, 'z_rho')

print('loaded')

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

#s_rho_wed = z_rho_wed.coord('S-coordinate at RHO-points').points
#s_rho_ros = z_rho_ros.coord('S-coordinate at RHO-points').points
#s_rho_pry = z_rho_pry.coord('S-coordinate at RHO-points').points
#s_rho_neu = z_rho_neu.coord('S-coordinate at RHO-points').points
#s_rho_mai = z_rho_mai.coord('S-coordinate at RHO-points').points
#s_rho_mer = z_rho_mer.coord('S-coordinate at RHO-points').points
#s_rho_dar = z_rho_dar.coord('S-coordinate at RHO-points').points
#s_rho_amu = z_rho_amu.coord('S-coordinate at RHO-points').points

depth_wed = z_rho_wed.coord('bathymetry at RHO-points').points
depth_ros = z_rho_ros.coord('bathymetry at RHO-points').points
depth_pry = z_rho_pry.coord('bathymetry at RHO-points').points
depth_neu = z_rho_neu.coord('bathymetry at RHO-points').points
depth_mai = z_rho_mai.coord('bathymetry at RHO-points').points
depth_mer = z_rho_mer.coord('bathymetry at RHO-points').points
depth_dar = z_rho_dar.coord('bathymetry at RHO-points').points
depth_amu = z_rho_amu.coord('bathymetry at RHO-points').points

#print(z_rho_wed)

print(temp_wed.shape, salt_wed.shape)

fig_path='/users/boeiradi/COLD_project/postprocessing/figs/TS_diagrams'


# plot ts diagram: All sections
fig = plt.figure(figsize=(20, 10))

plt.subplot(2,4,1)
plt.gcf().subplots_adjust(bottom=0.15)
print(salt_wed.shape, temp_wed.shape,z_rho_wed.shape)
for s, t in iris.iterate.izip(salt_wed, temp_wed, coords='bathymetry at RHO-points'):
    scat = iplt.scatter(s,t, c=depth_wed, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=200, vmax=4000)
ax = plt.gca()
ax.set_xlim([33.,35.])
ax.set_ylabel('Temperature ($^{\circ}$C)', fontsize=16)
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
ax.set_xlabel('Salinity', fontsize=16)
ax.set_xlim([33.,35.])
ax.set_xticks([33, 33.5, 34, 34.5, 35])
ax.set_ylim([-3,5])
plt.title('Amundsen Sea', fontsize=14)
ax.tick_params(labelsize=16)
CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(33.5,38,.2),linestyles='dashed', colors='k')
cbw = plt.contour(Si,Ti,dens.transpose(), levels=[36.9],linestyles='dashed', colors='m',linewidth=2)
plt.clabel(CS, fontsize=12, inline=1, fmt='%1.1f')
#name_fig="waom10_amudzBay_TSdiag_depth.png"
#plt.savefig(fig_path + name_fig)

cbar_axim = fig.add_axes([0.92, 0.15, 0.015, 0.7])
cb = fig.colorbar(scat, cax=cbar_axim)
cb.ax.tick_params(labelsize=12)
cb.set_label('Depth (m)',fontsize=14)

name_fig="waom10_shflim_S_ORAS5em_0.25Q_AllDSWsections_TSdiag_depth_yr15.png"
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

name_fig="waom10_shflim_S_ORAS5_0.25Q_DSW+sections_TSdiag_depth_yr15.png"
plt.savefig(fig_path + name_fig, dpi=300)



