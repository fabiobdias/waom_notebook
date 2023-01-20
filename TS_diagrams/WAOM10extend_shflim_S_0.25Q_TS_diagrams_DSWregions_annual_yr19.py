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
ds = xr.open_dataset("/scratch/gh9/fbd581/ROMS/waom10extend_shflim_S_0.25Q/output_01-20yr/ocean_avg_0019.nc")

#print(output.variables.keys()) # get all variable names

temp = ds.variables["temp"]
salt = ds.variables["salt"]
print(temp.shape)
m_len = temp[:,0,0,0].size
k_len = temp[0,:,0,0].size
i_len = temp[0,0,:,0].size
j_len = temp[0,0,0,:].size

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

dg = xr.open_dataset("/g/data3/hh5/tmp/access-om/fbd581/ROMS/waom10_frc/waom10extend_grd.nc")

lat_rho = dg.variables["lat_rho"]
lon_rho = dg.variables["lon_rho"]
ds.coords['lat_rho']=lat_rho.transpose() # put lat_rho into ds dataset
ds.coords['lon_rho']=lon_rho.transpose() # put lon_rho into ds dataset

sigma_t = gsw.rho(salt,temp,2000) - 1000
# mask out values ~1e37

# max,min temp,salt|
print(np.max(temp))
print(np.min(temp))
print(np.max(salt))
print(np.min(salt))


print(z_rho.shape)
print(ds.temp.shape)
print(lon_rho[310:480,200].shape)
print(np.ones((170,31)).shape)

## for new extended grid, change in eta_rho indices (+30):
# wed: eta_rho=slice(280,450)) -> eta_rho=slice(310,480))
# ros: eta_rho=slice(20,250)) -> eta_rho=slice(50,280))
# pry: eta_rho=338 -> eta_rho=368
# neu/mai: eta_rho=slice(420,520)) -> eta_rho=slice(450,550))
# mer: eta_rho=slice(10,80)) -> eta_rho=slice(40,110))
# dar:  eta_rho=350 ->  eta_rho=380
# amu: eta_rho=195 -> eta_rho=225


# plot some sections: Salt
# Weddell Sea
section_rho_wed = np.nanmean(np.squeeze(sigma_t[:,:,310:480,170]), axis=0)
section_salt_wed = np.nanmean(ds.salt.isel(xi_rho=170, eta_rho=slice(310,480)), axis=0)
section_temp_wed = np.nanmean(ds.temp.isel(xi_rho=170, eta_rho=slice(310,480)), axis=0)
section_lat_wed_tmp = ds.lat_rho.isel(xi_rho=170, eta_rho=slice(310,480))
section_z_wed = np.nanmean(ds.z_rho.isel(xi_rho=170, eta_rho=slice(310,480)), axis=2)
section_z_wed_mask = ma.array(section_z_wed,mask=np.isnan(section_z_wed))
section_lat_wed = np.ones((31,170))
for ii in np.arange(0,31):
    section_lat_wed[ii,:] = section_lat_wed_tmp
section_lat_wed_mask = ma.array(section_lat_wed,mask=np.isnan(section_z_wed))

print(section_lat_wed.shape, section_z_wed.shape, section_salt_wed.shape)

# Ross Sea
section_rho_ros = np.nanmean(np.squeeze(sigma_t[:,:,50:280,310]), axis=0)
section_salt_ros = np.nanmean(ds.salt.isel(xi_rho=310, eta_rho=slice(50,280)), axis=0)
section_temp_ros = np.nanmean(ds.temp.isel(xi_rho=310, eta_rho=slice(50,280)), axis=0)
section_lat_ros_tmp = ds.lat_rho.isel(xi_rho=310, eta_rho=slice(50,280))
section_z_ros = np.nanmean(ds.z_rho.isel(xi_rho=310, eta_rho=slice(50,280)), axis=2)
section_z_ros_mask = ma.array(section_z_ros,mask=np.isnan(section_z_ros))
section_lat_ros = np.ones((31,230))
for ii in np.arange(0,31):
    section_lat_ros[ii,:] = section_lat_ros_tmp
section_lat_ros_mask = ma.array(section_lat_ros,mask=np.isnan(section_z_ros))

# Prydz Bay
section_rho_pry = np.nanmean(np.squeeze(sigma_t[:,:,368,420:600]), axis=0)
section_salt_pry = np.nanmean(ds.salt.isel(xi_rho=slice(420,600), eta_rho=368), axis=0)
section_temp_pry = np.nanmean(ds.temp.isel(xi_rho=slice(420,600), eta_rho=368), axis=0)
section_lon_pry_tmp = ds.lon_rho.isel(xi_rho=slice(420,600), eta_rho=368)
section_z_pry = np.nanmean(ds.z_rho.isel(xi_rho=slice(420,600), eta_rho=368), axis=2)
section_z_pry_mask = ma.array(section_z_pry,mask=np.isnan(section_z_pry))
section_lon_pry = np.ones((31,180))
for ii in np.arange(0,31):
    section_lon_pry[ii,:] = section_lon_pry_tmp
section_lon_pry_mask = ma.array(section_lon_pry,mask=np.isnan(section_z_pry))

# plot some sections: Salt
# Maud Land1
section_rho_neu = np.nanmean(np.squeeze(sigma_t[:,:,450:550,280]), axis=0)
section_salt_neu = np.nanmean(ds.salt.isel(xi_rho=280, eta_rho=slice(450,550)), axis=0)
section_temp_neu = np.nanmean(ds.temp.isel(xi_rho=280, eta_rho=slice(450,550)), axis=0)
section_lat_neu_tmp = ds.lat_rho.isel(xi_rho=280, eta_rho=slice(450,550))
section_z_neu = np.nanmean(ds.z_rho.isel(xi_rho=280, eta_rho=slice(450,550)), axis=2)
section_z_neu_mask = ma.array(section_z_neu,mask=np.isnan(section_z_neu))
section_lat_neu = np.ones((31,100))
for ii in np.arange(0,31):
    section_lat_neu[ii,:] = section_lat_neu_tmp
section_lat_neu_mask = ma.array(section_lat_neu,mask=np.isnan(section_z_neu))

print(section_lat_neu.shape, section_z_neu.shape, section_salt_neu.shape)

section_rho_mai = np.nanmean(np.squeeze(sigma_t[:,:,450:550,320]), axis=0)
section_salt_mai = np.nanmean(ds.salt.isel(xi_rho=320, eta_rho=slice(450,550)), axis=0)
section_temp_mai = np.nanmean(ds.temp.isel(xi_rho=320, eta_rho=slice(450,550)), axis=0)
section_lat_mai_tmp = ds.lat_rho.isel(xi_rho=320, eta_rho=slice(450,550))
section_z_mai = np.nanmean(ds.z_rho.isel(xi_rho=320, eta_rho=slice(450,550)), axis=2)
section_z_mai_mask = ma.array(section_z_mai,mask=np.isnan(section_z_mai))
section_lat_mai = np.ones((31,100))
for ii in np.arange(0,31):
    section_lat_mai[ii,:] = section_lat_mai_tmp
section_lat_mai_mask = ma.array(section_lat_mai,mask=np.isnan(section_z_mai))

print(section_lat_mai.shape, section_z_mai.shape, section_salt_mai.shape)

# Merz Glacier
section_rho_mer = np.nanmean(np.squeeze(sigma_t[:,:,40:110,435]), axis=0)
section_salt_mer = np.nanmean(ds.salt.isel(xi_rho=435, eta_rho=slice(40,110)), axis=0)
section_temp_mer = np.nanmean(ds.temp.isel(xi_rho=435, eta_rho=slice(40,110)), axis=0)
section_lat_mer_tmp = ds.lat_rho.isel(xi_rho=435, eta_rho=slice(40,110))
section_z_mer = np.nanmean(ds.z_rho.isel(xi_rho=435, eta_rho=slice(40,110)), axis=2)
section_z_mer_mask = ma.array(section_z_mer,mask=np.isnan(section_z_mer))
section_lat_mer = np.ones((31,70))
for ii in np.arange(0,31):
    section_lat_mer[ii,:] = section_lat_mer_tmp
section_lat_mer_mask = ma.array(section_lat_mer,mask=np.isnan(section_z_mer))

# Darnley Polynya
section_rho_dar = np.nanmean(np.squeeze(sigma_t[:,:,380,460:615]), axis=0)
section_salt_dar = np.nanmean(ds.salt.isel(xi_rho=slice(460,615), eta_rho=380), axis=0)
section_temp_dar = np.nanmean(ds.temp.isel(xi_rho=slice(460,615), eta_rho=380), axis=0)
section_lat_dar_tmp = ds.lat_rho.isel(xi_rho=slice(460,615), eta_rho=380)
section_z_dar = np.nanmean(ds.z_rho.isel(xi_rho=slice(460,615), eta_rho=380), axis=2)
section_z_dar_mask = ma.array(section_z_dar,mask=np.isnan(section_z_dar))
section_lat_dar = np.ones((31,155))
for ii in np.arange(0,31):
    section_lat_dar[ii,:] = section_lat_dar_tmp
section_lat_dar_mask = ma.array(section_lat_dar,mask=np.isnan(section_z_dar))

# Amundsen Sea
section_rho_amu = np.nanmean(np.squeeze(sigma_t[:,:,50:200,225]), axis=0)
section_salt_amu = np.nanmean(ds.salt.isel(eta_rho=225, xi_rho=slice(50,200)), axis=0)
section_temp_amu = np.nanmean(ds.temp.isel(eta_rho=225, xi_rho=slice(50,200)), axis=0)
section_lat_amu_tmp = ds.lat_rho.isel(eta_rho=225, xi_rho=slice(50,200))
section_z_amu = np.nanmean(ds.z_rho.isel(eta_rho=225, xi_rho=slice(50,200)), axis=2)
section_z_amu_mask = ma.array(section_z_amu,mask=np.isnan(section_z_amu))
section_lat_amu = np.ones((31,150))
for ii in np.arange(0,31):
    section_lat_amu[ii,:] = section_lat_amu_tmp
section_lat_amu_mask = ma.array(section_lat_amu,mask=np.isnan(section_z_amu))

# save temp,salt sections
ncdftmp_path = '/g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/ncdf_tmp/'

variables = ['temp','salt','zeta','z_rho']
ds[variables].isel(ocean_time=slice(0, 12), #wed
                  xi_rho=170, eta_rho=slice(310,480)).to_netcdf(ncdftmp_path + 'WAOM10extend_shflim_S_0.25Q_sec_wed_ann_yr19.nc', mode='w')

ds[variables].isel(ocean_time=slice(0, 12), # ros
                  xi_rho=310, eta_rho=slice(50,280)).to_netcdf(ncdftmp_path + 'WAOM10extend_shflim_S_0.25Q_sec_ros_ann_yr19.nc', mode='w')

ds[variables].isel(ocean_time=slice(0, 12), # pry
                  xi_rho=slice(420,600), eta_rho=368).to_netcdf(ncdftmp_path + 'WAOM10extend_shflim_S_0.25Q_sec_pry_ann_yr19.nc', mode='w')

ds[variables].isel(ocean_time=slice(0, 12), # neu
                  xi_rho=280, eta_rho=slice(450,550)).to_netcdf(ncdftmp_path + 'WAOM10extend_shflim_S_0.25Q_sec_neu_ann_yr19.nc', mode='w')

ds[variables].isel(ocean_time=slice(0, 12), # mai
                  xi_rho=320, eta_rho=slice(450,550)).to_netcdf(ncdftmp_path + 'WAOM10extend_shflim_S_0.25Q_sec_mai_ann_yr19.nc', mode='w')

ds[variables].isel(ocean_time=slice(0, 12), # dar
                  eta_rho=380, xi_rho=slice(460,615)).to_netcdf(ncdftmp_path + 'WAOM10extend_shflim_S_0.25Q_sec_dar_ann_yr19.nc', mode='w')

ds[variables].isel(ocean_time=slice(0, 12), # mer
                  xi_rho=435, eta_rho=slice(40,110)).to_netcdf(ncdftmp_path + 'WAOM10extend_shflim_S_0.25Q_sec_mer_ann_yr19.nc', mode='w')

ds[variables].isel(ocean_time=slice(0, 12), # amu
                  eta_rho=225, xi_rho=slice(50,200)).to_netcdf(ncdftmp_path + 'WAOM10extend_shflim_S_0.25Q_sec_amu_ann_yr19.nc', mode='w')

print('saved')

# load vars for TS plot
wed=ncdftmp_path + 'WAOM10extend_shflim_S_0.25Q_sec_wed_ann_yr19.nc'
temp_wed = iris.load_cube(wed, 'time-averaged potential temperature')
salt_wed = iris.load_cube(wed, 'time-averaged salinity')
z_rho_wed = iris.load_cube(wed, 'z_rho')

ros=ncdftmp_path + 'WAOM10extend_shflim_S_0.25Q_sec_ros_ann_yr19.nc'
temp_ros = iris.load_cube(ros, 'time-averaged potential temperature')
salt_ros = iris.load_cube(ros, 'time-averaged salinity')
z_rho_ros = iris.load_cube(ros, 'z_rho')

pry=ncdftmp_path + 'WAOM10extend_shflim_S_0.25Q_sec_pry_ann_yr19.nc'
temp_pry = iris.load_cube(pry, 'time-averaged potential temperature')
salt_pry = iris.load_cube(pry, 'time-averaged salinity')
z_rho_pry = iris.load_cube(pry, 'z_rho')

neu=ncdftmp_path + 'WAOM10extend_shflim_S_0.25Q_sec_neu_ann_yr19.nc'
temp_neu = iris.load_cube(neu, 'time-averaged potential temperature')
salt_neu = iris.load_cube(neu, 'time-averaged salinity')
z_rho_neu = iris.load_cube(neu, 'z_rho')

mai=ncdftmp_path + 'WAOM10extend_shflim_S_0.25Q_sec_mai_ann_yr19.nc'
temp_mai = iris.load_cube(mai, 'time-averaged potential temperature')
salt_mai = iris.load_cube(mai, 'time-averaged salinity')
z_rho_mai = iris.load_cube(mai, 'z_rho')

mer=ncdftmp_path + 'WAOM10extend_shflim_S_0.25Q_sec_mer_ann_yr19.nc'
temp_mer = iris.load_cube(mer, 'time-averaged potential temperature')
salt_mer = iris.load_cube(mer, 'time-averaged salinity')
z_rho_mer = iris.load_cube(mer, 'z_rho')

dar=ncdftmp_path + 'WAOM10extend_shflim_S_0.25Q_sec_dar_ann_yr19.nc'
temp_dar = iris.load_cube(dar, 'time-averaged potential temperature')
salt_dar = iris.load_cube(dar, 'time-averaged salinity')
z_rho_dar = iris.load_cube(dar, 'z_rho')

amu=ncdftmp_path + 'WAOM10extend_shflim_S_0.25Q_sec_amu_ann_yr19.nc'
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

fig_path='/g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/figs/TS_diagrams/'


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

name_fig="waom10extend_shflim_S_0.25Q_AllDSWsections_TSdiag_depth_yr19.png"
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

name_fig="waom10extend_shflim_S_0.25Q_DSW+sections_TSdiag_depth_yr19.png"
plt.savefig(fig_path + name_fig, dpi=300)



