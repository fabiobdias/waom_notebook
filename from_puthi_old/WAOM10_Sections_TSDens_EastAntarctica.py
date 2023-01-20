# plot sections monthly:

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
ds = xr.open_dataset("/scratch/project_2000339/boeiradi/waom10_diag/output_17_float/ocean_avg_0017.nc")
#print(output.variables.keys()) # get all variable names

temp = ds.variables["temp"]
salt = ds.variables["salt"]

dens=gsw.rho(salt,temp,2000)
        # Substract 1000 to convert to sigma-t
dens = dens - 1000
print(temp.shape,dens.shape)
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

dg = xr.open_dataset("/scratch/project_2000339/boeiradi/waom10_frc/waom10_grd.nc")

lat_rho = dg.variables["lat_rho"]
lon_rho = dg.variables["lon_rho"]
ds.coords['lat_rho']=lat_rho.transpose() # put lat_rho into ds dataset
ds.coords['lon_rho']=lon_rho.transpose() # put lon_rho into ds dataset

# Period for storms in the winter of 2007:
events_ind = np.array([31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 53, 54])
date_storms = [('03.Jun.2007'),('08.Jun.2007'),('13.Jun.2007'),('18.Jun.2007'),('23.Jun.2007'),('28.Jun.2007'), \
               ('03.Jul.2007'),('08.Jul.2007'),('13.Jul.2007'),('18.Jul.2007'),('23.Jul.2007'),('28.Jul.2007'), \
               ('02.Aug.2007'),('07.Aug.2007'),('12.Aug.2007'),('17.Aug.2007'),('22.Aug.2007'), \
               ('21.Sep.2007'),('26.Sep.2007')]

# Darnley Polynya
# 1 - slice section for salt/temp/z_rho
section_salt_dar = np.empty((19,31,155)) # allocating
section_temp_dar = np.empty((19,31,155))
section_rho_dar = np.empty((19,31,155))
section_z_dar = np.empty((19,31,155))
section_z_dar_mask = np.empty((19,31,155))
for ee in np.arange(0,19):
    section_salt_dar[ee,:,:] = ds.salt.isel(xi_rho=slice(460,615), eta_rho=350, ocean_time=ee)
    section_temp_dar[ee,:,:] = ds.temp.isel(xi_rho=slice(460,615), eta_rho=350, ocean_time=ee)
    section_z_dar[ee,:,:] = ds.z_rho.isel(xi_rho=slice(460,615), eta_rho=350, ocean_time=ee)
    # 1.1 - mask land values in z_rho slice                                           
    section_z_dar_mask[...] = ma.array(section_z_dar[ee,:,:],mask=np.isnan(section_z_dar[ee,:,:]))
    section_rho_dar[ee,:,:] = np.squeeze(dens[ee,:,350,460:615])
    
# 2 - slide section for lat or lat
section_lat_dar_tmp = ds.lat_rho.isel(xi_rho=slice(460,615), eta_rho=350)
# 2.1 - mask land values for lat/lat needs loop for repeat through vert layers
section_lat_dar = np.ones((31,155))                                                
for ii in np.arange(0,31):
    section_lat_dar[ii,:] = section_lat_dar_tmp
section_lat_dar_mask = ma.array(section_lat_dar,mask=np.isnan(section_z_dar[0,:,:]))

print(section_lat_dar.shape, section_z_dar.shape, section_salt_dar.shape)

section_z_dar_mask_ann = np.squeeze(np.nanmean(section_z_dar_mask[:,:,:], axis=0))
section_temp_dar_ann = np.squeeze(np.nanmean(section_temp_dar[:,:,:], axis=0))
section_salt_dar_ann = np.squeeze(np.nanmean(section_salt_dar[:,:,:], axis=0))
section_rho_dar_ann = np.squeeze(np.nanmean(section_rho_dar[:,:,:], axis=0))

# Prydz Bay (Mackenzi Polynya)
# 1 - slice section for salt/temp/z_rho
section_salt_pry = np.empty((19,31,155)) # allocating
section_temp_pry = np.empty((19,31,155))
section_rho_pry = np.empty((19,31,155))
section_z_pry = np.empty((19,31,155))
section_z_pry_mask = np.empty((19,31,155))
for ee in np.arange(0,19):
    section_salt_pry[ee,:,:] = ds.salt.isel(xi_rho=slice(460,615), eta_rho=338, ocean_time=ee)
    section_temp_pry[ee,:,:] = ds.temp.isel(xi_rho=slice(460,615), eta_rho=338, ocean_time=ee)
    section_z_pry[ee,:,:] = ds.z_rho.isel(xi_rho=slice(460,615), eta_rho=338, ocean_time=ee)
    # 1.1 - mask land values in z_rho slice                                           
    section_z_pry_mask[...] = ma.array(section_z_pry[ee,:,:],mask=np.isnan(section_z_pry[ee,:,:]))
    section_rho_pry[ee,:,:] = np.squeeze(dens[ee,:,338,460:615])

# 2 - slide section for lat or lat
section_lat_pry_tmp = ds.lat_rho.isel(xi_rho=slice(460,615), eta_rho=338)
# 2.1 - mask land values for lat/lat needs loop for repeat through vert layers
section_lat_pry = np.ones((31,155))                                                
for ii in np.arange(0,31):
    section_lat_pry[ii,:] = section_lat_pry_tmp
section_lat_pry_mask = ma.array(section_lat_pry,mask=np.isnan(section_z_pry[0,:,:]))

print(section_lat_pry.shape, section_z_pry.shape, section_salt_pry.shape)

section_z_pry_mask_ann = np.squeeze(np.nanmean(section_z_pry_mask[:,:,:], axis=0))
section_temp_pry_ann = np.squeeze(np.nanmean(section_temp_pry[:,:,:], axis=0))
section_salt_pry_ann = np.squeeze(np.nanmean(section_salt_pry[:,:,:], axis=0))
section_rho_pry_ann = np.squeeze(np.nanmean(section_rho_pry[:,:,:], axis=0))

# Barrier Polynya
# 1 - slice section for salt/temp/z_rho
section_salt_bar = np.empty((19,31,120)) # allocating
section_temp_bar = np.empty((19,31,120))
section_rho_bar = np.empty((19,31,120))
section_z_bar = np.empty((19,31,120))
section_z_bar_mask = np.empty((19,31,120))
for ee in np.arange(0,19):
    section_salt_bar[ee,:,:] = ds.salt.isel(xi_rho=slice(500,620), eta_rho=315, ocean_time=ee)
    section_temp_bar[ee,:,:] = ds.temp.isel(xi_rho=slice(500,620), eta_rho=315, ocean_time=ee)
    section_z_bar[ee,:,:] = ds.z_rho.isel(xi_rho=slice(500,620), eta_rho=315, ocean_time=ee)
    # 1.1 - mask land values in z_rho slice                                           
    section_z_bar_mask[...] = ma.array(section_z_bar[ee,:,:],mask=np.isnan(section_z_bar[ee,:,:]))
    section_rho_bar[ee,:,:] = np.squeeze(dens[ee,:,315,500:620])

# 2 - slide section for lat or lat
section_lat_bar_tmp = ds.lat_rho.isel(xi_rho=slice(500,620), eta_rho=315)
# 2.1 - mask land values for lat/lat needs loop for repeat through vert layers
section_lat_bar = np.ones((31,120))                                                
for ii in np.arange(0,31):
    section_lat_bar[ii,:] = section_lat_bar_tmp
section_lat_bar_mask = ma.array(section_lat_bar,mask=np.isnan(section_z_bar[0,:,:]))

print(section_lat_bar.shape, section_z_bar.shape, section_salt_bar.shape)

section_z_bar_mask_ann = np.squeeze(np.nanmean(section_z_bar_mask[:,:,:], axis=0))
section_temp_bar_ann = np.squeeze(np.nanmean(section_temp_bar[:,:,:], axis=0))
section_salt_bar_ann = np.squeeze(np.nanmean(section_salt_bar[:,:,:], axis=0))
section_rho_bar_ann = np.squeeze(np.nanmean(section_rho_bar[:,:,:], axis=0))

#Shackleton
# 1 - slice section for salt/temp/z_rho
section_salt_sha = np.empty((19,31,100)) # allocating
section_temp_sha = np.empty((19,31,100))
section_rho_sha = np.empty((19,31,100))
section_z_sha = np.empty((19,31,100))
section_z_sha_mask = np.empty((19,31,100))
for ee in np.arange(0,19):
    section_salt_sha[ee,:,:] = ds.salt.isel(xi_rho=slice(520,620), eta_rho=240, ocean_time=ee)
    section_temp_sha[ee,:,:] = ds.temp.isel(xi_rho=slice(520,620), eta_rho=240, ocean_time=ee)
    section_z_sha[ee,:,:] = ds.z_rho.isel(xi_rho=slice(520,620), eta_rho=240, ocean_time=ee)
    # 1.1 - mask land values in z_rho slice                                           
    section_z_sha_mask[...] = ma.array(section_z_sha[ee,:,:],mask=np.isnan(section_z_sha[ee,:,:]))
    section_rho_sha[ee,:,:] = np.squeeze(dens[ee,:,240,520:620])

# 2 - slide section for lat or lat
section_lat_sha_tmp = ds.lat_rho.isel(xi_rho=slice(520,620), eta_rho=240)
# 2.1 - mask land values for lat/lat needs loop for repeat through vert layers
section_lat_sha = np.ones((31,100))                                                
for ii in np.arange(0,31):
    section_lat_sha[ii,:] = section_lat_sha_tmp
section_lat_sha_mask = ma.array(section_lat_sha,mask=np.isnan(section_z_sha[0,:,:]))

print(section_lat_sha.shape, section_z_sha.shape, section_salt_sha.shape)

section_z_sha_mask_ann = np.squeeze(np.nanmean(section_z_sha_mask[:,:,:], axis=0))
section_temp_sha_ann = np.squeeze(np.nanmean(section_temp_sha[:,:,:], axis=0))
section_salt_sha_ann = np.squeeze(np.nanmean(section_salt_sha[:,:,:], axis=0))
section_rho_sha_ann = np.squeeze(np.nanmean(section_rho_sha[:,:,:], axis=0))

#Vincennes
# 1 - slice section for salt/temp/z_rho
section_salt_vin = np.empty((19,31,80)) # allocating
section_temp_vin = np.empty((19,31,80))
section_rho_vin = np.empty((19,31,80))
section_z_vin = np.empty((19,31,80))
section_z_vin_mask = np.empty((19,31,80))
for ee in np.arange(0,19):
    section_salt_vin[ee,:,:] = ds.salt.isel(xi_rho=slice(520,600), eta_rho=170, ocean_time=ee)
    section_temp_vin[ee,:,:] = ds.temp.isel(xi_rho=slice(520,600), eta_rho=170, ocean_time=ee)
    section_z_vin[ee,:,:] = ds.z_rho.isel(xi_rho=slice(520,600), eta_rho=170, ocean_time=ee)
    # 1.1 - mask land values in z_rho slice                                           
    section_z_vin_mask[...] = ma.array(section_z_vin[ee,:,:],mask=np.isnan(section_z_vin[ee,:,:]))
    section_rho_vin[ee,:,:] = np.squeeze(dens[ee,:,170,520:600])

# 2 - slide section for lat or lat
section_lat_vin_tmp = ds.lat_rho.isel(xi_rho=slice(520,600), eta_rho=170)
# 2.1 - mask land values for lat/lat needs loop for repeat through vert layers
section_lat_vin = np.ones((31,80))                                                
for ii in np.arange(0,31):
    section_lat_vin[ii,:] = section_lat_vin_tmp
section_lat_vin_mask = ma.array(section_lat_vin,mask=np.isnan(section_z_vin[0,:,:]))

print(section_lat_vin.shape, section_z_vin.shape, section_salt_vin.shape)

section_z_vin_mask_ann = np.squeeze(np.nanmean(section_z_vin_mask[:,:,:], axis=0))
section_temp_vin_ann = np.squeeze(np.nanmean(section_temp_vin[:,:,:], axis=0))
section_salt_vin_ann = np.squeeze(np.nanmean(section_salt_vin[:,:,:], axis=0))
section_rho_vin_ann = np.squeeze(np.nanmean(section_rho_vin[:,:,:], axis=0))

# Merz Glacier
# 1 - slice section for salt/temp/z_rho
section_salt_mer = np.empty((19,31,70)) # allocating
section_temp_mer = np.empty((19,31,70))
section_rho_mer = np.empty((19,31,70))
section_z_mer = np.empty((19,31,70))
section_z_mer_mask = np.empty((19,31,70))
for ee in np.arange(0,19):
    section_salt_mer[ee,:,:] = ds.salt.isel(xi_rho=430, eta_rho=slice(10,80), ocean_time=ee)
    section_temp_mer[ee,:,:] = ds.temp.isel(xi_rho=430, eta_rho=slice(10,80), ocean_time=ee)
    section_z_mer[ee,:,:] = ds.z_rho.isel(xi_rho=430, eta_rho=slice(10,80), ocean_time=ee)
    # 1.1 - mask land values in z_rho slice                                           
    section_z_mer_mask[ee,:,:] = ma.array(section_z_mer[ee,:,:],mask=np.isnan(section_z_mer[ee,:,:]))
    section_rho_mer[ee,:,:] = np.squeeze(dens[ee,:,10:80,430])
    
# 2 - slide section for lon or lat
section_lat_mer_tmp = ds.lat_rho.isel(xi_rho=430, eta_rho=slice(10,80))
# 2.1 - mask land values for lon/lat needs loop for repeat through vert layers
section_lat_mer = np.ones((31,70))                                                
for ii in np.arange(0,31):
    section_lat_mer[ii,:] = section_lat_mer_tmp
section_lat_mer_mask = ma.array(section_lat_mer,mask=np.isnan(section_z_mer[0,:,:]))

print(section_lat_mer.shape, section_z_mer.shape, section_salt_mer.shape)

section_z_mer_mask_ann = np.squeeze(np.nanmean(section_z_mer_mask[:,:,:], axis=0))
section_temp_mer_ann = np.squeeze(np.nanmean(section_temp_mer[:,:,:], axis=0))
section_salt_mer_ann = np.squeeze(np.nanmean(section_salt_mer[:,:,:], axis=0))
section_rho_mer_ann = np.squeeze(np.nanmean(section_rho_mer[:,:,:], axis=0))

# Plots transects: Monthly
fig_path='/users/boeiradi/COLD_project/postprocessing/figs/WAOM10/'

levelsT = np.arange(-2.25,2.26,.25)
levelsS = np.arange(33.2,35.1,.15)
levelsR = np.arange(36.2,37.3,.1)

# Plots transects: loop through storm event dates:

for ee in np.arange(0,19):
    fig = plt.figure(figsize=(20,26))
    ax7 = fig.add_subplot(6,3,1)
    ct = plt.contourf(section_lat_dar_mask, np.squeeze(section_z_dar_mask[ee,:,:]), np.squeeze(section_temp_dar[ee,:,:]), levels=levelsT, cmap=plt.cm.coolwarm)
    plt.contour(section_lat_dar_mask, np.squeeze(section_z_dar_mask[ee,:,:]), np.squeeze(section_temp_dar[ee,:,:]), colors='k', levels=levelsT) #[-2.3,-2.1,-1.9,-1.5,-1.0,-.5,.0])
#    plt.colorbar(ct, extend='both')
    plt.xlim([-70,-61])
    plt.ylim([-1000,0])
    plt.title('Darnley Polynya ' + date_storms[ee])
    plt.ylabel('Depth (m)')

    ax8 = fig.add_subplot(6,3,2)
    cs = plt.contourf(section_lat_dar_mask, np.squeeze(section_z_dar_mask[ee,:,:]), np.squeeze(section_salt_dar[ee,:,:]), levels=levelsS, linestyle='dash', cmap=plt.cm.coolwarm)
    plt.contour(section_lat_dar_mask, np.squeeze(section_z_dar_mask[ee,:,:]), np.squeeze(section_salt_dar[ee,:,:]), colors='k', levels=levelsS, linestyle='dashed') #[34.3,34.46,34.6,34.56,34.6])
#    plt.colorbar(cs, extend='both')
    plt.xlim([-70,-61])
    plt.ylim([-1000,0])

    ax9 = fig.add_subplot(6,3,3)
    cr = plt.contourf(section_lat_dar_mask, np.squeeze(section_z_dar_mask[ee,:,:]), np.squeeze(section_rho_dar[ee,:,:]), levels=levelsR, linestyle='dash', cmap=plt.cm.viridis)
    plt.contour(section_lat_dar_mask, np.squeeze(section_z_dar_mask[ee,:,:]), np.squeeze(section_rho_dar[ee,:,:]), colors='k', levels=levelsR, linestyle='dashed') #[34.3,34.46,34.6,34.56,34.6])
#    plt.colorbar(cs, extend='both')
    plt.xlim([-70,-61])
    plt.ylim([-1000,0])
    
    ax4 = fig.add_subplot(6,3,4)
    ct = plt.contourf(section_lat_pry_mask, np.squeeze(section_z_pry_mask[ee,:,:]), np.squeeze(section_temp_pry[ee,:,:]), levels=levelsT, cmap=plt.cm.coolwarm)
    plt.contour(section_lat_pry_mask, np.squeeze(section_z_pry_mask[ee,:,:]), np.squeeze(section_temp_pry[ee,:,:]), colors='k', levels=levelsT) #[-2.3,-2.1,-1.9,-1.5,-1.0,-.5,.0])
#    plt.colorbar(ct, extend='both')
    plt.xlim([-70,-61])
    plt.ylim([-1000,0])
    plt.title('Mackenzi Polynya ' + date_storms[ee])
    plt.ylabel('Depth (m)')

    ax5 = fig.add_subplot(6,3,5)
    cs = plt.contourf(section_lat_pry_mask, np.squeeze(section_z_pry_mask[ee,:,:]), np.squeeze(section_salt_pry[ee,:,:]), levels=levelsS, linestyle='dash', cmap=plt.cm.coolwarm)
    plt.contour(section_lat_pry_mask, np.squeeze(section_z_pry_mask[ee,:,:]), np.squeeze(section_salt_pry[ee,:,:]), colors='k', levels=levelsS, linestyle='dashed') #[34.3,34.46,34.6,34.56,34.6])
#    plt.colorbar(cs, extend='both')
    plt.xlim([-70,-61])
    plt.ylim([-1000,0])
    
    ax6 = fig.add_subplot(6,3,6)
    cr = plt.contourf(section_lat_pry_mask, np.squeeze(section_z_pry_mask[ee,:,:]), np.squeeze(section_rho_pry[ee,:,:]), levels=levelsR, linestyle='dash', cmap=plt.cm.viridis)
    plt.contour(section_lat_pry_mask, np.squeeze(section_z_pry_mask[ee,:,:]), np.squeeze(section_rho_pry[ee,:,:]), colors='k', levels=levelsR, linestyle='dashed') #[34.3,34.46,34.6,34.56,34.6])
#    plt.colorbar(cs, extend='both')
    plt.xlim([-70,-61])
    plt.ylim([-1000,0])
    
    ax1 = fig.add_subplot(6,3,7)
    ct = plt.contourf(section_lat_bar_mask, np.squeeze(section_z_bar_mask[ee,:,:]), np.squeeze(section_temp_bar[ee,:,:]), levels=levelsT, cmap=plt.cm.coolwarm)
    plt.contour(section_lat_bar_mask, np.squeeze(section_z_bar_mask[ee,:,:]), np.squeeze(section_temp_bar[ee,:,:]), colors='k', levels=levelsT) #[-2.3,-2.1,-1.9,-1.5,-1.0,-.5,.0])
#    plt.colorbar(ct, extend='both')
    plt.xlim([-70,-61])
    plt.ylim([-1000,0])
    plt.title('Barrier Polynya ' + date_storms[ee])
    plt.ylabel('Depth (m)')

    ax2 = fig.add_subplot(6,3,8)
    cs = plt.contourf(section_lat_bar_mask, np.squeeze(section_z_bar_mask[ee,:,:]), np.squeeze(section_salt_bar[ee,:,:]), levels=levelsS, linestyle='dash', cmap=plt.cm.coolwarm)
    plt.contour(section_lat_bar_mask, np.squeeze(section_z_bar_mask[ee,:,:]), np.squeeze(section_salt_bar[ee,:,:]), colors='k', levels=levelsS, linestyle='dashed') #[34.3,34.46,34.6,34.56,34.6])
#    plt.colorbar(cs, extend='both')
    plt.xlim([-70,-61])
    plt.ylim([-1000,0])
    
    ax3 = fig.add_subplot(6,3,9)
    cr = plt.contourf(section_lat_bar_mask, np.squeeze(section_z_bar_mask[ee,:,:]), np.squeeze(section_rho_bar[ee,:,:]), levels=levelsR, linestyle='dash', cmap=plt.cm.viridis)
    plt.contour(section_lat_bar_mask, np.squeeze(section_z_bar_mask[ee,:,:]), np.squeeze(section_rho_bar[ee,:,:]), colors='k', levels=levelsR, linestyle='dashed') #[34.3,34.46,34.6,34.56,34.6])
#    plt.colorbar(cs, extend='both')
    plt.xlim([-70,-61])
    plt.ylim([-1000,0])
    
    ax16 = fig.add_subplot(6,3,10)
    ct = plt.contourf(section_lat_sha_mask, np.squeeze(section_z_sha_mask[ee,:,:]), np.squeeze(section_temp_sha[ee,:,:]), levels=levelsT, cmap=plt.cm.coolwarm)
    plt.contour(section_lat_sha_mask, np.squeeze(section_z_sha_mask[ee,:,:]), np.squeeze(section_temp_sha[ee,:,:]), colors='k', levels=levelsT) #[-2.3,-2.1,-1.9,-1.5,-1.0,-.5,.0])
#    plt.colorbar(ct, extend='both')
    plt.xlim([-66,-61.5])
    plt.ylim([-1000,0])
    plt.title('Shackleton Polynya ' + date_storms[ee])
    plt.ylabel('Depth (m)')

    ax17 = fig.add_subplot(6,3,11)
    cs = plt.contourf(section_lat_sha_mask, np.squeeze(section_z_sha_mask[ee,:,:]), np.squeeze(section_salt_sha[ee,:,:]), levels=levelsS, linestyle='dash', cmap=plt.cm.coolwarm)
    plt.contour(section_lat_sha_mask, np.squeeze(section_z_sha_mask[ee,:,:]), np.squeeze(section_salt_sha[ee,:,:]), colors='k', levels=levelsS, linestyle='dashed') #[34.3,34.46,34.6,34.56,34.6])
#    plt.colorbar(cs, extend='both')
    plt.xlim([-66,-61.5])
    plt.ylim([-1000,0])
    
    ax18 = fig.add_subplot(6,3,12)
    cr = plt.contourf(section_lat_sha_mask, np.squeeze(section_z_sha_mask[ee,:,:]), np.squeeze(section_rho_sha[ee,:,:]), levels=levelsR, linestyle='dash', cmap=plt.cm.viridis)
    plt.contour(section_lat_sha_mask, np.squeeze(section_z_sha_mask[ee,:,:]), np.squeeze(section_rho_sha[ee,:,:]), colors='k', levels=levelsR, linestyle='dashed') #[34.3,34.46,34.6,34.56,34.6])
#    plt.colorbar(cs, extend='both')
    plt.xlim([-66,-61.5])
    plt.ylim([-1000,0])

    ax13 = fig.add_subplot(6,3,13)
    ct = plt.contourf(section_lat_vin_mask, np.squeeze(section_z_vin_mask[ee,:,:]), np.squeeze(section_temp_vin[ee,:,:]), levels=levelsT, cmap=plt.cm.coolwarm)
    plt.contour(section_lat_vin_mask, np.squeeze(section_z_vin_mask[ee,:,:]), np.squeeze(section_temp_vin[ee,:,:]), colors='k', levels=levelsT) #[-2.3,-2.1,-1.9,-1.5,-1.0,-.5,.0])
#    plt.colorbar(ct, extend='both')
    plt.xlim([-66,-61.8])
    plt.ylim([-1000,0])
    plt.title('Vincennes Polynya ' + date_storms[ee])
    plt.ylabel('Depth (m)')

    ax14 = fig.add_subplot(6,3,14)
    cs = plt.contourf(section_lat_vin_mask, np.squeeze(section_z_vin_mask[ee,:,:]), np.squeeze(section_salt_vin[ee,:,:]), levels=levelsS, linestyle='dash', cmap=plt.cm.coolwarm)
    plt.contour(section_lat_vin_mask, np.squeeze(section_z_vin_mask[ee,:,:]), np.squeeze(section_salt_vin[ee,:,:]), colors='k', levels=levelsS, linestyle='dashed') #[34.3,34.46,34.6,34.56,34.6])
#    plt.colorbar(cs, extend='both')
    plt.xlim([-66,-61.8])
    plt.ylim([-1000,0])
    
    ax15 = fig.add_subplot(6,3,15)
    cr = plt.contourf(section_lat_vin_mask, np.squeeze(section_z_vin_mask[ee,:,:]), np.squeeze(section_rho_vin[ee,:,:]), levels=levelsR, linestyle='dash', cmap=plt.cm.viridis)
    plt.contour(section_lat_vin_mask, np.squeeze(section_z_vin_mask[ee,:,:]), np.squeeze(section_rho_vin[ee,:,:]), colors='k', levels=levelsR, linestyle='dashed') #[34.3,34.46,34.6,34.56,34.6])
#    plt.colorbar(cs, extend='both')
    plt.xlim([-66,-61.6])
    plt.ylim([-1000,0])
        
    ax10 = fig.add_subplot(6,3,16)
    ct = plt.contourf(section_lat_mer_mask, np.squeeze(section_z_mer_mask[ee,:,:]), np.squeeze(section_temp_mer[ee,:,:]), levels=levelsT, cmap=plt.cm.coolwarm)
    plt.contour(section_lat_mer_mask, np.squeeze(section_z_mer_mask[ee,:,:]), np.squeeze(section_temp_mer[ee,:,:]), colors='k', levels=levelsT) #[-2.3,-2.1,-1.9,-1.5,-1.0,-.5,.0])
#    plt.colorbar(ct, extend='both')
    plt.xlim([-70,-63.7])
    plt.ylim([-1000,0])
    plt.title('Merz Glacier ' + date_storms[ee])
    plt.ylabel('Depth (m)')
    plt.xlabel('Latitude')

    ax11 = fig.add_subplot(6,3,17)
    cs = plt.contourf(section_lat_mer_mask, np.squeeze(section_z_mer_mask[ee,:,:]), np.squeeze(section_salt_mer[ee,:,:]), levels=levelsS, linestyle='dash', cmap=plt.cm.coolwarm)
    plt.contour(section_lat_mer_mask, np.squeeze(section_z_mer_mask[ee,:,:]), np.squeeze(section_salt_mer[ee,:,:]), colors='k', levels=levelsS, linestyle='dashed') #[34.3,34.46,34.6,34.56,34.6])
#    plt.colorbar(cs, extend='both')
    plt.xlim([-70,-63.7])
    plt.ylim([-1000,0])
    plt.xlabel('Latitude')

    ax12 = fig.add_subplot(6,3,18)
    cr = plt.contourf(section_lat_mer_mask, np.squeeze(section_z_mer_mask[ee,:,:]), np.squeeze(section_rho_mer[ee,:,:]), levels=levelsR, linestyle='dash', cmap=plt.cm.viridis)
    plt.contour(section_lat_mer_mask, np.squeeze(section_z_mer_mask[ee,:,:]), np.squeeze(section_rho_mer[ee,:,:]), colors='k', levels=levelsR, linestyle='dashed') #[34.3,34.46,34.6,34.56,34.6])
#    plt.colorbar(cs, extend='both')
    plt.xlim([-70,-63.7])
    plt.ylim([-1000,0])
    plt.xlabel('Latitude')
    
    # temp colorbar
    cbarT_axim = fig.add_axes([0.125, 0.09, 0.23, 0.01])
    cbarT =fig.colorbar(ct, cax=cbarT_axim, extend='both', orientation='horizontal')
    cbarT.ax.tick_params(labelsize=14)
    plt.xlabel('$^{\circ}$C',fontsize=14)
    # salt colorbar
    cbarS_axim = fig.add_axes([0.395, 0.09, 0.23, 0.01])
    cbarS =fig.colorbar(cs, cax=cbarS_axim, extend='both', orientation='horizontal')
    cbarS.ax.tick_params(labelsize=14)
    plt.xlabel('psu',fontsize=14)
    # rho colorbar
    cbarR_axim = fig.add_axes([0.67, 0.09, 0.23, 0.01])
    cbarR =fig.colorbar(cr, cax=cbarR_axim, extend='both', orientation='horizontal')
    cbarR.ax.tick_params(labelsize=14)
    plt.xlabel('kg m$^{-3}$',fontsize=14)
    
    name_fig="waom10_All_EastAntarct_sections_event=" + date_storms[ee] + ".png"
    plt.savefig(fig_path + name_fig)

