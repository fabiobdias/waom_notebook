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
import matplotlib as mpl
mpl.use('Agg')

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


# dar
# slice(460,615), eta_rho=380
# pry
# xi_rho=slice(460,615), eta_rho=368
# bar
# xi_rho=slice(500,620), eta_rho=345
# sha
# xi_rho=slice(520,620), eta_rho=270
# vin
# xi_rho=slice(520,600), eta_rho=200
# mer
# xi_rho=435, eta_rho=slice(40,110)

# index for eta var sections (mer)
xi_pt = [435]
eta_sec_ini = [40]
eta_sec_end = [110]
# index for xi var sections (dar, pry, bar, sha, vin)
xi_sec_ini = [460, 460, 500, 520, 520]
xi_sec_end = [615, 615, 620, 620, 600]
#eta_pt = [380, 368, 345, 270, 200]
eta_pt = [388, 368, 342, 270, 200] # testing Darnley further West
# length of the section
sec_len = [155, 155, 120, 100, 80, 70]

def read_section_roms(sec_name,ind_ij, ind_len):

    # allocate sections arrays:
    section_salt = np.empty((73,31,sec_len[ind_len]))
    section_temp = np.empty((73,31,sec_len[ind_len]))
    section_rho = np.empty((73,31,sec_len[ind_len]))
    section_z = np.empty((73,31,sec_len[ind_len]))
    section_z_mask = np.empty((73,31,sec_len[ind_len]))

    monthly_ind_ini = [0, 7, 13, 19, 25, 31, 37, 43, 49, 55, 61, 67]
    monthly_ind_end = [6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 73]

# load ROMS avg output
    count = 0
    for mm in ['01','02','03','04','05','06','07','08','09','10','11','12']:
        ds = xr.open_dataset("/scratch/project_2000789/boeiradi/waom10extend_shflim_S_0.25Q/output_21yr/ocean_avg_00" + mm + ".nc")
        print('size temp and time length: ', ds.temp.shape, len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
        if sec_name == 'mer':
            print(sec_name, 'xi =', xi_pt[ind_ij], 'eta =', eta_sec_ini[ind_ij], eta_sec_end[ind_ij])
            section_salt_avg = ds.salt.isel(xi_rho=xi_pt[ind_ij], eta_rho=slice(eta_sec_ini[ind_ij],eta_sec_end[ind_ij]))
            section_temp_avg = ds.temp.isel(xi_rho=xi_pt[ind_ij], eta_rho=slice(eta_sec_ini[ind_ij],eta_sec_end[ind_ij]))
            section_zeta = ds.zeta.isel(xi_rho=xi_pt[ind_ij], eta_rho=slice(eta_sec_ini[ind_ij],eta_sec_end[ind_ij]))
            dens_avg = gsw.rho(section_salt_avg,section_temp_avg,2000) - 1000
        elif sec_name == 'dar' or sec_name ==  'pry' or sec_name ==  'bar' or sec_name ==  'sha' or sec_name ==  'vin':
            print(sec_name, 'eta =', eta_pt[ind_ij], 'xi =', xi_sec_ini[ind_ij], xi_sec_end[ind_ij])
            section_salt_avg = ds.salt.isel(xi_rho=slice(xi_sec_ini[ind_ij],xi_sec_end[ind_ij]), eta_rho=eta_pt[ind_ij])
            section_temp_avg = ds.temp.isel(xi_rho=slice(xi_sec_ini[ind_ij],xi_sec_end[ind_ij]), eta_rho=eta_pt[ind_ij])
            section_zeta = ds.zeta.isel(xi_rho=slice(xi_sec_ini[ind_ij],xi_sec_end[ind_ij]), eta_rho=eta_pt[ind_ij])
            dens_avg = gsw.rho(section_salt_avg,section_temp_avg,2000) - 1000

        if mm == '01':
            ds = ds.set_coords(['Cs_r', 'Cs_w', 'hc', 'Vtransform'])
            if sec_name ==  'mer':
                h = ds.h.isel(xi_rho=xi_pt[ind_ij], eta_rho=slice(eta_sec_ini[ind_ij],eta_sec_end[ind_ij]))
            elif sec_name == 'dar' or sec_name ==  'pry' or sec_name ==  'bar' or sec_name ==  'sha' or sec_name ==  'vin':
                h = ds.h.isel(xi_rho=slice(xi_sec_ini[ind_ij],xi_sec_end[ind_ij]), eta_rho=eta_pt[ind_ij])

        print("Vtransform=2")
        Zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * h) / (ds.hc + h)
        print('Zo_rho, zeta, and h sizes: ', Zo_rho.shape, section_zeta.shape, h.shape)
        z_rho_tmp = section_zeta + (section_zeta + h) * Zo_rho
        print('z_rho_tmp size: ', z_rho_tmp.shape)
        del section_zeta, Zo_rho
        print('temp/salt/z_rho_tmp size = ', section_temp_avg.shape, section_salt_avg.shape, z_rho_tmp.shape)
        # assign to monthly sections:
        section_salt[monthly_ind_ini[count]:monthly_ind_end[count]+1,:] = section_salt_avg
        section_temp[monthly_ind_ini[count]:monthly_ind_end[count]+1,:] = section_temp_avg
        if sec_name ==  'mer':
            z_rho_tmp = z_rho_tmp.transpose('ocean_time', 's_rho', 'eta_rho')
        elif sec_name == 'dar' or sec_name ==  'pry' or sec_name ==  'bar' or sec_name ==  'sha' or sec_name ==  'vin':
            z_rho_tmp = z_rho_tmp.transpose('ocean_time', 's_rho', 'xi_rho')
        section_z[monthly_ind_ini[count]:monthly_ind_end[count]+1,:] = z_rho_tmp
        section_rho[monthly_ind_ini[count]:monthly_ind_end[count]+1,:] = dens_avg
        # mask
        section_z_mask[monthly_ind_ini[count]:monthly_ind_end[count]+1,:] = ma.array(z_rho_tmp, mask=np.isnan(z_rho_tmp))
        del section_salt_avg, section_temp_avg, dens_avg, z_rho_tmp
        count = count+1
    del h
    return section_salt, section_temp, section_z, section_rho, section_z_mask

section_salt_dar, section_temp_dar, section_z_dar, section_rho_dar, section_z_dar_mask = read_section_roms('dar',0, 0)
section_salt_pry, section_temp_pry, section_z_pry, section_rho_pry, section_z_pry_mask = read_section_roms('pry',1, 1)
section_salt_bar, section_temp_bar, section_z_bar, section_rho_bar, section_z_bar_mask = read_section_roms('bar',2, 2)
section_salt_sha, section_temp_sha, section_z_sha, section_rho_sha, section_z_sha_mask = read_section_roms('sha',3, 3)
section_salt_vin, section_temp_vin, section_z_vin, section_rho_vin, section_z_vin_mask = read_section_roms('vin',4, 4)
section_salt_mer, section_temp_mer, section_z_mer, section_rho_mer, section_z_mer_mask = read_section_roms('mer',0, 5)



dg = xr.open_dataset("/scratch/project_2000339/boeiradi/waom10_frc/waom10extend_grd.nc")

lat_rho = dg.variables["lat_rho"]
lon_rho = dg.variables["lon_rho"]
#dx.coords['lat_rho']=lat_rho.transpose() # put lat_rho into ds dataset
#dx.coords['lon_rho']=lon_rho.transpose() # put lon_rho into ds dataset

# Period for storms in the winter of 2007:
events_ind = np.array([31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 52, 53])
date_storms = [('03.Jun.2007'),('08.Jun.2007'),('13.Jun.2007'),('18.Jun.2007'),('23.Jun.2007'),('28.Jun.2007'), \
               ('03.Jul.2007'),('08.Jul.2007'),('13.Jul.2007'),('18.Jul.2007'),('23.Jul.2007'),('28.Jul.2007'), \
               ('02.Aug.2007'),('07.Aug.2007'),('12.Aug.2007'),('17.Aug.2007'),('22.Aug.2007'), \
               ('22.Sep.2007'),('27.Sep.2007')]

# Darnley
ind_ij = 0
ind_len = 0
# 2 - slide section for lon or lat
section_lat_dar_tmp = dg.lat_rho.isel(xi_rho=slice(xi_sec_ini[ind_ij],xi_sec_end[ind_ij]), eta_rho=eta_pt[ind_ij])
# 2.1 - mask land values for lon/lat needs loop for repeat through vert layers
section_lat_dar = np.ones((31,sec_len[ind_len]))
for ii in np.arange(0,31):
    section_lat_dar[ii,:] = section_lat_dar_tmp
section_lat_dar_mask = ma.array(section_lat_dar,mask=np.isnan(section_z_dar[0,:,:]))
print(section_lat_dar.shape, section_z_dar.shape, section_salt_dar.shape)
# annual mean
section_z_dar_mask_ann = np.squeeze(np.nanmean(section_z_dar_mask[:,:,:], axis=0))
section_temp_dar_ann = np.squeeze(np.nanmean(section_temp_dar[:,:,:], axis=0))
section_salt_dar_ann = np.squeeze(np.nanmean(section_salt_dar[:,:,:], axis=0))
section_rho_dar_ann = np.squeeze(np.nanmean(section_rho_dar[:,:,:], axis=0))

# Mackenzi
ind_ij = 1
ind_len = 1
# 2 - slide section for lon or lat
section_lat_pry_tmp = dg.lat_rho.isel(xi_rho=slice(xi_sec_ini[ind_ij],xi_sec_end[ind_ij]), eta_rho=eta_pt[ind_ij])
# 2.1 - mask land values for lon/lat needs loop for repeat through vert layers
section_lat_pry = np.ones((31,sec_len[ind_len]))
for ii in np.arange(0,31):
    section_lat_pry[ii,:] = section_lat_pry_tmp
section_lat_pry_mask = ma.array(section_lat_pry,mask=np.isnan(section_z_pry[0,:,:]))
print(section_lat_pry.shape, section_z_pry.shape, section_salt_pry.shape)
# annual mean
section_z_pry_mask_ann = np.squeeze(np.nanmean(section_z_pry_mask[:,:,:], axis=0))
section_temp_pry_ann = np.squeeze(np.nanmean(section_temp_pry[:,:,:], axis=0))
section_salt_pry_ann = np.squeeze(np.nanmean(section_salt_pry[:,:,:], axis=0))
section_rho_pry_ann = np.squeeze(np.nanmean(section_rho_pry[:,:,:], axis=0))

# Barrier
ind_ij = 2
ind_len = 2
# 2 - slide section for lon or lat
section_lat_bar_tmp = dg.lat_rho.isel(xi_rho=slice(xi_sec_ini[ind_ij],xi_sec_end[ind_ij]), eta_rho=eta_pt[ind_ij])
# 2.1 - mask land values for lon/lat needs loop for repeat through vert layers
section_lat_bar = np.ones((31,sec_len[ind_len]))
for ii in np.arange(0,31):
    section_lat_bar[ii,:] = section_lat_bar_tmp
section_lat_bar_mask = ma.array(section_lat_bar,mask=np.isnan(section_z_bar[0,:,:]))
print(section_lat_bar.shape, section_z_bar.shape, section_salt_bar.shape)
# annual mean
section_z_bar_mask_ann = np.squeeze(np.nanmean(section_z_bar_mask[:,:,:], axis=0))
section_temp_bar_ann = np.squeeze(np.nanmean(section_temp_bar[:,:,:], axis=0))
section_salt_bar_ann = np.squeeze(np.nanmean(section_salt_bar[:,:,:], axis=0))
section_rho_bar_ann = np.squeeze(np.nanmean(section_rho_bar[:,:,:], axis=0))

# Shackleton
ind_ij = 3
ind_len = 3
# 2 - slide section for lon or lat
section_lat_sha_tmp = dg.lat_rho.isel(xi_rho=slice(xi_sec_ini[ind_ij],xi_sec_end[ind_ij]), eta_rho=eta_pt[ind_ij])
# 2.1 - mask land values for lon/lat needs loop for repeat through vert layers
section_lat_sha = np.ones((31,sec_len[ind_len]))
for ii in np.arange(0,31):
    section_lat_sha[ii,:] = section_lat_sha_tmp
section_lat_sha_mask = ma.array(section_lat_sha,mask=np.isnan(section_z_sha[0,:,:]))
print(section_lat_sha.shape, section_z_sha.shape, section_salt_sha.shape)
# annual mean
section_z_sha_mask_ann = np.squeeze(np.nanmean(section_z_sha_mask[:,:,:], axis=0))
section_temp_sha_ann = np.squeeze(np.nanmean(section_temp_sha[:,:,:], axis=0))
section_salt_sha_ann = np.squeeze(np.nanmean(section_salt_sha[:,:,:], axis=0))
section_rho_sha_ann = np.squeeze(np.nanmean(section_rho_sha[:,:,:], axis=0))

# Vincennes
ind_ij = 4
ind_len = 4
# 2 - slide section for lon or lat
section_lat_vin_tmp = dg.lat_rho.isel(xi_rho=slice(xi_sec_ini[ind_ij],xi_sec_end[ind_ij]), eta_rho=eta_pt[ind_ij])
# 2.1 - mask land values for lon/lat needs loop for repeat through vert layers
section_lat_vin = np.ones((31,sec_len[ind_len]))
for ii in np.arange(0,31):
    section_lat_vin[ii,:] = section_lat_vin_tmp
section_lat_vin_mask = ma.array(section_lat_vin,mask=np.isnan(section_z_vin[0,:,:]))
print(section_lat_vin.shape, section_z_vin.shape, section_salt_vin.shape)
# annual mean
section_z_vin_mask_ann = np.squeeze(np.nanmean(section_z_vin_mask[:,:,:], axis=0))
section_temp_vin_ann = np.squeeze(np.nanmean(section_temp_vin[:,:,:], axis=0))
section_salt_vin_ann = np.squeeze(np.nanmean(section_salt_vin[:,:,:], axis=0))
section_rho_vin_ann = np.squeeze(np.nanmean(section_rho_vin[:,:,:], axis=0))

# Merz
ind_ij = 0
ind_len = 5
# 2 - slide section for lon or lat
section_lat_mer_tmp = dg.lat_rho.isel(eta_rho=slice(eta_sec_ini[ind_ij],eta_sec_end[ind_ij]), xi_rho=xi_pt[ind_ij])
# 2.1 - mask land values for lon/lat needs loop for repeat through vert layers
section_lat_mer = np.ones((31,sec_len[ind_len]))
for ii in np.arange(0,31):
    section_lat_mer[ii,:] = section_lat_mer_tmp
section_lat_mer_mask = ma.array(section_lat_mer,mask=np.isnan(section_z_mer[0,:,:]))
print(section_lat_mer.shape, section_z_mer.shape, section_salt_mer.shape)
# annual mean
section_z_mer_mask_ann = np.squeeze(np.nanmean(section_z_mer_mask[:,:,:], axis=0))
section_temp_mer_ann = np.squeeze(np.nanmean(section_temp_mer[:,:,:], axis=0))
section_salt_mer_ann = np.squeeze(np.nanmean(section_salt_mer[:,:,:], axis=0))
section_rho_mer_ann = np.squeeze(np.nanmean(section_rho_mer[:,:,:], axis=0))



# Plots transects: Monthly
fig_path='/users/boeiradi/COLD_project/postprocessing/figs/Storm_analyses/'

levelsT = np.arange(-2.25,2.26,.25)
levelsS = np.arange(34.,34.8,.1)
levelsR = np.arange(36.2,37.3,.1)

# Plots transects: loop through storm event dates:

for ee in np.arange(0,19):
    fig = plt.figure(figsize=(20,26))
    ax7 = fig.add_subplot(6,3,1)
    ct = plt.contourf(section_lat_dar_mask, np.squeeze(section_z_dar_mask[ee,:,:]), np.squeeze(section_temp_dar[ee,:,:]), levels=levelsT, cmap=plt.cm.coolwarm, extend='both')
    plt.contour(section_lat_dar_mask, np.squeeze(section_z_dar_mask[ee,:,:]), np.squeeze(section_temp_dar[ee,:,:]), colors='k', levels=levelsT) #[-2.3,-2.1,-1.9,-1.5,-1.0,-.5,.0])
#    plt.colorbar(ct, extend='both')
    plt.xlim([-70,-61])
    plt.ylim([-600,0])
    plt.title('Darnley Polynya ' + date_storms[ee])
    plt.ylabel('Depth (m)')

    ax8 = fig.add_subplot(6,3,2)
    cs = plt.contourf(section_lat_dar_mask, np.squeeze(section_z_dar_mask[ee,:,:]), np.squeeze(section_salt_dar[ee,:,:]), levels=levelsS, linestyle='dash', cmap=plt.cm.coolwarm, extend='both')
    plt.contour(section_lat_dar_mask, np.squeeze(section_z_dar_mask[ee,:,:]), np.squeeze(section_salt_dar[ee,:,:]), colors='k', levels=levelsS, linestyle='dashed') #[34.3,34.46,34.6,34.56,34.6])
#    plt.colorbar(cs, extend='both')
    plt.xlim([-70,-61])
    plt.ylim([-600,0])

    ax9 = fig.add_subplot(6,3,3)
    cr = plt.contourf(section_lat_dar_mask, np.squeeze(section_z_dar_mask[ee,:,:]), np.squeeze(section_rho_dar[ee,:,:]), levels=levelsR, linestyle='dash', cmap=plt.cm.viridis, extend='both')
    plt.contour(section_lat_dar_mask, np.squeeze(section_z_dar_mask[ee,:,:]), np.squeeze(section_rho_dar[ee,:,:]), colors='k', levels=levelsR, linestyle='dashed') #[34.3,34.46,34.6,34.56,34.6])
#    plt.colorbar(cs, extend='both')
    plt.xlim([-70,-61])
    plt.ylim([-600,0])
    
    ax4 = fig.add_subplot(6,3,4)
    ct = plt.contourf(section_lat_pry_mask, np.squeeze(section_z_pry_mask[ee,:,:]), np.squeeze(section_temp_pry[ee,:,:]), levels=levelsT, cmap=plt.cm.coolwarm, extend='both')
    plt.contour(section_lat_pry_mask, np.squeeze(section_z_pry_mask[ee,:,:]), np.squeeze(section_temp_pry[ee,:,:]), colors='k', levels=levelsT) #[-2.3,-2.1,-1.9,-1.5,-1.0,-.5,.0])
#    plt.colorbar(ct, extend='both')
    plt.xlim([-70,-61])
    plt.ylim([-600,0])
    plt.title('Mackenzi Polynya ' + date_storms[ee])
    plt.ylabel('Depth (m)')

    ax5 = fig.add_subplot(6,3,5)
    cs = plt.contourf(section_lat_pry_mask, np.squeeze(section_z_pry_mask[ee,:,:]), np.squeeze(section_salt_pry[ee,:,:]), levels=levelsS, linestyle='dash', cmap=plt.cm.coolwarm, extend='both')
    plt.contour(section_lat_pry_mask, np.squeeze(section_z_pry_mask[ee,:,:]), np.squeeze(section_salt_pry[ee,:,:]), colors='k', levels=levelsS, linestyle='dashed') #[34.3,34.46,34.6,34.56,34.6])
#    plt.colorbar(cs, extend='both')
    plt.xlim([-70,-61])
    plt.ylim([-600,0])
    
    ax6 = fig.add_subplot(6,3,6)
    cr = plt.contourf(section_lat_pry_mask, np.squeeze(section_z_pry_mask[ee,:,:]), np.squeeze(section_rho_pry[ee,:,:]), levels=levelsR, linestyle='dash', cmap=plt.cm.viridis, extend='both')
    plt.contour(section_lat_pry_mask, np.squeeze(section_z_pry_mask[ee,:,:]), np.squeeze(section_rho_pry[ee,:,:]), colors='k', levels=levelsR, linestyle='dashed') #[34.3,34.46,34.6,34.56,34.6])
#    plt.colorbar(cs, extend='both')
    plt.xlim([-70,-61])
    plt.ylim([-600,0])
    
    ax1 = fig.add_subplot(6,3,7)
    ct = plt.contourf(section_lat_bar_mask, np.squeeze(section_z_bar_mask[ee,:,:]), np.squeeze(section_temp_bar[ee,:,:]), levels=levelsT, cmap=plt.cm.coolwarm, extend='both')
    plt.contour(section_lat_bar_mask, np.squeeze(section_z_bar_mask[ee,:,:]), np.squeeze(section_temp_bar[ee,:,:]), colors='k', levels=levelsT) #[-2.3,-2.1,-1.9,-1.5,-1.0,-.5,.0])
#    plt.colorbar(ct, extend='both')
    plt.xlim([-70,-61])
    plt.ylim([-600,0])
    plt.title('Barrier Polynya ' + date_storms[ee])
    plt.ylabel('Depth (m)')

    ax2 = fig.add_subplot(6,3,8)
    cs = plt.contourf(section_lat_bar_mask, np.squeeze(section_z_bar_mask[ee,:,:]), np.squeeze(section_salt_bar[ee,:,:]), levels=levelsS, linestyle='dash', cmap=plt.cm.coolwarm, extend='both')
    plt.contour(section_lat_bar_mask, np.squeeze(section_z_bar_mask[ee,:,:]), np.squeeze(section_salt_bar[ee,:,:]), colors='k', levels=levelsS, linestyle='dashed') #[34.3,34.46,34.6,34.56,34.6])
#    plt.colorbar(cs, extend='both')
    plt.xlim([-70,-61])
    plt.ylim([-600,0])
    
    ax3 = fig.add_subplot(6,3,9)
    cr = plt.contourf(section_lat_bar_mask, np.squeeze(section_z_bar_mask[ee,:,:]), np.squeeze(section_rho_bar[ee,:,:]), levels=levelsR, linestyle='dash', cmap=plt.cm.viridis, extend='both')
    plt.contour(section_lat_bar_mask, np.squeeze(section_z_bar_mask[ee,:,:]), np.squeeze(section_rho_bar[ee,:,:]), colors='k', levels=levelsR, linestyle='dashed') #[34.3,34.46,34.6,34.56,34.6])
#    plt.colorbar(cs, extend='both')
    plt.xlim([-70,-61])
    plt.ylim([-600,0])
    
    ax16 = fig.add_subplot(6,3,10)
    ct = plt.contourf(section_lat_sha_mask, np.squeeze(section_z_sha_mask[ee,:,:]), np.squeeze(section_temp_sha[ee,:,:]), levels=levelsT, cmap=plt.cm.coolwarm, extend='both')
    plt.contour(section_lat_sha_mask, np.squeeze(section_z_sha_mask[ee,:,:]), np.squeeze(section_temp_sha[ee,:,:]), colors='k', levels=levelsT) #[-2.3,-2.1,-1.9,-1.5,-1.0,-.5,.0])
#    plt.colorbar(ct, extend='both')
    plt.xlim([-66,-61.5])
    plt.ylim([-600,0])
    plt.title('Shackleton Polynya ' + date_storms[ee])
    plt.ylabel('Depth (m)')

    ax17 = fig.add_subplot(6,3,11)
    cs = plt.contourf(section_lat_sha_mask, np.squeeze(section_z_sha_mask[ee,:,:]), np.squeeze(section_salt_sha[ee,:,:]), levels=levelsS, linestyle='dash', cmap=plt.cm.coolwarm, extend='both')
    plt.contour(section_lat_sha_mask, np.squeeze(section_z_sha_mask[ee,:,:]), np.squeeze(section_salt_sha[ee,:,:]), colors='k', levels=levelsS, linestyle='dashed') #[34.3,34.46,34.6,34.56,34.6])
#    plt.colorbar(cs, extend='both')
    plt.xlim([-66,-61.5])
    plt.ylim([-600,0])
    
    ax18 = fig.add_subplot(6,3,12)
    cr = plt.contourf(section_lat_sha_mask, np.squeeze(section_z_sha_mask[ee,:,:]), np.squeeze(section_rho_sha[ee,:,:]), levels=levelsR, linestyle='dash', cmap=plt.cm.viridis, extend='both')
    plt.contour(section_lat_sha_mask, np.squeeze(section_z_sha_mask[ee,:,:]), np.squeeze(section_rho_sha[ee,:,:]), colors='k', levels=levelsR, linestyle='dashed') #[34.3,34.46,34.6,34.56,34.6])
#    plt.colorbar(cs, extend='both')
    plt.xlim([-66,-61.5])
    plt.ylim([-600,0])

    ax13 = fig.add_subplot(6,3,13)
    ct = plt.contourf(section_lat_vin_mask, np.squeeze(section_z_vin_mask[ee,:,:]), np.squeeze(section_temp_vin[ee,:,:]), levels=levelsT, cmap=plt.cm.coolwarm, extend='both')
    plt.contour(section_lat_vin_mask, np.squeeze(section_z_vin_mask[ee,:,:]), np.squeeze(section_temp_vin[ee,:,:]), colors='k', levels=levelsT) #[-2.3,-2.1,-1.9,-1.5,-1.0,-.5,.0])
#    plt.colorbar(ct, extend='both')
    plt.xlim([-66,-61.8])
    plt.ylim([-600,0])
    plt.title('Vincennes Polynya ' + date_storms[ee])
    plt.ylabel('Depth (m)')

    ax14 = fig.add_subplot(6,3,14)
    cs = plt.contourf(section_lat_vin_mask, np.squeeze(section_z_vin_mask[ee,:,:]), np.squeeze(section_salt_vin[ee,:,:]), levels=levelsS, linestyle='dash', cmap=plt.cm.coolwarm, extend='both')
    plt.contour(section_lat_vin_mask, np.squeeze(section_z_vin_mask[ee,:,:]), np.squeeze(section_salt_vin[ee,:,:]), colors='k', levels=levelsS, linestyle='dashed') #[34.3,34.46,34.6,34.56,34.6])
#    plt.colorbar(cs, extend='both')
    plt.xlim([-66,-61.8])
    plt.ylim([-600,0])
    
    ax15 = fig.add_subplot(6,3,15)
    cr = plt.contourf(section_lat_vin_mask, np.squeeze(section_z_vin_mask[ee,:,:]), np.squeeze(section_rho_vin[ee,:,:]), levels=levelsR, linestyle='dash', cmap=plt.cm.viridis, extend='both')
    plt.contour(section_lat_vin_mask, np.squeeze(section_z_vin_mask[ee,:,:]), np.squeeze(section_rho_vin[ee,:,:]), colors='k', levels=levelsR, linestyle='dashed') #[34.3,34.46,34.6,34.56,34.6])
#    plt.colorbar(cs, extend='both')
    plt.xlim([-66,-61.6])
    plt.ylim([-600,0])
        
    ax10 = fig.add_subplot(6,3,16)
    ct = plt.contourf(section_lat_mer_mask, np.squeeze(section_z_mer_mask[ee,:,:]), np.squeeze(section_temp_mer[ee,:,:]), levels=levelsT, cmap=plt.cm.coolwarm, extend='both')
    plt.contour(section_lat_mer_mask, np.squeeze(section_z_mer_mask[ee,:,:]), np.squeeze(section_temp_mer[ee,:,:]), colors='k', levels=levelsT) #[-2.3,-2.1,-1.9,-1.5,-1.0,-.5,.0])
#    plt.colorbar(ct, extend='both')
    plt.xlim([-70,-63.7])
    plt.ylim([-600,0])
    plt.title('Merz Glacier ' + date_storms[ee])
    plt.ylabel('Depth (m)')
    plt.xlabel('Latitude')

    ax11 = fig.add_subplot(6,3,17)
    cs = plt.contourf(section_lat_mer_mask, np.squeeze(section_z_mer_mask[ee,:,:]), np.squeeze(section_salt_mer[ee,:,:]), levels=levelsS, linestyle='dash', cmap=plt.cm.coolwarm, extend='both')
    plt.contour(section_lat_mer_mask, np.squeeze(section_z_mer_mask[ee,:,:]), np.squeeze(section_salt_mer[ee,:,:]), colors='k', levels=levelsS, linestyle='dashed') #[34.3,34.46,34.6,34.56,34.6])
#    plt.colorbar(cs, extend='both')
    plt.xlim([-70,-63.7])
    plt.ylim([-600,0])
    plt.xlabel('Latitude')

    ax12 = fig.add_subplot(6,3,18)
    cr = plt.contourf(section_lat_mer_mask, np.squeeze(section_z_mer_mask[ee,:,:]), np.squeeze(section_rho_mer[ee,:,:]), levels=levelsR, linestyle='dash', cmap=plt.cm.viridis, extend='both')
    plt.contour(section_lat_mer_mask, np.squeeze(section_z_mer_mask[ee,:,:]), np.squeeze(section_rho_mer[ee,:,:]), colors='k', levels=levelsR, linestyle='dashed') #[34.3,34.46,34.6,34.56,34.6])
#    plt.colorbar(cs, extend='both')
    plt.xlim([-70,-63.7])
    plt.ylim([-600,0])
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
    
    name_fig="waom10extend_All_EastAntarct_sections_event=" + date_storms[ee] + ".png"
    plt.savefig(fig_path + name_fig)

