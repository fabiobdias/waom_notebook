# plot sections monthly:

import xarray as xr
import pandas as p
import numpy as np
import numpy.ma as ma
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

import gsw

# index for eta var sections (wed, ros, neu, mai, mer)
xi_pt = [425, 775, 700, 875, 1075]
eta_sec_ini = [700, 125, 1050, 1050, 25]
eta_sec_end = [1075, 500, 1300, 1300, 200]
# index for xi var sections (pry, dar, amu)
xi_sec_ini = [1050, 1050, 125]
xi_sec_end = [1500, 1500, 500]
eta_pt = [845, 875, 512]
# length of the section
sec_len = [375, 375, 450, 250, 250, 450, 175, 375]

def read_section_roms(sec_name,ind_ij, ind_len): # 'wed',0,0
#if sec_name == 'wed' or sec_name ==  'ros' or sec_name ==  'neu' or sec_name ==  'mai' or sec_name ==  'mer':
#elif sec_name == 'pry' or sec_name ==  'dar' or sec_name ==  'amu':

    # allocate sections arrays:
    section_salt = np.empty((12,31,sec_len[ind_len]))
    section_temp = np.empty((12,31,sec_len[ind_len]))
    section_rho = np.empty((12,31,sec_len[ind_len]))
    section_z = np.empty((12,31,sec_len[ind_len]))
    section_z_mask = np.empty((12,31,sec_len[ind_len]))

# load ROMS avg output
    count = 0
    for mm in ['01','02','03','04','05','06','07','08','09','10','11','12']:
        ds = xr.open_dataset("/scratch/project_2000339/boeiradi/waom4_mahti/output_6yr_fabio_inifile/ocean_avg_00" + mm + ".nc")
        print('size temp and time length: ', ds.temp.shape, len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
        if sec_name == 'wed' or sec_name ==  'ros' or sec_name ==  'neu' or sec_name ==  'mai' or sec_name ==  'mer':
            print(sec_name, 'xi =', xi_pt[ind_ij], 'eta =', eta_sec_ini[ind_ij], eta_sec_end[ind_ij])
            section_salt_avg = np.divide(np.nansum(ds.salt.isel(xi_rho=xi_pt[ind_ij], eta_rho=slice(eta_sec_ini[ind_ij],eta_sec_end[ind_ij])), axis=0), len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
            section_temp_avg = np.divide(np.nansum(ds.temp.isel(xi_rho=xi_pt[ind_ij], eta_rho=slice(eta_sec_ini[ind_ij],eta_sec_end[ind_ij])), axis=0), len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
            section_zeta = ds.zeta.isel(xi_rho=xi_pt[ind_ij], eta_rho=slice(eta_sec_ini[ind_ij],eta_sec_end[ind_ij]))
            dens_avg = gsw.rho(section_salt_avg,section_temp_avg,2000) - 1000
        elif sec_name == 'pry' or sec_name ==  'dar' or sec_name ==  'amu':
            print(sec_name, 'eta =', eta_pt[ind_ij], 'xi =', xi_sec_ini[ind_ij], xi_sec_end[ind_ij])
            section_salt_avg = np.divide(np.nansum(ds.salt.isel(xi_rho=slice(xi_sec_ini[ind_ij],xi_sec_end[ind_ij]), eta_rho=eta_pt[ind_ij]), axis=0), len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
            section_temp_avg = np.divide(np.nansum(ds.temp.isel(xi_rho=slice(xi_sec_ini[ind_ij],xi_sec_end[ind_ij]), eta_rho=eta_pt[ind_ij]), axis=0), len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
            section_zeta = ds.zeta.isel(xi_rho=slice(xi_sec_ini[ind_ij],xi_sec_end[ind_ij]), eta_rho=eta_pt[ind_ij])
            dens_avg = gsw.rho(section_salt_avg,section_temp_avg,2000) - 1000

        if mm == '01':
            ds = ds.set_coords(['Cs_r', 'Cs_w', 'hc', 'Vtransform'])
            if sec_name == 'wed' or sec_name ==  'ros' or sec_name ==  'neu' or sec_name ==  'mai' or sec_name ==  'mer':
                h = ds.h.isel(xi_rho=xi_pt[ind_ij], eta_rho=slice(eta_sec_ini[ind_ij],eta_sec_end[ind_ij]))   
            elif sec_name == 'pry' or sec_name ==  'dar' or sec_name ==  'amu':
                h = ds.h.isel(xi_rho=slice(xi_sec_ini[ind_ij],xi_sec_end[ind_ij]), eta_rho=eta_pt[ind_ij])
 
        print("Vtransform=2")
        Zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * h) / (ds.hc + h)
        print('Zo_rho, zeta, and h sizes: ', Zo_rho.shape, section_zeta.shape, h.shape)
        z_rho_tmp = section_zeta + (section_zeta + h) * Zo_rho
        print('z_rho_tmp size: ', z_rho_tmp.shape)
        del section_zeta, Zo_rho
        z_rho_tmp_avg = np.divide(np.nansum(z_rho_tmp, axis=0), len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
        print('temp/salt/z_rho_tmp_avg size = ', section_temp_avg.shape, section_salt_avg.shape, z_rho_tmp_avg.shape)
        # assign to monthly sections:
        section_salt[count,:] = section_salt_avg
        section_temp[count,:] = section_temp_avg
        z_rho_tmp_avg = z_rho_tmp_avg.transpose(1,0)
        section_z[count,:] = z_rho_tmp_avg
        section_rho[count,:] = dens_avg
        # mask
        section_z_mask[count,:] = ma.array(z_rho_tmp_avg, mask=np.isnan(z_rho_tmp_avg))
        del section_salt_avg, section_temp_avg, z_rho_tmp_avg, dens_avg, z_rho_tmp
        count = count+1

    return section_salt, section_temp, section_z, section_rho, section_z_mask

section_salt_wed, section_temp_wed, section_z_wed, section_rho_wed, section_z_wed_mask = read_section_roms('wed',0, 0)
section_salt_ros, section_temp_ros, section_z_ros, section_rho_ros, section_z_ros_mask = read_section_roms('ros',1, 1)
section_salt_pry, section_temp_pry, section_z_pry, section_rho_pry, section_z_pry_mask = read_section_roms('pry',0, 2)
section_salt_neu, section_temp_neu, section_z_neu, section_rho_neu, section_z_neu_mask = read_section_roms('neu',2, 3)
section_salt_mai, section_temp_mai, section_z_mai, section_rho_mai, section_z_mai_mask = read_section_roms('mai',3, 4)
section_salt_dar, section_temp_dar, section_z_dar, section_rho_dar, section_z_dar_mask = read_section_roms('dar',1, 5)
section_salt_mer, section_temp_mer, section_z_mer, section_rho_mer, section_z_mer_mask = read_section_roms('mer',4, 6)
section_salt_amu, section_temp_amu, section_z_amu, section_rho_amu, section_z_amu_mask = read_section_roms('amu',2, 7)

dg = xr.open_dataset("/scratch/project_2000339/boeiradi/waom4_frc/waom4_grd.nc")

lat_rho = dg.variables["lat_rho"]
lon_rho = dg.variables["lon_rho"]
#ds.coords['lat_rho']=lat_rho.transpose() # put lat_rho into ds dataset
#ds.coords['lon_rho']=lon_rho.transpose() # put lon_rho into ds dataset

print('lon/lat shape =', lon_rho.shape, lat_rho.shape)

# Weddell Sea
# 2 - slide section for lon or lat
section_lat_wed_tmp = dg.lat_rho.isel(xi_rho=425, eta_rho=slice(700,1075))
# 2.1 - mask land values for lon/lat needs loop for repeat through vert layers
section_lat_wed = np.ones((31,375))                                                
for ii in np.arange(0,31):
    section_lat_wed[ii,:] = section_lat_wed_tmp
section_lat_wed_mask = ma.array(section_lat_wed,mask=np.isnan(section_z_wed[0,:,:]))
print(section_lat_wed.shape, section_z_wed.shape, section_salt_wed.shape)
# annual mean
section_z_wed_mask_ann = np.squeeze(np.nanmean(section_z_wed_mask[:,:,:], axis=0))
section_temp_wed_ann = np.squeeze(np.nanmean(section_temp_wed[:,:,:], axis=0))
section_salt_wed_ann = np.squeeze(np.nanmean(section_salt_wed[:,:,:], axis=0))
section_rho_wed_ann = np.squeeze(np.nanmean(section_rho_wed[:,:,:], axis=0))

# Ross Sea
# 2 - slide section for lon or lat
section_lat_ros_tmp = dg.lat_rho.isel(xi_rho=775, eta_rho=slice(125,500))
# 2.1 - mask land values for lon/lat needs loop for repeat through vert layers
section_lat_ros = np.ones((31,375))                                                
for ii in np.arange(0,31):
    section_lat_ros[ii,:] = section_lat_ros_tmp
section_lat_ros_mask = ma.array(section_lat_ros,mask=np.isnan(section_z_ros[0,:,:]))
print(section_lat_ros.shape, section_z_ros.shape, section_salt_ros.shape)
# annual mean
section_z_ros_mask_ann = np.squeeze(np.nanmean(section_z_ros_mask[:,:,:], axis=0))
section_temp_ros_ann = np.squeeze(np.nanmean(section_temp_ros[:,:,:], axis=0))
section_salt_ros_ann = np.squeeze(np.nanmean(section_salt_ros[:,:,:], axis=0))
section_rho_ros_ann = np.squeeze(np.nanmean(section_rho_ros[:,:,:], axis=0))

# Prydz Bay
# 2 - slide section for lon or lat
section_lon_pry_tmp = dg.lon_rho.isel(xi_rho=slice(1050,1500), eta_rho=845)
# 2.1 - mask land values for lon/lat needs loop for repeat through vert layers
section_lon_pry = np.ones((31,450))                                                
for ii in np.arange(0,31):
    section_lon_pry[ii,:] = section_lon_pry_tmp
section_lon_pry_mask = ma.array(section_lon_pry,mask=np.isnan(section_z_pry[0,:,:]))
print(section_lon_pry.shape, section_z_pry.shape, section_salt_pry.shape)
# annual mean
section_z_pry_mask_ann = np.squeeze(np.nanmean(section_z_pry_mask[:,:,:], axis=0))
section_temp_pry_ann = np.squeeze(np.nanmean(section_temp_pry[:,:,:], axis=0))
section_salt_pry_ann = np.squeeze(np.nanmean(section_salt_pry[:,:,:], axis=0))
section_rho_pry_ann = np.squeeze(np.nanmean(section_rho_pry[:,:,:], axis=0))

# Maud Land 1:
# 2 - slide section for lon or lat
section_lat_neu_tmp = dg.lat_rho.isel(xi_rho=700, eta_rho=slice(1050,1300))
# 2.1 - mask land values for lon/lat needs loop for repeat through vert layers
section_lat_neu = np.ones((31,250))                                                
for ii in np.arange(0,31):
    section_lat_neu[ii,:] = section_lat_neu_tmp
section_lat_neu_mask = ma.array(section_lat_neu,mask=np.isnan(section_z_neu[0,:,:]))
print(section_lat_neu.shape, section_z_neu.shape, section_salt_neu.shape)
# annual mean
section_z_neu_mask_ann = np.squeeze(np.nanmean(section_z_neu_mask[:,:,:], axis=0))
section_temp_neu_ann = np.squeeze(np.nanmean(section_temp_neu[:,:,:], axis=0))
section_salt_neu_ann = np.squeeze(np.nanmean(section_salt_neu[:,:,:], axis=0))
section_rho_neu_ann = np.squeeze(np.nanmean(section_rho_neu[:,:,:], axis=0))

# Maud Land 2:
# 2 - slide section for lon or lat
section_lat_mai_tmp = dg.lat_rho.isel(xi_rho=875, eta_rho=slice(1050,1300))
# 2.1 - mask land values for lon/lat needs loop for repeat through vert layers
section_lat_mai = np.ones((31,250))                                                
for ii in np.arange(0,31):
    section_lat_mai[ii,:] = section_lat_mai_tmp
section_lat_mai_mask = ma.array(section_lat_mai,mask=np.isnan(section_z_mai[0,:,:]))
print(section_lat_mai.shape, section_z_mai.shape, section_salt_mai.shape)
# annual mean
section_z_mai_mask_ann = np.squeeze(np.nanmean(section_z_mai_mask[:,:,:], axis=0))
section_temp_mai_ann = np.squeeze(np.nanmean(section_temp_mai[:,:,:], axis=0))
section_salt_mai_ann = np.squeeze(np.nanmean(section_salt_mai[:,:,:], axis=0))
section_rho_mai_ann = np.squeeze(np.nanmean(section_rho_mai[:,:,:], axis=0))

# Darnley Polynya
# 2 - slide section for lon or lat
section_lon_dar_tmp = dg.lon_rho.isel(xi_rho=slice(1050,1500), eta_rho=875)
# 2.1 - mask land values for lon/lat needs loop for repeat through vert layers
section_lon_dar = np.ones((31,450))                                                
for ii in np.arange(0,31):
    section_lon_dar[ii,:] = section_lon_dar_tmp
section_lon_dar_mask = ma.array(section_lon_dar,mask=np.isnan(section_z_dar[0,:,:]))
print(section_lon_dar.shape, section_z_dar.shape, section_salt_dar.shape)
# annual mean
section_z_dar_mask_ann = np.squeeze(np.nanmean(section_z_dar_mask[:,:,:], axis=0))
section_temp_dar_ann = np.squeeze(np.nanmean(section_temp_dar[:,:,:], axis=0))
section_salt_dar_ann = np.squeeze(np.nanmean(section_salt_dar[:,:,:], axis=0))
section_rho_dar_ann = np.squeeze(np.nanmean(section_rho_dar[:,:,:], axis=0))

# Merz Glacier
# 2 - slide section for lon or lat
section_lat_mer_tmp = dg.lat_rho.isel(xi_rho=1075, eta_rho=slice(25,200))
# 2.1 - mask land values for lon/lat needs loop for repeat through vert layers
section_lat_mer = np.ones((31,175))                                                
for ii in np.arange(0,31):
    section_lat_mer[ii,:] = section_lat_mer_tmp
section_lat_mer_mask = ma.array(section_lat_mer,mask=np.isnan(section_z_mer[0,:,:]))
print(section_lat_mer.shape, section_z_mer.shape, section_salt_mer.shape)
# annnual mean
section_z_mer_mask_ann = np.squeeze(np.nanmean(section_z_mer_mask[:,:,:], axis=0))
section_temp_mer_ann = np.squeeze(np.nanmean(section_temp_mer[:,:,:], axis=0))
section_salt_mer_ann = np.squeeze(np.nanmean(section_salt_mer[:,:,:], axis=0))
section_rho_mer_ann = np.squeeze(np.nanmean(section_rho_mer[:,:,:], axis=0))

# Amundsen Sea embayment
# 2 - slide section for lon or lat
section_lon_amu_tmp = dg.lon_rho.isel(xi_rho=slice(125,500), eta_rho=512)
# 2.1 - mask land values for lon/lat needs loop for repeat through vert layers
section_lon_amu = np.ones((31,375))                                                
for ii in np.arange(0,31):
    section_lon_amu[ii,:] = section_lon_amu_tmp
section_lon_amu_mask = ma.array(section_lon_amu,mask=np.isnan(section_z_amu[0,:,:]))
print(section_lon_amu.shape, section_z_amu.shape, section_salt_amu.shape)
# annual mean
section_z_amu_mask_ann = np.squeeze(np.nanmean(section_z_amu_mask[:,:,:], axis=0))
section_temp_amu_ann = np.squeeze(np.nanmean(section_temp_amu[:,:,:], axis=0))
section_salt_amu_ann = np.squeeze(np.nanmean(section_salt_amu[:,:,:], axis=0))
section_rho_amu_ann = np.squeeze(np.nanmean(section_rho_amu[:,:,:], axis=0))

# Plot transects: Annual
levelsT = np.arange(-2.5,2.5,.25)
levelsS = np.arange(33.1,35.2,.1)
levelsR = np.arange(36.4,37.2,.1)

# Amundsen Sea
fig = plt.figure(figsize=(16,18))

ax1 = fig.add_subplot(4,2,1)
ct = plt.contourf(section_lon_amu_mask, section_z_amu_mask_ann, section_temp_amu_ann, cmap=plt.cm.coolwarm, levels=levelsT)
plt.contour(section_lon_amu_mask, section_z_amu_mask_ann, section_temp_amu_ann, colors='k', levels=levelsT)
plt.contour(section_lon_amu_mask, section_z_amu_mask_ann, section_rho_amu_ann, colors='white', levels=[36.9],linewidth=2)
plt.xlim([-115,-106])
plt.ylim([-4200,0])
plt.title('Amundsen Sea',fontsize=12)
plt.ylabel('Depth (m)',fontsize=12)

# Weddell
ax2 = fig.add_subplot(4,2,2)
ct = plt.contourf(section_lat_wed_mask, section_z_wed_mask_ann, section_temp_wed_ann, cmap=plt.cm.coolwarm, levels=levelsT)
plt.contour(section_lat_wed_mask, section_z_wed_mask_ann, section_temp_wed_ann, colors='k', levels=levelsT)
plt.contour(section_lat_wed_mask, section_z_wed_mask_ann, section_rho_wed_ann, colors='white', levels=[36.9],linewidth=2)
plt.xlim([-81,-71])
plt.ylim([-4200,0])
plt.title('Weddell Sea',fontsize=12)

# Maud Land 1
ax3 = fig.add_subplot(4,2,4)
ct = plt.contourf(section_lat_neu_mask, section_z_neu_mask_ann, section_temp_neu_ann, cmap=plt.cm.coolwarm, levels=levelsT)
plt.contour(section_lat_neu_mask, section_z_neu_mask_ann, section_temp_neu_ann, colors='k', levels=levelsT)
plt.contour(section_lat_neu_mask, section_z_neu_mask_ann, section_rho_neu_ann, colors='white', levels=[36.9],linewidth=2)
plt.xlim([-72,-67])
plt.ylim([-4200,0])
plt.title('Maud Land 1',fontsize=12)

# Maud Land 2
ax4 = fig.add_subplot(4,2,6)
ct = plt.contourf(section_lat_mai_mask, section_z_mai_mask_ann, section_temp_mai_ann, cmap=plt.cm.coolwarm, levels=levelsT)
plt.contour(section_lat_mai_mask, section_z_mai_mask_ann, section_temp_mai_ann, colors='k', levels=levelsT)
plt.contour(section_lat_mai_mask, section_z_mai_mask_ann, section_rho_mai_ann, colors='white', levels=[36.9],linewidth=2)
plt.xlim([-72,-67])
plt.ylim([-4200,0])
plt.title('Maud Land 2',fontsize=12)

# Ross
ax5 = fig.add_subplot(4,2,3)
ct = plt.contourf(section_lat_ros_mask, section_z_ros_mask_ann, section_temp_ros_ann, cmap=plt.cm.coolwarm, levels=levelsT)
plt.contour(section_lat_ros_mask, section_z_ros_mask_ann, section_temp_ros_ann, colors='k', levels=levelsT)
plt.contour(section_lat_ros_mask, section_z_ros_mask_ann, section_rho_ros_ann, colors='white', levels=[36.9],linewidth=2)
plt.xlim([-86,-69.5])
plt.ylim([-4200,0])
plt.title('Ross Sea',fontsize=12)
plt.ylabel('Depth (m)',fontsize=12)

# Merz
ax6 = fig.add_subplot(4,2,8)
ct = plt.contourf(section_lat_mer_mask, section_z_mer_mask_ann, section_temp_mer_ann, cmap=plt.cm.coolwarm, levels=levelsT)
plt.contour(section_lat_mer_mask, section_z_mer_mask_ann, section_temp_mer_ann, colors='k', levels=levelsT)
plt.contour(section_lat_mer_mask, section_z_mer_mask_ann, section_rho_mer_ann, colors='white', levels=[36.9],linewidth=2)
plt.xlim([-70,-63.7])
plt.ylim([-4200,0])
plt.title('Merz Glacier',fontsize=12)
#plt.ylabel('Depth (m)')
plt.xlabel('Latitude',fontsize=12)

# prydz bay
ax8 = fig.add_subplot(4,2,7)
ct = plt.contourf(section_lon_pry_mask, section_z_pry_mask_ann, section_temp_pry_ann, cmap=plt.cm.coolwarm, levels=levelsT)
plt.contour(section_lon_pry_mask, section_z_pry_mask_ann, section_temp_pry_ann, colors='k', levels=levelsT)
plt.contour(section_lon_pry_mask, section_z_pry_mask_ann, section_rho_pry_ann, colors='white', levels=[36.9],linewidth=2)
plt.xlim([69,75.5])
plt.ylim([-4200,0])
plt.title('Prydz Bay')
plt.ylabel('Depth (m)',fontsize=12)
plt.xlabel('Longitude',fontsize=12)

# Darnley Polynya
ax7 = fig.add_subplot(4,2,5)
ct = plt.contourf(section_lon_dar_mask, section_z_dar_mask_ann, section_temp_dar_ann, cmap=plt.cm.coolwarm, levels=levelsT)
plt.contour(section_lon_dar_mask, section_z_dar_mask_ann, section_temp_dar_ann, colors='k', levels=levelsT)
plt.xlim([69,75])
plt.ylim([-4200,0])
plt.title('Darnley Polynya',fontsize=12)
plt.ylabel('Depth (m)',fontsize=12)

cbar_axim = fig.add_axes([0.92, 0.15, 0.015, 0.7])
cbar =fig.colorbar(ct, cax=cbar_axim, extend='both', format='%.1f')
cbar.ax.tick_params(labelsize=14)
plt.ylabel('Temperature ($^{\circ}$C)',fontsize=14)

#plt.tight_layout()
plt.savefig('figs/WAOM4/waom4_all_sections_temp_annual_new.png')
plt.close(fig)

# Plot transects: Annual
# Amundsen Sea
fig = plt.figure(figsize=(16,18))

ax1 = fig.add_subplot(4,2,1)
ct = plt.contourf(section_lon_amu_mask, section_z_amu_mask_ann, section_salt_amu_ann, cmap=plt.cm.coolwarm, levels=levelsS)
plt.contour(section_lon_amu_mask, section_z_amu_mask_ann, section_salt_amu_ann, colors='k', levels=levelsS)
plt.contour(section_lon_amu_mask, section_z_amu_mask_ann, section_rho_amu_ann, colors='white', levels=[36.9],linewidth=2)
plt.xlim([-115,-106])
plt.ylim([-4200,0])
plt.title('Amundsen Sea',fontsize=12)
plt.ylabel('Depth (m)',fontsize=12)

# Weddell
ax2 = fig.add_subplot(4,2,2)
ct = plt.contourf(section_lat_wed_mask, section_z_wed_mask_ann, section_salt_wed_ann, cmap=plt.cm.coolwarm, levels=levelsS)
plt.contour(section_lat_wed_mask, section_z_wed_mask_ann, section_salt_wed_ann, colors='k', levels=levelsS)
plt.contour(section_lat_wed_mask, section_z_wed_mask_ann, section_rho_wed_ann, colors='white', levels=[36.9],linewidth=2)
plt.xlim([-81,-71])
plt.ylim([-4200,0])
plt.title('Weddell Sea',fontsize=12)

# Maud Land 1
ax3 = fig.add_subplot(4,2,4)
ct = plt.contourf(section_lat_neu_mask, section_z_neu_mask_ann, section_salt_neu_ann, cmap=plt.cm.coolwarm, levels=levelsS)
plt.contour(section_lat_neu_mask, section_z_neu_mask_ann, section_salt_neu_ann, colors='k', levels=levelsS)
plt.contour(section_lat_neu_mask, section_z_neu_mask_ann, section_rho_neu_ann, colors='white', levels=[36.9],linewidth=2)
plt.xlim([-72,-67])
plt.ylim([-4200,0])
plt.title('Maud Land 1',fontsize=12)

# Maud Land 2
ax4 = fig.add_subplot(4,2,6)
ct = plt.contourf(section_lat_mai_mask, section_z_mai_mask_ann, section_salt_mai_ann, cmap=plt.cm.coolwarm, levels=levelsS)
plt.contour(section_lat_mai_mask, section_z_mai_mask_ann, section_salt_mai_ann, colors='k', levels=levelsS)
plt.contour(section_lat_mai_mask, section_z_mai_mask_ann, section_rho_mai_ann, colors='white', levels=[36.9],linewidth=2)
plt.xlim([-72,-67])
plt.ylim([-4200,0])
plt.title('Maud Land 2',fontsize=12)

# Ross
ax5 = fig.add_subplot(4,2,3)
ct = plt.contourf(section_lat_ros_mask, section_z_ros_mask_ann, section_salt_ros_ann, cmap=plt.cm.coolwarm, levels=levelsS)
plt.contour(section_lat_ros_mask, section_z_ros_mask_ann, section_salt_ros_ann, colors='k', levels=levelsS)
plt.contour(section_lat_ros_mask, section_z_ros_mask_ann, section_rho_ros_ann, colors='white', levels=[36.9],linewidth=2)
plt.xlim([-86,-69.5])
plt.ylim([-4200,0])
plt.title('Ross Sea',fontsize=12)
plt.ylabel('Depth (m)',fontsize=12)

# Merz
ax6 = fig.add_subplot(4,2,8)
ct = plt.contourf(section_lat_mer_mask, section_z_mer_mask_ann, section_salt_mer_ann, cmap=plt.cm.coolwarm, levels=levelsS)
plt.contour(section_lat_mer_mask, section_z_mer_mask_ann, section_salt_mer_ann, colors='k', levels=levelsS)
plt.contour(section_lat_mer_mask, section_z_mer_mask_ann, section_rho_mer_ann, colors='white', levels=[36.9],linewidth=2)
plt.xlim([-70,-63.7])
plt.ylim([-4200,0])
plt.title('Merz Glacier',fontsize=12)
#plt.ylabel('Depth (m)')
plt.xlabel('Latitude',fontsize=12)

# prydz bay
ax8 = fig.add_subplot(4,2,7)
ct = plt.contourf(section_lon_pry_mask, section_z_pry_mask_ann, section_salt_pry_ann, cmap=plt.cm.coolwarm, levels=levelsS)
plt.contour(section_lon_pry_mask, section_z_pry_mask_ann, section_salt_pry_ann, colors='k', levels=levelsS)
plt.contour(section_lon_pry_mask, section_z_pry_mask_ann, section_rho_pry_ann, colors='white', levels=[36.9],linewidth=2)
plt.xlim([69,75.5])
plt.ylim([-4200,0])
plt.title('Prydz Bay')
plt.ylabel('Depth (m)',fontsize=12)
plt.xlabel('Longitude',fontsize=12)

# Darnley Polynya
ax7 = fig.add_subplot(4,2,5)
ct = plt.contourf(section_lon_dar_mask, section_z_dar_mask_ann, section_salt_dar_ann, cmap=plt.cm.coolwarm, levels=levelsS)
plt.contour(section_lon_dar_mask, section_z_dar_mask_ann, section_salt_dar_ann, colors='k', levels=levelsS)
plt.contour(section_lon_dar_mask, section_z_dar_mask_ann, section_rho_dar_ann, colors='white', levels=[36.9],linewidth=2)
plt.xlim([69,75])
plt.ylim([-4200,0])
plt.title('Darnley Polynya',fontsize=12)
plt.ylabel('Depth (m)',fontsize=12)

cbar_axim = fig.add_axes([0.92, 0.15, 0.015, 0.7])
cbar =fig.colorbar(ct, cax=cbar_axim, extend='both', format='%.1f')
cbar.ax.tick_params(labelsize=14)
plt.ylabel('Salinity',fontsize=14)

#plt.tight_layout()
plt.savefig('figs/WAOM4/waom4_all_sections_salt_annual_new.png')
plt.close(fig)

# Plot transects: Annual
# Amundsen Sea
fig = plt.figure(figsize=(16,18))

ax1 = fig.add_subplot(4,2,1)
ct = plt.contourf(section_lon_amu_mask, section_z_amu_mask_ann, section_rho_amu_ann, cmap=plt.cm.viridis, levels=levelsR)
plt.contour(section_lon_amu_mask, section_z_amu_mask_ann, section_rho_amu_ann, colors='k', levels=levelsR)
plt.contour(section_lon_amu_mask, section_z_amu_mask_ann, section_rho_amu_ann, colors='white', levels=[36.9],linewidth=2)
plt.xlim([-115,-106])
plt.ylim([-4200,0])
plt.title('Amundsen Sea',fontsize=12)
plt.ylabel('Depth (m)',fontsize=12)

# Weddell
ax2 = fig.add_subplot(4,2,2)
ct = plt.contourf(section_lat_wed_mask, section_z_wed_mask_ann, section_rho_wed_ann, cmap=plt.cm.viridis, levels=levelsR)
plt.contour(section_lat_wed_mask, section_z_wed_mask_ann, section_rho_wed_ann, colors='k', levels=levelsR)
plt.contour(section_lat_wed_mask, section_z_wed_mask_ann, section_rho_wed_ann, colors='white', levels=[36.9],linewidth=2)
plt.xlim([-81,-71])
plt.ylim([-4200,0])
plt.title('Weddell Sea',fontsize=12)

# Maud Land 1
ax3 = fig.add_subplot(4,2,4)
ct = plt.contourf(section_lat_neu_mask, section_z_neu_mask_ann, section_rho_neu_ann, cmap=plt.cm.viridis, levels=levelsR)
plt.contour(section_lat_neu_mask, section_z_neu_mask_ann, section_rho_neu_ann, colors='k', levels=levelsR)
plt.contour(section_lat_neu_mask, section_z_neu_mask_ann, section_rho_neu_ann, colors='white', levels=[36.9],linewidth=2)
plt.xlim([-72,-67])
plt.ylim([-4200,0])
plt.title('Maud Land 1',fontsize=12)

# Maud Land 2
ax4 = fig.add_subplot(4,2,6)
ct = plt.contourf(section_lat_mai_mask, section_z_mai_mask_ann, section_rho_mai_ann, cmap=plt.cm.viridis, levels=levelsR)
plt.contour(section_lat_mai_mask, section_z_mai_mask_ann, section_rho_mai_ann, colors='k', levels=levelsR)
plt.contour(section_lat_mai_mask, section_z_mai_mask_ann, section_rho_mai_ann, colors='white', levels=[36.9],linewidth=2)
plt.xlim([-72,-67])
plt.ylim([-4200,0])
plt.title('Maud Land 2',fontsize=12)

# Ross
ax5 = fig.add_subplot(4,2,3)
ct = plt.contourf(section_lat_ros_mask, section_z_ros_mask_ann, section_rho_ros_ann, cmap=plt.cm.viridis, levels=levelsR)
plt.contour(section_lat_ros_mask, section_z_ros_mask_ann, section_rho_ros_ann, colors='k', levels=levelsR)
plt.contour(section_lat_ros_mask, section_z_ros_mask_ann, section_rho_ros_ann, colors='white', levels=[36.9],linewidth=2)
plt.xlim([-86,-69.5])
plt.ylim([-4200,0])
plt.title('Ross Sea',fontsize=12)
plt.ylabel('Depth (m)',fontsize=12)

# Merz
ax6 = fig.add_subplot(4,2,8)
ct = plt.contourf(section_lat_mer_mask, section_z_mer_mask_ann, section_rho_mer_ann, cmap=plt.cm.viridis, levels=levelsR)
plt.contour(section_lat_mer_mask, section_z_mer_mask_ann, section_rho_mer_ann, colors='k', levels=levelsR)
plt.contour(section_lat_mer_mask, section_z_mer_mask_ann, section_rho_mer_ann, colors='white', levels=[36.9],linewidth=2)
plt.xlim([-70,-63.7])
plt.ylim([-4200,0])
plt.title('Merz Glacier',fontsize=12)
#plt.ylabel('Depth (m)')
plt.xlabel('Latitude',fontsize=12)

# prydz bay
ax8 = fig.add_subplot(4,2,7)
ct = plt.contourf(section_lon_pry_mask, section_z_pry_mask_ann, section_rho_pry_ann, cmap=plt.cm.viridis, levels=levelsR)
plt.contour(section_lon_pry_mask, section_z_pry_mask_ann, section_rho_pry_ann, colors='k', levels=levelsR)
plt.contour(section_lon_pry_mask, section_z_pry_mask_ann, section_rho_pry_ann, colors='white', levels=[36.9],linewidth=2)
plt.xlim([69,75.5])
plt.ylim([-4200,0])
plt.title('Prydz Bay')
plt.ylabel('Depth (m)',fontsize=12)
plt.xlabel('Longitude',fontsize=12)

# Darnley Polynya
ax7 = fig.add_subplot(4,2,5)
ct = plt.contourf(section_lon_dar_mask, section_z_dar_mask_ann, section_rho_dar_ann, cmap=plt.cm.viridis, levels=levelsR)
plt.contour(section_lon_dar_mask, section_z_dar_mask_ann, section_rho_dar_ann, colors='k', levels=levelsR)
plt.contour(section_lon_dar_mask, section_z_dar_mask_ann, section_rho_dar_ann, colors='white', levels=[36.9],linewidth=2)
plt.xlim([69,75])
plt.ylim([-4200,0])
plt.title('Darnley Polynya',fontsize=12)
plt.ylabel('Depth (m)',fontsize=12)

cbar_axim = fig.add_axes([0.92, 0.15, 0.015, 0.7])
cbar =fig.colorbar(ct, cax=cbar_axim, extend='both', format='%.1f')
cbar.ax.tick_params(labelsize=14)
plt.ylabel('Potential density ($\sigma_2$)',fontsize=14)

#plt.tight_layout()
plt.savefig('figs/WAOM4/waom4_all_sections_rho_annual_new.png', dpi=300)
plt.close(fig)

# Plot transects: Annual
fig = plt.figure(figsize=(16,18))

# Weddell
ax2 = fig.add_subplot(2,2,1)
ct = plt.contourf(section_lat_wed_mask, section_z_wed_mask_ann, section_rho_wed_ann, cmap=plt.cm.viridis, levels=levelsR)
plt.contour(section_lat_wed_mask, section_z_wed_mask_ann, section_rho_wed_ann, colors='k', levels=levelsR)
plt.contour(section_lat_wed_mask, section_z_wed_mask_ann, section_rho_wed_ann, colors='white', levels=[36.9],linewidth=2)
plt.xlim([-81,-71])
plt.ylim([-4200,0])
plt.title('Weddell Sea',fontsize=12)
plt.xlabel('Latitude',fontsize=12)
plt.ylabel('Depth (m)',fontsize=12)

# Ross
ax5 = fig.add_subplot(2,2,2)
ct = plt.contourf(section_lat_ros_mask, section_z_ros_mask_ann, section_rho_ros_ann, cmap=plt.cm.viridis, levels=levelsR)
plt.contour(section_lat_ros_mask, section_z_ros_mask_ann, section_rho_ros_ann, colors='k', levels=levelsR)
plt.contour(section_lat_ros_mask, section_z_ros_mask_ann, section_rho_ros_ann, colors='white', levels=[36.9],linewidth=2)
plt.xlim([-86,-69.5])
plt.ylim([-4200,0])
plt.title('Ross Sea',fontsize=12)
#plt.ylabel('Depth (m)',fontsize=12)
plt.xlabel('Latitude',fontsize=12)

# Merz
ax6 = fig.add_subplot(2,2,3)
ct = plt.contourf(section_lat_mer_mask, section_z_mer_mask_ann, section_rho_mer_ann, cmap=plt.cm.viridis, levels=levelsR)
plt.contour(section_lat_mer_mask, section_z_mer_mask_ann, section_rho_mer_ann, colors='k', levels=levelsR)
plt.contour(section_lat_mer_mask, section_z_mer_mask_ann, section_rho_mer_ann, colors='white', levels=[36.9],linewidth=2)
plt.xlim([-70,-63.7])
plt.ylim([-4200,0])
plt.title('Merz Glacier',fontsize=12)
plt.ylabel('Depth (m)')
plt.xlabel('Latitude',fontsize=12)

# Darnley Polynya
ax7 = fig.add_subplot(2,2,4)
ct = plt.contourf(section_lon_dar_mask, section_z_dar_mask_ann, section_rho_dar_ann, cmap=plt.cm.viridis, levels=levelsR)
plt.contour(section_lon_dar_mask, section_z_dar_mask_ann, section_rho_dar_ann, colors='k', levels=levelsR)
plt.contour(section_lon_dar_mask, section_z_dar_mask_ann, section_rho_dar_ann, colors='white', levels=[36.9],linewidth=2)
plt.xlim([69,75])
plt.ylim([-4200,0])
plt.title('Darnley Polynya',fontsize=12)
#plt.ylabel('Depth (m)',fontsize=12)
plt.xlabel('Longitude',fontsize=12)

cbar_axim = fig.add_axes([0.92, 0.15, 0.015, 0.7])
cbar =fig.colorbar(ct, cax=cbar_axim, extend='both', format='%.1f')
cbar.ax.tick_params(labelsize=14)
plt.ylabel('Potential density ($\sigma_2$)',fontsize=14)

#plt.tight_layout()

plt.savefig('figs/WAOM4/waom4_DSW_sections_rho_annual_new.png')
plt.close(fig)
