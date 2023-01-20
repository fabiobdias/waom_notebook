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

# load ROMS avg output
temp = np.empty((12,31,1325,1575))
salt = np.empty((12,31,1325,1575))
z_rho = np.empty((12,31,1325,1575))
count = 0
for mm in ['01','02','03','04','05','06','07','08','09','10','11','12']:
    ds = xr.open_dataset("/scratch/project_2000339/boeiradi/waom4_mahti/output_6yr_fabio_inifile/ocean_avg_00" + mm + ".nc")
    temp_tmp = ds.variables["temp"]
    salt_tmp = ds.variables["salt"]
    zeta = ds.variables["zeta"]
    print('Reading ROMS avg month = ', mm)
    temp_tmp_avg = np.divide(np.nansum(temp_tmp, axis=0), len(zeta[:,1,1]))
    salt_tmp_avg = np.divide(np.nansum(salt_tmp, axis=0), len(zeta[:,1,1]))
    del temp_tmp, salt_tmp
    if mm == '01':
        ds = ds.set_coords(['Cs_r', 'Cs_w', 'hc', 'h', 'Vtransform'])

    if ds.Vtransform == 1:
        Zo_rho = ds.hc * (ds.s_rho - ds.Cs_r) + ds.Cs_r * ds.h
        z_rho = Zo_rho + ds.zeta * (1 + Zo_rho/ds.h)
        print("Vtransform=1")
    elif ds.Vtransform == 2:
        Zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * ds.h) / (ds.hc + ds.h)
        print('Zo_rho, zeta, zeta_avg, and h sizes: ', Zo_rho.shape, zeta.shape, ds.h.shape)
        z_rho_tmp = zeta + (zeta + ds.h) * Zo_rho
        print("Vtransform=2")
    z_rho_tmp = z_rho_tmp.transpose('ocean_time', 's_rho', 'eta_rho', 'xi_rho')
    z_rho_tmp_avg = np.divide(np.nansum(z_rho_tmp, axis=0), len(zeta[:,1,1]))
    del z_rho_tmp, zeta
    print('temp/salt/z_rho_tmp_avg size = ', temp_tmp_avg.shape, salt_tmp_avg.shape, z_rho_tmp_avg.shape)
    temp[count,:] = temp_tmp_avg
    salt[count,:] = salt_tmp_avg
    z_rho[count,:] = z_rho_tmp_avg
    del temp_tmp_avg, salt_tmp_avg
    del z_rho_tmp_avg
    count = count+1

dens=gsw.rho(salt,temp,2000)
        # Substract 1000 to convert to sigma-t
dens = dens - 1000
print(temp.shape,dens.shape)
m_len = temp[:,0,0,0].size
k_len = temp[0,:,0,0].size
i_len = temp[0,0,:,0].size
j_len = temp[0,0,0,:].size

dg = xr.open_dataset("/scratch/project_2000339/boeiradi/waom4_frc/waom4_grd.nc")

lat_rho = dg.variables["lat_rho"]
lon_rho = dg.variables["lon_rho"]
ds.coords['lat_rho']=lat_rho.transpose() # put lat_rho into ds dataset
ds.coords['lon_rho']=lon_rho.transpose() # put lon_rho into ds dataset

# Weddell Sea
# 1 - slice section for salt/temp/z_rho
section_salt_wed = np.empty((12,31,375)) # allocating
section_temp_wed = np.empty((12,31,375))
section_rho_wed = np.empty((12,31,375))
section_z_wed = np.empty((12,31,375))
section_z_wed_mask = np.empty((12,31,375))
for mm in np.arange(0,12):
    section_salt_wed[mm,:,:] = np.squeeze(salt[mm,:,700:1075,425]) #ds.salt.isel(xi_rho=425, eta_rho=slice(700,1075), ocean_time=mm)
    section_temp_wed[mm,:,:] = np.squeeze(temp[mm,:,700:1075,425]) #ds.temp.isel(xi_rho=425, eta_rho=slice(700,1075), ocean_time=mm)
    section_z_wed[mm,:,:] = np.squeeze(z_rho[mm,:,700:1075,425]) #ds.z_rho.isel(xi_rho=425, eta_rho=slice(700,1075), ocean_time=mm)
    # 1.1 - mask land values in z_rho slice                                           
    section_z_wed_mask[mm,:,:] = ma.array(section_z_wed[mm,:,:],mask=np.isnan(section_z_wed[mm,:,:]))
    section_rho_wed[mm,:,:] = np.squeeze(dens[mm,:,700:1075,425])
    
# 2 - slide section for lon or lat
section_lat_wed_tmp = ds.lat_rho.isel(xi_rho=425, eta_rho=slice(700,1075))
# 2.1 - mask land values for lon/lat needs loop for repeat through vert layers
section_lat_wed = np.ones((31,375))                                                
for ii in np.arange(0,31):
    section_lat_wed[ii,:] = section_lat_wed_tmp
section_lat_wed_mask = ma.array(section_lat_wed,mask=np.isnan(section_z_wed[0,:,:]))

print(section_lat_wed.shape, section_z_wed.shape, section_salt_wed.shape)

## season averages:
#section_z_wed_mask_JFM = np.squeeze(np.nanmean(section_z_wed_mask[0:18,:,:], axis=0))
#section_z_wed_mask_AMJ = np.squeeze(np.nanmean(section_z_wed_mask[19:36,:,:], axis=0))
#section_z_wed_mask_JAS = np.squeeze(np.nanmean(section_z_wed_mask[37:54,:,:], axis=0))
#section_z_wed_mask_OND = np.squeeze(np.nanmean(section_z_wed_mask[55:72,:,:], axis=0))
#
#section_temp_wed_JFM = np.squeeze(np.nanmean(section_temp_wed[0:18,:,:], axis=0))
#section_temp_wed_AMJ = np.squeeze(np.nanmean(section_temp_wed[19:36,:,:], axis=0))
#section_temp_wed_JAS = np.squeeze(np.nanmean(section_temp_wed[37:54,:,:], axis=0))
#section_temp_wed_OND = np.squeeze(np.nanmean(section_temp_wed[55:72,:,:], axis=0))
#
#section_salt_wed_JFM = np.squeeze(np.nanmean(section_salt_wed[0:18,:,:], axis=0))
#section_salt_wed_AMJ = np.squeeze(np.nanmean(section_salt_wed[19:36,:,:], axis=0))
#section_salt_wed_JAS = np.squeeze(np.nanmean(section_salt_wed[37:54,:,:], axis=0))
#section_salt_wed_OND = np.squeeze(np.nanmean(section_salt_wed[55:72,:,:], axis=0))
#
#section_rho_wed_JFM = np.squeeze(np.nanmean(section_rho_wed[0:18,:,:], axis=0))
#section_rho_wed_AMJ = np.squeeze(np.nanmean(section_rho_wed[19:36,:,:], axis=0))
#section_rho_wed_JAS = np.squeeze(np.nanmean(section_rho_wed[37:54,:,:], axis=0))
#section_rho_wed_OND = np.squeeze(np.nanmean(section_rho_wed[55:72,:,:], axis=0))
#
section_z_wed_mask_ann = np.squeeze(np.nanmean(section_z_wed_mask[:,:,:], axis=0))
section_temp_wed_ann = np.squeeze(np.nanmean(section_temp_wed[:,:,:], axis=0))
section_salt_wed_ann = np.squeeze(np.nanmean(section_salt_wed[:,:,:], axis=0))
section_rho_wed_ann = np.squeeze(np.nanmean(section_rho_wed[:,:,:], axis=0))
#
#print(section_lat_wed.shape, section_z_wed_mask_JFM.shape, section_salt_wed_JFM.shape)

# Ross Sea
# 1 - slice section for salt/temp/z_rho
section_salt_ros = np.empty((12,31,375)) # allocating
section_temp_ros = np.empty((12,31,375))
section_rho_ros = np.empty((12,31,375))
section_z_ros = np.empty((12,31,375))
section_z_ros_mask = np.empty((12,31,375))
for mm in np.arange(0,12):
    section_salt_ros[mm,:,:] = np.squeeze(salt[mm,:,125:500,775]) #ds.salt.isel(xi_rho=775, eta_rho=slice(125,500), ocean_time=mm)
    section_temp_ros[mm,:,:] = np.squeeze(temp[mm,:,125:500,775]) #ds.temp.isel(xi_rho=775, eta_rho=slice(125,500), ocean_time=mm)
    section_z_ros[mm,:,:] = np.squeeze(z_rho[mm,:,125:500,775]) #ds.z_rho.isel(xi_rho=775, eta_rho=slice(125,500), ocean_time=mm)
    # 1.1 - mask land values in z_rho slice                                           
    section_z_ros_mask[...] = ma.array(section_z_ros[mm,:,:],mask=np.isnan(section_z_ros[mm,:,:]))
    section_rho_ros[mm,:,:] = np.squeeze(dens[mm,:,125:500,775])

# 2 - slide section for lon or lat
section_lat_ros_tmp = ds.lat_rho.isel(xi_rho=775, eta_rho=slice(125,500))
# 2.1 - mask land values for lon/lat needs loop for repeat through vert layers
section_lat_ros = np.ones((31,375))                                                
for ii in np.arange(0,31):
    section_lat_ros[ii,:] = section_lat_ros_tmp
section_lat_ros_mask = ma.array(section_lat_ros,mask=np.isnan(section_z_ros[0,:,:]))

print(section_lat_ros.shape, section_z_ros.shape, section_salt_ros.shape)

## season averages:
#section_z_ros_mask_JFM = np.squeeze(np.nanmean(section_z_ros_mask[0:18,:,:], axis=0))
#section_z_ros_mask_AMJ = np.squeeze(np.nanmean(section_z_ros_mask[19:36,:,:], axis=0))
#section_z_ros_mask_JAS = np.squeeze(np.nanmean(section_z_ros_mask[37:54,:,:], axis=0))
#section_z_ros_mask_OND = np.squeeze(np.nanmean(section_z_ros_mask[55:72,:,:], axis=0))
#
#section_temp_ros_JFM = np.squeeze(np.nanmean(section_temp_ros[0:18,:,:], axis=0))
#section_temp_ros_AMJ = np.squeeze(np.nanmean(section_temp_ros[19:36,:,:], axis=0))
#section_temp_ros_JAS = np.squeeze(np.nanmean(section_temp_ros[37:54,:,:], axis=0))
#section_temp_ros_OND = np.squeeze(np.nanmean(section_temp_ros[55:72,:,:], axis=0))
#
#section_salt_ros_JFM = np.squeeze(np.nanmean(section_salt_ros[0:18,:,:], axis=0))
#section_salt_ros_AMJ = np.squeeze(np.nanmean(section_salt_ros[19:36,:,:], axis=0))
#section_salt_ros_JAS = np.squeeze(np.nanmean(section_salt_ros[37:54,:,:], axis=0))
#section_salt_ros_OND = np.squeeze(np.nanmean(section_salt_ros[55:72,:,:], axis=0))
#
#section_rho_ros_JFM = np.squeeze(np.nanmean(section_rho_ros[0:18,:,:], axis=0))
#section_rho_ros_AMJ = np.squeeze(np.nanmean(section_rho_ros[19:36,:,:], axis=0))
#section_rho_ros_JAS = np.squeeze(np.nanmean(section_rho_ros[37:54,:,:], axis=0))
#section_rho_ros_OND = np.squeeze(np.nanmean(section_rho_ros[55:72,:,:], axis=0))
#
section_z_ros_mask_ann = np.squeeze(np.nanmean(section_z_ros_mask[:,:,:], axis=0))
section_temp_ros_ann = np.squeeze(np.nanmean(section_temp_ros[:,:,:], axis=0))
section_salt_ros_ann = np.squeeze(np.nanmean(section_salt_ros[:,:,:], axis=0))
section_rho_ros_ann = np.squeeze(np.nanmean(section_rho_ros[:,:,:], axis=0))
#
#print(section_lat_ros.shape, section_z_ros_mask_JFM.shape, section_salt_ros_JFM.shape)

# Prydz Bay
# 1 - slice section for salt/temp/z_rho
section_salt_pry = np.empty((12,31,450)) # allocating
section_temp_pry = np.empty((12,31,450))
section_rho_pry = np.empty((12,31,450))
section_z_pry = np.empty((12,31,450))
section_z_pry_mask = np.empty((12,31,450))
for mm in np.arange(0,12):
    section_salt_pry[mm,:,:] = np.squeeze(salt[mm,:,845,1050:1500]) #ds.salt.isel(xi_rho=slice(1050,1500), eta_rho=845, ocean_time=mm)
    section_temp_pry[mm,:,:] = np.squeeze(temp[mm,:,845,1050:1500]) #ds.temp.isel(xi_rho=slice(1050,1500), eta_rho=845, ocean_time=mm)
    section_z_pry[mm,:,:] = np.squeeze(z_rho[mm,:,845,1050:1500]) #ds.z_rho.isel(xi_rho=slice(1050,1500), eta_rho=845, ocean_time=mm)
    # 1.1 - mask land values in z_rho slice                                           
    section_z_pry_mask[...] = ma.array(section_z_pry[mm,:,:],mask=np.isnan(section_z_pry[mm,:,:]))
    section_rho_pry[mm,:,:] = np.squeeze(dens[mm,:,845,1050:1500])

# 2 - slide section for lon or lat
section_lon_pry_tmp = ds.lon_rho.isel(xi_rho=slice(1050,1500), eta_rho=845)
# 2.1 - mask land values for lon/lat needs loop for repeat through vert layers
section_lon_pry = np.ones((31,450))                                                
for ii in np.arange(0,31):
    section_lon_pry[ii,:] = section_lon_pry_tmp
section_lon_pry_mask = ma.array(section_lon_pry,mask=np.isnan(section_z_pry[0,:,:]))

print(section_lon_pry.shape, section_z_pry.shape, section_salt_pry.shape)

## season averages:
#section_z_pry_mask_JFM = np.squeeze(np.nanmean(section_z_pry_mask[0:18,:,:], axis=0))
#section_z_pry_mask_AMJ = np.squeeze(np.nanmean(section_z_pry_mask[19:36,:,:], axis=0))
#section_z_pry_mask_JAS = np.squeeze(np.nanmean(section_z_pry_mask[37:54,:,:], axis=0))
#section_z_pry_mask_OND = np.squeeze(np.nanmean(section_z_pry_mask[55:72,:,:], axis=0))
#
#section_temp_pry_JFM = np.squeeze(np.nanmean(section_temp_pry[0:18,:,:], axis=0))
#section_temp_pry_AMJ = np.squeeze(np.nanmean(section_temp_pry[19:36,:,:], axis=0))
#section_temp_pry_JAS = np.squeeze(np.nanmean(section_temp_pry[37:54,:,:], axis=0))
#section_temp_pry_OND = np.squeeze(np.nanmean(section_temp_pry[55:72,:,:], axis=0))
#
#section_salt_pry_JFM = np.squeeze(np.nanmean(section_salt_pry[0:18,:,:], axis=0))
#section_salt_pry_AMJ = np.squeeze(np.nanmean(section_salt_pry[19:36,:,:], axis=0))
#section_salt_pry_JAS = np.squeeze(np.nanmean(section_salt_pry[37:54,:,:], axis=0))
#section_salt_pry_OND = np.squeeze(np.nanmean(section_salt_pry[55:72,:,:], axis=0))
#
#section_rho_pry_JFM = np.squeeze(np.nanmean(section_rho_pry[0:18,:,:], axis=0))
#section_rho_pry_AMJ = np.squeeze(np.nanmean(section_rho_pry[19:36,:,:], axis=0))
#section_rho_pry_JAS = np.squeeze(np.nanmean(section_rho_pry[37:54,:,:], axis=0))
#section_rho_pry_OND = np.squeeze(np.nanmean(section_rho_pry[55:72,:,:], axis=0))
#
section_z_pry_mask_ann = np.squeeze(np.nanmean(section_z_pry_mask[:,:,:], axis=0))
section_temp_pry_ann = np.squeeze(np.nanmean(section_temp_pry[:,:,:], axis=0))
section_salt_pry_ann = np.squeeze(np.nanmean(section_salt_pry[:,:,:], axis=0))
section_rho_pry_ann = np.squeeze(np.nanmean(section_rho_pry[:,:,:], axis=0))
#
#print(section_lon_pry.shape, section_z_pry_mask_JFM.shape, section_salt_pry_JFM.shape)

# Maud Land 1:
# 1 - slice section for salt/temp/z_rho
section_salt_neu = np.empty((12,31,250)) # allocating
section_temp_neu = np.empty((12,31,250))
section_rho_neu = np.empty((12,31,250))
section_z_neu = np.empty((12,31,250))
section_z_neu_mask = np.empty((12,31,250))
for mm in np.arange(0,12):
    section_salt_neu[mm,:,:] = np.squeeze(salt[mm,:,1050:1300,700]) # ds.salt.isel(xi_rho=700, eta_rho=slice(1050,1300), ocean_time=mm)
    section_temp_neu[mm,:,:] = np.squeeze(temp[mm,:,1050:1300,700]) #ds.temp.isel(xi_rho=700, eta_rho=slice(1050,1300), ocean_time=mm)
    section_z_neu[mm,:,:] = np.squeeze(z_rho[mm,:,1050:1300,700]) #ds.z_rho.isel(xi_rho=700, eta_rho=slice(1050,1300), ocean_time=mm)
    # 1.1 - mask land values in z_rho slice                                           
    section_z_neu_mask[mm,:,:] = ma.array(section_z_neu[mm,:,:],mask=np.isnan(section_z_neu[mm,:,:]))
    section_rho_neu[mm,:,:] = np.squeeze(dens[mm,:,1050:1300,700])

# 2 - slide section for lon or lat
section_lat_neu_tmp = ds.lat_rho.isel(xi_rho=700, eta_rho=slice(1050,1300))
# 2.1 - mask land values for lon/lat needs loop for repeat through vert layers
section_lat_neu = np.ones((31,250))                                                
for ii in np.arange(0,31):
    section_lat_neu[ii,:] = section_lat_neu_tmp
section_lat_neu_mask = ma.array(section_lat_neu,mask=np.isnan(section_z_neu[0,:,:]))

print(section_lat_neu.shape, section_z_neu.shape, section_salt_neu.shape)

## season averages:
#section_z_neu_mask_JFM = np.squeeze(np.nanmean(section_z_neu_mask[0:18,:,:], axis=0))
#section_z_neu_mask_AMJ = np.squeeze(np.nanmean(section_z_neu_mask[19:36,:,:], axis=0))
#section_z_neu_mask_JAS = np.squeeze(np.nanmean(section_z_neu_mask[37:54,:,:], axis=0))
#section_z_neu_mask_OND = np.squeeze(np.nanmean(section_z_neu_mask[55:72,:,:], axis=0))
#
#section_temp_neu_JFM = np.squeeze(np.nanmean(section_temp_neu[0:18,:,:], axis=0))
#section_temp_neu_AMJ = np.squeeze(np.nanmean(section_temp_neu[19:36,:,:], axis=0))
#section_temp_neu_JAS = np.squeeze(np.nanmean(section_temp_neu[37:54,:,:], axis=0))
#section_temp_neu_OND = np.squeeze(np.nanmean(section_temp_neu[55:72,:,:], axis=0))
#
#section_salt_neu_JFM = np.squeeze(np.nanmean(section_salt_neu[0:18,:,:], axis=0))
#section_salt_neu_AMJ = np.squeeze(np.nanmean(section_salt_neu[19:36,:,:], axis=0))
#section_salt_neu_JAS = np.squeeze(np.nanmean(section_salt_neu[37:54,:,:], axis=0))
#section_salt_neu_OND = np.squeeze(np.nanmean(section_salt_neu[55:72,:,:], axis=0))
#
#section_rho_neu_JFM = np.squeeze(np.nanmean(section_rho_neu[0:18,:,:], axis=0))
#section_rho_neu_AMJ = np.squeeze(np.nanmean(section_rho_neu[19:36,:,:], axis=0))
#section_rho_neu_JAS = np.squeeze(np.nanmean(section_rho_neu[37:54,:,:], axis=0))
#section_rho_neu_OND = np.squeeze(np.nanmean(section_rho_neu[55:72,:,:], axis=0))
#
section_z_neu_mask_ann = np.squeeze(np.nanmean(section_z_neu_mask[:,:,:], axis=0))
section_temp_neu_ann = np.squeeze(np.nanmean(section_temp_neu[:,:,:], axis=0))
section_salt_neu_ann = np.squeeze(np.nanmean(section_salt_neu[:,:,:], axis=0))
section_rho_neu_ann = np.squeeze(np.nanmean(section_rho_neu[:,:,:], axis=0))
#
#print(section_lat_neu.shape, section_z_neu_mask_JFM.shape, section_salt_neu_JFM.shape)

# Maud Land 2:
# 1 - slice section for salt/temp/z_rho
section_salt_mai = np.empty((12,31,250)) # allocating
section_temp_mai = np.empty((12,31,250))
section_rho_mai = np.empty((12,31,250))
section_z_mai = np.empty((12,31,250))
section_z_mai_mask = np.empty((12,31,250))
for mm in np.arange(0,12):
    section_salt_mai[mm,:,:] = np.squeeze(salt[mm,:,1050:1300,875]) #ds.salt.isel(xi_rho=875, eta_rho=slice(1050,1300), ocean_time=mm)
    section_temp_mai[mm,:,:] = np.squeeze(temp[mm,:,1050:1300,875]) #ds.temp.isel(xi_rho=875, eta_rho=slice(1050,1300), ocean_time=mm)
    section_z_mai[mm,:,:] = np.squeeze(z_rho[mm,:,1050:1300,875]) #ds.z_rho.isel(xi_rho=875, eta_rho=slice(1050,1300), ocean_time=mm)
    # 1.1 - mask land values in z_rho slice                                           
    section_z_mai_mask[mm,:,:] = ma.array(section_z_mai[mm,:,:],mask=np.isnan(section_z_mai[mm,:,:]))
    section_rho_mai[mm,:,:] = np.squeeze(dens[mm,:,1050:1300,875])
    
# 2 - slide section for lon or lat
section_lat_mai_tmp = ds.lat_rho.isel(xi_rho=875, eta_rho=slice(1050,1300))
# 2.1 - mask land values for lon/lat needs loop for repeat through vert layers
section_lat_mai = np.ones((31,250))                                                
for ii in np.arange(0,31):
    section_lat_mai[ii,:] = section_lat_mai_tmp
section_lat_mai_mask = ma.array(section_lat_mai,mask=np.isnan(section_z_mai[0,:,:]))

print(section_lat_mai.shape, section_z_mai.shape, section_salt_mai.shape)

## season averages:
#section_z_mai_mask_JFM = np.squeeze(np.nanmean(section_z_mai_mask[0:18,:,:], axis=0))
#section_z_mai_mask_AMJ = np.squeeze(np.nanmean(section_z_mai_mask[19:36,:,:], axis=0))
#section_z_mai_mask_JAS = np.squeeze(np.nanmean(section_z_mai_mask[37:54,:,:], axis=0))
#section_z_mai_mask_OND = np.squeeze(np.nanmean(section_z_mai_mask[55:72,:,:], axis=0))
#
#section_temp_mai_JFM = np.squeeze(np.nanmean(section_temp_mai[0:18,:,:], axis=0))
#section_temp_mai_AMJ = np.squeeze(np.nanmean(section_temp_mai[19:36,:,:], axis=0))
#section_temp_mai_JAS = np.squeeze(np.nanmean(section_temp_mai[37:54,:,:], axis=0))
#section_temp_mai_OND = np.squeeze(np.nanmean(section_temp_mai[55:72,:,:], axis=0))
#
#section_salt_mai_JFM = np.squeeze(np.nanmean(section_salt_mai[0:18,:,:], axis=0))
#section_salt_mai_AMJ = np.squeeze(np.nanmean(section_salt_mai[19:36,:,:], axis=0))
#section_salt_mai_JAS = np.squeeze(np.nanmean(section_salt_mai[37:54,:,:], axis=0))
#section_salt_mai_OND = np.squeeze(np.nanmean(section_salt_mai[55:72,:,:], axis=0))
#
#section_rho_mai_JFM = np.squeeze(np.nanmean(section_rho_mai[0:18,:,:], axis=0))
#section_rho_mai_AMJ = np.squeeze(np.nanmean(section_rho_mai[19:36,:,:], axis=0))
#section_rho_mai_JAS = np.squeeze(np.nanmean(section_rho_mai[37:54,:,:], axis=0))
#section_rho_mai_OND = np.squeeze(np.nanmean(section_rho_mai[55:72,:,:], axis=0))
#
section_z_mai_mask_ann = np.squeeze(np.nanmean(section_z_mai_mask[:,:,:], axis=0))
section_temp_mai_ann = np.squeeze(np.nanmean(section_temp_mai[:,:,:], axis=0))
section_salt_mai_ann = np.squeeze(np.nanmean(section_salt_mai[:,:,:], axis=0))
section_rho_mai_ann = np.squeeze(np.nanmean(section_rho_mai[:,:,:], axis=0))
#
#print(section_lat_mai.shape, section_z_mai_mask_JFM.shape, section_salt_mai_JFM.shape)

# Darnley Polynya
# 1 - slice section for salt/temp/z_rho
section_salt_dar = np.empty((12,31,450)) # allocating
section_temp_dar = np.empty((12,31,450))
section_rho_dar = np.empty((12,31,450))
section_z_dar = np.empty((12,31,450))
section_z_dar_mask = np.empty((12,31,450))
for mm in np.arange(0,12):
    section_salt_dar[mm,:,:] = np.squeeze(salt[mm,:,875,1050:1500]) #ds.salt.isel(xi_rho=slice(1050,1500), eta_rho=875, ocean_time=mm)
    section_temp_dar[mm,:,:] = np.squeeze(temp[mm,:,875,1050:1500]) #ds.temp.isel(xi_rho=slice(1050,1500), eta_rho=875, ocean_time=mm)
    section_z_dar[mm,:,:] = np.squeeze(z_rho[mm,:,875,1050:1500]) #ds.z_rho.isel(xi_rho=slice(1050,1500), eta_rho=875, ocean_time=mm)
    # 1.1 - mask land values in z_rho slice                                           
    section_z_dar_mask[...] = ma.array(section_z_dar[mm,:,:],mask=np.isnan(section_z_dar[mm,:,:]))
    section_rho_dar[mm,:,:] = np.squeeze(dens[mm,:,875,1050:1500])
    
# 2 - slide section for lon or lat
section_lon_dar_tmp = ds.lon_rho.isel(xi_rho=slice(1050,1500), eta_rho=875)
# 2.1 - mask land values for lon/lat needs loop for repeat through vert layers
section_lon_dar = np.ones((31,450))                                                
for ii in np.arange(0,31):
    section_lon_dar[ii,:] = section_lon_dar_tmp
section_lon_dar_mask = ma.array(section_lon_dar,mask=np.isnan(section_z_dar[0,:,:]))

print(section_lon_dar.shape, section_z_dar.shape, section_salt_dar.shape)

## season averages:
#section_z_dar_mask_JFM = np.squeeze(np.nanmean(section_z_dar_mask[0:18,:,:], axis=0))
#section_z_dar_mask_AMJ = np.squeeze(np.nanmean(section_z_dar_mask[19:36,:,:], axis=0))
#section_z_dar_mask_JAS = np.squeeze(np.nanmean(section_z_dar_mask[37:54,:,:], axis=0))
#section_z_dar_mask_OND = np.squeeze(np.nanmean(section_z_dar_mask[55:72,:,:], axis=0))
#
#section_temp_dar_JFM = np.squeeze(np.nanmean(section_temp_dar[0:18,:,:], axis=0))
#section_temp_dar_AMJ = np.squeeze(np.nanmean(section_temp_dar[19:36,:,:], axis=0))
#section_temp_dar_JAS = np.squeeze(np.nanmean(section_temp_dar[37:54,:,:], axis=0))
#section_temp_dar_OND = np.squeeze(np.nanmean(section_temp_dar[55:72,:,:], axis=0))
#
#section_salt_dar_JFM = np.squeeze(np.nanmean(section_salt_dar[0:18,:,:], axis=0))
#section_salt_dar_AMJ = np.squeeze(np.nanmean(section_salt_dar[19:36,:,:], axis=0))
#section_salt_dar_JAS = np.squeeze(np.nanmean(section_salt_dar[37:54,:,:], axis=0))
#section_salt_dar_OND = np.squeeze(np.nanmean(section_salt_dar[55:72,:,:], axis=0))
#
#section_rho_dar_JFM = np.squeeze(np.nanmean(section_rho_dar[0:18,:,:], axis=0))
#section_rho_dar_AMJ = np.squeeze(np.nanmean(section_rho_dar[19:36,:,:], axis=0))
#section_rho_dar_JAS = np.squeeze(np.nanmean(section_rho_dar[37:54,:,:], axis=0))
#section_rho_dar_OND = np.squeeze(np.nanmean(section_rho_dar[55:72,:,:], axis=0))
#
section_z_dar_mask_ann = np.squeeze(np.nanmean(section_z_dar_mask[:,:,:], axis=0))
section_temp_dar_ann = np.squeeze(np.nanmean(section_temp_dar[:,:,:], axis=0))
section_salt_dar_ann = np.squeeze(np.nanmean(section_salt_dar[:,:,:], axis=0))
section_rho_dar_ann = np.squeeze(np.nanmean(section_rho_dar[:,:,:], axis=0))
#
#
#print(section_lon_dar.shape, section_z_dar_mask_JFM.shape, section_salt_dar_JFM.shape)

# Merz Glacier
# 1 - slice section for salt/temp/z_rho
section_salt_mer = np.empty((12,31,175)) # allocating
section_temp_mer = np.empty((12,31,175))
section_rho_mer = np.empty((12,31,175))
section_z_mer = np.empty((12,31,175))
section_z_mer_mask = np.empty((12,31,175))
for mm in np.arange(0,12):
    section_salt_mer[mm,:,:] = np.squeeze(salt[mm,:,25:200,1075]) #ds.salt.isel(xi_rho=1075, eta_rho=slice(25,200), ocean_time=mm)
    section_temp_mer[mm,:,:] = np.squeeze(temp[mm,:,25:200,1075]) #ds.temp.isel(xi_rho=1075, eta_rho=slice(25,200), ocean_time=mm)
    section_z_mer[mm,:,:] = np.squeeze(z_rho[mm,:,25:200,1075]) #ds.z_rho.isel(xi_rho=1075, eta_rho=slice(25,200), ocean_time=mm)
    # 1.1 - mask land values in z_rho slice                                           
    section_z_mer_mask[mm,:,:] = ma.array(section_z_mer[mm,:,:],mask=np.isnan(section_z_mer[mm,:,:]))
    section_rho_mer[mm,:,:] = np.squeeze(dens[mm,:,25:200,1075])
    
# 2 - slide section for lon or lat
section_lat_mer_tmp = ds.lat_rho.isel(xi_rho=1075, eta_rho=slice(25,200))
# 2.1 - mask land values for lon/lat needs loop for repeat through vert layers
section_lat_mer = np.ones((31,175))                                                
for ii in np.arange(0,31):
    section_lat_mer[ii,:] = section_lat_mer_tmp
section_lat_mer_mask = ma.array(section_lat_mer,mask=np.isnan(section_z_mer[0,:,:]))

print(section_lat_mer.shape, section_z_mer.shape, section_salt_mer.shape)

## season averages:
#section_z_mer_mask_JFM = np.squeeze(np.nanmean(section_z_mer_mask[0:18,:,:], axis=0))
#section_z_mer_mask_AMJ = np.squeeze(np.nanmean(section_z_mer_mask[19:36,:,:], axis=0))
#section_z_mer_mask_JAS = np.squeeze(np.nanmean(section_z_mer_mask[37:54,:,:], axis=0))
#section_z_mer_mask_OND = np.squeeze(np.nanmean(section_z_mer_mask[55:72,:,:], axis=0))
#
#section_temp_mer_JFM = np.squeeze(np.nanmean(section_temp_mer[0:18,:,:], axis=0))
#section_temp_mer_AMJ = np.squeeze(np.nanmean(section_temp_mer[19:36,:,:], axis=0))
#section_temp_mer_JAS = np.squeeze(np.nanmean(section_temp_mer[37:54,:,:], axis=0))
#section_temp_mer_OND = np.squeeze(np.nanmean(section_temp_mer[55:72,:,:], axis=0))
#
#section_salt_mer_JFM = np.squeeze(np.nanmean(section_salt_mer[0:18,:,:], axis=0))
#section_salt_mer_AMJ = np.squeeze(np.nanmean(section_salt_mer[19:36,:,:], axis=0))
#section_salt_mer_JAS = np.squeeze(np.nanmean(section_salt_mer[37:54,:,:], axis=0))
#section_salt_mer_OND = np.squeeze(np.nanmean(section_salt_mer[55:72,:,:], axis=0))
#
#section_rho_mer_JFM = np.squeeze(np.nanmean(section_rho_mer[0:18,:,:], axis=0))
#section_rho_mer_AMJ = np.squeeze(np.nanmean(section_rho_mer[19:36,:,:], axis=0))
#section_rho_mer_JAS = np.squeeze(np.nanmean(section_rho_mer[37:54,:,:], axis=0))
#section_rho_mer_OND = np.squeeze(np.nanmean(section_rho_mer[55:72,:,:], axis=0))
#
section_z_mer_mask_ann = np.squeeze(np.nanmean(section_z_mer_mask[:,:,:], axis=0))
section_temp_mer_ann = np.squeeze(np.nanmean(section_temp_mer[:,:,:], axis=0))
section_salt_mer_ann = np.squeeze(np.nanmean(section_salt_mer[:,:,:], axis=0))
section_rho_mer_ann = np.squeeze(np.nanmean(section_rho_mer[:,:,:], axis=0))
#
#print(section_lat_mer.shape, section_z_mer_mask_JFM.shape, section_salt_mer_JFM.shape)

# Amundsen Sea embayment
# 1 - slice section for salt/temp/z_rho
section_salt_amu = np.empty((12,31,375)) # allocating
section_temp_amu = np.empty((12,31,375))
section_rho_amu = np.empty((12,31,375))
section_z_amu = np.empty((12,31,375))
section_z_amu_mask = np.empty((12,31,375))
for mm in np.arange(0,12):
    section_salt_amu[mm,:,:] = np.squeeze(salt[mm,:,512,125:500]) #ds.salt.isel(xi_rho=slice(125,500), eta_rho=512, ocean_time=mm)
    section_temp_amu[mm,:,:] = np.squeeze(temp[mm,:,512,125:500]) #ds.temp.isel(xi_rho=slice(125,500), eta_rho=512, ocean_time=mm)
    section_z_amu[mm,:,:] = np.squeeze(z_rho[mm,:,512,125:500]) #ds.z_rho.isel(xi_rho=slice(125,500), eta_rho=512, ocean_time=mm)
    # 1.1 - mask land values in z_rho slice                                           
    section_z_amu_mask[...] = ma.array(section_z_amu[mm,:,:],mask=np.isnan(section_z_amu[mm,:,:]))
    section_rho_amu[mm,:,:] = np.squeeze(dens[mm,:,512,125:500])
   
# 2 - slide section for lon or lat
section_lon_amu_tmp = ds.lon_rho.isel(xi_rho=slice(125,500), eta_rho=512)
# 2.1 - mask land values for lon/lat needs loop for repeat through vert layers
section_lon_amu = np.ones((31,375))                                                
for ii in np.arange(0,31):
    section_lon_amu[ii,:] = section_lon_amu_tmp
section_lon_amu_mask = ma.array(section_lon_amu,mask=np.isnan(section_z_amu[0,:,:]))

print(section_lon_amu.shape, section_z_amu.shape, section_salt_amu.shape)

## season averages:
#section_z_amu_mask_JFM = np.squeeze(np.nanmean(section_z_amu_mask[0:18,:,:], axis=0))
#section_z_amu_mask_AMJ = np.squeeze(np.nanmean(section_z_amu_mask[19:36,:,:], axis=0))
#section_z_amu_mask_JAS = np.squeeze(np.nanmean(section_z_amu_mask[37:54,:,:], axis=0))
#section_z_amu_mask_OND = np.squeeze(np.nanmean(section_z_amu_mask[55:72,:,:], axis=0))
#
#section_temp_amu_JFM = np.squeeze(np.nanmean(section_temp_amu[0:18,:,:], axis=0))
#section_temp_amu_AMJ = np.squeeze(np.nanmean(section_temp_amu[19:36,:,:], axis=0))
#section_temp_amu_JAS = np.squeeze(np.nanmean(section_temp_amu[37:54,:,:], axis=0))
#section_temp_amu_OND = np.squeeze(np.nanmean(section_temp_amu[55:72,:,:], axis=0))
#
#section_salt_amu_JFM = np.squeeze(np.nanmean(section_salt_amu[0:18,:,:], axis=0))
#section_salt_amu_AMJ = np.squeeze(np.nanmean(section_salt_amu[19:36,:,:], axis=0))
#section_salt_amu_JAS = np.squeeze(np.nanmean(section_salt_amu[37:54,:,:], axis=0))
#section_salt_amu_OND = np.squeeze(np.nanmean(section_salt_amu[55:72,:,:], axis=0))
#
#section_rho_amu_JFM = np.squeeze(np.nanmean(section_rho_amu[0:18,:,:], axis=0))
#section_rho_amu_AMJ = np.squeeze(np.nanmean(section_rho_amu[19:36,:,:], axis=0))
#section_rho_amu_JAS = np.squeeze(np.nanmean(section_rho_amu[37:54,:,:], axis=0))
#section_rho_amu_OND = np.squeeze(np.nanmean(section_rho_amu[55:72,:,:], axis=0))
#
section_z_amu_mask_ann = np.squeeze(np.nanmean(section_z_amu_mask[:,:,:], axis=0))
section_temp_amu_ann = np.squeeze(np.nanmean(section_temp_amu[:,:,:], axis=0))
section_salt_amu_ann = np.squeeze(np.nanmean(section_salt_amu[:,:,:], axis=0))
section_rho_amu_ann = np.squeeze(np.nanmean(section_rho_amu[:,:,:], axis=0))
#
#print(section_lon_amu.shape, section_z_amu_mask_JFM.shape, section_salt_amu_JFM.shape)

# Plot transects: Annual
levelsT = np.arange(-2.5,2.5,.25)
levelsS = np.arange(33.1,35.2,.1)
#levelsR = [24.5,26.6,27.,27.5,27.6,27.7,27.8]#np.arange(26.8,27.9,.1)
levelsR = np.arange(36.4,37.2,.1)

# Amundsen Sea
fig = plt.figure(figsize=(16,18))

ax1 = fig.add_subplot(4,2,1)
ct = plt.contourf(section_lon_amu_mask, section_z_amu_mask_ann, section_temp_amu_ann, levels=levelsT, cmap=plt.cm.coolwarm)
plt.contour(section_lon_amu_mask, section_z_amu_mask_ann, section_temp_amu_ann, colors='k', levels=levelsT)
plt.contour(section_lon_amu_mask, section_z_amu_mask_ann, section_rho_amu_ann, colors='white', levels=[36.9],linewidth=2)
plt.xlim([-115,-106])
plt.ylim([-4200,0])
plt.title('Amundsen Sea',fontsize=12)
plt.ylabel('Depth (m)',fontsize=12)

# Weddell
ax2 = fig.add_subplot(4,2,2)
ct = plt.contourf(section_lat_wed_mask, section_z_wed_mask_ann, section_temp_wed_ann, levels=levelsT, cmap=plt.cm.coolwarm)
plt.contour(section_lat_wed_mask, section_z_wed_mask_ann, section_temp_wed_ann, colors='k', levels=levelsT)
plt.contour(section_lat_wed_mask, section_z_wed_mask_ann, section_rho_wed_ann, colors='white', levels=[36.9],linewidth=2)
plt.xlim([-81,-71])
plt.ylim([-4200,0])
plt.title('Weddell Sea',fontsize=12)

# Maud Land 1
ax3 = fig.add_subplot(4,2,4)
ct = plt.contourf(section_lat_neu_mask, section_z_neu_mask_ann, section_temp_neu_ann, levels=levelsT, cmap=plt.cm.coolwarm)
plt.contour(section_lat_neu_mask, section_z_neu_mask_ann, section_temp_neu_ann, colors='k', levels=levelsT)
plt.contour(section_lat_neu_mask, section_z_neu_mask_ann, section_rho_neu_ann, colors='white', levels=[36.9],linewidth=2)
plt.xlim([-72,-67])
plt.ylim([-4200,0])
plt.title('Maud Land 1',fontsize=12)

# Maud Land 2
ax4 = fig.add_subplot(4,2,6)
ct = plt.contourf(section_lat_mai_mask, section_z_mai_mask_ann, section_temp_mai_ann, levels=levelsT, cmap=plt.cm.coolwarm)
plt.contour(section_lat_mai_mask, section_z_mai_mask_ann, section_temp_mai_ann, colors='k', levels=levelsT)
plt.contour(section_lat_mai_mask, section_z_mai_mask_ann, section_rho_mai_ann, colors='white', levels=[36.9],linewidth=2)
plt.xlim([-72,-67])
plt.ylim([-4200,0])
plt.title('Maud Land 2',fontsize=12)

# Ross
ax5 = fig.add_subplot(4,2,3)
ct = plt.contourf(section_lat_ros_mask, section_z_ros_mask_ann, section_temp_ros_ann, levels=levelsT, cmap=plt.cm.coolwarm)
plt.contour(section_lat_ros_mask, section_z_ros_mask_ann, section_temp_ros_ann, colors='k', levels=levelsT)
plt.contour(section_lat_ros_mask, section_z_ros_mask_ann, section_rho_ros_ann, colors='white', levels=[36.9],linewidth=2)
plt.xlim([-86,-69.5])
plt.ylim([-4200,0])
plt.title('Ross Sea',fontsize=12)
plt.ylabel('Depth (m)',fontsize=12)

# Merz
ax6 = fig.add_subplot(4,2,8)
ct = plt.contourf(section_lat_mer_mask, section_z_mer_mask_ann, section_temp_mer_ann, levels=levelsT, cmap=plt.cm.coolwarm)
plt.contour(section_lat_mer_mask, section_z_mer_mask_ann, section_temp_mer_ann, colors='k', levels=levelsT)
plt.contour(section_lat_mer_mask, section_z_mer_mask_ann, section_rho_mer_ann, colors='white', levels=[36.9],linewidth=2)
plt.xlim([-70,-63.7])
plt.ylim([-4200,0])
plt.title('Merz Glacier',fontsize=12)
#plt.ylabel('Depth (m)')
plt.xlabel('Latitude',fontsize=12)

# prydz bay
ax8 = fig.add_subplot(4,2,7)
ct = plt.contourf(section_lon_pry_mask, section_z_pry_mask_ann, section_temp_pry_ann, levels=levelsT, cmap=plt.cm.coolwarm)
plt.contour(section_lon_pry_mask, section_z_pry_mask_ann, section_temp_pry_ann, colors='k', levels=levelsT)
plt.contour(section_lon_pry_mask, section_z_pry_mask_ann, section_rho_pry_ann, colors='white', levels=[36.9],linewidth=2)
plt.xlim([69,75.5])
plt.ylim([-4200,0])
plt.title('Prydz Bay')
plt.ylabel('Depth (m)',fontsize=12)
plt.xlabel('Longitude',fontsize=12)

# Darnley Polynya
ax7 = fig.add_subplot(4,2,5)
ct = plt.contourf(section_lon_dar_mask, section_z_dar_mask_ann, section_temp_dar_ann, levels=levelsT, cmap=plt.cm.coolwarm)
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
plt.savefig('figs/WAOM4/waom4_all_sections_temp_annual.png')
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
plt.savefig('figs/WAOM4/waom4_all_sections_salt_annual.png')
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
plt.savefig('figs/WAOM4/waom4_all_sections_rho_annual.png', dpi=300)
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

plt.savefig('figs/WAOM4/waom4_DSW_sections_rho_annual.png')
plt.close(fig)
