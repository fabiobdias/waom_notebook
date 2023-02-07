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

# allocate sections arrays:
section_salt_wed = np.empty((12,31,375)) # Weddell
section_temp_wed = np.empty((12,31,375))
section_rho_wed = np.empty((12,31,375))
section_z_wed = np.empty((12,31,375))
section_z_wed_mask = np.empty((12,31,375))

section_salt_ros = np.empty((12,31,375)) # Ross
section_temp_ros = np.empty((12,31,375))
section_rho_ros = np.empty((12,31,375))
section_z_ros = np.empty((12,31,375))
section_z_ros_mask = np.empty((12,31,375))

section_salt_pry = np.empty((12,31,450)) # Prydz
section_temp_pry = np.empty((12,31,450))
section_rho_pry = np.empty((12,31,450))
section_z_pry = np.empty((12,31,450))
section_z_pry_mask = np.empty((12,31,450))

section_salt_neu = np.empty((12,31,250)) # Maud 1
section_temp_neu = np.empty((12,31,250))
section_rho_neu = np.empty((12,31,250))
section_z_neu = np.empty((12,31,250))
section_z_neu_mask = np.empty((12,31,250))

section_salt_mai = np.empty((12,31,250)) # Maud 2
section_temp_mai = np.empty((12,31,250))
section_rho_mai = np.empty((12,31,250))
section_z_mai = np.empty((12,31,250))
section_z_mai_mask = np.empty((12,31,250))

section_salt_dar = np.empty((12,31,450)) # Darnley
section_temp_dar = np.empty((12,31,450))
section_rho_dar = np.empty((12,31,450))
section_z_dar = np.empty((12,31,450))
section_z_dar_mask = np.empty((12,31,450))

section_salt_mer = np.empty((12,31,175)) # Merz
section_temp_mer = np.empty((12,31,175))
section_rho_mer = np.empty((12,31,175))
section_z_mer = np.empty((12,31,175))
section_z_mer_mask = np.empty((12,31,175))

section_salt_amu = np.empty((12,31,375)) # Amundsen
section_temp_amu = np.empty((12,31,375))
section_rho_amu = np.empty((12,31,375))
section_z_amu = np.empty((12,31,375))
section_z_amu_mask = np.empty((12,31,375))

# load ROMS avg output
temp = np.empty((12,31,1325,1575))
salt = np.empty((12,31,1325,1575))
z_rho = np.empty((12,31,1325,1575))
count = 0
for mm in ['01','02','03','04','05','06','07','08','09','10','11','12']:
    ds = xr.open_dataset("/scratch/project_2000339/boeiradi/waom4_mahti/output_6yr_fabio_inifile/ocean_avg_00" + mm + ".nc")
# salt
    section_salt_wed_avg = np.divide(np.nansum(ds.salt.isel(xi_rho=425, eta_rho=slice(700,1075)), axis=0), len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
    section_salt_ros_avg = np.divide(np.nansum(ds.salt.isel(xi_rho=775, eta_rho=slice(125,500)), axis=0), len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
    section_salt_pry_avg = np.divide(np.nansum(ds.salt.isel(xi_rho=slice(1050,1500), eta_rho=845), axis=0), len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
    section_salt_neu_avg = np.divide(np.nansum(ds.salt.isel(xi_rho=700, eta_rho=slice(1050,1300)), axis=0), len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
    section_salt_mai_avg = np.divide(np.nansum(ds.salt.isel(xi_rho=875, eta_rho=slice(1050,1300)), axis=0), len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
    section_salt_dar_avg = np.divide(np.nansum(ds.salt.isel(xi_rho=slice(1050,1500), eta_rho=875), axis=0), len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
    section_salt_mer_avg = np.divide(np.nansum(ds.salt.isel(xi_rho=1075, eta_rho=slice(25,200)), axis=0), len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
    section_salt_amu_avg = np.divide(np.nansum(ds.salt.isel(xi_rho=slice(125,500), eta_rho=512), axis=0), len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
# temp
    section_temp_wed_avg = np.divide(np.nansum(ds.temp.isel(xi_rho=425, eta_rho=slice(700,1075)), axis=0), len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
    section_temp_ros_avg = np.divide(np.nansum(ds.temp.isel(xi_rho=775, eta_rho=slice(125,500)), axis=0), len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
    section_temp_pry_avg = np.divide(np.nansum(ds.temp.isel(xi_rho=slice(1050,1500), eta_rho=845), axis=0), len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
    section_temp_neu_avg = np.divide(np.nansum(ds.temp.isel(xi_rho=700, eta_rho=slice(1050,1300)), axis=0), len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
    section_temp_mai_avg = np.divide(np.nansum(ds.temp.isel(xi_rho=875, eta_rho=slice(1050,1300)), axis=0), len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
    section_temp_dar_avg = np.divide(np.nansum(ds.temp.isel(xi_rho=slice(1050,1500), eta_rho=875), axis=0), len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
    section_temp_mer_avg = np.divide(np.nansum(ds.temp.isel(xi_rho=1075, eta_rho=slice(25,200)), axis=0), len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
    section_temp_amu_avg = np.divide(np.nansum(ds.temp.isel(xi_rho=slice(125,500), eta_rho=512), axis=0), len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
# zeta
    section_zeta_wed = ds.zeta.isel(xi_rho=425, eta_rho=slice(700,1075))
    section_zeta_ros = ds.zeta.isel(xi_rho=775, eta_rho=slice(125,500))
    section_zeta_pry = ds.zeta.isel(xi_rho=slice(1050,1500), eta_rho=845)
    section_zeta_neu = ds.zeta.isel(xi_rho=700, eta_rho=slice(1050,1300))
    section_zeta_mai = ds.zeta.isel(xi_rho=875, eta_rho=slice(1050,1300))
    section_zeta_dar = ds.zeta.isel(xi_rho=slice(1050,1500), eta_rho=875)
    section_zeta_mer = ds.zeta.isel(xi_rho=1075, eta_rho=slice(25,200))
    section_zeta_amu = ds.zeta.isel(xi_rho=slice(125,500), eta_rho=512)
# rho/dens
    dens_wed_avg = gsw.rho(section_salt_wed_avg,section_temp_wed_avg,2000) - 1000
    dens_ros_avg = gsw.rho(section_salt_ros_avg,section_temp_ros_avg,2000) - 1000
    dens_pry_avg = gsw.rho(section_salt_pry_avg,section_temp_pry_avg,2000) - 1000
    dens_neu_avg = gsw.rho(section_salt_neu_avg,section_temp_neu_avg,2000) - 1000
    dens_mai_avg = gsw.rho(section_salt_mai_avg,section_temp_mai_avg,2000) - 1000
    dens_dar_avg = gsw.rho(section_salt_dar_avg,section_temp_dar_avg,2000) - 1000
    dens_mer_avg = gsw.rho(section_salt_mer_avg,section_temp_mer_avg,2000) - 1000
    dens_amu_avg = gsw.rho(section_salt_amu_avg,section_temp_amu_avg,2000) - 1000

    if mm == '01':
        ds = ds.set_coords(['Cs_r', 'Cs_w', 'hc', 'Vtransform'])
        h_wed = ds.h.isel(xi_rho=425, eta_rho=slice(700,1075))   
        h_ros = ds.h.isel(xi_rho=775, eta_rho=slice(125,500))
        h_pry = ds.h.isel(xi_rho=slice(1050,1500), eta_rho=845)
        h_neu = ds.h.isel(xi_rho=700, eta_rho=slice(1050,1300))
        h_mai = ds.h.isel(xi_rho=875, eta_rho=slice(1050,1300))
        h_dar = ds.h.isel(xi_rho=slice(1050,1500), eta_rho=875)
        h_mer = ds.h.isel(xi_rho=1075, eta_rho=slice(25,200))
        h_amu = ds.h.isel(xi_rho=slice(125,500), eta_rho=512)
 
    print("Vtransform=2")
#    Zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * ds.h) / (ds.hc + ds.h)
    Zo_rho_wed = (ds.hc * ds.s_rho + ds.Cs_r * h_wed) / (ds.hc + h_wed)
    Zo_rho_ros = (ds.hc * ds.s_rho + ds.Cs_r * h_ros) / (ds.hc + h_ros)
    Zo_rho_pry = (ds.hc * ds.s_rho + ds.Cs_r * h_pry) / (ds.hc + h_pry)
    Zo_rho_neu = (ds.hc * ds.s_rho + ds.Cs_r * h_neu) / (ds.hc + h_neu)
    Zo_rho_mai = (ds.hc * ds.s_rho + ds.Cs_r * h_mai) / (ds.hc + h_mai)
    Zo_rho_dar = (ds.hc * ds.s_rho + ds.Cs_r * h_dar) / (ds.hc + h_dar)
    Zo_rho_mer = (ds.hc * ds.s_rho + ds.Cs_r * h_mer) / (ds.hc + h_mer)
    Zo_rho_amu = (ds.hc * ds.s_rho + ds.Cs_r * h_amu) / (ds.hc + h_amu)

    print('Zo_rho, zeta, zeta_avg, and h sizes: ', Zo_rho_wed.shape, section_zeta_wed.shape, h_wed.shape)
    z_rho_tmp_wed = section_zeta_wed + (section_zeta_wed + h_wed) * Zo_rho_wed
    z_rho_tmp_ros = section_zeta_ros + (section_zeta_ros + h_ros) * Zo_rho_ros
    z_rho_tmp_pry = section_zeta_pry + (section_zeta_pry + h_pry) * Zo_rho_pry
    z_rho_tmp_neu = section_zeta_neu + (section_zeta_neu + h_neu) * Zo_rho_neu
    z_rho_tmp_mai = section_zeta_mai + (section_zeta_mai + h_mai) * Zo_rho_mai
    z_rho_tmp_dar = section_zeta_dar + (section_zeta_dar + h_dar) * Zo_rho_dar
    z_rho_tmp_mer = section_zeta_mer + (section_zeta_mer + h_mer) * Zo_rho_mer
    z_rho_tmp_amu = section_zeta_amu + (section_zeta_amu + h_amu) * Zo_rho_amu
    print('z_rho_tmp_wed size: ', z_rho_tmp_wed.shape)
    del section_zeta_wed, section_zeta_ros, section_zeta_pry, section_zeta_neu, section_zeta_mai, section_zeta_dar, section_zeta_mer, section_zeta_amu
#    z_rho_tmp = z_rho_tmp.transpose('ocean_time', 's_rho', 'eta_rho', 'xi_rho')
    z_rho_tmp_avg_wed = np.divide(np.nansum(z_rho_tmp_wed, axis=0), len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
    z_rho_tmp_avg_ros = np.divide(np.nansum(z_rho_tmp_ros, axis=0), len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
    z_rho_tmp_avg_pry = np.divide(np.nansum(z_rho_tmp_pry, axis=0), len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
    z_rho_tmp_avg_neu = np.divide(np.nansum(z_rho_tmp_neu, axis=0), len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
    z_rho_tmp_avg_mai = np.divide(np.nansum(z_rho_tmp_mai, axis=0), len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
    z_rho_tmp_avg_dar = np.divide(np.nansum(z_rho_tmp_dar, axis=0), len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
    z_rho_tmp_avg_mer = np.divide(np.nansum(z_rho_tmp_mer, axis=0), len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
    z_rho_tmp_avg_amu = np.divide(np.nansum(z_rho_tmp_amu, axis=0), len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
    del z_rho_tmp_wed, z_rho_tmp_ros, z_rho_tmp_pry, z_rho_tmp_neu, z_rho_tmp_mai, z_rho_tmp_dar, z_rho_tmp_mer, z_rho_tmp_amu
    print('temp/salt/z_rho_tmp_avg size = ', section_temp_wed_avg.shape, section_salt_wed_avg.shape, z_rho_tmp_avg_wed.shape)
# assign to monthly sections:
    section_salt_wed[count,:] = section_salt_wed_avg
    section_salt_ros[count,:] = section_salt_ros_avg
    section_salt_pry[count,:] = section_salt_pry_avg
    section_salt_neu[count,:] = section_salt_neu_avg
    section_salt_mai[count,:] = section_salt_mai_avg
    section_salt_dar[count,:] = section_salt_dar_avg
    section_salt_mer[count,:] = section_salt_mer_avg
    section_salt_amu[count,:] = section_salt_amu_avg
    section_temp_wed[count,:] = section_temp_wed_avg
    section_temp_ros[count,:] = section_temp_ros_avg
    section_temp_pry[count,:] = section_temp_pry_avg
    section_temp_neu[count,:] = section_temp_neu_avg
    section_temp_mai[count,:] = section_temp_mai_avg
    section_temp_dar[count,:] = section_temp_dar_avg
    section_temp_mer[count,:] = section_temp_mer_avg
    section_temp_amu[count,:] = section_temp_amu_avg
# (375,31) into shape (31,375)
    z_rho_tmp_avg_wed = z_rho_tmp_avg_wed.transpose(1,0)
    z_rho_tmp_avg_ros = z_rho_tmp_avg_ros.transpose(1,0)
    z_rho_tmp_avg_pry = z_rho_tmp_avg_pry.transpose(1,0)
    z_rho_tmp_avg_neu = z_rho_tmp_avg_neu.transpose(1,0)
    z_rho_tmp_avg_mai = z_rho_tmp_avg_mai.transpose(1,0)
    z_rho_tmp_avg_dar = z_rho_tmp_avg_dar.transpose(1,0)
    z_rho_tmp_avg_mer = z_rho_tmp_avg_mer.transpose(1,0)
    z_rho_tmp_avg_amu = z_rho_tmp_avg_amu.transpose(1,0)
    section_z_wed[count,:] = z_rho_tmp_avg_wed
    section_z_ros[count,:] = z_rho_tmp_avg_ros
    section_z_pry[count,:] = z_rho_tmp_avg_pry
    section_z_neu[count,:] = z_rho_tmp_avg_neu
    section_z_mai[count,:] = z_rho_tmp_avg_mai
    section_z_dar[count,:] = z_rho_tmp_avg_dar
    section_z_mer[count,:] = z_rho_tmp_avg_mer
    section_z_amu[count,:] = z_rho_tmp_avg_amu
    section_rho_wed[count,:] = dens_wed_avg
    section_rho_ros[count,:] = dens_ros_avg
    section_rho_pry[count,:] = dens_pry_avg
    section_rho_neu[count,:] = dens_neu_avg
    section_rho_mai[count,:] = dens_mai_avg
    section_rho_dar[count,:] = dens_dar_avg
    section_rho_mer[count,:] = dens_mer_avg
    section_rho_amu[count,:] = dens_amu_avg
# mask
    section_z_wed_mask[count,:] = ma.array(z_rho_tmp_avg_wed, mask=np.isnan(z_rho_tmp_avg_wed))
    section_z_ros_mask[count,:] = ma.array(z_rho_tmp_avg_ros, mask=np.isnan(z_rho_tmp_avg_ros))
    section_z_pry_mask[count,:] = ma.array(z_rho_tmp_avg_pry, mask=np.isnan(z_rho_tmp_avg_pry))
    section_z_neu_mask[count,:] = ma.array(z_rho_tmp_avg_neu, mask=np.isnan(z_rho_tmp_avg_neu))
    section_z_mai_mask[count,:] = ma.array(z_rho_tmp_avg_mai, mask=np.isnan(z_rho_tmp_avg_mai))
    section_z_dar_mask[count,:] = ma.array(z_rho_tmp_avg_dar, mask=np.isnan(z_rho_tmp_avg_dar))
    section_z_mer_mask[count,:] = ma.array(z_rho_tmp_avg_mer, mask=np.isnan(z_rho_tmp_avg_mer))
    section_z_amu_mask[count,:] = ma.array(z_rho_tmp_avg_amu, mask=np.isnan(z_rho_tmp_avg_amu))
    del section_salt_wed_avg, section_temp_wed_avg, z_rho_tmp_avg_wed, dens_wed_avg
    del section_salt_ros_avg, section_temp_ros_avg, z_rho_tmp_avg_ros, dens_ros_avg
    del section_salt_pry_avg, section_temp_pry_avg, z_rho_tmp_avg_pry, dens_pry_avg
    del section_salt_neu_avg, section_temp_neu_avg, z_rho_tmp_avg_neu, dens_neu_avg
    del section_salt_mai_avg, section_temp_mai_avg, z_rho_tmp_avg_mai, dens_mai_avg
    del section_salt_dar_avg, section_temp_dar_avg, z_rho_tmp_avg_dar, dens_dar_avg
    del section_salt_mer_avg, section_temp_mer_avg, z_rho_tmp_avg_mer, dens_mer_avg
    del section_salt_amu_avg, section_temp_amu_avg, z_rho_tmp_avg_amu, dens_amu_avg
    count = count+1

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
# 2 - slide section for lon or lat
section_lat_wed_tmp = ds.lat_rho.isel(xi_rho=425, eta_rho=slice(700,1075))
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
section_lat_ros_tmp = ds.lat_rho.isel(xi_rho=775, eta_rho=slice(125,500))
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
section_lon_pry_tmp = ds.lon_rho.isel(xi_rho=slice(1050,1500), eta_rho=845)
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
section_lat_neu_tmp = ds.lat_rho.isel(xi_rho=700, eta_rho=slice(1050,1300))
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
section_lat_mai_tmp = ds.lat_rho.isel(xi_rho=875, eta_rho=slice(1050,1300))
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
section_lon_dar_tmp = ds.lon_rho.isel(xi_rho=slice(1050,1500), eta_rho=875)
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
section_lat_mer_tmp = ds.lat_rho.isel(xi_rho=1075, eta_rho=slice(25,200))
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
section_lon_amu_tmp = ds.lon_rho.isel(xi_rho=slice(125,500), eta_rho=512)
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
