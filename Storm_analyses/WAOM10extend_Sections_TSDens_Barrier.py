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

section_salt_bar, section_temp_bar, section_z_bar, section_rho_bar, section_z_bar_mask = read_section_roms('bar',2, 2)


dg = xr.open_dataset("/scratch/project_2000339/boeiradi/waom10_frc/waom10extend_grd.nc")

lat_rho = dg.variables["lat_rho"]
lon_rho = dg.variables["lon_rho"]

# Period for storms in the winter of 2007:
events_ind = np.array([37, 38, 42, 43, 44, 45])
date_storms = [('03.Jul.2007'),('08.Jul.2007'),('28.Jul.2007'), \
               ('02.Aug.2007'),('07.Aug.2007'),('12.Aug.2007')]

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

# Plots transects: Monthly
fig_path='/users/boeiradi/COLD_project/postprocessing/figs/Storm_analyses/Barrier_polynya/'

levelsT = np.arange(-2.25,2.26,.25)
levelsS = np.arange(34.,34.8,.1)
levelsR = np.arange(36.2,37.3,.1)

levelsTa = np.arange(-.6,.61,.2)
levelsSa = np.arange(-.3,.31,.1)
levelsRa = np.arange(-.3,.31,.1)

# Plots transects: loop through storm event dates:

for ee in np.arange(0,6,2):
    fig = plt.figure(figsize=(20,12))

    ax1 = fig.add_subplot(2,3,1)
    ct = plt.contourf(section_lat_bar_mask, np.squeeze(section_z_bar_mask[ee,:,:]), np.squeeze(section_temp_bar[ee,:,:]), levels=levelsT, cmap=plt.cm.coolwarm, extend='both')
    plt.contour(section_lat_bar_mask, np.squeeze(section_z_bar_mask[ee,:,:]), np.squeeze(section_temp_bar[ee,:,:]), colors='k', levels=levelsT) #[-2.3,-2.1,-1.9,-1.5,-1.0,-.5,.0])
    plt.xlim([-70,-61])
    plt.ylim([-100,0])
    plt.title(date_storms[ee])
    plt.ylabel('Depth (m)')
    cbarT =fig.colorbar(ct, extend='both')
    cbarT.ax.tick_params(labelsize=14)
    cbarT.set_label('$^{\circ}$C',fontsize=12)

    ax2 = fig.add_subplot(2,3,2)
    cs = plt.contourf(section_lat_bar_mask, np.squeeze(section_z_bar_mask[ee,:,:]), np.squeeze(section_salt_bar[ee,:,:]), levels=levelsS, linestyle='dash', cmap=plt.cm.coolwarm, extend='both')
    plt.contour(section_lat_bar_mask, np.squeeze(section_z_bar_mask[ee,:,:]), np.squeeze(section_salt_bar[ee,:,:]), colors='k', levels=levelsS, linestyle='dashed') #[34.3,34.46,34.6,34.56,34.6])
    plt.xlim([-70,-61])
    plt.ylim([-100,0])
    cbarS =fig.colorbar(cs, extend='both')
    cbarS.ax.tick_params(labelsize=14)
    cbarS.set_label('psu',fontsize=12)
    
    ax3 = fig.add_subplot(2,3,3)
    cr = plt.contourf(section_lat_bar_mask, np.squeeze(section_z_bar_mask[ee,:,:]), np.squeeze(section_rho_bar[ee,:,:]), levels=levelsR, linestyle='dash', cmap=plt.cm.viridis, extend='both')
    plt.contour(section_lat_bar_mask, np.squeeze(section_z_bar_mask[ee,:,:]), np.squeeze(section_rho_bar[ee,:,:]), colors='k', levels=levelsR, linestyle='dashed') #[34.3,34.46,34.6,34.56,34.6])
#    plt.colorbar(cs, extend='both')
    plt.xlim([-70,-61])
    plt.ylim([-100,0])
#    plt.ylabel('Depth (m)')
    cbarR =fig.colorbar(cr, extend='both')
    cbarR.ax.tick_params(labelsize=14)
    cbarR.set_label('kg m$^{-3}$',fontsize=12)

    # ANOMALY:
    ax4 = fig.add_subplot(2,3,4)
    cta = plt.contourf(section_lat_bar_mask, np.squeeze(section_z_bar_mask[ee,:,:]), np.squeeze(section_temp_bar[ee+1,:,:] - section_temp_bar[ee,:,:]), levels=levelsTa, cmap=plt.cm.coolwarm, extend='both')
    plt.contour(section_lat_bar_mask, np.squeeze(section_z_bar_mask[ee,:,:]), np.squeeze(section_temp_bar[ee+1,:,:] - section_temp_bar[ee,:,:]), colors='k', levels=levelsTa)
    plt.xlim([-70,-61])
    plt.ylim([-100,0])
    plt.title(date_storms[ee+1] + ' - ' + date_storms[ee])
    plt.ylabel('Depth (m)')
    plt.xlabel('Latitude')
    cbarTa =fig.colorbar(cta, extend='both')
    cbarTa.ax.tick_params(labelsize=14)
    cbarTa.set_label('$^{\circ}$C',fontsize=12)

    ax5 = fig.add_subplot(2,3,5)
    csa = plt.contourf(section_lat_bar_mask, np.squeeze(section_z_bar_mask[ee,:,:]), np.squeeze(section_salt_bar[ee+1,:,:] - section_salt_bar[ee,:,:]), levels=levelsSa, linestyle='dash', cmap=plt.cm.coolwarm, extend='both')
    plt.contour(section_lat_bar_mask, np.squeeze(section_z_bar_mask[ee,:,:]), np.squeeze(section_salt_bar[ee+1,:,:] -section_salt_bar[ee,:,:]), colors='k', levels=levelsSa, linestyle='dashed')
    plt.xlim([-70,-61])
    plt.ylim([-100,0])
    cbarSa =fig.colorbar(csa, extend='both')
    cbarSa.ax.tick_params(labelsize=14)
    cbarSa.set_label('psu',fontsize=12)
    plt.xlabel('Latitude')

    ax6 = fig.add_subplot(2,3,6)
    cra = plt.contourf(section_lat_bar_mask, np.squeeze(section_z_bar_mask[ee,:,:]), np.squeeze(section_rho_bar[ee+1,:,:] - section_rho_bar[ee,:,:]), levels=levelsRa, linestyle='dash', cmap=plt.cm.coolwarm, extend='both')
    plt.contour(section_lat_bar_mask, np.squeeze(section_z_bar_mask[ee,:,:]), np.squeeze(section_rho_bar[ee+1,:,:] - section_rho_bar[ee,:,:]), colors='k', levels=levelsRa, linestyle='dashed')
#    plt.colorbar(cs, extend='both')
    plt.xlim([-70,-61])
    plt.ylim([-100,0])
#    plt.ylabel('Depth (m)')
    cbarRa =fig.colorbar(cra, extend='both')
    cbarRa.ax.tick_params(labelsize=14)
    cbarRa.set_label('kg m$^{-3}$',fontsize=12)
    plt.xlabel('Latitude')

    name_fig="waom10extend_BarrierPolynya_sections_event=" + date_storms[ee] + ".png"
    plt.savefig(fig_path + name_fig)

