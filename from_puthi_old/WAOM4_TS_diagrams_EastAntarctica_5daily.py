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
xi_pt = [1075]
eta_sec_ini = [25]
eta_sec_end = [200]
# index for xi var sections (dar, pry, bar, sha, vin)
xi_sec_ini = [1150, 1150, 1250, 1300, 1300]
xi_sec_end = [1538, 1538, 1550, 1550, 1500]
eta_pt = [875, 845, 788, 600, 425]
# length of the section
sec_len = [388, 388, 300, 250, 200, 175]

#def read_section_roms(sec_name,ind_ij, ind_len):
#
#    # allocate sections arrays:
#    section_salt = np.empty((73,31,sec_len[ind_len]))
#    section_temp = np.empty((73,31,sec_len[ind_len]))
#    #section_z = np.empty((73,31,sec_len[ind_len]))
#    section_depth = np.empty((73,31,sec_len[ind_len]))
#    # limits for time index
#    monthly_ind_ini = [0, 7, 13, 19, 25, 31, 37, 43, 49, 55, 61, 67]
#    monthly_ind_end = [6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 73]
#    # vars to export
#    variables = ['temp','salt','zeta','z_rho']
#
## load ROMS avg output
#    count = 0
#    for mm in ['01','02','03','04','05','06','07','08','09','10','11','12']:
#        print('Saving month =', mm)
#        ds = xr.open_dataset("/scratch/project_2000339/boeiradi/waom4_mahti/output_6yr_fabio_inifile/ocean_avg_00" + mm + ".nc")
## ----- # first calculate z_rho spatially and put into ds xarray
#        ds = ds.set_coords(['Cs_r', 'Cs_w', 'hc', 'h', 'Vtransform'])
#        Zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * ds.h) / (ds.hc + ds.h)
#        z_rho = ds.zeta + (ds.zeta + ds.h) * Zo_rho
#        print("Vtransform=2, z_rho shape", z_rho.shape)
#        ds.coords['z_rho'] = z_rho.transpose('ocean_time', 's_rho', 'eta_rho', 'xi_rho') 
#
#        print('size temp and time length: ', ds.temp.shape, len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
#        if sec_name == 'mer':
#            print(sec_name, 'xi =', xi_pt[ind_ij], 'eta =', eta_sec_ini[ind_ij], eta_sec_end[ind_ij])
##            section_salt_avg = ds.salt.isel(xi_rho=xi_pt[ind_ij], eta_rho=slice(eta_sec_ini[ind_ij],eta_sec_end[ind_ij]))
##            section_temp_avg = ds.temp.isel(xi_rho=xi_pt[ind_ij], eta_rho=slice(eta_sec_ini[ind_ij],eta_sec_end[ind_ij]))
#            # saving all variables including temp/salt
#            ds[variables].isel(ocean_time=slice(0, 8),xi_rho=xi_pt[ind_ij], eta_rho=slice(eta_sec_ini[ind_ij],eta_sec_end[ind_ij])).to_netcdf('/scratch/project_2000339/boeiradi/postprocessing/ncdf_tmp/WAOM4_sec_' + sec_name + '_5days_mm' + mm + '.nc', mode='w')
#        elif sec_name == 'dar' or sec_name ==  'pry' or sec_name ==  'bar' or sec_name ==  'sha' or sec_name ==  'vin':
#            print(sec_name, 'eta =', eta_pt[ind_ij], 'xi =', xi_sec_ini[ind_ij], xi_sec_end[ind_ij])
#            section_salt_avg = ds.salt.isel(xi_rho=slice(xi_sec_ini[ind_ij],xi_sec_end[ind_ij]), eta_rho=eta_pt[ind_ij])
#            section_temp_avg = ds.temp.isel(xi_rho=slice(xi_sec_ini[ind_ij],xi_sec_end[ind_ij]), eta_rho=eta_pt[ind_ij])
#            # saving all variables including temp/salt
#            ds[variables].isel(ocean_time=slice(0, 8),xi_rho=slice(xi_sec_ini[ind_ij],xi_sec_end[ind_ij]), eta_rho=eta_pt[ind_ij]).to_netcdf('/scratch/project_2000339/boeiradi/postprocessing/ncdf_tmp/WAOM4_sec_' + sec_name + '_5days_mm' + mm + '.nc', mode='w')
#
##        section_salt[monthly_ind_ini[count]:monthly_ind_end[count]+1,:] = section_salt_avg
##        section_temp[monthly_ind_ini[count]:monthly_ind_end[count]+1,:] = section_temp_avg
##        del section_salt_avg, section_temp_avg
##        count = count+1
#
#    return #section_salt, section_temp
#
##section_salt_dar, section_temp_dar = read_section_roms('dar',0, 0)
##section_salt_pry, section_temp_pry = read_section_roms('pry',1, 1)
##section_salt_bar, section_temp_bar = read_section_roms('bar',2, 2)
##section_salt_sha, section_temp_sha = read_section_roms('sha',3, 3)
##section_salt_vin, section_temp_vin = read_section_roms('vin',4, 4)
##section_salt_mer, section_temp_mer = read_section_roms('mer',0, 5)
#read_section_roms('dar',0, 0)
#read_section_roms('pry',1, 1)
#read_section_roms('bar',2, 2)
#read_section_roms('sha',3, 3)
#read_section_roms('vin',4, 4)
#read_section_roms('mer',0, 5)

# Period for storms in the winter of 2007:
events_ind = np.array([31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 53, 54])
date_storms = [('03.Jun.2007'),('08.Jun.2007'),('13.Jun.2007'),('18.Jun.2007'),('23.Jun.2007'),('28.Jun.2007'), \
               ('03.Jul.2007'),('08.Jul.2007'),('13.Jul.2007'),('18.Jul.2007'),('23.Jul.2007'),('28.Jul.2007'), \
               ('02.Aug.2007'),('07.Aug.2007'),('12.Aug.2007'),('17.Aug.2007'),('22.Aug.2007'), \
               ('21.Sep.2007'),('26.Sep.2007')]

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

section_salt_dar = np.empty((73,31,sec_len[0]))
section_temp_dar = np.empty((73,31,sec_len[0]))
z_rho_dar = np.empty((73,31,sec_len[0]))
section_salt_pry = np.empty((73,31,sec_len[1]))
section_temp_pry = np.empty((73,31,sec_len[1]))
z_rho_pry = np.empty((73,31,sec_len[1]))
section_salt_bar = np.empty((73,31,sec_len[2]))
section_temp_bar = np.empty((73,31,sec_len[2]))
z_rho_bar = np.empty((73,31,sec_len[2]))
section_salt_sha = np.empty((73,31,sec_len[3]))
section_temp_sha = np.empty((73,31,sec_len[3]))
z_rho_sha = np.empty((73,31,sec_len[3]))
section_salt_vin = np.empty((73,31,sec_len[4]))
section_temp_vin = np.empty((73,31,sec_len[4]))
z_rho_vin = np.empty((73,31,sec_len[4]))
section_salt_mer = np.empty((73,31,sec_len[5]))
section_temp_mer = np.empty((73,31,sec_len[5]))
z_rho_mer = np.empty((73,31,sec_len[5]))

for mm in ['01','02','03','04','05','06','07','08','09','10','11','12']:
    print('Reading month =', mm)
    dar=iris.load('/scratch/project_2000339/boeiradi/postprocessing/ncdf_tmp/WAOM4_sec_dar_5days_mm' + mm + '.nc')
    print(dar)
    tmp_z_rho_dar = dar[1]
    tmp_section_salt_dar = dar[0]
    tmp_section_temp_dar = dar[3]

    pry=iris.load('/scratch/project_2000339/boeiradi/postprocessing/ncdf_tmp/WAOM4_sec_pry_5days_mm' + mm + '.nc')
    print(pry)
    tmp_z_rho_pry = pry[1]
    tmp_section_salt_pry = pry[0]
    tmp_section_temp_pry = pry[3]

    bar=iris.load('/scratch/project_2000339/boeiradi/postprocessing/ncdf_tmp/WAOM4_sec_bar_5days_mm' + mm + '.nc')
    print(bar)
    tmp_z_rho_bar = bar[1]
    tmp_section_salt_bar = bar[0]
    tmp_section_temp_bar = bar[3]

    sha=iris.load('/scratch/project_2000339/boeiradi/postprocessing/ncdf_tmp/WAOM4_sec_sha_5days_mm' + mm + '.nc')
    print(sha)
    tmp_z_rho_sha = sha[1]
    tmp_section_salt_sha= sha[0]
    tmp_section_temp_sha = sha[3]
   
    vin=iris.load('/scratch/project_2000339/boeiradi/postprocessing/ncdf_tmp/WAOM4_sec_vin_5days_mm' + mm + '.nc')
    print(vin)
    tmp_z_rho_vin = vin[0]
    tmp_section_salt_vin = vin[2]
    tmp_section_temp_vin = vin[3]
  
    mer=iris.load('/scratch/project_2000339/boeiradi/postprocessing/ncdf_tmp/WAOM4_sec_mer_5days_mm' + mm + '.nc')
    print(mer)
    tmp_z_rho_mer = mer[0]
    tmp_section_salt_mer = mer[2]
    tmp_section_temp_mer = mer[3]

    z_rho_dar[monthly_ind_ini[count]:monthly_ind_end[count]+1,:] = tmp_z_rho_dar
    z_rho_pry[monthly_ind_ini[count]:monthly_ind_end[count]+1,:] = tmp_z_rho_pry
    z_rho_bar[monthly_ind_ini[count]:monthly_ind_end[count]+1,:] = tmp_z_rho_bar
    z_rho_sha[monthly_ind_ini[count]:monthly_ind_end[count]+1,:] = tmp_z_rho_sha
    z_rho_vin[monthly_ind_ini[count]:monthly_ind_end[count]+1,:] = tmp_z_rho_vin
    z_rho_mer[monthly_ind_ini[count]:monthly_ind_end[count]+1,:] = tmp_z_rho_mer
    section_salt_dar[monthly_ind_ini[count]:monthly_ind_end[count]+1,:] = tmp_section_salt_dar
    section_salt_pry[monthly_ind_ini[count]:monthly_ind_end[count]+1,:] = tmp_section_salt_pry
    section_salt_bar[monthly_ind_ini[count]:monthly_ind_end[count]+1,:] = tmp_section_salt_bar
    section_salt_sha[monthly_ind_ini[count]:monthly_ind_end[count]+1,:] = tmp_section_salt_sha
    section_salt_vin[monthly_ind_ini[count]:monthly_ind_end[count]+1,:] = tmp_section_salt_vin
    section_salt_mer[monthly_ind_ini[count]:monthly_ind_end[count]+1,:] = tmp_section_salt_mer
    section_temp_dar[monthly_ind_ini[count]:monthly_ind_end[count]+1,:] = tmp_section_temp_dar
    section_temp_pry[monthly_ind_ini[count]:monthly_ind_end[count]+1,:] = tmp_section_temp_pry
    section_temp_bar[monthly_ind_ini[count]:monthly_ind_end[count]+1,:] = tmp_section_temp_bar
    section_temp_sha[monthly_ind_ini[count]:monthly_ind_end[count]+1,:] = tmp_section_temp_sha
    section_temp_vin[monthly_ind_ini[count]:monthly_ind_end[count]+1,:] = tmp_section_temp_vin
    section_temp_mer[monthly_ind_ini[count]:monthly_ind_end[count]+1,:] = tmp_section_temp_mer

    del tmp_z_rho_dar, tmp_z_rho_pry, tmp_z_rho_bar, tmp_z_rho_sha, tmp_z_rho_vin, tmp_z_rho_mer 
    del tmp_section_salt_dar, tmp_section_salt_pry, tmp_section_salt_bar, tmp_section_salt_sha, tmp_section_salt_vin, tmp_section_salt_mer
    del tmp_section_temp_dar, tmp_section_temp_pry, tmp_section_temp_bar, tmp_section_temp_shatmp_section_temp_vin, tmp_section_temp_mer
   
    count = count + 1

depth_dar = z_rho_dar.coord('bathymetry at RHO-points').points
depth_pry = z_rho_pry.coord('bathymetry at RHO-points').points
depth_bar = z_rho_bar.coord('bathymetry at RHO-points').points
depth_sha = z_rho_sha.coord('bathymetry at RHO-points').points
depth_vin = z_rho_vin.coord('bathymetry at RHO-points').points
depth_mer = z_rho_mer.coord('bathymetry at RHO-points').points

fig_path='/scratch/project_2000339/boeiradi/postprocessing/figs/WAOM4/'

# TS plots for 5-day averages:
for ee in np.arange(0,19):
    print('Timestep = ', events_ind[ee])
    fig = plt.figure(figsize=(20, 10))

    plt.subplot(2,4,1)
    plt.gcf().subplots_adjust(bottom=0.15)
    #print(section_salt_dar.shape, section_temp_dar.shape)#,section_z_dar.shape)
    for s, t in iris.iterate.izip(section_salt_dar[events_ind[ee],:,:], section_temp_dar[events_ind[ee],:,:], coords='bathymetry at RHO-points'):
        scat = iplt.scatter(s,t, c=depth_dar, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=200, vmax=4000)
    ax = plt.gca()
    ax.set_ylabel('Temperature ($^{\circ}$C)', fontsize=16)
    #ax.set_xlabel('Salinity', fontsize=16)
    ax.set_xlim([33.,35.])
    ax.set_xticks([33, 33.5, 34, 34.5, 35])
    ax.set_ylim([-3,5])
    plt.title('Darnley Polynya', fontsize=14)
    ax.tick_params(labelsize=16)
    CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(33.5,38,.2),linestyles='dashed', colors='k')
    cbw = plt.contour(Si,Ti,dens.transpose(), levels=[36.9],linestyles='dashed', colors='m',linewidth=2)
    plt.clabel(CS, fontsize=12, inline=1, fmt='%1.1f')

    plt.subplot(2,4,2)
    plt.gcf().subplots_adjust(bottom=0.15)
    #print(section_salt_pry.shape, section_temp_pry.shape)#,section_z_pry.shape)
    for s, t in iris.iterate.izip(section_salt_pry[events_ind[ee],:,:], section_temp_pry[events_ind[ee],:,:], coords='bathymetry at RHO-points'):
        scat = iplt.scatter(s,t, c=depth_pry, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=200, vmax=4000)
    ax = plt.gca()
    ax.set_xlim([33.,35.])
    ax.set_xticks([33, 33.5, 34, 34.5, 35])
    ax.set_ylim([-3,5])
    plt.title('Mackenzi Polynya', fontsize=14)
    ax.tick_params(labelsize=16)
    CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(33.5,38,.2),linestyles='dashed', colors='k')
    cbw = plt.contour(Si,Ti,dens.transpose(), levels=[36.9],linestyles='dashed', colors='m',linewidth=2)
    plt.clabel(CS, fontsize=12, inline=1, fmt='%1.1f')

    plt.subplot(2,4,3)
    plt.gcf().subplots_adjust(bottom=0.15)
    for s, t in iris.iterate.izip(section_salt_bar[events_ind[ee],:,:], section_temp_bar[events_ind[ee],:,:], coords='bathymetry at RHO-points'):
        scat = iplt.scatter(s,t, c=depth_bar, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=200, vmax=4000)
    ax = plt.gca()
    ax.set_xlim([33.,35.])
    ax.set_xticks([33, 33.5, 34, 34.5, 35])
    ax.set_ylim([-3,5])
    plt.title('Barrier Polynya', fontsize=14)
    ax.tick_params(labelsize=16)
    CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(33.5,38,.2),linestyles='dashed', colors='k')
    cbw = plt.contour(Si,Ti,dens.transpose(), levels=[36.9],linestyles='dashed', colors='m',linewidth=2)
    plt.clabel(CS, fontsize=12, inline=1, fmt='%1.1f')

    plt.subplot(2,4,4)
    plt.gcf().subplots_adjust(bottom=0.15)
    for s, t in iris.iterate.izip(section_salt_sha[events_ind[ee],:,:], section_temp_sha[events_ind[ee],:,:], coords='bathymetry at RHO-points'):
        scat = iplt.scatter(s,t, c=depth_sha, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=200, vmax=4000)
    ax = plt.gca()
    ax.set_xlabel('Salinity', fontsize=16)
    ax.set_xlim([33.,35.])
    ax.set_xticks([33, 33.5, 34, 34.5, 35])
    ax.set_ylim([-3,5])
    plt.title('Shackleton Polynya', fontsize=14)
    ax.tick_params(labelsize=16)
    CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(33.5,38,.2),linestyles='dashed', colors='k')
    cbw = plt.contour(Si,Ti,dens.transpose(), levels=[36.9],linestyles='dashed', colors='m',linewidth=2)
    plt.clabel(CS, fontsize=12, inline=1, fmt='%1.1f')

    plt.subplot(2,4,5)
    plt.gcf().subplots_adjust(bottom=0.15)
    for s, t in iris.iterate.izip(section_salt_vin[events_ind[ee],:,:], section_temp_vin[events_ind[ee],:,:], coords='bathymetry at RHO-points'):
        scat = iplt.scatter(s,t, c=depth_vin, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=200, vmax=4000)
    ax = plt.gca()
    ax.set_ylabel('Temperature ($^{\circ}$C)', fontsize=16)
    ax.set_xlabel('Salinity', fontsize=16)
    ax.set_xlim([33.,35.])
    ax.set_xticks([33, 33.5, 34, 34.5, 35])
    ax.set_ylim([-3,5])
    plt.title('Vincennes Polynya', fontsize=14)
    ax.tick_params(labelsize=16)
    CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(33.5,38,.2),linestyles='dashed', colors='k')
    cbw = plt.contour(Si,Ti,dens.transpose(), levels=[36.9],linestyles='dashed', colors='m',linewidth=2)
    plt.clabel(CS, fontsize=12, inline=1, fmt='%1.1f')

    plt.subplot(2,4,6)
    plt.gcf().subplots_adjust(bottom=0.15)
    print(section_salt_vin.shape, section_temp_vin.shape)
    for s, t in iris.iterate.izip(section_salt_mer[events_ind[ee],:,:], section_temp_mer[events_ind[ee],:,:], coords='bathymetry at RHO-points'):
        scat = iplt.scatter(s,t, c=depth_mer, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=200, vmax=4000)
    ax = plt.gca()
    ax.set_xlabel('Salinity', fontsize=16)
    ax.set_xlim([33.,35.])
    ax.set_xticks([33, 33.5, 34, 34.5, 35])
    ax.set_ylim([-3,5])
    plt.title('Merz Polynya', fontsize=14)
    ax.tick_params(labelsize=16)
    CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(33.5,38,.2),linestyles='dashed', colors='k')
    cbw = plt.contour(Si,Ti,dens.transpose(), levels=[36.9],linestyles='dashed', colors='m',linewidth=2)
    plt.clabel(CS, fontsize=12, inline=1, fmt='%1.1f')

    cbar_axim = fig.add_axes([0.92, 0.15, 0.015, 0.7])
    cb = fig.colorbar(scat, cax=cbar_axim)
    cb.ax.tick_params(labelsize=12)
    cb.set_label('Depth (m)',fontsize=14)

    name_fig="waom4_EastAntarct_sections_TSdiag_" + date_storms[ee] + ".png"
    plt.savefig(fig_path + name_fig, dpi=300)
    plt.close()
