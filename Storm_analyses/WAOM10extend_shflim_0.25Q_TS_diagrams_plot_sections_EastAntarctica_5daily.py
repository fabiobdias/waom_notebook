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
dx = xr.open_dataset("/scratch/project_2000339/boeiradi/waom10_frc/waom10extend_grd.nc")
lat_rho = dx.variables["lat_rho"]
lon_rho = dx.variables["lon_rho"]
#ds.coords['lat_rho']=lat_rho.transpose() # put lat_rho into ds dataset
#ds.coords['lon_rho']=lon_rho.transpose() # put lon_rho into ds dataset

## load section vars:
#bar = '/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM10extend_shflim_S_sec_bar_5days_annual.nc'
#z_rho_bar_tmp = iris.load_cube(bar, 'z_rho')
#temp_bar_tmp = iris.load_cube(bar, 'time-averaged potential temperature') 
#salt_bar_tmp = iris.load_cube(bar, 'time-averaged salinity') 
#depth_bar = z_rho_bar.coord('bathymetry at RHO-points').points
#print(temp_bar_tmp.shape,salt_bar_tmp.shape,depth_bar.shape, temp_bar.shape, z_rho_bar.shape, depth_bar.shape)
#rho_bar = gsw.rho(xr.DataArray.from_iris(salt_bar),xr.DataArray.from_iris(temp_bar),2000) - 1000
#lat_bar_tmp = dx.lat_rho.isel(xi_rho=slice(1250,1550), eta_rho=863)
#z_rho_bar2 = xr.DataArray.from_iris(z_rho_bar)
#z_bar_mask = ma.array(z_rho_bar2,mask=np.isnan(z_rho_bar2))
#lat_bar = np.ones((31,300))
#for ii in np.arange(0,31):
#    lat_bar[ii,:] = lat_bar_tmp
#lat_bar_mask =  ma.array(lat_bar,mask=np.isnan(lat_bar))
#
#vin = '/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM10extend_shflim_S_sec_vin_5days_annual.nc'
#print(vin)
#z_rho_vin_tmp = iris.load_cube(vin, 'z_rho')
#temp_vin_tmp = iris.load_cube(vin, 'time-averaged potential temperature')
#salt_vin_tmp = iris.load_cube(vin, 'time-averaged salinity')
#depth_vin = z_rho_vin.coord('bathymetry at RHO-points').points
#print(temp_vin_tmp.shape,salt_vin_tmp.shape,depth_vin.shape, temp_vin.shape, z_rho_vin.shape, depth_vin.shape)
#rho_vin = gsw.rho(xr.DataArray.from_iris(salt_vin),xr.DataArray.from_iris(temp_vin),2000) - 1000
#lat_vin_tmp = dx.lat_rho.isel(xi_rho=slice(1300,1500), eta_rho=500)
#z_vin_mask = ma.array(xr.DataArray.from_iris(z_rho_vin),mask=np.isnan(xr.DataArray.from_iris(z_rho_vin)))
#lat_vin = np.ones((31,200))
#for ii in np.arange(0,31):
#    lat_vin[ii,:] = lat_vin_tmp
#lat_vin_mask =  ma.array(lat_vin,mask=np.isnan(lat_vin))
#
#pry = '/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM10extend_shflim_S_sec_pry_5days_annual.nc'
##sec = iris.load('/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM10extend_shflim_S_ORAS5em_0.25Q_sec_pry_5days_annual.nc')
#print(pry)
#z_rho_pry_tmp = iris.load_cube(pry, 'z_rho')
#temp_pry_tmp = iris.load_cube(pry, 'time-averaged potential temperature')
#salt_pry_tmp = iris.load_cube(pry, 'time-averaged salinity')
#depth_pry = z_rho_pry.coord('bathymetry at RHO-points').points
#print(temp_pry_tmp.shape,salt_pry_tmp.shape,depth_pry.shape, temp_pry.shape, z_rho_pry.shape, depth_pry.shape)
#rho_pry = gsw.rho(xr.DataArray.from_iris(salt_pry),xr.DataArray.from_iris(temp_pry),2000) - 1000
#lat_pry_tmp = dx.lat_rho.isel(eta_rho=920, xi_rho=slice(1050,1538))
#z_pry_mask = ma.array(xr.DataArray.from_iris(z_rho_pry),mask=np.isnan(xr.DataArray.from_iris(z_rho_pry)))
#lat_pry = np.ones((31,388))
#for ii in np.arange(0,31):
#    lat_pry[ii,:] = lat_pry_tmp
#lat_pry_mask =  ma.array(lat_pry,mask=np.isnan(lat_pry))
#
#sha = '/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM10extend_shflim_S_sec_sha_5days_annual.nc'
##sec = iris.load('/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM10extend_shflim_S_ORAS5em_0.25Q_sec_sha_5days_annual.nc')
#z_rho_sha_tmp = iris.load_cube(sha, 'z_rho')
#temp_sha_tmp = iris.load_cube(sha, 'time-averaged potential temperature')
#salt_sha_tmp = iris.load_cube(sha, 'time-averaged salinity')
#depth_sha = z_rho_sha.coord('bathymetry at RHO-points').points
#print(temp_sha_tmp.shape,salt_sha_tmp.shape,depth_sha.shape, temp_sha.shape, z_rho_sha.shape, depth_sha.shape)
#rho_sha = gsw.rho(xr.DataArray.from_iris(salt_sha),xr.DataArray.from_iris(temp_sha),2000) - 1000
#lat_sha_tmp = dx.lat_rho.isel(xi_rho=slice(1300,1550), eta_rho=675)
#z_sha_mask = ma.array(xr.DataArray.from_iris(z_rho_sha),mask=np.isnan(xr.DataArray.from_iris(z_rho_sha)))
#lat_sha = np.ones((31,250))
#for ii in np.arange(0,31):
#    lat_sha[ii,:] = lat_sha_tmp
#lat_sha_mask =  ma.array(lat_sha,mask=np.isnan(lat_sha))
#
#dar = '/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM10extend_shflim_S_sec_dar_5days_annual.nc'
##sec = iris.load('/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM10extend_shflim_S_ORAS5em_0.25Q_sec_dar_5days_annual.nc')
#print(dar)
#z_rho_dar_tmp = iris.load_cube(dar, 'z_rho')
#temp_dar_tmp = iris.load_cube(dar, 'time-averaged potential temperature')
#salt_dar_tmp = iris.load_cube(dar, 'time-averaged salinity')
#depth_dar = z_rho_dar.coord('bathymetry at RHO-points').points
#print(temp_dar_tmp.shape,salt_dar_tmp.shape,depth_dar.shape, temp_dar.shape, z_rho_dar.shape, depth_dar.shape)
#rho_dar = gsw.rho(xr.DataArray.from_iris(salt_dar),xr.DataArray.from_iris(temp_dar),2000) - 1000
#lat_dar_tmp = dx.lat_rho.isel(eta_rho=950, xi_rho=slice(1150,1538))
#z_dar_mask = ma.array(xr.DataArray.from_iris(z_rho_dar),mask=np.isnan(xr.DataArray.from_iris(z_rho_dar)))
#lat_dar = np.ones((31,388))
#for ii in np.arange(0,31):
#    lat_dar[ii,:] = lat_dar_tmp
#lat_dar_mask =  ma.array(lat_dar,mask=np.isnan(lat_dar))
#
#mer = '/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM10extend_shflim_S_sec_mer_5days_annual.nc'
##sec = iris.load('/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM10extend_shflim_S_ORAS5em_0.25Q_sec_mer_5days_annual.nc')
#print(mer)
#z_rho_mer_tmp = iris.load_cube(mer, 'z_rho')
#temp_mer_tmp = iris.load_cube(mer, 'time-averaged potential temperature')
#salt_mer_tmp = iris.load_cube(mer, 'time-averaged salinity')
#depth_mer = z_rho_mer.coord('bathymetry at RHO-points').points
#print(temp_mer_tmp.shape,salt_mer_tmp.shape,depth_mer.shape, temp_mer.shape, z_rho_mer.shape, depth_mer.shape)
#rho_mer = gsw.rho(xr.DataArray.from_iris(salt_mer),xr.DataArray.from_iris(temp_mer),2000) - 1000
#lat_mer_tmp = dx.lat_rho.isel(xi_rho=1088, eta_rho=slice(100,275))
#z_mer_mask = ma.array(xr.DataArray.from_iris(z_rho_mer),mask=np.isnan(xr.DataArray.from_iris(z_rho_mer)))
#lat_mer = np.ones((31,175))
#for ii in np.arange(0,31):
#    lat_mer[ii,:] = lat_mer_tmp
#lat_mer_mask =  ma.array(lat_mer,mask=np.isnan(lat_mer))

## plots:

fig_path='/users/boeiradi/COLD_project/postprocessing/figs/Storm_analyses/'

# Period for storms in the winter of 2007:
#events_ind = np.array([31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 53, 54]) # annual file
events_ind = np.array([5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 3, 4]) 
date_storms = [('03.Jun.2007'),('08.Jun.2007'),('13.Jun.2007'),('18.Jun.2007'),('23.Jun.2007'),('28.Jun.2007'), \
               ('03.Jul.2007'),('08.Jul.2007'),('13.Jul.2007'),('18.Jul.2007'),('23.Jul.2007'),('28.Jul.2007'), \
               ('02.Aug.2007'),('07.Aug.2007'),('12.Aug.2007'),('17.Aug.2007'),('22.Aug.2007'),('27.Aug.2007'), \
               ('22.Sep.2007'),('27.Sep.2007')]

# loop through events
for ee in np.arange(0,20):

    if ee == 0:
        mm='05'
    elif ee >=1 and ee <= 6:
        mm='06'
    elif ee >= 7 and ee <= 12:
        mm='07'
    elif ee >= 13 and ee <= 17:
        mm='08'
    elif ee >= 18:
        mm='09'

    print(mm, events_ind[ee])

    time_ind = iris.Constraint(averaged_time_since_initialization=events_ind[ee])

    bar = '/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM10extend_shflim_S_sec_bar_5days_mm' + mm + '.nc'
    z_rho_bar = iris.load_cube(bar, 'z_rho')
    print(z_rho_bar[events_ind[ee]])
    temp_bar = iris.load_cube(bar, 'time-averaged potential temperature') 
    salt_bar = iris.load_cube(bar, 'time-averaged salinity')
    print(temp_bar[events_ind[ee]], salt_bar[events_ind[ee]]) 

    sha = '/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM10extend_shflim_S_sec_sha_5days_mm' + mm + '.nc'
    z_rho_sha = iris.load_cube(sha, 'z_rho') 
    temp_sha = iris.load_cube(sha, 'time-averaged potential temperature') 
    salt_sha = iris.load_cube(sha, 'time-averaged salinity') 

    vin = '/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM10extend_shflim_S_sec_vin_5days_mm' + mm + '.nc'
    z_rho_vin = iris.load_cube(vin, 'z_rho') 
    temp_vin = iris.load_cube(vin, 'time-averaged potential temperature') 
    salt_vin = iris.load_cube(vin, 'time-averaged salinity') 

    pry = '/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM10extend_shflim_S_sec_pry_5days_mm' + mm + '.nc'
    z_rho_pry = iris.load_cube(pry, 'z_rho') 
    temp_pry = iris.load_cube(pry, 'time-averaged potential temperature') 
    salt_pry = iris.load_cube(pry, 'time-averaged salinity') 

    dar = '/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM10extend_shflim_S_sec_dar_5days_mm' + mm + '.nc'
    z_rho_dar = iris.load_cube(dar, 'z_rho') 
    temp_dar = iris.load_cube(dar, 'time-averaged potential temperature') 
    salt_dar = iris.load_cube(dar, 'time-averaged salinity') 

    mer = '/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM10extend_shflim_S_sec_mer_5days_mm' + mm + '.nc'
    z_rho_mer = iris.load_cube(mer, 'z_rho') 
    temp_mer = iris.load_cube(mer, 'time-averaged potential temperature') 
    salt_mer = iris.load_cube(mer, 'time-averaged salinity') 

    # plot ts diagram
    fig = plt.figure(figsize=(20, 10))
    
    plt.subplot(2,3,3)
    plt.gcf().subplots_adjust(bottom=0.15)
    for s, t, z in iris.iterate.izip(salt_bar[events_ind[ee]], temp_bar[events_ind[ee]], z_rho_bar[events_ind[ee]], coords='S-coordinate at RHO-points'):
        sr_tmp = z.coord('S-coordinate at RHO-points').points
        dep_tmp = z.coord('bathymetry at RHO-points').points
        scat = iplt.scatter(s,t, c=-sr_tmp*dep_tmp, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=0, vmax=4000)
    ax = plt.gca()
    ax.set_xlim([33.5,35.])
    ax.set_xticks([33.5, 34, 34.5, 35])
    ax.set_ylim([-3,2])
    plt.title('Barrier Polynya', fontsize=14)
    ax.tick_params(labelsize=16)
    CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(33.5,38,.2),linestyles='dashed', colors='k')
    cbw = plt.contour(Si,Ti,dens.transpose(), levels=[36.9],linestyles='dashed', colors='m',linewidth=2)
    plt.clabel(CS, fontsize=12, inline=1, fmt='%1.1f')
    
    plt.subplot(2,3,4)
    plt.gcf().subplots_adjust(bottom=0.15)
    for s, t, z in iris.iterate.izip(salt_sha[events_ind[ee]], temp_sha[events_ind[ee]], z_rho_sha[events_ind[ee]], coords='S-coordinate at RHO-points'):
        sr_tmp = z.coord('S-coordinate at RHO-points').points
        dep_tmp = z.coord('bathymetry at RHO-points').points
        scat = iplt.scatter(s,t, c=-sr_tmp*dep_tmp, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=0, vmax=4000)
    ax = plt.gca()
    ax.set_xlim([33.5,35.])
    ax.set_xticks([33.5, 34, 34.5, 35])
    ax.set_ylim([-3,2])
    plt.title('Shackleton Polynya', fontsize=14)
    ax.set_xlabel('Salinity', fontsize=16)
    ax.set_ylabel('Temperature ($^{\circ}$C)', fontsize=16)
    ax.tick_params(labelsize=16)
    CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(33.5,38,.2),linestyles='dashed', colors='k')
    cbw = plt.contour(Si,Ti,dens.transpose(), levels=[36.9],linestyles='dashed', colors='m',linewidth=2)
    plt.clabel(CS, fontsize=12, inline=1, fmt='%1.1f')
    
    plt.subplot(2,3,5)
    plt.gcf().subplots_adjust(bottom=0.15)
    for s, t, z in iris.iterate.izip(salt_vin[events_ind[ee]], temp_vin[events_ind[ee]], z_rho_vin[events_ind[ee]], coords='S-coordinate at RHO-points'):
        sr_tmp = z.coord('S-coordinate at RHO-points').points
        dep_tmp = z.coord('bathymetry at RHO-points').points
        scat = iplt.scatter(s,t, c=-sr_tmp*dep_tmp, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=0, vmax=4000)
    ax = plt.gca()
    ax.set_xlim([33.5,35.])
    ax.set_xticks([33.5, 34, 34.5, 35])
    ax.set_ylim([-3,2])
    ax.set_xlabel('Salinity', fontsize=16)
    plt.title('Vincennes Polynya', fontsize=14)
    ax.tick_params(labelsize=16)
    CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(33.5,38,.2),linestyles='dashed', colors='k')
    cbw = plt.contour(Si,Ti,dens.transpose(), levels=[36.9],linestyles='dashed', colors='m',linewidth=2)
    plt.clabel(CS, fontsize=12, inline=1, fmt='%1.1f')
    
    plt.subplot(2,3,1)
    plt.gcf().subplots_adjust(bottom=0.15)
    for s, t, z in iris.iterate.izip(salt_dar[events_ind[ee]], temp_dar[events_ind[ee]], z_rho_dar[events_ind[ee]], coords='S-coordinate at RHO-points'):
        sr_tmp = z.coord('S-coordinate at RHO-points').points
        dep_tmp = z.coord('bathymetry at RHO-points').points
        scat = iplt.scatter(s,t, c=-sr_tmp*dep_tmp, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=0, vmax=4000)
    ax = plt.gca()
    #ax.set_xlabel('Salinity', fontsize=16)
    ax.set_ylabel('Temperature ($^{\circ}$C)', fontsize=16)
    ax.set_xlim([33.5,35.])
    ax.set_xticks([33.5, 34, 34.5, 35])
    ax.set_ylim([-3,2])
    plt.title('Darnley Polynya', fontsize=14)
    ax.tick_params(labelsize=16)
    CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(33.5,38,.2),linestyles='dashed', colors='k')
    cbw = plt.contour(Si,Ti,dens.transpose(), levels=[36.9],linestyles='dashed', colors='m',linewidth=2)
    plt.clabel(CS, fontsize=12, inline=1, fmt='%1.1f')
   
    plt.subplot(2,3,2)
    plt.gcf().subplots_adjust(bottom=0.15)
    for s, t, z in iris.iterate.izip(salt_pry[events_ind[ee]], temp_pry[events_ind[ee]], z_rho_pry[events_ind[ee]], coords='S-coordinate at RHO-points'):
        sr_tmp = z.coord('S-coordinate at RHO-points').points
        dep_tmp = z.coord('bathymetry at RHO-points').points
        scat = iplt.scatter(s,t, c=-sr_tmp*dep_tmp, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=0, vmax=4000)
    ax = plt.gca()
    ax.set_xlim([33.5,35.])
    ax.set_xticks([33.5, 34, 34.5, 35])
    ax.set_ylim([-3,2])
    plt.title('Mackenzie Polynya', fontsize=14)
    ax.tick_params(labelsize=16)
    CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(33.5,38,.2),linestyles='dashed', colors='k')
    cbw = plt.contour(Si,Ti,dens.transpose(), levels=[36.9],linestyles='dashed', colors='m',linewidth=2)
    plt.clabel(CS, fontsize=12, inline=1, fmt='%1.1f')
 
    plt.subplot(2,3,6)
    plt.gcf().subplots_adjust(bottom=0.15)
    for s, t, z in iris.iterate.izip(salt_mer[events_ind[ee]], temp_mer[events_ind[ee]], z_rho_mer[events_ind[ee]], coords='S-coordinate at RHO-points'):
        sr_tmp = z.coord('S-coordinate at RHO-points').points
        dep_tmp = z.coord('bathymetry at RHO-points').points
        scat = iplt.scatter(s,t, c=-sr_tmp*dep_tmp, marker='o', edgecolor='none' ,cmap=plt.cm.gist_ncar, vmin=0, vmax=4000)
    ax = plt.gca()
    #ax.set_ylabel('Temperature ($^{\circ}$C)', fontsize=16)
    ax.set_xlabel('Salinity', fontsize=16)
    ax.set_xlim([33.5,35.])
    ax.set_xticks([33.5, 34, 34.5, 35])
    ax.set_ylim([-3,2])
    plt.title('Mertz Polynyar', fontsize=14)
    ax.tick_params(labelsize=16)
    CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(33.5,38,.2),linestyles='dashed', colors='k')
    cbw = plt.contour(Si,Ti,dens.transpose(), levels=[36.9],linestyles='dashed', colors='m',linewidth=2)
    plt.clabel(CS, fontsize=12, inline=1, fmt='%1.1f')
    
    cbar_axim = fig.add_axes([0.92, 0.15, 0.015, 0.7])
    cb = fig.colorbar(scat, cax=cbar_axim)
    cb.ax.tick_params(labelsize=12)
    cb.set_label('Depth (m)',fontsize=14)
    
    name_fig="waom10extend_shflim_S_0.25Q_EastAntarct_sections_TSdiag_yr20_ee=" + date_storms[ee] + ".png"
    plt.savefig(fig_path + name_fig, dpi=300)

