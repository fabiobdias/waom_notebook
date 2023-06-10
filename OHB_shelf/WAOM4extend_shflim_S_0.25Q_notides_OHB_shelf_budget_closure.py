# Ocean Heat Budget Analyses in the Antarctica continental shelf (WAOM)

# Fabio B Dias - 3 May 2023
# Description:
#     this script compares methods of calculate the OHC and surface heat flux integrated over the continental shelf
#     - update to xarray open_mfdataset from open_dataset brought new functions more efficient to calculate model layer thickness (dz)
# and ocean heat content tendencies (rho0*Cp*dT/dt). 

# read nc output from WAOM 10km run

import xarray as xr
# import pandas as p
import numpy as np
import numpy.ma as ma
import cartopy.crs as ccrs
import cartopy.feature as cfeature
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

import gsw

import pyresample

from dask.distributed import Client

import warnings
warnings.filterwarnings('ignore')

fig_path = '/g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/figs/Contour_isobath/'

# using xr.open_mfdataset

vars2drop = ["ubar","vbar","w","Hsbl","Hbbl","swrad"]

ds = xr.open_mfdataset(paths="/g/data3/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom4extend_shflim_S_0.25Q/output_yr10_notides_diag/ocean_avg_00*.nc" , chunks={'eta_rho': '200MB'}, parallel=bool, drop_variables=vars2drop, decode_times=False) # , concat_dim="ocean_time"

#- preserving 5-days avgs
temp = ds.variables["temp"]
salt = ds.variables["salt"]
shflux = ds.variables["shflux"]
ssflux = ds.variables["ssflux"]
m = ds.variables["m"]
# HvomT = ds.variables["Hvom_temp"]       ## !!! Huon_temp/Hvom_temp were not saved in the original run
# HuonT = ds.variables["Huon_temp"]       ## now it's running here: /scratch/gi0/fbd581/waom4extend_shflim_S_0.25Q/output_yr10_diag
# Hvom = ds.variables["Hvom"]
# Huon = ds.variables["Huon"]

ds = ds.set_coords(['Cs_r', 'Cs_w', 'hc', 'h', 'Vtransform'])

Zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * ds.h) / (ds.hc + ds.h)
z_rho = ds.zeta + (ds.zeta + ds.h) * Zo_rho + ds.zice
print("Vtransform=2")
Zo_w = (ds.hc * ds.s_w + ds.Cs_w * ds.h) / (ds.hc + ds.h)
z_w = ds.zeta + (ds.zeta + ds.h) * Zo_w + ds.zice

ds.close()

ds = xr.open_mfdataset(paths="/g/data3/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom4extend_shflim_S_0.25Q/output_yr10_notides_diag/ocean_his_00*.nc" , chunks={'eta_rho': '200MB'}, parallel=bool, decode_times=False) #, chuncks="auto", concat_dim="ocean_time"

#- preserving 5-days avgs
temp_snap = ds.variables["temp"]

ds.close()

# calculate dT/dt by differentiating temp_snap:
temp_Rate = np.empty(temp_snap.shape)
dT = np.empty(temp_snap.shape)

# needs the initial conditions:
ds = xr.open_dataset('/g/data3/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom4extend_shflim_S_0.25Q/output_yr9_notides/ocean_rst.nc')
temp_ini = ds.variables["temp"].isel(ocean_time=-1, two=0) # 5-days mean
ds.close()


tlen = len(temp[:,0,0,0])

# transform to DataArray
temp_snap = xr.DataArray(temp_snap)
temp_ini = xr.DataArray(temp_ini)

# - append temp_ini to first time index in temp_snap and then do diff
temp_snap = xr.concat([temp_ini,temp_snap], 'ocean_time')
dT = temp_snap.diff('ocean_time')
print(dT.shape)

#
dt = 5*86400 # 5-days in seconds
temp_Rate = np.divide(dT, dt)
temp_Rate = xr.DataArray(temp_Rate)
# temp_Rate=temp_Rate.rename({'dim_0':'ocean_time','dim_1':'s_rho','dim_2':'eta_rho','dim_3':'xi_rho'})

print(dT.shape)

temp_Rate = temp_Rate.transpose('ocean_time','s_rho','eta_rho','xi_rho')

print(temp_Rate.shape)

# calculate surface sigma_theta (potential density)
sigma_t_sfc = gsw.rho(salt[:,-1,:,:],temp[:,-1,:,:],0) - 1000

# load ice draft to create masks
di = xr.open_dataset('/g/data3/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom4extend_shflim_S_0.25Q/output_yr10_notides_diag/ocean_avg_0001.nc')
ice_draft = di.variables["zice"]

mask_zice = ma.masked_where(ice_draft < 0, np.ones(ice_draft.shape))
mask_outice = ma.masked_where(ice_draft >= 0, np.ones(ice_draft.shape))
di.close()

dg = xr.open_dataset("/g/data3/hh5/tmp/access-om/fbd581/ROMS/waom4_frc/waom4extend_grd.nc")

lat_rho = dg.variables["lat_rho"]
lon_rho = dg.variables["lon_rho"]
lat_u = dg.variables["lat_u"]
lon_u = dg.variables["lon_u"]
lat_v = dg.variables["lat_v"]
lon_v = dg.variables["lon_v"]
pm = dg.variables["pm"]
pn = dg.variables["pn"]
h = dg.variables["h"]

ds.coords['lat_rho']=lat_rho.transpose() # put lat_rho into ds dataset
ds.coords['lon_rho']=lon_rho.transpose() # put lon_rho into ds dataset

area=np.divide(1,pm*pn)

# method 2 to calculate Dz:
# z_w=z_w.chunks(chunks={'eta_rho': '200MB'}) # couldn't change chunks.

Z_w = z_w.transpose('ocean_time','s_w','eta_rho','xi_rho')
print(z_w.shape, Z_w.shape)
dz = np.diff(Z_w,axis=1)

# this mask needs to be replaced by mask_shelf2 that uses the 1500-m contour

mask_shelf = ma.masked_where(h > 1500, np.ones(h.shape))

# determine constants:
rho0 = 1025 # kg. m-3
Cp = 3989.245 # J.kg-1.degC-1
Tf = -1.95 # degC

# use -1000 mask to compute integral of surface heat fluxes and ocean heat content tendency:
# temp_Rate=xr.DataArray(temp_Rate)
temp_rate = temp_Rate.transpose('ocean_time','s_rho','eta_rho','xi_rho')
dT = dT.transpose('ocean_time','s_rho','eta_rho','xi_rho')

## Integrated heat tendency and sfc heat flux terms to check heat budget closure:

# 1. area-integral surface heat flux
tlen = len(temp_rate[:,0,0,0])
area_sum =  np.nansum(np.nansum(area,axis=1), axis=0)

shflux_int = np.empty((tlen))
for mm in np.arange(0,tlen):
    ###shflux_area = shflux[mm,:]*area*mask_shelf*mask_coast
    # needs xarray masking:
    shflux_area = shflux[mm,:].where(h < 1500)*area
    shflux_int[mm] = np.nansum(np.nansum(shflux_area,axis=1), axis=0)
    # del shflux_area
    
# 2. volume-integral heat tendency
temp_rate_int = np.empty((tlen))
temp_rate_vol = np.empty(np.squeeze(temp_Rate[:,0,:,:]).shape)
for mm in np.arange(0,tlen):
# - multplying by dz:
    temp_rate_dz = temp_Rate[mm,:].where(h < 1500)*dz[mm,:]
    temp_rate_vint = np.nansum(temp_rate_dz, axis=0)
    temp_rate_vol[mm,:] = temp_rate_vint*area*mask_shelf#*mask_coast
    del temp_rate_vint
    temp_rate_int[mm] = np.nansum(np.nansum(temp_rate_vol[mm,:],axis=1), axis=0)*Cp*rho0

print(mask_shelf.shape, temp_rate_vol.shape)

# mask of the coast
mask_coast = ma.masked_where(np.isnan(shflux[0,:,:]), np.ones(shflux[0,:,:].shape))

# plot maps dOHC/dt and shflux:
fig, ax = plt.subplots(nrows=1, ncols=2, figsize = (10,5))

ax[0].title.set_text('WAOM4-NOTIDE \n Heat content tendencies')
aa=ax[0].pcolormesh(np.nanmean(temp_rate_vol[1:-1], axis=0)*mask_shelf*mask_coast, vmin=-100, vmax=100, cmap='coolwarm')
# plt.colorbar(aa)

ax[1].title.set_text('WAOM4-NOTIDE \n Surface heat flux')
bb=ax[1].pcolormesh(np.nanmean(shflux[1:-1], axis=0)*mask_shelf*mask_coast, vmin=-100, vmax=100, cmap='coolwarm')

cax2 = plt.axes([0.92, 0.12, 0.01, 0.75])
cb = plt.colorbar(bb, cax=cax2, orientation='vertical')
cb.ax.set_ylabel('W.m$^{-2}$')
name_fig="waom4-NOTIDE_OHC+shflux_shelf_maps_raw.png"
plt.savefig(fig_path + name_fig, dpi=300)
plt.close()

# plot integrated values (5-daily) dOHC/dt and shflux:
ohc_tend = temp_rate_int - np.nanmean(temp_rate_int)
shflux_tend = shflux_int - np.nanmean(shflux_int)

months=np.arange(0,73)*(5/30.41667)

fig = plt.figure(figsize=(10,6))
plt.plot(months,shflux_tend, label='Sfc heat flux')
plt.plot(months,ohc_tend, label='OHC tendency')
plt.plot(months,shflux_int*0,'--k')
plt.plot(months,ohc_tend - shflux_tend, '--r', label='residual',linewidth=0.5)
plt.ylim([-1.5e14,1.5e14])

plt.legend()
plt.ylabel('Heat transport (W)')
plt.xlabel('Time (months)')
plt.title('WAOM4-NOTIDE')

name_fig="waom4-NOTIDE_OHC+shflux_shelf_integ_5daily.png"
plt.savefig(fig_path + name_fig, dpi=300)
plt.close()

## plot maps dOHC/dt and shflux: attempt to use cartopy
proj = ccrs.SouthPolarStereo()

fig = plt.figure(figsize=(12,6))

ax1 = fig.add_subplot(121, projection=proj)
ct1=plt.pcolormesh(lon_rho, lat_rho, np.nanmean(temp_rate_vol[1:-1], axis=0)*mask_shelf*mask_coast, vmin=-100, vmax=100, cmap='coolwarm', transform=ccrs.PlateCarree())
plt.title('WAOM4-NOTIDE \n Heat content tendencies')
# ax1.gridlines() # draw_labels=True,linewidth=
ax1.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')

ax2 = fig.add_subplot(122, projection=proj)
ct1=plt.pcolormesh(lon_rho, lat_rho, np.nanmean(shflux[1:-1], axis=0)*mask_shelf*mask_coast, vmin=-100, vmax=100, cmap='coolwarm', transform=ccrs.PlateCarree())
plt.title('WAOM4-NOTIDE \n Surface heat flux')
# ax2.gridlines() # draw_labels=True,linewidth=
ax2.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')

cax2 = plt.axes([0.92, 0.12, 0.01, 0.75])
cb = plt.colorbar(ct1, cax=cax2, orientation='vertical')
cb.ax.set_ylabel('W.m$^{-2}$')

name_fig="waom4-NOTIDE_OHC+shflux_shelf_maps_annual.png"
plt.savefig(fig_path + name_fig, dpi=300)
plt.close()
