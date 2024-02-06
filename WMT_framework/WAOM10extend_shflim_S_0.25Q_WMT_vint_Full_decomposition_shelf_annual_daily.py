# read nc output from WAOM 10km run

import xarray as xr
import pandas as p
import numpy as np
import numpy.ma as ma
import cartopy.crs as ccrs
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

from dask.distributed import Client
import logging
import warnings
warnings.filterwarnings('ignore')

client = Client(threads_per_worker=1, memory_limit=0, silence_logs=logging.ERROR)
print(client)

ds = xr.open_mfdataset(paths="/g/data/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom10extend_shflim_S_0.25Q/output_20yr_diag_daily/ocean_avg_00*.nc" , chunks={'eta_rho': '200MB'}, parallel=True, decode_times=False) # , concat_dim="ocean_time"

#- preserving 5-days avgs
temp = ds.variables["temp"]
salt = ds.variables["salt"]
shflux = ds.variables["shflux"]
ssflux = ds.variables["ssflux"]
m = ds.variables["m"]

ds = ds.set_coords(['Cs_r', 'Cs_w', 'hc', 'h', 'Vtransform'])

Zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * ds.h) / (ds.hc + ds.h)
z_rho = ds.zeta + (ds.zeta + ds.h) * Zo_rho + ds.zice
print("Vtransform=2")
Zo_w = (ds.hc * ds.s_w + ds.Cs_w * ds.h) / (ds.hc + ds.h)
z_w = ds.zeta + (ds.zeta + ds.h) * Zo_w + ds.zice

ds.close()

sigma_t_sfc = gsw.rho(salt.isel(s_rho=-1),temp.isel(s_rho=-1),0) - 1000

di = xr.open_mfdataset(paths='/g/data/hh5/tmp/access-om/fbd581/ROMS/output/waom10extend_shflim_S_0.25Q/output_20yr_diag_daily/ocean_avg_0001.nc')
ice_draft = di.variables["zice"]

mask_zice = ma.masked_where(ice_draft < 0, np.ones(ice_draft.shape))

di.close()

# calculate dz following:
dz = np.empty((365,560,630,31))
dz_inv = np.empty((365,560,630,31))


for tt in np.arange(0,365):
    z_w_sorted = -1*z_w.i[tt,:,:,::-1]
    print(z_w_sorted.shape)
    dz_inv[tt,:,:,:] = np.diff(z_w_sorted,axis=2)
    dz[tt,:,:,:] = dz_inv[tt,:,:,::-1]

dg = xr.open_dataset("/g/data/hh5/tmp/access-om/fbd581/ROMS/waom10_frc/waom10extend_grd.nc")

lat_rho = dg.variables["lat_rho"]
lon_rho = dg.variables["lon_rho"]
pm = dg.variables["pm"]
pn = dg.variables["pn"]

ds.coords['lat_rho']=lat_rho.transpose() # put lat_rho into ds dataset
ds.coords['lon_rho']=lon_rho.transpose() # put lon_rho into ds dataset

# shelf/open-ocean masks:
mask_open = ma.masked_where(dg.h <= 2000, np.ones(dg.h.shape))
mask_shelf = ma.masked_where(dg.h > 2000, np.ones(dg.h.shape))

fig = plt.figure(figsize=(10,12))
ax1 = fig.add_subplot(221)#, projection=proj)
cy=plt.pcolor(mask_shelf)#, transform=ccrs.PlateCarree())
plt.colorbar(cy)
plt.clim(0.,1.)
ax2 = fig.add_subplot(222)#, projection=proj)
cy=plt.pcolor(mask_open)#, transform=ccrs.PlateCarree())
plt.colorbar(cy)
plt.clim(0.,1.)
ax3 = fig.add_subplot(223)#, projection=proj)
cy=plt.pcolor(mask_zice)#, transform=ccrs.PlateCarree())
plt.colorbar(cy)
plt.clim(0.,1.)

dx = xr.open_dataset('/g/data/hh5/tmp/access-om/fbd581/ROMS/output/waom10extend_shflim_S_0.25Q/output_20yr_diag_daily/Full_vint_vars_for_WMT_m.s-1.nc')

# - variables integrated throughout the ML; multiply by -1 b/c dz is negative.
temp_vdia_diff_full_vint = dx.variables["temp_vdia_diff_full_vint"]
salt_vdia_diff_full_vint = dx.variables["salt_vdia_diff_full_vint"]
temp_hdia_diff_full_vint = dx.variables["temp_hdia_diff_full_vint"]
salt_hdia_diff_full_vint = dx.variables["salt_hdia_diff_full_vint"]
temp_vdia_adv_full_vint = dx.variables["temp_vdia_adv_full_vint"]
salt_vdia_adv_full_vint = dx.variables["salt_vdia_adv_full_vint"]
temp_hdia_adv_full_vint = dx.variables["temp_hdia_adv_full_vint"]
salt_hdia_adv_full_vint = dx.variables["salt_hdia_adv_full_vint"]
temp_tend_full_vint = dx.variables["temp_tend_avg_full_vint"]
salt_tend_full_vint = dx.variables["salt_tend_avg_full_vint"]

sigma_t = gsw.rho(salt[:,-1,:,:],temp[:,-1,:,:],0) - 1000

dx.close()

# obtain thermal expansion (alpha) & salinity contraction (beta) coefficients:
SA = np.empty(salt.shape)
# neet Absolute Salinity, converting from Pratical Salinity:
print('salt and z_rho shape:', np.squeeze(salt[0,0,:,:]).shape,np.squeeze(z_rho[0,:,:,0].shape))
for mm in np.arange(0,12):
    for kk in np.arange(0,31):
        SA_tmp =gsw.SA_from_SP(np.squeeze(salt[mm,kk,:,:]),np.squeeze(z_rho[mm,:,:,kk]),lon_rho,lat_rho)
        SA[mm,kk,:,:] = SA_tmp
        del SA_tmp

# gsw.alpha/gsw.beta
#alpha = gsw_alpha(SA,CT,p)
[specvol, alpha, beta] = gsw.specvol_alpha_beta(SA,temp,z_rho.transpose(0,3,1,2))

print('Alpha/beta shapes:', alpha.shape, beta.shape)

# calculate the LHS term in Pellichero et al (2018):
# ps: Diffusion (R_s, R_t) terms already include the sfc fluxes

# heat (eqn 5)
rho0 = 1025 #1000
Cp = 3985

# total diffusion terms:
R_s_vint = (salt_hdia_diff_full_vint + salt_vdia_diff_full_vint)
R_t_vint = (temp_hdia_diff_full_vint + temp_vdia_diff_full_vint)

# surface flux terms:
salt_sfc = beta[:,-1,:,:]*(ssflux)
temp_sfc = alpha[:,-1,:,:]*(np.divide(shflux, rho0*Cp))

# advection terms:
salt_adv_full_vint = (salt_hdia_adv_full_vint + salt_vdia_adv_full_vint)
temp_adv_full_vint = (temp_hdia_adv_full_vint + temp_vdia_adv_full_vint)

# net tendencies
salt_net_full_vint = salt_tend_full_vint
temp_net_full_vint = temp_tend_full_vint

#  Function to calculate Water Mass Transformation (in m3/s):

# rho grid for binning:
#rho_grid=np.arange(35.5,37.4,0.1) # for sigma-2
#rho_grid=np.arange(24.4,29.1,0.1) # for sigma-0
rho_grid=np.arange(26.,28.,0.05) # for sigma-0

len_rho_grid=len(rho_grid)

dx = np.divide(1,pm)
dy = np.divide(1,pn)
dt = 86400#30#/12 #why divide by 12?

def wmt(var_int, dx, dy,var_type):
    # var_type: 'budget' or 'sfc_frc'
    
    F_rate_var_vint = np.empty(var_int.shape)
    
    for mm in np.arange(0,365):
        if var_type == 'budget':
            F_rate_var_vint[mm,:] = dx*dy*var_int.isel(times=mm)
        elif var_type == 'sfc_frc':
            F_rate_var_vint[mm,:] = dx*dy*var_int[mm,:,:]

    print(F_rate_var_vint.shape)

    F_rate_delta_var_vint_mm = np.empty((365,len_rho_grid,560,630))

    for mm in np.arange(0,365):
        sigma_tmp = sigma_t[mm,:,:]
    
        #print(mm)
        for irho in np.arange(0,len_rho_grid):
    
            #print(irho)
            F_rate_tmp = ma.masked_where(np.logical_or(sigma_tmp <= (rho_grid[irho]-(0.05/2)),sigma_tmp > (rho_grid[irho]+(0.05/2))), F_rate_var_vint[mm,:,:])

            if irho == 0:
                F_rate_delta = F_rate_tmp.copy()
                F_rate_delta[np.logical_or(sigma_tmp <= (rho_grid[irho]-(0.05/2)),sigma_tmp > (rho_grid[irho]+(0.05/2)))] = np.nan
            elif irho == 1:
                F_rate_tmp[np.logical_or(sigma_tmp <= (rho_grid[irho]-(0.05/2)),sigma_tmp > (rho_grid[irho]+(0.05/2)))] = np.nan
                F_rate_delta = np.stack((F_rate_delta,F_rate_tmp), axis=0)
            else:
                F_rate_tmp[np.logical_or(sigma_tmp <= (rho_grid[irho]-(0.05/2)),sigma_tmp > (rho_grid[irho]+(0.05/2)))] = np.nan
                F_rate_extradim = np.expand_dims(F_rate_tmp, axis=0)
                F_rate_delta = np.concatenate((F_rate_delta,F_rate_extradim), axis=0)
            del F_rate_tmp
        
        F_rate_delta_var_vint_mm[mm,:] = F_rate_delta
    
    print('completed, size: ', F_rate_delta_var_vint_mm.shape)
    
    return F_rate_delta_var_vint_mm

# -- Units --
# Heat: m.degC/s -> m3.degC/s
# Salt: m/s -> m3/s
# Fwf: Kg.m-2.s-1 = Kg/s

# Shelf only: excluding open ocean

Fs_rate_delta_adv_vint_shelf_mm = wmt(salt_adv_full_vint*mask_shelf, dx, dy,'budget')
Fs_rate_delta_diff_vint_shelf_mm = wmt(R_s_vint*mask_shelf, dx, dy,'budget')
Fs_rate_delta_net_vint_shelf_mm = wmt(salt_net_full_vint*mask_shelf, dx, dy,'budget')
Fs_rate_delta_sfc_shelf_mm = wmt(salt_sfc*mask_shelf, dx, dy,'sfc_frc')

# decomponsing vert/horiz components:
Fs_rate_delta_vadv_vint_shelf_mm = wmt(salt_vdia_adv_full_vint*mask_shelf, dx, dy,'budget')
Fs_rate_delta_hadv_vint_shelf_mm = wmt(salt_hdia_adv_full_vint*mask_shelf, dx, dy,'budget')
Fs_rate_delta_vdiff_vint_shelf_mm = wmt(salt_vdia_diff_full_vint*mask_shelf, dx, dy,'budget')
Fs_rate_delta_hdiff_vint_shelf_mm = wmt(salt_hdia_diff_full_vint*mask_shelf, dx, dy,'budget')

Fh_rate_delta_adv_vint_shelf_mm = wmt(temp_adv_full_vint*mask_shelf, dx, dy,'budget')
Fh_rate_delta_diff_vint_shelf_mm = wmt(R_t_vint*mask_shelf, dx, dy,'budget')
Fh_rate_delta_net_vint_shelf_mm = wmt(temp_net_full_vint*mask_shelf, dx, dy,'budget')
Fh_rate_delta_sfc_shelf_mm = wmt(temp_sfc*mask_shelf, dx, dy,'sfc_frc')

# decomponsing vert/horiz components:
Fh_rate_delta_vadv_vint_shelf_mm = wmt(temp_vdia_adv_full_vint*mask_shelf, dx, dy,'budget')
Fh_rate_delta_hadv_vint_shelf_mm = wmt(temp_hdia_adv_full_vint*mask_shelf, dx, dy,'budget')
Fh_rate_delta_vdiff_vint_shelf_mm = wmt(temp_vdia_diff_full_vint*mask_shelf, dx, dy,'budget')
Fh_rate_delta_hdiff_vint_shelf_mm = wmt(temp_hdia_diff_full_vint*mask_shelf, dx, dy,'budget')

# development mask for density classes:
sigma_sept = sigma_t[8,0:50,0:50]

for irho in np.arange(0,len_rho_grid):
    print(irho, rho_grid[irho]-(0.05/2), rho_grid[irho]+(0.05/2), sigma_sept.shape)    
    icritst_tmp = ma.masked_where(np.logical_or(sigma_sept <= (rho_grid[irho]-(0.05/2)),sigma_sept > (rho_grid[irho]+(0.05/2))), sigma_sept)    
    
    if irho == 0:
        icritst = icritst_tmp.copy()
        icritst[np.logical_or(sigma_sept <= (rho_grid[irho]-(0.05/2)),sigma_sept > (rho_grid[irho]+(0.05/2)))] = np.nan
    elif irho == 1:
        icritst_tmp[np.logical_or(sigma_sept <= (rho_grid[irho]-(0.05/2)),sigma_sept > (rho_grid[irho]+(0.05/2)))] = np.nan
        icritst = np.stack((icritst,icritst_tmp), axis=0)
    else:
        icritst_tmp[np.logical_or(sigma_sept <= (rho_grid[irho]-(0.05/2)),sigma_sept > (rho_grid[irho]+(0.05/2)))] = np.nan
        icritst_extradim = np.expand_dims(icritst_tmp, axis=0)
        icritst = np.concatenate((icritst,icritst_extradim), axis=0)
    del icritst_tmp
    
print(icritst.shape)

# integrated over x, y directions: only continental shelf

Fs_rate_sfc_shelf_mm_int = np.empty((len(rho_grid),365))
Fh_rate_sfc_shelf_mm_int = np.empty((len(rho_grid),365))

Fs_rate_adv_vint_shelf_mm_int = np.empty((len(rho_grid),365))
Fh_rate_adv_vint_shelf_mm_int = np.empty((len(rho_grid),365))
Fs_rate_diff_vint_shelf_mm_int = np.empty((len(rho_grid),365))
Fh_rate_diff_vint_shelf_mm_int = np.empty((len(rho_grid),365))
Fs_rate_net_vint_shelf_mm_int = np.empty((len(rho_grid),365))
Fh_rate_net_vint_shelf_mm_int = np.empty((len(rho_grid),365))

Fs_rate_vadv_vint_shelf_mm_int = np.empty((len(rho_grid),365))
Fh_rate_vadv_vint_shelf_mm_int = np.empty((len(rho_grid),365))
Fs_rate_vdiff_vint_shelf_mm_int = np.empty((len(rho_grid),365))
Fh_rate_vdiff_vint_shelf_mm_int = np.empty((len(rho_grid),365))
Fs_rate_hadv_vint_shelf_mm_int = np.empty((len(rho_grid),365))
Fh_rate_hadv_vint_shelf_mm_int = np.empty((len(rho_grid),365))
Fs_rate_hdiff_vint_shelf_mm_int = np.empty((len(rho_grid),365))
Fh_rate_hdiff_vint_shelf_mm_int = np.empty((len(rho_grid),365))

for irho in np.arange(0,len(rho_grid)):   
    for mm in np.arange(0,365):
        
        Fs_rate_sfc_shelf_mm_int[irho,mm] = np.nansum(np.nansum(Fs_rate_delta_sfc_shelf_mm[mm,irho,:], axis=1), axis=0)
        Fh_rate_sfc_shelf_mm_int[irho,mm] = np.nansum(np.nansum(Fh_rate_delta_sfc_shelf_mm[mm,irho,:], axis=1), axis=0)
        
        Fs_rate_adv_vint_shelf_mm_int[irho,mm] = np.nansum(np.nansum(Fs_rate_delta_adv_vint_shelf_mm[mm,irho,:], axis=1), axis=0)
        Fh_rate_adv_vint_shelf_mm_int[irho,mm] = np.nansum(np.nansum(Fh_rate_delta_adv_vint_shelf_mm[mm,irho,:], axis=1), axis=0)

        Fs_rate_vadv_vint_shelf_mm_int[irho,mm] = np.nansum(np.nansum(Fs_rate_delta_vadv_vint_shelf_mm[mm,irho,:], axis=1), axis=0)
        Fh_rate_vadv_vint_shelf_mm_int[irho,mm] = np.nansum(np.nansum(Fh_rate_delta_vadv_vint_shelf_mm[mm,irho,:], axis=1), axis=0)
        Fs_rate_hadv_vint_shelf_mm_int[irho,mm] = np.nansum(np.nansum(Fs_rate_delta_hadv_vint_shelf_mm[mm,irho,:], axis=1), axis=0)
        Fh_rate_hadv_vint_shelf_mm_int[irho,mm] = np.nansum(np.nansum(Fh_rate_delta_hadv_vint_shelf_mm[mm,irho,:], axis=1), axis=0)
        
        Fs_rate_diff_vint_shelf_mm_int[irho,mm] = np.nansum(np.nansum(Fs_rate_delta_diff_vint_shelf_mm[mm,irho,:], axis=1), axis=0)
        Fh_rate_diff_vint_shelf_mm_int[irho,mm] = np.nansum(np.nansum(Fh_rate_delta_diff_vint_shelf_mm[mm,irho,:], axis=1), axis=0)

        Fs_rate_vdiff_vint_shelf_mm_int[irho,mm] = np.nansum(np.nansum(Fs_rate_delta_vdiff_vint_shelf_mm[mm,irho,:], axis=1), axis=0)
        Fh_rate_vdiff_vint_shelf_mm_int[irho,mm] = np.nansum(np.nansum(Fh_rate_delta_vdiff_vint_shelf_mm[mm,irho,:], axis=1), axis=0)
        Fs_rate_hdiff_vint_shelf_mm_int[irho,mm] = np.nansum(np.nansum(Fs_rate_delta_hdiff_vint_shelf_mm[mm,irho,:], axis=1), axis=0)
        Fh_rate_hdiff_vint_shelf_mm_int[irho,mm] = np.nansum(np.nansum(Fh_rate_delta_hdiff_vint_shelf_mm[mm,irho,:], axis=1), axis=0)
        
        Fs_rate_net_vint_shelf_mm_int[irho,mm] = np.nansum(np.nansum(Fs_rate_delta_net_vint_shelf_mm[mm,irho,:], axis=1), axis=0)
        Fh_rate_net_vint_shelf_mm_int[irho,mm] = np.nansum(np.nansum(Fh_rate_delta_net_vint_shelf_mm[mm,irho,:], axis=1), axis=0)

# figures
fig_path = '/users/boeiradi/COLD_project/postprocessing/figs/WMT/'

# plot with bars
width=.023

# convert to rate per year:
Dt = 1000/0.05

# divide for Dt (sec in 1 year) and 10^6 to obtain rate in Sv year-1 

# SFC: continental shelf
Fs_sig_sfc_shelf_mm = -Fs_rate_sfc_shelf_mm_int*Dt/1e6
Fh_sig_sfc_shelf_mm = -Fh_rate_sfc_shelf_mm_int*Dt/1e6 
F_sig_sfc_shelf_mm = -Fs_sig_sfc_shelf_mm + Fh_sig_sfc_shelf_mm
# - calculate the anual average of the monthly ars:
Fs_sig_sfc_shelf =  np.nanmean(Fs_sig_sfc_shelf_mm, axis=1)
Fh_sig_sfc_shelf =  np.nanmean(Fh_sig_sfc_shelf_mm, axis=1)
F_sig_sfc_shelf = -Fs_sig_sfc_shelf + Fh_sig_sfc_shelf

# plots:


fig = plt.figure(figsize=(12,10))


ax3 = fig.add_subplot(223)
plt.title('Continental shelf')
cf=plt.plot(rho_grid,-Fs_sig_sfc_shelf,'--b',label='WMT due to sfc salt flux')
ch=plt.plot(rho_grid,Fh_sig_sfc_shelf,'r',label='WMT due to sfc heat flux')
ct=plt.plot(rho_grid,F_sig_sfc_shelf,'k',label='WMT due to sfc total flux')
plt.legend()
plt.xlim(26.5,28),plt.ylim(-5,5)
plt.grid(True)
plt.text(26.6,1.5,'Buoyancy loss')
plt.text(26.6,-1.5,'Buoyancy gain')

name_fig="waom10extend_shflim_S_0.25Q_WMT_Full_sfc_fluxes_annual_daily_shelf.png"
plt.savefig(fig_path + name_fig, dpi=300)

# ADV (vint)
Fs_sig_adv_vint_shelf_mm = -Fs_rate_adv_vint_shelf_mm_int*Dt/1e6 
Fh_sig_adv_vint_shelf_mm = -Fh_rate_adv_vint_shelf_mm_int*Dt/1e6 
F_sig_adv_vint_shelf_mm = -Fs_sig_adv_vint_shelf_mm + Fh_sig_adv_vint_shelf_mm
# - calculate the anual average of the monthly ars:
Fs_sig_adv_vint_shelf =  np.nanmean(Fs_sig_adv_vint_shelf_mm, axis=1)
Fh_sig_adv_vint_shelf =  np.nanmean(Fh_sig_adv_vint_shelf_mm, axis=1)
F_sig_adv_vint_shelf = -Fs_sig_adv_vint_shelf + Fh_sig_adv_vint_shelf

# Vert adv
Fs_sig_vadv_vint_shelf_mm = -Fs_rate_vadv_vint_shelf_mm_int*Dt/1e6
Fh_sig_vadv_vint_shelf_mm = -Fh_rate_vadv_vint_shelf_mm_int*Dt/1e6
# - calculate the anual average of the monthly ars:
Fs_sig_vadv_vint_shelf =  np.nanmean(Fs_sig_vadv_vint_shelf_mm, axis=1)
Fh_sig_vadv_vint_shelf =  np.nanmean(Fh_sig_vadv_vint_shelf_mm, axis=1)

# Horz adv
Fs_sig_hadv_vint_shelf_mm = -Fs_rate_hadv_vint_shelf_mm_int*Dt/1e6
Fh_sig_hadv_vint_shelf_mm = -Fh_rate_hadv_vint_shelf_mm_int*Dt/1e6
# - calculate the anual average of the monthly ars:
Fs_sig_hadv_vint_shelf =  np.nanmean(Fs_sig_hadv_vint_shelf_mm, axis=1)
Fh_sig_hadv_vint_shelf =  np.nanmean(Fh_sig_hadv_vint_shelf_mm, axis=1)

# DIFF
Fs_sig_diff_vint_shelf_mm = -Fs_rate_diff_vint_shelf_mm_int*Dt/1e6 
Fh_sig_diff_vint_shelf_mm = -Fh_rate_diff_vint_shelf_mm_int*Dt/1e6 
F_sig_diff_vint_shelf_mm = -Fs_sig_diff_vint_shelf_mm + Fh_sig_diff_vint_shelf_mm
# - calculate the anual average of the monthly ars:
Fs_sig_diff_vint_shelf =  np.nanmean(Fs_sig_diff_vint_shelf_mm, axis=1)
Fh_sig_diff_vint_shelf =  np.nanmean(Fh_sig_diff_vint_shelf_mm, axis=1)
F_sig_diff_vint_shelf = -Fs_sig_diff_vint_shelf + Fh_sig_diff_vint_shelf

# vert diff
Fs_sig_vdiff_vint_shelf_mm = -Fs_rate_vdiff_vint_shelf_mm_int*Dt/1e6
Fh_sig_vdiff_vint_shelf_mm = -Fh_rate_vdiff_vint_shelf_mm_int*Dt/1e6
# - calculate the anual average of the monthly ars:
Fs_sig_vdiff_vint_shelf =  np.nanmean(Fs_sig_vdiff_vint_shelf_mm, axis=1)
Fh_sig_vdiff_vint_shelf =  np.nanmean(Fh_sig_vdiff_vint_shelf_mm, axis=1)

# horz diff
Fs_sig_hdiff_vint_shelf_mm = -Fs_rate_hdiff_vint_shelf_mm_int*Dt/1e6
Fh_sig_hdiff_vint_shelf_mm = -Fh_rate_hdiff_vint_shelf_mm_int*Dt/1e6
# - calculate the anual average of the monthly ars:
Fs_sig_hdiff_vint_shelf =  np.nanmean(Fs_sig_hdiff_vint_shelf_mm, axis=1)
Fh_sig_hdiff_vint_shelf =  np.nanmean(Fh_sig_hdiff_vint_shelf_mm, axis=1)

# NET
Fs_sig_net_vint_shelf_mm = -Fs_rate_net_vint_shelf_mm_int*Dt/1e6 
Fh_sig_net_vint_shelf_mm = -Fh_rate_net_vint_shelf_mm_int*Dt/1e6 
F_sig_net_vint_shelf_mm = -Fs_sig_net_vint_shelf_mm + Fh_sig_net_vint_shelf_mm
# - calculate the anual average of the monthly ars:
Fs_sig_net_vint_shelf =  np.nanmean(Fs_sig_net_vint_shelf_mm, axis=1)
Fh_sig_net_vint_shelf =  np.nanmean(Fh_sig_net_vint_shelf_mm, axis=1)
F_sig_net_vint_shelf = -Fs_sig_net_vint_shelf + Fh_sig_net_vint_shelf

# obs (16.6.21): salt budget has opposite effect to buoyancy than fw budget,
#  so term are inverted (x-1) to plot the transformation rates. This seems to agree
#  better with the surface fluxes.
# obs (21.6.21): only inverted sign is not correct. Need to do full conversion as done from 
#  salt flux to fw flux, i.e. multiply by -1000 (fw density) and divide by salinity (avg in the MLD)

fig = plt.figure(figsize=(16,10))
ax1 = fig.add_subplot(222)
plt.title('Salt budget')
cs=plt.plot(rho_grid,-Fs_sig_net_vint_shelf,'b',label='WMT due to salt tend')
cs=plt.plot(rho_grid,-Fs_sig_diff_vint_shelf,'--b',label='WMT due to salt diff')
cs=plt.plot(rho_grid,-Fs_sig_adv_vint_shelf,':b',label='WMT due to salt adv')
plt.legend()
plt.ylabel('Transformation rate (Sv)')
plt.xlim(26.5,28),plt.ylim(-5,5)
plt.grid(True)
plt.text(26.6,1.5,'Buoyancy loss')
plt.text(26.6,-1.5,'Buoyancy gain')
ax1 = fig.add_subplot(223)
plt.title('Heat budget')
ch=plt.plot(rho_grid,Fh_sig_net_vint_shelf,'r',label='WMT due to heat tend')
ch=plt.plot(rho_grid,Fh_sig_diff_vint_shelf,'--r',label='WMT due to heat diff')
ch=plt.plot(rho_grid,Fh_sig_adv_vint_shelf,':r',label='WMT due to heat adv')
plt.legend()
plt.ylabel('Transformation rate (Sv)')
plt.xlim(26.5,28),plt.ylim(-5,5)
plt.grid(True)

ax2 = fig.add_subplot(221)
plt.title('Surface fluxes')
cs=plt.plot(rho_grid,-Fs_sig_sfc_shelf,'b',label='WMT due to salt sfc')
ch=plt.plot(rho_grid,Fh_sig_sfc_shelf,'r',label='WMT due to heat sfc')
ct=plt.plot(rho_grid,F_sig_sfc_shelf,'k',label='WMT due to total sfc')
plt.legend()
plt.ylabel('Transformation rate (Sv)')
plt.xlim(26.5,28),plt.ylim(-5,5)
plt.grid(True)

ax2 = fig.add_subplot(224)
plt.title('Net tendencies')
cs=plt.plot(rho_grid,-Fs_sig_net_vint_shelf,'-b',label='WMT tendency due to salt sfc')
ch=plt.plot(rho_grid,Fh_sig_net_vint_shelf,'-r',label='WMT tendency due to heat sfc')
ct=plt.plot(rho_grid,Fh_sig_net_vint_shelf-Fs_sig_net_vint_shelf,'k',label='WMT tendency total')
plt.legend()
plt.ylabel('Transformation rate (Sv)')
plt.xlim(26.5,28),plt.ylim(-5,5)
plt.grid(True)
plt.title('Total tendencies')

name_fig="waom10extend_shflim_S_0.25Q_WMT_Full_heat-salt_vint_annual_daily_shelf.png"
plt.savefig(fig_path + name_fig, dpi=300)

### plot some maps
import matplotlib.path as mpath
import cartopy.feature as cfeature

def lonlat_labels(ax):
    # latitude labels
    ax.text(120,-80,'80$^{\circ}$S',transform=ccrs.PlateCarree(),color='gray')
    ax.text(120,-70,'70$^{\circ}$S',transform=ccrs.PlateCarree(),color='gray')
    # longitude labels
    ax.text(0,-66,'0$^{\circ}$',transform=ccrs.PlateCarree(),color='gray')
    #ax.text(60,-53,'60$^{\circ}$E',transform=ccrs.PlateCarree(),color='gray')
    #ax.text(120,-53,'120$^{\circ}$E',transform=ccrs.PlateCarree(),color='gray')
    ax.text(-60,-48,'60$^{\circ}$W',transform=ccrs.PlateCarree(),color='gray')
    ax.text(-120,-48,'120$^{\circ}$W',transform=ccrs.PlateCarree(),color='gray')
    ax.text(180,-60,'180$^{\circ}$',transform=ccrs.PlateCarree(),color='gray')
    return

proj = ccrs.SouthPolarStereo()

# plot maps for 27.5 kg.m-3 isopycnal
print(Fs_rate_delta_net_vint_shelf_mm.shape)

for irho in np.arange(29,35):#15,37):

# call cartopy projection
    proj = ccrs.SouthPolarStereo()
    fig = plt.figure(figsize=(16,10))
    ax1 = fig.add_subplot(221, projection=proj)
    plt.title('')
    cy=plt.pcolormesh(lon_rho,lat_rho,sigma_t[4,:,:], transform=ccrs.PlateCarree())
    plt.colorbar(cy)
    plt.clim(25.,28.)
    ax1.gridlines() # draw_labels=True,linewidth=2, color='white', alpha=0.5, linestyle='--')
    lonlat_labels(ax1)
    ax1.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
    ax1.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white') 

    ax2 = fig.add_subplot(222, projection=proj)
    plt.title('FWF, $\sigma_{\Theta}$ = ' + str(rho_grid[irho]))
    cy=plt.pcolormesh(lon_rho,lat_rho,np.nanmean(Fs_rate_delta_net_vint_shelf_mm[:,irho,:], axis=0)*Dt/(dx*dy), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm)
    plt.colorbar(cy)
    ax2.gridlines() # draw_labels=True,linewidth=2, color='white', alpha=0.5, linestyle='--')
    lonlat_labels(ax2)
    ax2.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
    ax2.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white') 
    plt.clim(-5e-5,5e-5)

    ax3 = fig.add_subplot(223, projection=proj)
    plt.title('HF, $\sigma_{\Theta}$ = ' + str(rho_grid[irho]))
    cy=plt.pcolormesh(lon_rho,lat_rho,np.nanmean(Fh_rate_delta_net_vint_shelf_mm[:,irho,:], axis=0)*Dt/(dx*dy), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm)
    plt.colorbar(cy)
    ax3.gridlines() # draw_labels=True,linewidth=2, color='white', alpha=0.5, linestyle='--')
    lonlat_labels(ax3)
    ax3.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
    ax3.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white') 
    plt.clim(-5e-5,5e-5)
                                                
    ax4 = fig.add_subplot(224, projection=proj)
    plt.title('Total, $\sigma_{\Theta}$ = ' + str(rho_grid[irho]))
    cy=plt.pcolormesh(lon_rho,lat_rho,np.nanmean(Fh_rate_delta_net_vint_shelf_mm[:,irho,:], axis=0)*Dt/(dx*dy)+np.nanmean(Fs_rate_delta_net_vint_shelf_mm[:,irho,:], axis=0)*Dt/(dx*dy), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm)
    plt.colorbar(cy)
    ax4.gridlines() # draw_labels=True,linewidth=2, color='white', alpha=0.5, linestyle='--')
    lonlat_labels(ax4)
    ax4.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
    ax4.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white') 
    plt.clim(-5e-5,5e-5)
                                                
    name_fig="waom10extend_shflim_S_0.25Q_WMTmaps_Full_annual_yr20_l" + str(irho)  + "_daily_shelf.png"
    plt.savefig(fig_path + name_fig, dpi=300)
    plt.close()

proj = ccrs.SouthPolarStereo()
fig = plt.figure(figsize=(16,10))

irho=30
ax1 = fig.add_subplot(221, projection=proj)
plt.title('Total, $\sigma_{\Theta}$ = ' + str(rho_grid[irho]))
cy=plt.pcolormesh(lon_rho,lat_rho,np.nanmean(Fh_rate_delta_net_vint_shelf_mm[:,irho,:], axis=0)*Dt/(dx*dy)+np.nanmean(Fs_rate_delta_net_vint_shelf_mm[:,irho,:], axis=0)*Dt/(dx*dy), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm)
plt.colorbar(cy, extend='both')
ax1.gridlines() # draw_labels=True,linewidth=2, color='white', alpha=0.5, linestyle='--')
lonlat_labels(ax1)
ax1.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
ax1.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white') 
plt.clim(-5e-5,5e-5)

irho=31
ax2 = fig.add_subplot(222, projection=proj)
plt.title('Total, $\sigma_{\Theta}$ = ' + str(rho_grid[irho]))
cy=plt.pcolormesh(lon_rho,lat_rho,np.nanmean(Fh_rate_delta_net_vint_shelf_mm[:,irho,:], axis=0)*Dt/(dx*dy)+np.nanmean(Fs_rate_delta_net_vint_shelf_mm[:,irho,:], axis=0)*Dt/(dx*dy), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm)
plt.colorbar(cy, extend='both')
ax2.gridlines() # draw_labels=True,linewidth=2, color='white', alpha=0.5, linestyle='--')
lonlat_labels(ax2)
ax2.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
ax2.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white') 
plt.clim(-5e-5,5e-5)

irho=32
ax3 = fig.add_subplot(223, projection=proj)
plt.title('Total, $\sigma_{\Theta}$ = ' + str(rho_grid[irho]))
cy=plt.pcolormesh(lon_rho,lat_rho,np.nanmean(Fh_rate_delta_net_vint_shelf_mm[:,irho,:], axis=0)*Dt/(dx*dy)+np.nanmean(Fs_rate_delta_net_vint_shelf_mm[:,irho,:], axis=0)*Dt/(dx*dy), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm)
plt.colorbar(cy, extend='both')
ax3.gridlines() # draw_labels=True,linewidth=2, color='white', alpha=0.5, linestyle='--')
lonlat_labels(ax3)
ax3.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
ax3.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white') 
plt.clim(-5e-5,5e-5)

irho=33
ax4 = fig.add_subplot(224, projection=proj)
plt.title('Total, $\sigma_{\Theta}$ = ' + str(rho_grid[irho]))
cy=plt.pcolormesh(lon_rho,lat_rho,np.nanmean(Fh_rate_delta_net_vint_shelf_mm[:,irho,:], axis=0)*Dt/(dx*dy)+np.nanmean(Fs_rate_delta_net_vint_shelf_mm[:,irho,:], axis=0)*Dt/(dx*dy), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm)
plt.colorbar(cy, extend='both')
ax4.gridlines() # draw_labels=True,linewidth=2, color='white', alpha=0.5, linestyle='--')
lonlat_labels(ax4)
ax4.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
ax4.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white') 
plt.clim(-5e-5,5e-5)
                                                
name_fig="waom10extend_shflim_S_0.25Q_WMTmaps_Full_annual_daily_yr20_shelf.png"
plt.savefig(fig_path + name_fig, dpi=300)
plt.close()

# --- Save transformation arrays:
npy_path = '/users/boeiradi/COLD_project/postprocessing/tmp_files/'

print('Saving files .....')

np.savez(npy_path + 'WAOM10extend_daily_Full/gc_WMT_shelf_sfc_total',rho_grid = rho_grid, F_sig_sfc_shelf = F_sig_sfc_shelf) 
np.savez(npy_path + 'WAOM10extend_daily_Full/gc_WMT_shelf_sfc_salt',rho_grid = rho_grid, Fs_sig_sfc_shelf = Fs_sig_sfc_shelf) 
np.savez(npy_path + 'WAOM10extend_daily_Full/gc_WMT_shelf_sfc_heat',rho_grid = rho_grid, Fh_sig_sfc_shelf = Fh_sig_sfc_shelf) 

np.savez(npy_path + 'WAOM10extend_daily_Full/gc_WMT_shelf_vint_salt_net',rho_grid = rho_grid, Fs_sig_net_vint_shelf = Fs_sig_net_vint_shelf) 
np.savez(npy_path + 'WAOM10extend_daily_Full/gc_WMT_shelf_vint_salt_diff',rho_grid = rho_grid, Fs_sig_diff_vint_shelf = Fs_sig_diff_vint_shelf,\
        Fs_sig_vdiff_vint_shelf = Fs_sig_vdiff_vint_shelf, Fs_sig_hdiff_vint_shelf = Fs_sig_hdiff_vint_shelf)
np.savez(npy_path + 'WAOM10extend_daily_Full/gc_WMT_shelf_vint_salt_adv',rho_grid = rho_grid, Fs_sig_adv_vint_shelf = Fs_sig_adv_vint_shelf,\
        Fs_sig_vadv_vint_shelf = Fs_sig_vadv_vint_shelf, Fs_sig_hadv_vint_shelf = Fs_sig_hadv_vint_shelf)

np.savez(npy_path + 'WAOM10extend_daily_Full/gc_WMT_shelf_vint_heat_net', rho_grid = rho_grid, Fh_sig_net_vint_shelf = Fh_sig_net_vint_shelf)
np.savez(npy_path + 'WAOM10extend_daily_Full/gc_WMT_shelf_vint_heat_diff', rho_grid = rho_grid, Fh_sig_diff_vint_shelf = Fh_sig_diff_vint_shelf,\
        Fh_sig_vdiff_vint_shelf = Fh_sig_vdiff_vint_shelf, Fh_sig_hdiff_vint_shelf = Fh_sig_hdiff_vint_shelf)
np.savez(npy_path + 'WAOM10extend_daily_Full/gc_WMT_shelf_vint_heat_adv', rho_grid = rho_grid, Fh_sig_adv_vint_shelf = Fh_sig_adv_vint_shelf,\
        Fh_sig_vadv_vint_shelf = Fh_sig_vadv_vint_shelf, Fh_sig_hadv_vint_shelf = Fh_sig_hadv_vint_shelf)
