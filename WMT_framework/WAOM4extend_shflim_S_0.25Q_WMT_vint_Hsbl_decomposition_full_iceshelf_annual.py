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

#import iris
#import iris.iterate
#import iris.coords
#import iris.plot as iplt
import gsw

# load ROMS avg output
for mm  in ['01','02','03','04','05','06','07','08','09','10','11','12']:
    ds = xr.open_dataset('/scratch/project_2000789/boeiradi/waom4extend_shflim_S_0.25Q/output_yr10_diag/ocean_avg_00' + mm + '.nc')
    print(ds.variables["temp"].shape)
    temp_tmp = np.nanmean(ds.variables["temp"], axis=0)
    salt_tmp = np.nanmean(ds.variables["salt"], axis=0)
    shflux_tmp = np.nanmean(ds.variables["shflux"], axis=0)
    ssflux_tmp = np.nanmean(ds.variables["ssflux"], axis=0)
    m_tmp = np.nanmean(ds.variables["m"], axis=0)
    
    ds = ds.set_coords(['Cs_r', 'Cs_w', 'hc', 'h', 'Vtransform'])
    if ds.Vtransform == 1:
        Zo_rho = ds.hc * (ds.s_rho - ds.Cs_r) + ds.Cs_r * ds.h
        z_rho_tmp = Zo_rho + ds.zeta * (1 + Zo_rho/ds.h)
        print("Vtransform=1")
    elif ds.Vtransform == 2:
        Zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * ds.h) / (ds.hc + ds.h)
        z_rho_tmp = ds.zeta + (ds.zeta + ds.h) * Zo_rho + ds.zice
        print("Vtransform=2")
        Zo_w = (ds.hc * ds.s_w + ds.Cs_w * ds.h) / (ds.hc + ds.h)
        z_w_tmp = ds.zeta + (ds.zeta + ds.h) * Zo_w + ds.zice
        
    z_rho_avg = np.nanmean(z_rho_tmp, axis=0)
    z_w_avg = np.nanmean(z_w_tmp,axis=0)
    
    # concatanate monthly avgs into a yearly variable
    if mm == '01':
        temp = temp_tmp
        salt = salt_tmp
        shflux = shflux_tmp
        ssflux = ssflux_tmp
        m = m_tmp
        z_rho = z_rho_avg
        z_w = z_w_avg
    elif mm == '02':
        temp = np.stack((temp,temp_tmp), axis=0)
        salt = np.stack((salt,salt_tmp), axis=0)
        shflux = np.stack((shflux,shflux_tmp), axis=0)
        ssflux = np.stack((ssflux,ssflux_tmp), axis=0)
        m = np.stack((m,m_tmp), axis=0)
        z_rho = np.stack((z_rho,z_rho_avg), axis=0)
        z_w = np.stack((z_w,z_w_avg), axis=0)
    else:
        temp_tmp_4thdim = np.expand_dims(temp_tmp, axis=0)
        temp = np.concatenate((temp,temp_tmp_4thdim), axis=0)
        salt_tmp_4thdim = np.expand_dims(salt_tmp, axis=0)
        salt = np.concatenate((salt,salt_tmp_4thdim), axis=0)
        shflux_tmp_4thdim = np.expand_dims(shflux_tmp, axis=0)
        shflux = np.concatenate((shflux,shflux_tmp_4thdim), axis=0)    
        ssflux_tmp_4thdim = np.expand_dims(ssflux_tmp, axis=0)
        ssflux = np.concatenate((ssflux,ssflux_tmp_4thdim), axis=0) 
        m_tmp_4thdim = np.expand_dims(m_tmp, axis=0)
        m = np.concatenate((m,m_tmp_4thdim), axis=0) 
        z_rho_tmp_4thdim = np.expand_dims(z_rho_avg, axis=0)
        z_rho = np.concatenate((z_rho,z_rho_tmp_4thdim), axis=0) 
        z_w_tmp_4thdim = np.expand_dims(z_w_avg, axis=0)
        z_w = np.concatenate((z_w,z_w_tmp_4thdim), axis=0) 
#ds.coords['flux'] = flux#.transpose() # put flux into ds dataset

    ds.close()
    
sigma_t_sfc = gsw.rho(salt[:,-1,:,:],temp[:,-1,:,:],0) - 1000

di = xr.open_dataset('/scratch/project_2000789/boeiradi/waom4extend_shflim_S_0.25Q/output_yr10_diag/ocean_avg_0001.nc')
ice_draft = di.variables["zice"]

mask_zice = ma.masked_where(ice_draft < 0, np.ones(ice_draft.shape))

di.close()

# calculate dz following:
dz = np.empty((12,1400,1575,31))
dz_inv = np.empty((12,1400,1575,31))


for tt in np.arange(0,12):
    z_w_sorted = -1*z_w[tt,:,:,::-1]
    print(z_w_sorted.shape)
    dz_inv[tt,:,:,:] = np.diff(z_w_sorted,axis=2)
    dz[tt,:,:,:] = dz_inv[tt,:,:,::-1]

dg = xr.open_dataset("/scratch/project_2000339/boeiradi/waom4_frc/waom4extend_grd.nc")

lat_rho = dg.variables["lat_rho"]
lon_rho = dg.variables["lon_rho"]
pm = dg.variables["pm"]
pn = dg.variables["pn"]

ds.coords['lat_rho']=lat_rho.transpose() # put lat_rho into ds dataset
ds.coords['lon_rho']=lon_rho.transpose() # put lon_rho into ds dataset

# shelf/open-ocean masks:
mask_open = ma.masked_where(dg.h <= 2000, np.ones(dg.h.shape))
mask_shelf = ma.masked_where(dg.h > 2000, np.ones(dg.h.shape))


#mask_shelf = np.empty((dg.h.shape))
#mask_open = np.empty((dg.h.shape))

#open_ind=ds.h.where(dg.h > 1000)
#shelf_ind=ds.h.where(dg.h <= 1000)
#print(open_ind)

#mask_shelf = np.divide(shelf_ind,shelf_ind)
#mask_open = np.divide(open_ind,open_ind)

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
ax4 = fig.add_subplot(224)#, projection=proj)
cy=plt.pcolor(np.squeeze(dz[0,:,:,-1]))#, transform=ccrs.PlateCarree())
plt.colorbar(cy)

dx = xr.open_dataset('/scratch/project_2000789/boeiradi/waom4extend_shflim_S_0.25Q/output_yr10_diag/MLD_vint_vars_for_WMT.nc')
print(dx.variables["sfc_mld"].shape)
sfc_mld = dx.variables["sfc_mld"]
# - variables across ML base
temp_avg_mld = dx.variables["temp_avg_mld"]
salt_avg_mld = dx.variables["salt_avg_mld"]

# - variables integrated throughout the ML; multiply by -1 b/c dz is negative.
temp_vdia_diff_mld_vint = dx.variables["temp_vdia_diff_mld_vint"]
salt_vdia_diff_mld_vint = dx.variables["salt_vdia_diff_mld_vint"]
temp_hdia_diff_mld_vint = dx.variables["temp_hdia_diff_mld_vint"]
salt_hdia_diff_mld_vint = dx.variables["salt_hdia_diff_mld_vint"]
temp_vdia_adv_mld_vint = dx.variables["temp_vdia_adv_mld_vint"]
salt_vdia_adv_mld_vint = dx.variables["salt_vdia_adv_mld_vint"]
temp_hdia_adv_mld_vint = dx.variables["temp_hdia_adv_mld_vint"]
salt_hdia_adv_mld_vint = dx.variables["salt_hdia_adv_mld_vint"]
temp_tend_mld_vint = dx.variables["temp_tend_avg_mld_vint"]
salt_tend_mld_vint = dx.variables["salt_tend_avg_mld_vint"]

sigma_t = gsw.rho(salt_avg_mld,temp_avg_mld,0) - 1000

dx.close()

# obtain thermal expansion (alpha) & salinity contraction (beta) coefficients:
SA_avg_mld = np.empty(salt_avg_mld.shape)
# neet Absolute Salinity, converting from Pratical Salinity:
for mm in np.arange(0,12):
    SA_tmp =gsw.SA_from_SP(np.squeeze(salt_avg_mld[mm,:,:]),0,lon_rho,lat_rho)
    SA_avg_mld[mm,:,:] = SA_tmp
    del SA_tmp

# gsw.alpha/gsw.beta
#alpha = gsw_alpha(SA,CT,p)
[specvol_mld, alpha_mld, beta_mld] = gsw.specvol_alpha_beta(SA_avg_mld,temp_avg_mld,0)

print(alpha_mld.shape)

fig = plt.figure(figsize=(14,12))
ax1 = fig.add_subplot(221)#, projection=proj)
cy=plt.pcolor(np.nanmean(specvol_mld, axis=0))
plt.colorbar(cy)
#plt.clim(0.,1.)
ax2 = fig.add_subplot(222)#, projection=proj)
cy=plt.pcolor(np.nanmean(alpha_mld, axis=0))
plt.colorbar(cy)
#plt.clim(0.,1.)
ax3 = fig.add_subplot(223)#, projection=proj)
cy=plt.pcolor(np.nanmean(beta_mld, axis=0))
plt.colorbar(cy)
#plt.clim(0.,1.)
ax4 = fig.add_subplot(224)#, projection=proj)
cy=plt.pcolor(np.nanmean(beta_mld*specvol_mld, axis=0))
plt.colorbar(cy)

# calculate the LHS term in Pellichero et al (2018):
# ps: Diffusion (R_s, R_t) terms already include the sfc fluxes

# heat (eqn 5)
rho0 = 1025 #1000
Cp = 3985

# convert salt convergenvces to equivalent FW
#fwf_sfc = (np.divide(ssflux*rho0, -salt[:,-1,:,:]))*specvol_mld
fwf_sfc = (np.divide(beta_mld*ssflux*rho0*dz[:,:,:,-1], -salt[:,-1,:,:]))*specvol_mld

fwf_hdia_adv_mld_vint = np.divide(salt_hdia_adv_mld_vint*rho0, -salt_avg_mld[:,:,:])*specvol_mld
fwf_vdia_adv_mld_vint = np.divide(salt_vdia_adv_mld_vint*rho0, -salt_avg_mld[:,:,:])*specvol_mld
fwf_hdia_diff_mld_vint = np.divide(salt_hdia_diff_mld_vint*rho0, -salt_avg_mld[:,:,:])*specvol_mld
fwf_vdia_diff_mld_vint = np.divide(salt_vdia_diff_mld_vint*rho0, -salt_avg_mld[:,:,:])*specvol_mld

# total diffusion terms:
R_s_vint = beta_mld*(salt_hdia_diff_mld_vint + salt_vdia_diff_mld_vint)
R_t_vint = alpha_mld*(temp_hdia_diff_mld_vint + temp_vdia_diff_mld_vint)
R_f_vint = fwf_hdia_diff_mld_vint + fwf_vdia_diff_mld_vint

# surface flux terms:
salt_sfc = beta_mld*(ssflux)
temp_sfc = alpha_mld*(np.divide(shflux, rho0*Cp))

# advection terms:
salt_adv_mld_vint = beta_mld*(salt_hdia_adv_mld_vint + salt_vdia_adv_mld_vint)
temp_adv_mld_vint = alpha_mld*(temp_hdia_adv_mld_vint + temp_vdia_adv_mld_vint)
fwf_adv_mld_vint = fwf_hdia_adv_mld_vint + fwf_vdia_adv_mld_vint

# net tendencies
salt_net_mld_vint = beta_mld*salt_tend_mld_vint
temp_net_mld_vint = alpha_mld*temp_tend_mld_vint
fwf_net_mld_vint = (np.divide(salt_tend_mld_vint*rho0, -salt_avg_mld[:,:,:]))*specvol_mld

# ML budget equation: 
# salt:
# salt_sfc - R_s = -salt_tend_mld + salt_hdia_adv_mld + salt_vdia_adv_mld

# temp:
# temp_sfc - R_t = -temp_tend_mld + temp_hdia_adv_mld + temp_vdia_adv_mld

#  Function to calculate Water Mass Transformation (in m3/s):

# rho grid for binning:
#rho_grid=np.arange(35.5,37.4,0.1) # for sigma-2
rho_grid=np.arange(24.4,29.1,0.1) # for sigma-0
len_rho_grid=len(rho_grid)

dx = np.divide(1,pm)
dy = np.divide(1,pn)
dt = 86400#30#/12 #why divide by 12?

def wmt(var_int, dx, dy,var_type):
    # var_type: 'budget' or 'sfc_frc'
    
    F_rate_var_vint = np.empty(var_int.shape)
    
    for mm in np.arange(0,12):
        if var_type == 'budget':
            F_rate_var_vint[mm,:] = dx*dy*var_int.isel(times=mm)
        elif var_type == 'sfc_frc':
            F_rate_var_vint[mm,:] = dx*dy*var_int[mm,:,:]

    print(F_rate_var_vint.shape)

    F_rate_delta_var_vint_mm = np.empty((12,len_rho_grid,1400,1575))

    for mm in np.arange(0,12):
        sigma_tmp = sigma_t[mm,:,:]
    
        #print(mm)
        for irho in np.arange(0,len_rho_grid):
    
            #print(irho)
            F_rate_tmp = ma.masked_where(np.logical_or(sigma_tmp <= (rho_grid[irho]-(0.1/2)),sigma_tmp > (rho_grid[irho]+(0.1/2))), F_rate_var_vint[mm,:,:])

            if irho == 0:
                F_rate_delta = F_rate_tmp.copy()
                F_rate_delta[np.logical_or(sigma_tmp <= (rho_grid[irho]-(0.1/2)),sigma_tmp > (rho_grid[irho]+(0.1/2)))] = np.nan
            elif irho == 1:
                F_rate_tmp[np.logical_or(sigma_tmp <= (rho_grid[irho]-(0.1/2)),sigma_tmp > (rho_grid[irho]+(0.1/2)))] = np.nan
                F_rate_delta = np.stack((F_rate_delta,F_rate_tmp), axis=0)
            else:
                F_rate_tmp[np.logical_or(sigma_tmp <= (rho_grid[irho]-(0.1/2)),sigma_tmp > (rho_grid[irho]+(0.1/2)))] = np.nan
                F_rate_extradim = np.expand_dims(F_rate_tmp, axis=0)
                F_rate_delta = np.concatenate((F_rate_delta,F_rate_extradim), axis=0)
            del F_rate_tmp
        
        F_rate_delta_var_vint_mm[mm,:] = F_rate_delta
    
    print('completed, size: ', F_rate_delta_var_vint_mm.shape)
    
    return F_rate_delta_var_vint_mm

# -- Units --
# Heat: m/s -> m3/s
# Salt: m/s -> m3/s
# Fwf: Kg.m-2.s-1 = Kg/s

# Whole domain:
Fs_rate_delta_adv_vint_mm = wmt(salt_adv_mld_vint, dx, dy,'budget')
Fs_rate_delta_diff_vint_mm = wmt(R_s_vint, dx, dy,'budget')
Fs_rate_delta_net_vint_mm = wmt(salt_net_mld_vint, dx, dy,'budget')
Fs_rate_delta_sfc_mm = wmt(salt_sfc, dx, dy,'sfc_frc')
    
Fh_rate_delta_adv_vint_mm = wmt(temp_adv_mld_vint, dx, dy,'budget')
Fh_rate_delta_diff_vint_mm = wmt(R_t_vint, dx, dy,'budget')
Fh_rate_delta_net_vint_mm = wmt(temp_net_mld_vint, dx, dy,'budget')
Fh_rate_delta_sfc_mm = wmt(temp_sfc, dx, dy,'sfc_frc')

Ff_rate_delta_adv_vint_mm = wmt(fwf_adv_mld_vint, dx, dy,'budget')
Ff_rate_delta_diff_vint_mm = wmt(R_f_vint, dx, dy,'budget')
Ff_rate_delta_net_vint_mm = wmt(fwf_net_mld_vint, dx, dy,'budget')
Ff_rate_delta_sfc_mm = wmt(fwf_sfc, dx, dy,'sfc_frc')
Ffm_rate_delta_sfc_mm = wmt(fwf_sfc+m, dx, dy,'sfc_frc')
Fm_rate_delta_sfc_mm = wmt(m, dx, dy,'sfc_frc')

# Excluding ice shelf:

Fs_rate_delta_adv_vint_noice_mm = wmt(salt_adv_mld_vint*mask_zice, dx, dy,'budget')
Fs_rate_delta_diff_vint_noice_mm = wmt(R_s_vint*mask_zice, dx, dy,'budget')
Fs_rate_delta_net_vint_noice_mm = wmt(salt_net_mld_vint*mask_zice, dx, dy,'budget')
Fs_rate_delta_sfc_noice_mm = wmt(salt_sfc*mask_zice, dx, dy,'sfc_frc')
    
Ff_rate_delta_adv_vint_noice_mm = wmt(fwf_adv_mld_vint*mask_zice, dx, dy,'budget')
Ff_rate_delta_diff_vint_noice_mm = wmt(R_f_vint*mask_zice, dx, dy,'budget')
Ff_rate_delta_net_vint_noice_mm = wmt(fwf_net_mld_vint*mask_zice, dx, dy,'budget')
Ff_rate_delta_sfc_noice_mm = wmt(fwf_sfc*mask_zice, dx, dy,'sfc_frc')

Fh_rate_delta_adv_vint_noice_mm = wmt(temp_adv_mld_vint*mask_zice, dx, dy,'budget')
Fh_rate_delta_diff_vint_noice_mm = wmt(R_t_vint*mask_zice, dx, dy,'budget')
Fh_rate_delta_net_vint_noice_mm = wmt(temp_net_mld_vint*mask_zice, dx, dy,'budget')
Fh_rate_delta_sfc_noice_mm = wmt(temp_sfc*mask_zice, dx, dy,'sfc_frc')

# development mask for density classes:
sigma_sept = sigma_t[8,0:50,0:50]

for irho in np.arange(0,len_rho_grid):
    print(irho, rho_grid[irho]-(0.1/2), rho_grid[irho]+(0.1/2), sigma_sept.shape)    
    icritst_tmp = ma.masked_where(np.logical_or(sigma_sept <= (rho_grid[irho]-(0.1/2)),sigma_sept > (rho_grid[irho]+(0.1/2))), sigma_sept)    
    
    if irho == 0:
        icritst = icritst_tmp.copy()
        icritst[np.logical_or(sigma_sept <= (rho_grid[irho]-(0.1/2)),sigma_sept > (rho_grid[irho]+(0.1/2)))] = np.nan
    elif irho == 1:
        icritst_tmp[np.logical_or(sigma_sept <= (rho_grid[irho]-(0.1/2)),sigma_sept > (rho_grid[irho]+(0.1/2)))] = np.nan
        icritst = np.stack((icritst,icritst_tmp), axis=0)
    else:
        icritst_tmp[np.logical_or(sigma_sept <= (rho_grid[irho]-(0.1/2)),sigma_sept > (rho_grid[irho]+(0.1/2)))] = np.nan
        icritst_extradim = np.expand_dims(icritst_tmp, axis=0)
        icritst = np.concatenate((icritst,icritst_extradim), axis=0)
    del icritst_tmp
    
print(icritst.shape)

# integrated over x, y directions: full domain

Fs_rate_sfc_mm_int = np.empty((len(rho_grid),12))
Ff_rate_sfc_mm_int = np.empty((len(rho_grid),12))
Fm_rate_sfc_mm_int = np.empty((len(rho_grid),12))
Ffm_rate_sfc_mm_int = np.empty((len(rho_grid),12))
Fh_rate_sfc_mm_int = np.empty((len(rho_grid),12))

Fs_rate_adv_vint_mm_int = np.empty((len(rho_grid),12))
Ff_rate_adv_vint_mm_int = np.empty((len(rho_grid),12))
Fh_rate_adv_vint_mm_int = np.empty((len(rho_grid),12))
Fs_rate_diff_vint_mm_int = np.empty((len(rho_grid),12))
Ff_rate_diff_vint_mm_int = np.empty((len(rho_grid),12))
Fh_rate_diff_vint_mm_int = np.empty((len(rho_grid),12))
Fs_rate_net_vint_mm_int = np.empty((len(rho_grid),12))
Ff_rate_net_vint_mm_int = np.empty((len(rho_grid),12))
Fh_rate_net_vint_mm_int = np.empty((len(rho_grid),12))


for irho in np.arange(0,len(rho_grid)):   
    for mm in np.arange(0,12):
        
        Fs_rate_sfc_mm_int[irho,mm] = np.nansum(np.nansum(Fs_rate_delta_sfc_mm[mm,irho,:], axis=1), axis=0)
        Ff_rate_sfc_mm_int[irho,mm] = np.nansum(np.nansum(Ff_rate_delta_sfc_mm[mm,irho,:], axis=1), axis=0)
        Fm_rate_sfc_mm_int[irho,mm] = np.nansum(np.nansum(Fm_rate_delta_sfc_mm[mm,irho,:], axis=1), axis=0)
        Ffm_rate_sfc_mm_int[irho,mm] = np.nansum(np.nansum(Ffm_rate_delta_sfc_mm[mm,irho,:], axis=1), axis=0)
        Fh_rate_sfc_mm_int[irho,mm] = np.nansum(np.nansum(Fh_rate_delta_sfc_mm[mm,irho,:], axis=1), axis=0)
    
        
        Fs_rate_adv_vint_mm_int[irho,mm] = np.nansum(np.nansum(Fs_rate_delta_adv_vint_mm[mm,irho,:], axis=1), axis=0)
        Ff_rate_adv_vint_mm_int[irho,mm] = np.nansum(np.nansum(Ff_rate_delta_adv_vint_mm[mm,irho,:], axis=1), axis=0)
        Fh_rate_adv_vint_mm_int[irho,mm] = np.nansum(np.nansum(Fh_rate_delta_adv_vint_mm[mm,irho,:], axis=1), axis=0)
        
        Fs_rate_diff_vint_mm_int[irho,mm] = np.nansum(np.nansum(Fs_rate_delta_diff_vint_mm[mm,irho,:], axis=1), axis=0)
        Ff_rate_diff_vint_mm_int[irho,mm] = np.nansum(np.nansum(Ff_rate_delta_diff_vint_mm[mm,irho,:], axis=1), axis=0)
        Fh_rate_diff_vint_mm_int[irho,mm] = np.nansum(np.nansum(Fh_rate_delta_diff_vint_mm[mm,irho,:], axis=1), axis=0)
        
        Fs_rate_net_vint_mm_int[irho,mm] = np.nansum(np.nansum(Fs_rate_delta_net_vint_mm[mm,irho,:], axis=1), axis=0)
        Ff_rate_net_vint_mm_int[irho,mm] = np.nansum(np.nansum(Ff_rate_delta_net_vint_mm[mm,irho,:], axis=1), axis=0)
        Fh_rate_net_vint_mm_int[irho,mm] = np.nansum(np.nansum(Fh_rate_delta_net_vint_mm[mm,irho,:], axis=1), axis=0)

# integrated over x, y directions: no ice shelf
Fs_rate_sfc_noice_mm_int = np.empty((len(rho_grid),12))
Ff_rate_sfc_noice_mm_int = np.empty((len(rho_grid),12))
Fh_rate_sfc_noice_mm_int = np.empty((len(rho_grid),12))

Fs_rate_adv_vint_noice_mm_int = np.empty((len(rho_grid),12))
Ff_rate_adv_vint_noice_mm_int = np.empty((len(rho_grid),12))
Fh_rate_adv_vint_noice_mm_int = np.empty((len(rho_grid),12))
Fs_rate_diff_vint_noice_mm_int = np.empty((len(rho_grid),12))
Ff_rate_diff_vint_noice_mm_int = np.empty((len(rho_grid),12))
Fh_rate_diff_vint_noice_mm_int = np.empty((len(rho_grid),12))
Fs_rate_net_vint_noice_mm_int = np.empty((len(rho_grid),12))
Ff_rate_net_vint_noice_mm_int = np.empty((len(rho_grid),12))
Fh_rate_net_vint_noice_mm_int = np.empty((len(rho_grid),12))


for irho in np.arange(0,len(rho_grid)):   
    for mm in np.arange(0,12):
        
        Fs_rate_sfc_noice_mm_int[irho,mm] = np.nansum(np.nansum(Fs_rate_delta_sfc_noice_mm[mm,irho,:], axis=1), axis=0)
        Ff_rate_sfc_noice_mm_int[irho,mm] = np.nansum(np.nansum(Ff_rate_delta_sfc_noice_mm[mm,irho,:], axis=1), axis=0)
        Fh_rate_sfc_noice_mm_int[irho,mm] = np.nansum(np.nansum(Fh_rate_delta_sfc_noice_mm[mm,irho,:], axis=1), axis=0)
        
        Fs_rate_adv_vint_noice_mm_int[irho,mm] = np.nansum(np.nansum(Fs_rate_delta_adv_vint_noice_mm[mm,irho,:], axis=1), axis=0)
        Ff_rate_adv_vint_noice_mm_int[irho,mm] = np.nansum(np.nansum(Ff_rate_delta_adv_vint_noice_mm[mm,irho,:], axis=1), axis=0)
        Fh_rate_adv_vint_noice_mm_int[irho,mm] = np.nansum(np.nansum(Fh_rate_delta_adv_vint_noice_mm[mm,irho,:], axis=1), axis=0)
        
        Fs_rate_diff_vint_noice_mm_int[irho,mm] = np.nansum(np.nansum(Fs_rate_delta_diff_vint_noice_mm[mm,irho,:], axis=1), axis=0)
        Ff_rate_diff_vint_noice_mm_int[irho,mm] = np.nansum(np.nansum(Ff_rate_delta_diff_vint_noice_mm[mm,irho,:], axis=1), axis=0)
        Fh_rate_diff_vint_noice_mm_int[irho,mm] = np.nansum(np.nansum(Fh_rate_delta_diff_vint_noice_mm[mm,irho,:], axis=1), axis=0)
        
        Fs_rate_net_vint_noice_mm_int[irho,mm] = np.nansum(np.nansum(Fs_rate_delta_net_vint_noice_mm[mm,irho,:], axis=1), axis=0)
        Ff_rate_net_vint_noice_mm_int[irho,mm] = np.nansum(np.nansum(Ff_rate_delta_net_vint_noice_mm[mm,irho,:], axis=1), axis=0)
        Fh_rate_net_vint_noice_mm_int[irho,mm] = np.nansum(np.nansum(Fh_rate_delta_net_vint_noice_mm[mm,irho,:], axis=1), axis=0)


# figures
fig_path = '/users/boeiradi/COLD_project/postprocessing/figs/WMT/'

# plot with bars
width=.023

# con/ert to rate per year: not doing it, keep seconds to convert to Sv.
Dt = rho0

# divide for Dt (sec in 1 year) and 10^6 to obtain rate in Sv year-1 

# SFC: whole domain
Fs_sig_sfc_mm = -Fs_rate_sfc_mm_int*rho0/1e6
Fm_sig_sfc_mm = -Fm_rate_sfc_mm_int*rho0/1e6
Fh_sig_sfc_mm = -Fh_rate_sfc_mm_int*rho0/1e6 
F_sig_sfc_mm = -Fs_sig_sfc_mm + Fh_sig_sfc_mm + Fm_sig_sfc_mm
# - calculate the anual average of the monthly ars:
Fs_sig_sfc =  np.nanmean(Fs_sig_sfc_mm, axis=1)
Fm_sig_sfc =  np.nanmean(Fm_sig_sfc_mm, axis=1)
Fh_sig_sfc =  np.nanmean(Fh_sig_sfc_mm, axis=1)
F_sig_sfc = -Fs_sig_sfc + Fh_sig_sfc + Fm_sig_sfc

# SFC: no ice shelf
Fs_sig_sfc_noice_mm = -Fs_rate_sfc_noice_mm_int*rho0/1e6
Fh_sig_sfc_noice_mm = -Fh_rate_sfc_noice_mm_int*rho0/1e6 
F_sig_sfc_noice_mm = -Fs_sig_sfc_noice_mm + Fh_sig_sfc_noice_mm
# - calculate the anual average of the monthly ars:
Fs_sig_sfc_noice =  np.nanmean(Fs_sig_sfc_noice_mm, axis=1)
Fh_sig_sfc_noice =  np.nanmean(Fh_sig_sfc_noice_mm, axis=1)
F_sig_sfc_noice = -Fs_sig_sfc_noice + Fh_sig_sfc_noice

# plots:


fig = plt.figure(figsize=(12,10))

ax1 = fig.add_subplot(221)
plt.title('Full domain')
cs=plt.plot(rho_grid,-Fs_sig_sfc,'b',label='WMT due to sfc salt flux')
cf=plt.plot(rho_grid,Fm_sig_sfc,'-g',label='WMT due to sfc ice melt')
ch=plt.plot(rho_grid,Fh_sig_sfc,'r',label='WMT due to sfc heat flux')
ct=plt.plot(rho_grid,F_sig_sfc,'k',label='WMT due to sfc total flux')
plt.legend()
plt.ylabel('Transformation rate sfc (Sv)')
plt.xlim(26.5,28),plt.ylim(-16,2.5)
plt.grid(True)
plt.text(26.6,1.5,'Buoyancy loss')
plt.text(26.6,-1.5,'Buoyancy gain')

ax2 = fig.add_subplot(222)
plt.title('Full domain minus ice shelfs')
cs=plt.plot(rho_grid,-Fs_sig_sfc_noice,'b',label='WMT due to sfc salt flux')
#cf=plt.plot(rho_grid,Ff_sig_sfc_noice,'--b',label='WMT due to sfc fw flux')
ch=plt.plot(rho_grid,Fh_sig_sfc_noice,'r',label='WMT due to sfc heat flux')
ct=plt.plot(rho_grid,F_sig_sfc_noice,'k',label='WMT due to sfc total flux')
plt.legend()
#plt.ylabel('Transformation rate sfc (Sv)')
plt.xlim(26.5,28),plt.ylim(-16,2.5)
plt.grid(True)

name_fig="waom4extend_shflim_S_0.25Q_WMT_sfc_fluxes_annual_iceshelf_comparison.png"
plt.savefig(fig_path + name_fig, dpi=300)

# ADV (vint)
Fs_sig_adv_vint_mm = -Fs_rate_adv_vint_mm_int*rho0/1e6 
Fh_sig_adv_vint_mm = -Fh_rate_adv_vint_mm_int*rho0/1e6 
F_sig_adv_vint_mm = -Fs_sig_adv_vint_mm + Fh_sig_adv_vint_mm
# - calculate the anual average of the monthly ars:
Fs_sig_adv_vint =  np.nanmean(Fs_sig_adv_vint_mm, axis=1)
Fh_sig_adv_vint =  np.nanmean(Fh_sig_adv_vint_mm, axis=1)
F_sig_adv_vint = -Fs_sig_adv_vint + Fh_sig_adv_vint 

# DIFF
Fs_sig_diff_vint_mm = -Fs_rate_diff_vint_mm_int*rho0/1e6 
Fh_sig_diff_vint_mm = -Fh_rate_diff_vint_mm_int*rho0/1e6 
F_sig_diff_vint_mm = Fs_sig_diff_vint_mm + Fh_sig_diff_vint_mm
# - calculate the anual average of the monthly ars:
Fs_sig_diff_vint =  np.nanmean(Fs_sig_diff_vint_mm, axis=1)
Fh_sig_diff_vint =  np.nanmean(Fh_sig_diff_vint_mm, axis=1)
F_sig_diff_vint = -Fs_sig_diff_vint + Fh_sig_diff_vint

# NET
Fs_sig_net_vint_mm = -Fs_rate_net_vint_mm_int*rho0/1e6 
Fh_sig_net_vint_mm = -Fh_rate_net_vint_mm_int*rho0/1e6 
F_sig_net_vint_mm = -Fs_sig_net_vint_mm + Fh_sig_net_vint_mm
# - calculate the anual average of the monthly ars:
Fs_sig_net_vint =  np.nanmean(Fs_sig_net_vint_mm, axis=1)
Fh_sig_net_vint =  np.nanmean(Fh_sig_net_vint_mm, axis=1)
F_sig_net_vint = -Fs_sig_net_vint + Fh_sig_net_vint

# obs (16.6.21): salt budget has opposite effect to buoyancy than fw budget,
#  so term are inverted (x-1) to plot the transformation rates. This seems to agree
#  better with the surface fluxes.
# obs (21.6.21): only inverted sign is not correct. Need to do full conversion as done from 
#  salt flux to salt flux, i.e. multiply by -1000 (fw density) and divide by salinity (avg in the MLD)

fig = plt.figure(figsize=(16,10))
ax1 = fig.add_subplot(222)
plt.title('Salt budget')
cs=plt.plot(rho_grid,-Fs_sig_net_vint,'b',label='WMT due to salt tend')
cs=plt.plot(rho_grid,-Fs_sig_diff_vint,'--b',label='WMT due to salt diff')
cs=plt.plot(rho_grid,-Fs_sig_adv_vint,':b',label='WMT due to salt adv')
plt.legend()
plt.ylabel('Transformation rate (Sv)')
plt.xlim(26.5,28),plt.ylim(-16,2.5)
plt.grid(True)
plt.text(26.6,1.5,'Buoyancy loss')
plt.text(26.6,-1.5,'Buoyancy gain')
ax1 = fig.add_subplot(223)
plt.title('Heat budget')
ch=plt.plot(rho_grid,Fh_sig_net_vint,'r',label='WMT due to heat tend')
ch=plt.plot(rho_grid,Fh_sig_diff_vint,'--r',label='WMT due to heat diff')
ch=plt.plot(rho_grid,Fh_sig_adv_vint,':r',label='WMT due to heat adv')
plt.legend()
plt.ylabel('Transformation rate (Sv)')
plt.xlim(26.5,28),plt.ylim(-16,2.5)
plt.grid(True)

ax2 = fig.add_subplot(221)
plt.title('Surface Fluxes budget')
cs=plt.plot(rho_grid,-Fs_sig_sfc,'b',label='WMT due to salt sfc')
cf=plt.plot(rho_grid,Fm_sig_sfc,'-g',label='WMT due to sfc ice melt')
ch=plt.plot(rho_grid,Fh_sig_sfc,'r',label='WMT due to heat sfc')
ct=plt.plot(rho_grid,F_sig_sfc,'k',label='WMT due to total sfc')
plt.legend()
plt.ylabel('Transformation rate (Sv)')
plt.xlim(26.5,28),plt.ylim(-16,2.5)
plt.grid(True)

ax2 = fig.add_subplot(224)
plt.title('Net tendencies')
cs=plt.plot(rho_grid,-Fs_sig_net_vint,'-b',label='WMT tendency due to salt sfc')
ch=plt.plot(rho_grid,Fh_sig_net_vint,'-r',label='WMT tendency due to heat sfc')
ct=plt.plot(rho_grid,Fh_sig_net_vint-Fs_sig_net_vint,'k',label='WMT tendency total')
plt.legend()
plt.ylabel('Transformation rate (Sv)')
plt.xlim(26.5,28),plt.ylim(-16,2.5)
plt.grid(True)

name_fig="waom4extend_shflim_S_0.25Q_WMT_heat-salt_vint_annual_fulldomain.png"
plt.savefig(fig_path + name_fig, dpi=300)

# ADV (vint)
Fs_sig_adv_vint_noice_mm = -Fs_rate_adv_vint_noice_mm_int*rho0/1e6 
Fh_sig_adv_vint_noice_mm = -Fh_rate_adv_vint_noice_mm_int*rho0/1e6 
F_sig_adv_vint_noice_mm = -Fs_sig_adv_vint_noice_mm + Fh_sig_adv_vint_noice_mm
# - calculate the anual average of the monthly ars:
Fs_sig_adv_vint_noice =  np.nanmean(Fs_sig_adv_vint_noice_mm, axis=1)
Fh_sig_adv_vint_noice =  np.nanmean(Fh_sig_adv_vint_noice_mm, axis=1)
F_sig_adv_vint_noice = -Fs_sig_adv_vint_noice + Fh_sig_adv_vint_noice

# DIFF
Fs_sig_diff_vint_noice_mm = -Fs_rate_diff_vint_noice_mm_int*rho0/1e6 
Fh_sig_diff_vint_noice_mm = -Fh_rate_diff_vint_noice_mm_int*rho0/1e6 
F_sig_diff_vint_noice_mm = -Fs_sig_diff_vint_noice_mm + Fh_sig_diff_vint_noice_mm
# - calculate the anual average of the monthly ars:
Fs_sig_diff_vint_noice =  np.nanmean(Fs_sig_diff_vint_noice_mm, axis=1)
Fh_sig_diff_vint_noice =  np.nanmean(Fh_sig_diff_vint_noice_mm, axis=1)
F_sig_diff_vin_noicet = -Fs_sig_diff_vint_noice + Fh_sig_diff_vint_noice

# NET
Fs_sig_net_vint_noice_mm = -Fs_rate_net_vint_noice_mm_int*rho0/1e6 
Fh_sig_net_vint_noice_mm = -Fh_rate_net_vint_noice_mm_int*rho0/1e6 
F_sig_net_vint_noice_mm = -Fs_sig_net_vint_noice_mm + Fh_sig_net_vint_noice_mm
# - calculate the anual average of the monthly ars:
Fs_sig_net_vint_noice =  np.nanmean(Fs_sig_net_vint_noice_mm, axis=1)
Fh_sig_net_vint_noice =  np.nanmean(Fh_sig_net_vint_noice_mm, axis=1)
F_sig_net_vint_noice = -Fs_sig_net_vint_noice + Fh_sig_net_vint_noice

# obs (16.6.21): salt budget has opposite effect to buoyancy than fw budget,
#  so term are inverted (x-1) to plot the transformation rates. This seems to agree
#  better with the surface fluxes.
# obs (21.6.21): only inverted sign is not correct. Need to do full conversion as done from 
#  salt flux to fw flux, i.e. multiply by -1000 (fw density) and divide by salinity (avg in the MLD)

fig = plt.figure(figsize=(16,10))
ax1 = fig.add_subplot(222)
plt.title('Salt budget')
cs=plt.plot(rho_grid,-Fs_sig_net_vint_noice,'b',label='WMT due to salt tend')
cs=plt.plot(rho_grid,-Fs_sig_diff_vint_noice,'--b',label='WMT due to salt diff')
cs=plt.plot(rho_grid,-Fs_sig_adv_vint_noice,':b',label='WMT due to salt adv')
plt.legend()
plt.ylabel('Transformation rate (Sv)')
plt.xlim(26.5,28),plt.ylim(-2.5,2.5)
plt.grid(True)
#plt.text(26.6,25,'Buoyancy loss')
#plt.text(26.6,-25,'Buoyancy gain')
ax1 = fig.add_subplot(223)
plt.title('Heat budget')
ch=plt.plot(rho_grid,Fh_sig_net_vint_noice,'r',label='WMT due to heat tend')
ch=plt.plot(rho_grid,Fh_sig_diff_vint_noice,'--r',label='WMT due to heat diff')
ch=plt.plot(rho_grid,Fh_sig_adv_vint_noice,':r',label='WMT due to heat adv')
plt.legend()
plt.ylabel('Transformation rate (Sv)')
plt.xlim(26.5,28),plt.ylim(-2.5,2.5)
plt.grid(True)
#plt.text(26.6,25,'Buoyancy loss')
#plt.text(26.6,-25,'Buoyancy gain')
ax2 = fig.add_subplot(221)
plt.title('Surface fluxes')
cs=plt.plot(rho_grid,-Fs_sig_sfc_noice,'b',label='WMT due to salt sfc')
ch=plt.plot(rho_grid,Fh_sig_sfc_noice,'r',label='WMT due to heat sfc')
ct=plt.plot(rho_grid,F_sig_sfc_noice,'k',label='WMT due to total sfc')
plt.legend()
plt.ylabel('Transformation rate (Sv)')
plt.xlim(26.5,28),plt.ylim(-2.5,2.5)
plt.grid(True)

ax2 = fig.add_subplot(224)
plt.title('Net tendencies')
cs=plt.plot(rho_grid,-Fs_sig_net_vint_noice,'-b',label='WMT tendency due to salt sfc')
ch=plt.plot(rho_grid,Fh_sig_net_vint_noice,'-r',label='WMT tendency due to heat sfc')
ct=plt.plot(rho_grid,Fh_sig_net_vint_noice-Fs_sig_net_vint_noice,'k',label='WMT tendency total')
plt.legend()
plt.ylabel('Transformation rate (Sv)')
plt.xlim(26.5,28),plt.ylim(-2.5,2.5)
plt.grid(True)
plt.title('Total tendencies')

name_fig="waom4extend_shflim_S_0.25Q_WMT_heat-salt_vint_annual_noice.png"
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
#Fs_rate_mm_int[irho,mm] = np.nansum(np.nansum(Fs_rate_delta_mm[mm,irho,:], axis=1), axis=0)*rho0/0.1
#Fh_rate_mm_int[irho,mm] = np.nansum(np.nansum(Fh_rate_delta_mm[mm,irho,:], axis=1), axis=0)*rho0/0.1
print(Fs_rate_delta_net_vint_mm.shape)

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
    cy=plt.pcolormesh(lon_rho,lat_rho,np.nanmean(Fs_rate_delta_net_vint_mm[:,irho,:], axis=0)*rho0/(dx*dy), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm)
    plt.colorbar(cy)
    ax2.gridlines() # draw_labels=True,linewidth=2, color='white', alpha=0.5, linestyle='--')
    lonlat_labels(ax2)
    ax2.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
    ax2.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white') 
    plt.clim(-5e-6,5e-6)

    ax3 = fig.add_subplot(223, projection=proj)
    plt.title('HF, $\sigma_{\Theta}$ = ' + str(rho_grid[irho]))
    cy=plt.pcolormesh(lon_rho,lat_rho,np.nanmean(Fh_rate_delta_net_vint_mm[:,irho,:], axis=0)*rho0/(dx*dy), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm)
    plt.colorbar(cy)
    ax3.gridlines() # draw_labels=True,linewidth=2, color='white', alpha=0.5, linestyle='--')
    lonlat_labels(ax3)
    ax3.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
    ax3.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white') 
    plt.clim(-5e-6,5e-6)
                                                
    ax4 = fig.add_subplot(224, projection=proj)
    plt.title('Total, $\sigma_{\Theta}$ = ' + str(rho_grid[irho]))
    cy=plt.pcolormesh(lon_rho,lat_rho,np.nanmean(Fh_rate_delta_net_vint_mm[:,irho,:], axis=0)*rho0/(dx*dy)+np.nanmean(Fs_rate_delta_net_vint_mm[:,irho,:], axis=0)*rho0/(dx*dy), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm)
    plt.colorbar(cy)
    ax4.gridlines() # draw_labels=True,linewidth=2, color='white', alpha=0.5, linestyle='--')
    lonlat_labels(ax4)
    ax4.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
    ax4.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white') 
    plt.clim(-5e-6,5e-6)
                                                
    name_fig="waom4extend_shflim_S_0.25Q_WMTmaps_annual_yr10_l" + str(irho)  + "_fulldomain.png"
    plt.savefig(fig_path + name_fig, dpi=300)
    plt.close()

# plot maps for 27.5 kg.m-3 isopycnal
#Fs_rate_mm_int[irho,mm] = np.nansum(np.nansum(Fs_rate_delta_mm[mm,irho,:], axis=1), axis=0)*rho0/0.1
#Fh_rate_mm_int[irho,mm] = np.nansum(np.nansum(Fh_rate_delta_mm[mm,irho,:], axis=1), axis=0)*rho0/0.1
print(Fs_rate_delta_net_vint_noice_mm.shape)

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
    cy=plt.pcolormesh(lon_rho,lat_rho,np.nanmean(Fs_rate_delta_net_vint_noice_mm[:,irho,:], axis=0)*rho0/(dx*dy), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm)
    plt.colorbar(cy)
    ax2.gridlines() # draw_labels=True,linewidth=2, color='white', alpha=0.5, linestyle='--')
    lonlat_labels(ax2)
    ax2.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
    ax2.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white') 
    plt.clim(-5e-6,5e-6)

    ax3 = fig.add_subplot(223, projection=proj)
    plt.title('HF, $\sigma_{\Theta}$ = ' + str(rho_grid[irho]))
    cy=plt.pcolormesh(lon_rho,lat_rho,np.nanmean(Fh_rate_delta_net_vint_noice_mm[:,irho,:], axis=0)*rho0/(dx*dy), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm)
    plt.colorbar(cy)
    ax3.gridlines() # draw_labels=True,linewidth=2, color='white', alpha=0.5, linestyle='--')
    lonlat_labels(ax3)
    ax3.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
    ax3.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white') 
    plt.clim(-5e-6,5e-6)
                                                
    ax4 = fig.add_subplot(224, projection=proj)
    plt.title('Total, $\sigma_{\Theta}$ = ' + str(rho_grid[irho]))
    cy=plt.pcolormesh(lon_rho,lat_rho,np.nanmean(Fh_rate_delta_net_vint_noice_mm[:,irho,:], axis=0)*rho0/(dx*dy)+np.nanmean(Fs_rate_delta_net_vint_noice_mm[:,irho,:], axis=0)*rho0/(dx*dy), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm)
    plt.colorbar(cy)
    ax4.gridlines() # draw_labels=True,linewidth=2, color='white', alpha=0.5, linestyle='--')
    lonlat_labels(ax4)
    ax4.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
    ax4.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white') 
    plt.clim(-5e-6,5e-6)
                                                
    name_fig="waom4extend_shflim_S_0.25Q_WMTmaps_annual_yr10_l" + str(irho)  + "_noice.png"
    plt.savefig(fig_path + name_fig, dpi=300)
    plt.close()

proj = ccrs.SouthPolarStereo()
fig = plt.figure(figsize=(16,10))

irho=30
ax1 = fig.add_subplot(221, projection=proj)
plt.title('Total, $\sigma_{\Theta}$ = ' + str(rho_grid[irho]))
cy=plt.pcolormesh(lon_rho,lat_rho,np.nanmean(Fh_rate_delta_net_vint_mm[:,irho,:], axis=0)*rho0/(dx*dy)+np.nanmean(Fs_rate_delta_net_vint_mm[:,irho,:], axis=0)*rho0/(dx*dy), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm)
plt.colorbar(cy, extend='both')
ax1.gridlines() # draw_labels=True,linewidth=2, color='white', alpha=0.5, linestyle='--')
lonlat_labels(ax1)
ax1.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
ax1.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white') 
plt.clim(-5e-6,5e-6)

irho=31
ax2 = fig.add_subplot(222, projection=proj)
plt.title('Total, $\sigma_{\Theta}$ = ' + str(rho_grid[irho]))
cy=plt.pcolormesh(lon_rho,lat_rho,np.nanmean(Fh_rate_delta_net_vint_mm[:,irho,:], axis=0)*rho0/(dx*dy)+np.nanmean(Fs_rate_delta_net_vint_mm[:,irho,:], axis=0)*rho0/(dx*dy), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm)
plt.colorbar(cy, extend='both')
ax2.gridlines() # draw_labels=True,linewidth=2, color='white', alpha=0.5, linestyle='--')
lonlat_labels(ax2)
ax2.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
ax2.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white') 
plt.clim(-5e-6,5e-6)

irho=32
ax3 = fig.add_subplot(223, projection=proj)
plt.title('Total, $\sigma_{\Theta}$ = ' + str(rho_grid[irho]))
cy=plt.pcolormesh(lon_rho,lat_rho,np.nanmean(Fh_rate_delta_net_vint_mm[:,irho,:], axis=0)*rho0/(dx*dy)+np.nanmean(Fs_rate_delta_net_vint_mm[:,irho,:], axis=0)*rho0/(dx*dy), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm)
plt.colorbar(cy, extend='both')
ax3.gridlines() # draw_labels=True,linewidth=2, color='white', alpha=0.5, linestyle='--')
lonlat_labels(ax3)
ax3.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
ax3.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white') 
plt.clim(-5e-6,5e-6)

irho=33
ax4 = fig.add_subplot(224, projection=proj)
plt.title('Total, $\sigma_{\Theta}$ = ' + str(rho_grid[irho]))
cy=plt.pcolormesh(lon_rho,lat_rho,np.nanmean(Fh_rate_delta_net_vint_mm[:,irho,:], axis=0)*rho0/(dx*dy)+np.nanmean(Fs_rate_delta_net_vint_mm[:,irho,:], axis=0)*rho0/(dx*dy), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm)
plt.colorbar(cy, extend='both')
ax4.gridlines() # draw_labels=True,linewidth=2, color='white', alpha=0.5, linestyle='--')
lonlat_labels(ax4)
ax4.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
ax4.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white') 
plt.clim(-5e-6,5e-6)
                                                
name_fig="waom4extend_shflim_S_0.25Q_WMTmaps_annual_yr10_fulldomain.png"
plt.savefig(fig_path + name_fig, dpi=300)
plt.close()

proj = ccrs.SouthPolarStereo()
fig = plt.figure(figsize=(16,10))

irho=30
ax1 = fig.add_subplot(221, projection=proj)
plt.title('Total, $\sigma_{\Theta}$ = ' + str(rho_grid[irho]))
cy=plt.pcolormesh(lon_rho,lat_rho,np.nanmean(Fh_rate_delta_net_vint_noice_mm[:,irho,:], axis=0)*rho0/(dx*dy)+np.nanmean(Fs_rate_delta_net_vint_noice_mm[:,irho,:], axis=0)*rho0/(dx*dy), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm)
plt.colorbar(cy, extend='both')
ax1.gridlines() # draw_labels=True,linewidth=2, color='white', alpha=0.5, linestyle='--')
lonlat_labels(ax1)
ax1.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
ax1.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white') 
plt.clim(-5e-6,5e-6)

irho=31
ax2 = fig.add_subplot(222, projection=proj)
plt.title('Total, $\sigma_{\Theta}$ = ' + str(rho_grid[irho]))
cy=plt.pcolormesh(lon_rho,lat_rho,np.nanmean(Fh_rate_delta_net_vint_noice_mm[:,irho,:], axis=0)*rho0/(dx*dy)+np.nanmean(Fs_rate_delta_net_vint_noice_mm[:,irho,:], axis=0)*rho0/(dx*dy), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm)
plt.colorbar(cy, extend='both')
ax2.gridlines() # draw_labels=True,linewidth=2, color='white', alpha=0.5, linestyle='--')
lonlat_labels(ax2)
ax2.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
ax2.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white') 
plt.clim(-5e-6,5e-6)

irho=32
ax3 = fig.add_subplot(223, projection=proj)
plt.title('Total, $\sigma_{\Theta}$ = ' + str(rho_grid[irho]))
cy=plt.pcolormesh(lon_rho,lat_rho,np.nanmean(Fh_rate_delta_net_vint_noice_mm[:,irho,:], axis=0)*rho0/(dx*dy)+np.nanmean(Fs_rate_delta_net_vint_noice_mm[:,irho,:], axis=0)*rho0/(dx*dy), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm)
plt.colorbar(cy, extend='both')
ax3.gridlines() # draw_labels=True,linewidth=2, color='white', alpha=0.5, linestyle='--')
lonlat_labels(ax3)
ax3.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
ax3.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white') 
plt.clim(-5e-6,5e-6)

irho=33
ax4 = fig.add_subplot(224, projection=proj)
plt.title('Total, $\sigma_{\Theta}$ = ' + str(rho_grid[irho]))
cy=plt.pcolormesh(lon_rho,lat_rho,np.nanmean(Fh_rate_delta_net_vint_noice_mm[:,irho,:], axis=0)*rho0/(dx*dy)+np.nanmean(Fs_rate_delta_net_vint_noice_mm[:,irho,:], axis=0)*rho0/(dx*dy), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm)
plt.colorbar(cy, extend='both')
ax4.gridlines() # draw_labels=True,linewidth=2, color='white', alpha=0.5, linestyle='--')
lonlat_labels(ax4)
ax4.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
ax4.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white') 
plt.clim(-5e-6,5e-6)
                                                
name_fig="waom4extend_shflim_S_0.25Q_WMTmaps_annual_yr10_noice.png"
plt.savefig(fig_path + name_fig, dpi=300)
plt.close()


