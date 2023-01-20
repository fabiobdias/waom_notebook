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
    ds = xr.open_dataset('/g/data3/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom4extend_shflim_S_0.25Q/output_yr10_10km-bathy_diag/ocean_avg_00' + mm + '.nc')
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

di = xr.open_dataset('/g/data3/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom4extend_shflim_S_0.25Q/output_yr10_10km-bathy_diag/ocean_avg_0001.nc')
ice_draft = di.variables["zice"]

mask_zice = ma.masked_where(ice_draft < 0, np.ones(ice_draft.shape))
mask_outice = ma.masked_where(ice_draft >= 0, np.ones(ice_draft.shape))

di.close()

# calculate dz following:
dz = np.empty((12,1400,1575,31))
dz_inv = np.empty((12,1400,1575,31))


for tt in np.arange(0,12):
    z_w_sorted = -1*z_w[tt,:,:,::-1]
    print(z_w_sorted.shape)
    dz_inv[tt,:,:,:] = np.diff(z_w_sorted,axis=2)
    dz[tt,:,:,:] = dz_inv[tt,:,:,::-1]

dg = xr.open_dataset("/g/data3/hh5/tmp/access-om/fbd581/ROMS/waom4_frc/waom4extend_grd.nc")

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

dx = xr.open_dataset('/g/data3/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom4extend_shflim_S_0.25Q/output_yr10_10km-bathy_diag/Full_vint_vars_for_WMT_m.s-1.nc')

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

# ML budget equation: 
# salt:
# salt_sfc - R_s = -salt_tend_mld + salt_hdia_adv_mld + salt_vdia_adv_mld

# temp:
# temp_sfc - R_t = -temp_tend_mld + temp_hdia_adv_mld + temp_vdia_adv_mld

#  Function to calculate Water Mass Transformation (in m3/s):

# rho grid for binning:
#rho_grid=np.arange(35.5,37.4,0.1) # for sigma-2
#rho_grid=np.arange(24.4,29.1,0.1) # for sigma-0
rho_grid=np.arange(26.,28.,0.05) # for sigma-0

len_rho_grid=len(rho_grid)

dx = np.divide(1,pm)
dy = np.divide(1,pn)

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

# create mask per longitude: Weddell, Maud Land, East Antarctica, Ross, West Antarctica
import matplotlib.path as mpath
import cartopy.feature as cfeature

# 1) Mask Weddell (90W : 25W, lat < -75, 65W : 25W, lat <= -75)
mask_Wed1lon = ma.masked_where(lon_rho <= -90, np.ones(lon_rho.shape)) # West limit for lat<-75
mask_Wed1lat = ma.masked_where(lat_rho > -75, np.ones(lat_rho.shape))
mask_Wed1 = mask_Wed1lon*mask_Wed1lat

mask_Wed2lon = ma.masked_where(lon_rho <= -65, np.ones(lon_rho.shape)) # West limit for lat<-75
mask_Wed2lat = ma.masked_where(lat_rho <= -75, np.ones(lat_rho.shape))
mask_Wed2 = mask_Wed2lon*mask_Wed2lat

mask_Wed12 = np.ma.array(mask_Wed1.filled(1) * mask_Wed2.filled(1), mask=(mask_Wed1.mask * mask_Wed2.mask))

mask_Wed3= ma.masked_where(lon_rho > -25, np.ones(lon_rho.shape)) # East limit for any latitude

mask_Wed = mask_Wed12*mask_Wed3

# 2) Mask Maud Land (25W : 60E, lat < -60)

mask_Mau1lon = ma.masked_where(lon_rho <= -25, np.ones(lon_rho.shape)) # East limit for any latitude
mask_Mau1lat = ma.masked_where(lat_rho > -60, np.ones(lat_rho.shape))
mask_Mau1 = mask_Mau1lon*mask_Mau1lat

mask_Mau2lon = ma.masked_where(lon_rho >= 60, np.ones(lon_rho.shape)) # East limit for any latitude
mask_Mau2lat = ma.masked_where(lat_rho > -60, np.ones(lat_rho.shape))
mask_Mau2 = mask_Mau2lon*mask_Mau2lat

mask_Mau = mask_Mau1*mask_Mau2

# 3) Mask East Antarctica (60E : 160E, lat < -60)

mask_EAnt1lon = ma.masked_where(lon_rho < 60, np.ones(lon_rho.shape)) # East limit for any latitude
mask_EAnt1lat = ma.masked_where(lat_rho > -60, np.ones(lat_rho.shape))
mask_EAnt1 = mask_EAnt1lon*mask_EAnt1lat

mask_EAnt2lon = ma.masked_where(lon_rho >= 160, np.ones(lon_rho.shape)) # East limit for any latitude
mask_EAnt2lat = ma.masked_where(lat_rho > -60, np.ones(lat_rho.shape))
mask_EAnt2 = mask_EAnt2lon*mask_EAnt2lat

mask_EAnt = mask_EAnt1*mask_EAnt2

# 4) Mask Ross (140W : 160E, lat < -60)

mask_Ros1lon = ma.masked_where(lon_rho < 160, np.ones(lon_rho.shape)) # East limit for any latitude
mask_Ros1lat = ma.masked_where(lat_rho > -60, np.ones(lat_rho.shape))
mask_Ros1 = mask_Ros1lon*mask_Ros1lat

mask_Ros2lon = ma.masked_where(lon_rho >= -150, np.ones(lon_rho.shape)) # East limit for any latitude
mask_Ros2lat = ma.masked_where(lat_rho > -60, np.ones(lat_rho.shape))
mask_Ros2 = mask_Ros2lon*mask_Ros2lat

mask_Ros = np.ma.array(mask_Ros1.filled(1) * mask_Ros2.filled(1), mask=(mask_Ros1.mask * mask_Ros2.mask))

# 5) Mask West Antarctica (150W : 90W, lat < -65)

mask_WAnt1 = ma.masked_where(lon_rho < -150, np.ones(lon_rho.shape)) # West limit for any latitude

mask_WAnt2lon = ma.masked_where(lon_rho > -90, np.ones(lon_rho.shape)) # East limit for lat <-73
mask_WAnt2lat = ma.masked_where(lat_rho > -75, np.ones(lat_rho.shape))
mask_WAnt2 = mask_WAnt2lon*mask_WAnt2lat

mask_WAnt3lon = ma.masked_where(lon_rho > -65, np.ones(lon_rho.shape)) # East limit for lat >-73
mask_WAnt3lat = ma.masked_where(lat_rho <= -75, np.ones(lat_rho.shape))
mask_WAnt3 = mask_WAnt3lon*mask_WAnt3lat

mask_WAnt23 = np.ma.array(mask_WAnt2.filled(1) * mask_WAnt3.filled(1), mask=(mask_WAnt2.mask * mask_WAnt3.mask))

mask_WAnt = mask_WAnt1*mask_WAnt23

# Shelf only: Weddell

Fs_rate_delta_adv_vint_iceshelf_Wed_mm = wmt(salt_adv_full_vint*mask_outice*mask_Wed*mask_Wed, dx, dy,'budget')
Fs_rate_delta_diff_vint_iceshelf_Wed_mm = wmt(R_s_vint*mask_outice*mask_Wed, dx, dy,'budget')
Fs_rate_delta_net_vint_iceshelf_Wed_mm = wmt(salt_net_full_vint*mask_outice*mask_Wed, dx, dy,'budget')
Fs_rate_delta_sfc_iceshelf_Wed_mm = wmt(salt_sfc*mask_outice*mask_Wed, dx, dy,'sfc_frc')

Fh_rate_delta_adv_vint_iceshelf_Wed_mm = wmt(temp_adv_full_vint*mask_outice*mask_Wed, dx, dy,'budget')
Fh_rate_delta_diff_vint_iceshelf_Wed_mm = wmt(R_t_vint*mask_outice*mask_Wed, dx, dy,'budget')
Fh_rate_delta_net_vint_iceshelf_Wed_mm = wmt(temp_net_full_vint*mask_outice*mask_Wed, dx, dy,'budget')
Fh_rate_delta_sfc_iceshelf_Wed_mm = wmt(temp_sfc*mask_outice*mask_Wed, dx, dy,'sfc_frc')


# integrated over x, y directions: only continental shelf

Fs_rate_sfc_iceshelf_Wed_mm_int = np.empty((len(rho_grid),12))
Fh_rate_sfc_iceshelf_Wed_mm_int = np.empty((len(rho_grid),12))

Fs_rate_adv_vint_iceshelf_Wed_mm_int = np.empty((len(rho_grid),12))
Fh_rate_adv_vint_iceshelf_Wed_mm_int = np.empty((len(rho_grid),12))
Fs_rate_diff_vint_iceshelf_Wed_mm_int = np.empty((len(rho_grid),12))
Fh_rate_diff_vint_iceshelf_Wed_mm_int = np.empty((len(rho_grid),12))
Fs_rate_net_vint_iceshelf_Wed_mm_int = np.empty((len(rho_grid),12))
Fh_rate_net_vint_iceshelf_Wed_mm_int = np.empty((len(rho_grid),12))

for irho in np.arange(0,len(rho_grid)):
    for mm in np.arange(0,12):

        Fs_rate_sfc_iceshelf_Wed_mm_int[irho,mm] = np.nansum(np.nansum(Fs_rate_delta_sfc_iceshelf_Wed_mm[mm,irho,:], axis=1), axis=0)
        Fh_rate_sfc_iceshelf_Wed_mm_int[irho,mm] = np.nansum(np.nansum(Fh_rate_delta_sfc_iceshelf_Wed_mm[mm,irho,:], axis=1), axis=0)
        Fs_rate_adv_vint_iceshelf_Wed_mm_int[irho,mm] = np.nansum(np.nansum(Fs_rate_delta_adv_vint_iceshelf_Wed_mm[mm,irho,:], axis=1), axis=0)
        Fh_rate_adv_vint_iceshelf_Wed_mm_int[irho,mm] = np.nansum(np.nansum(Fh_rate_delta_adv_vint_iceshelf_Wed_mm[mm,irho,:], axis=1), axis=0)
        Fs_rate_diff_vint_iceshelf_Wed_mm_int[irho,mm] = np.nansum(np.nansum(Fs_rate_delta_diff_vint_iceshelf_Wed_mm[mm,irho,:], axis=1), axis=0)
        Fh_rate_diff_vint_iceshelf_Wed_mm_int[irho,mm] = np.nansum(np.nansum(Fh_rate_delta_diff_vint_iceshelf_Wed_mm[mm,irho,:], axis=1), axis=0)
        Fs_rate_net_vint_iceshelf_Wed_mm_int[irho,mm] = np.nansum(np.nansum(Fs_rate_delta_net_vint_iceshelf_Wed_mm[mm,irho,:], axis=1), axis=0)
        Fh_rate_net_vint_iceshelf_Wed_mm_int[irho,mm] = np.nansum(np.nansum(Fh_rate_delta_net_vint_iceshelf_Wed_mm[mm,irho,:], axis=1), axis=0)

# Shelf only: Maud Land

Fs_rate_delta_adv_vint_iceshelf_Mau_mm = wmt(salt_adv_full_vint*mask_outice*mask_Mau*mask_Mau, dx, dy,'budget')
Fs_rate_delta_diff_vint_iceshelf_Mau_mm = wmt(R_s_vint*mask_outice*mask_Mau, dx, dy,'budget')
Fs_rate_delta_net_vint_iceshelf_Mau_mm = wmt(salt_net_full_vint*mask_outice*mask_Mau, dx, dy,'budget')
Fs_rate_delta_sfc_iceshelf_Mau_mm = wmt(salt_sfc*mask_outice*mask_Mau, dx, dy,'sfc_frc')

Fh_rate_delta_adv_vint_iceshelf_Mau_mm = wmt(temp_adv_full_vint*mask_outice*mask_Mau, dx, dy,'budget')
Fh_rate_delta_diff_vint_iceshelf_Mau_mm = wmt(R_t_vint*mask_outice*mask_Mau, dx, dy,'budget')
Fh_rate_delta_net_vint_iceshelf_Mau_mm = wmt(temp_net_full_vint*mask_outice*mask_Mau, dx, dy,'budget')
Fh_rate_delta_sfc_iceshelf_Mau_mm = wmt(temp_sfc*mask_outice*mask_Mau, dx, dy,'sfc_frc')


# integrated over x, y directions: only continental shelf

Fs_rate_sfc_iceshelf_Mau_mm_int = np.empty((len(rho_grid),12))
Fh_rate_sfc_iceshelf_Mau_mm_int = np.empty((len(rho_grid),12))

Fs_rate_adv_vint_iceshelf_Mau_mm_int = np.empty((len(rho_grid),12))
Fh_rate_adv_vint_iceshelf_Mau_mm_int = np.empty((len(rho_grid),12))
Fs_rate_diff_vint_iceshelf_Mau_mm_int = np.empty((len(rho_grid),12))
Fh_rate_diff_vint_iceshelf_Mau_mm_int = np.empty((len(rho_grid),12))
Fs_rate_net_vint_iceshelf_Mau_mm_int = np.empty((len(rho_grid),12))
Fh_rate_net_vint_iceshelf_Mau_mm_int = np.empty((len(rho_grid),12))

for irho in np.arange(0,len(rho_grid)):
    for mm in np.arange(0,12):

        Fs_rate_sfc_iceshelf_Mau_mm_int[irho,mm] = np.nansum(np.nansum(Fs_rate_delta_sfc_iceshelf_Mau_mm[mm,irho,:], axis=1), axis=0)
        Fh_rate_sfc_iceshelf_Mau_mm_int[irho,mm] = np.nansum(np.nansum(Fh_rate_delta_sfc_iceshelf_Mau_mm[mm,irho,:], axis=1), axis=0)
        Fs_rate_adv_vint_iceshelf_Mau_mm_int[irho,mm] = np.nansum(np.nansum(Fs_rate_delta_adv_vint_iceshelf_Mau_mm[mm,irho,:], axis=1), axis=0)
        Fh_rate_adv_vint_iceshelf_Mau_mm_int[irho,mm] = np.nansum(np.nansum(Fh_rate_delta_adv_vint_iceshelf_Mau_mm[mm,irho,:], axis=1), axis=0)
        Fs_rate_diff_vint_iceshelf_Mau_mm_int[irho,mm] = np.nansum(np.nansum(Fs_rate_delta_diff_vint_iceshelf_Mau_mm[mm,irho,:], axis=1), axis=0)
        Fh_rate_diff_vint_iceshelf_Mau_mm_int[irho,mm] = np.nansum(np.nansum(Fh_rate_delta_diff_vint_iceshelf_Mau_mm[mm,irho,:], axis=1), axis=0)
        Fs_rate_net_vint_iceshelf_Mau_mm_int[irho,mm] = np.nansum(np.nansum(Fs_rate_delta_net_vint_iceshelf_Mau_mm[mm,irho,:], axis=1), axis=0)
        Fh_rate_net_vint_iceshelf_Mau_mm_int[irho,mm] = np.nansum(np.nansum(Fh_rate_delta_net_vint_iceshelf_Mau_mm[mm,irho,:], axis=1), axis=0)

# Shelf only: East Antarctica

Fs_rate_delta_adv_vint_iceshelf_EAnt_mm = wmt(salt_adv_full_vint*mask_outice*mask_EAnt*mask_EAnt, dx, dy,'budget')
Fs_rate_delta_diff_vint_iceshelf_EAnt_mm = wmt(R_s_vint*mask_outice*mask_EAnt, dx, dy,'budget')
Fs_rate_delta_net_vint_iceshelf_EAnt_mm = wmt(salt_net_full_vint*mask_outice*mask_EAnt, dx, dy,'budget')
Fs_rate_delta_sfc_iceshelf_EAnt_mm = wmt(salt_sfc*mask_outice*mask_EAnt, dx, dy,'sfc_frc')

Fh_rate_delta_adv_vint_iceshelf_EAnt_mm = wmt(temp_adv_full_vint*mask_outice*mask_EAnt, dx, dy,'budget')
Fh_rate_delta_diff_vint_iceshelf_EAnt_mm = wmt(R_t_vint*mask_outice*mask_EAnt, dx, dy,'budget')
Fh_rate_delta_net_vint_iceshelf_EAnt_mm = wmt(temp_net_full_vint*mask_outice*mask_EAnt, dx, dy,'budget')
Fh_rate_delta_sfc_iceshelf_EAnt_mm = wmt(temp_sfc*mask_outice*mask_EAnt, dx, dy,'sfc_frc')


# integrated over x, y directions: only continental shelf

Fs_rate_sfc_iceshelf_EAnt_mm_int = np.empty((len(rho_grid),12))
Fh_rate_sfc_iceshelf_EAnt_mm_int = np.empty((len(rho_grid),12))

Fs_rate_adv_vint_iceshelf_EAnt_mm_int = np.empty((len(rho_grid),12))
Fh_rate_adv_vint_iceshelf_EAnt_mm_int = np.empty((len(rho_grid),12))
Fs_rate_diff_vint_iceshelf_EAnt_mm_int = np.empty((len(rho_grid),12))
Fh_rate_diff_vint_iceshelf_EAnt_mm_int = np.empty((len(rho_grid),12))
Fs_rate_net_vint_iceshelf_EAnt_mm_int = np.empty((len(rho_grid),12))
Fh_rate_net_vint_iceshelf_EAnt_mm_int = np.empty((len(rho_grid),12))


for irho in np.arange(0,len(rho_grid)):
    for mm in np.arange(0,12):

        Fs_rate_sfc_iceshelf_EAnt_mm_int[irho,mm] = np.nansum(np.nansum(Fs_rate_delta_sfc_iceshelf_EAnt_mm[mm,irho,:], axis=1), axis=0)
        Fh_rate_sfc_iceshelf_EAnt_mm_int[irho,mm] = np.nansum(np.nansum(Fh_rate_delta_sfc_iceshelf_EAnt_mm[mm,irho,:], axis=1), axis=0)
        Fs_rate_adv_vint_iceshelf_EAnt_mm_int[irho,mm] = np.nansum(np.nansum(Fs_rate_delta_adv_vint_iceshelf_EAnt_mm[mm,irho,:], axis=1), axis=0)
        Fh_rate_adv_vint_iceshelf_EAnt_mm_int[irho,mm] = np.nansum(np.nansum(Fh_rate_delta_adv_vint_iceshelf_EAnt_mm[mm,irho,:], axis=1), axis=0)
        Fs_rate_diff_vint_iceshelf_EAnt_mm_int[irho,mm] = np.nansum(np.nansum(Fs_rate_delta_diff_vint_iceshelf_EAnt_mm[mm,irho,:], axis=1), axis=0)
        Fh_rate_diff_vint_iceshelf_EAnt_mm_int[irho,mm] = np.nansum(np.nansum(Fh_rate_delta_diff_vint_iceshelf_EAnt_mm[mm,irho,:], axis=1), axis=0)
        Fs_rate_net_vint_iceshelf_EAnt_mm_int[irho,mm] = np.nansum(np.nansum(Fs_rate_delta_net_vint_iceshelf_EAnt_mm[mm,irho,:], axis=1), axis=0)
        Fh_rate_net_vint_iceshelf_EAnt_mm_int[irho,mm] = np.nansum(np.nansum(Fh_rate_delta_net_vint_iceshelf_EAnt_mm[mm,irho,:], axis=1), axis=0)

# Shelf only: Ross

Fs_rate_delta_adv_vint_iceshelf_Ros_mm = wmt(salt_adv_full_vint*mask_outice*mask_Ros*mask_Ros, dx, dy,'budget')
Fs_rate_delta_diff_vint_iceshelf_Ros_mm = wmt(R_s_vint*mask_outice*mask_Ros, dx, dy,'budget')
Fs_rate_delta_net_vint_iceshelf_Ros_mm = wmt(salt_net_full_vint*mask_outice*mask_Ros, dx, dy,'budget')
Fs_rate_delta_sfc_iceshelf_Ros_mm = wmt(salt_sfc*mask_outice*mask_Ros, dx, dy,'sfc_frc')

Fh_rate_delta_adv_vint_iceshelf_Ros_mm = wmt(temp_adv_full_vint*mask_outice*mask_Ros, dx, dy,'budget')
Fh_rate_delta_diff_vint_iceshelf_Ros_mm = wmt(R_t_vint*mask_outice*mask_Ros, dx, dy,'budget')
Fh_rate_delta_net_vint_iceshelf_Ros_mm = wmt(temp_net_full_vint*mask_outice*mask_Ros, dx, dy,'budget')
Fh_rate_delta_sfc_iceshelf_Ros_mm = wmt(temp_sfc*mask_outice*mask_Ros, dx, dy,'sfc_frc')


# integrated over x, y directions: only continental shelf

Fs_rate_sfc_iceshelf_Ros_mm_int = np.empty((len(rho_grid),12))
Fh_rate_sfc_iceshelf_Ros_mm_int = np.empty((len(rho_grid),12))

Fs_rate_adv_vint_iceshelf_Ros_mm_int = np.empty((len(rho_grid),12))
Fh_rate_adv_vint_iceshelf_Ros_mm_int = np.empty((len(rho_grid),12))
Fs_rate_diff_vint_iceshelf_Ros_mm_int = np.empty((len(rho_grid),12))
Fh_rate_diff_vint_iceshelf_Ros_mm_int = np.empty((len(rho_grid),12))
Fs_rate_net_vint_iceshelf_Ros_mm_int = np.empty((len(rho_grid),12))
Fh_rate_net_vint_iceshelf_Ros_mm_int = np.empty((len(rho_grid),12))


for irho in np.arange(0,len(rho_grid)):
    for mm in np.arange(0,12):

        Fs_rate_sfc_iceshelf_Ros_mm_int[irho,mm] = np.nansum(np.nansum(Fs_rate_delta_sfc_iceshelf_Ros_mm[mm,irho,:], axis=1), axis=0)
        Fh_rate_sfc_iceshelf_Ros_mm_int[irho,mm] = np.nansum(np.nansum(Fh_rate_delta_sfc_iceshelf_Ros_mm[mm,irho,:], axis=1), axis=0)
        Fs_rate_adv_vint_iceshelf_Ros_mm_int[irho,mm] = np.nansum(np.nansum(Fs_rate_delta_adv_vint_iceshelf_Ros_mm[mm,irho,:], axis=1), axis=0)
        Fh_rate_adv_vint_iceshelf_Ros_mm_int[irho,mm] = np.nansum(np.nansum(Fh_rate_delta_adv_vint_iceshelf_Ros_mm[mm,irho,:], axis=1), axis=0)
        Fs_rate_diff_vint_iceshelf_Ros_mm_int[irho,mm] = np.nansum(np.nansum(Fs_rate_delta_diff_vint_iceshelf_Ros_mm[mm,irho,:], axis=1), axis=0)
        Fh_rate_diff_vint_iceshelf_Ros_mm_int[irho,mm] = np.nansum(np.nansum(Fh_rate_delta_diff_vint_iceshelf_Ros_mm[mm,irho,:], axis=1), axis=0)
        Fs_rate_net_vint_iceshelf_Ros_mm_int[irho,mm] = np.nansum(np.nansum(Fs_rate_delta_net_vint_iceshelf_Ros_mm[mm,irho,:], axis=1), axis=0)
        Fh_rate_net_vint_iceshelf_Ros_mm_int[irho,mm] = np.nansum(np.nansum(Fh_rate_delta_net_vint_iceshelf_Ros_mm[mm,irho,:], axis=1), axis=0)

# Shelf only: West Antarctica

Fs_rate_delta_adv_vint_iceshelf_WAnt_mm = wmt(salt_adv_full_vint*mask_outice*mask_WAnt*mask_WAnt, dx, dy,'budget')
Fs_rate_delta_diff_vint_iceshelf_WAnt_mm = wmt(R_s_vint*mask_outice*mask_WAnt, dx, dy,'budget')
Fs_rate_delta_net_vint_iceshelf_WAnt_mm = wmt(salt_net_full_vint*mask_outice*mask_WAnt, dx, dy,'budget')
Fs_rate_delta_sfc_iceshelf_WAnt_mm = wmt(salt_sfc*mask_outice*mask_WAnt, dx, dy,'sfc_frc')

Fh_rate_delta_adv_vint_iceshelf_WAnt_mm = wmt(temp_adv_full_vint*mask_outice*mask_WAnt, dx, dy,'budget')
Fh_rate_delta_diff_vint_iceshelf_WAnt_mm = wmt(R_t_vint*mask_outice*mask_WAnt, dx, dy,'budget')
Fh_rate_delta_net_vint_iceshelf_WAnt_mm = wmt(temp_net_full_vint*mask_outice*mask_WAnt, dx, dy,'budget')
Fh_rate_delta_sfc_iceshelf_WAnt_mm = wmt(temp_sfc*mask_outice*mask_WAnt, dx, dy,'sfc_frc')


# integrated over x, y directions: only continental shelf

Fs_rate_sfc_iceshelf_WAnt_mm_int = np.empty((len(rho_grid),12))
Fh_rate_sfc_iceshelf_WAnt_mm_int = np.empty((len(rho_grid),12))

Fs_rate_adv_vint_iceshelf_WAnt_mm_int = np.empty((len(rho_grid),12))
Fh_rate_adv_vint_iceshelf_WAnt_mm_int = np.empty((len(rho_grid),12))
Fs_rate_diff_vint_iceshelf_WAnt_mm_int = np.empty((len(rho_grid),12))
Fh_rate_diff_vint_iceshelf_WAnt_mm_int = np.empty((len(rho_grid),12))
Fs_rate_net_vint_iceshelf_WAnt_mm_int = np.empty((len(rho_grid),12))
Fh_rate_net_vint_iceshelf_WAnt_mm_int = np.empty((len(rho_grid),12))


for irho in np.arange(0,len(rho_grid)):
    for mm in np.arange(0,12):

        Fs_rate_sfc_iceshelf_WAnt_mm_int[irho,mm] = np.nansum(np.nansum(Fs_rate_delta_sfc_iceshelf_WAnt_mm[mm,irho,:], axis=1), axis=0)
        Fh_rate_sfc_iceshelf_WAnt_mm_int[irho,mm] = np.nansum(np.nansum(Fh_rate_delta_sfc_iceshelf_WAnt_mm[mm,irho,:], axis=1), axis=0)
        Fs_rate_adv_vint_iceshelf_WAnt_mm_int[irho,mm] = np.nansum(np.nansum(Fs_rate_delta_adv_vint_iceshelf_WAnt_mm[mm,irho,:], axis=1), axis=0)
        Fh_rate_adv_vint_iceshelf_WAnt_mm_int[irho,mm] = np.nansum(np.nansum(Fh_rate_delta_adv_vint_iceshelf_WAnt_mm[mm,irho,:], axis=1), axis=0)
        Fs_rate_diff_vint_iceshelf_WAnt_mm_int[irho,mm] = np.nansum(np.nansum(Fs_rate_delta_diff_vint_iceshelf_WAnt_mm[mm,irho,:], axis=1), axis=0)
        Fh_rate_diff_vint_iceshelf_WAnt_mm_int[irho,mm] = np.nansum(np.nansum(Fh_rate_delta_diff_vint_iceshelf_WAnt_mm[mm,irho,:], axis=1), axis=0)
        Fs_rate_net_vint_iceshelf_WAnt_mm_int[irho,mm] = np.nansum(np.nansum(Fs_rate_delta_net_vint_iceshelf_WAnt_mm[mm,irho,:], axis=1), axis=0)
        Fh_rate_net_vint_iceshelf_WAnt_mm_int[irho,mm] = np.nansum(np.nansum(Fh_rate_delta_net_vint_iceshelf_WAnt_mm[mm,irho,:], axis=1), axis=0)

# figures
fig_path = '/users/boeiradi/COLD_project/postprocessing/figs/WMT/'

# plot with bars
width=.023

# convert to rate per year:
Dt = 1000/0.05

# divide for Dt (sec in 1 year) and 10^6 to obtain rate in Sv year-1

#Weddell

# SFC: continental shelf
Fs_sig_sfc_iceshelf_Wed_mm = -Fs_rate_sfc_iceshelf_Wed_mm_int*Dt/1e6
Fh_sig_sfc_iceshelf_Wed_mm = -Fh_rate_sfc_iceshelf_Wed_mm_int*Dt/1e6
F_sig_sfc_iceshelf_Wed_mm = -Fs_sig_sfc_iceshelf_Wed_mm + Fh_sig_sfc_iceshelf_Wed_mm
# - calculate the anual average of the monthly ars:
Fs_sig_sfc_iceshelf_Wed =  np.nanmean(Fs_sig_sfc_iceshelf_Wed_mm, axis=1)
Fh_sig_sfc_iceshelf_Wed =  np.nanmean(Fh_sig_sfc_iceshelf_Wed_mm, axis=1)
F_sig_sfc_iceshelf_Wed = -Fs_sig_sfc_iceshelf_Wed+ Fh_sig_sfc_iceshelf_Wed

# ADV (vint)
Fs_sig_adv_vint_iceshelf_Wed_mm = -Fs_rate_adv_vint_iceshelf_Wed_mm_int*Dt/1e6
Fh_sig_adv_vint_iceshelf_Wed_mm = -Fh_rate_adv_vint_iceshelf_Wed_mm_int*Dt/1e6
F_sig_adv_vint_iceshelf_Wed_mm = -Fs_sig_adv_vint_iceshelf_Wed_mm + Fh_sig_adv_vint_iceshelf_Wed_mm
# - calculate the anual average of the monthly ars:
Fs_sig_adv_vint_iceshelf_Wed =  np.nanmean(Fs_sig_adv_vint_iceshelf_Wed_mm, axis=1)
Fh_sig_adv_vint_iceshelf_Wed =  np.nanmean(Fh_sig_adv_vint_iceshelf_Wed_mm, axis=1)
F_sig_adv_vint_iceshelf_Wed = -Fs_sig_adv_vint_iceshelf_Wed+ Fh_sig_adv_vint_iceshelf_Wed

# DIFF
Fs_sig_diff_vint_iceshelf_Wed_mm = -Fs_rate_diff_vint_iceshelf_Wed_mm_int*Dt/1e6
Fh_sig_diff_vint_iceshelf_Wed_mm = -Fh_rate_diff_vint_iceshelf_Wed_mm_int*Dt/1e6
F_sig_diff_vint_iceshelf_Wed_mm = -Fs_sig_diff_vint_iceshelf_Wed_mm + Fh_sig_diff_vint_iceshelf_Wed_mm
# - calculate the anual average of the monthly ars:
Fs_sig_diff_vint_iceshelf_Wed =  np.nanmean(Fs_sig_diff_vint_iceshelf_Wed_mm, axis=1)
Fh_sig_diff_vint_iceshelf_Wed =  np.nanmean(Fh_sig_diff_vint_iceshelf_Wed_mm, axis=1)
F_sig_diff_vin_iceshelf_Wed = -Fs_sig_diff_vint_iceshelf_Wed+ Fh_sig_diff_vint_iceshelf_Wed

# NET
Fs_sig_net_vint_iceshelf_Wed_mm = -Fs_rate_net_vint_iceshelf_Wed_mm_int*Dt/1e6
Fh_sig_net_vint_iceshelf_Wed_mm = -Fh_rate_net_vint_iceshelf_Wed_mm_int*Dt/1e6
F_sig_net_vint_iceshelf_Wed_mm = -Fs_sig_net_vint_iceshelf_Wed_mm + Fh_sig_net_vint_iceshelf_Wed_mm
# - calculate the anual average of the monthly ars:
Fs_sig_net_vint_iceshelf_Wed =  np.nanmean(Fs_sig_net_vint_iceshelf_Wed_mm, axis=1)
Fh_sig_net_vint_iceshelf_Wed =  np.nanmean(Fh_sig_net_vint_iceshelf_Wed_mm, axis=1)
F_sig_net_vint_iceshelf_Wed = -Fs_sig_net_vint_iceshelf_Wed+ Fh_sig_net_vint_iceshelf_Wed

# MAud Land

# SFC: continental shelf
Fs_sig_sfc_iceshelf_Mau_mm = -Fs_rate_sfc_iceshelf_Mau_mm_int*Dt/1e6
Fh_sig_sfc_iceshelf_Mau_mm = -Fh_rate_sfc_iceshelf_Mau_mm_int*Dt/1e6
F_sig_sfc_iceshelf_Mau_mm = -Fs_sig_sfc_iceshelf_Mau_mm + Fh_sig_sfc_iceshelf_Mau_mm
# - calculate the anual average of the monthly ars:
Fs_sig_sfc_iceshelf_Mau =  np.nanmean(Fs_sig_sfc_iceshelf_Mau_mm, axis=1)
Fh_sig_sfc_iceshelf_Mau =  np.nanmean(Fh_sig_sfc_iceshelf_Mau_mm, axis=1)
F_sig_sfc_iceshelf_Mau = -Fs_sig_sfc_iceshelf_Mau+ Fh_sig_sfc_iceshelf_Mau

# ADV (vint)
Fs_sig_adv_vint_iceshelf_Mau_mm = -Fs_rate_adv_vint_iceshelf_Mau_mm_int*Dt/1e6
Fh_sig_adv_vint_iceshelf_Mau_mm = -Fh_rate_adv_vint_iceshelf_Mau_mm_int*Dt/1e6
F_sig_adv_vint_iceshelf_Mau_mm = -Fs_sig_adv_vint_iceshelf_Mau_mm + Fh_sig_adv_vint_iceshelf_Mau_mm
# - calculate the anual average of the monthly ars:
Fs_sig_adv_vint_iceshelf_Mau =  np.nanmean(Fs_sig_adv_vint_iceshelf_Mau_mm, axis=1)
Fh_sig_adv_vint_iceshelf_Mau =  np.nanmean(Fh_sig_adv_vint_iceshelf_Mau_mm, axis=1)
F_sig_adv_vint_iceshelf_Mau = -Fs_sig_adv_vint_iceshelf_Mau+ Fh_sig_adv_vint_iceshelf_Mau

# DIFF
Fs_sig_diff_vint_iceshelf_Mau_mm = -Fs_rate_diff_vint_iceshelf_Mau_mm_int*Dt/1e6
Fh_sig_diff_vint_iceshelf_Mau_mm = -Fh_rate_diff_vint_iceshelf_Mau_mm_int*Dt/1e6
F_sig_diff_vint_iceshelf_Mau_mm = -Fs_sig_diff_vint_iceshelf_Mau_mm + Fh_sig_diff_vint_iceshelf_Mau_mm
# - calculate the anual average of the monthly ars:
Fs_sig_diff_vint_iceshelf_Mau =  np.nanmean(Fs_sig_diff_vint_iceshelf_Mau_mm, axis=1)
Fh_sig_diff_vint_iceshelf_Mau =  np.nanmean(Fh_sig_diff_vint_iceshelf_Mau_mm, axis=1)
F_sig_diff_vin_iceshelf_Mau = -Fs_sig_diff_vint_iceshelf_Mau+ Fh_sig_diff_vint_iceshelf_Mau

# NET
Fs_sig_net_vint_iceshelf_Mau_mm = -Fs_rate_net_vint_iceshelf_Mau_mm_int*Dt/1e6
Fh_sig_net_vint_iceshelf_Mau_mm = -Fh_rate_net_vint_iceshelf_Mau_mm_int*Dt/1e6
F_sig_net_vint_iceshelf_Mau_mm = -Fs_sig_net_vint_iceshelf_Mau_mm + Fh_sig_net_vint_iceshelf_Mau_mm
# - calculate the anual average of the monthly ars:
Fs_sig_net_vint_iceshelf_Mau =  np.nanmean(Fs_sig_net_vint_iceshelf_Mau_mm, axis=1)
Fh_sig_net_vint_iceshelf_Mau =  np.nanmean(Fh_sig_net_vint_iceshelf_Mau_mm, axis=1)
F_sig_net_vint_iceshelf_Mau = -Fs_sig_net_vint_iceshelf_Mau+ Fh_sig_net_vint_iceshelf_Mau

# East Antarctic

# SFC: continental shelf
Fs_sig_sfc_iceshelf_EAnt_mm = -Fs_rate_sfc_iceshelf_EAnt_mm_int*Dt/1e6
Fh_sig_sfc_iceshelf_EAnt_mm = -Fh_rate_sfc_iceshelf_EAnt_mm_int*Dt/1e6
F_sig_sfc_iceshelf_EAnt_mm = -Fs_sig_sfc_iceshelf_EAnt_mm + Fh_sig_sfc_iceshelf_EAnt_mm
# - calculate the anual average of the monthly ars:
Fs_sig_sfc_iceshelf_EAnt =  np.nanmean(Fs_sig_sfc_iceshelf_EAnt_mm, axis=1)
Fh_sig_sfc_iceshelf_EAnt =  np.nanmean(Fh_sig_sfc_iceshelf_EAnt_mm, axis=1)
F_sig_sfc_iceshelf_EAnt = -Fs_sig_sfc_iceshelf_EAnt+ Fh_sig_sfc_iceshelf_EAnt

# ADV (vint)
Fs_sig_adv_vint_iceshelf_EAnt_mm = -Fs_rate_adv_vint_iceshelf_EAnt_mm_int*Dt/1e6
Fh_sig_adv_vint_iceshelf_EAnt_mm = -Fh_rate_adv_vint_iceshelf_EAnt_mm_int*Dt/1e6
F_sig_adv_vint_iceshelf_EAnt_mm = -Fs_sig_adv_vint_iceshelf_EAnt_mm + Fh_sig_adv_vint_iceshelf_EAnt_mm
# - calculate the anual average of the monthly ars:
Fs_sig_adv_vint_iceshelf_EAnt =  np.nanmean(Fs_sig_adv_vint_iceshelf_EAnt_mm, axis=1)
Fh_sig_adv_vint_iceshelf_EAnt =  np.nanmean(Fh_sig_adv_vint_iceshelf_EAnt_mm, axis=1)
F_sig_adv_vint_iceshelf_EAnt = -Fs_sig_adv_vint_iceshelf_EAnt+ Fh_sig_adv_vint_iceshelf_EAnt

# DIFF
Fs_sig_diff_vint_iceshelf_EAnt_mm = -Fs_rate_diff_vint_iceshelf_EAnt_mm_int*Dt/1e6
Fh_sig_diff_vint_iceshelf_EAnt_mm = -Fh_rate_diff_vint_iceshelf_EAnt_mm_int*Dt/1e6
F_sig_diff_vint_iceshelf_EAnt_mm = -Fs_sig_diff_vint_iceshelf_EAnt_mm + Fh_sig_diff_vint_iceshelf_EAnt_mm
# - calculate the anual average of the monthly ars:
Fs_sig_diff_vint_iceshelf_EAnt =  np.nanmean(Fs_sig_diff_vint_iceshelf_EAnt_mm, axis=1)
Fh_sig_diff_vint_iceshelf_EAnt =  np.nanmean(Fh_sig_diff_vint_iceshelf_EAnt_mm, axis=1)
F_sig_diff_vin_iceshelf_EAnt = -Fs_sig_diff_vint_iceshelf_EAnt+ Fh_sig_diff_vint_iceshelf_EAnt

# NET
Fs_sig_net_vint_iceshelf_EAnt_mm = -Fs_rate_net_vint_iceshelf_EAnt_mm_int*Dt/1e6
Fh_sig_net_vint_iceshelf_EAnt_mm = -Fh_rate_net_vint_iceshelf_EAnt_mm_int*Dt/1e6
F_sig_net_vint_iceshelf_EAnt_mm = -Fs_sig_net_vint_iceshelf_EAnt_mm + Fh_sig_net_vint_iceshelf_EAnt_mm
# - calculate the anual average of the monthly ars:
Fs_sig_net_vint_iceshelf_EAnt =  np.nanmean(Fs_sig_net_vint_iceshelf_EAnt_mm, axis=1)
Fh_sig_net_vint_iceshelf_EAnt =  np.nanmean(Fh_sig_net_vint_iceshelf_EAnt_mm, axis=1)
F_sig_net_vint_iceshelf_EAnt = -Fs_sig_net_vint_iceshelf_EAnt+ Fh_sig_net_vint_iceshelf_EAnt

# Ross

# SFC: continental shelf
Fs_sig_sfc_iceshelf_Ros_mm = -Fs_rate_sfc_iceshelf_Ros_mm_int*Dt/1e6
Fh_sig_sfc_iceshelf_Ros_mm = -Fh_rate_sfc_iceshelf_Ros_mm_int*Dt/1e6
F_sig_sfc_iceshelf_Ros_mm = -Fs_sig_sfc_iceshelf_Ros_mm + Fh_sig_sfc_iceshelf_Ros_mm
# - calculate the anual average of the monthly ars:
Fs_sig_sfc_iceshelf_Ros =  np.nanmean(Fs_sig_sfc_iceshelf_Ros_mm, axis=1)
Fh_sig_sfc_iceshelf_Ros =  np.nanmean(Fh_sig_sfc_iceshelf_Ros_mm, axis=1)
F_sig_sfc_iceshelf_Ros = -Fs_sig_sfc_iceshelf_Ros+ Fh_sig_sfc_iceshelf_Ros

# ADV (vint)
Fs_sig_adv_vint_iceshelf_Ros_mm = -Fs_rate_adv_vint_iceshelf_Ros_mm_int*Dt/1e6
Fh_sig_adv_vint_iceshelf_Ros_mm = -Fh_rate_adv_vint_iceshelf_Ros_mm_int*Dt/1e6
F_sig_adv_vint_iceshelf_Ros_mm = -Fs_sig_adv_vint_iceshelf_Ros_mm + Fh_sig_adv_vint_iceshelf_Ros_mm
# - calculate the anual average of the monthly ars:
Fs_sig_adv_vint_iceshelf_Ros =  np.nanmean(Fs_sig_adv_vint_iceshelf_Ros_mm, axis=1)
Fh_sig_adv_vint_iceshelf_Ros =  np.nanmean(Fh_sig_adv_vint_iceshelf_Ros_mm, axis=1)
F_sig_adv_vint_iceshelf_Ros = -Fs_sig_adv_vint_iceshelf_Ros+ Fh_sig_adv_vint_iceshelf_Ros

# DIFF
Fs_sig_diff_vint_iceshelf_Ros_mm = -Fs_rate_diff_vint_iceshelf_Ros_mm_int*Dt/1e6
Fh_sig_diff_vint_iceshelf_Ros_mm = -Fh_rate_diff_vint_iceshelf_Ros_mm_int*Dt/1e6
F_sig_diff_vint_iceshelf_Ros_mm = -Fs_sig_diff_vint_iceshelf_Ros_mm + Fh_sig_diff_vint_iceshelf_Ros_mm
# - calculate the anual average of the monthly ars:
Fs_sig_diff_vint_iceshelf_Ros =  np.nanmean(Fs_sig_diff_vint_iceshelf_Ros_mm, axis=1)
Fh_sig_diff_vint_iceshelf_Ros =  np.nanmean(Fh_sig_diff_vint_iceshelf_Ros_mm, axis=1)
F_sig_diff_vin_iceshelf_Ros = -Fs_sig_diff_vint_iceshelf_Ros+ Fh_sig_diff_vint_iceshelf_Ros

# NET
Fs_sig_net_vint_iceshelf_Ros_mm = -Fs_rate_net_vint_iceshelf_Ros_mm_int*Dt/1e6
Fh_sig_net_vint_iceshelf_Ros_mm = -Fh_rate_net_vint_iceshelf_Ros_mm_int*Dt/1e6
F_sig_net_vint_iceshelf_Ros_mm = -Fs_sig_net_vint_iceshelf_Ros_mm + Fh_sig_net_vint_iceshelf_Ros_mm
# - calculate the anual average of the monthly ars:
Fs_sig_net_vint_iceshelf_Ros =  np.nanmean(Fs_sig_net_vint_iceshelf_Ros_mm, axis=1)
Fh_sig_net_vint_iceshelf_Ros =  np.nanmean(Fh_sig_net_vint_iceshelf_Ros_mm, axis=1)
F_sig_net_vint_iceshelf_Ros = -Fs_sig_net_vint_iceshelf_Ros+ Fh_sig_net_vint_iceshelf_Ros

# West Antarctic


# SFC: continental shelf
Fs_sig_sfc_iceshelf_WAnt_mm = -Fs_rate_sfc_iceshelf_WAnt_mm_int*Dt/1e6
Fh_sig_sfc_iceshelf_WAnt_mm = -Fh_rate_sfc_iceshelf_WAnt_mm_int*Dt/1e6
F_sig_sfc_iceshelf_WAnt_mm = -Fs_sig_sfc_iceshelf_WAnt_mm + Fh_sig_sfc_iceshelf_WAnt_mm
# - calculate the anual average of the monthly ars:
Fs_sig_sfc_iceshelf_WAnt =  np.nanmean(Fs_sig_sfc_iceshelf_WAnt_mm, axis=1)
Fh_sig_sfc_iceshelf_WAnt =  np.nanmean(Fh_sig_sfc_iceshelf_WAnt_mm, axis=1)
F_sig_sfc_iceshelf_WAnt = -Fs_sig_sfc_iceshelf_WAnt+ Fh_sig_sfc_iceshelf_WAnt

# ADV (vint)
Fs_sig_adv_vint_iceshelf_WAnt_mm = -Fs_rate_adv_vint_iceshelf_WAnt_mm_int*Dt/1e6
Fh_sig_adv_vint_iceshelf_WAnt_mm = -Fh_rate_adv_vint_iceshelf_WAnt_mm_int*Dt/1e6
F_sig_adv_vint_iceshelf_WAnt_mm = -Fs_sig_adv_vint_iceshelf_WAnt_mm + Fh_sig_adv_vint_iceshelf_WAnt_mm
# - calculate the anual average of the monthly ars:
Fs_sig_adv_vint_iceshelf_WAnt =  np.nanmean(Fs_sig_adv_vint_iceshelf_WAnt_mm, axis=1)
Fh_sig_adv_vint_iceshelf_WAnt =  np.nanmean(Fh_sig_adv_vint_iceshelf_WAnt_mm, axis=1)
F_sig_adv_vint_iceshelf_WAnt = -Fs_sig_adv_vint_iceshelf_WAnt+ Fh_sig_adv_vint_iceshelf_WAnt

# DIFF
Fs_sig_diff_vint_iceshelf_WAnt_mm = -Fs_rate_diff_vint_iceshelf_WAnt_mm_int*Dt/1e6
Fh_sig_diff_vint_iceshelf_WAnt_mm = -Fh_rate_diff_vint_iceshelf_WAnt_mm_int*Dt/1e6
F_sig_diff_vint_iceshelf_WAnt_mm = -Fs_sig_diff_vint_iceshelf_WAnt_mm + Fh_sig_diff_vint_iceshelf_WAnt_mm
# - calculate the anual average of the monthly ars:
Fs_sig_diff_vint_iceshelf_WAnt =  np.nanmean(Fs_sig_diff_vint_iceshelf_WAnt_mm, axis=1)
Fh_sig_diff_vint_iceshelf_WAnt =  np.nanmean(Fh_sig_diff_vint_iceshelf_WAnt_mm, axis=1)
F_sig_diff_vin_iceshelf_WAnt = -Fs_sig_diff_vint_iceshelf_WAnt+ Fh_sig_diff_vint_iceshelf_WAnt

# NET
Fs_sig_net_vint_iceshelf_WAnt_mm = -Fs_rate_net_vint_iceshelf_WAnt_mm_int*Dt/1e6
Fh_sig_net_vint_iceshelf_WAnt_mm = -Fh_rate_net_vint_iceshelf_WAnt_mm_int*Dt/1e6
F_sig_net_vint_iceshelf_WAnt_mm = -Fs_sig_net_vint_iceshelf_WAnt_mm + Fh_sig_net_vint_iceshelf_WAnt_mm
# - calculate the anual average of the monthly ars:
Fs_sig_net_vint_iceshelf_WAnt =  np.nanmean(Fs_sig_net_vint_iceshelf_WAnt_mm, axis=1)
Fh_sig_net_vint_iceshelf_WAnt =  np.nanmean(Fh_sig_net_vint_iceshelf_WAnt_mm, axis=1)
F_sig_net_vint_iceshelf_WAnt = -Fs_sig_net_vint_iceshelf_WAnt+ Fh_sig_net_vint_iceshelf_WAnt

# WMT figures:
fig_path = '/g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/figs/WMT/'

fig = plt.figure(figsize=(15,12))

ax1 = fig.add_subplot(321)
plt.title('Weddell Sea')
cs=plt.plot(rho_grid,-Fs_sig_net_vint_iceshelf_Wed,'b',label='WMT due to salt tend')
cs=plt.plot(rho_grid,-Fs_sig_diff_vint_iceshelf_Wed,'--b',label='WMT due to salt diff')
cs=plt.plot(rho_grid,-Fs_sig_adv_vint_iceshelf_Wed,':b',label='WMT due to salt adv')
plt.legend()
plt.ylabel('Transformation rate (Sv)')
plt.xlim(26.5,28),plt.ylim(-5,5)
plt.grid(True)
plt.text(26.6,.15,'Buoyancy loss')
plt.text(26.6,-.15,'Buoyancy gain')

ax2 = fig.add_subplot(322)
plt.title('Maud Land')
cs=plt.plot(rho_grid,-Fs_sig_net_vint_iceshelf_Mau,'b',label='WMT due to salt tend')
cs=plt.plot(rho_grid,-Fs_sig_diff_vint_iceshelf_Mau,'--b',label='WMT due to salt diff')
cs=plt.plot(rho_grid,-Fs_sig_adv_vint_iceshelf_Mau,':b',label='WMT due to salt adv')
plt.legend()
#plt.ylabel('Transformation rate (Sv)')
plt.xlim(26.5,28),plt.ylim(-5,5)
plt.grid(True)

ax3 = fig.add_subplot(323)
plt.title('West Antarctica')
cs=plt.plot(rho_grid,-Fs_sig_net_vint_iceshelf_EAnt,'b',label='WMT due to salt tend')
cs=plt.plot(rho_grid,-Fs_sig_diff_vint_iceshelf_EAnt,'--b',label='WMT due to salt diff')
cs=plt.plot(rho_grid,-Fs_sig_adv_vint_iceshelf_EAnt,':b',label='WMT due to salt adv')
plt.legend()
plt.ylabel('Transformation rate (Sv)')
plt.xlim(26.5,28),plt.ylim(-5,5)
plt.grid(True)

ax4 = fig.add_subplot(324)
plt.title('Ross Sea')
cs=plt.plot(rho_grid,-Fs_sig_net_vint_iceshelf_Ros,'b',label='WMT due to salt tend')
cs=plt.plot(rho_grid,-Fs_sig_diff_vint_iceshelf_Ros,'--b',label='WMT due to salt diff')
cs=plt.plot(rho_grid,-Fs_sig_adv_vint_iceshelf_Ros,':b',label='WMT due to salt adv')
plt.legend()
#plt.ylabel('Transformation rate (Sv)')
plt.xlim(26.5,28),plt.ylim(-5,5)
plt.grid(True)

ax5 = fig.add_subplot(325)
plt.title('West Antarctica')
cs=plt.plot(rho_grid,-Fs_sig_net_vint_iceshelf_WAnt,'b',label='WMT due to salt tend')
cs=plt.plot(rho_grid,-Fs_sig_diff_vint_iceshelf_WAnt,'--b',label='WMT due to salt diff')
cs=plt.plot(rho_grid,-Fs_sig_adv_vint_iceshelf_WAnt,':b',label='WMT due to salt adv')
plt.legend()
plt.ylabel('Transformation rate (Sv)')
plt.xlim(26.5,28),plt.ylim(-5,5)
plt.grid(True)

name_fig="waom4extend_shflim_S_0.25Q_10km-bathy_WMT_Full_salt_vint_annual_iceshelf_subregions.png"
plt.savefig(fig_path + name_fig, dpi=300)

fig = plt.figure(figsize=(15,12))

ax1 = fig.add_subplot(321)
plt.title('Weddell Sea')
cs=plt.plot(rho_grid,-Fh_sig_net_vint_iceshelf_Wed,'r',label='WMT due to heat tend')
cs=plt.plot(rho_grid,-Fh_sig_diff_vint_iceshelf_Wed,'--r',label='WMT due to heat diff')
cs=plt.plot(rho_grid,-Fh_sig_adv_vint_iceshelf_Wed,':r',label='WMT due to heat adv')
plt.legend()
plt.ylabel('Transformation rate (Sv)')
plt.xlim(26.5,28),plt.ylim(-5,5)
plt.grid(True)
plt.text(26.6,.15,'Buoyancy loss')
plt.text(26.6,-.15,'Buoyancy gain')

ax2 = fig.add_subplot(322)
plt.title('Maud Land')
cs=plt.plot(rho_grid,-Fh_sig_net_vint_iceshelf_Mau,'r',label='WMT due to heat tend')
cs=plt.plot(rho_grid,-Fh_sig_diff_vint_iceshelf_Mau,'--r',label='WMT due to heat diff')
cs=plt.plot(rho_grid,-Fh_sig_adv_vint_iceshelf_Mau,':r',label='WMT due to heat adv')
plt.legend()
#plt.ylabel('Transformation rate (Sv)')
plt.xlim(26.5,28),plt.ylim(-5,5)
plt.grid(True)

ax3 = fig.add_subplot(323)
plt.title('East Antarctica')
cs=plt.plot(rho_grid,-Fh_sig_net_vint_iceshelf_EAnt,'r',label='WMT due to heat tend')
cs=plt.plot(rho_grid,-Fh_sig_diff_vint_iceshelf_EAnt,'--r',label='WMT due to heat diff')
cs=plt.plot(rho_grid,-Fh_sig_adv_vint_iceshelf_EAnt,':r',label='WMT due to heat adv')
plt.legend()
plt.ylabel('Transformation rate (Sv)')
plt.xlim(26.5,28),plt.ylim(-5,5)
plt.grid(True)

ax4 = fig.add_subplot(324)
plt.title('Ross Sea')
cs=plt.plot(rho_grid,-Fh_sig_net_vint_iceshelf_Ros,'r',label='WMT due to heat tend')
cs=plt.plot(rho_grid,-Fh_sig_diff_vint_iceshelf_Ros,'--r',label='WMT due to heat diff')
cs=plt.plot(rho_grid,-Fh_sig_adv_vint_iceshelf_Ros,':r',label='WMT due to heat adv')
plt.legend()
#plt.ylabel('Transformation rate (Sv)')
plt.xlim(26.5,28),plt.ylim(-5,5)
plt.grid(True)

ax5 = fig.add_subplot(325)
plt.title('West Antarctica')
cs=plt.plot(rho_grid,-Fh_sig_net_vint_iceshelf_WAnt,'r',label='WMT due to heat tend')
cs=plt.plot(rho_grid,-Fh_sig_diff_vint_iceshelf_WAnt,'--r',label='WMT due to heat diff')
cs=plt.plot(rho_grid,-Fh_sig_adv_vint_iceshelf_WAnt,':r',label='WMT due to heat adv')
plt.legend()
plt.ylabel('Transformation rate (Sv)')
plt.xlim(26.5,28),plt.ylim(-5,5)
plt.grid(True)

name_fig="waom4extend_shflim_S_0.25Q_10km-bathy_WMT_Full_heat_vint_annual_iceshelf_subregions.png"
plt.savefig(fig_path + name_fig, dpi=300)

# --- Save transformation arrays
npy_path = '/g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/tmp_files/'

print('Saving files .....')

np.savez(npy_path + 'WAOM4extend_10km-bathy_Full_WMT_iceshelf_sfc_total_subregions',rho_grid = rho_grid, F_sig_sfc_iceshelf_Wed = F_sig_sfc_iceshelf_Wed, F_sig_sfc_iceshelf_Mau = F_sig_sfc_iceshelf_Mau, \
        F_sig_sfc_iceshelf_EAnt = F_sig_sfc_iceshelf_EAnt, F_sig_sfc_iceshelf_Ros = F_sig_sfc_iceshelf_Ros, F_sig_sfc_iceshelf_WAnt = F_sig_sfc_iceshelf_WAnt)
np.savez(npy_path + 'WAOM4extend_10km-bathy_Full_WMT_iceshelf_sfc_salt_subregions',rho_grid = rho_grid, Fs_sig_sfc_iceshelf_Wed = Fs_sig_sfc_iceshelf_Wed, Fs_sig_sfc_iceshelf_Mau = Fs_sig_sfc_iceshelf_Mau, \
        Fs_sig_sfc_iceshelf_EAnt = Fs_sig_sfc_iceshelf_EAnt, Fs_sig_sfc_iceshelf_Ros = Fs_sig_sfc_iceshelf_Ros, Fs_sig_sfc_iceshelf_WAnt = Fs_sig_sfc_iceshelf_WAnt)
np.savez(npy_path + 'WAOM4extend_10km-bathy_Full_WMT_iceshelf_sfc_heat_subregions',rho_grid = rho_grid, Fh_sig_sfc_iceshelf_Wed = Fh_sig_sfc_iceshelf_Wed, Fh_sig_sfc_iceshelf_Mau = Fh_sig_sfc_iceshelf_Mau, \
        Fh_sig_sfc_iceshelf_EAnt = Fh_sig_sfc_iceshelf_EAnt, Fh_sig_sfc_iceshelf_Ros = Fh_sig_sfc_iceshelf_Ros, Fh_sig_sfc_iceshelf_WAnt = Fh_sig_sfc_iceshelf_WAnt)

np.savez(npy_path + 'WAOM4extend_10km-bathy_Full_WMT_iceshelf_vint_salt_net_subregions',rho_grid = rho_grid, Fs_sig_net_vint_iceshelf_Wed = Fs_sig_net_vint_iceshelf_Wed, Fs_sig_net_vint_iceshelf_Mau = Fs_sig_net_vint_iceshelf_Mau, \
        Fs_sig_net_vint_iceshelf_EAnt = Fs_sig_net_vint_iceshelf_EAnt, Fs_sig_net_vint_iceshelf_Ros = Fs_sig_net_vint_iceshelf_Ros, Fs_sig_net_vint_iceshelf_WAnt = Fs_sig_net_vint_iceshelf_WAnt)
np.savez(npy_path + 'WAOM4extend_10km-bathy_Full_WMT_iceshelf_vint_salt_diff_subregions',rho_grid = rho_grid, Fs_sig_diff_vint_iceshelf_Wed = Fs_sig_diff_vint_iceshelf_Wed, Fs_sig_diff_vint_iceshelf_Mau = Fs_sig_diff_vint_iceshelf_Mau, \
        Fs_sig_diff_vint_iceshelf_EAnt = Fs_sig_diff_vint_iceshelf_EAnt, Fs_sig_diff_vint_iceshelf_Ros = Fs_sig_diff_vint_iceshelf_Ros, Fs_sig_diff_vint_iceshelf_WAnt = Fs_sig_diff_vint_iceshelf_WAnt)
np.savez(npy_path + 'WAOM4extend_10km-bathy_Full_WMT_iceshelf_vint_salt_adv_subregions',rho_grid = rho_grid, Fs_sig_adv_vint_iceshelf_Wed = Fs_sig_adv_vint_iceshelf_Wed, Fs_sig_adv_vint_iceshelf_Mau = Fs_sig_adv_vint_iceshelf_Mau, \
        Fs_sig_adv_vint_iceshelf_EAnt = Fs_sig_adv_vint_iceshelf_EAnt, Fs_sig_adv_vint_iceshelf_Ros = Fs_sig_adv_vint_iceshelf_Ros, Fs_sig_adv_vint_iceshelf_WAnt = Fs_sig_adv_vint_iceshelf_WAnt)

np.savez(npy_path + 'WAOM4extend_10km-bathy_Full_WMT_iceshelf_vint_heat_net_subregions',rho_grid = rho_grid, Fh_sig_net_vint_iceshelf_Wed = Fh_sig_net_vint_iceshelf_Wed, Fh_sig_net_vint_iceshelf_Mau = Fh_sig_net_vint_iceshelf_Mau, \
        Fh_sig_net_vint_iceshelf_EAnt = Fh_sig_net_vint_iceshelf_EAnt, Fh_sig_net_vint_iceshelf_Ros = Fh_sig_net_vint_iceshelf_Ros, Fh_sig_net_vint_iceshelf_WAnt = Fh_sig_net_vint_iceshelf_WAnt)
np.savez(npy_path + 'WAOM4extend_10km-bathy_Full_WMT_iceshelf_vint_heat_diff_subregions',rho_grid = rho_grid, Fh_sig_diff_vint_iceshelf_Wed = Fh_sig_diff_vint_iceshelf_Wed, Fh_sig_diff_vint_iceshelf_Mau = Fh_sig_diff_vint_iceshelf_Mau, \
        Fh_sig_diff_vint_iceshelf_EAnt = Fh_sig_diff_vint_iceshelf_EAnt, Fh_sig_diff_vint_iceshelf_Ros = Fh_sig_diff_vint_iceshelf_Ros, Fh_sig_diff_vint_iceshelf_WAnt = Fh_sig_diff_vint_iceshelf_WAnt)
np.savez(npy_path + 'WAOM4extend_10km-bathy_Full_WMT_iceshelf_vint_heat_adv_subregions',rho_grid = rho_grid, Fh_sig_adv_vint_iceshelf_Wed = Fh_sig_adv_vint_iceshelf_Wed, Fh_sig_adv_vint_iceshelf_Mau = Fh_sig_adv_vint_iceshelf_Mau, \
        Fh_sig_adv_vint_iceshelf_EAnt = Fh_sig_adv_vint_iceshelf_EAnt, Fh_sig_adv_vint_iceshelf_Ros = Fh_sig_adv_vint_iceshelf_Ros, Fh_sig_adv_vint_iceshelf_WAnt = Fh_sig_adv_vint_iceshelf_WAnt)

