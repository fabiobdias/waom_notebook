# read nc output from WAOM 10km run

import xarray as xr
import pandas as p
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

#import iris
#import iris.iterate
#import iris.coords
#import iris.plot as iplt
import gsw

# load ROMS avg output
for mm  in ['01','02','03','04','05','06','07','08','09','10','11','12']:
    ds = xr.open_dataset('/g/data3/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom2extend_shflim_S_0.25Q/output_yr5_diag/ocean_avg_00' + mm + '.nc')
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

di = xr.open_dataset('/g/data3/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom2extend_shflim_S_0.25Q/output_yr5_diag/ocean_avg_0001.nc')
ice_draft = di.variables["zice"]
h = di.variables["h"]
mask_zice = ma.masked_where(ice_draft < 0, np.ones(ice_draft.shape))
mask_outice = ma.masked_where(ice_draft >= 0, np.ones(ice_draft.shape))


di.close()

# calculate dz following:
dz = np.empty((12,2800,3150,31))
dz_inv = np.empty((12,2800,3150,31))


for tt in np.arange(0,12):
    z_w_sorted = -1*z_w[tt,:,:,::-1]
    print(z_w_sorted.shape)
    dz_inv[tt,:,:,:] = np.diff(z_w_sorted,axis=2)
    dz[tt,:,:,:] = dz_inv[tt,:,:,::-1]

dg = xr.open_dataset("/g/data3/hh5/tmp/access-om/fbd581/ROMS/waom2_frc/waom2extend_grd.nc")

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

dx = xr.open_dataset('/g/data3/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom2extend_shflim_S_0.25Q/output_yr5_diag/Full_vint_vars_for_WMT_m.s-1.nc')

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

    for mm in np.arange(0,12):
        if var_type == 'budget':
            # Zika advices  not to multiply by area
            #F_rate_var_vint[mm,:] = dx*dy*var_int.isel(times=mm)
            F_rate_var_vint[mm,:] = var_int.isel(times=mm)
        elif var_type == 'sfc_frc':
            #F_rate_var_vint[mm,:] = dx*dy*var_int[mm,:,:]
            F_rate_var_vint[mm,:] = dx*dy*var_int[mm,:,:]

    print(F_rate_var_vint.shape)

    F_rate_delta_var_vint_mm = np.empty((12,len_rho_grid,2800,3150))

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

# -- Units --
# Heat: m.degC/s -> m3.degC/s
# Salt: m/s -> m3/s
# Fwf: Kg.m-2.s-1 = Kg/s

# Shelf only: excluding open ocean

Fs_rate_delta_adv_vint_shelf_mm = wmt(salt_adv_full_vint*mask_shelf*mask_zice, dx, dy,'budget')
Fs_rate_delta_diff_vint_shelf_mm = wmt(R_s_vint*mask_shelf*mask_zice, dx, dy,'budget')
Fs_rate_delta_net_vint_shelf_mm = wmt(salt_net_full_vint*mask_shelf*mask_zice, dx, dy,'budget')
Fs_rate_delta_sfc_shelf_mm = wmt(salt_sfc*mask_shelf*mask_zice, dx, dy,'sfc_frc')

# decomposition advection and diffusion into vertical and horizontal components
Fs_rate_delta_vadv_vint_shelf_mm = wmt(salt_vdia_adv_full_vint*mask_shelf*mask_zice, dx, dy,'budget')
Fs_rate_delta_hadv_vint_shelf_mm = wmt(salt_hdia_adv_full_vint*mask_shelf*mask_zice, dx, dy,'budget')
Fs_rate_delta_vdiff_vint_shelf_mm = wmt(salt_vdia_diff_full_vint*mask_shelf*mask_zice, dx, dy,'budget')
Fs_rate_delta_hdiff_vint_shelf_mm = wmt(salt_hdia_diff_full_vint*mask_shelf*mask_zice, dx, dy,'budget')

#Fh_rate_delta_adv_vint_shelf_mm = wmt(temp_adv_full_vint*mask_shelf*mask_zice, dx, dy,'budget')
#Fh_rate_delta_diff_vint_shelf_mm = wmt(R_t_vint*mask_shelf*mask_zice, dx, dy,'budget')
#Fh_rate_delta_net_vint_shelf_mm = wmt(temp_net_full_vint*mask_shelf*mask_zice, dx, dy,'budget')
#Fh_rate_delta_sfc_shelf_mm = wmt(temp_sfc*mask_shelf*mask_zice, dx, dy,'sfc_frc')

# figures
fig_path = '/g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/figs/WMT/'

# plot with bars
width=.023

# convert to rate per year:
Dt = 1000/0.05

for irho in np.arange(19,40): # 27.0:27.9 kg m-3 for sigma interval 0.05 (to save maps)

#    salt_net_tmp = np.nanmean(Fs_rate_delta_net_vint_shelf_mm[:,irho,:], axis=0)*Dt/(dx*dy)
    salt_net_tmp = np.empty((12,2800,3150))
    for mm in np.arange(0,12):
        salt_net_tmp[mm,:] = Fs_rate_delta_net_vint_shelf_mm[mm,irho,:]*Dt/(dx*dy)
    salt_net_tmp = np.expand_dims(salt_net_tmp, axis=0)
    if irho >= 20:
        salt_net = np.concatenate((salt_net,salt_net_tmp), axis=0)
    else:
        salt_net = salt_net_tmp
    del salt_net_tmp

    #salt_adv_tmp = Fs_rate_delta_adv_vint_shelf_mm[:,irho,:]
    salt_adv_tmp = np.empty((12,2800,3150))
    for mm in np.arange(0,12):
        salt_adv_tmp[mm,:] = Fs_rate_delta_adv_vint_shelf_mm[mm,irho,:]*Dt/(dx*dy)
    salt_adv_tmp = np.expand_dims(salt_adv_tmp, axis=0)
    if irho >= 20:
        salt_adv = np.concatenate((salt_adv,salt_adv_tmp), axis=0)
    else:
        salt_adv = salt_adv_tmp
    del salt_adv_tmp

    #salt_diff_tmp = Fs_rate_delta_diff_vint_shelf_mm[:,irho,:]
    salt_diff_tmp = np.empty((12,2800,3150))
    for mm in np.arange(0,12):
        salt_diff_tmp[mm,:] = Fs_rate_delta_diff_vint_shelf_mm[mm,irho,:]*Dt/(dx*dy)
    salt_diff_tmp = np.expand_dims(salt_diff_tmp, axis=0)
    if irho >= 20:
        salt_diff = np.concatenate((salt_diff,salt_diff_tmp), axis=0)
    else:
        salt_diff = salt_diff_tmp
    del salt_diff_tmp

# VDIFF
    salt_vdiff_tmp = np.empty((12,2800,3150))
    for mm in np.arange(0,12):
        salt_vdiff_tmp[mm,:] = Fs_rate_delta_vdiff_vint_shelf_mm[mm,irho,:]*Dt/(dx*dy)
    salt_vdiff_tmp = np.expand_dims(salt_vdiff_tmp, axis=0)
    if irho >= 20:
        salt_vdiff = np.concatenate((salt_vdiff,salt_vdiff_tmp), axis=0)
    else:
        salt_vdiff = salt_vdiff_tmp
    del salt_vdiff_tmp

# HDIFF
    salt_hdiff_tmp = np.empty((12,2800,3150))
    for mm in np.arange(0,12):
        salt_hdiff_tmp[mm,:] = Fs_rate_delta_hdiff_vint_shelf_mm[mm,irho,:]*Dt/(dx*dy)
    salt_hdiff_tmp = np.expand_dims(salt_hdiff_tmp, axis=0)
    if irho >= 20:
        salt_hdiff = np.concatenate((salt_hdiff,salt_hdiff_tmp), axis=0)
    else:
        salt_hdiff = salt_hdiff_tmp
    del salt_hdiff_tmp

    #salt_sfc_tmp = Fs_rate_delta_sfc_shelf_mm[:,irho,:]
    salt_sfc_tmp = np.empty((12,2800,3150))
    for mm in np.arange(0,12):
        salt_sfc_tmp[mm,:] = Fs_rate_delta_sfc_shelf_mm[mm,irho,:]*Dt/(dx*dy)
    salt_sfc_tmp = np.expand_dims(salt_sfc_tmp, axis=0)
    if irho >= 20:
        salt_sfc = np.concatenate((salt_sfc,salt_sfc_tmp), axis=0)
    else:
        salt_sfc = salt_sfc_tmp
    del salt_sfc_tmp


print('salt_net size, after concatenating = ',salt_net.shape)

## saving annual mean transformation maps for each density classes between 27.0 and 27.9 kg m-3:

fn = '/g/data3/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom2extend_shflim_S_0.25Q/output_yr5_diag/Full_vint_WMTmaps_shelf.nc'
#fn = 'Full_vint_WMTmaps.nc'
dx = Dataset(fn, 'w', format='NETCDF4')

times = dx.createDimension('times', 1)
rho = dx.createDimension('rho', 21)
xi_rho = dx.createDimension('xi_rho', 3150)
eta_rho = dx.createDimension('eta_rho', 2800)

ocean_times = dx.createVariable('times', 'f4', ('times',))
ocean_times.units = 'seconds of the year'
irho = dx.createVariable('rho', 'f4', ('rho',))
irho.units = 'kg m-3'
xi = dx.createVariable('xi_rho', 'f4', ('xi_rho',))
eta = dx.createVariable('eta_rho', 'f4', ('eta_rho',))

# transformation maps variables:
map_salt_net = dx.createVariable('map_salt_net', 'f4', ('rho','times', 'eta_rho', 'xi_rho',))
map_salt_net.units = 'transformation maps salt tendency full-depth integrated (m^3 s)'

map_salt_adv = dx.createVariable('map_salt_adv', 'f4', ('rho','times', 'eta_rho', 'xi_rho',))
map_salt_adv.units = 'transformation maps salt total advection full-depth integrated (m^3 s)'

map_salt_diff = dx.createVariable('map_salt_diff', 'f4', ('rho','times', 'eta_rho', 'xi_rho',))
map_salt_diff.units = 'transformation maps salt total diffusion full-depth integrated (m^3 s)'

map_salt_hdiff = dx.createVariable('map_salt_hdiff', 'f4', ('rho','times', 'eta_rho', 'xi_rho',))
map_salt_hdiff.units = 'transformation maps salt horizontal diffusion full-depth integrated (m^3 s)'

map_salt_vdiff = dx.createVariable('map_salt_vdiff', 'f4', ('rho','times', 'eta_rho', 'xi_rho',))
map_salt_vdiff.units = 'transformation maps salt vertical diffusion full-depth integrated (m^3 s)'

map_salt_sfc = dx.createVariable('map_salt_sfc', 'f4', ('rho','times', 'eta_rho', 'xi_rho',))
map_salt_sfc.units = 'transformation maps surface salt flux integrated (m^3 s)'

irho[:] = rho_grid[19:40]
ocean_times[:] = 0 #annual   #np.arange(1328000,31536000,2628000) # monthly
xi[:] = np.arange(0,3150)
eta[:] = np.arange(0,2800)

map_salt_net[:] = np.nanmean(salt_net, axis=1)
map_salt_adv[:] = np.nanmean(salt_adv, axis=1)
map_salt_diff[:] = np.nanmean(salt_diff, axis=1)
map_salt_hdiff[:] = np.nanmean(salt_hdiff, axis=1)
map_salt_vdiff[:] = np.nanmean(salt_vdiff, axis=1)
map_salt_sfc[:] = np.nanmean(salt_sfc, axis=1)

dx.close()
