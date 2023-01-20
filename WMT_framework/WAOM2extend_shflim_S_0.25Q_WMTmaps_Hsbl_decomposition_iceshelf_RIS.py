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

dx = xr.open_dataset('/g/data3/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom2extend_shflim_S_0.25Q/output_yr5_diag/MLD_vint_vars_for_WMT.nc')
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

# calculate the LHS term in Pellichero et al (2018):
# ps: Diffusion (R_s, R_t) terms already include the sfc fluxes

# heat (eqn 5)
rho0 = 1025 #1000
Cp = 3985

# convert salt convergenvces to equivalent FW
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

    F_rate_delta_var_vint_mm = np.empty((12,len_rho_grid,2800,3150))

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
# Heat: m.degC/s -> m3.degC/s
# Salt: m/s -> m3/s
# Fwf: Kg.m-2.s-1 = Kg/s

# Shelf only: excluding open ocean

Fs_rate_delta_adv_vint_iceshelf_mm = wmt(salt_adv_mld_vint*mask_outice, dx, dy,'budget')
Fs_rate_delta_diff_vint_iceshelf_mm = wmt(R_s_vint*mask_outice, dx, dy,'budget')
Fs_rate_delta_net_vint_iceshelf_mm = wmt(salt_net_mld_vint*mask_outice, dx, dy,'budget')
Fs_rate_delta_sfc_iceshelf_mm = wmt(salt_sfc*mask_outice, dx, dy,'sfc_frc')

Fm_rate_delta_sfc_iceshelf_mm = wmt(beta_mld*m*mask_outice, dx, dy,'sfc_frc')
    
#Fh_rate_delta_adv_vint_iceshelf_mm = wmt(temp_adv_mld_vint*mask_outice, dx, dy,'budget')
#Fh_rate_delta_diff_vint_iceshelf_mm = wmt(R_t_vint*mask_outice, dx, dy,'budget')
#Fh_rate_delta_net_vint_iceshelf_mm = wmt(temp_net_mld_vint*mask_outice, dx, dy,'budget')
#Fh_rate_delta_sfc_iceshelf_mm = wmt(temp_sfc*mask_outice, dx, dy,'sfc_frc')

# figures
fig_path = '/g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/figs/WMT/'

# plot with bars
width=.023

# convert to rate per year:
Dt = rho0

### plot maps:
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

# test new plot: sum irho from 28-34 (27.2 - 27.7 kg.m-3)
# only salt

# 1) annual mean

print
for irho in np.arange(28,34):
#    salt_net_tmp = np.nanmean(Fs_rate_delta_net_vint_iceshelf_mm[:,irho,:], axis=0)*rho0/(dx*dy)
    #salt_net_tmp = Fs_rate_delta_net_vint_iceshelf_mm[:,irho,:]
    salt_net_tmp = np.empty((12,2800,3150))
    for mm in np.arange(0,12):
        salt_net_tmp[mm,:] = Fs_rate_delta_net_vint_iceshelf_mm[mm,irho,:]*rho0/(dx*dy)
    salt_net_tmp = np.expand_dims(salt_net_tmp, axis=0)
    if irho >= 29:
        salt_net = np.concatenate((salt_net,salt_net_tmp), axis=0)
    else:
        salt_net = salt_net_tmp
    del salt_net_tmp

    #salt_adv_tmp = Fs_rate_delta_adv_vint_iceshelf_mm[:,irho,:]
    salt_adv_tmp = np.empty((12,2800,3150))
    for mm in np.arange(0,12):
        salt_adv_tmp[mm,:] = Fs_rate_delta_adv_vint_iceshelf_mm[mm,irho,:]*rho0/(dx*dy)
    salt_adv_tmp = np.expand_dims(salt_adv_tmp, axis=0)
    if irho >= 29:
        salt_adv = np.concatenate((salt_adv,salt_adv_tmp), axis=0)
    else:
        salt_adv = salt_adv_tmp
    del salt_adv_tmp

    #salt_diff_tmp = Fs_rate_delta_diff_vint_iceshelf_mm[:,irho,:]
    salt_diff_tmp = np.empty((12,2800,3150))
    for mm in np.arange(0,12):
        salt_diff_tmp[mm,:] = Fs_rate_delta_diff_vint_iceshelf_mm[mm,irho,:]*rho0/(dx*dy)
    salt_diff_tmp = np.expand_dims(salt_diff_tmp, axis=0)
    if irho >= 29:
        salt_diff = np.concatenate((salt_diff,salt_diff_tmp), axis=0)
    else:
        salt_diff = salt_diff_tmp
    del salt_diff_tmp

    #salt_sfc_tmp = Fs_rate_delta_sfc_iceshelf_mm[:,irho,:]
    salt_sfc_tmp = np.empty((12,2800,3150))
    for mm in np.arange(0,12):
        salt_sfc_tmp[mm,:] = Fs_rate_delta_sfc_iceshelf_mm[mm,irho,:]*rho0/(dx*dy)
    salt_sfc_tmp = np.expand_dims(salt_sfc_tmp, axis=0)
    if irho >= 29:
        salt_sfc = np.concatenate((salt_sfc,salt_sfc_tmp), axis=0)
    else:
        salt_sfc = salt_sfc_tmp
    del salt_sfc_tmp

    m_sfc_tmp = np.empty((12,2800,3150))
    for mm in np.arange(0,12):
        m_sfc_tmp[mm,:] = Fm_rate_delta_sfc_iceshelf_mm[mm,irho,:]*rho0/(dx*dy)
    m_sfc_tmp = np.expand_dims(m_sfc_tmp, axis=0)
    if irho >= 29:
        m_sfc = np.concatenate((m_sfc,m_sfc_tmp), axis=0)
    else:
        m_sfc = m_sfc_tmp
    del m_sfc_tmp

print('salt_net size, after concatenating = ',salt_net.shape)

salt_net_sum = np.nansum(salt_net, axis=0)
salt_adv_sum = np.nansum(salt_adv, axis=0)
salt_diff_sum = np.nansum(salt_diff, axis=0)
salt_sfc_sum = np.nansum(salt_sfc, axis=0)
m_sfc_sum = np.nansum(m_sfc, axis=0)

proj = ccrs.SouthPolarStereo()
fig = plt.figure(figsize=(5,10))

ax1 = fig.add_subplot(311, projection=proj)
plt.title('Salt tendency, $\sigma_{\Theta}$ = ' + str(np.around(rho_grid[28],decimals=2)) + ':'  + str(np.around(rho_grid[34],decimals=2)), fontsize=10)
cy=plt.pcolormesh(lon_rho,lat_rho,np.nanmean(salt_net_sum, axis=0), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm)
plt.colorbar(cy, extend='both')
ax1.gridlines()
lonlat_labels(ax1)
ax1.set_extent([130, 280, -85, -70], crs=ccrs.PlateCarree())
ax1.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white')
plt.clim(-5e-6,5e-6)

ax2 = fig.add_subplot(312, projection=proj)
plt.title('Salt advection, $\sigma_{\Theta}$ = ' + str(np.around(rho_grid[28],decimals=2)) + ':'  + str(np.around(rho_grid[34],decimals=2)), fontsize=10)
cy=plt.pcolormesh(lon_rho,lat_rho,np.nanmean(salt_adv_sum, axis=0), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm)
plt.colorbar(cy, extend='both')
ax2.gridlines()
lonlat_labels(ax2)
ax2.set_extent([130, 280, -85, -70], crs=ccrs.PlateCarree())
ax2.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white')
plt.clim(-5e-6,5e-6)

ax3 = fig.add_subplot(313, projection=proj)
plt.title('Salt diffusion, $\sigma_{\Theta}$ = ' + str(np.around(rho_grid[28],decimals=2)) + ':'  + str(np.around(rho_grid[34],decimals=2)), fontsize=10)
cy=plt.pcolormesh(lon_rho,lat_rho,np.nanmean(salt_diff_sum, axis=0), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm)
plt.colorbar(cy, extend='both')
ax3.gridlines()
lonlat_labels(ax3)
ax3.set_extent([130, 280, -85, -70], crs=ccrs.PlateCarree())
ax3.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white')
plt.clim(-5e-6,5e-6)

#ax4 = fig.add_subplot(314, projection=proj)
#plt.title('Salt sfc flux, $\sigma_{\Theta}$ = ' + str(np.around(rho_grid[28],decimals=2)) + ':'  + str(np.around(rho_grid[34],decimals=2)), fontsize=10)
#cy=plt.pcolormesh(lon_rho,lat_rho,np.nanmean(salt_sfc_sum, axis=0), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm)
#plt.colorbar(cy, extend='both')
#ax4.gridlines()
#lonlat_labels(ax4)
#ax4.set_extent([130, 280, -85, -70], crs=ccrs.PlateCarree())
#ax4.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white')
#plt.clim(-5e-6,5e-6)

name_fig="waom2extend_shflim_S_0.25Q_WMTmaps_salt_budget_annual_yr10_iceshelf_RIS.png"
plt.savefig(fig_path + name_fig, dpi=300)
plt.close()

proj = ccrs.SouthPolarStereo()
fig = plt.figure(figsize=(5,10))

ax1 = fig.add_subplot(311, projection=proj)
plt.title('Basal melt, $\sigma_{\Theta}$ = ' + str(np.around(rho_grid[28],decimals=2)) + ':'  + str(np.around(rho_grid[34],decimals=2)), fontsize=10)
cy=plt.pcolormesh(lon_rho,lat_rho,-np.nanmean(m_sfc_sum, axis=0), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm)
plt.colorbar(cy, extend='both')
ax1.gridlines()
lonlat_labels(ax1)
ax1.set_extent([130, 280, -85, -70], crs=ccrs.PlateCarree())
ax1.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white')
plt.clim(-5e-6,5e-6)

name_fig="waom2extend_shflim_S_0.25Q_WMTmaps_salt_budget_annual_yr20_basalmelt_iceshelf_RIS.png"
plt.savefig(fig_path + name_fig, dpi=300)
plt.close()


# monthly plots:
for mm in np.arange(0,12):
    salt_net_sum_mm = salt_net_sum[mm,:,:]
    salt_adv_sum_mm = salt_adv_sum[mm,:,:]
    salt_diff_sum_mm = salt_diff_sum[mm,:,:]
    salt_sfc_sum_mm = salt_sfc_sum[mm,:,:]

    proj = ccrs.SouthPolarStereo()
    fig = plt.figure(figsize=(5,10))

    ax1 = fig.add_subplot(311, projection=proj)
    plt.title('Salt tendency, $\sigma_{\Theta}$ = ' + str(np.around(rho_grid[28],decimals=2)) + ':'  + str(np.around(rho_grid[34],decimals=2)), fontsize=10)
    cy=plt.pcolormesh(lon_rho,lat_rho,salt_net_sum_mm, transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm)
    plt.colorbar(cy, extend='both')
    ax1.gridlines()
    lonlat_labels(ax1)
    ax1.set_extent([130, 280, -85, -70], crs=ccrs.PlateCarree())
    ax1.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white')
    plt.clim(-5e-6,5e-6)

    ax2 = fig.add_subplot(312, projection=proj)
    plt.title('Salt advection, $\sigma_{\Theta}$ = ' + str(np.around(rho_grid[28],decimals=2)) + ':'  + str(np.around(rho_grid[34],decimals=2)), fontsize=10)
    cy=plt.pcolormesh(lon_rho,lat_rho,salt_adv_sum_mm, transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm)
    plt.colorbar(cy, extend='both')
    ax2.gridlines()
    lonlat_labels(ax2)
    ax2.set_extent([130, 280, -85, -70], crs=ccrs.PlateCarree())
    ax2.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white')
    plt.clim(-5e-6,5e-6)

    ax3 = fig.add_subplot(313, projection=proj)
    plt.title('Salt diffusion, $\sigma_{\Theta}$ = ' + str(np.around(rho_grid[28],decimals=2)) + ':'  + str(np.around(rho_grid[34],decimals=2)), fontsize=10)
    cy=plt.pcolormesh(lon_rho,lat_rho,salt_diff_sum_mm, transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm)
    plt.colorbar(cy, extend='both')
    ax3.gridlines()
    lonlat_labels(ax3)
    ax3.set_extent([130, 280, -85, -70], crs=ccrs.PlateCarree())
    ax3.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white')
    plt.clim(-5e-6,5e-6)

#    ax4 = fig.add_subplot(314, projection=proj)
#    plt.title('Salt sfc flux, $\sigma_{\Theta}$ = ' + str(np.around(rho_grid[28],decimals=2)) + ':'  + str(np.around(rho_grid[34],decimals=2)), fontsize=10)
#    cy=plt.pcolormesh(lon_rho,lat_rho,salt_sfc_sum_mm, transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm)
#    plt.colorbar(cy, extend='both')
#    ax4.gridlines()
#    lonlat_labels(ax4)
#    ax4.set_extent([130, 280, -85, -70], crs=ccrs.PlateCarree())
#    ax4.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='white')
#    plt.clim(-5e-6,5e-6)

    name_fig="waom2extend_shflim_S_0.25Q_WMTmaps_salt_budget_annual_yr10_iceshelf_mm=" + str(mm) + "_RIS.png"
    plt.savefig(fig_path + name_fig, dpi=300)
    plt.close()
