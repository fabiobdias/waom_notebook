# Vertically-integrated heat/salt budgets over the MLD
# multiplied by heat expansion/salinty contraction coeffs = units [m.s-1]

#import venv

import xarray as xr
import pandas as p
import numpy as np
import numpy.ma as ma
#import cartopy.crs as ccrs
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

# load ROMS avg output
for mm  in ['01','02','03','04','05','06','07','08','09','10','11','12']:
    ds = xr.open_dataset('/scratch/project_2000789/boeiradi/waom10extend_shflim_S_0.25Q/output_20yr_diag/ocean_avg_00' + mm + '.nc')
    print(ds.variables["temp"].shape)
    temp_tmp = np.nanmean(ds.variables["temp"], axis=0)
    salt_tmp = np.nanmean(ds.variables["salt"], axis=0)
    Hsbl_tmp = np.nanmean(ds.variables["Hsbl"], axis=0)
    
    ds = ds.set_coords(['Cs_r', 'Cs_w', 'hc', 'h', 'Vtransform'])
    if ds.Vtransform == 1:
        Zo_rho = ds.hc * (ds.s_rho - ds.Cs_r) + ds.Cs_r * ds.h
        z_rho_tmp = Zo_rho + ds.zeta * (1 + Zo_rho/ds.h) + ds.zice
        print("Vtransform=1")
    elif ds.Vtransform == 2:
        Zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * ds.h) / (ds.hc + ds.h)
        z_rho_tmp = ds.zeta + (ds.zeta + ds.h) * Zo_rho   + ds.zice
        print("Vtransform=2")
        Zo_w = (ds.hc * ds.s_w + ds.Cs_w * ds.h) / (ds.hc + ds.h)
        z_w_tmp = ds.zeta + (ds.zeta + ds.h) * Zo_w   + ds.zice
        
    z_rho_avg = np.nanmean(z_rho_tmp, axis=0)
    z_w_avg = np.nanmean(z_w_tmp, axis=0)
    
    # concatanate monthly avgs into a yearly variable
    if mm == '01':
        temp = temp_tmp
        salt = salt_tmp
        Hsbl = Hsbl_tmp
        z_rho = z_rho_avg
        z_w = z_w_avg
    elif mm == '02':
        temp = np.stack((temp,temp_tmp), axis=0)
        salt = np.stack((salt,salt_tmp), axis=0)
        Hsbl = np.stack((Hsbl,Hsbl_tmp), axis=0)
        z_rho = np.stack((z_rho,z_rho_avg), axis=0)
        z_w = np.stack((z_w,z_w_avg), axis=0)
    else:
        temp_tmp_4thdim = np.expand_dims(temp_tmp, axis=0)
        temp = np.concatenate((temp,temp_tmp_4thdim), axis=0)
        salt_tmp_4thdim = np.expand_dims(salt_tmp, axis=0)
        salt = np.concatenate((salt,salt_tmp_4thdim), axis=0)
        Hsbl_tmp_4thdim = np.expand_dims(Hsbl_tmp, axis=0)
        Hsbl = np.concatenate((Hsbl,Hsbl_tmp_4thdim), axis=0)    
        z_rho_tmp_4thdim = np.expand_dims(z_rho_avg, axis=0)
        z_rho = np.concatenate((z_rho,z_rho_tmp_4thdim), axis=0)    
        z_w_tmp_4thdim = np.expand_dims(z_w_avg, axis=0)
        z_w = np.concatenate((z_w,z_w_tmp_4thdim), axis=0)
            
#ds.coords['z_rho'] = z_rho#.transpose() # put z_rho into ds dataset


dg = xr.open_dataset("/scratch/project_2000339/boeiradi/waom10_frc/waom10extend_grd.nc")

lat_rho = dg.variables["lat_rho"]
lon_rho = dg.variables["lon_rho"]
ds.coords['lat_rho']=lat_rho.transpose() # put lat_rho into ds dataset
ds.coords['lon_rho']=lon_rho.transpose() # put lon_rho into ds dataset

sigma_t = gsw.rho(salt,temp,0) - 1000

print(sigma_t.shape)

diff_sigma = np.empty(sigma_t.shape)
for mm in np.arange(0,12):
    for kk in np.arange(0,30):
       # print(kk)
        diff_sigma[mm,kk,:,:] = np.squeeze(sigma_t[mm,kk,:,:]) - np.squeeze(sigma_t[mm,30,:,:])
        
## read diffusion/advection of temperature - horizontal and vertical terms:

for mm  in ['01','02','03','04','05','06','07','08','09','10','11','12']:
    dd = xr.open_dataset('/scratch/project_2000789/boeiradi/waom10extend_shflim_S_0.25Q/output_20yr_diag/ocean_dia_00' + mm + '.nc')
    print(dd.variables["temp_vdiff"].shape)
    temp_hdiff_tmp = np.nanmean(dd.variables["temp_hdiff"], axis=0)
    temp_vdiff_tmp = np.nanmean(dd.variables["temp_vdiff"], axis=0)
    temp_hadv_tmp = np.nanmean(dd.variables["temp_hadv"], axis=0)
    temp_vadv_tmp = np.nanmean(dd.variables["temp_vadv"], axis=0)
    
     # concatenate monthly avgs into a yearly variable
    if mm == '01':
        temp_hdiff = temp_hdiff_tmp
        temp_vdiff = temp_vdiff_tmp
        temp_hadv = temp_hadv_tmp
        temp_vadv = temp_vadv_tmp
    elif mm == '02':
        temp_hdiff = np.stack((temp_hdiff,temp_hdiff_tmp), axis=0)
        temp_vdiff = np.stack((temp_vdiff,temp_vdiff_tmp), axis=0)
        temp_hadv = np.stack((temp_hadv,temp_hadv_tmp), axis=0)
        temp_vadv = np.stack((temp_vadv,temp_vadv_tmp), axis=0)
    else:
        temp_hdiff_tmp_4thdim = np.expand_dims(temp_hdiff_tmp, axis=0)
        temp_hdiff = np.concatenate((temp_hdiff,temp_hdiff_tmp_4thdim), axis=0)
        temp_vdiff_tmp_4thdim = np.expand_dims(temp_vdiff_tmp, axis=0)
        temp_vdiff = np.concatenate((temp_vdiff,temp_vdiff_tmp_4thdim), axis=0) 
        temp_hadv_tmp_4thdim = np.expand_dims(temp_hadv_tmp, axis=0)
        temp_hadv= np.concatenate((temp_hadv,temp_hadv_tmp_4thdim), axis=0)
        temp_vadv_tmp_4thdim = np.expand_dims(temp_vadv_tmp, axis=0)
        temp_vadv = np.concatenate((temp_vadv,temp_vadv_tmp_4thdim), axis=0) 
        
    dd.close()

## read diffusion/advection of temperature - horizontal and vertical terms:

for mm  in ['01','02','03','04','05','06','07','08','09','10','11','12']:
    dd = xr.open_dataset('/scratch/project_2000789/boeiradi/waom10extend_shflim_S_0.25Q/output_20yr_diag/ocean_dia_00' + mm + '.nc')
    print(dd.variables["salt_vdiff"].shape)
    salt_hdiff_tmp = np.nanmean(dd.variables["salt_hdiff"], axis=0)
    salt_vdiff_tmp = np.nanmean(dd.variables["salt_vdiff"], axis=0)
    salt_hadv_tmp = np.nanmean(dd.variables["salt_hadv"], axis=0)
    salt_vadv_tmp = np.nanmean(dd.variables["salt_vadv"], axis=0)
    
     # concatenate monthly avgs into a yearly variable
    if mm == '01':
        salt_hdiff = salt_hdiff_tmp
        salt_vdiff = salt_vdiff_tmp
        salt_hadv = salt_hadv_tmp
        salt_vadv = salt_vadv_tmp
    elif mm == '02':
        salt_hdiff = np.stack((salt_hdiff,salt_hdiff_tmp), axis=0)
        salt_vdiff = np.stack((salt_vdiff,salt_vdiff_tmp), axis=0)
        salt_hadv = np.stack((salt_hadv,salt_hadv_tmp), axis=0)
        salt_vadv = np.stack((salt_vadv,salt_vadv_tmp), axis=0)
    else:
        salt_hdiff_tmp_4thdim = np.expand_dims(salt_hdiff_tmp, axis=0)
        salt_hdiff = np.concatenate((salt_hdiff,salt_hdiff_tmp_4thdim), axis=0)
        salt_vdiff_tmp_4thdim = np.expand_dims(salt_vdiff_tmp, axis=0)
        salt_vdiff = np.concatenate((salt_vdiff,salt_vdiff_tmp_4thdim), axis=0) 
        salt_hadv_tmp_4thdim = np.expand_dims(salt_hadv_tmp, axis=0)
        salt_hadv= np.concatenate((salt_hadv,salt_hadv_tmp_4thdim), axis=0)
        salt_vadv_tmp_4thdim = np.expand_dims(salt_vadv_tmp, axis=0)
        salt_vadv = np.concatenate((salt_vadv,salt_vadv_tmp_4thdim), axis=0) 
        
    dd.close()

## read tendency terms for temperature and salinity:

for mm  in ['01','02','03','04','05','06','07','08','09','10','11','12']:
    dd = xr.open_dataset('/scratch/project_2000789/boeiradi/waom10extend_shflim_S_0.25Q/output_20yr_diag/ocean_dia_00' + mm + '.nc')
    salt_rate_tmp = np.nanmean(dd.variables["salt_rate"], axis=0)
    temp_rate_tmp = np.nanmean(dd.variables["temp_rate"], axis=0)
    
    # concatanate monthly avgs into a yearly variable
    if mm == '01':
        salt_rate = salt_rate_tmp
        temp_rate = temp_rate_tmp
    elif mm == '02':
        salt_rate = np.stack((salt_rate,salt_rate_tmp), axis=0)
        temp_rate = np.stack((temp_rate,temp_rate_tmp), axis=0)
    else:
        salt_rate_tmp_4thdim = np.expand_dims(salt_rate_tmp, axis=0)
        salt_rate = np.concatenate((salt_rate,salt_rate_tmp_4thdim), axis=0)
        temp_rate_tmp_4thdim = np.expand_dims(temp_rate_tmp, axis=0)
        temp_rate = np.concatenate((temp_rate,temp_rate_tmp_4thdim), axis=0)
        
    dd.close()

# calculate model thickness for vertical integrals:
# use z_w_avg (monthly)

# calculate dz following:
dz = np.empty((12,560,630,31))
dz_inv = np.empty((12,560,630,31))


for tt in np.arange(0,12):
    z_w_sorted = -1*z_w[tt,:,:,::-1]
    print(z_w_sorted.shape)
    dz_inv[tt,:,:,:] = np.diff(z_w_sorted,axis=2)
    dz[tt,:,:,:] = dz_inv[tt,:,:,::-1]
print(dz.shape, temp_hdiff.shape)

# transpose dz to the same dim order as temp_xdiff
dz_trsp = dz.transpose(0,3,1,2)

print(dz_trsp.shape)

# obtain thermal expansion (alpha) & salinity contraction (beta) coefficients:
SA = np.empty(salt.shape)
# neet Absolute Salinity, converting from Pratical Salinity:
print('salt and z_rho shape:', np.squeeze(salt[0,0,:,:]).shape,np.squeeze(z_rho[0,:,:,0].shape))
for mm in np.arange(0,12):
    for kk in np.arange(0,31):
        SA_tmp = gsw.SA_from_SP(np.squeeze(salt[mm,kk,:,:]),np.squeeze(z_rho[mm,:,:,kk]),lon_rho,lat_rho)
        SA[mm,kk,:,:] = SA_tmp
        del SA_tmp

[specvol, alpha, beta] = gsw.specvol_alpha_beta(SA,temp,z_rho.transpose(0,3,1,2))

print('Alpha/beta shapes:', alpha.shape, beta.shape)

# multiply heat/salt budget diagnostics by heat expansion/salinity contraction coeffs.:
Atemp_hdiff = temp_hdiff*alpha
Atemp_vdiff = temp_vdiff*alpha
Atemp_hadv = temp_hadv*alpha
Atemp_vadv = temp_vadv*alpha
Atemp_rate = temp_rate*alpha
Bsalt_hdiff = salt_hdiff*beta
Bsalt_vdiff = salt_vdiff*beta
Bsalt_hadv = salt_hadv*beta
Bsalt_vadv = salt_vadv*beta
Bsalt_rate = salt_rate*beta

# allocate depth-integrated variables:
temp_hdiff_full_vint = np.empty(salt[:,0,:,:].shape)
temp_vdiff_full_vint = np.empty(salt[:,0,:,:].shape)
salt_hdiff_full_vint = np.empty(salt[:,0,:,:].shape)
salt_vdiff_full_vint = np.empty(salt[:,0,:,:].shape)

temp_hadv_full_vint = np.empty(salt[:,0,:,:].shape)
temp_vadv_full_vint = np.empty(salt[:,0,:,:].shape)
salt_hadv_full_vint = np.empty(salt[:,0,:,:].shape)
salt_vadv_full_vint = np.empty(salt[:,0,:,:].shape)

temp_rate_full_vint = np.empty(salt[:,0,:,:].shape)
salt_rate_full_vint = np.empty(salt[:,0,:,:].shape)

for mm in np.arange(0,12):
    depth = np.squeeze(z_rho[mm,:,:,:])
    print(mm)
    
    for ii in np.arange(0,560):
        for jj in np.arange(0,630):
            
           # for vertically-integrated vars: calculate var*dz_trsp
            # Diffusion
            temp_hdiff_full_dz_trsp = Atemp_hdiff[mm,0:,ii,jj]*dz_trsp[mm,0:,ii,jj]
            temp_vdiff_full_dz_trsp = Atemp_vdiff[mm,0:,ii,jj]*dz_trsp[mm,0:,ii,jj]
            salt_hdiff_full_dz_trsp = Bsalt_hdiff[mm,0:,ii,jj]*dz_trsp[mm,0:,ii,jj]
            salt_vdiff_full_dz_trsp = Bsalt_vdiff[mm,0:,ii,jj]*dz_trsp[mm,0:,ii,jj]
            #print(temp_hdiff_full_dz_trsp.shape)
            # sum in the vertical
            temp_hdiff_full_vint[mm,ii,jj] = np.nansum(temp_hdiff_full_dz_trsp, axis=0)
            temp_vdiff_full_vint[mm,ii,jj] = np.nansum(temp_vdiff_full_dz_trsp, axis=0)
            salt_hdiff_full_vint[mm,ii,jj] = np.nansum(salt_hdiff_full_dz_trsp, axis=0)
            salt_vdiff_full_vint[mm,ii,jj] = np.nansum(salt_vdiff_full_dz_trsp, axis=0)
            del temp_hdiff_full_dz_trsp, temp_vdiff_full_dz_trsp
            del salt_hdiff_full_dz_trsp, salt_vdiff_full_dz_trsp

            # Advection
            temp_hadv_full_dz_trsp = Atemp_hadv[mm,0:,ii,jj]*dz_trsp[mm,0:,ii,jj]
            temp_vadv_full_dz_trsp = Atemp_vadv[mm,0:,ii,jj]*dz_trsp[mm,0:,ii,jj]
            salt_hadv_full_dz_trsp = Bsalt_hadv[mm,0:,ii,jj]*dz_trsp[mm,0:,ii,jj]
            salt_vadv_full_dz_trsp = Bsalt_vadv[mm,0:,ii,jj]*dz_trsp[mm,0:,ii,jj]
            # sum in the vertical
            temp_hadv_full_vint[mm,ii,jj] = np.nansum(temp_hadv_full_dz_trsp, axis=0)
            temp_vadv_full_vint[mm,ii,jj] = np.nansum(temp_vadv_full_dz_trsp, axis=0)
            salt_hadv_full_vint[mm,ii,jj] = np.nansum(salt_hadv_full_dz_trsp, axis=0)
            salt_vadv_full_vint[mm,ii,jj] = np.nansum(salt_vadv_full_dz_trsp, axis=0)
            del temp_hadv_full_dz_trsp, temp_vadv_full_dz_trsp
            del salt_hadv_full_dz_trsp, salt_vadv_full_dz_trsp
                
            # Tendency
            temp_rate_full_dz_trsp = Atemp_rate[mm,0:,ii,jj]*dz_trsp[mm,0:,ii,jj]
            salt_rate_full_dz_trsp = Bsalt_rate[mm,0:,ii,jj]*dz_trsp[mm,0:,ii,jj]
            # sum in the vertical
            temp_rate_full_vint[mm,ii,jj] = np.nansum(temp_rate_full_dz_trsp, axis=0)
            salt_rate_full_vint[mm,ii,jj] = np.nansum(salt_rate_full_dz_trsp, axis=0)
            del temp_rate_full_dz_trsp, salt_rate_full_dz_trsp
print(temp_hdiff_full_vint.shape)
# fix mask crazy values in temp_hdiff var (not sure what cause it, iteractive run doesn't happen)
for mm in np.arange(0,12):
    temp_hdiff_full_vint[mm,:,:] = temp_hdiff_full_vint[mm,:,:]*dg.mask_rho

mm=8

fig = plt.figure(figsize=(16,14))
ax1 = fig.add_subplot(321)
plt.title('hdiff vint')
cy=plt.pcolor(temp_hdiff_full_vint[mm,:,:])
plt.colorbar(cy)
plt.clim(-.0001,.0001)
ax2 = fig.add_subplot(322)
plt.title('vdiff vint')
cy=plt.pcolor(temp_vdiff_full_vint[mm,:,:])
plt.colorbar(cy)
plt.clim(-.0001,.0001)
ax3 = fig.add_subplot(323)
plt.title('hadv vint')
cy=plt.pcolor(temp_hadv_full_vint[mm,:,:])
plt.colorbar(cy)
plt.clim(-.0001,.0001)
ax4 = fig.add_subplot(324)
plt.title('vadv vint')
cy=plt.pcolor(temp_vadv_full_vint[mm,:,:])
plt.colorbar(cy)
plt.clim(-.0001,.0001)
ax5 = fig.add_subplot(325)
plt.title('tend vint')
cy=plt.pcolor(temp_rate_full_vint[mm,:,:])
plt.colorbar(cy)
plt.clim(-.0001,.0001)
ax6 = fig.add_subplot(326)
rhs = temp_hdiff_full_vint[mm,:,:]+temp_vdiff_full_vint[mm,:,:]+temp_hadv_full_vint[mm,:,:]+temp_vadv_full_vint[mm,:,:]
plt.title('residue')
cy=plt.pcolor(temp_rate_full_vint[mm,:,:]-rhs)
plt.colorbar(cy)
#plt.clim(-.0001,.0001)

fig_path = '/users/boeiradi/COLD_project/postprocessing/figs/WMT/'
name_fig='waom10extend_shflim_S_0.25Q_Full_OHB_Sep_yr20.png'
plt.savefig(fig_path + name_fig, dpi=300)
plt.close()

# export netCDF file containing:
# - MLD (Hsbl)
# - temp and salt averaged over MLD [temp_MLD, salt_MLD] 
# - temp/salt diffusion components (xz and yz, across MLD)

fn = '/scratch/project_2000339/boeiradi/waom10extend_shflim_S_0.25Q/output_20yr_diag/Full_vint_vars_for_WMT_m.s-1.nc'
dx = Dataset(fn, 'w', format='NETCDF4')

times = dx.createDimension('times', 12)
xi_rho = dx.createDimension('xi_rho', 630,)
eta_rho = dx.createDimension('eta_rho', 560)

ocean_times = dx.createVariable('times', 'f4', ('times',))
ocean_times.units = 'seconds of the year'
xi = dx.createVariable('xi_rho', 'f4', ('xi_rho',))
eta = dx.createVariable('eta_rho', 'f4', ('eta_rho',))

# vertically integrated terms
# diffusion terms
temp_vdia_diff_full_vint = dx.createVariable('temp_vdia_diff_full_vint', 'f4', ('times', 'eta_rho', 'xi_rho',))
salt_vdia_diff_full_vint = dx.createVariable('salt_vdia_diff_full_vint', 'f4', ('times', 'eta_rho', 'xi_rho',))
temp_vdia_diff_full_vint.units = 'vertical Full-depth component  of temp_diff full depth-integrated (degC s-1)'
salt_vdia_diff_full_vint.units = 'vertical Full-depth component  of salt_diff full depth-integrated (s-1)'
temp_hdia_diff_full_vint = dx.createVariable('temp_hdia_diff_full_vint', 'f4', ('times', 'eta_rho', 'xi_rho',))
salt_hdia_diff_full_vint = dx.createVariable('salt_hdia_diff_full_vint', 'f4', ('times', 'eta_rho', 'xi_rho',))
temp_hdia_diff_full_vint.units = 'horizontal Full-depth component of temp_diff full depth-integrated (degC s-1)'
salt_hdia_diff_full_vint.units = 'horizontal Full-depth component of salt_diff full depth-integrated (s-1)'
# advective terms
temp_vdia_adv_full_vint = dx.createVariable('temp_vdia_adv_full_vint', 'f4', ('times', 'eta_rho', 'xi_rho',))
salt_vdia_adv_full_vint = dx.createVariable('salt_vdia_adv_full_vint', 'f4', ('times', 'eta_rho', 'xi_rho',))
temp_vdia_adv_full_vint.units = 'vertical Full-depth component  of temp_adv to the base of the MLD (degC s-1)'
salt_vdia_adv_full_vint.units = 'vertical Full-depth component  of salt_adv to the base of the MLD (s-1)'
temp_hdia_adv_full_vint = dx.createVariable('temp_hdia_adv_full_vint', 'f4', ('times', 'eta_rho', 'xi_rho',))
salt_hdia_adv_full_vint = dx.createVariable('salt_hdia_adv_full_vint', 'f4', ('times', 'eta_rho', 'xi_rho',))
temp_hdia_adv_full_vint.units = 'horizontal Full-depth component of temp_adv to the base of the MLD (degC s-1)'
salt_hdia_adv_full_vint.units = 'horizontal Full-depth component of salt_adv to the base of the MLD (s-1)'
# tendency terms
temp_tend_full_vint = dx.createVariable('temp_tend_avg_full_vint', 'f4', ('times', 'eta_rho', 'xi_rho',))
salt_tend_full_vint = dx.createVariable('salt_tend_avg_full_vint', 'f4', ('times', 'eta_rho', 'xi_rho',))
temp_tend_full_vint.units = 'temp_rate full depth-integrated (degC s-1)'
salt_tend_full_vint.units = 'salt_rate full depth-integrated (s-1)'

ocean_times[:] = np.arange(1314000,31536000,2628000) # monthly 
xi[:] = np.arange(0,630)
eta[:] = np.arange(0,560)

temp_vdia_diff_full_vint[:,:,:] = temp_vdiff_full_vint
salt_vdia_diff_full_vint[:,:,:] = salt_vdiff_full_vint
temp_hdia_diff_full_vint[:,:,:] = temp_hdiff_full_vint
salt_hdia_diff_full_vint[:,:,:] = salt_hdiff_full_vint

temp_vdia_adv_full_vint[:,:,:] = temp_vadv_full_vint
salt_vdia_adv_full_vint[:,:,:] = salt_vadv_full_vint
temp_hdia_adv_full_vint[:,:,:] = temp_hadv_full_vint
salt_hdia_adv_full_vint[:,:,:] = salt_hadv_full_vint

temp_tend_full_vint[:,:,:] = temp_rate_full_vint
salt_tend_full_vint[:,:,:] = salt_rate_full_vint

dx.close()
