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

from xhistogram.xarray import histogram
from xmovie import Movie
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
    Hsbl_tmp = np.nanmean(ds.variables["Hsbl"], axis=0)

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
        Hsbl = Hsbl_tmp
    elif mm == '02':
        temp = np.stack((temp,temp_tmp), axis=0)
        salt = np.stack((salt,salt_tmp), axis=0)
        shflux = np.stack((shflux,shflux_tmp), axis=0)
        ssflux = np.stack((ssflux,ssflux_tmp), axis=0)
        m = np.stack((m,m_tmp), axis=0)
        z_rho = np.stack((z_rho,z_rho_avg), axis=0)
        z_w = np.stack((z_w,z_w_avg), axis=0)
        Hsbl = np.stack((Hsbl,Hsbl_tmp), axis=0)
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
        Hsbl_tmp_4thdim = np.expand_dims(Hsbl_tmp, axis=0)
        Hsbl = np.concatenate((Hsbl,Hsbl_tmp_4thdim), axis=0)
#ds.coords['flux'] = flux#.transpose() # put flux into ds dataset

    ds.close()

sigma_t_sfc = gsw.rho(salt[:,-1,:,:],temp[:,-1,:,:],0) - 1000


di = xr.open_dataset('/scratch/project_2000789/boeiradi/waom4extend_shflim_S_0.25Q/output_yr10_diag/ocean_avg_0001.nc')
ice_draft = di.variables["zice"]
h = di.variables["h"]

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

dg = xr.open_dataset("/scratch/project_2000789/boeiradi/waom4_frc/waom4extend_grd.nc")

lat_rho = dg.variables["lat_rho"]
lon_rho = dg.variables["lon_rho"]
pm = dg.variables["pm"]
pn = dg.variables["pn"]

ds.coords['lat_rho']=lat_rho.transpose() # put lat_rho into ds dataset
ds.coords['lon_rho']=lon_rho.transpose() # put lon_rho into ds dataset

# shelf/open-ocean masks:
mask_open = ma.masked_where(dg.h <= 2000, np.ones(dg.h.shape))
mask_shelf = ma.masked_where(dg.h > 2000, np.ones(dg.h.shape))


dx = xr.open_dataset('/scratch/project_2000789/boeiradi/waom4extend_shflim_S_0.25Q/output_yr10_diag/MLD_vint_vars_for_WMT_m.s-1_iteractive_copy.nc')
print(dx.variables["sfc_mld"].shape)
sfc_mld = dx.variables["sfc_mld"]
# - variables across ML base

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

sigma_t = gsw.rho(salt[:,-1,:,:],temp[:,-1,:,:],0) - 1000

dx.close()

# calculate salt/temp averaged for ML and below:

salt_mld = np.empty(salt[:,0,:,:].shape)
temp_mld = np.empty(salt[:,0,:,:].shape)

for mm in np.arange(0,12):
    depth = np.squeeze(z_rho[mm,:,:,:])
    depth_sort = depth.transpose(2,0,1) #'s_rho','eta_rho','xi_rho')
#    print(depth_sort.shape)
    print(mm)

    salt_mld_tmp = ma.masked_where(-depth_sort > -np.squeeze(Hsbl[mm,:,:]), salt[mm,:,:,:])
    temp_mld_tmp = ma.masked_where(-depth_sort > -np.squeeze(Hsbl[mm,:,:]), temp[mm,:,:,:])
    print(salt_mld_tmp.shape)

    # put SSS/SST into salt/temp_mld var independently to get rid of masked values when Hsbl < first layer of depth_sort
    salt_mld_tmp[-1,:,:] = salt[mm,-1,:,:]
    temp_mld_tmp[-1,:,:] = temp[mm,-1,:,:]

    salt_mld[mm,:,:] = np.nanmean(salt_mld_tmp, axis=0)
    temp_mld[mm,:,:] = np.nanmean(temp_mld_tmp, axis=0)
    del depth, depth_sort #, salt_mld_tmp, temp_mld_tmp, salt_deep_tmp, temp_deep_tmp

# calculate WMT from budget vars:

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
rho0 = 1025
Cp = 3985

# total diffusion terms:
R_s_vint = (salt_hdia_diff_mld_vint + salt_vdia_diff_mld_vint)
R_t_vint = (temp_hdia_diff_mld_vint + temp_vdia_diff_mld_vint)

# surface flux terms:
salt_sfc = beta[:,-1,:,:]*(ssflux)
temp_sfc = alpha[:,-1,:,:]*(np.divide(shflux, rho0*Cp))

# advection terms:
salt_adv_mld_vint = (salt_hdia_adv_mld_vint + salt_vdia_adv_mld_vint)
temp_adv_mld_vint = (temp_hdia_adv_mld_vint + temp_vdia_adv_mld_vint)

# net tendencies
salt_net_mld_vint = salt_tend_mld_vint
temp_net_mld_vint = temp_tend_mld_vint

# gather sections only:
dx = np.divide(1,pm)
dy = np.divide(1,pn)
dt = 86400#30#/12 #why divide by 12?

# new section, too far east
# West Weddell = (xi_rho=438, eta_rho=slice(320,389))
# East Weddell = (xi_rho=500, eta_rho=slice(330,394))
# SO, new west one::
# West Weddell = (xi_rho=450, eta_rho=slice(320,389))
# East Weddell = (xi_rho=500, eta_rho=slice(330,394))

salt_net_WWed = salt_net_mld_vint.isel(xi_rho=slice(413,437), eta_rho=slice(800,973))
salt_adv_WWed = salt_adv_mld_vint.isel(xi_rho=slice(413,437), eta_rho=slice(800,973))
salt_diff_WWed = R_s_vint.isel(xi_rho=slice(413,437), eta_rho=slice(800,973))
salt_sfc_WWed = salt_sfc[:,800:973,413:437]
salt_mld_WWed = salt_mld[:,800:973,413:437]
temp_mld_WWed = temp_mld[:,800:973,413:437]
dx_WWed = dx[800:973,413:437]
dy_WWed = dy[800:973,413:437]

salt_net_EWed = salt_net_mld_vint.isel(xi_rho=slice(488,512), eta_rho=slice(825,985))
salt_adv_EWed = salt_adv_mld_vint.isel(xi_rho=slice(488,512), eta_rho=slice(825,985))
salt_diff_EWed = R_s_vint.isel(xi_rho=slice(488,512), eta_rho=slice(825,985))
salt_sfc_EWed = salt_sfc[:,825:985,488:512]
salt_mld_EWed = salt_mld[:,825:985,488:512]
temp_mld_EWed = temp_mld[:,825:985,488:512]
dx_EWed = dx[825:985,488:512]
dy_EWed = dy[825:985,488:512]

print(Hsbl.shape, z_rho.shape)
# (12, 560, 630) (12, 31, 560, 630)

# fix z_rho dimensions order
z_rho = z_rho.transpose(0,3,1,2)
print(z_rho.shape)


# 2nd attempt TS-space; combined T and S


sbins = np.arange(31,35, 0.01) # aim for S interval = 0.025
tbins = np.arange(-3, 8, 0.05)+273.15 # aim for T interval = 0.1 in Kelvin!
len_tbin = len(tbins)
len_sbin = len(sbins)

def wmt_TS(var_int, temp, salt, dx, dy,var_type):
    # var_type: 'budget' or 'sfc_frc'

    F_rate_var_vint = np.empty(var_int.shape)

    for mm in np.arange(0,12):
        if var_type == 'budget':
            F_rate_var_vint[mm,:] = dx*dy*var_int.isel(times=mm)
        elif var_type == 'sfc_frc':
            F_rate_var_vint[mm,:] = dx*dy*var_int[mm,:]

    print(F_rate_var_vint.shape)

    F_TS_rate_delta_var_vint_mm = np.empty((12,len_tbin,len_sbin,len(var_int[0,:]),24))
#     F_TS_rate_delta_var_vint_mm = np.empty((1,len_tbin,len_sbin,len(var_int[0,:])))
    
    for mm in np.arange(0,12):
        temp_tmp = temp[mm,:]+273.15 # convert to Kelvin
        salt_tmp = salt[mm,:]
    
        F_rate_delta = np.empty((len_tbin,len_sbin,len(var_int[0,:]),24))

        #print(mm)
        for itemp in np.arange(0,len_tbin):
            for isalt in np.arange(0,len_sbin):
                # create mask=True for values outside Tbin/Sbin, and F_rate_var_vint inside Tbin/Sbin
                F_rate_tmp = ma.masked_where( \
                    np.logical_or( \
                    np.logical_or(temp_tmp <= (tbins[itemp]-(0.05/2)),temp_tmp > (tbins[itemp]+(0.05/2))), \
                    np.logical_or(salt_tmp <= (sbins[isalt]-(0.01/2)),salt_tmp > (sbins[isalt]+(0.01/2)))), \
                    F_rate_var_vint[mm,:])
                # define NaN for values outside Tbin/Sbin
                F_rate_tmp[np.logical_or( \
                    np.logical_or(temp_tmp <= (tbins[itemp]-(0.05/2)),temp_tmp > (tbins[itemp]+(0.05/2))), \
                    np.logical_or(salt_tmp <= (sbins[isalt]-(0.01/2)),salt_tmp > (sbins[isalt]+(0.01/2))))] = np.nan

                F_rate_delta[itemp,isalt,:] = F_rate_tmp # allocate xi_rho, eta_rho array

                del F_rate_tmp
        print(F_rate_delta.shape)

        F_TS_rate_delta_var_vint_mm[mm,:] = F_rate_delta

    print('completed, size: ', F_TS_rate_delta_var_vint_mm.shape)

    return F_TS_rate_delta_var_vint_mm

print(len(salt_net_WWed[0,:]), dx_WWed.shape)

# TS-space
F_TS_rate_delta_net_WWed_mm = wmt_TS(salt_net_WWed, temp_mld_WWed, salt_mld_WWed, dx_WWed, dy_WWed,'budget')

F_TS_rate_delta_net_EWed_mm = wmt_TS(salt_net_EWed, temp_mld_EWed, salt_mld_EWed, dx_EWed, dy_EWed,'budget')

F_TS_rate_delta_net_WWed_mm.shape

# TS-space WMT: integrating horizontally

# TS-space
Fs_TS_rate_net_WWed_mm_int = np.nansum(np.nansum(F_TS_rate_delta_net_WWed_mm, axis=4), axis=3)
Fs_TS_rate_net_EWed_mm_int = np.nansum(np.nansum(F_TS_rate_delta_net_EWed_mm, axis=4), axis=3)


# make grid for density contours
smin = 34 - (0.01 * 34)    #salt_ctrl_subregR.min - (0.01 * salt_ctrl_subregR.min)
smax = 35. + (0.01 * 35.)    #salt_ctrl_subregR.max + (0.01 * salt_ctrl_subregR.max)
tmin = -3. + (0.1 * -3.)       #temp_ctrl_subregR.min - (0.1 * temp_ctrl_subregR.max)
tmax = 1 + (0.1 * 1.)       #temp_ctrl_subregR.max + (0.1 * temp_ctrl_subregR.max)
print('tmin, tmax, smin, smax sizes=,', tmin, tmax, smin, smax)
# Calculate how many gridcells we need in the x and y dimensions
xdim = 40
ydim = 20
# Create empty grid of zeros
dens = np.zeros((ydim,xdim))
# Create temp and salt vectors of appropiate dimensions
ti = np.linspace(-3,1,ydim)
si = np.linspace(34,35,xdim)

Si, Ti = np.meshgrid(si, ti, sparse=False, indexing='ij')
# Loop to fill in grid with densities
for j in range(0,int(ydim)):
    for i in range(0, int(xdim)):
        dens[j,i]=gsw.rho(si[i],ti[j],0)
        # Substract 1000 to convert to sigma-0
dens = dens - 1000

# freezing point temperature:

Freez_tempW = gsw.CT_freezing(si,500,0)
Freez_tempE = gsw.CT_freezing(si,400,0)

print(Freez_tempW.shape, si.shape)

# Gade line: WDW
# Twdw = -.75
# Swdw = 34.5
Twdw = -1.
Swdw = 34.5
Lf = 334 # kJ/kg
#Cp = 3.97 #kJ/kg/Kelvin
Cp = gsw.cp_t_exact(Swdw,Twdw,0)/1000
print(Cp)

Tgade = Twdw + (Lf/Cp)*(1 - (Swdw/si))

# Gade line: WDW
Thssw = -2.2
Shssw = 34.4
Lf = 334 # kJ/kg
Cp2 = gsw.cp_t_exact(Shssw,Thssw,0)/1000
print(Cp2)

Tgade2 = Thssw + (Lf/Cp2)*(1 - (Shssw/si))

Dt = 1000/(0.05*0.01)

[Sbins,Tbins] = np.meshgrid(sbins,tbins)
print(Tbins.shape, tbins.shape)

fig = plt.figure(figsize=(12,16))

plt.subplot(3,2,3)
# cw = plt.pcolormesh(Sbins,Tbins-273.15,np.nanmean(Fs_TS_rate_net_WWed_mm_int*Dt/1e6, axis=0), vmin=-1, vmax=1, cmap=plt.cm.seismic)
cw = plt.contourf(Sbins,Tbins-273.15,np.nanmean(Fs_TS_rate_net_WWed_mm_int*Dt/1e6, axis=0), levels=np.arange(-5,5.1,.1), cmap=plt.cm.coolwarm, extend='both')
#plt.colorbar()
plt.title('C) WAOM4 - West FRIS')
ax = plt.gca()
#plt.xlabel('Salinity')
plt.ylabel('Temperature ($^{\circ}$C)')
# ax.set_xlabel('Salinity',fontsize=12)
ax.set_xlim([34,34.7])
ax.set_ylim([-2.65,-.5])
CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(26.,28.2,.05),linestyles='dotted', colors='gray',linewidths=.5)
plt.clabel(CS, fontsize=11, inline=1, fmt='%1.2f')
plt.plot(si,Freez_tempE,'-',color='k', linewidth=2)
plt.plot(si,Tgade2,'-',color='green')

plt.subplot(3,2,4)
#cw = plt.pcolormesh(Sbins,Tbins-273.15,np.nanmean(Fs_TS_rate_net_EWed_mm_int*Dt/1e6, axis=0), vmin=-1, vmax=1, cmap=plt.cm.seismic)
cw = plt.contourf(Sbins,Tbins-273.15,np.nanmean(Fs_TS_rate_net_EWed_mm_int*Dt/1e6, axis=0), levels=np.arange(-5,5.1,.1), cmap=plt.cm.coolwarm, extend='both')
#plt.colorbar()
plt.title('D) WAOM4 - Central FRIS')
ax = plt.gca()
# ax.set_ylabel('Temperature ($^{\circ}$C)')
# ax.set_xlabel('Salinity',fontsize=12)
ax.set_xlim([34,34.7])
ax.set_ylim([-2.65,-.5])
CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(26.,28.2,.05),linestyles='dotted', colors='gray',linewidths=.5)
plt.clabel(CS, fontsize=11, inline=1, fmt='%1.2f')
plt.plot(si,Freez_tempE,'-',color='k', linewidth=2)
plt.plot(si,Tgade2,'-',color='green')

cbar_ax1 = fig.add_axes([0.92, 0.1, 0.015, 0.8])
fig.colorbar(cw, cax=cbar_ax1, orientation='vertical')
cbar_ax1.set_ylabel('Annual mean WMT (Sv)')

fig_path='/users/boeiradi/COLD_project/postprocessing/figs/WMT/'
name_fig="waom4extend_shflim_S_0.25Q_WMT_TS_RFISsections.png"
plt.savefig(fig_path + name_fig, dpi=300)

# print(Fs_TS_rate_net_shelf_mm_int*Dt/1e6)

# Movie:

Fs_TS_rate_net_WWed_mov = xr.DataArray(Fs_TS_rate_net_WWed_mm_int, name='Fs_salt_net', coords= {'time':np.arange(0,12),'T-bins':tbins,'S-bins':sbins})
Fs_TS_rate_net_EWed_mov = xr.DataArray(Fs_TS_rate_net_EWed_mm_int, name='Fs_salt_net', coords= {'time':np.arange(0,12),'T-bins':tbins,'S-bins':sbins})

Fs_TS_rate_net_mov = xr.Dataset(
    {
        'salt_WWed': (['time','eta_rho','xi_rho'], Fs_TS_rate_net_WWed_mm_int),
        'salt_EWed': (['time','eta_rho','xi_rho'], Fs_TS_rate_net_EWed_mm_int),

    },
    coords= {'time':np.arange(0,12),'T-bins':tbins,'S-bins':sbins}
    )

fig = plt.figure(figsize=(12,8))

def custom_plotfunc(salt, fig, tt, *args, **kwargs):
    
    plt.subplot(1,2,1)
    cw = plt.contourf(Sbins,Tbins-273.15,salt['salt_WWed'].isel(time=tt)*Dt/1e6, levels=np.arange(-15,15.1,.1), cmap=plt.cm.seismic, extend='both')
    plt.title('A) WAOM4 - West FRIS \n month = ' + str (tt))
    plt.xlabel('Salinity'),plt.ylabel('Temperature ($^{\circ}$C)')
    CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(26.,28.2,.05),linestyles='dotted', colors='gray',linewidths=.5)
    plt.clabel(CS,fontsize=8)
    plt.xlim([34,34.7])
    plt.ylim([-2.65,-.5])
    plt.plot(si,Freez_tempE,'-',color='k', linewidth=2)
    plt.plot(si,Tgade2,'-',color='green')

    plt.subplot(1,2,2)
    cw = plt.contourf(Sbins,Tbins-273.15,salt['salt_EWed'].isel(time=tt)*Dt/1e6, levels=np.arange(-15,15.1,.1), cmap=plt.cm.seismic, extend='both')
    plt.title('B) WAOM4 - Central FRIS \n month = ' + str (tt))
    plt.xlabel('Salinity'),plt.ylabel('Temperature ($^{\circ}$C)')
    CS = plt.contour(Si,Ti,dens.transpose(), levels=np.arange(26.,28.2,.05),linestyles='dotted', colors='gray',linewidths=.5)
    plt.clabel(CS,fontsize=8)
    plt.xlim([34,34.7])
    plt.ylim([-2.65,-.5])
    plt.plot(si,Freez_tempE,'-',color='k', linewidth=2)
    plt.plot(si,Tgade2,'-',color='green')

    cbar_ax1 = fig.add_axes([0.92, 0.1, 0.015, 0.8])
    fig.colorbar(cw, cax=cbar_ax1, orientation='vertical')
    cbar_ax1.set_ylabel('Annual mean WMT (Sv)')

    return None, None

mov_custom = Movie(Fs_TS_rate_net_mov, custom_plotfunc, input_check=False)
mov_custom.save('/users/boeiradi/COLD_project/postprocessing/figs/WMT/movie/WAOM4extend_WMT_TS-space_FRIS.mp4', overwrite_existing=True, progress=True, framerate=1)

