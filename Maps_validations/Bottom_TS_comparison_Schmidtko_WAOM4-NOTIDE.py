
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
from matplotlib import colors
import matplotlib.path as mpath

from datetime import datetime, timedelta

from netCDF4 import Dataset
from netCDF4 import num2date, date2num

import os

import gsw
import pyresample
from xgcm import Grid

# read grid file for lon/lat coordinates
dg = xr.open_dataset("/g/data3/hh5/tmp/access-om/fbd581/ROMS/waom4_frc/waom4extend_grd.nc")
lat_rho_4km= dg.variables["lat_rho"]
lon_rho_4km = dg.variables["lon_rho"]
lat_u_4km= dg.variables["lat_u"]
lon_u_4km = dg.variables["lon_u"]
lat_v_4km= dg.variables["lat_v"]
lon_v_4km = dg.variables["lon_v"]
cor_4km = dg.variables["f"]
pm_4km = dg.variables["pm"]
pn_4km = dg.variables["pn"]
zice_4km = dg.variables["zice"]
h_4km = dg.variables["h"]
dg.close()
print('Print lon/lat_rho shapes',lon_rho_4km.shape, lat_rho_4km.shape)
print('Print lon/lat_rho shapes',lon_rho_4km[0:-1,0:-1].shape, lat_rho_4km[0:-1,0:-1].shape)

#### load ROMS avg output

def read_roms_ts_4km(exp_path,year):
    for yr  in [year]:
        ds = xr.open_dataset(exp_path + 'ocean_avg_00' + yr + '.nc')
        print('size temp and time length: ', ds.temp.shape, len(ds.salt.isel(xi_rho=20, eta_rho=100, s_rho=0)))
        temp_tmp = ds.variables["temp"]
        salt_tmp = ds.variables["salt"]
        zeta_tmp = ds.variables["zeta"]
        melt_tmp = ds.variables["m"]
        temp_tmp_ann = np.nanmean(temp_tmp, axis=0)
        salt_tmp_ann = np.nanmean(salt_tmp, axis=0)
        melt_tmp_ann = np.nanmean(melt_tmp, axis=0)

        print('size temp_tmp_ann = ', temp_tmp_ann.shape)

        ds = ds.set_coords(['Cs_r', 'Cs_w', 'hc', 'h', 'Vtransform'])

        if ds.Vtransform == 1:
            Zo_rho = ds.hc * (ds.s_rho - ds.Cs_r) + ds.Cs_r * ds.h
            z_rho_tmp = Zo_rho + ds.zeta * (1 + Zo_rho/ds.h)
            print("Vtransform=1")
        elif ds.Vtransform == 2:
            Zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * ds.h) / (ds.hc + ds.h)
            z_rho_tmp = ds.zeta + (ds.zeta + ds.h) * Zo_rho
            print("Vtransform=2")
        z_rho_tmp_ann = np.nanmean(z_rho_tmp, axis=0)

        Zo_w = (ds.hc * ds.s_w + ds.Cs_w * ds.h) / (ds.hc + ds.h)
        z_w_tmp = ds.zeta + (ds.zeta + ds.h) * Zo_w
        z_w_tmp_ann = np.nanmean(z_w_tmp, axis=0)

        # Handle interpolation from u and v grid to rho points:
        ds = ds.rename({'eta_u': 'eta_rho', 'xi_v': 'xi_rho', 'xi_psi': 'xi_u', 'eta_psi': 'eta_v'})

        coords={'X':{'center':'xi_rho', 'inner':'xi_u'},
            'Y':{'center':'eta_rho', 'inner':'eta_v'},
            'Z':{'center':'s_rho', 'outer':'s_w'}}

        grid = Grid(ds, coords=coords, periodic=[])

        Zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * ds.h) / (ds.hc + ds.h)
        z_rho = Zo_rho * (ds.zeta + ds.h) + ds.zeta + ds.zice
        Zo_w = (ds.hc * ds.s_w + ds.Cs_w * ds.h) / (ds.hc + ds.h)
        z_w = Zo_w * (ds.zeta + ds.h) + ds.zeta + ds.zice

        ds.coords['z_w'] = z_w.where(ds.mask_rho, 0).transpose('ocean_time', 's_w', 'eta_rho', 'xi_rho')
        ds.coords['z_rho'] = z_rho.where(ds.mask_rho, 0).transpose('ocean_time', 's_rho', 'eta_rho', 'xi_rho')

        ds['pm_v'] = grid.interp(ds.pm, 'Y')
        ds['pn_u'] = grid.interp(ds.pn, 'X')
        ds['pm_u'] = grid.interp(ds.pm, 'X')
        ds['pn_v'] = grid.interp(ds.pn, 'Y')
        ds['pm_psi'] = grid.interp(grid.interp(ds.pm, 'Y'),  'X') # at psi points (eta_v, xi_u)
        ds['pn_psi'] = grid.interp(grid.interp(ds.pn, 'X'),  'Y') # at psi points (eta_v, xi_u)

        ds['dx'] = 1/ds.pm
        ds['dx_u'] = 1/ds.pm_u
        ds['dx_v'] = 1/ds.pm_v
        ds['dx_psi'] = 1/ds.pm_psi

        ds['dy'] = 1/ds.pn
        ds['dy_u'] = 1/ds.pn_u
        ds['dy_v'] = 1/ds.pn_v
        ds['dy_psi'] = 1/ds.pn_psi

        ds['dz'] = grid.diff(ds.z_w, 'Z', boundary='fill')
        ds['dz_w'] = grid.diff(ds.z_rho, 'Z', boundary='fill')
        ds['dz_u'] = grid.interp(ds.dz, 'X')
        ds['dz_w_u'] = grid.interp(ds.dz_w, 'X')
        ds['dz_v'] = grid.interp(ds.dz, 'Y')
        ds['dz_w_v'] = grid.interp(ds.dz_w, 'Y')

        ds['dA'] = ds.dx * ds.dy

        metrics = {
            ('X',): ['dx', 'dx_u', 'dx_v', 'dx_psi'], # X distances
            ('Y',): ['dy', 'dy_u', 'dy_v', 'dy_psi'], # Y distances
            ('Z',): ['dz', 'dz_u', 'dz_v', 'dz_w', 'dz_w_u', 'dz_w_v'], # Z distances
            ('X', 'Y'): ['dA'] # Areas
        }
        grid = Grid(ds, coords=coords, metrics=metrics, periodic=False)

        u4_rho_sfc = np.zeros((12,1400,1575))
        v4_rho_sfc = np.zeros((12,1400,1575))
        u4_rho_bot = np.zeros((12,1400,1575))
        v4_rho_bot = np.zeros((12,1400,1575))

        for mm in np.arange(0,12):
            #interpoate u, v to rho grid:
            u4_interp = grid.interp(ds.u.isel(s_rho=0,ocean_time=mm), 'X',boundary='fill')
            v4_interp = grid.interp(ds.v.isel(s_rho=0,ocean_time=mm), 'Y',boundary='fill')
            u4_rho_bot[mm,:,:]=u4_interp
            v4_rho_bot[mm,:,:]=v4_interp
            del u4_interp,v4_interp
            u4_interp = grid.interp(ds.u.isel(s_rho=-1,ocean_time=mm), 'X',boundary='fill')
            v4_interp = grid.interp(ds.v.isel(s_rho=-1,ocean_time=mm), 'Y',boundary='fill')
            u4_rho_sfc[mm,:,:]=u4_interp
            v4_rho_sfc[mm,:,:]=v4_interp
            del u4_interp,v4_interp

        u4_rho_bot_ann = np.nanmean(u4_rho_bot, axis=0)
        v4_rho_bot_ann = np.nanmean(v4_rho_bot, axis=0)
        u4_rho_sfc_ann = np.nanmean(u4_rho_sfc, axis=0)
        v4_rho_sfc_ann = np.nanmean(v4_rho_sfc, axis=0)
        # concantenate annual averaged temp/salt
        if yr == '10':
            temp_ann = temp_tmp_ann
            salt_ann = salt_tmp_ann
            melt_ann = melt_tmp_ann
            u4_bot = u4_rho_bot_ann
            v4_bot = v4_rho_bot_ann
            u4_sfc = u4_rho_sfc_ann
            v4_sfc = v4_rho_sfc_ann

    return temp_ann, salt_ann, u4_sfc, v4_sfc, u4_bot, v4_bot, melt_ann


path_ECCO2_4km = '/g/data3/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom4extend_shflim_S_0.25Q/output_01-20yr/'
path_ECCO2_4kmNT = '/g/data3/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom4extend_shflim_S_0.25Q/output_01-10yr_notides/'


temp_ann_4km, salt_ann_4km, u4_sfc, v4_sfc, u4_bot, v4_bot, melt_ann_4km = read_roms_ts_4km(path_ECCO2_4km,'10')

temp_ann_4kmNT, salt_ann_4kmNT, u4NT_sfc, v4NT_sfc, u4_bot, v10_bot, melt_ann_4kmNT = read_roms_ts_4km(path_ECCO2_4kmNT,'10')

mask_zice_4km = ma.masked_where(zice_4km < 0, np.ones(zice_4km.shape))
mask_outice_4km = ma.masked_where(zice_4km >= 0, np.ones(zice_4km.shape))

mask_shelf_4km = ma.masked_where(h_4km > 2000, np.ones(zice_4km.shape))

fig_path='/g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/figs/Maps/'

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


# call cartopy projection
proj = ccrs.SouthPolarStereo()

# ========== subregions plots:
# for 2000m isoline
bathym = cfeature.NaturalEarthFeature(name='bathymetry_J_1000', scale='10m', category='physical')
# limits for contour of ice front (Ronne-Filchner IS):
# xlimit = np.arange(300,500,1)
# ylimit = np.arange(100,300,1)

# --- bottom

tmin = -2.7
tmax = 1.5
smin = 34.1
smax = 35.

fig = plt.figure(figsize=(12,8))

ax1 = fig.add_subplot(221, projection=proj)
ct1=plt.pcolormesh(lon_rho_4km,lat_rho_4km,np.squeeze(temp_ann_4km[0,:,:])*mask_shelf_4km, transform=ccrs.PlateCarree(), cmap=plt.cm.plasma, vmin=tmin, vmax=tmax)
# plt.title('Bottom temperature \n WAOM10')
plt.title('WAOM4')
ax1.gridlines() # draw_labels=True,linewidth=
ax1.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')

ax2 = fig.add_subplot(222, projection=proj)
cs1=plt.pcolormesh(lon_rho_4km,lat_rho_4km,np.squeeze(salt_ann_4km[0,:,:])*mask_shelf_4km, transform=ccrs.PlateCarree(), cmap=plt.cm.plasma, vmin=smin, vmax=smax)
# plt.title('Bottom salinity \n WAOM10')
plt.title('WAOM4')
ax2.gridlines()
ax2.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')

ax3 = fig.add_subplot(223, projection=proj)
ct1a=plt.pcolormesh(lon_rho_4km,lat_rho_4km,(np.squeeze(temp_ann_4kmNT[0,:,:])-np.squeeze(temp_ann_4km[0,:,:]))*mask_shelf_4km, transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=-0.6, vmax=0.6)
# plt.title('Bottom temperature \n WAOM10')
plt.title('WAOM4-NOTIDE minus WAOM4')
ax3.gridlines() # draw_labels=True,linewidth=
# ax1.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
ratio = .9
ax3.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')

ax4 = fig.add_subplot(224, projection=proj)
cs1a=plt.pcolormesh(lon_rho_4km,lat_rho_4km,(np.squeeze(salt_ann_4kmNT[0,:,:])-np.squeeze(salt_ann_4km[0,:,:]))*mask_shelf_4km, transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=-.25, vmax=.25)
# plt.title('Bottom salinity \n WAOM10')
plt.title('WAOM4-NOTIDE minus WAOM4')
ax4.gridlines()
ax4.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')

cbar_ax1 = fig.add_axes([0.08, 0.55, 0.01, 0.35])
fig.colorbar(ct1, cax=cbar_ax1, orientation='vertical')
cbar_ax1.set_ylabel('$^{\circ}$C')#, labelpad=-35)

cbar_ax2 = fig.add_axes([0.9, 0.55, 0.01, 0.35])
fig.colorbar(cs1, cax=cbar_ax2, orientation='vertical')
cbar_ax2.set_xlabel('')#, labelpad=-35)

cbar_ax3 = fig.add_axes([0.08, 0.1, 0.01, 0.35])
fig.colorbar(ct1a, cax=cbar_ax3, orientation='vertical')
cbar_ax3.set_ylabel('$^{\circ}$C')#, labelpad=-35)

cbar_ax4 = fig.add_axes([0.9, 0.1, 0.01, 0.35])
fig.colorbar(cs1a, cax=cbar_ax4, orientation='vertical')
cbar_ax4.set_xlabel('')#, labelpad=-35)

name_fig='waom4x4-NOTIDE_bottomTSmaps_yr10.png'
plt.savefig(fig_path + name_fig, dpi=300)
plt.close()
