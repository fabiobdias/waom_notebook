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
from matplotlib import colors
import matplotlib.path as mpath

from datetime import datetime, timedelta

import os
import sys

from netCDF4 import Dataset
from netCDF4 import num2date, date2num

import gsw
import pyresample
from xgcm import Grid
import xesmf as xe

# read grid file for lon/lat coordinates
dg = xr.open_dataset("/scratch/project_2000339/boeiradi/waom10_frc/waom10extend_grd.nc")
lat_rho_10km= dg.variables["lat_rho"]
lon_rho_10km = dg.variables["lon_rho"]
lat_u_10km= dg.variables["lat_u"]
lon_u_10km = dg.variables["lon_u"]
lat_v_10km= dg.variables["lat_v"]
lon_v_10km = dg.variables["lon_v"]
cor_10km = dg.variables["f"]
pm_10km = dg.variables["pm"]
pn_10km = dg.variables["pn"]
zice_10km = dg.variables["zice"]
h_10km = dg.variables["h"]
dg.close()
print('Print lon/lat_rho shapes',lon_rho_10km.shape, lat_rho_10km.shape)
print('Print lon/lat_rho shapes',lon_rho_10km[0:-1,0:-1].shape, lat_rho_10km[0:-1,0:-1].shape)

dg4 = xr.open_dataset("/scratch/project_2000789/boeiradi/waom4_frc/waom4extend_grd.nc")
lat_rho_4km = dg4.variables["lat_rho"]
lon_rho_4km = dg4.variables["lon_rho"]
lat_u_4km= dg4.variables["lat_u"]
lon_u_4km = dg4.variables["lon_u"]
lat_v_4km= dg4.variables["lat_v"]
lon_v_4km = dg4.variables["lon_v"]
cor_4km = dg4.variables["f"]
pm_4km = dg4.variables["pm"]
pn_4km = dg4.variables["pn"]
zice_4km = dg4.variables["zice"]
dg4.close()

dg2 = xr.open_dataset("/scratch/project_2000339/boeiradi/waom2_frc/waom2extend_grd.nc")
lat_rho_2km = dg2.variables["lat_rho"]
lon_rho_2km = dg2.variables["lon_rho"]
lat_u_2km= dg2.variables["lat_u"]
lon_u_2km = dg2.variables["lon_u"]
lat_v_2km= dg2.variables["lat_v"]
lon_v_2km = dg2.variables["lon_v"]
cor_2km = dg2.variables["f"]
pm_2km = dg2.variables["pm"]
pn_2km = dg2.variables["pn"]
zice_2km = dg2.variables["zice"]
h_2km = dg2.variables["h"]
dg2.close()

# load ROMS avg output

def read_roms_ts_10km(exp_path):
    for yr  in ['20']:
        ds = xr.open_dataset(exp_path + 'ocean_avg_00' + yr + '.nc')
        print('size temp and time length: ', ds.temp.shape, len(ds.salt.isel(xi_rho=20, eta_rho=100, s_rho=0)))
        temp_tmp = ds.variables["temp"]
        salt_tmp = ds.variables["salt"]
        zeta_tmp = ds.variables["zeta"]
        temp_tmp_ann = np.nanmean(temp_tmp, axis=0)
        salt_tmp_ann = np.nanmean(salt_tmp, axis=0)
        print('size temp_tmp_ann = ', temp_tmp_ann.shape)

        # calculate MKE and EKE using bottom velocity:
        # 1. load bottom vel
        mmax = 12
        usfc = ds.variables["u"].isel(s_rho=-1)
        vsfc = ds.variables["v"].isel(s_rho=-1)
        # 2. interpolate to rho grid
        w10_def = pyresample.geometry.SwathDefinition(lons=lon_rho_10km,lats=lat_rho_10km)
        w10u_def = pyresample.geometry.SwathDefinition(lons=lon_u_10km,lats=lat_u_10km)
        w10v_def = pyresample.geometry.SwathDefinition(lons=lon_v_10km,lats=lat_v_10km)
        wf = lambda r: 1/r
        usfci = np.empty(ds.salt.isel(s_rho=0).shape)
        vsfci = np.empty(ds.salt.isel(s_rho=0).shape)
        for mm in np.arange(0,mmax):
            usfci[mm,:,:] = pyresample.kd_tree.resample_custom(w10u_def,usfc[mm,:,:].values,w10_def,\
                                             radius_of_influence=5000,neighbours=1,weight_funcs=wf)
            vsfci[mm,:,:] = pyresample.kd_tree.resample_custom(w10v_def,vsfc[mm,:,:].values,w10_def,\
                                             radius_of_influence=5000,neighbours=1,weight_funcs=wf)

        usfc_ann_tmp = np.nanmean(usfci, axis=0)
        vsfc_ann_tmp = np.nanmean(vsfci, axis=0)
        V_ann_tmp = np.sqrt(usfc_ann_tmp*usfc_ann_tmp + vsfc_ann_tmp*vsfc_ann_tmp)
        MKE_tmp = np.divide(V_ann_tmp*V_ann_tmp,2)
        print("Size MKE_tmp", MKE_tmp.shape)
        
        usfc_anom_tmp = usfci - usfc_ann_tmp
        vsfc_anom_tmp = vsfci - vsfc_ann_tmp
        V_anom_tmp = np.sqrt(usfc_anom_tmp*usfc_anom_tmp + vsfc_anom_tmp*vsfc_anom_tmp)
        V_anom_ann = np.nanmean(V_anom_tmp, axis=0)
        EKE_tmp = np.divide(V_anom_ann*V_anom_ann,2)
        print("Size EKE_tmp", EKE_tmp.shape)

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

        # concantenate annual averaged temp/salt
        if yr == '20':
            temp_ann = temp_tmp_ann
            salt_ann = salt_tmp_ann
            z_w_ann = z_w_tmp_ann
            z_rho_ann = z_rho_tmp_ann
            MKE = MKE_tmp
            EKE = EKE_tmp
#        elif yr == '02':
#            temp_ann = np.stack((temp_ann,temp_tmp_ann), axis=0)
#            salt_ann = np.stack((salt_ann,salt_tmp_ann), axis=0)
#            z_w_ann = np.stack((z_w_ann,z_w_tmp_ann), axis=0)
#            z_rho_ann = np.stack((z_rho_ann,z_rho_tmp_ann), axis=0)
#            MKE = np.stack((MKE,MKE_tmp), axis=0)
#            EKE = np.stack((EKE,EKE_tmp), axis=0)
#        else:
#            temp_tmp_ann_4thdim = np.expand_dims(temp_tmp_ann, axis=0)
#            temp_ann = np.concatenate((temp_ann,temp_tmp_ann_4thdim), axis=0)
#            salt_tmp_ann_4thdim = np.expand_dims(salt_tmp_ann, axis=0)
#            salt_ann = np.concatenate((salt_ann,salt_tmp_ann_4thdim), axis=0)
#            z_w_tmp_ann_4thdim = np.expand_dims(z_w_tmp_ann, axis=0)
#            z_w_ann = np.concatenate((z_w_ann,z_w_tmp_ann_4thdim), axis=0)
#            z_rho_tmp_ann_4thdim = np.expand_dims(z_rho_tmp_ann, axis=0)
#            z_rho_ann = np.concatenate((z_rho_ann,z_rho_tmp_ann_4thdim), axis=0)
#            MKE_tmp_4thdim = np.expand_dims(MKE_tmp, axis=0)
#            MKE = np.concatenate((MKE,MKE_tmp_4thdim), axis=0)
#            EKE_tmp_4thdim = np.expand_dims(EKE_tmp, axis=0)
#            EKE = np.concatenate((EKE,EKE_tmp_4thdim), axis=0)
#            del temp_tmp_ann_4thdim, salt_tmp_ann_4thdim, z_w_tmp_ann_4thdim, z_rho_tmp_ann_4thdim, MKE_tmp_4thdim, EKE_tmp_4thdim

    return temp_ann, salt_ann, MKE, EKE

# read ROMS 4km:
def read_roms_ts_4km(exp_path,yr_file):
    for yr in yr_file:
        ds = xr.open_dataset(exp_path + 'ocean_avg_00' + yr_file + '.nc')
        print('size temp and time length: ', ds.temp.shape, len(ds.salt.isel(xi_rho=20, eta_rho=100, s_rho=0)))
        temp_tmp = ds.variables["temp"]
        salt_tmp = ds.variables["salt"]
        zeta_tmp = ds.variables["zeta"]
        temp_tmp_ann = np.nanmean(temp_tmp, axis=0)
        salt_tmp_ann = np.nanmean(salt_tmp, axis=0)
        print('size temp_tmp_ann = ', temp_tmp_ann.shape)
        del temp_tmp, salt_tmp

        # calculate MKE and EKE using bottom velocity:
        # 1. load bottom vel
        mmax = 12
        usfc = ds.variables["u"].isel(s_rho=-1)
        vsfc = ds.variables["v"].isel(s_rho=-1)
        # 2. interpolate to rho grid
        w4_def = pyresample.geometry.SwathDefinition(lons=lon_rho_4km,lats=lat_rho_4km)
        w4u_def = pyresample.geometry.SwathDefinition(lons=lon_u_4km,lats=lat_u_4km)
        w4v_def = pyresample.geometry.SwathDefinition(lons=lon_v_4km,lats=lat_v_4km)
        wf = lambda r: 1/r
        usfci = np.empty(ds.salt.isel(s_rho=0).shape)
        vsfci = np.empty(ds.salt.isel(s_rho=0).shape)
        for mm in np.arange(0,mmax):
            usfci[mm,:,:] = pyresample.kd_tree.resample_custom(w4u_def,usfc[mm,:,:].values,w4_def,\
                                             radius_of_influence=5000,neighbours=1,weight_funcs=wf)
            vsfci[mm,:,:] = pyresample.kd_tree.resample_custom(w4v_def,vsfc[mm,:,:].values,w4_def,\
                                             radius_of_influence=5000,neighbours=1,weight_funcs=wf)

        usfc_ann_tmp = np.nanmean(usfci, axis=0)
        vsfc_ann_tmp = np.nanmean(vsfci, axis=0)
        V_ann_tmp = np.sqrt(usfc_ann_tmp*usfc_ann_tmp + vsfc_ann_tmp*vsfc_ann_tmp)
        MKE_tmp = np.divide(V_ann_tmp*V_ann_tmp,2)
        print("Size MKE_tmp", MKE_tmp.shape)

        usfc_anom_tmp = usfci - usfc_ann_tmp
        vsfc_anom_tmp = vsfci - vsfc_ann_tmp
        V_anom_tmp = np.sqrt(usfc_anom_tmp*usfc_anom_tmp + vsfc_anom_tmp*vsfc_anom_tmp)
        V_anom_ann = np.nanmean(V_anom_tmp, axis=0)
        EKE_tmp = np.divide(V_anom_ann*V_anom_ann,2)
        print("Size EKE_tmp", EKE_tmp.shape)

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

        # concantenate annual averaged temp/salt
        if yr_file == '10' or yr_file == '01':
            temp_ann = temp_tmp_ann
            salt_ann = salt_tmp_ann
            z_w_ann = z_w_tmp_ann
            z_rho_ann = z_rho_tmp_ann
            MKE_ann = MKE_tmp
            EKE_ann = EKE_tmp
#        elif yr == '02':
#            temp_ann = np.stack((temp_ann,temp_tmp_ann), axis=0)
#            salt_ann = np.stack((salt_ann,salt_tmp_ann), axis=0)
#            z_w_ann = np.stack((z_w_ann,z_w_tmp_ann), axis=0)
#            z_rho_ann = np.stack((z_rho_ann,z_rho_tmp_ann), axis=0)
#            MKE_ann = np.stack((MKE_ann,MKE_tmp), axis=0)
#            EKE_ann = np.stack((EKE_ann,EKE_tmp), axis=0)
#        else:
#            temp_tmp_ann_4thdim = np.expand_dims(temp_tmp_ann, axis=0)
#            temp_ann = np.concatenate((temp_ann,temp_tmp_ann_4thdim), axis=0)
#            salt_tmp_ann_4thdim = np.expand_dims(salt_tmp_ann, axis=0)
#            salt_ann = np.concatenate((salt_ann,salt_tmp_ann_4thdim), axis=0)
#            z_w_tmp_ann_4thdim = np.expand_dims(z_w_tmp_ann, axis=0)
#            z_w_ann = np.concatenate((z_w_ann,z_w_tmp_ann_4thdim), axis=0)
#            z_rho_tmp_ann_4thdim = np.expand_dims(z_rho_tmp_ann, axis=0)
#            z_rho_ann = np.concatenate((z_rho_ann,z_rho_tmp_ann_4thdim), axis=0)
#            MKE_tmp_4thdim = np.expand_dims(MKE_tmp, axis=0)
#            MKE_ann = np.concatenate((MKE_ann,MKE_tmp_4thdim), axis=0)
#            EKE_tmp_4thdim = np.expand_dims(EKE_tmp, axis=0)
#            EKE_ann = np.concatenate((EKE_ann,EKE_tmp_4thdim), axis=0)
#            del temp_tmp_ann_4thdim, salt_tmp_ann_4thdim, z_w_tmp_ann_4thdim, z_rho_tmp_ann_4thdim, MKE_tmp_4thdim, EKE_tmp_4thdim

        print('Annual temp and annual tmp temp sizes = ', temp_ann.shape, temp_tmp_ann.shape)
        print('Annual z_w and annual tmp z_w sizes = ', z_w_ann.shape, z_w_tmp_ann.shape)

        del temp_tmp_ann, salt_tmp_ann, z_w_tmp_ann, z_rho_tmp_ann

    print('Annual temp, salt, z_w, z_rho sizes = ', temp_ann.shape, salt_ann.shape, z_w_ann.shape, z_rho_ann.shape)

    return temp_ann, salt_ann, MKE_ann, EKE_ann


# read ROMS 2km:
def read_roms_ts_2km(exp_path):
    for yr in ['09','10']:
        ds = xr.open_dataset(exp_path + 'ocean_avg_00' + yr + '.nc')
        print('size temp and time length: ', ds.temp.shape, len(ds.salt.isel(xi_rho=20, eta_rho=100, s_rho=0)))
        temp_tmp = ds.variables["temp"]
        salt_tmp = ds.variables["salt"]
        zeta_tmp = ds.variables["zeta"]
        temp_tmp_ann = np.nanmean(temp_tmp, axis=0)
        salt_tmp_ann = np.nanmean(salt_tmp, axis=0)
        print('size temp_tmp_ann = ', temp_tmp_ann.shape)
        del temp_tmp, salt_tmp

        # calculate MKE and EKE using bottom velocity:
        # 1. load bottom vel
        mmax = len(zeta_tmp)
        usfc = ds.variables["u"].isel(s_rho=-1)
        vsfc = ds.variables["v"].isel(s_rho=-1)
        # 2. interpolate to rho grid
        w2_def = pyresample.geometry.SwathDefinition(lons=lon_rho_2km,lats=lat_rho_2km)
        w2u_def = pyresample.geometry.SwathDefinition(lons=lon_u_2km,lats=lat_u_2km)
        w2v_def = pyresample.geometry.SwathDefinition(lons=lon_v_2km,lats=lat_v_2km)
        wf = lambda r: 1/r
        usfci = np.empty(ds.salt.isel(s_rho=0).shape)
        vsfci = np.empty(ds.salt.isel(s_rho=0).shape)
        for mm in np.arange(0,mmax):
            usfci[mm,:,:] = pyresample.kd_tree.resample_custom(w2u_def,usfc[mm,:,:].values,w2_def,\
                                             radius_of_influence=5000,neighbours=1,weight_funcs=wf)
            vsfci[mm,:,:] = pyresample.kd_tree.resample_custom(w2v_def,vsfc[mm,:,:].values,w2_def,\
                                             radius_of_influence=5000,neighbours=1,weight_funcs=wf)
        
        usfc_ann_tmp = np.nanmean(usfci, axis=0)
        vsfc_ann_tmp = np.nanmean(vsfci, axis=0)
        V_ann_tmp = np.sqrt(usfc_ann_tmp*usfc_ann_tmp + vsfc_ann_tmp*vsfc_ann_tmp)
        MKE_tmp = np.divide(V_ann_tmp*V_ann_tmp,2)
        print("Size MKE_tmp", MKE_tmp.shape)

        usfc_anom_tmp = usfci - usfc_ann_tmp
        vsfc_anom_tmp = vsfci - vsfc_ann_tmp
        V_anom_tmp = np.sqrt(usfc_anom_tmp*usfc_anom_tmp + vsfc_anom_tmp*vsfc_anom_tmp)
        V_anom_ann = np.nanmean(V_anom_tmp, axis=0)
        EKE_tmp = np.divide(V_anom_ann*V_anom_ann,2)
        print("Size EKE_tmp", EKE_tmp.shape)

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

        # concantenate annual averaged temp/salt
        if yr == '09':
            temp_ann = np.expand_dims(temp_tmp_ann, axis=0)
            salt_ann = np.expand_dims(salt_tmp_ann, axis=0)
            z_w_ann = np.expand_dims(z_w_tmp_ann, axis=0)
            z_rho_ann = np.expand_dims(z_rho_tmp_ann, axis=0)
            MKE_ann = np.expand_dims(MKE_tmp, axis=0)
            EKE_ann = np.expand_dims(EKE_tmp, axis=0)
        elif yr == '10':
            temp_ann = np.squeeze(np.nanmean(np.stack((temp_ann,np.expand_dims(temp_tmp_ann, axis=0)), axis=0), axis=0))
            salt_ann = np.squeeze(np.nanmean(np.stack((salt_ann,np.expand_dims(salt_tmp_ann, axis=0)), axis=0), axis=0))
            z_w_ann = np.squeeze(np.nanmean(np.stack((z_w_ann,np.expand_dims(z_w_tmp_ann, axis=0)), axis=0), axis=0))
            z_rho_ann = np.squeeze(np.nanmean(np.stack((z_rho_ann,np.expand_dims(z_rho_tmp_ann, axis=0)), axis=0), axis=0))
            MKE_ann = np.squeeze(np.nanmean(np.stack((MKE_ann,np.expand_dims(MKE_tmp, axis=0)), axis=0), axis=0))
            EKE_ann = np.squeeze(np.nanmean(np.stack((EKE_ann,np.expand_dims(EKE_tmp, axis=0)), axis=0), axis=0))

        print('Annual temp and annual tmp temp sizes = ', temp_ann.shape, temp_tmp_ann.shape)

    return temp_ann, salt_ann, MKE_ann, EKE_ann

path_ECCO2_10km = '/scratch/project_2000339/boeiradi/waom10extend_shflim_S_0.25Q/output_01-20yr/'
path_ECCO2_4km = '/scratch/project_2000789/boeiradi/waom4extend_shflim_S_0.25Q/output_01-20yr/'
path_ECCO2_4kmB = '/scratch/project_2000789/boeiradi/waom4extend_shflim_S_0.25Q/output_yr10_10km-bathy/'
path_ECCO2_2km = '/scratch/project_2000339/boeiradi/waom2extend_shflim_S_0.25Q/output_01-05yr/'

temp_ann_10km, salt_ann_10km, MKE_10km, EKE_10km = read_roms_ts_10km(path_ECCO2_10km)
temp_ann_4km, salt_ann_4km, MKE_4km, EKE_4km = read_roms_ts_4km(path_ECCO2_4km,'10')
temp_ann_4kmB, salt_ann_4kmB, MKE_4kmB, EKE_4kmB = read_roms_ts_4km(path_ECCO2_4kmB,'01')
temp_ann_2km, salt_ann_2km, MKE_2km, EKE_2km = read_roms_ts_2km(path_ECCO2_2km)

# define mask zice:
mask_zice_10km = ma.masked_where(zice_10km < 0, np.ones(zice_10km.shape))
mask_zice_4km = ma.masked_where(zice_4km < 0, np.ones(zice_4km.shape))
mask_zice_2km = ma.masked_where(zice_2km < 0, np.ones(zice_2km.shape))
print('Print temp, MKE, EKE shapes =',temp_ann_10km.shape, MKE_10km.shape, EKE_10km.shape)
mask_shelf_10km = ma.masked_where(h_10km > 1500, np.ones(zice_10km.shape))

fig_path='/users/boeiradi/COLD_project/postprocessing/figs/Maps_validations/'

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

# print array shapes before plotting:
print('lat_rho_2km shape', lat_rho_2km.shape)
print('lon_rho_2km shape', lon_rho_2km.shape)
print('temp_ann_2km (squeezed) shape', np.squeeze(temp_ann_2km[0,:,:]).shape)
print('temp_ann_2km shape', temp_ann_2km.shape)

ratio = .9
## MKE/EKE
MEmin = 0.
MEmax = 0.01
EKmin = 0.
EKmax = 0.005
# call cartopy projection
proj = ccrs.SouthPolarStereo()
fig = plt.figure(figsize=(8,10))

ax1 = fig.add_subplot(321, projection=proj)
ct=plt.pcolormesh(lon_rho_10km,lat_rho_10km,MKE_10km*mask_zice_10km, transform=ccrs.PlateCarree(), cmap=plt.cm.plasma, vmin=MEmin, vmax=MEmax)
#cbart =fig.colorbar(ct)# , extend='both')
plt.title('MKE \n WAOM10')
ax1.gridlines() # draw_labels=True,linewidth=
#ax1.coastlines(resolution='110m')
#lonlat_labels(ax1)
x_left, x_right = ax1.get_xlim()
y_low, y_high = ax1.get_ylim()
ax1.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax1.set_extent([-180, 180, -90, -62], crs=ccrs.PlateCarree())
ax1.add_feature(cfeature.LAND, zorder=3, facecolor='white')

ax2 = fig.add_subplot(322, projection=proj)
cs=plt.pcolormesh(lon_rho_10km,lat_rho_10km,EKE_10km*mask_zice_10km, transform=ccrs.PlateCarree(), cmap=plt.cm.plasma, vmin=EKmin, vmax=EKmax)
#cbars =fig.colorbar(cs)# , extend='both')
plt.title('EKE \n WAOM10')
#ax2.coastlines(resolution='110m')
ax2.gridlines()
#lonlat_labels(ax2)
x_left, x_right = ax2.get_xlim()
y_low, y_high = ax2.get_ylim()
ax2.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax2.set_extent([-180, 180, -90, -62], crs=ccrs.PlateCarree())
ax2.add_feature(cfeature.LAND, zorder=3, facecolor='white')

ax3 = fig.add_subplot(323, projection=proj)
ct=plt.pcolormesh(lon_rho_4km,lat_rho_4km,MKE_4km*mask_zice_4km, transform=ccrs.PlateCarree(), cmap=plt.cm.plasma, vmin=MEmin, vmax=MEmax)
#cbart =fig.colorbar(ct)# , extend='both')
plt.title('WAOM4')
ax3.gridlines() # draw_labels=True,linewidth=
#ax3.coastlines(resolution='110m')
#lonlat_labels(ax3)
x_left, x_right = ax3.get_xlim()
y_low, y_high = ax3.get_ylim()
ax3.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax3.set_extent([-180, 180, -90, -62], crs=ccrs.PlateCarree())
ax3.add_feature(cfeature.LAND, zorder=3, facecolor='white')

ax4 = fig.add_subplot(324, projection=proj)
cs=plt.pcolormesh(lon_rho_4km,lat_rho_4km,EKE_4km*mask_zice_4km, transform=ccrs.PlateCarree(), cmap=plt.cm.plasma, vmin=EKmin, vmax=EKmax)
#cbars =fig.colorbar(cs)# , extend='both')
plt.title('WAOM4')
#ax4.coastlines(resolution='110m')
ax4.gridlines()
#lonlat_labels(ax4)
x_left, x_right = ax4.get_xlim()
y_low, y_high = ax4.get_ylim()
ax4.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax4.set_extent([-180, 180, -90, -62], crs=ccrs.PlateCarree())
ax4.add_feature(cfeature.LAND, zorder=3, facecolor='white')

ax5 = fig.add_subplot(325, projection=proj)
ct=plt.pcolormesh(lon_rho_2km,lat_rho_2km,MKE_2km*mask_zice_2km, transform=ccrs.PlateCarree(), cmap=plt.cm.plasma, vmin=MEmin, vmax=MEmax)
#cbart =fig.colorbar(ct)# , extend='both')
plt.title('WAOM2')
ax5.gridlines() # draw_labels=True,linewidth=
#ax5.coastlines(resolution='110m')
#lonlat_labels(ax5)
x_left, x_right = ax5.get_xlim()
y_low, y_high = ax5.get_ylim()
ax5.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax5.set_extent([-180, 180, -90, -62], crs=ccrs.PlateCarree())
ax5.add_feature(cfeature.LAND, zorder=3, facecolor='white')

ax6 = fig.add_subplot(326, projection=proj)
cs=plt.pcolormesh(lon_rho_2km,lat_rho_2km,EKE_2km*mask_zice_2km, transform=ccrs.PlateCarree(), cmap=plt.cm.plasma, vmin=EKmin, vmax=EKmax)
#cbars =fig.colorbar(cs)# , extend='both')
plt.title('WAOM2')
#ax6.coastlines(resolution='110m')
ax6.gridlines()
#lonlat_labels(ax6)
x_left, x_right = ax6.get_xlim()
y_low, y_high = ax6.get_ylim()
ax6.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax6.set_extent([-180, 180, -90, -62], crs=ccrs.PlateCarree())
ax6.add_feature(cfeature.LAND, zorder=3, facecolor='white')

cbar_ax1 = fig.add_axes([0.12, 0.05, 0.35, 0.01])
fig.colorbar(ct, cax=cbar_ax1, orientation='horizontal')
cbar_ax1.set_xlabel('MKE (m$^2$ s$^{-2}$)')#, labelpad=-35)

cbar_ax2 = fig.add_axes([0.55, 0.05, 0.35, 0.01])
fig.colorbar(cs, cax=cbar_ax2, orientation='horizontal')
cbar_ax2.set_xlabel('EKE (m$^2$ s$^{-2}$)')#, labelpad=-35)

name_fig="waom10x4x2extend_shflim_S_0.25Q_MKExEKE_surface_maps.png"
plt.savefig(fig_path + name_fig, dpi=300)


## ==========================================
# interpolate WAOM4 and WAOM2 field to the WAOM10 grid for maps of difference:

w10_def = pyresample.geometry.SwathDefinition(lons=lon_rho_10km,lats=lat_rho_10km)
w4_def = pyresample.geometry.SwathDefinition(lons=lon_rho_4km,lats=lat_rho_4km)
w2_def = pyresample.geometry.SwathDefinition(lons=lon_rho_2km,lats=lat_rho_2km)
 
wf = lambda r: 1/r

#w10_on_waom4 = pyresample.kd_tree.resample_custom(w10_def,waom10.h,w4_def,\
#                                         radius_of_influence=30000,neighbours=4,weight_funcs=wf)
#temp_ann_4km[-1,:,:]
#temp_ann_2km[-1,:,:]
#salt_ann_4km[-1,:,:]
#salt_ann_2km[-1,:,:]

sfc_temp_ann_4km_interp = pyresample.kd_tree.resample_custom(w4_def,temp_ann_4km[-1,:,:],w10_def,\
                                         radius_of_influence=30000,neighbours=4,weight_funcs=wf)
sfc_salt_ann_4km_interp = pyresample.kd_tree.resample_custom(w4_def,salt_ann_4km[-1,:,:],w10_def,\
                                         radius_of_influence=30000,neighbours=4,weight_funcs=wf)
bot_temp_ann_4km_interp = pyresample.kd_tree.resample_custom(w4_def,temp_ann_4km[0,:,:],w10_def,\
                                         radius_of_influence=30000,neighbours=4,weight_funcs=wf)
bot_salt_ann_4km_interp = pyresample.kd_tree.resample_custom(w4_def,salt_ann_4km[0,:,:],w10_def,\
                                         radius_of_influence=30000,neighbours=4,weight_funcs=wf)
MKE_4km_interp = pyresample.kd_tree.resample_custom(w4_def,MKE_4km,w10_def,\
                                         radius_of_influence=30000,neighbours=4,weight_funcs=wf)
EKE_4km_interp = pyresample.kd_tree.resample_custom(w4_def,EKE_4km,w10_def,\
                                         radius_of_influence=30000,neighbours=4,weight_funcs=wf)
MKE_4kmB_interp = pyresample.kd_tree.resample_custom(w4_def,MKE_4kmB,w10_def,\
                                         radius_of_influence=30000,neighbours=4,weight_funcs=wf)
EKE_4kmB_interp = pyresample.kd_tree.resample_custom(w4_def,EKE_4kmB,w10_def,\
                                         radius_of_influence=30000,neighbours=4,weight_funcs=wf)
sfc_temp_ann_2km_interp = pyresample.kd_tree.resample_custom(w2_def,temp_ann_2km[-1,:,:],w10_def,\
                                         radius_of_influence=5000,neighbours=1,weight_funcs=wf)
sfc_salt_ann_2km_interp = pyresample.kd_tree.resample_custom(w2_def,salt_ann_2km[-1,:,:],w10_def,\
                                         radius_of_influence=5000,neighbours=1,weight_funcs=wf)
bot_temp_ann_2km_interp = pyresample.kd_tree.resample_custom(w2_def,temp_ann_2km[0,:,:],w10_def,\
                                         radius_of_influence=5000,neighbours=1,weight_funcs=wf)
bot_salt_ann_2km_interp = pyresample.kd_tree.resample_custom(w2_def,salt_ann_2km[0,:,:],w10_def,\
                                         radius_of_influence=5000,neighbours=1,weight_funcs=wf)
MKE_2km_interp = pyresample.kd_tree.resample_custom(w2_def,MKE_2km,w10_def,\
                                         radius_of_influence=5000,neighbours=1,weight_funcs=wf)
EKE_2km_interp = pyresample.kd_tree.resample_custom(w2_def,EKE_2km,w10_def,\
                                         radius_of_influence=5000,neighbours=1,weight_funcs=wf)

# for bathymetric isoline of 1500
bathym = cfeature.NaturalEarthFeature(name='bathymetry_J_1000', scale='10m', category='physical')

# call cartopy projection
proj = ccrs.SouthPolarStereo()

# MKE/EKE for RFIS
# limits for contour of ice front (Ronne-Filchner IS):
xlimit = np.arange(1500,2500,1)
ylimit = np.arange(500,1500,1)

## MKE/EKE
MEmin = 0.
MEmax = 0.01
EKmin = 0.
EKmax = 0.001
# call cartopy projection
proj = ccrs.SouthPolarStereo()
fig = plt.figure(figsize=(8,10))


ax1 = fig.add_subplot(321, projection=proj)
ct1=plt.pcolormesh(lon_rho_10km,lat_rho_10km,np.squeeze(MKE_10km[:,:]), transform=ccrs.PlateCarree(), cmap=plt.cm.plasma, vmin=MEmin, vmax=MEmax)
plt.title('MKE \n WAOM10')
ax1.gridlines() # draw_labels=True,linewidth=
ax1.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
ratio = .9
x_left, x_right = ax1.get_xlim()
y_low, y_high = ax1.get_ylim()
#ax1.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax1.add_feature(cfeature.LAND, zorder=3, facecolor='white')
plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],h_2km[xlimit,ylimit],levels=(2000,2001), transform=ccrs.PlateCarree(), colors='k', linewidths=1)
plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],zice_2km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],h_2km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='grey')

#cbar_ax3 = fig.add_axes([0.04, 0.65, 0.01, 0.25])
#fig.colorbar(ct1, cax=cbar_ax3, orientation='vertical')
#cbar_ax3.set_xlabel('$^{\circ}$C')#, labelpad=-35)

ax2 = fig.add_subplot(322, projection=proj)
cs1=plt.pcolormesh(lon_rho_10km,lat_rho_10km,np.squeeze(EKE_10km[:,:]), transform=ccrs.PlateCarree(), cmap=plt.cm.plasma, vmin=EKmin, vmax=EKmax)
plt.title('EKE \n WAOM10')
ax2.gridlines()
ax2.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
x_left, x_right = ax2.get_xlim()
y_low, y_high = ax2.get_ylim()
#ax2.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax2.add_feature(cfeature.LAND, zorder=3, facecolor='white')
plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],h_2km[xlimit,ylimit],levels=(2000,2001), transform=ccrs.PlateCarree(), colors='k', linewidths=1)
plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],zice_2km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],h_2km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='grey')

#cbar_ax4 = fig.add_axes([0.92, 0.65, 0.01, 0.25])
#fig.colorbar(cs1, cax=cbar_ax4, orientation='vertical')
#cbar_ax4.set_xlabel('')#, labelpad=-35)

ax3 = fig.add_subplot(323, projection=proj)
ct2=plt.pcolormesh(lon_rho_10km,lat_rho_10km,MKE_4km_interp, transform=ccrs.PlateCarree(), cmap=plt.cm.plasma, vmin=MEmin, vmax=MEmax)
plt.title('WAOM4')
ax3.gridlines() # draw_labels=True,linewidth=
ax3.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
x_left, x_right = ax3.get_xlim()
y_low, y_high = ax3.get_ylim()
#ax3.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax3.add_feature(cfeature.LAND, zorder=3, facecolor='white')
plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],h_2km[xlimit,ylimit],levels=(2000,2001), transform=ccrs.PlateCarree(), colors='k', linewidths=1)
plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],zice_2km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],h_2km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='grey')

ax4 = fig.add_subplot(324, projection=proj)
cs2=plt.pcolormesh(lon_rho_10km,lat_rho_10km,EKE_4km_interp, transform=ccrs.PlateCarree(), cmap=plt.cm.plasma, vmin=EKmin, vmax=EKmax)
plt.title('WAOM4')
ax4.gridlines()
ax4.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
x_left, x_right = ax4.get_xlim()
y_low, y_high = ax4.get_ylim()
#ax4.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax4.add_feature(cfeature.LAND, zorder=3, facecolor='white')
plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],h_2km[xlimit,ylimit],levels=(2000,2001), transform=ccrs.PlateCarree(), colors='k', linewidths=1)
plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],zice_2km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],h_2km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='grey')

ax5 = fig.add_subplot(325, projection=proj)
ct2=plt.pcolormesh(lon_rho_10km,lat_rho_10km,MKE_2km_interp, transform=ccrs.PlateCarree(), cmap=plt.cm.plasma, vmin=MEmin, vmax=MEmax)
plt.title('WAOM2')
ax5.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
x_left, x_right = ax5.get_xlim()
y_low, y_high = ax5.get_ylim()
#ax5.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax5.add_feature(cfeature.LAND, zorder=3, facecolor='white')
plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],h_2km[xlimit,ylimit],levels=(2000,2001), transform=ccrs.PlateCarree(), colors='k', linewidths=1)
plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],zice_2km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],h_2km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='grey')

ax6 = fig.add_subplot(326, projection=proj)
cs2=plt.pcolormesh(lon_rho_10km,lat_rho_10km,EKE_2km_interp, transform=ccrs.PlateCarree(), cmap=plt.cm.plasma, vmin=EKmin, vmax=EKmax)
plt.title('WAOM2')
ax6.gridlines()
ax6.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
x_left, x_right = ax6.get_xlim()
y_low, y_high = ax6.get_ylim()
#ax6.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax6.add_feature(cfeature.LAND, zorder=3, facecolor='white')
plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],h_2km[xlimit,ylimit],levels=(2000,2001), transform=ccrs.PlateCarree(), colors='k', linewidths=1)
plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],zice_2km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],h_2km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='grey')

cbar_ax1 = fig.add_axes([0.12, 0.08, 0.35, 0.01])
fig.colorbar(ct2, cax=cbar_ax1, orientation='horizontal')
cbar_ax1.set_xlabel('$^{\circ}$C')#, labelpad=-35)

cbar_ax2 = fig.add_axes([0.55, 0.08, 0.35, 0.01])
fig.colorbar(cs2, cax=cbar_ax2, orientation='horizontal')
cbar_ax2.set_xlabel('')#, labelpad=-35)

name_fig="waom10x4x2extend_shflim_S_0.25Q_MKE-EKE_surface_RFIS_anomalies.png"
plt.savefig(fig_path + name_fig, dpi=300)
plt.close()

# --  East Antarctic MKExEKE
# limits for contour of ice front (East Antarctica):
xlimit = np.arange(250,2500,1)
ylimit = np.arange(2000,3000,1)

MEmin = 0.
MEmax = 0.01
EKmin = 0.
EKmax = 0.001
# call cartopy projection
proj = ccrs.SouthPolarStereo()
fig = plt.figure(figsize=(10,8))


ax1 = fig.add_subplot(231, projection=proj)
ct1=plt.pcolormesh(lon_rho_10km,lat_rho_10km,np.squeeze(MKE_10km[:,:]), transform=ccrs.PlateCarree(), cmap=plt.cm.plasma, vmin=MEmin, vmax=MEmax)
plt.title('MKE \n WAOM10')
ax1.gridlines() # draw_labels=True,linewidth=
ax1.set_extent([65, 150, -70, -63.5], crs=ccrs.PlateCarree())
ratio = .9
x_left, x_right = ax1.get_xlim()
y_low, y_high = ax1.get_ylim()
#ax1.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax1.add_feature(cfeature.LAND, zorder=3, facecolor='white')
plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],zice_2km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
#plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],h_2km[xlimit,ylimit],levels=(2000,2001), transform=ccrs.PlateCarree(), colors='k', linewidths=1)
#plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],h_2km[xlimit,ylimit],levels=(300,6000), transform=ccrs.PlateCarree(), colors='grey', linewidths=1)

#cbar_ax3 = fig.add_axes([0.04, 0.65, 0.01, 0.25])
#fig.colorbar(ct1, cax=cbar_ax3, orientation='vertical')
#cbar_ax3.set_xlabel('$^{\circ}$C')#, labelpad=-35)

ax2 = fig.add_subplot(234, projection=proj)
cs1=plt.pcolormesh(lon_rho_10km,lat_rho_10km,np.squeeze(EKE_10km[:,:]), transform=ccrs.PlateCarree(), cmap=plt.cm.plasma, vmin=EKmin, vmax=EKmax)
plt.title('EKE \n WAOM10')
ax2.gridlines()
ax2.set_extent([65, 150, -70, -63.5], crs=ccrs.PlateCarree())
x_left, x_right = ax2.get_xlim()
y_low, y_high = ax2.get_ylim()
#ax2.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax2.add_feature(cfeature.LAND, zorder=3, facecolor='white')
plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],zice_2km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
#plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],h_2km[xlimit,ylimit],levels=(2000,2001), transform=ccrs.PlateCarree(), colors='k', linewidths=1)
#plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],h_2km[xlimit,ylimit],levels=(300,6000), transform=ccrs.PlateCarree(), colors='grey', linewidths=1)

#cbar_ax4 = fig.add_axes([0.92, 0.65, 0.01, 0.25])
#fig.colorbar(cs1, cax=cbar_ax4, orientation='vertical')
#cbar_ax4.set_xlabel('')#, labelpad=-35)

ax3 = fig.add_subplot(232, projection=proj)
ct2=plt.pcolormesh(lon_rho_10km,lat_rho_10km,MKE_4km_interp, transform=ccrs.PlateCarree(), cmap=plt.cm.plasma, vmin=MEmin, vmax=MEmax)
plt.title('MKE \n WAOM4')
ax3.gridlines() # draw_labels=True,linewidth=
ax3.set_extent([65, 150, -70, -63.5], crs=ccrs.PlateCarree())
x_left, x_right = ax3.get_xlim()
y_low, y_high = ax3.get_ylim()
#ax3.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax3.add_feature(cfeature.LAND, zorder=3, facecolor='white')
plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],zice_2km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
#plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],h_2km[xlimit,ylimit],levels=(2000,2001), transform=ccrs.PlateCarree(), colors='k', linewidths=1)
#plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],h_2km[xlimit,ylimit],levels=(300,6000), transform=ccrs.PlateCarree(), colors='grey', linewidths=1)

ax4 = fig.add_subplot(235, projection=proj)
cs2=plt.pcolormesh(lon_rho_10km,lat_rho_10km,EKE_4km_interp, transform=ccrs.PlateCarree(), cmap=plt.cm.plasma, vmin=EKmin, vmax=EKmax)
plt.title('EKE \n WAOM4')
ax4.gridlines()
ax4.set_extent([65, 150, -70, -63.5], crs=ccrs.PlateCarree())
x_left, x_right = ax4.get_xlim()
y_low, y_high = ax4.get_ylim()
#ax4.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax4.add_feature(cfeature.LAND, zorder=3, facecolor='white')
plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],zice_2km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
#plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],h_2km[xlimit,ylimit],levels=(2000,2001), transform=ccrs.PlateCarree(), colors='k', linewidths=1)
#plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],h_2km[xlimit,ylimit],levels=(300,6000), transform=ccrs.PlateCarree(), colors='grey', linewidths=1)

ax5 = fig.add_subplot(233, projection=proj)
ct2=plt.pcolormesh(lon_rho_10km,lat_rho_10km,MKE_4kmB_interp, transform=ccrs.PlateCarree(), cmap=plt.cm.plasma, vmin=MEmin, vmax=MEmax)
plt.title('MKE \n WAOM4-COARSE')
ax5.set_extent([65, 150, -70, -63.5], crs=ccrs.PlateCarree())
x_left, x_right = ax5.get_xlim()
y_low, y_high = ax5.get_ylim()
#ax5.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax5.add_feature(cfeature.LAND, zorder=3, facecolor='white')
plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],zice_2km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
#plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],h_2km[xlimit,ylimit],levels=(2000,2001), transform=ccrs.PlateCarree(), colors='k', linewidths=1)
#plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],h_2km[xlimit,ylimit],levels=(300,6000), transform=ccrs.PlateCarree(), colors='grey', linewidths=1)

ax6 = fig.add_subplot(236, projection=proj)
cs2=plt.pcolormesh(lon_rho_10km,lat_rho_10km,EKE_4kmB_interp, transform=ccrs.PlateCarree(), cmap=plt.cm.plasma, vmin=EKmin, vmax=EKmax)
plt.title('EKE \n WAOM4-COARSE')
ax6.gridlines()
ax6.set_extent([65, 150, -70, -63.5], crs=ccrs.PlateCarree())
x_left, x_right = ax6.get_xlim()
y_low, y_high = ax6.get_ylim()
#ax6.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax6.add_feature(cfeature.LAND, zorder=3, facecolor='white')
plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],zice_2km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
#plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],h_2km[xlimit,ylimit],levels=(2000,2001), transform=ccrs.PlateCarree(), colors='k', linewidths=1)
#plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],h_2km[xlimit,ylimit],levels=(300,6000), transform=ccrs.PlateCarree(), colors='grey', linewidths=1)

cbar_ax1 = fig.add_axes([0.12, 0.08, 0.35, 0.01])
fig.colorbar(ct2, cax=cbar_ax1, orientation='horizontal')
cbar_ax1.set_xlabel('$^{\circ}$C')#, labelpad=-35)

cbar_ax2 = fig.add_axes([0.55, 0.08, 0.35, 0.01])
fig.colorbar(cs2, cax=cbar_ax2, orientation='horizontal')
cbar_ax2.set_xlabel('')#, labelpad=-35)

name_fig="waom10x4x4cextend_shflim_S_0.25Q_MKE-EKE_surface_EAnt_anomalies.png"
plt.savefig(fig_path + name_fig, dpi=300)
plt.close()

# West Antarctica MKExEKE

# --  East Antarctic MKExEKE
# limits for contour of ice front (Wast Antarctica):
xlimit = np.arange(500,2800,1)
ylimit = np.arange(250,1250,1)

MEmin = 0.
MEmax = 0.01
EKmin = 0.
EKmax = 0.001
# call cartopy projection
proj = ccrs.SouthPolarStereo()
fig = plt.figure(figsize=(10,8))


ax1 = fig.add_subplot(231, projection=proj)
ct1=plt.pcolormesh(lon_rho_10km,lat_rho_10km,np.squeeze(MKE_10km[:,:]), transform=ccrs.PlateCarree(), cmap=plt.cm.plasma, vmin=MEmin, vmax=MEmax)
plt.title('MKE \n WAOM10')
ax1.gridlines() # draw_labels=True,linewidth=
ax1.set_extent([-145, -68, -70, -67], crs=ccrs.PlateCarree())
ratio = .9
x_left, x_right = ax1.get_xlim()
y_low, y_high = ax1.get_ylim()
#ax1.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax1.add_feature(cfeature.LAND, zorder=3, facecolor='white')
plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],zice_2km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
#plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],h_2km[xlimit,ylimit],levels=(2000,2001), transform=ccrs.PlateCarree(), colors='k', linewidths=1)
#plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],h_2km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='grey')

#cbar_ax3 = fig.add_axes([0.04, 0.65, 0.01, 0.25])
#fig.colorbar(ct1, cax=cbar_ax3, orientation='vertical')
#cbar_ax3.set_xlabel('$^{\circ}$C')#, labelpad=-35)

ax2 = fig.add_subplot(234, projection=proj)
cs1=plt.pcolormesh(lon_rho_10km,lat_rho_10km,np.squeeze(EKE_10km[:,:]), transform=ccrs.PlateCarree(), cmap=plt.cm.plasma, vmin=EKmin, vmax=EKmax)
plt.title('EKE \n WAOM10')
ax2.gridlines()
ax2.set_extent([-145, -68, -70, -67], crs=ccrs.PlateCarree())
x_left, x_right = ax2.get_xlim()
y_low, y_high = ax2.get_ylim()
#ax2.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax2.add_feature(cfeature.LAND, zorder=3, facecolor='white')
plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],zice_2km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
#plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],h_2km[xlimit,ylimit],levels=(2000,2001), transform=ccrs.PlateCarree(), colors='k', linewidths=1)
#plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],h_2km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='grey')

#cbar_ax4 = fig.add_axes([0.92, 0.65, 0.01, 0.25])
#fig.colorbar(cs1, cax=cbar_ax4, orientation='vertical')
#cbar_ax4.set_xlabel('')#, labelpad=-35)

ax3 = fig.add_subplot(232, projection=proj)
ct2=plt.pcolormesh(lon_rho_10km,lat_rho_10km,MKE_4km_interp, transform=ccrs.PlateCarree(), cmap=plt.cm.plasma, vmin=MEmin, vmax=MEmax)
plt.title('MKE \n WAOM4')
ax3.gridlines() # draw_labels=True,linewidth=
ax3.set_extent([-145, -68, -70, -67], crs=ccrs.PlateCarree())
x_left, x_right = ax3.get_xlim()
y_low, y_high = ax3.get_ylim()
#ax3.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax3.add_feature(cfeature.LAND, zorder=3, facecolor='white')
plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],zice_2km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
#plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],h_2km[xlimit,ylimit],levels=(2000,2001), transform=ccrs.PlateCarree(), colors='k', linewidths=1)
#plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],h_2km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='grey')

ax4 = fig.add_subplot(235, projection=proj)
cs2=plt.pcolormesh(lon_rho_10km,lat_rho_10km,EKE_4km_interp, transform=ccrs.PlateCarree(), cmap=plt.cm.plasma, vmin=EKmin, vmax=EKmax)
plt.title('EKE \n WAOM4')
ax4.gridlines()
ax4.set_extent([-145, -68, -70, -67], crs=ccrs.PlateCarree())
x_left, x_right = ax4.get_xlim()
y_low, y_high = ax4.get_ylim()
#ax4.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax4.add_feature(cfeature.LAND, zorder=3, facecolor='white')
plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],zice_2km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
#plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],h_2km[xlimit,ylimit],levels=(2000,2001), transform=ccrs.PlateCarree(), colors='k', linewidths=1)
#plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],h_2km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='grey')

ax5 = fig.add_subplot(233, projection=proj)
ct2=plt.pcolormesh(lon_rho_10km,lat_rho_10km,MKE_2km_interp, transform=ccrs.PlateCarree(), cmap=plt.cm.plasma, vmin=MEmin, vmax=MEmax)
plt.title('MKE \n WAOM2')
ax5.set_extent([-145, -68, -70, -67], crs=ccrs.PlateCarree())
x_left, x_right = ax5.get_xlim()
y_low, y_high = ax5.get_ylim()
#ax5.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax5.add_feature(cfeature.LAND, zorder=3, facecolor='white')
plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],zice_2km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
#plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],h_2km[xlimit,ylimit],levels=(2000,2001), transform=ccrs.PlateCarree(), colors='k', linewidths=1)
#plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],h_2km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='grey')

ax6 = fig.add_subplot(236, projection=proj)
cs2=plt.pcolormesh(lon_rho_10km,lat_rho_10km,EKE_2km_interp, transform=ccrs.PlateCarree(), cmap=plt.cm.plasma, vmin=EKmin, vmax=EKmax)
plt.title('EKE \n WAOM2')
ax6.gridlines()
ax6.set_extent([-145, -68, -70, -67], crs=ccrs.PlateCarree())
x_left, x_right = ax6.get_xlim()
y_low, y_high = ax6.get_ylim()
#ax6.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax6.add_feature(cfeature.LAND, zorder=3, facecolor='white')
plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],zice_2km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
#plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],h_2km[xlimit,ylimit],levels=(2000,2001), transform=ccrs.PlateCarree(), colors='k', linewidths=1)
#plt.contour(lon_rho_2km[xlimit,ylimit], lat_rho_2km[xlimit,ylimit],h_2km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='grey')

cbar_ax1 = fig.add_axes([0.12, 0.08, 0.35, 0.01])
fig.colorbar(ct2, cax=cbar_ax1, orientation='horizontal')
cbar_ax1.set_xlabel('$^{\circ}$C')#, labelpad=-35)

cbar_ax2 = fig.add_axes([0.55, 0.08, 0.35, 0.01])
fig.colorbar(cs2, cax=cbar_ax2, orientation='horizontal')
cbar_ax2.set_xlabel('')#, labelpad=-35)

name_fig="waom10x4x2extend_shflim_S_0.25Q_MKE-EKE_surface_WAnt_anomalies.png"
plt.savefig(fig_path + name_fig, dpi=300)
plt.close()

