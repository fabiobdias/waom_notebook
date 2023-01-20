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

from netCDF4 import Dataset
from netCDF4 import num2date, date2num

import iris
import iris.iterate
import iris.coords
import iris.plot as iplt
import gsw


# read grid file for lon/lat coordinates
dg = xr.open_dataset("/scratch/project_2000339/boeiradi/waom10_frc/waom10extend_grd.nc")
lat_rho_10km= dg.variables["lat_rho"]
lon_rho_10km = dg.variables["lon_rho"]
dg.close()

dg4 = xr.open_dataset("/scratch/project_2000339/boeiradi/waom4_frc/waom4extend_grd.nc")
lat_rho_4km = dg4.variables["lat_rho"]
lon_rho_4km = dg4.variables["lon_rho"]
dg4.close()

dg2 = xr.open_dataset("/scratch/project_2000339/boeiradi/waom2_frc/waom2extend_grd.nc")
lat_rho_2km = dg2.variables["lat_rho"]
lon_rho_2km = dg2.variables["lon_rho"]
dg2.close()

# load ROMS avg output

def read_roms_ts_10km(exp_path):
    for yr  in ['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20']:
        ds = xr.open_dataset(exp_path + 'ocean_avg_00' + yr + '.nc')
        print('size temp and time length: ', ds.temp.shape, len(ds.salt.isel(xi_rho=20, eta_rho=100, s_rho=0)))
        temp_tmp = ds.variables["temp"]
        salt_tmp = ds.variables["salt"]
        zeta_tmp = ds.variables["zeta"]
        temp_tmp_ann = np.nanmean(temp_tmp, axis=0)
        salt_tmp_ann = np.nanmean(salt_tmp, axis=0)
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
    
        # concantenate annual averaged temp/salt
        if yr == '01':
            temp_ann = temp_tmp_ann
            salt_ann = salt_tmp_ann
            z_w_ann = z_w_tmp_ann
            z_rho_ann = z_rho_tmp_ann
        elif yr == '02':
            temp_ann = np.stack((temp_ann,temp_tmp_ann), axis=0)
            salt_ann = np.stack((salt_ann,salt_tmp_ann), axis=0)
            z_w_ann = np.stack((z_w_ann,z_w_tmp_ann), axis=0)
            z_rho_ann = np.stack((z_rho_ann,z_rho_tmp_ann), axis=0)
        else:
            temp_tmp_ann_4thdim = np.expand_dims(temp_tmp_ann, axis=0)
            temp_ann = np.concatenate((temp_ann,temp_tmp_ann_4thdim), axis=0)
            salt_tmp_ann_4thdim = np.expand_dims(salt_tmp_ann, axis=0)
            salt_ann = np.concatenate((salt_ann,salt_tmp_ann_4thdim), axis=0)
            z_w_tmp_ann_4thdim = np.expand_dims(z_w_tmp_ann, axis=0)
            z_w_ann = np.concatenate((z_w_ann,z_w_tmp_ann_4thdim), axis=0)
            z_rho_tmp_ann_4thdim = np.expand_dims(z_rho_tmp_ann, axis=0)
            z_rho_ann = np.concatenate((z_rho_ann,z_rho_tmp_ann_4thdim), axis=0)
            del temp_tmp_ann_4thdim, salt_tmp_ann_4thdim, z_w_tmp_ann_4thdim, z_rho_tmp_ann_4thdim
            
        print('Annual temp and annual tmp temp sizes = ', temp_ann.shape, temp_tmp_ann.shape)
        print('Annual z_w and annual tmp z_w sizes = ', z_w_ann.shape, z_w_tmp_ann.shape)
        
        del temp_tmp_ann, salt_tmp_ann, z_w_tmp_ann, z_rho_tmp_ann

    print('Annual temp, salt, z_w, z_rho sizes = ', temp_ann.shape, salt_ann.shape, z_w_ann.shape, z_rho_ann.shape)

    # shelf/open-ocean masks:
    mask_shelf = np.empty((ds.h.shape))

    shelf_ind=ds.h.where(ds.h < 2000)
    #print(open_ind)

    mask_shelf = np.divide(shelf_ind,shelf_ind)

    # calculate dz following:
    dz = np.empty((20,560,630,31))
    dz_inv = np.empty((20,560,630,31))

    for tt in np.arange(0,20):
        z_w_sorted = -1*z_w_ann[tt,:,::-1]
        #print(z_w_sorted.shape)
        dz_inv[tt,:,:,:] = np.diff(z_w_sorted,axis=2)
        dz[tt,:,:,:] = dz_inv[tt,:,:,::-1]

    print('size dz = ', dz.shape)

    # first, vertical integral:
    ohc_dz = np.empty((20,31,560,630))
    salt_dz = np.empty((20,31,560,630))

    cp_rho = 3989.245*1035 # J/kg/degC
    ohc_ann = cp_rho*temp_ann

    for tt in np.arange(0,20):
        for zz in np.arange(30,-1,-1):
            ohc_dz[tt,zz,:,:] = ohc_ann[tt,zz,:,:]*dz[tt,:,:,zz]
            salt_dz[tt,zz,:,:] = salt_ann[tt,zz,:,:]*dz[tt,:,:,zz]

    ohc_vertint = np.nansum(ohc_dz[:,::-1,:,:], axis=1)
    salt_vertint = np.nansum(salt_dz[:,::-1,:,:], axis=1)
    print('size ohc_ann, ohc_dz, ohc_vertint = ', ohc_ann.shape, ohc_dz.shape, ohc_vertint.shape)

    # horizontal integral
    p_area = ds.pm.isel(eta_rho=slice(0,560), xi_rho=slice(0,630))*ds.pn.isel(eta_rho=slice(0,560), xi_rho=slice(0,630))
    area = 1/p_area # area in meters

    print('size area = ', area.shape)

    ohc_area_shelf = np.empty((20,560,630))
    salt_area_shelf = np.empty((20,560,630))
    ohc_hint_shelf = np.empty((20,))
    salt_hint_shelf = np.empty((20,))
    ohc_hintn_shelf = np.empty((20,))
    salt_hintn_shelf = np.empty((20,))

    for tt in np.arange(0,20,1):
        ohc_area_shelf[tt,:,:] = ohc_vertint[tt,:,:]*area*mask_shelf
        salt_area_shelf[tt,:,:] =salt_vertint[tt,:,:]*area*mask_shelf
        ohc_hint_shelf[tt] = np.nansum(np.nansum((ohc_area_shelf[tt,:,:]), axis=1),axis=0)
        salt_hint_shelf[tt] = np.nansum(np.nansum((salt_area_shelf[tt,:,:]), axis=1),axis=0)
        area_sum_shelf = np.nansum(np.nansum(area*mask_shelf, axis=1),axis=0)

        ohc_hintn_shelf[tt] = np.divide(ohc_hint_shelf[tt],area_sum_shelf)
        salt_hintn_shelf[tt] = np.divide(salt_hint_shelf[tt],area_sum_shelf)

    SST_havg_shelf = np.empty((20,))
    SSS_havg_shelf = np.empty((20,))


    for tt in np.arange(0,20,1):
        SST_havg_shelf[tt] = np.nanmean(np.nanmean((temp_ann[tt,-1,:,:]*mask_shelf), axis=1),axis=0)
        SSS_havg_shelf[tt] = np.nanmean(np.nanmean((salt_ann[tt,-1,:,:]*mask_shelf), axis=1),axis=0)


    return SST_havg_shelf, SSS_havg_shelf, ohc_hintn_shelf, salt_hintn_shelf

# read ROMS 4km:
def read_roms_ts_4km(exp_path):
    for yr  in ['01','02','03','04','05','06','07','08','09','10']:
        ds = xr.open_dataset(exp_path + 'ocean_avg_00' + yr + '.nc')
        print('size temp and time length: ', ds.temp.shape, len(ds.salt.isel(xi_rho=20, eta_rho=100, s_rho=0)))
        temp_tmp = ds.variables["temp"]
        salt_tmp = ds.variables["salt"]
        zeta_tmp = ds.variables["zeta"]
        temp_tmp_ann = np.nanmean(temp_tmp, axis=0)
        salt_tmp_ann = np.nanmean(salt_tmp, axis=0)
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

        # concantenate annual averaged temp/salt
        if yr == '01':
            temp_ann = temp_tmp_ann
            salt_ann = salt_tmp_ann
            z_w_ann = z_w_tmp_ann
            z_rho_ann = z_rho_tmp_ann
        elif yr == '02':
            temp_ann = np.stack((temp_ann,temp_tmp_ann), axis=0)
            salt_ann = np.stack((salt_ann,salt_tmp_ann), axis=0)
            z_w_ann = np.stack((z_w_ann,z_w_tmp_ann), axis=0)
            z_rho_ann = np.stack((z_rho_ann,z_rho_tmp_ann), axis=0)
        else:
            temp_tmp_ann_4thdim = np.expand_dims(temp_tmp_ann, axis=0)
            temp_ann = np.concatenate((temp_ann,temp_tmp_ann_4thdim), axis=0)
            salt_tmp_ann_4thdim = np.expand_dims(salt_tmp_ann, axis=0)
            salt_ann = np.concatenate((salt_ann,salt_tmp_ann_4thdim), axis=0)
            z_w_tmp_ann_4thdim = np.expand_dims(z_w_tmp_ann, axis=0)
            z_w_ann = np.concatenate((z_w_ann,z_w_tmp_ann_4thdim), axis=0)
            z_rho_tmp_ann_4thdim = np.expand_dims(z_rho_tmp_ann, axis=0)
            z_rho_ann = np.concatenate((z_rho_ann,z_rho_tmp_ann_4thdim), axis=0)
            del temp_tmp_ann_4thdim, salt_tmp_ann_4thdim, z_w_tmp_ann_4thdim, z_rho_tmp_ann_4thdim

        print('Annual temp and annual tmp temp sizes = ', temp_ann.shape, temp_tmp_ann.shape)
        print('Annual z_w and annual tmp z_w sizes = ', z_w_ann.shape, z_w_tmp_ann.shape)

        del temp_tmp_ann, salt_tmp_ann, z_w_tmp_ann, z_rho_tmp_ann

    print('Annual temp, salt, z_w, z_rho sizes = ', temp_ann.shape, salt_ann.shape, z_w_ann.shape, z_rho_ann.shape)

    # shelf/open-ocean masks:
    mask_shelf = np.empty((ds.h.shape))

    shelf_ind=ds.h.where(ds.h < 2000)
    #print(open_ind)

    mask_shelf = np.divide(shelf_ind,shelf_ind)

    # calculate dz following:
    dz = np.empty((10,1400,1575,31))
    dz_inv = np.empty((10,1400,1575,31))

    for tt in np.arange(0,10):
        z_w_sorted = -1*z_w_ann[tt,:,::-1]
        dz_inv[tt,:,:,:] = np.diff(z_w_sorted,axis=2)
        dz[tt,:,:,:] = dz_inv[tt,:,:,::-1]

    print('size dz = ', dz.shape)

    # first, vertical integral:
    ohc_dz = np.empty((10,31,1400,1575))
    salt_dz = np.empty((10,31,1400,1575))

    cp_rho = 3989.245*1035 # J/kg/degC
    ohc_ann = cp_rho*temp_ann

    for tt in np.arange(0,10):
        for zz in np.arange(30,-1,-1):
            ohc_dz[tt,zz,:,:] = ohc_ann[tt,zz,:,:]*dz[tt,:,:,zz]
            salt_dz[tt,zz,:,:] = salt_ann[tt,zz,:,:]*dz[tt,:,:,zz]

    ohc_vertint = np.nansum(ohc_dz[:,::-1,:,:], axis=1)
    salt_vertint = np.nansum(salt_dz[:,::-1,:,:], axis=1)

    print('size ohc_ann, ohc_dz, ohc_vertint = ', ohc_ann.shape, ohc_dz.shape, ohc_vertint.shape)

    # horizontal integral
    p_area = ds.pm.isel(eta_rho=slice(0,1400), xi_rho=slice(0,1575))*ds.pn.isel(eta_rho=slice(0,1400), xi_rho=slice(0,1575))
    area = 1/p_area # area in meters

    print('size area = ', area.shape)

    ohc_area_shelf = np.empty((10,1400,1575))
    salt_area_shelf = np.empty((10,1400,1575))
    ohc_hint_shelf = np.empty((10,))
    salt_hint_shelf = np.empty((10,))
    ohc_hintn_shelf = np.empty((10,))
    salt_hintn_shelf = np.empty((10,))

    for tt in np.arange(0,10,1):
        ohc_area_shelf[tt,:,:] = ohc_vertint[tt,:,:]*area*mask_shelf
        salt_area_shelf[tt,:,:] =salt_vertint[tt,:,:]*area*mask_shelf

        ohc_hint_shelf[tt] = np.nansum(np.nansum((ohc_area_shelf[tt,:,:]), axis=1),axis=0)
        salt_hint_shelf[tt] = np.nansum(np.nansum((salt_area_shelf[tt,:,:]), axis=1),axis=0)
        area_sum_shelf = np.nansum(np.nansum(area*mask_shelf, axis=1),axis=0)

        ohc_hintn_shelf[tt] = np.divide(ohc_hint_shelf[tt],area_sum_shelf)
        salt_hintn_shelf[tt] = np.divide(salt_hint_shelf[tt],area_sum_shelf)

    SST_havg_shelf = np.empty((10,))
    SSS_havg_shelf = np.empty((10,))

    for tt in np.arange(0,10,1):
        SST_havg_shelf[tt] = np.nanmean(np.nanmean((temp_ann[tt,-1,:,:]*mask_shelf), axis=1),axis=0)
        SSS_havg_shelf[tt] = np.nanmean(np.nanmean((salt_ann[tt,-1,:,:]*mask_shelf), axis=1),axis=0)

    return SST_havg_shelf, SSS_havg_shelf, ohc_hintn_shelf, salt_hintn_shelf

# read ROMS 2km:
def read_roms_ts_2km(exp_path):
    for yr in ['01','02','03','04','05','06','07','08','09','10']:
        ds = xr.open_dataset(exp_path + 'ocean_avg_00' + yr + '.nc')
        print('size temp and time length: ', ds.temp.shape, len(ds.salt.isel(xi_rho=20, eta_rho=100, s_rho=0)))
        temp_tmp = ds.variables["temp"]
        salt_tmp = ds.variables["salt"]
        zeta_tmp = ds.variables["zeta"]
        temp_tmp_ann = np.nanmean(temp_tmp, axis=0)
        salt_tmp_ann = np.nanmean(salt_tmp, axis=0)
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

        # concantenate annual averaged temp/salt 
        if yr == '01':
            temp_ann = np.expand_dims(temp_tmp_ann, axis=0)
            salt_ann = np.expand_dims(salt_tmp_ann, axis=0)
            z_w_ann = np.expand_dims(z_w_tmp_ann, axis=0)
            z_rho_ann = np.expand_dims(z_rho_tmp_ann, axis=0)
        elif yr == '02':
            temp_ann = np.nanmean(np.stack((temp_ann,np.expand_dims(temp_tmp_ann, axis=0)), axis=0), axis=0)
            salt_ann = np.nanmean(np.stack((salt_ann,np.expand_dims(salt_tmp_ann, axis=0)), axis=0), axis=0)
            z_w_ann = np.nanmean(np.stack((z_w_ann,np.expand_dims(z_w_tmp_ann, axis=0)), axis=0), axis=0)
            z_rho_ann = np.nanmean(np.stack((z_rho_ann,np.expand_dims(z_rho_tmp_ann, axis=0)), axis=0), axis=0)
        elif yr == '03' or yr == '05' or yr == '07' or yr == '09':
            temp_tmp_a = temp_tmp_ann
            salt_tmp_a = salt_tmp_ann
            z_w_tmp_a = z_w_tmp_ann
            z_rho_tmp_a = z_rho_tmp_ann
        elif yr == '04' or yr == '06' or yr == '08' or yr == '10':
            temp_tmp_a_4thdim = np.expand_dims(temp_tmp_a, axis=0)
            temp_tmp_b_4thdim = np.expand_dims(temp_tmp_ann, axis=0)
            print('array sizes before concatenating (files 03/05/07 + 04/06/08):', temp_tmp_a_4thdim.shape, temp_tmp_b_4thdim.shape)
            temp_tmp_ann_4thdim = np.nanmean(np.stack((temp_tmp_a_4thdim,temp_tmp_b_4thdim), axis=0), axis=0)
            print(temp_ann.shape, temp_tmp_ann_4thdim.shape)
            temp_ann = np.concatenate((temp_ann,temp_tmp_ann_4thdim), axis=0)

            salt_tmp_a_4thdim = np.expand_dims(salt_tmp_a, axis=0)
            salt_tmp_b_4thdim = np.expand_dims(salt_tmp_ann, axis=0)
            salt_tmp_ann_4thdim = np.nanmean(np.stack((salt_tmp_a_4thdim,salt_tmp_b_4thdim), axis=0), axis=0)
            salt_ann = np.concatenate((salt_ann,salt_tmp_ann_4thdim), axis=0)

            z_w_tmp_a_4thdim = np.expand_dims(z_w_tmp_a, axis=0)
            z_w_tmp_b_4thdim = np.expand_dims(z_w_tmp_ann, axis=0)
            z_w_tmp_ann_4thdim = np.nanmean(np.stack((z_w_tmp_a_4thdim,z_w_tmp_b_4thdim), axis=0), axis=0)
            z_w_ann = np.concatenate((z_w_ann,z_w_tmp_ann_4thdim), axis=0)

            z_rho_tmp_a_4thdim = np.expand_dims(z_rho_tmp_a, axis=0)
            z_rho_tmp_b_4thdim = np.expand_dims(z_rho_tmp_ann, axis=0)
            z_rho_tmp_ann_4thdim = np.nanmean(np.stack((z_rho_tmp_a_4thdim,z_rho_tmp_b_4thdim), axis=0), axis=0)
            z_rho_ann = np.concatenate((z_rho_ann,z_rho_tmp_ann_4thdim), axis=0)

            del temp_tmp_a, salt_tmp_a, z_w_tmp_a, z_rho_tmp_a
            del temp_tmp_a_4thdim, temp_tmp_b_4thdim, salt_tmp_a_4thdim, salt_tmp_b_4thdim, z_w_tmp_a_4thdim, z_w_tmp_b_4thdim, z_rho_tmp_a_4thdim, z_rho_tmp_b_4thdim
            del temp_tmp_ann_4thdim, salt_tmp_ann_4thdim, z_w_tmp_ann_4thdim, z_rho_tmp_ann_4thdim

        print('Annual temp and annual tmp temp sizes = ', temp_ann.shape, temp_tmp_ann.shape)
        print('Annual z_w and annual tmp z_w sizes = ', z_w_ann.shape, z_w_tmp_ann.shape)

        del temp_tmp_ann, salt_tmp_ann, z_w_tmp_ann, z_rho_tmp_ann

    print('Annual temp, salt, z_w, z_rho sizes = ', temp_ann.shape, salt_ann.shape, z_w_ann.shape, z_rho_ann.shape)

    # shelf/open-ocean masks:
    mask_shelf = np.empty((ds.h.shape))

    shelf_ind=ds.h.where(ds.h < 2000)

    mask_shelf = np.divide(shelf_ind,shelf_ind)

    # calculate dz following:
    dz = np.empty((5,2800,3150,31))
    dz_inv = np.empty((5,2800,3150,31))

    for tt in np.arange(0,5):
        z_w_sorted = -1*z_w_ann[tt,:,::-1]
        #print(z_w_sorted.shape)
        dz_inv[tt,:,:,:] = np.diff(z_w_sorted,axis=2)
        dz[tt,:,:,:] = dz_inv[tt,:,:,::-1]

    print('size dz = ', dz.shape)

    # first, vertical integral:
    ohc_dz = np.empty((5,31,2800,3150))
    salt_dz = np.empty((5,31,2800,3150))

    cp_rho = 3989.245*1035 # J/kg/degC
    ohc_ann = cp_rho*temp_ann

    for tt in np.arange(0,5):
        for zz in np.arange(30,-1,-1):
            ohc_dz[tt,zz,:,:] = ohc_ann[tt,zz,:,:]*dz[tt,:,:,zz]
            salt_dz[tt,zz,:,:] = salt_ann[tt,zz,:,:]*dz[tt,:,:,zz]

    ohc_vertint = np.nansum(ohc_dz[:,::-1,:,:], axis=1)
    salt_vertint = np.nansum(salt_dz[:,::-1,:,:], axis=1)

    print('size ohc_ann, ohc_dz, ohc_vertint = ', ohc_ann.shape, ohc_dz.shape, ohc_vertint.shape)

    # horizontal integral
    p_area = ds.pm.isel(eta_rho=slice(0,2800), xi_rho=slice(0,3150))*ds.pn.isel(eta_rho=slice(0,2800), xi_rho=slice(0,3150))
    area = 1/p_area # area in meters

    print('size area = ', area.shape)

    ohc_area_shelf = np.empty((5,2800,3150))
    salt_area_shelf = np.empty((5,2800,3150))
    ohc_hint = np.empty((5,))
    ohc_hint_shelf = np.empty((5,))
    salt_hint_shelf = np.empty((5,))
    ohc_hintn_shelf = np.empty((5,))
    salt_hintn_shelf = np.empty((5,))

    for tt in np.arange(0,5,1):
        ohc_area_shelf[tt,:,:] = ohc_vertint[tt,:,:]*area*mask_shelf
        salt_area_shelf[tt,:,:] =salt_vertint[tt,:,:]*area*mask_shelf

        ohc_hint_shelf[tt] = np.nansum(np.nansum((ohc_area_shelf[tt,:,:]), axis=1),axis=0)
        salt_hint_shelf[tt] = np.nansum(np.nansum((salt_area_shelf[tt,:,:]), axis=1),axis=0)
        area_sum_shelf = np.nansum(np.nansum(area*mask_shelf, axis=1),axis=0)

        ohc_hintn_shelf[tt] = np.divide(ohc_hint_shelf[tt],area_sum_shelf)
        salt_hintn_shelf[tt] = np.divide(salt_hint_shelf[tt],area_sum_shelf)

    SST_havg_shelf = np.empty((5,))
    SSS_havg_shelf = np.empty((5,))

    for tt in np.arange(0,5,1):
        SST_havg_shelf[tt] = np.nanmean(np.nanmean((temp_ann[tt,-1,:,:]), axis=1),axis=0)
        SSS_havg_shelf[tt] = np.nanmean(np.nanmean((salt_ann[tt,-1,:,:]), axis=1),axis=0)

    return SST_havg_shelf, SSS_havg_shelf, ohc_hintn_shelf, salt_hintn_shelf


path_ECCO2_10km = '/scratch/project_2000789/boeiradi/waom10extend_shflim_S_0.25Q/output_01-20yr/'
path_ECCO2_4km = '/scratch/project_2000789/boeiradi/waom4extend_shflim_S_0.25Q/output_01-10yr/'
path_ECCO2_2km = '/scratch/project_2000339/boeiradi/waom2extend_shflim_S_0.25Q/output_01-05yr/'

E10_SST_havg_shelf, E10_SSS_havg_shelf, E10_ohc_hintn_shelf, E10_salt_hintn_shelf = read_roms_ts_10km(path_ECCO2_10km)
E4_SST_havg_shelf, E4_SSS_havg_shelf, E4_ohc_hintn_shelf, E4_salt_hintn_shelf = read_roms_ts_4km(path_ECCO2_4km)
E2_SST_havg_shelf, E2_SSS_havg_shelf, E2_ohc_hintn_shelf, E2_salt_hintn_shelf = read_roms_ts_2km(path_ECCO2_2km)

fig_path='/users/boeiradi/COLD_project/postprocessing/figs/Grid_avg/'

fig = plt.figure(figsize=(12,8))
ax1 = fig.add_subplot(411)
plt.plot(np.arange(0,20)+1,E10_SST_havg_shelf,'--r',label='10km')
plt.plot(np.arange(10,20)+1,E4_SST_havg_shelf,'--b',label='4km')
plt.plot(np.arange(15,20)+1,E2_SST_havg_shelf,'--g',label='2km')
l1 = plt.legend(loc='lower left', borderaxespad=0)
plt.ylabel('SST')
ax2 = fig.add_subplot(412)
plt.plot(np.arange(0,20)+1,E10_ohc_hintn_shelf,'--r',label='10km')
plt.plot(np.arange(10,20)+1,E4_ohc_hintn_shelf,'--b',label='4km')
plt.plot(np.arange(15,20)+1,E2_ohc_hintn_shelf,'--g',label='2km')
plt.ylabel('heat content (J)')
ax3 = fig.add_subplot(413)
plt.plot(np.arange(0,20)+1,E10_SSS_havg_shelf,'--r',label='10km')
plt.plot(np.arange(10,20)+1,E4_SSS_havg_shelf,'--b',label='4km')
plt.plot(np.arange(15,20)+1,E2_SSS_havg_shelf,'--g',label='2km')
plt.ylabel('SSS')
ax4 = fig.add_subplot(414)
plt.plot(np.arange(0,20)+1,-E10_salt_hintn_shelf,'--r',label='10km')
plt.plot(np.arange(10,20)+1,-E4_salt_hintn_shelf,'--b',label='4km')
plt.plot(np.arange(15,20)+1,-E2_salt_hintn_shelf,'--g',label='2km')
plt.ylabel('salt content')
plt.xlabel('years')

name_fig="waom10x4x2extend_shflim_S_0.25Q_OHC+Salt+SST+SSS_evolution_shelf.png"
plt.savefig(fig_path + name_fig, dpi=300)
plt.close()

