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
dg = xr.open_dataset("/scratch/project_2000339/boeiradi/waom10_frc/waom10_grd.nc")

lat_rho_10km= dg.variables["lat_rho"]
lon_rho_10km = dg.variables["lon_rho"]
#ds.coords['lat_rho']=lat_rho#.transpose() # put lat_rho into ds dataset
#ds.coords['lon_rho']=lon_rho#.transpose() # put lon_rho into ds dataset

dg4 = xr.open_dataset("/scratch/project_2000339/boeiradi/waom4_frc/waom4_grd.nc")

lat_rho_4km = dg4.variables["lat_rho"]
lon_rho_4km = dg4.variables["lon_rho"]


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
    
        # calculate geostrophic velocities from ssh:
        cor = 2*(2*np.pi/86400)*np.sin(lat_rho_10km[0:-1,0:-1])
        g = 9.81
        dzeta_dy = np.diff(zeta_tmp, axis=1)
        dzeta_dx = np.diff(zeta_tmp, axis=2)
        #print("Sizes of cor, dzeta_dy, dzeta_dx = ", cor.shape, dzeta_dy[:,:,0:-1].shape, dzeta_dx[:,0:-1,:].shape )
    
        ugeo = np.empty((dzeta_dy[:,:,0:-1].shape))
        vgeo = np.empty((dzeta_dy[:,:,0:-1].shape))
        mmax = 12
        for mm in np.arange(0,mmax,1):
            ugeo[mm,:,:] = -g/cor*dzeta_dy[mm,:,0:-1]
            vgeo[mm,:,:] = g/cor*dzeta_dx[mm,0:-1,:]
        ugeo_ann_tmp = np.nanmean(ugeo, axis=0)
        vgeo_ann_tmp = np.nanmean(vgeo, axis=0)
        #ugeo_anom_tmp = ugeo - ugeo_ann
        #vgeo_anom_tmp = vgeo - vgeo_ann
        V_ann_tmp = np.sqrt(ugeo_ann_tmp*ugeo_ann_tmp + vgeo_ann_tmp*vgeo_ann_tmp)
        MKE_tmp = np.divide(V_ann_tmp*V_ann_tmp,2)
        print("Size MKE_tmp", MKE_tmp.shape)
        
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
            MKE = MKE_tmp
        elif yr == '02':
            temp_ann = np.stack((temp_ann,temp_tmp_ann), axis=0)
            salt_ann = np.stack((salt_ann,salt_tmp_ann), axis=0)
            z_w_ann = np.stack((z_w_ann,z_w_tmp_ann), axis=0)
            z_rho_ann = np.stack((z_rho_ann,z_rho_tmp_ann), axis=0)
            MKE = np.stack((MKE,MKE_tmp), axis=0)
        else:
            temp_tmp_ann_4thdim = np.expand_dims(temp_tmp_ann, axis=0)
            temp_ann = np.concatenate((temp_ann,temp_tmp_ann_4thdim), axis=0)
            salt_tmp_ann_4thdim = np.expand_dims(salt_tmp_ann, axis=0)
            salt_ann = np.concatenate((salt_ann,salt_tmp_ann_4thdim), axis=0)
            z_w_tmp_ann_4thdim = np.expand_dims(z_w_tmp_ann, axis=0)
            z_w_ann = np.concatenate((z_w_ann,z_w_tmp_ann_4thdim), axis=0)
            z_rho_tmp_ann_4thdim = np.expand_dims(z_rho_tmp_ann, axis=0)
            z_rho_ann = np.concatenate((z_rho_ann,z_rho_tmp_ann_4thdim), axis=0)
            MKE_tmp_4thdim = np.expand_dims(MKE_tmp, axis=0)
            MKE = np.concatenate((MKE,MKE_tmp_4thdim), axis=0)
            del temp_tmp_ann_4thdim, salt_tmp_ann_4thdim, z_w_tmp_ann_4thdim, z_rho_tmp_ann_4thdim, MKE_tmp_4thdim
            
        print('Annual temp and annual tmp temp sizes = ', temp_ann.shape, temp_tmp_ann.shape)
        print('Annual z_w and annual tmp z_w sizes = ', z_w_ann.shape, z_w_tmp_ann.shape)
        
        del temp_tmp_ann, salt_tmp_ann, z_w_tmp_ann, z_rho_tmp_ann

    print('Annual temp, salt, z_w, z_rho, MKE sizes = ', temp_ann.shape, salt_ann.shape, z_w_ann.shape, z_rho_ann.shape, MKE.shape)

    # shelf/open-ocean masks:
    mask_shelf = np.empty((ds.h.shape))
    mask_open = np.empty((ds.h.shape))

    open_ind=ds.h.where(ds.h >= 1000)
    shelf_ind=ds.h.where(ds.h < 1000)
    #print(open_ind)

    mask_shelf = np.divide(shelf_ind,shelf_ind)
    mask_open = np.divide(open_ind,open_ind)

    # calculate dz following:
    dz = np.empty((20,530,630,31))
    dz_inv = np.empty((20,530,630,31))

    for tt in np.arange(0,20):
        z_w_sorted = -1*z_w_ann[tt,:,::-1]
        #print(z_w_sorted.shape)
        dz_inv[tt,:,:,:] = np.diff(z_w_sorted,axis=2)
        dz[tt,:,:,:] = dz_inv[tt,:,:,::-1]

    print('size dz = ', dz.shape)

    # first, vertical integral:
    ohc_dz = np.empty((20,31,530,630))
    salt_dz = np.empty((20,31,530,630))

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
    p_area = ds.pm.isel(eta_rho=slice(0,530), xi_rho=slice(0,630))*ds.pn.isel(eta_rho=slice(0,530), xi_rho=slice(0,630))
    area = 1/p_area # area in meters

    print('size area = ', area.shape)

    ohc_area = np.empty((20,530,630))
    salt_area = np.empty((20,530,630))
    ohc_area_shelf = np.empty((20,530,630))
    salt_area_shelf = np.empty((20,530,630))
    ohc_area_open = np.empty((20,530,630))
    salt_area_open = np.empty((20,530,630))
    ohc_hint = np.empty((20,))
    salt_hint = np.empty((20,))
    ohc_hintn = np.empty((20,))
    salt_hintn = np.empty((20,))
    ohc_hint_shelf = np.empty((20,))
    salt_hint_shelf = np.empty((20,))
    ohc_hintn_shelf = np.empty((20,))
    salt_hintn_shelf = np.empty((20,))
    ohc_hint_open = np.empty((20,))
    salt_hint_open = np.empty((20,))
    ohc_hintn_open = np.empty((20,))
    salt_hintn_open = np.empty((20,))

    for tt in np.arange(0,20,1):
        ohc_area[tt,:,:] = ohc_vertint[tt,:,:]*area
        salt_area[tt,:,:] =salt_vertint[tt,:,:]*area
        print('size ohc_area = ', ohc_area[tt,:,:].shape)
        ohc_area_shelf[tt,:,:] = ohc_vertint[tt,:,:]*area*mask_shelf
        salt_area_shelf[tt,:,:] =salt_vertint[tt,:,:]*area*mask_shelf
        ohc_area_open[tt,:,:] = ohc_vertint[tt,:,:]*area*mask_open
        salt_area_open[tt,:,:] =salt_vertint[tt,:,:]*area*mask_open

        ohc_hint[tt] = np.nansum(np.nansum((ohc_area[tt,:,:]), axis=1),axis=0)
        salt_hint[tt] = np.nansum(np.nansum((salt_area[tt,:,:]), axis=1),axis=0)
        area_sum = np.nansum(np.nansum(area, axis=1),axis=0)

        ohc_hint_shelf[tt] = np.nansum(np.nansum((ohc_area_shelf[tt,:,:]), axis=1),axis=0)
        salt_hint_shelf[tt] = np.nansum(np.nansum((salt_area_shelf[tt,:,:]), axis=1),axis=0)
        area_sum_shelf = np.nansum(np.nansum(area*mask_shelf, axis=1),axis=0)

        ohc_hint_open[tt] = np.nansum(np.nansum((ohc_area_open[tt,:,:]), axis=1),axis=0)
        salt_hint_open[tt] = np.nansum(np.nansum((salt_area_open[tt,:,:]), axis=1),axis=0)
        area_sum_open = np.nansum(np.nansum(area*mask_open, axis=1),axis=0)

        ohc_hintn[tt] = np.divide(ohc_hint[tt],area_sum)
        salt_hintn[tt] = np.divide(salt_hint[tt],area_sum)

        ohc_hintn_shelf[tt] = np.divide(ohc_hint_shelf[tt],area_sum_shelf)
        salt_hintn_shelf[tt] = np.divide(salt_hint_shelf[tt],area_sum_shelf)

        ohc_hintn_open[tt] = np.divide(ohc_hint_open[tt],area_sum_open)
        salt_hintn_open[tt] = np.divide(salt_hint_open[tt],area_sum_open)

    # horizontal-integral of the surface fluxes

    # from forcing file:
    #shflx_ann = dshf.shflux.mean(axis=0)
    #sfwflx_ann = dfw.swflux.mean(axis=0)

    # from model diagnostic:
    shflx_ann = ds.shflux.mean(axis=0)
    sfwflx_ann = ds.ssflux.mean(axis=0)

    shflx_area = np.empty((530,630))
    sfwflx_area = np.empty((530,630))

    shflx_area[:,:] = shflx_ann[:,:]*area
    sfwflx_area[:,:] =sfwflx_ann[:,:]*area

    shflx_hint = np.nansum(np.nansum((shflx_area[:,:]), axis=1),axis=0)
    sfwflx_hint = np.nansum(np.nansum((sfwflx_area[:,:]), axis=1),axis=0)

    shflx_hintn = np.divide(shflx_hint,area_sum)
    sfwflx_hintn = np.divide(sfwflx_hint,area_sum)

    print('Horizontally-averaged surface heat flux = ', shflx_hintn)
    print('Horizontally-averaged surface FW flux = ', sfwflx_hintn)

    # horizontally-average Mean Kinect Energy (MKE)

    MKE_havg = np.empty((20,))
    SST_havg = np.empty((20,))
    SSS_havg = np.empty((20,))

    print(MKE.shape)

    for tt in np.arange(0,20,1):
        MKE_havg[tt] = np.nanmean(np.nanmean((MKE[tt,:,:]), axis=1),axis=0)
        SST_havg[tt] = np.nanmean(np.nanmean((temp_ann[tt,-1,:,:]), axis=1),axis=0)
        SSS_havg[tt] = np.nanmean(np.nanmean((salt_ann[tt,-1,:,:]), axis=1),axis=0)

    print(ohc_hintn, SST_havg)

    return ohc_hintn, salt_hintn, shflx_hintn, sfwflx_hintn, MKE_havg, SST_havg, SSS_havg, ohc_hintn_shelf, ohc_hintn_open, salt_hintn_shelf, salt_hintn_open

# read ROMS 4km:


def read_roms_ts_4km(exp_path):
    for yr in ['output_1yr','output_2yr','output_3yr','output_4yr','output_5yr','output_6yr','output_7yr','output_8yr','output_9yr','output_10yr']:
        for mm  in ['01','02','03','04','05','06','07','08','09','10','11','12']:
            ds = xr.open_dataset(exp_path + yr + '/ocean_avg_00' + mm + '.nc')
            print('MM: size temp and time length: ', yr, mm, ds.temp.shape, len(ds.salt.isel(xi_rho=20, eta_rho=100, s_rho=0)))
            temp_tmp = ds.variables["temp"]
            salt_tmp = ds.variables["salt"]
            zeta_tmp = ds.variables["zeta"]
            temp_mavg = np.nanmean(temp_tmp, axis=0) # monthly average
            salt_mavg = np.nanmean(salt_tmp, axis=0)
            zeta_mavg = np.nanmean(zeta_tmp, axis=0)
            print('size temp_mavg = ', temp_mavg.shape)
            del temp_tmp, salt_tmp, zeta_tmp

            # calculate geostrophic velocities from ssh:
            cor = 2*(2*np.pi/86400)*np.sin(lat_rho_4km[0:-1,0:-1])
            g = 9.81
            dzeta_dy = np.diff(zeta_mavg, axis=0)
            dzeta_dx = np.diff(zeta_mavg, axis=1)
            print("Sizes of cor, dzeta_dy, dzeta_dx = ", cor.shape, dzeta_dy[:,0:-1].shape, dzeta_dx[0:-1,:].shape )
            ugeo_mavg = -g/cor*dzeta_dy[:,0:-1]
            vgeo_mavg = g/cor*dzeta_dx[0:-1,:]
            V_mavg = np.sqrt(ugeo_mavg*ugeo_mavg + vgeo_mavg*vgeo_mavg)
            MKE_mavg = np.divide(V_mavg*V_mavg,2)
            print("Size MKE_mavg", MKE_mavg.shape)

            ds = ds.set_coords(['Cs_r', 'Cs_w', 'hc', 'h', 'Vtransform'])

            if ds.Vtransform == 1:
                Zo_rho = ds.hc * (ds.s_rho - ds.Cs_r) + ds.Cs_r * ds.h
                z_rho_tmp = Zo_rho + ds.zeta * (1 + Zo_rho/ds.h)
                print("Vtransform=1")
            elif ds.Vtransform == 2:
                Zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * ds.h) / (ds.hc + ds.h)
                z_rho_tmp = ds.zeta + (ds.zeta + ds.h) * Zo_rho
                print("Vtransform=2")
            z_rho_mavg = np.nanmean(z_rho_tmp, axis=0)

            Zo_w = (ds.hc * ds.s_w + ds.Cs_w * ds.h) / (ds.hc + ds.h)
            z_w_tmp = ds.zeta + (ds.zeta + ds.h) * Zo_w
            z_w_mavg = np.nanmean(z_w_tmp, axis=0)

	# concatenate monthly averages into annual vars:
            if mm == '01':
                temp_ann = temp_mavg
                salt_ann = salt_mavg
                z_w_ann = z_w_mavg
                z_rho_ann = z_rho_mavg
                MKE_ann = MKE_mavg
            elif mm == '02':
                temp_ann = np.stack((temp_ann,temp_mavg), axis=0)
                salt_ann = np.stack((salt_ann,salt_mavg), axis=0)
                z_w_ann = np.stack((z_w_ann,z_w_mavg), axis=0)
                z_rho_ann = np.stack((z_rho_ann,z_rho_mavg), axis=0)
                MKE_ann = np.stack((MKE_ann,MKE_mavg), axis=0)
            else:
                temp_mavg_4thdim = np.expand_dims(temp_mavg, axis=0)
                temp_ann = np.concatenate((temp_ann,temp_mavg_4thdim), axis=0)
                salt_mavg_4thdim = np.expand_dims(salt_mavg, axis=0)
                salt_ann = np.concatenate((salt_ann,salt_mavg_4thdim), axis=0)
                z_w_mavg_4thdim = np.expand_dims(z_w_mavg, axis=0)
                z_w_ann = np.concatenate((z_w_ann,z_w_mavg_4thdim), axis=0)
                z_rho_mavg_4thdim = np.expand_dims(z_rho_mavg, axis=0)
                z_rho_ann = np.concatenate((z_rho_ann,z_rho_mavg_4thdim), axis=0)
                MKE_tmp_4thdim = np.expand_dims(MKE_mavg, axis=0)
                MKE_ann = np.concatenate((MKE_ann,MKE_tmp_4thdim), axis=0)
                del temp_mavg_4thdim, salt_mavg_4thdim, z_w_mavg_4thdim, z_rho_mavg_4thdim, MKE_tmp_4thdim
            del temp_mavg, salt_mavg, z_w_mavg, z_rho_mavg, MKE_mavg

        # do annual mean of annual vars:
        temp_yavg = np.nanmean(temp_ann, axis=0)
        salt_yavg = np.nanmean(salt_ann, axis=0)
        z_w_yavg = np.nanmean(z_w_ann, axis=0)
        z_rho_yavg = np.nanmean(z_rho_ann, axis=0)
        MKE_yavg = np.nanmean(MKE_ann, axis=0)

        print('Annual-avg temp, annual tmp temp sizes = ', temp_yavg.shape, temp_ann.shape)

        del temp_ann, salt_ann, z_w_ann, z_rho_ann, MKE_ann

        # concantenate annual averaged temp/salt
        if yr == 'output_1yr':
            temp_yrly = temp_yavg
            salt_yrly = salt_yavg
            z_w_yrly = z_w_yavg
            z_rho_yrly = z_rho_yavg
            MKE_yrly = MKE_yavg
        elif yr == 'output_2yr':
            temp_yrly = np.stack((temp_yrly,temp_yavg), axis=0)
            salt_yrly = np.stack((salt_yrly,salt_yavg), axis=0)
            z_w_yrly = np.stack((z_w_yrly,z_w_yavg), axis=0)
            z_rho_yrly = np.stack((z_rho_yrly,z_rho_yavg), axis=0)
            MKE_yrly = np.stack((MKE_yrly,MKE_yavg), axis=0)
        else:
            temp_tmp_ann_4thdim = np.expand_dims(temp_yavg, axis=0)
            temp_yrly = np.concatenate((temp_yrly,temp_tmp_ann_4thdim), axis=0)
            salt_tmp_ann_4thdim = np.expand_dims(salt_yavg, axis=0)
            salt_yrly = np.concatenate((salt_yrly,salt_tmp_ann_4thdim), axis=0)
            z_w_tmp_ann_4thdim = np.expand_dims(z_w_yavg, axis=0)
            z_w_yrly = np.concatenate((z_w_yrly,z_w_tmp_ann_4thdim), axis=0)
            z_rho_tmp_ann_4thdim = np.expand_dims(z_rho_yavg, axis=0)
            z_rho_yrly = np.concatenate((z_rho_yrly,z_rho_tmp_ann_4thdim), axis=0)
            MKE_tmp_4thdim = np.expand_dims(MKE_yavg, axis=0)
            MKE_yrly = np.concatenate((MKE_yrly,MKE_tmp_4thdim), axis=0)
            del temp_tmp_ann_4thdim, salt_tmp_ann_4thdim, z_w_tmp_ann_4thdim, z_rho_tmp_ann_4thdim, MKE_tmp_4thdim

        print('Annual temp and annual tmp temp sizes = ', temp_yrly.shape, temp_yavg.shape)
        print('Annual z_w and annual tmp z_w sizes = ', z_w_yrly.shape, z_w_yavg.shape)

        del temp_yavg, salt_yavg, z_w_yavg, z_rho_yavg, MKE_yavg

    print('Annual temp, salt, z_w, z_rho, MKE sizes = ', temp_yrly.shape, salt_yrly.shape, z_w_yrly.shape, z_rho_yrly.shape, MKE_yrly.shape)

    # shelf/open-ocean masks:
    mask_shelf = np.empty((ds.h.shape))
    mask_open = np.empty((ds.h.shape))

    open_ind=ds.h.where(ds.h >= 1000)
    shelf_ind=ds.h.where(ds.h < 1000)
    #print(open_ind)

    mask_shelf = np.divide(shelf_ind,shelf_ind)
    mask_open = np.divide(open_ind,open_ind)

    # calculate dz following:
    dz = np.empty((10,1325,1575,31))
    dz_inv = np.empty((10,1325,1575,31))

    for tt in np.arange(0,10):
        z_w_sorted = -1*z_w_yrly[tt,:,::-1]
        #print(z_w_sorted.shape)
        dz_inv[tt,:,:,:] = np.diff(z_w_sorted,axis=2)
        dz[tt,:,:,:] = dz_inv[tt,:,:,::-1]

    print('size dz = ', dz.shape)

    # first, vertical integral:
    ohc_dz = np.empty((10,31,1325,1575))
    salt_dz = np.empty((10,31,1325,1575))

    cp_rho = 3989.245*1035 # J/kg/degC
    ohc_yrly = cp_rho*temp_yrly

    for tt in np.arange(0,10):
        for zz in np.arange(30,-1,-1):
            ohc_dz[tt,zz,:,:] = ohc_yrly[tt,zz,:,:]*dz[tt,:,:,zz]
            salt_dz[tt,zz,:,:] = salt_yrly[tt,zz,:,:]*dz[tt,:,:,zz]

    ohc_vertint = np.nansum(ohc_dz[:,::-1,:,:], axis=1)
    salt_vertint = np.nansum(salt_dz[:,::-1,:,:], axis=1)

    print('size ohc_ann, ohc_dz, ohc_vertint = ', ohc_yrly.shape, ohc_dz.shape, ohc_vertint.shape)

    # horizontal integral
    p_area = ds.pm.isel(eta_rho=slice(0,1325), xi_rho=slice(0,1575))*ds.pn.isel(eta_rho=slice(0,1325), xi_rho=slice(0,1575))
    area = 1/p_area # area in meters

    print('size area = ', area.shape)

    ohc_area = np.empty((10,1325,1575))
    salt_area = np.empty((10,1325,1575))
    ohc_area_shelf = np.empty((10,1325,1575))
    salt_area_shelf = np.empty((10,1325,1575))
    ohc_area_open = np.empty((10,1325,1575))
    salt_area_open = np.empty((10,1325,1575))
    ohc_hint = np.empty((10,))
    salt_hint = np.empty((10,))
    ohc_hintn = np.empty((10,))
    salt_hintn = np.empty((10,))
    ohc_hint_shelf = np.empty((10,))
    salt_hint_shelf = np.empty((10,))
    ohc_hintn_shelf = np.empty((10,))
    salt_hintn_shelf = np.empty((10,))
    ohc_hint_open = np.empty((10,))
    salt_hint_open = np.empty((10,))
    ohc_hintn_open = np.empty((10,))
    salt_hintn_open = np.empty((10,))

    for tt in np.arange(0,10,1):
        ohc_area[tt,:,:] = ohc_vertint[tt,:,:]*area
        salt_area[tt,:,:] =salt_vertint[tt,:,:]*area
        print('size ohc_area = ', ohc_area[tt,:,:].shape)
        ohc_area_shelf[tt,:,:] = ohc_vertint[tt,:,:]*area*mask_shelf
        salt_area_shelf[tt,:,:] =salt_vertint[tt,:,:]*area*mask_shelf
        ohc_area_open[tt,:,:] = ohc_vertint[tt,:,:]*area*mask_open
        salt_area_open[tt,:,:] =salt_vertint[tt,:,:]*area*mask_open

        ohc_hint[tt] = np.nansum(np.nansum((ohc_area[tt,:,:]), axis=1),axis=0)
        salt_hint[tt] = np.nansum(np.nansum((salt_area[tt,:,:]), axis=1),axis=0)
        area_sum = np.nansum(np.nansum(area, axis=1),axis=0)

        ohc_hint_shelf[tt] = np.nansum(np.nansum((ohc_area_shelf[tt,:,:]), axis=1),axis=0)
        salt_hint_shelf[tt] = np.nansum(np.nansum((salt_area_shelf[tt,:,:]), axis=1),axis=0)
        area_sum_shelf = np.nansum(np.nansum(area*mask_shelf, axis=1),axis=0)

        ohc_hint_open[tt] = np.nansum(np.nansum((ohc_area_open[tt,:,:]), axis=1),axis=0)
        salt_hint_open[tt] = np.nansum(np.nansum((salt_area_open[tt,:,:]), axis=1),axis=0)
        area_sum_open = np.nansum(np.nansum(area*mask_open, axis=1),axis=0)

        ohc_hintn[tt] = np.divide(ohc_hint[tt],area_sum)
        salt_hintn[tt] = np.divide(salt_hint[tt],area_sum)

        ohc_hintn_shelf[tt] = np.divide(ohc_hint_shelf[tt],area_sum_shelf)
        salt_hintn_shelf[tt] = np.divide(salt_hint_shelf[tt],area_sum_shelf)

        ohc_hintn_open[tt] = np.divide(ohc_hint_open[tt],area_sum_open)
        salt_hintn_open[tt] = np.divide(salt_hint_open[tt],area_sum_open)

    # horizontal-integral of the surface fluxes

    # from forcing file:
    #shflx_ann = dshf.shflux.mean(axis=0)
    #sfwflx_ann = dfw.swflux.mean(axis=0)

    # from model diagnostic:
    shflx_ann = ds.shflux.mean(axis=0)
    sfwflx_ann = ds.ssflux.mean(axis=0)

    shflx_area = np.empty((1325,1575))
    sfwflx_area = np.empty((1325,1575))

    shflx_area[:,:] = shflx_ann[:,:]*area
    sfwflx_area[:,:] =sfwflx_ann[:,:]*area

    shflx_hint = np.nansum(np.nansum((shflx_area[:,:]), axis=1),axis=0)
    sfwflx_hint = np.nansum(np.nansum((sfwflx_area[:,:]), axis=1),axis=0)

    shflx_hintn = np.divide(shflx_hint,area_sum)
    sfwflx_hintn = np.divide(sfwflx_hint,area_sum)

    print('Horizontally-averaged surface heat flux = ', shflx_hintn)
    print('Horizontally-averaged surface FW flux = ', sfwflx_hintn)

    # horizontally-average Mean Kinect Energy (MKE)

    MKE_havg = np.empty((10,))
    SST_havg = np.empty((10,))
    SSS_havg = np.empty((10,))

    print(MKE_yrly.shape)

    for tt in np.arange(0,10,1):
        MKE_havg[tt] = np.nanmean(np.nanmean((MKE_yrly[tt,:,:]), axis=1),axis=0)
        SST_havg[tt] = np.nanmean(np.nanmean((temp_yrly[tt,-1,:,:]), axis=1),axis=0)
        SSS_havg[tt] = np.nanmean(np.nanmean((salt_yrly[tt,-1,:,:]), axis=1),axis=0)

    print(ohc_hintn, SST_havg)

    return ohc_hintn, salt_hintn, shflx_hintn, sfwflx_hintn, MKE_havg, SST_havg, SSS_havg, ohc_hintn_shelf, ohc_hintn_open, salt_hintn_shelf, salt_hintn_open

path_ORAS5_10km = '/scratch/project_2000339/boeiradi/waom10_shflim_S_ORAS5em_0.25Q/output_01-20yr/'
path_ECCO2_10km = '/scratch/project_2000339/boeiradi/waom10_shflim_S_0.25Q/output_01-20yr/'
path_ORAS5_4km = '/scratch/project_2000789/boeiradi/waom4_shflim_S_ORAS5em_0.25Q/'

O10_ohc_hintn, O10_salt_hintn, O10_shflx_hintn, O10_sfwflx_hintn, O10_MKE_havg, O10_SST_havg, O10_SSS_havg, O10_ohc_hintn_shelf, O10_ohc_hintn_open, O10_salt_hintn_shelf, O10_salt_hintn_open = read_roms_ts_10km(path_ORAS5_10km)
E10_ohc_hintn, E10_salt_hintn, E10_shflx_hintn, E10_sfwflx_hintn, E10_MKE_havg, E10_SST_havg, E10_SSS_havg, E10_ohc_hintn_shelf, E10_ohc_hintn_open, E10_salt_hintn_shelf, E10_salt_hintn_open = read_roms_ts_10km(path_ECCO2_10km)
O4_ohc_hintn, O4_salt_hintn, O4_shflx_hintn, O4_sfwflx_hintn, O4_MKE_havg, O4_SST_havg, O4_SSS_havg, O4_ohc_hintn_shelf, O4_ohc_hintn_open, O4_salt_hintn_shelf, O4_salt_hintn_open = read_roms_ts_4km(path_ORAS5_4km)

fig_path='/users/boeiradi/COLD_project/postprocessing/figs/Grid_avg/'

fig = plt.figure(figsize=(12,8))
ax1 = fig.add_subplot(311)
plt.plot(np.arange(0,20)+1,E10_ohc_hintn,'r',label='ECCO2 10km')
plt.plot(np.arange(0,20)+1,O10_ohc_hintn,'b',label='ORAS5 10km')
plt.plot(np.arange(10,20)+1,O4_ohc_hintn,'--b',label='ORAS5 4km')
plt.text(15.,-6.05e9,'$\overline{Q}_{ECCO2-10km}$ = ' +  str(np.around(E10_shflx_hintn, decimals=0)) +' W m$^{-2}$')
plt.text(15.,-6.4e9,'$\overline{Q}_{ORAS5-10km}$ = ' +  str(np.around(O10_shflx_hintn, decimals=0)) +' W m$^{-2}$')
plt.text(15.,-6.8e9,'$\overline{Q}_{ORAS5-4km}$ = ' +  str(np.around(O4_shflx_hintn, decimals=0)) +' W m$^{-2}$')
#l1 = plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)
l1 = plt.legend(loc='lower left', borderaxespad=0)
plt.ylabel('heat content (J)')
ax2 = fig.add_subplot(312)
plt.plot(np.arange(0,20)+1,E10_salt_hintn,'r',label='ECCO2 10km')
plt.plot(np.arange(0,20)+1,O10_salt_hintn,'b',label='ORAS5 10km')
plt.plot(np.arange(10,20)+1,O4_salt_hintn,'--b',label='ORAS5 4km')
plt.text(.5,-59950,'$\overline{FWF}_{ORAS5-10km}$ = ' +  str(np.around(O10_sfwflx_hintn, decimals=0)) +' m s$^{-1}$')
plt.text(.5,-59990,'$\overline{FWF}_{ECCO2-10km}$ = ' +  str(np.around(E10_sfwflx_hintn, decimals=0)) +' m s$^{-1}$')
plt.text(.5,-60030,'$\overline{FWF}_{ORAS5-4km}$ = ' +  str(np.around(O4_sfwflx_hintn, decimals=0)) +' m s$^{-1}$')
#l2 = plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)
plt.ylabel('salt content')
ax3 = fig.add_subplot(313)
plt.plot(np.arange(0,20)+1,E10_MKE_havg,'r',label='ECCO2 10km')
plt.plot(np.arange(0,20)+1,O10_MKE_havg,'b',label='ORAS5 10km')
plt.plot(np.arange(10,20)+1,O4_MKE_havg,'--b',label='ORAS5 4km')
#l3 = plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)
plt.ylabel('MKE (m$^{2}$ s$^{-2}$)')
plt.xlabel('years')

name_fig="waom10x4_shflim_S_ORAS5emXECCO2_0.25Q_OHC+Salt+MKE_evolution.png"
plt.savefig(fig_path + name_fig, dpi=300)
plt.close()

fig = plt.figure(figsize=(12,8))
ax1 = fig.add_subplot(211)
plt.plot(np.arange(0,20)+1,E10_SST_havg,'r',label='ECCO2 10km')
plt.plot(np.arange(0,20)+1,O10_SST_havg,'b',label='ORAS5 10km')
plt.plot(np.arange(10,20)+1,O4_SST_havg,'--b',label='ORAS5 4km')
#l3 = plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)
l3 = plt.legend(loc="upper left", borderaxespad=0)
plt.ylabel('SST')
ax2 = fig.add_subplot(212)
plt.plot(np.arange(0,20)+1,E10_SSS_havg,'r',label='ECCO2 10km')
plt.plot(np.arange(0,20)+1,O10_SSS_havg,'b',label='ORAS5 10km')
plt.plot(np.arange(10,20)+1,O4_SSS_havg,'--b',label='ORAS5 4km')
#l3 = plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)
plt.ylabel('SSS')

name_fig="waom10x4_shflim_S_ORAS5emXECCO2_0.25Q_SST+SSS_evolution.png"
plt.savefig(fig_path + name_fig, dpi=300)
plt.close()

#fig = plt.figure(figsize=(14,20))
#ax1 = fig.add_subplot(211)
#plt.plot(np.arange(0,20)+1,ohc_hintn,'b',label='10km total')
#plt.plot(np.arange(0,20)+1,ohc_hintn_shelf,':b',label='10km shelf')
#plt.plot(np.arange(0,20)+1,ohc_hintn_open,'-*b',label='10km open-ocean')
#l1 = plt.legend(bbox_to_anchor=(0.7,0.72), loc="center left", borderaxespad=0)
#plt.ylabel('heat content (J)')
#ax3 = fig.add_subplot(212)
#plt.plot(np.arange(0,20)+1,salt_hintn,'r',label='10km total')
#plt.plot(np.arange(0,20)+1,salt_hintn_shelf,':r',label='10km shelf')
#plt.plot(np.arange(0,20)+1,salt_hintn_open,'-*r',label='10km open-ocean')
#l3 = plt.legend(bbox_to_anchor=(0.7,0.72), loc="center left", borderaxespad=0)
#
#name_fig="waom10_shflim_S_ORAS5em_0.25Q_OHC+Salt_subregions_evolution.png"
#plt.savefig(fig_path + name_fig, dpi=300)
#plt.close()
