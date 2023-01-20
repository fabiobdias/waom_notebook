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

import gsw
import pyresample
from xgcm import Grid

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

        ds = ds.set_coords(['Cs_r', 'Cs_w', 'hc', 'h', 'Vtransform'])

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

        u10_rho_sfc = np.zeros((12,560,630))
        v10_rho_sfc = np.zeros((12,560,630))
        u10_rho_bot = np.zeros((12,560,630))
        v10_rho_bot = np.zeros((12,560,630))

        for mm in np.arange(0,12):
            #interpoate u, v to rho grid:
            u10_interp = grid.interp(ds.u.isel(s_rho=0,ocean_time=mm), 'X',boundary='fill')
            v10_interp = grid.interp(ds.v.isel(s_rho=0,ocean_time=mm), 'Y',boundary='fill')
            u10_rho_bot[mm,:,:]=u10_interp
            v10_rho_bot[mm,:,:]=v10_interp
            del u10_interp,v10_interp
            u10_interp = grid.interp(ds.u.isel(s_rho=-1,ocean_time=mm), 'X',boundary='fill')
            v10_interp = grid.interp(ds.v.isel(s_rho=-1,ocean_time=mm), 'Y',boundary='fill')
            u10_rho_sfc[mm,:,:]=u10_interp
            v10_rho_sfc[mm,:,:]=v10_interp
            del u10_interp,v10_interp

        u10_rho_bot_ann = np.nanmean(u10_rho_bot, axis=0)
        v10_rho_bot_ann = np.nanmean(v10_rho_bot, axis=0)
        u10_rho_sfc_ann = np.nanmean(u10_rho_sfc, axis=0)
        v10_rho_sfc_ann = np.nanmean(v10_rho_sfc, axis=0)
        # concantenate annual averaged temp/salt
        if yr == '20':
            temp_ann = temp_tmp_ann
            salt_ann = salt_tmp_ann
            u10_bot = u10_rho_bot_ann
            v10_bot = v10_rho_bot_ann
            u10_sfc = u10_rho_sfc_ann
            v10_sfc = v10_rho_sfc_ann

        # calculated freezing temperature:
        Tf_bot = gsw.CT_freezing(salt_ann[0,:,:], np.nanmean(-z_rho[0,:,:,:],axis=2), 0)
        Tf_sfc = gsw.CT_freezing(salt_ann[-1,:,:], np.nanmean(-z_rho[-1,:,:,:],axis=2), 0)

        print('z_rho shape: ', z_rho.shape)
        print('Tf_bot shape: ', Tf_bot.shape)

    return temp_ann, salt_ann, u10_sfc, v10_sfc, u10_bot, v10_bot, Tf_sfc, Tf_bot

# read ROMS 4km:
def read_roms_ts_4km(exp_path,year):
    for yr  in [year]:
        ds = xr.open_dataset(exp_path + 'ocean_avg_00' + yr + '.nc')
        print('size temp and time length: ', ds.temp.shape, len(ds.salt.isel(xi_rho=20, eta_rho=100, s_rho=0)))
        temp_tmp = ds.variables["temp"]
        salt_tmp = ds.variables["salt"]
        zeta_tmp = ds.variables["zeta"]
        temp_tmp_ann = np.nanmean(temp_tmp, axis=0)
        salt_tmp_ann = np.nanmean(salt_tmp, axis=0)
        print('size temp_tmp_ann = ', temp_tmp_ann.shape)
        del temp_tmp, salt_tmp

        ds = ds.set_coords(['Cs_r', 'Cs_w', 'hc', 'h', 'Vtransform'])

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
        temp_ann = temp_tmp_ann
        salt_ann = salt_tmp_ann
        u4_bot = u4_rho_bot_ann
        v4_bot = v4_rho_bot_ann
        u4_sfc = u4_rho_sfc_ann
        v4_sfc = v4_rho_sfc_ann

        # calculated freezing temperature:
        Tf_bot = gsw.CT_freezing(salt_ann[0,:,:], np.nanmean(-z_rho[0,:,:,:],axis=2), 0)
        Tf_sfc = gsw.CT_freezing(salt_ann[-1,:,:], np.nanmean(-z_rho[-1,:,:,:],axis=2), 0)

    return temp_ann, salt_ann, u4_sfc, v4_sfc, u4_bot, v4_bot, Tf_sfc, Tf_bot

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

        ds = ds.set_coords(['Cs_r', 'Cs_w', 'hc', 'h', 'Vtransform'])

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

        u2_rho_sfc = np.zeros((12,2800,3150))
        v2_rho_sfc = np.zeros((12,2800,3150))
        u2_rho_bot = np.zeros((12,2800,3150))
        v2_rho_bot = np.zeros((12,2800,3150))

        for mm in np.arange(0,6):
            #interpoate u, v to rho grid:
            u2_interp = grid.interp(ds.u.isel(s_rho=0,ocean_time=mm), 'X',boundary='fill')
            v2_interp = grid.interp(ds.v.isel(s_rho=0,ocean_time=mm), 'Y',boundary='fill')
            u2_rho_bot[mm,:,:]=u2_interp
            v2_rho_bot[mm,:,:]=v2_interp
            del u2_interp,v2_interp
            u2_interp = grid.interp(ds.u.isel(s_rho=-1,ocean_time=mm), 'X',boundary='fill')
            v2_interp = grid.interp(ds.v.isel(s_rho=-1,ocean_time=mm), 'Y',boundary='fill')
            u2_rho_sfc[mm,:,:]=u2_interp
            v2_rho_sfc[mm,:,:]=v2_interp
            del u2_interp,v2_interp

        u2_rho_bot_ann = np.nanmean(u2_rho_bot, axis=0)
        v2_rho_bot_ann = np.nanmean(v2_rho_bot, axis=0)
        u2_rho_sfc_ann = np.nanmean(u2_rho_sfc, axis=0)
        v2_rho_sfc_ann = np.nanmean(v2_rho_sfc, axis=0)

        # concantenate annual averaged temp/salt
        if yr == '09':
            temp_ann = np.expand_dims(temp_tmp_ann, axis=0)
            salt_ann = np.expand_dims(salt_tmp_ann, axis=0)
            u2_sfc_ann = np.expand_dims(u2_rho_sfc_ann, axis=0)
            v2_sfc_ann = np.expand_dims(v2_rho_sfc_ann, axis=0)
            u2_bot_ann = np.expand_dims(u2_rho_bot_ann, axis=0)
            v2_bot_ann = np.expand_dims(v2_rho_bot_ann, axis=0)
        elif yr == '10':
            temp_ann = np.squeeze(np.nanmean(np.stack((temp_ann,np.expand_dims(temp_tmp_ann, axis=0)), axis=0), axis=0))
            salt_ann = np.squeeze(np.nanmean(np.stack((salt_ann,np.expand_dims(salt_tmp_ann, axis=0)), axis=0), axis=0))
            u2_sfc_ann = np.squeeze(np.nanmean(np.stack((u2_sfc_ann,np.expand_dims(u2_rho_sfc_ann, axis=0)), axis=0), axis=0))
            v2_sfc_ann = np.squeeze(np.nanmean(np.stack((v2_sfc_ann,np.expand_dims(v2_rho_sfc_ann, axis=0)), axis=0), axis=0))
            u2_bot_ann = np.squeeze(np.nanmean(np.stack((u2_bot_ann,np.expand_dims(u2_rho_bot_ann, axis=0)), axis=0), axis=0))
            v2_bot_ann = np.squeeze(np.nanmean(np.stack((v2_bot_ann,np.expand_dims(v2_rho_bot_ann, axis=0)), axis=0), axis=0))

        print('Annual temp and annual tmp temp sizes = ', temp_ann.shape, temp_tmp_ann.shape)

# calculated freezing temperature:
        Tf_bot = gsw.CT_freezing(salt_ann[0,:,:], np.nanmean(-z_rho[0,:,:,:],axis=2), 0)
        Tf_sfc = gsw.CT_freezing(salt_ann[-1,:,:], np.nanmean(-z_rho[-1,:,:,:],axis=2), 0)

    return temp_ann, salt_ann, u2_sfc_ann, v2_sfc_ann, u2_bot_ann, v2_bot_ann, Tf_sfc, Tf_bot

path_ECCO2_10km = '/scratch/project_2000339/boeiradi/waom10extend_shflim_S_0.25Q/output_01-20yr/'
path_ECCO2_4km = '/scratch/project_2000789/boeiradi/waom4extend_shflim_S_0.25Q/output_01-20yr/'
path_ECCO2_2km = '/scratch/project_2000339/boeiradi/waom2extend_shflim_S_0.25Q/output_01-05yr/'
# sensitivity experiments:
path_ECCO2_4km_notide = '/scratch/project_2000789/boeiradi/waom4extend_shflim_S_0.25Q/output_yr10_notides/'
path_ECCO2_4km_lowbathy = '/scratch/project_2000789/boeiradi/waom4extend_shflim_S_0.25Q/output_yr10_10km-bathy/'

temp_ann_10km, salt_ann_10km, u10_sfc, v10_sfc, u10_bot, v10_bot, Tf_sfc_10, Tf_bot_10 = read_roms_ts_10km(path_ECCO2_10km)
temp_ann_4km, salt_ann_4km, u4_sfc, v4_sfc, u4_bot, v4_bot, Tf_sfc_4, Tf_bot_4 = read_roms_ts_4km(path_ECCO2_4km,'10')
temp_ann_2km, salt_ann_2km, u2_sfc, v2_sfc, u2_bot, v2_bot, Tf_sfc_2, Tf_bot_2 = read_roms_ts_2km(path_ECCO2_2km)
temp_ann_4kmA, salt_ann_4kmA, u4_sfcA, v4_sfcA, u4_botA, v4_botA, Tf_sfc_4A, Tf_bot_4A = read_roms_ts_4km(path_ECCO2_4km_notide,'01')
temp_ann_4kmB, salt_ann_4kmB, u4_sfcB, v4_sfcB, u4_botB, v4_botB, Tf_sfc_4B, Tf_bot_4B = read_roms_ts_4km(path_ECCO2_4km_lowbathy,'01')

# define mask zice:
mask_zice_10km = ma.masked_where(zice_10km < 0, np.ones(zice_10km.shape))
mask_zice_4km = ma.masked_where(zice_4km < 0, np.ones(zice_4km.shape))
mask_zice_2km = ma.masked_where(zice_2km < 0, np.ones(zice_2km.shape))
mask_outice_10km = ma.masked_where(zice_10km >= 0, np.ones(zice_10km.shape))
mask_outice_4km = ma.masked_where(zice_4km >= 0, np.ones(zice_4km.shape))
mask_outice_2km = ma.masked_where(zice_2km >= 0, np.ones(zice_2km.shape))

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

# define sections to investigate Ross-Filchner and Ross IS dynamics:
xi_pt = [850, 1000, 1275, 1550]
##xi_pt = [875, 1000, 1275, 1550] # new sections, too far east
eta_sec_ini = [1600, 1650,  750,  500]
eta_sec_end = [1945, 1970, 1075, 1196]

# print array shapes before plotting:
print('lat_rho_2km shape', lat_rho_2km.shape)
print('lon_rho_2km shape', lon_rho_2km.shape)
print('temp_ann_2km (squeezed) shape', np.squeeze(temp_ann_2km[0,:,:]).shape)
print('temp_ann_2km shape', temp_ann_2km.shape)


## ==========================================
# interpolate WAOM4 and WAOM2 field to the WAOM10 grid for maps of difference:

w10_def = pyresample.geometry.SwathDefinition(lons=lon_rho_10km,lats=lat_rho_10km)
w4_def = pyresample.geometry.SwathDefinition(lons=lon_rho_4km,lats=lat_rho_4km)
w2_def = pyresample.geometry.SwathDefinition(lons=lon_rho_2km,lats=lat_rho_2km)
 
wf = lambda r: 1/r

sfc_temp_ann_4km_interp = pyresample.kd_tree.resample_custom(w4_def,temp_ann_4km[-1,:,:],w10_def,\
                                         radius_of_influence=30000,neighbours=4,weight_funcs=wf)
sfc_salt_ann_4km_interp = pyresample.kd_tree.resample_custom(w4_def,salt_ann_4km[-1,:,:],w10_def,\
                                         radius_of_influence=30000,neighbours=4,weight_funcs=wf)
bot_temp_ann_4km_interp = pyresample.kd_tree.resample_custom(w4_def,temp_ann_4km[0,:,:],w10_def,\
                                         radius_of_influence=30000,neighbours=4,weight_funcs=wf)
bot_salt_ann_4km_interp = pyresample.kd_tree.resample_custom(w4_def,salt_ann_4km[0,:,:],w10_def,\
                                         radius_of_influence=30000,neighbours=4,weight_funcs=wf)
sfc_temp_ann_2km_interp = pyresample.kd_tree.resample_custom(w2_def,temp_ann_2km[-1,:,:],w10_def,\
                                         radius_of_influence=5000,neighbours=1,weight_funcs=wf)
sfc_salt_ann_2km_interp = pyresample.kd_tree.resample_custom(w2_def,salt_ann_2km[-1,:,:],w10_def,\
                                         radius_of_influence=5000,neighbours=1,weight_funcs=wf)
bot_temp_ann_2km_interp = pyresample.kd_tree.resample_custom(w2_def,temp_ann_2km[0,:,:],w10_def,\
                                         radius_of_influence=5000,neighbours=1,weight_funcs=wf)
bot_salt_ann_2km_interp = pyresample.kd_tree.resample_custom(w2_def,salt_ann_2km[0,:,:],w10_def,\
                                         radius_of_influence=5000,neighbours=1,weight_funcs=wf)

# interpolate sensitive expts too: A=notide, B=coarse
bot_temp_ann_4kmA_interp = pyresample.kd_tree.resample_custom(w4_def,temp_ann_4kmA[0,:,:],w10_def,\
                                         radius_of_influence=30000,neighbours=4,weight_funcs=wf)
bot_salt_ann_4kmA_interp = pyresample.kd_tree.resample_custom(w4_def,salt_ann_4kmA[0,:,:],w10_def,\
                                         radius_of_influence=30000,neighbours=4,weight_funcs=wf)

bot_temp_ann_4kmB_interp = pyresample.kd_tree.resample_custom(w4_def,temp_ann_4kmB[0,:,:],w10_def,\
                                         radius_of_influence=30000,neighbours=4,weight_funcs=wf)
bot_salt_ann_4kmB_interp = pyresample.kd_tree.resample_custom(w4_def,salt_ann_4kmB[0,:,:],w10_def,\
                                         radius_of_influence=30000,neighbours=4,weight_funcs=wf)

sfc_temp_ann_4kmA_interp = pyresample.kd_tree.resample_custom(w4_def,temp_ann_4kmA[-1,:,:],w10_def,\
                                         radius_of_influence=30000,neighbours=4,weight_funcs=wf)
sfc_salt_ann_4kmA_interp = pyresample.kd_tree.resample_custom(w4_def,salt_ann_4kmA[-1,:,:],w10_def,\
                                         radius_of_influence=30000,neighbours=4,weight_funcs=wf)

sfc_temp_ann_4kmB_interp = pyresample.kd_tree.resample_custom(w4_def,temp_ann_4kmB[-1,:,:],w10_def,\
                                         radius_of_influence=30000,neighbours=4,weight_funcs=wf)
sfc_salt_ann_4kmB_interp = pyresample.kd_tree.resample_custom(w4_def,salt_ann_4kmB[-1,:,:],w10_def,\
                                         radius_of_influence=30000,neighbours=4,weight_funcs=wf)


# regridding u/v for plotting:
src = ccrs.SouthPolarStereo()
nx10, ny10 = 630, 560
xscale= [-.308e7, .345e7]
yscale= [-.308e7, .267e7]
xs10 = np.linspace(xscale[0], xscale[1], nx10)
ys10 = np.linspace(yscale[0], yscale[1], ny10)
xs10, ys10 = np.meshgrid(xs10, ys10)

nx4, ny4 = 1575, 1400
xs4 = np.linspace(xscale[0], xscale[1], nx4)
ys4 = np.linspace(yscale[0], yscale[1], ny4)
xs4, ys4 = np.meshgrid(xs4, ys4)

nx2, ny2 = 3150, 2800
xs2 = np.linspace(xscale[0], xscale[1], nx2)
ys2 = np.linspace(yscale[0], yscale[1], ny2)
xs2, ys2 = np.meshgrid(xs2, ys2)

# call cartopy projection
proj = ccrs.SouthPolarStereo()

# =============================================================================================================================
## Saving bottom temperature and currents in netcedf file for plotting within WMT maps:

fn = '/scratch/project_2000339/boeiradi/waom10extend_shflim_S_0.25Q/output_20yr_diag/WAOM10x4x2_Annual_bottom_temp_UV_maps.nc'
#fn = 'Full_vint_WMTmaps.nc'
dx = Dataset(fn, 'w', format='NETCDF4')

times = dx.createDimension('times', 1)
xi_rho_10 = dx.createDimension('xi_rho_10', 630)
eta_rho_10 = dx.createDimension('eta_rho_10', 560)
xi_rho_4 = dx.createDimension('xi_rho_4', 1575)
eta_rho_4 = dx.createDimension('eta_rho_4', 1400)
xi_rho_2 = dx.createDimension('xi_rho_2', 3150)
eta_rho_2 = dx.createDimension('eta_rho_2', 2800)

ocean_times = dx.createVariable('times', 'f4', ('times',))
ocean_times.units = 'seconds of the year'
xi = dx.createVariable('xi_rho', 'f4', ('xi_rho_10',))
eta = dx.createVariable('eta_rho', 'f4', ('eta_rho_10',))
xs_10 = dx.createVariable('xs_10', 'f4', ('eta_rho_10','xi_rho_10',))
ys_10 = dx.createVariable('ys_10', 'f4', ('eta_rho_10','xi_rho_10',))
Tf_10_sfc = dx.createVariable('Tf_10_sfc', 'f4', ('eta_rho_10','xi_rho_10',))
Tf_10_bot = dx.createVariable('Tf_10_bot', 'f4', ('eta_rho_10','xi_rho_10',))
xs_4 = dx.createVariable('xs_4', 'f4', ('eta_rho_4','xi_rho_4',))
ys_4 = dx.createVariable('ys_4', 'f4', ('eta_rho_4','xi_rho_4',))
Tf_4_sfc = dx.createVariable('Tf_4_sfc', 'f4', ('eta_rho_4','xi_rho_4',))
Tf_4_bot = dx.createVariable('Tf_4_bot', 'f4', ('eta_rho_4','xi_rho_4',))
xs_2 = dx.createVariable('xs_2', 'f4', ('eta_rho_2','xi_rho_2',))
ys_2 = dx.createVariable('ys_2', 'f4', ('eta_rho_2','xi_rho_2',))
Tf_2_sfc = dx.createVariable('Tf_2_sfc', 'f4', ('eta_rho_2','xi_rho_2',))
Tf_2_bot = dx.createVariable('Tf_2_bot', 'f4', ('eta_rho_2','xi_rho_2',))
Tf_4A_sfc = dx.createVariable('Tf_4A_sfc', 'f4', ('eta_rho_4','xi_rho_4',))
Tf_4A_bot = dx.createVariable('Tf_4A_bot', 'f4', ('eta_rho_4','xi_rho_4',))
Tf_4B_sfc = dx.createVariable('Tf_4B_sfc', 'f4', ('eta_rho_4','xi_rho_4',))
Tf_4B_bot = dx.createVariable('Tf_4B_bot', 'f4', ('eta_rho_4','xi_rho_4',))

# transformation maps variables:
# WAOM10:
temp_ann_bot_10 = dx.createVariable('temp_ann_bot_10', 'f4', ('times', 'eta_rho_10', 'xi_rho_10',))
temp_ann_bot_10.units = '10km, Annual map of bottom temperature (degC)'
salt_ann_bot_10 = dx.createVariable('salt_ann_bot_10', 'f4', ('times', 'eta_rho_10', 'xi_rho_10',))
salt_ann_bot_10.units = '10km, Annual map of bottom salinity (psu)'
u_10_bot = dx.createVariable('u_10_bot', 'f4', ('times', 'eta_rho_10', 'xi_rho_10',))
u_10_bot.units = '10km Annual map of bottom U-currents (m/s)'
v_10_bot = dx.createVariable('v_10_bot', 'f4', ('times', 'eta_rho_10', 'xi_rho_10',))
v_10_bot.units = '10km Annual map of bottom V-currents (m/s)'

temp_ann_sfc_10 = dx.createVariable('temp_ann_sfc_10', 'f4', ('times', 'eta_rho_10', 'xi_rho_10',))
temp_ann_sfc_10.units = '10km, Annual map of surface temperature (degC)'
salt_ann_sfc_10 = dx.createVariable('salt_ann_sfc_10', 'f4', ('times', 'eta_rho_10', 'xi_rho_10',))
salt_ann_sfc_10.units = '10km, Annual map of surfcae salinity (psu)'
u_10_sfc = dx.createVariable('u_10_sfc', 'f4', ('times', 'eta_rho_10', 'xi_rho_10',))
u_10_sfc.units = '10km Annual map of surface U-currents (m/s)'
v_10_sfc = dx.createVariable('v_10_sfc', 'f4', ('times', 'eta_rho_10', 'xi_rho_10',))
v_10_sfc.units = '10km Annual map of surfcae V-currents (m/s)'

# WAOM4:
temp_ann_bot_4 = dx.createVariable('temp_ann_bot_4', 'f4', ('times', 'eta_rho_10', 'xi_rho_10',))
temp_ann_bot_4.units = '4km, Annual map of bottom temperature (degC)'
salt_ann_bot_4 = dx.createVariable('salt_ann_bot_4', 'f4', ('times', 'eta_rho_10', 'xi_rho_10',))
salt_ann_bot_4.units = '4km, Annual map of bottom salinity (psu)'
u_4_bot = dx.createVariable('u_4_bot', 'f4', ('times', 'eta_rho_4', 'xi_rho_4',))
u_4_bot.units = '4km Annual map of bottom U-currents (m/s)'
v_4_bot = dx.createVariable('v_4_bot', 'f4', ('times', 'eta_rho_4', 'xi_rho_4',))
v_4_bot.units = '4km Annual map of bottom V-currents (m/s)'

temp_ann_sfc_4 = dx.createVariable('temp_ann_sfc_4', 'f4', ('times', 'eta_rho_10', 'xi_rho_10',))
temp_ann_sfc_4.units = '4km, Annual map of surfcace temperature (degC)'
salt_ann_sfc_4 = dx.createVariable('salt_ann_sfc_4', 'f4', ('times', 'eta_rho_10', 'xi_rho_10',))
salt_ann_sfc_4.units = '4km, Annual map of surface salinity (psu)'
u_4_sfc = dx.createVariable('u_4_sfc', 'f4', ('times', 'eta_rho_4', 'xi_rho_4',))
u_4_sfc.units = '4km Annual map of surface U-currents (m/s)'
v_4_sfc = dx.createVariable('v_4_sfc', 'f4', ('times', 'eta_rho_4', 'xi_rho_4',))
v_4_sfc.units = '4km Annual map of surface V-currents (m/s)'

# WAOM2:
temp_ann_bot_2 = dx.createVariable('temp_ann_bot_2', 'f4', ('times', 'eta_rho_10', 'xi_rho_10',))
temp_ann_bot_2.units = '2km, Annual map of bottom temperature (degC)'
salt_ann_bot_2 = dx.createVariable('salt_ann_bot_2', 'f4', ('times', 'eta_rho_10', 'xi_rho_10',))
salt_ann_bot_2.units = '2km, Annual map of bottom salinity (psu)'
u_2_bot = dx.createVariable('u_2_bot', 'f4', ('times', 'eta_rho_2', 'xi_rho_2',))
u_2_bot.units = '2km Annual map of bottom U-currents (m/s)'
v_2_bot = dx.createVariable('v_2_bot', 'f4', ('times', 'eta_rho_2', 'xi_rho_2',))
v_2_bot.units = '2km Annual map of bottom V-currents (m/s)'

temp_ann_sfc_2 = dx.createVariable('temp_anns_fc_2', 'f4', ('times', 'eta_rho_10', 'xi_rho_10',))
temp_ann_sfc_2.units = '2km, Annual map of surface temperature (degC)'
salt_ann_sfc_2 = dx.createVariable('salt_ann_sfc_2', 'f4', ('times', 'eta_rho_10', 'xi_rho_10',))
salt_ann_sfc_2.units = '2km, Annual map of surface salinity (psu)'
u_2_sfc = dx.createVariable('u_2_sfc', 'f4', ('times', 'eta_rho_2', 'xi_rho_2',))
u_2_sfc.units = '2km Annual map of surface U-currents (m/s)'
v_2_sfc = dx.createVariable('v_2_sfc', 'f4', ('times', 'eta_rho_2', 'xi_rho_2',))
v_2_sfc.units = '2km Annual map of surface V-currents (m/s)'

# WAOM4-NOTIDE:
temp_ann_bot_4nt = dx.createVariable('temp_ann_bot_4nt', 'f4', ('times', 'eta_rho_10', 'xi_rho_10',))
temp_ann_bot_4nt.units = '4km notide, Annual map of bottom temperature (degC)'
salt_ann_bot_4nt = dx.createVariable('salt_ann_bot_4nt', 'f4', ('times', 'eta_rho_10', 'xi_rho_10',))
salt_ann_bot_4nt.units = '4km notide, Annual map of bottom salinity (psu)'
u_4nt_bot = dx.createVariable('u_4nt_bot', 'f4', ('times', 'eta_rho_4', 'xi_rho_4',))
u_4nt_bot.units = '4km notide Annual map of bottom U-currents (m/s)'
v_4nt_bot = dx.createVariable('v_4nt_bot', 'f4', ('times', 'eta_rho_4', 'xi_rho_4',))
v_4nt_bot.units = '4km notide Annual map of bottom V-currents (m/s)'

temp_ann_sfc_4nt = dx.createVariable('temp_ann_sfc_4nt', 'f4', ('times', 'eta_rho_10', 'xi_rho_10',))
temp_ann_sfc_4nt.units = '4km notide, Annual map of surface temperature (degC)'
salt_ann_sfc_4nt = dx.createVariable('salt_ann_sfc_4nt', 'f4', ('times', 'eta_rho_10', 'xi_rho_10',))
salt_ann_sfc_4nt.units = '4km notide, Annual map of surface salinity (psu)'
u_4nt_sfc = dx.createVariable('u_4nt_sfc', 'f4', ('times', 'eta_rho_4', 'xi_rho_4',))
u_4nt_sfc.units = '4km notide Annual map of surface U-currents (m/s)'
v_4nt_sfc = dx.createVariable('v_4nt_sfc', 'f4', ('times', 'eta_rho_4', 'xi_rho_4',))
v_4nt_sfc.units = '4km notide Annual map of surface V-currents (m/s)'

# WAOM4-COARSE:
temp_ann_bot_4c = dx.createVariable('temp_ann_bot_4c', 'f4', ('times', 'eta_rho_10', 'xi_rho_10',))
temp_ann_bot_4c.units = '4km coarse, Annual map of bottom temperature (degC)'
salt_ann_bot_4c = dx.createVariable('salt_ann_bot_4c', 'f4', ('times', 'eta_rho_10', 'xi_rho_10',))
salt_ann_bot_4c.units = '4km coarse, Annual map of bottom salinity (psu)'
u_4c_bot = dx.createVariable('u_4c_bot', 'f4', ('times', 'eta_rho_4', 'xi_rho_4',))
u_4c_bot.units = '4km coarse Annual map of bottom U-currents (m/s)'
v_4c_bot = dx.createVariable('v_4c_bot', 'f4', ('times', 'eta_rho_4', 'xi_rho_4',))
v_4c_bot.units = '4km coarse Annual map of bottom V-currents (m/s)'

temp_ann_sfc_4c = dx.createVariable('temp_ann_sfc_4c', 'f4', ('times', 'eta_rho_10', 'xi_rho_10',))
temp_ann_sfc_4c.units = '4km coarse, Annual map of surface temperature (degC)'
salt_ann_sfc_4c = dx.createVariable('salt_ann_sfc_4c', 'f4', ('times', 'eta_rho_10', 'xi_rho_10',))
salt_ann_sfc_4c.units = '4km coarse, Annual map of surface salinity (psu)'
u_4c_sfc = dx.createVariable('u_4c_sfc', 'f4', ('times', 'eta_rho_4', 'xi_rho_4',))
u_4c_sfc.units = '4km coarse Annual map of surface U-currents (m/s)'
v_4c_sfc = dx.createVariable('v_4c_sfc', 'f4', ('times', 'eta_rho_4', 'xi_rho_4',))
v_4c_sfc.units = '4km coarse Annual map of surface V-currents (m/s)'


ocean_times[:] = 0 #annual   #np.arange(1314000,31536000,2628000) # monthly
xi[:] = np.arange(0,630)
eta[:] = np.arange(0,560)
xs_10[:] = xs10
ys_10[:] = ys10
Tf_10_sfc[:] = Tf_sfc_10
Tf_10_bot[:] = Tf_bot_10
xs_4[:] = xs4
ys_4[:] = ys4
Tf_4_sfc[:] = Tf_sfc_4
Tf_4_bot[:] = Tf_bot_4
xs_2[:] = xs2
ys_2[:] = ys2
Tf_2_sfc[:] = Tf_sfc_2
Tf_2_bot[:] = Tf_bot_2
Tf_4A_sfc[:] = Tf_sfc_4A
Tf_4A_bot[:] = Tf_bot_4A
Tf_4B_sfc[:] = Tf_sfc_4B
Tf_4B_bot[:] = Tf_bot_4B

# WAOM10
temp_ann_bot_10[:] = np.squeeze(temp_ann_10km[0,:,:])
salt_ann_bot_10[:] = np.squeeze(salt_ann_10km[0,:,:])
u_10_bot[:] = u10_bot
v_10_bot[:] = v10_bot

temp_ann_sfc_10[:] = np.squeeze(temp_ann_10km[-1,:,:])
salt_ann_sfc_10[:] = np.squeeze(salt_ann_10km[-1,:,:])
u_10_sfc[:] = u10_sfc
v_10_sfc[:] = v10_sfc

# WAOM4
temp_ann_bot_4[:] = bot_temp_ann_4km_interp
salt_ann_bot_4[:] = bot_salt_ann_4km_interp
u_4_bot[:] = u4_bot
v_4_bot[:] = v4_bot

temp_ann_sfc_4[:] = sfc_temp_ann_4km_interp
salt_ann_sfc_4[:] = sfc_salt_ann_4km_interp
u_4_sfc[:] = u4_sfc
v_4_sfc[:] = v4_sfc

# WAOM2
temp_ann_bot_2[:] = bot_temp_ann_2km_interp
salt_ann_bot_2[:] = bot_salt_ann_2km_interp
u_2_bot[:] = u2_bot
v_2_bot[:] = v2_bot

temp_ann_sfc_2[:] = sfc_temp_ann_2km_interp
salt_ann_sfc_2[:] = sfc_salt_ann_2km_interp
u_2_sfc[:] = u2_sfc
v_2_sfc[:] = v2_sfc

# WAOM4-NOTIDE
temp_ann_bot_4nt[:] = bot_temp_ann_4kmA_interp
salt_ann_bot_4nt[:] = bot_salt_ann_4kmA_interp
u_4nt_bot[:] = u4_botA
v_4nt_bot[:] = v4_botA

temp_ann_sfc_4nt[:] = sfc_temp_ann_4kmA_interp
salt_ann_sfc_4nt[:] = sfc_salt_ann_4kmA_interp
u_4nt_sfc[:] = u4_sfcA
v_4nt_sfc[:] = v4_sfcA

# WAOM4-COARSE
temp_ann_bot_4c[:] = bot_temp_ann_4kmB_interp
salt_ann_bot_4c[:] = bot_salt_ann_4kmB_interp
u_4c_bot[:] = u4_botB
v_4c_bot[:] = v4_botB

temp_ann_sfc_4c[:] = sfc_temp_ann_4kmB_interp
salt_ann_sfc_4c[:] = sfc_salt_ann_4kmB_interp
u_4c_sfc[:] = u4_sfcB
v_4c_sfc[:] = v4_sfcB

dx.close()

# =============================================================================================================================


# ========== subregions plots:

# for 2000m isoline
bathym = cfeature.NaturalEarthFeature(name='bathymetry_J_1000', scale='10m', category='physical')
# limits for contour of ice front (Ronne-Filchner IS):
xlimit = np.arange(300,500,1)
ylimit = np.arange(100,300,1)

tmin = -2.7
tmax = 0.5
smin = 34.1
smax = 34.5

fig = plt.figure(figsize=(8,10))

ax1 = fig.add_subplot(321, projection=proj)
ct1=plt.pcolormesh(lon_rho_10km,lat_rho_10km,np.squeeze(temp_ann_10km[-1,:,:]), transform=ccrs.PlateCarree(), cmap=plt.cm.plasma, vmin=tmin, vmax=tmax)
cv=ax1.quiver(xs10[::3,::3],ys10[::3,::3],u10_sfc[::3,::3],v10_sfc[::3,::3], color='k', transform=src, scale=1.)
plt.title('A) Surface temperature \n WAOM10')
ax1.gridlines() # draw_labels=True,linewidth=
ax1.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
ratio = .9
x_left, x_right = ax1.get_xlim()
y_low, y_high = ax1.get_ylim()
#ax1.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax1.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
# Add limits subregions:
plt.plot(-25*np.ones(len(np.arange(-76,-64,.1))),np.arange(-76,-64,.1), linestyle='dashed', color='w', transform=ccrs.PlateCarree(), linewidth=2)
plt.plot(-65*np.ones(len(np.arange(-67,-60.2,.1))),np.arange(-67,-60.2,.1), linestyle='dashed', color='w', transform=ccrs.PlateCarree(), linewidth=2)
plt.text(-61,-65,'Ronne- \n Filchner IS',weight='bold', color='w', transform=ccrs.PlateCarree(), fontsize=9)
plt.plot(60*np.ones(len(np.arange(-67.5,-55.7,.1))),np.arange(-67.5,-55.7,.1), linestyle='dashed', color='w', transform=ccrs.PlateCarree(), linewidth=2)
plt.text(15,-69,'Maud Land',weight='bold', color='w', transform=ccrs.PlateCarree(), fontsize=9)
plt.plot(160*np.ones(len(np.arange(-70,-61,.1))),np.arange(-70,-61,.1), linestyle='dashed', color='w', transform=ccrs.PlateCarree(), linewidth=2)
plt.text(155,-62,'East \n Antarctic',weight='bold', color='w', transform=ccrs.PlateCarree(), fontsize=9)
plt.plot(-140*np.ones(len(np.arange(-75.5,-55,.1))),np.arange(-75.5,-55,.1), linestyle='dashed', color='w', transform=ccrs.PlateCarree(), linewidth=2)
plt.text(-150,-63,'Ross IS',weight='bold', color='w', transform=ccrs.PlateCarree(), fontsize=9)
plt.text(-115,-60,'West \n Antarctic',weight='bold', color='w', transform=ccrs.PlateCarree(), fontsize=9)
ax1.add_feature(bathym, zorder=3, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='gray', alpha=.65)

cbar_ax3 = fig.add_axes([0.04, 0.65, 0.01, 0.25])
fig.colorbar(ct1, cax=cbar_ax3, orientation='vertical')
cbar_ax3.set_xlabel('$^{\circ}$C')#, labelpad=-35)

ax2 = fig.add_subplot(322, projection=proj)
cs1=plt.pcolormesh(lon_rho_10km,lat_rho_10km,np.squeeze(salt_ann_10km[-1,:,:]), transform=ccrs.PlateCarree(), cmap=plt.cm.plasma, vmin=smin, vmax=smax)
cv=ax2.quiver(xs10[::3,::3],ys10[::3,::3],u10_sfc[::3,::3],v10_sfc[::3,::3], color='k', transform=src, scale=1.)
plt.title('B) Surface salinity \n WAOM10')
ax2.gridlines()
x_left, x_right = ax2.get_xlim()
y_low, y_high = ax2.get_ylim()
#ax2.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax2.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
ax2.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
ax2.add_feature(bathym, zorder=3, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='gray', alpha=.65)

cbar_ax4 = fig.add_axes([0.92, 0.65, 0.01, 0.25])
fig.colorbar(cs1, cax=cbar_ax4, orientation='vertical')
cbar_ax4.set_xlabel('')#, labelpad=-35)

# 4/2km anomalies
Atmin=-.15
Atmax=.15
Asmin=-0.15
Asmax=0.15

ax3 = fig.add_subplot(323, projection=proj)
ct2=plt.pcolormesh(lon_rho_10km,lat_rho_10km,(sfc_temp_ann_4km_interp-np.squeeze(temp_ann_10km[-1,:,:])), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=Atmin, vmax=Atmax)
cv=ax3.quiver(xs4[::6,::6],ys4[::6,::6],u4_sfc[::6,::6],v4_sfc[::6,::6], color='k', transform=src, scale=1.)
plt.title('C) WAOM4 - WAOM10')
ax3.gridlines() # draw_labels=True,linewidth=
x_left, x_right = ax3.get_xlim()
y_low, y_high = ax3.get_ylim()
#ax3.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax3.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
ax3.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
ax3.add_feature(bathym, zorder=3, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='gray', alpha=.65)

ax4 = fig.add_subplot(324, projection=proj)
cs2=plt.pcolormesh(lon_rho_10km,lat_rho_10km,(sfc_salt_ann_4km_interp-np.squeeze(salt_ann_10km[-1,:,:])), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=Asmin, vmax=Asmax)
cv=ax4.quiver(xs4[::6,::6],ys4[::6,::6],u4_sfc[::6,::6],v4_sfc[::6,::6], color='k', transform=src, scale=1.)
plt.title('D) WAOM4 - WAOM10')
ax4.gridlines()
x_left, x_right = ax4.get_xlim()
y_low, y_high = ax4.get_ylim()
#ax4.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax4.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
ax4.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
ax4.add_feature(bathym, zorder=3, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='gray', alpha=.65)

ax5 = fig.add_subplot(325, projection=proj)
ct2=plt.pcolormesh(lon_rho_10km,lat_rho_10km,(sfc_temp_ann_2km_interp-np.squeeze(temp_ann_10km[-1,:,:])), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=Atmin, vmax=Atmax)
cv=ax5.quiver(xs2[::12,::12],ys2[::12,::12],u2_sfc[::12,::12],v2_sfc[::12,::12], color='k', transform=src, scale=1.)
plt.title('E) WAOM2 - WAOM10')
ax5.gridlines() # draw_labels=True,linewidth=
x_left, x_right = ax5.get_xlim()
y_low, y_high = ax5.get_ylim()
#ax5.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax5.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
ax5.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
ax5.add_feature(bathym, zorder=3, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='gray', alpha=.65)
ii=0
plt.plot([lon_rho_2km[eta_sec_ini[ii],xi_pt[ii]], lon_rho_2km[eta_sec_end[ii],xi_pt[ii]]],[lat_rho_2km[eta_sec_ini[ii],xi_pt[ii]], lat_rho_2km[eta_sec_end[ii],xi_pt[ii]]],color='yellow',transform=ccrs.PlateCarree(),linewidth=1.5, linestyle='dashed')
ii=1
plt.plot([lon_rho_2km[eta_sec_ini[ii],xi_pt[ii]], lon_rho_2km[eta_sec_end[ii],xi_pt[ii]]],[lat_rho_2km[eta_sec_ini[ii],xi_pt[ii]], lat_rho_2km[eta_sec_end[ii],xi_pt[ii]]],color='yellow',transform=ccrs.PlateCarree(),linewidth=1.5, linestyle='dashed')

ax6 = fig.add_subplot(326, projection=proj)
cs2=plt.pcolormesh(lon_rho_10km,lat_rho_10km,(sfc_salt_ann_2km_interp-np.squeeze(salt_ann_10km[-1,:,:])), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=Asmin, vmax=Asmax)
cv=ax6.quiver(xs2[::12,::12],ys2[::12,::12],u2_sfc[::12,::12],v2_sfc[::12,::12], color='k', transform=src, scale=1.)
plt.title('F) WAOM2 - WAOM10')
ax6.gridlines()
x_left, x_right = ax6.get_xlim()
y_low, y_high = ax6.get_ylim()
#ax6.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax6.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
ax6.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
ax6.add_feature(bathym, zorder=3, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='gray', alpha=.65)
ii=0
plt.plot([lon_rho_2km[eta_sec_ini[ii],xi_pt[ii]], lon_rho_2km[eta_sec_end[ii],xi_pt[ii]]],[lat_rho_2km[eta_sec_ini[ii],xi_pt[ii]], lat_rho_2km[eta_sec_end[ii],xi_pt[ii]]],color='yellow',transform=ccrs.PlateCarree(),linewidth=1.5, linestyle='dashed')
ii=1
plt.plot([lon_rho_2km[eta_sec_ini[ii],xi_pt[ii]], lon_rho_2km[eta_sec_end[ii],xi_pt[ii]]],[lat_rho_2km[eta_sec_ini[ii],xi_pt[ii]], lat_rho_2km[eta_sec_end[ii],xi_pt[ii]]],color='yellow',transform=ccrs.PlateCarree(),linewidth=1.5, linestyle='dashed')

cbar_ax1 = fig.add_axes([0.12, 0.08, 0.35, 0.01])
fig.colorbar(ct2, cax=cbar_ax1, orientation='horizontal')
cbar_ax1.set_xlabel('$^{\circ}$C')#, labelpad=-35)

cbar_ax2 = fig.add_axes([0.55, 0.08, 0.35, 0.01])
fig.colorbar(cs2, cax=cbar_ax2, orientation='horizontal')
cbar_ax2.set_xlabel('')#, labelpad=-35)

name_fig="waom10x4x2extend_shflim_S_0.25Q_sfc_TSVel_RFIS_anomalies.png"
plt.savefig(fig_path + name_fig, dpi=300)
plt.close()


# --- bottom

tmin = -2.7
tmax = 0.5
smin = 34.1
smax = 34.5

fig = plt.figure(figsize=(8,10))

ax1 = fig.add_subplot(321, projection=proj)
ct1=plt.pcolormesh(lon_rho_10km,lat_rho_10km,np.squeeze(temp_ann_10km[0,:,:]), transform=ccrs.PlateCarree(), cmap=plt.cm.plasma, vmin=tmin, vmax=tmax)
cv=ax1.quiver(xs10[::3,::3],ys10[::3,::3],u10_bot[::3,::3],v10_bot[::3,::3], color='k', transform=src, scale=1.)
plt.title('Bottom temperature \n WAOM10')
ax1.gridlines() # draw_labels=True,linewidth=
ax1.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
ratio = .9
x_left, x_right = ax1.get_xlim()
y_low, y_high = ax1.get_ylim()
#ax1.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax1.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
ax1.add_feature(bathym, zorder=3, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='gray', alpha=.65)

cbar_ax3 = fig.add_axes([0.04, 0.65, 0.01, 0.25])
fig.colorbar(ct1, cax=cbar_ax3, orientation='vertical')
cbar_ax3.set_xlabel('$^{\circ}$C')#, labelpad=-35)

ax2 = fig.add_subplot(322, projection=proj)
cs1=plt.pcolormesh(lon_rho_10km,lat_rho_10km,np.squeeze(salt_ann_10km[0,:,:]), transform=ccrs.PlateCarree(), cmap=plt.cm.plasma, vmin=smin, vmax=smax)
cv=ax2.quiver(xs10[::3,::3],ys10[::3,::3],u10_bot[::3,::3],v10_bot[::3,::3], color='k', transform=src, scale=1.)
plt.title('Bottom salinity \n WAOM10')
ax2.gridlines()
x_left, x_right = ax2.get_xlim()
y_low, y_high = ax2.get_ylim()
#ax2.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax2.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
ax2.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
ax2.add_feature(bathym, zorder=3, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='gray', alpha=.65)

cbar_ax4 = fig.add_axes([0.92, 0.65, 0.01, 0.25])
fig.colorbar(cs1, cax=cbar_ax4, orientation='vertical')
cbar_ax4.set_xlabel('')#, labelpad=-35)

# 4/2km anomalies
Atmin=-.6
Atmax=.6
Asmin=-0.1
Asmax=0.1

ax3 = fig.add_subplot(323, projection=proj)
ct2=plt.pcolormesh(lon_rho_10km,lat_rho_10km,(bot_temp_ann_4km_interp-np.squeeze(temp_ann_10km[0,:,:])), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=Atmin, vmax=Atmax)
cv=ax3.quiver(xs4[::6,::6],ys4[::6,::6],u4_bot[::6,::6],v4_bot[::6,::6], color='k', transform=src, scale=1.)
plt.title('WAOM4 - WAOM10')
ax3.gridlines() # draw_labels=True,linewidth=
x_left, x_right = ax3.get_xlim()
y_low, y_high = ax3.get_ylim()
#ax3.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax3.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
ax3.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
ax3.add_feature(bathym, zorder=3, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='gray', alpha=.65)

ax4 = fig.add_subplot(324, projection=proj)
cs2=plt.pcolormesh(lon_rho_10km,lat_rho_10km,(bot_salt_ann_4km_interp-np.squeeze(salt_ann_10km[0,:,:])), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=Asmin, vmax=Asmax)
cv=ax4.quiver(xs4[::6,::6],ys4[::6,::6],u4_bot[::6,::6],v4_bot[::6,::6], color='k', transform=src, scale=1.)
plt.title('WAOM4 - WAOM10')
ax4.gridlines()
x_left, x_right = ax4.get_xlim()
y_low, y_high = ax4.get_ylim()
#ax4.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax4.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
ax4.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
ax4.add_feature(bathym, zorder=3, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='gray', alpha=.65)

ax5 = fig.add_subplot(325, projection=proj)
ct2=plt.pcolormesh(lon_rho_10km,lat_rho_10km,(bot_temp_ann_2km_interp-np.squeeze(temp_ann_10km[0,:,:])), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=Atmin, vmax=Atmax)
cv=ax5.quiver(xs2[::12,::12],ys2[::12,::12],u2_bot[::12,::12],v2_bot[::12,::12], color='k', transform=src, scale=1.)
plt.title('WAOM2 - WAOM10')
ax5.gridlines() # draw_labels=True,linewidth=
x_left, x_right = ax5.get_xlim()
y_low, y_high = ax5.get_ylim()
#ax5.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax5.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
ax5.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
ax5.add_feature(bathym, zorder=3, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='gray', alpha=.65)
ii=0
plt.plot([lon_rho_2km[eta_sec_ini[ii],xi_pt[ii]], lon_rho_2km[eta_sec_end[ii],xi_pt[ii]]],[lat_rho_2km[eta_sec_ini[ii],xi_pt[ii]], lat_rho_2km[eta_sec_end[ii],xi_pt[ii]]],color='yellow',transform=ccrs.PlateCarree(),linewidth=1.5, linestyle='dashed')
ii=1
plt.plot([lon_rho_2km[eta_sec_ini[ii],xi_pt[ii]], lon_rho_2km[eta_sec_end[ii],xi_pt[ii]]],[lat_rho_2km[eta_sec_ini[ii],xi_pt[ii]], lat_rho_2km[eta_sec_end[ii],xi_pt[ii]]],color='yellow',transform=ccrs.PlateCarree(),linewidth=1.5, linestyle='dashed')

ax6 = fig.add_subplot(326, projection=proj)
cs2=plt.pcolormesh(lon_rho_10km,lat_rho_10km,(bot_salt_ann_2km_interp-np.squeeze(salt_ann_10km[0,:,:])), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=Asmin, vmax=Asmax)
cv=ax6.quiver(xs2[::12,::12],ys2[::12,::12],u2_bot[::12,::12],v2_bot[::12,::12], color='k', transform=src, scale=1.)
plt.title('WAOM2 - WAOM10')
ax6.gridlines()
x_left, x_right = ax6.get_xlim()
y_low, y_high = ax6.get_ylim()
#ax6.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax6.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
ax6.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
ax6.add_feature(bathym, zorder=3, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='gray', alpha=.65)
ii=0
plt.plot([lon_rho_2km[eta_sec_ini[ii],xi_pt[ii]], lon_rho_2km[eta_sec_end[ii],xi_pt[ii]]],[lat_rho_2km[eta_sec_ini[ii],xi_pt[ii]], lat_rho_2km[eta_sec_end[ii],xi_pt[ii]]],color='yellow',transform=ccrs.PlateCarree(),linewidth=1.5, linestyle='dashed')
ii=1
plt.plot([lon_rho_2km[eta_sec_ini[ii],xi_pt[ii]], lon_rho_2km[eta_sec_end[ii],xi_pt[ii]]],[lat_rho_2km[eta_sec_ini[ii],xi_pt[ii]], lat_rho_2km[eta_sec_end[ii],xi_pt[ii]]],color='yellow',transform=ccrs.PlateCarree(),linewidth=1.5, linestyle='dashed')

cbar_ax1 = fig.add_axes([0.12, 0.08, 0.35, 0.01])
fig.colorbar(ct2, cax=cbar_ax1, orientation='horizontal')
cbar_ax1.set_xlabel('$^{\circ}$C')#, labelpad=-35)

cbar_ax2 = fig.add_axes([0.55, 0.08, 0.35, 0.01])
fig.colorbar(cs2, cax=cbar_ax2, orientation='horizontal')
cbar_ax2.set_xlabel('')#, labelpad=-35)

name_fig="waom10x4x2extend_shflim_S_0.25Q_bottom_TSVel_RFIS_anomalies.png"
plt.savefig(fig_path + name_fig, dpi=300)
plt.close()

# ---- plot for paper: left column = temp bottom; right column = surface salimity

tmin = -2.7
tmax = 0.5
smin = 34.1
smax = 34.5

fig = plt.figure(figsize=(8,10))

ax1 = fig.add_subplot(321, projection=proj)
ct1=plt.pcolormesh(lon_rho_10km,lat_rho_10km,np.squeeze(temp_ann_10km[0,:,:]), transform=ccrs.PlateCarree(), cmap=plt.cm.plasma, vmin=tmin, vmax=tmax)
cv=ax1.quiver(xs10[::3,::3],ys10[::3,::3],u10_bot[::3,::3]*mask_zice_10km[::3,::3],v10_bot[::3,::3]*mask_zice_10km[::3,::3], color='k', transform=src, scale=2.)
cv=ax1.quiver(xs10[::3,::3],ys10[::3,::3],u10_bot[::3,::3]*mask_outice_10km[::3,::3],v10_bot[::3,::3]*mask_outice_10km[::3,::3], color='k', transform=src, scale=1.)
ax1.quiver(-.7e6,1.2e6,.2,0,color='w', transform=src, scale=2., zorder=4)
ax1.text(-.7e6,1.25e6,'off-cavity \n 0.2 m/s',color='w', fontsize=8, transform=src, zorder=4)
ax1.quiver(-.8e6,.1e6,.1,0,color='w', transform=src, scale=1., zorder=4)
ax1.text(-.8e6,.15e6,'sub-cavity \n 0.1 m/s',color='w', fontsize=8, transform=src, zorder=4)
plt.title('A) Bottom temperature \n WAOM10')
ax1.gridlines() # draw_labels=True,linewidth=
ax1.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
ratio = .9
x_left, x_right = ax1.get_xlim()
y_low, y_high = ax1.get_ylim()
#ax1.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='gray', alpha=.65)
ax1.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
ax1.add_feature(bathym, zorder=3, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)

cbar_ax3 = fig.add_axes([0.05, 0.65, 0.01, 0.25])
fig.colorbar(ct1, cax=cbar_ax3, orientation='vertical')
cbar_ax3.set_ylabel('$^{\circ}$C')#, labelpad=-35)

ax2 = fig.add_subplot(322, projection=proj)
cs1=plt.pcolormesh(lon_rho_10km,lat_rho_10km,np.squeeze(salt_ann_10km[-1,:,:]), transform=ccrs.PlateCarree(), cmap=plt.cm.plasma, vmin=smin, vmax=smax)
cv=ax2.quiver(xs10[::3,::3],ys10[::3,::3],u10_sfc[::3,::3]*mask_zice_10km[::3,::3],v10_sfc[::3,::3]*mask_zice_10km[::3,::3], color='k', transform=src, scale=2.)
cv=ax2.quiver(xs10[::3,::3],ys10[::3,::3],u10_sfc[::3,::3]*mask_outice_10km[::3,::3],v10_sfc[::3,::3]*mask_outice_10km[::3,::3], color='k', transform=src, scale=1.)
ax2.quiver(-.7e6,1.2e6,.2,0,color='w', transform=src, scale=2., zorder=4)
ax2.text(-.7e6,1.25e6,'off-cavity \n 0.2 m/s',color='w', fontsize=8, transform=src, zorder=4)
ax2.quiver(-.8e6,.1e6,.1,0,color='w', transform=src, scale=1., zorder=4)
ax2.text(-.8e6,.15e6,'sub-cavity \n 0.1 m/s',color='w', fontsize=8, transform=src, zorder=4)
plt.title('B) Surface salinity \n WAOM10')
ax2.gridlines()
x_left, x_right = ax2.get_xlim()
y_low, y_high = ax2.get_ylim()
#ax2.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax2.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='gray', alpha=.65)
ax2.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
ax2.add_feature(bathym, zorder=3, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)


cbar_ax4 = fig.add_axes([0.90, 0.65, 0.01, 0.25])
fig.colorbar(cs1, cax=cbar_ax4, orientation='vertical')
cbar_ax4.set_xlabel('')#, labelpad=-35)

# 4/2km anomalies
Atmin=-.6
Atmax=.6
Asmin=-0.15
Asmax=0.15

ax3 = fig.add_subplot(323, projection=proj)
ct2=plt.pcolormesh(lon_rho_10km,lat_rho_10km,(bot_temp_ann_4km_interp-np.squeeze(temp_ann_10km[0,:,:])), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=Atmin, vmax=Atmax)
cv=ax3.quiver(xs4[::6,::6],ys4[::6,::6],u4_bot[::6,::6]*mask_zice_4km[::6,::6],v4_bot[::6,::6]*mask_zice_4km[::6,::6], color='k', transform=src, scale=2.)
cv=ax3.quiver(xs4[::6,::6],ys4[::6,::6],u4_bot[::6,::6]*mask_outice_4km[::6,::6],v4_bot[::6,::6]*mask_outice_4km[::6,::6], color='k', transform=src, scale=1.)
ax3.quiver(-.7e6,1.2e6,.2,0,color='w', transform=src, scale=2., zorder=4)
ax3.text(-.7e6,1.25e6,'off-cavity \n 0.2 m/s',color='w', fontsize=8, transform=src, zorder=4)
ax3.quiver(-.8e6,.1e6,.1,0,color='w', transform=src, scale=1., zorder=4)
ax3.text(-.8e6,.15e6,'sub-cavity \n 0.1 m/s',color='w', fontsize=8, transform=src, zorder=4)
plt.title('C) WAOM4 - WAOM10')
ax3.gridlines() # draw_labels=True,linewidth=
x_left, x_right = ax3.get_xlim()
y_low, y_high = ax3.get_ylim()
#ax3.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax3.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='gray', alpha=.65)
ax3.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
ax3.add_feature(bathym, zorder=3, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)

ax4 = fig.add_subplot(324, projection=proj)
cs2=plt.pcolormesh(lon_rho_10km,lat_rho_10km,(sfc_salt_ann_4km_interp-np.squeeze(salt_ann_10km[-1,:,:])), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=Asmin, vmax=Asmax)
cv=ax4.quiver(xs4[::6,::6],ys4[::6,::6],u4_sfc[::6,::6]*mask_zice_4km[::6,::6],v4_sfc[::6,::6]*mask_zice_4km[::6,::6], color='k', transform=src, scale=2.)
cv=ax4.quiver(xs4[::6,::6],ys4[::6,::6],u4_sfc[::6,::6]*mask_outice_4km[::6,::6],v4_sfc[::6,::6]*mask_outice_4km[::6,::6], color='k', transform=src, scale=1.)
ax4.quiver(-.7e6,1.2e6,.2,0,color='w', transform=src, scale=2., zorder=4)
ax4.text(-.7e6,1.25e6,'off-cavity \n 0.2 m/s',color='w', fontsize=8, transform=src, zorder=4)
ax4.quiver(-.8e6,.1e6,.1,0,color='w', transform=src, scale=1., zorder=4)
ax4.text(-.8e6,.15e6,'sub-cavity \n 0.1 m/s',color='w', fontsize=8, transform=src, zorder=4)
plt.title('D) WAOM4 - WAOM10')
ax4.gridlines()
x_left, x_right = ax4.get_xlim()
y_low, y_high = ax4.get_ylim()
#ax4.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax4.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='gray', alpha=.65)
ax4.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
ax4.add_feature(bathym, zorder=3, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)

ax5 = fig.add_subplot(325, projection=proj)
ct2=plt.pcolormesh(lon_rho_10km,lat_rho_10km,(bot_temp_ann_2km_interp-np.squeeze(temp_ann_10km[0,:,:])), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=Atmin, vmax=Atmax)
cv=ax5.quiver(xs2[::12,::12],ys2[::12,::12],u2_bot[::12,::12]*mask_zice_2km[::12,::12],v2_bot[::12,::12]*mask_zice_2km[::12,::12], color='k', transform=src, scale=2.)
cv=ax5.quiver(xs2[::12,::12],ys2[::12,::12],u2_bot[::12,::12]*mask_outice_2km[::12,::12],v2_bot[::12,::12]*mask_outice_2km[::12,::12], color='k', transform=src, scale=1.)
ax5.quiver(-.7e6,1.2e6,.2,0,color='w', transform=src, scale=2., zorder=4)
ax5.text(-.7e6,1.25e6,'off-cavity \n 0.2 m/s',color='w', fontsize=8, transform=src, zorder=4)
ax5.quiver(-.8e6,.1e6,.1,0,color='w', transform=src, scale=1., zorder=4)
ax5.text(-.8e6,.15e6,'sub-cavity \n 0.1 m/s',color='w', fontsize=8, transform=src, zorder=4)
plt.title('E) WAOM2 - WAOM10')
ax5.gridlines() # draw_labels=True,linewidth=
x_left, x_right = ax5.get_xlim()
y_low, y_high = ax5.get_ylim()
#ax5.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax5.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='gray', alpha=.65)
ax5.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
ax5.add_feature(bathym, zorder=3, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
ii=0
plt.plot([lon_rho_2km[eta_sec_ini[ii],xi_pt[ii]], lon_rho_2km[eta_sec_end[ii],xi_pt[ii]]],[lat_rho_2km[eta_sec_ini[ii],xi_pt[ii]], lat_rho_2km[eta_sec_end[ii],xi_pt[ii]]],color='yellow',transform=ccrs.PlateCarree(),linewidth=1.5, linestyle='solid')
# yellow dashed line represents the outer limit for WMT TS and TS-diag.:
plt.plot([lon_rho_2km[eta_sec_ini[ii],xi_pt[ii]-25], lon_rho_2km[eta_sec_end[ii],xi_pt[ii]-25]],[lat_rho_2km[eta_sec_ini[ii],xi_pt[ii]-25], lat_rho_2km[eta_sec_end[ii],xi_pt[ii]-25]],color='yellow',transform=ccrs.PlateCarree(),linewidth=1.5, linestyle='dashed')
plt.plot([lon_rho_2km[eta_sec_ini[ii],xi_pt[ii]+25], lon_rho_2km[eta_sec_end[ii],xi_pt[ii]+25]],[lat_rho_2km[eta_sec_ini[ii],xi_pt[ii]+25], lat_rho_2km[eta_sec_end[ii],xi_pt[ii]+25]],color='yellow',transform=ccrs.PlateCarree(),linewidth=1.5, linestyle='dashed')
ii=1
plt.plot([lon_rho_2km[eta_sec_ini[ii],xi_pt[ii]], lon_rho_2km[eta_sec_end[ii],xi_pt[ii]]],[lat_rho_2km[eta_sec_ini[ii],xi_pt[ii]], lat_rho_2km[eta_sec_end[ii],xi_pt[ii]]],color='yellow',transform=ccrs.PlateCarree(),linewidth=1.5, linestyle='solid')
plt.plot([lon_rho_2km[eta_sec_ini[ii],xi_pt[ii]-25], lon_rho_2km[eta_sec_end[ii],xi_pt[ii]-25]],[lat_rho_2km[eta_sec_ini[ii],xi_pt[ii]-25], lat_rho_2km[eta_sec_end[ii],xi_pt[ii]-25]],color='yellow',transform=ccrs.PlateCarree(),linewidth=1.5, linestyle='dashed')
plt.plot([lon_rho_2km[eta_sec_ini[ii],xi_pt[ii]+25], lon_rho_2km[eta_sec_end[ii],xi_pt[ii]+25]],[lat_rho_2km[eta_sec_ini[ii],xi_pt[ii]+25], lat_rho_2km[eta_sec_end[ii],xi_pt[ii]+25]],color='yellow',transform=ccrs.PlateCarree(),linewidth=1.5, linestyle='dashed')

ax6 = fig.add_subplot(326, projection=proj)
cs2=plt.pcolormesh(lon_rho_10km,lat_rho_10km,(sfc_salt_ann_2km_interp-np.squeeze(salt_ann_10km[-1,:,:])), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=Asmin, vmax=Asmax)
cv=ax6.quiver(xs2[::12,::12],ys2[::12,::12],u2_sfc[::12,::12]*mask_zice_2km[::12,::12],v2_sfc[::12,::12]*mask_zice_2km[::12,::12], color='k', transform=src, scale=2.)
cv=ax6.quiver(xs2[::12,::12],ys2[::12,::12],u2_sfc[::12,::12]*mask_outice_2km[::12,::12],v2_sfc[::12,::12]*mask_outice_2km[::12,::12], color='k', transform=src, scale=1.)
ax6.quiver(-.7e6,1.2e6,.2,0,color='w', transform=src, scale=2., zorder=4)
ax6.text(-.7e6,1.25e6,'off-cavity \n 0.2 m/s',color='w', fontsize=8, transform=src, zorder=4)
ax6.quiver(-.8e6,.1e6,.1,0,color='w', transform=src, scale=1., zorder=4)
ax6.text(-.8e6,.15e6,'sub-cavity \n 0.1 m/s',color='w', fontsize=8, transform=src, zorder=4)
plt.title('F) WAOM2 - WAOM10')
ax6.gridlines()
x_left, x_right = ax6.get_xlim()
y_low, y_high = ax6.get_ylim()
#ax6.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax6.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='gray', alpha=.65)
ax6.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
ax6.add_feature(bathym, zorder=3, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
ii=0
plt.plot([lon_rho_2km[eta_sec_ini[ii],xi_pt[ii]], lon_rho_2km[eta_sec_end[ii],xi_pt[ii]]],[lat_rho_2km[eta_sec_ini[ii],xi_pt[ii]], lat_rho_2km[eta_sec_end[ii],xi_pt[ii]]],color='yellow',transform=ccrs.PlateCarree(),linewidth=1.5, linestyle='solid')
# yellow dashed line represents the outer limit for WMT TS and TS-diag.:
plt.plot([lon_rho_2km[eta_sec_ini[ii],xi_pt[ii]-25], lon_rho_2km[eta_sec_end[ii],xi_pt[ii]-25]],[lat_rho_2km[eta_sec_ini[ii],xi_pt[ii]-25], lat_rho_2km[eta_sec_end[ii],xi_pt[ii]-25]],color='yellow',transform=ccrs.PlateCarree(),linewidth=1.5, linestyle='dashed')
plt.plot([lon_rho_2km[eta_sec_ini[ii],xi_pt[ii]+25], lon_rho_2km[eta_sec_end[ii],xi_pt[ii]+25]],[lat_rho_2km[eta_sec_ini[ii],xi_pt[ii]+25], lat_rho_2km[eta_sec_end[ii],xi_pt[ii]+25]],color='yellow',transform=ccrs.PlateCarree(),linewidth=1.5, linestyle='dashed')
ii=1
plt.plot([lon_rho_2km[eta_sec_ini[ii],xi_pt[ii]], lon_rho_2km[eta_sec_end[ii],xi_pt[ii]]],[lat_rho_2km[eta_sec_ini[ii],xi_pt[ii]], lat_rho_2km[eta_sec_end[ii],xi_pt[ii]]],color='yellow',transform=ccrs.PlateCarree(),linewidth=1.5, linestyle='solid')
plt.plot([lon_rho_2km[eta_sec_ini[ii],xi_pt[ii]-25], lon_rho_2km[eta_sec_end[ii],xi_pt[ii]-25]],[lat_rho_2km[eta_sec_ini[ii],xi_pt[ii]-25], lat_rho_2km[eta_sec_end[ii],xi_pt[ii]-25]],color='yellow',transform=ccrs.PlateCarree(),linewidth=1.5, linestyle='dashed')
plt.plot([lon_rho_2km[eta_sec_ini[ii],xi_pt[ii]+25], lon_rho_2km[eta_sec_end[ii],xi_pt[ii]+25]],[lat_rho_2km[eta_sec_ini[ii],xi_pt[ii]+25], lat_rho_2km[eta_sec_end[ii],xi_pt[ii]+25]],color='yellow',transform=ccrs.PlateCarree(),linewidth=1.5, linestyle='dashed')

cbar_ax1 = fig.add_axes([0.12, 0.08, 0.35, 0.01])
fig.colorbar(ct2, cax=cbar_ax1, orientation='horizontal')
cbar_ax1.set_xlabel('$^{\circ}$C')#, labelpad=-35)

cbar_ax2 = fig.add_axes([0.55, 0.08, 0.35, 0.01])
fig.colorbar(cs2, cax=cbar_ax2, orientation='horizontal')
cbar_ax2.set_xlabel('')#, labelpad=-35)

name_fig="waom10x4x2extend_shflim_S_0.25Q_bottomT_sfcS+Vel_RFIS_anomalies.png"
plt.savefig(fig_path + name_fig, dpi=300)
plt.close()

# only velocity fields:
# velocities magnitude
V10_sfc = np.sqrt(u10_sfc*u10_sfc + v10_sfc*v10_sfc)
V10_bot = np.sqrt(u10_bot*u10_bot + v10_bot*v10_bot)
V4_sfc = np.sqrt(u4_sfc*u4_sfc + v4_sfc*v4_sfc)
V4_bot = np.sqrt(u4_bot*u4_bot + v4_bot*v4_bot)
V2_sfc = np.sqrt(u2_sfc*u2_sfc + v2_sfc*v2_sfc)
V2_bot = np.sqrt(u2_bot*u2_bot + v2_bot*v2_bot)
# sensitivity expt.:
V4_sfcA = np.sqrt(u4_sfcA*u4_sfcA + v4_sfcA*v4_sfcA)
V4_botA = np.sqrt(u4_botA*u4_botA + v4_botA*v4_botA)
V4_sfcB = np.sqrt(u4_sfcB*u4_sfcB + v4_sfcB*v4_sfcB)
V4_botB = np.sqrt(u4_botB*u4_botB + v4_botB*v4_botB)



fig = plt.figure(figsize=(10,8))

ax1 = fig.add_subplot(231, projection=proj)
ct1=plt.pcolormesh(lon_rho_10km,lat_rho_10km,V10_sfc, transform=ccrs.PlateCarree(), cmap=plt.cm.GnBu, vmin=0, vmax=.5)
cv=ax1.quiver(xs10[::3,::3],ys10[::3,::3],u10_sfc[::3,::3],v10_sfc[::3,::3], color='k', transform=src, scale=1.)
plt.title('Surface velocities \n WAOM10')
ax1.gridlines() # draw_labels=True,linewidth=
ax1.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
ax1.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
ax1.add_feature(bathym, zorder=3, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='gray', alpha=.65)

ax2 = fig.add_subplot(232, projection=proj)
ct1=plt.pcolormesh(lon_rho_4km,lat_rho_4km,V4_sfc, transform=ccrs.PlateCarree(), cmap=plt.cm.GnBu, vmin=0, vmax=.5)
cv=ax2.quiver(xs4[::6,::6],ys4[::6,::6],u4_sfc[::6,::6],v4_sfc[::6,::6], color='k', transform=src, scale=1.)
plt.title('Surface velocities \n WAOM4')
ax2.gridlines() # draw_labels=True,linewidth=
ax2.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
ax2.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
ax2.add_feature(bathym, zorder=3, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='gray', alpha=.65)

ax3 = fig.add_subplot(233, projection=proj)
ct1=plt.pcolormesh(lon_rho_2km,lat_rho_2km,V2_sfc, transform=ccrs.PlateCarree(), cmap=plt.cm.GnBu, vmin=0, vmax=.5)
cv=ax3.quiver(xs2[::12,::12],ys2[::12,::12],u2_sfc[::12,::12],v2_sfc[::12,::12], color='k', transform=src, scale=1.)
plt.title('Surface velocities \n WAOM2')
ax3.gridlines() # draw_labels=True,linewidth=
ax3.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
ax3.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
ax3.add_feature(bathym, zorder=3, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='gray', alpha=.65)

ax4 = fig.add_subplot(234, projection=proj)
ct1=plt.pcolormesh(lon_rho_10km,lat_rho_10km,V10_bot, transform=ccrs.PlateCarree(), cmap=plt.cm.GnBu, vmin=0, vmax=.5)
cv=ax4.quiver(xs10[::3,::3],ys10[::3,::3],u10_bot[::3,::3],v10_bot[::3,::3], color='k', transform=src, scale=1.)
plt.title('Bottom velocities \n WAOM10')
ax4.gridlines() # draw_labels=True,linewidth=
ax4.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
ax4.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
ax4.add_feature(bathym, zorder=3, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='gray', alpha=.65)

ax5 = fig.add_subplot(235, projection=proj)
ct1=plt.pcolormesh(lon_rho_4km,lat_rho_4km,V4_bot, transform=ccrs.PlateCarree(), cmap=plt.cm.GnBu, vmin=0, vmax=.5)
cv=ax5.quiver(xs4[::6,::6],ys4[::6,::6],u4_bot[::6,::6],v4_bot[::6,::6], color='k', transform=src, scale=1.)
plt.title('Bottom velocities \n WAOM4')
ax5.gridlines() # draw_labels=True,linewidth=
ax5.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
ax5.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
ax5.add_feature(bathym, zorder=3, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='gray', alpha=.65)

ax6 = fig.add_subplot(236, projection=proj)
ct1=plt.pcolormesh(lon_rho_2km,lat_rho_2km,V2_bot, transform=ccrs.PlateCarree(), cmap=plt.cm.GnBu, vmin=0, vmax=.5)
cv=ax6.quiver(xs2[::12,::12],ys2[::12,::12],u2_bot[::12,::12],v2_bot[::12,::12], color='k', transform=src, scale=1.)
plt.title('Bottom velocities \n WAOM2')
ax6.gridlines() # draw_labels=True,linewidth=
ax6.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
ax6.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
ax6.add_feature(bathym, zorder=3, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='gray', alpha=.65)

cbar_ax1 = fig.add_axes([0.12, 0.08, 0.78, 0.01])
fig.colorbar(ct1, cax=cbar_ax1, orientation='horizontal')
cbar_ax1.set_xlabel('m/s')#, labelpad=-35)

name_fig="waom10x4x2extend_shflim_S_0.25Q_bottom_Vel_RFIS.png"
plt.savefig(fig_path + name_fig, dpi=300)
plt.close()


### PLOTS WITH SENSITIVITY EXPTS: A) NOTIDE AND B) 10KM-BATHY

# --- bottom

tmin = -2.7
tmax = 0.5
smin = 34.1
smax = 34.5

fig = plt.figure(figsize=(8,10))

ax1 = fig.add_subplot(321, projection=proj)
ct1=plt.pcolormesh(lon_rho_4km,lat_rho_4km,np.squeeze(temp_ann_4km[0,:,:]), transform=ccrs.PlateCarree(), cmap=plt.cm.plasma, vmin=tmin, vmax=tmax)
cv=ax1.quiver(xs4[::6,::6],ys4[::6,::6],u4_bot[::6,::6]*mask_zice_4km[::6,::6],v4_bot[::6,::6]*mask_zice_4km[::6,::6], color='k', transform=src, scale=2.)
cv=ax1.quiver(xs4[::6,::6],ys4[::6,::6],u4_bot[::6,::6]*mask_outice_4km[::6,::6],v4_bot[::6,::6]*mask_outice_4km[::6,::6], color='k', transform=src, scale=1.)
ax1.quiver(-.7e6,1.2e6,.2,0,color='w', transform=src, scale=2., zorder=4)
ax1.text(-.7e6,1.25e6,'off-cavity \n 0.2 m/s',color='w', fontsize=8, transform=src, zorder=4)
ax1.quiver(-.8e6,.1e6,.1,0,color='w', transform=src, scale=1., zorder=4)
ax1.text(-.8e6,.15e6,'sub-cavity \n 0.1 m/s',color='w', fontsize=8, transform=src, zorder=4)
plt.title('A) Bottom temperature \n WAOM4')
ax1.gridlines() # draw_labels=True,linewidth=
ax1.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
ratio = .9
x_left, x_right = ax1.get_xlim()
y_low, y_high = ax1.get_ylim()
#ax1.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax1.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
ax1.add_feature(bathym, zorder=3, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='gray', alpha=.65)

cbar_ax3 = fig.add_axes([0.04, 0.65, 0.01, 0.25])
fig.colorbar(ct1, cax=cbar_ax3, orientation='vertical')
cbar_ax3.set_ylabel('$^{\circ}$C')#, labelpad=-35)

ax2 = fig.add_subplot(322, projection=proj)
cs1=plt.pcolormesh(lon_rho_4km,lat_rho_4km,np.squeeze(salt_ann_4km[-1,:,:]), transform=ccrs.PlateCarree(), cmap=plt.cm.plasma, vmin=smin, vmax=smax)
cv=ax2.quiver(xs4[::6,::6],ys4[::6,::6],u4_sfc[::6,::6]*mask_zice_4km[::6,::6],v4_sfc[::6,::6]*mask_zice_4km[::6,::6], color='k', transform=src, scale=2.)
cv=ax2.quiver(xs4[::6,::6],ys4[::6,::6],u4_sfc[::6,::6]*mask_outice_4km[::6,::6],v4_sfc[::6,::6]*mask_outice_4km[::6,::6], color='k', transform=src, scale=1.)
ax2.quiver(-.7e6,1.2e6,.2,0,color='w', transform=src, scale=2., zorder=4)
ax2.text(-.7e6,1.25e6,'off-cavity \n 0.2 m/s',color='w', fontsize=8, transform=src, zorder=4)
ax2.quiver(-.8e6,.1e6,.1,0,color='w', transform=src, scale=1., zorder=4)
ax2.text(-.8e6,.15e6,'sub-cavity \n 0.1 m/s',color='w', fontsize=8, transform=src, zorder=4)
plt.title('B) Surface salinity \n WAOM4')
ax2.gridlines()
x_left, x_right = ax2.get_xlim()
y_low, y_high = ax2.get_ylim()
#ax2.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax2.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
ax2.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
ax2.add_feature(bathym, zorder=3, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='gray', alpha=.65)

cbar_ax4 = fig.add_axes([0.92, 0.65, 0.01, 0.25])
fig.colorbar(cs1, cax=cbar_ax4, orientation='vertical')
cbar_ax4.set_xlabel('')#, labelpad=-35)

# 4/2km anomalies
Atmin=-.6
Atmax=.6
Asmin=-0.15
Asmax=0.15

ax3 = fig.add_subplot(323, projection=proj)
ct2=plt.pcolormesh(lon_rho_4km,lat_rho_4km,np.squeeze(temp_ann_4kmA[0,:,:])-np.squeeze(temp_ann_4km[0,:,:]), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=Atmin, vmax=Atmax)
cv=ax3.quiver(xs4[::6,::6],ys4[::6,::6],u4_botA[::6,::6]*mask_zice_4km[::6,::6],v4_botA[::6,::6]*mask_zice_4km[::6,::6], color='k', transform=src, scale=2.)
cv=ax3.quiver(xs4[::6,::6],ys4[::6,::6],u4_botA[::6,::6]*mask_outice_4km[::6,::6],v4_botA[::6,::6]*mask_outice_4km[::6,::6], color='k', transform=src, scale=1.)
ax3.quiver(-.7e6,1.2e6,.2,0,color='w', transform=src, scale=2., zorder=4)
ax3.text(-.7e6,1.25e6,'off-cavity \n 0.2 m/s',color='w', fontsize=8, transform=src, zorder=4)
ax3.quiver(-.8e6,.1e6,.1,0,color='w', transform=src, scale=1., zorder=4)
ax3.text(-.8e6,.15e6,'sub-cavity \n 0.1 m/s',color='w', fontsize=8, transform=src, zorder=4)
plt.title('C) WAOM4-NOTIDE - WAOM4')
ax3.gridlines() # draw_labels=True,linewidth=
x_left, x_right = ax3.get_xlim()
y_low, y_high = ax3.get_ylim()
#ax3.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax3.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
ax3.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
ax3.add_feature(bathym, zorder=3, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='gray', alpha=.65)

ax4 = fig.add_subplot(324, projection=proj)
cs2=plt.pcolormesh(lon_rho_4km,lat_rho_4km,np.squeeze(salt_ann_4kmA[-1,:,:])-np.squeeze(salt_ann_4km[-1,:,:]), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=Asmin, vmax=Asmax)
cv=ax4.quiver(xs4[::6,::6],ys4[::6,::6],u4_sfcA[::6,::6]*mask_zice_4km[::6,::6],v4_sfcA[::6,::6]*mask_zice_4km[::6,::6], color='k', transform=src, scale=2.)
cv=ax4.quiver(xs4[::6,::6],ys4[::6,::6],u4_sfcA[::6,::6]*mask_outice_4km[::6,::6],v4_sfcA[::6,::6]*mask_outice_4km[::6,::6], color='k', transform=src, scale=1.)
ax4.quiver(-.7e6,1.2e6,.2,0,color='w', transform=src, scale=2., zorder=4)
ax4.text(-.7e6,1.25e6,'off-cavity \n 0.2 m/s',color='w', fontsize=8, transform=src, zorder=4)
ax4.quiver(-.8e6,.1e6,.1,0,color='w', transform=src, scale=1., zorder=4)
ax4.text(-.8e6,.15e6,'sub-cavity \n 0.1 m/s',color='w', fontsize=8, transform=src, zorder=4)
plt.title('D) WAOM4-NOTIDE - WAOM4')
ax4.gridlines()
x_left, x_right = ax4.get_xlim()
y_low, y_high = ax4.get_ylim()
#ax4.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax4.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
ax4.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
ax4.add_feature(bathym, zorder=3, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='gray', alpha=.65)

ax5 = fig.add_subplot(325, projection=proj)
ct2=plt.pcolormesh(lon_rho_4km,lat_rho_4km,np.squeeze(temp_ann_4kmB[0,:,:])-np.squeeze(temp_ann_4km[0,:,:]), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=Atmin, vmax=Atmax)
cv=ax5.quiver(xs4[::6,::6],ys4[::6,::6],u4_botB[::6,::6]*mask_zice_4km[::6,::6],v4_botB[::6,::6]*mask_zice_4km[::6,::6], color='k', transform=src, scale=2.)
cv=ax5.quiver(xs4[::6,::6],ys4[::6,::6],u4_botB[::6,::6]*mask_outice_4km[::6,::6],v4_botB[::6,::6]*mask_outice_4km[::6,::6], color='k', transform=src, scale=1.)
ax5.quiver(-.7e6,1.2e6,.2,0,color='w', transform=src, scale=2., zorder=4)
ax5.text(-.7e6,1.25e6,'off-cavity \n 0.2 m/s',color='w', fontsize=8, transform=src, zorder=4)
ax5.quiver(-.8e6,.1e6,.1,0,color='w', transform=src, scale=1., zorder=4)
ax5.text(-.8e6,.15e6,'sub-cavity \n 0.1 m/s',color='w', fontsize=8, transform=src, zorder=4)
plt.title('E) WAOM4-COARSE - WAOM4')
ax5.gridlines() # draw_labels=True,linewidth=
x_left, x_right = ax5.get_xlim()
y_low, y_high = ax5.get_ylim()
#ax5.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax5.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
ax5.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
ax5.add_feature(bathym, zorder=3, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='gray', alpha=.65)

ax6 = fig.add_subplot(326, projection=proj)
cs2=plt.pcolormesh(lon_rho_4km,lat_rho_4km,np.squeeze(salt_ann_4kmB[-1,:,:])-np.squeeze(salt_ann_4km[-1,:,:]), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=Asmin, vmax=Asmax)
cv=ax6.quiver(xs4[::6,::6],ys4[::6,::6],u4_sfcB[::6,::6]*mask_zice_4km[::6,::6],v4_sfcB[::6,::6]*mask_zice_4km[::6,::6], color='k', transform=src, scale=2.)
cv=ax6.quiver(xs4[::6,::6],ys4[::6,::6],u4_sfcB[::6,::6]*mask_outice_4km[::6,::6],v4_sfcB[::6,::6]*mask_outice_4km[::6,::6], color='k', transform=src, scale=1.)
ax6.quiver(-.7e6,1.2e6,.2,0,color='w', transform=src, scale=2., zorder=4)
ax6.text(-.7e6,1.25e6,'off-cavity \n 0.2 m/s',color='w', fontsize=8, transform=src, zorder=4)
ax6.quiver(-.8e6,.1e6,.1,0,color='w', transform=src, scale=1., zorder=4)
ax6.text(-.8e6,.15e6,'sub-cavity \n 0.1 m/s',color='w', fontsize=8, transform=src, zorder=4)
plt.title('F) WAOM4-COARSE - WAOM4')
ax6.gridlines()
x_left, x_right = ax6.get_xlim()
y_low, y_high = ax6.get_ylim()
#ax6.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
ax6.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
ax6.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
ax6.add_feature(bathym, zorder=3, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='gray', alpha=.65)

cbar_ax1 = fig.add_axes([0.12, 0.08, 0.35, 0.01])
fig.colorbar(ct2, cax=cbar_ax1, orientation='horizontal')
cbar_ax1.set_xlabel('$^{\circ}$C')#, labelpad=-35)

cbar_ax2 = fig.add_axes([0.55, 0.08, 0.35, 0.01])
fig.colorbar(cs2, cax=cbar_ax2, orientation='horizontal')
cbar_ax2.set_xlabel('')#, labelpad=-35)

name_fig="waom4extend_Sensitivity_bottom_TSVel_RFIS_anomalies.png"
plt.savefig(fig_path + name_fig, dpi=300)
plt.close()
