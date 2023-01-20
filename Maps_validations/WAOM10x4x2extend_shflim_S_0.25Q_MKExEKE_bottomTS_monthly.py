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
        print('size temp_tmp = ', temp_tmp.shape)

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

    return temp_tmp, salt_tmp, u10_rho_sfc, v10_rho_sfc, u10_rho_bot, v10_rho_bot

# read ROMS 4km:
def read_roms_ts_4km(exp_path):
    for yr  in ['10']:
        ds = xr.open_dataset(exp_path + 'ocean_avg_00' + yr + '.nc')
        print('size temp and time length: ', ds.temp.shape, len(ds.salt.isel(xi_rho=20, eta_rho=100, s_rho=0)))
        temp_tmp = ds.variables["temp"]
        salt_tmp = ds.variables["salt"]
        print('size temp_tmp = ', temp_tmp.shape)

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

    return temp_tmp, salt_tmp, u4_rho_sfc, v4_rho_sfc, u4_rho_bot, v4_rho_bot


# read ROMS 2km:
def read_roms_ts_2km(exp_path):
    for yr in ['09','10']:
        ds = xr.open_dataset(exp_path + 'ocean_avg_00' + yr + '.nc')
        print('size temp and time length: ', ds.temp.shape, len(ds.salt.isel(xi_rho=20, eta_rho=100, s_rho=0)))
        temp_tmp = ds.variables["temp"]
        salt_tmp = ds.variables["salt"]
        u_tmp = ds.variables["u"]
        v_tmp = ds.variables["v"]

        ds = ds.set_coords(['Cs_r', 'Cs_w', 'hc', 'h', 'Vtransform'])

        if ds.Vtransform == 1:
            Zo_rho = ds.hc * (ds.s_rho - ds.Cs_r) + ds.Cs_r * ds.h
            z_rho_tmp = Zo_rho + ds.zeta * (1 + Zo_rho/ds.h)
            print("Vtransform=1")
        elif ds.Vtransform == 2:
            Zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * ds.h) / (ds.hc + ds.h)
            z_rho_tmp = ds.zeta + (ds.zeta + ds.h) * Zo_rho
            print("Vtransform=2")

        Zo_w = (ds.hc * ds.s_w + ds.Cs_w * ds.h) / (ds.hc + ds.h)
        z_w_tmp = ds.zeta + (ds.zeta + ds.h) * Zo_w


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

        # concantenate annual averaged temp/salt
        if yr == '09':
            temp = temp_tmp 
            salt = salt_tmp
            u2_sfc = u2_rho_sfc
            v2_sfc = v2_rho_sfc
            u2_bot = u2_rho_bot
            v2_bot = v2_rho_bot
        elif yr == '10':
            temp = np.concatenate((temp,temp_tmp), axis=0)
            salt = np.concatenate((salt,salt_tmp), axis=0)
            u2_sfc = np.concatenate((u2_sfc,u2_rho_sfc), axis=0)
            v2_sfc = np.concatenate((v2_sfc,v2_rho_sfc), axis=0)
            u2_bot = np.concatenate((u2_bot,u2_rho_bot), axis=0)
            v2_bot = np.concatenate((v2_bot,v2_rho_bot), axis=0)

    return temp, salt, u2_sfc, v2_sfc, u2_bot, v2_bot

path_ECCO2_10km = '/scratch/project_2000789/boeiradi/waom10extend_shflim_S_0.25Q/output_01-20yr/'
path_ECCO2_4km = '/scratch/project_2000789/boeiradi/waom4extend_shflim_S_0.25Q/output_01-10yr/'
path_ECCO2_2km = '/scratch/project_2000339/boeiradi/waom2extend_shflim_S_0.25Q/output_01-05yr/'

temp_10km, salt_10km, u10_sfc, v10_sfc, u10_bot, v10_bot = read_roms_ts_10km(path_ECCO2_10km)
temp_4km, salt_4km, u4_sfc, v4_sfc, u4_sfc, v4_sfc = read_roms_ts_4km(path_ECCO2_4km)
temp_2km, salt_2km, u2_sfc, v2_sfc, u2_sfc, v2_sfc = read_roms_ts_2km(path_ECCO2_2km)

# define mask zice:
mask_zice_10km = ma.masked_where(zice_10km < 0, np.ones(zice_10km.shape))
mask_zice_4km = ma.masked_where(zice_4km < 0, np.ones(zice_4km.shape))
mask_zice_2km = ma.masked_where(zice_2km < 0, np.ones(zice_2km.shape))
print('Print temp, u, v shapes =',temp_10km.shape, u10_sfc.shape, v10_sfc.shape)

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
print('temp_2km shape', temp_2km.shape)

## ==========================================
# interpolate WAOM4 and WAOM2 field to the WAOM10 grid for maps of difference:

w10_def  = pyresample.geometry.SwathDefinition(lons=lon_rho_10km,lats=lat_rho_10km)
w4_def   = pyresample.geometry.SwathDefinition(lons=lon_rho_4km,lats=lat_rho_4km)
w2_def   = pyresample.geometry.SwathDefinition(lons=lon_rho_2km,lats=lat_rho_2km)
 
wf = lambda r: 1/r

sfc_temp_4km_interp = np.empty((12,560, 630))
bot_temp_4km_interp = np.empty((12,560, 630))
sfc_salt_4km_interp = np.empty((12,560, 630))
bot_salt_4km_interp = np.empty((12,560, 630))

sfc_temp_2km_interp = np.empty((12,560, 630))
bot_temp_2km_interp = np.empty((12,560, 630))
sfc_salt_2km_interp = np.empty((12,560, 630))
bot_salt_2km_interp = np.empty((12,560, 630))

print('Before interpolation, sizes: w4_def,np.array(temp_4km[mm,-1,:,:]),w10_def =', w4_def.shape, np.squeeze(temp_4km[0,-1,:,:]).shape, w10_def.shape)

for mm in np.arange(0,12,1):
    sfc_temp_4km_interp[mm,:,:] = pyresample.kd_tree.resample_custom(w4_def,np.array(temp_4km[mm,-1,:,:]),w10_def,\
                                             radius_of_influence=30000,neighbours=4,weight_funcs=wf)
    sfc_salt_4km_interp[mm,:,:] = pyresample.kd_tree.resample_custom(w4_def,np.array(salt_4km[mm,-1,:,:]),w10_def,\
                                             radius_of_influence=30000,neighbours=4,weight_funcs=wf)
    bot_temp_4km_interp[mm,:,:] = pyresample.kd_tree.resample_custom(w4_def,np.array(temp_4km[mm,0,:,:]),w10_def,\
                                             radius_of_influence=30000,neighbours=4,weight_funcs=wf)
    bot_salt_4km_interp[mm,:,:] = pyresample.kd_tree.resample_custom(w4_def,np.array(salt_4km[mm,0,:,:]),w10_def,\
                                             radius_of_influence=30000,neighbours=4,weight_funcs=wf)

    sfc_temp_2km_interp[mm,:,:] = pyresample.kd_tree.resample_custom(w2_def,np.array(temp_2km[mm,-1,:,:]),w10_def,\
                                             radius_of_influence=5000,neighbours=1,weight_funcs=wf)
    sfc_salt_2km_interp[mm,:,:] = pyresample.kd_tree.resample_custom(w2_def,np.array(salt_2km[mm,-1,:,:]),w10_def,\
                                             radius_of_influence=5000,neighbours=1,weight_funcs=wf)
    bot_temp_2km_interp[mm,:,:] = pyresample.kd_tree.resample_custom(w2_def,np.array(temp_2km[mm,0,:,:]),w10_def,\
                                             radius_of_influence=5000,neighbours=1,weight_funcs=wf)
    bot_salt_2km_interp[mm,:,:] = pyresample.kd_tree.resample_custom(w2_def,np.array(salt_2km[mm,0,:,:]),w10_def,\
                                             radius_of_influence=5000,neighbours=1,weight_funcs=wf)

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

# ========== subregions plots:

# for 2000m isoline
bathym = cfeature.NaturalEarthFeature(name='bathymetry_J_1000', scale='10m', category='physical')
# limits for contour of ice front (Ronne-Filchner IS):
xlimit = np.arange(300,500,1)
ylimit = np.arange(100,300,1)

tmin = -2.4
tmax = 0.
smin = 33.8
smax = 34.5

proj = ccrs.SouthPolarStereo()

for mm in np.arange(0,12,1):

    fig = plt.figure(figsize=(8,10))
    
    ax1 = fig.add_subplot(321, projection=proj)
    ct1=plt.pcolormesh(lon_rho_10km,lat_rho_10km,np.squeeze(temp_10km[mm,-1,:,:]), transform=ccrs.PlateCarree(), cmap=plt.cm.inferno, vmin=tmin, vmax=tmax)
    cv=ax1.quiver(xs10[::6,::6],ys10[::6,::6],u10_sfc[::6,::6],v10_sfc[::6,::6], color='gray', transform=src)
    plt.title('Surface temperature \n WAOM10')
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
    ax1.add_feature(bathym, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
    plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
    plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='grey')
    
    cbar_ax3 = fig.add_axes([0.04, 0.65, 0.01, 0.25])
    fig.colorbar(ct1, cax=cbar_ax3, orientation='vertical')
    cbar_ax3.set_xlabel('$^{\circ}$C')#, labelpad=-35)
    
    ax2 = fig.add_subplot(322, projection=proj)
    cs1=plt.pcolormesh(lon_rho_10km,lat_rho_10km,np.squeeze(salt_10km[mm,-1,:,:]), transform=ccrs.PlateCarree(), cmap=plt.cm.inferno, vmin=smin, vmax=smax)
    cv=ax2.quiver(xs10[::6,::6],ys10[::6,::6],u10_sfc[::6,::6],v10_sfc[::6,::6], color='gray', transform=src)
    plt.title('Surface salinity \n WAOM10')
    ax2.gridlines()
    x_left, x_right = ax2.get_xlim()
    y_low, y_high = ax2.get_ylim()
    #ax2.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
    ax2.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
    ax2.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
    ax2.add_feature(bathym, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
    plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
    plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='grey')
    
    cbar_ax4 = fig.add_axes([0.92, 0.65, 0.01, 0.25])
    fig.colorbar(cs1, cax=cbar_ax4, orientation='vertical')
    cbar_ax4.set_xlabel('')#, labelpad=-35)
    
    # 4/2km anomalies
    Atmin=-.5
    Atmax=.5
    Asmin=-0.3
    Asmax=0.3
    
    ax3 = fig.add_subplot(323, projection=proj)
    ct2=plt.pcolormesh(lon_rho_10km,lat_rho_10km,(sfc_temp_4km_interp[mm,:,:]-np.squeeze(temp_10km[mm,-1,:,:])), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=Atmin, vmax=Atmax)
    cv=ax3.quiver(xs4[::15,::15],ys4[::15,::15],u4_sfc[::15,::15],v4_sfc[::15,::15], color='gray', transform=src)
    plt.title('WAOM4')
    ax3.gridlines() # draw_labels=True,linewidth=
    x_left, x_right = ax3.get_xlim()
    y_low, y_high = ax3.get_ylim()
    #ax3.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
    ax3.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
    ax3.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
    ax3.add_feature(bathym, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
    plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
    plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='grey')
    
    ax4 = fig.add_subplot(324, projection=proj)
    cs2=plt.pcolormesh(lon_rho_10km,lat_rho_10km,(sfc_salt_4km_interp[mm,:,:]-np.squeeze(salt_10km[mm,-1,:,:])), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=Asmin, vmax=Asmax)
    cv=ax4.quiver(xs4[::15,::15],ys4[::15,::15],u4_sfc[::15,::15],v4_sfc[::15,::15], color='gray', transform=src)
    plt.title('WAOM4')
    ax4.gridlines()
    x_left, x_right = ax4.get_xlim()
    y_low, y_high = ax4.get_ylim()
    #ax4.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
    ax4.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
    ax4.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
    ax4.add_feature(bathym, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
    plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
    plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='grey')
    
    ax5 = fig.add_subplot(325, projection=proj)
    ct2=plt.pcolormesh(lon_rho_10km,lat_rho_10km,(sfc_temp_2km_interp[mm,:,:]-np.squeeze(temp_10km[mm,-1,:,:])), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=Atmin, vmax=Atmax)
    cv=ax5.quiver(xs2[::30,::30],ys2[::30,::30],u2_sfc[::30,::30],v2_sfc[::30,::30], color='gray', transform=src)
    plt.title('WAOM2')
    ax5.gridlines() # draw_labels=True,linewidth=
    x_left, x_right = ax5.get_xlim()
    y_low, y_high = ax5.get_ylim()
    #ax5.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
    ax5.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
    ax5.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
    ax5.add_feature(bathym, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
    plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
    plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='grey')
    
    ax6 = fig.add_subplot(326, projection=proj)
    cs2=plt.pcolormesh(lon_rho_10km,lat_rho_10km,(sfc_salt_2km_interp[mm,:,:]-np.squeeze(salt_10km[mm,-1,:,:])), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=Asmin, vmax=Asmax)
    cv=ax6.quiver(xs2[::30,::30],ys2[::30,::30],u2_sfc[::30,::30],v2_sfc[::30,::30], color='gray', transform=src)
    plt.title('WAOM2')
    ax6.gridlines()
    x_left, x_right = ax6.get_xlim()
    y_low, y_high = ax6.get_ylim()
    #ax6.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
    ax6.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
    ax6.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
    ax6.add_feature(bathym, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
    plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
    plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='grey')
    
    cbar_ax1 = fig.add_axes([0.12, 0.08, 0.35, 0.01])
    fig.colorbar(ct2, cax=cbar_ax1, orientation='horizontal')
    cbar_ax1.set_xlabel('$^{\circ}$C')#, labelpad=-35)
    
    cbar_ax2 = fig.add_axes([0.55, 0.08, 0.35, 0.01])
    fig.colorbar(cs2, cax=cbar_ax2, orientation='horizontal')
    cbar_ax2.set_xlabel('')#, labelpad=-35)
    
    name_fig="waom10x4x2extend_shflim_S_0.25Q_sfc_TS_mm=" + str(mm) + "_RFIS_anomalies.png"
    plt.savefig(fig_path + name_fig, dpi=300)
    plt.close()
    
    
    # --- bottom
    
    tmin = -2.4
    tmax = 0.
    smin = 33.8
    smax = 34.5
    
    fig = plt.figure(figsize=(8,10))
    
    ax1 = fig.add_subplot(321, projection=proj)
    ct1=plt.pcolormesh(lon_rho_10km,lat_rho_10km,np.squeeze(temp_10km[mm,0,:,:]), transform=ccrs.PlateCarree(), cmap=plt.cm.inferno, vmin=tmin, vmax=tmax)
    cv=ax1.quiver(xs10[::6,::6],ys10[::6,::6],u10_bot[::6,::6],v10_bot[::6,::6], color='gray', transform=src)
    plt.title('Bottom temperature \n WAOM10')
    ax1.gridlines() # draw_labels=True,linewidth=
    ax1.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
    ratio = .9
    x_left, x_right = ax1.get_xlim()
    y_low, y_high = ax1.get_ylim()
    #ax1.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
    ax1.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
    ax1.add_feature(bathym, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
    plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
    plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='grey')
    
    cbar_ax3 = fig.add_axes([0.04, 0.65, 0.01, 0.25])
    fig.colorbar(ct1, cax=cbar_ax3, orientation='vertical')
    cbar_ax3.set_xlabel('$^{\circ}$C')#, labelpad=-35)
    
    ax2 = fig.add_subplot(322, projection=proj)
    cs1=plt.pcolormesh(lon_rho_10km,lat_rho_10km,np.squeeze(salt_10km[mm,0,:,:]), transform=ccrs.PlateCarree(), cmap=plt.cm.inferno, vmin=smin, vmax=smax)
    cv=ax2.quiver(xs10[::6,::6],ys10[::6,::6],u10_bot[::6,::6],v10_bot[::6,::6], color='gray', transform=src)
    plt.title('Bottom salinity \n WAOM10')
    ax2.gridlines()
    x_left, x_right = ax2.get_xlim()
    y_low, y_high = ax2.get_ylim()
    #ax2.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
    ax2.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
    ax2.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
    ax2.add_feature(bathym, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
    plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
    plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='grey')
    
    cbar_ax4 = fig.add_axes([0.92, 0.65, 0.01, 0.25])
    fig.colorbar(cs1, cax=cbar_ax4, orientation='vertical')
    cbar_ax4.set_xlabel('')#, labelpad=-35)
    
    # 4/2km anomalies
    Atmin=-.6
    Atmax=.6
    Asmin=-0.2
    Asmax=0.2
    
    ax3 = fig.add_subplot(323, projection=proj)
    ct2=plt.pcolormesh(lon_rho_10km,lat_rho_10km,(bot_temp_4km_interp[mm,:,:]-np.squeeze(temp_10km[mm,0,:,:])), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=Atmin, vmax=Atmax)
    cv=ax3.quiver(xs4[::15,::15],ys4[::15,::15],u4_bot[::15,::15],v4_bot[::15,::15], color='gray', transform=src)
    plt.title('WAOM4')
    ax3.gridlines() # draw_labels=True,linewidth=
    x_left, x_right = ax3.get_xlim()
    y_low, y_high = ax3.get_ylim()
    #ax3.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
    ax3.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
    ax3.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
    ax3.add_feature(bathym, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
    plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
    plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='grey')
    
    ax4 = fig.add_subplot(324, projection=proj)
    cs2=plt.pcolormesh(lon_rho_10km,lat_rho_10km,(bot_salt_4km_interp[mm,:,:]-np.squeeze(salt_10km[mm,0,:,:])), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=Asmin, vmax=Asmax)
    cv=ax4.quiver(xs4[::15,::15],ys4[::15,::15],u4_bot[::15,::15],v4_bot[::15,::15], color='gray', transform=src)
    plt.title('WAOM4')
    ax4.gridlines()
    x_left, x_right = ax4.get_xlim()
    y_low, y_high = ax4.get_ylim()
    #ax4.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
    ax4.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
    ax4.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
    ax4.add_feature(bathym, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
    plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
    plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='grey')
    
    ax5 = fig.add_subplot(325, projection=proj)
    ct2=plt.pcolormesh(lon_rho_10km,lat_rho_10km,(bot_temp_2km_interp[mm,:,:]-np.squeeze(temp_10km[mm,0,:,:])), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=Atmin, vmax=Atmax)
    cv=ax5.quiver(xs2[::30,::30],ys2[::30,::30],u2_bot[::30,::30],v2_bot[::30,::30], color='gray', transform=src)
    plt.title('WAOM2')
    ax5.gridlines() # draw_labels=True,linewidth=
    x_left, x_right = ax5.get_xlim()
    y_low, y_high = ax5.get_ylim()
    #ax5.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
    ax5.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
    ax5.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
    ax5.add_feature(bathym, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
    plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
    plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='grey')
    
    ax6 = fig.add_subplot(326, projection=proj)
    cs2=plt.pcolormesh(lon_rho_10km,lat_rho_10km,(bot_salt_2km_interp[mm,:,:]-np.squeeze(salt_10km[mm,0,:,:])), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, vmin=Asmin, vmax=Asmax)
    cv=ax6.quiver(xs2[::30,::30],ys2[::30,::30],u2_bot[::30,::30],v2_bot[::30,::30], color='gray', transform=src)
    plt.title('WAOM2')
    ax6.gridlines()
    x_left, x_right = ax6.get_xlim()
    y_low, y_high = ax6.get_ylim()
    #ax6.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
    ax6.set_extent([-85, -30, -84, -74], crs=ccrs.PlateCarree())
    ax6.add_feature(cfeature.LAND, zorder=3, edgecolor='white', facecolor='gray')
    ax6.add_feature(bathym, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
    plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],zice_10km[xlimit,ylimit],levels=[-.1],linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
    plt.contour(lon_rho_10km[xlimit,ylimit], lat_rho_10km[xlimit,ylimit],h_10km[xlimit,ylimit],levels=(200,400,600,800,1000), transform=ccrs.PlateCarree(), colors='grey')
    
    cbar_ax1 = fig.add_axes([0.12, 0.08, 0.35, 0.01])
    fig.colorbar(ct2, cax=cbar_ax1, orientation='horizontal')
    cbar_ax1.set_xlabel('$^{\circ}$C')#, labelpad=-35)
    
    cbar_ax2 = fig.add_axes([0.55, 0.08, 0.35, 0.01])
    fig.colorbar(cs2, cax=cbar_ax2, orientation='horizontal')
    cbar_ax2.set_xlabel('')#, labelpad=-35)
    
    name_fig="waom10x4x2extend_shflim_S_0.25Q_bottom_TS_mm=" + str(mm) + "_RFIS_anomalies.png"
    plt.savefig(fig_path + name_fig, dpi=300)
    plt.close()

