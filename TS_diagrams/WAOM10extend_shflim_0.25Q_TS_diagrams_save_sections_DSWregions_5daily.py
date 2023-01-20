import xarray as xr
import pandas as p
import numpy as np
import numpy.ma as ma
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.dates as dates
from matplotlib.ticker import NullFormatter
from matplotlib.colors import LinearSegmentedColormap   # for custom colormaps
import matplotlib as mpl
mpl.use('Agg')

from datetime import datetime, timedelta

import netCDF4
from matplotlib.ticker import NullFormatter
from matplotlib.colors import LinearSegmentedColormap   # for custom colormaps

import iris
import iris.iterate
import iris.coords
import iris.plot as iplt
import gsw

# !!! updated sections following the ones for WAOM10_shflim_S_ORAS5em_0.25Q:

# waom10extend sections DSW:
# index for eta var sections (wed, ros, neu, mai, mer)
xi_pt = [170, 310, 280, 320, 435]
eta_sec_ini = [310, 50, 450, 450, 40]
eta_sec_end = [480, 280, 550, 550, 110]
# index for xi var sections (dar, pry, amu)
xi_sec_ini = [460, 420, 50]
xi_sec_end = [615, 600, 200]
eta_pt = [380, 368, 225]
# length of the section (wed, ros, neu, mai, mer, dar, pry, amu)
sec_len = [170, 230, 100, 100, 70, 155, 180, 150]

def read_section_roms(sec_name,ind_ij, ind_len): # 'wed',0,0
#if sec_name == 'wed' or sec_name ==  'ros' or sec_name ==  'neu' or sec_name ==  'mai' or sec_name ==  'mer':
#elif sec_name == 'pry' or sec_name ==  'dar' or sec_name ==  'amu':

# vars to export
    variables = ['temp','salt','zeta','z_rho']

# load ROMS avg output
    count = 0
    for mm in ['01','02','03','04','05','06','07','08','09','10','11','12']:
        ds = xr.open_dataset("/scratch/project_2000789/boeiradi/waom10extend_shflim_S_0.25Q/output_20yr_diag/ocean_avg_00" + mm + ".nc")
        print('size temp and time length: ', ds.temp.shape, len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
        # ----- # first calculate z_rho spatially and put into ds xarray
        ds = ds.set_coords(['Cs_r', 'Cs_w', 'hc', 'h', 'Vtransform'])
        Zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * ds.h) / (ds.hc + ds.h)
        z_rho = ds.zeta + (ds.zeta + ds.h) * Zo_rho
        print("Vtransform=2, z_rho shape", z_rho.shape)
        ds.coords['z_rho'] = z_rho.transpose('ocean_time', 's_rho', 'eta_rho', 'xi_rho')
        
        print('size temp and time length: ', ds.temp.shape, len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
        if sec_name == 'wed' or sec_name ==  'ros' or sec_name ==  'neu' or sec_name ==  'mai' or sec_name ==  'mer':
            print(sec_name, 'xi =', xi_pt[ind_ij], 'eta =', eta_sec_ini[ind_ij], eta_sec_end[ind_ij])
            ds[variables].isel(ocean_time=slice(0, 8),xi_rho=xi_pt[ind_ij], eta_rho=slice(eta_sec_ini[ind_ij],eta_sec_end[ind_ij])).to_netcdf('/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM10extend_shflim_S_sec_'+ sec_name + '_5days_mm' + mm + '.nc', mode='w', format='NETCDF4')
        elif sec_name == 'pry' or sec_name ==  'dar' or sec_name ==  'amu':
            print(sec_name, 'eta =', eta_pt[ind_ij], 'xi =', xi_sec_ini[ind_ij], xi_sec_end[ind_ij])
            ds[variables].isel(ocean_time=slice(0, 8),xi_rho=slice(xi_sec_ini[ind_ij],xi_sec_end[ind_ij]), eta_rho=eta_pt[ind_ij]).to_netcdf('/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM10extend_shflim_S_sec_'+ sec_name + '_5days_mm' + mm + '.nc', mode='w', format='NETCDF4')

    return 

read_section_roms('wed',0, 0)

read_section_roms('ros',1, 1)

read_section_roms('neu',2, 2)

read_section_roms('mai',3, 3)

read_section_roms('mer',4, 4)

read_section_roms('dar',0, 5)

read_section_roms('pry',1, 6)

read_section_roms('amu',2, 7)
