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

# !!! sections for Ronne and Ross ice shelves:

# waom2extend sections :
# index for eta var sections (WWed, EWed, ERos, WRos)
xi_pt = [875, 1000, 1275, 1550]
eta_sec_ini = [1600, 1650,  750,  500]
eta_sec_end = [2200, 2250, 1075, 1196]
sec_len = [600, 600, 325, 696]

def read_section_roms(sec_name,ind_ij, ind_len): # 'wed',0,0
#if sec_name == 'WWed' or sec_name ==  'EWed' or sec_name ==  'ERos' or sec_name ==  'Wros':

# vars to export
    variables = ['temp','salt','zeta','z_rho','Hsbl']

# load ROMS avg output
    count = 0
    for mm in ['01','02','03','04','05','06','07','08','09','10','11','12']:
        ds = xr.open_dataset("/g/data3/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom2extend_shflim_S_0.25Q/output_yr5_diag/ocean_avg_00" + mm + ".nc")
        print('size temp and time length: ', ds.temp.shape, len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
        # ----- # first calculate z_rho spatially and put into ds xarray
        ds = ds.set_coords(['Cs_r', 'Cs_w', 'hc', 'h', 'Vtransform'])
        Zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * ds.h) / (ds.hc + ds.h)
        z_rho = ds.zeta + (ds.zeta + ds.h) * Zo_rho + ds.zice
        print("Vtransform=2, z_rho shape", z_rho.shape)
        ds.coords['z_rho'] = z_rho.transpose('ocean_time', 's_rho', 'eta_rho', 'xi_rho')
        
        print('size temp and time length: ', ds.temp.shape, len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
        print(sec_name, 'xi =', xi_pt[ind_ij], 'eta =', eta_sec_ini[ind_ij], eta_sec_end[ind_ij])
        ds[variables].isel(ocean_time=slice(0, 8),xi_rho=xi_pt[ind_ij], eta_rho=slice(eta_sec_ini[ind_ij],eta_sec_end[ind_ij])).to_netcdf('/g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/ncdf_tmp/WAOM2extend_shflim_S_sec_'+ sec_name + '_5days_mm' + mm + '.nc', mode='w', format='NETCDF4')
        
    return 

read_section_roms('WWed',0, 0)

read_section_roms('EWed',1, 1)

#read_section_roms('ERos',2, 2)

#read_section_roms('WRos',3, 3)
