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

# index for eta var sections (mer)
xi_pt = [1088]
eta_sec_ini = [100]
eta_sec_end = [275]
# index for xi var sections (dar, pry, bar, sha, vin)
xi_sec_ini = [1150, 1150, 1250, 1300, 1300]
xi_sec_end = [1538, 1538, 1550, 1550, 1500]
#eta_pt = [950, 920, 863, 675, 500]
eta_pt = [970, 920, 860, 675, 500] # testing Darnley further West
# length of the section
sec_len = [388, 388, 300, 250, 200, 175]


def read_section_roms(sec_name,ind_ij, ind_len): # 'wed',0,0

# vars to export
    variables = ['temp','salt','zeta','z_rho']

# load ROMS avg output
    count = 0
    for mm in ['01','02','03','04','05','06','07','08','09','10','11','12']:
        ds = xr.open_dataset("/scratch/project_2000789/boeiradi/waom4extend_shflim_S_0.25Q/output_yr10_diag/ocean_avg_00" + mm + ".nc")
        print('size temp and time length: ', ds.temp.shape, len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
        # ----- # first calculate z_rho spatially and put into ds xarray
        ds = ds.set_coords(['Cs_r', 'Cs_w', 'hc', 'h', 'Vtransform'])
        Zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * ds.h) / (ds.hc + ds.h)
        z_rho = ds.zeta + (ds.zeta + ds.h) * Zo_rho
        print("Vtransform=2, z_rho shape", z_rho.shape)
        ds.coords['z_rho'] = z_rho.transpose('ocean_time', 's_rho', 'eta_rho', 'xi_rho')
        
        print('size temp and time length: ', ds.temp.shape, len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
        if sec_name ==  'mer':
            print(sec_name, 'xi =', xi_pt[ind_ij], 'eta =', eta_sec_ini[ind_ij], eta_sec_end[ind_ij])
            ds[variables].isel(ocean_time=slice(0, 8),xi_rho=xi_pt[ind_ij], eta_rho=slice(eta_sec_ini[ind_ij],eta_sec_end[ind_ij])).to_netcdf('/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM4extend_shflim_S_sec_'+ sec_name + '_5days_mm' + mm + '.nc', mode='w', format='NETCDF4')
        elif sec_name == 'pry' or sec_name ==  'dar' or sec_name ==  'bar' or sec_name ==  'sha' or sec_name ==  'vin':
            print(sec_name, 'eta =', eta_pt[ind_ij], 'xi =', xi_sec_ini[ind_ij], xi_sec_end[ind_ij])
            ds[variables].isel(ocean_time=slice(0, 8),xi_rho=slice(xi_sec_ini[ind_ij],xi_sec_end[ind_ij]), eta_rho=eta_pt[ind_ij]).to_netcdf('/users/boeiradi/COLD_project/postprocessing/ncdf_tmp/WAOM4extend_shflim_S_sec_'+ sec_name + '_5days_mm' + mm + '.nc', mode='w', format='NETCDF4')

    return 

read_section_roms('dar',0, 0)
read_section_roms('pry',1, 1)
read_section_roms('bar',2, 2)
read_section_roms('sha',3, 3)
read_section_roms('vin',4, 4)
read_section_roms('mer',0, 5)
