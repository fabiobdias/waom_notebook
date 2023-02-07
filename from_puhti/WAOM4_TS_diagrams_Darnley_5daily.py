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

from netCDF4 import Dataset
from netCDF4 import num2date, date2num
from matplotlib.ticker import NullFormatter
from matplotlib.colors import LinearSegmentedColormap   # for custom colormaps

import iris
import iris.iterate
import iris.coords
import iris.plot as iplt
import gsw

# index for eta var sections (mer)
xi_pt = [1075]
eta_sec_ini = [25]
eta_sec_end = [200]
# index for xi var sections (dar, pry, bar, sha, vin)
xi_sec_ini = [1150, 1150, 1250, 1300, 1300]
xi_sec_end = [1538, 1538, 1550, 1550, 1500]
eta_pt = [875, 845, 788, 600, 425]
# length of the section
sec_len = [388, 388, 300, 250, 200, 175]

def read_section_roms(sec_name,ind_ij, ind_len):

    # limits for time index
    monthly_ind_ini = [0, 7, 13, 19, 25, 31, 37, 43, 49, 55, 61, 67]
    monthly_ind_end = [6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 73]
    # vars to export
    variables = ['temp','salt','zeta','z_rho']

# load ROMS avg output
    count = 0
    for mm in ['01','02','03','04','05','06','07','08','09','10','11','12']:
        print('Saving month =', mm)
        ds = xr.open_dataset("/scratch/project_2000339/boeiradi/waom4_mahti/output_6yr_fabio_inifile/ocean_avg_00" + mm + ".nc")
# ----- # first calculate z_rho spatially and put into ds xarray
        ds = ds.set_coords(['Cs_r', 'Cs_w', 'hc', 'h', 'Vtransform'])
        Zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * ds.h) / (ds.hc + ds.h)
        z_rho = ds.zeta + (ds.zeta + ds.h) * Zo_rho
        print("Vtransform=2, z_rho shape", z_rho.shape)
        ds.coords['z_rho'] = z_rho.transpose('ocean_time', 's_rho', 'eta_rho', 'xi_rho')

        print('size temp and time length: ', ds.temp.shape, len(ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
        if sec_name == 'mer':
            print(sec_name, 'xi =', xi_pt[ind_ij], 'eta =', eta_sec_ini[ind_ij], eta_sec_end[ind_ij])
# saving all variables including temp/salt
            ds[variables].isel(ocean_time=slice(0, 8),xi_rho=xi_pt[ind_ij], eta_rho=slice(eta_sec_ini[ind_ij],eta_sec_end[ind_ij])).to_netcdf('/scratch/project_2000339/boeiradi/postprocessing/ncdf_tmp/WAOM4_sec_' + sec_name + '_5days_mm' + mm + '.nc', mode='w', format='NETCDF3_64BIT')
        elif sec_name == 'dar' or sec_name ==  'pry' or sec_name ==  'bar' or sec_name ==  'sha' or sec_name ==  'vin':
            print(sec_name, 'eta =', eta_pt[ind_ij], 'xi =', xi_sec_ini[ind_ij], xi_sec_end[ind_ij])
            # saving all variables including temp/salt
            ds[variables].isel(ocean_time=slice(0, 8),xi_rho=slice(xi_sec_ini[ind_ij],xi_sec_end[ind_ij]), eta_rho=eta_pt[ind_ij]).to_netcdf('/scratch/project_2000339/boeiradi/postprocessing/ncdf_tmp/WAOM4_sec_' + sec_name + '_5days_mm' + mm + '.nc', mode='w', format='NETCDF3_64BIT')

    return

read_section_roms('dar',0, 0)
