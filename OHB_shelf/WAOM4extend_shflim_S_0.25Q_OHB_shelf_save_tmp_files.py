# read nc output from WAOM 10km run

import xarray as xr
import pandas as p
import numpy as np
import numpy.ma as ma
import cartopy.crs as ccrs
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as dates
from matplotlib.ticker import NullFormatter
from matplotlib.colors import LinearSegmentedColormap   # for custom colormaps

from datetime import datetime, timedelta

from netCDF4 import Dataset
from netCDF4 import num2date, date2num
from matplotlib.ticker import NullFormatter
from matplotlib.colors import LinearSegmentedColormap   # for custom colormaps

import gsw

import pyresample

from dask.distributed import Client

fig_path = '/g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/figs/OHB_shelf/'

# starting dask

client = Client()
client

# load ROMS avg output
for mm  in ['01','02','03','04','05','06','07','08','09','10','11','12']:
    ds = xr.open_dataset('/g/data3/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom4extend_shflim_S_0.25Q/output_yr10_diag/ocean_avg_00' + mm + '.nc')
    print(ds.variables["temp"].shape)
    ##- preserving 5-days avgs
    temp_tmp = ds.variables["temp"] # X,31,560,630
    salt_tmp = ds.variables["salt"]
    shflux_tmp = ds.variables["shflux"] # X,560,630
    ssflux_tmp = ds.variables["ssflux"]
    m_tmp = ds.variables["m"]
    # HvomT_tmp = ds.variables["Hvom_temp"] # X,31,559,630      ## !!! Huon_temp/Hvom_temp were not saved in the original run
    # HuonT_tmp = ds.variables["Huon_temp"] # X,31,560,629      ## now it's running here: /scratch/gi0/fbd581/waom4extend_shflim_S_0.25Q/output_yr10_diag
    Hvom_tmp = ds.variables["Hvom"] # X,31,559,630
    Huon_tmp = ds.variables["Huon"] # X,31,560,629

    ds = ds.set_coords(['Cs_r', 'Cs_w', 'hc', 'h', 'Vtransform'])
    if ds.Vtransform == 1:
        Zo_rho = ds.hc * (ds.s_rho - ds.Cs_r) + ds.Cs_r * ds.h
        z_rho_tmp = Zo_rho + ds.zeta * (1 + Zo_rho/ds.h)
        print("Vtransform=1")
    elif ds.Vtransform == 2:
        Zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * ds.h) / (ds.hc + ds.h)
        z_rho_tmp = ds.zeta + (ds.zeta + ds.h) * Zo_rho + ds.zice
        print("Vtransform=2")
        Zo_w = (ds.hc * ds.s_w + ds.Cs_w * ds.h) / (ds.hc + ds.h)
        z_w_tmp = ds.zeta + (ds.zeta + ds.h) * Zo_w + ds.zice

    ##- for 5-days avg:
    z_rho_avg = z_rho_tmp
    z_w_avg = z_w_tmp

    # concatanate monthly avgs into a yearly variable
    if mm == '01':
        temp = temp_tmp # 31,560,630
        salt = salt_tmp
        shflux = shflux_tmp
        ssflux = ssflux_tmp
        m = m_tmp
        z_rho = z_rho_avg
        z_w = z_w_avg
        # HvomT = HvomT_tmp
        # HuonT = HuonT_tmp
        Hvom = Hvom_tmp
        Huon = Huon_tmp
    else:
        temp = np.concatenate((temp,temp_tmp), axis=0) # 2,31,560,630
        salt = np.concatenate((salt,salt_tmp), axis=0)
        shflux = np.concatenate((shflux,shflux_tmp), axis=0)
        ssflux = np.concatenate((ssflux,ssflux_tmp), axis=0)
        m = np.concatenate((m,m_tmp), axis=0)
        z_rho = np.concatenate((z_rho,z_rho_avg), axis=0)
        z_w = np.concatenate((z_w,z_w_avg), axis=0)
        # HvomT = np.concatenate((HvomT,HvomT_tmp), axis=0)
        # HuonT = np.concatenate((HuonT,HuonT_tmp), axis=0)
        Hvom = np.concatenate((Hvom,Hvom_tmp), axis=0)
        Huon = np.concatenate((Huon,Huon_tmp), axis=0)

    ds.close()

# load ROMS avg output
for mm  in ['01','02','03','04','05','06','07','08','09','10','11','12']:
    ds = xr.open_dataset('/g/data3/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom4extend_shflim_S_0.25Q/output_yr10_diag/ocean_his_00' + mm + '.nc')
    print(ds.variables["temp"].shape)
    #- preserving 5-days avgs
    temp_tmp = ds.variables["temp"] # X,31,560,630

    # concatanate monthly avgs into a yearly variable
    if mm == '01':
        temp_snap = temp_tmp # 31,560,630
    else: # only used for monthly averages:
        temp_snap = np.concatenate((temp_snap,temp_tmp), axis=0) # 2,31,560,630

    ds.close()


# --- Save transformation arrays:
npy_path = '/g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/tmp_files/'

print('Saving files .....')

np.savez(npy_path + 'WAOM4extend_5-days_output4OHB_shelf_analyses',temp=temp, salt=salt, shflux=shflux, ssflux=ssflux, m=m, z_rho=z_rho, z_w=z_w, Hvom=Hvom, Huon=Huon, temp_snap=temp_snap)
