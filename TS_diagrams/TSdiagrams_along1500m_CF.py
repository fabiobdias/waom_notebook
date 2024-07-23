# read nc output from WAOM 10km_CF run

import xarray as xr
import pandas as p
import numpy as  np
import numpy.ma as ma
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as dates
from matplotlib.ticker import NullFormatter
from matplotlib.colors import LinearSegmentedColormap   # for custom colormaps

from datetime import datetime, timedelta

from netCDF4 import Dataset
from netCDF4 import num2date, date2num
from matplotlib.ticker import NullFormatter
from matplotlib.colors import LinearSegmentedColormap   # for custom colormaps
import cmocean

import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

import pyresample

import gsw

from dask.distributed import Client
import logging
import warnings
warnings.filterwarnings('ignore')

if __name__== '__main__':

    client = Client(threads_per_worker=1, memory_limit=0, silence_logs=logging.ERROR)
    print(client)

    tmp_files_dir = '/g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/cross_contour_tmp/'
    # determine constants:
    # rho0 = 1025 # kg. m-3
    # Cp = 3989.245 # J.kg-1.degC-1
    # Tf = -1.95 # degC

    # use same values as access-om2-01
    rho0 = 1035 # kg. m-3
    Cp = 3992.1 # J.kg-1.degC-1

    # load ice draft to create masks
    expt = 'WAOM10'

    ds = xr.open_dataset(tmp_files_dir + expt + '_vol_trans_1500m_daily_v3')
    vol_transport_10km = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    # ds = xr.open_dataset(tmp_files_dir + expt + '_heat_trans_1500m_daily_v3')
    ds = xr.open_dataset(tmp_files_dir + expt + '_heat_trans_1500m_daily_v4')
    heat_transport_10km = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_temp_1500m_daily_v3')
    temp_10km = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_salt_1500m_daily_v3')
    salt_10km = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_z_rho_1500m_daily_v4')
    z_rho_10km = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_dist_along_1500m_v3')
    dist_along_10km = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_lon_along_1500m_v3')
    lon_along_10km = ds.variables["one"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_lat_along_1500m_v3')
    lat_along_10km = ds.variables["two"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_dist_indices_1500m_v3')
    distance_indices_10km = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    expt = 'WAOM4'

    ds = xr.open_dataset(tmp_files_dir + expt + '_vol_trans_1500m_daily_v3')
    vol_transport_4km = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    # ds = xr.open_dataset(tmp_files_dir + expt + '_heat_trans_1500m_daily_v3')
    ds = xr.open_dataset(tmp_files_dir + expt + '_heat_trans_1500m_daily_v4')
    heat_transport_4km = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_temp_1500m_daily_v3')
    temp_4km = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_salt_1500m_daily_v3')
    salt_4km = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_z_rho_1500m_daily_v4')
    z_rho_4km = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_dist_along_1500m_v3')
    dist_along_4km = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_lon_along_1500m_v3')
    lon_along_4km = ds.variables["one"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_lat_along_1500m_v3')
    lat_along_4km = ds.variables["two"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_dist_indices_1500m_v3')
    distance_indices_4km = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    expt = 'WAOM4_notides'

    ds = xr.open_dataset(tmp_files_dir + expt + '_vol_trans_1500m_daily_v3')
    vol_transport_4kmNT = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    #ds = xr.open_dataset(tmp_files_dir + expt + '_heat_trans_1500m_daily_v3')
    ds = xr.open_dataset(tmp_files_dir + expt + '_heat_trans_1500m_daily_v4')
    heat_transport_4kmNT = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_temp_1500m_daily_v3')
    temp_4kmNT = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_salt_1500m_daily_v3')
    salt_4kmNT = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_z_rho_1500m_daily_v4') # replace for z_rho_1500m_5daily when it finishes (9/8/23)
    z_rho_4kmNT = ds.variables["__xarray_dataarray_variable__"]
    ds.close()


    # load ice draft to create masks
    expt = 'WAOM10'

    ds = xr.open_dataset(tmp_files_dir + expt + '_vol_trans_CF_daily_v3')
    vol_transport_10km_CF = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    # ds = xr.open_dataset(tmp_files_dir + expt + '_heat_trans_CF_daily_v3')
    ds = xr.open_dataset(tmp_files_dir + expt + '_heat_trans_CF_daily_v4')
    heat_transport_10km_CF = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_temp_CF_daily_v3')
    temp_10km_CF = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_salt_CF_daily_v3')
    salt_10km_CF = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_z_rho_CF_daily_v4')
    z_rho_10km_CF = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_dist_along_CF_v3')
    dist_along_10km_CF = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_lon_along_CF_v3')
    lon_along_10km_CF = ds.variables["one"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_lat_along_CF_v3')
    lat_along_10km_CF = ds.variables["two"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_dist_indices_CF_v3')
    distance_indices_10km_CF = ds.variables["__xarray_dataarray_variable__"]
    ds.close()


    expt = 'WAOM4'

    ds = xr.open_dataset(tmp_files_dir + expt + '_vol_trans_CF_daily_v3')
    vol_transport_4km_CF = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    # ds = xr.open_dataset(tmp_files_dir + expt + '_heat_trans_CF_daily_v3')
    ds = xr.open_dataset(tmp_files_dir + expt + '_heat_trans_CF_daily_v4')
    heat_transport_4km_CF = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_temp_CF_daily_v3')
    temp_4km_CF = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_salt_CF_daily_v3')
    salt_4km_CF = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_z_rho_CF_daily_v4')
    z_rho_4km_CF = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_dist_along_CF_v3')
    dist_along_4km_CF = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_lon_along_CF_v3')
    lon_along_4km_CF = ds.variables["one"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_lat_along_CF_v3')
    lat_along_4km_CF = ds.variables["two"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_dist_indices_CF_v3')
    distance_indices_4km_CF = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    expt = 'WAOM4_notides'

    ds = xr.open_dataset(tmp_files_dir + expt + '_vol_trans_CF_daily_v3')
    vol_transport_4kmNT_CF = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    # ds = xr.open_dataset(tmp_files_dir + expt + '_heat_trans_CF_daily_v3')
    ds = xr.open_dataset(tmp_files_dir + expt + '_heat_trans_CF_daily_v4')
    heat_transport_4kmNT_CF = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_temp_CF_daily_v3')
    temp_4kmNT_CF = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_salt_CF_daily_v3')
    salt_4kmNT_CF = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_z_rho_CF_daily_v4')
    z_rho_4kmNT_CF = ds.variables["__xarray_dataarray_variable__"]
    ds.close()


    ## obtain area/distance:

    [dist_along_axis_10km, Zaxis_10km] = np.meshgrid(dist_along_10km, np.arange(0,31))
    [dist_along_axis_4km, Zaxis_4km] = np.meshgrid(dist_along_4km, np.arange(0,31))

    [dist_along_axis_10km_CF, Zaxis_10k_CFm] = np.meshgrid(dist_along_10km_CF, np.arange(0,31))
    [dist_along_axis_4km_CF, Zaxis_4km_CF] = np.meshgrid(dist_along_4km_CF, np.arange(0,31))

    # define fig path:
    fig_path = '/g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/figs/OHB_shelf/'

    # define density range:

    # rho grid for binning:
    rho_grid2=np.arange(36.2,37.4,0.025) # for sigma-2
    len_rho_grid2=len(rho_grid2)

    # # use salt and temp along the contour to calculate sigma:
    sigma_2_10km = gsw.rho(salt_10km[:,:],temp_10km[:,:],2000) - 1000
    sigma_2_4km = gsw.rho(salt_4km[:,:],temp_4km[:,:],2000) - 1000
    sigma_2_4kmNT = gsw.rho(salt_4kmNT[:,:],temp_4kmNT[:,:],2000) - 1000

    # same but for calving front (use salt and temp along the contour to calculate sigma):
    sigma_2_10km_CF = gsw.rho(salt_10km_CF[:,:],temp_10km_CF[:,:],2000) - 1000
    sigma_2_4km_CF = gsw.rho(salt_4km_CF[:,:],temp_4km_CF[:,:],2000) - 1000
    sigma_2_4kmNT_CF = gsw.rho(salt_4kmNT_CF[:,:],temp_4kmNT_CF[:,:],2000) - 1000

    # DEFINE coldest temp along contour:

    Tf_10km = temp_10km.min().values
    Tf_4km = temp_4km.min().values
    Tf_4kmNT = temp_4kmNT.min().values

    Tf_10km_CF = temp_10km_CF.min().values
    Tf_4km_CF = temp_4km_CF.min().values
    Tf_4kmNT_CF = temp_4kmNT_CF.min().values

    print(Tf_10km, Tf_4km, Tf_4kmNT)
    print(Tf_10km_CF, Tf_4km_CF, Tf_4kmNT_CF)

    Tf = Tf_4kmNT # lowest of them all
    print(Tf)

    ## TS diagram to choose isoypcnals:
    # print(sigma_2_10km.shape, temp_10km.shape, salt_10km.shape)


    # make grid for density contours
    smin = 30 - (0.01 * 30)    #salt_ctrl_subregR.min - (0.01 * salt_ctrl_subregR.min)
    smax = 36. + (0.01 * 36.)    #salt_ctrl_subregR.max + (0.01 * salt_ctrl_subregR.max)
    tmin = -4. + (0.1 * -4.)       #temp_ctrl_subregR.min - (0.1 * temp_ctrl_subregR.max)
    tmax = 5 + (0.1 * 5.)       #temp_ctrl_subregR.max + (0.1 * temp_ctrl_subregR.max)
    print('tmin, tmax, smin, smax sizes=,', tmin, tmax, smin, smax)
    # Calculate how many gridcells we need in the x and y dimensions
    xdim = 30
    ydim = 20
    # Create empty grid of zeros
    dens = np.zeros((ydim,xdim))
    # Create temp and salt vectors of appropiate dimensions
    ti = np.linspace(-4,5,ydim)
    si = np.linspace(30,36,xdim)

    Si, Ti = np.meshgrid(si, ti, sparse=False, indexing='ij')
    # Loop to fill in grid with densities
    for j in range(0,int(ydim)):
        for i in range(0, int(xdim)):
            dens[j,i]=gsw.rho(si[i],ti[j],2000) # sigma-2000
    # Substract 1000 to convert to sigma-2
    dens = dens - 1000

    salt_10km.load()
    temp_10km.load()
    z_rho_10km.load()

    salt_4km.load()
    temp_4km.load()
    z_rho_4km.load()

    salt_4kmNT.load()
    temp_4kmNT.load()
    z_rho_4kmNT.load()

    fig, ax = plt.subplots(nrows=2, ncols=3, figsize = (12, 6))
    # normal isopycnals
    for aa in range(0,2):
        for bb in range(0,3):
            CS1 = ax[aa,bb].contour(Si,Ti,dens.transpose(), levels=np.arange(35.5,38,.1),linestyles='solid', colors=[(.8,0.8,0.8)], linewidth=0.1)
            ax[aa,bb].clabel(CS1, CS1.levels, inline=True, fontsize=10)
            CS2 = ax[aa,bb].contour(Si,Ti,dens.transpose(), levels=np.arange(36.7,36.71),linestyles='dashed', colors='cyan', linewidth=0.1)
            CS3 = ax[aa,bb].contour(Si,Ti,dens.transpose(), levels=np.arange(36.9,36.91),linestyles='dashed', colors='deepskyblue', linewidth=0.1)
            CS4 = ax[aa,bb].contour(Si,Ti,dens.transpose(), levels=np.arange(37.,37.01),linestyles='dashed', colors='blue', linewidth=0.1)
            CS5 = ax[aa,bb].contour(Si,Ti,dens.transpose(), levels=np.arange(37.3,37.31),linestyles='dashed', colors='navy', linewidth=0.1)
            ax[aa,bb].set_xlim([33.2,35])
            ax[aa,bb].set_ylim([-3,2])
            ax[aa,bb].set_ylabel('Temperature')
            ax[aa,bb].set_xlabel('Salinity')

    c = ax[0,0].scatter(salt_10km,temp_10km,
                     c=-z_rho_10km,marker="p", s=(72./fig.dpi)**2,
                     label='WAOM10-1500m', cmap=cmocean.cm.deep, vmin=0, vmax=2000)
    c = ax[0,1].scatter(salt_4km,temp_4km,
                     c=-z_rho_4km,marker="p", s=(72./fig.dpi)**2,
                     label='WAOM4-1500m', cmap=cmocean.cm.deep, vmin=0, vmax=2000)
    c = ax[0,2].scatter(salt_4kmNT,temp_4kmNT,
                     c=-z_rho_4kmNT,marker="p", s=(72./fig.dpi)**2,
                     label='WAOM4NT-1500m', cmap=cmocean.cm.deep, vmin=0, vmax=2000)
    c = ax[1,0].scatter(salt_10km_CF,temp_10km_CF,
                     c=-z_rho_10km_CF,marker="p", s=(72./fig.dpi)**2,
                     label='WAOM10-CF', cmap=cmocean.cm.deep, vmin=0, vmax=2000)
    c = ax[1,1].scatter(salt_4km_CF,temp_4km_CF,
                     c=-z_rho_4km_CF,marker="p", s=(72./fig.dpi)**2,
                     label='WAOM4-CF', cmap=cmocean.cm.deep, vmin=0, vmax=2000)
    c = ax[1,2].scatter(salt_4kmNT_CF,temp_4kmNT_CF,
                     c=-z_rho_4kmNT_CF,marker="p", s=(72./fig.dpi)**2,
                     label='WAOM4NT-CF', cmap=cmocean.cm.deep, vmin=0, vmax=2000)

    cbar_ax1 = fig.add_axes([0.915, 0.12,  0.01, 0.75])
    fig.colorbar(c, cax=cbar_ax1, orientation='vertical')
    cbar_ax1.set_ylabel('Depth (m)')#, labelpad=-35)

    #
    plt.savefig(fig_path + 'WAOM10x4x4NT_Along-1500m_CF_TSdiags_daily.png', bbox_inches='tight', dpi=300)
