# Ocean Heat Budget Analyses in the Antarctica continental shelf (WAOM)

# Fabio B Dias - 3 May 2023
# Description:
#     this script compares methods of calculate the OHC and surface heat flux integrated over the continental shelf
#     - update to xarray open_mfdataset from open_dataset brought new functions more efficient to calculate model layer thickness (dz)
# and ocean heat content tendencies (rho0*Cp*dT/dt). 

# read nc output from WAOM 10km run

import xarray as xr
# import pandas as p
import numpy as np
import numpy.ma as ma
import cartopy.crs as ccrs
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

import gsw

import pyresample

from dask.distributed import Client
import logging
import warnings
warnings.filterwarnings('ignore')

if __name__== '__main__':

    client = Client(threads_per_worker=1, memory_limit=0, silence_logs=logging.ERROR)
    print(client)

    fig_path = '/g/data/hh5/tmp/access-om/fbd581/ROMS/postprocessing/figs/Contour_isobath/'
    
    # using xr.open_mfdataset
    
    vars2drop = ["ubar","vbar","w","Hsbl","Hbbl","swrad"]
    
    ds = xr.open_mfdataset(paths="/g/data/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom10extend_shflim_S_0.25Q/output_20yr_diag/ocean_avg_00*.nc" , chunks={'eta_rho': '200MB'}, parallel=True, drop_variables=vars2drop, decode_times=False) # , concat_dim="ocean_time"
    
    #- preserving 5-days avgs
    temp = ds.variables["temp"]
    salt = ds.variables["salt"]
    shflux = ds.variables["shflux"]
    ssflux = ds.variables["ssflux"]
    m = ds.variables["m"]
    HvomT = ds.variables["Hvom_temp"]       ## !!! Huon_temp/Hvom_temp were not saved in the original run
    HuonT = ds.variables["Huon_temp"]       ## now it's running here: /scratch/gi0/fbd581/waom4extend_shflim_S_0.25Q/output_yr10_diag
    Hvom = ds.variables["Hvom"]
    Huon = ds.variables["Huon"]
    
    ds = ds.set_coords(['Cs_r', 'Cs_w', 'hc', 'h', 'Vtransform'])
    
    Zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * ds.h) / (ds.hc + ds.h)
    z_rho = ds.zeta + (ds.zeta + ds.h) * Zo_rho + ds.zice
    print("Vtransform=2")
    Zo_w = (ds.hc * ds.s_w + ds.Cs_w * ds.h) / (ds.hc + ds.h)
    z_w = ds.zeta + (ds.zeta + ds.h) * Zo_w + ds.zice
    
    ds.close()
    
    import pickle as pk# load mask of south of the 1500m isobath:
    files_path = '/g/data/hh5/tmp/access-om/fbd581/ROMS/postprocessing/cross_contour_tmp/'
    
    # aaa = np.load(files_path + 'WAOM10_mask_shelf2_1500m.npy', allow_pickle=True)
    
    bbb = np.load(files_path + 'WAOM10_mask_shelf2.mask_1500m', allow_pickle=True)
    
    mask_shelf2 = ma.masked_where((bbb == True), np.ones(bbb.shape))
    ##mask_shelf_1500m = -(bbb-1)
    
    contour_masked_above = np.load(files_path + 'WAOM10_contour_masked_above_1500m', allow_pickle=True)
    contour_masked_above = np.load(files_path + 'WAOM10_contour_masked_above_1500m', allow_pickle=True)
    plt.pcolormesh(ma.masked_where(contour_masked_above == -1000, temp.isel(ocean_time=0, s_rho=0)))
    plt.colorbar()
    
    ds = xr.open_mfdataset(paths="/g/data/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom10extend_shflim_S_0.25Q/output_20yr_diag/ocean_his_00*.nc" , chunks={'eta_rho': '200MB'}, parallel=True, decode_times=False) #, chuncks="auto", concat_dim="ocean_time"
    
    #- preserving 5-days avgs
    temp_snap = ds.variables["temp"] ##+273.15 (changing to Kelvin didn't change any results)
    
    ds.close()
    
    # calculate dT/dt by differentiating temp_snap:
    temp_Rate = np.empty(temp_snap.shape)
    dT = np.empty(temp_snap.shape)
    
    # needs the initial conditions:
    ds = xr.open_dataset('/g/data/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom10extend_shflim_S_0.25Q/output_11-20yr/ocean_rst.nc')
    temp_ini = ds.variables["temp"].isel(ocean_time=8, two=0) ##+273.15 (changing to Kelvin didn't change any results) # 5-days mean
    ds.close()
    
    
    tlen = len(temp[:,0,0,0])
    
    
    # transform to DataArray
    temp_snap = xr.DataArray(temp_snap)
    temp_ini = xr.DataArray(temp_ini)
    
    # - append temp_ini to first time index in temp_snap and then do diff
    temp_snap = xr.concat([temp_ini,temp_snap], 'ocean_time')
    dT = temp_snap.diff('ocean_time')
    print(dT.shape)
    
    #
    dt = 5*86400 # 5-days in seconds
    temp_Rate = np.divide(dT, dt)
    temp_Rate = xr.DataArray(temp_Rate)
    # temp_Rate=temp_Rate.rename({'dim_0':'ocean_time','dim_1':'s_rho','dim_2':'eta_rho','dim_3':'xi_rho'})
    
    print(dT.shape)
    
    temp_Rate = temp_Rate.transpose('ocean_time','s_rho','eta_rho','xi_rho')
    
    print(temp_Rate.shape)
    
    # calculate surface sigma_theta (potential density)
    sigma_t_sfc = gsw.rho(salt[:,-1,:,:],temp[:,-1,:,:],0) - 1000
    
    # load ice draft to create masks
    di = xr.open_dataset('/g/data/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom10extend_shflim_S_0.25Q/output_20yr_diag/ocean_avg_0001.nc')
    ice_draft = di.variables["zice"]
    
    mask_zice = ma.masked_where(ice_draft < 0, np.ones(ice_draft.shape))
    mask_outice = ma.masked_where(ice_draft >= 0, np.ones(ice_draft.shape))
    di.close()
    
    dg = xr.open_dataset("/g/data/hh5/tmp/access-om/fbd581/ROMS/waom10_frc/waom10extend_grd.nc")
    
    lat_rho = dg.variables["lat_rho"]
    lon_rho = dg.variables["lon_rho"]
    lat_u = dg.variables["lat_u"]
    lon_u = dg.variables["lon_u"]
    lat_v = dg.variables["lat_v"]
    lon_v = dg.variables["lon_v"]
    pm = dg.variables["pm"]
    pn = dg.variables["pn"]
    h = dg.variables["h"]
    
    ds.coords['lat_rho']=lat_rho.transpose() # put lat_rho into ds dataset
    ds.coords['lon_rho']=lon_rho.transpose() # put lon_rho into ds dataset
    
    area=np.divide(1,pm*pn)
    
    # method 2 to calculate Dz:
    # z_w=z_w.chunks(chunks={'eta_rho': '200MB'}) # couldn't change chunks.
    
    Z_w = z_w.transpose('ocean_time','s_w','eta_rho','xi_rho')
    print(z_w.shape, Z_w.shape)
    dz = np.diff(Z_w,axis=1)
    
    # determine constants:
    rho0 = 1025 # kg. m-3
    Cp = 3989.245 # J.kg-1.degC-1
    Tf = -1.95 # degC
    
    # use -1000 mask to compute integral of surface heat fluxes and ocean heat content tendency:
    # temp_Rate=xr.DataArray(temp_Rate)
    temp_rate = temp_Rate.transpose('ocean_time','s_rho','eta_rho','xi_rho')
    dT = dT.transpose('ocean_time','s_rho','eta_rho','xi_rho')
    
    ## Integrated heat tendency and sfc heat flux terms to check heat budget closure*1
    # 1. area-integral surface heat flux
    tlen = len(temp_rate[:,0,0,0])
    area_sum =  np.nansum(np.nansum(area,axis=1), axis=0)
    
    shflux_int = np.empty((tlen))
    for mm in np.arange(0,tlen):
        # 20/07/2023
        shflux_masked = ma.masked_where(contour_masked_above == -1000, shflux[mm,:]) # -1000 is just the mask for values south of 1500m isobath; not depths.
        shflux_area = shflux_masked*area
        shflux_int[mm] = np.nansum(np.nansum(shflux_area,axis=1), axis=0)
        del shflux_area
    
    # 2. volume-integral heat tendency
    temp_rate_int = np.empty((tlen))
    temp_rate_vol = np.empty(np.squeeze(temp_Rate[:,0,:,:]).shape)
    for mm in np.arange(0,tlen):
    # - multplying by dz:
        # 20/7/2023:
        temp_rate_dz = temp_Rate[mm,:]*dz[mm,:]
        temp_rate_vint = np.nansum(temp_rate_dz, axis=0)
        temp_rate_vint_masked =  ma.masked_where(contour_masked_above == -1000, temp_rate_vint)
        temp_rate_vol[mm,:] = temp_rate_vint_masked*area
    
        del temp_rate_vint
        temp_rate_int[mm] = np.nansum(np.nansum(ma.masked_where(contour_masked_above == -1000,temp_rate_vol[mm,:]),axis=1), axis=0)*Cp*rho0
    
    # calculate horizontal-integrated area:
    area_masked = ma.masked_where(contour_masked_above == -1000, area)
    area_sum = np.nansum(np.nansum(area_masked,axis=1), axis=0)
    
    print(area_sum*1e-15)
    
    print(temp_Rate.shape, dz.shape)
    
    #print(shflux_area.shape, shflux_int.shape)
    print(np.mean(shflux_int))#/area_sum)
    
    print(temp_rate_vol.shape, temp_rate_int.shape)
    print(np.mean(temp_rate_int))#/area_sum)
    
    mm=35
    #plt.pcolormesh(ma.masked_where(contour_masked_above == -1000,temp_rate_vol[mm,:]), vmin=-2e4, vmax=2e4, cmap='coolwarm')
    plt.colorbar()
    
    #plt.pcolormesh(temp_rate_vol[-1,:], vmin=-2e4, vmax=2e4, cmap='coolwarm')
    #plt.colorbar()
    
    #plt.pcolormesh(shflux_area, cmap='coolwarm')
    #plt.colorbar()
    
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize = (10,5))
    
    ax[0].title.set_text('WAOM10 \n Heat content tendencies')
    # aa=ax[0].pcolormesh(np.nanmean(temp_rate_vol[1:-1], axis=0)*mask_shelf*mask_coast, vmin=-100, vmax=100, cmap='coolwarm')
    # divide by 2.5 to be comparable to 4km maps (normalized by area instead?)
    # aa=ax[0].pcolormesh(np.divide(np.nanmean(temp_rate_vol[1:-1], axis=0),2.5)*mask_shelf*mask_coast, vmin=-100, vmax=100, cmap='coolwarm')
    # divided by area
    aa=ax[0].pcolormesh(np.divide(np.nanmean(temp_rate_vol[1:-1,:], axis=0),area_sum), vmin=-1e-10, vmax=1e-10, cmap='coolwarm')
    # plt.colorbar(aa)
    
    #  define 5-daily masked shflux:
    shflux_masked = np.empty(shflux.shape)
    for mm in np.arange(0,tlen):
        shflux_masked[mm,:,:] = ma.masked_where(contour_masked_above == -1000, shflux[mm,:]) # -1000 is just the mask for values south of 1500m isobath; not depths.
    
    
    ax[1].title.set_text('WAOM10 \n Surface heat flux')
    # bb=ax[1].pcolormesh(np.nanmean(shflux[1:-1], axis=0)*mask_shelf, vmin=-100, vmax=100, cmap='coolwarm')
    bb=ax[1].pcolormesh(np.divide(np.nanmean(shflux_masked[1:-1,:], axis=0),area_sum), vmin=-1e-9, vmax=1e-9, cmap='coolwarm') ##vmin=-100, vmax=100, cmap='coolwarm')
    
    cax2 = plt.axes([0.92, 0.11, 0.01, 0.77])
    cb = plt.colorbar(bb, cax=cax2, orientation='vertical')
    cb.ax.set_ylabel('W.m$^{-2}$', fontsize=12)
    name_fig='WAOM10_OHB_south1500m_temp_rate_vint+shflux_maps_annual.png'
    plt.savefig(fig_path + name_fig, dpi=300)
    
    # OHB integrated
    print('OHC tendency annual avg: ',np.mean(temp_rate_int)*1e-15)
    # print(np.mean(temp_rate2_int)*1e-15)
    print('Net sfc heat flux annual avg: ',np.mean(shflux_int)*1e-15)
    print('Residue (OHC - shflux): ',(np.mean(temp_rate_int)-np.mean(shflux_int))*1e-15)
    
    months=np.arange(0,73)*(5/30.41667)
    
    # save to netcdf file:
    coordinatesC=dict(ocean_time=months)
    
    temp_rate_int_xr = xr.DataArray(temp_rate_int, coords = coordinatesC, dims = ['ocean_time'])
    shflux_int_xr = xr.DataArray(shflux_int, coords = coordinatesC, dims = ['ocean_time'])
    
    files_path = '/g/data/hh5/tmp/access-om/fbd581/ROMS/postprocessing/cross_contour_tmp/'
    temp_rate_int_xr.to_netcdf(files_path + 'WAOM10_OHB_1500m_temp_rate_vint_5daily', mode='w', format="NETCDF4")
    shflux_int_xr.to_netcdf(files_path + 'WAOM10_OHB_1500m_shflux_vint_5daily', mode='w', format="NETCDF4")
    
    fig, ax = plt.subplots(ncols=1, figsize = (10, 6))
    plt.plot(months,shflux_int, label='Sfc heat flux')
    plt.plot(months,temp_rate_int, label='OHC tendency')
    # plt.plot(months,ohc_tend2,'--', label='OHC tendency old')
    plt.plot(months,shflux_int*0,'--k')
    plt.plot(months,temp_rate_int - shflux_int, '--r', label='residual',linewidth=0.5)
    plt.ylim([-2e14,2e14])
    # print annual avg values:
    plt.text(10.2,-.5e14,str(np.round(np.mean(shflux_int)*1e-15,decimals=4)) + 'PW', color='b')
    plt.text(10.2,1.65e14,str(np.round(np.mean(temp_rate_int)*1e-15,decimals=4)) + 'PW', color='darkorange')
    plt.text(9.2,1.1e14,str(np.round(np.mean(temp_rate_int - shflux_int)*1e-15,decimals=4)) + 'PW', color='red')
    
    plt.grid()
    plt.legend()
    plt.ylabel('Heat transport (W)')
    plt.xlabel('Time (months)')
    plt.title('Absolute values; WAOM10')
    name_fig='WAOM10_OHB_south1500m_vint_annual.png'
    plt.savefig(fig_path + name_fig, dpi=300)





