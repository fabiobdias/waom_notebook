# read nc output from WAOM 10km run

import xarray as xr
# import pandas as p
import numpy as np
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

import gsw
import cmocean

import pyresample

from dask.distributed import Client
import logging
import warnings
warnings.filterwarnings('ignore')

if __name__== '__main__':

    client = Client(threads_per_worker=1, memory_limit=0, silence_logs=logging.ERROR)
    print(client)

# using xr.open_mfdataset
    
    vars2drop = ["ubar","vbar","w","Hsbl","Hbbl","swrad"]
    
    ds = xr.open_mfdataset(paths="/g/data/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom4extend_shflim_S_0.25Q/output_yr10_diag_daily/ocean_avg_00*.nc" , chunks={'eta_rho': '200MB'}, parallel=bool, drop_variables=vars2drop, decode_times=False) # , concat_dim="ocean_time"

    m = ds.variables["m"]
    time_avg = ds.variables["ocean_time"] 
    ice_draft = ds.variables["zice"]
    
    mask_zice = ma.masked_where(ice_draft < 0, np.ones(ice_draft.shape))
    mask_outice = ma.masked_where(ice_draft >= 0, np.ones(ice_draft.shape))
   
    # define masks for cavities less_equal than 300m & more than 300m ice draft depth (to separate modes into shallow and deep):

    msk_ice_deep = ma.masked_where(ice_draft>=-300, np.ones(ice_draft.shape))
    msk_ice_shallow = ma.masked_where(ice_draft<-300, np.ones(ice_draft.shape))

    print("Vtransform=2")
    #  New formulation (Vtransform(ng) = 2):
    #
    #         z_w(x,y,s,t) = zeta(x,y,t) + [zeta(x,y,t)+ h(x,y)] * Zo_w
    #                 Zo_w = [hc * s(k) + C(k) * h(x,y)] / [hc + h(x,y)]
    hwater = ds.h- abs(ds.zice) # replace ds.h for hwater below
    Zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * hwater) / (ds.hc + hwater)
    z_rho = ds.zeta + (ds.zeta + hwater) * Zo_rho 
    
    Zo_w = (ds.hc * ds.s_w + ds.Cs_w * hwater) / (ds.hc + hwater)
    z_w = ds.zeta + (ds.zeta + hwater) * Zo_w 
    
    dg = xr.open_dataset("/g/data/hh5/tmp/access-om/fbd581/ROMS/waom4_frc/waom4extend_grd.nc")
    
    lat_rho = dg.variables["lat_rho"]
    lon_rho = dg.variables["lon_rho"]
    lat_u = dg.variables["lat_u"]
    lon_u = dg.variables["lon_u"]
    lat_v = dg.variables["lat_v"]
    lon_v = dg.variables["lon_v"]
    pm = dg.variables["pm"]
    pn = dg.variables["pn"]
    h = dg.variables["h"]
    zice = dg.variables["zice"]
    
    ds.coords['lat_rho']=lat_rho.transpose() # put lat_rho into ds dataset
    ds.coords['lon_rho']=lon_rho.transpose() # put lon_rho into ds dataset
    ds = ds.set_coords(['Cs_r', 'Cs_w', 'hc', 'h', 'Vtransform'])
    
    ds.close()
    dg.close()
    
    # calculate surface sigma_theta (potential density)
    # sigma_t_sfc = gsw.rho(salt[:,-1,:,:],temp[:,-1,:,:],0) - 1000
    area=np.divide(1,pm*pn)

    # define constants:
    rho0=1035
    Cp=3992.1
    # Tf = -1.95 # degC
    Tf =  -3.534879684448242

    # load contour (1500m/calving-front):
    files_path = '/g/data/hh5/tmp/access-om/fbd581/ROMS/postprocessing/cross_contour_tmp/'
    contour_masked_above = np.load(files_path + 'WAOM4_contour_masked_above_1500m', allow_pickle=True)
    contour_masked_above_CF = np.load(files_path + 'WAOM4_contour_masked_above_CF', allow_pickle=True)
    mask_shelf = ma.masked_where(contour_masked_above==-1000, np.ones(h.shape))
    mask_iceshelf = ma.masked_where(contour_masked_above_CF!=-1000, np.ones(h.shape))
    mask_outiceshelf = ma.masked_where(contour_masked_above_CF==-1000, np.ones(h.shape))
    mask_land = ma.masked_where(h<=40, np.ones(h.shape))

    # load mask lon bins:
    dm = xr.open_dataset('/g/data/hh5/tmp/access-om/fbd581/ROMS/postprocessing/cross_contour_tmp/WAOM4_masked_lonlatBins_v2.nc')
    
    mask_LonBins = dm.variables['__xarray_dataarray_variable__']
    dm.close()

    melt_bins_shallow = np.empty((365,141))
    melt_bins_deep = np.empty((365,141))
    for ll in np.arange(141):
        print('------- Debug, figuring out files with new ice draft masks:')
        print('mask_LonBins[ll,:]*mask_shelf*mask_land*msk_ice_shallow shapes =')
        print(mask_LonBins[ll,:].shape, mask_shelf.shape, mask_land.shape, msk_ice_shallow[0,:].shape)

    
        comb_masks_shallow = mask_LonBins[ll,:]*mask_shelf*mask_land*msk_ice_shallow[0,:]
        comb_masks_deep = mask_LonBins[ll,:]*mask_shelf*mask_land*msk_ice_deep[0,:]
        #comb_masks_shallow = mask_LonBins[ll,:]*msk_ice_shallow[0,:]
        #comb_masks_deep = mask_LonBins[ll,:]*msk_ice_deep[0,:]

        condition1_s = comb_masks_shallow != 1
        # apply mask for area_sum too:
        area_masked_s = ma.masked_where(condition1_s, area)
        area_sum_masked_s = np.nansum(np.nansum(area_masked_s,axis=1), axis=0)


        condition1_d = comb_masks_deep != 1
        # apply mask for area_sum too:
        area_masked_d = ma.masked_where(condition1_d, area)
        area_sum_masked_d = np.nansum(np.nansum(area_masked_d,axis=1), axis=0)

        for tt in np.arange(365):
            print('ll, tt =', ll, tt)
            melt_masked_s =  ma.masked_where(condition1_s, m.isel(ocean_time=tt))
            melt_area_s = (melt_masked_s*area)*86400*365.25 # Convert to m /yr
            melt_area_avg_s = np.nansum(np.nansum(melt_area_s,axis=1), axis=0)
            melt_bins_shallow[tt,ll] = melt_area_avg_s  #np.divide(melt_area_avg_s,area_sum_masked_s)
            del melt_masked_s, melt_area_s, melt_area_avg_s

            melt_masked_d =  ma.masked_where(condition1_d, m.isel(ocean_time=tt))
            melt_area_d = (melt_masked_d*area)*86400*365.25 # Convert to m /yr
            melt_area_avg_d = np.nansum(np.nansum(melt_area_d,axis=1), axis=0)
            melt_bins_deep[tt,ll] = melt_area_avg_d    #np.divide(melt_area_avg_d,area_sum_masked_d)
            del melt_masked_d, melt_area_d, melt_area_avg_d

    print('===========> Finished loop days/longitudinal bins!!!!')

    months=np.arange(0,365)*(1/30.41667)

    lon_bin_midpts = np.arange(-178.5,180,3)
    lat_bin_midpts = np.arange(lon_bin_midpts[41],lon_bin_midpts[43],0.6667)
    lon_bin_midpts_WRoss = np.arange(170,172,.5)

    # define lonlat
    lonlat_bin_midpts = np.empty((141))
    lonlat_bin_midpts[0:41] = lon_bin_midpts[0:41]
    lonlat_bin_midpts[41:50] = lat_bin_midpts
    lonlat_bin_midpts[50:60] = np.arange(-50.1664, -49.5, .07)
    lonlat_bin_midpts[60:134] = lon_bin_midpts[43:117]
    lonlat_bin_midpts[134:138] = lon_bin_midpts_WRoss
    lonlat_bin_midpts[138:] = lon_bin_midpts[117:]

    coordinatesC=dict(lon_bins=lonlat_bin_midpts, ocean_time=months)
    melt_bins_shallow_xr = xr.DataArray(melt_bins_shallow, coords = coordinatesC, dims = ['ocean_time','lon_bins'])
    melt_bins_deep_xr = xr.DataArray(melt_bins_deep, coords = coordinatesC, dims = ['ocean_time','lon_bins'])

    files_path = '/g/data/hh5/tmp/access-om/fbd581/ROMS/postprocessing/cross_contour_tmp/'
    melt_bins_shallow_xr.to_netcdf(files_path + 'WAOM4_OHB_lonlatbins-v2_melt_shallow_daily_v2_300m.nc', mode='w', format="NETCDF4")
    melt_bins_deep_xr.to_netcdf(files_path + 'WAOM4_OHB_lonlatbins-v2_melt_deep_daily_v2_300m.nc', mode='w', format="NETCDF4")
    print('===========> Finished saving melt file!!!!')

