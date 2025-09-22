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

    du = xr.open_mfdataset(paths="/g/data/hh5/tmp/access-om/fbd581/ROMS/waom4_frc/waom4extend_sustr.nc" , chunks={'eta_rho': '200MB'}, parallel=bool, decode_times=False)
    sustr = du.sustr
    du.close()

    dv = xr.open_mfdataset(paths="/g/data/hh5/tmp/access-om/fbd581/ROMS/waom4_frc/waom4extend_svstr.nc" , chunks={'eta_rho': '200MB'}, parallel=bool, decode_times=False)
    svstr = dv.svstr
    dv.close()

# using xr.open_mfdataset
    
    vars2drop = ["ubar","vbar","w","Hsbl","Hbbl","swrad"]
    
    ds = xr.open_mfdataset(paths="/g/data/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom4extend_shflim_S_0.25Q/output_yr10_diag_daily/ocean_avg_00*.nc" , chunks={'eta_rho': '200MB'}, parallel=bool, drop_variables=vars2drop, decode_times=False) # , concat_dim="ocean_time"

    time_avg = ds.variables["ocean_time"] 
    ice_draft = ds.variables["zice"]
    
    mask_zice = ma.masked_where(ice_draft < 0, np.ones(ice_draft.shape))
    mask_outice = ma.masked_where(ice_draft >= 0, np.ones(ice_draft.shape))
    
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
    dm = xr.open_dataset('/g/data/hh5/tmp/access-om/fbd581/ROMS/postprocessing/cross_contour_tmp/WAOM4_masked_lonlatBins.nc')
    
    mask_LonBins = dm.variables['__xarray_dataarray_variable__']
    dm.close()

    print('Finished to read wind and mask datas!')

    # sustr/svstr # waom4
    Ustr = sustr.rename({'eta_u': 'eta','xi_u': 'xi'})*mask_shelf[:,:-1]
    Vstr = svstr.rename({'eta_v': 'eta','xi_v': 'xi'})*mask_shelf[:-1,:]
    Usqr = (Ustr*Ustr)
    Vsqr = (Vstr*Vstr)
    print(Usqr[:,:-1,:].shape, Vsqr[:,:,:-1].shape)
    WSmgn = np.sqrt(Usqr[:,:-1,:]+Vsqr[:,:,:-1])

    # define lon/lat bin mid-points:
    lon_bin_midpts = np.arange(-178.5,180,3)
    lat_bin_midpts = np.arange(lon_bin_midpts[41],lon_bin_midpts[43],0.6667)

    # define lonlat
    lonlat_bin_midpts = np.empty((127))
    lonlat_bin_midpts[0:41] = lon_bin_midpts[0:41]
    lonlat_bin_midpts[41:50] = lat_bin_midpts
    lonlat_bin_midpts[50:] = lon_bin_midpts[43:]

    ds = xr.open_mfdataset(paths="/g/data/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom4extend_shflim_S_0.25Q/output_yr10_diag_daily/ocean_dia_00*.nc" , chunks={'eta_rho': '200MB'}, parallel=bool, drop_variables=vars2drop, decode_times=False)
    temp_hdiff = ds.temp_hdiff
    ds.close()

    ## convert mask_LonBins_4km to dataArray:
    coordinatesC=dict(lon_bins=lonlat_bin_midpts, eta=temp_hdiff.eta_rho, xi=temp_hdiff.xi_rho)
    mask_LonBins_xr = xr.DataArray(mask_LonBins, coords = coordinatesC, dims = ['lon_bins','eta_rho','xi_rho'])

    mask_LonBins_xr=mask_LonBins_xr.rename({'eta_rho': 'eta','xi_rho': 'xi'})

    # define Wind Stress x mask_LonLat_Bins
    WS_masked = (mask_LonBins_xr.isel(eta=slice(0,-1),xi=slice(0,-1))*WSmgn)

    wind_str_bins = np.empty((1,127))
    for ll in np.arange(127):
        comb_masks = mask_LonBins[ll,:-1,:-1]*mask_shelf[:-1,:-1]*mask_land[:-1,:-1]
        print('Bug tracking - comb_masks shape = ', comb_masks.shape)
        condition1 = comb_masks != 1
        # apply mask for area_sum too:
        comb_masks_f = mask_LonBins[ll,:,:]*mask_shelf*mask_land
        condition1_f = comb_masks_f != 1
        area_masked = ma.masked_where(condition1_f, area)
        area_sum_masked = np.nansum(np.nansum(area_masked,axis=1), axis=0)

        #for tt in np.arange(365):
        tt=0
        print('ll, tt =', ll, tt)
        print('Bug tracking - condition1 & WS_masked.isel(sms_time=tt) shapes = ', condition1.shape, WS_masked.isel(sms_time=tt, lon_bins=ll).shape)
        wind_str_masked =  ma.masked_where(condition1, WS_masked.isel(sms_time=tt, lon_bins=ll))
        wind_str_area = (wind_str_masked*area[:-1,:-1])
        wind_str_area_avg = np.nansum(np.nansum(wind_str_area,axis=1), axis=0)
        wind_str_bins[tt,ll] = np.divide(wind_str_area_avg,area_sum_masked)
        del wind_str_masked, wind_str_area, wind_str_area_avg

    print('===========> Finished loop days/longitudinal bins!!!!')

    ##months=np.arange(0,365)*(1/30.41667)
    months=tt*(1/30.41667)
    lon_bin_midpts = np.arange(-178.5,180,3)

    coordinatesC=dict(lon_bins=lonlat_bin_midpts, ocean_time=months)
    wind_str_bins_xr = xr.DataArray(wind_str_bins, coords = coordinatesC, dims = ['ocean_time','lon_bins'])

    files_path = '/g/data/hh5/tmp/access-om/fbd581/ROMS/postprocessing/cross_contour_tmp/'
    wind_str_bins_xr.to_netcdf(files_path + 'WAOM4_OHB_lonlatbins_wind_str_daily_d' + str(tt) + '.nc', mode='w', format="NETCDF4")
    print('===========> Finished saving melt file!!!!')

