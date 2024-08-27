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
    
    ds = xr.open_mfdataset(paths="/g/data/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom4extend_shflim_S_0.25Q/output_yr10_notides_diag_daily/ocean_avg_00*.nc" , chunks={'eta_rho': '200MB'}, parallel=bool, drop_variables=vars2drop, decode_times=False) # , concat_dim="ocean_time"

    m = ds.variables["m"]
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
    dm = xr.open_dataset('/g/data/hh5/tmp/access-om/fbd581/ROMS/postprocessing/cross_contour_tmp/WAOM4_masked_lonBins.nc')
    mask_LonBins = dm.variables['__xarray_dataarray_variable__']
    dm.close()

    melt_bins = np.empty((365,120))
    for ll in np.arange(120):
        comb_masks = mask_LonBins[ll,:]*mask_shelf*mask_land
        condition1 = comb_masks != 1

        for tt in np.arange(365):
            print('ll, tt =', ll, tt)
            melt_masked =  ma.masked_where(condition1, m.isel(ocean_time=tt))
            melt_area = (melt_masked*area)*86400*365.25 # Convert to m /yr
            melt_area_avg = np.nansum(np.nansum(melt_area,axis=1), axis=0)
            melt_bins[tt,ll] = melt_area_avg
            del melt_masked, melt_area, melt_area_avg

    print('===========> Finished loop days/longitudinal bins!!!!')

    months=np.arange(0,365)*(1/30.41667)
    lon_bin_midpts = np.arange(-178.5,180,3)
    coordinatesC=dict(ocean_time=months, lon_bins=lon_bin_midpts)
    melt_bins_xr = xr.DataArray(melt_bins, coords = coordinatesC, dims = ['ocean_time','lon_bins'])

    files_path = '/g/data/hh5/tmp/access-om/fbd581/ROMS/postprocessing/cross_contour_tmp/'
    melt_bins_xr.to_netcdf(files_path + 'WAOM4_notides_OHB_lonbins_melt_daily', mode='w', format="NETCDF4")
    print('===========> Finished saving melt file!!!!')

