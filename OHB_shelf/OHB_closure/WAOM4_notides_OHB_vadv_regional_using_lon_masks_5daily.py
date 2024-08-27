# Ocean Heat Budget Analyses in the Antarctica continental shelf (WAOM)

# Fabio B Dias - 19 July 2024

# using mask of longitudinal-bins created by 
# - WAOM4extend_shflim_S_0.25Q_OHB_shelf_budget_closure.ipynb
# to integrated over regions.
# Using mask for cont. shelf, land and ice shelf cavities, we
# can separate the heat tendency in each region due to individual
# processes (e.g. advection, diffusion, sfc fluxes)


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
    

    fig_path = '/g/data/hh5/tmp/access-om/fbd581/ROMS/postprocessing/figs/OHB_shelf/'

    # using xr.open_mfdataset
    
    vars2drop = ["ubar","vbar","w","Hsbl","Hbbl","swrad"]
    
    ds = xr.open_mfdataset(paths="/g/data/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom4extend_shflim_S_0.25Q/output_yr10_notides_diag/ocean_avg_00*.nc" , chunks={'eta_rho': '200MB'}, parallel=bool, drop_variables=vars2drop, decode_times=False) # , concat_dim="ocean_time"

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

        # load horizontal diffusion of heat calculated online:
    	# float temp_hdiff(ocean_time, s_rho, eta_rho, xi_rho) ;

    ds = xr.open_mfdataset(paths="/g/data/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom4extend_shflim_S_0.25Q/output_yr10_notides_diag/ocean_dia_00*.nc" , chunks={'eta_rho': '200MB'}, parallel=bool, drop_variables=vars2drop, decode_times=False)
    temp_hdiff = ds.temp_hdiff
    temp_vdiff = ds.temp_vdiff
    temp_hadv = ds.temp_hadv
    temp_vadv = ds.temp_vadv
    ds.close()

    # method 2 to calculate Dz:
    # z_w=z_w.chunks(chunks={'eta_rho': '200MB'}) # couldn't change chunks.

    Z_w = z_w.transpose('ocean_time','s_w','eta_rho','xi_rho')
    print(z_w.shape, Z_w.shape)
    dz = np.diff(Z_w,axis=1)

        # convert dz to xarray:
    months=np.arange(0,73)*(5/30.41667)
    # save to netcdf file:
    coordinatesC=dict(ocean_time=months, s_rho=(['s_rho'], np.arange(0,31)), eta_rho=(['eta_rho'], np.arange(0,1400)), xi_rho=(['xi_rho'], np.arange(0,1575)))

    dz_xr = xr.DataArray(dz, coords = coordinatesC, dims = ['ocean_time','s_rho','eta_rho','xi_rho'])

    heat_vadv_Adz = temp_vadv*dz_xr*area
    heat_hadv_Adz = temp_hadv*dz_xr*area
    heat_vdiff_Adz = temp_vdiff*dz_xr*area
    heat_hdiff_Adz = temp_hdiff*dz_xr*area

    RHS_budget = heat_hadv_Adz + heat_vadv_Adz + heat_hdiff_Adz + heat_vdiff_Adz

    # load contour (1500m/calving-front):
    files_path = '/g/data/hh5/tmp/access-om/fbd581/ROMS/postprocessing/cross_contour_tmp/'
    contour_masked_above = np.load(files_path + 'WAOM4_contour_masked_above_1500m', allow_pickle=True)
    #contour_masked_above_CF = np.load(files_path + 'WAOM4_contour_masked_above_CF', allow_pickle=True)

    mask_shelf = ma.masked_where(contour_masked_above==-1000, np.ones(h.shape))
    #mask_iceshelf = ma.masked_where(contour_masked_above_CF!=-1000, np.ones(h.shape))
    #mask_outiceshelf = ma.masked_where(contour_masked_above_CF==-1000, np.ones(h.shape))

    mask_land = ma.masked_where(h<=40, np.ones(h.shape))

    # define constants:
    rho0=1035
    Cp=3992.1
    # Tf = -1.95 # degC
    Tf =  -3.534879684448242

    # load masks

    dm = xr.open_dataset('/g/data/hh5/tmp/access-om/fbd581/ROMS/postprocessing/cross_contour_tmp/WAOM4_masked_lonBins.nc')
    mask_LonBins = dm.variables['__xarray_dataarray_variable__']
    dm.close()

    # Define conditions to mask out (1) off the shelf, (2) land (see plots above)
    # PS: condtions used for the OHB closure (ice shelf front contour doesn't matter)
    #condition1 = contour_masked_above == -1000
    #condition2 = (zice+h) <= 1

    # Combine conditions with logical OR
    #combined_condition = condition1 | condition2

    # here I need the ice shelf mask obtained from the contour:
    # define condition:
    #lonbin = 0
    # 1. cont.shelf + ice shelf cavities:
    #comb_masks = mask_LonBins[lonbin]*mask_shelf*mask_land
    #condition1 = comb_masks != 1
    # 2. only cont. shelf:
    #comb_masks_noice = mask_LonBins[lonbin]*mask_shelf*mask_iceshelf*mask_land
    #condition2 = comb_masks_noice != 1


    tlen=73
    # heat tendency due to horizontal advection:
    temp_vadv_int = np.empty((len(mask_LonBins),tlen))

    temp_vadv.load()
    #dz.load()

    for ll in np.arange(0,len(mask_LonBins)):
        # 1. cont.shelf + ice shelf cavities:
        comb_masks = mask_LonBins[ll]*mask_shelf*mask_land
        condition1 = comb_masks != 1

        for mm in np.arange(0,tlen):
        # - multplying by dz:
            temp_vadv_dz = temp_vadv[mm,:]*dz[mm,:]
            temp_vadv_vint = np.nansum(temp_vadv_dz, axis=0)
            temp_vadv_vint_masked =  ma.masked_where(condition1, temp_vadv_vint)
            temp_vadv_vol = temp_vadv_vint_masked*area
            temp_vadv_int[ll,mm] = np.nansum(np.nansum(temp_vadv_vol,axis=1), axis=0)*Cp*rho0
            del temp_vadv_dz, temp_vadv_vint, temp_vadv_vint_masked, temp_vadv_vol


    # save into a netcdf file:
    months=np.arange(0,73)*(5/30.41667)
    lon_bin_midpts = np.arange(-178.5,180,3)

    coordinatesC=dict(lon_bins=lon_bin_midpts, ocean_time=months)
    temp_vadv_int_xr = xr.DataArray(temp_vadv_int, coords = coordinatesC, dims = ['lon_bins','ocean_time'])

    files_path = '/g/data/hh5/tmp/access-om/fbd581/ROMS/postprocessing/cross_contour_tmp/'
    temp_vadv_int_xr.to_netcdf(files_path + 'WAOM4_notides_OHB_lonbins_temp_vadv_vint_5daily', mode='w', format="NETCDF4")

