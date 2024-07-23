#!/usr/bin/env python3
# Ocean Heat Budget Analyses in the Antarctica continental shelf (WAOM)

# Fabio B Dias - 28 June 2023
# Description:
#     this script obtain and save the 1500m isobath contour variables, which is used for the 
#     cross-shelf heat transport estimates

# read nc output from WAOM 10km run

import xarray as xr
import pandas as p
import numpy as np
import numpy.ma as ma
import cartopy.crs as ccrs
import matplotlib as mpl
#mpl.use('Agg')
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
    
    # load ice draft to create masks
    di = xr.open_dataset('/g/data/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom10extend_shflim_S_0.25Q/output_20yr_diag_daily/ocean_avg_0001.nc')
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
    pm = dg.pm
    pn = dg.pn
    h = dg.variables["h"]
    
    # ds.coords['lat_rho']=lat_rho.transpose() # put lat_rho into ds dataset
    # ds.coords['lon_rho']=lon_rho.transpose() # put lon_rho into ds dataset
    
    area=np.divide(1,pm*pn)
    
    ## creating the contour, such as a isobath, and extracting the coordinates using matplotlib's Path class
    # based on https://github.com/COSIMA/cosima-recipes/blob/master/DocumentedExamples/Cross-contour_transport.ipynb
    
    h = dg.h.load()
    
    h = h*mask_zice
    
    # Fill in land with zeros:
    h = h.fillna(0)
    
    contour_depth = 1500.
    
    ## Choose whether you want your contour on the u or t grid.
    grid_sel = 't'
    if grid_sel == 'u':
        x_var = lon_u
        y_var = lat_u
    elif grid_sel == 't':
        x_var = lon_rho
        y_var = lat_rho

# load contour generated from WAOM10_Extract_1500m_contour.ipynb
    fig_path = '/g/data/hh5/tmp/access-om/fbd581/ROMS/postprocessing/figs/Contour_isobath/'
    xcon_np=np.loadtxt("/g/data/hh5/tmp/access-om/fbd581/ROMS/postprocessing/figs/Contour_isobath/WAOM10_1500m_x_contour.csv")
    x_contour = xcon_np.tolist()
    ycon_np=np.loadtxt("/g/data/hh5/tmp/access-om/fbd581/ROMS/postprocessing/figs/Contour_isobath/WAOM10_1500m_y_contour.csv")
    y_contour = ycon_np.tolist()

    ## SHOULD I SMOOTH IT? PROBABLY YES!
    
    # Difference between two neighbouring indices
    diff_x_contour = np.diff(x_contour)
    diff_y_contour = np.diff(y_contour)
    
    # Get a list with the indices of duplicates
    diff_ind = []
    for ii in range(len(diff_x_contour)):
        if (diff_x_contour[ii]==0) and (diff_y_contour[ii]==0):
            diff_ind.append(ii)
    
    # Now remove the indices (start from the end so the indices don't shift)
    for ii in range(len(diff_ind)):
        index = diff_ind[::-1][ii]
        del x_contour[index]
        del y_contour[index]
    
    h_contour = np.zeros(len(x_contour))
    
    for ii in range(len(h_contour)):
        h_contour[ii] = h[int(y_contour[ii]), int(x_contour[ii])]
    
    # Get lat/lon along the contour
    
    # Choose whether you want your contour on the u or t grid.
    grid_sel = 't'
    
    if grid_sel == 'u':
        x_var = lon_u
        y_var = lat_u
    elif grid_sel == 'v':
        x_var = lon_v
        y_var = lat_v
    elif grid_sel == 't':
        x_var = lon_rho
        y_var = lat_rho
    
    lat_along_contour = np.zeros((len(x_contour)))
    lon_along_contour = np.zeros((len(x_contour)))
    
    for ii in range(len(h_contour)):
        lon_along_contour[ii] = x_var[int(y_contour[ii-1]),int(x_contour[ii-1])]
        lat_along_contour[ii] = y_var[int(y_contour[ii-1]),int(x_contour[ii-1])]
    
    # Repeat the leftmost point at the end of the array.
    # (Required for masking contour above and below)
    
    lat_along_contour = np.append(lat_along_contour, lat_along_contour[0])
    lon_along_contour = np.append(lon_along_contour, lon_along_contour[0])
    
    # Number of grid points on the contour
    num_points = len(lat_along_contour)
    
    # Now we number the points along the contour
    contour_mask_numbered = np.zeros_like(lon_along_contour)
    
    for ii in range(num_points-1):
        contour_mask_numbered[ii] = ii
    
    contour_mask = h*0
    
    for ii in range(num_points-1):
        contour_mask[int(y_contour[ii]), int(x_contour[ii])] = contour_mask_numbered[ii]+1
    mask_value = -1000
    contour_mask_numbered = contour_mask
    
    # fill in points to south of contour:
    contour_masked_above = np.copy(contour_mask_numbered)
    contour_masked_above[-1, 0] = mask_value
    
    #Create mask
    #Now we create a mask below contour sothat the direction of the contour can be determined
    
    #Remark on computational inefficiency:
    #Note that creating masks with nested for loops is very inefficient. We should probably use boolean masks (just compare the entire array with mask_value), and DataArray.shift() or DataArray.roll() from each of the directions to generate the masks without using loops.
    #See discussion in: https://github.com/COSIMA/cosima-recipes/issues/179
    
    print(contour_masked_above.shape, contour_mask_numbered.shape)
    print(contour_masked_above[-20:-1, 0])
    
    # fill in points to south of contour:
    contour_masked_above = np.copy(contour_mask_numbered)
    contour_masked_above[-1, 0] = mask_value
    
    # from top left:
    for ii in range(len(contour_mask[0,:])-1): #x: len(x-axis) - 1
        for jj in range(len(contour_mask[:,0]))[::-1][:-1]: #y: len(y-axis)[from end to start, inverse order][from first to (end-1)]
            if contour_masked_above[jj, ii] == mask_value: # if north of contour line
                if contour_masked_above[jj-1, ii] == 0: # if previous cell in Y-dir is zero (= south of contour line)
                    contour_masked_above[jj-1, ii] = mask_value # make it -1000
                if contour_masked_above[jj, ii+1] == 0: # if next cell in X-dir is zero
                    contour_masked_above[jj, ii+1] = mask_value # make it -1000
    
    #from top right:
    for ii in range(len(contour_mask[0,:]))[::-1][:-1]:
        for jj in range(len(contour_mask[:,0]))[::-1][:-1]:
            if contour_masked_above[jj, ii] == mask_value:
                if contour_masked_above[jj-1, ii] == 0: # if previous cell in Y-dir is zero
                    contour_masked_above[jj-1, ii] = mask_value
                if contour_masked_above[jj, ii-1] == 0: # if previous cell in X-dir is zero
                    contour_masked_above[jj, ii-1] = mask_value
    
    # from bottom right:
    for ii in range(len(contour_mask[0,:]))[::-1][:-1]:
        for jj in range(len(contour_mask[:,0])-1):
            if contour_masked_above[jj, ii] == mask_value:
                if contour_masked_above[jj+1, ii] == 0: # if next cell in Y-dir is zero
                    contour_masked_above[jj+1, ii] = mask_value
                if contour_masked_above[jj, ii-1] == 0: # if previous cell in X-dir is zero
                    contour_masked_above[jj, ii-1] = mask_value
    
    #from bottom left:
    for ii in range(len(contour_mask[0,:])-1):
        for jj in range(len(contour_mask[:,0])-1):
            if contour_masked_above[jj, ii] == mask_value:
                if contour_masked_above[jj+1, ii] == 0: # if next cell in Y-dir is zero
                    contour_masked_above[jj+1, ii] = mask_value
                if contour_masked_above[jj, ii+1] == 0: # if next cell in X-dir is zero
                    contour_masked_above[jj, ii+1] = mask_value
    
    mask_shelf2 = ma.masked_where(contour_masked_above == -1000, np.ones(h.shape))
    
    # Direction of cross-contour transport
    
    mask_x_transport = np.zeros_like(lon_u)
    mask_y_transport = np.zeros_like(lon_v)
    
    mask_x_transport_numbered = np.zeros_like(lon_u)
    mask_y_transport_numbered = np.zeros_like(lon_v)
    
    shape = contour_masked_above.shape
    
    new_number_count = 1
    for mask_loc in range(1, int(np.max(contour_mask_numbered))+1):
        #if mask_loc%100 == 0:
        #    print('mask for x/y transport at point '+str(mask_loc))
        index_i = np.where(contour_mask_numbered==mask_loc)[1]
        index_j = np.where(contour_mask_numbered==mask_loc)[0]
        # if point above is towards Antarctica and point below is away from Antarctica:
        # take transport grid point to north of t grid:
        if (contour_masked_above[index_j+1, index_i]==0) and (contour_masked_above[index_j-1, index_i]!=0):
            mask_y_transport[index_j, index_i] = -1
            # important to do
            mask_y_transport_numbered[index_j, index_i] = new_number_count
            new_number_count += 1
        # if point below is towards Antarctica and point above is away from Antarctica:
        # take transport grid point to south of t grid:
        elif (contour_masked_above[index_j-1, index_i]==0) and (contour_masked_above[index_j+1, index_i]!=0):
            mask_y_transport[index_j-1, index_i] = 1
            mask_y_transport_numbered[index_j-1, index_i] = new_number_count
            new_number_count += 1
        # if point below and point above are BOTH towards Antarctica:
        # take transport grid point to south of t grid:
        elif (contour_masked_above[index_j-1, index_i]==0) and (contour_masked_above[index_j+1, index_i]==0):
            mask_y_transport[index_j-1, index_i] = 1
            mask_y_transport[index_j, index_i] = -1
            mask_y_transport_numbered[index_j-1, index_i] = new_number_count
            mask_y_transport_numbered[index_j, index_i] = new_number_count+1
            new_number_count += 2
        # if point to right is towards Antarctica and point to left is away from Antarctica:
        # take transport grid point on right of t grid:
        if (contour_masked_above[index_j, index_i+1]==0) and (contour_masked_above[index_j, index_i-1]!=0):
            mask_x_transport[index_j, index_i] = -1
            mask_x_transport_numbered[index_j, index_i] = new_number_count
            new_number_count += 1
        # if point to left is towards Antarctica and point to right is away from Antarctica:
        # take transport grid point on left of t grid:
        elif (contour_masked_above[index_j, index_i-1]==0) and (contour_masked_above[index_j, index_i+1]!=0):
            mask_x_transport[index_j, index_i-1] = 1
            mask_x_transport_numbered[index_j, index_i-1] = new_number_count
            new_number_count += 1
        # if point to left and right BOTH toward Antarctica
        elif (contour_masked_above[index_j, index_i-1]==0) and (contour_masked_above[index_j, index_i+1]==0):
            mask_x_transport[index_j, index_i-1] = 1
            mask_x_transport[index_j, index_i] = -1
            mask_x_transport_numbered[index_j, index_i-1] = new_number_count
            mask_x_transport_numbered[index_j, index_i] = new_number_count+1
            new_number_count += 2
    
    # We now have the coordinates of the contours, and whether the x or y transport is needed to calculate cross-contour transport.
    
    print(pm.shape, pn.shape)
    
    # Convert contour masks to data arrays, so we can multiply them later.
    # We need to ensure the lat lon coordinates correspond to the actual data location:
    #       The y masks are used for ty_trans, so like vhrho this should have dimensions (yu_ocean, xt_ocean).
    #       The x masks are used for tx_trans, so like uhrho this should have dimensions (yt_ocean, xu_ocean).
    #       However the actual name will always be simply y_ocean/x_ocean irrespective of the variable
    #       to make concatenation of transports in both direction and sorting possible.
    coordinates=dict(one=lon_rho, two=lat_rho)
    coordinatesU=dict(one=lon_u, two=lat_u)
    coordinatesV=dict(one=lon_v, two=lat_v)

    mask_x_transport = xr.DataArray(mask_x_transport, coords = coordinatesU, dims = ['eta_u', 'xi_u'])
    mask_y_transport = xr.DataArray(mask_y_transport, coords = coordinatesV, dims = ['eta_v', 'xi_v'])
    mask_x_transport_numbered = xr.DataArray(mask_x_transport_numbered, coords = coordinatesU, dims = ['eta_u', 'xi_u'])
    mask_y_transport_numbered = xr.DataArray(mask_y_transport_numbered, coords = coordinatesV, dims = ['eta_v', 'xi_v'])

    # rename dimensions as simply eta/xi
    mask_x_transport = mask_x_transport.rename({'eta_u': 'eta','xi_u': 'xi'})
    mask_y_transport = mask_y_transport.rename({'eta_v': 'eta','xi_v': 'xi'})
    mask_x_transport_numbered = mask_x_transport_numbered.rename({'eta_u': 'eta','xi_u': 'xi'})
    mask_y_transport_numbered = mask_y_transport_numbered.rename({'eta_v': 'eta','xi_v': 'xi'})

    # Create the contour order data-array. Note that in this procedure the x-grid counts have x-grid
    #   dimensions and the y-grid counts have y-grid dimensions, but these are implicit, the dimension
    #   *names* are kept general across the counts, the generic y_ocean, x_ocean, so that concatening works
    #   but we dont double up with numerous counts for one lat/lon point.
    
    # stack contour data into 1d:
    mask_x_numbered_1d = mask_x_transport_numbered.stack(contour_index = ['eta', 'xi'])
    mask_x_numbered_1d = mask_x_numbered_1d.where(mask_x_numbered_1d > 0, drop = True)
    
    mask_y_numbered_1d = mask_y_transport_numbered.stack(contour_index = ['eta', 'xi'])
    mask_y_numbered_1d = mask_y_numbered_1d.where(mask_y_numbered_1d > 0, drop = True)
    
    contour_ordering = xr.concat((mask_x_numbered_1d, mask_y_numbered_1d), dim = 'contour_index', data_vars="all")
    contour_ordering = contour_ordering.sortby(contour_ordering)
    contour_index_array = np.arange(1, len(contour_ordering)+1)
    
    # using xr.open_mfdataset
    
    vars2drop = ["ubar","vbar","w","Hsbl","Hbbl","swrad"]
    
    ds = xr.open_mfdataset(paths="/g/data/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom10extend_shflim_S_0.25Q/output_20yr_diag_daily/ocean_avg_00*.nc" , chunks={'eta_rho': '200MB'}, parallel=True, drop_variables=vars2drop, decode_times=False)
    
    #- preserving 5-days avgs
    temp = ds.temp
    #shflux = ds.variables["shflux"]
    #ssflux = ds.variables["ssflux"]
    #m = ds.variables["m"]
    HvomT = ds.variables["Hvom_temp"]       ## !!! Huon_temp/Hvom_temp were not saved in the original run
    HuonT = ds.variables["Huon_temp"]       ## now it's running here: /scratch/gi0/fbd581/waom4extend_shflim_S_0.25Q/output_yr10_diag
    Hvom = ds.variables["Hvom"]
    Huon = ds.variables["Huon"]
    
    print("Vtransform=2")
    hwater = ds.h - abs(ds.zice)
    Zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * hwater) / (ds.hc + hwater)
    z_rho = ds.zeta + (ds.zeta + hwater) * Zo_rho - abs(ds.zice)
    Zo_w = (ds.hc * ds.s_w + ds.Cs_w * hwater) / (ds.hc + hwater)
    z_w = ds.zeta + (ds.zeta + hwater) * Zo_w - abs(ds.zice)

    ds.close()
    
    # subtract Tf heat transport from abs heat transport in the original grid (like done in ACCESS-OM2)
    # rho0 = 1025 # kg. m-3
    # Cp = 3989.245 # J.kg-1.degC-1
    # use same values as in access-om2
    rho0=1035
    Cp=3992.1
    # Tf = -1.95 # degC
    Tf =  -3.534879684448242 # coldest temp along 1500m among all three WAOM expts (10km, 4km, 4km-notide)

## ======================== section below if for the Advective horizontal heat transport =========================
    # --- calculate horizontal diffusion term:
    #  Dx = d/dx [(Ni * Hz * pm)/pn * dT/dx]
    #  Dy = d/dy [(Ni * Hz * pn)/pm * dT/dy]
    #

    # parameters: pm/pn already load
    ni_10km = 50 #m2/s
    ni_4km = 20 #m2/s
    # ni_2km = 10 #m2/s

    # 1st correct temp with Tf and do 1st derivative
    temp_Tf = temp - Tf
    dTdx = temp_Tf.diff('xi_rho',label='lower').rename({'eta_rho': 'eta','xi_rho': 'xi'})
    dTdy = temp_Tf.diff('eta_rho',label='lower').rename({'eta_rho': 'eta','xi_rho': 'xi'})
    print('HDiff calc; finish 1st derivatives dTdx,dTdy shapes: ',dTdx.shape, dTdy.shape)

    # compute layer thickness (on rho-grid):
    Z_w = z_w.transpose('ocean_time','s_w','eta_rho','xi_rho')
    Dz = Z_w.diff('s_w') # thickness w/ 31 layers
    Dz.rename({'s_w': 's_rho'})

    # compute m/n, n/m as from roms code (ROMS/Utility/metrics.F)

    # !-----------------------------------------------------------------------
    # !  Compute m/n, 1/m, and 1/n at horizontal U-points.
    # !-----------------------------------------------------------------------
    # !
    #       DO j=JstrT,JendT
    #         DO i=IstrP,IendT
    #           pmon_u(i,j)=(pm(i-1,j)+pm(i,j))/(pn(i-1,j)+pn(i,j))
    #           pnom_u(i,j)=(pn(i-1,j)+pn(i,j))/(pm(i-1,j)+pm(i,j))
    #           om_u(i,j)=2.0_r8/(pm(i-1,j)+pm(i,j))
    #           on_u(i,j)=2.0_r8/(pn(i-1,j)+pn(i,j))
    #         END DO
    #       END DO

    # attempt to do it with xarray:
    # - follows the order of the eqn in https://www.myroms.org/wiki/Horizontal_Mixing#Horizontal_Diffusion)
    # (i.e., derivativates are done in order)

    pm_u = pm.rolling({'xi_rho': 2}).mean() # first column is nan (so actual dimension is 629x560
    pn_v = pn.rolling({'eta_rho': 2}).mean() # first row is nan (so actual dimension is 630x559
    pmon_u = pm_u*pn_v
    pnom_v = pn_v*pm_u

    pmon_u=pmon_u.rename({'eta_rho': 'eta','xi_rho': 'xi'})
    pnom_v=pnom_v.rename({'eta_rho': 'eta','xi_rho': 'xi'})

    Dz_u = Dz.rolling({'xi_rho': 2}).mean().rename({'eta_rho': 'eta','xi_rho': 'xi','s_w':'s_rho'}) # first row is nan (so actual dimension is 630x559
    Dz_v = Dz.rolling({'eta_rho': 2}).mean().rename({'eta_rho': 'eta','xi_rho': 'xi','s_w':'s_rho'}) # first row is nan (so actual dimension is 630x559

    # print(pmon_u)
    pmon_u = pmon_u.drop_sel({'xi': pmon_u['xi'][0]}) # drop first line
    pnom_v = pnom_v.drop_sel({'eta': pnom_v['eta'][0]})
    Dz_u = Dz_u.drop_sel({'xi': Dz_u['xi'][0]})
    Dz_v = Dz_v.drop_sel({'eta': Dz_v['eta'][0]})

    # # replace s_rho in one variable to match:
    Dz_u['s_rho'] = dTdx.s_rho
    Dz_v['s_rho'] = dTdy.s_rho

    # calculate Hz*dT/dx & Hz*dT/dy
    Hzmon_dTdx = np.zeros(dTdx.shape)
    Hznom_dTdy = np.zeros(dTdy.shape)
    for tt in range(365): # 1 day all z-levels take ~10-11min
        for kk in range(31):
            Hzmon_dTdx[tt,kk,:] = dTdx.isel(ocean_time=tt, s_rho=kk)*Dz_u.isel(ocean_time=tt, s_rho=kk)*pmon_u*ni_10km
            Hznom_dTdy[tt,kk,:] = dTdy.isel(ocean_time=tt, s_rho=kk)*Dz_v.isel(ocean_time=tt, s_rho=kk)*pnom_v*ni_10km

    # convert Hzmon_dTdx/Hznom_dTdy back to DataArray:
    new_coords = dTdx.coords
    Hzmon_dTdx_da = xr.DataArray(Hzmon_dTdx, dims=dTdx.dims, coords=new_coords)
    new_coords = dTdy.coords
    Hznom_dTdy_da = xr.DataArray(Hznom_dTdy, dims=dTdx.dims, coords=new_coords)
    ## now calculate the 2nd derivative:
    dTdx2 = Hzmon_dTdx_da.diff('xi',label='lower')
    dTdy2 = Hznom_dTdy_da.diff('eta',label='lower')

    # apply mask_x/y_transport:
    #dTdx2_m=(dTdx2*mask_x_transport.isel(xi=slice(1,629))) # here I'm throwing out the 1st and last Xi (x) indexes
    #dTdy2_m=(dTdy2*mask_y_transport.isel(eta=slice(1,559))) # same but for Eta (y) indexes
    dTdx2_m=(dTdx2*mask_x_transport.isel(xi=slice(0,628))) # here I'm throwing out the 1st and last Xi (x) indexes
    dTdy2_m=(dTdy2*mask_y_transport.isel(eta=slice(0,558))) # same but for Eta (y) indexes

    # units: until here, the units are in: m.degC/sec
    # - apply rho0*Cp to convert to W/m2:
    hdiff_HTu = dTdx2_m*rho0*Cp
    hdiff_HTv = dTdy2_m*rho0*Cp

    
## ======================= section below if for the Advective horizontal heat transport =========================
#    # 1) multiply rho0*Cp:
#    HTv = HvomT*rho0*Cp
#    HTu = HuonT*rho0*Cp
#    
#    # 2) calculate mean freezing point heat transport:
#    #VTv_avg = Hvom.mean('ocean_time')
#    #VTu_avg = Huon.mean('ocean_time')
#    #HTf_v_avg = VTv_avg*Tf*rho0*Cp
#    #HTf_u_avg = VTu_avg*Tf*rho0*Cp
#    HTf_v = Hvom*Tf*rho0*Cp
#    HTf_u = Huon*Tf*rho0*Cp
#
#    # 3) subtract Tf HT from abs HT:
#    HTv = HTv - HTf_v
#    HTu = HTu - HTf_u
#    
#    # 4) rename dimensions before multiplying for mask_x/y_transport:
#    # convert to DataArray first:
#    HTv = xr.DataArray(HTv)
#    HTu = xr.DataArray(HTu)
#    HTv = HTv.rename({'eta_v': 'eta','xi_v': 'xi'})
#    HTu = HTu.rename({'eta_u': 'eta','xi_u': 'xi'})
#
#    # 5) multiply mask_x/y_transport:
#    HTvm = HTv*mask_y_transport
#    HTum = HTu*mask_x_transport
#
#    HTv.load()
#    HTu.load()
#    # ===========================================================================================================

    # convert temp to DataArray to extract values along contour: Not necessary in v4!
    months=np.arange(0,365)*(1/30.41667)
    
    # Convert heat transport to data arrays:
    coordinates3Du = dict(ocean_time=months, s_rho=(['s_rho'], np.arange(0,31)),
                        eta_u=(['eta_u'], np.arange(0,560)), xi_u=(['xi_u'], np.arange(0,629)))
    coordinates3Dv = dict(ocean_time=months, s_rho=(['s_rho'], np.arange(0,31)),
                        eta_v=(['eta_v'], np.arange(0,559)), xi_v=(['xi_v'], np.arange(0,630)))
    
    # - handling x/y transports (Hvom, Huon [m3.s-1]) to calculate heat transport
    #HuonT_xr = xr.DataArray(HTvm, coords = coordinates3Du, dims = ['ocean_time','s_rho','eta_u', 'xi_u'])
    #HvomT_xr = xr.DataArray(HTum, coords = coordinates3Dv, dims = ['ocean_time','s_rho','eta_v', 'xi_v'])
    
    # rename dimensions as simply eta/xi
    #HuonT_xr = HuonT_xr.rename({'eta_u': 'eta','xi_u': 'xi'})
    #HvomT_xr = HvomT_xr.rename({'eta_v': 'eta','xi_v': 'xi'})

    # define a function to extract any 4D var along the contour line
    

    def extract_transp_across_contour(var_x, var_y):   # var:4D [time,eta_rho,xi_rho]
        tlen = len(var_x[:,0,0,0])
        zlen = len(var_x[0,:,0,0])
        transp_across_contour = np.empty((tlen,zlen,len(contour_ordering)))
        for tt in range(0,tlen):
            for zz in range(0,zlen): # loop through z-levels
                var_x_tmp = var_x[zz,:,:]
                var_y_tmp = var_y[zz,:,:]

                # stack transports into 1d and drop any points not on contour:
                x_var_1d_tmp = var_x_tmp.stack(contour_index = ['eta', 'xi'])
                x_var_1d_tmp = x_var_1d_tmp.where(mask_x_numbered_1d>0, drop = True)
                y_var_1d_tmp = var_y_tmp.stack(contour_index = ['eta', 'xi'])
                y_var_1d_tmp = y_var_1d_tmp.where(mask_y_numbered_1d>0, drop = True)

                # combine all points on contour:
                transp_across_contour_tmp = xr.concat((x_var_1d_tmp, y_var_1d_tmp), dim = 'contour_index')
                transp_across_contour_tmp = transp_across_contour_tmp.reset_index('contour_index') # added by fabio, otherwise it crashes due to duplicated indices
                transp_across_contour_tmp = transp_across_contour_tmp.sortby(contour_ordering)
                transp_across_contour_tmp.coords['contour_index'] = contour_index_array
                transp_across_contour_tmp = transp_across_contour_tmp.load()

                transp_across_contour[zz,:] = transp_across_contour_tmp
                del transp_across_contour_tmp

        return transp_across_contour

    # extract variables:
    # 1. vol transp
    diff_heat_trans_across_contour = extract_transp_across_contour(hdiff_HTu, hdiff_HTv)
    
    # save to netcdf file:
    coordinatesC=dict(ocean_time=months, s_rho=(['s_rho'], np.arange(0,31)),
                        contour_index_array=(['contour_index_array'], np.arange(0,len(contour_index_array))))
    
    diff_heat_trans_across_contour_xr = xr.DataArray(diff_heat_trans_across_contour, coords = coordinatesC, dims = ['ocean_time','s_rho','contour_index_array'])
    files_path = '/g/data/hh5/tmp/access-om/fbd581/ROMS/postprocessing/cross_contour_tmp/'
    diff_heat_trans_across_contour_xr.to_netcdf(files_path + 'WAOM10_hdiff_heat_trans_1500m_daily_v4', mode='w', format="NETCDF4")
    # -> PS: V4 already accounts for (-Tf HT), contrary to V3.
