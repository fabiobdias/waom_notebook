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
    di = xr.open_dataset('/scratch/gi0/fbd581/waom10extend_testcode/output_21yr_diag_daily/ocean_avg_0001.nc')
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
        lat1 = lat_along_contour[ii]
        lat2 = lat_along_contour[ii+1]
        lon1 = lon_along_contour[ii]
        lon2 = lon_along_contour[ii+1]
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
    
    mask_x_transport = np.zeros_like(contour_mask_numbered)
    mask_y_transport = np.zeros_like(contour_mask_numbered)
    
    mask_y_transport_numbered = np.zeros_like(contour_mask_numbered)
    mask_x_transport_numbered = np.zeros_like(contour_mask_numbered)
    
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
    
    
    mask_x_transport = xr.DataArray(mask_x_transport, coords = coordinates, dims = ['eta_rho', 'xi_rho'])
    mask_y_transport = xr.DataArray(mask_y_transport, coords = coordinates, dims = ['eta_rho', 'xi_rho'])
    mask_x_transport_numbered = xr.DataArray(mask_x_transport_numbered, coords = coordinates, dims = ['eta_rho', 'xi_rho'])
    mask_y_transport_numbered = xr.DataArray(mask_y_transport_numbered, coords = coordinates, dims = ['eta_rho', 'xi_rho'])
    
    # rename dimensions as simply eta/xi
    mask_x_transport = mask_x_transport.rename({'eta_rho': 'eta','xi_rho': 'xi'})
    mask_y_transport = mask_y_transport.rename({'eta_rho': 'eta','xi_rho': 'xi'})
    mask_x_transport_numbered = mask_x_transport_numbered.rename({'eta_rho': 'eta','xi_rho': 'xi'})
    mask_y_transport_numbered = mask_y_transport_numbered.rename({'eta_rho': 'eta','xi_rho': 'xi'})
    
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
    
    ds = xr.open_mfdataset(paths="/scratch/gi0/fbd581/waom10extend_testcode/output_21yr_diag_daily/ocean_avg_0*.nc" , chunks={'eta_rho': '200MB'}, parallel=True, drop_variables=vars2drop, decode_times=False)
    
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
    
    # convert temp to DataArray to extract values along contour:
    months=np.arange(0,365)*(1/30.41667)
    
    coordinatesT=dict(ocean_time=months, s_rho=(['s_rho'], np.arange(0,31)),
                        eta_rho=(['eta_rho'], np.arange(0,560)), xi_rho=(['xi_rho'], np.arange(0,630)))
    temp_xr = xr.DataArray(temp, coords = coordinatesT, dims = ['ocean_time','s_rho','eta_rho', 'xi_rho'])
   
    # rename dimensions as simply eta/xi
    temp_xr = temp_xr.rename({'eta_rho': 'eta','xi_rho': 'xi'})

    # define a function to extract any 4D var along the contour line
    
    def extract_var_along_contour(var):   # var:4D [time,eta_rho,xi_rho]

        zlen = len(temp[0,:,0,0])
        tlen = len(temp[:,0,0,0])
        print(tlen,zlen)
        var_along_contour = np.empty((tlen,zlen,len(contour_ordering)))

        for tt in range(0,tlen): # loop through time
            for zz in range(0,zlen): # loop through z-levels
                var_tmp = var[tt,zz,:,:]
                
                # stack transports into 1d and drop any points not on contour:
                x_var_1d_tmp = var_tmp.stack(contour_index = ['eta', 'xi'])
                x_var_1d_tmp = x_var_1d_tmp.where(mask_x_numbered_1d>0, drop = True)
                y_var_1d_tmp = var_tmp.stack(contour_index = ['eta', 'xi'])
                y_var_1d_tmp = y_var_1d_tmp.where(mask_y_numbered_1d>0, drop = True)

                # combine all points on contour:
                var_along_contour_tmp = xr.concat((x_var_1d_tmp, y_var_1d_tmp), dim = 'contour_index')
                var_along_contour_tmp = var_along_contour_tmp.reset_index('contour_index') # added by fabio, otherwise it crashes due to duplicated indices
                var_along_contour_tmp = var_along_contour_tmp.sortby(contour_ordering)
                var_along_contour_tmp.coords['contour_index'] = contour_index_array
                var_along_contour_tmp = var_along_contour_tmp.load()

                print(tt, zz, var_along_contour_tmp.shape)
                var_along_contour[tt,zz,:] = var_along_contour_tmp
                del var_along_contour_tmp

        return var_along_contour

    # extract variables:
    # 1. temp
    temp_along_contour = extract_var_along_contour(temp_xr)

    # save to netcdf file:
    coordinatesC=dict(ocean_time=months, s_rho=(['s_rho'], np.arange(0,31)),
                        contour_index_array=(['contour_index_array'], np.arange(0,len(contour_index_array))))
    
    temp_along_contour_xr = xr.DataArray(temp_along_contour, coords = coordinatesC, dims = ['ocean_time','s_rho','contour_index_array'])
    files_path = '/g/data/hh5/tmp/access-om/fbd581/ROMS/postprocessing/cross_contour_tmp/'
    temp_along_contour_xr.to_netcdf(files_path + 'WAOM10_newcode_temp_1500m_daily_v3', mode='w', format="NETCDF4")
    
    
