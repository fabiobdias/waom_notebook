#!/usr/bin/env python3
# Ocean Heat Budget Analyses in the Antarctica continental shelf (WAOM)

# Fabio B Dias - 28 June 2023
# Description:
#     this script obtain and save the CF isobath contour variables, which is used for the 
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
    di = xr.open_dataset('/g/data/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom4extend_shflim_S_0.25Q/output_yr10_diag_daily/ocean_avg_0001.nc')
    ice_draft = di.variables["zice"]
    
    mask_zice = ma.masked_where(ice_draft < 0, np.ones(ice_draft.shape))
    mask_outice = ma.masked_where(ice_draft >= 0, np.ones(ice_draft.shape))
    mask_zice_1000 = ma.masked_where(ice_draft < -1000, np.ones(ice_draft.shape))
    di.close()
    
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
    
    # ds.coords['lat_rho']=lat_rho.transpose() # put lat_rho into ds dataset
    # ds.coords['lon_rho']=lon_rho.transpose() # put lon_rho into ds dataset
    
    area=np.divide(1,pm*pn)
    
    ## creating the contour, such as a isobath, and extracting the coordinates using matplotlib's Path class
    # based on https://github.com/COSIMA/cosima-recipes/blob/master/DocumentedExamples/Cross-contour_transport.ipynb
    zice = dg.zice.load()
    zice = zice*mask_zice_1000
    
    # Fill in land with zeros:
    zice = zice.fillna(0)
    
    contour_depth = -.01
    
    ## Choose whether you want your contour on the u or t grid.
    grid_sel = 't'
    if grid_sel == 'u':
        x_var = lon_u
        y_var = lat_u
    elif grid_sel == 't':
        x_var = lon_rho
        y_var = lat_rho

# load contour generated from WAOM4_Extract_CF_contour.ipynb
    fig_path = '/g/data/hh5/tmp/access-om/fbd581/ROMS/postprocessing/figs/Contour_isobath/'
    xcon_np=np.loadtxt("/g/data/hh5/tmp/access-om/fbd581/ROMS/postprocessing/figs/Contour_isobath/WAOM4_CF_x_contour.csv")
    x_contour = xcon_np.tolist()
    ycon_np=np.loadtxt("/g/data/hh5/tmp/access-om/fbd581/ROMS/postprocessing/figs/Contour_isobath/WAOM4_CF_y_contour.csv")
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
    
    # make halos: add 2 extra columns with the value of the last/first columns of the original
    shape = contour_masked_above.shape
    contour_masked_above_halo = np.zeros((shape[0], shape[1]+2))
    contour_masked_above_halo[:, 0] = contour_masked_above[:, -1]
    contour_masked_above_halo[:, 1:-1] = contour_masked_above
    contour_masked_above_halo[:, -1] = contour_masked_above[:, 0]
    
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
        # zonal indices increased by 1 due to halos
        # take transport grid point on right of t grid:
        if (contour_masked_above_halo[index_j, index_i+2]==0) and (contour_masked_above_halo[index_j, index_i]!=0):
            mask_x_transport[index_j, index_i] = -1
            mask_x_transport_numbered[index_j, index_i] = new_number_count
            new_number_count += 1
        # if point to left is towards Antarctica and point to right is away from Antarctica:
        # take transport grid point on left of t grid:
        elif (contour_masked_above_halo[index_j, index_i]==0) and (contour_masked_above_halo[index_j, index_i+2]!=0):
            mask_x_transport[index_j, index_i-1] = 1
            mask_x_transport_numbered[index_j, index_i-1] = new_number_count
            new_number_count += 1
        # if point to left and right BOTH toward Antarctica
        elif (contour_masked_above_halo[index_j, index_i]==0) and (contour_masked_above_halo[index_j, index_i+2]==0):
            mask_x_transport[index_j, index_i-1] = 1
            mask_x_transport[index_j, index_i] = -1
            mask_x_transport_numbered[index_j, index_i-1] = new_number_count
            mask_x_transport_numbered[index_j, index_i] = new_number_count+1
            new_number_count += 2
    
    # We now have the coordinates of the contours, and whether the x or y transport is needed to calculate cross-contour transport.
    
    print(pm.shape, pn.shape)
    
    # Now we need to interpolate the mask_x/y_transport to the corresponding lon/lat_u/v grids
    # so we can multiply by the U/V transport
    
    # re-grid from rho to u/v-grids
    rho_def = pyresample.geometry.SwathDefinition(lons=lon_rho,lats=lat_rho)
    u_def = pyresample.geometry.SwathDefinition(lons=lon_u,lats=lat_u)
    v_def = pyresample.geometry.SwathDefinition(lons=lon_v,lats=lat_v)
    
    wf = lambda r: 1/r
    
    mask_x_transport_Ugrd = pyresample.kd_tree.resample_custom(rho_def,mask_x_transport,u_def,\
                                             radius_of_influence=100000,neighbours=1,weight_funcs=wf)
    mask_y_transport_Vgrd = pyresample.kd_tree.resample_custom(rho_def,mask_y_transport,v_def,\
                                             radius_of_influence=100000,neighbours=1,weight_funcs=wf)
    mask_x_transport_numbered_Ugrd = pyresample.kd_tree.resample_custom(rho_def,mask_x_transport_numbered,u_def,\
                                             radius_of_influence=100000,neighbours=1,weight_funcs=wf)
    mask_y_transport_numbered_Vgrd = pyresample.kd_tree.resample_custom(rho_def,mask_y_transport_numbered,v_def,\
                                             radius_of_influence=100000,neighbours=1,weight_funcs=wf)
    
    # Convert contour masks to data arrays, so we can multiply them later.
    # We need to ensure the lat lon coordinates correspond to the actual data location:
    #       The y masks are used for ty_trans, so like vhrho this should have dimensions (yu_ocean, xt_ocean).
    #       The x masks are used for tx_trans, so like uhrho this should have dimensions (yt_ocean, xu_ocean).
    #       However the actual name will always be simply y_ocean/x_ocean irrespective of the variable
    #       to make concatenation of transports in both direction and sorting possible.
    coordinates=dict(one=lon_rho, two=lat_rho)
    coordinatesU=dict(one=lon_u, two=lat_u)
    coordinatesV=dict(one=lon_v, two=lat_v)
    
    
    mask_x_transport_Ugrd = xr.DataArray(mask_x_transport_Ugrd, coords = coordinatesU, dims = ['eta_u', 'xi_u'])
    mask_y_transport_Vgrd = xr.DataArray(mask_y_transport_Vgrd, coords = coordinatesV, dims = ['eta_v', 'xi_v'])
    mask_x_transport_numbered_Ugrd = xr.DataArray(mask_x_transport_numbered_Ugrd, coords = coordinatesU, dims = ['eta_u', 'xi_u'])
    mask_y_transport_numbered_Vgrd = xr.DataArray(mask_y_transport_numbered_Vgrd, coords = coordinatesV, dims = ['eta_v', 'xi_v'])
    
    # rename dimensions as simply eta/xi
    mask_x_transport_Ugrd = mask_x_transport_Ugrd.rename({'eta_u': 'eta','xi_u': 'xi'})
    mask_y_transport_Vgrd = mask_y_transport_Vgrd.rename({'eta_v': 'eta','xi_v': 'xi'})
    mask_x_transport_numbered_Ugrd = mask_x_transport_numbered_Ugrd.rename({'eta_u': 'eta','xi_u': 'xi'})
    mask_y_transport_numbered_Vgrd = mask_y_transport_numbered_Vgrd.rename({'eta_v': 'eta','xi_v': 'xi'})
    
    # Create the contour order data-array. Note that in this procedure the x-grid counts have x-grid
    #   dimensions and the y-grid counts have y-grid dimensions, but these are implicit, the dimension
    #   *names* are kept general across the counts, the generic y_ocean, x_ocean, so that concatening works
    #   but we dont double up with numerous counts for one lat/lon point.
    
    # stack contour data into 1d:
    mask_x_numbered_1d = mask_x_transport_numbered_Ugrd.stack(contour_index = ['eta', 'xi'])
    mask_x_numbered_1d = mask_x_numbered_1d.where(mask_x_numbered_1d > 0, drop = True)
    
    mask_y_numbered_1d = mask_y_transport_numbered_Vgrd.stack(contour_index = ['eta', 'xi'])
    mask_y_numbered_1d = mask_y_numbered_1d.where(mask_y_numbered_1d > 0, drop = True)
    
    contour_ordering = xr.concat((mask_x_numbered_1d, mask_y_numbered_1d), dim = 'contour_index', data_vars="all")
    contour_ordering = contour_ordering.sortby(contour_ordering)

    # get lat and lon along contour, useful for plotting later:
    lat_along_contour = contour_ordering.two
    lon_along_contour = contour_ordering.one

    contour_index_array = np.arange(1, len(contour_ordering)+1)

    # don't need the multi-index anymore, replace with contour count and save
    lat_along_contour.coords['contour_index'] = contour_index_array
    lon_along_contour.coords['contour_index'] = contour_index_array
    
    # using xr.open_mfdataset
    
    vars2drop = ["ubar","vbar","w","Hsbl","Hbbl","swrad"]

    ### Code to extract distance in between contour coordinates, using length of diagonal if there is a bend.
    #Loop through the contour, determining if diagonal is required or not, and save the distance along each segment. Then, cumulatively sum the distances along each segment to get the distance from the first point.
    #If there is a bend in the contour, then half the diagonal distance is added to each side to avoid artifically inflating the along-contour distance metric, according to this diagram:

    num_points = len(lat_along_contour)

    # if there is a bend in the contour, add the distance using the half-length of the diagonal
    # instead of the sum of 2 edges, to be more representative.
    distance_along_contour = np.zeros((num_points))

    x_indices = np.sort(mask_x_transport_numbered[mask_x_transport_numbered>0])
    y_indices = np.sort(mask_y_transport_numbered[mask_y_transport_numbered>0])

    skip = False
    # note dxu and dyt do not vary in x, so we can just take the first value (as long as there is no land there,
    # which for this latitude range there is not. If using a different latitude range, choose an x value that is
    # not a nan/land for the entire latitude range
    dxu = np.divide(1,pm)
    dyt = np.divide(1,pn)

    for count in range(1, num_points):
        if skip == True:
            skip = False
            continue
        if count in y_indices:
            if count + 1 in y_indices:
                # note dxu and dyt do not vary in x:
                jj = np.where(mask_y_transport_numbered==count)#[0]
                # print(count, jj)
                distance_along_contour[count-1] = (dxu[jj])#[0]
            else:
                jj0 = np.where(mask_y_transport_numbered==count)#[0]
                jj1 = np.where(mask_x_transport_numbered==count+1)#[0]
                half_diagonal_distance = 0.5 * np.sqrt((dxu[jj0])**2 + (dyt[jj1])**2)
                distance_along_contour[count-1] = half_diagonal_distance
                distance_along_contour[count] = half_diagonal_distance
                # skip to next count:
                skip = True

        # count in x_indices:
        else:
            if count + 1 in x_indices:
                jj = np.where(mask_x_transport_numbered==count)#[0]
                distance_along_contour[count-1] = (dyt[jj])#[0]
            else:
                jj0 = np.where(mask_x_transport_numbered==count)#[0]
                jj1 = np.where(mask_y_transport_numbered==count+1)#[0]
                half_diagonal_distance = 0.5 * np.sqrt((dyt[jj0])**2 + (dxu[jj1])**2)
                distance_along_contour[count-1] = half_diagonal_distance
                distance_along_contour[count] = half_diagonal_distance
                # skip to next count:
                skip = True

    # fix last value:
    if distance_along_contour[-1] == 0:
        count = count + 1
        if count in y_indices:
            jj = np.where(mask_y_transport_numbered==count)#[0]
            distance_along_contour[-1] = (dxu[jj])#[0]
        else:
            jj = np.where(mask_x_transport_numbered==count)#[0]
            distance_along_contour[-1] = (dyt[jj])#[0]

    # units are 10^3 km:
    distance_along_contour = np.cumsum(distance_along_contour) / 1e3 / 1e3


    # --- Select the indices for axis labels of specific longitudes,
    # so we can plot transport vs distance but have longitude labels instead of length

    distance_indices = np.zeros(9)

    for i in np.arange(100, len(lon_along_contour.one)):
        # print(i)
        if (distance_indices[1]==0):
            if (lon_along_contour.one[i]>-160 and lon_along_contour.one[i]<-159):
                distance_indices[1] = lon_along_contour.contour_index.values[i]
        if (distance_indices[2]==0):
            if (lon_along_contour.one[i]>-120 and lon_along_contour.one[i]<-119):
                distance_indices[2] = lon_along_contour.contour_index.values[i]
        if (distance_indices[3]==0):
            if (lon_along_contour.one[i]>-80  and lon_along_contour.one[i]<-79):
                distance_indices[3] = lon_along_contour.contour_index.values[i]
        if (distance_indices[4]==0):
            if (lon_along_contour.one[i]>-40 and lon_along_contour.one[i]<-39):
                distance_indices[4] = lon_along_contour.contour_index.values[i]
        if (distance_indices[5]==0):
            if (lon_along_contour.one[i]>0 and lon_along_contour.one[i]<1):
                distance_indices[5] = lon_along_contour.contour_index.values[i]
        if (distance_indices[6]==0):
            if (lon_along_contour.one[i]>40 and lon_along_contour.one[i]<41):
                distance_indices[6] = lon_along_contour.contour_index.values[i]
        if (distance_indices[7]==0):
            if (lon_along_contour.one[i]>80 and lon_along_contour.one[i]<81):
                distance_indices[7] = lon_along_contour.contour_index.values[i]
        if (distance_indices[8]==0):
            if (lon_along_contour.one[i]>120 and lon_along_contour.one[i]<121):
                distance_indices[8] = lon_along_contour.contour_index.values[i]

    #distance_indices[13] = len(lon_along_contour.contour_index.values)-1

    ## save coordinates along contour to use for plotting etc:
    files_path = '/g/data/hh5/tmp/access-om/fbd581/ROMS/postprocessing/cross_contour_tmp/'

    coordinatesD=dict(contour_index_array=(['contour_index_array'], np.arange(0,len(contour_index_array))))

    #lon_along_contour
    lon_along_contour_xr = xr.DataArray(lon_along_contour, coords = coordinatesD, dims = ['contour_index_array'])
    lon_along_contour_xr.to_netcdf(files_path + 'WAOM4_lon_along_CF_v2', mode='w', format="NETCDF4")

    #lat_along_contour
    lat_along_contour_xr = xr.DataArray(lat_along_contour, coords = coordinatesD, dims = ['contour_index_array'])
    lat_along_contour_xr.to_netcdf(files_path + 'WAOM4_lat_along_CF_v2', mode='w', format="NETCDF4")

    ## distance
    distance_along_contour_xr = xr.DataArray(distance_along_contour, coords = coordinatesD, dims = ['contour_index_array'])
    distance_along_contour_xr.to_netcdf(files_path + 'WAOM4_dist_along_CF_v2', mode='w', format="NETCDF4")

    coordinatesI=dict(lon_index=(['lon_index'], np.arange(0,9)))
    distance_indices_xr = xr.DataArray(distance_indices, coords = coordinatesI, dims = ['lon_index'])
    distance_indices_xr.to_netcdf(files_path + 'WAOM4_dist_indices_CF_v2', mode='w', format="NETCDF4")




