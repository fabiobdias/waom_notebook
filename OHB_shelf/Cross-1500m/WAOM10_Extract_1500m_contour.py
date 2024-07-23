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
import cartopy.feature as cfeature
import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as dates
from matplotlib.ticker import NullFormatter
from matplotlib.colors import LinearSegmentedColormap   # for custom colormaps
import cmocean

from datetime import datetime, timedelta

from netCDF4 import Dataset
from netCDF4 import num2date, date2num
from matplotlib.ticker import NullFormatter
from matplotlib.colors import LinearSegmentedColormap   # for custom colormaps

import gsw

import pyresample

import scipy.io

from dask.distributed import Client
import logging
import warnings
warnings.filterwarnings('ignore')

client = Client(n_workers=56,threads_per_worker=1, memory_limit=0, silence_logs=logging.ERROR)
client

    %%time

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

    dg.close()

    # to re-load
    xcon_np=np.loadtxt("/g/data/hh5/tmp/access-om/fbd581/ROMS/postprocessing/figs/Contour_isobath/WAOM10_1500m_x_contour.csv")
    x_contour = xcon_np.tolist()
    ycon_np=np.loadtxt("/g/data/hh5/tmp/access-om/fbd581/ROMS/postprocessing/figs/Contour_isobath/WAOM10_1500m_y_contour.csv")
    y_contour = ycon_np.tolist()

    # Difference between two neighbouring indices
    diff_x_contour = np.diff(x_contour)
    diff_y_contour = np.diff(y_contour)

    # Get a list with the indices of duplicates
    diff_ind = []
    for ii in range(len(diff_x_contour)):
        if (diff_x_contour[ii]==0) and (diff_y_contour[ii]==0):
            diff_ind.append(ii)
    # print(diff_ind)
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
        # lat1 = lat_along_contour[ii]
        # lat2 = lat_along_contour[ii+1]
        # lon1 = lon_along_contour[ii]
        # lon2 = lon_along_contour[ii+1]
        contour_mask_numbered[ii] = ii

    contour_mask = h*0

    #Create mask
    #Now we create a mask below contour so that the direction of the contour can be determined

    #Remark on computational inefficiency:
    #Note that creating masks with nested for loops is very inefficient. We should probably use boolean masks (just compare the entire array with mask_value), and DataArray.shift() or DataArray.roll() from each of the directions to generate the masks without using loops.
    #See discussion in: https://github.com/COSIMA/cosima-recipes/issues/179

    for ii in range(num_points-1):
        contour_mask[int(y_contour[ii]), int(x_contour[ii])] = contour_mask_numbered[ii]+1
    mask_value = -1000
    contour_mask_numbered = contour_mask

    # fill in points to south of contour:
    contour_masked_above = np.copy(contour_mask_numbered)
    contour_masked_above[-1, 0] = mask_value # this makes pt x=0, y=-1 equal to -1000
    # then this will be the first point to enter in the fiirst IF-case in the next cell;

    #from top left:
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

        # 1A) if point above is towards Antarctica and point below is away from Antarctica:
        # take transport grid point to north of t grid: (same as t-grid index)
        if (contour_masked_above[index_j+1, index_i]==0) and (contour_masked_above[index_j-1, index_i]!=0):
            mask_y_transport[index_j, index_i] = -1
            # important to do
            mask_y_transport_numbered[index_j, index_i] = new_number_count
            new_number_count += 1

        # 1B) if point below is towards Antarctica and point above is away from Antarctica:
        # take transport grid point to south of t grid: (j=-1 from t-grid index)
        elif (contour_masked_above[index_j-1, index_i]==0) and (contour_masked_above[index_j+1, index_i]!=0):
            mask_y_transport[index_j-1, index_i] = 1
            mask_y_transport_numbered[index_j-1, index_i] = new_number_count
            new_number_count += 1

        # 1C) if point below and point above are BOTH towards Antarctica:
        # take transport grid point to south of t grid: (Embayment case; where t-grid is surronded by both points towards Antarctica)
        elif (contour_masked_above[index_j-1, index_i]==0) and (contour_masked_above[index_j+1, index_i]==0):
            mask_y_transport[index_j-1, index_i] = 1
            mask_y_transport[index_j, index_i] = -1
            mask_y_transport_numbered[index_j-1, index_i] = new_number_count
            mask_y_transport_numbered[index_j, index_i] = new_number_count+1
            new_number_count += 2

        # 2A) if point to right is towards Antarctica and point to left is away from Antarctica:
        # take transport grid point on right of t grid: (same as t-grid index)
        if (contour_masked_above[index_j, index_i+1]==0) and (contour_masked_above[index_j, index_i-1]!=0):
            mask_x_transport[index_j, index_i] = -1
            mask_x_transport_numbered[index_j, index_i] = new_number_count
            new_number_count += 1

        # 2B) if point to left is towards Antarctica and point to right is away from Antarctica:
        # take transport grid point on left of t grid: (i=-1 from t-grid index)
        elif (contour_masked_above[index_j, index_i-1]==0) and (contour_masked_above[index_j, index_i+1]!=0):
            mask_x_transport[index_j, index_i-1] = 1
            mask_x_transport_numbered[index_j, index_i-1] = new_number_count
            new_number_count += 1

        # 2C) if point to left and right BOTH toward Antarctica
        elif (contour_masked_above[index_j, index_i-1]==0) and (contour_masked_above[index_j, index_i+1]==0):
            mask_x_transport[index_j, index_i-1] = 1
            mask_x_transport[index_j, index_i] = -1
            mask_x_transport_numbered[index_j, index_i-1] = new_number_count
            mask_x_transport_numbered[index_j, index_i] = new_number_count+1
            new_number_count += 2

    dg = xr.open_dataset("/g/data/hh5/tmp/access-om/fbd581/ROMS/waom10_frc/waom10extend_grd.nc")

    # save mask_x/y_transport for inspection;
    coordinates=dict(lon=(['lon_rho','lat_rho'],dg.lon_rho.data), lat=(['lon_rho','lat_rho'], dg.lat_rho.data))
    coordinatesU=dict(lon=(['lon_u','lat_u'],dg.lon_u.data), lat=(['lon_u','lat_u'], dg.lat_u.data))
    coordinatesV=dict(lon=(['lon_v','lat_v'],dg.lon_v.data), lat=(['lon_v','lat_v'], dg.lat_v.data))

    mask_x_transport_xr = xr.DataArray(mask_x_transport, coords = coordinatesU, dims = ['lon_u','lat_u'])
    mask_y_transport_xr = xr.DataArray(mask_y_transport, coords = coordinatesV, dims = ['lon_v','lat_v'])

    # rename vars:
    mask_x_transport_xr.name = 'mask_x_transport'
    mask_y_transport_xr.name = 'mask_y_transport'

    files_path = '/g/data/hh5/tmp/access-om/fbd581/ROMS/postprocessing/cross-1500m_transp/'

    mask_x_transport_xr.to_netcdf(files_path + 'mask_x_transport.nc', mode='w', format="NETCDF4")
    mask_y_transport_xr.to_netcdf(files_path + 'mask_y_transport.nc', mode='w', format="NETCDF4")

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

    # get lat and lon along contour, useful for plotting later:
    lat_along_contour = contour_ordering.two
    lon_along_contour = contour_ordering.one

    contour_index_array = np.arange(1, len(contour_ordering)+1)

    # don't need the multi-index anymore, replace with contour count and save
    lat_along_contour.coords['contour_index'] = contour_index_array
    lon_along_contour.coords['contour_index'] = contour_index_array

    # using xr.open_mfdataset

    vars2drop = ["ubar","vbar","w","Hsbl","Hbbl","swrad"]

    ds = xr.open_mfdataset(paths="/g/data/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom10extend_shflim_S_0.25Q/output_20yr_diag_daily/ocean_avg_00*.nc" , chunks={'eta_rho': '200MB'}, parallel=True, drop_variables=vars2drop, decode_times=False)

    #- preserving daily avgs
    temp = ds.temp
    # salt = ds.variables["salt"]
    # shflux = ds.variables["shflux"]
    # ssflux = ds.variables["ssflux"]
    # m = ds.variables["m"]
    HvomT = ds.variables["Hvom_temp"]       ## !!! Huon_temp/Hvom_temp were not saved in the original run
    HuonT = ds.variables["Huon_temp"]       ## now it's running here: /scratch/gi0/fbd581/waom4extend_shflim_S_0.25Q/output_yr10_diag
    Hvom = ds.variables["Hvom"]
    Huon = ds.variables["Huon"]
    ssh = ds.variables["zeta"]

    ds = ds.set_coords(['Cs_r', 'Cs_w', 'hc', 'h', 'Vtransform'])

    print("Vtransform=2")
    hwater = ds.h - abs(ds.zice)
    Zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * hwater) / (ds.hc + hwater)
    z_rho = ds.zeta + (ds.zeta + hwater) * Zo_rho - abs(ds.zice)
    Zo_w = (ds.hc * ds.s_w + ds.Cs_w * hwater) / (ds.hc + hwater)
    z_w = ds.zeta + (ds.zeta + hwater) * Zo_w - abs(ds.zice)


    ds.close()

    # compute layer thickness (on rho-grid):
    Z_w = z_w.transpose('ocean_time','s_w','eta_rho','xi_rho')

    print(Z_w.shape)
    Dz = Z_w.diff('s_w') # thickness w/ 31 layers

    Dz = Dz.rename({'s_w': 's_rho'})

    ds = xr.open_mfdataset(paths="/g/data/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom10extend_shflim_S_0.25Q/output_20yr_diag_daily/ocean_dia_00*.nc" , chunks={'eta_rho': '200MB'}, parallel=True, drop_variables=vars2drop, decode_times=False)

    #- preserving daily avgs
    temp_hadv = ds.temp_hadv

    ds.close()

    ## testing on the grid-point budget closure between heat and temperature diagnostics:

    heat_hadv = temp_hadv*Dz

    xlen=560
    ylen=630

    HADV=np.empty((ylen,xlen))
    HUON_VOM_T=np.empty((ylen,xlen))
    for I in range(xlen-1):
        for J in range(ylen-1):

            tadv = temp_hadv[0,:,J,I].values # degC s-1
            t = temp[0,:,J,I].values # degC
            dz = Dz[0,:,J,I].values # m
            m=1/pm[J,I].values #m
            n=1/pn[J,I].values #m
            HADV[J,I] = np.sum(dz*tadv)*m*n # m3 degC s-1


            if I==0:
                huonT = HuonT[0,:,J,I].values # m3 degC s-1
                hvomT = HvomT[0,:,J,I].values
                hvomT_ = HvomT[0,:,J-1,I].values
                HUON_VOM_T[J,I] = np.sum((hvomT-hvomT_)+(huonT))
            elif J==0:
                huonT = HuonT[0,:,J,I].values # m3 degC s-1
                hvomT = HvomT[0,:,J,I].values
                huonT_ = HuonT[0,:,J,I-1].values
                HUON_VOM_T[J,I] = np.sum((hvomT)+(huonT-huonT_))
            elif I==0 and J==0:
                huonT = HuonT[0,:,J,I].values # m3 degC s-1
                hvomT = HvomT[0,:,J,I].values
                HUON_VOM_T[J,I] = np.sum((hvomT)+(huonT))
            elif I==(xlen-2):
                huonT = HuonT[0,:,J,I].values # m3 degC s-1
                huonT_ = HuonT[0,:,J,I-1].values
                hvomT_ = HvomT[0,:,J-1,I].values
                HUON_VOM_T[J,I] = np.sum((-hvomT_)+(huonT-huonT_))
            elif J==(ylen-2):
                hvomT = HvomT[0,:,J,I].values
                huonT_ = HuonT[0,:,J,I-1].values
                hvomT_ = HvomT[0,:,J-1,I].values
                HUON_VOM_T[J,I] = np.sum((hvomT-hvomT_)+(-huonT_))
            elif I==(xlen-2) and J==(ylen-2):
                huonT_ = HuonT[0,:,J,I-1].values
                hvomT_ = HvomT[0,:,J-1,I].values
                HUON_VOM_T[J,I] = np.sum((-hvomT_)+(-huonT_))
            elif I==0 and J==(ylen-1):
                huonT_ = HuonT[0,:,J,I].values
                hvomT = HvomT[0,:,J-1,I].values
                HUON_VOM_T[J,I] = np.sum((hvomT)+(huonT))
            elif J==0 and I==(xlen-1):
                huonT = HuonT[0,:,J,I-1].values
                hvomT = HvomT[0,:,J,I].values
                HUON_VOM_T[J,I] = np.sum((hvomT)+(huonT))
            else:
                huonT = HuonT[0,:,J,I].values # m3 degC s-1
                hvomT = HvomT[0,:,J,I].values
                huonT_ = HuonT[0,:,J,I-1].values
                hvomT_ = HvomT[0,:,J-1,I].values
                HUON_VOM_T[J,I] = np.sum((hvomT-hvomT_)+(huonT-huonT_))

    scipy.io.savemat('HADV.mat', {'HADV': HADV})
    scipy.io.savemat('HUON_VOM_T.mat', {'HUON_VOM_T': HUON_VOM_T})
