# read nc output from WAOM 10km run

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
import sys

from dask.distributed import Client
import logging
import warnings
warnings.filterwarnings('ignore')

if __name__== '__main__':

    # Get the arguments passed to the script
    Bin_index = int(sys.argv[1] )# if len(sys.argv) > 1 else 0

    client = Client(threads_per_worker=1, memory_limit=0, silence_logs=logging.ERROR)
    print(client)

    tmp_files_dir = '/g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/cross_contour_tmp/'

    # determine constants:
    rho0 = 1025 # kg. m-3
    Cp = 3989.245 # J.kg-1.degC-1
    Tf = -1.95 # degC

    ds = xr.open_mfdataset(paths='/g/data/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom10extend_shflim_S_0.25Q/output_01-20yr/ocean_avg_0010.nc')#, chunks={'eta_rho': '200MB'}, parallel=True, decode_times=False)
    temp3d_10km= ds.variables["temp"]
    ds = ds.set_coords(['Cs_r', 'Cs_w', 'hc', 'h', 'Vtransform'])
    Zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * ds.h) / (ds.hc + ds.h)
    z_rho3d_10km = ds.zeta + (ds.zeta + ds.h) * Zo_rho + ds.zice
    del Zo_rho
    ds.close()

    dg = xr.open_dataset("/g/data3/hh5/tmp/access-om/fbd581/ROMS/waom10_frc/waom10extend_grd.nc")
    lat_rho_10km= dg.variables["lat_rho"]
    lon_rho_10km = dg.variables["lon_rho"]
    lat_u_10km= dg.variables["lat_u"]
    lon_u_10km = dg.variables["lon_u"]
    lat_v_10km= dg.variables["lat_v"]
    lon_v_10km = dg.variables["lon_v"]
    pm_10km = dg.variables["pm"]
    pn_10km = dg.variables["pn"]
    zice_10km = dg.variables["zice"]
    h_10km = dg.variables["h"]
    dg.close()

    temp3d_10km.load()
    z_rho3d_10km.load()

    # load ice draft to create masks
    expt = 'WAOM10'

    ds = xr.open_dataset(tmp_files_dir + expt + '_vol_trans_1500m_daily_v3')
    vol_transport_10km = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_heat_trans_1500m_daily_v3')
    heat_transport_10km = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_temp_1500m_daily_v3')
    temp_10km = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_salt_1500m_daily_v3')
    salt_10km = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_z_rho_1500m_daily_v3')
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


    [dist_along_axis_10km, Zaxis_10km] = np.meshgrid(dist_along_10km, np.arange(0,31))

    fig_path = '/g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/figs/OHB_shelf/'

    Tf = -3.534879684448242

    # WAOM10

    Tf_heat_transp_10km = Tf*vol_transport_10km.mean('ocean_time') # it doesn't matter if the mean is done now or later; Tf is cte.

    heat_transp_10km = heat_transport_10km - Tf_heat_transp_10km

    # Total HT (online diag) annual-avg:
    annual_transp_10km = heat_transp_10km.mean('ocean_time')

    # calculate Mean as v_bar*T_bar, vert-discretised version:
    Mean_transp_10km = vol_transport_10km.mean('ocean_time')*(temp_10km.mean('ocean_time')-Tf)
    # temp_tf_10km=temp_10km-Tf
    # Mean_transp_10km = vol_transport_10km.mean('ocean_time')*temp_tf_10km.mean('ocean_time')

    # Eddy heat transport: total - Mean
    Eddy_transp_10km = np.empty((365,31,2690))
    for mm in np.arange(0,365):
        # Eddy_transp_10km[mm,:,:] = (heat_transport_10km[mm,:,:]) - Mean_transp_10km
        Eddy_transp_10km[mm,:,:] = (heat_transp_10km[mm,:,:]) - Mean_transp_10km #- Tf_heat_transp_10km



    # Zonal Convergence calculation

    # Re-index my contour so the longitude is mostly monotonic:

    shift_index1 = 191# np.argmax(lon_along_10km.values < 0) # = 192
    # use first point after stoping going back and forth ~179E/~179W.
    shift_index2=226
    shift_index3=192
    shift_index4=215
    shift_index5=216
    shift_index6=225

    # Split the longitude array at the shift index and concatenate in reverse order
    lon_reindexed = np.concatenate((lon_along_10km[shift_index3:shift_index4],
                                    lon_along_10km[shift_index2:],
                                    lon_along_10km[:shift_index1],
                                    lon_along_10km[shift_index5:shift_index6]))
    lat_reindexed = np.concatenate((lat_along_10km[shift_index3:shift_index4],
                                    lat_along_10km[shift_index2:],
                                    lat_along_10km[:shift_index1],
                                    lat_along_10km[shift_index5:shift_index6]))

    # Create a new xarray.DataArray with the reindexed longitudes
    lon_da_10km = xr.DataArray(lon_reindexed, dims=['contour_index_array'])
    lat_da_10km = xr.DataArray(lat_reindexed, dims=['contour_index_array'])

    print(lon_da_10km.shape, lon_da_10km)
    # Update the longitudes in your dataset with the reindexed values
    # -> use lon_da_10km/lat_da_10km instead of lon_along_10km/lat_along_10km

    # re-index heat transport and other variables:
    heat_transp_10km_reindexed = np.concatenate((heat_transp_10km[:,:,shift_index3:shift_index4],
                                    heat_transp_10km[:,:,shift_index2:],
                                    heat_transp_10km[:,:,:shift_index1],
                                    heat_transp_10km[:,:,shift_index5:shift_index6]), axis=-1)

    # Create a new xarray.DataArray with the reindexed longitudes
    heat_transp_da_10km = xr.DataArray(heat_transp_10km_reindexed, dims=['ocean_time','s_rho','contour_index_array'])


    heat_transp_10km_reindexed = np.concatenate((heat_transp_10km[:,:,shift_index3:shift_index4],
                                    heat_transp_10km[:,:,shift_index2:],
                                    heat_transp_10km[:,:,:shift_index1],
                                    heat_transp_10km[:,:,shift_index5:shift_index6]), axis=-1)

    # Create a new xarray.DataArray with the reindexed longitudes
    heat_transp_da_10km = xr.DataArray(heat_transp_10km_reindexed, dims=['ocean_time','s_rho','contour_index_array'])

    # 1) bin the along contour data into larger longitude intervals:

    # convert to longitude coordinate and average into 3 degree longitude bins:
    # in degrees:
    bin_width = 3
    bin_spacing = 0.25
    lon_west = -180
    lon_east = 180

    # new coordinate and midpoints of longitude bins:
    full_lon_coord = np.arange(lon_west,lon_east+bin_spacing,bin_spacing)
    lon_bin_midpoints = np.arange(lon_west+bin_width/2,lon_east-bin_width/2,bin_spacing)
    n_bin_edges = len(full_lon_coord)

    # sum into longitude bins:
    # need to be very careful of loops, we can't just mask over longitude values, but instead pick indices
    # on the isobath contour and sum continously along contour between defined indices.
    # (i.e. lon_along_contour is not monotonic)
    # find points on contour to define edges of longitude bins:
    bin_edge_indices = np.zeros(n_bin_edges)
    for lon_bin in range(n_bin_edges-1):
        # print(lon_bin)
            # find first isobath point that has the right longitude:
        first_point = np.where(lon_da_10km>=full_lon_coord[lon_bin])[0][0]
            # then find all other isobath points with the same longitude as that first point:
        same_lon_points = np.where(lon_da_10km==lon_da_10km[first_point])[0]

        # print(lat_along_10km.shape, same_lon_points)
        # we want the most southerly of these points on the same longitude line:
        bin_edge_indices[lon_bin] = same_lon_points[np.argmin(lat_da_10km[same_lon_points].values)]

    # define east/west edges:
    bin_edge_indices = bin_edge_indices.astype(int)
    bin_edge_indices_west = bin_edge_indices[:-int(bin_width/bin_spacing)-1]
    bin_edge_indices_east = bin_edge_indices[int(bin_width/bin_spacing):-1]
    n_bins = len(bin_edge_indices_west)

    lat_bin_midpoints = np.zeros(n_bins)
    for lon_bin in range(n_bins):
        # find nearest isobath point:
        lon_index = np.where(lon_da_10km>=lon_bin_midpoints[lon_bin])[0][0]
        lat_bin_midpoints[lon_bin] = lat_da_10km[lon_index]

    # sum heat transport from isobath coord into new longitude-binned coord:
    cross_slope_heat_trans = np.zeros([365,31,n_bins])
    for lon_bin in range(n_bins):
            heat_trans_this_bin0 = heat_transp_da_10km[:,:,bin_edge_indices_west[lon_bin]:bin_edge_indices_east[lon_bin]]
            cross_slope_heat_trans[:,:,lon_bin] = np.sum(heat_trans_this_bin0,axis=2)

    # Create meridional transect, interpolate variables (T,U/V,DZT), and rotate velocities to get ONLY zonal component
    # 1) Load gridded t, u, v, dzt
    ds = xr.open_mfdataset(paths='/g/data/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom10extend_shflim_S_0.25Q/output_20yr_diag_daily/ocean_avg_00*.nc')#, chunks={'eta_rho': '200MB'}, parallel=True, decode_times=False)
    temp3d_10km= ds.variables["temp"]
    u3d_10km= ds.variables["u"]
    v3d_10km= ds.variables["v"]
    ds = ds.set_coords(['Cs_r', 'Cs_w', 'hc', 'h', 'Vtransform'])
    Zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * ds.h) / (ds.hc + ds.h)
    z_rho3d_10km = ds.zeta + (ds.zeta + ds.h) * Zo_rho + ds.zice
    del Zo_rho
    ds.close()

    mask_land = ma.masked_where(zice_10km<=-.1, np.ones(h_10km.shape)) # mask out where larger than zero (i.e. everywhere)
    mask_shelf = ma.masked_where(h_10km>=1000, np.ones(h_10km.shape)) # masked out where h deeper than 1000m
    mask_cshelf = mask_land*mask_shelf


    # 2) create vector (using lat.min:lat.max; then apply mask to remove values north of the 1500m isobath)
    # vector should have same resolution as the original data (to minimise errors due to non-exact zonal transect)

    minlat = -90
    maxlat = -60

    lat_maxpts = np.ceil(30/0.091985)

    lat_vector = np.linspace(minlat,maxlat,num=int(lat_maxpts))

    ## HERE START THE LOOP through bin edges... ==================================
    for ibin in range(Bin_index, Bin_index+1):#len(full_lon_coord)):

        lon_vector = np.linspace(full_lon_coord[ibin],full_lon_coord[ibin],num=int(lat_maxpts))  # replace 0 -> lon_bin looping

        # Create a meshgrid from the longitude and latitude vectors
        lon_mesh, lat_mesh = np.meshgrid(lon_vector, lat_vector)

        # Create a new DataArray for the vector
        vector_da = xr.DataArray(
            np.zeros_like(lon_mesh),  # Placeholder values (will be overwritten)
            dims=['lat_rho', 'lon_rho'],
            coords={'lat_rho': lat_vector, 'lon_rho': lon_vector}
        )

        # re-grid high-res zice/h to 10km grid:
        w10_def = pyresample.geometry.SwathDefinition(lons=lon_rho_10km,lats=lat_rho_10km)
        w10u_def = pyresample.geometry.SwathDefinition(lons=lon_u_10km,lats=lat_u_10km)
        w10v_def = pyresample.geometry.SwathDefinition(lons=lon_v_10km,lats=lat_v_10km)
        transect_def = pyresample.geometry.SwathDefinition(lons=vector_da.lon_rho,lats=vector_da.lat_rho)
        wf = lambda r: 1/r

        h_merid_transect = pyresample.kd_tree.resample_custom(w10_def,h_10km.values,transect_def,\
                                                 radius_of_influence=10000,neighbours=4,weight_funcs=wf)
        ## Find the index of the last point before depth >= 1500m -> to plot/calc only on the shelf
        last_pt_ind = np.where(h_merid_transect >= 1500)[0][0] - 1
        last_pt_dep = h_merid_transect[last_pt_ind]
        # The first point before depth > 500m
        first_pt_ind = np.where(h_merid_transect > 50)[0][0] -1
        first_pt_dep = h_merid_transect[first_pt_ind]

        temp_merid_transect = np.empty((365,31,327))
        Hz_merid_transect = np.empty((365,31,327))
        u_merid_transect = np.empty((365,31,327))
        v_merid_transect = np.empty((365,31,327))
        for tt in range(365):
            for zz in range(0,31):
                temp_merid_transect[tt,zz,:] = pyresample.kd_tree.resample_custom(w10_def,temp3d_10km.isel(s_rho=zz,ocean_time=tt).values,transect_def,\
                                                 radius_of_influence=10000,neighbours=4,weight_funcs=wf)
                Hz_merid_transect[tt,zz,:] = pyresample.kd_tree.resample_custom(w10_def,z_rho3d_10km.isel(s_rho=zz,ocean_time=tt).values,transect_def,\
                                                 radius_of_influence=10000,neighbours=4,weight_funcs=wf)
                u_merid_transect[tt,zz,:] = pyresample.kd_tree.resample_custom(w10u_def,u3d_10km.isel(s_rho=zz,ocean_time=tt).values,transect_def,\
                                                 radius_of_influence=10000,neighbours=4,weight_funcs=wf)
                v_merid_transect[tt,zz,:] = pyresample.kd_tree.resample_custom(w10v_def,v3d_10km.isel(s_rho=zz,ocean_time=tt).values,transect_def,\
                                                 radius_of_influence=10000,neighbours=4,weight_funcs=wf)

        temp_merid_transect = temp_merid_transect[:,:,first_pt_ind:last_pt_ind]
        u_merid_transect = u_merid_transect[:,:,first_pt_ind:last_pt_ind]
        v_merid_transect = v_merid_transect[:,:,first_pt_ind:last_pt_ind]
        Hz_merid_transect = Hz_merid_transect[:,:,first_pt_ind:last_pt_ind]

        lat_31lev = np.empty(Hz_merid_transect[0,:,:].shape)
        for zz in range(31):
            lat_31lev[zz,:] = vector_da.lat_rho[first_pt_ind:last_pt_ind]

        # rotate the velocities based on the angle of the meridional transect:
        #
        lon=full_lon_coord[ibin]
        # Function to calculate rotation matrix for each longitude
        def rotation_matrix(lon):
            angle = np.radians(90-lon)  # Calculate the angle for each longitude
            # [90-lon] gives the angle to rotate conterclockwise the XY-coordinate system
            # to be truly zonal-meridional coordinates.
            return np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])

        # Calculate the rotation matrix for each point along the longitude
        rotation_matrices = rotation_matrix(lon)

        U_rotated = np.empty(u_merid_transect.shape)
        V_rotated = np.empty(v_merid_transect.shape)
        for tt in range(365):
            for zz in range(0,31):
                U = u_merid_transect[tt,zz,:]
                V = v_merid_transect[tt,zz,:]
                # Rotate the velocity vector components
                U_rotated[tt,zz,:] = U * rotation_matrices[0, 0] - V * rotation_matrices[1, 0]
                #V_rotated[tt,zz,:] = U * rotation_matrices[0, 1] - V * rotation_matrices[1, 1]

        # Calculate the zonal component (U_zonal) by setting the meridional component (V_rotated) to zero
        U_zonal = U_rotated
        
        # calculate heat transport across meridional transect:
        heat_transp_transect = U_zonal*(temp_merid_transect-Tf)*-Hz_merid_transect*10000*rho0*Cp # 10k = dxt/dyt approximately

        # testing save heat transport:
        months=np.arange(0,365)*(1/30.41667)
        coordinatesC=dict(ocean_time=months, s_rho=(['s_rho'], np.arange(0,31)),
                          transect_lat=(['transect_lat'], np.arange(0,len(temp_merid_transect[0,0,:]))))
                          
        heat_transp_transect_xr = xr.DataArray(heat_transp_transect, coords = coordinatesC, dims = ['ocean_time','s_rho','transect_lat'])
        u_zonal_transect_xr = xr.DataArray(U_zonal, coords = coordinatesC, dims = ['ocean_time','s_rho','transect_lat'])
        temp_transect_xr = xr.DataArray(temp_merid_transect, coords = coordinatesC, dims = ['ocean_time','s_rho','transect_lat'])
        dzt_transect_xr = xr.DataArray(Hz_merid_transect, coords = coordinatesC, dims = ['ocean_time','s_rho','transect_lat'])

        # rename vars:
        heat_transp_transect_xr.name = 'heat_transp'
        temp_transect_xr.name = 'temp'
        u_zonal_transect_xr.name = 'u_zonal'
        dzt_transect_xr.name = 'dzt'

        files_path = '/g/data/hh5/tmp/access-om/fbd581/ROMS/postprocessing/cross-1500m_transp/'

        heat_transp_transect_xr.to_netcdf(files_path + 'heat_transp_transect_bin=' + str(ibin) + '.nc', mode='w', format="NETCDF4")
        u_zonal_transect_xr.to_netcdf(files_path + 'u_zonal_transect_bin=' + str(ibin) + '.nc', mode='w', format="NETCDF4")
        temp_transect_xr.to_netcdf(files_path + 'temp_transect_bin=' + str(ibin) + '.nc', mode='w', format="NETCDF4")
        dzt_transect_xr.to_netcdf(files_path + 'dzt_transect_bin=' + str(ibin) + '.nc', mode='w', format="NETCDF4")

### obs:
# - should I also saved interpolated temp, u_zonal, thickness? Also salt/depth, to be able to calculate density.
#
#

