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

    # load waom4 3D temp field to plot some maps
    ds = xr.open_mfdataset(paths='/g/data/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom4extend_shflim_S_0.25Q/output_01-20yr/ocean_avg_0010.nc')

    temp3d_4km= ds.temp

    ds = ds.set_coords(['Cs_r', 'Cs_w', 'hc', 'h', 'Vtransform'])
    # calc dz:
    hwater = ds.h -abs(ds.zice) # replace ds.h for hwater below
    Zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * hwater) / (ds.hc + hwater)
    z_rho3d_4km = ds.zeta + (ds.zeta + hwater) * Zo_rho -abs(ds.zice)
    del Zo_rho
    ds.close()

     # read grid file for lon/lat coordinates
    dg = xr.open_dataset("/g/data3/hh5/tmp/access-om/fbd581/ROMS/waom4_frc/waom4extend_grd.nc")
    angle_4km= dg.variables["angle"]
    lat_rho_4km= dg.variables["lat_rho"]
    lon_rho_4km = dg.variables["lon_rho"]
    lat_u_4km= dg.variables["lat_u"]
    lon_u_4km = dg.variables["lon_u"]
    lat_v_4km= dg.variables["lat_v"]
    lon_v_4km = dg.variables["lon_v"]
    pm_4km = dg.variables["pm"]
    pn_4km = dg.variables["pn"]
    zice_4km = dg.variables["zice"]
    h_4km = dg.variables["h"]
    dg.close()
    print('Print lon/lat_rho shapes',lon_rho_4km.shape, lat_rho_4km.shape)
    print('Print lon/lat_rho shapes',lon_rho_4km[0:-1,0:-1].shape, lat_rho_4km[0:-1,0:-1].shape)


    # Assuming temp_data is your temperature xarray.DataArray
    # and depth_data is your depth xarray.DataArray


    depth_target = -500  # The depth you want to interpolate to (500 meters)

    # Create an empty array to hold the interpolated temperature values
#    interp_temp = np.full((temp3d_4km.shape[0],temp3d_4km[0,:].shape[1], temp3d_4km[0,:].shape[2]), np.nan)

### this takes forever!
#    for otime in range(365):
#
#        # Loop through each grid point and interpolate
#        for i in range(temp3d_4km.shape[1]):
#            for j in range(temp3d_4km.shape[2]):
#                print(i,j)
#                depth_profile = z_rho3d_4km[otime, i, j, :]
#                temp_profile = temp3d_4km[otime, :, i, j]
#
#                if np.any(np.isnan(depth_profile)) or np.any(np.isnan(temp_profile)):
#                    print('This i,j is nan!')
#                    continue  # Skip if there are NaNs in the profile
#
#                # Find the indices of the two layers surrounding the target depth
#                idx_above = np.where(depth_profile <= depth_target)[0].max()
#                idx_below = np.where(depth_profile > depth_target)[0].min()
#
#                # Linear interpolation
#                if idx_above == idx_below:
#                    interp_temp[i, j] = temp_profile[idx_above]
#                else:
#                    depth_above = depth_profile[idx_above]
#                    depth_below = depth_profile[idx_below]
#                    temp_above = temp_profile[idx_above]
#                    temp_below = temp_profile[idx_below]
#
#                    interp_temp[otime, i, j] = temp_above + (temp_below - depth_above) * (temp_below - temp_above) / (depth_below - depth_above)
#
#    # # Create a DataArray for the interpolated temperature data
#    interp_temp3d_4km = xr.DataArray(interp_temp, coords=[temp3d_4km.coords['ocean_time'], temp3d_4km.coords['eta_rho'], temp3d_4km.coords['xi_rho']], dims=['ocean_time', 'eta_rho', 'xi_rho'])
#

from scipy.interpolate import RegularGridInterpolator

def interpolate_to_depth_3d(data_2d, depth_levels, target_depth):
    """
    Interpolates 2D arrays with terrain-following vertical levels to a specific depth using 3D interpolation.

    Parameters:
    data_2d (numpy.ndarray): 3D array with shape (levels, lat, lon)
    depth_levels (numpy.ndarray): 3D array with shape (levels, lat, lon) representing the depths of each level
    target_depth (float): The specific depth to interpolate to

    Returns:
    numpy.ndarray: 2D array with shape (lat, lon) interpolated to the target depth
    """
    levels, lat, lon = data_2d.shape

    # Create a meshgrid for the original coordinates
    level_coords = np.arange(levels)
    lat_coords = np.arange(lat)
    lon_coords = np.arange(lon)

    # Create the interpolator function
    interpolator = RegularGridInterpolator((level_coords, lat_coords, lon_coords), data_2d, bounds_error=False, fill_value=None)

    # Create a meshgrid for the target depth
    target_depths = np.full((lat, lon), target_depth)

    # Flatten the arrays for interpolation
    points = np.array([depth_levels.flatten(), np.repeat(lat_coords, lon), np.tile(lon_coords, lat)]).T

    # Interpolate the data
    interpolated_data_flat = interpolator(points)

    # Reshape the interpolated data back to 2D
    interpolated_data = interpolated_data_flat.reshape((lat, lon))

    return interpolated_data

    interp_temp_annual = np.full((temp3d_4km[0,:].shape[1], temp3d_4km[0,:].shape[2]), np.nan)
    interp_temp_feb = np.full((temp3d_4km[0,:].shape[1], temp3d_4km[0,:].shape[2]), np.nan)
    interp_temp_aug = np.full((temp3d_4km[0,:].shape[1], temp3d_4km[0,:].shape[2]), np.nan)
    # Loop through each grid point and interpolate
#    for i in range(temp3d_4km.shape[1]):
#        for j in range(temp3d_4km.shape[2]):
            print(i,j)
    depth_annual = z_rho3d_4km.mean('ocean_time')
    temp_annual = temp3d_4km.mean('ocean_time')
    depth_feb = z_rho3d_4km[31:60, :].mean('ocean_time')
    temp_feb = temp3d_4km[31:60, :].mean('ocean_time')
    depth_aug = z_rho3d_4km[213:244, :].mean('ocean_time')
    temp_aug = temp3d_4km[213:244, :].mean('ocean_time')

## annual
    depth_profile = depth_annual
    temp_profile = temp_annual


## february
            depth_profile = depth_feb
            temp_profile = temp_feb
            print('Tracking crash; print depth_profile shape before finding indices of the two layers above/below:')
            print(depth_profile.shape)
            print('========================================================')

            # Find the indices of the two layers surrounding the target depth
            idx_above = np.where(depth_profile <= depth_target)[0].max()
            idx_below = np.where(depth_profile > depth_target)[0].min()
            # Linear interpolation
            if idx_above == idx_below:
                interp_temp[i, j] = temp_profile[idx_above]
            else:
                depth_above = depth_profile[idx_above]
                depth_below = depth_profile[idx_below]
                temp_above = temp_profile[idx_above]
                temp_below = temp_profile[idx_below]
                interp_temp_feb[i, j] = temp_above + (temp_below - depth_above) * (temp_below - temp_above) / (depth_below - depth_above)
            del depth_profile, temp_profile

## august
            depth_profile = depth_aug
            temp_profile = temp_aug
            # Find the indices of the two layers surrounding the target depth
            idx_above = np.where(depth_profile <= depth_target)[0].max()
            idx_below = np.where(depth_profile > depth_target)[0].min()
            # Linear interpolation
            if idx_above == idx_below:
                interp_temp[i, j] = temp_profile[idx_above]
            else:
                depth_above = depth_profile[idx_above]
                depth_below = depth_profile[idx_below]
                temp_above = temp_profile[idx_above]
                temp_below = temp_profile[idx_below]
                interp_temp_aug[i, j] = temp_above + (temp_below - depth_above) * (temp_below - temp_above) / (depth_below - depth_above)
            del depth_profile, temp_profile

    # # Create a DataArray for the interpolated temperature data
    interp_temp3d_4km_annual = xr.DataArray(interp_temp_annual, coords=[temp3d_4km.coords['eta_rho'], temp3d_4km.coords['xi_rho']], dims=['eta_rho', 'xi_rho'])
    interp_temp3d_4km_feb = xr.DataArray(interp_temp_feb, coords=[temp3d_4km.coords['eta_rho'], temp3d_4km.coords['xi_rho']], dims=['eta_rho', 'xi_rho'])
    interp_temp3d_4km_aug = xr.DataArray(interp_temp_aug, coords=[temp3d_4km.coords['eta_rho'], temp3d_4km.coords['xi_rho']], dims=['eta_rho', 'xi_rho'])


    # Add the interpolated data to the original dataset (optional)
    #dataset = temp3d_4km.to_dataset(name='temp')
    dataset['interp_temp_at_500m_annual'] = interp_temp3d_4km_annual
    dataset['interp_temp_at_500m_feb'] = interp_temp3d_4km_feb
    dataset['interp_temp_at_500m_aug'] = interp_temp3d_4km_aug

    files_path = '/g/data/hh5/tmp/access-om/fbd581/ROMS/postprocessing/cross_contour_tmp/'
##    dataset.to_netcdf(files_path + 'WAOM4_temp_500m_daily.nc', mode='w', format="NETCDF4")
    dataset.to_netcdf(files_path + 'WAOM4_temp_500m_avg.nc', mode='w', format="NETCDF4")

