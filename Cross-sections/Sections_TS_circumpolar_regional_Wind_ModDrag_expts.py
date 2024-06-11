# read nc output from WAOM 4km run
import venv

import xarray as xr
import pandas as p
import numpy as np
import numpy.ma as ma
#import cartopy.crs as ccrs
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as dates
from matplotlib.ticker import NullFormatter
from matplotlib.colors import LinearSegmentedColormap   # for custom colormaps

from datetime import datetime, timedelta

from netCDF4 import Dataset
from netCDF4 import num2date, date2num
from matplotlib.ticker import NullFormatter
from matplotlib.colors import LinearSegmentedColormap   # for custom colormaps

#import iris
#import iris.iterate
#import iris.coords
#import iris.plot as iplt
from xgcm import Grid

import gsw

import pyresample

from dask.distributed import Client
import logging
import warnings
warnings.filterwarnings('ignore')

import sys

longitude = float(sys.argv[1])
print('Processing longitude = ', longitude)

if __name__== '__main__':

    client = Client(threads_per_worker=1, memory_limit=0, silence_logs=logging.ERROR)
    print(client)

    # =====================================================================
    # LOADING MODEL OUTPUT:
    #
    # read grid file for lon/lat coordinates
    dg = xr.open_dataset("/g/data/hh5/tmp/access-om/fbd581/ROMS/waom10_frc/waom10extend_grd.nc")
    lat_rho_10km= dg.variables["lat_rho"]
    lon_rho_10km = dg.variables["lon_rho"]
    lat_u_10km= dg.variables["lat_u"]
    lon_u_10km = dg.variables["lon_u"]
    lat_v_10km= dg.variables["lat_v"]
    lon_v_10km = dg.variables["lon_v"]
    cor_10km = dg.variables["f"]
    pm_10km = dg.variables["pm"]
    pn_10km = dg.variables["pn"]
    zice_10km = dg.variables["zice"]
    h_10km = dg.variables["h"]
    dg.close()

    ds = xr.open_mfdataset(paths='/g/data/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom10extend_shflim_S_0.25Q/output_01-20yr/ocean_avg_0020.nc')#, chunks={'eta_rho': '200MB'}, parallel=True, decode_times=False)
    temp= ds.variables["temp"]
    salt= ds.variables["salt"]
    zeta= ds.variables["zeta"]
    melt= ds.variables["m"]
    ds = ds.set_coords(['Cs_r', 'Cs_w', 'hc', 'h', 'Vtransform'])
    # calc dz:
    hwater = ds.h- abs(ds.zice) # replace ds.h for hwater below
    Zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * hwater) / (ds.hc + hwater)
    z_rho = ds.zeta + (ds.zeta + hwater) * Zo_rho - abs(ds.zice)
    Zo_w = (ds.hc * ds.s_w + ds.Cs_w * hwater) / (ds.hc + hwater)
    z_w = ds.zeta + (ds.zeta + hwater) * Zo_w - abs(ds.zice)
    del Zo_rho, Zo_w

    ds.close()

    ds = xr.open_mfdataset(paths='/g/data/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom10extend_shflim_S_0.25Q/output_01-20yr_ModDragA/ocean_avg_0020.nc')#, chunks={'eta_rho': '200MB'}, parallel=True, decode_times=False)
    temp_mdA= ds.variables["temp"]
    salt_mdA= ds.variables["salt"]
    zeta_mdA= ds.variables["zeta"]
    melt_mdA= ds.variables["m"]
    ds = ds.set_coords(['Cs_r', 'Cs_w', 'hc', 'h', 'Vtransform'])
    # calc dz:
    hwater = ds.h- abs(ds.zice) # replace ds.h for hwater below
    Zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * hwater) / (ds.hc + hwater)
    z_rho_mdA = ds.zeta + (ds.zeta + hwater) * Zo_rho - abs(ds.zice)
    Zo_w = (ds.hc * ds.s_w + ds.Cs_w * hwater) / (ds.hc + hwater)
    z_w_mdA = ds.zeta + (ds.zeta + hwater) * Zo_w - abs(ds.zice)
    del Zo_rho, Zo_w
    ds.close()


    ds = xr.open_mfdataset(paths='/g/data/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom10extend_shflim_S_0.25Q/output_01-20yr_ModDragB/ocean_avg_0020.nc')#, chunks={'eta_rho': '200MB'}, parallel=True, decode_times=False)
    temp_mdB= ds.variables["temp"]
    salt_mdB= ds.variables["salt"]
    zeta_mdB= ds.variables["zeta"]
    melt_mdB= ds.variables["m"]
    # calc dz:
    hwater = ds.h- abs(ds.zice) # replace ds.h for hwater below
    Zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * hwater) / (ds.hc + hwater)
    z_rho_mdB = ds.zeta + (ds.zeta + hwater) * Zo_rho - abs(ds.zice)
    Zo_w = (ds.hc * ds.s_w + ds.Cs_w * hwater) / (ds.hc + hwater)
    z_w_mdB = ds.zeta + (ds.zeta + hwater) * Zo_w - abs(ds.zice)
    del Zo_rho, Zo_w
    ds.close()

    ds = xr.open_mfdataset(paths='/g/data/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom10extend_shflim_S_0.25Q/output_01-20yr_ModDragC/ocean_avg_0020.nc')#, chunks={'eta_rho': '200MB'}, parallel=True, decode_times=False)
    # ds = xr.open_mfdataset(path_ECCO2_mdC + "ocean_avg_00*.nc")
    temp_mdC= ds.variables["temp"]
    salt_mdC= ds.variables["salt"]
    zeta_mdC= ds.variables["zeta"]
    melt_mdC= ds.variables["m"]
    ds = ds.set_coords(['Cs_r', 'Cs_w', 'hc', 'h', 'Vtransform'])
    # calc dz:
    hwater = ds.h- abs(ds.zice) # replace ds.h for hwater below
    Zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * hwater) / (ds.hc + hwater)
    z_rho_mdC = ds.zeta + (ds.zeta + hwater) * Zo_rho - abs(ds.zice)
    Zo_w = (ds.hc * ds.s_w + ds.Cs_w * hwater) / (ds.hc + hwater)
    z_w_mdC = ds.zeta + (ds.zeta + hwater) * Zo_w - abs(ds.zice)
    del Zo_rho, Zo_w
    ds.close()

    ds = xr.open_mfdataset(paths='/g/data/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom10extend_shflim_S_0.25Q/output_01-20yr_ModDragD/ocean_avg_0020.nc')#, chunks={'eta_rho': '200MB'}, parallel=True, decode_times=False)
    # ds = xr.open_mfdataset(path_ECCO2_mdD + "ocean_avg_00*.nc")
    temp_mdD= ds.variables["temp"]
    salt_mdD= ds.variables["salt"]
    zeta_mdD= ds.variables["zeta"]
    melt_mdD= ds.variables["m"]
    ds = ds.set_coords(['Cs_r', 'Cs_w', 'hc', 'h', 'Vtransform'])
    # calc dz:
    hwater = ds.h- abs(ds.zice) # replace ds.h for hwater below
    Zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * hwater) / (ds.hc + hwater)
    z_rho_mdD = ds.zeta + (ds.zeta + hwater) * Zo_rho - abs(ds.zice)
    Zo_w = (ds.hc * ds.s_w + ds.Cs_w * hwater) / (ds.hc + hwater)
    z_w_mdD = ds.zeta + (ds.zeta + hwater) * Zo_w - abs(ds.zice)
    del Zo_rho, Zo_w
    ds.close()

    mask_zice_10km = ma.masked_where(zice_10km < 0, np.ones(zice_10km.shape))
    mask_outice_10km = ma.masked_where(zice_10km >= 0, np.ones(zice_10km.shape))
    mask_shelf_10km = ma.masked_where(h_10km > 2000, np.ones(zice_10km.shape))

    # load 1500 and calving front contours:
    tmp_files_dir = '/g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/cross_contour_tmp/'
    expt = 'WAOM10'
    ds = xr.open_dataset(tmp_files_dir + expt + '_lon_along_1500m')
    lon_along_10km_shelf = ds.variables["one"]
    ds.close()
    ds = xr.open_dataset(tmp_files_dir + expt + '_lat_along_1500m')
    lat_along_10km_shelf = ds.variables["two"]
    ds.close()
    ds = xr.open_dataset(tmp_files_dir + expt + '_lon_along_CalvingFront')
    lon_along_10km_CF = ds.variables["one"]
    ds.close()
    ds = xr.open_dataset(tmp_files_dir + expt + '_lat_along_CalvingFront')
    lat_along_10km_CF = ds.variables["two"]
    ds.close()

    ## load 4D vars:

    temp.load()
    temp_mdA.load()
    temp_mdB.load()
    temp_mdC.load()
    temp_mdD.load()
    
    print('WAOM output loading completed.')
    # =====================================================================


    # =====================================================================
    # Do cross-sections using interpolated sections (per longitude):

    def extract_merid_section(longitude,tracer_var,z_var): # longitude could be between -180:180

        # 1) create vector (using lat.min:lat.max; then apply mask to remove values north of the 1500m isobath)
        # vector should have same resolution as the original data (to minimise errors due to non-exact zonal transect)

        minlat = -90
        maxlat = -60

        lat_maxpts = np.ceil(30/0.091985) # took 3.5min
        ##lat_maxpts = np.ceil(30/0.5) # took
        lat_vector = np.linspace(minlat,maxlat,num=int(lat_maxpts))

        lon=longitude # west of PIT trough
        lon_vector = np.linspace(lon,lon,num=int(lat_maxpts))  # replace 0 -> lon_bin looping

        # Create a meshgrid from the longitude and latitude vectors
        lon_mesh, lat_mesh = np.meshgrid(lon_vector, lat_vector)


        # Create a new DataArray for the vector
        vector_da = xr.DataArray(
            np.zeros_like(lon_mesh),  # Placeholder values (will be overwritten)
            dims=['lat_rho', 'lon_rho'],
            coords={'lat_rho': lat_vector, 'lon_rho': lon_vector}
        )

        ##re-grid high-res zice/h to 10km grid:
        w10_def = pyresample.geometry.SwathDefinition(lons=lon_rho_10km,lats=lat_rho_10km)
        w10u_def = pyresample.geometry.SwathDefinition(lons=lon_u_10km,lats=lat_u_10km)
        w10v_def = pyresample.geometry.SwathDefinition(lons=lon_v_10km,lats=lat_v_10km)
        transect_def = pyresample.geometry.SwathDefinition(lons=vector_da.lon_rho,lats=vector_da.lat_rho)
        wf = lambda r: 1/r

        h_merid_transect = pyresample.kd_tree.resample_custom(w10_def,h_10km.values,transect_def,\
                                                 radius_of_influence=10000,neighbours=4,weight_funcs=wf)

        temp_merid_transect = np.empty((12,31,len(lon_vector)))
        Hz_merid_transect = np.empty((12,31,len(lon_vector)))

        lat_31lev = np.empty(Hz_merid_transect[0,:,:].shape)
        for zz in range(31):
            lat_31lev[zz,:] = vector_da.lat_rho

        for tt in range(0,12):
            for zz in range(0,31):
                temp_merid_transect[tt,zz,:] = pyresample.kd_tree.resample_custom(w10_def,tracer_var.isel(s_rho=zz,ocean_time=tt).values,transect_def,\
                                                 radius_of_influence=10000,neighbours=4,weight_funcs=wf)
                Hz_merid_transect[tt,zz,:] = pyresample.kd_tree.resample_custom(w10_def,z_var.isel(s_rho=zz,ocean_time=tt).values,transect_def,\
                                                 radius_of_influence=10000,neighbours=4,weight_funcs=wf)

        return lat_31lev, temp_merid_transect, Hz_merid_transect

    # =====================================================================

    # =====================================================================
    # Call function and run individual sections/experiments.
    #
    # jupyter notebook using normal queue/XXLargeMem times:
    # 1 (run control): 4min23s
    # 2 (4 ModDrag expts): 56min27s (~14min per expt)
    # similar queue should take ~1h

    # - now longitude is set as input:
    #longitude=-165
    #longitude = input("Enter something: ")

    lat_31lev, temp_merid_transect, Hz_merid_transect = extract_merid_section(longitude,temp,z_rho)
    lat_31lev, temp_mdA_merid_transect, Hz_mdA_merid_transect = extract_merid_section(longitude,temp_mdA,z_rho_mdA)
    lat_31lev, temp_mdB_merid_transect, Hz_mdB_merid_transect = extract_merid_section(longitude,temp_mdB,z_rho_mdB)
    lat_31lev, temp_mdC_merid_transect, Hz_mdC_merid_transect = extract_merid_section(longitude,temp_mdC,z_rho_mdC)
    lat_31lev, temp_mdD_merid_transect, Hz_mdD_merid_transect = extract_merid_section(longitude,temp_mdD,z_rho_mdD)

    print('Section interpolation completed.')
    # =====================================================================
    # convert to xarray and save to a netcdf file:

    # path:
    files_path='/g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/cross-sections/'

    # first convert do xarray dataset:
    coordinatesSec=dict(s_layer=(['s_layer'], np.arange(0,31,1)),
            lat_pts=(['lat_pts'], np.arange(0,len(temp_merid_transect[0,0,:]),1)))

    coordinatesSecT=dict(time_months=(['time_months'],np.arange(0,12,1)),
            s_layer=(['s_layer'], np.arange(0,31,1)),
            lat_pts=(['lat_pts'], np.arange(0,len(temp_merid_transect[0,0,:]),1)))

    temp_merid_transect_xr = xr.DataArray(temp_merid_transect, coords = coordinatesSecT, dims = ['time_months','s_layer','lat_pts'])
    temp_mdA_merid_transect_xr = xr.DataArray(temp_mdA_merid_transect, coords = coordinatesSecT, dims = ['time_months','s_layer','lat_pts'])
    temp_mdB_merid_transect_xr = xr.DataArray(temp_mdB_merid_transect, coords = coordinatesSecT, dims = ['time_months','s_layer','lat_pts'])
    temp_mdC_merid_transect_xr = xr.DataArray(temp_mdC_merid_transect, coords = coordinatesSecT, dims = ['time_months','s_layer','lat_pts'])
    temp_mdD_merid_transect_xr = xr.DataArray(temp_mdD_merid_transect, coords = coordinatesSecT, dims = ['time_months','s_layer','lat_pts'])

    Hz_merid_transect_xr = xr.DataArray(Hz_merid_transect, coords = coordinatesSecT, dims = ['time_months','s_layer','lat_pts'])
    Hz_mdA_merid_transect_xr = xr.DataArray(Hz_mdA_merid_transect, coords = coordinatesSecT, dims = ['time_months','s_layer','lat_pts'])
    Hz_mdB_merid_transect_xr = xr.DataArray(Hz_mdB_merid_transect, coords = coordinatesSecT, dims = ['time_months','s_layer','lat_pts'])
    Hz_mdC_merid_transect_xr = xr.DataArray(Hz_mdC_merid_transect, coords = coordinatesSecT, dims = ['time_months','s_layer','lat_pts'])
    Hz_mdD_merid_transect_xr = xr.DataArray(Hz_mdD_merid_transect, coords = coordinatesSecT, dims = ['time_months','s_layer','lat_pts'])

    lat_31lev_xr = xr.DataArray(lat_31lev, coords = coordinatesSec, dims = ['s_layer','lat_pts'])

    print('conversion to xarray dataset completed.')
    # rename vars:
    temp_merid_transect_xr.name = 'temp'
    temp_mdA_merid_transect_xr.name = 'temp'
    temp_mdB_merid_transect_xr.name = 'temp'
    temp_mdC_merid_transect_xr.name = 'temp'
    temp_mdD_merid_transect_xr.name = 'temp'

    Hz_merid_transect_xr.name = 'Hz'
    Hz_mdA_merid_transect_xr.name = 'Hz'
    Hz_mdB_merid_transect_xr.name = 'Hz'
    Hz_mdC_merid_transect_xr.name = 'Hz'
    Hz_mdD_merid_transect_xr.name = 'Hz'

    lat_31lev_xr.name = 'latitude_rho'
    print('renaming vars completed')

    # save to netcdf
    lat_31lev_xr.to_netcdf(files_path + 'Lat_cross-sec_' + str(longitude) + '_WAOM10.nc', mode='w', format="NETCDF4")
    temp_merid_transect_xr.to_netcdf(files_path + 'Temp_cross-sec_' + str(longitude) + '_WAOM10.nc', mode='w', format="NETCDF4")
    Hz_merid_transect_xr.to_netcdf(files_path + 'Hz_cross-sec_' + str(longitude) + '_WAOM10.nc', mode='w', format="NETCDF4")

    temp_mdA_merid_transect_xr.to_netcdf(files_path + 'Temp_cross-sec_' + str(longitude) + '_WAOM10_ModDragA.nc', mode='w', format="NETCDF4")
    Hz_mdA_merid_transect_xr.to_netcdf(files_path + 'Hz_cross-sec_' + str(longitude) + '_WAOM10_ModDragA.nc', mode='w', format="NETCDF4")

    temp_mdB_merid_transect_xr.to_netcdf(files_path + 'Temp_cross-sec_' + str(longitude) + '_WAOM10_ModDragB.nc', mode='w', format="NETCDF4")
    Hz_mdB_merid_transect_xr.to_netcdf(files_path + 'Hz_cross-sec_' + str(longitude) + '_WAOM10_ModDragB.nc', mode='w', format="NETCDF4")

    temp_mdC_merid_transect_xr.to_netcdf(files_path + 'Temp_cross-sec_' + str(longitude) + '_WAOM10_ModDragC.nc', mode='w', format="NETCDF4")
    Hz_mdC_merid_transect_xr.to_netcdf(files_path + 'Hz_cross-sec_' + str(longitude) + '_WAOM10_ModDragC.nc', mode='w', format="NETCDF4")

    temp_mdD_merid_transect_xr.to_netcdf(files_path + 'Temp_cross-sec_' + str(longitude) + '_WAOM10_ModDragD.nc', mode='w', format="NETCDF4")
    Hz_mdD_merid_transect_xr.to_netcdf(files_path + 'Hz_cross-sec_' + str(longitude) + '_WAOM10_ModDragD.nc', mode='w', format="NETCDF4")

    print('Saving processed sections completed.')


