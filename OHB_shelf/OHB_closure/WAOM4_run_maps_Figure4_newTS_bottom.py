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

import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

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

    # load waom4 3D temp field to plot some maps
    # path_ECCO2_4km = '/g/data/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom4extend_shflim_S_0.25Q/output_01-20yr/ocean_avg_00*.nc'

    # ds = xr.open_mfdataset(paths='/g/data/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom4extend_shflim_S_0.25Q/output_yr10_diag_daily/ocean_avg_00*.nc')#, chunks={'eta_rho': '200MB'}, parallel=True, decode_times=False)
    # ds = xr.open_mfdataset(paths='/g/data/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom4extend_shflim_S_0.25Q/output_01-20yr/ocean_avg_00*.nc~')#, chunks={'eta_rho': '200MB'}, parallel=True, decode_times=False)
    ds = xr.open_mfdataset(paths='/g/data/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom4extend_shflim_S_0.25Q/output_01-20yr/ocean_avg_0010.nc')#, chunks={'eta_rho': '200MB'}, parallel=True, decode_times=False)
    salt3d_4km= ds.salt
    temp3d_4km= ds.temp
    u3d_4km= ds.variables["u"]
    v3d_4km= ds.variables["v"]
    melt_4km= ds.m
    shflux_4km= ds.shflux

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

    area_4km=np.divide(1,pm_4km*pn_4km)
    area_sum_4km = area_4km.sum('eta_rho').sum('xi_rho')

    # load along 1500m and CF WAOM4 variables:
    tmp_files_dir = '/g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/cross_contour_tmp/'

    expt = 'WAOM4'

    ds = xr.open_dataset(tmp_files_dir + expt + '_temp_1500m_daily_v3')
    temp_4km = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_salt_1500m_daily_v3')
    salt_4km = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_z_rho_1500m_daily_v4')
    z_rho_4km = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_dist_along_1500m_v3')
    dist_along_4km = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_lon_along_1500m_v3')
    lon_along_4km = ds.variables["one"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_lat_along_1500m_v3')
    lat_along_4km = ds.variables["two"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_temp_CF_daily_v3')
    temp_4km_CF = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_salt_CF_daily_v3')
    salt_4km_CF = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_z_rho_CF_daily_v4')
    z_rho_4km_CF = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_dist_along_CF_v3')
    dist_along_4km_CF = ds.variables["__xarray_dataarray_variable__"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_lon_along_CF_v3')
    lon_along_4km_CF = ds.variables["one"]
    ds.close()

    ds = xr.open_dataset(tmp_files_dir + expt + '_lat_along_CF_v3')
    lat_along_4km_CF = ds.variables["two"]
    ds.close()

    # define along_contour axis with same length for both resolutions:
    [dist_along_axis_4km, Zaxis_4km] = np.meshgrid(dist_along_4km, np.arange(0,31))

    [dist_along_axis_4km_CF, Zaxis_4km_CF] = np.meshgrid(dist_along_4km_CF, np.arange(0,31))

    ## load contours for plotting on top of Maps:
    files_path = '/g/data/hh5/tmp/access-om/fbd581/ROMS/postprocessing/cross_contour_tmp/'
    contour_masked_above_4km = np.load(files_path + 'WAOM4_contour_masked_above_1500m', allow_pickle=True)
    contour_masked_above_CF_4km = np.load(files_path + 'WAOM4_contour_masked_above_CF', allow_pickle=True)

    # 2) calculate OHC depth-integrated:
    Tf = -3.534879684448242
    rho0 = 1035 # kg. m-3
    Cp = 3992.1 # J.kg-1.degC-1
    temp3d_4km.load()
    temp3d_4km_K = temp3d_4km-Tf
    temp3d_vint_4km_annual = temp3d_4km_K.sum('s_rho').mean('ocean_time')*rho0*Cp
    print(temp3d_vint_4km_annual.shape)

    # 3) melt fields
    melt_4km.load()

    # define masks for cavities less_equal than 500m & more than 500m ice draft depth (to separate modes into shallow and deep):

    msk_ice_shallow = ma.masked_where(zice_4km<-500, np.ones(zice_4km.shape))
    msk_ice_deep = ma.masked_where(zice_4km>=-500, np.ones(zice_4km.shape))

    w4_mask_shelf = ma.masked_where(contour_masked_above_4km==-1000, np.ones(h_4km.shape))
    w4_mask_iceshelf = ma.masked_where(contour_masked_above_CF_4km!=-1000, np.ones(h_4km.shape))
    w4_mask_outiceshelf = ma.masked_where(contour_masked_above_CF_4km==-1000, np.ones(h_4km.shape))

    w4_mask_land = ma.masked_where(h_4km<=40, np.ones(h_4km.shape))
    w4_mask_zice = ma.masked_where(zice_4km < 0, np.ones(zice_4km.shape))

    # define constants:
    rho0=1035
    Cp=3992.1
    # Tf = -1.95 # degC
    Tf =  -3.534879684448242

    dm = xr.open_dataset('/g/data/hh5/tmp/access-om/fbd581/ROMS/postprocessing/cross_contour_tmp/WAOM4_masked_lonlatBins_v2.nc')
    mask_LonBins_4km = dm.variables['__xarray_dataarray_variable__']
    dm.close()

    file = files_path + "WAOM4_OHB_lonlatbins-v2_temp_hadv_vint_daily.nc"
    ds = xr.open_mfdataset(paths=file, chunks={'eta_rho': '200MB'}, parallel=bool)
    w4_hadv = ds.variables['__xarray_dataarray_variable__']
    w4_ocean_time = ds.ocean_time
    w4_lon_bins = ds.lon_bins
    ds.close()

    file = files_path + "WAOM4_OHB_lonlatbins-v2_temp_hdiff_vint_daily.nc"
    ds = xr.open_mfdataset(paths=file, chunks={'eta_rho': '200MB'}, parallel=bool)
    w4_hdiff = ds.variables['__xarray_dataarray_variable__']
    ds.close()

    file = files_path + "WAOM4_OHB_lonlatbins-v2_temp_vadv_vint_daily.nc"
    ds = xr.open_mfdataset(paths=file, chunks={'eta_rho': '200MB'}, parallel=bool)
    w4_vadv = ds.variables['__xarray_dataarray_variable__']
    ds.close()

    file = files_path + "WAOM4_OHB_lonlatbins-v2_temp_vdiff_vint_daily.nc"
    ds = xr.open_mfdataset(paths=file, chunks={'eta_rho': '200MB'}, parallel=bool)
    w4_vdiff = ds.variables['__xarray_dataarray_variable__']
    ds.close()

    file = files_path + "WAOM4_OHB_lonlatbins-v2_temp_hadv_vint_daily_iceshelf.nc"
    ds = xr.open_mfdataset(paths=file, chunks={'eta_rho': '200MB'}, parallel=bool)
    w4_hadv_is = ds.variables['__xarray_dataarray_variable__']
    ds.close()

    file = files_path + "WAOM4_OHB_lonlatbins-v2_temp_hdiff_vint_daily_iceshelf.nc"
    ds = xr.open_mfdataset(paths=file, chunks={'eta_rho': '200MB'}, parallel=bool)
    w4_hdiff_is = ds.variables['__xarray_dataarray_variable__']
    ds.close()

    file = files_path + "WAOM4_OHB_lonlatbins-v2_temp_vadv_vint_daily_iceshelf.nc"
    ds = xr.open_mfdataset(paths=file, chunks={'eta_rho': '200MB'}, parallel=bool)
    w4_vadv_is = ds.variables['__xarray_dataarray_variable__']
    ds.close()

    file = files_path + "WAOM4_OHB_lonlatbins-v2_temp_vdiff_vint_daily_iceshelf.nc"
    ds = xr.open_mfdataset(paths=file, chunks={'eta_rho': '200MB'}, parallel=bool)
    w4_vdiff_is = ds.variables['__xarray_dataarray_variable__']
    ds.close()

    file = files_path + "WAOM4_OHB_lonlatbins-v2_melt_daily_v2.nc"
    ds = xr.open_mfdataset(paths=file, chunks={'eta_rho': '200MB'}, parallel=bool)
    w4_bmelt = ds.variables['__xarray_dataarray_variable__']
    ds.close()

    # basal melting divided into shallow vs deep (only done for WAOM4 so far):

    file = files_path + "WAOM4_OHB_lonlatbins-v2_melt_shallow_daily_v2.nc"
    ds = xr.open_mfdataset(paths=file, chunks={'eta_rho': '200MB'}, parallel=bool)
    w4_bmelt_shallow = ds.variables['__xarray_dataarray_variable__']
    lon_bins = ds.variables['lon_bins']

    ds.close()

    file = files_path + "WAOM4_OHB_lonlatbins-v2_melt_deep_daily_v2.nc"
    ds = xr.open_mfdataset(paths=file, chunks={'eta_rho': '200MB'}, parallel=bool)
    w4_bmelt_deep = ds.variables['__xarray_dataarray_variable__']
    ds.close()

    w4_net = w4_hadv + w4_vadv + w4_hdiff + w4_vdiff
    w4_net_is = w4_hadv_is + w4_vadv_is + w4_hdiff_is + w4_vdiff_is

    # define intervals:

    # win1=90
    # win2=349
    # sum1=350
    # sum2=89

    win1=121 #Apr-Nov
    win2=335
    sum1=336 #Dec-Mar
    sum2=120

    # plt.plot(lon_bins,net.isel(ocean_time=slice(0,sum2)).mean('ocean_time')*1e-12,'--',color='r')

    # define summer net (combining end of the year + beginning of the year):
    # winter is defined while plotting...

    # WAOM4:
    w4_net_sumA=w4_net.isel(ocean_time=slice(0,sum2))
    w4_net_sumB=w4_net.isel(ocean_time=slice(sum1,-1))
    w4_net_SUM_A=xr.DataArray(w4_net_sumA)
    w4_net_SUM_B=xr.DataArray(w4_net_sumB)
    w4_net_summer = xr.concat((w4_net_SUM_A,w4_net_SUM_B), 'ocean_time')

    w4_net_is_sumA=w4_net_is.isel(ocean_time=slice(0,sum2))
    w4_net_is_sumB=w4_net_is.isel(ocean_time=slice(sum1,-1))
    w4_net_is_SUM_A=xr.DataArray(w4_net_is_sumA)
    w4_net_is_SUM_B=xr.DataArray(w4_net_is_sumB)
    w4_net_is_summer = xr.concat((w4_net_is_SUM_A,w4_net_is_SUM_B), 'ocean_time')

    w4_bmelt_sumA=w4_bmelt.isel(ocean_time=slice(0,sum2))
    w4_bmelt_sumB=w4_bmelt.isel(ocean_time=slice(sum1,-1))
    w4_bmelt_SUM_A=xr.DataArray(w4_bmelt_sumA)
    w4_bmelt_SUM_B=xr.DataArray(w4_bmelt_sumB)
    w4_bmelt_summer = xr.concat((w4_bmelt_SUM_A,w4_bmelt_SUM_B), 'ocean_time')


    w4_bmeltS_sumA=w4_bmelt_shallow.isel(ocean_time=slice(0,sum2))
    w4_bmeltS_sumB=w4_bmelt_shallow.isel(ocean_time=slice(sum1,-1))
    w4_bmeltS_SUM_A=xr.DataArray(w4_bmeltS_sumA)
    w4_bmeltS_SUM_B=xr.DataArray(w4_bmeltS_sumB)
    w4_bmeltS_summer = xr.concat((w4_bmeltS_SUM_A,w4_bmeltS_SUM_B), 'ocean_time')


    w4_bmeltD_sumA=w4_bmelt_deep.isel(ocean_time=slice(0,sum2))
    w4_bmeltD_sumB=w4_bmelt_deep.isel(ocean_time=slice(sum1,-1))
    w4_bmeltD_SUM_A=xr.DataArray(w4_bmeltD_sumA)
    w4_bmeltD_SUM_B=xr.DataArray(w4_bmeltD_sumB)
    w4_bmeltD_summer = xr.concat((w4_bmeltD_SUM_A,w4_bmeltD_SUM_B), 'ocean_time')

    # WAOM4:
    w4_hadv_sumA=w4_hadv.isel(ocean_time=slice(0,sum2))
    w4_hadv_sumB=w4_hadv.isel(ocean_time=slice(sum1,-1))
    w4_hadv_SUM_A=xr.DataArray(w4_hadv_sumA)
    w4_hadv_SUM_B=xr.DataArray(w4_hadv_sumB)
    w4_hadv_summer = xr.concat((w4_hadv_SUM_A,w4_hadv_SUM_B), 'ocean_time')

    w4_hadv_is_sumA=w4_hadv_is.isel(ocean_time=slice(0,sum2))
    w4_hadv_is_sumB=w4_hadv_is.isel(ocean_time=slice(sum1,-1))
    w4_hadv_is_SUM_A=xr.DataArray(w4_hadv_is_sumA)
    w4_hadv_is_SUM_B=xr.DataArray(w4_hadv_is_sumB)
    w4_hadv_is_summer = xr.concat((w4_hadv_is_SUM_A,w4_hadv_is_SUM_B), 'ocean_time')


    # WAOM4:
    w4_hdiff_sumA=w4_hdiff.isel(ocean_time=slice(0,sum2))
    w4_hdiff_sumB=w4_hdiff.isel(ocean_time=slice(sum1,-1))
    w4_hdiff_SUM_A=xr.DataArray(w4_hdiff_sumA)
    w4_hdiff_SUM_B=xr.DataArray(w4_hdiff_sumB)
    w4_hdiff_summer = xr.concat((w4_hdiff_SUM_A,w4_hdiff_SUM_B), 'ocean_time')

    w4_hdiff_is_sumA=w4_hdiff_is.isel(ocean_time=slice(0,sum2))
    w4_hdiff_is_sumB=w4_hdiff_is.isel(ocean_time=slice(sum1,-1))
    w4_hdiff_is_SUM_A=xr.DataArray(w4_hdiff_is_sumA)
    w4_hdiff_is_SUM_B=xr.DataArray(w4_hdiff_is_sumB)
    w4_hdiff_is_summer = xr.concat((w4_hdiff_is_SUM_A,w4_hdiff_is_SUM_B), 'ocean_time')


    # WAOM4:
    w4_vdiff_sumA=w4_vdiff.isel(ocean_time=slice(0,sum2))
    w4_vdiff_sumB=w4_vdiff.isel(ocean_time=slice(sum1,-1))
    w4_vdiff_SUM_A=xr.DataArray(w4_vdiff_sumA)
    w4_vdiff_SUM_B=xr.DataArray(w4_vdiff_sumB)
    w4_vdiff_summer = xr.concat((w4_vdiff_SUM_A,w4_vdiff_SUM_B), 'ocean_time')

    w4_vdiff_is_sumA=w4_vdiff_is.isel(ocean_time=slice(0,sum2))
    w4_vdiff_is_sumB=w4_vdiff_is.isel(ocean_time=slice(sum1,-1))
    w4_vdiff_is_SUM_A=xr.DataArray(w4_vdiff_is_sumA)
    w4_vdiff_is_SUM_B=xr.DataArray(w4_vdiff_is_sumB)
    w4_vdiff_is_summer = xr.concat((w4_vdiff_is_SUM_A,w4_vdiff_is_SUM_B), 'ocean_time')

    # treat velocities; interpolate u/v to rho-grid


    w4r_def = pyresample.geometry.SwathDefinition(lons=lon_rho_4km[:,:],lats=lat_rho_4km[:,:])
    w4u_def = pyresample.geometry.SwathDefinition(lons=lon_u_4km[:,:],lats=lat_u_4km[:,:])
    w4v_def = pyresample.geometry.SwathDefinition(lons=lon_v_4km[:,:],lats=lat_v_4km[:,:])

    wf = lambda r: 1/r

    u_rho_i = np.empty(temp3d_4km.shape)
    v_rho_i = np.empty(temp3d_4km.shape)
    for tt in np.arange(0,12):
        for kk in np.arange(0,31):
            u_rho_tmp = pyresample.kd_tree.resample_custom(w4u_def,u3d_4km[tt,kk,:,:].values,w4r_def,\
                                             radius_of_influence=30000,neighbours=1,weight_funcs=wf)
            v_rho_tmp = pyresample.kd_tree.resample_custom(w4v_def,v3d_4km[tt,kk,:,:].values,w4r_def,\
                                             radius_of_influence=30000,neighbours=1,weight_funcs=wf)
            u_rho_i[tt,kk,:] = u_rho_tmp
            v_rho_i[tt,kk,:] = v_rho_tmp
            del u_rho_tmp, v_rho_tmp

    ## Re-map u/v for plotting:
    # call cartopy projectione
    proj = ccrs.PlateCarree()

    # re-gridding
    src = ccrs.PlateCarree()
    nx, ny = 1575,1400

    scale = 1e7
    xs = np.linspace(-scale, scale, nx)
    ys = np.linspace(-scale, scale, ny)

    xs, ys = np.meshgrid(xs, ys)

    # - new (lonlatbins-v2):
    # sind_subreg = [7, 14, 24, 34, 60, 75, 104, 110, 114, 124, 135]
    # eind_subreg = [13, 19, 27, 41, 66, 79, 108, 113, 118, 129, 139]
    # - fix after inspecting each sub-regions (10/10/24):
    sind_subreg = [8, 14, 22, 32, 49, 75, 99, 108, 114, 124, 135, 0]
    eind_subreg = [12, 21, 27, 37, 67, 80, 105, 113, 118, 129, 140, 6]

    max_BM = [3, 2, 1, 2.5, 1, 1, 1,
              3.8, 2, 2, 1, 1]        # max ylim for basal melting plot (vary b/w regions)

    names_subreg = ['Sulzberger','Getz','Thwaites-PIG','Bellingshausen','FRIS',
                   'Fimbul','Amery','Shackleton','Totten','Mertz','Ross1','Ross2']
    max_BM = [3, 2, 1, 2.5, 1, 1, 1,
              3.8, 2, 2, 1, 1]        # max ylim for basal melting plot (vary b/w regions)

    Melt_fac=365.25*86400


    dayofmonth = [slice(0,32), slice(32,61), slice(61,93), slice(93,123), slice(123,155), slice(155,184), slice(184,215),
                  slice(215,247), slice(247,277), slice(277,304), slice(304,334), slice(334,365)]

    # define limits for maps of each region:

    # order is: LonMin, LonMax, LatMin, LatMax
    map_limits=[[-165, -135, -79, -74], # Sulzberger
    [-140, -110, -71, -76],             # Getz
    [-120, -90, -76, -69],              # PIG-Thwaites
    [-100, -65, -75, -67],            # Bellingshausen
    [-85, -27, -84, -71],            # FRIS
    [-15, 15, -73.5, -69],              # Fimbul
    [65, 85, -74, -65],                 # Amery
    [80, 110, -68.1, -63],              # Shackleton
    [110, 135, -68, -63.5],             # Totten
    [130,160, -70, -64],              # Mertz
    [-195,-165, -85, -69]]            # Ross
    ###[-215,-160, -84, -67]]              # Ross

    ## this piece do the maps for Figure 4 (containing annual temp, anom 1month in summer, anom 1 month in winter) + TS-diag for 2 regions (totten x PIG) at SURFACE, MID-DEPTH, BOTTOM sigma layers

    slayer=[15, -1, 0]
    zname =['Mid-depth','Surface','Bottom']
    pname =['mid-depth','surface','bottom']
    anomlim=[0.5, 2, 0.5]

    for dd in np.arange(2,3):
        # dd=2

        mm=1 # So sum_MM=February (2nd month) and win_MM=August (8th month)
        sum_MM=mm
        win_MM=mm+6

        # define subregion:
        # names_subreg = ['Sulzberger','Getz','Thwaites-PIG','Bellingshausen','FRIS',
                       # 'Fimbul','Amery','Shackleton','Totten','Mertz','Ross1','Ross2']

        for rr in range(10,11):
            Pname = names_subreg[rr]

            # define 2d Map vars for plottins: annual, summer month, winter month
            annual_T = temp3d_4km.mean('ocean_time')
            annual_S = salt3d_4km.mean('ocean_time')
            annual_Z = z_rho3d_4km.mean('ocean_time')
            annual_melt = melt_4km.mean('ocean_time')*Melt_fac

            sum_Ta = (temp3d_4km.isel(ocean_time=sum_MM)+273.15) - (temp3d_4km.mean('ocean_time')+273.15) #  Jan - Annual (+ means summer is warmer than mean)
            sum_T = (temp3d_4km.isel(ocean_time=sum_MM))
            sum_S = (salt3d_4km.isel(ocean_time=sum_MM))
            sum_Z = z_rho3d_4km.isel(ocean_time=sum_MM)

            sum_U = (u_rho_i[sum_MM,:,:])
            sum_V = (v_rho_i[sum_MM,:,:])
            sum_melt = melt_4km.isel(ocean_time=sum_MM)*Melt_fac

            win_Ta = (temp3d_4km.isel(ocean_time=win_MM)+273.15) - (temp3d_4km.mean('ocean_time')+273.15) # Jul - Anual (- mean winter is colder...)
            win_T = (temp3d_4km.isel(ocean_time=win_MM))
            win_S = (salt3d_4km.isel(ocean_time=win_MM))
            win_Z = z_rho3d_4km.isel(ocean_time=win_MM)
            win_U = (u_rho_i[win_MM,:,:])
            win_V = (v_rho_i[win_MM,:,:])
            win_melt = melt_4km.isel(ocean_time=win_MM)*Melt_fac

            # shflux
            annual_shf = shflux_4km.mean('ocean_time')
            sum_shf = shflux_4km.isel(ocean_time=sum_MM)
            win_shf = shflux_4km.isel(ocean_time=win_MM)

            # get map limits:
            LonMin=map_limits[rr][0]
            LonMax=map_limits[rr][1]
            LatMin=map_limits[rr][2]
            LatMax=map_limits[rr][3]

            #### TS-diag related definitions (along the calving front):
            # 2 - sum individual bins:
            mask_LonBins_tmp = mask_LonBins_4km[sind_subreg[rr]:eind_subreg[rr],:]
            mask_LonBins_tmp = np.nansum(mask_LonBins_tmp, axis=0)

            # 3 - multiply ice shelves mask (masked outside it) and Lon Bins:
            newmask_bin = w4_mask_outiceshelf*mask_LonBins_tmp
            # 3b - mask out values that are 0 
            mask = (newmask_bin == 1)
            newmask_bin = np.where(mask, newmask_bin, np.nan)  # Replace masked values with NaN

            # 4 - define T/S/Z variables for the TS-plot:
            temp_masked_ann = newmask_bin*annual_T
            salt_masked_ann = newmask_bin*annual_S
            Z_masked_ann = newmask_bin*annual_Z.transpose('s_rho','eta_rho','xi_rho')
            # 5 -filter & stack example:
            # 5a - Flatten the 3D array to 2D while removing NaNs
            flat_temp_ann = temp_masked_ann.stack(z=('s_rho', 'eta_rho', 'xi_rho')).values
            flat_salt_ann = salt_masked_ann.stack(z=('s_rho', 'eta_rho', 'xi_rho')).values
            flat_z_rho_ann = Z_masked_ann.stack(z=('s_rho', 'eta_rho', 'xi_rho')).values
            non_nan_mask = ~np.isnan(flat_temp_ann)
            # 5b -Filter out NaN values
            filt_temp_ann = flat_temp_ann[non_nan_mask]
            filt_salt_ann = flat_salt_ann[non_nan_mask]
            filt_z_rho_ann = flat_z_rho_ann[non_nan_mask]

## Summer:
            # 4 - define T/S/Z variables for the TS-plot:
            temp_masked_sum = newmask_bin*sum_T
            salt_masked_sum = newmask_bin*sum_S
            Z_masked_sum = newmask_bin*sum_Z.transpose('s_rho','eta_rho','xi_rho')
            # 5 -filter & stack example:
            # 5a - Flatten the 3D array to 2D while removing NaNs
            flat_temp_sum = temp_masked_sum.stack(z=('s_rho', 'eta_rho', 'xi_rho')).values
            flat_salt_sum = salt_masked_sum.stack(z=('s_rho', 'eta_rho', 'xi_rho')).values
            flat_z_rho_sum = Z_masked_sum.stack(z=('s_rho', 'eta_rho', 'xi_rho')).values
            non_nan_mask = ~np.isnan(flat_temp_sum)
            # 5b -Filter out NaN values
            filt_temp_sum = flat_temp_sum[non_nan_mask]
            filt_salt_sum = flat_salt_sum[non_nan_mask]
            filt_z_rho_sum = flat_z_rho_sum[non_nan_mask]
## Winter:
            # 4 - define T/S/Z variables for the TS-plot:
            temp_masked_win = newmask_bin*win_T
            salt_masked_win = newmask_bin*win_S
            Z_masked_win = newmask_bin*win_Z.transpose('s_rho','eta_rho','xi_rho')
            # 5 -filter & stack example:
            # 5a - Flatten the 3D array to 2D while removing NaN
            flat_temp_win = temp_masked_win.stack(z=('s_rho', 'eta_rho', 'xi_rho')).values
            flat_salt_win = salt_masked_win.stack(z=('s_rho', 'eta_rho', 'xi_rho')).values
            flat_z_rho_win = Z_masked_win.stack(z=('s_rho', 'eta_rho', 'xi_rho')).values
            non_nan_mask = ~np.isnan(flat_temp_win)
            # 5b -Filter out NaN values
            filt_temp_win = flat_temp_win[non_nan_mask]
            filt_salt_win = flat_salt_win[non_nan_mask]
            filt_z_rho_win = flat_z_rho_win[non_nan_mask]


## calling projections:
            proj_v = ccrs.PlateCarree()
            proj = ccrs.PlateCarree(central_longitude=0)#-24)
            proj_v2 = ccrs.SouthPolarStereo()
            proj2 = ccrs.SouthPolarStereo()

            # fig = plt.figure(figsize=(15,8))
            fig = plt.figure(figsize=(14,10))


            # if Ross, need to use proj2!
            if rr>=10:
                PROJ=proj2
                PROJV=proj_v2
            else:
                PROJ=proj
                PROJV=proj_v
            PROJN=proj_v

            ax2 = fig.add_subplot(3,4,5, projection=PROJ) # Anual
            ct=ax2.pcolormesh(lon_rho_4km,lat_rho_4km,annual_T.isel(s_rho=slayer[dd]).values,
                           vmin=-3, vmax=3 ,cmap=cmocean.cm.thermal, transform=PROJN)

            cm=ax2.pcolormesh(lon_rho_4km,lat_rho_4km,
                           ma.masked_where(contour_masked_above_CF_4km == -1000, annual_melt),
                           transform=PROJN, vmin=-8, vmax=8, cmap='bwr') #, cmap=cmocean.cm.balance

            ax2.quiver(lon_rho_4km[::5,::5],lat_rho_4km[::5,::5],
                       np.mean(u_rho_i[:,slayer[dd],::5,::5],axis=0),
                       np.mean(v_rho_i[:,slayer[dd],::5,::5],axis=0),scale=5,transform=PROJV)
            ax2.set_title(names_subreg[rr] + '\n' + zname[dd] + ': Annual mean')
            ax2.set_extent([LonMin, LonMax, LatMin, LatMax], crs=PROJN) # Getz
            x_left, x_right = ax2.get_xlim()
            y_low, y_high = ax2.get_ylim()
            ratio = .8
            ax2.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
            gl=ax2.gridlines(draw_labels=True,linewidth=.1)
            gl.top_labels=False   # suppress top labels
            gl.right_labels=False   # suppress right labels
            gl.left_labels=False
            gl.bottom_labels=False
            plt.scatter(lon_along_4km,lat_along_4km, s=.8, alpha=0.2, color='gold',label='1500m isobath', transform=PROJN);
            ax2.add_feature(cfeature.LAND, zorder=3, facecolor='darkgray',edgecolor='white')

            ax2 = fig.add_subplot(2,4,2, projection=PROJ) # jan, TEMP ANOM
            cta=ax2.pcolormesh(lon_rho_4km,lat_rho_4km,sum_Ta.isel(s_rho=slayer[dd]),
                           vmin=-anomlim[dd], vmax=anomlim[dd] ,cmap='twilight_shifted', transform=PROJN)
            cma=ax2.pcolormesh(lon_rho_4km,lat_rho_4km,
                           ma.masked_where(contour_masked_above_CF_4km == -1000, sum_melt),
                           transform=PROJN, vmin=-8, vmax=8, cmap='bwr') #, cmap=cmocean.cm.balance
            ax2.quiver(lon_rho_4km[::5,::5],lat_rho_4km[::5,::5],sum_U[slayer[dd],::5,::5],sum_V[slayer[dd],::5,::5],scale=5,transform=PROJV)
            ax2.set_title('Temperature anomaly \n Summer')
            ax2.set_extent([LonMin, LonMax, LatMin, LatMax], crs=PROJN) # Getz
            x_left, x_right = ax2.get_xlim()
            y_low, y_high = ax2.get_ylim()
            ratio = .8
            ax2.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
            gl=ax2.gridlines(draw_labels=True,linewidth=.1)
            gl.top_labels=False   # suppress top labels
            gl.right_labels=False   # suppress right labels
            gl.bottom_labels=False
            plt.scatter(lon_along_4km,lat_along_4km, s=.8, alpha=0.2, color='gold',label='1500m isobath', transform=PROJN);
            ax2.add_feature(cfeature.LAND, zorder=3, facecolor='darkgray',edgecolor='white')

            ax2 = fig.add_subplot(2,4,3, projection=PROJ) # SHFLUX
            cts=ax2.pcolormesh(lon_rho_4km,lat_rho_4km,sum_shf,
                           vmin=-100, vmax=100 ,cmap=cmocean.cm.curl, transform=PROJN)
            cma=ax2.pcolormesh(lon_rho_4km,lat_rho_4km,
                           ma.masked_where(contour_masked_above_CF_4km == -1000, sum_melt),
                           transform=PROJN, vmin=-8, vmax=8, cmap='bwr') #, cmap=cmocean.cm.balance
            ax2.quiver(lon_rho_4km[::5,::5],lat_rho_4km[::5,::5],sum_U[slayer[dd],::5,::5],sum_V[slayer[dd],::5,::5],scale=5,transform=PROJV)
            ax2.set_title('Sfc heat flux anomaly \n Summer')
            # ax2.set_extent([-135, -75, -76, -69], crs=ccrs.PlateCarree()) # PIG-Thwaites region
            ax2.set_extent([LonMin, LonMax, LatMin, LatMax], crs=PROJN) # Getz
            x_left, x_right = ax2.get_xlim()
            y_low, y_high = ax2.get_ylim()
            ratio = .8
            ax2.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
            gl=ax2.gridlines(draw_labels=True,linewidth=.1)
            gl.top_labels=False   # suppress top labels
            gl.right_labels=False   # suppress right labels
            gl.bottom_labels=False
            plt.scatter(lon_along_4km,lat_along_4km, s=.8, alpha=0.2, color='gold',label='1500m isobath', transform=PROJN);
            ax2.add_feature(cfeature.LAND, zorder=3, facecolor='darkgray',edgecolor='white')

            ax2 = fig.add_subplot(2,4,6, projection=PROJ) # july temp anom
            cta=ax2.pcolormesh(lon_rho_4km,lat_rho_4km,win_Ta.isel(s_rho=slayer[dd]),
                           vmin=-anomlim[dd], vmax=anomlim[dd] ,cmap='twilight_shifted', transform=PROJN)
            cma=ax2.pcolormesh(lon_rho_4km,lat_rho_4km,
                           ma.masked_where(contour_masked_above_CF_4km == -1000, win_melt),
                           transform=PROJN, vmin=-8, vmax=8, cmap='bwr') #, cmap=cmocean.cm.balance
            ax2.quiver(lon_rho_4km[::5,::5],lat_rho_4km[::5,::5],
                       win_U[slayer[dd],::5,::5],
                       win_V[slayer[dd],::5,::5],scale=5,transform=PROJV)
            ax2.set_title('Temperature anomaly \n Winter')
            # ax2.set_extent([-135, -75, -76, -69], crs=ccrs.PlateCarree()) # PIG-Thwaites region
            ax2.set_extent([LonMin, LonMax, LatMin, LatMax], crs=PROJN) # Getz
            x_left, x_right = ax2.get_xlim()
            y_low, y_high = ax2.get_ylim()
            ratio = .8
            ax2.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
            gl=ax2.gridlines(draw_labels=True,linewidth=.1)
            gl.top_labels=False   # suppress top labels
            gl.right_labels=False   # suppress right labels
            gl.bottom_labels=False
            plt.scatter(lon_along_4km,lat_along_4km, s=.8, alpha=0.2, color='gold',label='1500m isobath', transform=PROJN);
            ax2.add_feature(cfeature.LAND, zorder=3, facecolor='darkgray',edgecolor='white')

            ax2 = fig.add_subplot(2,4,7, projection=PROJ) # SHFLX
            cts=ax2.pcolormesh(lon_rho_4km,lat_rho_4km,win_shf,
                           vmin=-100, vmax=100 ,cmap=cmocean.cm.curl, transform=PROJN)
            cma=ax2.pcolormesh(lon_rho_4km,lat_rho_4km,
                           ma.masked_where(contour_masked_above_CF_4km == -1000, win_melt),
                           transform=PROJN, vmin=-8, vmax=8, cmap='bwr') #, cmap=cmocean.cm.balance
            ax2.quiver(lon_rho_4km[::5,::5],lat_rho_4km[::5,::5],
                       win_U[slayer[dd],::5,::5],
                       win_V[slayer[dd],::5,::5],scale=5,transform=PROJV)
            ax2.set_title('Sfc heat flux anomaly \n Winter')
            ax2.set_extent([LonMin, LonMax, LatMin, LatMax], crs=PROJN) # Getz
            x_left, x_right = ax2.get_xlim()
            y_low, y_high = ax2.get_ylim()
            ratio = .8
            ax2.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
            gl=ax2.gridlines(draw_labels=True,linewidth=.1)
            gl.top_labels=False   # suppress top labels
            gl.right_labels=False   # suppress right labels
            gl.bottom_labels=False
            plt.scatter(lon_along_4km,lat_along_4km, s=.8, alpha=0.2, color='gold',label='1500m isobath', transform=PROJN);
            ax2.add_feature(cfeature.LAND, zorder=3, facecolor='darkgray',edgecolor='white')


            ax3 = fig.add_subplot(2,4,4)
            cd = ax3.scatter(filt_salt_ann, filt_temp_ann, c='lightgray', s=(22/fig.dpi)**2)
            c = ax3.scatter(filt_salt_sum, filt_temp_sum, c=-filt_z_rho_sum, s=(22/fig.dpi)**2,
                                    marker="o", cmap=cmocean.cm.deep, vmin=0, vmax=900)
            ax3.set_title('TS-diagram - Summer')
            ax3.set_xlim([34,35])
            ax3.set_ylim([-2.5,.5])
            # ax3.legend()
            # ax3.set_ylabel('Temperature ($^\circ$C)')
            ax3.set_aspect(0.32)

            ax3 = fig.add_subplot(2,4,8)
            cd = ax3.scatter(filt_salt_ann, filt_temp_ann, c='lightgray', s=(22/fig.dpi)**2)
            c = ax3.scatter(filt_salt_win, filt_temp_win, c=-filt_z_rho_win, s=(22/fig.dpi)**2,
                                    marker="o", cmap=cmocean.cm.deep, vmin=0, vmax=900)
            ax3.set_title('TS-diagram - Winter')
            ax3.set_xlim([34,35])
            ax3.set_ylim([-2.5,.5])
            # ax3.legend()
            # ax3.set_ylabel('Temperature ($^\circ$C)')
            ax3.set_aspect(0.32)


            cbar_ax1 = fig.add_axes([0.07, 0.42,  0.006, 0.18])
            fig.colorbar(ct, cax=cbar_ax1, orientation='vertical')
            cbar_ax1.set_ylabel('Temperature ($^{\circ}$C)')#, labelpad=-35)

            cbar_ax12 = fig.add_axes([0.13, 0.365,  0.16, 0.01])
            fig.colorbar(cm, cax=cbar_ax12, orientation='horizontal')
            cbar_ax12.set_xlabel('Melt rate (m/yr)')#, labelpad=-35)

            cbar_ax2 = fig.add_axes([0.33, 0.525,  0.16, 0.01])
            fig.colorbar(cta, cax=cbar_ax2, orientation='horizontal')
            cbar_ax2.set_xlabel('Temperature anomaly ($^{\circ}$C)')#, labelpad=-35)

            cbar_ax3 = fig.add_axes([0.53, 0.525,  0.16, 0.01])
            fig.colorbar(cts, cax=cbar_ax3, orientation='horizontal')
            cbar_ax3.set_xlabel('Surface heat flux (W m$^{-2}$)')#, labelpad=-35)

            #cbar_ax1 = fig.add_axes([0.912, 0.2,  0.01, 0.6])
            #fig.colorbar(c, cax=cbar_ax1, orientation='vertical')
            #cbar_ax1.set_ylabel('Temperature-Salinity diagram; Depth (m)')#, labelpad=-35)
            cbar_ax1 = fig.add_axes([0.735, 0.525,  0.16, 0.01])
            fig.colorbar(c, cax=cbar_ax1, orientation='horizontal')
            cbar_ax1.set_xlabel('Temperature-Salinity diagram; Depth (m)')

            name_fig='waom4_MapsT+SHF_season_' + pname[dd] + '_' + Pname +'_newTS.png'
            plt.savefig(fig_path + name_fig, dpi=600, bbox_inches='tight')
