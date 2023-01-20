import xarray as xr
import pandas as p
import numpy as np
import numpy.ma as ma
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.dates as dates
from matplotlib.ticker import NullFormatter
from matplotlib.colors import LinearSegmentedColormap   # for custom colormaps
import matplotlib as mpl
mpl.use('Agg')

from datetime import datetime, timedelta

import os
from pyresample import kd_tree, geometry, utils # for interpolation
from scipy import interpolate

import netCDF4
from matplotlib.ticker import NullFormatter
from matplotlib.colors import LinearSegmentedColormap   # for custom colormaps

data_dir = os.path.join('/scratch','project_2000339','boeiradi','waom2extend_shflim_S_0.25Q','output_yr5')
out_path = os.path.join(data_dir,'waom2extend_to10km_grd_monthly_avg.nc')

print(out_path)

def high_to_low(hr_da,hr_grd,lr_grd,gt,dim,fill_value=0.0):
    
    print('set up empty lr data array')
    if dim == 2:
    
        dummy = np.zeros(lr_grd['lon_'+gt].shape)
        x = lr_grd['xi_'+gt]
        y = lr_grd['eta_'+gt]
        lr_da = xr.DataArray(dummy,coords=[y,x],dims=['eta_'+gt,'xi_'+gt])
        
    elif dim == 3:
        
        N = hr_da.s_rho.size
        dummy = np.tile(np.zeros(lr_grd['lon_'+gt].shape),(N,1,1))
        x = lr_grd['xi_'+gt]
        y = lr_grd['eta_'+gt]
        z = hr_da['s_rho']
        lr_da = xr.DataArray(dummy,coords=[z,y,x],dims=['s_rho','eta_'+gt,'xi_'+gt])
    
    
    # Fill the mask of low resolution data with nearest neibghours and fill in known values on high res grid.
    if dim == 2:
        
        print('Fill in the mask of lr data')
        data = hr_da.values

        valid_mask = ~np.isnan(data)
        coords = np.array(np.nonzero(valid_mask)).T
        values = data[valid_mask]

        it = interpolate.NearestNDInterpolator(coords,values)

        filled = it(list(np.ndindex(data.shape))).reshape(data.shape)
        
        print('Resample to high resolution grid')
        lr_da[:,:] = kd_tree.get_sample_from_neighbour_info('custom', lr_def[gt].shape, filled,\
                                             valid_input_index[gt],\
                                             valid_output_index[gt],index_array[gt],distance_array[gt],wf)
        
        # Fill with zeros where mask is present
        #print('fill hr mask areas with fill value: ',fill_value)
        #hr_da.values[hr_grd['mask_'+gt].values == 0] = fill_value
            
    if dim == 3:
        
        print('Fill in the mask of lr data and resample to high resolution grid')
        for k in np.arange(N):
            
            print('processing depth level: ',k)
            data = hr_da[k].values

            valid_mask = ~np.isnan(data)
            coords = np.array(np.nonzero(valid_mask)).T
            values = data[valid_mask]

            it = interpolate.NearestNDInterpolator(coords,values)

            filled = it(list(np.ndindex(data.shape))).reshape(data.shape)
    
            # Fill in known values on high res grid
            lr_da[k] = kd_tree.get_sample_from_neighbour_info('custom', lr_def[gt].shape,filled,\
                                             valid_input_index[gt],\
                                             valid_output_index[gt],index_array[gt],distance_array[gt],wf)
            
            # Fill with zeros where mask is present
            #print('fill hr mask areas with fill value: ',fill_value)
            #hr_da[k].values[hr_grd['mask_'+gt].values == 0] = fill_value
            
    return lr_da

def interp_roms_2to10km(): # 'wed',0,0

    # load ROMS avg output
    count = 0
    for mm in ['01','02','03','04','05','06','07','08','09','10','11','12']:
        hr_ds = xr.open_dataset("/scratch/project_2000339/boeiradi/waom2extend_shflim_S_0.25Q/output_yr5/ocean_avg_00" + mm + ".nc")
        print('size temp and time length: ', hr_ds.temp.shape, len(hr_ds.salt.isel(xi_rho=10, eta_rho=100, s_rho=0)))
        # ----- # first calculate z_rho spatially and put into hr_ds xarray
        hr_ds = hr_ds.set_coords(['Cs_r', 'Cs_w', 'hc', 'h', 'Vtransform'])
        Zo_rho = (hr_ds.hc * hr_ds.s_rho + hr_ds.Cs_r * hr_ds.h) / (hr_ds.hc + hr_ds.h)
        z_rho = hr_ds.zeta + (hr_ds.zeta + hr_ds.h) * Zo_rho
        #print("Vtransform=2, z_rho shape", z_rho.shape)
        hr_ds.coords['z_rho'] = z_rho.transpose('ocean_time', 's_rho', 'eta_rho', 'xi_rho')
        print('size z_rho:', hr_ds.z_rho.shape)
        
        for var,gt,dim in zip(['temp','salt','zeta','z_rho'],['rho','rho','rho','rho'],[3,3,2,3]):
        #for var,gt,dim in zip(['temp'],['rho'],[3]):
            if mm == '01':
                print('processing: ',var, mm)
                lr_jan[var] = high_to_low(hr_ds[var].mean(dim='ocean_time', skipna=True),hr_grd,lr_grd,gt,dim)
            elif mm == '02':
                print('processing: ',var, mm)
                lr_feb[var] = high_to_low(hr_ds[var].mean(dim='ocean_time', skipna=True),hr_grd,lr_grd,gt,dim)    
            elif mm == '03':
                print('processing: ',var, mm)
                lr_mar[var] = high_to_low(hr_ds[var].mean(dim='ocean_time', skipna=True),hr_grd,lr_grd,gt,dim)  
            elif mm == '04':
                print('processing: ',var, mm)
                lr_apr[var] = high_to_low(hr_ds[var].mean(dim='ocean_time', skipna=True),hr_grd,lr_grd,gt,dim)  
            elif mm == '05':
                print('processing: ',var, mm)
                lr_mai[var] = high_to_low(hr_ds[var].mean(dim='ocean_time', skipna=True),hr_grd,lr_grd,gt,dim)                  
            elif mm == '06':
                print('processing: ',var, mm)
                lr_jun[var] = high_to_low(hr_ds[var].mean(dim='ocean_time', skipna=True),hr_grd,lr_grd,gt,dim)  
            elif mm == '07':
                print('processing: ',var, mm)
                lr_jul[var] = high_to_low(hr_ds[var].mean(dim='ocean_time', skipna=True),hr_grd,lr_grd,gt,dim)  
            elif mm == '08':
                print('processing: ',var, mm)
                lr_aug[var] = high_to_low(hr_ds[var].mean(dim='ocean_time', skipna=True),hr_grd,lr_grd,gt,dim)  
            elif mm == '09':
                print('processing: ',var, mm)
                lr_sep[var] = high_to_low(hr_ds[var].mean(dim='ocean_time', skipna=True),hr_grd,lr_grd,gt,dim)  
            elif mm == '10':
                print('processing: ',var, mm)
                lr_oct[var] = high_to_low(hr_ds[var].mean(dim='ocean_time', skipna=True),hr_grd,lr_grd,gt,dim)  
            elif mm == '11':
                print('processing: ',var, mm)
                lr_nov[var] = high_to_low(hr_ds[var].mean(dim='ocean_time', skipna=True),hr_grd,lr_grd,gt,dim)  
            elif mm == '12':
                print('processing: ',var, mm)
                lr_dec[var] = high_to_low(hr_ds[var].mean(dim='ocean_time', skipna=True),hr_grd,lr_grd,gt,dim)                                                 
    lr_annual = xr.concat([lr_jan, lr_feb, lr_mar, lr_apr, lr_mai, lr_jun, lr_jul,lr_aug, lr_sep, lr_oct, lr_nov, lr_dec], dim='ocean_time')
    lr_annual = lr_annual.mean(dim='ocean_time', skipna=True)
    print('shapes lr_01==jan and lr_annual: ', lr_jan.temp.shape, lr_annual.temp.shape)
    #del lr_jan, lr_feb, lr_mar, lr_apr, lr_mai, lr_jun, lr_jul,lr_aug, lr_sep, lr_oct, lr_nov, lr_dec
    lr_annual.to_netcdf(out_path,unlimited_dims='ocean_time')
       
    return 

hr_grd = xr.open_dataset("/scratch/project_2000339/boeiradi/waom2_frc/waom2extend_grd.nc")
lr_grd = xr.open_dataset("/scratch/project_2000339/boeiradi/waom10_frc/waom10extend_grd.nc")
#lr_ini = xr.open_dataset("/scratch/project_2000339/boeiradi/waom10_frc/waom10_ini.nc")
lr_avg =  xr.open_dataset("/scratch/project_2000789/boeiradi/waom10extend_shflim_S_0.25Q/output_01-20yr/ocean_avg_0001.nc") #/scratch/project_2000339/boeiradi/waom10_shflim_S/output_11-15yr/ocean_avg_0003.nc")
    
lr_def = {}
hr_def= {}

valid_input_index = {}
valid_output_index = {}
index_array = {}
distance_array = {}

for gt in ['rho','u','v']:

    lr_def[gt] = geometry.SwathDefinition(lons=lr_grd['lon_'+gt].values,lats=lr_grd['lat_'+gt].values)
    hr_def[gt] = geometry.SwathDefinition(lons=hr_grd['lon_'+gt].values,lats=hr_grd['lat_'+gt].values)

    valid_input_index[gt], valid_output_index[gt], index_array[gt], distance_array[gt] = \
    kd_tree.get_neighbour_info(hr_def[gt],lr_def[gt], 20000,neighbours=4,nprocs=40)


wf = lambda r: 1/r

lr_jan = lr_avg.drop(['u','v','ubar','vbar','salt','temp','zeta','ocean_time'])
lr_feb = lr_avg.drop(['u','v','ubar','vbar','salt','temp','zeta','ocean_time'])
lr_mar = lr_avg.drop(['u','v','ubar','vbar','salt','temp','zeta','ocean_time'])
lr_apr = lr_avg.drop(['u','v','ubar','vbar','salt','temp','zeta','ocean_time'])
lr_mai = lr_avg.drop(['u','v','ubar','vbar','salt','temp','zeta','ocean_time'])
lr_jun = lr_avg.drop(['u','v','ubar','vbar','salt','temp','zeta','ocean_time'])
lr_jul = lr_avg.drop(['u','v','ubar','vbar','salt','temp','zeta','ocean_time'])
lr_aug = lr_avg.drop(['u','v','ubar','vbar','salt','temp','zeta','ocean_time'])
lr_sep = lr_avg.drop(['u','v','ubar','vbar','salt','temp','zeta','ocean_time'])
lr_oct = lr_avg.drop(['u','v','ubar','vbar','salt','temp','zeta','ocean_time'])
lr_nov = lr_avg.drop(['u','v','ubar','vbar','salt','temp','zeta','ocean_time'])
lr_dec = lr_avg.drop(['u','v','ubar','vbar','salt','temp','zeta','ocean_time'])

### ======== run

interp_roms_2to10km()
