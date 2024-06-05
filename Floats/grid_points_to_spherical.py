
import numpy as np
#import cartopy.crs as ccrs
#import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib
import time

from matplotlib.axes import Axes
#from cartopy.mpl.geoaxes import GeoAxes
#GeoAxes._pcolormesh_patched = Axes.pcolormesh
from scipy.interpolate import interpn

#from mpl_toolkits.basemap import Basemap
import math
import warnings
#%matplotlib inline
import xarray as xr
warnings.filterwarnings('ignore')
warnings.filterwarnings(action='ignore', message='All-NaN slice encountered')

from matplotlib import cm
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
ds = xr.open_dataset('/scratch/project_2000789/boeiradi/waom4extend_shflim_S_0.25Q/output_yr10a_floats/ocean_flt.nc')
dg = xr.open_dataset('/scratch/project_2000789/boeiradi/waom4_frc/waom4extend_grd.nc')
# %load latlon_functions.py
def interpolate(xarr,yarr,ds):
    """
    inputs:
        ds: the dataset with lat and long values
        xarr: the current row of the Xgrid array
        yarr: the current row of the Ygrid array
    returns:
        points: the lat/lon values of the current array
    """
    lat_values = ds.lat_rho.values
    long_values = ds.lon_rho.values
    lat = np.array([])
    long = np.array([])
    
    x_shape1 = np.arange(ds.lat_rho.shape[0])
    y_shape1 = np.arange(ds.lat_rho.shape[1])
    x_shape2 = np.arange(ds.lon_rho.shape[0])
    y_shape2 = np.arange(ds.lon_rho.shape[1])
    
    for i,j in zip(xarr,yarr):
        #print(len(i))
        interp_x=i[~np.isnan(i)]
        interp_y = j[~np.isnan(j)]
        
        if len(interp_x)>0 and len(interp_y)>0:
            # latitude points
            interp_mesh = np.array(np.meshgrid(interp_x, interp_y))
            interp_points = np.rollaxis(interp_mesh, 0, 3).reshape((1, 2))
            # Perform the interpolation
            interp_arr1 = interpn((x_shape1, y_shape1), lat_values, interp_points)
            
            # true values with interp_arr1 values
            #np.maximum.accumulate(idx, out=idx)
            #print(out)
            
            lat = np.append(lat,interp_arr1[0])

            # interpolate in longitude
            # Note the following two lines that are used to set up the
            interp_x = j[~np.isnan(j)]
            interp_y = i[~np.isnan(i)]
            interp_mesh = np.array(np.meshgrid(interp_x, interp_y))
            interp_points = np.rollaxis(interp_mesh, 0, 3).reshape((1, 2))
            # Perform the interpolation
            interp_arr2 = interpn((x_shape2, y_shape2), long_values, interp_points)
            long = np.append(long,interp_arr2[0])
        else:
            lat = np.append(lat,np.nan)
            long = np.append(long,np.nan)
    
    return (long,lat)


    # write a function to get the long and lat from the grid points
def grid_to_spherical(x,y,ds):
    """
    x: array of grid x position values
    y: array of grid y position values
    rho: the used xi_rho value 
    """
    start_time = time.time()
    x_long = np.empty((0,x.shape[1]))
    y_lat = np.empty((0,y.shape[1]))
    
    # first create the points in the other array
    # get the maximum and minimum values of x and y
    for i in range(len(x)):
        if (i%500)==0:
            print(i,'/',len(x),'lines at time:',(time.time()-start_time)/60,'min')
        # get the min and max of x and y
        minx = (np.nanmin(x[i]))
        maxx = (np.nanmax(x[i]))
        if not (math.isnan(minx)) and not (math.isnan(maxx)):
            #x[i]=np.nan_to_num(x[i], nan=0.0)
            #y[i]=np.nan_to_num(y[i], nan=0.0)
            x_i,y_i = interpolate(x[i],y[i],ds)
        else:
            x_i,y_i = x[i],y[i]
        x_long = np.append(x_long,np.array([x_i]),axis=0)
        y_lat = np.append(y_lat,np.array([y_i]),axis=0)
        
    #x_long = np.fliplr(x_long)
    #y_lat = np.fliplr(y_lat)
    return (x_long,y_lat)

# x = (ds.variables['Xgrid'][5800:5810].values)
# y = (ds.variables['Ygrid'][5800:5810].values)
# #print(x)
# long,lat = grid_to_spherical(x,y,dg)
# print(x[0][36])
# print(long)
# create the latitude and longitude points from the spherical grids
x = (ds.variables['Xgrid'][5761:].values)
y = (ds.variables['Ygrid'][5761:].values)

#firstlong = np.array([(row[20]) for s,row in enumerate(x) if not math.isnan(row[20])])
long,lat = grid_to_spherical(x,y,dg)

# write the new latitude and longitude points to a file
lines = []
for i in (long):
    line = ""
    for x in i:
        mystr = str(x)
        line += mystr
        line += ","
    # create the line to add to lines
    lines.append(line)
with open('lonpoints531.txt', 'w') as f:
    for line in lines:
        f.write(line)
        f.write('\n')
f.close()

lines=[]
for i in (lat):
    line = ""
    for x in i:
        mystr = str(x)
        line += mystr
        line += ","
    # create the line to add to lines
    lines.append(line)
with open('latpoints531.txt', 'w') as f:
    for line in lines:
        f.write(line)
        f.write('\n')
f.close()
# TESTING
ds1 = xr.open_dataset('/scratch/project_2000789/muramarg/run_5_31/output_WAOM_check/ocean_avg_0003.nc')
ds1 = ds1.chunk({'ocean_time': 1})
ds1.coords["lon_rho"] = dg.lon_rho
ds1.coords["lat_rho"] = dg.lat_rho
x = (ds.variables['Xgrid'].values)
y = (ds.variables['Ygrid'].values)

proj = ccrs.SouthPolarStereo(central_longitude=0.0, true_scale_latitude=None, globe=None)
fig = plt.figure(figsize=(20,10))
ax = fig.add_subplot(1, 1, 1, projection=proj)
ax.coastlines(zorder=7,facecolor='white',edgecolor='black')
#plt.pcolormesh(dg.lon_rho[300:560,200:630],dg.lat_rho[300:560,200:630],ds1.temp.isel(xi_rho=slice(200, 630),eta_rho=slice(300,560)),transform=ccrs.PlateCarree())

xlimit = np.arange(300,560,1)
ylimit = np.arange(350,630,1)

#myds = ds1.temp[xlimit,ylimit]
ds1.temp.isel(s_rho=-1, ocean_time=0,xi_rho=slice(350, 630),eta_rho=slice(300,630)).plot(x="lon_rho", y="lat_rho",transform=ccrs.PlateCarree())
plt.contour(dg.lon_rho[xlimit,ylimit], dg.lat_rho[xlimit,ylimit],ds1.zice[xlimit,ylimit],levels=[-0.2],zorder=5,linestyles='dashed', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
plt.contour(dg.lon_rho[xlimit,ylimit], dg.lat_rho[xlimit,ylimit],ds1.h[xlimit,ylimit],levels=[200,400,600],zorder=5, linestyles='solid', transform=ccrs.PlateCarree(), cmap=plt.cm.binary)
bathym = cfeature.NaturalEarthFeature(name='bathymetry_J_1000', scale='10m', category='physical')
ax.add_feature(bathym, zorder=3, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=1)
ax.add_feature(cfeature.LAND,facecolor='#c9c9c9',zorder=6)
plt.scatter(long,lat,marker = '.', linewidths=0.02,color='g',alpha=0.1,transform=ccrs.PlateCarree())

plt.savefig("point_contours.png")
plt.show()
 
