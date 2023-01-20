# read nc output from WAOM 10km run

import xarray as xr
import pandas as p
import numpy as np
import numpy.ma as ma
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.dates as dates
from matplotlib.ticker import NullFormatter
from matplotlib.colors import LinearSegmentedColormap   # for custom colormaps

from datetime import datetime, timedelta

from netCDF4 import Dataset
from netCDF4 import num2date, date2num
from matplotlib.ticker import NullFormatter
from matplotlib.colors import LinearSegmentedColormap   # for custom colormaps

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import gsw

#/scratch/gh9/fbd581/ROMS/waom4extend_shflim_S_0.25Q/output_yr10_notides_diag/

# load ROMS avg output: 4km
for mm  in ['01','02','03','04','05','06','07','08','09','10','11','12']:
    ds = xr.open_dataset('/scratch/gh9/fbd581/ROMS/waom4extend_shflim_S_0.25Q/output_yr10_notides_diag/ocean_avg_00' + mm + '.nc')
    print(ds.variables["temp"].shape)
    m_tmp = np.nanmean(ds.variables["m"], axis=0)

    if mm == '01':
        m_4km = m_tmp

    elif mm == '02':
        m_4km = np.stack((m_4km,m_tmp), axis=0)

    else:
        m_tmp_4thdim = np.expand_dims(m_tmp, axis=0)
        m_4km = np.concatenate((m_4km,m_tmp_4thdim), axis=0)

    ds.close()

dg = xr.open_dataset("/g/data3/hh5/tmp/access-om/fbd581/ROMS/waom4_frc/waom4extend_grd.nc")

lat_rho_4km = dg.variables["lat_rho"]
lon_rho_4km = dg.variables["lon_rho"]

ds.coords['lat_rho']=lat_rho_4km.transpose() # put lat_rho into ds dataset
ds.coords['lon_rho']=lon_rho_4km.transpose() # put lon_rho into ds dataset

### plot some maps
import matplotlib.path as mpath
import cartopy.feature as cfeature
import matplotlib.colors as colors

def lonlat_labels(ax):
    # latitude labels
    ax.text(120,-80,'80$^{\circ}$S',transform=ccrs.PlateCarree(),color='gray')
    ax.text(120,-70,'70$^{\circ}$S',transform=ccrs.PlateCarree(),color='gray')
    # longitude labels
    ax.text(0,-66,'0$^{\circ}$',transform=ccrs.PlateCarree(),color='gray')
    #ax.text(60,-53,'60$^{\circ}$E',transform=ccrs.PlateCarree(),color='gray')
    #ax.text(120,-53,'120$^{\circ}$E',transform=ccrs.PlateCarree(),color='gray')
    ax.text(-60,-48,'60$^{\circ}$W',transform=ccrs.PlateCarree(),color='gray')
    ax.text(-120,-48,'120$^{\circ}$W',transform=ccrs.PlateCarree(),color='gray')
    ax.text(180,-60,'180$^{\circ}$',transform=ccrs.PlateCarree(),color='gray')
    return

proj = ccrs.SouthPolarStereo()

# plotting surface fluxes:
fig_path = '/g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/figs/Melt_rates/'

# define your scale, with white at zero
vmin = -2
vmax = 6
norm = colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

proj = ccrs.SouthPolarStereo()
fig = plt.figure(figsize=(16,10))

ax2 = fig.add_subplot(222, projection=proj)
plt.title('Annual melt rates (m/yr): 4km no tides')
cy=plt.pcolormesh(lon_rho_4km,lat_rho_4km,np.nanmean(m_4km, axis=0)*86400*365, transform=ccrs.PlateCarree(), cmap=plt.cm.bwr, vmin=vmin, vmax=vmax, norm=norm)
plt.colorbar(cy, extend='both')
ax2.gridlines() # draw_labels=True,linewidth=2, color='white', alpha=0.5, linestyle='--')
lonlat_labels(ax2)
ax2.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
ax2.add_feature(cfeature.COASTLINE, zorder=4, edgecolor='black', alpha=0.2)
plt.clim(-2,6)

name_fig="waom4extend_shflim_S_0.25Q_notides_melt_rates_maps_annual_yr20.png"
plt.savefig(fig_path + name_fig, dpi=300)

fig = plt.figure(figsize=(12,5))
 
ax2 = fig.add_subplot(132, projection=proj)
plt.title('Annual melt rates (m/yr): 4km no tides')
cy=plt.pcolormesh(lon_rho_4km,lat_rho_4km,np.nanmean(m_4km, axis=0)*86400*365, transform=ccrs.PlateCarree(), cmap=plt.cm.bwr, vmin=vmin, vmax=vmax, norm=norm) 
ax2.gridlines() # draw_labels=True,linewidth=2, color='white', alpha=0.5, linestyle='--')
#lonlat_labels(ax2)
ax2.set_extent([-80, 10, -80, -65], ccrs.PlateCarree())
ax2.add_feature(cfeature.COASTLINE, zorder=4, edgecolor='black')#, alpha=0.2) 
plt.clim(-2,6)

axins = inset_axes(ax2,
                   width="5%",  # width = 5% of parent_bbox width
                   height="100%",  # height : 50%
                   loc='lower left',
                   bbox_to_anchor=(1.05, 0., 1, 1),
                   bbox_transform=ax2.transAxes,
                   borderpad=0,
                   )
fig.colorbar(cy, cax=axins, orientation="vertical")
                                                
name_fig="waom4extend_shflim_S_0.25Q_notides_melt_rates_maps_annual_yr20_RFIS.png"
plt.savefig(fig_path + name_fig, dpi=300)
#plt.close()

fig = plt.figure(figsize=(12,5))

ax2 = fig.add_subplot(132, projection=proj)
plt.title('Annual melt rates (m/yr): 4km no tides')
cy=plt.pcolormesh(lon_rho_4km,lat_rho_4km,np.nanmean(m_4km, axis=0)*86400*365, transform=ccrs.PlateCarree(), cmap=plt.cm.bwr, vmin=vmin, vmax=vmax, norm=norm)
ax2.gridlines() # draw_labels=True,linewidth=2, color='white', alpha=0.5, linestyle='--')
#lonlat_labels(ax2)
ax2.set_extent([130, 280, -85, -70], ccrs.PlateCarree())
ax2.add_feature(cfeature.COASTLINE, zorder=4, edgecolor='black')#, alpha=0.2)
plt.clim(-2,6)
axins = inset_axes(ax2,
                   width="5%",  # width = 5% of parent_bbox width
                   height="100%",  # height : 50%
                   loc='lower left',
                   bbox_to_anchor=(1.05, 0., 1, 1),
                   bbox_transform=ax2.transAxes,
                   borderpad=0,
                   )
fig.colorbar(cy, cax=axins, orientation="vertical")

name_fig="waom4extend_shflim_S_0.25Q_notides_melt_rates_maps_annual_yr20_RIS.png"
plt.savefig(fig_path + name_fig, dpi=300)
#plt.close()
