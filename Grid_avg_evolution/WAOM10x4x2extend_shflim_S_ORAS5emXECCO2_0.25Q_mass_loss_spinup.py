import os
import sys
import xarray as xr
import dask
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import matplotlib

roms_dir = os.path.join('/users','boeiradi','COLD_project')
proj_dir = os.path.join(roms_dir,'postprocessing','figs','Melt_rates')

#from dask.distributed import Client
#C = Client()
#C

def total_bmb(avg,grd,vostock,time_slice):

    s2a = 3600*24*365.25
    rhoi = 916

    mask = (grd.mask_rho == 1 )& (grd.zice < 0)
    mask[vostock[0],vostock[1]] = False

    dA = (1/(grd.pm*grd.pn)).where(mask)
    weights = dA/dA.sum()

    ismr = (avg.m.sel(ocean_time=time_slice).where(mask)*weights*s2a).sum(['xi_rho','eta_rho'])
    bmb = (avg.m.sel(ocean_time=time_slice).where(mask)*dA*rhoi*(10**-12)*s2a).sum(['xi_rho','eta_rho'])

    return bmb,ismr,mask

grd10_path = os.path.join('/scratch','project_2000339','boeiradi','waom10_frc','waom10extend_grd.nc')
waom10_path = os.path.join('/scratch','project_2000789','boeiradi','waom10extend_shflim_S_0.25Q','output_01-20yr','ocean_avg_00??.nc')
grd4_path = os.path.join('/scratch','project_2000339','boeiradi','waom4_frc','waom4extend_grd.nc')
waom4_path = os.path.join('/scratch','project_2000789','boeiradi','waom4extend_shflim_S_0.25Q','output_01-10yr','ocean_avg_00??.nc')
grd2_path = os.path.join('/scratch','project_2000339','boeiradi','waom2_frc','waom2extend_grd.nc')
waom2_path = os.path.join('/scratch','project_2000339','boeiradi','waom2extend_shflim_S_0.25Q','output_01-05yr','ocean_avg_00??.nc')

vostock10 = [np.arange(200,260),np.arange(400,480)]
vostock4 = [np.arange(500,650),np.arange(1000,1200)]
vostock2 = [np.arange(1000,1300),np.arange(2000,2400)]

time_10km = slice('2007','2026')
time_4km = slice('2007','2016')
time_2km = slice('2007','2011')

avg10 = xr.open_mfdataset(waom10_path)
print('WAOM10: ', avg10)
grd10 = xr.open_dataset(grd10_path)
bmb10,ismr10,mask10 = total_bmb(avg10,grd10,vostock10,time_10km)
avg10.close()
grd10.close()

avg4 = xr.open_mfdataset(waom4_path)
print('WAOM4: ', avg4)
grd4 = xr.open_dataset(grd4_path)
bmb4,ismr4,mask4 = total_bmb(avg4,grd4,vostock4,time_4km)
avg4.close()
grd4.close()

avg2 = xr.open_mfdataset(waom2_path)
print('WAOM2: ', avg2)
grd2 = xr.open_dataset(grd2_path)
bmb2,ismr2,mask2 = total_bmb(avg2,grd2,vostock2,time_2km)
avg2.close()
grd2.close()

out_path = os.path.join(proj_dir,'waom10x4x2extend_mass_loss_spinup.png')
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
rcParams.update({'font.size': 18})
plt.close()

fig,ax = plt.subplots(nrows=2,figsize=(15,12))
ismr10.plot(label='10 km',ax=ax[0])
ismr4.plot(label='4 km',ax=ax[0])
ismr2.plot(label='2 km',ax=ax[0])
ax[0].set_ylim(0.0,2.2)
ax[0].set_title('Area average melt rate',fontsize=24)
ax[0].set_ylabel('(m/yr)',fontsize=20)
ax[0].set_xlabel('',fontsize=20)
ax[0].legend()
ax[0].grid()
#ax[0].set_xticklabels(range(20))
ax[0].legend(markerscale=1,fontsize=20)

bmb10.plot(label='10 km', ax=ax[1])
bmb4.plot(label='4 km', ax=ax[1])
bmb2.plot(label='2 km', ax=ax[1])
ax[1].set_ylim(500,4000)
#ax[1].set_title('Area average melt rate (Gt/yr)',fontsize=24)
ax[1].set_ylabel('(Gt/yr)',fontsize=20)
ax[1].grid()
ax[1].set_xlabel('Years after initialization',fontsize=20)
#ax[1].set_xticklabels(range(20))
plt.text(2008,3700,'10 km = ' + '{0:.2f}'.format(np.nanmean(bmb10[-12:-1])))
plt.text(2008,3300,'4 km = ' + '{0:.2f}'.format(np.nanmean(bmb4[-12:-1])))
plt.text(2008,2900,'2 km = ' + '{0:.2f}'.format(np.nanmean(bmb2[-12:-1])))
plt.savefig(out_path,dpi=300)
#plt.show()
plt.close()

# plot using my setup:
fig = plt.figure(figsize=(12,8))
ax1 = fig.add_subplot(211)
plt.plot(np.arange(0,20)+1,ismr10,'--r',label='10km')
plt.plot(np.arange(10,20)+1,ismr4,'--b',label='4km')
plt.plot(np.arange(15,20)+1,ismr2,'--g',label='2km')
#l3 = plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)
l3 = plt.legend(loc="upper left", borderaxespad=0)
plt.title('Area average melt rate')
plt.ylabel('(m/yr)')
plt.ylim(0,2.2)
ax2 = fig.add_subplot(212)
plt.plot(np.arange(0,20)+1,bmb10,'r',label='10km')
plt.plot(np.arange(10,20)+1,bmb4,'--b',label='4km')
plt.plot(np.arange(15,20)+1,bmb2,'--g',label='2km')
#l3 = plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)
plt.ylabel('(Gt/yr)')
plt.ylim(500,4000)

out_path = os.path.join(proj_dir,'waom10x4x2extend_mass_loss_spinup2.png')
plt.savefig(fig_path + name_fig, dpi=300)
plt.close()
