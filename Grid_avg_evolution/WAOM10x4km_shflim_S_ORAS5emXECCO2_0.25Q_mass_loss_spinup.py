import os
import sys
import xarray as xr
import dask
import matplotlib.pyplot as plt
import numpy as np
import matplotlib

roms_dir = os.path.join('/g','data3','hh5','tmp','access-om','fbd581', 'ROMS')
proj_dir = os.path.join(roms_dir,'postprocessing','figs','Grid_avg')
data_dir = os.path.join(roms_dir,'OUTPUT')

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

grd10_path = os.path.join(roms_dir,'waom10_frc','waom10_grd.nc')
waom10_ORAS5_path = os.path.join(data_dir,'waom10_shflim_S_ORAS5em_0.25Q','output_01-20yr','ocean_avg_00??.nc')
waom10_ECCO2_path = os.path.join(data_dir,'waom10_shflim_S_0.25Q','output_01-20yr','ocean_avg_00??.nc')
vostock10 = [np.arange(200,260),np.arange(400,480)]
time = slice('2007','2026')

avg10_ORAS5 = xr.open_mfdataset(waom10_ORAS5_path)
avg10_ECCO2 = xr.open_mfdataset(waom10_ECCO2_path)
grd10 = xr.open_dataset(grd10_path)

bmb10_ORAS5,ismr10_ORAS5,mask10_ORAS5 = total_bmb(avg10_ORAS5,grd10,vostock10,time)
bmb10_ECCO2,ismr10_ECCO2,mask10_ECCO2 = total_bmb(avg10_ECCO2,grd10,vostock10,time)


out_path = os.path.join(proj_dir,'mass_loss_spinup.png')
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
rcParams.update({'font.size': 18})

plt.close()
fig,ax = plt.subplots(nrows=2,figsize=(15,12))
ismr10_ORAS5.plot(label='10 km ORAS5',ax=ax[0])
ismr10_ECCO2.plot(label='10 km ECCO2',ax=ax[0])
#ismr4.plot(label='4 km',ax=ax)
#ismr2.plot(label='2 km',ax=ax
ax[0].set_ylim(0.6,2.2)
ax[0].set_title('Area average melt rate',fontsize=24)
ax[0].set_ylabel('(m/yr)',fontsize=20)
ax[0].set_xlabel('',fontsize=20)
ax[0].legend()
ax[0].grid()
#ax[0].set_xticklabels(range(20))
ax[0].legend(markerscale=1,fontsize=20)

bmb10_ORAS5.plot(label='10 km', ax=ax[1])
bmb10_ECCO2.plot(label='10 km', ax=ax[1])
ax[1].set_ylim(1000,4000)
#ax[1].set_title('Area average melt rate (Gt/yr)',fontsize=24)
ax[1].set_ylabel('(Gt/yr)',fontsize=20)
ax[1].grid()
ax[1].set_xlabel('Years after initialization',fontsize=20)
#ax[1].set_xticklabels(range(20))
plt.text(2008,3700,'ORAS5 = ' + '{0:.2f}'.format(np.nanmean(bmb10_ORAS5[-12:-1])))
plt.text(2008,3300,'ECCO2 = ' + '{0:.2f}'.format(np.nanmean(bmb10_ECCO2[-12:-1])))
plt.savefig(out_path,dpi=300)

plt.show()


