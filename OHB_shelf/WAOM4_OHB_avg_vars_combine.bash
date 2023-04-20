# This script concatenate the variables from ocean_avg_00*.nc files that will be used
# for the cross-shelf heat transport calculations:

# from WAOM10extend_shflim_S_0.25Q_OHB_shelf.ipynb:

#    temp_tmp = ds.variables["temp"] # X,31,560,630
#    salt_tmp = ds.variables["salt"]
#    shflux_tmp = ds.variables["shflux"] # X,560,630
#    ssflux_tmp = ds.variables["ssflux"]
#    m_tmp = ds.variables["m"]
#    HvomT_tmp = ds.variables["Hvom_temp"] # X,31,559,630
#    HuonT_tmp = ds.variables["Huon_temp"] # X,31,560,629
#    Hvom_tmp = ds.variables["Hvom"] # X,31,559,630
#    Huon_tmp = ds.variables["Huon"] # X,31,560,629


# using nco:
module load nco

# waom4extend_shflim_S_0.25Q outputs:
cd /g/data3/hh5/tmp/access-om/fbd581/ROMS/OUTPUT/waom4extend_shflim_S_0.25Q/output_yr10_diag/

for n in 01 02 03 04 05 06 07 08 09 10 11 12;
do
   echo ${n}
#   ncks -v temp,salt,shflux,ssflux,Hvom,Huon ocean_avg_00${n}.nc ocean_avg_00${n}_sel.nc
done

ncks -A ocean_avg_0001_sel.nc ocean_avg_0002_sel.nc
ncks -A ocean_avg_0002_sel.nc ocean_avg_0003_sel.nc
ncks -A ocean_avg_0003_sel.nc ocean_avg_0004_sel.nc
ncks -A ocean_avg_0004_sel.nc ocean_avg_0005_sel.nc
ncks -A ocean_avg_0005_sel.nc ocean_avg_0006_sel.nc
ncks -A ocean_avg_0006_sel.nc ocean_avg_0007_sel.nc
ncks -A ocean_avg_0007_sel.nc ocean_avg_0008_sel.nc
ncks -A ocean_avg_0008_sel.nc ocean_avg_0009_sel.nc
ncks -A ocean_avg_0009_sel.nc ocean_avg_0010_sel.nc
ncks -A ocean_avg_0010_sel.nc ocean_avg_0011_sel.nc
ncks -A ocean_avg_0011_sel.nc ocean_avg_0012_sel.nc
