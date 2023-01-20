#/usr/bin/bash

module load nco

#mv ocean_rst.nc ocean_rst_sedkg.nc
cp ocean_rst_mm10.nc ocean_rst_mm10_sedkg.nc

ncap2 -A -s 'sand_01=temp*0' ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc
ncap2 -A -s 'sand_02=temp*0' ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc
ncap2 -A -s 'sand_03=temp*0' ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc
ncap2 -A -s 'sand_04=temp*0' ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc
ncap2 -A -s 'sand_05=temp*0' ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc
ncap2 -A -s 'sand_06=temp*0' ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc
ncap2 -A -s 'sand_07=temp*0' ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc
ncap2 -A -s 'sand_08=temp*0' ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc

# Modify attributes:
ncatted -O -a field,sand_01,m,c,'sand_01, scalar, series' ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc
ncatted -O -a long_name,sand_01,m,c,'time-averaged suspended noncohesive sediment, size class 01' ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc
ncatted -O -a units,sand_01,m,c,'kilogram meter-3' ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc

ncatted -O -a field,sand_02,m,c,'sand_02, scalar, series' ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc
ncatted -O -a long_name,sand_02,m,c,'time-averaged suspended noncohesive sediment, size class 02' ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc
ncatted -O -a units,sand_02,m,c,'kilogram meter-3' ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc

ncatted -O -a field,sand_03,m,c,'sand_03, scalar, series' ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc
ncatted -O -a long_name,sand_03,m,c,'time-averaged suspended noncohesive sediment, size class 03' ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc
ncatted -O -a units,sand_03,m,c,'kilogram meter-3' ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc

ncatted -O -a field,sand_04,m,c,'sand_04, scalar, series' ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc
ncatted -O -a long_name,sand_04,m,c,'time-averaged suspended noncohesive sediment, size class 04' ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc
ncatted -O -a units,sand_04,m,c,'kilogram meter-3' ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc

ncatted -O -a field,sand_05,m,c,'sand_05, scalar, series' ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc
ncatted -O -a long_name,sand_05,m,c,'time-averaged suspended noncohesive sediment, size class 05' ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc
ncatted -O -a units,sand_05,m,c,'kilogram meter-3' ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc

ncatted -O -a field,sand_06,m,c,'sand_06, scalar, series' ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc
ncatted -O -a long_name,sand_06,m,c,'time-averaged suspended noncohesive sediment, size class 06' ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc
ncatted -O -a units,sand_06,m,c,'kilogram meter-3' ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc

ncatted -O -a field,sand_07,m,c,'sand_07, scalar, series' ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc
ncatted -O -a long_name,sand_07,m,c,'time-averaged suspended noncohesive sediment, size class 07' ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc
ncatted -O -a units,sand_07,m,c,'kilogram meter-3' ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc

ncatted -O -a field,sand_08,m,c,'sand_08, scalar, series' ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc
ncatted -O -a long_name,sand_08,m,c,'time-averaged suspended noncohesive sediment, size class 08' ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc
ncatted -O -a units,sand_08,m,c,'kilogram meter-3' ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc

# change sand vars name (sandfrac_XX):
#ncrename -O -v sand_01,sandfrac_01 ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc
#ncrename -O -v sand_02,sandfrac_02 ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc
#ncrename -O -v sand_03,sandfrac_03 ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc
#ncrename -O -v sand_04,sandfrac_04 ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc
#ncrename -O -v sand_05,sandfrac_05 ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc
#ncrename -O -v sand_06,sandfrac_06 ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc
#ncrename -O -v sand_07,sandfrac_07 ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc
#ncrename -O -v sand_08,sandfrac_08 ocean_rst_mm10_sedkg.nc ocean_rst_mm10_sedkg.nc

