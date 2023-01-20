#!/bin/bash
#PBS -P v45
#PBS -q hugemem 
#PBS -l ncpus=10
#PBS -l mem=1000GB
#PBS -l walltime=10:00:00
#PBS -l software=python
#PBS -l wd
#PBS -l storage=gdata/hh5+scratch/gh9
###PBS -l storage=scratch/gh9

 
# Run Python applications
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_notides_save_MLD_vint_Hsbl_puhti.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_10km-bathy_save_MLD_vint_Hsbl_puhti.py > $PBS_JOBID.log

# for full-depth integ.:
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_notides_save_Full_vint_m.s-1.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_10km-bathy_save_Full_vint_m.s-1.py > $PBS_JOBID.log

# extended WAOM4:
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_yr18_save_Full_vint_m.s-1.py > $PBS_JOBID.log
python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_yr20_save_Full_vint_m.s-1.py > $PBS_JOBID.log

~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
"run_WAOM4_WMT.bash" 23L, 1126C                                                                                                                                                                                                                               1,1           All
