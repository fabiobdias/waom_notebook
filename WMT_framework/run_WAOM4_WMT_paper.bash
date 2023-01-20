#!/bin/bash
#PBS -P v45
#PBS -q hugemem
#PBS -l ncpus=4
#PBS -l mem=1000GB
#PBS -l walltime=10:00:00
#PBS -l software=python
#PBS -l wd
#PBS -l storage=gdata/hh5+scratch/gh9
###PBS -l storage=scratch/gh9

 
# Run Python applications
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_WMT_vint_Hsbl_decomposition_shelf+iceshelf_annual.py > $PBS_JOBID.log

#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_WMT_vint_Hsbl_decomposition_shelf_annual.py > $PBS_JOBID.log

#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_WMT_vint_Hsbl_decomposition_iceshelf_annual.py > $PBS_JOBID.log

# sensitive expts:
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_notides_WMT_vint_Hsbl_decomposition_shelf+iceshelf_annual.py > $PBS_JOBID.log 
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_notides_WMT_vint_Hsbl_decomposition_shelf_annual.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_notides_WMT_vint_Hsbl_decomposition_iceshelf_annual.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_notides_WMT_vint_Hsbl_decomposition_iceshelf_annual_subregions.py > $PBS_JOBID.log

#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_10km-bathy_WMT_vint_Hsbl_decomposition_shelf+iceshelf_annual.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_10km-bathy_WMT_vint_Hsbl_decomposition_shelf_annual.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_10km-bathy_WMT_vint_Hsbl_decomposition_iceshelf_annual.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_10km-bathy_WMT_vint_Hsbl_decomposition_iceshelf_annual_subregions.py > $PBS_JOBID.log

# sensitive expts full-depth integ.:
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_notides_WMT_vint_Full_decomposition_shelf+iceshelf_annual.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_notides_WMT_vint_Full_decomposition_shelf_annual.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_notides_WMT_vint_Full_decomposition_iceshelf_annual.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_notides_WMT_vint_Full_decomposition_iceshelf_annual_subregions.py > $PBS_JOBID.log

#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_10km-bathy_WMT_vint_Full_decomposition_shelf+iceshelf_annual.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_10km-bathy_WMT_vint_Full_decomposition_shelf_annual.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_10km-bathy_WMT_vint_Full_decomposition_iceshelf_annual.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_10km-bathy_WMT_vint_Full_decomposition_iceshelf_annual_subregions.py > $PBS_JOBID.log

#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_WMT_vint_Full_decomposition_shelf+iceshelf_annual.py > $PBS_JOBID.log
python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_WMT_vint_Full_decomposition_shelf_annual.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_WMT_vint_Full_decomposition_iceshelf_annual.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_WMT_vint_Full_decomposition_iceshelf_annual_subregions.py > $PBS_JOBID.log
