#!/bin/bash
#PBS -P v45
#PBS -q normal 
#PBS -l ncpus=1
#PBS -l mem=190GB
#PBS -l walltime=02:00:00
#PBS -l software=python
#PBS -l wd
#PBS -l storage=gdata/hh5+scratch/gh9
###PBS -l storage=scratch/gh9

 
# Run Python applications
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_WMTmaps_Hsbl_decomposition_shelf_RFIS.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_notides_WMTmaps_Hsbl_decomposition_shelf_RFIS.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_10km-bathy_WMTmaps_Hsbl_decomposition_shelf_RFIS.py > $PBS_JOBID.log

# maps for full-depth integ. and new extended 20yr waom4:
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_WMTmaps_Full_decomposition_shelf_RFIS.py > $PBS_JOBID.log

## save WMT maps for each density class (27.0:27.9):
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_notides_WMTmaps_Full_decomposition_iceshelf_save_maps.py > $PBS_JOBID.log
python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_10km-bathy_WMTmaps_Full_decomposition_iceshelf_save_maps.py > $PBS_JOBID.log
