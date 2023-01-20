#!/bin/bash
#PBS -P gi0
#PBS -q hugemem
#PBS -l ncpus=8
#PBS -l mem=1000GB 
#PBS -l walltime=10:00:00
#PBS -l software=python
#PBS -l wd
#PBS -l storage=gdata/hh5+scratch/gh9
###PBS -l storage=scratch/gh9

 
# Run Python applications
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM2extend_shflim_S_0.25Q_WMTmaps_Hsbl_decomposition_shelf.py > $PBS_JOBID.log

#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM2extend_shflim_S_0.25Q_WMTmaps_Hsbl_decomposition_iceshelf.py > $PBS_JOBID.log

# maps for RFIS and RIS
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM2extend_shflim_S_0.25Q_WMTmaps_Hsbl_decomposition_shelf_RFIS.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM2extend_shflim_S_0.25Q_WMTmaps_Hsbl_decomposition_shelf_RIS.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM2extend_shflim_S_0.25Q_WMTmaps_Hsbl_decomposition_shelf_WAnt.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM2extend_shflim_S_0.25Q_WMTmaps_Hsbl_decomposition_shelf_EAnt.py > $PBS_JOBID.log

#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM2extend_shflim_S_0.25Q_WMTmaps_Hsbl_decomposition_iceshelf_RFIS.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM2extend_shflim_S_0.25Q_WMTmaps_Hsbl_decomposition_iceshelf_RIS.py > $PBS_JOBID.log

# movie with monthly WMT
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM2extend_shflim_S_0.25Q_WMTmaps_Hsbl_decomposition_shelf_RFIS_monthly_movie.py

# maps for RFIS and RIS: full-dept integ.
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM2extend_shflim_S_0.25Q_WMTmaps_Full_decomposition_shelf_RFIS.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM2extend_shflim_S_0.25Q_WMTmaps_Full_decomposition_iceshelf_RFIS.py > $PBS_JOBID.log

## saving maps for each density class (27.0:27.9):
python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM2extend_shflim_S_0.25Q_WMTmaps_Full_decomposition_shelf_save_maps.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM2extend_shflim_S_0.25Q_WMTmaps_Full_decomposition_iceshelf_save_maps.py > $PBS_JOBID.log
