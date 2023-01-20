#!/bin/bash
#PBS -P v45
#PBS -q hugemem
#PBS -l ncpus=8
#PBS -l mem=1500GB
#PBS -l walltime=10:00:00
#PBS -l software=python
#PBS -l wd
#PBS -l storage=gdata/hh5+scratch/gh9
###PBS -l storage=scratch/gh9

 
# Run Python applications
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM2extend_shflim_S_0.25Q_save_MLD_vint_Hsbl_gadi.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM2extend_shflim_S_0.25Q_save_MLD_vint_Hsbl_m.s-1.py > $PBS_JOBID.log

# full-depth integrated:
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM2extend_shflim_S_0.25Q_save_Full_vint_m.s-1.py > $PBS_JOBID.log
# attempt w/o the bottom layer 
python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM2extend_shflim_S_0.25Q_save_SemiFull_vint_m.s-1.py > $PBS_JOBID.log

# WMT rates per density, shelf+iceshelf:
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM2extend_shflim_S_0.25Q_WMT_vint_Hsbl_decomposition_shelf+iceshelf_annual.py > $PBS_JOBID.log

#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM2extend_shflim_S_0.25Q_WMT_vint_Hsbl_decomposition_fulldomain_annual.py > $PBS_JOBID.log

#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM2extend_shflim_S_0.25Q_WMT_vint_Hsbl_decomposition_noice_annual.py > $PBS_JOBID.log
