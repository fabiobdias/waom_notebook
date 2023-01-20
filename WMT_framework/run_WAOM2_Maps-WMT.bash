#!/bin/bash
#PBS -P gh9
#PBS -q hugemem
#PBS -l ncpus=8
#PBS -l mem=1000GB
#PBS -l walltime=10:00:00
#PBS -l software=python
#PBS -l wd
#PBS -l storage=gdata/hh5+scratch/gh9
###PBS -l storage=scratch/gh9

 
# Run Python applications
python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM2extend_shflim_S_0.25Q_WMT-maps_vint_Hsbl_decomposition_shelf+iceshelf_open_annual.py > $PBS_JOBID.log

