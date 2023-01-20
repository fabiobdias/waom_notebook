#!/bin/bash
#PBS -P jk72
#PBS -q normal
#PBS -l ncpus=1
#PBS -l mem=16GB
#PBS -l walltime=01:00:00
#PBS -l software=python
#PBS -l wd
#PBS -l storage=gdata/hh5+scratch/gh9
###PBS -l storage=scratch/gh9

 
python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/WMT_framework/WAOM2extend_shflim_S_0.25Q_WMTmaps_TEST.py > $PBS_JOBID.log
