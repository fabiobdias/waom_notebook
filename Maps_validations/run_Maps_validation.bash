#!/bin/bash
#PBS -P v45
#PBS -q normal 
#PBS -l ncpus=1
#PBS -l mem=16GB
#PBS -l walltime=1:00:00
#PBS -l software=python
#PBS -l wd
#PBS -l storage=gdata/hh5

 
# Run Python applications
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/Maps_validations/WAOM10_shflim_S_ORAS5em_0.25Q_EN4_validation+sea-ice_maps.py > $PBS_JOBID.log

python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/Maps_validations/WAOM10_shflim_S_0.25Q_EN4_validation+sea-ice_maps.py > $PBS_JOBID.log
