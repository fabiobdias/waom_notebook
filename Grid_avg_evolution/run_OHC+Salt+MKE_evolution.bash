#!/bin/bash
#PBS -P v45
#PBS -q express
#PBS -l ncpus=1
#PBS -l mem=64GB
#PBS -l walltime=1:00:00
#PBS -l software=python
#PBS -l wd
#PBS -l storage=gdata/hh5
###PBS -l storage=scratch/gh9

 
# Run Python applications
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/Grid_avg_evolution/WAOM10x4km_shflim_S_ORAS5em_0.25Q_OHC+Salt+MKE_evolution.py > $PBS_JOBID.log

python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/Grid_avg_evolution/WAOM10x4km_shflim_S_ORAS5emXECCO2_0.25Q_OHC+Salt+MKE_evolution.py > $PBS_JOBID.log


