#!/bin/bash
#PBS -P e14
#PBS -q normal
#PBS -l ncpus=8
#PBS -l mem=128GB
#PBS -l walltime=03:00:00
#PBS -l software=python
#PBS -l wd
#PBS -l storage=gdata/hh5+scratch/gh9
###PBS -l storage=scratch/gh9

 
# Run Python applications
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/Grid_avg_evolution/WAOM10x4km_shflim_S_ORAS5em_0.25Q_OHC+Salt+MKE_evolution.py > $PBS_JOBID.log

#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/Grid_avg_evolution/WAOM10x4km_shflim_S_ORAS5emXECCO2_0.25Q_OHC+Salt+MKE_evolution.py > $PBS_JOBID.log

#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/Grid_avg_evolution/WAOM10extend_shflim_S_0.25Q_OHC+Salt+MKE_evolution.py > $PBS_JOBID.log

#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/Grid_avg_evolution/WAOM10extend_x4km_shflim_S_ORAS5emXECCO2_0.25Q_OHC+Salt+MKE_evolution.py > $PBS_JOBID.log

#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/Grid_avg_evolution/WAOM10extendXshrink_shflim_S_ORAS5emXECCO2_0.25Q_OHC+Salt+MKE_evolution.py > $PBS_JOBID.log

# WAOM4 x WAOM4-NOTIDE
python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/Grid_avg_evolution/WAOM4extend_tidesXnotides_OHC+OSC+MKE_evolution_shelf.py > $PBS_JOBID.log
