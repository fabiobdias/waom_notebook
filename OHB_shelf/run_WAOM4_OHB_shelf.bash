#!/bin/bash
#PBS -P gv90
#PBS -q hugemem
#PBS -l ncpus=8
#PBS -l mem=1500GB
#PBS -l walltime=10:00:00
#PBS -l software=python
#PBS -l wd
#PBS -l storage=gdata/hh5+scratch/gh9+scratch/gi0
###PBS -l storage=scratch/gh9

module use /g/data/hh5/public/modules
module load conda/analysis3-unstable

# Run Python applications
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/WAOM4extend_shflim_S_0.25Q_Cross-1500m-Transport.py > $PBS_JOBID.log
python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/WAOM4extend_shflim_S_0.25Q_notides_Cross-1500m-Transport.py > $PBS_JOBID.log
