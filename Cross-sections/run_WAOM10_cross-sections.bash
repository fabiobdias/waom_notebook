#!/bin/bash
#PBS -P e14
###PBS -q hugemem
#PBS -q normal
#PBS -l ncpus=48
###PBS -l mem=1500GB
#PBS -l mem=190GB
#PBS -l walltime=24:00:00
#PBS -l software=python
#PBS -l wd
#PBS -l storage=gdata/hh5+scratch/gh9+scratch/gi0

module use /g/data/hh5/public/modules
#module load conda/analysis3-unstable
module load conda/analysis3-24.01

# run OHB closure:

# run cross-1500m transport/vars to save on tmp file:
python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/Cross-sections/Sections_TS_circumpolar_regional_Wind_ModDrag_expts.py > $PBS_JOBID.log
