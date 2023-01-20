#!/bin/bash
#PBS -P v45
#PBS -q express 
#PBS -l ncpus=1
#PBS -l mem=16GB
#PBS -l walltime=1:00:00
#PBS -l software=python
#PBS -l wd
#PBS -l storage=gdata/hh5

 
# Run Python applications

# ORAS5em
python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/TS_diagrams/WAOM10_shflim_S_ORAS5em_0.25Q_TS_diagrams_DSWregions_annual_yr10.py > $PBS_JOBID.log
python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/TS_diagrams/WAOM10_shflim_S_ORAS5em_0.25Q_TS_diagrams_DSWregions_annual_yr20.py > $PBS_JOBID.log

# ECCO2
python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/TS_diagrams/WAOM10_shflim_S_0.25Q_TS_diagrams_DSWregions_annual_yr10.py > $PBS_JOBID.log
python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/TS_diagrams/WAOM10_shflim_S_0.25Q_TS_diagrams_DSWregions_annual_yr20.py > $PBS_JOBID.log

