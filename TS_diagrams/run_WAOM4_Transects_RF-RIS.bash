#!/bin/bash
#PBS -P gi0
#PBS -q normal  
#PBS -l ncpus=1
#PBS -l mem=190GB
#PBS -l walltime=3:00:00
#PBS -l software=python
#PBS -l wd
#PBS -l storage=gdata/hh5+scratch/gh9

 
# Run Python applications

python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/TS_diagrams/WAOM4extend_notides_Transects_RF-RIS_gridline.py > $PBS_JOBID.log
