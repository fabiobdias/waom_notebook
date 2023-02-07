#!/bin/bash
#PBS -q normal 
#PBS -l ncpus=48
#PBS -l mem=128GB
#PBS -l jobfs=400GB
#PBS -l walltime=10:00:00
#PBS -l software=python
#PBS -l wd
 
# Load modules.
module unload intel-fc intel-cc
module load python3/3.8.5
module load hdf5/1.10.5
module load netcdf/4.7.1

export PYTHONPATH=$PYTHONPATH:/g/data3/hh5/tmp/access-om/fbd581/postprocessing/py3libs/lib/python3.8/site-packages/

#Set number of OMP threads
 
#export OMP_NUM_THREADS=$PBS_NCPUS
 
 
# Run Python applications
python3 mPpythonScript.py > $PBS_JOBID.log
