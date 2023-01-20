#!/bin/bash
#SBATCH --job-name=postproc
#SBATCH --account=Project_2000339
#SBATCH --time=04:00:00
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=64G

#module load python-data/3.7.6-1
export HDF5_USE_FILE_LOCKING="FALSE"

srun python /users/boeiradi/COLD_project/postprocessing/WAOM4_TS_diagrams_Barrier_5daily.py
