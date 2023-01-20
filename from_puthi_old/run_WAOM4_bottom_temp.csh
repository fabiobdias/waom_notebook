#!/bin/bash
#SBATCH --job-name=postproc
#SBATCH --account=Project_2000339
#SBATCH --time=5:00:00
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=64G

#module load python-data/3.7.6-1

srun python /users/boeiradi/COLD_project/postprocessing/WAOM4_bottom_temp_maps.py
