#!/bin/bash
#SBATCH --job-name=postproc
#SBATCH --account=Project_2000789
#SBATCH --time=1:00:00
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16G

#module load python-data/3.7.6-1

srun python /users/boeiradi/COLD_project/postprocessing/WAOM10_EN4_validation+sea-ice_maps.py 

