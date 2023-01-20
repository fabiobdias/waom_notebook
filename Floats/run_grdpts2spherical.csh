#!/bin/bash
#SBATCH --job-name=postproc
#SBATCH --account=Project_2000339
#SBATCH --time=15:00:00
#SBATCH --partition=hugemem ###small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1024G ###256G

srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/Floats/grid_points_to_spherical.py
