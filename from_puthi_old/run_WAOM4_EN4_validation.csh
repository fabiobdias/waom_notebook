#!/bin/bash
#SBATCH --job-name=postproc
#SBATCH --account=Project_2000789
#SBATCH --time=5:00:00
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=128G

#module load python-data/3.7.6-1

#srun python /users/boeiradi/COLD_project/postprocessing/WAOM4_EN4_validation+sea-ice_maps.py
srun python /users/boeiradi/COLD_project/postprocessing/WAOM4_shflim_S_EN4_validation+sea-ice_maps.py
