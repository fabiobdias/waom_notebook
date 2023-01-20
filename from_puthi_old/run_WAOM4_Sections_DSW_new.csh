#!/bin/bash
#SBATCH --job-name=postproc
#SBATCH --account=Project_2000339
#SBATCH --time=04:00:00
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=64G

#module load python-data/3.7.6-1

#srun python /users/boeiradi/COLD_project/postprocessing/WAOM4_Sections_TSDens_annual_new.py
#srun python /users/boeiradi/COLD_project/postprocessing/WAOM4_Sections_TSDens_EastAntarctica.py

srun python /users/boeiradi/COLD_project/postprocessing/WAOM4_shflim_S_Sections_TSDens_annual_new.py
