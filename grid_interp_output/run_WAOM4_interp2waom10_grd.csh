#!/bin/bash
#SBATCH --job-name=postproc
#SBATCH --account=Project_2000339
#SBATCH --time=15:00:00
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=256G

#module load python-data/3.7.6-1

#srun python /users/boeiradi/COLD_project/postprocessing/WAOM4_interp2waom10_grd.py

#srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/grid_interp_output/WAOM4_shflim_S_ORAS5em_0.25Q_interp2waom10_grd.py

#srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/grid_interp_output/WAOM4extend_shflim_S_0.25Q_interp2waom10_grd.py

srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/grid_interp_output/WAOM2extend_shflim_S_0.25Q_interp2waom10_grd.py
