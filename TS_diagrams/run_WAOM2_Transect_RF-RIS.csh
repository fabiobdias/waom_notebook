#!/bin/bash
#SBATCH --job-name=postproc
#SBATCH --account=Project_2000339
#SBATCH --time=3:00:00
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=256G

srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/TS_diagrams/WAOM2extend_Transects_RF-RIS_gridline.py
