#!/bin/bash
#SBATCH --job-name=postproc
#SBATCH --account=Project_2000789
#SBATCH --time=15:00:00
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=356G

#module load python-data/3.7.6-1

srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/Grid_avg_evolution/WAOM10x4x2extend_shflim_S_ORAS5emXECCO2_0.25Q_mass_loss_spinup.py 

