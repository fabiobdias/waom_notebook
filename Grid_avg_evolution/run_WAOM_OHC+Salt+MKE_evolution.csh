#!/bin/bash
#SBATCH --job-name=postproc
#SBATCH --account=Project_2000789
#SBATCH --time=5:00:00
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=256G

#module load python-data/3.7.6-1

#srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/Grid_avg_evolution/WAOM10x4km_shflim_S_ORAS5emXECCO2_0.25Q_OHC+Salt+MKE_evolution.py
#srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/Grid_avg_evolution/WAOM10x4extend_shflim_S_0.25Q_OHC+Salt+MKE_evolution.py
#srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/Grid_avg_evolution/WAOM10x4x2extend_shflim_S_0.25Q_OHC+Salt+MKE_evolution.py

#srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/Grid_avg_evolution/WAOM10x4x2extend_shflim_S_0.25Q_OHC+OSC_evolution_shelf.py
#srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/Grid_avg_evolution/WAOM10x4x2extend_shflim_S_0.25Q_OHC+OSC+MKE_evolution_shelf.py

# OSC evolution with subregions
#srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/Grid_avg_evolution/WAOM4extend_shflim_S_0.25Q_SSS+OSC_evolution_shelf.py
srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/Grid_avg_evolution/WAOM4x2extend_shflim_S_0.25Q_SSS+OSC_evolution_shelf.py
