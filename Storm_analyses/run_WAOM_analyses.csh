#!/bin/bash
#SBATCH --job-name=postproc
#SBATCH --account=Project_2000789
#SBATCH --time=5:00:00
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=128G

#module load python-data/3.7.6-1

# - Surface fluxes: using jupyter-notebook
# WAOM4extend_yr10_StormEvents_Sfc_fluxes.ipynb
# WAOM10extend_yr21_StormEvents_Sfc_fluxes.ipynb

# Cross-sections:
#srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/Storm_analyses/WAOM10extend_Sections_TSDens_EastAntarctica.py
#srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/Storm_analyses/WAOM4extend_Sections_TSDens_EastAntarctica.py
#srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/Storm_analyses/WAOM2extend_Sections_TSDens_EastAntarctica.py
srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/Storm_analyses/WAOM10extend_Sections_TSDens_Barrier.py
srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/Storm_analyses/WAOM4extend_Sections_TSDens_Barrier.py
srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/Storm_analyses/WAOM2extend_Sections_TSDens_Barrier.py

# - TS-diagrams
# waom10extend
#srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/Storm_analyses/WAOM10extend_shflim_0.25Q_TS_diagrams_save_sections_EastAntarctica_5daily.py
#srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/Storm_analyses/WAOM10extend_shflim_0.25Q_TS_diagrams_plot_sections_EastAntarctica_5daily.py

# waom4extend:
#srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/Storm_analyses/WAOM4extend_shflim_0.25Q_TS_diagrams_save_sections_EastAntarctica_5daily.py
#srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/Storm_analyses/WAOM4extend_shflim_0.25Q_TS_diagrams_plot_sections_EastAntarctica_5daily.py

# waom2extend:
#srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/Storm_analyses/WAOM2extend_shflim_0.25Q_TS_diagrams_save_sections_EastAntarctica_5daily.py
#srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/Storm_analyses/WAOM2extend_shflim_0.25Q_TS_diagrams_plot_sections_EastAntarctica_5daily.py
