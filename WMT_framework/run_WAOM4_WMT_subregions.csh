#!/bin/bash
#SBATCH --job-name=postproc
#SBATCH --account=Project_2000339
#SBATCH --time=05:00:00
#SBATCH --partition=hugemem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1024G

# mld integ.:
#srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_WMT_vint_Hsbl_decomposition_shelf_annual_subregions.py
#srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_WMT_vint_Hsbl_decomposition_iceshelf_annual_subregions.py
#srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_WMT_vint_Hsbl_decomposition_shelf+iceshelf_annual_subregions.py

# full-depth integ.:
srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_WMT_vint_Full_decomposition_shelf_annual_subregions.py
#srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_WMT_vint_Full_decomposition_iceshelf_annual_subregions.py
#srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_WMT_vint_Full_decomposition_shelf+iceshelf_annual_subregions.py
