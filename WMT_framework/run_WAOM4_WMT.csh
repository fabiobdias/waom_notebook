#!/bin/bash
#SBATCH --job-name=postproc
#SBATCH --account=Project_2000339
#SBATCH --time=15:00:00
#SBATCH --partition=hugemem ###small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1024G ###256G

#srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_save_MLD_vint_Hsbl_puhti.py
#srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_save_MLD_vint_Hsbl_m.s-1.py
#srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_save_MLD_vint_Hsbl_m.s-1_iteractive_copy.py

# full-depth integrated:
srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_save_Full_vint_m.s-1.py

#srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_WMT_vint_Hsbl_decomposition_shelf+iceshelf_open_annual.py

#srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_WMT_vint_Hsbl_decomposition_iceshelf_annual.py

#srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_WMT_vint_Hsbl_decomposition_shelf_open_annual.py
