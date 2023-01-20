#!/bin/bash
#SBATCH --job-name=postproc
#SBATCH --account=Project_2000339
#SBATCH --time=05:00:00
#SBATCH --partition=hugemem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=512G

#srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_WMTmaps_Hsbl_decomposition_shelf.py

#srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_WMTmaps_Hsbl_decomposition_iceshelf.py

# maps for RFIS and RIS: MLD integ.
#srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_WMTmaps_Hsbl_decomposition_shelf_RFIS.py
#srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_WMTmaps_Hsbl_decomposition_shelf_RIS.py
#srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_WMTmaps_Hsbl_decomposition_shelf_WAnt.py
#srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_WMTmaps_Hsbl_decomposition_shelf_EAnt.py

#srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_WMTmaps_Hsbl_decomposition_iceshelf_RFIS.py
#srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_WMTmaps_Hsbl_decomposition_iceshelf_RIS.py

# movie with monthly WMT:
#srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_WMTmaps_Hsbl_decomposition_shelf_RFIS_monthly_movie.py

# maps for RFIS and RIS: full-depth integ.
#srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_WMTmaps_Full_decomposition_shelf_RFIS.py
#srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_WMTmaps_Full_decomposition_iceshelf_RFIS.py

## --- saving maps for each density class (27.0:27.9)
#srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_WMTmaps_Full_decomposition_iceshelf_save_maps.py
srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/WMT_framework/WAOM4extend_shflim_S_0.25Q_WMTmaps_Full_decomposition_shelf_save_maps.py
