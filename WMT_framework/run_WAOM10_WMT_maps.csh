#!/bin/bash
#SBATCH --job-name=postproc
#SBATCH --account=Project_2000339
#SBATCH --time=05:00:00
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=256G

#srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/WMT_framework/WAOM10extend_shflim_S_0.25Q_WMTmaps_Hsbl_decomposition_shelf.py

#srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/WMT_framework/WAOM10extend_shflim_S_0.25Q_WMTmaps_Hsbl_decomposition_iceshelf.py

# maps for RFIS, RIS, WAnt, EAnt: MLD integ.
#srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/WMT_framework/WAOM10extend_shflim_S_0.25Q_WMTmaps_Hsbl_decomposition_shelf_RFIS.py
#srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/WMT_framework/WAOM10extend_shflim_S_0.25Q_WMTmaps_Hsbl_decomposition_shelf_RIS.py
#srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/WMT_framework/WAOM10extend_shflim_S_0.25Q_WMTmaps_Hsbl_decomposition_shelf_WAnt.py
#srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/WMT_framework/WAOM10extend_shflim_S_0.25Q_WMTmaps_Hsbl_decomposition_shelf_EAnt.py

#srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/WMT_framework/WAOM10extend_shflim_S_0.25Q_WMTmaps_Hsbl_decomposition_iceshelf_RFIS.py
#srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/WMT_framework/WAOM10extend_shflim_S_0.25Q_WMTmaps_Hsbl_decomposition_iceshelf_RIS.py

# maps for RFIS, RIS, WAnt, EAnt: full-depth integ.
srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/WMT_framework/WAOM10extend_shflim_S_0.25Q_WMTmaps_Full_decomposition_shelf_RFIS.py
# #srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/WMT_framework/WAOM10extend_shflim_S_0.25Q_WMTmaps_Full_decomposition_shelf_RIS.py
# #srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/WMT_framework/WAOM10extend_shflim_S_0.25Q_WMTmaps_Full_decomposition_shelf_WAnt.py
# #srun python /users/boeiradi/COLD_project/postprocessing/waom_notebook/WMT_framework/WAOM10extend_shflim_S_0.25Q_WMTmaps_Full_decomposition_shelf_EAnt.py
#
