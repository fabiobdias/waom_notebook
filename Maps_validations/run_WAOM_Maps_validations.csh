#!/bin/bash
#SBATCH --job-name=postproc
#SBATCH --account=Project_2000339
#SBATCH --time=15:00:00
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=256G

#module load python-data/3.7.6-1

#srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/Maps_validations/WAOM4_shflim_S_ORAS5em_0.25Q_EN4_validation+sea-ice_maps.py

#srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/Maps_validations/WAOM4extend_shflim_S_0.25Q_EN4_validation+sea-ice_maps.py

#srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/Maps_validations/WAOM2extend_shflim_S_0.25Q_EN4_validation+sea-ice_maps.py

# for paper, TS bottom and MKE/EKE figures:
srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/Maps_validations/WAOM10x4x2extend_shflim_S_0.25Q_MKExEKE_bottomTS.py
# only MKE/EKE using bottom velocities:
#srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/Maps_validations/WAOM10x4x2extend_shflim_S_0.25Q_MKExEKE_bottom.py
# only MKE/EKE using surface velocities:
#srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/Maps_validations/WAOM10x4x2extend_shflim_S_0.25Q_MKExEKE_surface.py
#srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/Maps_validations/WAOM10x4x2extend_shflim_S_0.25Q_TSVel_botXsfc.py
# test WAOM2 using yr4:
#srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/Maps_validations/WAOM10x4x2extend_shflim_S_0.25Q_MKExEKE_bottomTS_testWAOM2diffYear.py


# monthly maps:
#srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/Maps_validations/WAOM10x4x2extend_shflim_S_0.25Q_MKExEKE_bottomTS_monthly.py
#srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/Maps_validations/WAOM10x4x2extend_shflim_S_0.25Q_TSVel_botXsfc_monthly.py

# Susheel melt rate X WAOM
#srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/Maps_validations/Meltrate_comparison_Susheel.py

