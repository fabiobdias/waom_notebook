#!/bin/bash
#SBATCH --job-name=postproc
#SBATCH --account=Project_2000339
#SBATCH --time=5:00:00
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=128G

#module load python-data/3.7.6-1

PYTHONNOUSERSITE=true python
module load nco

#srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/TS_diagrams/WAOM10extend_shflim_S_0.25Q_TS_diagrams_save_sections_RF-RIS_5daily.py
#srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/TS_diagrams/WAOM10extend_shflim_S_0.25Q_TS_diagrams_save_sections_RF-RIS_5daily_iceshelf.py
# extended version to match with WMT TS-space analysis:
srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/TS_diagrams/WAOM10extend_shflim_S_0.25Q_TS_diagrams_save_sections_ext_RF-RIS_5daily_iceshelf.py

#cd ~/COLD_project/postprocessing/ncdf_tmp/

#./concat_WAOM10extend_shflim_S_sec_files.bash

#srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/TS_diagrams/WAOM10extend_shflim_S_0.25Q_TS_diagrams_plot_sections_RF-RIS_5daily.py
