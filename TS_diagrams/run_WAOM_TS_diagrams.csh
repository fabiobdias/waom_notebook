#!/bin/bash
#SBATCH --job-name=postproc
#SBATCH --account=Project_2000789
#SBATCH --time=5:00:00
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=128G

#module load python-data/3.7.6-1


#srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/TS_diagrams/WAOM4_shflim_S_ORAS5em_0.25Q_TS_diagrams_DSWregions_annual_yr15.py
#srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/TS_diagrams/WAOM4_gridded_shflim_S_ORAS5em_0.25Q_TS_diagrams_DSWregions_annual_yr15.py

## using waom4 native grid:
#srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/TS_diagrams/WAOM4_shflim_S_ORAS5em_0.25Q_TS_diagrams_save_sections_DSWregions_5daily.py
#srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/TS_diagrams/WAOM4_shflim_S_ORAS5em_0.25Q_TS_diagrams_plot_sections_DSWregions_5daily.py

## using interpolated fields to waom10km: grid_interp_output/run_WAOM4_interp2waom10_grd.csh
#srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/TS_diagrams/WAOM4_shflim_S_ORAS5em_0.25Q_TS_diagrams_DSWregions_from10km_grd_annual.py

# waom4extend:

#srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/TS_diagrams/WAOM4extend_shflim_0.25Q_TS_diagrams_save_sections_DSWregions_5daily.py
#srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/TS_diagrams/WAOM4extend_shflim_0.25Q_TS_diagrams_plot_sections_DSWregions_5daily.py

#srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/TS_diagrams/WAOM2extend_shflim_0.25Q_TS_diagrams_save_sections_DSWregions_5daily.py
#srun python -s /users/boeiradi/COLD_project/postprocessing/waom_notebook/TS_diagrams/WAOM2extend_shflim_0.25Q_TS_diagrams_plot_sections_DSWregions_5daily.py
