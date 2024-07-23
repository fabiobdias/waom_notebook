#!/bin/bash
#PBS -P e14
#PBS -q normal
#PBS -l ncpus=48
#PBS -l mem=190GB
#PBS -l walltime=05:00:00
#PBS -l software=python
#PBS -l wd
#PBS -l storage=gdata/hh5+scratch/gh9+scratch/gi0
###PBS -l storage=scratch/gh9

module use /g/data/hh5/public/modules
module load conda/analysis3-unstable

# run OHB closure:

# run cross-1500m transport/vars to save on tmp file:
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM10extend_shflim_S_0.25Q_OHB_shelf_budget_closure.py > $PBS_JOBID.log

#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM10_OHB_hadv_regional_using_lon_masks.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM10_OHB_vadv_regional_using_lon_masks.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM10_OHB_hdiff_regional_using_lon_masks.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM10_OHB_vdiff_regional_using_lon_masks.py > $PBS_JOBID.log

#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM10_OHB_hadv_regional_using_lon_masks_iceshelf.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM10_OHB_vadv_regional_using_lon_masks_iceshelf.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM10_OHB_hdiff_regional_using_lon_masks_iceshelf.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM10_OHB_vdiff_regional_using_lon_masks_iceshelf.py > $PBS_JOBID.log

# Basal melt (m/yr) averaged over the longitudinal bins:
python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM10_basal_melt_avg_regional.py > $PBS_JOBID.log



