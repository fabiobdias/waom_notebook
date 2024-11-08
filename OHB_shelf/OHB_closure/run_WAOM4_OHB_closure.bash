#!/bin/bash
#PBS -P jk72
#PBS -q normal
###hugemem 
###megamem
#PBS -l ncpus=48   
###48
#PBS -l mem=190GB
###1400GB     
###2900GB
#PBS -l walltime=48:00:00 
###48:00:00
#PBS -l software=python
#PBS -l wd
#PBS -l storage=gdata/hh5+scratch/gh9+scratch/gi0
###PBS -l storage=scratch/gh9

module use /g/data/hh5/public/modules
module load conda/analysis3-unstable

# run OHB closure:

# run cross-1500m transport/vars to save on tmp file:
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4extend_shflim_S_0.25Q_OHB_shelf_budget_closure.py > $PBS_JOBID.log

#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_OHB_hadv_regional_using_lon_masks.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_OHB_vadv_regional_using_lon_masks.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_OHB_hdiff_regional_using_lon_masks.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_OHB_vdiff_regional_using_lon_masks.py > $PBS_JOBID.log

#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_OHB_hadv_regional_using_lon_masks_iceshelf.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_OHB_vadv_regional_using_lon_masks_iceshelf.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_OHB_hdiff_regional_using_lon_masks_iceshelf.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_OHB_vdiff_regional_using_lon_masks_iceshelf.py > $PBS_JOBID.log

# 5-daily:
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_OHB_hadv_regional_using_lon_masks_5daily.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_OHB_vadv_regional_using_lon_masks_5daily.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_OHB_hdiff_regional_using_lon_masks_5daily.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_OHB_vdiff_regional_using_lon_masks_5daily.py > $PBS_JOBID.log

#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_OHB_hadv_regional_using_lon_masks_iceshelf_5daily.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_OHB_vadv_regional_using_lon_masks_iceshelf_5daily.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_OHB_hdiff_regional_using_lon_masks_iceshelf_5daily.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_OHB_vdiff_regional_using_lon_masks_iceshelf_5daily.py > $PBS_JOBID.log

# Basal melt (m/yr) averaged over the longitudinal bins:
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_basal_melt_avg_regional.py > $PBS_JOBID.log

# --- notides:
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4extend_shflim_S_0.25Q_notides_OHB_shelf_budget_closure.py > $PBS_JOBID.log

#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_notides_OHB_hadv_regional_using_lon_masks.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_notides_OHB_vadv_regional_using_lon_masks.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_notides_OHB_hdiff_regional_using_lon_masks.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_notides_OHB_vdiff_regional_using_lon_masks.py > $PBS_JOBID.log

#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_notides_OHB_hadv_regional_using_lon_masks_iceshelf.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_notides_OHB_vadv_regional_using_lon_masks_iceshelf.py > $PBS_JOBID.log
# 2 belows are done:
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_notides_OHB_hdiff_regional_using_lon_masks_iceshelf.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_notides_OHB_vdiff_regional_using_lon_masks_iceshelf.py > $PBS_JOBID.log

# Basal melt (m/yr) averaged over the longitudinal bins:
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_notides_basal_melt_avg_regional.py > $PBS_JOBID.log

# 5-daily:
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_notides_OHB_hadv_regional_using_lon_masks_5daily.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_notides_OHB_vadv_regional_using_lon_masks_5daily.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_notides_OHB_hdiff_regional_using_lon_masks_5daily.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_notides_OHB_vdiff_regional_using_lon_masks_5daily.py > $PBS_JOBID.log

#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_notides_OHB_hadv_regional_using_lon_masks_iceshelf_5daily.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_notides_OHB_vadv_regional_using_lon_masks_iceshelf_5daily.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_notides_OHB_hdiff_regional_using_lon_masks_iceshelf_5daily.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_notides_OHB_vdiff_regional_using_lon_masks_iceshelf_5daily.py > $PBS_JOBID.log

## -- using LonLat mask with adjusted Antarctic Peninsula

#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_OHB_hadv_regional_using_lonlat_masks.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_OHB_vadv_regional_using_lonlat_masks.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_OHB_hdiff_regional_using_lonlat_masks.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_OHB_vdiff_regional_using_lonlat_masks.py > $PBS_JOBID.log

#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_OHB_hadv_regional_using_lonlat_masks_iceshelf.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_OHB_vadv_regional_using_lonlat_masks_iceshelf.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_OHB_hdiff_regional_using_lonlat_masks_iceshelf.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_OHB_vdiff_regional_using_lonlat_masks_iceshelf.py > $PBS_JOBID.log
#

# no tide ->>> NEED TO RUN THIS (19/8/24)
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_notides_OHB_hadv_regional_using_lonlat_masks.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_notides_OHB_vadv_regional_using_lonlat_masks.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_notides_OHB_hdiff_regional_using_lonlat_masks.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_notides_OHB_vdiff_regional_using_lonlat_masks.py > $PBS_JOBID.log

#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_notides_OHB_hadv_regional_using_lonlat_masks_iceshelf.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_notides_OHB_vadv_regional_using_lonlat_masks_iceshelf.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_notides_OHB_hdiff_regional_using_lonlat_masks_iceshelf.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_notides_OHB_vdiff_regional_using_lonlat_masks_iceshelf.py > $PBS_JOBID.log



## -- using LonLat mask with adjusted Antarctic Peninsula v2 (also Ross fix)
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_OHB_hadv_regional_using_lonlat_masks_v2.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_OHB_vadv_regional_using_lonlat_masks_v2.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_OHB_hdiff_regional_using_lonlat_masks_v2.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_OHB_vdiff_regional_using_lonlat_masks_v2.py > $PBS_JOBID.log

#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_OHB_hadv_regional_using_lonlat_masks_iceshelf_v2.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_OHB_vadv_regional_using_lonlat_masks_iceshelf_v2.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_OHB_hdiff_regional_using_lonlat_masks_iceshelf_v2.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_OHB_vdiff_regional_using_lonlat_masks_iceshelf_v2.py > $PBS_JOBID.log

# Basal melt (m/yr) averaged over the new lon-lat bins:
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_basal_melt_avg_regional_lonlatBins.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_basal_melt_avg_regional_lonlatBins_v2.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_notides_basal_melt_avg_regional_lonlatBins.py > $PBS_JOBID.log

# basal melting divided into deep/shallow melting:
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_basal_melt_Decomp_avg_regional_lonlatBins.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_basal_melt_Decomp_avg_regional_lonlatBins_v2.py > $PBS_JOBID.log
# 400m and 300m split tests:
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_basal_melt_Decomp_avg_regional_lonlatBins_v2_400m.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_basal_melt_Decomp_avg_regional_lonlatBins_v2_300m.py > $PBS_JOBID.log

# Sea ice:
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_SeaIce_avg_regional_lonlatBins.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_SeaIce_avg_regional_lonlatBins_v2.py > $PBS_JOBID.log

# Wind stress:
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_WindStress_avg_regional_lonlatBins.py > $PBS_JOBID.log
#- test with only 30days:
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_WindStress_avg_regional_lonlatBins_tst.py > $PBS_JOBID.log
# doing daily files:
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_WindStress_avg_regional_lonlatBins_eachDay.py > $PBS_JOBID.log

# Interpolate temperature field to 500m:
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_interp_temp-500m.py > $PBS_JOBID.log

# Run Fig 4 (maps + TS-diag)
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_run_maps_Figure4.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_run_maps_Figure4_newTS.py > $PBS_JOBID.log
python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/OHB_closure/WAOM4_run_maps_Figure4_newTS_bottom.py > $PBS_JOBID.log
