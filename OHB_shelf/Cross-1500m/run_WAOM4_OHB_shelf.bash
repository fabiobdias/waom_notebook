#!/bin/bash
#PBS -P e14
#PBS -q normal
#PBS -l ncpus=48
#PBS -l mem=190GB
#PBS -l walltime=48:00:00
#PBS -l software=python
#PBS -l wd
#PBS -l storage=gdata/hh5+scratch/gh9+scratch/gi0
###PBS -l storage=scratch/gh9

module use /g/data/hh5/public/modules
module load conda/analysis3-unstable

## WAOM4
# run OHB closure:

# run cross-1500m transport/vars to save on tmp file:
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM4_Cross-1500m_temp_daily.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM4_Cross-1500m_salt_daily.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM4_Cross-1500m_z_rho_daily.py > $PBS_JOBID.log

#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM4_Cross-1500m_vol_transp_daily.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM4_Cross-1500m_heat_transp_daily.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM4_Cross-1500m_Tf_heat_transp_daily.py > $PBS_JOBID.log


## WAOM4_NOTIDE
# run cross-1500m transport/vars to save on tmp file:
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM4_notides_Cross-1500m_temp_daily.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM4_notides_Cross-1500m_salt_daily.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM4_notides_Cross-1500m_z_rho_daily.py > $PBS_JOBID.log

#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM4_notides_Cross-1500m_vol_transp_daily.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM4_notides_Cross-1500m_heat_transp_daily.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM4_notides_Cross-1500m_Tf_heat_transp_daily.py > $PBS_JOBID.log

# V2 with pre-selected contour:
## WAOM4
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM4_Cross-1500m_temp_daily_v2.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM4_Cross-1500m_salt_daily_v2.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM4_Cross-1500m_z_rho_daily_v2.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM4_Cross-1500m_vol_transp_daily_v2.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM4_Cross-1500m_Coordinates_daily_v2.py > $PBS_JOBID.log

# v3
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM4_Cross-1500m_temp_daily_v3.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM4_Cross-1500m_salt_daily_v3.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM4_Cross-1500m_z_rho_daily_v3.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM4_Cross-1500m_vol_transp_daily_v3.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM4_Cross-1500m_Coordinates_daily_v3.py > $PBS_JOBID.log
python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM4_Cross-1500m_heat_transp_daily_v3.py > $PBS_JOBID.log

## WAOM4_NOTIDE
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM4_notides_Cross-1500m_temp_daily_v2.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM4_notides_Cross-1500m_salt_daily_v2.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM4_notides_Cross-1500m_z_rho_daily_v2.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM4_notides_Cross-1500m_vol_transp_daily_v2.py > $PBS_JOBID.log

# v3
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM4_notides_Cross-1500m_temp_daily_v3.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM4_notides_Cross-1500m_salt_daily_v3.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM4_notides_Cross-1500m_z_rho_daily_v3.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM4_notides_Cross-1500m_vol_transp_daily_v3.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM4_notides_Cross-1500m_heat_transp_daily_v3.py > $PBS_JOBID.log
