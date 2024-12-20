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

# run cross-1000m transport/vars to save on tmp file:
# daily outputs:
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1000m/WAOM10_Cross-1000m_temp_daily.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1000m/WAOM10_Cross-1000m_salt_daily.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1000m/WAOM10_Cross-1000m_z_rho_daily.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1000m/WAOM10_Cross-1000m_vol_transp_daily.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1000m/WAOM10_Cross-1000m_heat_transp_daily.py > $PBS_JOBID.log
python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1000m/WAOM10_Cross-1000m_Tf_heat_transp_daily.py > $PBS_JOBID.log



