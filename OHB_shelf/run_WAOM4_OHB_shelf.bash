#!/bin/bash
#PBS -P v45
#PBS -l ncpus=48
#PBS -q megamem
#PBS -l mem=2990GB
###PBS -q hugemem
###PBS -l mem=1460GB
###PBS -q normal
###PBS -l mem=190GB
#PBS -l walltime=48:00:00
#PBS -l software=python
#PBS -l wd
#PBS -l storage=gdata/hh5+scratch/gh9+scratch/gi0
###PBS -l storage=scratch/gh9

module use /g/data/hh5/public/modules
module load conda/analysis3-unstable

# run OHB closure:
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/WAOM4extend_shflim_S_0.25Q_notides_OHB_shelf_budget_closure.py > $PBS_JOBID.log

# run cross-1500m transport:
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/WAOM4extend_shflim_S_0.25Q_Cross-1500m-Transport.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/WAOM4extend_shflim_S_0.25Q_Cross-1500m_temp_daily.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/WAOM4extend_shflim_S_0.25Q_Cross-1500m_salt_daily.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/WAOM4extend_shflim_S_0.25Q_Cross-1500m_z_rho_daily.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/WAOM4extend_shflim_S_0.25Q_Cross-1500m_Tf_heat_transp_daily.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/WAOM4extend_shflim_S_0.25Q_Cross-1500m_heat_transp_daily.py > $PBS_JOBID.log
python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/WAOM4extend_shflim_S_0.25Q_Cross-1500m_vol_transp_daily.py > $PBS_JOBID.log

# run cross-CalvingFront transport:
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/WAOM4extend_shflim_S_0.25Q_Cross-CalvingFront-Transport.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/WAOM4extend_shflim_S_0.25Q_notides_Cross-CalvingFront-Transport.py > $PBS_JOBID.log
