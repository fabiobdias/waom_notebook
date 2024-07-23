#!/bin/bash
#PBS -P e14
#PBS -q hugemem
#PBS -l ncpus=1
#PBS -l mem=1500GB
#PBS -l walltime=5:00:00
#PBS -l software=python
#PBS -l wd
#PBS -l storage=gdata/hh5+scratch/gh9

module use /g/data/hh5/public/modules
#module load conda/analysis3-unstable
module load conda/analysis3-24.01
 
# Run Python applications

# ORAS5em
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/TS_diagrams/WAOM10_shflim_S_ORAS5em_0.25Q_TS_diagrams_DSWregions_annual_yr10.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/TS_diagrams/WAOM10_shflim_S_ORAS5em_0.25Q_TS_diagrams_DSWregions_annual_yr20.py > $PBS_JOBID.log

# ECCO2
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/TS_diagrams/WAOM10_shflim_S_0.25Q_TS_diagrams_DSWregions_annual_yr10.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/TS_diagrams/WAOM10_shflim_S_0.25Q_TS_diagrams_DSWregions_annual_yr20.py > $PBS_JOBID.log

# new waom10extend
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/TS_diagrams/WAOM10extend_shflim_S_0.25Q_TS_diagrams_DSWregions_annual_yr10.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/TS_diagrams/WAOM10extend_shflim_S_0.25Q_TS_diagrams_DSWregions_annual_yr20.py > $PBS_JOBID.log

#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/TS_diagrams/WAOM10extend_shflim_S_0.25Q_ORAS5em_TS_diagrams_DSWregions_annual_yr10.py > $PBS_JOBID.log

#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/TS_diagrams/WAOM10shrink_shflim_S_0.25Q_ORAS5em_TS_diagrams_DSWregions_annual_yr10.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/TS_diagrams/WAOM10shrink_shflim_S_0.25Q_ORAS5em_TS_diagrams_DSWregions_annual_yr20.py > $PBS_JOBID.log

# TS diags for along-contour analyses
python3 /g/data/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/TS_diagrams/TSdiagrams_along1500m_CF.py > $PBS_JOBID.log
