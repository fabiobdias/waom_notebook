#!/bin/bash
#PBS -P e14
#PBS -q express
#PBS -l ncpus=8
#PBS -l mem=128GB
#PBS -l walltime=15:00:00
#PBS -l software=python
#PBS -l wd
#PBS -l storage=gdata/hh5

module use /g/data/hh5/public/modules
module load conda/analysis3-unstable

# Run Python applications
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/Maps_validations/WAOM10_shflim_S_ORAS5em_0.25Q_EN4_validation+sea-ice_maps.py > $PBS_JOBID.log

#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/Maps_validations/WAOM10_shflim_S_0.25Q_EN4_validation+sea-ice_maps.py > $PBS_JOBID.log

#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/Maps_validations/WAOM10extend_shflim_S_0.25Q_EN4_validation+sea-ice_maps.py > $PBS_JOBID.log

#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/Maps_validations/WAOM10extend_shflim_S_0.25Q_ORAS5em_EN4_validation+sea-ice_maps.py > $PBS_JOBID.log

#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/Maps_validations/WAOM10shrink_shflim_S_0.25Q_ORAS5em_EN4_validation+sea-ice_maps.py > $PBS_JOBID.log

# 4x4-NOTIDE vs Susheel
python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/Maps_validations/Meltrate_comparison_Susheel_WAOM4-NOTIDE.py > $PBS_JOBID.log

# 4x4-NOTIDE bottom TS
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/Maps_validations/Bottom_TS_comparison_Schmidtko_WAOM4-NOTIDE.py > $PBS_JOBID.log
