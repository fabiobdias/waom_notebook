#!/bin/bash
#PBS -P e14
#PBS -q hugemem
#PBS -l ncpus=48
#PBS -l mem=1500GB
#PBS -l walltime=24:00:00
#PBS -l software=python
#PBS -l wd
#PBS -l storage=gdata/hh5+scratch/gh9+scratch/gi0
###PBS -l storage=scratch/gh9

module use /g/data/hh5/public/modules
#module load conda/analysis3-unstable
module load conda/analysis3-24.01

# run OHB closure:

# run cross-1500m transport/vars to save on tmp file:
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_Cross-1500m_temp_5daily.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_Cross-1500m_salt_5daily.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_Cross-1500m_z_rho_5daily.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_Cross-1500m_vol_transp_5daily.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_Cross-1500m_heat_transp_5daily.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_Cross-1500m_Tf_heat_transp_5daily.py > $PBS_JOBID.log

# daily outputs:
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_Cross-1500m_temp_daily.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_Cross-1500m_salt_daily.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_Cross-1500m_z_rho_daily.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_Cross-1500m_vol_transp_daily.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_Cross-1500m_heat_transp_daily.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_Cross-1500m_Tf_heat_transp_daily.py > $PBS_JOBID.log


## 	TESTS WITH NEW CODE+BUG FIX: results seems better, longer contour with more details but overall not much difference from previous contour.
# https://github.com/COSIMA/cosima-recipes/pull/286

#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_Cross-1500m_vol_transp_daily_oldCode.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_Cross-1500m_vol_transp_daily_newCode.py > $PBS_JOBID.log

# Repeating daily analyses with new (v2) code for contour:
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_Cross-1500m_temp_daily_v2.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_Cross-1500m_salt_daily_v2.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_Cross-1500m_z_rho_daily_v2.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_Cross-1500m_vol_transp_daily_v2.py > $PBS_JOBID.log

#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_Cross-1500m_Coordinates_daily_v2.py > $PBS_JOBID.log

# fixing Vol. transport not closing to zero when integrated circum-polarly; error where in the interpolation of masks x/y transports; fixed now in v3.0:
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_Cross-1500m_temp_daily_v3.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_Cross-1500m_salt_daily_v3.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_Cross-1500m_z_rho_daily_v3.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_Cross-1500m_vol_transp_daily_v3.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_Cross-1500m_heat_transp_daily_v3.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_Cross-1500m_Coordinates_daily_v3.py > $PBS_JOBID.log

# calculate of heat transport Zonal Convergence:
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_Cross-1500m_zonal_convergence.py > $PBS_JOBID.log

# attempt to input Bin index for ZC calculation:
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_Cross-1500m_zonal_convergence_inputBin.py 6 > $PBS_JOBID.log

# v4 (only for heat transport): subtract Tf HT before extracting along/contour (i.e., in the original model grid) = Mimicing access-om2 script (Wilton's):
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_Cross-1500m_heat_transp_daily_v4.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_newcode_Cross-1500m_heat_transp_daily_v4.py > $PBS_JOBID.log
# corrected z_rho calculation:
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_Cross-1500m_z_rho_daily_v4.py > $PBS_JOBID.log
python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_Cross-1500m_temp_daily_v4.py > $PBS_JOBID.log

## testing with model run using T0=-3.6 -> NO difference in the Huon_temp/Hvom_temp diagnostics.
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10-T0_Cross-1500m_temp_daily_v3.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10-T0_Cross-1500m_salt_daily_v3.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10-T0_Cross-1500m_z_rho_daily_v3.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10-T0_Cross-1500m_vol_transp_daily_v3.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10-T0_Cross-1500m_heat_transp_daily_v3.py > $PBS_JOBID.log

## TESTS with WAOM_V2 with new ROMSIceShelf_devel code ROMS 4.1 (1Q/2024).
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_newcode_Cross-1500m_temp_daily_v3.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_newcode_Cross-1500m_salt_daily_v3.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_newcode_Cross-1500m_z_rho_daily_v3.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_newcode_Cross-1500m_vol_transp_daily_v3.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_newcode_Cross-1500m_heat_transp_daily_v3.py > $PBS_JOBID.log

# year 41:
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_newcode_yr41_Cross-1500m_temp_daily_v3.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_newcode_yr41_Cross-1500m_salt_daily_v3.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_newcode_yr41_Cross-1500m_z_rho_daily_v3.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_newcode_yr41_Cross-1500m_vol_transp_daily_v3.py > $PBS_JOBID.log
#python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/Cross-1500m/WAOM10_newcode_yr41_Cross-1500m_heat_transp_daily_v3.py > $PBS_JOBID.log
