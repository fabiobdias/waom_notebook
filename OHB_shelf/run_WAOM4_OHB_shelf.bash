#!/bin/bash
#PBS -P e14
#PBS -q normal
#PBS -l ncpus=48
#PBS -l mem=190GB
#PBS -l walltime=48:00:00
#PBS -l software=python
#PBS -l wd
#PBS -l storage=gdata/hh5+scratch/gh9
###PBS -l storage=scratch/gh9

module use /g/data/hh5/public/modules
module load conda/analysis3-unstable

# Run Python applications
python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/WAOM4extend_shflim_S_0.25Q_OHB_shelf.py > $PBS_JOBID.log
##python3 /g/data3/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook/OHB_shelf/WAOM4extend_shflim_S_0.25Q_OHB_shelf_save_tmp_files.py > $PBS_JOBID.log


~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
"run_WAOM4_WMT.bash" 23L, 1126C                                                                                                                                                                                                                               1,1           All
