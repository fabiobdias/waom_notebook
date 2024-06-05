#PBS -P e14
#PBS -q hugemem
#PBS -l ncpus=48
#PBS -l mem=1500GB
#PBS -l walltime=06:00:00
#PBS -l software=python
#PBS -l wd
#PBS -l storage=gdata/hh5+scratch/gh9+scratch/gi0

module use /g/data/hh5/public/modules
module load conda/analysis3-unstable

cd /g/data/hh5/tmp/access-om/fbd581/ROMS/postprocessing/waom_notebook

# Run Python applications
python3 WAOM2_save_dzt_yr5_monthly.py > $PBS_JOBID.log
