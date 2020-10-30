#!/bin/bash


#$ -cwd
# #$ -j y

#$ -o /lustre/scratch/astro/im281/fitIR/log/
#$ -e /lustre/scratch/astro/im281/fitIR/log/

# sumbit with the following command
# qsub -t 1-10 -l m_mem_free=20G run_xidplus_elais_n1.sh
echo starting script
. /etc/profile.d/modules.sh
module load sge

module load easybuild/software
#module load python/intelpython3/3.5.3
module load Anaconda3/4.0.0
echo conda_activate
#source activate /its/home/im281/.conda/envs/herschelhelp
source activate herschelhelp

cd /lustre/scratch/astro/im281/fitIR/examples
# export PATH="/its/home/im281/.conda/envs/herschelhelp/bin/python":$PATH
echo $SGE_TASK_ID
export HDF5_USE_FILE_LOCKING='FALSE'
echo running_python
python fit_greybody_HELP_z_delta.py #> /lustre/scratch/astro/im281/fitIR/log/$SGE_TASK_ID.txt
#/its/home/im281/.conda/envs/herschelhelp/bin/python3.6 fit_greybody_HELP.py
# /its/home/im281/.conda/envs/herschelhelp/bin/python3.6 XID+_run_ELAIS-N1_PACS.py
echo finished

