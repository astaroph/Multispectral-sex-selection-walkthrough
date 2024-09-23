#!/bin/bash
#SBATCH -p  serial_requeue # Partition to submit to
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-16:00 # Runtime in days-hours:minutes
#SBATCH --mem 16000 # Memory in MB
#SBATCH -o /path/to/Generalized_color_metrics_py%A.out # File to which standard out will be written
#SBATCH -e /path/to/Generalized_color_metrics_py%A.err # File to which standard err will be written
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=EMAIL # Email to which notifications will be sent

# #arg1
# matrices_dir_path='/path/to/dir/'
# #arg2
# matrices_file='/path/to/specimen_matrices_list.txt'
# #arg3
# body_part='1'
# #arg4
# channel_start=0
# #arg5
# channel_stop=6
# #arg6
# startprop=0
# #arg7
# endprop=1
# #arg8
# stepprop=0.05
##Note, see the script for description of the body parts and channels
## As of 5-4-21, we need to run this separately for body parts 1-6 (all), and for the following channel ranges
#0-6 (UV,B,G,R,740,940)-
#1-4 (B,G,R)-
#0-4 (UV,B,G,R)-
#0-5 (UV,B,G,R,740)-
#5-6 (940)

module load python
source activate /path/to/environment/mypython3_R8
module load gcc/12.2.0-fasrc01
cd ~/workingdir/

python generalized_color_metrics_alltaxa.py -d ${1} \
-m ${2} \
-b ${3} \
-f ${4} \
-l ${5} \
-s 0 \
-e 1 \
-t 0.05 \
-T ${6}

