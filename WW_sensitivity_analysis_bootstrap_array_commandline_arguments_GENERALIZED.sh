#!/bin/bash
#SBATCH -p  serial_requeue # Partition to submit to
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-00:10 # Runtime in days-hours:minutes
#SBATCH --mem 1800 # Memory in MB
#SBATCH -o /path/to/WW_sensitivity_analysis_bootstrap_array_commandline_arguments_GENERALIZED%A_%a.out # File to which standard out will be written
#SBATCH -e /path/to/WW_sensitivity_analysis_bootstrap_array_commandline_arguments_GENERALIZED%A_%a.err # File to which standard err will be written
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=EMAIL # Email to which notifications will be sent
module load R/4.0.2-fasrc01
module load gcc/12.1.0-fasrc01
export R_LIBS_USER=/path/to/RLIBRARY/
cd ~/workingdir/

#supply 2 arguments following the sbatch command, a single-quoted
#taxon name (one of the options given in the R script) as the first: 
#choices are 'Butterflies' 'Moths' 'Geometrids' 
#Second argument is an arbitrary offset, usually left as zero unless you wish to run very large numbers of bootstraps, in which case you need to submit
##separate array jobs with unique offets to go beyond the number of jobs you can have in a single array (9999). In my case I ran an array of 1000 with an offset of 0, and an 
##array of 9000 with an offset of 1000, for 10000 unique bootstraps. Note, offsets should be chosen so that all output files have a unique value of "boot" or they will
##overwrite each other.
num=${SLURM_ARRAY_TASK_ID}
boot=$((num+${2}))
R CMD BATCH --quiet --no-save --no-restore --no-echo "--args boot=${boot} Taxon_group='${1}'" \
/path/to/Bootstrap_sensitivity_analysis_butterflies_and_moths_ARRAY_GENERALIZED.R \
/path/to/${1}_sensitivity_analysis_R_outputfile${boot}.out
