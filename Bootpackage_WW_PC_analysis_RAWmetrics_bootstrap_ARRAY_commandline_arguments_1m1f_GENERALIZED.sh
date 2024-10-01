#!/bin/bash
#SBATCH -p  serial_requeue # Partition to submit to
#SBATCH -n 6 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-00:55 # Runtime in days-hours:minutes
#SBATCH --mem 48000 # Memory in MB
#SBATCH -o /PATH/TO/Bootpackage_WW_PC_analysis_RAWmetrics_bootstrap_ARRAY_commandline_arguments_1m1f_GENERALIZED%A_%a.out # File to which standard out will be written
#SBATCH -e /PATH/TO/Bootpackage_WW_PC_analysis_RAWmetrics_bootstrap_ARRAY_commandline_arguments_1m1f_GENERALIZED%A_%a.err # File to which standard err will be written
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=EMAIL # Email to which notifications will be sent

module load R/4.2.2-fasrc01
module load gcc/12.2.0-fasrc01
export R_LIBS_USER=/PATH/TO/R_4.2.2/
cd ~/workingdir/

#supply arguments following the sbatch command, a single-quoted
#taxon name (one of the options given in the R script) as the first: 
#choices are 'Butterflies' 'Moths' and, experimentally 'Geometrids'

# the categorical variable column name exactly as it is in the dataset as the second
# a numeric batch size variable, which is the number of bootstraps per array task
# and the full path to an Rdata file which contains all of the initial part of the function up to the bootstrap for-loop.
# and the minimum number of specimens per species threshold supplied as a single unquoted int
#NOTE: if you don't have an Rdata file to supply, just put 'No_Data' which will instruct the program to first do all of the initial steps to generate the PC rates and per specimen
##datasets

R CMD BATCH --quiet --no-save --no-restore --no-echo "--args metric_num=${SLURM_ARRAY_TASK_ID} Taxon_group='${1}' categorical_variable='${2}' boot_reps=${3} rdata='${4}' specimen_threshold=${5} n_cores=${6}" \
/n/home06/astaroph/WingsAndWavelengths/scripts/Boot_package_redo_Moths_AND_butterflies_RAW_metrics_M_minus_F_contrasts_rates_MD_script_1m1f_datatable_cleanversion_sexISdistnull_GENERALIZED.R \
/n/holyscratch01/pierce_lab/astaroph/WW_R_output/Raw_metrics_array/boot_package_permetric/${SLURM_ARRAY_TASK_ID}_${1}_${2}_${5}m${5}f_RAWmetrics_bootpkg_R_outputfile_${3}BS.out
