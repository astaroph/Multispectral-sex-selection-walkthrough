#!/bin/bash
#SBATCH -p  shared # Partition to submit to
#SBATCH -n 12 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-02:45 # Runtime in days-hours:minutes
#SBATCH --mem 50000 # Memory in MB
#SBATCH -o /path/to/Bootpackage_WW_allPCs_PC_analysis_overall_highboot_bootstrap_ARRAY_commandline_arguments_1m1f_singleband_GENERALIZED%A_%a.out # File to which standard out will be written
#SBATCH -e /path/to/Bootpackage_WW_allPCs_PC_analysis_overall_highboot_bootstrap_ARRAY_commandline_arguments_1m1f_singleband_GENERALIZED%A_%a.err # File to which standard err will be written
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=EMAIL # Email to which notifications will be sent
module load R/4.1.0-fasrc01
module load gcc/9.2.0-fasrc01
export R_LIBS_USER=/path/to/Rlibrary/R_4.1.0_cluster/:$R_LIBS_USER
cd ~/workingdir/

#supply arguments following the sbatch command, a single-quoted
#taxon name (one of the options given in the R script) as the first: 
#choices are 'Butterflies' 'Moths' 'Moths_and_butterflies' but 'Eumaeini' still needs its own code

# the categorical variable column name exactly as it is in the dataset as the second
# a numeric batch size variable, which is the number of bootstraps per array task
# and the full path to an Rdata file which contains all of the initial part of the function up to the bootstrap for-loop.
# and the minimum number of specimens per species threshold supplied as a single unquoted int
#NOTE: if you don't have an Rdata file to supply, just put 'No_Data' which will instruct the program to first do all of the initial steps to generate the PC rates and per specimen
##datasets

R CMD BATCH --quiet --no-save --no-restore --no-echo "--args Taxon_group='${1}' categorical_variable='${2}' boot_reps=${3} rdata='${4}' specimen_threshold=${5} n_cores=${6} rep=${7} wband='${8}'" \
/path/to/Boot_package_Moths_AND_butterflies_PCmetrics_M_minus_F_contrasts_rates_1m1f_datatable_commandline_allPCs_singleband_GENERALIZED.R \
/path/to/outputfiles/AllPCs${1}_${2}_${5}m${5}f_OverallPC_PCmetrics_bootpkg_R_outputfile_${3}BS${7}_${8}.out
