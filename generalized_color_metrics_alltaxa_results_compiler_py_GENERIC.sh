#!/bin/bash
#SBATCH -p  serial_requeue # Partition to submit to
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-01:00 # Runtime in days-hours:minutes
#SBATCH --mem 30000 # Memory in MB
#SBATCH -o /path/to/generalized_color_metrics_alltaxa_results_compiler_py%A.out # File to which standard out will be written
#SBATCH -e /path/to/generalized_color_metrics_alltaxa_results_compiler_py%A.err # File to which standard err will be written
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=EMAIL # Email to which notifications will be sent

source new-modules.sh
module load python/3.6.3-fasrc02
source activate mypython3
cd ~/workingdir/

python generalized_color_metrics_alltaxa_results_compiler.py -d '/path/to/Color_results_output/' \
-m 'color_palette_metrics_compilation_input_variable_dataframe.csv' 
