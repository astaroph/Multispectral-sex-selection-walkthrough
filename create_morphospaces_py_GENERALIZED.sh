#!/bin/bash
#SBATCH -p  serial_requeue # Partition to submit to
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-01:00 # Runtime in days-hours:minutes
#SBATCH --mem 30000 # Memory in MB
#SBATCH -o /path/to/create_morphospaces_py%A.out # File to which standard out will be written
#SBATCH -e /path/to/create_morphospaces_py%A.err # File to which standard err will be written
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=EMAIL # Email to which notifications will be sent

module load python
source activate /path/to/mypython3_ENV
module load gcc/12.2.0-fasrc01
cd ~/workingdir/

# parser.add_option("-d", "--dataset_str", dest='dataset_str', type="string",
#                  help="input dataset filename")
# parser.add_option("-t", "--Taxon", dest='Taxon', type="string",default='Butterflies,'
#                  help="Taxon group corresponding to dataset, can be 'Butterflies' or 'Moths' currently")
# parser.add_option("-c", "--cat_var", dest='cat_var', type="string",
#                  help="Categorical variable to include in final PC datasets and to optionally filter by")
# parser.add_option("-f", "--filter_str", dest='filter_str', type="string",
#                  help="Restrict dataset to those matching this string in cat_var")
# parser.add_option("-o", "--output_prefix", dest='output_prefix', type="string",
#                  help="unique prefix to append to filenames")
# parser.add_option("-b", "--traits_begin", dest='traits_begin', type="int",default=20,
#                  help="where traits begin in the trait dataset")
# parser.add_option("-i", "--ID_var", dest='ID_var', type="string",default='Phy_ID',
#                  help="string name of the unique species ID column")
# parser.add_option("-s", "--specimen_cutoff", dest='specimen_cutoff', type="int",default=1,
#                  help="Minimum number of male and female specimens a species must have to be retained before PCA")
# parser.add_option("-F", "--output_mf_filt", dest='output_mf_filt', action="store_true",default=True,
#                  help="Boolean. whether to save the raw trait dataset that has been filtered by the specimen threshold. Defaults to True")
# parser.add_option("-A", "--output_All_PCs", dest='output_All_PCs',action="store_true", default=True,
#                  help="Boolean, whether to save the PC dataset that contains all PCs. Defaults to True")
# parser.add_option("-N", "--output_topN_PCs", dest='output_topN_PCs',action="store_true", default=True,
#                  help="Boolean, whether to save the PC dataset that contains only the PCs explaining the top -p proportion of variance. True by default ")
# parser.add_option("-p", "--topN_PC_prop", dest='topN_PC_prop', type="float",default=.90,
#                  help="The threshold proportion of variance explained by the top N PCs")
# parser.add_option("-P", "--plot_PCA", dest='plot_PCA',action="store_true",default=True,
#                  help="Boolean, whether to plot the PC plot. Defaults to true")
# parser.add_option("-x", "--PC1", dest='PC1', type="int",default=1,
#                  help="PC to set on the x axis")
# parser.add_option("-y", "--PC2", dest='PC2', type="int",default=2,
#                  help="PC to set on the y axis")
# parser.add_option("-V", "--variable", dest='variable', type="string",default='Family',
#                  help="Grouping variable for PC plot")
# parser.add_option("-S", "--shapevar", dest='shapevar', type="string",
#                  help="Point shape variable for PC plot")
# parser.add_option("-O", "--open_closed_var", dest='open_closed_var', type="string",
#                  help="Point open or closed variable for PC plot")
# parser.add_option("-C", "--create_singleband_morphospaces", dest='create_singleband_morphospaces', action="store_true",default=True,
#                  help="Whether to create separate morphospaces for the reflectance based single band metrics datasets")

###Note, do not supply the following lines of code, and they will default to None type
##If you supply None as a value it fails
# --cat_var None \
# --filter_str None \


python create_morphospaces_GENERIC.py --dataset_str 'Master_Moths_trait_dataset_9_2_2023_sensitivitycut.csv' \
--Taxon 'Moths' \
--cat_var 'Diel_behavior' \
--filter_str 'Nocturnal' \
--output_prefix 'Nocturnal_Moths_1m1f_sensitivitycut_Final' \
--traits_begin 23 \
--ID_var 'Moth_ID' \
--specimen_cutoff 1 \
--output_mf_filt \
--output_All_PCs \
--output_topN_PCs \
--topN_PC_prop .90 \
--plot_PCA \
--PC1 1 \
--PC2 2 \
--variable 'Family' \
--shapevar 'Sex_Side_Wing' \
--open_closed_var 'Sex' \
--create_singleband_morphospaces

python create_morphospaces.py --dataset_str 'Master_Butterflies_trait_dataset_8_2_2023_sensitivitycut.csv' \
--Taxon 'Butterflies' \
--output_prefix 'Butterflies_1m1f_sensitivitycut_Final' \
--traits_begin 20 \
--ID_var 'Phy_ID' \
--specimen_cutoff 1 \
--output_mf_filt \
--output_All_PCs \
--output_topN_PCs \
--topN_PC_prop .90 \
--plot_PCA \
--PC1 1 \
--PC2 2 \
--variable 'Family' \
--shapevar 'Sex_Side_Wing' \
--open_closed_var 'Sex' \
--create_singleband_morphospaces
