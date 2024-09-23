import pandas as pd
import numpy as np
import scipy as sp
import pickle
import itertools
import logging
from optparse import OptionParser
import sys
from time import time
from sklearn.pipeline import make_pipeline
from itertools import product
from time import localtime
import ast


#############################################################################################################################
#############################################################################################################################
start=time()
usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)
parser.add_option("-c", "--color_palette_pickle", dest='color_palette_pickle', type="string",
                 help="filename of the color palette .p file output by generalized_color_metrics_alltaxa_results_compiler.py assumes current directory")
parser.add_option("-d", "--metrics", dest='metrics', type="string",
                 help="filename of the .csv file in the current directory that lists the complete filenames of the metric datasets to be joined together, \
                 under a column titled 'Metric_datasets'")
parser.add_option("-m", "--metadata", dest='metadata', type="string",
                 help="filename of the .csv file in the current directory that lists the complete filenames of the metadata datasets that will be intersected\
                 via 'Specimen' values with the compiled data. Requires 'Metadata_datasets' 'Grouping_variables' 'ID_variable' 'Taxon_group' columns.\
                 which contain metadata filenames, a list of user supplied grouping variables in the format ['var1','var2',...], the unique species ID variable name,\
                 and the taxon group name, respectively")
parser.add_option("-s", "--save_overall_traits", dest='save_overall_traits', type="string",
                 help="Save the overall trait database before splitting by metadata (can be large)")

(options, args) = parser.parse_args()


###Inputs
color_palette_pickle=options.color_palette_pickle
metrics=options.metrics
metadata=options.metadata
save_overall_traits=options.save_overall_traits


color_compiled_wide=pickle.load(open('%s'%(color_palette_pickle),"rb"))

metrics_datasets=pd.read_csv('%s' %(metrics))

master_metrics=color_compiled_wide
for i in range(0,len(metrics_datasets.Metric_datasets)):
    print(metrics_datasets.Metric_datasets[i])
    Dataset=pd.read_csv(metrics_datasets.Metric_datasets[i])
    Dataset.Specimen=Dataset.Specimen.str.upper()
    Dataset.insert(1,'UID',Dataset.Specimen+'_'+Dataset.Wing+'_'+Dataset.Side)
    Dataset.insert(1,'UID2',Dataset.Specimen+'_'+Dataset.body_part+'_'+Dataset.Side)
    Dataset=Dataset.set_index('UID2')
    master_metrics=master_metrics.join(Dataset.iloc[:,6:])
master_metrics=master_metrics.fillna(0)

##Making the average reflectance ("melanism") across all bands and for RGB only
def melanism(dataframe,thresh,UV=False,B=False,G=False,R=False,FR=False,nIR=False):
    if UV==True:
        UV=dataframe['%s_UV_refl_cm2' %(thresh)]
    else:
        UV=pd.Series()
    if B==True:
        B=dataframe['%s_B_refl_cm2' %(thresh)]
    else:
        B=pd.Series()
    if G==True:
        G=dataframe['%s_G_refl_cm2' %(thresh)]
    else:
        G=pd.Series()
    if R==True:
        R=dataframe['%s_R_refl_cm2' %(thresh)]
    else:
        R=pd.Series()
    if FR==True:
        FR=dataframe['%s_740_refl_cm2' %(thresh)]
    else:
        FR=pd.Series()
    if nIR==True:
        nIR=dataframe['%s_940_refl_cm2' %(thresh)]
    else:
        nIR=pd.Series()
    l=[UV,B,G,R,FR,nIR]
    join=pd.concat([i for i in l if len(i)==len(dataframe)],axis=1)
    return join.mean(axis=1)

master_metrics['1_avg_reflectance_UV940']=melanism(master_metrics,1,UV=True,B=True,G=True,R=True,FR=True,nIR=True)
master_metrics['5_avg_reflectance_UV940']=melanism(master_metrics,5,UV=True,B=True,G=True,R=True,FR=True,nIR=True)
master_metrics['10_avg_reflectance_UV940']=melanism(master_metrics,10,UV=True,B=True,G=True,R=True,FR=True,nIR=True)
master_metrics['20_avg_reflectance_UV940']=melanism(master_metrics,20,UV=True,B=True,G=True,R=True,FR=True,nIR=True)

master_metrics['1_avg_reflectance_UV740']=melanism(master_metrics,1,UV=True,B=True,G=True,R=True,FR=True)
master_metrics['5_avg_reflectance_UV740']=melanism(master_metrics,5,UV=True,B=True,G=True,R=True,FR=True)
master_metrics['10_avg_reflectance_UV740']=melanism(master_metrics,10,UV=True,B=True,G=True,R=True,FR=True)
master_metrics['20_avg_reflectance_UV740']=melanism(master_metrics,20,UV=True,B=True,G=True,R=True,FR=True)

master_metrics['1_avg_reflectance_UVR']=melanism(master_metrics,1,UV=True,B=True,G=True,R=True)
master_metrics['5_avg_reflectance_UVR']=melanism(master_metrics,5,UV=True,B=True,G=True,R=True)
master_metrics['10_avg_reflectance_UVR']=melanism(master_metrics,10,UV=True,B=True,G=True,R=True)
master_metrics['20_avg_reflectance_UVR']=melanism(master_metrics,20,UV=True,B=True,G=True,R=True)

master_metrics['1_avg_reflectance_BR']=melanism(master_metrics,1,B=True,G=True,R=True)
master_metrics['5_avg_reflectance_BR']=melanism(master_metrics,5,B=True,G=True,R=True)
master_metrics['10_avg_reflectance_BR']=melanism(master_metrics,10,B=True,G=True,R=True)
master_metrics['20_avg_reflectance_BR']=melanism(master_metrics,20,B=True,G=True,R=True)

cols=[col for col in master_metrics.columns if ('Color_richness' in col) or ('Color_contrast' in col) ]
for i in cols:
    print(i)
    master_metrics['%s_SC'%(i)]=(master_metrics['%s'%(i)]/master_metrics.Area_Mask_cm2)*master_metrics.Area_Mask_cm2.mean()

time_struct=localtime()

if save_overall_traits=='Yes':
    pickle.dump(master_metrics, open("compiled_reflectance_color_complexity_metrics_%s_%s_%s.p" %(time_struct[2],time_struct[1],time_struct[0]), "wb"))
    master_metrics.to_csv("compiled_reflectance_color_complexity_metrics_%s_%s_%s.csv" %(time_struct[2],time_struct[1],time_struct[0]))
else:
    print('overall compiled traits are not being saved')

#########Now splitting by metadata
print('Now splitting by metadata')
master_metrics.insert(0,'UID3',master_metrics.Specimen+'_'+master_metrics.Wing+'_'+master_metrics.Side)
master_metrics=master_metrics.set_index('UID3',drop=False)

import ast
metadata_datasets=pd.read_csv('%s' %(metadata),converters={"Grouping_variables": ast.literal_eval})

lista=['UID3','Specimen','Sex','Side','Wing']
def metadata_join_filter_group(overall_dataset,metadata_dataset, sp_id_str,fixed_variables,user_variables):
    overall_dataset1=overall_dataset.set_index('Specimen')
    overall_dataset2=metadata_dataset.join(overall_dataset1,how='inner')
    overall_dataset2=overall_dataset2.set_index('UID3',drop=False)
    overall_dataset3=overall_dataset2[overall_dataset2.body_part!='body']
    overall_dataset3=overall_dataset3[overall_dataset3.body_part!='antenna']
    overall_dataset_grouped=overall_dataset3.reset_index(drop=True).groupby(fixed_variables+['%s' %(sp_id_str)]+user_variables,axis=0).mean(numeric_only=True).reset_index().set_index('UID3',drop=False)
    overall_dataset_grouped['%s' %(sp_id_str)]=overall_dataset_grouped['%s' %(sp_id_str)].fillna('missing')
    overall_dataset_grouped=overall_dataset_grouped.fillna(0)
    return(overall_dataset_grouped)

if len(metadata_datasets.Metadata_datasets)>0:
    metadata_A=pd.read_csv(metadata_datasets.Metadata_datasets[0])
    metadata_A=metadata_A.drop_duplicates('Specimen').set_index('Specimen',drop=False).dropna(subset=[metadata_datasets.ID_variable[0]])
    metadata_A.Specimen=metadata_A.Specimen.str.upper()
    master_metrics_A=metadata_join_filter_group(master_metrics,metadata_A,'%s' %(metadata_datasets.ID_variable[0]),lista,metadata_datasets.Grouping_variables[0])
    pickle.dump(master_metrics_A, open("Master_%s_trait_dataset_%s_%s_%s.p" %(metadata_datasets.Taxon_group[0],time_struct[2],time_struct[1],time_struct[0]), "wb"))
    master_metrics_A.to_csv("Master_%s_trait_dataset_%s_%s_%s.csv" %(metadata_datasets.Taxon_group[0],time_struct[2],time_struct[1],time_struct[0]),index=False)
    if len(metadata_datasets.Metadata_datasets)>1:
        metadata_B=pd.read_csv(metadata_datasets.Metadata_datasets[1])
        metadata_B=metadata_B.drop_duplicates('Specimen').set_index('Specimen',drop=False).dropna(subset=[metadata_datasets.ID_variable[1]])
        metadata_B.Specimen=metadata_B.Specimen.str.upper()
        master_metrics_B=metadata_join_filter_group(master_metrics,metadata_B,'%s' %(metadata_datasets.ID_variable[1]),lista,metadata_datasets.Grouping_variables[1])
        pickle.dump(master_metrics_B, open("Master_%s_trait_dataset_%s_%s_%s.p" %(metadata_datasets.Taxon_group[1],time_struct[2],time_struct[1],time_struct[0]), "wb"))
        master_metrics_B.to_csv("Master_%s_trait_dataset_%s_%s_%s.csv" %(metadata_datasets.Taxon_group[1],time_struct[2],time_struct[1],time_struct[0]),index=False)
        if len(metadata_datasets.Metadata_datasets)>2:
            metadata_C=pd.read_csv(metadata_datasets.Metadata_datasets[2])
            metadata_C=metadata_C.drop_duplicates('Specimen').set_index('Specimen',drop=False).dropna(subset=[metadata_datasets.ID_variable[2]])
            metadata_C.Specimen=metadata_C.Specimen.str.upper()
            master_metrics_C=metadata_join_filter_group(master_metrics,metadata_C,'%s' %(metadata_datasets.ID_variable[2]),lista,metadata_datasets.Grouping_variables[2])
            pickle.dump(master_metrics_C, open("Master_%s_trait_dataset_%s_%s_%s.p" %(metadata_datasets.Taxon_group[2],time_struct[2],time_struct[1],time_struct[0]), "wb"))
            master_metrics_C.to_csv("Master_%s_trait_dataset_%s_%s_%s.csv" %(metadata_datasets.Taxon_group[2],time_struct[2],time_struct[1],time_struct[0]),index=False)

stop=time()
print("elapsed time:",stop-start)
print("")
