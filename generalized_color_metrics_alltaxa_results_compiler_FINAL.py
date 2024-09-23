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


#############################################################################################################################
#############################################################################################################################
start=time()
usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)
parser.add_option("-d", "--directory", dest='directory', type="string",
                 help="Directory containing color palette .p files")
parser.add_option("-m", "--matrix_var_file", dest='matrix_var_file', type="string",
                 help="The following should read in a user supplied csv file with the column names 'Taxon_groups' 'Body_parts' 'Channel_groups' in that order. All are read as strings though body parts are the numeric codes for the different body parts detected by our pipeline. Channel groups must be specified as 'C[x]_to_X[y]'' where x is the start channel and y is the stop channel output by the color palette script (0=UV to 6=nIR)")

(options, args) = parser.parse_args()


###Inputs
filepath=options.directory
filename=options.matrix_var_file
###The following should read in a user supplied csv file with the column names "Taxon_groups" "Body_parts" "Channel_groups" in that order. All are read as strings though body parts 
##are the numeric codes for the different body parts detected by our pipeline. Channel groups must be specified as "C[x]_to_X[y]" where x is the start channel and y is the stop channel
##output by the color palette script (0=UV to 6=nIR)
iteration_frame=pd.read_csv(filename,dtype=str)
def compile_color_results(Taxon=None,compiled=pd.DataFrame(),picklename=None,body_part=None,channel_col=None,custom_column1_name=None,custom_column1_value=None,new_column_index=0):
    helptext="""Usage: compile_color_results(compiled,picklename,body_part,channel_col,custom_column1_name,custom_column1_value,new_column_index)

    Description:
    Compile color results from different body parts and channel ranges into one dataframe with descriptive metadata columns
    User provides a blank pandas dataframe instance,the pickle .p results file output from the generalized color palette pipeline, 
    the body part and channel range as strings, one custom column name and the value to fill it with as strings,
    and an index value at which the new metadata columns will be added.
    Several identical .p files can be combined by calling the function on the output of this function, to join them together, provided
    that the new file does not have identical values for the metadata or channel range columns.

    Parameters
    ----------
    compiled: a pandas DataFrame, OPTIONAL
        defaults to a blank instance of a pandas DataFrame to populate with the results.p files, but the output of a 
        compile_color_results function call can be supplied by the user to which the function will join new results.
    
    pickle: a pickle file, REQUIRED
        the full path and filename to the raw pickle files output by the generalized_color_metrics.py file, as a string
    
    body_part: string, REQUIRED
        a user provided string describing the body part of the results.p being imported

    channel_col: string, REQUIRED
        a user provided string describing the channel range of the results.p being imported

    custom_column1_name: string, OPTIONAL
        an optional user provided string name for a custom descriptive metadata column. defaults to None

    custom_column1_value: any value, OPTIONAL, required if custom_column1_name is passed
        the value with which to fill the user metadata column

    new_column_index: int, OPTIONAL
        the index integer at which to insert new metadata columns. defaults to 0 (the beginning of the dataframe) 
        if not supplied by user
    Returns
    -------
    pandas DataFrame
    """
    if picklename!=None:
        trial=pickle.load(open(picklename, "rb"))
    else:
        print(helptext)
        return
    if channel_col!=None:
        trial.insert(new_column_index,'Channel_range',channel_col)
    else:
        print(helptext)
        return
    if body_part!=None:
        trial.insert(new_column_index+1,'body_part',body_part)
    else:
        print(helptext)
        return
    if custom_column1_name!=None:
        trial.insert(new_column_index+2,custom_column1_name,custom_column1_value)
    UID_str=trial.UID+'_'+trial.body_part+'_'+trial.Channel_range
    trial.insert(new_column_index+3,'UID2',UID_str)
    trial['UID2']=trial.UID2.str.replace('-','_')
    trial=trial.set_index('UID2',drop=False)
    ##Explicitly coding the numeric columns as numeric. for some reason this isn't automatic and they are loading as objects
    trial[['1cm_scale_pixels','Primary_color_8bit_prop','Secondary_color_8bit_prop',
       'Tertiary_color_8bit_prop', 'Primary_color_3bit_prop', 'Secondary_color_3bit_prop', 'Tertiary_color_3bit_prop',
   'MDI_8bit', 'MDI_3bit','Color_richness_8bit', 'Color_richness_3bit','Norm_perchannel_Color_richness_8bit',
       'Norm_perchannel_Color_richness_3bit', 'SW_entropy_8bit',
       'SW_entropy_3bit', 'Perchannel_effective_num_colors_8bit',
       'Perchannel_effective_num_colors_3bit',
       'Norm_actual_to_effective_colors_8bit',
       'Norm_actual_to_effective_colors_3bit', 'SW_evenness_8bit',
       'SW_evenness_3bit', 'Color_contrast_8bit', 'Color_contrast_3bit']]= trial[['1cm_scale_pixels','Primary_color_8bit_prop','Secondary_color_8bit_prop',
       'Tertiary_color_8bit_prop', 'Primary_color_3bit_prop', 'Secondary_color_3bit_prop', 'Tertiary_color_3bit_prop',
   'MDI_8bit', 'MDI_3bit','Color_richness_8bit', 'Color_richness_3bit','Norm_perchannel_Color_richness_8bit',
       'Norm_perchannel_Color_richness_3bit', 'SW_entropy_8bit',
       'SW_entropy_3bit', 'Perchannel_effective_num_colors_8bit',
       'Perchannel_effective_num_colors_3bit',
       'Norm_actual_to_effective_colors_8bit',
       'Norm_actual_to_effective_colors_3bit', 'SW_evenness_8bit',
       'SW_evenness_3bit', 'Color_contrast_8bit', 'Color_contrast_3bit']].apply(pd.to_numeric)
    trial.insert(new_column_index+4,'Taxon_group',Taxon)
    if len(compiled)>0:
        compiled=pd.concat([compiled,trial],axis=0)
    else:
        compiled=trial.join(compiled)
    return compiled
color_compiled=pd.DataFrame()
def body_call_name(str):
    bod_dict={'1':'LFW','2':'RFW','3':'LHW','4':'RHW','5':'body','6':'antenna'}
    return(bod_dict[str])

def channel_call_name(str):
    channel_dict={'C0_to_C1': 'UV', 'C0_to_C2': 'UV-B', 'C0_to_C3': 'UV-G', 'C0_to_C4': 'UV-R', 'C0_to_C5': 'UV-740', 'C0_to_C6': 'UV-940', 'C1_to_C2': 'B', 'C1_to_C3': 'B-G', 'C1_to_C4': 'B-R',
 'C1_to_C5': 'B-740', 'C1_to_C6': 'B-940', 'C2_to_C3': 'G', 'C2_to_C4': 'G-R', 'C2_to_C5': 'G-740', 'C2_to_C6': 'G-940', 'C3_to_C4': 'R', 'C3_to_C5': 'R-740', 'C3_to_C6': 'R-940',
 'C4_to_C5': '740', 'C4_to_C6': '740-940', 'C5_to_C6': '940'}
    return(channel_dict[str])

for ind in iteration_frame.index:
    call=iteration_frame.loc[ind,]
    channel_col=channel_call_name(call[2])
    body_part=body_call_name(call[1])
    custom_column1_name='Wing'
    new_column_index=4
    if body_part=='LFW' or body_part=='RFW':
        custom_column1_value='FW'
    elif body_part=='LHW' or body_part=='RHW':
        custom_column1_value='HW'
    else:
        custom_column1_value='none'
    picklename='%s%s_All_in_one_NClist_colors_MDI_streamlined_bit_color_results_BodyPart%s_%s.p' %(filepath,call[0],call[1],call[2])
    try:
        color_compiled=compile_color_results(Taxon=call[0],compiled=color_compiled,picklename=picklename,body_part=body_part,channel_col=channel_col,custom_column1_name=custom_column1_name,custom_column1_value=custom_column1_value,new_column_index=new_column_index)
    except:
        print('%s NOT PRESENT' %(picklename))
        continue

        
color_compiled.Channel_range=color_compiled.Channel_range.str.replace('-','_')
color_compiled2=color_compiled.drop('UID2',axis=1)
color_compiled_wide=color_compiled2.set_index(['ID','UID','Taxon_group','Specimen','Side','body_part','Wing','Channel_range']).unstack()
color_compiled_wide.columns = color_compiled_wide.columns.map(lambda x: '{}_{}'.format(x[0], x[1]))
color_compiled_wide = color_compiled_wide.reset_index()
color_compiled_wide.Specimen=color_compiled_wide.Specimen.str.upper()
color_compiled_wide.insert(0,'UID2',color_compiled_wide.Specimen+'_'+color_compiled_wide.body_part+'_'+color_compiled_wide.Side)
color_compiled_wide=color_compiled_wide.set_index('UID2')
# color_compiled_wide=color_compiled_wide.drop(['1cm_scale_pixels_940','1cm_scale_pixels_B_R','1cm_scale_pixels_UV_740','1cm_scale_pixels_UV_R'],axis=1)

time_struct=localtime()
pickle.dump(color_compiled_wide, open("Compiled_color_palette_metrics_%s_%s__%s_%s_%s.p" %(time_struct[3],time_struct[4],time_struct[2],time_struct[1],time_struct[0]), "wb"))
color_compiled_wide.to_csv("Compiled_color_palette_metrics_%s_%s__%s_%s_%s.csv" %(time_struct[3],time_struct[4],time_struct[2],time_struct[1],time_struct[0]))


stop=time()
print("elapsed time:",stop-start)
print("")
