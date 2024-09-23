import pandas as pd
import numpy as np
import scipy as sp
import itertools
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import Normalizer
import logging
from optparse import OptionParser
import sys
from time import time
from scipy.io import loadmat
from collections import Counter
import pickle
import math
#############################################################################################################################
#############################################################################################################################

usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)
parser.add_option("-d", "--directory", dest='matrices_dir_path', type="string",
                 help="Directory containing .mat files")
parser.add_option("-m", "--mattxtfile", dest='matrices_file', type="string",
                 help="full path to list of .mat filenames, one per line")
parser.add_option("-b", "--body_part", dest='body_part',type="string",
                 help="string number defining the body part of interest")
parser.add_option("-f", "--channel_start", dest='channel_start',type="int",
                 help="The start channel of a continuous range of image channels to evaluate (inclusive)")
parser.add_option("-l", "--channel_stop", dest='channel_stop',type="int",
                 help="the stop channel of a continuous range of image channels to evaluate (exclusive)")
parser.add_option("-s", "--startprop", dest='startprop',type="float",
                 help="starting proportion for color frequency histograms")
parser.add_option("-e", "--endprop", dest='endprop',type="float",
                 help="ending proportion for color frequency histograms")
parser.add_option("-t", "--stepprop", dest='stepprop',type="float",
                 help="interval proportion for color frequency histograms")
parser.add_option("-T", "--taxon_group", dest='taxon_group', type="string",
                 help="Name of taxon group to label output files")
parser.add_option("-o", "--output_dir", dest='output_dir', type="string",
                 help="Path to directory where output files should be written. Must end with '/'")

(options, args) = parser.parse_args()

taxon_group=options.taxon_group

##matrices_dir is the full path, passed externally, to the directory where the .mat matrix files are contained.
### MATRIX_EXTERNAL_FILENAME_AND_PATH is the full path and filename for a simple list with the filenames of the .mat matrix files 
# which are contained one per line
#it will be read in below and processed into a simple dataframe with a single column 'name.' each of these will be processed jointly
##as a for loop, resulting in a single dataframe for each body part and combination of channels, across all mat files in this text file
matrices_dir_path=options.matrices_dir_path
matrices_file=options.matrices_file
output_dir=options.output_dir
#body part is a number passed as a string that is the last number in the string colors supplied by Wei-Ping Chan's pipeline
#currently:
#0 is background
#1 is left fore wing from dorsal side
#2 is right fore wing from dorsal side
#3 is left hind wing from dorsal side
#4 is right hind wing from dorsal side
#5 is body
#6 is antenna
#7 is uncertain categories
body_part=options.body_part

#Channel refers to the different imaging wavelength bands. This pipeline is designed to take a CONTIGUOUS range of these, from
##0-13, currently supplied as an int
# 1: UV,# 2: B,# 3: G,# 4: R,# 5: 740,# 6: 940,# 7: F (B),# 8: F (G),# 9: F (R),# 10: Pol (B),# 11: Pol (G),# 12: Pol (R),# 13: mask
channel_start=options.channel_start
channel_stop=options.channel_stop

##These variables refer to the start and end proportions (currently from 0 to 1) over which to calculate the color frequency histograms
##these histograms are the number of colors that cover the given proportion of the total body area.
##increasing the startprop increases the starting proportional frequency, whereas decreasing the endprop changes the high proportional frequency
#changing the stepprop changes the bin size of this histogram, with smaller values resulting in finer resolution. currently set to
#5% increments.
startprop=options.startprop
endprop=options.endprop
stepprop=options.stepprop

######Defining functions
##Many of these functions refer to counter object, which contains of all the color labels and their frequencies in a given bodypart

##This function returns the number of colors that account for the given wing area frequency threshold
def color_prop2(counter,threshold):
    labels, values =zip(*counter.most_common()[1:])
    prop_values=[float(i)/np.sum(values) for i in values]
    cumu_values=np.cumsum(prop_values)
    cumu_values1=[x for x in cumu_values if x<threshold]
    return len(cumu_values1)

##This function also takes the same color frequency counter object and returns the color label name and frequency for the supplied
#index, which is 0-indexed from the most common color (0) to the next most common color (1) and so on and so forth.
def top_color(counter,index,threshold=.99):
    labels, values =zip(*counter.most_common()[1:])
    prop_values=[float(i)/np.sum(values) for i in values]
    cumu_values=np.cumsum(prop_values)
    cumu_values1=[x for x in cumu_values if x<threshold]
    labels1 = [x for x, y in zip(labels, cumu_values) if y <threshold]
    prop_values1=[x for x, y in zip(prop_values, cumu_values) if y <threshold]
    if len(labels1)>=index+1:
        return labels1[index],prop_values1[index]
    else:
        return (np.nan,np.nan)

#This function takes the original input string from Wei-Ping Chan's pipeline and if the string matches the given body part
#chops it to the specified channels. If it isnt, it makes it a none value. Meant to be passed to the forloop below to restrict the
##analysis to the given channels and bodyparts.
##NOTE, body_part, channel_start, and channel_stop are required components of this function, but because I can't use apply
##element-wise to the whole dataframe, I have defined these variables globally and continued to use applymap
def paredown_general(string):
    if string[-1:]==body_part:
        return string[channel_start:channel_stop]
    else:
        return None

##Defining a function to relabel 8 bit colors into "3-bit" low high none colors, where "1" is 1, 2-5 are 2 and 6-8 are 3
##This is to get a coarser view of the way colors are distributed.
def threshold_relabel(string):
    if pd.isnull(string)==False:
        inputstr=str(string)
        if len(inputstr)>=1:
            trialstr=list(inputstr)
            for i in range(len(trialstr)):
                if trialstr[i]!='1':
                    if trialstr[i]=="2" or trialstr[i]=="3" or trialstr[i]=="4" or trialstr[i]=="5":
                        trialstr[i]="2"
                    else:
                        trialstr[i]="3"
                else:
                    trialstr[i]="1"
            outputstr="".join(trialstr)
        else:
            outputstr=string
    else:
        outputstr=string
    return outputstr

###This function takes the same counter and makes a color frequency histogram given the supplied values for starting and stopping proportions
def prop_color_histogram3(counter,thresh_start,thresh_stop,step):
    num_colors=[]
    for i in np.arange(thresh_start,thresh_stop,step):
        num_colors+=[color_prop2(counter,i)]
    return num_colors,np.arange(thresh_start,thresh_stop,step)

#This function is a convenience function for finding the closest value in an array to a supplied value, used in the MDI function below
def find_nearest(array, value): array = np.asarray(array); idx = (np.abs(array - value)).argmin(); return idx

##This function calculates the Moment distance index statistic, which is a way of describing the shape of histograms.
#As applied in this pipeline, this differentiates butterflies who have differently shaped color-frequency histograms
def Moment_Distance_Index(x_list,y_list,xleft,xright):
    MDlp=[]
    MDrp=[]
    left_idx=find_nearest(x_list,xleft)
    right_idx=find_nearest(x_list,xright)
    for i,j in zip(y_list[left_idx:right_idx+1],x_list[left_idx:right_idx+1]):
        MDlp+=[(i**2+(j-x_list[left_idx])**2)**0.5]
        MDrp+=[(i**2+(x_list[right_idx]-j)**2)**0.5]
    return sum(MDrp)-sum(MDlp)

#This function calculates the Shannon-Weaver entropy and evenness values for the given colors and frequencies in the counter object
def SW_and_evenness(counter):
    labels, values =zip(*counter.most_common()[1:])
    if len(values)>0:
        prop_values=[float(i)/np.sum(values) for i in values]
        hcol = -sum([prop*math.log(prop) for prop in prop_values]) #Shannon-Wiener index of diversity (also known as image entropy)
        evenness=hcol/math.log(len(values)); #Pielou's evenness index
        return hcol, evenness
    else:
        return None,None
##this function turns a string name into a list of ints for manipulation by other functions.
def listify(string):
    if pd.isnull(string)==False:
        trialist=[int(x) for x in list(string)]
        return trialist
    else:
        return None

##Importing the list of matrix files to process
matrices=pd.read_csv(matrices_file,header=None)
matrices.columns=["name"]
matrices=matrices.set_index('name',drop=False)
##Setting up a results frame container to add values to
Results_frame_base=pd.DataFrame()
###The ID column is the only one of these metadata columns to be invoked in the code below. Specimen and side could be removed for making this more general
Results_frame_base['ID']=matrices.name
Results_frame_base['UID']=matrices.name.str.split('_',expand=True)[0]+'-'+matrices.name.str.split('_',expand=True)[1]
Results_frame_base['Specimen']=matrices.name.str.split('_',expand=True)[0]
Results_frame_base['Side']=matrices.name.str.split('_',expand=True)[1]
#the following metrics are the basic metric supplied by this color palette analysis. Currently these are separately calculated for 
##the full "8bit" palette with 8 brightness steps + a zero, and a none, low, high palette (3bit) that provide different levels of resolution
##Note, thhe 3 bit values may not be very informative for single channels with low variation, and may even break the code if there is no variation in colors.
Results_frame_base['1cm_scale_pixels']=None
Results_frame_base['Primary_color_8bit']=None
Results_frame_base['Primary_color_8bit']=None
Results_frame_base['Primary_color_8bit_prop']=None
Results_frame_base['Secondary_color_8bit']=None
Results_frame_base['Secondary_color_8bit_prop']=None
Results_frame_base['Tertiary_color_8bit']=None
Results_frame_base['Tertiary_color_8bit_prop']=None
Results_frame_base['Primary_color_3bit']=None
Results_frame_base['Primary_color_3bit_prop']=None
Results_frame_base['Secondary_color_3bit']=None
Results_frame_base['Secondary_color_3bit_prop']=None
Results_frame_base['Tertiary_color_3bit']=None
Results_frame_base['Tertiary_color_3bit_prop']=None
Results_frame_base['MDI_8bit']=None
Results_frame_base['MDI_3bit']=None
Results_frame_base['Color_richness_8bit']=None
Results_frame_base['Color_richness_3bit']=None
Results_frame_base['Norm_perchannel_Color_richness_8bit']=None
Results_frame_base['Norm_perchannel_Color_richness_3bit']=None
Results_frame_base['SW_entropy_8bit']=None
Results_frame_base['SW_entropy_3bit']=None
Results_frame_base['Perchannel_effective_num_colors_8bit']=None
Results_frame_base['Perchannel_effective_num_colors_3bit']=None
Results_frame_base['Norm_actual_to_effective_colors_8bit']=None
Results_frame_base['Norm_actual_to_effective_colors_3bit']=None
Results_frame_base['SW_evenness_8bit']=None
Results_frame_base['SW_evenness_3bit']=None
Results_frame_base['Color_centroid_8bit']=None
Results_frame_base['Color_centroid_3bit']=None
Results_frame_base['Color_contrast_8bit']=None
Results_frame_base['Color_contrast_3bit']=None
Results_frame_base=Results_frame_base.set_index('ID',drop=False)
Results_frame_bodypart=Results_frame_base.copy()

start=time()
count=0
names_8bit=[]
list_o_lists_8bit=[]
names_3bit=[]
list_o_lists_3bit=[]

num_channels=channel_stop-channel_start
def forloop_frameOPS(pixel_frame,i,Results_frame,names_8bit,list_o_lists_8bit,names_3bit,list_o_lists_3bit):
    #'A function that takes a raw pixel frame array, the specified body part string number, the starting and ending channel numbers, as well as start stop and step proportion values and populates with values an empty results frame and empty list containers')
    trial=pixel_frame.applymap(paredown_general)
    aggregate_counter=Counter(trial.values.flatten())
    names_8bit+=[i]
    list_o_lists_8bit+=[prop_color_histogram3(aggregate_counter,startprop,endprop,stepprop)[0]]
    topC,topP=top_color(aggregate_counter,0)
    Results_frame.at[i,'Primary_color_8bit']=topC
    Results_frame.at[i,'Primary_color_8bit_prop']=topP
    topC,topP=top_color(aggregate_counter,1)
    Results_frame.at[i,'Secondary_color_8bit']=topC
    Results_frame.at[i,'Secondary_color_8bit_prop']=top_color(aggregate_counter,1)[1]
    topC,topP=top_color(aggregate_counter,2)
    Results_frame.at[i,'Tertiary_color_8bit']=topC
    Results_frame.at[i,'Tertiary_color_8bit_prop']=topP
    colors,props=prop_color_histogram3(aggregate_counter,startprop,endprop,stepprop)
    MDI=Moment_Distance_Index(props,colors,startprop,endprop)
    Results_frame.at[i,'MDI_8bit']=MDI
    richness=color_prop2(aggregate_counter,0.99)
    Results_frame.at[i,'Color_richness_8bit']=richness
    Results_frame.at[i,'Norm_perchannel_Color_richness_8bit']=richness**(1/num_channels)
    SW,evenness=SW_and_evenness(aggregate_counter)
    Results_frame.at[i,'SW_entropy_8bit']=SW
    Results_frame.at[i,'Perchannel_effective_num_colors_8bit']=(math.exp(SW))**(1/num_channels)
    Results_frame.at[i,'Norm_actual_to_effective_colors_8bit']=(richness**(1/num_channels))/((math.exp(SW))**(1/num_channels))
    Results_frame.at[i,'SW_evenness_8bit']=evenness
    ###now for color contrast calculations
    ##First remove nones from the color counter
    #then make a list of lists by converting each str color labels into a list of ints, corresponding to the brightness in each wavelength band.
    #next create a weighting list, which is the occurrence of each color in the image in number of pixels
    #multiply the array created by your list of color int-lists by their respective weight to obtain an array of color ints weighted by their occurence
    ##taking the sum of these weighted color ints divided by the total pixels in the dataset provides a list of ints of the "average" unweighted int for each wavelength band
    ##Creating a weighted centroid dataframe, which is the centroid color multiplied by their ocurrences, allows us to simply manipulate
    ##these arrays together, essentially comparing an array of weighted actual colors against an array of the centroid color weighted by the original color ocurrence data
    del aggregate_counter[None]
    color_lists=[listify(x) for x in list(aggregate_counter.keys())]
    weight_mat=list(aggregate_counter.values())
    weighted_colors=np.array([np.array(x).dot(y) for x, y in zip(color_lists, weight_mat)])
    color_centroid=np.sum(weighted_colors,axis=0)/np.sum(weight_mat)
    Results_frame.at[i,'Color_centroid_8bit']=list(color_centroid)
    weighted_centroid=np.array([[x] for x in color_centroid]).dot(np.array([[x] for x in weight_mat]).T).T
    distance = math.sqrt(np.sum([(a - b) ** 2 for a, b in zip(weighted_colors, weighted_centroid)]))
    Results_frame.at[i,'Color_contrast_8bit']=distance
    #################################################################
    ###Now the 3bit color data
    retrial=trial.applymap(threshold_relabel)
    aggregate_counter=Counter(retrial.values.flatten())
    names_3bit+=[i]
    list_o_lists_3bit+=[prop_color_histogram3(aggregate_counter,startprop,endprop,stepprop)[0]]
    topC,topP=top_color(aggregate_counter,0)
    Results_frame.at[i,'Primary_color_3bit']=topC
    Results_frame.at[i,'Primary_color_3bit_prop']=topP
    topC,topP=top_color(aggregate_counter,1)
    Results_frame.at[i,'Secondary_color_3bit']=topC
    Results_frame.at[i,'Secondary_color_3bit_prop']=top_color(aggregate_counter,1)[1]
    topC,topP=top_color(aggregate_counter,2)
    Results_frame.at[i,'Tertiary_color_3bit']=topC
    Results_frame.at[i,'Tertiary_color_3bit_prop']=topP
    colors,props=prop_color_histogram3(aggregate_counter,startprop,endprop,stepprop)
    MDI=Moment_Distance_Index(props,colors,startprop,endprop)
    Results_frame.at[i,'MDI_3bit']=MDI
    richness=color_prop2(aggregate_counter,0.99)
    Results_frame.at[i,'Color_richness_3bit']=richness
    Results_frame.at[i,'Norm_perchannel_Color_richness_3bit']=richness**(1/num_channels)
    SW,evenness=SW_and_evenness(aggregate_counter)
    Results_frame.at[i,'SW_entropy_3bit']=SW
    Results_frame.at[i,'Perchannel_effective_num_colors_3bit']=(math.exp(SW))**(1/num_channels)
    Results_frame.at[i,'Norm_actual_to_effective_colors_3bit']=(richness**(1/num_channels))/((math.exp(SW))**(1/num_channels))
    Results_frame.at[i,'SW_evenness_3bit']=evenness
    ###now for color contrast calculations
    del aggregate_counter[None]
    color_lists=[listify(x) for x in list(aggregate_counter.keys())]
    weight_mat=list(aggregate_counter.values())
    weighted_colors=np.array([np.array(x).dot(y) for x, y in zip(color_lists, weight_mat)])
    color_centroid=np.sum(weighted_colors,axis=0)/np.sum(weight_mat)
    Results_frame.at[i,'Color_centroid_3bit']=list(color_centroid)
    weighted_centroid=np.array([[x] for x in color_centroid]).dot(np.array([[x] for x in weight_mat]).T).T
    distance = math.sqrt(np.sum([(a - b) ** 2 for a, b in zip(weighted_colors, weighted_centroid)]))
    Results_frame.at[i,'Color_contrast_3bit']=distance

###This for loop iterates over every matrix file in 'matrices' and executes the above 'forloop_frameOPS' function, which
##automatically adds values to the master frame
for i in matrices.name:
    count+=1
    x = loadmat(matrices_dir_path+'%s' % (i))
    ##Loading the scale
    scale=x['transspp'][0][1][0][0]
    Results_frame_bodypart.at[i,'1cm_scale_pixels']=scale
    #loading the array
    array=x['transspp'][0][0]
    pixel_frame=pd.DataFrame(array,dtype='int')
    pixel_frame=pd.DataFrame(pixel_frame,dtype='str')
    try: 
        forloop_frameOPS(pixel_frame,i,Results_frame_bodypart,names_8bit,list_o_lists_8bit,names_3bit,list_o_lists_3bit)
        print(i,'body_part_%s:done'%(body_part))
    except:
        print(i,'body_part_%s:ABSENT'%(body_part))
    
cols=np.round(np.arange(startprop,endprop,stepprop),2)

NC_frame_8bit=pd.concat([pd.DataFrame(names_8bit,columns=["names_8bit"]),pd.DataFrame(list_o_lists_8bit,columns=cols)],axis=1)
NC_frame_3bit=pd.concat([pd.DataFrame(names_3bit,columns=["names_3bit"]),pd.DataFrame(list_o_lists_3bit,columns=cols)],axis=1)

###This outputs the main results needed for future steps in the color palette ppeline
Results_frame_bodypart.to_csv('%s%s_All_in_one_NClist_colors_MDI_streamlined_bit_color_results_BodyPart%s_C%s_to_C%s.csv'%(output_dir,taxon_group,body_part,channel_start,channel_stop))
pickle.dump(Results_frame_bodypart, open("%s%s_All_in_one_NClist_colors_MDI_streamlined_bit_color_results_BodyPart%s_C%s_to_C%s.p"%(output_dir,taxon_group,body_part,channel_start,channel_stop), "wb"))

###the following two blocks outputs color frequency histograms (8bit or 3bit) with specimens as rows and percentages from 0-0.95 as columns,
#with the values being integers that describe ##how many of the most frequent colors cumulatively account for that percentage of the wing area. 
#e.g. If the top 3 colors account for 60,20 and 18% of the wing area then the histogram would contain 0 values until column 0.6, then 1 values until column 0.8 and so on.
##commented out for now
# NC_frame_8bit.to_csv('%s%s_All_in_one_NClist_colors_MDI_streamlined_bit_color_NC_frame_8bit_BodyPart%s_C%s_to_C%s.csv'%(output_dir,taxon_groupbody_part,channel_start,channel_stop))
# pickle.dump(NC_frame_8bit, open("%s%s_All_in_one_NClist_colors_MDI_streamlined_bit_color_NC_frame_8bit_BodyPart%s_C%s_to_C%s.p"%(output_dir,taxon_group,body_part,channel_start,channel_stop), "wb"))

# NC_frame_3bit.to_csv('%s%s_All_in_one_NClist_colors_MDI_streamlined_bit_color_NC_frame_3bit_BodyPart%s_C%s_to_C%s.csv'%(output_dir,taxon_group,body_part,channel_start,channel_stop))
# pickle.dump(NC_frame_3bit, open("%s%s_All_in_one_NClist_colors_MDI_streamlined_bit_color_NC_frame_3bit_BodyP,art%s_C%s_to_C%s.p"%(output_dir,taxon_group,body_part,channel_start,channel_stop), "wb"))

stop=time()
print("elapsed time:",stop-start)
print("")
