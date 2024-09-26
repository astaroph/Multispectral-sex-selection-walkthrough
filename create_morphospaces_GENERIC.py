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
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
from matplotlib.lines import Line2D
import matplotlib as mpl
def confidence_ellipse(x, y, ax, n_std=2.0, facecolor='none', **kwargs):
    """
    Create a plot of the covariance confidence ellipse of *x* and *y*.

    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    Returns
    -------
    matplotlib.patches.Ellipse

    Other parameters
    ----------------
    kwargs : `~matplotlib.patches.Patch` properties
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensionl dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0),
        width=ell_radius_x * 2,
        height=ell_radius_y * 2,lw=2,
        facecolor=facecolor,
        **kwargs)

    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)

#############################################################################################################################
#############################################################################################################################
start=time()
usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)
parser.add_option("-d", "--dataset_str", dest='dataset_str', type="string",
                 help="input dataset filename")
parser.add_option("-t", "--Taxon", dest='Taxon', type="string",default='Butterflies',
                 help="Taxon group corresponding to dataset, can be 'Butterflies' or 'Moths' currently")
parser.add_option("-c", "--cat_var", dest='cat_var', type="string",
                 help="Categorical variable to include in final PC datasets and to optionally filter by")
parser.add_option("-f", "--filter_str", dest='filter_str', type="string",
                 help="Restrict dataset to those matching this string in cat_var")
parser.add_option("-o", "--output_prefix", dest='output_prefix', type="string",
                 help="unique prefix to append to filenames")
parser.add_option("-b", "--traits_begin", dest='traits_begin', type="int",default=20,
                 help="where traits begin in the trait dataset")
parser.add_option("-i", "--ID_var", dest='ID_var', type="string",default='Phy_ID',
                 help="string name of the unique species ID column")
parser.add_option("-s", "--specimen_cutoff", dest='specimen_cutoff', type="int",default=1,
                 help="Minimum number of male and female specimens a species must have to be retained before PCA")
parser.add_option("-F", "--output_mf_filt", dest='output_mf_filt', action="store_true",default=True,
                 help="Boolean. whether to save the raw trait dataset that has been filtered by the specimen threshold. Defaults to True")
parser.add_option("-A", "--output_All_PCs", dest='output_All_PCs',action="store_true", default=True,
                 help="Boolean, whether to save the PC dataset that contains all PCs. Defaults to True")
parser.add_option("-N", "--output_topN_PCs", dest='output_topN_PCs',action="store_true", default=True,
                 help="Boolean, whether to save the PC dataset that contains only the PCs explaining the top -p proportion of variance. True by default ")
parser.add_option("-p", "--topN_PC_prop", dest='topN_PC_prop', type="float",default=.90,
                 help="The threshold proportion of variance explained by the top N PCs")
parser.add_option("-P", "--plot_PCA", dest='plot_PCA',action="store_true",default=True,
                 help="Boolean, whether to plot the PC plot. Defaults to true")
parser.add_option("-x", "--PC1", dest='PC1', type="int",default=1,
                 help="PC to set on the x axis")
parser.add_option("-y", "--PC2", dest='PC2', type="int",default=2,
                 help="PC to set on the y axis")
parser.add_option("-V", "--variable", dest='variable', type="string",default='Family',
                 help="Grouping variable for PC plot")
parser.add_option("-S", "--shapevar", dest='shapevar', type="string",
                 help="Point shape variable for PC plot")
parser.add_option("-O", "--open_closed_var", dest='open_closed_var', type="string",
                 help="Point open or closed variable for PC plot")
parser.add_option("-C", "--create_singleband_morphospaces", dest='create_singleband_morphospaces', action="store_true",default=True,
                 help="Whether to create separate morphospaces for the reflectance based single band metrics datasets")


(options, args) = parser.parse_args()

##Specifying the inputs
dataset_str=options.dataset_str
#Taxon: "Butterflies" or "Moths"
Taxon=options.Taxon
#specify the datasets and the unique categorical grouping variables of interest that you wish to include in the final PC dataframe (note, "Family" "Sex", "Side", "Wing" and the ID variable 
#are included by default)
#If no categorical_variable is needed, put cat_var=None
cat_var=options.cat_var
##The value of "cat_var" to which the dataset will be restricted, if desired
filter_str=options.filter_str
##Unique string descriptor to add to filenames
output_prefix=options.output_prefix
#where the traits begin in the original supplied trait dataset
traits_begin=options.traits_begin
##The name of the species ID variable
ID_var=options.ID_var

##How many males and females needed for a species to be retained
specimen_cutoff=options.specimen_cutoff
##Whether to save the male-female cutoff filtered datasets
output_mf_filt=options.output_mf_filt
##Whether to output the full PC dataframe without filtering for the PCs explaining the top N proportion variance
output_All_PCs=options.output_All_PCs
##Whether to output the PC dataframe  filtering for the PCs explaining the top N proportion variance
output_topN_PCs=options.output_topN_PCs
##The proportion of variance explained for filtering the topN PCs
topN_PC_prop=options.topN_PC_prop
# ##Whether to save the raw 3 males and 3 females filtered datasets
# output3m3f=False
##Whether to plot the PC plots with the given PCs
plot_PCA=options.plot_PCA
##Which PCs to plot in the overall PC data
PC1=options.PC1
PC2=options.PC2
##Grouping variable for PC plot
variable=options.variable
##Point shape variable for PC plot
shapevar=options.shapevar
##Open-closed point variable for PC plot
open_closed_var=options.open_closed_var
###Whether to create the single band morphospaces
create_singleband_morphospaces=options.create_singleband_morphospaces

###############################################################################################################################################
###Loading ALL of the MCZ imaged specimens, including those that are unsexed or don't have enough males and females
master_wings2=pd.read_csv(dataset_str)
master_wings2=master_wings2.set_index('%s'%(ID_var),drop=False)
master_wings2.insert(0,'UID',master_wings2.Specimen+"_"+master_wings2['%s'%(ID_var)].astype('str')+"_"+master_wings2.Sex+"_"+master_wings2.Side+"_"+master_wings2.Wing)
master_wings2.insert(0,'Sex_Side_Wing',master_wings2.Sex+"_"+master_wings2.Side+"_"+master_wings2.Wing)
master_wings2.insert(0,'Family_Side_Wing',master_wings2.Family+"_"+"_"+master_wings2.Side+"_"+master_wings2.Wing)
master_wings2=master_wings2[(master_wings2.Sex!='unknown') & (master_wings2.Sex!='Unknown')]
##filtering steps
# master_wings3=master_wings2[master_wings2.Diel_behavior!='All']
# master_wings3=master_wings3[master_wings3.Diel_behavior!='Crepuscular']
master_wings3=master_wings2
#Now making a 3m3f cutoff and 1m1f cutoff dataset from the overall image data
sexcounts=master_wings3.set_index('UID',drop=False).groupby(['%s'%(ID_var),'Family','Sex','Side','Wing'],axis=0).count().reset_index().iloc[:,:7]
sexcounts=sexcounts.rename(columns={'Family_Side_Wing':'Count'}).set_index('%s'%(ID_var),drop=False)
# sexcounts_M=sexcounts[sexcounts.Sex=='male']
# sexcounts_F=sexcounts[sexcounts.Sex=='female']
# m1f1=set(sexcounts_M['%s'%(ID_var)]).intersection(set(sexcounts_F['%s'%(ID_var)]))
sexcounts_thresh=sexcounts[sexcounts.Count>=specimen_cutoff]
sexcounts_M=sexcounts_thresh[sexcounts_thresh.Sex=='male']
sexcounts_F=sexcounts_thresh[sexcounts_thresh.Sex=='female']
mf_filt=list(set(sexcounts_M['%s'%(ID_var)]).intersection(set(sexcounts_F['%s'%(ID_var)])))
print(len(mf_filt))
##Outputing these filtered datasets
master_wings2_mf=master_wings3.loc[mf_filt]
# master_wings2_3m3f=master_wings3.loc[m3f3]


# time_struct=localtime()

# if output_mf_filt==True:
#     ####Outputing the raw datasets filtered to the number of males and females matching the minimum specimen_threshold without any specimens of unknown sex
#     pickle.dump(master_wings2_mf, open("Master_%s_trait_dataset_%sMF_%s_%s_%s.p" %(output_prefix,specimen_cutoff,time_struct[2],time_struct[1],time_struct[0]), "wb"))
#     master_wings2_mf.to_csv("Master_%s_trait_dataset_%sMF_%s_%s_%s.csv" %(output_prefix,specimen_cutoff,time_struct[2],time_struct[1],time_struct[0]),index=False)

if output_mf_filt==True:
    ####Outputing the raw datasets filtered to the number of males and females matching the minimum specimen_threshold without any specimens of unknown sex
    pickle.dump(master_wings2_mf, open("Master_%s_trait_dataset_%sMF.p" %(output_prefix,specimen_cutoff), "wb"))
    master_wings2_mf.to_csv("Master_%s_trait_dataset_%sMF.csv" %(output_prefix,specimen_cutoff),index=False)


######Now proceeding to make the PCA morphospace
entry_dataset=master_wings2_mf  

if 'UID' in entry_dataset.columns:
    subset=entry_dataset.set_index('UID',drop=False)
else:
    entry_dataset.insert(0,'UID',entry_dataset.Specimen+"_"+entry_dataset['%s'%(ID_var)].astype('str')+"_"+entry_dataset.Sex+"_"+entry_dataset.Side+"_"+entry_dataset.Wing)
    subset=entry_dataset.set_index('UID',drop=False)

if cat_var==None:
    catnames=['Specimen','Family','Genus','Species','%s'%(ID_var),'Sex','Side','Wing']
    catnames=[element for element in catnames if element in subset.columns]
else:
    catnames=['Specimen','Family','Genus','Species','%s'%(ID_var),'Sex','Side','Wing']
    catnames=[element for element in catnames if element in subset.columns]
    catnames+=[cat_var]

if filter_str!=None:
    subset=subset[subset[cat_var]==filter_str]

#use this to select only metrics which are less than 90% 0s and have a variance greater than 0
subset2=subset.iloc[:,traits_begin:].loc[:,subset.iloc[:,traits_begin:].quantile(.9)>0].loc[:,subset.var(numeric_only=True)>0]
#use this to not perform any quantile based subselection
# subset2=subset.iloc[:,9:]
full_dataset=subset
metrics_dataset=subset2
X = metrics_dataset
scaler = StandardScaler()
scaler.fit(X)
X=scaler.transform(X)    
pca = PCA()
x_pca=pca.fit_transform(X)
x_new = pd.DataFrame(x_pca)

PCA_scoreframe=x_new.set_index(full_dataset['UID'])
loadings = pd.DataFrame(pca.components_.T,index=subset2.columns)

loadings_subset=loadings.iloc[:,:pca.explained_variance_ratio_[np.cumsum(pca.explained_variance_ratio_)<topN_PC_prop].shape[0]]
loading_score=np.sum(abs(loadings_subset*pca.explained_variance_ratio_[np.cumsum(pca.explained_variance_ratio_)<topN_PC_prop]),axis=1)
loadings_subset=loadings_subset.loc[loading_score.sort_values(ascending=False).index]

##making a PCA_scores frame based on only the PCs that explain 90% of the variation, as above with the loadings frame
PCA_scores=full_dataset.loc[:,catnames].join(PCA_scoreframe.iloc[:,:pca.explained_variance_ratio_[np.cumsum(pca.explained_variance_ratio_)<topN_PC_prop].shape[0]])
PCA_scores_full=full_dataset.loc[:,catnames].join(PCA_scoreframe.iloc[:,:pca.explained_variance_ratio_[np.cumsum(pca.explained_variance_ratio_)<1].shape[0]])

# ###The above code automatically updates based on which taxon is selected, but categorical variables of interest need to be manually input here
elements_to_remove=['Specimen','Species','Genus','Genus_species']
catnames=[element for element in catnames if element not in elements_to_remove]
# catnames.remove('Specimen')
# catnames.remove('Species')
# catnames.remove('Genus')
PCA_scores_grouped=PCA_scores_full.groupby(catnames,axis=0).median(numeric_only=True).reset_index()
PCA_var_start=len(catnames)

PCA_scores_grouped.insert(0,'Phy_Side_Wing',PCA_scores_grouped['%s'%(ID_var)].astype('str')+"_"+PCA_scores_grouped.Side+"_"+PCA_scores_grouped.Wing)
PCA_F=PCA_scores_grouped[PCA_scores_grouped.Sex=='female']
PCA_F=PCA_F.set_index('Phy_Side_Wing')
PCA_M=PCA_scores_grouped[PCA_scores_grouped.Sex=='male']
PCA_M=PCA_M.set_index('Phy_Side_Wing')
a=set(PCA_M.index)
b=set(PCA_F.index)
indices=list(a.intersection(b))
PCA_M1=PCA_M.loc[indices]
PCA_F1=PCA_F.loc[indices]
PCA_scores1=PCA_scores_grouped.set_index('Phy_Side_Wing').loc[indices]

distframe=PCA_M1.iloc[:,:PCA_var_start]
distframe['MF_PC_distance']=np.linalg.norm(PCA_M1.iloc[:,PCA_var_start:]-PCA_F1.iloc[:,PCA_var_start:],axis=1)
distframe['Sex']='Male-Female distance'
distframe1=pd.concat([distframe,PCA_M1.iloc[:,PCA_var_start:]-PCA_F1.iloc[:,PCA_var_start:]],axis=1)

if output_All_PCs==True:
    loadings.to_csv('%s_raw_loadings_sensitivity_curated_%sMF.csv'  %(output_prefix,specimen_cutoff))
    PCA_scores_full.to_csv('%s_PCA_scoreframe_sensitivity_curated_%sMF_AllPCs.csv' %(output_prefix,specimen_cutoff))
    pickle.dump(PCA_scores_full,open("%s_PCA_scoreframe_sensitivity_curated_%sMF_AllPCs.p"  %(output_prefix,specimen_cutoff), "wb"))
    
if output_topN_PCs==True:
    loadings_subset.to_csv('%s_PCA_sorted_90percent_loadings_sensitivity_curated_%sMF.csv' %(output_prefix,specimen_cutoff))
    PCA_scores.to_csv('%s_PCA_scoreframe_sensitivity_curated_%sMF.csv' %(output_prefix,specimen_cutoff))
    pickle.dump(PCA_scores,open("%s_PCA_scoreframe_sensitivity_curated_%sMF.p"  %(output_prefix,specimen_cutoff), "wb"))

distframe.to_csv('%s_M_minus_F_perspecimenPCA_medians_scoreframe_euclidean_distance_vector_%sMF_allPCs.csv'%(output_prefix,specimen_cutoff))
distframe1.to_csv('%s_M_minus_F_perspecimenPCA_medians_scoreframe_euclidean_distance_vector_%sMF_allPCs_with_M_minus_F_PCs.csv'%(output_prefix,specimen_cutoff))
#######

if plot_PCA==True:
    cumsum=np.cumsum(pca.explained_variance_ratio_[np.cumsum(pca.explained_variance_ratio_)<0.95])
    plt.plot(range(1,len(cumsum)+1),cumsum)
    plt.xlabel('number of components')
    plt.ylabel('cumulative explained variance');
    plt.savefig('%s_M_minus_F_perspecimenPCA_medians_%sMF_var_explained.png'%(output_prefix,specimen_cutoff),format='png',dpi=200)
    plt.savefig('%s_M_minus_F_perspecimenPCA_medians_%sMF_var_explained.pdf'%(output_prefix,specimen_cutoff),format='pdf',dpi=200)

    # PCA=pd.read_csv('Butterflies_Mar2022_nofluorpol_BS_sensitivity_redo_PCA_scoreframe_sensitivity_curated_AllPCs.csv')
    PCA_frame=PCA_scores
    PCA_frame.insert(0,'Sex_Side_Wing',PCA_frame.Sex+'_'+PCA_frame.Side+'_'+PCA_frame.Wing)

    # grouped=PCA_frame.groupby(['Family','%s'%(ID_var),'Sex_Side_Wing','Sex','Side','Wing'],axis=0).median(numeric_only=True).reset_index()
    grouped=PCA_frame.groupby(catnames,axis=0).median(numeric_only=True).reset_index()

    grouped.insert(0,'UID2',grouped['%s'%(ID_var)].astype(str) + '_' + grouped.Sex + '_' + grouped.Side + '_' + grouped.Wing)
    grouped=grouped.set_index('UID2',drop=False)
    # grouped_DFW=grouped[(grouped.Side=='dorsal') & (grouped.Wing=='FW')]
    # grouped_VFW=grouped[(grouped.Side=='ventral') & (grouped.Wing=='FW')]
    # grouped_DHW=grouped[(grouped.Side=='dorsal') & (grouped.Wing=='HW')]
    # grouped_VHW=grouped[(grouped.Side=='ventral') & (grouped.Wing=='HW')]
    grouped_D=grouped[grouped.Side=='dorsal']
    grouped_V=grouped[grouped.Side=='ventral']

    ####Now plotting the 
    ##I had to define a custom scatter function to allow a list of markers to be given to matplotlib
    #modified it from https://stackoverflow.com/questions/51810492/how-can-i-add-a-list-of-marker-styles-in-matplotlib
    def mscatter(x,y, ax=None, m=None, **kw):
        import matplotlib.markers as mmarkers
        ax = ax or plt.gca()
        sc = ax.scatter(x,y,**kw)
        if (m is not None) and (len(m)==len(x)):
            paths = []
            for marker in m:
                if isinstance(marker, mmarkers.MarkerStyle):
                    marker_obj = marker
                else:
                    marker_obj = mmarkers.MarkerStyle(marker)
                path = marker_obj.get_path().transformed(
                            marker_obj.get_transform())
                paths.append(path)
            sc.set_paths(paths)
        return sc

    def myplot(titlestring,entry_dataset=grouped_D,variable='Family',shapevar='Sex_Side_Wing',open_closed_var='Sex',PC1=1,PC2=2):
        #for the grouped per species data:
        subset=entry_dataset
        subset2=subset.iloc[:,7:]
        full_dataset=subset
        full_dataset.insert(0,'Sex_Side_Wing',PCA_frame.Sex+'_'+PCA_frame.Side+'_'+PCA_frame.Wing)
        metrics_dataset=subset2

        #Making a callable dataframe of color lists.
        Wing_colors1=['Brown','Blue']
        Sex_colors1=['brown','Green','Red']
        Side_colors1=['Orange','Purple']
        # Family_colors1=['blue','magenta','green','red','orange','purple','black']
        if Taxon=='Butterflies':
            Family_colors1=['black','magenta','blue','green','purple','orange','red']
        else:
            kinds=len(grouped.loc[:,variable].unique())
            hsv = plt.get_cmap('tab20')
            colors = hsv(np.linspace(0, 1.0, kinds))
            Family_colors1=[mpl.colors.to_hex(i) for i in colors ]
            diel_colors=[mpl.colors.to_hex(i) for i in colors]
            # Family_colors1=['#201923','#ffffff','#fcff5d','#7dfc00','#0ec434','#228c68','#8ad8e8','#235b54','#29bdab','#3998f5','#37294f','#277da7','#3750db','#f22020','#991919','#ffcba5','#e68f66','#c56133','#96341c','#632819','#ffc413','#f47a22','#2f2aa0','#b732cc','#772b9d','#f07cab','#d30b94','#edeff3','#c3a5b4','#946aa2','#5d4c86']
        # Family_colors1=['green','blue','red','yellow','magenta','black','purple']
        mothbutt_colors1=['black','orange']
        # diel_colors=['black','orange']

        colorframe=pd.DataFrame()
        colorframe['Variable']=['Wing','Sex','Side','Family','Nocturnality','Diurnality_habit']
        colorframe['Color_List']=[Wing_colors1,Sex_colors1,Side_colors1,Family_colors1,mothbutt_colors1,diel_colors]
        colorframe=colorframe.set_index('Variable',drop=False)

        groups=full_dataset['%s' %(variable)].unique()
        print(groups)
        colors=colorframe.Color_List.loc['%s' %(variable)]
        print(colors)

        PCA_scoreframe=metrics_dataset.set_index(full_dataset['UID2'])
        x_frame=metrics_dataset.set_index(full_dataset['%s' %(variable)])

        #the following code makes a dictionary of each shape var and marker type, and incorporates those into a callable shape dictionary
        ##Which will replace the actual dataset of the 'variable' indexed full dataset with the shape markers so they can be called in
        #the pca plotting for loop below
        Side_shapes={"marker":["o","s"],"var":full_dataset.Side.unique()}
        Wing_shapes={"marker":["o","s"],"var":full_dataset.Wing.unique()}
        Sex_shapes={"marker":["o","s",'d'],"var":full_dataset.Sex.unique()}
        # Sex_side_shapes={"marker":["d","X",'D','P'],"var":full_dataset.Sex_Side.unique()}
        # Sex_side_wing_shapes={"marker":["2","4",'1','3','^','>','v','<'],"var":full_dataset.Sex_Side_Wing.unique()}
        Sex_side_wing_shapes={"marker":["^","o",'','','.','+','x','d'],"var":full_dataset.Sex_Side_Wing.unique()}

        shape_dict={'Side':Side_shapes,'Wing':Wing_shapes,'Sex':Sex_shapes,'Sex_Side_Wing':Sex_side_wing_shapes}
        # shape_dict={'Side':Side_shapes,'Sex':Sex_shapes}
        # shape_dict={'Side':Side_shapes}

        Side_open_close={"facecolors":["none",np.nan],"var":full_dataset.Side.unique()}
        Wing_open_close={"facecolors":["none",np.nan],"var":full_dataset.Wing.unique()}
        Sex_open_close={"facecolors":["none",np.nan],"var":full_dataset.Sex.unique()}
        open_close_dict={'Side':Side_open_close,'Wing':Wing_open_close,'Sex':Sex_open_close}

        trial_dataset=full_dataset.set_index('%s' %(variable)).loc[:,['%s'%(shapevar),'%s'%(open_closed_var)]]
        for i in range(len(full_dataset['%s'%(shapevar)].unique())):
            trial_dataset=trial_dataset.replace(shape_dict['%s' %(shapevar)]['var'][i],shape_dict['%s' %(shapevar)]['marker'][i])
        print(shape_dict['%s' %(shapevar)])       
        for i in range(len(full_dataset['%s'%(open_closed_var)].unique())):
            trial_dataset=trial_dataset.replace(open_close_dict['%s' %(open_closed_var)]['var'][i],open_close_dict['%s' %(open_closed_var)]['facecolors'][i])
        for i in range(len(full_dataset['%s'%(variable)].unique())):
            trial_dataset.loc[groups[i]]=trial_dataset.loc[groups[i]].replace(np.nan,colors[i])
        print(open_close_dict['%s' %(open_closed_var)])

        ####Adding male female connectors to my centroids
        x_new1= metrics_dataset.set_index(full_dataset.UID2)
        x_new1.insert(0,'Sex',full_dataset.Sex)
        x_new1.insert(0,'Side',full_dataset.Side)
        x_new1.insert(0,'Wing',full_dataset.Wing)
        x_new1.insert(0,variable,full_dataset[variable])
        print(variable)
        groupingvars=['Sex','Side','Wing']
        groupingvars+=[variable]            
        x_new2=x_new1.groupby(groupingvars,axis=0).median(numeric_only=True).reset_index()

        males=x_new2[x_new2.Sex=='male']
        females=x_new2[x_new2.Sex=='female']

        score=x_frame
        xs = score.iloc[:,PC1-1]
        ys = score.iloc[:,PC2-1]
        scalex = 1.0/(xs.max() - xs.min())
        scaley = 1.0/(ys.max() - ys.min())
    #     scalex = 1
    #     scaley=1
        Fig=plt.figure(figsize=(40,35),dpi=80)
        ax1=Fig.add_subplot(111)
        for i,j in zip(groups,colors):
            print(i,j)
            xs = score.loc[i].iloc[:,PC1-1]
            ys = score.loc[i].iloc[:,PC2-1]
    #         scalex = 1.0/(xs.max() - xs.min())
    #         scaley = 1.0/(ys.max() - ys.min())
            mscatter(xs * scalex,ys * scaley,ax=ax1, c = j,s=600,alpha=0.1,label=i,m=trial_dataset.loc[i,'%s'%(shapevar)], facecolors=trial_dataset.loc[i,'%s'%(open_closed_var)])
            ##Adding centroid values for each shapevar split by open closed var
            for x in trial_dataset.loc[i,'%s' %(shapevar)].unique():
                pos=trial_dataset.loc[i,'%s' %(shapevar)]==x
                # print(pos)
                facecolored=trial_dataset.loc[i][pos]['%s'%(open_closed_var)]
                xs_new=xs[pos]
                ys_new=ys[pos]
            #             scalex_new = 1.0/(xs_new.max() - xs_new.min())
            #             scaley_new = 1.0/(ys_new.max() - ys_new.min())
                cent_x=np.median(xs_new * scalex)
                cent_y=np.median(ys_new * scaley)
                ax1.scatter(cent_x,cent_y,c=j,s=2500,marker=x,facecolors=facecolored)
            confidence_ellipse(xs * scalex,ys * scaley,ax1,2,edgecolor=j,alpha=0.7)
        X_coords= np.array([males[PC1-1], 
                               females[PC1-1]])*scalex
        Y_coords=np.array([males[PC2-1], 
                               females[PC2-1]])*scaley
        ax1.plot(X_coords, 
             Y_coords, 
             color='black',ls='solid',lw=5,alpha=0.8)
        ax1.set_xlabel("PC{}".format(PC1),fontsize=64)
        ax1.set_ylabel("PC{}".format(PC2),fontsize=64)
        ax1.tick_params(axis='x',labelsize=54)
        ax1.tick_params(axis='y',labelsize=54)
        ax1.grid()
        ax1.grid()
        handles, labels = ax1.get_legend_handles_labels()
        leg1=ax1.legend(handles, labels,loc='upper right',fontsize=36,markerscale=2.0)
        ax1.set_title(titlestring,fontsize=48)
        #to add a dorsal/ventral legend
        legend_elements = []
        for i in range(len(trial_dataset['%s'%(shapevar)].unique())):
            legend_elements+=[Line2D([0], [0], marker=shape_dict['%s'%(shapevar)]['marker'][i], color='black', label=shape_dict['%s'%(shapevar)]['var'][i], markersize=30)]
        leg2=ax1.legend(handles=legend_elements, loc='lower right',fontsize=36)
        ax1.add_artist(leg1)
        return Fig

    ##for color trait plot only
    # ax=myplot(x_frame,1,2,titlestring,subset2.columns)
    ax_D=myplot(titlestring='Per species PCA of spectral and pattern traits by %s'%(variable),entry_dataset=grouped_D,variable=variable,shapevar=shapevar,open_closed_var=open_closed_var,PC1=PC1,PC2=PC2)
    ax_V=myplot(titlestring='Per species PCA of spectral and pattern traits by %s'%(variable),entry_dataset=grouped_V,variable=variable,shapevar=shapevar,open_closed_var=open_closed_var,PC1=PC1,PC2=PC2)
    # plt.show()
    ax_D.savefig('%s_M_minus_F_perspecimenPCA_medians_%sMF_PCA_plot_by_%s_%s_%s_dorsal.png'%(output_prefix,specimen_cutoff,variable,shapevar,open_closed_var),format='png',dpi=200)
    ax_D.savefig('%s_M_minus_F_perspecimenPCA_medians_%sMF_PCA_plot_by_%s_%s_%s_dorsal.pdf'%(output_prefix,specimen_cutoff,variable,shapevar,open_closed_var),format='pdf',dpi=200)
    ax_V.savefig('%s_M_minus_F_perspecimenPCA_medians_%sMF_PCA_plot_by_%s_%s_%s_ventral.png'%(output_prefix,specimen_cutoff,variable,shapevar,open_closed_var),format='png',dpi=200)
    ax_V.savefig('%s_M_minus_F_perspecimenPCA_medians_%sMF_PCA_plot_by_%s_%s_%s_ventral.pdf'%(output_prefix,specimen_cutoff,variable,shapevar,open_closed_var),format='pdf',dpi=200)

if create_singleband_morphospaces==True:
    #######Now doing creating the single band morphospaces
    bands=['UV','B','G','R','740','940']
    # Taxon='Moths'
    # cat_var=['Diel_behavior']
    # filter_str='Nocturnal'
    # output_prefix='Nocturnal_Moths_1m1f_sensitivitycut_Final'

    # Taxon='Butterflies'
    # cat_var=None
    # filter_str=None
    # output_prefix='Butterflies_1m1f_sensitivitycut_Final'

    for i in bands:
        band_cols=[x for x in master_wings2_mf.columns if '_%s_refl' %(i) in x or '_%s_pct' %(i) in x or '%s_png' %(i) == x]
        entry_dataset=master_wings2_mf
        if Taxon=='Moths':
            entry_dataset=master_wings2_mf
            traits_begin=23
            ID_var='Moth_ID'
        elif Taxon=='Butterflies':
            entry_dataset=master_wings2_mf
            # entry_dataset=entry_dataset.rename(columns={"Specimen_Barcode": "Barcode"})
            traits_begin=20
            ID_var='Phy_ID'
        if 'UID' in entry_dataset.columns:
            subset=entry_dataset.set_index('UID',drop=False)
        else:
            entry_dataset.insert(0,'UID',entry_dataset.Specimen+"_"+entry_dataset['%s'%(ID_var)].astype('str')+"_"+entry_dataset.Sex+"_"+entry_dataset.Side+"_"+entry_dataset.Wing)
            subset=entry_dataset.set_index('UID',drop=False)
        if cat_var==None:
            catnames=['Specimen','Family','Genus','Species','%s'%(ID_var),'Sex','Side','Wing']
            catnames=[element for element in catnames if element in subset.columns]
        else:
            catnames=['Specimen','Family','Genus','Species','%s'%(ID_var),'Sex','Side','Wing']
            catnames=[element for element in catnames if element in subset.columns]
            catnames+=[cat_var]

        if filter_str!=None:
            subset=subset[subset[cat_var]==filter_str]
        metric_cols=subset.loc[:,band_cols]
        #use this to select only metrics which are less than 90% 0s and have a variance greater than 0
        subset2=metric_cols.loc[:,metric_cols.quantile(.9)>0].loc[:,metric_cols.var(numeric_only=True)>0]
        full_dataset=subset
        metrics_dataset=subset2
        X = metrics_dataset
        scaler = StandardScaler()
        scaler.fit(X)
        X=scaler.transform(X)    
        pca = PCA()
        x_pca=pca.fit_transform(X)
        x_new = pd.DataFrame(x_pca)

        PCA_scoreframe=x_new.set_index(full_dataset['UID'])
        loadings = pd.DataFrame(pca.components_.T,index=subset2.columns)

        loadings_subset=loadings.iloc[:,:pca.explained_variance_ratio_[np.cumsum(pca.explained_variance_ratio_)<topN_PC_prop].shape[0]]
        loading_score=np.sum(abs(loadings_subset*pca.explained_variance_ratio_[np.cumsum(pca.explained_variance_ratio_)<topN_PC_prop]),axis=1)
        loadings_subset=loadings_subset.loc[loading_score.sort_values(ascending=False).index]

        ##making a PCA_scores frame based on only the PCs that explain 90% of the variation, as above with the loadings frame
        PCA_scores=full_dataset.loc[:,catnames].join(PCA_scoreframe.iloc[:,:pca.explained_variance_ratio_[np.cumsum(pca.explained_variance_ratio_)<topN_PC_prop].shape[0]])
        PCA_scores_full=full_dataset.loc[:,catnames].join(PCA_scoreframe.iloc[:,:pca.explained_variance_ratio_[np.cumsum(pca.explained_variance_ratio_)<1].shape[0]])

        ###The above code automatically updates based on which taxon is selected, but categorical variables of interest need to be manually input here
        elements_to_remove=['Specimen','Species','Genus','Genus_species']
        catnames=[element for element in catnames if element not in elements_to_remove]
        # catnames.remove('Specimen')
        # catnames.remove('Species')
        # catnames.remove('Genus')
        PCA_scores_grouped=PCA_scores_full.groupby(catnames,axis=0).median(numeric_only=True).reset_index()
        PCA_var_start=len(catnames)
        PCA_scores_grouped.insert(0,'Phy_Side_Wing',PCA_scores_grouped['%s'%(ID_var)].astype('str')+"_"+PCA_scores_grouped.Side+"_"+PCA_scores_grouped.Wing)
        PCA_F=PCA_scores_grouped[PCA_scores_grouped.Sex=='female']
        PCA_F=PCA_F.set_index('Phy_Side_Wing')
        PCA_M=PCA_scores_grouped[PCA_scores_grouped.Sex=='male']
        PCA_M=PCA_M.set_index('Phy_Side_Wing')
        a=set(PCA_M.index)
        b=set(PCA_F.index)
        indices=list(a.intersection(b))
        PCA_M1=PCA_M.loc[indices]
        PCA_F1=PCA_F.loc[indices]
        PCA_scores1=PCA_scores_grouped.set_index('Phy_Side_Wing').loc[indices]
        distframe=PCA_M1.iloc[:,:PCA_var_start]
        distframe['MF_PC_distance']=np.linalg.norm(PCA_M1.iloc[:,PCA_var_start:]-PCA_F1.iloc[:,PCA_var_start:],axis=1)
        distframe['Sex']='Male-Female distance'
        distframe1=pd.concat([distframe,PCA_M1.iloc[:,PCA_var_start:]-PCA_F1.iloc[:,PCA_var_start:]],axis=1)
        # print(len(PCA_scores_full.index))
        PCA_scores_full.to_csv('%s_PCA_scoreframe_sensitivity_curated_%sMF_AllPCs_%s.csv'%(output_prefix,specimen_cutoff,i))
        distframe.to_csv('%s_M_minus_F_perspecimenPCA_medians_scoreframe_euclidean_distance_vector_%sMF_AllPCs_%s.csv'%(output_prefix,specimen_cutoff,i))
       

        #######



stop=time()
print("elapsed time:",stop-start)
print("")
