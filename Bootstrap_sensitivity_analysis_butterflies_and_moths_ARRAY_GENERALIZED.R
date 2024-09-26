library(geiger)
library(plyr)
library(ape)
library(phytools)

make_csvs=FALSE
final_csvs=FALSE

val_func<-function(metric_dataset, metric_str, treefile){
  dtt.subtree<-dtt.to.df(calculate.feature.dtt(treefile, metric_dataset, metric_str, 1), treefile)
  vector<-dtt.subtree$dtt[2:length(treefile$tip.label)]
  dtt_val<-mean(vector)
  mean_val<-mean(metric_dataset[,metric_str])
  return(list(mean_val,dtt_val))
}


##MDI throught time
calculate.feature.dtt<-function(tree, df.phy.x.fore.p, colname, n.sim){
  tree.dat.x.fore.PC1<-unlist(df.phy.x.fore.p[,colname])
  names(tree.dat.x.fore.PC1)<-df.phy.x.fore.p$tipID
  tree.dat.x.fore.PC1<-tree.dat.x.fore.PC1[!is.na(tree.dat.x.fore.PC1)]
  #Drop unused tips
  tree.trimmed<-tree
  keep.tips<-names(tree.dat.x.fore.PC1)
  tree.trimmed<-drop.tip(tree.trimmed, tree.trimmed$tip.label[!tree.trimmed$tip.label %in% keep.tips])
  #remove unnecessary sub lists
  tree.trimmed$node.labels<-NULL
  tree.trimmed$node.label<-NULL
  
  shp.dtt.x.fore<-dtt(tree.trimmed, tree.dat.x.fore.PC1,index=c("avg.sq", "avg.manhattan", "num.states"),nsim=n.sim,plot=FALSE)
  return(shp.dtt.x.fore)
}

#a function to turn dtt result into data frame
dtt.to.df<-function(shp.dtt.all.fore, tree){
  dtt.df<-data.frame("time"=(1-shp.dtt.all.fore$time)*tree$root.time, "dtt"=shp.dtt.all.fore$dtt)
  return(dtt.df)
}

result.dir<-paste(getwd(),'/results',sep='')
indir=getwd()

load(file.path(indir,'key_examplar_variables_Dec2020.RData'))
######################set parameters here----
start_time<-Sys.time()
##First read in the arguments listed at the command line
args=(commandArgs(TRUE))
##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
  print("No arguments supplied.")
  ##supply default values
  boot = 1
  Taxon_group<-'Butterflies'
}else{
  for(i in 1:length(args)){
    print(args)
    eval(parse(text=args[[i]]))
  }
}
print(paste('bootstrapping ',boot))
set.seed(boot)
##Loading in the  metrics ----

if(Taxon_group=='Moths'){
  Metric_input<-read.csv('Master_Moths_trait_dataset_9_2_2023.csv')
  Metric_input<-Metric_input[Metric_input$Diel_behavior=='Nocturnal',]
  prefix_str='Moth_nocturnalONLY_Final'
  IDvar<-'Moth_ID'
  Metric_col_start=21
  Metric_col_end=263
}
if(Taxon_group=='Butterflies'){
  Metric_input<-read.csv('Master_Butterflies_trait_dataset_8_2_2023.csv')
  prefix_str='Butterflies_Final'
  IDvar<-'Phy_ID'
  Metric_col_start=18
  Metric_col_end=260
}

if(Taxon_group=='Geometrids'){
  Metric_input<-read.csv('Master_all_Geometrids_compiled_trait_dataset_7_1_2024_onlymatched_tips.csv')
  prefix_str='Geometrids_Final'
  IDvar<-'Geom_ID'
  Metric_col_start=18
  Metric_col_end=260
}

Metric_colnames<-colnames(Metric_input[,Metric_col_start:Metric_col_end])


###Add unnecessary variables to drop from the original dataframe here
if(Taxon_group=='Moths'){
  # drop=c('UID','Family_Side_Wing','Sex_Side_Wing','Diel_behavior_reference',
  #        'Genus','Species','X1cm_scale_pixels_UV_940','Moth_ID.1')
}


if(Taxon_group=='Butterflies'){
  # drop=c('UID','Family_Side_Wing','Sex_Side_Wing','Family_Side_Wing','Sex_Side_Wing',
  #        'Genus_valid','Identity_of_this_specimen','Phy_oldID','Swap_or_not','Loc_ID',
  #        'X1cm_scale_pixels_UV_940')
}

# if(Taxon_group=='Geometrids'){
#   drop=c('UID','Family_Side_Wing','Sex_Side_Wing','Diel_behavior_reference',
#          'Genus','Species','X1cm_scale_pixels_UV_940','Moth_ID.1')
# }


meanframe_MDF<-data.frame('Dataset'='MDF','metrics'=Metric_colnames)
meanframe_MVF<-data.frame('Dataset'='MVF','metrics'=Metric_colnames)
meanframe_MDH<-data.frame('Dataset'='MDH','metrics'=Metric_colnames)
meanframe_MVH<-data.frame('Dataset'='MVH','metrics'=Metric_colnames)
meanframe_FDF<-data.frame('Dataset'='FDF','metrics'=Metric_colnames)
meanframe_FVF<-data.frame('Dataset'='FVF','metrics'=Metric_colnames)
meanframe_FDH<-data.frame('Dataset'='FDH','metrics'=Metric_colnames)
meanframe_FVH<-data.frame('Dataset'='FVH','metrics'=Metric_colnames)

dttframe_MDF<-data.frame('Dataset'='MDF','metrics'=Metric_colnames)
dttframe_MVF<-data.frame('Dataset'='MVF','metrics'=Metric_colnames)
dttframe_MDH<-data.frame('Dataset'='MDH','metrics'=Metric_colnames)
dttframe_MVH<-data.frame('Dataset'='MVH','metrics'=Metric_colnames)
dttframe_FDF<-data.frame('Dataset'='FDF','metrics'=Metric_colnames)
dttframe_FVF<-data.frame('Dataset'='FVF','metrics'=Metric_colnames)
dttframe_FDH<-data.frame('Dataset'='FDH','metrics'=Metric_colnames)
dttframe_FVH<-data.frame('Dataset'='FVH','metrics'=Metric_colnames)


Metric_all<-Metric_input
names(Metric_all)[names(Metric_all) == paste(IDvar)] <- 'ID_var'

#####function continues below
#####making a boostraped dataset by separately sampling the specimens in males and females WITH REPLACEMENT, then recombining the new dataset
##before proceeding

all_metrics_MDF=subset(Metric_all,Sex=='male'& Side=='dorsal' & Wing=='FW')
all_metrics_MVF=subset(Metric_all,Sex=='male'& Side=='ventral' & Wing=='FW')
all_metrics_MDH=subset(Metric_all,Sex=='male'& Side=='dorsal' & Wing=='HW')
all_metrics_MVH=subset(Metric_all,Sex=='male'& Side=='ventral' & Wing=='HW')
all_metrics_FDF=subset(Metric_all,Sex=='female'& Side=='dorsal' & Wing=='FW')
all_metrics_FVF=subset(Metric_all,Sex=='female'& Side=='ventral' & Wing=='FW')
all_metrics_FDH=subset(Metric_all,Sex=='female'& Side=='dorsal' & Wing=='HW')
all_metrics_FVH=subset(Metric_all,Sex=='female'& Side=='ventral' & Wing=='HW')

dorsal_IDs=intersect(intersect(all_metrics_MDF$ID_var,all_metrics_FDF$ID_var),
                     intersect(all_metrics_MDH$ID_var,all_metrics_FDH$ID_var))
ventral_IDs=intersect(intersect(all_metrics_MVF$ID_var,all_metrics_FVF$ID_var),
                      intersect(all_metrics_MVH$ID_var,all_metrics_FVH$ID_var))
final_ids=intersect(dorsal_IDs,ventral_IDs)
all_metrics_MDF=subset(all_metrics_MDF,ID_var %in% final_ids)
all_metrics_MVF=subset(all_metrics_MVF,ID_var %in% final_ids)
all_metrics_FDF=subset(all_metrics_FDF,ID_var %in% final_ids)
all_metrics_FVF=subset(all_metrics_FVF,ID_var %in% final_ids)

all_metrics_MDH=subset(all_metrics_MDH,ID_var %in% final_ids)
all_metrics_MVH=subset(all_metrics_MVH,ID_var %in% final_ids)
all_metrics_FDH=subset(all_metrics_FDH,ID_var %in% final_ids)
all_metrics_FVH=subset(all_metrics_FVH,ID_var %in% final_ids)

#Use the following block of code to bootstrap such that each family has the same number of rows as the original data, but not necessarily the same number of species
# all_metrics_MDF <- ddply(all_metrics_MDF,.(Family),function(x) x[sample(nrow(x),nrow(x),replace=TRUE),])
# all_metrics_MVF <- ddply(all_metrics_MVF,.(Family),function(x) x[sample(nrow(x),nrow(x),replace=TRUE),])
# all_metrics_MDH <- ddply(all_metrics_MDH,.(Family),function(x) x[sample(nrow(x),nrow(x),replace=TRUE),])
# all_metrics_MVH <- ddply(all_metrics_MVH,.(Family),function(x) x[sample(nrow(x),nrow(x),replace=TRUE),])
# Metric_all_M<-rbind(all_metrics_MDF,all_metrics_MVF,all_metrics_MDH,all_metrics_MVH)

##Use the following block of code to bootstrap randomly across families
rownames(all_metrics_MDF)<-all_metrics_MDF$Specimen
rownames(all_metrics_MVF)<-all_metrics_MVF$Specimen
rownames(all_metrics_MDH)<-all_metrics_MDH$Specimen
rownames(all_metrics_MVH)<-all_metrics_MVH$Specimen
Specimens_M<-intersect(intersect(intersect(all_metrics_MDF$Specimen,all_metrics_MVF$Specimen),all_metrics_MDH$Specimen),all_metrics_MVH$Specimen)
Specimens_M1<-sample(Specimens_M,length(Specimens_M),replace=TRUE)
Metric_all_M<-rbind(all_metrics_MDF[Specimens_M1,],all_metrics_MVF[Specimens_M1,],all_metrics_MDH[Specimens_M1,],all_metrics_MVH[Specimens_M1,])



#Use the following block of code to bootstrap such that each family has the same number of rows as the original data, but not necessarily the same number of species
# all_metrics_FDF <- ddply(all_metrics_FDF,.(Family),function(x) x[sample(nrow(x),nrow(x),replace=TRUE),])
# all_metrics_FVF <- ddply(all_metrics_FVF,.(Family),function(x) x[sample(nrow(x),nrow(x),replace=TRUE),])
# all_metrics_FDH <- ddply(all_metrics_FDH,.(Family),function(x) x[sample(nrow(x),nrow(x),replace=TRUE),])
# all_metrics_FVH <- ddply(all_metrics_FVH,.(Family),function(x) x[sample(nrow(x),nrow(x),replace=TRUE),])
# Metric_all_F<-rbind(all_metrics_FDF,all_metrics_FVF,all_metrics_FDH,all_metrics_FVH)

##Use the following block of code to bootstrap randomly across families
rownames(all_metrics_FDF)<-all_metrics_FDF$Specimen
rownames(all_metrics_FVF)<-all_metrics_FVF$Specimen
rownames(all_metrics_FDH)<-all_metrics_FDH$Specimen
rownames(all_metrics_FVH)<-all_metrics_FVH$Specimen
Specimens_F<-intersect(intersect(intersect(all_metrics_FDF$Specimen,all_metrics_FVF$Specimen),all_metrics_FDH$Specimen),all_metrics_FVH$Specimen)
Specimens_F1<-sample(Specimens_F,length(Specimens_F),replace=TRUE)
Metric_all_F<-rbind(all_metrics_FDF[Specimens_F1,],all_metrics_FVF[Specimens_F1,],all_metrics_FDH[Specimens_F1,],all_metrics_FVH[Specimens_F1,])

Metric_all<-rbind(Metric_all_M,Metric_all_F)
rownames(Metric_all)<-NULL

names(Metric_all)[names(Metric_all) == paste(IDvar)] <- 'ID_var'

#####Drop column names from the dataset, specified in the 'drop' variable----
# Metric_all=Metric_all[,!names(Metric_all)%in%drop]
Phy_Side_Wing<-paste(Metric_all$ID_var,Metric_all$Side,Metric_all$Wing,sep='_')
Metric_all<-cbind(Phy_Side_Wing,Metric_all)


all_metrics<-Metric_all
species<-as.character(all_metrics$ID_var)
tipID<-as.character(all_metrics$ID_var)

all_metrics1<-cbind(tipID,species,all_metrics)

all_metrics_MDF=subset(all_metrics1,Sex=='male'& Side=='dorsal' & Wing=='FW')
all_metrics_MVF=subset(all_metrics1,Sex=='male'& Side=='ventral' & Wing=='FW')
all_metrics_FDF=subset(all_metrics1,Sex=='female'& Side=='dorsal' & Wing=='FW')
all_metrics_FVF=subset(all_metrics1,Sex=='female'& Side=='ventral' & Wing=='FW')

all_metrics_MDH=subset(all_metrics1,Sex=='male'& Side=='dorsal' & Wing=='HW')
all_metrics_MVH=subset(all_metrics1,Sex=='male'& Side=='ventral' & Wing=='HW')
all_metrics_FDH=subset(all_metrics1,Sex=='female'& Side=='dorsal' & Wing=='HW')
all_metrics_FVH=subset(all_metrics1,Sex=='female'& Side=='ventral' & Wing=='HW')


dorsal_IDs=intersect(intersect(all_metrics_MDF$ID_var,all_metrics_FDF$ID_var),
                     intersect(all_metrics_MDH$ID_var,all_metrics_FDH$ID_var))
ventral_IDs=intersect(intersect(all_metrics_MVF$ID_var,all_metrics_FVF$ID_var),
                      intersect(all_metrics_MVH$ID_var,all_metrics_FVH$ID_var))
final_ids=intersect(dorsal_IDs,ventral_IDs)
all_metrics_MDF1=subset(all_metrics_MDF,ID_var %in% final_ids)
all_metrics_MVF1=subset(all_metrics_MVF,ID_var %in% final_ids)
all_metrics_FDF1=subset(all_metrics_FDF,ID_var %in% final_ids)
all_metrics_FVF1=subset(all_metrics_FVF,ID_var %in% final_ids)

all_metrics_MDH1=subset(all_metrics_MDH,ID_var %in% final_ids)
all_metrics_MVH1=subset(all_metrics_MVH,ID_var %in% final_ids)
all_metrics_FDH1=subset(all_metrics_FDH,ID_var %in% final_ids)
all_metrics_FVH1=subset(all_metrics_FVH,ID_var %in% final_ids)

all_metrics_MDF2=ddply(all_metrics_MDF1,.(species,tipID,Family,ID_var,Sex,Side,Wing),colwise(median))
all_metrics_MVF2=ddply(all_metrics_MVF1,.(species,tipID,Family,ID_var,Sex,Side,Wing),colwise(median))
all_metrics_FDF2=ddply(all_metrics_FDF1,.(species,tipID,Family,ID_var,Sex,Side,Wing),colwise(median))
all_metrics_FVF2=ddply(all_metrics_FVF1,.(species,tipID,Family,ID_var,Sex,Side,Wing),colwise(median))

all_metrics_MDH2=ddply(all_metrics_MDH1,.(species,tipID,Family,ID_var,Sex,Side,Wing),colwise(median))
all_metrics_MVH2=ddply(all_metrics_MVH1,.(species,tipID,Family,ID_var,Sex,Side,Wing),colwise(median))
all_metrics_FDH2=ddply(all_metrics_FDH1,.(species,tipID,Family,ID_var,Sex,Side,Wing),colwise(median))
all_metrics_FVH2=ddply(all_metrics_FVH1,.(species,tipID,Family,ID_var,Sex,Side,Wing),colwise(median))


if(Taxon_group=='Moths'){
  tree.pglmm=read.nexus('Kawahara_2019_moth_tree_edited')
  ##Reading in a key file that maps the modified tree tip names from Kawahara et al. 2019 to my
  ##"ID_var" variable. It is already sorted and can be directly assigned to the tree
  name_keys=read.csv('Kawahara_renamed_moth_tree_species_keyfile.csv')
  names(name_keys)[names(name_keys) == paste(IDvar)] <- 'ID_var'
  tip.label.pglmm<-name_keys$ID_var
  tree.pglmm$tip.label<-as.character(tip.label.pglmm) #assign tip labels to the tree
  #removing trichopteran outgroups
  no_outgroups=name_keys[name_keys[,"ID_var"]>10,"ID_var"]
  tree.pglmm.no_outgroups<-drop.tip(tree.pglmm,tree.pglmm$tip.label[-match(as.character(no_outgroups), tree.pglmm$tip.label)])
  ##adding a root age for lepidoptera from Kawahara et al. 2019
  tree.pglmm.no_outgroups$root.time=299.5
  #Remove unused tips  ##note that it needs to interpret tip labels as CHARACTERS, NOT INTS
  keep.tips.d<-unique(all_metrics_MDF2[,"ID_var"])
  tree.pglmm.trimmed.d<-drop.tip(tree.pglmm.no_outgroups,tree.pglmm.no_outgroups$tip.label[-match(as.character(keep.tips.d), tree.pglmm.no_outgroups$tip.label)])
  phylopar_tree=force.ultrametric(tree.pglmm.trimmed.d)
  ###Rescaling the phylopar tree to have a height of 299.5 million years
  phylopar_tree_copy<-phylopar_tree
  phylopar_tree_copy$tip.label<-paste('tip',phylopar_tree_copy$tip.label,sep='_')
  scale=299.5
  phylopar_tree_copy$edge.length<-
    phylopar_tree_copy$edge.length/max(nodeHeights(phylopar_tree_copy)[,2])*scale
  edge_lengths=reorder(phylopar_tree_copy,order='cladewise')$edge.length
}


if(Taxon_group=='Butterflies'){
  tree_copy<-tree
  tip.label.pglmm<-as.character(unlist(c(1:Ntip(tree_copy)))) #create tip labels (1 to 195)
  tree_copy$tip.label<-tip.label.pglmm #assign tip labels to the tree
  keep.tips.d<-unique(all_metrics_MDF2[paste('ID_var')][[1]])
  tree.pglmm.trimmed.d<-drop.tip(tree_copy,tree_copy$tip.label[-match(keep.tips.d, tree_copy$tip.label)])
  H<-nodeHeights(tree.pglmm.trimmed.d)
  phylopar_tree=force.ultrametric(tree.pglmm.trimmed.d)
  phylopar_tree_copy<-phylopar_tree
  phylopar_tree_copy$tip.label<-paste('tip',phylopar_tree_copy$tip.label,sep='_')
  edge_lengths=reorder(phylopar_tree_copy,order='cladewise')$edge.length
}



if(Taxon_group=='Geometrids'){
  tree.pglmm=read.nexus('Geometridae_08_04_2019_rooted_uraniidae_sematuridae')
  ##Reading in a key file that maps the modified tree tip names from Murillo Ramos et al. 2019 to my
  ##"Geom_ID" variable. It is already sorted and can be directly assigned to the tree
  name_keys=read.csv('Geometridae_tips_from_Murillo_Ramos_etal_2019_Geometridae_08_04_2019_rooted_uraniidea_sematuridae_treefile_keyfile.csv')
  names(name_keys)[names(name_keys) == 'Geom_ID'] <- 'tip_ID'
  tip.label.pglmm<-name_keys$tip_ID
  tree.pglmm$tip.label<-as.character(tip.label.pglmm) #assign tip labels to the tree
  ##adding a root age for geometroidea (Geometridae+Uraniidae+Sematuridae) from run 3 of Kawahara et al. 2019 divergence time estimates table S11
  tree.pglmm$root.time=82.6
  #Remove unused tips  ##note that it needs to interpret tip labels as CHARACTERS, NOT INTS
  keep.tips.d<-unique(all_metrics_MDF2[,"tipID"])
  
  keep.tips.d
  tree.pglmm.trimmed.d<-drop.tip(tree.pglmm,tree.pglmm$tip.label[-match(as.character(keep.tips.d), tree.pglmm$tip.label)])
  phylopar_tree=force.ultrametric(tree.pglmm.trimmed.d)
  ###Rescaling the phylopar tree to have a height of 82.6 million years
  phylopar_tree_copy<-phylopar_tree
  phylopar_tree_copy$tip.label<-paste('tip',phylopar_tree_copy$tip.label,sep='_')
  scale=82.6
  phylopar_tree_copy$edge.length<-
    phylopar_tree_copy$edge.length/max(nodeHeights(phylopar_tree_copy)[,2])*scale
  edge_lengths=reorder(phylopar_tree_copy,order='cladewise')$edge.length
  root_time<-max(nodeHeights(phylopar_tree_copy))
}


mean_dtt<-function(input_metrics){
  metrics=input_metrics
  mean_vals <- vector(mode = "list", length = 0)
  DTT_stat_mean <- vector(mode = "list", length = 0)
  metrics2=metrics[,Metric_colnames]
  
  for (i in colnames(metrics2)){
    # print(i)
    try(index<-val_func(metrics,i,tree.pglmm.trimmed.d))
    try(mean_vals[[length(mean_vals) + 1]] <- index[[1]])
    try(DTT_stat_mean[[length(DTT_stat_mean) + 1]] <- index[[2]])
  }
  return(list(mean_vals,DTT_stat_mean))
}
results_MDF<-mean_dtt(all_metrics_MDF2)
results_MVF<-mean_dtt(all_metrics_MVF2)
results_MDH<-mean_dtt(all_metrics_MDH2)
results_MVH<-mean_dtt(all_metrics_MVH2)
results_FDF<-mean_dtt(all_metrics_FDF2)
results_FVF<-mean_dtt(all_metrics_FVF2)
results_FDH<-mean_dtt(all_metrics_FDH2)
results_FVH<-mean_dtt(all_metrics_FVH2)

meanframe_MDF<-cbind(meanframe_MDF,'mean_vals'=unlist(results_MDF[[1]]))
meanframe_MVF<-cbind(meanframe_MVF,'mean_vals'=unlist(results_MVF[[1]]))
meanframe_MDH<-cbind(meanframe_MDH,'mean_vals'=unlist(results_MDH[[1]]))
meanframe_MVH<-cbind(meanframe_MVH,'mean_vals'=unlist(results_MVH[[1]]))
meanframe_FDF<-cbind(meanframe_FDF,'mean_vals'=unlist(results_FDF[[1]]))
meanframe_FVF<-cbind(meanframe_FVF,'mean_vals'=unlist(results_FVF[[1]]))
meanframe_FDH<-cbind(meanframe_FDH,'mean_vals'=unlist(results_FDH[[1]]))
meanframe_FVH<-cbind(meanframe_FVH,'mean_vals'=unlist(results_FVH[[1]]))

names(meanframe_MDF)[names(meanframe_MDF) == 'mean_vals'] <- paste('boot',boot,sep='_')
names(meanframe_MVF)[names(meanframe_MVF) == 'mean_vals'] <- paste('boot',boot,sep='_')
names(meanframe_MDH)[names(meanframe_MDH) == 'mean_vals'] <- paste('boot',boot,sep='_')
names(meanframe_MVH)[names(meanframe_MVH) == 'mean_vals'] <- paste('boot',boot,sep='_')
names(meanframe_FDF)[names(meanframe_FDF) == 'mean_vals'] <- paste('boot',boot,sep='_')
names(meanframe_FVF)[names(meanframe_FVF) == 'mean_vals'] <- paste('boot',boot,sep='_')
names(meanframe_FDH)[names(meanframe_FDH) == 'mean_vals'] <- paste('boot',boot,sep='_')
names(meanframe_FVH)[names(meanframe_FVH) == 'mean_vals'] <- paste('boot',boot,sep='_')

dttframe_MDF<-cbind(dttframe_MDF,'dtt_vals'=unlist(results_MDF[[2]]))
dttframe_MVF<-cbind(dttframe_MVF,'dtt_vals'=unlist(results_MVF[[2]]))
dttframe_MDH<-cbind(dttframe_MDH,'dtt_vals'=unlist(results_MDH[[2]]))
dttframe_MVH<-cbind(dttframe_MVH,'dtt_vals'=unlist(results_MVH[[2]]))
dttframe_FDF<-cbind(dttframe_FDF,'dtt_vals'=unlist(results_FDF[[2]]))
dttframe_FVF<-cbind(dttframe_FVF,'dtt_vals'=unlist(results_FVF[[2]]))
dttframe_FDH<-cbind(dttframe_FDH,'dtt_vals'=unlist(results_FDH[[2]]))
dttframe_FVH<-cbind(dttframe_FVH,'dtt_vals'=unlist(results_FVH[[2]]))

names(dttframe_MDF)[names(dttframe_MDF) == 'dtt_vals'] <- paste('boot',boot,sep='_')
names(dttframe_MVF)[names(dttframe_MVF) == 'dtt_vals'] <- paste('boot',boot,sep='_')
names(dttframe_MDH)[names(dttframe_MDH) == 'dtt_vals'] <- paste('boot',boot,sep='_')
names(dttframe_MVH)[names(dttframe_MVH) == 'dtt_vals'] <- paste('boot',boot,sep='_')
names(dttframe_FDF)[names(dttframe_FDF) == 'dtt_vals'] <- paste('boot',boot,sep='_')
names(dttframe_FVF)[names(dttframe_FVF) == 'dtt_vals'] <- paste('boot',boot,sep='_')
names(dttframe_FDH)[names(dttframe_FDH) == 'dtt_vals'] <- paste('boot',boot,sep='_')
names(dttframe_FVH)[names(dttframe_FVH) == 'dtt_vals'] <- paste('boot',boot,sep='_')


Value<-'Mean'
mean_frame_out<-rbind(meanframe_MDF,meanframe_MVF,meanframe_MDH,meanframe_MVH,
                      meanframe_FDF,meanframe_FVF,meanframe_FDH,meanframe_FVH)
mean_frame_out<-cbind(Value,mean_frame_out)
Value<-'DTT'
dtt_frame_out<-rbind(dttframe_MDF,dttframe_MVF,dttframe_MDH,dttframe_MVH,
                     dttframe_FDF,dttframe_FVF,dttframe_FDH,dttframe_FVH)
dtt_frame_out<-cbind(Value,dtt_frame_out)


combined_frame<-rbind(mean_frame_out,dtt_frame_out)
print(boot)
##REPLACE THE PATH TO THE OUTPUT TO DESIRED LOCATION
write.csv(combined_frame,paste('/PATH/TO/boot',boot,prefix_str,'BS_sensitivity_analysis_output.csv',sep='_'))
end_time<-Sys.time()

end_time-start_time


