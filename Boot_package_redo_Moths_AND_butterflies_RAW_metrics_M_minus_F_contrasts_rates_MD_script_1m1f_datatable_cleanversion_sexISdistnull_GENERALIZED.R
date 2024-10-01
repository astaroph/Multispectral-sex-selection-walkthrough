library('geiger')
library('plyr')
library('ape')
library('phytools')
library('tidyr')
library('data.table')
library('boot.pval')
library('boot')
library('Rphylopars')
library('Rcpp')
library('stringr')
library('dispRity')


start_time<-Sys.time()

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

##First read in the arguments listed at the command line
args=(commandArgs(TRUE))
##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
  print("No arguments supplied.")
  ##supply default values
  metric_num<-1
  boot_reps<-1000
  n_cores<-1
  Taxon_group<-'Butterflies'
  categorical_variable='All'
  rdata<-'FINAL_raw_values_bootstrapping_overall_data_Butterflies_All_1m1f.RData'
  specimen_threshold=1
}else{
  for(i in 1:length(args)){
    # print(args)
    eval(parse(text=args[[i]]))
  }
}

if(rdata!='No_Data'){
  load(rdata)
  args=(commandArgs(TRUE))
  for(i in 1:length(args)){
    print(args)
    eval(parse(text=args[[i]]))
  }
}else{
  
  #########################################################Making the overall raw metrics data from scratch----
  ######################set parameters here----
  # Taxon_group<-'Eumaeini'
  # categorical_variable='brush'
  #specimen_t
  ##Loading in the PC metrics ----
  ##Currently works for both butterflies and moths and the butterflies and moths from Kawahara 2019 as well as for both categorical variables
  ##insert the taxon group as the first argument supplied to the sbatch command
  #choices are 'Butterflies' 'Moths' and experimentally 'Geometrids' but 'Eumaeini' still needs its own code
  ##Change to your own IDvar, start and stop columns and input dataframes as needed.
  if(Taxon_group=='Moths'){
    PC_input<-read.csv('Master_Moths_1m1f_sensitivitycut_Final_trait_dataset_1MF.csv')
    prefix_str='Moth_nocturnality_BS_raw'
    if(categorical_variable=='All'){
      PC_input<-PC_input[PC_input$Diel_behavior=='Nocturnal',]
      prefix_str='Moth_nocturnaL_BS_raw_FINAL'
    }
    IDvar<-'Moth_ID'
    PC_col_start=24
    PC_col_end=248
  }
  if(Taxon_group=='Butterflies'){
    PC_input<-read.csv('Master_Butterflies_1m1f_sensitivitycut_Final_trait_dataset_1MF.csv')
    prefix_str='Butterflies_BS_raw_FINAL'
    IDvar<-'Phy_ID'
    PC_col_start=21
    PC_col_end=224
  }
  
  if(Taxon_group=='Geometrids'){
    PC_input<-read.csv('Master_Geometrids_1m1f_sens_cut_daynightOnly_Final_trait_dataset_1MF.csv')
    prefix_str='Geometrids_BS_raw_Final'
    IDvar<-'Geom_ID'
    PC_col_start=21
    PC_col_end=240
  }
  
  #############SET THESE BEFORE RUNNING
  # categorical_variable='Diel_behavior'
  # categorical_variable='Family'
  # categorical_variable='All'
  # categorical_variable='Moth_butterfly'
  make_csvs=FALSE
  final_csvs=FALSE
  IS_var_thresh=3
  
  ###Add unnecessary variables to drop from the original dataframe here
  if(Taxon_group=='Moths'){
    # drop=c('UID','UID.1','Unnamed..0','Family_Side_Wing','Sex_Side_Wing','Diel_behavior_source',
    #        'Genus','Species','X1cm_scale_pixels_UV_940')
    # drop=c('UID','UID.1','Family_Side_Wing','Sex_Side_Wing','Diel_behavior_reference',
    #        'Genus','Species','X1cm_scale_pixels_UV_940','Moth_ID.1')
    drop=c('UID','Family_Side_Wing','Sex_Side_Wing','Diel_behavior_reference',
           'Genus','Species','X1cm_scale_pixels_UV_940','Moth_ID.1')
  }
  
  if(Taxon_group=='Butterflies'){
    # drop=c('UID','Family_Side_Wing','Sex_Side_Wing','Family_Side_Wing','Sex_Side_Wing',
    #        'Genus_valid','Identity_of_this_specimen','Phy_oldID','Swap_or_not','Loc_ID',
    #        'X1cm_scale_pixels_UV_940')
    drop=c(NULL)
    names(PC_input)[names(PC_input) == 'Specimen_Barcode'] <- 'Specimen'
  }
  
  if(Taxon_group=='Geometrids'){
    drop=c(NULL)
  }
  PC_all<-PC_input
  
  # #####function continues below----
  # #####making a boostrapped dataset by separately sampling the specimens in males and females WITH REPLACEMENT, then recombining the new dataset
  # ##before proceeding
  # 
  # all_metrics_MDF=subset(PC_all,Sex=='male'& Side=='dorsal' & Wing=='FW')
  # all_metrics_MVF=subset(PC_all,Sex=='male'& Side=='ventral' & Wing=='FW')
  # all_metrics_MDH=subset(PC_all,Sex=='male'& Side=='dorsal' & Wing=='HW')
  # all_metrics_MVH=subset(PC_all,Sex=='male'& Side=='ventral' & Wing=='HW')
  # 
  # #Use the following block of code to bootstrap such that each family has the same number of rows as the original data, but not necessarily the same number of species
  # # all_metrics_MDF <- ddply(all_metrics_MDF,.(Family),function(x) x[sample(nrow(x),nrow(x),replace=TRUE),])
  # # all_metrics_MVF <- ddply(all_metrics_MVF,.(Family),function(x) x[sample(nrow(x),nrow(x),replace=TRUE),])
  # # all_metrics_MDH <- ddply(all_metrics_MDH,.(Family),function(x) x[sample(nrow(x),nrow(x),replace=TRUE),])
  # # all_metrics_MVH <- ddply(all_metrics_MVH,.(Family),function(x) x[sample(nrow(x),nrow(x),replace=TRUE),])
  # # PC_all_M<-rbind(all_metrics_MDF,all_metrics_MVF,all_metrics_MDH,all_metrics_MVH)
  # 
  # ##Use the following block of code to bootstrap randomly across families
  # rownames(all_metrics_MDF)<-all_metrics_MDF$Specimen
  # rownames(all_metrics_MVF)<-all_metrics_MVF$Specimen
  # rownames(all_metrics_MDH)<-all_metrics_MDH$Specimen
  # rownames(all_metrics_MVH)<-all_metrics_MVH$Specimen
  # Specimens_M<-intersect(intersect(intersect(all_metrics_MDF$Specimen,all_metrics_MVF$Specimen),all_metrics_MDH$Specimen),all_metrics_MVH$Specimen)
  # 
  # Specimens_M1<-sample(Specimens_M,length(Specimens_M),replace=TRUE)
  # PC_all_M<-rbind(all_metrics_MDF[Specimens_M1,],all_metrics_MVF[Specimens_M1,],all_metrics_MDH[Specimens_M1,],all_metrics_MVH[Specimens_M1,])
  # 
  # all_metrics_FDF=subset(PC_all,Sex=='female'& Side=='dorsal' & Wing=='FW')
  # all_metrics_FVF=subset(PC_all,Sex=='female'& Side=='ventral' & Wing=='FW')
  # all_metrics_FDH=subset(PC_all,Sex=='female'& Side=='dorsal' & Wing=='HW')
  # all_metrics_FVH=subset(PC_all,Sex=='female'& Side=='ventral' & Wing=='HW')
  # 
  # #Use the following block of code to bootstrap such that each family has the same number of rows as the original data, but not necessarily the same number of species
  # # all_metrics_FDF <- ddply(all_metrics_FDF,.(Family),function(x) x[sample(nrow(x),nrow(x),replace=TRUE),])
  # # all_metrics_FVF <- ddply(all_metrics_FVF,.(Family),function(x) x[sample(nrow(x),nrow(x),replace=TRUE),])
  # # all_metrics_FDH <- ddply(all_metrics_FDH,.(Family),function(x) x[sample(nrow(x),nrow(x),replace=TRUE),])
  # # all_metrics_FVH <- ddply(all_metrics_FVH,.(Family),function(x) x[sample(nrow(x),nrow(x),replace=TRUE),])
  # # PC_all_F<-rbind(all_metrics_FDF,all_metrics_FVF,all_metrics_FDH,all_metrics_FVH)
  # 
  # ##Use the following block of code to bootstrap randomly across families
  # rownames(all_metrics_FDF)<-all_metrics_FDF$Specimen
  # rownames(all_metrics_FVF)<-all_metrics_FVF$Specimen
  # rownames(all_metrics_FDH)<-all_metrics_FDH$Specimen
  # rownames(all_metrics_FVH)<-all_metrics_FVH$Specimen
  # Specimens_F<-intersect(intersect(intersect(all_metrics_FDF$Specimen,all_metrics_FVF$Specimen),all_metrics_FDH$Specimen),all_metrics_FVH$Specimen)
  # Specimens_F1<-sample(Specimens_F,length(Specimens_F),replace=TRUE)
  # PC_all_F<-rbind(all_metrics_FDF[Specimens_F1,],all_metrics_FVF[Specimens_F1,],all_metrics_FDH[Specimens_F1,],all_metrics_FVH[Specimens_F1,])
  # 
  # PC_all<-rbind(PC_all_M,PC_all_F)
  rownames(PC_all)<-NULL
  if(categorical_variable=='All'){
    Categorical_var<-rep('All',length(nrow(PC_all)))
  }else{
    Categorical_var<-PC_all[paste(categorical_variable)]
  }
  
  
  PC_colnames<-colnames(PC_all[,PC_col_start:PC_col_end])
  names(PC_all)[names(PC_all) == paste(IDvar)] <- 'ID_var'
  names(Categorical_var)<-'Categorical_var'
  PC_all<-cbind(Categorical_var,PC_all)
  
  # names(PC_all)[names(PC_all) == paste(categorical_variable)] <- 'Categorical_var'
  PC_all=PC_all[,!names(PC_all)%in%drop]
  Phy_Side_Wing<-paste(PC_all$ID_var,PC_all$Side,PC_all$Wing,sep='_')
  PC_all<-cbind(Phy_Side_Wing,PC_all)
  
  PC<-ddply(PC_all,.(Phy_Side_Wing,Family,ID_var,Categorical_var,Sex,Side,Wing),colwise(median))
  PC_catnames<-c('Phy_Side_Wing','Family','ID_var','Categorical_var','Sex','Side','Wing')
  
  #####If the specimen threshold is set higher than 1, this will filter the input datasets to match
  if(specimen_threshold>1){
    PC_IS_count<-ddply(PC_all,.(Phy_Side_Wing,ID_var,Family,Categorical_var,Sex,Side,Wing),colwise(length))
    PC_IS_med<-ddply(PC_all,.(Phy_Side_Wing,ID_var,Family,Categorical_var,Sex,Side,Wing),colwise(median))
    cat_count<-ddply(PC_IS_med,.(Categorical_var,Sex,Side,Wing),colwise(length))
    cat_count_filt<-subset(cat_count,ID_var>1)
    PC_IS_med<-PC_IS_med[PC_IS_med$Categorical_var %in% unique(cat_count_filt$Categorical_var),]
    PC_IS_count<-PC_IS_count[PC_IS_count$Categorical_var %in% unique(cat_count_filt$Categorical_var),]
    PC_IS_med<-PC_IS_med[,c(PC_catnames,PC_colnames)]
    count<-PC_IS_count$Specimen
    PC_sub<-cbind(count,PC_IS_med)
    PC_sub<-PC_sub[PC_sub$count>(specimen_threshold-1),]
    M_IDs<-unique(PC_sub[PC_sub$Sex=='male',]$ID_var)
    F_IDs<-unique(PC_sub[PC_sub$Sex=='female',]$ID_var)
    final_ids<-intersect(M_IDs,F_IDs)
    print(final_ids)
    PC_all<-PC_all[PC_all$ID_var %in% final_ids,]
    PC<-ddply(PC_all,.(Phy_Side_Wing,Family,ID_var,Categorical_var,Sex,Side,Wing),colwise(median))
  }
  # PC<-subset(PC,select= -c(Specimen))
  PC_categories<-PC[PC_catnames]
  PC_metrics<-PC[,PC_colnames]
  PC<-cbind(PC_categories,PC_metrics)
  
  PC_start<-as.numeric(length(PC)-length(PC_colnames)+1)
  PC_end<-as.numeric(length(PC))
  
  ###Removing any data in levels of the categorical variable
  ##that only have a single species level data point 
  cat_count<-ddply(PC,.(Categorical_var,Sex,Side,Wing),colwise(length))
  cat_count_filt<-subset(cat_count,ID_var>1)
  PC<-PC[PC$Categorical_var %in% unique(cat_count_filt$Categorical_var),]
  
  # PC_diurnal<-unique(PC[PC$Categorical_var=='Diurnal',]$ID_var)
  # PC_nocturnal<-unique(PC[PC$Categorical_var=='Nocturnal',]$Family)
  # sort(PC_nocturnal)
  # PC<-read.csv('Moth_perspecimenPC_medians_Sep2021_noFluor_Pol_PCA_scoreframe_sensitivity_curated_spp_w_males_females.csv')
  
  ##function starts below. The first part runs the ancestral state reconstruction for the entire
  ##dataset supplied. The next part produces values for the unique levels of a supplied
  ##categorical variable to create a final plotting dataframe of the final PC metrics split by
  ##that variable
  # names(PC)[names(PC) == paste(categorical_variable)] <- 'Categorical_var'
  
  #####Now making the M-F sex specific rates ----
  ####Rphylopars final ancestral state reconstruction workflow:
  #I'm using Eric Goolsby's Rphylopar package, which can handle both missing data
  ##as well as within species variation and alternative evolutionary models, but is 
  #much faster than phytools anc.ml
  # This allows me to evaluate 5 different models for each trait and keep the ancestral
  ##states from the ones with the best fits to the actual data on the phylogeny
  ##importing the per specimen data, as Rphylopars can handle
  #within species variation in its ancestral state reconstruction
  all_metrics<-PC
  species<-all_metrics$ID_var
  all_metrics1<-cbind(species,all_metrics)
  
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
    keep.tips.d<-unique(all_metrics_MDF1[,"ID_var"])
    tree.pglmm.trimmed.d<-drop.tip(tree.pglmm.no_outgroups,tree.pglmm.no_outgroups$tip.label[-match(as.character(keep.tips.d), tree.pglmm.no_outgroups$tip.label)])
    phylopar_tree=force.ultrametric(tree.pglmm.trimmed.d)
    # phylopar_tree$tip.label<-paste('tip',phylopar_tree$tip.label,sep='_')
    # writeNexus(phylopar_tree,'phylopar_tree_forced_filtered_MF_by_phyID_ultrametric.nexus')
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
    keep.tips.d<-unique(all_metrics_MDF1[paste('ID_var')][[1]])
    tree.pglmm.trimmed.d<-drop.tip(tree_copy,tree_copy$tip.label[-match(keep.tips.d, tree_copy$tip.label)])
    H<-nodeHeights(tree.pglmm.trimmed.d)
    # tree.pglmm.trimmed.d$tip.label<-paste("t",tree.pglmm.trimmed.d$tip.label,sep='_') #assign tip labels to the tree
    phylopar_tree=force.ultrametric(tree.pglmm.trimmed.d)
    phylopar_tree_copy<-phylopar_tree
    phylopar_tree_copy$tip.label<-paste('tip',phylopar_tree_copy$tip.label,sep='_')
    edge_lengths=reorder(phylopar_tree_copy,order='cladewise')$edge.length
  }
  
  if(Taxon_group=='Geometrids'){
    # tree.pglmm=read.nexus('Geometridae_08_04_2019_rooted_uraniidae_sematuridae')
    # ##Reading in a key file that maps the modified tree tip names from Murillo Ramos et al. 2019 to my
    # ##"Geom_ID" variable. It is already sorted and can be directly assigned to the tree
    # name_keys=read.csv('Geometridae_tips_from_Murillo_Ramos_etal_2019_Geometridae_08_04_2019_rooted_uraniidea_sematuridae_treefile_keyfile.csv')
    # names(name_keys)[names(name_keys) == 'Geom_ID'] <- 'tip_ID'
    # tip.label.pglmm<-name_keys$tip_ID
    # tree.pglmm$tip.label<-as.character(tip.label.pglmm) #assign tip labels to the tree
    # ##adding a root age for geometroidea (Geometridae+Uraniidae+Sematuridae) from run 3 of Kawahara et al. 2019 divergence time estimates table S11
    # tree.pglmm$root.time=82.6
    # #Remove unused tips  ##note that it needs to interpret tip labels as CHARACTERS, NOT INTS
    # keep.tips.d<-unique(all_metrics_MDF1[paste('ID_var')][[1]])
    # tree.pglmm.trimmed.d<-drop.tip(tree.pglmm,tree.pglmm$tip.label[-match(as.character(keep.tips.d), tree.pglmm$tip.label)])
    # phylopar_tree=force.ultrametric(tree.pglmm.trimmed.d)
    # ###Rescaling the phylopar tree to have a height of 78.2 million years
    # phylopar_tree_copy<-phylopar_tree
    # phylopar_tree_copy$tip.label<-paste('tip',phylopar_tree_copy$tip.label,sep='_')
    # scale=82.6
    # phylopar_tree_copy$edge.length<-
    #   phylopar_tree_copy$edge.length/max(nodeHeights(phylopar_tree_copy)[,2])*scale
    # edge_lengths=reorder(phylopar_tree_copy,order='cladewise')$edge.length
    # root_time<-max(nodeHeights(phylopar_tree_copy))
    name_keys=read.csv('Geometridae_tips_from_Murillo_Ramos_etal_2019_Geometridae_08_04_2019_rooted_uraniidea_sematuridae_treefile_keyfile.csv')
    names(name_keys)[names(name_keys) == 'Geom_ID'] <- 'tip_ID'
    tree.pglmm<-read.tree("Geometridae_Murillo_Ramos2019_Geometridae_08_04_2019_rooted_uraniidae_sematuridae_timetree_relaxed_lambda0_1.newick")
    keep.tips.d<-unique(all_metrics_MDF1[paste('ID_var')][[1]])
    tree.pglmm.trimmed.d<-drop.tip(tree.pglmm,tree.pglmm$tip.label[-match(as.character(keep.tips.d), tree.pglmm$tip.label)])
    phylopar_tree<-tree.pglmm.trimmed.d
    phylopar_tree_copy<-tree.pglmm.trimmed.d
    phylopar_tree_copy$tip.label<-paste('tip',phylopar_tree_copy$tip.label,sep='_')
    edge_lengths=reorder(phylopar_tree_copy,order='cladewise')$edge.length
    root_time<-max(nodeHeights(phylopar_tree_copy))
  }
  #####
  #testcode
  # trialmat<-all_metrics_MDF1[,10:24]
  # species<-paste('tip',all_metrics_MDF1$ID_var,sep='_')
  # trialmat1<-cbind(species,trialmat)
  # A=cbind(trialmat1[1],trialmat1[2])
  # trial=phylopars(A,phylopar_tree_copy,model='BM')
  # edges=trial$anc_recon
  # rownames(edges)<-as.character(rownames(edges))
  # edgematrix=data.frame('edgename'=rownames(edges))
  #####
  phylopars_multi_anc<-function(fullmetrics,treefile,contextstr,metric_start,metrics_end){
    trialmat<-fullmetrics[,metric_start:metrics_end]
    species<-paste('tip',fullmetrics$ID_var,sep='_')
    trialmat1<-cbind(species,trialmat)
    A=cbind(trialmat1[1],trialmat1[2])
    trial=phylopars(A,treefile,model='BM')
    edges=trial$anc_recon
    rownames(edges)<-as.character(rownames(edges))
    edgematrix=data.frame('edgename'=rownames(edges))
    metriclist=vector()
    modellist=vector()
    for(i in seq(2,length(trialmat1))){
      # print(names(trialmat1[i]))
      metriclist[[length(metriclist)+1]]<-names(trialmat1[i])
      A=cbind(trialmat1[1],trialmat1[i])
      if(var(A[[names(trialmat1[i])]])==0){
        modellist[[length(modellist)+1]]<-'No model: VAR=0'
        edgematrix[names(trialmat1[i])]<-mean(A[[names(trialmat1[i])]])
        # print('No model: VAR=0')
      } else {
        BM=phylopars(A,treefile,model='BM')
        OU=phylopars(A,treefile,model='OU')
        EB=phylopars(A,treefile,model='EB')
        lambda=phylopars(A,treefile,model='lambda')
        kappa=phylopars(A,treefile,model='kappa')
        combine=data.frame('BM'=BM$anc_recon,'OU'=OU$anc_recon,'EB'=EB$anc_recon,
                           'lambda'=lambda$anc_recon,'kappa'=kappa$anc_recon)
        akaike=AIC(BM,OU,EB,lambda,kappa)
        ind=which.min(akaike[,2])
        model=rownames(akaike)[ind]
        modellist[[length(modellist)+1]]<-model
        edgematrix[names(trialmat1[i])]<-combine[ind]  
      }
    }
    modelframe=data.frame('Context'=contextstr,'Metric'=metriclist,"Model"=modellist)
    return(list(modelframe,edgematrix))
  }
  # phylopars_multi_anc(all_metrics_FVH1,phylopar_tree,'FVH',(PC_start+1),(PC_end+1))
  results_MDF<-phylopars_multi_anc(all_metrics_MDF1,phylopar_tree_copy,'MDF',PC_start+1,PC_end+1)
  print('Done with Anc_state MDF')
  results_MVF<-phylopars_multi_anc(all_metrics_MVF1,phylopar_tree_copy,'MVF',PC_start+1,PC_end+1)
  print('Done with Anc_state MVF')
  results_FDF<-phylopars_multi_anc(all_metrics_FDF1,phylopar_tree_copy,'FDF',PC_start+1,PC_end+1)
  print('Done with Anc_state FDF')
  results_FVF<-phylopars_multi_anc(all_metrics_FVF1,phylopar_tree_copy,'FVF',PC_start+1,PC_end+1)
  print('Done with Anc_state FVF')
  results_MDH<-phylopars_multi_anc(all_metrics_MDH1,phylopar_tree_copy,'MDH',PC_start+1,PC_end+1)
  print('Done with Anc_state MDH')
  results_MVH<-phylopars_multi_anc(all_metrics_MVH1,phylopar_tree_copy,'MVH',PC_start+1,PC_end+1)
  print('Done with Anc_state MVH')
  results_FDH<-phylopars_multi_anc(all_metrics_FDH1,phylopar_tree_copy,'FDH',PC_start+1,PC_end+1)
  print('Done with Anc_state FDH')
  results_FVH<-phylopars_multi_anc(all_metrics_FVH1,phylopar_tree_copy,'FVH',PC_start+1,PC_end+1)
  print('Done with Anc_state FVH')
  modelframe_MDF=results_MDF[[1]]
  modelframe_MVF=results_MVF[[1]]
  modelframe_FDF=results_FDF[[1]]
  modelframe_FVF=results_FVF[[1]]
  modelframe_MDH=results_MDH[[1]]
  modelframe_MVH=results_MVH[[1]]
  modelframe_FDH=results_FDH[[1]]
  modelframe_FVH=results_FVH[[1]]
  
  edgematrix_MDF=results_MDF[[2]]
  edgematrix_MVF=results_MVF[[2]]
  edgematrix_FDF=results_FDF[[2]]
  edgematrix_FVF=results_FVF[[2]]
  edgematrix_MDH=results_MDH[[2]]
  edgematrix_MVH=results_MVH[[2]]
  edgematrix_FDH=results_FDH[[2]]
  edgematrix_FVH=results_FVH[[2]]
  
  edgematrix_MDF1<-edgematrix_MDF
  edgematrix_MVF1<-edgematrix_MVF
  edgematrix_FDF1<-edgematrix_FDF
  edgematrix_FVF1<-edgematrix_FVF
  edgematrix_MDH1<-edgematrix_MDH
  edgematrix_MVH1<-edgematrix_MVH
  edgematrix_FDH1<-edgematrix_FDH
  edgematrix_FVH1<-edgematrix_FVH
  
  if(Taxon_group=='Moths'){
    tree_copy<-tree.pglmm
    # tree_copy$tip.label<-paste('t',tree_copy$tip.label,sep='_')
  }
  if(Taxon_group=='Moths_and_butterflies'){
    tree_copy<-tree.pglmm
    # tree_copy$tip.label<-paste('t',tree_copy$tip.label,sep='_')
  }
  if(Taxon_group=='Butterflies'){
    tree_copy<-tree
    # tree_copy$tip.label<-paste('t',tree_copy$tip.label,sep='_')
    # tipmatch=matchLabels(tree_copy,phylopar_tree)
    # tipclean=na.omit(tipmatch)
    # pruned_ids=tipclean[,1]
    # names(pruned_ids)<-NULL
    # name_matrix=reorder(phylopar_tree,order='cladewise')$edge
    # name_matrix[,1]=name_matrix[,1]+10
    # name_matrix[,2][name_matrix[,2]>185]<-name_matrix[,2][name_matrix[,2]>185]+10
    # name_matrix[,2][name_matrix[,2]<186]<-paste('t',pruned_ids,sep='_')
    # name_matrix
  }
  if(Taxon_group=='Eumaeini'){
    tree_copy<-tree.pglmm
    # tree_copy$tip.label<-paste('t',tree_copy$tip.label,sep='_')
  }
  
  if(Taxon_group=='Geometrids'){
    tree_copy<-tree.pglmm
  }
  
  tipmatch=matchLabels(tree_copy,phylopar_tree)
  tipclean=na.omit(tipmatch)
  rownames(tipclean)<-as.character(rownames(tipclean))
  pruned_ids=tipclean[,1]
  names(pruned_ids)<-NULL
  name_matrix=reorder(phylopar_tree,order='cladewise')$edge
  rownames(name_matrix)<-name_matrix[,2]
  tipnum=length(phylopar_tree$tip.label)
  # name_matrix[,2][(as.numeric(rownames(name_matrix))>tipnum)]<-paste('anc',name_matrix[,2][(as.numeric(rownames(name_matrix))>tipnum)],sep='_')
  name_matrix[,2][(as.numeric(rownames(name_matrix))<=tipnum)]<-paste('tip',pruned_ids,sep='_')
  
  #####
  #testcode
  # name_matrix[,1]=name_matrix[,1]+10
  # name_matrix[,2][name_matrix[,2]>16]<-name_matrix[,2][name_matrix[,2]>185]+10
  # name_matrix[,2][name_matrix[,2]<tipnum+1]<-pruned_ids
  # name_matrix[,2][name_matrix[,2]>tipnum]<-paste('anc',name_matrix[,2][name_matrix[,2]>tipnum],sep='_')
  #####
  
  if(make_csvs==TRUE){
    write.csv(name_matrix,paste(prefix_str,'phylopar_tree_node_tip_name_matrix.csv',sep='_'))
  }
  rownames(edgematrix_MDF1)<-edgematrix_MDF1$edgename
  rownames(edgematrix_MVF1)<-edgematrix_MVF1$edgename
  rownames(edgematrix_FDF1)<-edgematrix_FDF1$edgename
  rownames(edgematrix_FVF1)<-edgematrix_FVF1$edgename
  
  rownames(edgematrix_MDH1)<-edgematrix_MDH1$edgename
  rownames(edgematrix_MVH1)<-edgematrix_MVH1$edgename
  rownames(edgematrix_FDH1)<-edgematrix_FDH1$edgename
  rownames(edgematrix_FVH1)<-edgematrix_FVH1$edgename
  
  #####Sanity check code- hide----
  # edgelengths<-as.data.frame(name_matrix)
  # edgelengths$edgel<-edge_lengths
  # colnames(edgelengths)<-c('Ancestors','Descendants','Edge_lengths')
  # # edgelengths
  # rownames(edgelengths)<-edgelengths$Descendants
  # descendants=sort(edgelengths$Descendants)[15:30]
  # # sort(edgelengths$Descendants)[15:30]
  # #
  # edgelengths_tips<-edgelengths[descendants,]
  # edgelengths[descendants,]
  # #
  # #
  # edgemat_m=edgematrix_MVH1
  # edgemat_f=edgematrix_FVH1
  # metricsM=all_metrics_MVH1
  # metricsF=all_metrics_FVH1
  # name_matrix[,2]
  # edgemat_m[name_matrix[,2],'X0']
  # #generating a 2 column matrix of ancestral and descendant trait values for each edge for males and females
  # desc=edgemat_m[name_matrix[,2],'X0']
  # names(desc)<-edgemat_m[name_matrix[,2],'edgename']
  # anc=edgemat_m[name_matrix[,1],'X0']
  # names(anc)<-edgemat_m[name_matrix[,1],'edgename']
  # T_m=cbind(anc,desc)
  # rownames(T_m)<-names(desc)
  # dT_m=T_m[,2]-T_m[,1]
  # rate_m=dT_m/edge_lengths
  # desc=edgemat_f[name_matrix[,2],'X0']
  # names(desc)<-edgemat_f[name_matrix[,2],'edgename']
  # anc=edgemat_f[name_matrix[,1],'X0']
  # names(anc)<-edgemat_f[name_matrix[,1],'edgename']
  # T_f=cbind(anc,desc)
  # rownames(T_f)<-names(desc)
  # dT_f=T_f[,2]-T_f[,1]
  # rate_f=dT_f/edge_lengths
  # 
  # ################Creating sex specific 2-column edge dimorphism matrices for each ancestor or descendant
  # D_m=T_m-T_f
  # D_f=T_f-T_m
  # 
  # ##Making sex specific delta D matrices that describe the change in dimorphism
  # dD_m=D_m[,2]-D_m[,1]
  # dD_f=D_f[,2]-D_f[,1]
  # 
  # ####Creating "alignment" dataframes, which are positive values when changes in traits align with
  # ##the direction (positive or negative) of dimorphism, and negative when not aligned
  # dT_D_m=D_m*dT_m
  # dT_D_f=D_f*dT_f
  # ##Creating filter lists that contain node or tip IDs to call for assigning final contribution values
  # ##aligned IDs
  # D_plus_m=names(dT_D_m[,2][dT_D_m[,2]>0])
  # D_plus_f=names(dT_D_f[,2][dT_D_f[,2]>0])
  # 
  # ##unaligned IDs
  # D_minus_m=names(dT_D_m[,2][dT_D_m[,2]<0])
  # D_minus_f=names(dT_D_f[,2][dT_D_f[,2]<0])
  # ##IDs where D descendant==0
  # D_zero_m=names(dT_D_m[,2][dT_D_m[,2]==0])
  # D_zero_f=names(dT_D_f[,2][dT_D_f[,2]==0])
  # ##IDs where D ancestral==0
  # Da_zero_m=names(dT_D_m[,2][dT_D_m[,1]==0])
  # Da_zero_f=names(dT_D_f[,2][dT_D_f[,1]==0])
  # ###IDs where D descendant==0 and D ancestral is not zero
  # Da_m=setdiff(D_zero_m,Da_zero_m)
  # Da_f=setdiff(D_zero_f,Da_zero_f)
  # ##Ids where dT and ancestral dimorphism align
  # Da_plus_m=intersect(names(dT_D_m[,2][dT_D_m[,1]>0]),Da_m)
  # Da_plus_f=intersect(names(dT_D_f[,2][dT_D_f[,1]>0]),Da_f)
  # ##Ids where dT and ancestral dimorphism dont align
  # Da_minus_m=intersect(names(dT_D_m[,2][dT_D_m[,1]<0]),Da_m)
  # Da_minus_f=intersect(names(dT_D_f[,2][dT_D_f[,1]<0]),Da_m)
  # ###Ids where both D ancestral and D descendant are zero
  # all_zero_m=intersect(D_zero_m,Da_zero_m)
  # all_zero_f=intersect(D_zero_f,Da_zero_f)
  # 
  # #####Making the final sex-specific contributions vector, one value per descendant
  # ##Starting with the raw change in trait values
  # C_frame_m<-rate_m
  # C_frame_f<-rate_f
  # ##making aligned contributions equal to absolute values of trait changes
  # C_frame_m[D_plus_m]<-abs(C_frame_m[D_plus_m])
  # C_frame_f[D_plus_f]<-abs(C_frame_f[D_plus_f])
  # ##making unaligned contributions equal to negative absolute values of trait changes
  # C_frame_m[D_minus_m]<--abs(C_frame_m[D_minus_m])
  # C_frame_f[D_minus_f]<--abs(C_frame_f[D_minus_f])
  # ##For cases where descendant dimorphism is zero, making aligned contributions equal to absolute values of trait changes
  # C_frame_m[Da_plus_m]<-abs(C_frame_m[Da_plus_m])
  # C_frame_f[Da_plus_f]<-abs(C_frame_f[Da_plus_f])
  # ##For cases where descendant dimorphism is zero, making unaligned contributions equal to negtaive absolute values of trait changes
  # C_frame_m[Da_minus_m]<--abs(C_frame_m[Da_minus_m])
  # C_frame_f[Da_minus_f]<--abs(C_frame_f[Da_minus_f])
  # ##For cases where ancestral and descendant dimorphism is zero, contribution is zero
  # C_frame_m[all_zero_m]<-0
  # C_frame_f[all_zero_f]<-0
  # 
  # # ##this plot makes a little bit of sense:
  # # #A plot of the
  # plot(abs(D_m[,2]),log(abs(rate_m)/abs(rate_f)))
  # # plot(dD_m,rate_m/rate_f)
  #####
  
  #####The main rates, condtribution and change in dimorphism function----
  contributions<-function(M_metrics,F_metrics,metric_str,M_edgemat,F_edgemat,name_mat,edge_l){
    ##generating a 2 column matrix of ancestral and descendant trait values for each edge for males and females 
    desc=M_edgemat[name_mat[,2],metric_str]
    if(length(desc[!is.na(desc)])==0){
      len<-length(M_edgemat[name_mat[,2],metric_str])
      C_frame_m<-rep(0,len)
      C_frame_f<-rep(0,len)
      dD_m<-rep(0,len)
      dD_f<-rep(0,len)
      rate_m<-rep(0,len)
      names(rate_m)<-M_edgemat[name_mat[,2],'edgename']
      rate_f<-rep(0,len)
      names(rate_f)<-F_edgemat[name_mat[,2],'edgename']
    } else{
      names(desc)<-M_edgemat[name_mat[,2],'edgename']
      anc=M_edgemat[name_mat[,1],metric_str]
      names(anc)<-M_edgemat[name_mat[,1],'edgename']
      T_m=cbind(anc,desc)
      rownames(T_m)<-names(desc)
      dT_m=T_m[,2]-T_m[,1]
      rate_m=dT_m/edge_l
      ###
      desc=F_edgemat[name_mat[,2],metric_str]
      if(length(desc[!is.na(desc)])==0){
        len<-length(F_edgemat[name_mat[,2],metric_str])
        C_frame_m<-rep(0,len)
        C_frame_f<-rep(0,len)
        dD_m<-rep(0,len)
        dD_f<-rep(0,len)
        rate_m<-rate_m
        rate_f<-rep(0,len)
        names(rate_f)<-F_edgemat[name_mat[,2],'edgename']
      } else{
        names(desc)<-F_edgemat[name_mat[,2],'edgename']
        anc=F_edgemat[name_mat[,1],metric_str]
        names(anc)<-F_edgemat[name_mat[,1],'edgename']
        T_f=cbind(anc,desc)
        rownames(T_f)<-names(desc)
        dT_f=T_f[,2]-T_f[,1]
        rate_f=dT_f/edge_l
        ################Creating sex specific 2-column edge dimorphism matrices for each ancestor or descendant
        D_m=T_m-T_f
        D_f=T_f-T_m
        dD_m=(abs(D_m[,2])-abs(D_m[,1]))/edge_l
        dD_f=(abs(D_f[,2])-abs(D_f[,1]))/edge_l
        ####Creating "alignment" dataframes, which are positive values when changes in traits align with
        ##the direction (positive or negative) of dimorphism, and negative when not aligned
        dT_D_m=D_m*dT_m
        dT_D_f=D_f*dT_f
        ##Creating filter lists that contain node or tip IDs to call for assigning final contribution values
        ##aligned IDs
        D_plus_m=names(dT_D_m[,2][dT_D_m[,2]>0])
        D_plus_f=names(dT_D_f[,2][dT_D_f[,2]>0])
        ##unaligned IDs
        D_minus_m=names(dT_D_m[,2][dT_D_m[,2]<0])
        D_minus_f=names(dT_D_f[,2][dT_D_f[,2]<0])
        ##IDs where D descendant==0
        D_zero_m=names(dT_D_m[,2][dT_D_m[,2]==0])
        D_zero_f=names(dT_D_f[,2][dT_D_f[,2]==0])
        ##IDs where D ancestral==0
        Da_zero_m=names(dT_D_m[,2][dT_D_m[,1]==0])
        Da_zero_f=names(dT_D_f[,2][dT_D_f[,1]==0])
        ###IDs where D descendant==0 and D ancestral is not zero
        Da_m=setdiff(D_zero_m,Da_zero_m)
        Da_f=setdiff(D_zero_f,Da_zero_f)
        ##Ids where dT and ancestral dimorphism align
        Da_plus_m=intersect(names(dT_D_m[,2][dT_D_m[,1]>0]),Da_m)
        Da_plus_f=intersect(names(dT_D_f[,2][dT_D_f[,1]>0]),Da_f)
        ##Ids where dT and ancestral dimorphism dont align
        Da_minus_m=intersect(names(dT_D_m[,2][dT_D_m[,1]<0]),Da_m)
        Da_minus_f=intersect(names(dT_D_f[,2][dT_D_f[,1]<0]),Da_m)
        ###Ids where both D ancestral and D descendant are zero
        all_zero_m=intersect(D_zero_m,Da_zero_m)
        all_zero_f=intersect(D_zero_f,Da_zero_f)
        #####Making the final sex-specific contributions vector, one value per descendant 
        ##Starting with the raw change in trait values
        C_frame_m<-rate_m
        C_frame_f<-rate_f
        ##making aligned contributions equal to absolute values of trait changes 
        C_frame_m[D_plus_m]<-abs(C_frame_m[D_plus_m])
        C_frame_f[D_plus_f]<-abs(C_frame_f[D_plus_f])
        ##making unaligned contributions equal to negtaive absolute values of trait changes 
        C_frame_m[D_minus_m]<--abs(C_frame_m[D_minus_m])
        C_frame_f[D_minus_f]<--abs(C_frame_f[D_minus_f])
        ##For cases where descendant dimorphism is zero, making aligned contributions equal to absolute values of trait changes 
        C_frame_m[Da_plus_m]<-abs(C_frame_m[Da_plus_m])
        C_frame_f[Da_plus_f]<-abs(C_frame_f[Da_plus_f])
        ##For cases where descendant dimorphism is zero, making unaligned contributions equal to negtaive absolute values of trait changes 
        C_frame_m[Da_minus_m]<--abs(C_frame_m[Da_minus_m])
        C_frame_f[Da_minus_f]<--abs(C_frame_f[Da_minus_f])
        ##For cases where ancestral and descendant dimorphism is zero, contribution is zero
        C_frame_m[all_zero_m]<-0
        C_frame_f[all_zero_f]<-0
      }
    }
    return(list(C_frame_m,C_frame_f,dD_m,dD_f,rate_m,rate_f))
  }
  
  # trial=contributions(all_metrics_MDF1,all_metrics_FDF1,'X0',edgematrix_MDF1,edgematrix_FDF1,name_matrix,edge_lengths)
  
  c_dataframe_MDF=data.frame('Ancestral_node'=name_matrix[,1],'Descendant_node'=name_matrix[,2],'Sex'='male','Side'='dorsal','Wing'='FW')
  c_dataframe_MVF=data.frame('Ancestral_node'=name_matrix[,1],'Descendant_node'=name_matrix[,2],'Sex'='male','Side'='ventral','Wing'='FW')
  c_dataframe_FDF=data.frame('Ancestral_node'=name_matrix[,1],'Descendant_node'=name_matrix[,2],'Sex'='female','Side'='dorsal','Wing'='FW')
  c_dataframe_FVF=data.frame('Ancestral_node'=name_matrix[,1],'Descendant_node'=name_matrix[,2],'Sex'='female','Side'='ventral','Wing'='FW')
  c_dataframe_MDH=data.frame('Ancestral_node'=name_matrix[,1],'Descendant_node'=name_matrix[,2],'Sex'='male','Side'='dorsal','Wing'='HW')
  c_dataframe_MVH=data.frame('Ancestral_node'=name_matrix[,1],'Descendant_node'=name_matrix[,2],'Sex'='male','Side'='ventral','Wing'='HW')
  c_dataframe_FDH=data.frame('Ancestral_node'=name_matrix[,1],'Descendant_node'=name_matrix[,2],'Sex'='female','Side'='dorsal','Wing'='HW')
  c_dataframe_FVH=data.frame('Ancestral_node'=name_matrix[,1],'Descendant_node'=name_matrix[,2],'Sex'='female','Side'='ventral','Wing'='HW')
  
  dD_dataframe_MDF=data.frame('Ancestral_node'=name_matrix[,1],'Descendant_node'=name_matrix[,2],'Sex'='male','Side'='dorsal','Wing'='FW')
  dD_dataframe_MVF=data.frame('Ancestral_node'=name_matrix[,1],'Descendant_node'=name_matrix[,2],'Sex'='male','Side'='ventral','Wing'='FW')
  dD_dataframe_FDF=data.frame('Ancestral_node'=name_matrix[,1],'Descendant_node'=name_matrix[,2],'Sex'='female','Side'='dorsal','Wing'='FW')
  dD_dataframe_FVF=data.frame('Ancestral_node'=name_matrix[,1],'Descendant_node'=name_matrix[,2],'Sex'='female','Side'='ventral','Wing'='FW')
  dD_dataframe_MDH=data.frame('Ancestral_node'=name_matrix[,1],'Descendant_node'=name_matrix[,2],'Sex'='male','Side'='dorsal','Wing'='HW')
  dD_dataframe_MVH=data.frame('Ancestral_node'=name_matrix[,1],'Descendant_node'=name_matrix[,2],'Sex'='male','Side'='ventral','Wing'='HW')
  dD_dataframe_FDH=data.frame('Ancestral_node'=name_matrix[,1],'Descendant_node'=name_matrix[,2],'Sex'='female','Side'='dorsal','Wing'='HW')
  dD_dataframe_FVH=data.frame('Ancestral_node'=name_matrix[,1],'Descendant_node'=name_matrix[,2],'Sex'='female','Side'='ventral','Wing'='HW')
  
  r_dataframe_MDF=data.frame('Ancestral_node'=name_matrix[,1],'Descendant_node'=name_matrix[,2],'Sex'='male','Side'='dorsal','Wing'='FW')
  r_dataframe_MVF=data.frame('Ancestral_node'=name_matrix[,1],'Descendant_node'=name_matrix[,2],'Sex'='male','Side'='ventral','Wing'='FW')
  r_dataframe_FDF=data.frame('Ancestral_node'=name_matrix[,1],'Descendant_node'=name_matrix[,2],'Sex'='female','Side'='dorsal','Wing'='FW')
  r_dataframe_FVF=data.frame('Ancestral_node'=name_matrix[,1],'Descendant_node'=name_matrix[,2],'Sex'='female','Side'='ventral','Wing'='FW')
  r_dataframe_MDH=data.frame('Ancestral_node'=name_matrix[,1],'Descendant_node'=name_matrix[,2],'Sex'='male','Side'='dorsal','Wing'='HW')
  r_dataframe_MVH=data.frame('Ancestral_node'=name_matrix[,1],'Descendant_node'=name_matrix[,2],'Sex'='male','Side'='ventral','Wing'='HW')
  r_dataframe_FDH=data.frame('Ancestral_node'=name_matrix[,1],'Descendant_node'=name_matrix[,2],'Sex'='female','Side'='dorsal','Wing'='HW')
  r_dataframe_FVH=data.frame('Ancestral_node'=name_matrix[,1],'Descendant_node'=name_matrix[,2],'Sex'='female','Side'='ventral','Wing'='HW')
  
  for (i in PC_colnames){
    # print(i)
    trial_DF=contributions(all_metrics_MDF1,all_metrics_FDF1,i,edgematrix_MDF1,edgematrix_FDF1,name_matrix,edge_lengths)
    trial_VF=contributions(all_metrics_MVF1,all_metrics_FVF1,i,edgematrix_MVF1,edgematrix_FVF1,name_matrix,edge_lengths)
    trial_DH=contributions(all_metrics_MDH1,all_metrics_FDH1,i,edgematrix_MDH1,edgematrix_FDH1,name_matrix,edge_lengths)
    trial_VH=contributions(all_metrics_MVH1,all_metrics_FVH1,i,edgematrix_MVH1,edgematrix_FVH1,name_matrix,edge_lengths)
    c_dataframe_MDF[paste(i)]<-trial_DF[[1]]
    c_dataframe_FDF[paste(i)]<-trial_DF[[2]]
    c_dataframe_MVF[paste(i)]<-trial_VF[[1]]
    c_dataframe_FVF[paste(i)]<-trial_VF[[2]]
    c_dataframe_MDH[paste(i)]<-trial_DH[[1]]
    c_dataframe_FDH[paste(i)]<-trial_DH[[2]]
    c_dataframe_MVH[paste(i)]<-trial_VH[[1]]
    c_dataframe_FVH[paste(i)]<-trial_VH[[2]]
    dD_dataframe_MDF[paste(i)]<-trial_DF[3]
    dD_dataframe_FDF[paste(i)]<-trial_DF[4]
    dD_dataframe_MVF[paste(i)]<-trial_VF[3]
    dD_dataframe_FVF[paste(i)]<-trial_VF[4]
    dD_dataframe_MDH[paste(i)]<-trial_DH[3]
    dD_dataframe_FDH[paste(i)]<-trial_DH[4]
    dD_dataframe_MVH[paste(i)]<-trial_VH[3]
    dD_dataframe_FVH[paste(i)]<-trial_VH[4]
    r_dataframe_MDF[paste(i)]<-trial_DF[5]
    r_dataframe_FDF[paste(i)]<-trial_DF[6]
    r_dataframe_MVF[paste(i)]<-trial_VF[5]
    r_dataframe_FVF[paste(i)]<-trial_VF[6]
    r_dataframe_MDH[paste(i)]<-trial_DH[5]
    r_dataframe_FDH[paste(i)]<-trial_DH[6]
    r_dataframe_MVH[paste(i)]<-trial_VH[5]
    r_dataframe_FVH[paste(i)]<-trial_VH[6]
  }
  
  contributions_combined=rbind(c_dataframe_MDF,c_dataframe_FDF,c_dataframe_MVF,c_dataframe_FVF,
                               c_dataframe_MDH,c_dataframe_FDH,c_dataframe_MVH,c_dataframe_FVH)
  
  rates_combined=rbind(r_dataframe_MDF,r_dataframe_FDF,r_dataframe_MVF,r_dataframe_FVF,
                       r_dataframe_MDH,r_dataframe_FDH,r_dataframe_MVH,r_dataframe_FVH)
  
  deltaD_combined=rbind(dD_dataframe_MDF,dD_dataframe_FDF,dD_dataframe_MVF,dD_dataframe_FVF,
                        dD_dataframe_MDH,dD_dataframe_FDH,dD_dataframe_MVH,dD_dataframe_FVH)
  if(make_csvs==TRUE){
    write.csv(contributions_combined,paste(prefix_str,'Sex_specific_contributions_dataset_perspecimenPCmedians_combined.csv',sep='_'))
    write.csv(deltaD_combined,paste(prefix_str,'Sex_specific_deltaD_dataset_perspecimenPCmedians_combined.csv',sep='_'))
    write.csv(rates_combined,paste(prefix_str,'Sex_specific_rates_dataset_perspecimenPCmedians_combined.csv',sep='_'))
  }
  # write.csv(contributions_combined,'MOTHS_Sex_specific_contributions_dataset_perspecimenPCmedians_combined.csv')
  # write.csv(deltaD_combined,'MOTHS_Sex_specific_deltaD_dataset_perspecimenPCmedians_combined.csv')
  # write.csv(rates_combined,'MOTHS_Sex_specific_rates_dataset_perspecimenPCmedians_combined.csv')
  
  ##Now subsetting to descendants_only
  tips_only<-function(inputframe){
    rownames(inputframe)<-as.character(inputframe$Descendant_node)
    # descendants=as.character(sort(as.numeric(inputframe$Descendant_node))[15:30])
    # outframe<-inputframe[descendants,]
    outframe<-inputframe[str_detect(inputframe$Descendant_node,'tip'),]
    return(outframe)
  }
  
  c_dataframe_MDF_d<-tips_only(c_dataframe_MDF)
  c_dataframe_MVF_d<-tips_only(c_dataframe_MVF)
  c_dataframe_FDF_d<-tips_only(c_dataframe_FDF)
  c_dataframe_FVF_d<-tips_only(c_dataframe_FVF)
  c_dataframe_MDH_d<-tips_only(c_dataframe_MDH)
  c_dataframe_MVH_d<-tips_only(c_dataframe_MVH)
  c_dataframe_FDH_d<-tips_only(c_dataframe_FDH)
  c_dataframe_FVH_d<-tips_only(c_dataframe_FVH)
  
  dD_dataframe_MDF_d<-tips_only(dD_dataframe_MDF)
  dD_dataframe_MVF_d<-tips_only(dD_dataframe_MVF)
  dD_dataframe_FDF_d<-tips_only(dD_dataframe_FDF)
  dD_dataframe_FVF_d<-tips_only(dD_dataframe_FVF)
  dD_dataframe_MDH_d<-tips_only(dD_dataframe_MDH)
  dD_dataframe_MVH_d<-tips_only(dD_dataframe_MVH)
  dD_dataframe_FDH_d<-tips_only(dD_dataframe_FDH)
  dD_dataframe_FVH_d<-tips_only(dD_dataframe_FVH)
  
  r_dataframe_MDF_d<-tips_only(r_dataframe_MDF)
  r_dataframe_MVF_d<-tips_only(r_dataframe_MVF)
  r_dataframe_FDF_d<-tips_only(r_dataframe_FDF)
  r_dataframe_FVF_d<-tips_only(r_dataframe_FVF)
  r_dataframe_MDH_d<-tips_only(r_dataframe_MDH)
  r_dataframe_MVH_d<-tips_only(r_dataframe_MVH)
  r_dataframe_FDH_d<-tips_only(r_dataframe_FDH)
  r_dataframe_FVH_d<-tips_only(r_dataframe_FVH)
  
  
  
  contributions_combined_tips=rbind(c_dataframe_MDF_d,c_dataframe_FDF_d,c_dataframe_MVF_d,c_dataframe_FVF_d,
                                    c_dataframe_MDH_d,c_dataframe_FDH_d,c_dataframe_MVH_d,c_dataframe_FVH_d)
  
  rates_combined_tips=rbind(r_dataframe_MDF_d,r_dataframe_FDF_d,r_dataframe_MVF_d,r_dataframe_FVF_d,
                            r_dataframe_MDH_d,r_dataframe_FDH_d,r_dataframe_MVH_d,r_dataframe_FVH_d)
  
  deltaD_combined_tips=rbind(dD_dataframe_MDF_d,dD_dataframe_FDF_d,dD_dataframe_MVF_d,dD_dataframe_FVF_d,
                             dD_dataframe_MDH_d,dD_dataframe_FDH_d,dD_dataframe_MVH_d,dD_dataframe_FVH_d)
  rownames(contributions_combined_tips)<-NULL
  rownames(rates_combined_tips)<-NULL
  rownames(deltaD_combined_tips)<-NULL
  ID_var<-str_remove(contributions_combined_tips$Descendant_node,'tip_')
  contributions_combined_tips1<-cbind(ID_var,contributions_combined_tips)
  ID_var<-str_remove(rates_combined_tips$Descendant_node,'tip_')
  rates_combined_tips1<-cbind(ID_var,rates_combined_tips)
  ID_var<-str_remove(deltaD_combined_tips$Descendant_node,'tip_')
  deltaD_combined_tips1<-cbind(ID_var,deltaD_combined_tips)
  
  if(make_csvs==TRUE){
    write.csv(contributions_combined_tips1,paste(prefix_str,'MOTHS_Sex_specific_contributions_dataset_perspecimenPCmedians_combined_tips.csv',sep='_'))
    write.csv(deltaD_combined_tips1,paste(prefix_str,'MOTHS_Sex_specific_deltaD_dataset_perspecimenPCmedians_combined_tips.csv',sep='_'))
    write.csv(rates_combined_tips1,paste(prefix_str,'MOTHS_Sex_specific_rates_dataset_perspecimenPCmedians_combined_tips.csv',sep='_'))
  }
  #####Setting up the input dataframes for the rates resampling dataframes----
  metadatafile<-data.frame(unique(cbind(PC$ID_var,PC$Family,PC$Categorical_var)))
  colnames(metadatafile)<-c('ID_var','Family','Categorical_var')
  # metadatafile<-read.csv('Imaged_all_species_moths_metadata_file.csv')
  rownames(metadatafile)<-metadatafile$ID_var
  rates<-merge(metadatafile,rates_combined_tips1,by='ID_var')
  rates_MDF=subset(rates,Sex=='male'& Side=='dorsal' & Wing=='FW')
  rates_MVF=subset(rates,Sex=='male'& Side=='ventral' & Wing=='FW')
  rates_FDF=subset(rates,Sex=='female'& Side=='dorsal' & Wing=='FW')
  rates_FVF=subset(rates,Sex=='female'& Side=='ventral' & Wing=='FW')
  rates_MDH=subset(rates,Sex=='male'& Side=='dorsal' & Wing=='HW')
  rates_MVH=subset(rates,Sex=='male'& Side=='ventral' & Wing=='HW')
  rates_FDH=subset(rates,Sex=='female'& Side=='dorsal' & Wing=='HW')
  rates_FVH=subset(rates,Sex=='female'& Side=='ventral' & Wing=='HW')
  
  #####Setting up the input dataframes for the IS-variance resampling----
  all_spec_metrics<-PC_all
  # species<-all_spec_metrics$ID_var
  # all_spec_metrics1<-cbind(species,all_spec_metrics)
  
  all_spec_metrics_MDF=subset(all_spec_metrics,Sex=='male'& Side=='dorsal' & Wing=='FW')
  all_spec_metrics_MVF=subset(all_spec_metrics,Sex=='male'& Side=='ventral' & Wing=='FW')
  all_spec_metrics_FDF=subset(all_spec_metrics,Sex=='female'& Side=='dorsal' & Wing=='FW')
  all_spec_metrics_FVF=subset(all_spec_metrics,Sex=='female'& Side=='ventral' & Wing=='FW')
  all_spec_metrics_MDH=subset(all_spec_metrics,Sex=='male'& Side=='dorsal' & Wing=='HW')
  all_spec_metrics_MVH=subset(all_spec_metrics,Sex=='male'& Side=='ventral' & Wing=='HW')
  all_spec_metrics_FDH=subset(all_spec_metrics,Sex=='female'& Side=='dorsal' & Wing=='HW')
  all_spec_metrics_FVH=subset(all_spec_metrics,Sex=='female'& Side=='ventral' & Wing=='HW')
  
  dorsal_IDs=intersect(intersect(all_spec_metrics_MDF$ID_var,all_spec_metrics_FDF$ID_var),
                       intersect(all_spec_metrics_MDH$ID_var,all_spec_metrics_FDH$ID_var))
  ventral_IDs=intersect(intersect(all_spec_metrics_MVF$ID_var,all_spec_metrics_FVF$ID_var),
                        intersect(all_spec_metrics_MVH$ID_var,all_spec_metrics_FVH$ID_var))
  final_ids=intersect(dorsal_IDs,ventral_IDs)
  all_spec_metrics_MDF1=subset(all_spec_metrics_MDF,ID_var %in% final_ids)
  all_spec_metrics_MVF1=subset(all_spec_metrics_MVF,ID_var %in% final_ids)
  all_spec_metrics_FDF1=subset(all_spec_metrics_FDF,ID_var %in% final_ids)
  all_spec_metrics_FVF1=subset(all_spec_metrics_FVF,ID_var %in% final_ids)
  
  all_spec_metrics_MDH1=subset(all_spec_metrics_MDH,ID_var %in% final_ids)
  all_spec_metrics_MVH1=subset(all_spec_metrics_MVH,ID_var %in% final_ids)
  all_spec_metrics_FDH1=subset(all_spec_metrics_FDH,ID_var %in% final_ids)
  all_spec_metrics_FVH1=subset(all_spec_metrics_FVH,ID_var %in% final_ids)
  
  # rm(eco.facet.05d.map.full.ext)
}
######Beginning the bootstrapping below----
metric<-PC_colnames[metric_num]
PC<-setDT(PC)
all_spec_metrics<-setDT(all_spec_metrics)
PC_all<-setDT(PC_all)

all_metrics_MDF2<-setDT(all_metrics_MDF1)
all_metrics_MVF2<-setDT(all_metrics_MVF1)
all_metrics_MDH2<-setDT(all_metrics_MDH1)
all_metrics_MVH2<-setDT(all_metrics_MVH1)
all_metrics_FDF2<-setDT(all_metrics_FDF1)
all_metrics_FVF2<-setDT(all_metrics_FVF1)
all_metrics_FDH2<-setDT(all_metrics_FDH1)
all_metrics_FVH2<-setDT(all_metrics_FVH1)

PC2<-rbind(all_metrics_MDF2,all_metrics_MVF2,all_metrics_MDH2,all_metrics_MVH2,
           all_metrics_FDF2,all_metrics_FVF2,all_metrics_FDH2,all_metrics_FVH2)

rates_MDF1<-setDT(rates_MDF)
rates_MVF1<-setDT(rates_MVF)
rates_FDF1<-setDT(rates_FDF)
rates_FVF1<-setDT(rates_FVF)
rates_MDH1<-setDT(rates_MDH)
rates_MVH1<-setDT(rates_MVH)
rates_FDH1<-setDT(rates_FDH)
rates_FVH1<-setDT(rates_FVH)
rates1<-rbind(rates_MDF1,rates_MVF1,rates_FDF1,rates_FVF1,
              rates_MDH1,rates_MVH1,rates_FDH1,rates_FVH1)
Phy_Side_Wing<-paste(rates1$ID_var,rates1$Side,rates1$Wing,sep='_')

rates1<-cbind(Phy_Side_Wing,rates1)
#####Now making the per-specimen level resampled IS-var dataframe
all_spec_metrics_MDF2<-setDT(all_spec_metrics_MDF1)
all_spec_metrics_MVF2<-setDT(all_spec_metrics_MVF1)
all_spec_metrics_FDF2<-setDT(all_spec_metrics_FDF1)
all_spec_metrics_FVF2<-setDT(all_spec_metrics_FVF1)
all_spec_metrics_MDH2<-setDT(all_spec_metrics_MDH1)
all_spec_metrics_MVH2<-setDT(all_spec_metrics_MVH1)
all_spec_metrics_FDH2<-setDT(all_spec_metrics_FDH1)
all_spec_metrics_FVH2<-setDT(all_spec_metrics_FVH1)
all_spec_metrics1<-rbind(all_spec_metrics_MDF2,all_spec_metrics_MVF2,all_spec_metrics_FDF2,all_spec_metrics_FVF2,
                         all_spec_metrics_MDH2,all_spec_metrics_MVH2,all_spec_metrics_FDH2,all_spec_metrics_FVH2)
Phy_Side_Wing<-paste(all_spec_metrics1$ID_var,all_spec_metrics1$Side,all_spec_metrics1$Wing,sep='_')
all_spec_metrics1<-cbind(Phy_Side_Wing,all_spec_metrics1)

###Editing the boot.pval function from the boot.pval package. For some reason it throws an error when the minimum potential
##p-value is set to 1e-16, but runs just fine if I change this to 1e-15. Probably something to do with rounding and precision
##In any case, the important point is that it is identical in every way but that small change in precision to the boot.pval
##function, but importantly, will no longer throw an error when calculating the p-values associated with BCa intervals.
boot.pval_fixed <- function(boot_res,
                            type = "perc",
                            theta_null = 0,
                            pval_precision = NULL,
                            ...)
{
  if(is.null(pval_precision)) { pval_precision = 1/boot_res$R }
  # Create a sequence of alphas:
  # alpha_seq <- seq(1e-15, 1-1e-15, pval_precision)
  # alpha_seq <- seq(1e-06, 1-1e-06, 1e-06)
  alpha_seq <- seq(pval_precision, 1-pval_precision, pval_precision)
  # alpha_seq<-sort(unique(c(2:10 %o% 10^-seq(1,14,0.001))))
  # Compute the 1-alpha confidence intervals, and extract
  # their bounds:
  ci <- suppressWarnings(boot::boot.ci(boot_res,
                                       conf = 1- alpha_seq,
                                       type = type,
                                       ...))
  bounds <- switch(type,
                   norm = ci$normal[,2:3],
                   basic = ci$basic[,4:5],
                   stud = ci$student[,4:5],
                   perc = ci$percent[,4:5],
                   bca = ci$bca[,4:5])
  
  # Find the smallest alpha such that theta_null is not contained in the 1-alpha
  # confidence interval:
  alpha <- alpha_seq[which.min(theta_null >= bounds[,1] & theta_null <= bounds[,2])]
  
  # Return the p-value:
  return(alpha)
}


########Sex contribution to interspecimen distance (AKA sexual distance dimorphism)----
mean_dist<-function(x){
  return(median(dist(x)))
}
scaled_all<-data.table(scale(all_spec_metrics1[,paste(PC_colnames),with=FALSE]))
all_spec_metrics2<-cbind(all_spec_metrics1[,!paste(PC_colnames),with=FALSE],scaled_all)

# Mean centering by species and sex
mean_cent_sp_sex <- all_spec_metrics2[, lapply(.SD, function(x) x - mean(x)), by = .(Phy_Side_Wing, Sex), .SDcols = PC_colnames]
dt_mean_cent_sp_sex <- cbind(all_spec_metrics1[, !paste(PC_colnames), with = FALSE], mean_cent_sp_sex)

##now mean centering just by species, leaving sex differences intact
mean_cent_sp <- all_spec_metrics2[, lapply(.SD, function(x) x - mean(x)), by = .(Phy_Side_Wing), .SDcols = PC_colnames]
dt_mean_cent_sp <- cbind(all_spec_metrics1[, !paste(PC_colnames), with = FALSE], mean_cent_sp)

##Calculating the mean interspecimen pairwise distance per species with sex differences removed
dt_mean_cent_sp_sex_dist<-setDT(ddply(dt_mean_cent_sp_sex,.(Phy_Side_Wing,ID_var,Family,Categorical_var,Side,Wing),colwise(mean_dist)))

##Calculating the mean  distance per species with sex differences retained
dt_mean_cent_sp_dist<-setDT(ddply(dt_mean_cent_sp,.(Phy_Side_Wing,ID_var,Family,Categorical_var,Side,Wing),colwise(mean_dist)))

# sanity checks, skip
# median(abs(PC_IS_med3_M[,paste(i),with=FALSE]-PC_IS_med3_F[,paste(i),with=FALSE])[[1]])
# 
# 
# dt_mean_cent_sp_dist_sub<-subset(dt_mean_cent_sp_dist,Wing=='HW'&Side=='ventral')
# dt_mean_cent_sp_sex_dist_sub<-subset(dt_mean_cent_sp_sex_dist,Wing=='HW'&Side=='ventral')
# 

# 
# median(dt_mean_cent_sp_sex_dist_sub$Area_Mask_cm2)
# median(dt_mean_cent_sp_dist_sub$Area_Mask_cm2)
# 
# sum(dt_mean_cent_sp_sex_dist_sub$Area_Mask_cm2 >= median(dt_mean_cent_sp_dist_sub$Area_Mask_cm2)) / length(dt_mean_cent_sp_sex_dist_sub$Area_Mask_cm2)
# 

###Sanity check shows that the median interspecimen distance is greater for all metrics when sex differences are retained
##this indicates sex differences.
##bootstrapping with random species resampling for the difference between these and testing whether these differences
##differ significantly from zero should test for sexual dimorphism no matter whether it varies in direction across species
# for(i in PC_colnames){
#   xx<-dt_mean_cent_sp_dist[,paste(i),with=FALSE]-dt_mean_cent_sp_sex_dist[,paste(i),with=FALSE]
#   print(median(xx[[1]]))
# }

median_Sex_IS_dist_boot<-function(indata1,indata2,variable,group,cores,metric,nboot){
  # print(metric)
  indata1<-setDT(indata1)
  indata2<-setDT(indata2)
  dat_out<-data.table()
  n_tests<-as.numeric(length(PC_colnames)*length(unique(indata1$Side))*length(unique(indata1$Wing)))
  Bonferroni_conf_level<-1-0.05/n_tests
  DF1<-subset(indata1,Wing=='FW'& Side=='dorsal'  & Categorical_var ==group)
  VF1<-subset(indata1,Wing=='FW'& Side=='ventral'  & Categorical_var ==group)
  DH1<-subset(indata1,Wing=='HW'& Side=='dorsal'  & Categorical_var ==group)
  VH1<-subset(indata1,Wing=='HW'& Side=='ventral'  & Categorical_var ==group)
  DF2<-subset(indata2,Wing=='FW'& Side=='dorsal'  & Categorical_var ==group)
  VF2<-subset(indata2,Wing=='FW'& Side=='ventral'  & Categorical_var ==group)
  DH2<-subset(indata2,Wing=='HW'& Side=='dorsal'  & Categorical_var ==group)
  VH2<-subset(indata2,Wing=='HW'& Side=='ventral'  & Categorical_var ==group)
  diff_df<-na.omit(DF1[,paste(metric),with=FALSE]-DF2[,paste(metric),with=FALSE])
  diff_dh<-na.omit(DH1[,paste(metric),with=FALSE]-DH2[,paste(metric),with=FALSE])
  diff_vf<-na.omit(VF1[,paste(metric),with=FALSE]-VF2[,paste(metric),with=FALSE])
  diff_vh<-na.omit(VH1[,paste(metric),with=FALSE]-VH2[,paste(metric),with=FALSE])
  median_func<-function(data,IDvec){
    med<-median(data[IDvec])
    return(med)
  }
  evaluate_boot_CI<-function(overall_value,indata1){
    if(is.na(overall_value)==TRUE){
      CI_lo<-NA
      CI_hi<-NA
      P_val<-NA
      CI_type=NA
    }else if(overall_value==0){
      boot_out<-boot(indata1,median_func,R=nboot,parallel='multicore',ncpus=cores)
      CI<-boot.ci(boot_out,conf=Bonferroni_conf_level,type='perc')
      CI_lo<-CI$percent[4]
      CI_hi<-CI$percent[5]
      P_val<-boot.pval_fixed(boot_out,theta_null=0,type='perc')
      CI_type='percentile'
    }else{
      boot_out<-boot(indata1,median_func,R=nboot,parallel='multicore',ncpus=cores)
      CI<-boot.ci(boot_out,conf=Bonferroni_conf_level,type='bca')
      CI_lo<-CI$bca[4]
      CI_hi<-CI$bca[5]
      P_val<-boot.pval_fixed(boot_out,theta_null=0,type='bca')
      CI_type='bca'
    } 
    return(list(CI_lo,CI_hi,P_val,CI_type))
  }
  overall_med_df<-median(diff_df[[1]])
  outs<-evaluate_boot_CI(overall_med_df,diff_df[[1]])
  CI_df_lo<-outs[1]
  CI_df_hi<-outs[2]
  P_val_df<-outs[3]
  CI_type_df<-outs[4]
  overall_med_dh<-median(diff_dh[[1]])
  outs<-evaluate_boot_CI(overall_med_dh,diff_dh[[1]])
  CI_dh_lo<-outs[1]
  CI_dh_hi<-outs[2]
  P_val_dh<-outs[3]
  CI_type_dh<-outs[4]
  overall_med_vf<-median(diff_vf[[1]])
  outs<-evaluate_boot_CI(overall_med_vf,diff_vf[[1]])
  CI_vf_lo<-outs[1]
  CI_vf_hi<-outs[2]
  P_val_vf<-outs[3]
  CI_type_vf<-outs[4]
  overall_med_vh<-median(diff_vh[[1]])
  outs<-evaluate_boot_CI(overall_med_vh,diff_vh[[1]])
  CI_vh_lo<-outs[1]
  CI_vh_hi<-outs[2]
  P_val_vh<-outs[3]
  CI_type_vh<-outs[4]
  overall_med<-c(overall_med_df,overall_med_dh,overall_med_vf,overall_med_vh)
  CI_lo<-c(CI_df_lo,CI_dh_lo,CI_vf_lo,CI_vh_lo)
  CI_hi<-c(CI_df_hi,CI_dh_hi,CI_vf_hi,CI_vh_hi)
  P_val<-c(P_val_df,P_val_dh,P_val_vf,P_val_vh)
  CI_type<-c(CI_type_df,CI_type_dh,CI_type_vf,CI_type_vh)
  Dataset<-c('DFW','DHW','VFW','VHW')
  dat<-data.table('Taxon'=Taxon_group,'Categorical_variable'=group,'Metrics'=metric,'Dataset'=Dataset,'median'=overall_med,'CI_lo'=CI_lo,'CI_hi'=CI_hi,'P_val'=P_val,'CI_type'=CI_type)
  names(dat)[names(dat) == 'median'] <- paste(variable,'med',sep='_')
  names(dat)[names(dat) == 'CI_hi'] <- paste(variable,'CI_hi',sep='_')
  names(dat)[names(dat) == 'CI_lo'] <- paste(variable,'CI_lo',sep='_')
  names(dat)[names(dat) == 'P_val'] <- paste(variable,'P_val',sep='_')
  names(dat)[names(dat) == 'CI_type'] <- paste(variable,'CI_type',sep='_')
  dat_out<-rbind(dat_out,dat)
  return(dat_out)
}


# MF_difference<-median_sexdiff_boot(PC3,'MF_difference',categorical_variable,n_cores,metric,boot_reps)
t1<-Sys.time()
Sex_IS_dist<-data.table()
for(i in unique(PC$Categorical_var)){
  boot_out<-median_Sex_IS_dist_boot(dt_mean_cent_sp_dist,dt_mean_cent_sp_sex_dist,'Sex_IS_dist_effect',i,n_cores,metric,boot_reps)
  Sex_IS_dist<-rbind(Sex_IS_dist,boot_out)
}
print('Finished Sex-effect IS dist dimorphism')
t2<-Sys.time()
t2-t1



#####Absolute rate proportion sex bias----

median_sexbias_boot<-function(indata,variable,group,cores,metric,nboot){
  # print(metric)
  indata<-setDT(indata)
  dat_out<-data.table()
  n_tests<-as.numeric(length(PC_colnames)*length(unique(indata$Side))*length(unique(indata$Wing)))
  Bonferroni_conf_level<-1-0.05/n_tests
  MDF<-subset(indata,Wing=='FW'& Side=='dorsal' & Sex=='male' & Categorical_var ==group)
  MVF<-subset(indata,Wing=='FW'& Side=='ventral' & Sex=='male' & Categorical_var ==group)
  MDH<-subset(indata,Wing=='HW'& Side=='dorsal' & Sex=='male' & Categorical_var ==group)
  MVH<-subset(indata,Wing=='HW'& Side=='ventral' & Sex=='male' & Categorical_var ==group)
  FDF<-subset(indata,Wing=='FW'& Side=='dorsal' & Sex=='female' & Categorical_var ==group)
  FVF<-subset(indata,Wing=='FW'& Side=='ventral' & Sex=='female' & Categorical_var ==group)
  FDH<-subset(indata,Wing=='HW'& Side=='dorsal' & Sex=='female' & Categorical_var ==group)
  FVH<-subset(indata,Wing=='HW'& Side=='ventral' & Sex=='female' & Categorical_var ==group)
  prop_df<-na.omit(abs(MDF[,paste(metric),with=FALSE])/(abs(MDF[,paste(metric),with=FALSE])+abs(FDF[,paste(metric),with=FALSE])))
  prop_dh<-na.omit(abs(MDH[,paste(metric),with=FALSE])/(abs(MDH[,paste(metric),with=FALSE])+abs(FDH[,paste(metric),with=FALSE])))
  prop_vf<-na.omit(abs(MVF[,paste(metric),with=FALSE])/(abs(MVF[,paste(metric),with=FALSE])+abs(FVF[,paste(metric),with=FALSE])))
  prop_vh<-na.omit(abs(MVH[,paste(metric),with=FALSE])/(abs(MVH[,paste(metric),with=FALSE])+abs(FVH[,paste(metric),with=FALSE])))
  median_func<-function(data,IDvec){
    med<-median(data[IDvec])
    return(med)
  }
  evaluate_boot_CI<-function(overall_value,indata){
    if(is.na(overall_value)==TRUE){
      CI_lo<-NA
      CI_hi<-NA
      P_val<-NA
      CI_type=NA
    }else if(overall_value==0){
      boot_out<-boot(indata,median_func,R=nboot,parallel='multicore',ncpus=cores)
      CI<-boot.ci(boot_out,conf=Bonferroni_conf_level,type='perc')
      CI_lo<-CI$percent[4]
      CI_hi<-CI$percent[5]
      P_val<-boot.pval_fixed(boot_out,theta_null=0.5,type='perc')
      CI_type='percentile'
    }else{
      boot_out<-boot(indata,median_func,R=nboot,parallel='multicore',ncpus=cores)
      CI<-boot.ci(boot_out,conf=Bonferroni_conf_level,type='bca')
      CI_lo<-CI$bca[4]
      CI_hi<-CI$bca[5]
      P_val<-boot.pval_fixed(boot_out,theta_null=0.5,type='bca')
      CI_type='bca'
    } 
    return(list(CI_lo,CI_hi,P_val,CI_type))
  }
  overall_med_df<-median(prop_df[[1]])
  outs<-evaluate_boot_CI(overall_med_df,prop_df[[1]])
  CI_df_lo<-outs[1]
  CI_df_hi<-outs[2]
  P_val_df<-outs[3]
  CI_type_df<-outs[4]
  overall_med_dh<-median(prop_dh[[1]])
  outs<-evaluate_boot_CI(overall_med_dh,prop_dh[[1]])
  CI_dh_lo<-outs[1]
  CI_dh_hi<-outs[2]
  P_val_dh<-outs[3]
  CI_type_dh<-outs[4]
  overall_med_vf<-median(prop_vf[[1]])
  outs<-evaluate_boot_CI(overall_med_vf,prop_vf[[1]])
  CI_vf_lo<-outs[1]
  CI_vf_hi<-outs[2]
  P_val_vf<-outs[3]
  CI_type_vf<-outs[4]
  overall_med_vh<-median(prop_vh[[1]])
  outs<-evaluate_boot_CI(overall_med_vh,prop_vh[[1]])
  CI_vh_lo<-outs[1]
  CI_vh_hi<-outs[2]
  P_val_vh<-outs[3]
  CI_type_vh<-outs[4]
  overall_med<-c(overall_med_df,overall_med_dh,overall_med_vf,overall_med_vh)
  CI_lo<-c(CI_df_lo,CI_dh_lo,CI_vf_lo,CI_vh_lo)
  CI_hi<-c(CI_df_hi,CI_dh_hi,CI_vf_hi,CI_vh_hi)
  P_val<-c(P_val_df,P_val_dh,P_val_vf,P_val_vh)
  CI_type<-c(CI_type_df,CI_type_dh,CI_type_vf,CI_type_vh)
  Dataset<-c('DFW','DHW','VFW','VHW')
  dat<-data.table('Taxon'=Taxon_group,'Categorical_variable'=group,'Metrics'=metric,'Dataset'=Dataset,'median'=overall_med,'CI_lo'=CI_lo,'CI_hi'=CI_hi,'P_val'=P_val,'CI_type'=CI_type)
  names(dat)[names(dat) == 'median'] <- paste(variable,'med',sep='_')
  names(dat)[names(dat) == 'CI_hi'] <- paste(variable,'CI_hi',sep='_')
  names(dat)[names(dat) == 'CI_lo'] <- paste(variable,'CI_lo',sep='_')
  names(dat)[names(dat) == 'P_val'] <- paste(variable,'P_val',sep='_')
  names(dat)[names(dat) == 'CI_type'] <- paste(variable,'CI_type',sep='_')
  dat_out<-rbind(dat_out,dat)
  return(dat_out)
}
abs_rate_sexbias<-data.table()
for(i in unique(PC$Categorical_var)){
  boot_out<-median_sexbias_boot(rates1,'Abs_rate_prop_sexbias',i,n_cores,metric,boot_reps)
  abs_rate_sexbias<-rbind(abs_rate_sexbias,boot_out)
}

print('Finished absolute rate sexbias')
#################MD proportion sexbias----

MD_sexbias_boot<-function(indata,variable,group,cores,metric,nboot){
  dat_out<-data.table()
  MDF<-subset(indata,Wing=='FW'& Side=='dorsal' & Sex=='male' & Categorical_var ==group)
  MVF<-subset(indata,Wing=='FW'& Side=='ventral' & Sex=='male' & Categorical_var ==group)
  MDH<-subset(indata,Wing=='HW'& Side=='dorsal' & Sex=='male' & Categorical_var ==group)
  MVH<-subset(indata,Wing=='HW'& Side=='ventral' & Sex=='male' & Categorical_var ==group)
  FDF<-subset(indata,Wing=='FW'& Side=='dorsal' & Sex=='female' & Categorical_var ==group)
  FVF<-subset(indata,Wing=='FW'& Side=='ventral' & Sex=='female' & Categorical_var ==group)
  FDH<-subset(indata,Wing=='HW'& Side=='dorsal' & Sex=='female' & Categorical_var ==group)
  FVH<-subset(indata,Wing=='HW'& Side=='ventral' & Sex=='female' & Categorical_var ==group)
  n_tests<-as.numeric(length(PC_colnames)*length(unique(indata$Side))*length(unique(indata$Wing)))
  Bonferroni_conf_level<-1-0.05/n_tests
  disp_func<-function(data,IDvec){
    dats<-data[IDvec,]
    tm<-dats[,1]
    tf<-dats[,2]
    trial<-c(tm,tf)
    lenm<-as.numeric(length(tm))
    lenf<-as.numeric(length(tf))
    # try(disp<-dispRity.per.group(bootstraps=10,data=matrix(trial),list("Male"=c(1:lenm),"Female"=c((lenm+1):(lenm+lenf))),
    # metric=c(variances)),silent=FALSE)
    try(disp<-dispRity(boot.matrix(custom.subsets(matrix(trial), list("Male"=c(1:lenm),"Female"=c((lenm+1):(lenm+lenf)))),
                                   bootstraps = 1),
                       metric = c(variances)))
    male_results<-disp$disparity$Male$elements
    female_results<-disp$disparity$Female$elements
    prop_df<-male_results/(male_results+female_results)
    return(prop_df)
  }
  evaluate_boot_CI<-function(overall_value,indata){
    if(is.na(overall_value)==TRUE){
      CI_lo<-NA
      CI_hi<-NA
      P_val<-NA
      CI_type=NA
    }else if(overall_value==0){
      boot_out<-boot(indata,disp_func,R=nboot,parallel='multicore',ncpus=cores)
      CI<-boot.ci(boot_out,conf=Bonferroni_conf_level,type='perc')
      CI_lo<-CI$percent[4]
      CI_hi<-CI$percent[5]
      P_val<-boot.pval_fixed(boot_out,theta_null=0.5,type='perc')
      CI_type='percentile'
    }else{
      boot_out<-boot(indata,disp_func,R=nboot,parallel='multicore',ncpus=cores)
      CI<-boot.ci(boot_out,conf=Bonferroni_conf_level,type='bca')
      CI_lo<-CI$bca[4]
      CI_hi<-CI$bca[5]
      P_val<-boot.pval_fixed(boot_out,theta_null=0.5,type='bca')
      CI_type='bca'
    } 
    return(list(CI_lo,CI_hi,P_val,CI_type))
  }
  tt<-cbind(MDF[[paste(metric)]],FDF[[paste(metric)]])
  overall_prop_df<-disp_func(tt,(1:length(tt[,1])))
  outs<-evaluate_boot_CI(overall_prop_df,tt)
  CI_df_lo<-outs[1]
  CI_df_hi<-outs[2]
  P_val_df<-outs[3]
  CI_type_df<-outs[4]
  tt<-cbind(MDH[[paste(metric)]],FDH[[paste(metric)]])
  overall_prop_dh<-disp_func(tt,(1:length(tt[,1])))
  outs<-evaluate_boot_CI(overall_prop_dh,tt)
  CI_dh_lo<-outs[1]
  CI_dh_hi<-outs[2]
  P_val_dh<-outs[3]
  CI_type_dh<-outs[4]
  tt<-cbind(MVF[[paste(metric)]],FVF[[paste(metric)]])
  overall_prop_vf<-disp_func(tt,(1:length(tt[,1])))
  outs<-evaluate_boot_CI(overall_prop_vf,tt)
  CI_vf_lo<-outs[1]
  CI_vf_hi<-outs[2]
  P_val_vf<-outs[3]
  CI_type_vf<-outs[4]
  tt<-cbind(MVH[[paste(metric)]],FVH[[paste(metric)]])
  overall_prop_vh<-disp_func(tt,(1:length(tt[,1])))
  outs<-evaluate_boot_CI(overall_prop_vh,tt)
  CI_vh_lo<-outs[1]
  CI_vh_hi<-outs[2]
  P_val_vh<-outs[3]
  CI_type_vh<-outs[4]
  overall_prop<-c(overall_prop_df,overall_prop_dh,overall_prop_vf,overall_prop_vh)
  CI_lo<-c(CI_df_lo,CI_dh_lo,CI_vf_lo,CI_vh_lo)
  CI_hi<-c(CI_df_hi,CI_dh_hi,CI_vf_hi,CI_vh_hi)
  P_val<-c(P_val_df,P_val_dh,P_val_vf,P_val_vh)
  CI_type<-c(CI_type_df,CI_type_dh,CI_type_vf,CI_type_vh)
  Dataset<-c('DFW','DHW','VFW','VHW')
  dat<-data.table('Taxon'=Taxon_group,'Categorical_variable'=group,'Metrics'=metric,'Dataset'=Dataset,'prop'=overall_prop,'CI_lo'=CI_lo,'CI_hi'=CI_hi,'P_val'=P_val,'CI_type'=CI_type)
  names(dat)[names(dat) == 'prop'] <- paste(variable,'prop',sep='_')
  names(dat)[names(dat) == 'CI_hi'] <- paste(variable,'CI_hi',sep='_')
  names(dat)[names(dat) == 'CI_lo'] <- paste(variable,'CI_lo',sep='_')
  names(dat)[names(dat) == 'P_val'] <- paste(variable,'P_val',sep='_')
  names(dat)[names(dat) == 'CI_type'] <- paste(variable,'CI_type',sep='_')
  dat_out<-rbind(dat_out,dat)
  return(dat_out)
}
# 
# MD_prop_sexbias<-MD_sexbias_boot(PC2,'Extant_var_MD_prop_sexbias',categorical_variable,n_cores,metric,boot_reps)

t1<-Sys.time()
MD_prop_sexbias<-data.table()
for(i in unique(PC$Categorical_var)){
  boot_out<-MD_sexbias_boot(PC2,'Extant_var_MD_prop_sexbias',i,n_cores,metric,boot_reps)
  MD_prop_sexbias<-rbind(MD_prop_sexbias,boot_out)
}
print('Finished MD sexbias')
t2<-Sys.time()
t2-t1


########M-F sexual dimorphism----
scaled<-data.table(scale(PC2[,paste(PC_colnames),with=FALSE]))
PC3<-cbind(PC2[,!paste(PC_colnames),with=FALSE],scaled)
median_sexdiff_boot<-function(indata,variable,group,cores,metric,nboot){
  # print(metric)
  indata1<-setDT(indata)
  dat_out<-data.table()
  n_tests<-as.numeric(length(PC_colnames)*length(unique(indata1$Side))*length(unique(indata1$Wing)))
  Bonferroni_conf_level<-1-0.05/n_tests
  MDF<-subset(indata1,Wing=='FW'& Side=='dorsal' & Sex=='male' & Categorical_var ==group)
  MVF<-subset(indata1,Wing=='FW'& Side=='ventral' & Sex=='male' & Categorical_var ==group)
  MDH<-subset(indata1,Wing=='HW'& Side=='dorsal' & Sex=='male' & Categorical_var ==group)
  MVH<-subset(indata1,Wing=='HW'& Side=='ventral' & Sex=='male' & Categorical_var ==group)
  FDF<-subset(indata1,Wing=='FW'& Side=='dorsal' & Sex=='female' & Categorical_var ==group)
  FVF<-subset(indata1,Wing=='FW'& Side=='ventral' & Sex=='female' & Categorical_var ==group)
  FDH<-subset(indata1,Wing=='HW'& Side=='dorsal' & Sex=='female' & Categorical_var ==group)
  FVH<-subset(indata1,Wing=='HW'& Side=='ventral' & Sex=='female' & Categorical_var ==group)
  diff_df<-na.omit(MDF[,paste(metric),with=FALSE]-FDF[,paste(metric),with=FALSE])
  diff_dh<-na.omit(MDH[,paste(metric),with=FALSE]-FDH[,paste(metric),with=FALSE])
  diff_vf<-na.omit(MVF[,paste(metric),with=FALSE]-FVF[,paste(metric),with=FALSE])
  diff_vh<-na.omit(MVH[,paste(metric),with=FALSE]-FVH[,paste(metric),with=FALSE])
  median_func<-function(data,IDvec){
    med<-median(data[IDvec])
    return(med)
  }
  evaluate_boot_CI<-function(overall_value,indata1){
    if(is.na(overall_value)==TRUE){
      CI_lo<-NA
      CI_hi<-NA
      P_val<-NA
      CI_type=NA
    }else if(overall_value==0){
      boot_out<-boot(indata1,median_func,R=nboot,parallel='multicore',ncpus=cores)
      CI<-boot.ci(boot_out,conf=Bonferroni_conf_level,type='perc')
      CI_lo<-CI$percent[4]
      CI_hi<-CI$percent[5]
      P_val<-boot.pval_fixed(boot_out,theta_null=0,type='perc')
      CI_type='percentile'
    }else{
      boot_out<-boot(indata1,median_func,R=nboot,parallel='multicore',ncpus=cores)
      CI<-boot.ci(boot_out,conf=Bonferroni_conf_level,type='bca')
      CI_lo<-CI$bca[4]
      CI_hi<-CI$bca[5]
      P_val<-boot.pval_fixed(boot_out,theta_null=0,type='bca')
      CI_type='bca'
    } 
    return(list(CI_lo,CI_hi,P_val,CI_type))
  }
  overall_med_df<-median(diff_df[[1]])
  outs<-evaluate_boot_CI(overall_med_df,diff_df[[1]])
  CI_df_lo<-outs[1]
  CI_df_hi<-outs[2]
  P_val_df<-outs[3]
  CI_type_df<-outs[4]
  overall_med_dh<-median(diff_dh[[1]])
  outs<-evaluate_boot_CI(overall_med_dh,diff_dh[[1]])
  CI_dh_lo<-outs[1]
  CI_dh_hi<-outs[2]
  P_val_dh<-outs[3]
  CI_type_dh<-outs[4]
  overall_med_vf<-median(diff_vf[[1]])
  outs<-evaluate_boot_CI(overall_med_vf,diff_vf[[1]])
  CI_vf_lo<-outs[1]
  CI_vf_hi<-outs[2]
  P_val_vf<-outs[3]
  CI_type_vf<-outs[4]
  overall_med_vh<-median(diff_vh[[1]])
  outs<-evaluate_boot_CI(overall_med_vh,diff_vh[[1]])
  CI_vh_lo<-outs[1]
  CI_vh_hi<-outs[2]
  P_val_vh<-outs[3]
  CI_type_vh<-outs[4]
  overall_med<-c(overall_med_df,overall_med_dh,overall_med_vf,overall_med_vh)
  CI_lo<-c(CI_df_lo,CI_dh_lo,CI_vf_lo,CI_vh_lo)
  CI_hi<-c(CI_df_hi,CI_dh_hi,CI_vf_hi,CI_vh_hi)
  P_val<-c(P_val_df,P_val_dh,P_val_vf,P_val_vh)
  CI_type<-c(CI_type_df,CI_type_dh,CI_type_vf,CI_type_vh)
  Dataset<-c('DFW','DHW','VFW','VHW')
  dat<-data.table('Taxon'=Taxon_group,'Categorical_variable'=group,'Metrics'=metric,'Dataset'=Dataset,'median'=overall_med,'CI_lo'=CI_lo,'CI_hi'=CI_hi,'P_val'=P_val,'CI_type'=CI_type)
  names(dat)[names(dat) == 'median'] <- paste(variable,'med',sep='_')
  names(dat)[names(dat) == 'CI_hi'] <- paste(variable,'CI_hi',sep='_')
  names(dat)[names(dat) == 'CI_lo'] <- paste(variable,'CI_lo',sep='_')
  names(dat)[names(dat) == 'P_val'] <- paste(variable,'P_val',sep='_')
  names(dat)[names(dat) == 'CI_type'] <- paste(variable,'CI_type',sep='_')
  dat_out<-rbind(dat_out,dat)
  return(dat_out)
}


# MF_difference<-median_sexdiff_boot(PC3,'MF_difference',categorical_variable,n_cores,metric,boot_reps)

MF_difference<-data.table()
for(i in unique(PC$Categorical_var)){
  boot_out<-median_sexdiff_boot(PC3,'MF_difference',i,n_cores,metric,boot_reps)
  MF_difference<-rbind(MF_difference,boot_out)
}
print('Finished MF diff dimorphism')

######M-F absolute rate differences ----
scaled<-data.table(scale(rates1[,paste(PC_colnames),with=FALSE]))
rates_scaled<-cbind(rates1[,!paste(PC_colnames),with=FALSE],scaled)

# MF_abs_rates<-median_sexdiff_boot(rates_scaled,'MF_abs_rates',categorical_variable,n_cores,metric,boot_reps)

MF_abs_rates<-data.table()
for(i in unique(PC$Categorical_var)){
  boot_out<-median_sexdiff_boot(rates_scaled,'MF_abs_rates',i,n_cores,metric,boot_reps)
  MF_abs_rates<-rbind(MF_abs_rates,boot_out)
}

print('Finished MF abs rates difference')


#############M-F extant var difference----
MD_sexdiff_boot<-function(indata,variable,group,cores,metric,nboot){
  dat_out<-data.table()
  MDF<-subset(indata,Wing=='FW'& Side=='dorsal' & Sex=='male' & Categorical_var ==group)
  MVF<-subset(indata,Wing=='FW'& Side=='ventral' & Sex=='male' & Categorical_var ==group)
  MDH<-subset(indata,Wing=='HW'& Side=='dorsal' & Sex=='male' & Categorical_var ==group)
  MVH<-subset(indata,Wing=='HW'& Side=='ventral' & Sex=='male' & Categorical_var ==group)
  FDF<-subset(indata,Wing=='FW'& Side=='dorsal' & Sex=='female' & Categorical_var ==group)
  FVF<-subset(indata,Wing=='FW'& Side=='ventral' & Sex=='female' & Categorical_var ==group)
  FDH<-subset(indata,Wing=='HW'& Side=='dorsal' & Sex=='female' & Categorical_var ==group)
  FVH<-subset(indata,Wing=='HW'& Side=='ventral' & Sex=='female' & Categorical_var ==group)
  n_tests<-as.numeric(length(PC_colnames)*length(unique(indata$Side))*length(unique(indata$Wing)))
  Bonferroni_conf_level<-1-0.05/n_tests
  disp_func<-function(data,IDvec){
    dats<-data[IDvec,]
    tm<-dats[,1]
    tf<-dats[,2]
    trial<-c(tm,tf)
    lenm<-as.numeric(length(tm))
    lenf<-as.numeric(length(tf))
    # try(disp<-dispRity.per.group(data=matrix(trial),list("Male"=c(1:lenm),"Female"=c((lenm+1):(lenm+lenf))),
    #                              metric=c(variances)),silent=FALSE)
    try(disp<-dispRity(boot.matrix(custom.subsets(matrix(trial), list("Male"=c(1:lenm),"Female"=c((lenm+1):(lenm+lenf)))),
                                   bootstraps = 1),
                       metric = c(variances)))
    male_results<-disp$disparity$Male$elements
    female_results<-disp$disparity$Female$elements
    diff<-male_results-female_results
    return(diff)
  }
  evaluate_boot_CI<-function(overall_value,indata){
    if(is.na(overall_value)==TRUE){
      CI_lo<-NA
      CI_hi<-NA
      P_val<-NA
      CI_type=NA
    }else if(overall_value==0){
      boot_out<-boot(indata,disp_func,R=nboot,parallel='multicore',ncpus=cores)
      CI<-boot.ci(boot_out,conf=Bonferroni_conf_level,type='perc')
      CI_lo<-CI$percent[4]
      CI_hi<-CI$percent[5]
      P_val<-boot.pval_fixed(boot_out,theta_null=0,type='perc')
      CI_type='percentile'
    }else{
      boot_out<-boot(indata,disp_func,R=nboot,parallel='multicore',ncpus=cores)
      CI<-boot.ci(boot_out,conf=Bonferroni_conf_level,type='bca')
      CI_lo<-CI$bca[4]
      CI_hi<-CI$bca[5]
      P_val<-boot.pval_fixed(boot_out,theta_null=0,type='bca')
      CI_type='bca'
    } 
    return(list(CI_lo,CI_hi,P_val,CI_type))
  }
  tt<-cbind(MDF[[paste(metric)]],FDF[[paste(metric)]])
  overall_diff_df<-disp_func(tt,(1:length(tt[,1])))
  outs<-evaluate_boot_CI(overall_diff_df,tt)
  CI_df_lo<-outs[1]
  CI_df_hi<-outs[2]
  P_val_df<-outs[3]
  CI_type_df<-outs[4]
  tt<-cbind(MDH[[paste(metric)]],FDH[[paste(metric)]])
  overall_diff_dh<-disp_func(tt,(1:length(tt[,1])))
  outs<-evaluate_boot_CI(overall_diff_dh,tt)
  CI_dh_lo<-outs[1]
  CI_dh_hi<-outs[2]
  P_val_dh<-outs[3]
  CI_type_dh<-outs[4]
  tt<-cbind(MVF[[paste(metric)]],FVF[[paste(metric)]])
  overall_diff_vf<-disp_func(tt,(1:length(tt[,1])))
  outs<-evaluate_boot_CI(overall_diff_vf,tt)
  CI_vf_lo<-outs[1]
  CI_vf_hi<-outs[2]
  P_val_vf<-outs[3]
  CI_type_vf<-outs[4]
  tt<-cbind(MVH[[paste(metric)]],FVH[[paste(metric)]])
  overall_diff_vh<-disp_func(tt,(1:length(tt[,1])))
  outs<-evaluate_boot_CI(overall_diff_vh,tt)
  CI_vh_lo<-outs[1]
  CI_vh_hi<-outs[2]
  P_val_vh<-outs[3]
  CI_type_vh<-outs[4]
  overall_diff<-c(overall_diff_df,overall_diff_dh,overall_diff_vf,overall_diff_vh)
  CI_lo<-c(CI_df_lo,CI_dh_lo,CI_vf_lo,CI_vh_lo)
  CI_hi<-c(CI_df_hi,CI_dh_hi,CI_vf_hi,CI_vh_hi)
  P_val<-c(P_val_df,P_val_dh,P_val_vf,P_val_vh)
  CI_type<-c(CI_type_df,CI_type_dh,CI_type_vf,CI_type_vh)
  Dataset<-c('DFW','DHW','VFW','VHW')
  dat<-data.table('Taxon'=Taxon_group,'Categorical_variable'=group,'Metrics'=metric,'Dataset'=Dataset,'med'=overall_diff,'CI_lo'=CI_lo,'CI_hi'=CI_hi,'P_val'=P_val,'CI_type'=CI_type)
  names(dat)[names(dat) == 'med'] <- paste(variable,'med',sep='_')
  names(dat)[names(dat) == 'CI_hi'] <- paste(variable,'CI_hi',sep='_')
  names(dat)[names(dat) == 'CI_lo'] <- paste(variable,'CI_lo',sep='_')
  names(dat)[names(dat) == 'P_val'] <- paste(variable,'P_val',sep='_')
  names(dat)[names(dat) == 'CI_type'] <- paste(variable,'CI_type',sep='_')
  dat_out<-rbind(dat_out,dat)
  return(dat_out)
}

# MF_extant_var<-MD_sexdiff_boot(PC2,'MF_extant_var',categorical_variable,n_cores,metric,boot_reps)
MF_extant_var<-data.table()
for(i in unique(PC$Categorical_var)){
  boot_out<-MD_sexdiff_boot(PC2,'MF_extant_var',i,n_cores,metric,boot_reps)
  MF_extant_var<-rbind(MF_extant_var,boot_out)
}
print('Finished MF extant var difference')

##############SSQ metrics----
dimorph_ssq_boot<-function(indata,group,cores,metric,nboot){
  dat_out<-data.table()
  MDF<-subset(indata,Wing=='FW'& Side=='dorsal' & Sex=='male' & Categorical_var ==group)
  MVF<-subset(indata,Wing=='FW'& Side=='ventral' & Sex=='male' & Categorical_var ==group)
  MDH<-subset(indata,Wing=='HW'& Side=='dorsal' & Sex=='male' & Categorical_var ==group)
  MVH<-subset(indata,Wing=='HW'& Side=='ventral' & Sex=='male' & Categorical_var ==group)
  FDF<-subset(indata,Wing=='FW'& Side=='dorsal' & Sex=='female' & Categorical_var ==group)
  FVF<-subset(indata,Wing=='FW'& Side=='ventral' & Sex=='female' & Categorical_var ==group)
  FDH<-subset(indata,Wing=='HW'& Side=='dorsal' & Sex=='female' & Categorical_var ==group)
  FVH<-subset(indata,Wing=='HW'& Side=='ventral' & Sex=='female' & Categorical_var ==group)
  n_tests<-as.numeric(length(PC_colnames)*length(unique(indata$Side))*length(unique(indata$Wing)))
  Bonferroni_conf_level<-1-0.05/n_tests
  SSQ_func<-function(data,IDvec){
    dats<-data[IDvec,]
    tm<-dats[,1]
    tf<-dats[,2]
    trial<-c(tm,tf)
    MF<-abs(tm-tf)
    grand_mean<-mean(trial)
    mean_m<-mean(tm)
    mean_f<-mean(tf)
    ss_tm<-sum(((mean_m-grand_mean)^2)*length(tm))
    ss_tf<-sum(((mean_f-grand_mean)^2)*length(tf))
    ss_T<-ss_tm+ss_tf
    ss_tot<-sum((trial-grand_mean)^2)
    ssT_to_tot<-ss_T/ss_tot
    sec_moment_distance<-(sum((MF-0)^2))/length(MF)
    n<-as.numeric(length(trial))
    var_tot<-ss_tot/n
    return(c(ssT_to_tot,var_tot,sec_moment_distance))
  }
  evaluate_boot_CI<-function(overall_value,indata){
    dt_out<-data.table()
    for(i in 1:length(overall_value)){
      if(is.na(overall_value[i])==TRUE){
        CI_lo<-NA
        CI_hi<-NA
        CI_type=NA
      }else if(overall_value[i]==0){
        boot_out<-boot(indata,SSQ_func,R=nboot,parallel='multicore',ncpus=cores)
        CI<-boot.ci(boot_out,conf=Bonferroni_conf_level,type='perc',index=i)
        CI_lo<-CI$percent[4]
        CI_hi<-CI$percent[5]
        CI_type='percentile'
      }else{
        boot_out<-boot(indata,SSQ_func,R=nboot,parallel='multicore',ncpus=cores)
        CI<-boot.ci(boot_out,conf=Bonferroni_conf_level,type='bca',index=i)
        CI_lo<-CI$bca[4]
        CI_hi<-CI$bca[5]
        CI_type='bca'
      } 
      dt<-data.table('CI_lo'=CI_lo,'CI_hi'=CI_hi,'CI_type'=CI_type)
      dt_out<-rbind(dt_out,dt)
    }
    dt_out<-cbind(c('ssT_to_tot','var_tot','sec_moment_distance'),dt_out)
    colnames(dt_out)<-c('metric','CI_lo','CI_hi','CI_type')
    return(dt_out)
  }
  tt<-cbind(MDF[[paste(metric)]],FDF[[paste(metric)]])
  overall_diff_df<-SSQ_func(tt,(1:length(tt[,1])))
  outs<-evaluate_boot_CI(overall_diff_df,tt)
  CI_df_lo<-outs$CI_lo
  CI_df_hi<-outs$CI_hi
  CI_type_df<-outs$CI_type
  tt<-cbind(MDH[[paste(metric)]],FDH[[paste(metric)]])
  overall_diff_dh<-SSQ_func(tt,(1:length(tt[,1])))
  outs<-evaluate_boot_CI(overall_diff_dh,tt)
  CI_dh_lo<-outs$CI_lo
  CI_dh_hi<-outs$CI_hi
  CI_type_dh<-outs$CI_type
  tt<-cbind(MVF[[paste(metric)]],FVF[[paste(metric)]])
  overall_diff_vf<-SSQ_func(tt,(1:length(tt[,1])))
  outs<-evaluate_boot_CI(overall_diff_vf,tt)
  CI_vf_lo<-outs$CI_lo
  CI_vf_hi<-outs$CI_hi
  CI_type_vf<-outs$CI_type
  tt<-cbind(MVH[[paste(metric)]],FVH[[paste(metric)]])
  overall_diff_vh<-SSQ_func(tt,(1:length(tt[,1])))
  outs<-evaluate_boot_CI(overall_diff_vh,tt)
  CI_vh_lo<-outs$CI_lo
  CI_vh_hi<-outs$CI_hi
  CI_type_vh<-outs$CI_type
  final_out<-data.table('overall_diff'=c(overall_diff_df,overall_diff_dh,overall_diff_vf,overall_diff_vh),
                        'CI_lo'=c(CI_df_lo,CI_dh_lo,CI_vf_lo,CI_vh_lo), 'CI_hi'=c(CI_df_hi,CI_dh_hi,CI_vf_hi,CI_vh_hi),
                        'CI_type'=c(CI_type_df,CI_type_dh,CI_type_vf,CI_type_vh))
  final_out$Dataset<-c(rep('DFW',3),rep('DHW',3),rep('VFW',3),rep('VHW',3))
  final_out$value<-rep(c('ssT_to_tot','var_tot','sec_mom_dist'),4)
  dat<-data.table('Taxon'=Taxon_group,'Categorical_variable'=group,'Metrics'=metric,'Dataset'=c('DFW','DHW','VFW','VHW'),
                  'ssT_to_tot_med'=final_out[value=='ssT_to_tot']$overall_diff,'ssT_to_tot_CI_lo'=final_out[value=='ssT_to_tot']$CI_lo,'ssT_to_tot_CI_hi'=final_out[value=='ssT_to_tot']$CI_hi,'ssT_to_tot_CI_type'=final_out[value=='ssT_to_tot']$CI_type,
                  'tot_var_med'=final_out[value=='var_tot']$overall_diff,'tot_var_CI_lo'=final_out[value=='var_tot']$CI_lo,'tot_var_CI_hi'=final_out[value=='var_tot']$CI_hi,'tot_var_CI_type'=final_out[value=='var_tot']$CI_type,
                  'sec_mom_var_med'=final_out[value=='sec_mom_dist']$overall_diff,'sec_mom_var_CI_lo'=final_out[value=='sec_mom_dist']$CI_lo,'sec_mom_var_CI_hi'=final_out[value=='sec_mom_dist']$CI_hi,'sec_mom_var_CI_type'=final_out[value=='sec_mom_dist']$CI_type)
  dat_out<-rbind(dat_out,dat)
  return(dat_out)
}

# SSQ_frame<-dimorph_ssq_boot(PC3,categorical_variable,n_cores,metric,boot_reps)
SSQ_frame<-data.table()
for(i in unique(PC$Categorical_var)){
  boot_out<-dimorph_ssq_boot(PC3,i,n_cores,metric,boot_reps)
  SSQ_frame<-rbind(SSQ_frame,boot_out)
}
print('Finished SSQ values')

final_frame<-join_all(list(MF_difference,abs_rate_sexbias,MD_prop_sexbias,MF_abs_rates,MF_extant_var,SSQ_frame,Sex_IS_dist))

boot_reps2=format(boot_reps, scientific = FALSE)

fwrite(final_frame,paste('/path/to/results/',metric,prefix_str,categorical_variable,'Bootpkg_dimorph_rates_MD_sexbias_diff_ssq_bonferroniCI_rawpval',specimen_threshold,'mf',boot_reps2,'BS.csv',sep='_'))


end_time<-Sys.time()
end_time-start_time
