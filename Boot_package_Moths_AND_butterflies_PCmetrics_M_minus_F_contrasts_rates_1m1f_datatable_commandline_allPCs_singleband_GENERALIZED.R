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


# Taxon_group<-'Moths_and_butterflies'
# categorical_variable='Moth_butterfly_nocturnality'
# rdata<-'NEW_PC_values_bootstrapping_overall_data_Mothbutterfly_MothButtnocturnality_All_1m1f_noHemaris.RData'
# specimen_threshold=1
# n_cores=20
# boot_reps=100000
##First read in the arguments listed at the command line
args=(commandArgs(TRUE))
##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
  print("No arguments supplied.")
  ##supply default values
  boot_reps<-10000
  n_cores<-1
  Taxon_group<-'Moths_and_butterflies'
  categorical_variable='Moth_butterfly_nocturnality'
  # rdata<-'NEW_PC_individual_values_bootstrapping_overall_data_NocturnalONLY_moths_All_1m1f_redone.RData'
  rdata<-'NEW_PC_values_bootstrapping_overall_data_Mothbutterfly_MothButtnocturnality_All_1m1f_noHemaris.RData'
  specimen_threshold=1
  rep=1
  wband='UV'
}else{
  for(i in 1:length(args)){
    # print(args)
    eval(parse(text=args[[i]]))
  }
}

if(rdata!='No_Data'){
  load(rdata)
  Final_output_plotframe<-data.frame()
  Final_output_dimorphframe<-data.frame()
  #reparsing command line arguments in case loading the environment overwrote any of them
  args=(commandArgs(TRUE))
  print(args)
  for(i in 1:length(args)){
    print(args)
    eval(parse(text=args[[i]]))
  }
}else{
  #########################################################
  ##Loading in the PC metrics ----
  ##Currently works for both butterflies and moths with "All" as the categorical variable
  ##insert the taxon group as the first argument supplied to the sbatch command
  #choices are 'Butterflies' 'Moths'
  if(Taxon_group=='Moths'){
    if(categorical_variable=='All'){
      # PC_input<-PC_input[PC_input$Diel_behavior=='Nocturnal',]
      PC_input<-read.csv(paste('Nocturnal_Moths_1m1f_sensitivitycut_Final_PCA_scoreframe_sensitivity_curated_1MF_AllPCs_',wband,'.csv',sep=''))
      PC_col_start=11
      PC_col_end=length(colnames(PC_input))
      prefix_str='Moth_nocturnal_BS_FINAL_allPCs'
    }    
    IDvar<-'Moth_ID'
  }
  if(Taxon_group=='Butterflies'){
    PC_input<-read.csv(paste('Butterflies_1m1f_sensitivitycut_Final_PCA_scoreframe_sensitivity_curated_1MF_AllPCs_',wband,'.csv',sep=''))
    prefix_str='Butterflies_BS_FINAL_allPCs'
    IDvar<-'Phy_ID'
    PC_col_start=10
    PC_col_end=length(colnames(PC_input))
    # categorical_variable='All'
  }
  #############SET THESE BEFORE RUNNING
  # categorical_variable='Diel_behavior'
  # categorical_variable='Family'
  # categorical_variable='All'
  # categorical_variable='Moth_butterfly'
  make_csvs=FALSE
  final_csvs=FALSE
  # IS_var_thresh=3
  
  ###Add unnecessary variables to drop from the original dataframe here
  if(Taxon_group=='Moths'){
    # drop=c('UID','UID.1','Unnamed..0','Family_Side_Wing','Sex_Side_Wing','Diel_behavior_source',
    #        'Genus','Species','X1cm_scale_pixels_UV_940')
    # drop=c('UID','UID.1','Family_Side_Wing','Sex_Side_Wing','Diel_behavior_reference',
    #        'Genus','Species','X1cm_scale_pixels_UV_940','Moth_ID.1')
    drop=c(NULL)
  }
  
  if(Taxon_group=='Butterflies'){
    # drop=c('UID','Family_Side_Wing','Sex_Side_Wing','Family_Side_Wing','Sex_Side_Wing',
    #        'Genus_valid','Identity_of_this_specimen','Phy_oldID','Swap_or_not','Loc_ID',
    #        'X1cm_scale_pixels_UV_940')
    drop=c(NULL)
    # names(PC_input)[names(PC_input) == 'Specimen_Specimen'] <- 'Specimen'
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
  
  unique(PC[,"ID_var"])
  unique(all_metrics_FDF1[,"ID_var"])
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
  results_MVF<-phylopars_multi_anc(all_metrics_MVF1,phylopar_tree_copy,'MVF',PC_start+1,PC_end+1)
  results_FDF<-phylopars_multi_anc(all_metrics_FDF1,phylopar_tree_copy,'FDF',PC_start+1,PC_end+1)
  results_FVF<-phylopars_multi_anc(all_metrics_FVF1,phylopar_tree_copy,'FVF',PC_start+1,PC_end+1)
  results_MDH<-phylopars_multi_anc(all_metrics_MDH1,phylopar_tree_copy,'MDH',PC_start+1,PC_end+1)
  results_MVH<-phylopars_multi_anc(all_metrics_MVH1,phylopar_tree_copy,'MVH',PC_start+1,PC_end+1)
  results_FDH<-phylopars_multi_anc(all_metrics_FDH1,phylopar_tree_copy,'FDH',PC_start+1,PC_end+1)
  results_FVH<-phylopars_multi_anc(all_metrics_FVH1,phylopar_tree_copy,'FVH',PC_start+1,PC_end+1)
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
  
  if(Taxon_group=='Butterflies'){
    tree_copy<-tree
    
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

###################################################################################################
##################################################################################################

#####bootstrapping proceeds from here-----
# startnum<-30
# batchnum<-50
# boot1<-((startnum-1)*batchnum)
# boot2<-((startnum*batchnum)-1)
# 
# boot<-1
# print(boot)
#####Setting up the output dataframes----
# initial_plotframe<-data.frame(Boot=rep(boot,8),Dataset=rep(c('DFW','DHW','VFW','VHW'),2),
#                               Sex_bias=c(rep('female',4),rep('male',4)))
# levels=unique(PC$Categorical_var)
# Final_plotframe=as.data.frame(expand(initial_plotframe,levels,Boot,Sex_bias,Dataset))
# Final_plotframe[paste(categorical_variable)]<-Final_plotframe$levels
# Final_plotframe['UID']<-paste(Final_plotframe[paste(categorical_variable)][[1]],Final_plotframe$Sex_bias,Final_plotframe$Dataset,sep='_')
# rownames(Final_plotframe)<-Final_plotframe$UID
# 
# initial_dimorph_frame<-data.frame(Boot=rep(boot,4),Dataset=c('DFW','DHW','VFW','VHW'))
# Final_dimorph=as.data.frame(expand(initial_dimorph_frame,levels,Boot,Dataset))
# Final_dimorph[paste(categorical_variable)]<-Final_dimorph$levels
# Final_dimorph['UID']<-paste(Final_dimorph[paste(categorical_variable)][[1]],Final_dimorph$Dataset,sep='_')
# rownames(Final_dimorph)<-Final_dimorph$UID

set.seed(rep)
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


######The following code sets the P-value thresholds for the Confidence intervals calculated below. The n_tests parameter is the key variable, and assumes that these are the number of statistical inferences 
#that will be compared against each other simultaneously (i.e. how many confidence intervals are plotted in the same graph that could be interpreted as "significantly different" from each other or not)
#for the purposes of this work, I consider this to be 8 (Moths/Butterflies X dorsal/ventral x FW/HW)


###I'm not sure the following code is actually used any more, as I use BCA confidence intervals not bonferroni
##might have to remove
n_taxa=3
# n_tests<-as.numeric(n_taxa*length(unique(PC2$Side))*length(unique(PC2$Wing)))
n_tests<-1
Bonferroni_conf_level<-1-0.05/n_tests

rates_m<-subset(rates1,Sex=='male')
rates_f<-subset(rates1,Sex=='female')
rownames(rates_m)<-rates_m$Phy_Side_Wing
rownames(rates_f)<-rates_f$Phy_Side_Wing
final_ids=intersect(rates_m$Phy_Side_Wing,rates_f$Phy_Side_Wing)
setkey(rates_m,Phy_Side_Wing)
setkey(rates_f,Phy_Side_Wing)

rates_m1<-rates_m[final_ids,]
rates_f1<-rates_f[final_ids,]
MF<-abs(rates_m1[,PC_colnames,with=FALSE])-abs(rates_f1[,PC_colnames,with=FALSE])
MF1<-cbind(rates_m1[,!names(rates_m1) %in% PC_colnames,with=FALSE],MF)
MF1$Sex<-'M-F'


##
##This is my fixed boot.pval function from the boot.pval package. It wasn't outputting BCA intervals correctly so I 
##fixed it myself. Note that in order to do so the minimum p-value cannot be lower than the 1/#bootstraps
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
median_MF_sexbias_boot<-function(indata,variable,group,cores,nboot){
  # print(metric)
  indata1<-setDT(indata)
  dat_out<-data.table()
  diff_df<-subset(indata1,Side=='dorsal'& Wing=='FW' & Categorical_var ==group,select=PC_colnames)
  diff_vf<-subset(indata1,Side=='ventral'& Wing=='FW' & Categorical_var ==group,select=PC_colnames)
  diff_dh<-subset(indata1,Side=='dorsal'& Wing=='HW' & Categorical_var ==group,select=PC_colnames)
  diff_vh<-subset(indata1,Side=='ventral'& Wing=='HW'& Categorical_var ==group,select=PC_colnames)
  median_func<-function(data,IDvec){
    MF=data[IDvec,]
    meds<-apply(MF,2,median)
    MF_F<-abs(sum(replace(meds,meds > 0, 0)))
    MF_M<-abs(sum(replace(meds,meds < 0, 0)))
    MF_prop_M<-MF_M/(MF_M+MF_F)
    return(MF_prop_M)
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
      P_val<-boot.pval_fixed(boot_out,theta_null=0.5,type='perc')
      CI_type='percentile'
    }else{
      boot_out<-boot(indata1,median_func,R=nboot,parallel='multicore',ncpus=cores)
      CI<-boot.ci(boot_out,conf=Bonferroni_conf_level,type='bca')
      CI_lo<-CI$bca[4]
      CI_hi<-CI$bca[5]
      P_val<-boot.pval_fixed(boot_out,theta_null=0.5,type='bca')
      CI_type='bca'
    } 
    return(list(CI_lo,CI_hi,P_val,CI_type))
  }
  meds<-apply(diff_df,2,median)
  MF_F<-abs(sum(replace(meds,meds > 0, 0)))
  MF_M<-abs(sum(replace(meds,meds < 0, 0)))
  overall_med_df<-MF_M/(MF_M+MF_F)  
  outs<-evaluate_boot_CI(overall_med_df,diff_df)
  CI_df_lo<-outs[1]
  CI_df_hi<-outs[2]
  P_val_df<-outs[3]
  CI_type_df<-outs[4]
  meds<-apply(diff_dh,2,median)
  MF_F<-abs(sum(replace(meds,meds > 0, 0)))
  MF_M<-abs(sum(replace(meds,meds < 0, 0)))
  overall_med_dh<-MF_M/(MF_M+MF_F)  
  outs<-evaluate_boot_CI(overall_med_dh,diff_dh)
  CI_dh_lo<-outs[1]
  CI_dh_hi<-outs[2]
  P_val_dh<-outs[3]
  CI_type_dh<-outs[4]
  meds<-apply(diff_vf,2,median)
  MF_F<-abs(sum(replace(meds,meds > 0, 0)))
  MF_M<-abs(sum(replace(meds,meds < 0, 0)))
  overall_med_vf<-MF_M/(MF_M+MF_F)  
  outs<-evaluate_boot_CI(overall_med_vf,diff_vf)
  CI_vf_lo<-outs[1]
  CI_vf_hi<-outs[2]
  P_val_vf<-outs[3]
  CI_type_vf<-outs[4]
  meds<-apply(diff_dh,2,median)
  MF_F<-abs(sum(replace(meds,meds > 0, 0)))
  MF_M<-abs(sum(replace(meds,meds < 0, 0)))
  overall_med_vh<-MF_M/(MF_M+MF_F)  
  outs<-evaluate_boot_CI(overall_med_vh,diff_vh)
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
  dat<-data.table('Taxon'=Taxon_group,'Categorical_variable'=group,'Dataset'=Dataset,'median'=overall_med,'CI_lo'=CI_lo,'CI_hi'=CI_hi,'P_val'=P_val,'CI_type'=CI_type)
  names(dat)[names(dat) == 'median'] <- paste(variable,'med',sep='_')
  names(dat)[names(dat) == 'CI_hi'] <- paste(variable,'CI_hi',sep='_')
  names(dat)[names(dat) == 'CI_lo'] <- paste(variable,'CI_lo',sep='_')
  names(dat)[names(dat) == 'P_val'] <- paste(variable,'P_val',sep='_')
  names(dat)[names(dat) == 'CI_type'] <- paste(variable,'CI_type',sep='_')
  dat_out<-rbind(dat_out,dat)
  return(dat_out)
}
t1<-Sys.time()
MF_abs_rate_sexbias<-data.table()
for(i in unique(PC$Categorical_var)){
  boot_out<-median_MF_sexbias_boot(MF1,'MF_abs_rate_sexbias',i,n_cores,boot_reps)
  MF_abs_rate_sexbias<-rbind(MF_abs_rate_sexbias,boot_out)
}
print('Finished MF abs_rate prop sexbias')
t2<-Sys.time()
t2-t1


dimorph_ssq_boot<-function(indata,group,cores,nboot){
  dat_out<-data.table()
  MDF<-subset(indata,Wing=='FW'& Side=='dorsal' & Sex=='male' & Categorical_var ==group)
  MVF<-subset(indata,Wing=='FW'& Side=='ventral' & Sex=='male' & Categorical_var ==group)
  MDH<-subset(indata,Wing=='HW'& Side=='dorsal' & Sex=='male' & Categorical_var ==group)
  MVH<-subset(indata,Wing=='HW'& Side=='ventral' & Sex=='male' & Categorical_var ==group)
  FDF<-subset(indata,Wing=='FW'& Side=='dorsal' & Sex=='female' & Categorical_var ==group)
  FVF<-subset(indata,Wing=='FW'& Side=='ventral' & Sex=='female' & Categorical_var ==group)
  FDH<-subset(indata,Wing=='HW'& Side=='dorsal' & Sex=='female' & Categorical_var ==group)
  FVH<-subset(indata,Wing=='HW'& Side=='ventral' & Sex=='female' & Categorical_var ==group)
  SSQ_func<-function(data,IDvec){
    dats<-data[IDvec,]
    ss_T<-vector()
    ss_tot<-vector()
    n<-vector()
    sec_moment_distance<-vector()
    dats_m<-dats[,1:length(PC_colnames)]
    dats_f<-dats[,(1+length(PC_colnames)):(length(PC_colnames)+length(PC_colnames))]
    for(i in 1:length(PC_colnames)){
      tm<-dats_m[[i]]
      tf<-dats_f[[i]]
      trial<-c(tm,tf)
      MF<-abs(tm-tf)
      grand_mean<-mean(trial)
      mean_m<-mean(tm)
      mean_f<-mean(tf)
      ss_tm<-sum(((mean_m-grand_mean)^2)*length(tm))
      ss_tf<-sum(((mean_f-grand_mean)^2)*length(tf))
      ss_T<-append(ss_T,ss_tm+ss_tf)
      ss_tot<-append(ss_tot,sum((trial-grand_mean)^2))
      sec_moment_distance<-append(sec_moment_distance,(sum((MF-0)^2))/length(MF))
      n<-append(n,as.numeric(length(trial)))
    }
    SSQ_tot<-sum(ss_tot)
    ssT_to_tot<-sum(ss_T)/SSQ_tot
    var_tot<-SSQ_tot/mean(n)
    sec_moment_distance<-sum(sec_moment_distance)
    return(c(ssT_to_tot,var_tot,sec_moment_distance,SSQ_tot))
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
    dt_out<-cbind(c('ssT_to_tot','var_tot','sec_moment_distance','SSQ_tot'),dt_out)
    colnames(dt_out)<-c('metric','CI_lo','CI_hi','CI_type')
    return(dt_out)
  }
  tt<-cbind(MDF[,PC_colnames,with=FALSE],FDF[,PC_colnames,with=FALSE])
  overall_diff_df<-SSQ_func(tt,(1:nrow(tt)))
  outs<-evaluate_boot_CI(overall_diff_df,tt)
  CI_df_lo<-outs$CI_lo
  CI_df_hi<-outs$CI_hi
  CI_type_df<-outs$CI_type
  tt<-cbind(MDH[,PC_colnames,with=FALSE],FDH[,PC_colnames,with=FALSE])
  overall_diff_dh<-SSQ_func(tt,(1:nrow(tt)))
  outs<-evaluate_boot_CI(overall_diff_dh,tt)
  CI_dh_lo<-outs$CI_lo
  CI_dh_hi<-outs$CI_hi
  CI_type_dh<-outs$CI_type
  tt<-cbind(MVF[,PC_colnames,with=FALSE],FVF[,PC_colnames,with=FALSE])
  overall_diff_vf<-SSQ_func(tt,(1:nrow(tt)))
  outs<-evaluate_boot_CI(overall_diff_vf,tt)
  CI_vf_lo<-outs$CI_lo
  CI_vf_hi<-outs$CI_hi
  CI_type_vf<-outs$CI_type
  tt<-cbind(MVH[,PC_colnames,with=FALSE],FVH[,PC_colnames,with=FALSE])
  overall_diff_vh<-SSQ_func(tt,(1:nrow(tt)))
  outs<-evaluate_boot_CI(overall_diff_vh,tt)
  CI_vh_lo<-outs$CI_lo
  CI_vh_hi<-outs$CI_hi
  CI_type_vh<-outs$CI_type
  final_out<-data.table('overall_diff'=c(overall_diff_df,overall_diff_dh,overall_diff_vf,overall_diff_vh),
                        'CI_lo'=c(CI_df_lo,CI_dh_lo,CI_vf_lo,CI_vh_lo), 'CI_hi'=c(CI_df_hi,CI_dh_hi,CI_vf_hi,CI_vh_hi),
                        'CI_type'=c(CI_type_df,CI_type_dh,CI_type_vf,CI_type_vh))
  final_out$Dataset<-c(rep('DFW',4),rep('DHW',4),rep('VFW',4),rep('VHW',4))
  final_out$value<-rep(c('ssT_to_tot','var_tot','sec_mom_dist','SSQ_tot'),4)
  dat<-data.table('Taxon'=Taxon_group,'Categorical_variable'=group,'Dataset'=c('DFW','DHW','VFW','VHW'),
                  'ssT_to_tot_med'=final_out[value=='ssT_to_tot']$overall_diff,'ssT_to_tot_CI_lo'=final_out[value=='ssT_to_tot']$CI_lo,'ssT_to_tot_CI_hi'=final_out[value=='ssT_to_tot']$CI_hi,'ssT_to_tot_CI_type'=final_out[value=='ssT_to_tot']$CI_type,
                  'tot_var_med'=final_out[value=='var_tot']$overall_diff,'tot_var_CI_lo'=final_out[value=='var_tot']$CI_lo,'tot_var_CI_hi'=final_out[value=='var_tot']$CI_hi,'tot_var_CI_type'=final_out[value=='var_tot']$CI_type,
                  'sec_mom_var_med'=final_out[value=='sec_mom_dist']$overall_diff,'sec_mom_var_CI_lo'=final_out[value=='sec_mom_dist']$CI_lo,'sec_mom_var_CI_hi'=final_out[value=='sec_mom_dist']$CI_hi,'sec_mom_var_CI_type'=final_out[value=='sec_mom_dist']$CI_type,
                  'SSQ_tot_med'=final_out[value=='SSQ_tot']$overall_diff,'SSQ_tot_CI_lo'=final_out[value=='SSQ_tot']$CI_lo,'SSQ_tot_CI_hi'=final_out[value=='SSQ_tot']$CI_hi,'SSQ_tot_CI_type'=final_out[value=='SSQ_tot']$CI_type)
  dat_out<-rbind(dat_out,dat)
  return(dat_out)
}

# SSQ_frame<-dimorph_ssq_boot(PC3,categorical_variable,n_cores,metric,boot_reps)
t1<-Sys.time()
SSQ_frame<-data.table()
for(i in unique(PC$Categorical_var)){
  boot_out<-dimorph_ssq_boot(PC2,i,n_cores,boot_reps)
  SSQ_frame<-rbind(SSQ_frame,boot_out)
}
print('Finished SSQ values')
t2<-Sys.time()
t1-t2

#####Rate SSQ sexbias boot function----
rate_ssq_boot<-function(indata,variable,group,cores,nboot){
  dat_out<-data.table()
  MDF<-subset(indata,Wing=='FW'& Side=='dorsal' & Sex=='male' & Categorical_var ==group)
  MVF<-subset(indata,Wing=='FW'& Side=='ventral' & Sex=='male' & Categorical_var ==group)
  MDH<-subset(indata,Wing=='HW'& Side=='dorsal' & Sex=='male' & Categorical_var ==group)
  MVH<-subset(indata,Wing=='HW'& Side=='ventral' & Sex=='male' & Categorical_var ==group)
  FDF<-subset(indata,Wing=='FW'& Side=='dorsal' & Sex=='female' & Categorical_var ==group)
  FVF<-subset(indata,Wing=='FW'& Side=='ventral' & Sex=='female' & Categorical_var ==group)
  FDH<-subset(indata,Wing=='HW'& Side=='dorsal' & Sex=='female' & Categorical_var ==group)
  FVH<-subset(indata,Wing=='HW'& Side=='ventral' & Sex=='female' & Categorical_var ==group)
  SSQ_func<-function(data,IDvec){
    dats<-data[IDvec,]
    ss_T<-vector()
    ss_tot<-vector()
    n<-vector()
    sec_moment_distance<-vector()
    MF_abs_rate<-vector()
    dats_m<-dats[,1:length(PC_colnames)]
    dats_f<-dats[,(1+length(PC_colnames)):(length(PC_colnames)+length(PC_colnames))]
    for(i in 1:length(PC_colnames)){
      tm<-dats_m[[i]]
      tf<-dats_f[[i]]
      trial<-c(tm,tf)
      MF<-abs(tm-tf)
      PC_MF<-abs(tm)-abs(tf)
      MF_abs_rate<-append(MF_abs_rate,median(PC_MF))
      grand_mean<-mean(abs(trial))
      mean_m<-mean(abs(tm))
      mean_f<-mean(abs(tf))
      ss_tm<-sum(((mean_m-grand_mean)^2)*length(tm))
      ss_tf<-sum(((mean_f-grand_mean)^2)*length(tf))
      ss_T<-append(ss_T,ss_tm+ss_tf)
      ss_tot<-append(ss_tot,sum((trial-grand_mean)^2))
      sec_moment_distance<-append(sec_moment_distance,(sum((MF-0)^2))/length(MF))
      n<-append(n,as.numeric(length(trial)))
    }
    dat_out<-data.table('ss_T'=ss_T,'MF_abs_rate'=MF_abs_rate)
    MF_M<-sum(dat_out[MF_abs_rate>0]$ss_T)
    MF_F<-sum(dat_out[MF_abs_rate<0]$ss_T)
    MF_SSQ_rates_prop_sexbias<-MF_M/(MF_M+MF_F)
    var_tot<-sum(ss_tot)/mean(n)
    MF_SSQ<-sum(ss_T)
    sec_moment_distance<-sum(sec_moment_distance)
    return(c(MF_SSQ_rates_prop_sexbias,MF_SSQ))
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
    dt_out<-cbind(c('MF_SSQ_rates_prop_sexbias','MF_SSQ'),dt_out)
    colnames(dt_out)<-c('metric','CI_lo','CI_hi','CI_type')
    return(dt_out)
  }
  tt<-cbind(MDF[,PC_colnames,with=FALSE],FDF[,PC_colnames,with=FALSE])
  overall_diff_df<-SSQ_func(tt,(1:nrow(tt)))
  outs<-evaluate_boot_CI(overall_diff_df,tt)
  CI_df_lo<-outs$CI_lo
  CI_df_hi<-outs$CI_hi
  CI_type_df<-outs$CI_type
  tt<-cbind(MDH[,PC_colnames,with=FALSE],FDH[,PC_colnames,with=FALSE])
  overall_diff_dh<-SSQ_func(tt,(1:nrow(tt)))
  outs<-evaluate_boot_CI(overall_diff_dh,tt)
  CI_dh_lo<-outs$CI_lo
  CI_dh_hi<-outs$CI_hi
  CI_type_dh<-outs$CI_type
  tt<-cbind(MVF[,PC_colnames,with=FALSE],FVF[,PC_colnames,with=FALSE])
  overall_diff_vf<-SSQ_func(tt,(1:nrow(tt)))
  outs<-evaluate_boot_CI(overall_diff_vf,tt)
  CI_vf_lo<-outs$CI_lo
  CI_vf_hi<-outs$CI_hi
  CI_type_vf<-outs$CI_type
  tt<-cbind(MVH[,PC_colnames,with=FALSE],FVH[,PC_colnames,with=FALSE])
  overall_diff_vh<-SSQ_func(tt,(1:nrow(tt)))
  outs<-evaluate_boot_CI(overall_diff_vh,tt)
  CI_vh_lo<-outs$CI_lo
  CI_vh_hi<-outs$CI_hi
  CI_type_vh<-outs$CI_type
  final_out<-data.table('overall_diff'=c(overall_diff_df,overall_diff_dh,overall_diff_vf,overall_diff_vh),
                        'CI_lo'=c(CI_df_lo,CI_dh_lo,CI_vf_lo,CI_vh_lo), 'CI_hi'=c(CI_df_hi,CI_dh_hi,CI_vf_hi,CI_vh_hi),
                        'CI_type'=c(CI_type_df,CI_type_dh,CI_type_vf,CI_type_vh))
  final_out$Dataset<-c(rep('DFW',2),rep('DHW',2),rep('VFW',2),rep('VHW',2))
  final_out$value<-rep(c('MF_SSQ_rates_prop_sexbias','MF_SSQ'),4)
  dat<-data.table('Taxon'=Taxon_group,'Categorical_variable'=group,'Dataset'=c('DFW','DHW','VFW','VHW'),
                  'MF_SSQ_rates_prop_sexbias_med'=final_out[value=='MF_SSQ_rates_prop_sexbias']$overall_diff,'MF_SSQ_rates_prop_sexbias_CI_lo'=final_out[value=='MF_SSQ_rates_prop_sexbias']$CI_lo,'MF_SSQ_rates_prop_sexbias_CI_hi'=final_out[value=='MF_SSQ_rates_prop_sexbias']$CI_hi,'MF_SSQ_rates_prop_sexbias_CI_type'=final_out[value=='MF_SSQ_rates_prop_sexbias']$CI_type,
                  'MF_SSQ_med'=final_out[value=='MF_SSQ']$overall_diff,'MF_SSQ_CI_lo'=final_out[value=='MF_SSQ']$CI_lo,'MF_SSQ_CI_hi'=final_out[value=='MF_SSQ']$CI_hi,'MF_SSQ_CI_type'=final_out[value=='MF_SSQ']$CI_type)
  dat_out<-rbind(dat_out,dat)
}



# SSQ_frame<-dimorph_ssq_boot(PC3,categorical_variable,n_cores,metric,boot_reps)
t1<-Sys.time()
rate_SSQ_frame<-data.table()
for(i in unique(PC$Categorical_var)){
  boot_out<-rate_ssq_boot(rates1,'MF_SSQ_rates_prop_sexbias',i,n_cores,boot_reps)
  rate_SSQ_frame<-rbind(rate_SSQ_frame,boot_out)
}
print('Finished rate SSQ values')
t2<-Sys.time()
t2-t1



#####Morphological disparity bias bootstrap----
MF_MD_sexbias_boot<-function(indata,variable,group,cores,nboot){
  dat_out<-data.table()
  MDF<-subset(indata,Wing=='FW'& Side=='dorsal' & Sex=='male' & Categorical_var ==group)
  MVF<-subset(indata,Wing=='FW'& Side=='ventral' & Sex=='male' & Categorical_var ==group)
  MDH<-subset(indata,Wing=='HW'& Side=='dorsal' & Sex=='male' & Categorical_var ==group)
  MVH<-subset(indata,Wing=='HW'& Side=='ventral' & Sex=='male' & Categorical_var ==group)
  FDF<-subset(indata,Wing=='FW'& Side=='dorsal' & Sex=='female' & Categorical_var ==group)
  FVF<-subset(indata,Wing=='FW'& Side=='ventral' & Sex=='female' & Categorical_var ==group)
  FDH<-subset(indata,Wing=='HW'& Side=='dorsal' & Sex=='female' & Categorical_var ==group)
  FVH<-subset(indata,Wing=='HW'& Side=='ventral' & Sex=='female' & Categorical_var ==group)
  disp_func<-function(data,IDvec){
    dats<-data[IDvec,]
    dats_m<-dats[,1:length(PC_colnames)]
    dats_f<-dats[,(1+length(PC_colnames)):(length(PC_colnames)+length(PC_colnames))]
    trial<-rbind(dats_m,dats_f)
    lenm<-as.numeric(nrow(dats_m))
    lenf<-as.numeric(nrow(dats_f))
    # try(disp<-dispRity.per.group(bootstraps=10,data=matrix(trial),list("Male"=c(1:lenm),"Female"=c((lenm+1):(lenm+lenf))),
    # metric=c(variances)),silent=FALSE)
    try(disp<-dispRity(boot.matrix(custom.subsets(setDF(trial), list("Male"=c(1:lenm),"Female"=c((lenm+1):(lenm+lenf)))),
                                   bootstraps = 1),
                       metric = c(variances)))
    male_results<-disp$disparity$Male$elements
    female_results<-disp$disparity$Female$elements
    M_minus_F_disp<-male_results-female_results
    MF_F<-abs(sum(replace(M_minus_F_disp, M_minus_F_disp > 0, 0), na.rm = TRUE))
    MF_M<-abs(sum(replace(M_minus_F_disp, M_minus_F_disp < 0, 0), na.rm = TRUE))
    MF_extant_var_MD<-c(MF_F,MF_M)
    MF_extant_var_MD_prop_sexbias<-MF_M/(MF_M+MF_F)
    return(MF_extant_var_MD_prop_sexbias)
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
      P_val<-boot.pval_fixed(boot_out,theta_null=0.5,type='perc' )
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
  tt<-cbind(MDF[,PC_colnames,with=FALSE],FDF[,PC_colnames,with=FALSE])
  overall_prop_df<-disp_func(tt,(1:nrow(tt)))
  outs<-evaluate_boot_CI(overall_prop_df,tt)
  CI_df_lo<-outs[1]
  CI_df_hi<-outs[2]
  P_val_df<-outs[3]
  CI_type_df<-outs[4]
  tt<-cbind(MDH[,PC_colnames,with=FALSE],FDH[,PC_colnames,with=FALSE])
  overall_prop_dh<-disp_func(tt,(1:nrow(tt)))
  outs<-evaluate_boot_CI(overall_prop_dh,tt)
  CI_dh_lo<-outs[1]
  CI_dh_hi<-outs[2]
  P_val_dh<-outs[3]
  CI_type_dh<-outs[4]
  tt<-cbind(MVF[,PC_colnames,with=FALSE],FVF[,PC_colnames,with=FALSE])
  overall_prop_vf<-disp_func(tt,(1:nrow(tt)))
  outs<-evaluate_boot_CI(overall_prop_vf,tt)
  CI_vf_lo<-outs[1]
  CI_vf_hi<-outs[2]
  P_val_vf<-outs[3]
  CI_type_vf<-outs[4]
  tt<-cbind(MVH[,PC_colnames,with=FALSE],FVH[,PC_colnames,with=FALSE])
  overall_prop_vh<-disp_func(tt,(1:nrow(tt)))
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
  dat<-data.table('Taxon'=Taxon_group,'Categorical_variable'=group,'Dataset'=Dataset,'med'=overall_prop,'CI_lo'=CI_lo,'CI_hi'=CI_hi,'P_val'=P_val,'CI_type'=CI_type)
  names(dat)[names(dat) == 'med'] <- paste(variable,'med',sep='_')
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
MF_extant_var_MD_prop_sexbias<-data.table()
for(i in unique(PC$Categorical_var)){
  boot_out<-MF_MD_sexbias_boot(PC2,'MF_extant_var_MD_prop_sexbias',i,n_cores,boot_reps)
  MF_extant_var_MD_prop_sexbias<-rbind(MF_extant_var_MD_prop_sexbias,boot_out)
}
print('Finished MF MD sexbias')
t2<-Sys.time()
t2-t1

# ######Now making the 3m3f intraspecific variation stats----
# PC_IS_var<-ddply(all_spec_metrics1,.(Phy_Side_Wing,ID_var,Family,Categorical_var,Sex,Side,Wing),colwise(var))
# PC_IS_count<-ddply(all_spec_metrics1,.(Phy_Side_Wing,ID_var,Family,Categorical_var,Sex,Side,Wing),colwise(length))
# PC_IS_med<-ddply(all_spec_metrics1,.(Phy_Side_Wing,ID_var,Family,Categorical_var,Sex,Side,Wing),colwise(median))
# cat_count<-ddply(PC_IS_var,.(Categorical_var,Sex,Side,Wing),colwise(length))
# cat_count_filt<-subset(cat_count,ID_var>1)
# PC_IS_var<-PC_IS_var[PC_IS_var$Categorical_var %in% unique(cat_count_filt$Categorical_var),]
# PC_IS_med<-PC_IS_med[PC_IS_med$Categorical_var %in% unique(cat_count_filt$Categorical_var),]
# PC_IS_count<-PC_IS_count[PC_IS_count$Categorical_var %in% unique(cat_count_filt$Categorical_var),]
# 
# PC_IS_var<-PC_IS_var[,c(PC_catnames,PC_colnames)]
# PC_IS_med<-PC_IS_med[,c(PC_catnames,PC_colnames)]
# count<-PC_IS_count$Specimen
# PC_IS_var1<-cbind(count,PC_IS_var)
# PC_IS_var2<-PC_IS_var1[PC_IS_var1$count>(IS_var_thresh-1),]
# PC_IS_var2_M<-subset(PC_IS_var2,Sex=='male')
# PC_IS_var2_F<-subset(PC_IS_var2,Sex=='female')
# rownames(PC_IS_var2_M)<-PC_IS_var2_M$Phy_Side_Wing
# rownames(PC_IS_var2_F)<-PC_IS_var2_F$Phy_Side_Wing
# final_ids=intersect(PC_IS_var2_M$Phy_Side_Wing,PC_IS_var2_F$Phy_Side_Wing)
# PC_IS_var3_M<-PC_IS_var2_M[final_ids,]
# PC_IS_var3_F<-PC_IS_var2_F[final_ids,]
# 
# PC_IS_var3_MF<-cbind(PC_IS_var3_M[,!names(PC_IS_var3_M) %in% PC_colnames],PC_IS_var3_M[final_ids,PC_colnames]-PC_IS_var3_F[final_ids,PC_colnames])
# PC_IS_var3_MF$Sex<-'M-F'

# #####The main MF_ISvar_prop_sexbias bootstrap results----
# 
# # MF_difference<-median_sexdiff_boot(PC3,'MF_difference',categorical_variable,n_cores,metric,boot_reps)
# t1<-Sys.time()
# MF_IS_var_prop_sexbias<-data.table()
# for(i in unique(PC$Categorical_var)){
#   boot_out<-median_MF_sexbias_boot(PC_IS_var3_MF,'MF_IS_var_prop_sexbias',i,n_cores,boot_reps)
#   MF_IS_var_prop_sexbias<-rbind(MF_IS_var_prop_sexbias,boot_out)
# }
# print('Finished MF ISvar prop sexbias')
# 
# t2<-Sys.time()
# t2-t1

##########Now the absolute rate sexbias
median_sexbias_boot<-function(indata,variable,group,cores,nboot){
  # print(metric)
  indata1<-setDT(indata)
  dat_out<-data.table()
  MDF<-abs(subset(indata1,Wing=='FW'& Side=='dorsal' & Sex=='male' & Categorical_var ==group,select=PC_colnames))
  MVF<-abs(subset(indata1,Wing=='FW'& Side=='ventral' & Sex=='male' & Categorical_var ==group,select=PC_colnames))
  MDH<-abs(subset(indata1,Wing=='HW'& Side=='dorsal' & Sex=='male' & Categorical_var ==group,select=PC_colnames))
  MVH<-abs(subset(indata1,Wing=='HW'& Side=='ventral' & Sex=='male' & Categorical_var ==group,select=PC_colnames))
  FDF<-abs(subset(indata1,Wing=='FW'& Side=='dorsal' & Sex=='female' & Categorical_var ==group,select=PC_colnames))
  FVF<-abs(subset(indata1,Wing=='FW'& Side=='ventral' & Sex=='female' & Categorical_var ==group,select=PC_colnames))
  FDH<-abs(subset(indata1,Wing=='HW'& Side=='dorsal' & Sex=='female' & Categorical_var ==group,select=PC_colnames))
  FVH<-abs(subset(indata1,Wing=='HW'& Side=='ventral' & Sex=='female' & Categorical_var ==group,select=PC_colnames))
  median_bias_func<-function(data,IDvec){
    dats<-data[IDvec,]
    dats_m<-abs(dats[,1:length(PC_colnames)])
    dats_f<-abs(dats[,(1+length(PC_colnames)):(length(PC_colnames)+length(PC_colnames))])
    meds_F<-apply(dats_f,2,median)
    meds_M<-apply(dats_m,2,median)
    prop_M_bias<-sum(meds_M)/(sum(meds_M)+sum(meds_F))
    return(prop_M_bias)
  }
  evaluate_boot_CI<-function(overall_value,indata1){
    if(is.na(overall_value)==TRUE){
      CI_lo<-NA
      CI_hi<-NA
      P_val<-NA
      CI_type=NA
    }else if(overall_value==0){
      boot_out<-boot(indata1,median_bias_func,R=nboot,parallel='multicore',ncpus=cores)
      CI<-boot.ci(boot_out,conf=Bonferroni_conf_level,type='perc')
      CI_lo<-CI$percent[4]
      CI_hi<-CI$percent[5]
      P_val<-boot.pval_fixed(boot_out,theta_null=0.5,type='perc')
      CI_type='percentile'
    }else{
      boot_out<-boot(indata1,median_bias_func,R=nboot,parallel='multicore',ncpus=cores)
      CI<-boot.ci(boot_out,conf=Bonferroni_conf_level,type='bca')
      CI_lo<-CI$bca[4]
      CI_hi<-CI$bca[5]
      P_val<-boot.pval_fixed(boot_out,theta_null=0.5,type='bca')
      CI_type='bca'
    } 
    return(list(CI_lo,CI_hi,P_val,CI_type))
  }
  tt<-cbind(MDF[,PC_colnames,with=FALSE],FDF[,PC_colnames,with=FALSE])
  overall_prop_df<-median_bias_func(tt,(1:nrow(tt)))
  outs<-evaluate_boot_CI(overall_prop_df,tt)
  CI_df_lo<-outs[1]
  CI_df_hi<-outs[2]
  P_val_df<-outs[3]
  CI_type_df<-outs[4]
  tt<-cbind(MDH[,PC_colnames,with=FALSE],FDH[,PC_colnames,with=FALSE])
  overall_prop_dh<-median_bias_func(tt,(1:nrow(tt)))
  outs<-evaluate_boot_CI(overall_prop_dh,tt)
  CI_dh_lo<-outs[1]
  CI_dh_hi<-outs[2]
  P_val_dh<-outs[3]
  CI_type_dh<-outs[4]
  tt<-cbind(MVF[,PC_colnames,with=FALSE],FVF[,PC_colnames,with=FALSE])
  overall_prop_vf<-median_bias_func(tt,(1:nrow(tt)))
  outs<-evaluate_boot_CI(overall_prop_vf,tt)
  CI_vf_lo<-outs[1]
  CI_vf_hi<-outs[2]
  P_val_vf<-outs[3]
  CI_type_vf<-outs[4]
  tt<-cbind(MVH[,PC_colnames,with=FALSE],FVH[,PC_colnames,with=FALSE])
  overall_prop_vh<-median_bias_func(tt,(1:nrow(tt)))
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
  dat<-data.table('Taxon'=Taxon_group,'Categorical_variable'=group,'Dataset'=Dataset,'median'=overall_prop,'CI_lo'=CI_lo,'CI_hi'=CI_hi,'P_val'=P_val,'CI_type'=CI_type)
  names(dat)[names(dat) == 'median'] <- paste(variable,'med',sep='_')
  names(dat)[names(dat) == 'CI_hi'] <- paste(variable,'CI_hi',sep='_')
  names(dat)[names(dat) == 'CI_lo'] <- paste(variable,'CI_lo',sep='_')
  names(dat)[names(dat) == 'P_val'] <- paste(variable,'P_val',sep='_')
  names(dat)[names(dat) == 'CI_type'] <- paste(variable,'CI_type',sep='_')
  dat_out<-rbind(dat_out,dat)
  return(dat_out)
}
t1<-Sys.time()
Abs_rate_prop_sexbias<-data.table()
for(i in unique(PC$Categorical_var)){
  boot_out<-median_sexbias_boot(rates1,'Abs_rate_prop_sexbias',i,n_cores,boot_reps)
  Abs_rate_prop_sexbias<-rbind(Abs_rate_prop_sexbias,boot_out)
}
t2<-Sys.time()
t2-t1


# #####Now the ISvar sexbias
# PC_IS_var3<-rbind(PC_IS_var3_M,PC_IS_var3_F)
# t1<-Sys.time()
# IS_var_prop_sexbias<-data.table()
# for(i in unique(PC$Categorical_var)){
#   boot_out<-median_sexbias_boot(PC_IS_var3,'IS_var_prop_sexbias',i,n_cores,boot_reps)
#   IS_var_prop_sexbias<-rbind(IS_var_prop_sexbias,boot_out)
# }
# print('Finished ISvar prop sexbias values')
# 
# t2<-Sys.time()
# t2-t1

#####Now the MD sexbias
MD_sexbias_boot<-function(indata,variable,group,cores,nboot){
  dat_out<-data.table()
  indata1<-setDT(indata)
  MDF<-subset(indata1,Wing=='FW'& Side=='dorsal' & Sex=='male' & Categorical_var ==group)
  MVF<-subset(indata1,Wing=='FW'& Side=='ventral' & Sex=='male' & Categorical_var ==group)
  MDH<-subset(indata1,Wing=='HW'& Side=='dorsal' & Sex=='male' & Categorical_var ==group)
  MVH<-subset(indata1,Wing=='HW'& Side=='ventral' & Sex=='male' & Categorical_var ==group)
  FDF<-subset(indata1,Wing=='FW'& Side=='dorsal' & Sex=='female' & Categorical_var ==group)
  FVF<-subset(indata1,Wing=='FW'& Side=='ventral' & Sex=='female' & Categorical_var ==group)
  FDH<-subset(indata1,Wing=='HW'& Side=='dorsal' & Sex=='female' & Categorical_var ==group)
  FVH<-subset(indata1,Wing=='HW'& Side=='ventral' & Sex=='female' & Categorical_var ==group)
  disp_func<-function(data,IDvec){
    dats<-data[IDvec,]
    dats_m<-dats[,1:length(PC_colnames)]
    dats_f<-dats[,(1+length(PC_colnames)):(length(PC_colnames)+length(PC_colnames))]
    trial<-rbind(dats_m,dats_f)
    lenm<-as.numeric(nrow(dats_m))
    lenf<-as.numeric(nrow(dats_f))
    try(disp<-dispRity(boot.matrix(custom.subsets(setDF(trial), list("Male"=c(1:lenm),"Female"=c((lenm+1):(lenm+lenf)))),
                                   bootstraps = 1),
                       metric = c(variances)))
    male_results<-disp$disparity$Male$elements
    female_results<-disp$disparity$Female$elements
    sum_F<-abs(sum(female_results))
    sum_M<-abs(sum(male_results))
    extant_var_MD_prop_sexbias<-sum_M/(sum_M+sum_F)
    return(extant_var_MD_prop_sexbias)
  }
  evaluate_boot_CI<-function(overall_value,indata1){
    if(is.na(overall_value)==TRUE){
      CI_lo<-NA
      CI_hi<-NA
      P_val<-NA
      CI_type=NA
    }else if(overall_value==0){
      boot_out<-boot(indata1,disp_func,R=nboot,parallel='multicore',ncpus=cores)
      CI<-boot.ci(boot_out,conf=Bonferroni_conf_level,type='perc')
      CI_lo<-CI$percent[4]
      CI_hi<-CI$percent[5]
      P_val<-boot.pval_fixed(boot_out,theta_null=0.5,type='perc')
      CI_type='percentile'
    }else{
      boot_out<-boot(indata1,disp_func,R=nboot,parallel='multicore',ncpus=cores)
      CI<-boot.ci(boot_out,conf=Bonferroni_conf_level,type='bca')
      CI_lo<-CI$bca[4]
      CI_hi<-CI$bca[5]
      P_val<-boot.pval_fixed(boot_out,theta_null=0.5,type='bca')
      CI_type='bca'
    } 
    return(list(CI_lo,CI_hi,P_val,CI_type))
  }
  tt<-cbind(MDF[,PC_colnames,with=FALSE],FDF[,PC_colnames,with=FALSE])
  overall_prop_df<-disp_func(tt,(1:nrow(tt)))
  outs<-evaluate_boot_CI(overall_prop_df,tt)
  CI_df_lo<-outs[1]
  CI_df_hi<-outs[2]
  P_val_df<-outs[3]
  CI_type_df<-outs[4]
  tt<-cbind(MDH[,PC_colnames,with=FALSE],FDH[,PC_colnames,with=FALSE])
  overall_prop_dh<-disp_func(tt,(1:nrow(tt)))
  outs<-evaluate_boot_CI(overall_prop_dh,tt)
  CI_dh_lo<-outs[1]
  CI_dh_hi<-outs[2]
  P_val_dh<-outs[3]
  CI_type_dh<-outs[4]
  tt<-cbind(MVF[,PC_colnames,with=FALSE],FVF[,PC_colnames,with=FALSE])
  overall_prop_vf<-disp_func(tt,(1:nrow(tt)))
  outs<-evaluate_boot_CI(overall_prop_vf,tt)
  CI_vf_lo<-outs[1]
  CI_vf_hi<-outs[2]
  P_val_vf<-outs[3]
  CI_type_vf<-outs[4]
  tt<-cbind(MVH[,PC_colnames,with=FALSE],FVH[,PC_colnames,with=FALSE])
  overall_prop_vh<-disp_func(tt,(1:nrow(tt)))
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
  dat<-data.table('Taxon'=Taxon_group,'Categorical_variable'=group,'Dataset'=Dataset,'med'=overall_prop,'CI_lo'=CI_lo,'CI_hi'=CI_hi,'P_val'=P_val,'CI_type'=CI_type)
  names(dat)[names(dat) == 'med'] <- paste(variable,'med',sep='_')
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
Extant_var_MD_prop_sexbias<-data.table()
for(i in unique(PC$Categorical_var)){
  boot_out<-MD_sexbias_boot(PC2,'Extant_var_MD_prop_sexbias',i,n_cores,boot_reps)
  Extant_var_MD_prop_sexbias<-rbind(Extant_var_MD_prop_sexbias,boot_out)
}
print('Finished MD sexbias')
t2<-Sys.time()
t2-t1

#####Now making the IS_var MF dist proportion total bootstraps----
median_dist_boot<-function(indata,variable,group,cores,nboot){
  # print(metric)
  indata1<-setDT(indata)
  dat_out<-data.table()
  MDF<-abs(subset(indata1,Wing=='FW'& Side=='dorsal' & Sex=='male' & Categorical_var ==group,select=PC_colnames))
  MVF<-abs(subset(indata1,Wing=='FW'& Side=='ventral' & Sex=='male' & Categorical_var ==group,select=PC_colnames))
  MDH<-abs(subset(indata1,Wing=='HW'& Side=='dorsal' & Sex=='male' & Categorical_var ==group,select=PC_colnames))
  MVH<-abs(subset(indata1,Wing=='HW'& Side=='ventral' & Sex=='male' & Categorical_var ==group,select=PC_colnames))
  FDF<-abs(subset(indata1,Wing=='FW'& Side=='dorsal' & Sex=='female' & Categorical_var ==group,select=PC_colnames))
  FVF<-abs(subset(indata1,Wing=='FW'& Side=='ventral' & Sex=='female' & Categorical_var ==group,select=PC_colnames))
  FDH<-abs(subset(indata1,Wing=='HW'& Side=='dorsal' & Sex=='female' & Categorical_var ==group,select=PC_colnames))
  FVH<-abs(subset(indata1,Wing=='HW'& Side=='ventral' & Sex=='female' & Categorical_var ==group,select=PC_colnames))
  median_bias_func<-function(data,IDvec){
    dats<-data[IDvec,]
    dats_m<-abs(dats[,1:length(PC_colnames)])
    dats_f<-abs(dats[,(1+length(PC_colnames)):(length(PC_colnames)+length(PC_colnames))])
    meds_F<-apply(dats_f,2,median)
    meds_M<-apply(dats_m,2,median)
    dist<-abs(sum(meds_M)-sum(meds_F))
    return(dist)
  }
  evaluate_boot_CI<-function(overall_value,indata1){
    if(is.na(overall_value)==TRUE){
      CI_lo<-NA
      CI_hi<-NA
      CI_type=NA
    }else if(overall_value==0){
      boot_out<-boot(indata1,median_bias_func,R=nboot,parallel='multicore',ncpus=cores)
      CI<-boot.ci(boot_out,conf=Bonferroni_conf_level,type='perc')
      CI_lo<-CI$percent[4]
      CI_hi<-CI$percent[5]
      CI_type='percentile'
    }else{
      boot_out<-boot(indata1,median_bias_func,R=nboot,parallel='multicore',ncpus=cores)
      CI<-boot.ci(boot_out,conf=Bonferroni_conf_level,type='bca')
      CI_lo<-CI$bca[4]
      CI_hi<-CI$bca[5]
      CI_type='bca'
    } 
    return(list(CI_lo,CI_hi,CI_type))
  }
  tt<-cbind(MDF[,PC_colnames,with=FALSE],FDF[,PC_colnames,with=FALSE])
  overall_prop_df<-median_bias_func(tt,(1:nrow(tt)))
  outs<-evaluate_boot_CI(overall_prop_df,tt)
  CI_df_lo<-outs[1]
  CI_df_hi<-outs[2]
  CI_type_df<-outs[3]
  tt<-cbind(MDH[,PC_colnames,with=FALSE],FDH[,PC_colnames,with=FALSE])
  overall_prop_dh<-median_bias_func(tt,(1:nrow(tt)))
  outs<-evaluate_boot_CI(overall_prop_dh,tt)
  CI_dh_lo<-outs[1]
  CI_dh_hi<-outs[2]
  CI_type_dh<-outs[3]
  tt<-cbind(MVF[,PC_colnames,with=FALSE],FVF[,PC_colnames,with=FALSE])
  overall_prop_vf<-median_bias_func(tt,(1:nrow(tt)))
  outs<-evaluate_boot_CI(overall_prop_vf,tt)
  CI_vf_lo<-outs[1]
  CI_vf_hi<-outs[2]
  CI_type_vf<-outs[3]
  tt<-cbind(MVH[,PC_colnames,with=FALSE],FVH[,PC_colnames,with=FALSE])
  overall_prop_vh<-median_bias_func(tt,(1:nrow(tt)))
  outs<-evaluate_boot_CI(overall_prop_vh,tt)
  CI_vh_lo<-outs[1]
  CI_vh_hi<-outs[2]
  CI_type_vh<-outs[3]
  overall_prop<-c(overall_prop_df,overall_prop_dh,overall_prop_vf,overall_prop_vh)
  CI_lo<-c(CI_df_lo,CI_dh_lo,CI_vf_lo,CI_vh_lo)
  CI_hi<-c(CI_df_hi,CI_dh_hi,CI_vf_hi,CI_vh_hi)
  CI_type<-c(CI_type_df,CI_type_dh,CI_type_vf,CI_type_vh)
  Dataset<-c('DFW','DHW','VFW','VHW')
  dat<-data.table('Taxon'=Taxon_group,'Categorical_variable'=group,'Dataset'=Dataset,'median'=overall_prop,'CI_lo'=CI_lo,'CI_hi'=CI_hi,'CI_type'=CI_type)
  names(dat)[names(dat) == 'median'] <- paste(variable,'med',sep='_')
  names(dat)[names(dat) == 'CI_hi'] <- paste(variable,'CI_hi',sep='_')
  names(dat)[names(dat) == 'CI_lo'] <- paste(variable,'CI_lo',sep='_')
  names(dat)[names(dat) == 'CI_type'] <- paste(variable,'CI_type',sep='_')
  dat_out<-rbind(dat_out,dat)
  return(dat_out)
}
# t1<-Sys.time()
# MF_dist_Isvar_prop<-data.table()
# for(i in unique(PC$Categorical_var)){
#   boot_out<-median_dist_boot(PC_IS_var3,'MF_dist_Isvar_prop',i,n_cores,boot_reps)
#   MF_dist_Isvar_prop<-rbind(MF_dist_Isvar_prop,boot_out)
# }
# t2<-Sys.time()
# t2-t1


#####Now the MF distance in absolute rate
t1<-Sys.time()
MF_dist_abs_rate<-data.table()
for(i in unique(PC$Categorical_var)){
  boot_out<-median_dist_boot(rates1,'MF_dist_abs_rate',i,n_cores,boot_reps)
  MF_dist_abs_rate<-rbind(MF_dist_abs_rate,boot_out)
}
print('Finished MF dist absolute rate values')
t2<-Sys.time()
t2-t1

######Now the MF distance in morphological disparity

MF_MD_dist_boot<-function(indata,variable,group,cores,nboot){
  dat_out<-data.table()
  indata1<-setDT(indata)
  MDF<-subset(indata1,Wing=='FW'& Side=='dorsal' & Sex=='male' & Categorical_var ==group)
  MVF<-subset(indata1,Wing=='FW'& Side=='ventral' & Sex=='male' & Categorical_var ==group)
  MDH<-subset(indata1,Wing=='HW'& Side=='dorsal' & Sex=='male' & Categorical_var ==group)
  MVH<-subset(indata1,Wing=='HW'& Side=='ventral' & Sex=='male' & Categorical_var ==group)
  FDF<-subset(indata1,Wing=='FW'& Side=='dorsal' & Sex=='female' & Categorical_var ==group)
  FVF<-subset(indata1,Wing=='FW'& Side=='ventral' & Sex=='female' & Categorical_var ==group)
  FDH<-subset(indata1,Wing=='HW'& Side=='dorsal' & Sex=='female' & Categorical_var ==group)
  FVH<-subset(indata1,Wing=='HW'& Side=='ventral' & Sex=='female' & Categorical_var ==group)
  disp_func<-function(data,IDvec){
    dats<-data[IDvec,]
    dats_m<-dats[,1:length(PC_colnames)]
    dats_f<-dats[,(1+length(PC_colnames)):(length(PC_colnames)+length(PC_colnames))]
    trial<-rbind(dats_m,dats_f)
    lenm<-as.numeric(nrow(dats_m))
    lenf<-as.numeric(nrow(dats_f))
    try(disp<-dispRity(boot.matrix(custom.subsets(setDF(trial), list("Male"=c(1:lenm),"Female"=c((lenm+1):(lenm+lenf)))),
                                   bootstraps = 1),
                       metric = c(variances)))
    male_results<-disp$disparity$Male$elements
    female_results<-disp$disparity$Female$elements
    M_minus_F_dist<-sum(abs(male_results-female_results))
    return(M_minus_F_dist)
  }
  evaluate_boot_CI<-function(overall_value,indata1){
    if(is.na(overall_value)==TRUE){
      CI_lo<-NA
      CI_hi<-NA
      CI_type=NA
    }else if(overall_value==0){
      boot_out<-boot(indata1,disp_func,R=nboot,parallel='multicore',ncpus=cores)
      CI<-boot.ci(boot_out,conf=Bonferroni_conf_level,type='perc')
      CI_lo<-CI$percent[4]
      CI_hi<-CI$percent[5]
      CI_type='percentile'
    }else{
      boot_out<-boot(indata1,disp_func,R=nboot,parallel='multicore',ncpus=cores)
      CI<-boot.ci(boot_out,conf=Bonferroni_conf_level,type='bca')
      CI_lo<-CI$bca[4]
      CI_hi<-CI$bca[5]
      CI_type='bca'
    } 
    return(list(CI_lo,CI_hi,CI_type))
  }
  tt<-cbind(MDF[,PC_colnames,with=FALSE],FDF[,PC_colnames,with=FALSE])
  overall_prop_df<-disp_func(tt,(1:nrow(tt)))
  outs<-evaluate_boot_CI(overall_prop_df,tt)
  CI_df_lo<-outs[1]
  CI_df_hi<-outs[2]
  CI_type_df<-outs[3]
  tt<-cbind(MDH[,PC_colnames,with=FALSE],FDH[,PC_colnames,with=FALSE])
  overall_prop_dh<-disp_func(tt,(1:nrow(tt)))
  outs<-evaluate_boot_CI(overall_prop_dh,tt)
  CI_dh_lo<-outs[1]
  CI_dh_hi<-outs[2]
  CI_type_dh<-outs[3]
  tt<-cbind(MVF[,PC_colnames,with=FALSE],FVF[,PC_colnames,with=FALSE])
  overall_prop_vf<-disp_func(tt,(1:nrow(tt)))
  outs<-evaluate_boot_CI(overall_prop_vf,tt)
  CI_vf_lo<-outs[1]
  CI_vf_hi<-outs[2]
  CI_type_vf<-outs[3]
  tt<-cbind(MVH[,PC_colnames,with=FALSE],FVH[,PC_colnames,with=FALSE])
  overall_prop_vh<-disp_func(tt,(1:nrow(tt)))
  outs<-evaluate_boot_CI(overall_prop_vh,tt)
  CI_vh_lo<-outs[1]
  CI_vh_hi<-outs[2]
  CI_type_vh<-outs[3]
  overall_prop<-c(overall_prop_df,overall_prop_dh,overall_prop_vf,overall_prop_vh)
  CI_lo<-c(CI_df_lo,CI_dh_lo,CI_vf_lo,CI_vh_lo)
  CI_hi<-c(CI_df_hi,CI_dh_hi,CI_vf_hi,CI_vh_hi)
  CI_type<-c(CI_type_df,CI_type_dh,CI_type_vf,CI_type_vh)
  Dataset<-c('DFW','DHW','VFW','VHW')
  dat<-data.table('Taxon'=Taxon_group,'Categorical_variable'=group,'Dataset'=Dataset,'med'=overall_prop,'CI_lo'=CI_lo,'CI_hi'=CI_hi,'CI_type'=CI_type)
  names(dat)[names(dat) == 'med'] <- paste(variable,'med',sep='_')
  names(dat)[names(dat) == 'CI_hi'] <- paste(variable,'CI_hi',sep='_')
  names(dat)[names(dat) == 'CI_lo'] <- paste(variable,'CI_lo',sep='_')
  names(dat)[names(dat) == 'CI_type'] <- paste(variable,'CI_type',sep='_')
  dat_out<-rbind(dat_out,dat)
  return(dat_out)
}

t1<-Sys.time()
MF_dist_Extant_var_MD<-data.table()
for(i in unique(PC$Categorical_var)){
  boot_out<-MF_MD_dist_boot(PC2,'MF_dist_Extant_var_MD',i,n_cores,boot_reps)
  MF_dist_Extant_var_MD<-rbind(MF_dist_Extant_var_MD,boot_out)
}
print('Finished MD distance')
t2<-Sys.time()
t2-t1

#####Now making the ISvar SSQ dimorphism frame to use for the ISvar subset of the PC metrics----
# PC_sub<-cbind(count,PC_IS_med)
# PC_sub<-PC_sub[PC_sub$count>(IS_var_thresh-1),]
# PC_sub<-setDT(PC_sub)

# PC_sub=PC[Phy_Side_Wing %in% PC_IS_var3_F$Phy_Side_Wing]
# PC_sub=PC_sub[Phy_Side_Wing %in% PC_IS_var3_M$Phy_Side_Wing]

# 
# ISvar_ssq_boot<-function(indata,group,cores,nboot){
#   dat_out<-data.table()
#   indata1<-setDT(indata)
#   MDF<-subset(indata1,Wing=='FW'& Side=='dorsal' & Sex=='male' & Categorical_var ==group)
#   MVF<-subset(indata1,Wing=='FW'& Side=='ventral' & Sex=='male' & Categorical_var ==group)
#   MDH<-subset(indata1,Wing=='HW'& Side=='dorsal' & Sex=='male' & Categorical_var ==group)
#   MVH<-subset(indata1,Wing=='HW'& Side=='ventral' & Sex=='male' & Categorical_var ==group)
#   FDF<-subset(indata1,Wing=='FW'& Side=='dorsal' & Sex=='female' & Categorical_var ==group)
#   FVF<-subset(indata1,Wing=='FW'& Side=='ventral' & Sex=='female' & Categorical_var ==group)
#   FDH<-subset(indata1,Wing=='HW'& Side=='dorsal' & Sex=='female' & Categorical_var ==group)
#   FVH<-subset(indata1,Wing=='HW'& Side=='ventral' & Sex=='female' & Categorical_var ==group)
#   SSQ_func<-function(data,IDvec){
#     dats<-data[IDvec,]
#     ss_T<-vector()
#     ss_tot<-vector()
#     n<-vector()
#     sec_moment_distance<-vector()
#     dats_m<-dats[,1:length(PC_colnames)]
#     dats_f<-dats[,(1+length(PC_colnames)):(length(PC_colnames)+length(PC_colnames))]
#     for(i in 1:length(PC_colnames)){
#       tm<-dats_m[[i]]
#       tf<-dats_f[[i]]
#       trial<-c(tm,tf)
#       MF<-abs(tm-tf)
#       grand_mean<-mean(trial)
#       mean_m<-mean(tm)
#       mean_f<-mean(tf)
#       ss_tm<-sum(((mean_m-grand_mean)^2)*length(tm))
#       ss_tf<-sum(((mean_f-grand_mean)^2)*length(tf))
#       ss_T<-append(ss_T,ss_tm+ss_tf)
#       ss_tot<-append(ss_tot,sum((trial-grand_mean)^2))
#       sec_moment_distance<-append(sec_moment_distance,(sum((MF-0)^2))/length(MF))
#       n<-append(n,as.numeric(length(trial)))
#     }
#     SSQ_tot<-sum(ss_tot)
#     ssT_to_tot<-sum(ss_T)/SSQ_tot
#     var_tot<-SSQ_tot/mean(n)
#     sec_moment_distance<-sum(sec_moment_distance)
#     return(c(var_tot,SSQ_tot))
#   }
#   evaluate_boot_CI<-function(overall_value,indata1){
#     dt_out<-data.table()
#     for(i in 1:length(overall_value)){
#       if(is.na(overall_value[i])==TRUE){
#         CI_lo<-NA
#         CI_hi<-NA
#         CI_type=NA
#       }else if(overall_value[i]==0){
#         boot_out<-boot(indata1,SSQ_func,R=nboot,parallel='multicore',ncpus=cores)
#         CI<-boot.ci(boot_out,conf=Bonferroni_conf_level,type='perc',index=i)
#         CI_lo<-CI$percent[4]
#         CI_hi<-CI$percent[5]
#         CI_type='percentile'
#       }else{
#         boot_out<-boot(indata1,SSQ_func,R=nboot,parallel='multicore',ncpus=cores)
#         CI<-boot.ci(boot_out,conf=Bonferroni_conf_level,type='bca',index=i)
#         CI_lo<-CI$bca[4]
#         CI_hi<-CI$bca[5]
#         CI_type='bca'
#       } 
#       dt<-data.table('CI_lo'=CI_lo,'CI_hi'=CI_hi,'CI_type'=CI_type)
#       dt_out<-rbind(dt_out,dt)
#     }
#     dt_out<-cbind(c('var_tot','SSQ_tot'),dt_out)
#     colnames(dt_out)<-c('metric','CI_lo','CI_hi','CI_type')
#     return(dt_out)
#   }
#   tt<-cbind(MDF[,PC_colnames,with=FALSE],FDF[,PC_colnames,with=FALSE])
#   overall_diff_df<-SSQ_func(tt,(1:nrow(tt)))
#   outs<-evaluate_boot_CI(overall_diff_df,tt)
#   CI_df_lo<-outs$CI_lo
#   CI_df_hi<-outs$CI_hi
#   CI_type_df<-outs$CI_type
#   tt<-cbind(MDH[,PC_colnames,with=FALSE],FDH[,PC_colnames,with=FALSE])
#   overall_diff_dh<-SSQ_func(tt,(1:nrow(tt)))
#   outs<-evaluate_boot_CI(overall_diff_dh,tt)
#   CI_dh_lo<-outs$CI_lo
#   CI_dh_hi<-outs$CI_hi
#   CI_type_dh<-outs$CI_type
#   tt<-cbind(MVF[,PC_colnames,with=FALSE],FVF[,PC_colnames,with=FALSE])
#   overall_diff_vf<-SSQ_func(tt,(1:nrow(tt)))
#   outs<-evaluate_boot_CI(overall_diff_vf,tt)
#   CI_vf_lo<-outs$CI_lo
#   CI_vf_hi<-outs$CI_hi
#   CI_type_vf<-outs$CI_type
#   tt<-cbind(MVH[,PC_colnames,with=FALSE],FVH[,PC_colnames,with=FALSE])
#   overall_diff_vh<-SSQ_func(tt,(1:nrow(tt)))
#   outs<-evaluate_boot_CI(overall_diff_vh,tt)
#   CI_vh_lo<-outs$CI_lo
#   CI_vh_hi<-outs$CI_hi
#   CI_type_vh<-outs$CI_type
#   final_out<-data.table('overall_diff'=c(overall_diff_df,overall_diff_dh,overall_diff_vf,overall_diff_vh),
#                         'CI_lo'=c(CI_df_lo,CI_dh_lo,CI_vf_lo,CI_vh_lo), 'CI_hi'=c(CI_df_hi,CI_dh_hi,CI_vf_hi,CI_vh_hi),
#                         'CI_type'=c(CI_type_df,CI_type_dh,CI_type_vf,CI_type_vh))
#   final_out$Dataset<-c(rep('DFW',2),rep('DHW',2),rep('VFW',2),rep('VHW',2))
#   final_out$value<-rep(c('var_tot','SSQ_tot'),4)
#   dat<-data.table('Taxon'=Taxon_group,'Categorical_variable'=group,'Dataset'=c('DFW','DHW','VFW','VHW'),
#                   'ISvar_tot_var_med'=final_out[value=='var_tot']$overall_diff,'ISvar_tot_var_CI_lo'=final_out[value=='var_tot']$CI_lo,'ISvar_tot_var_CI_hi'=final_out[value=='var_tot']$CI_hi,'ISvar_tot_var_CI_type'=final_out[value=='var_tot']$CI_type,
#                   'ISvar_SSQ_tot_med'=final_out[value=='SSQ_tot']$overall_diff,'ISvar_SSQ_tot_CI_lo'=final_out[value=='SSQ_tot']$CI_lo,'ISvar_SSQ_tot_CI_hi'=final_out[value=='SSQ_tot']$CI_hi,'ISvar_SSQ_tot_CI_type'=final_out[value=='SSQ_tot']$CI_type)
#   dat_out<-rbind(dat_out,dat)
#   return(dat_out)
# }
# 
# 
# # SSQ_frame<-dimorph_ssq_boot(PC3,categorical_variable,n_cores,metric,boot_reps)
# t1<-Sys.time()
# IS_var_SSQ_frame<-data.table()
# for(i in unique(PC$Categorical_var)){
#   boot_out<-ISvar_ssq_boot(PC_sub,i,n_cores,boot_reps)
#   IS_var_SSQ_frame<-rbind(IS_var_SSQ_frame,boot_out)
# }
# print('Finished IS var SSQ values')
# t2<-Sys.time()
# t2-t1

final_frame<-join_all(list(SSQ_frame,rate_SSQ_frame,Abs_rate_prop_sexbias,MF_abs_rate_sexbias,MF_dist_abs_rate,MF_extant_var_MD_prop_sexbias,Extant_var_MD_prop_sexbias,MF_dist_Extant_var_MD))

# fwrite(final_frame,'Butterflies_PC_metrics_Bootpkg_redo_rate_dimorphismSSQ_MD_ISvar_100000BS_CI_pvals.csv')
boot_reps2=format(boot_reps, scientific = FALSE)

# fwrite(final_frame,paste(prefix_str,categorical_variable,'level_PC_metrics_Bootpkg_redo_rate_dimorphismSSQ_MD_ISvar_100000BS_CI_pvals',specimen_threshold,'M',specimen_threshold,'F.csv',sep='_'))
fwrite(final_frame,paste('/PATH/TO/OUTPUTFILES/noHemaris',rep,prefix_str,categorical_variable,'level_PC_metrics_Bootpkg_redo_rate_dimorphismSSQ_MD',boot_reps2,'BS_CI_pvals',specimen_threshold,'M',specimen_threshold,'F',wband,'.csv',sep='_'))

# Final_output_plotframe<-rbind(Final_output_plotframe,Final_plotframe4)
# Final_output_dimorphframe<-rbind(Final_output_dimorphframe,Final_dimorph1)
# print(paste('Done with',boot))


end_time<-Sys.time()

end_time-start_time
# write.csv(Final_output_plotframe,paste('/PATH/TO/OUTPUTFILES/boot',boot1,'to',boot2,prefix_str,'PerspecimenPC',categorical_variable,'level_rateSSQ_extant_dtt_var_IS_var_by_bias_R_ready',specimen_threshold,'males_and_females.csv',sep='_'))
# write.csv(Final_output_dimorphframe,paste('/PATH/TO/OUTPUTFILES/boot',boot1,'to',boot2,prefix_str,'PerspecimenPC',categorical_variable,'level_dimorphSSQ_R_ready',specimen_threshold,'males_and_females.csv',sep='_'))
# 


