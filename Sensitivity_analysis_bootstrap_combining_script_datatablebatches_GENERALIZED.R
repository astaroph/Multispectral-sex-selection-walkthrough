###The following is a helper script to take the results of the bootstrap parallel output and combine it into a combined output
##as well as a processed results output. 

##Add your unique string match below that matches the prefix_string used in the automated bootstrap script
# match_pattern<-'Butterflies_Final_BS'
match_pattern<-'Geometrids_Final_BS'


basename.matches <- list.files(path="/path/to/individualbootstrapoutput",
                               pattern=paste("*",match_pattern,"*",sep=""), recursive=TRUE,full.names = TRUE)

library(plyr)
library(dplyr)
library(data.table)


##Because the combining process slows down dramatically the longer each concatenated dataframe gets, I have split the process into
##five initial sub dataframes of 2000 bootstraps each, which are then combined with each other later. This dramatically speeds up the process
for(i in 1:2000){
  # print(i)
  if(i==1){
    outframe1<-fread(basename.matches[i])
  }else{
    trial<-fread(basename.matches[i])
    trial2<-trial[,5:length(colnames(trial))]
    outframe1<-cbind(outframe1,trial2)
  }
}
for(i in 2001:4000){
  # print(i)
  if(i==2001){
    outframe2<-fread(basename.matches[i])
  }else{
    trial<-fread(basename.matches[i])
    trial2<-trial[,5:length(colnames(trial))]
    outframe2<-cbind(outframe2,trial2)
  }
}
for(i in 4001:6000){
  # print(i)
  if(i==4001){
    outframe3<-fread(basename.matches[i])
  }else{
    trial<-fread(basename.matches[i])
    trial2<-trial[,5:length(colnames(trial))]
    outframe3<-cbind(outframe3,trial2)
  }
}
for(i in 6001:8000){
  # print(i)
  if(i==6001){
    outframe4<-fread(basename.matches[i])
  }else{
    trial<-fread(basename.matches[i])
    trial2<-trial[,5:length(colnames(trial))]
    outframe4<-cbind(outframe4,trial2)
  }
}
for(i in 8001:10000){
  # print(i)
  if(i==8001){
    outframe5<-fread(basename.matches[i])
  }else{
    trial<-fread(basename.matches[i])
    trial2<-trial[,5:length(colnames(trial))]
    outframe5<-cbind(outframe5,trial2)
  }
}

outframe6<-cbind(outframe1,outframe2[,5:2004],outframe3[,5:2004],outframe4[,5:2004],outframe5[,5:2004])

UID<-paste(outframe6$Value,outframe6$Dataset,outframe6$metrics,sep='_')
outframe7<-cbind(UID,outframe6)
rownames(outframe7)<-outframe7$UID
outframe8<-subset(outframe7,select=-c(V1,Value,Dataset,metrics))
trial<-data.frame(t(outframe8[,-1]))

quants <- c(0.025,.975)
quant_out<-apply(trial , 2 , quantile , probs = quants , na.rm = TRUE )
quant_range<-quant_out[2,]-quant_out[1,]
means<-colMeans(trial)
perc_error_range<-(quant_range/means)*100

outdata<-data.table(cbind(outframe7$Value,outframe7$Dataset,outframe7$metrics,means,quant_range,perc_error_range))
colnames(outdata)<-c('Value','Dataset','Metrics','Mean_boot_value','95_percent_CI','95_percent_CI_percentofmean')
size=length(trial[,1])

fwrite(outdata,paste(match_pattern,"sensitivity_analysis_results_error_range_mean",size,"calculated.csv",sep='_'))
